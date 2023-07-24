# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
#
#  Biomedical Sequencing Facility (BSF), part of the genomics core facility
#  of the Research Center for Molecular Medicine (CeMM) of the
#  Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
#  This file is part of BSF Python.
#
#  BSF Python is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BSF Python is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.
#
"""The :py:mod:`bsf.analyses.illumina_to_bam_tools` provides classes supporting the
`Illumina2Bam-Tools <https://github.com/wtsi-npg/illumina2bam>`_ package.
"""
import logging
import os
import warnings
from argparse import ArgumentParser
from subprocess import Popen
from typing import Callable, Optional

from Bio.Seq import reverse_complement

from bsf.analyses.illumina_run_folder import IlluminaRunFolderRestore
from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.argument import Argument
from bsf.connector import Connector
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.ngs import Collection, Sample, SampleAnnotationSheet
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, RunnableStep, RunnableStepChangeMode, RunnableStepJava, \
    RunnableStepMakeDirectory, RunnableStepLink, RunnableStepMove, RunnableStepPicard
from bsf.standards import Configuration, StandardFilePath, JavaArchive, JavaClassPath, Operator, VendorQualityFilter


class BamIndexDecoderSheet(AnnotationSheet):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoderSheet` class represents a
    :emphasis:`Tab-Separated Value` (TSV) table of library information for the
    :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` class.
    """

    _file_type = 'excel-tab'

    # The field names are defined in the IndexDecoder.java source file.
    # https://github.com/wtsi-npg/illumina2bam/blob/devel/src/uk/ac/sanger/npg/picard/IndexDecoder.java

    _field_name_list = [
        'barcode_sequence',
        'barcode_name',
        'library_name',
        'sample_name',
        'description',
    ]


class LibraryAnnotationSheet(AnnotationSheet):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet` class represents a
    :emphasis:`Comma-Separated Value` (CSV) table for de-multiplexing.
    """

    _file_type = 'excel'

    _field_name_list = [
        'lane',  # Lane number
        'barcode_mismatches',  # Number of mismatches allowed in index decoding
        'barcode_sequence_1',  # Index read sequence 1
        'barcode_sequence_2',  # Index read sequence 2
        'barcode_start',  # Start position of the barcode in the barcode read (numeric)
        'sample_name',  # Sample name (alphanumeric including underscore characters)
        'library_name',  # Library name (alphanumeric including underscore characters)
        'library_size',  # Library size (numeric)
        'read_structure',  # Picard tools Read Structure
    ]

    _check_dict: dict[str, list[Callable[[list[str], int, dict[str, str], str], None]]] = {
        'lane': [
            AnnotationSheet.check_alphanumeric_value,
        ],
        'barcode_mismatches': [
            AnnotationSheet.check_numeric,
        ],
        'barcode_sequence_1': [
            AnnotationSheet.check_ambiguous_sequence,
        ],
        'barcode_sequence_2': [
            AnnotationSheet.check_ambiguous_sequence,
        ],
        'barcode_start': [
            AnnotationSheet.check_numeric,
        ],
        'sample_name': [
            AnnotationSheet.check_alphanumeric_value,
            AnnotationSheet.check_underscore_leading,
            AnnotationSheet.check_underscore_trailing,
            AnnotationSheet.check_underscore_multiple,
        ],
        'library_name': [
            AnnotationSheet.check_alphanumeric_value,
            AnnotationSheet.check_underscore_leading,
            AnnotationSheet.check_underscore_trailing,
            AnnotationSheet.check_underscore_multiple,
        ],
        'read_structure': [
            AnnotationSheet.check_alphanumeric,
        ],
    }

    def get_annotation_dict(self) -> dict[int, list[dict[str, str]]]:
        """Get a Python :py:class:`dict` object of Python :py:class:`int` (lane) key and
        Python :py:class`list` objects of Python :py:class:`dict` (row dict) value objects.

        :return: A Python :py:class:`dict` object of Python :py:class:`int` (lane) key and
            Python :py:class:`list` objects of Python :py:class:`dict` (row dict) value objects.
        :rtype: dict[int, list[dict[str, str]]]
        """
        flow_cell_dict: dict[int, list[dict[str, str]]] = dict()

        for row_dict in self.row_dict_list:
            lane_int = int(row_dict['lane'])
            if lane_int not in flow_cell_dict:
                flow_cell_dict[lane_int] = list()
            flow_cell_dict[lane_int].append(row_dict)

        return flow_cell_dict

    def get_barcode_length_dict(self) -> dict[int, tuple[int, int]]:
        """Get a Python :py:class:`dict` object of :literal:`i7` and :literal:`i5` barcode sequence lengths by lane.

        :return: A Python :py:class:`dict` object of Python :py:class:`int` (lane) key and
            Python :py:class:`tuple` objects of Python :py:class:`int` (:emphasis:`i7`) and
            Python :py:class:`int` (:emphasis:`i5`) barcode sequence length value objects.
        :rtype: dict[int, (int, int)]
        """
        flow_cell_dict: dict[int, tuple[int, int]] = dict()

        for row_dict in self.row_dict_list:
            lane_int = int(row_dict['lane'])
            if lane_int not in flow_cell_dict:
                flow_cell_dict[lane_int] = (len(row_dict['barcode_sequence_1']), len(row_dict['barcode_sequence_2']))
                continue
            for index in range(0, 1 + 1):
                if flow_cell_dict[lane_int][index] != len(row_dict['barcode_sequence_' + str(index + 1)]):
                    warnings.warn(f"The length of 'barcode_sequence_{index + 1!s}' does not match previous records.\n"
                                  f"{row_dict!r}", UserWarning)

        return flow_cell_dict

    def validate(self, lanes: Optional[int] = None) -> list[str]:
        """
        Validate a :py:class:`bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet` object.

        If the number of lanes is provided, all lanes in the range are validated or else,
        the validation is just based on annotation for the lane variable.

        :param lanes: A number of lanes to validate.
        :type lanes: int | None
        :return: A Python :py:class:`list` object of Python :py:class:`str` (warning message) objects.
        :rtype: list[str]
        """
        message_list: list[str] = []

        # Check the header line via the pre-defined field names.

        for index in range(0, len(self._field_name_list)):
            if not self.field_name_list[index]:
                message_list.append(
                    f'Column with name {self._field_name_list[index]!r} is missing from the header line.')

            if self.field_name_list[index] != self._field_name_list[index]:
                message_list.append(
                    f'Column name {self.field_name_list[index]!r} in the header line '
                    f'does not match template {self._field_name_list[index]!r}.')

        # Validate the field values for alphanumeric or sequence grade in the context of the
        # AnnotationSheet super-class.

        message_list.extend(super(LibraryAnnotationSheet, self).validate())

        flow_cell_dict = {}

        for row_index, row_dict in enumerate(self.row_dict_list):
            row_number = row_index + 1
            lane_int = int(row_dict['lane'])

            # Check that all required fields are defined.

            if lane_int not in flow_cell_dict:
                flow_cell_dict[lane_int] = {
                    'barcode_dict': {},  # Barcode sequence (str) to row number (int).
                    'sample_dict': {},  # Sample name (str) to row number (int).
                    'library_name': row_dict['library_name'],
                }

                if 'barcode_mismatches' in row_dict:
                    flow_cell_dict[lane_int]['barcode_mismatches'] = row_dict['barcode_mismatches']

                if 'barcode_start' in row_dict:
                    flow_cell_dict[lane_int]['barcode_start'] = row_dict['barcode_start']

                if 'read_structure' in row_dict:
                    flow_cell_dict[lane_int]['read_structure'] = row_dict['read_structure']

            if 'barcode_sequence_1' in row_dict and row_dict['barcode_sequence_1']:
                barcode_sequence_1 = row_dict['barcode_sequence_1']
            else:
                barcode_sequence_1 = '-NoIndex-'

            if 'barcode_sequence_2' in row_dict and row_dict['barcode_sequence_2']:
                barcode_sequence_2 = row_dict['barcode_sequence_2']
            else:
                barcode_sequence_2 = '-NoIndex-'

            barcode_sequence = barcode_sequence_1 + barcode_sequence_2

            if barcode_sequence in flow_cell_dict[lane_int]['barcode_dict']:
                message_list.append(
                    f"Barcode sequence {barcode_sequence!r} from "
                    f"row {flow_cell_dict[lane_int]['barcode_dict'][barcode_sequence]!s} "
                    f"duplicated in row {row_number!s}.")
            else:
                flow_cell_dict[lane_int]['barcode_dict'][barcode_sequence] = row_number

            # Check for duplicate sample name values.
            if row_dict['sample_name'] in flow_cell_dict[lane_int]['sample_dict']:
                message_list.append(
                    f"Sample name {row_dict['sample_name']!s} from "
                    f"row {flow_cell_dict[lane_int]['sample_dict'][row_dict['sample_name']]!s} "
                    f"duplicated in row {row_number!s}.")
            else:
                flow_cell_dict[lane_int]['sample_dict'][row_dict['sample_name']] = row_number

            # Check for identical library name values.
            if flow_cell_dict[lane_int]['library_name'] != row_dict['library_name']:
                message_list.append(
                    f"Library name {row_dict['library_name']!r} in row {row_number!s} "
                    f"does not match previous name {flow_cell_dict[lane_int]['library_name']!r}.")

            # Check for identical barcode mismatches values.
            if 'barcode_mismatches' in row_dict:
                if flow_cell_dict[lane_int]['barcode_mismatches'] != row_dict['barcode_mismatches']:
                    message_list.append(
                        f"Barcode mismatches {row_dict['barcode_mismatches']!r} in row {row_number!s} "
                        f"does not match previous mismatches {flow_cell_dict[lane_int]['barcode_mismatches']!r}.")

            # Check for identical barcode start values.
            if 'barcode_start' in row_dict:
                if flow_cell_dict[lane_int]['barcode_start'] != row_dict['barcode_start']:
                    message_list.append(
                        f"Barcode start {row_dict['barcode_start']!r} in row {row_number!s} "
                        f"does not match previous start {flow_cell_dict[lane_int]['barcode_start']!r}.")

            # Check for identical read structure values.
            if 'read_structure' in row_dict:
                if flow_cell_dict[lane_int]['read_structure'] != row_dict['read_structure']:
                    message_list.append(
                        f"Read structure {row_dict['read_structure']!r} in row {row_number!s} "
                        f"does not match previous read structure {flow_cell_dict[lane_int]['read_structure']!r}.")

        if lanes is None:
            # If the number of lanes was not provided, validate just on the basis of the lane annotation.
            lane_list = sorted(flow_cell_dict)
        else:
            # If the number of lanes was provided, validate all lanes in the range.
            lane_list = [lane_int for lane_int in range(0 + 1, lanes + 1)]

        for lane_int in lane_list:
            # Check that all lanes have annotation.
            if lane_int not in flow_cell_dict:
                message_list.append(f'No annotation for lane number {lane_int!s}.')
                continue

            # Check that all or none of the rows have index sequence 1 or 2 populated.
            no_index_1 = 0
            no_index_2 = 0
            for key in flow_cell_dict[lane_int]['barcode_dict']:
                if key[:9] == '-NoIndex-':
                    no_index_1 += 1
                if key[-9:] == '-NoIndex-':
                    no_index_2 += 1

            if not (no_index_1 == 0 or no_index_1 == len(flow_cell_dict[lane_int]['barcode_dict'])):
                message_list.append(f'Some empty barcode_sequence_1 fields in lane {lane_int!s}.')
            if not (no_index_2 == 0 or no_index_2 == len(flow_cell_dict[lane_int]['barcode_dict'])):
                message_list.append(f'Some empty barcode_sequence_2 fields in lane {lane_int!s}.')

            # Check that all barcode sequences have the same length.
            # This test also finds cases of missing sequences tested for above.
            key_length: Optional[int] = None
            for key in flow_cell_dict[lane_int]['barcode_dict']:
                if key_length is None:
                    key_length = len(key)
                if len(key) != key_length:
                    message_list.append(f'Mismatching barcode sequence lengths in lane {lane_int!s}.')

        return message_list

    @classmethod
    def console_validate(
            cls,
            library_path: str,
            lanes: Optional[int] = None) -> int:
        """Console function to validate a :py:class:`bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet`.

        :param library_path: A library annotation sheet CSV file path.
        :type library_path: str
        :param lanes: A number of lanes to validate.
        :type lanes: int | None
        :return: A :py:class:`SystemExit` status value. Returns 0 upon successful validation, 1 otherwise.
        :rtype: int
        """
        library_annotation_sheet: LibraryAnnotationSheet = LibraryAnnotationSheet.from_file_path(
            file_path=library_path)

        message_list = library_annotation_sheet.validate(lanes=lanes)

        return_value = 0

        if message_list:
            # warnings.warn(f'Validation of library annotation sheet {library_path!r}:\n'
            #               '\n'.join(message_list), UserWarning)
            print(*message_list, sep='\n')
            return_value = 1

        return return_value

    @classmethod
    def entry_point_validate(cls) -> int:
        """Console entry point to validate a :py:class:`bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet`.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description='BSF Python utility script to validate a Library Annotation Sheet.')

        argument_parser.add_argument(
            '--logging-level',
            default='WARNING',
            choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
            help='logging level [WARNING]')

        argument_parser.add_argument(
            '--lanes',
            type=int,
            help='number of lanes to validate')

        argument_parser.add_argument(
            'library_path',
            help='library annotation sheet CSV file path')

        name_space = argument_parser.parse_args()

        if name_space.logging_level:
            logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
            logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

            logging.basicConfig(level=name_space.logging_level)

        return cls.console_validate(library_path=name_space.library_path, lanes=name_space.lanes)

    @classmethod
    def console_reverse(
            cls,
            input_path: str,
            output_path: str,
            lanes: Optional[str] = None,
            barcode_1: Optional[bool] = None,
            barcode_2: Optional[bool] = None) -> int:
        """Console function to reverse-complement barcodes in a
        :py:class:`bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet`.

        :param input_path: An input library annotation sheet path.
        :type input_path: str
        :param output_path: An output library annotation sheet path.
        :type output_path: str | None
        :param lanes: A comma-separated list of lanes to reverse-complement.
        :type lanes: str | None
        :param barcode_1: Request reverse-complementing barcode 1.
        :type barcode_1: bool | None
        :param barcode_2: Request reverse-complementing barcode 2.
        :type barcode_2: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if lanes:
            lane_list: Optional[list[str]] = lanes.split(',')
        else:
            lane_list: Optional[list[str]] = None

        library_annotation_sheet: LibraryAnnotationSheet = LibraryAnnotationSheet.from_file_path(file_path=input_path)

        for row_dict in library_annotation_sheet.row_dict_list:
            if not lane_list or row_dict['lane'] in lane_list:
                if row_dict['barcode_sequence_1'] and barcode_1:
                    row_dict['barcode_sequence_1'] = reverse_complement(sequence=row_dict['barcode_sequence_1'])

                if row_dict['barcode_sequence_2'] and barcode_2:
                    row_dict['barcode_sequence_2'] = reverse_complement(sequence=row_dict['barcode_sequence_2'])

        library_annotation_sheet.file_path = output_path
        library_annotation_sheet.to_file_path(adjust_field_names=True)

        return 0

    @classmethod
    def entry_point_reverse(cls) -> int:
        """Console entry point to reverse-complement barcodes in a
        :py:class:`bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet`.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description='BSF Python utility script to reverse-complement a Library Annotation Sheet.')

        argument_parser.add_argument(
            '--input-path',
            required=True,
            help='input library annotation sheet')

        argument_parser.add_argument(
            '--output-path',
            required=True,
            help='output library annotation sheet')

        argument_parser.add_argument(
            '--lanes',
            help='comma-separated list of lanes to reverse complement')

        argument_parser.add_argument(
            '--barcode-1',
            action='store_true',
            help='reverse complement barcode 1')

        argument_parser.add_argument(
            '--barcode-2',
            action='store_true',
            help='reverse complement barcode 2')

        name_space = argument_parser.parse_args()

        return cls.console_reverse(
            input_path=name_space.input_path,
            output_path=name_space.output_path,
            lanes=name_space.lanes,
            barcode_1=name_space.barcode_1,
            barcode_2=name_space.barcode_2)


class RunnableStepIlluminaToBam(RunnableStepJava):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.RunnableStepIlluminaToBam` class represents a
    :py:class:`bsf.process.RunnableStepJava` class specific to
    `Illumina2bam <https://github.com/wtsi-npg/illumina2bam>`_ tools, which use the old
    `Picard <http://broadinstitute.github.io/picard/>`_ interface where each algorithm is implemented as a separate
    :emphasis:`Java Archive` (JAR) file.
    """

    def __init__(
            self,
            name: Optional[str],
            program: Optional[str] = None,
            options: Optional[dict[str, list[Argument]]] = None,
            arguments: Optional[list[str]] = None,
            sub_command: Optional[Command] = None,
            stdin: Optional[Connector] = None,
            stdout: Optional[Connector] = None,
            stderr: Optional[Connector] = None,
            dependencies: Optional[list[str]] = None,
            hold: Optional[bool] = None,
            submit: bool = True,
            process_identifier: Optional[str] = None,
            process_name: Optional[str] = None,
            sub_process: Optional[Popen] = None,
            obsolete_file_path_list: Optional[list[str]] = None,
            java_temporary_path: Optional[str] = None,
            java_heap_maximum: Optional[str] = None,
            java_jar_path: Optional[str] = None,
            itb_classpath: Optional[str] = None,
            itb_command: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.process.RunnableStepIlluminaToBam` object.

        :param name: A name.
        :type name: str | None
        :param program: A program.
        :type program: str | None
        :param options: A Python :py:class:`dict` object of
            Python :py:class:`str` (:py:attr:`bsf.argument.Argument.key`) key and
            Python :py:class:`list` value objects of :py:class:`bsf.argument.Argument` objects.
        :type options: dict[str, list[Argument]] | None
        :param arguments: A Python :py:class:`list` object of Python :py:class:`str` (program argument) objects.
        :type arguments: list[str] | None
        :param sub_command: A subordinate :py:class:`bsf.process.Command` object.
        :type sub_command: Command | None
        :param stdin: A standard input :literal:`STDIN` :py:class:`bsf.connector.Connector` object.
        :type stdin: Connector | None
        :param stdout: A standard output :literal:`STDOUT` :py:class:`bsf.connector.Connector` object.
        :type stdout: Connector | None
        :param stderr: A standard error :literal:`STDERR` :py:class:`bsf.connector.Connector` object.
        :type stderr: Connector | None
        :param dependencies: A Python :py:class:`list` object of
            Python :py:class:`str` (:py:attr:`bsf.process.Executable.name`) objects
            in the context of :py:class:`bsf.analysis.Stage` dependencies.
        :type dependencies: list[str] | None
        :param hold: Request a hold on job scheduling.
        :type hold: bool | None
        :param submit: Request the submission via the :py:meth:`bsf.analysis.Stage.submit` method.
        :type submit: bool
        :param process_identifier: A process identifier.
        :type process_identifier: str | None
        :param process_name: A process name.
        :type process_name: str | None
        :param sub_process: A :py:class:`subprocess.Popen` object.
        :type sub_process: Popen | None
        :param obsolete_file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (file path) objects
            that can be removed after successfully completing the :py:meth:`bsf.process.RunnableStep.run` method.
        :type obsolete_file_path_list: list[str] | None
        :param java_temporary_path: A temporary directory path for the Java Virtual Machine.
        :type java_temporary_path: str | None
        :param java_heap_maximum: A Java heap maximum size (-Xmx option).
        :type java_heap_maximum: str | None
        :param java_jar_path: A :emphasis:`Java Archive` (JAR) file path.
        :type java_jar_path: str | None
        :param itb_classpath: An `Illumina2bam <https://github.com/wtsi-npg/illumina2bam>`_ Java Class Path directory.
        :type itb_classpath: str | None
        :param itb_command: An `Illumina2bam <https://github.com/wtsi-npg/illumina2bam>`_ command.
        :type itb_command: str | None
        """
        super(RunnableStepIlluminaToBam, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path,
            java_heap_maximum=java_heap_maximum,
            java_jar_path=java_jar_path)

        # Set the Illumina2bam Java Class Path and the command-specific Illumina2bam Java Archive (JAR) file path.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(itb_classpath, itb_command + '.jar'))

        # Introduce another empty sub-command to separate Illumina2Bam-tools options.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command()

        return

    def add_itb_option(self, key: str, value: str, override: bool = False) -> None:
        """Add a :py:class:`bsf.argument.OptionPair` object to a
        :py:class:`bsf.analyses.illumina_to_bam_tools.RunnableStepIlluminaToBam` object.

        :param key: An option key.
        :type key: str
        :param value: An option value.
        :type value: str
        :param override: Override an existing :py:class:`bsf.argument.Argument` object without warning.
        :type override: bool
        """
        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class FilePathIlluminaToBamLane(FilePath):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathIlluminaToBamLane` class models file names.

    :ivar unsorted_bam: An :emphasis:`unsorted` BAM file path.
    :type unsorted_bam: str
    :ivar unsorted_md5: An :emphasis:`unsorted` BAM file MD5 checksum file path.
    :type unsorted_md5: str
    :ivar sorted_bam: A :emphasis:`sorted` BAM file path.
    :type sorted_bam: str
    :ivar sorted_md5: A :emphasis:`sorted` BAM file MD5 checksum file path.
    :type sorted_md5: str
    :ivar archive_bam: An :emphasis:`archive` BAM file path.
    :type archive_bam: str
    :ivar archive_md5: An :emphasis:`archive` BAM file MD5 checksum file path.
    :type archive_md5: str
    """

    def __init__(self, prefix: str, project_name: str, lane: str, experiment_directory: str) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathIlluminaToBamLane` object.
        
        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param experiment_directory: An experiment-specific directory.
        :type experiment_directory: str
        """
        super(FilePathIlluminaToBamLane, self).__init__(prefix=prefix)

        self.unsorted_bam = prefix + '_unsorted.bam'
        self.unsorted_md5 = prefix + '_unsorted.bam.md5'
        self.sorted_bam = prefix + '_sorted.bam'
        self.sorted_md5 = prefix + '_sorted.bam.md5'
        # The final BAM and MD5 checksum files are non-standard in that they do not contain a prefix and
        # reside in the experiment_directory.
        self.archive_bam = os.path.join(experiment_directory, '_'.join((project_name, lane)) + '.bam')
        self.archive_md5 = os.path.join(experiment_directory, '_'.join((project_name, lane)) + '.bam.md5')

        return


class IlluminaToBam(Analysis):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` class represents the logic to
    convert Illumina BCL to a BAM or SAM files.

    :ivar run_directory: An :emphasis:`Illumina Run Folder` (IRF) directory path.
    :type run_directory: str | None
    :ivar intensity_directory: An :emphasis:`Illumina Run Folder` (IRF) :literal:`Intensities` directory path,
        defaults to :literal:`IRF/Data/Intensities`.
    :type intensity_directory: str | None
    :ivar basecalls_directory: An :emphasis:`Illumina Run Folder` (IRF) :literal:`BaseCalls` directory path,
        defaults to :literal:`IRF/Data/Intensities/BaseCalls`.
    :type basecalls_directory: str | None
    :ivar experiment_name: An experiment name (i.e., flow cell identifier) normally automatically read from
        :emphasis:`Illumina Run Folder` parameters.
    :type experiment_name: str | None
    :ivar sequencing_centre: A sequencing centre code.
    :type sequencing_centre: str | None
    :ivar sort_output: Request sorting of BAM files.
    :type sort_output: bool | None
    :ivar sequences_directory: A directory storing lane-specific unaligned BAM files.
    :type sequences_directory: str | None
    :ivar mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
    :type mode_directory: str | None
    :ivar mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
    :type mode_file: str | None
    :ivar java_classpath_illumina2bam: An Illumina2Bam tools Java Class Path directory path.
    :type java_classpath_illumina2bam: str | None
    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    :ivar vendor_quality_filter: Request vendor quality filtering.
    :type vendor_quality_filter: bool
    :ivar force: Request processing of incomplete :emphasis:`Illumina Run Folder` entities.
    :type force: bool | None
    """

    name = 'Illumina To Bam Analysis'
    prefix = 'illumina_to_bam'

    @classmethod
    def get_stage_name_cell(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cell'))

    @classmethod
    def get_stage_name_lane(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'lane'))

    @classmethod
    def get_prefix_cell(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_cell(), project_name))

    @classmethod
    def get_prefix_lane(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_lane(), project_name, lane))

    @classmethod
    def get_file_path_lane(cls, project_name: str, lane: str, experiment_directory: str) -> FilePathIlluminaToBamLane:
        """Get a :py:class:`bsf.analyses.illumina_to_bam.FilePathIlluminaToBamLane` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param experiment_directory: An experiment directory.
        :type experiment_directory: str
        :return: :py:class:`bsf.analyses.illumina_to_bam.FilePathIlluminaToBamLane` object.
        :rtype: FilePathIlluminaToBamLane
        """
        return FilePathIlluminaToBamLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane,
            experiment_directory=experiment_directory)

    def __init__(
            self,
            configuration: Optional[Configuration] = None,
            project_name: Optional[str] = None,
            genome_version: Optional[str] = None,
            input_directory: Optional[str] = None,
            output_directory: Optional[str] = None,
            project_directory: Optional[str] = None,
            genome_directory: Optional[str] = None,
            report_style_path: Optional[str] = None,
            report_header_path: Optional[str] = None,
            report_footer_path: Optional[str] = None,
            e_mail: Optional[str] = None,
            stage_list: Optional[list[Stage]] = None,
            collection: Optional[Collection] = None,
            sample_list: Optional[list[Sample]] = None,
            run_directory: Optional[str] = None,
            intensity_directory: Optional[str] = None,
            basecalls_directory: Optional[str] = None,
            experiment_name: Optional[str] = None,
            sequencing_centre: Optional[str] = None,
            sort_output: Optional[bool] = None,
            sequences_directory: Optional[str] = None,
            mode_directory: Optional[str] = None,
            mode_file: Optional[str] = None,
            java_classpath_illumina2bam: Optional[str] = None,
            java_archive_picard: Optional[str] = None,
            vendor_quality_filter: Optional[bool] = None,
            force: Optional[bool] = False) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` object.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration | None
        :param project_name: A project name.
        :type project_name: str | None
        :param genome_version: A genome assembly version.
        :type genome_version: str | None
        :param input_directory: An input directory path.
        :type input_directory: str | None
        :param output_directory: An output directory path.
        :type output_directory: str | None
        :param project_directory: A project directory path, normally under the output directory path.
        :type project_directory: str | None
        :param genome_directory: A genome directory path, normally under the project directory path.
        :type genome_directory: str | None
        :param report_style_path: A report style :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: A report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: A report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a :emphasis:`UCSC Genome Browser Track Hub`.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param run_directory: An :emphasis:`Illumina Run Folder` (IRF) directory path.
        :type run_directory: str | None
        :param intensity_directory: An :emphasis:`Illumina Run Folder` (IRF) :literal:`Intensities` directory path,
            defaults to :literal:`IRF/Data/Intensities`.
        :type intensity_directory: str | None
        :param basecalls_directory: File path to the :literal:`BaseCalls` directory,
            defaults to :literal:`IRF/Data/Intensities/BaseCalls`.
        :type basecalls_directory: str | None
        :param experiment_name: An experiment name (i.e., flow cell identifier) normally automatically read from
            :emphasis:`Illumina Run Folder` parameters.
        :type experiment_name: str | None
        :param sequencing_centre: A sequencing centre code.
        :type sequencing_centre: str | None
        :param sort_output: Request sorting of BAM files.
        :type sort_output: bool | None
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str | None
        :param mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
        :type mode_directory: str | None
        :param mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
        :type mode_file: str | None
        :param java_classpath_illumina2bam: An Illumina2Bam tools Java Class Path directory path.
        :type java_classpath_illumina2bam: str | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param vendor_quality_filter: Request vendor quality filtering.
        :type vendor_quality_filter: bool | None
        :param force: Request processing of incomplete :emphasis:`Illumina Run Folder` entities.
        :type force: bool | None
        """
        super(IlluminaToBam, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            report_style_path=report_style_path,
            report_header_path=report_header_path,
            report_footer_path=report_footer_path,
            e_mail=e_mail,
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list)

        # Sub-class specific ...

        self.run_directory = run_directory
        self.intensity_directory = intensity_directory
        self.basecalls_directory = basecalls_directory
        self.experiment_name = experiment_name
        self.sequencing_centre = sequencing_centre
        self.sort_output = sort_output
        self.sequences_directory = sequences_directory
        self.mode_directory = mode_directory
        self.mode_file = mode_file
        self.java_classpath_illumina2bam = java_classpath_illumina2bam
        self.java_archive_picard = java_archive_picard
        self.vendor_quality_filter = vendor_quality_filter
        self.force = force

        return

    @property
    def get_experiment_directory(self) -> Optional[str]:
        """Get the experiment directory.

        The experiment directory is a concatenation of the sequences directory and the project name.

        :return: The experiment directory.
        :rtype: str | None
        """
        if self.sequences_directory and self.project_name:
            return os.path.join(self.sequences_directory, self.project_name)

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(IlluminaToBam, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'illumina_run_folder'
        if configuration.config_parser.has_option(section=section, option=option):
            self.run_directory = configuration.config_parser.get(section=section, option=option)

        option = 'intensity_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.intensity_directory = configuration.config_parser.get(section=section, option=option)

        option = 'basecalls_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.basecalls_directory = configuration.config_parser.get(section=section, option=option)

        option = 'experiment_name'
        if configuration.config_parser.has_option(section=section, option=option):
            self.experiment_name = configuration.config_parser.get(section=section, option=option)

        option = 'sequencing_centre'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequencing_centre = configuration.config_parser.get(section=section, option=option)

        option = 'sort_output'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sort_output = configuration.config_parser.getboolean(section=section, option=option)

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        option = 'java_classpath_illumina2bam'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_classpath_illumina2bam = configuration.config_parser.get(section=section, option=option)

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        option = 'vendor_quality_filter'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vendor_quality_filter = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` object.

        Convert an Illumina flow cell into lane-specific archive BAM files.

        To convert an Illumina flow cell, :literal:`Illumina2bam` is run first, setting the SAM Read Group (@RG)
        library name (LB) and sample name (SM) to :literal:`<flow cell identifier>.<lane>`.
        The resulting archive BAM file is then sorted by query name with Picard SortSam.
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(f"An {self.name!s} requires a 'run_directory' configuration option.")

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=StandardFilePath.get_illumina_run(absolute=True))

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(f"The Illumina Run Folder ('run_directory') {self.run_directory!r} does not exist.")

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(f'The Illumina Run Folder {self.run_directory!r} is not complete.')

        # Define an 'Intensities' directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Illumina Run Folder path.

        if self.intensity_directory:
            intensity_directory = self.configuration.get_absolute_path(
                file_path=self.intensity_directory,
                default_path=self.run_directory)
        else:
            intensity_directory = os.path.join(self.run_directory, 'Data', 'Intensities')

        # Check that the "Intensities" directory exists.

        if not os.path.isdir(intensity_directory):
            raise Exception(f"The 'intensity_directory' {intensity_directory!r} does not exist.")

        # Define a 'BaseCalls' directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        # Check if an absolute path has been provided, if not,
        # automatically prepend the "Intensities" directory path.

        if self.basecalls_directory:
            basecalls_directory = self.configuration.get_absolute_path(
                file_path=self.basecalls_directory,
                default_path=intensity_directory)
        else:
            basecalls_directory = os.path.join(intensity_directory, 'BaseCalls')

        # Check that the BaseCalls directory exists.

        if not os.path.isdir(basecalls_directory):
            raise Exception(f"The 'basecalls_directory' {basecalls_directory!r} does not exist.")

        irf = RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g., BSF_0000) is used as the prefix for archive BAM files.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = irf.run_parameters.get_experiment_name

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, irf.run_information.flow_cell))

        # Get sequencing centre information.

        if not self.sequencing_centre:
            self.sequencing_centre = Operator.get_sequencing_centre()
            if not self.sequencing_centre:
                raise Exception(f"An {self.name!s} requires a 'sequencing_centre' configuration option.")

        # Define the "Sequences" directory in which to create the experiment directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        # An absolute path cannot be prepended.

        if self.sequences_directory:
            self.sequences_directory = self.configuration.get_absolute_path(
                file_path=self.sequences_directory)
        else:
            self.sequences_directory = StandardFilePath.get_sequences(absolute=True)

        # As a safety measure, to prevent creation of rogue directory paths, the sequences_directory has to exist.

        if not os.path.isdir(self.sequences_directory):
            raise Exception(f"The 'sequences_directory' {self.sequences_directory!r} does not exist.")

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the Illumina2Bam tools Java Class Path directory.

        if not self.java_classpath_illumina2bam:
            self.java_classpath_illumina2bam = JavaClassPath.get_illumina2bam()
            if not self.java_classpath_illumina2bam:
                raise Exception(f"An {self.name!s} requires a 'java_classpath_illumina2bam' configuration option.")

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"An {self.name!s} requires a 'java_archive_picard' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self.vendor_quality_filter is None:
            self.vendor_quality_filter = VendorQualityFilter.get_vendor_quality_filter(
                flow_cell_type=irf.run_parameters.get_flow_cell_type)

        # Call the run method of the super class after the project_name has been defined.

        super(IlluminaToBam, self).run()

        lane_dependency_str: Optional[str] = None

        cell_dependency_list: list[str] = list()

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        for lane_int in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            file_path_lane = self.get_file_path_lane(
                project_name=self.project_name,
                lane=lane_str,
                experiment_directory=experiment_directory)

            # NOTE: The bsf.procedure.Runnable.name has to match the Executable.name that gets submitted via the Stage.
            runnable_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))

            # TODO: The Runnable class could have dependencies just like the Executable class so that they could be
            #  passed on upon creation of the Executable from the Runnable via Executable.from_analysis_runnable().
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)
            executable_lane.dependencies.append(
                IlluminaRunFolderRestore.get_prefix_compress_base_calls(
                    project_name=self.project_name,
                    lane=lane_str))
            # For NextSeq instruments the number of open file handles can become really large.
            # Set dependencies to run Illumina2bam lane for lane.
            if irf.run_parameters.get_instrument_type in ('NextSeq',) and lane_dependency_str:
                executable_lane.dependencies.append(lane_dependency_str)

            lane_dependency_str = runnable_lane.name

            # Add the dependency for the cell-specific process.
            cell_dependency_list.append(executable_lane.name)

            # Only submit this Executable if the final result file does not exist.
            # absolute_sorted_md5 = os.path.join(experiment_directory, file_path_lane.sorted_md5)
            # if os.path.exists(absolute_sorted_md5) and os.path.getsize(absolute_sorted_md5):
            #     executable_lane.submit = False

            # Run Illumina2Bam tools Illumina2bam.

            # TODO: It would be good to allow for configuration of RunnableStep objects via configuration files.
            if self.sort_output:
                java_heap_maximum = 'Xmx16G'
            else:
                java_heap_maximum = 'Xmx4G'

            runnable_step = RunnableStepIlluminaToBam(
                name='illumina_to_bam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum=java_heap_maximum,
                itb_classpath=self.java_classpath_illumina2bam,
                itb_command='Illumina2bam')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            if self.intensity_directory:
                # RUN_FOLDER defaults to 'null'.
                # Only set the RUN_FOLDER option, if a separate 'Intensities' directory has been configured.
                # The default is to use the directory two up from the INTENSITY_DIR.
                runnable_step.add_itb_option(key='RUN_FOLDER', value=self.run_directory)
            # INTENSITY_DIR is required.
            runnable_step.add_itb_option(key='INTENSITY_DIR', value=intensity_directory)
            if self.basecalls_directory:
                # BASECALLS_DIR defaults to 'null'.
                # Only set the BASECALLS_DIR option, if a separate 'BaseCalls' directory has been configured.
                # The default is to use the 'BaseCalls' directory under the INTENSITY_DIR.
                runnable_step.add_itb_option(key='BASECALLS_DIR', value=basecalls_directory)
            # LANE is required.
            runnable_step.add_itb_option(key='LANE', value=lane_str)
            # OUTPUT is required.
            runnable_step.add_itb_option(key='OUTPUT', value=file_path_lane.unsorted_bam)
            # GENERATE_SECONDARY_BASE_CALLS defaults to 'false'.
            # PF_FILTER defaults to 'true'.
            if self.vendor_quality_filter:
                runnable_step.add_itb_option(key='PF_FILTER', value='true')
            else:
                runnable_step.add_itb_option(key='PF_FILTER', value='false')
            # READ_GROUP_ID defaults to '1'.
            runnable_step.add_itb_option(key='READ_GROUP_ID', value='_'.join((irf.run_information.flow_cell, lane_str)))
            # SAMPLE_ALIAS defaults to 'null', using LIBRARY_NAME.
            # LIBRARY_NAME defaults to 'unknown'.
            runnable_step.add_itb_option(key='LIBRARY_NAME', value='_'.join((irf.run_information.flow_cell, lane_str)))
            # STUDY_NAME defaults to 'null'.
            # PLATFORM_UNIT defaults to 'null', using run folder name plus lane number.
            # RUN_START_DATE defaults to 'null', using the configuration file value.
            # SEQUENCING_CENTER defaults to 'SC' for Sanger Center.
            runnable_step.add_itb_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)
            # PLATFORM defaults to 'ILLUMINA'.
            # FIRST_TILE defaults to 'null'.
            # TILE_LIMIT defaults to 'null'.
            # BARCODE_SEQUENCE_TAG_NAME defaults to 'BC'.
            # BARCODE_QUALITY_TAG_NAME defaults to 'QT'.
            # BC_READ defaults to null.
            # SECOND_BARCODE_SEQUENCE_TAG_NAME defaults to 'null'.
            # SECOND_BARCODE_QUALITY_TAG_NAME defaults to 'null'.
            # SEC_BC_READ defaults to 'null'.
            # FIRST_CYCLE
            # FINAL_CYCLE
            # FIRST_INDEX_CYCLE
            # FINAL_INDEX_CYCLE
            # ADD_CLUSTER_INDEX_TAG defaults to 'false'.
            # SORT_OUTPUT defaults to 'false'.
            if self.sort_output:
                runnable_step.add_itb_option(key='SORT_OUTPUT', value='true')
            else:
                runnable_step.add_itb_option(key='SORT_OUTPUT', value='false')
            # NEXTSEQ defaults to 'false'.
            if irf.run_parameters.get_instrument_type in ('NextSeq',):
                runnable_step.add_itb_option(key='NEXTSEQ', value='true')
            else:
                runnable_step.add_itb_option(key='NEXTSEQ', value='false')
            # TMP_DIR
            runnable_step.add_itb_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # VERBOSITY defaults to 'INFO'.
            runnable_step.add_itb_option(key='VERBOSITY', value='WARNING')
            # QUIET defaults to 'false'.
            # VALIDATION_STRINGENCY defaults to 'STRICT'.
            # COMPRESSION_LEVEL defaults to '5'.
            runnable_step.add_itb_option(key='COMPRESSION_LEVEL', value='9')
            # MAX_RECORDS_IN_RAM defaults to '500000'.
            runnable_step.add_itb_option(key='MAX_RECORDS_IN_RAM', value='2000000')
            # CREATE_INDEX defaults to 'false'.
            # CREATE_MD5_FILE defaults to 'false'.
            runnable_step.add_itb_option(key='CREATE_MD5_FILE', value='true')
            # OPTIONS_FILE

            # Run Picard SortSam, if not sorted already.

            if self.sort_output:
                runnable_step = RunnableStepMove(
                    name='move_unsorted_bam',
                    source_path=file_path_lane.unsorted_bam,
                    target_path=file_path_lane.sorted_bam)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step = RunnableStepMove(
                    name='move_unsorted_md5',
                    source_path=file_path_lane.unsorted_md5,
                    target_path=file_path_lane.sorted_md5)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)
            else:
                runnable_step = RunnableStepPicard(
                    name='picard_sort_sam',
                    obsolete_file_path_list=[
                        file_path_lane.unsorted_bam,
                        file_path_lane.unsorted_md5,
                    ],
                    java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx18G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='SortSam')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # INPUT is required.
                runnable_step.add_picard_option(key='INPUT', value=file_path_lane.unsorted_bam)
                # OUTPUT is required.
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.sorted_bam)
                # SORT_ORDER is required.
                runnable_step.add_picard_option(key='SORT_ORDER', value='queryname')
                # TMP_DIR
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_lane.temporary_directory_path(absolute=False))
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                # OPTIONS_FILE

            # Create the experiment directory if it does not exist already.

            runnable_step = RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move and rename the final, sorted BAM file.

            runnable_step = RunnableStepMove(
                name='move_sorted_bam',
                source_path=file_path_lane.sorted_bam,
                target_path=file_path_lane.archive_bam)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move and rename the checksum file.

            runnable_step = RunnableStepMove(
                name='move_sorted_md5',
                source_path=file_path_lane.sorted_md5,
                target_path=file_path_lane.archive_md5)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

        # Add another flow cell-specific Runnable to reset directory and file mode permissions if requested.

        runnable_cell = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                working_directory=self.project_directory))
        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        runnable_step = RunnableStepChangeMode(
            name='chmod',
            file_path=experiment_directory,
            mode_directory=self.mode_directory,
            mode_file=self.mode_file)
        runnable_cell.add_runnable_step(runnable_step=runnable_step)

        return

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            irf_path: Optional[str] = None,
            force: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` analysis.

        This analysis requires either a :literal:`configuration` argument or a :literal:`irf_path` argument.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param irf_path: An :emphasis:`Illumina Run Folder` path.
        :type irf_path: str | None
        :param force: Request processing incomplete :emphasis:`Illumina Run Folder` objects.
        :type force: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if configuration_path == Configuration.get_global_file_path():
            if not irf_path:
                raise Exception('The --irf argument is required if configuration is not set.')

        analysis: IlluminaToBam = cls.from_config_file_path(config_path=configuration_path)

        if irf_path:
            analysis.run_directory = irf_path

        if force:
            analysis.force = force

        analysis.run()
        analysis.check_state()
        analysis.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis.name)
        print('Project name:           ', analysis.project_name)
        print('Project directory:      ', analysis.project_directory)
        print('Illumina run directory: ', analysis.run_directory)
        print('Experiment directory:   ', analysis.get_experiment_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` analysis.

        This analysis requires either a positional :literal:`configuration` argument or a :literal:`--irf` argument.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description=cls.name + ' submission script.')

        argument_parser.add_argument(
            '--dry-run',
            action='store_false',
            help='dry run',
            dest='drms_submit')

        argument_parser.add_argument(
            '--logging-level',
            default='WARNING',
            choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
            help='logging level [WARNING]')

        argument_parser.add_argument(
            '--stage-name',
            help='limit job submission to a particular analysis stage')

        argument_parser.add_argument(
            '--irf',
            help='Illumina Run Folder name or file path',
            dest='irf_path')

        argument_parser.add_argument(
            '--force',
            action='store_true',
            help='force processing of an incomplete Illumina Run Folder')

        argument_parser.add_argument(
            'configuration',
            nargs='?',
            default=Configuration.get_global_file_path(),
            help=f'configuration (INI) file path [{Configuration.get_global_file_path()!s}]')

        name_space = argument_parser.parse_args()

        if name_space.logging_level:
            logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
            logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

            logging.basicConfig(level=name_space.logging_level)

        return cls.console_submit(
            configuration_path=name_space.configuration,
            stage_name=name_space.stage_name,
            drms_submit=name_space.drms_submit,
            irf_path=name_space.irf_path,
            force=name_space.force)


class FilePathBamIndexDecoderCell(FilePath):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell` class models
    flow cell-specific file paths.

    See also classes:
        - :py:class:`bsf.analyses.picard.FilePathExtractIlluminaCell`
        - :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell`

    :ivar prefix_cell: A non-standard, flow cell-specific (i.e., project_name) prefix.
    :type prefix_cell: str
    :ivar sample_annotation_sheet_csv: A Sample Annotation Sheet CSV file path.
    :type sample_annotation_sheet_csv: str
    """

    def __init__(self, prefix: str, project_name: str) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        """
        super(FilePathBamIndexDecoderCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the flow cell-specific prefix accordingly.

        self.prefix_cell = project_name

        self.sample_annotation_sheet_csv = self.prefix_cell + '_samples.csv'

        return


class FilePathBamIndexDecoderLane(FilePath):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane` class models
    lane-specific file paths.

    See also classes:
        - :py:class:`bsf.analyses.picard.FilePathExtractIlluminaLane`
        - :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane`

    :ivar prefix_lane: A non-standard, lane-specific (i.e., project_name and lane) prefix.
    :type prefix_lane: str
    :ivar project_barcode: A project-specific barcode CSV file path.
    :type project_barcode: str
    :ivar barcode_tsv: A lane-specific barcode TSV file path.
    :type barcode_tsv: str
    :ivar metrics_tsv: A lane-specific metrics TSV file path.
    :type metrics_tsv: str
    :ivar metrics_fraction_pdf: A lane-specific Illumina2bam :literal:`BamIndexDecoder` fraction metrics PDF file path.
    :type metrics_fraction_pdf: str
    :ivar metrics_fraction_png: A lane-specific Illumina2bam :literal:`BamIndexDecoder` fraction metrics PNG file path.
    :type metrics_fraction_png: str
    :ivar metrics_number_pdf: A lane-specific Illumina2bam :literal:`BamIndexDecoder` number metrics PDF file path.
    :type metrics_number_pdf: str
    :ivar metrics_number_png: A lane-specific Illumina2bam :literal:`BamIndexDecoder` number metrics PNG file path.
    :type metrics_number_png: str
    :ivar samples_directory: A directory storing sample-specific unaligned BAM files.
    :type samples_directory: str
    :ivar archive_bam: An :emphasis:`archive` BAM file path.
    :type archive_bam: str
    :ivar archive_md5: An :emphasis:`archive` BAM MD5 checksum file path.
    :type archive_md5: str
    """

    def __init__(self, prefix: str, project_name: str, lane: str, sequences_directory: str) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane` object.
        
        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str
        """
        super(FilePathBamIndexDecoderLane, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the lane-specific prefix accordingly.

        self.prefix_lane = '_'.join((project_name, lane))

        self.project_barcode = self.prefix_lane + '_barcode.tsv'
        self.barcode_tsv = self.prefix_lane + '_barcode.tsv'
        self.metrics_tsv = self.prefix_lane + '_metrics.tsv'
        self.metrics_fraction_pdf = self.prefix_lane + '_metrics_fraction.pdf'
        self.metrics_fraction_png = self.prefix_lane + '_metrics_fraction.png'
        self.metrics_number_pdf = self.prefix_lane + '_metrics_number.pdf'
        self.metrics_number_png = self.prefix_lane + '_metrics_number.png'
        self.samples_directory = self.prefix_lane + '_samples'
        # The input (sequence archive) files are non-standard in that they do not contain a prefix and
        # reside in the sequences_directory.
        self.archive_bam = os.path.join(sequences_directory, self.prefix_lane + '.bam')
        self.archive_md5 = os.path.join(sequences_directory, self.prefix_lane + '.bam.md5')

        return


class BamIndexDecoder(Analysis):
    """The :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` class represents the logic to
    decode sequence archive BAM files into sample-specific BAM files.

    :ivar hash_algorithm: Request using a BSF-specific hashing algorithm for demultiplexing.
    :type hash_algorithm: bool | None
    :ivar library_path: A library annotation file path.
    :type library_path: str | None
    :ivar sequences_directory: A directory storing lane-specific unaligned BAM files.
    :type sequences_directory: str | None
    :ivar samples_directory: A directory storing sample-specific unaligned BAM files.
    :type samples_directory: str | None
    :ivar mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
    :type mode_directory: str | None
    :ivar mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
    :type mode_file: str | None
    :ivar java_classpath_illumina2bam: An Illumina2Bam tools Java Class Path directory path.
    :type java_classpath_illumina2bam: str | None
    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    :ivar lanes: A number of lanes on the flow cell.
    :type lanes: int | None
    :ivar force: Request de-multiplexing with a Library Annotation sheet failing validation.
    :type force: bool | None
    """

    name = 'Bam Index Decoder Analysis'
    prefix = 'bam_index_decoder'

    @classmethod
    def get_stage_name_cell(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cell'))

    @classmethod
    def get_stage_name_lane(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'lane'))

    @classmethod
    def get_prefix_cell(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_cell(), project_name))

    @classmethod
    def get_prefix_lane(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_lane(), project_name, lane))

    @classmethod
    def get_file_path_cell(cls, project_name: str) -> FilePathBamIndexDecoderCell:
        """Get a :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell` object.

        :param project_name: A project name
        :type project_name: str
        :return: A :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell` object.
        :rtype: FilePathBamIndexDecoderCell
        """
        return FilePathBamIndexDecoderCell(
            prefix=cls.get_prefix_cell(project_name=project_name),
            project_name=project_name)

    @classmethod
    def get_file_path_lane(cls, project_name: str, lane: str, sequences_directory: str) -> FilePathBamIndexDecoderLane:
        """Get a :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str
        :return: A :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane` object.
        :rtype: FilePathBamIndexDecoderLane
        """
        return FilePathBamIndexDecoderLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane,
            sequences_directory=sequences_directory)

    @classmethod
    def get_sample_annotation_sheet(cls, file_path: str) -> SampleAnnotationSheet:
        """Get a :py:class:`bsf.ngs.SampleAnnotationSheet` object to annotate :py:class:`bsf.ngs.Sample` objects in a
        flow cell-specific :literal:`BSF samples` directory.

        :param file_path: A file path.
        :type file_path: str
        :return: A :py:class:`bsf.ngs.SampleAnnotationSheet` object.
        :rtype: SampleAnnotationSheet
        """
        return SampleAnnotationSheet(
            file_path=file_path,
            header=True,
            field_name_list=[
                'File Type',
                'ProcessedRunFolder Name',
                'Project Name',
                'Project Size',
                'Sample Name',
                'PairedReads Exclude',
                'PairedReads Flow Cell',
                'PairedReads Flow Cell Lane',
                'PairedReads Index 1',
                'PairedReads Index 2',
                'PairedReads Lane',
                'PairedReads ReadGroup',
                'PairedReads Structure',
                'Reads1 Name',
                'Reads1 File',
                'Reads2 Name',
                'Reads2 File',
            ])

    def __init__(
            self,
            configuration: Optional[Configuration] = None,
            project_name: Optional[str] = None,
            genome_version: Optional[str] = None,
            input_directory: Optional[str] = None,
            output_directory: Optional[str] = None,
            project_directory: Optional[str] = None,
            genome_directory: Optional[str] = None,
            report_style_path: Optional[str] = None,
            report_header_path: Optional[str] = None,
            report_footer_path: Optional[str] = None,
            e_mail: Optional[str] = None,
            stage_list: Optional[list[Stage]] = None,
            collection: Optional[Collection] = None,
            sample_list: Optional[list[Sample]] = None,
            hash_algorithm: Optional[bool] = None,
            library_path: Optional[str] = None,
            sequences_directory: Optional[str] = None,
            samples_directory: Optional[str] = None,
            mode_directory: Optional[str] = None,
            mode_file: Optional[str] = None,
            java_classpath_illumina2bam: Optional[str] = None,
            java_archive_picard: Optional[str] = None,
            lanes: Optional[int] = 8,
            force: Optional[bool] = False) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` object.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration | None
        :param project_name: A project name.
        :type project_name: str | None
        :param genome_version: A genome assembly version.
        :type genome_version: str | None
        :param input_directory: An input directory path.
        :type input_directory: str | None
        :param output_directory: An output directory path.
        :type output_directory: str | None
        :param project_directory: A project directory path, normally under the output directory path.
        :type project_directory: str | None
        :param genome_directory: A genome directory path, normally under the project directory path.
        :type genome_directory: str | None
        :param report_style_path: A report style :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: A report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: A report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a :emphasis:`UCSC Genome Browser Track Hub`.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param hash_algorithm: Request using a BSF-specific hashing algorithm for demultiplexing.
        :type hash_algorithm: bool | None
        :param library_path: A library annotation file path.
        :type library_path: str | None
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str | None
        :param samples_directory: A directory storing sample-specific unaligned BAM files.
        :type samples_directory: str | None
        :param mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
        :type mode_directory: str | None
        :param mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
        :type mode_file: str | None
        :param java_classpath_illumina2bam: An Illumina2Bam tools Java Class Path directory path.
        :type java_classpath_illumina2bam: str | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param lanes: A number of lanes on the flow cell.
        :type lanes: int | None
        :param force: Request de-multiplexing with a Library Annotation sheet failing validation.
        :type force: bool | None
        """
        super(BamIndexDecoder, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            report_style_path=report_style_path,
            report_header_path=report_header_path,
            report_footer_path=report_footer_path,
            e_mail=e_mail,
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list)

        # Sub-class specific ...

        self.hash_algorithm = hash_algorithm
        self.library_path = library_path
        self.sequences_directory = sequences_directory
        self.samples_directory = samples_directory
        self.mode_directory = mode_directory
        self.mode_file = mode_file
        self.java_classpath_illumina2bam = java_classpath_illumina2bam
        self.java_archive_picard = java_archive_picard
        self.lanes = lanes
        self.force = force

        return

    @property
    def get_experiment_directory(self) -> Optional[str]:
        """Get the experiment directory.

        The experiment directory is a concatenation of the samples directory and the project name.

        :return: An experiment directory.
        :rtype: str | None
        """
        if self.samples_directory and self.project_name:
            return os.path.join(self.samples_directory, self.project_name)

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(BamIndexDecoder, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'hash_algorithm'
        if configuration.config_parser.has_option(section=section, option=option):
            self.hash_algorithm = self.configuration.config_parser.getboolean(section=section, option=option)

        option = 'library_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_path = configuration.config_parser.get(section=section, option=option)

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        option = 'samples_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.samples_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        option = 'java_classpath_illumina2bam'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_classpath_illumina2bam = configuration.config_parser.get(section=section, option=option)

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        option = 'lanes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.lanes = configuration.config_parser.getint(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` object to
        decode an archive BAM file produced with Illumina2Bam tools into sample-specific BAM files.
        """

        def run_get_sample_file_name(sample_name: str) -> str:
            """Private function to format sample-specific BAM file names (i.e., :literal:`project_lane#sample.bam`).

            :param sample_name: A sample name.
            :type sample_name: str
            :return: A sample-specific BAM file name.
            :rtype: str
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        # a Picard tools *tab-separated* value (TSV) sample sheet.
        # lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        # barcode_sequence, barcode_name, library_name, sample_name, description

        super(BamIndexDecoder, self).run()

        # Load from the configuration file and override with the default if necessary.

        # Define the sequences and samples directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.sequences_directory:
            self.sequences_directory = self.project_name

        self.sequences_directory = self.configuration.get_absolute_path(
            file_path=self.sequences_directory,
            default_path=StandardFilePath.get_sequences(absolute=True))

        self.samples_directory = self.configuration.get_absolute_path(
            file_path=self.samples_directory,
            default_path=StandardFilePath.get_samples(absolute=True))

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(f"The 'samples_directory' {self.samples_directory!r} does not exist.")

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the library annotation sheet.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(self.library_path):
            raise Exception(f'The library annotation sheet {self.library_path!r} does not exist.')

        # Load the LibraryAnnotationSheet and validate.

        library_annotation_sheet: LibraryAnnotationSheet = LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)

        message_list = library_annotation_sheet.validate(lanes=self.lanes)

        if message_list:
            if self.force:
                warnings.warn(f'Validation of library annotation sheet {self.library_path!r}:\n'
                              '\n'.join(message_list), UserWarning)
            else:
                raise Exception(f'Validation of library annotation sheet {self.library_path!r}:\n'
                                '\n'.join(message_list))

        # Get the Illumina2Bam tools Java Class Path directory.

        if not self.java_classpath_illumina2bam:
            self.java_classpath_illumina2bam = JavaClassPath.get_illumina2bam()
            if not self.java_classpath_illumina2bam:
                raise Exception(f"A {self.name!s} requires a 'java_classpath_illumina2bam' configuration option.")

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        library_annotation_dict = library_annotation_sheet.get_annotation_dict()
        library_barcode_dict = library_annotation_sheet.get_barcode_length_dict()

        file_path_cell = self.get_file_path_cell(project_name=self.project_name)

        # Create a SampleAnnotationSheet in the project directory and
        # eventually transfer it into the experiment directory.
        sample_annotation_sheet = self.get_sample_annotation_sheet(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        cell_dependency_list: list[str] = list()

        for lane_int in sorted(library_annotation_dict):
            lane_str = str(lane_int)

            # Determine the read structure from the first line of the lane-specific library annotation.
            # Getting it from the run configuration via the bsf.illumina.RunInformation.get_picard_read_structure()
            # method would require a bsf.illumina.IlluminaRunFolder object.
            # See bsf.analyses.picard.ExtractIlluminaRunFolder.run() and
            # bsf.analyses.picard.IlluminaDemultiplexSam.run().

            if 'read_structure' in library_annotation_dict[lane_int][0] and \
                    library_annotation_dict[lane_int][0]['read_structure']:
                read_structure = library_annotation_dict[lane_int][0]['read_structure']
            else:
                # read_structure = self._irf.run_information.get_picard_read_structure
                read_structure = ''

            file_path_lane = self.get_file_path_lane(
                project_name=self.project_name,
                lane=lane_str,
                sequences_directory=self.sequences_directory)

            # Do not check whether the sorted input BAM file exists, because at the time of
            # BamIndexDecoder submission the IlluminaToBam analysis may not have finished.

            if library_barcode_dict[lane_int][1] > 0:
                barcode_number = 2
            elif library_barcode_dict[lane_int][0] > 0:
                barcode_number = 1
            else:
                barcode_number = 0

            bam_index_decoder_sheet = BamIndexDecoderSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.barcode_tsv))

            for row_dict in library_annotation_dict[lane_int]:
                # Add a row to the lane-specific tab-delimited Illumina2bam tools BamIndexDecoder barcode file.

                bam_index_decoder_sheet.row_dict_list.append({
                    'barcode_sequence': row_dict['barcode_sequence_1'] + row_dict['barcode_sequence_2'],
                    'barcode_name': row_dict['sample_name'],
                    'library_name': row_dict['library_name'],
                    'sample_name': row_dict['sample_name'],
                    'description': str(),
                })

                # Add a row to the flow cell-specific sample annotation sheet.

                sample_annotation_sheet.row_dict_list.append({
                    'File Type': '',
                    'ProcessedRunFolder Name': self.project_name,
                    'Project Name': row_dict['library_name'],
                    'Project Size': row_dict['library_size'],
                    'Sample Name': row_dict['sample_name'],
                    'PairedReads Exclude': 'FALSE',
                    'PairedReads Flow Cell': self.project_name,
                    'PairedReads Flow Cell Lane': '_'.join((self.project_name, lane_str)),
                    'PairedReads Index 1': row_dict['barcode_sequence_1'],
                    'PairedReads Index 2': row_dict['barcode_sequence_2'],
                    'PairedReads Lane': lane_str,
                    # TODO: It would be good to add a RunnableStep to populate the ReadGroup.
                    'PairedReads ReadGroup': '',
                    'PairedReads Structure': read_structure,
                    'Reads1 File': os.path.join(
                        os.path.basename(experiment_directory),
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    'Reads1 Name': '_'.join((self.project_name, lane_str, row_dict['sample_name'])),
                    'Reads2 File': '',
                    'Reads2 Name': '',
                })

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the Stage.
            runnable_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))

            # TODO: It would be good to extend the Runnable so that it holds dependencies on other Runnable objects
            #  and that it could be submitted to a Stage so that the Executable gets automatically created and
            #  submitted.
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)
            executable_lane.dependencies.append(
                IlluminaToBam.get_prefix_lane(
                    project_name=self.project_name,
                    lane=lane_str))

            cell_dependency_list.append(executable_lane.name)

            if executable_lane.submit:
                # Only if this Executable actually gets submitted ...
                # Write the lane-specific BamIndexDecoderSheet to the internal file path.

                bam_index_decoder_sheet.to_file_path()

            # Create a samples_directory in the project_directory.

            runnable_step = RunnableStepMakeDirectory(
                name='make_samples_directory',
                directory_path=file_path_lane.samples_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            if barcode_number:
                # Run the BamIndexDecoder if there is at least one line containing a barcode sequence.
                runnable_step = RunnableStepIlluminaToBam(
                    name='bam_index_decoder',
                    java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx4G',
                    itb_classpath=self.java_classpath_illumina2bam,
                    itb_command='BamIndexDecoder')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # INPUT is required
                runnable_step.add_itb_option(key='INPUT', value=file_path_lane.archive_bam)
                # OUTPUT is required, but cannot be used together with OUTPUT_FORMAT, OUTPUT_PREFIX and OUTPUT_DIR.
                # OUTPUT_DIR is required, but cannot be used with OUTPUT.
                runnable_step.add_itb_option(key='OUTPUT_DIR', value=file_path_lane.samples_directory)
                # OUTPUT_PREFIX is required, but cannot be used with OUTPUT.
                runnable_step.add_itb_option(key='OUTPUT_PREFIX', value='_'.join((self.project_name, lane_str)))
                # OUTPUT_FORMAT is required, but cannot be used with OUTPUT.
                runnable_step.add_itb_option(key='OUTPUT_FORMAT', value='bam')
                # BARCODE_TAG_NAME defaults to 'BC'.
                # BARCODE_QUALITY_TAG_NAME defaults to 'QT'.
                # BARCODE cannot be used with BARCODE_FILE.
                # BARCODE_FILE is required, but cannot be used with BARCODE.
                runnable_step.add_itb_option(key='BARCODE_FILE', value=file_path_lane.barcode_tsv)
                # METRICS_FILE is required.
                runnable_step.add_itb_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
                # MAX_MISMATCHES defaults to '1'.
                if 'barcode_mismatches' in library_annotation_dict[lane_int][0] and \
                        library_annotation_dict[lane_int][0]['barcode_mismatches']:
                    # If the barcode_mismatches field exists and has a meaningful value ...
                    runnable_step.add_itb_option(
                        key='MAX_MISMATCHES',
                        value=library_annotation_dict[lane_int][0]['barcode_mismatches'])
                elif barcode_number == 2:
                    runnable_step.add_itb_option(key='MAX_MISMATCHES', value='2')
                # MIN_MISMATCH_DELTA defaults to '1'.
                # MAX_NO_CALLS defaults to '2'.
                # CONVERT_LOW_QUALITY_TO_NO_CALL defaults to 'false'.
                # MAX_LOW_QUALITY_TO_CONVERT defaults to '15'.
                # USE_HASH_DEMULTIPLEXING defaults to 'false'.
                if self.hash_algorithm:
                    runnable_step.add_itb_option(key='USE_HASH_DEMULTIPLEXING', value='true')
                # BARCODE_STARTBASE_INDEX defaults to '0'.
                if 'barcode_start' in library_annotation_dict[lane_int][0] and \
                        library_annotation_dict[lane_int][0]['barcode_start']:
                    runnable_step.add_itb_option(
                        key='BARCODE_STARTBASE_INDEX',
                        value=library_annotation_dict[lane_int][0]['barcode_start'])
                # TMP_DIR
                runnable_step.add_itb_option(
                    key='TMP_DIR',
                    value=runnable_lane.temporary_directory_path(absolute=False))
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_itb_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                runnable_step.add_itb_option(key='COMPRESSION_LEVEL', value='9')
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                runnable_step.add_itb_option(key='CREATE_MD5_FILE', value='true')
                # OPTIONS_FILE

                # Plot the metrics file.

                runnable_step = RunnableStep(
                    name='plot_metrics',
                    program='bsf_illumina_demultiplex_sam.R')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='file-path', value=file_path_lane.metrics_tsv)

                # Create the experiment directory if it does not exist already.

                runnable_step = RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Move the sample directory into the experiment directory.

                runnable_step = RunnableStepMove(
                    name='move_samples_directory',
                    source_path=file_path_lane.samples_directory,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Move the metrics file into the experiment directory.

                runnable_step = RunnableStepMove(
                    name='move_metrics_tsv',
                    source_path=file_path_lane.metrics_tsv,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step = RunnableStepMove(
                    name='move_metrics_pdf',
                    source_path=file_path_lane.metrics_number_pdf,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step = RunnableStepMove(
                    name='move_metrics_png',
                    source_path=file_path_lane.metrics_number_png,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)
            else:
                # Run Picard CollectAlignmentSummaryMetrics if there is no line containing a barcode sequence.
                runnable_step = RunnableStepPicard(
                    name='picard_collect_alignment_summary_metrics',
                    obsolete_file_path_list=[
                        file_path_lane.barcode_tsv,
                    ],
                    java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx4G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='CollectAlignmentSummaryMetrics')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # MAX_INSERT_SIZE defaults to '100000'.
                # ADAPTER_SEQUENCE
                # METRIC_ACCUMULATION_LEVEL.
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP')
                # IS_BISULFITE_SEQUENCED defaults to 'false'.
                # INPUT is required.
                runnable_step.add_picard_option(key='INPUT', value=file_path_lane.archive_bam)
                # OUTPUT is required.
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.metrics_tsv)
                # REFERENCE_SEQUENCE defaults to 'null'.
                # ASSUME_SORTED defaults to 'true'.
                # STOP_AFTER defaults to '0'.
                # TMP_DIR
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_lane.temporary_directory_path(absolute=False))
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                # OPTIONS_FILE

                # TODO: It would be even better to run Picard AddOrReplaceReadGroups to get a correct SAM header.
                #  Unfortunately, the RG-LB and RG-SM fields cannot easily be overridden since all options need
                #  to be set. However, the RG-PU and other information is only available after the Illumina2bam stage
                #  has completed, which may be well after submission of this analysis. The solution could be another
                #  script to rad the BAM file and propagate RG information.

                # Create the experiment directory if it does not exist already.

                runnable_step = RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Move the sample directory into the experiment directory.

                runnable_step = RunnableStepMove(
                    name='move_samples_directory',
                    source_path=file_path_lane.samples_directory,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Add a symbolic link to the BSF Sequence Archive file within the sample directory.

                experiment_samples_directory = os.path.join(experiment_directory, file_path_lane.samples_directory)

                runnable_step = RunnableStepLink(
                    name='link_bam',
                    source_path=os.path.relpath(file_path_lane.archive_bam, experiment_samples_directory),
                    target_path=os.path.join(
                        experiment_samples_directory,
                        '{}_{}#{}.bam'.format(self.project_name, lane_str,
                                              library_annotation_dict[lane_int][0]['sample_name'])))
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step = RunnableStepLink(
                    name='link_md5',
                    source_path=os.path.relpath(file_path_lane.archive_md5, experiment_samples_directory),
                    target_path=os.path.join(
                        experiment_samples_directory,
                        '{}_{}#{}.bam.md5'.format(self.project_name, lane_str,
                                                  library_annotation_dict[lane_int][0]['sample_name'])))
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Move the metrics file into the experiment directory.

                runnable_step = RunnableStepMove(
                    name='move_metrics_tsv',
                    source_path=file_path_lane.metrics_tsv,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

        # Finally, write the flow cell-specific bsf.ngs.SampleAnnotationSheet to the internal file path.

        sample_annotation_sheet.to_file_path()

        # Create a flow-cell specific Runnable.

        runnable_cell = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                working_directory=self.project_directory))

        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        # Move the SampleAnnotationSheet from the project_directory to the experiment_directory.

        if os.path.exists(sample_annotation_sheet.file_path):
            runnable_step = RunnableStepMove(
                name='move_sample_annotation',
                source_path=file_path_cell.sample_annotation_sheet_csv,
                target_path=experiment_directory)
            runnable_cell.add_runnable_step(runnable_step=runnable_step)

        # Change directory and file access permissions.

        runnable_step = RunnableStepChangeMode(
            name='chmod',
            file_path=experiment_directory,
            mode_directory=self.mode_directory,
            mode_file=self.mode_file)
        runnable_cell.add_runnable_step(runnable_step=runnable_step)

        return

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            project_name: Optional[str] = None,
            library_path: Optional[str] = None,
            mode: Optional[str] = None,
            force: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` analysis.

        This analysis requires either a :literal:`configuration` argument or a :literal:`project_name`.
        The :literal:`library_path` argument can be automatically determined on the basis of the
        :literal:`project_name`.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param project_name: A project name.
        :type project_name: str | None
        :param library_path: A library path.
        :type library_path: str | None
        :param mode: A sequencer mode.
        :type mode: str | None
        :param force: Force processing of imperfect library annotation sheets.
        :type force: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if configuration_path == Configuration.get_global_file_path():
            if not project_name:
                raise Exception('The --project-name argument is required if configuration is not set.')

        analysis: BamIndexDecoder = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (INI) file.')

            analysis.project_name = project_name

        if mode:
            if mode == 'high':
                analysis.lanes = 8
            elif mode == 'rapid':
                analysis.lanes = 2
            elif mode == 'miseq':
                analysis.lanes = 1
            elif mode == 'nextseq':
                analysis.lanes = 4
            else:
                raise Exception(f'The --mode option {mode!r} is not supported.')

        if force:
            analysis.force = force

        if library_path:
            analysis.library_path = library_path

        # If a library file has not been defined so far, check,
        # if a standard library file i.e. PROJECT_NAME_libraries.csv exists in the current directory.

        if not analysis.library_path:
            library_path = '_'.join((analysis.project_name, 'libraries.csv'))
            if os.path.exists(library_path):
                analysis.library_path = library_path

        analysis.run()
        analysis.check_state()
        analysis.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis.name)
        print('Project name:         ', analysis.project_name)
        print('Project directory:    ', analysis.project_directory)
        print('Sequences directory:  ', analysis.sequences_directory)
        print('Experiment directory: ', analysis.get_experiment_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a :py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` analysis.

        This analysis requires either a positional :literal:`configuration` argument or a :literal:`--project-name`.
        The :literal:`--library-path` argument can be automatically determined on the basis of the
        :literal:`--project-name`.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description=cls.name + ' submission script.')

        argument_parser.add_argument(
            '--dry-run',
            action='store_false',
            help='dry run',
            dest='drms_submit')

        argument_parser.add_argument(
            '--logging-level',
            default='WARNING',
            choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
            help='logging level [WARNING]')

        argument_parser.add_argument(
            '--stage-name',
            help='limit job submission to a particular analysis stage')

        argument_parser.add_argument(
            '--project-name',
            help='project name i.e. flow cell identifier')

        argument_parser.add_argument(
            '--library-path',
            help='library annotation sheet (CSV) file path')

        argument_parser.add_argument(
            '--mode',
            help='HiSeq run mode i.e. high (high-output) or rapid (rapid run)')

        argument_parser.add_argument(
            '--force',
            action='store_true',
            help='force processing even if library annotation sheet validation fails')

        argument_parser.add_argument(
            'configuration',
            nargs='?',
            default=Configuration.get_global_file_path(),
            help=f'configuration (INI) file path [{Configuration.get_global_file_path()!s}]')

        name_space = argument_parser.parse_args()

        if name_space.logging_level:
            logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
            logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

            logging.basicConfig(level=name_space.logging_level)

        return cls.console_submit(
            configuration_path=name_space.configuration,
            stage_name=name_space.stage_name,
            drms_submit=name_space.drms_submit,
            project_name=name_space.project_name,
            library_path=name_space.library_path,
            mode=name_space.mode,
            force=name_space.force)
