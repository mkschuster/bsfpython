# -*- coding: utf-8 -*-
"""IlluminaToBam Tools Analysis module.

A package of classes and methods supporting analyses of the Illumina2Bam-Tools package.
"""
#  Copyright 2013 - 2019 Michael K. Schuster
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
import os
import warnings
from subprocess import Popen

from bsf.analyses.illumina_run_folder import IlluminaRunFolderRestore
from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.ngs import Collection, Sample, SampleAnnotationSheet
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, Executable, RunnableStep, RunnableStepChangeMode, RunnableStepJava, \
    RunnableStepMakeDirectory, RunnableStepLink, RunnableStepMove, RunnableStepPicard
from bsf.standards import Configuration, StandardFilePath, JavaArchive, JavaClassPath, Operator, VendorQualityFilter


class BamIndexDecoderSheet(AnnotationSheet):
    """The C{bsf.analyses.illumina_to_bam_tools.BamIndexDecoderSheet} class represents a
    Tab-Separated Value (TSV) table of library information for the
    C{bsf.analyses.illumina_to_bam_tools.BamIndexDecoder} C{bsf.analysis.Analysis}.

    Attributes:
    """

    _file_type = 'excel-tab'

    # The field names are defined in the IndexDecoder.java source file.
    # https://github.com/wtsi-npg/illumina2bam/blob/devel/src/uk/ac/sanger/npg/picard/IndexDecoder.java

    _header_line = True

    _field_names = [
        'barcode_sequence',
        'barcode_name',
        'library_name',
        'sample_name',
        'description',
    ]

    _test_methods = dict()


class LibraryAnnotationSheet(AnnotationSheet):
    """The C{LibraryAnnotationSheet} class represents a Comma-Separated Value (CSV) table for de-multiplexing.

    Attributes:
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'lane',  # Lane number
        'barcode_mismatches',  # Number of mismatches allowed in index decoding
        'barcode_sequence_1',  # Index read sequence 1
        'barcode_sequence_2',  # Index read sequence 2
        'barcode_start',  # Start position of the barcode in the barcode read (numeric)
        'sample_name',  # Sample name (alphanumeric including '_' characters)
        'library_name',  # Library name (alphanumeric including '_' characters)
        'library_size',  # Library size (numeric)
        'read_structure',  # Picard tools Read Structure
    ]

    _test_methods = {
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

    def get_annotation_dict(self):
        """Get a Python C{dict} of Python C{int} (lane) key and Python C{list} of Python C{dict} (row dict) values.

        @return: Python C{dict} of Python C{int} (lane) key and Python C{list} of Python C{dict} (row dict) values
        @rtype: dict[int, list[dict[str, str]]]
        """
        flow_cell_dict = dict()
        """ @type flow_cell_dict: dict[int, list[dict[str, str]]] """

        for row_dict in self.row_dicts:
            lane_int = int(row_dict['lane'])
            if lane_int not in flow_cell_dict:
                flow_cell_dict[lane_int] = list()
            flow_cell_dict[lane_int].append(row_dict)

        return flow_cell_dict

    def get_barcode_length_dict(self):
        """Get a Python C{dict} of i7 and i5 barcode sequence lengths by lane.

        @return: Python C{dict} of Python C{int} (lane) key and Python C{tuple} of Python C{int} (i7) and
            Python C{int} (i5) barcode sequence length value data
        @rtype: dict[int, (int, int)]
        """
        flow_cell_dict = dict()
        """ @type flow_cell_dict: dict[int, (int, int)] """

        for row_dict in self.row_dicts:
            lane_int = int(row_dict['lane'])
            if lane_int not in flow_cell_dict:
                flow_cell_dict[lane_int] = (len(row_dict['barcode_sequence_1']), len(row_dict['barcode_sequence_2']))
                continue
            for index in range(0, 1 + 1):
                if flow_cell_dict[lane_int][index] != len(row_dict['barcode_sequence_' + str(index + 1)]):
                    warnings.warn("The length of 'barcode_sequence_" + str(index + 1) +
                                  "' does not match previous records.\n" + repr(row_dict))

        return flow_cell_dict

    def validate(self, lanes=None):
        """
        Validate a C{LibraryAnnotationSheet}.

        If the number of lanes is provided, all lanes in the range are validated or else,
        the validation is just based on annotation for the lane variable.
        @param lanes: Number of lanes to validate
        @type lanes: int | None
        @return: Warning messages
        @rtype: str
        """
        messages = str()

        # Check the header line via the pre-defined field names.

        for index in range(0, len(self._field_names)):
            if not self.field_names[index]:
                messages += 'Column with name {!r} is missing from the header line.\n'. \
                    format(self._field_names[index])

            if self.field_names[index] != self._field_names[index]:
                messages += 'Column name {!r} in the header line does not match template {!r}.\n'. \
                    format(self.field_names[index], self._field_names[index])

        # Validate the field values for alphanumeric or sequence grade in the context of the
        # AnnotationSheet super-class.

        messages += super(LibraryAnnotationSheet, self).validate()

        flow_cell_dict = dict()
        # TODO: Since the Python dict holds a dict of complex types, another object would be justified.
        # See barcode_dict: {}, sample_dict: {} and library_name: ''

        row_number = 0

        for row_dict in self.row_dicts:
            row_number += 1

            # Check that all required fields are defined.

            if row_dict['lane'] not in flow_cell_dict:
                flow_cell_dict[row_dict['lane']] = {
                    'barcode_dict': {},  # Barcode sequence to row number.
                    'sample_dict': {},  # Sample name to row number.
                    'library_name': row_dict['library_name'],
                    # 'barcode_mismatches': row_dict['barcode_mismatches'],
                    # 'barcode_start': row_dict['barcode_start'],
                }

                if 'barcode_mismatches' in row_dict:
                    flow_cell_dict[row_dict['lane']]['barcode_mismatches'] = row_dict['barcode_mismatches']

                if 'barcode_start' in row_dict:
                    flow_cell_dict[row_dict['lane']]['barcode_start'] = row_dict['barcode_start']

                if 'read_structure' in row_dict:
                    flow_cell_dict[row_dict['lane']]['read_structure'] = row_dict['read_structure']

            barcode_sequence = str()

            if 'barcode_sequence_1' in row_dict and row_dict['barcode_sequence_1']:
                barcode_sequence += row_dict['barcode_sequence_1']
            else:
                barcode_sequence += '-NoIndex-'

            if 'barcode_sequence_2' in row_dict and row_dict['barcode_sequence_2']:
                barcode_sequence += row_dict['barcode_sequence_2']
            else:
                barcode_sequence += '-NoIndex-'

            if barcode_sequence in flow_cell_dict[row_dict['lane']]['barcode_dict']:
                messages += 'Barcode sequence {!r} from row {} duplicated in row {}.\n'. \
                    format(barcode_sequence,
                           flow_cell_dict[row_dict['lane']]['barcode_dict'][barcode_sequence],
                           row_number)
            else:
                flow_cell_dict[row_dict['lane']]['barcode_dict'][barcode_sequence] = row_number

            # Check for duplicate sample name values.
            if row_dict['sample_name'] in flow_cell_dict[row_dict['lane']]['sample_dict']:
                messages += 'Sample name {!r} from row {} duplicated in row {}.\n'. \
                    format(row_dict['sample_name'],
                           flow_cell_dict[row_dict['lane']]['sample_dict'][row_dict['sample_name']],
                           row_number)
            else:
                flow_cell_dict[row_dict['lane']]['sample_dict'][row_dict['sample_name']] = row_number

            # Check for identical library name values.
            if flow_cell_dict[row_dict['lane']]['library_name'] != row_dict['library_name']:
                messages += 'Library name {!r} in row {} does not match previous name {!r}.\n'. \
                    format(row_dict['library_name'],
                           row_number,
                           flow_cell_dict[row_dict['lane']]['library_name'])

            # Check for identical barcode mismatches values.
            if 'barcode_mismatches' in row_dict:
                if flow_cell_dict[row_dict['lane']]['barcode_mismatches'] != row_dict['barcode_mismatches']:
                    messages += 'Barcode mismatches {!r} in row {} does not match previous mismatches {!r}.\n'. \
                        format(row_dict['barcode_mismatches'],
                               row_number,
                               flow_cell_dict[row_dict['lane']]['barcode_mismatches'])

            # Check for identical barcode start values.
            if 'barcode_start' in row_dict:
                if flow_cell_dict[row_dict['lane']]['barcode_start'] != row_dict['barcode_start']:
                    messages += 'Barcode start {!r} in row {} does not match previous start {!r}.\n'. \
                        format(row_dict['barcode_start'],
                               row_number,
                               flow_cell_dict[row_dict['lane']]['barcode_start'])

            # Check for identical read structure values.
            if 'read_structure' in row_dict:
                if flow_cell_dict[row_dict['lane']]['read_structure'] != row_dict['read_structure']:
                    messages += 'Read structure {!r} in row {} does not match previous read structure {!r}.\n'. \
                        format(row_dict['read_structure'],
                               row_number,
                               flow_cell_dict[row_dict['lane']]['read_structure'])

        if lanes is None:
            # If the number of lanes was not provided, validate just on the basis if the lane annotation.
            lane_list = sorted(flow_cell_dict)
        else:
            # If the number of lanes was provided, validate all lanes in the range.
            lane_list = range(0 + 1, lanes + 1)

        for lane_int in lane_list:
            lane_str = str(lane_int)

            # Check that all lanes have annotation.
            if lane_str not in flow_cell_dict:
                messages += 'No annotation for lane number ' + repr(lane_int) + '.\n'
                continue

            # Check that all or none of the rows have index sequence 1 or 2 populated.
            no_index_1 = 0
            no_index_2 = 0
            for key in flow_cell_dict[lane_str]['barcode_dict']:
                if key[:9] == '-NoIndex-':
                    no_index_1 += 1
                if key[-9:] == '-NoIndex-':
                    no_index_2 += 1

            if not (no_index_1 == 0 or no_index_1 == len(flow_cell_dict[lane_str]['barcode_dict'])):
                messages += 'Some empty barcode_sequence_1 fields in lane ' + repr(lane_int) + '.\n'
            if not (no_index_2 == 0 or no_index_2 == len(flow_cell_dict[lane_str]['barcode_dict'])):
                messages += 'Some empty barcode_sequence_2 fields in lane ' + repr(lane_int) + '.\n'

            # Check that all barcode sequences have the same length.
            # This test also finds cases of missing sequences tested for above.
            key_length = None
            for key in flow_cell_dict[lane_str]['barcode_dict']:
                if key_length is None:
                    key_length = len(key)
                if len(key) != key_length:
                    messages += 'Mismatching barcode sequence lengths in lane ' + repr(lane_int) + '.\n'

        return messages


class RunnableStepIlluminaToBam(RunnableStepJava):
    """The C{bsf.analyses.illumina_to_bam_tools.RunnableStepIlluminaToBam} class represents a
    C{bsf.process.RunnableStepJava} specific to IlluminaToBam tools.

    IlluminaToBam tools use the old Picard tools interface where each algorithm is implemented as a separate
    Java Archive (JAR) file.

    Attributes:
    @ivar itb_classpath: IlluminaToBam Java Class Path directory
    @type itb_classpath: str | None
    @ivar itb_command: IlluminaToBam command
    @type itb_command: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_maximum=None,
            java_jar_path=None,
            itb_classpath=None,
            itb_command=None):
        """Initialise a C{bsf.process.RunnableStepIlluminaToBam} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: bsf.connector.Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: bsf.connector.Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: bsf.connector.Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} into the C{bsf.analysis.Stage}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java Archive (JAR) file path
        @type java_jar_path: str | None
        @param itb_classpath: IlluminaToBam Java Class Path directory
        @type itb_classpath: str | None
        @param itb_command: IlluminaToBam command
        @type itb_command: str | None
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

        # Set the IlluminaToBam Java Class Path and the command-specific IlluminaToBam Java Archive (JAR) file path.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(itb_classpath, itb_command + '.jar'))

        # Introduce another empty sub-command to separate Illumina2Bam-tools options.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command()

        return

    def add_itb_option(self, key, value, override=False):
        """Add a C{bsf.argument.OptionPair} to a C{bsf.analyses.illumina_to_bam_tools.RunnableStepIlluminaToBam}.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class FilePathIlluminaToBamLane(FilePath):
    """The C{bsf.analyses.illumina_to_bam_tools.FilePathIlluminaRunFolderRestore} models files in an archive directory.

    Attributes:
    @ivar unsorted_bam: Unsorted BAM file
    @type unsorted_bam: str
    @ivar unsorted_md5: Unsorted BAM file MD5 check sum
    @type unsorted_md5: str
    @ivar sorted_bam: Sorted BAM file
    @type sorted_bam: str
    @ivar sorted_md5: Sorted BAM file MD5 check sum
    @type sorted_md5: str
    @ivar archive_bam: Archive BAM file
    @type archive_bam: str
    @ivar archive_md5: Archive BAM file MD5 check sum
    @type archive_md5: str
    """

    def __init__(self, prefix, project_name, lane, experiment_directory):
        """Initialise a C{bsf.analyses.illumina_to_bam_tools.FilePathIlluminaToBamLane} object.
        
        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param experiment_directory: Experiment-specific directory
        @type experiment_directory: str
        """
        super(FilePathIlluminaToBamLane, self).__init__(prefix=prefix)

        self.unsorted_bam = prefix + '_unsorted.bam'
        self.unsorted_md5 = prefix + '_unsorted.bam.md5'
        self.sorted_bam = prefix + '_sorted.bam'
        self.sorted_md5 = prefix + '_sorted.bam.md5'
        # The final BAM and MD5 files are non-standard in that they do not contain a prefix and
        # reside in the experiment_directory.
        self.archive_bam = os.path.join(experiment_directory, '_'.join((project_name, lane)) + '.bam')
        self.archive_md5 = os.path.join(experiment_directory, '_'.join((project_name, lane)) + '.bam.md5')

        return


class IlluminaToBam(Analysis):
    """The C{bsf.analyses.illumina_to_bam_tools.IlluminaToBam} class represents the logic to
    convert Illumina BCL to a BAM or SAM files.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | None
    @ivar intensity_directory: File path to the I{Intensities} directory,
        defaults to I{illumina_run_folder/Data/Intensities}
    @type intensity_directory: str | None
    @ivar basecalls_directory: File path to the I{BaseCalls} directory,
        defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
    @type basecalls_directory: str | None
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str | None
    @ivar sequencing_centre: Sequencing centre
    @type sequencing_centre: str | None
    @ivar sort_output: Sort BAM files
    @type sort_output: bool | None
    @ivar sequences_directory: Sequences directory to store archive BAM files
    @type sequences_directory: str | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar java_classpath_illumina2bam: Illumina2Bam tools Java Class Path directory
    @type java_classpath_illumina2bam: str | None
    @ivar java_archive_picard: Picard tools Java Archive (JAR) file path
    @type java_archive_picard: str | None
    @ivar vendor_quality_filter: Vendor quality filter
    @type vendor_quality_filter: bool
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool | None
    """

    name = 'Illumina To Bam Analysis'
    prefix = 'illumina_to_bam'

    @classmethod
    def get_stage_name_cell(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cell'))

    @classmethod
    def get_stage_name_lane(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'lane'))

    @classmethod
    def get_prefix_cell(cls, project_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_cell(), project_name))

    @classmethod
    def get_prefix_lane(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_lane(), project_name, lane))

    @classmethod
    def get_file_path_lane(cls, project_name, lane, experiment_directory):
        """Get a C{FilePathIlluminaToBamLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param experiment_directory: BSF experiment directory
        @type experiment_directory: str
        @return: C{FilePathIlluminaToBamLane} object
        @rtype: FilePathIlluminaToBamLane
        """
        return FilePathIlluminaToBamLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane,
            experiment_directory=experiment_directory)

    def __init__(
            self,
            configuration=None,
            project_name=None,
            genome_version=None,
            input_directory=None,
            output_directory=None,
            project_directory=None,
            genome_directory=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            run_directory=None,
            intensity_directory=None,
            basecalls_directory=None,
            experiment_name=None,
            sequencing_centre=None,
            sort_output=None,
            sequences_directory=None,
            mode_directory=None,
            mode_file=None,
            java_classpath_illumina2bam=None,
            java_archive_picard=None,
            vendor_quality_filter=None,
            force=False):
        """Initialise a C{bsf.analyses.illumina_to_bam_tools.IlluminaToBam} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.analysis.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.analysis.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.analysis.Analysis}-wide project directory,
            normally under the C{bsf.analysis.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.analysis.Analysis}-wide genome directory,
            normally under the C{bsf.analysis.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.analysis.Stage} objects
        @type stage_list: list[Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[Sample]
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | None
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | None
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | None
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str | None
        @param sequencing_centre: Sequencing centre
        @type sequencing_centre: str | None
        @param sort_output: Sort BAM files
        @type sort_output: bool | None
        @param sequences_directory: Sequences directory to store archive BAM files
        @type sequences_directory: str | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param java_classpath_illumina2bam: Illumina2Bam tools Java Class Path directory
        @type java_classpath_illumina2bam: str | None
        @param java_archive_picard: Picard tools Java Archive (JAR) file path
        @type java_archive_picard: str | None
        @param vendor_quality_filter: Vendor quality filter
        @type vendor_quality_filter: bool | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
        """
        super(IlluminaToBam, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            e_mail=e_mail,
            debug=debug,
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
    def get_experiment_directory(self):
        """Get the experiment directory.

        The experiment directory is a concatenation of the sequences directory and the project name.
        @return: Experiment directory
        @rtype: str | None
        """
        if self.sequences_directory and self.project_name:
            return os.path.join(self.sequences_directory, self.project_name)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.illumina_to_bam_tools.IlluminaToBam} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
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

    def run(self):
        """Run this C{bsf.analyses.illumina_to_bam_tools.IlluminaToBam} C{bsf.analysis.Analysis}.

        Convert an Illumina flow cell into lane-specific archive BAM files.

        To convert an Illumina flow cell, Illumina2bam is run first, setting the SAM Read Group (@RG)
        library name (LB) and sample name (SM) to 'flow cell identifier.lane'.
        The resulting archive BAM file is then sorted by query name with Picard SortSam.
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception('An Illumina run directory or file path has not been defined.')

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=StandardFilePath.get_illumina_run(absolute=True))

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The Illumina run directory {!r} does not exist.'.format(self.run_directory))

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(
                'The Illumina run directory {!r} is not complete.'.format(self.run_directory))

        # Define an 'Intensities' directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Illumina Run Folder path.

        if self.intensity_directory:
            intensity_directory = self.configuration.get_absolute_path(
                file_path=self.intensity_directory,
                default_path=self.run_directory)
        else:
            intensity_directory = os.path.join(self.run_directory, 'Data', 'Intensities')

        # Check that the Intensities directory exists.

        if not os.path.isdir(intensity_directory):
            raise Exception(
                'The Intensity directory {!r} does not exist.'.format(intensity_directory))

        # Define a 'BaseCalls' directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Intensities directory path.

        if self.basecalls_directory:
            basecalls_directory = self.configuration.get_absolute_path(
                file_path=self.basecalls_directory,
                default_path=intensity_directory)
        else:
            basecalls_directory = os.path.join(intensity_directory, 'BaseCalls')

        # Check that the BaseCalls directory exists.

        if not os.path.isdir(basecalls_directory):
            raise Exception(
                'The BaseCalls directory {!r} does not exist.'.format(basecalls_directory))

        irf = RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g. BSF_0000) is used as the prefix for archive BAM files.
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
                raise Exception("An 'IlluminaToBam' analysis requires a "
                                "'sequencing_centre' configuration option.")

        # Define the sequences directory in which to create the experiment directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # An absolute path cannot be prepended.

        if self.sequences_directory:
            self.sequences_directory = self.configuration.get_absolute_path(
                file_path=self.sequences_directory)
        else:
            self.sequences_directory = StandardFilePath.get_sequences(absolute=True)

        # As a safety measure, to prevent creation of rogue directory paths, the sequences_directory has to exist.

        if not os.path.isdir(self.sequences_directory):
            raise Exception(
                'The IlluminaToBam sequences_directory {!r} does not exist.'.format(self.sequences_directory))

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the Illumina2Bam tools Java Class Path directory.

        if not self.java_classpath_illumina2bam:
            self.java_classpath_illumina2bam = JavaClassPath.get_illumina2bam()
            if not self.java_classpath_illumina2bam:
                raise Exception("An 'IlluminaToBam' analysis requires a "
                                "'java_classpath_illumina2bam' configuration option.")

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception("An 'IlluminaToBam' analysis requires a "
                                "'java_archive_picard' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self.vendor_quality_filter is None:
            self.vendor_quality_filter = VendorQualityFilter.get_vendor_quality_filter(
                flow_cell_type=irf.run_parameters.get_flow_cell_type)

        # Call the run method of the super class after the project_name has been defined.

        super(IlluminaToBam, self).run()

        lane_dependency_str = None
        """ @type lane_dependency_str: str | None """

        cell_dependency_list = list()
        """ @type cell_dependency_list: list[str] """

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        for lane_int in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            file_path_lane = self.get_file_path_lane(
                project_name=self.project_name,
                lane=lane_str,
                experiment_directory=experiment_directory)

            # NOTE: The bsf.procedure.Runnable.name has to match the Executable.name that gets submitted via the Stage.
            runnable_lane = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))

            # TODO: The Runnable class could have dependencies just like the Executable class so that they could be
            # passed on upon creation of the Executable from the Runnable via Executable.from_analysis_runnable().
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

        runnable_cell = self.add_runnable_consecutive(
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


class FilePathBamIndexDecoderCell(FilePath):
    """The C{bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell} models flow cell-specific file paths.

    Attributes:
    @ivar prefix_cell: Non-standard, flow cell-specific (i.e. project_name) prefix
    @type prefix_cell: str
    @ivar sample_annotation_sheet_csv: Sample Annotation Sheet CSV file
    @type sample_annotation_sheet_csv: str
    @see: FilePathIlluminaDemultiplexSamCell
    @see: FilePathExtractIlluminaCell
    """

    def __init__(self, prefix, project_name):
        """Initialise a C{FilePathBamIndexDecoderCell} object.

        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        """
        super(FilePathBamIndexDecoderCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the flow cell-specific prefix accordingly.

        self.prefix_cell = project_name

        self.sample_annotation_sheet_csv = self.prefix_cell + '_samples.csv'

        return


class FilePathBamIndexDecoderLane(FilePath):
    """The C{bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane} models lane-specific file paths.

    Attributes:
    @ivar prefix_lane: Non-standard, lane-specific (i.e. project_name and lane) prefix
    @type prefix_lane: str
    @ivar project_barcode: Project-specific barcode CSV file
    @type project_barcode: str
    @ivar barcode_tsv: Lane-specific barcode TSV file
    @type barcode_tsv: str
    @ivar metrics_tsv: Lane-specific metrics TSV file
    @type metrics_tsv: str
    @ivar metrics_fraction_pdf: Lane-specific Illumina2bam I{BamIndexDecoder} fraction metrics PDF file
    @type metrics_fraction_pdf: str
    @ivar metrics_fraction_png: Lane-specific Illumina2bam I{BamIndexDecoder} fraction metrics PNG file
    @type metrics_fraction_png: str
    @ivar metrics_number_pdf: Lane-specific Illumina2bam I{BamIndexDecoder} number metrics PDF file
    @type metrics_number_pdf: str
    @ivar metrics_number_png: Lane-specific Illumina2bam I{BamIndexDecoder} number metrics PNG file
    @type metrics_number_png: str
    @ivar samples_directory: Sample directory
    @type samples_directory: str
    @ivar archive_bam: Archive BAM file
    @type archive_bam: str
    @ivar archive_md5: Archive BAM MD5 check sum
    @type archive_md5: str
    @see: FilePathExtractIlluminaLane
    @see: FilePathIlluminaDemultiplexSamLane
    """

    def __init__(self, prefix, project_name, lane, sequences_directory):
        """Initialise a C{bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane} object.
        
        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str
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
    """The C{bsf.analyses.illumina_to_bam_tools.BamIndexDecoder} class represents the logic to
    decode sequence archive BAM files into sample-specific BAM files.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar hash_algorithm: Use a BSF-specific hashing algorithm for demultiplexing
    @type hash_algorithm: bool | None
    @ivar library_path: Library annotation file path
    @type library_path: str | None
    @ivar sequences_directory: BSF sequences directory
    @type sequences_directory: str | None
    @ivar samples_directory: BSF samples directory
    @type samples_directory: str | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar java_classpath_illumina2bam: Illumina2Bam tools Java Class Path directory
    @type java_classpath_illumina2bam: str | None
    @ivar java_archive_picard: Picard tools Java Archive (JAR) file path
    @type java_archive_picard: str | None
    @ivar lanes: Number of lanes on the flow cell
    @type lanes: int | None
    @ivar force: Force de-multiplexing with a Library Annotation sheet failing validation
    @type force: bool | None
    """

    name = 'Bam Index Decoder Analysis'
    prefix = 'bam_index_decoder'

    @classmethod
    def get_stage_name_cell(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cell'))

    @classmethod
    def get_stage_name_lane(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'lane'))

    @classmethod
    def get_prefix_cell(cls, project_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_cell(), project_name))

    @classmethod
    def get_prefix_lane(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_lane(), project_name, lane))

    @classmethod
    def get_file_path_cell(cls, project_name):
        """Get a C{FilePathBamIndexDecoderCell} object.

        @param project_name: Project name
        @type project_name: str
        @return: C{FilePathBamIndexDecoderCell} object
        @rtype: FilePathBamIndexDecoderCell
        """
        return FilePathBamIndexDecoderCell(
            prefix=cls.get_prefix_cell(project_name=project_name),
            project_name=project_name)

    @classmethod
    def get_file_path_lane(cls, project_name, lane, sequences_directory):
        """Get a C{FilePathBamIndexDecoderLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str
        @return: C{FilePathBamIndexDecoderLane} object
        @rtype: FilePathBamIndexDecoderLane
        """
        return FilePathBamIndexDecoderLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane,
            sequences_directory=sequences_directory)

    @classmethod
    def get_sample_annotation_sheet(cls, file_path):
        """Get a C{bsf.ngs.SampleAnnotationSheet} object to annotate C{bsf.ngs.Sample} objects in a
        flow cell-specific I{BSF samples} directory.

        @param file_path: File path
        @type file_path: str
        @return: C{bsf.ngs.SampleAnnotationSheet}
        @rtype: SampleAnnotationSheet
        """
        return SampleAnnotationSheet(
            file_path=file_path,
            header=True,
            field_names=[
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
            configuration=None,
            project_name=None,
            genome_version=None,
            input_directory=None,
            output_directory=None,
            project_directory=None,
            genome_directory=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            hash_algorithm=None,
            library_path=None,
            sequences_directory=None,
            samples_directory=None,
            mode_directory=None,
            mode_file=None,
            java_classpath_illumina2bam=None,
            java_archive_picard=None,
            lanes=8,
            force=False):
        """Initialise a C{BamIndexDecoder} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.analysis.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.analysis.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.analysis.Analysis}-wide project directory,
            normally under the C{bsf.analysis.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.analysis.Analysis}-wide genome directory,
            normally under the C{bsf.analysis.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.analysis.Stage} objects
        @type stage_list: list[Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[Sample]
        @param hash_algorithm: Use a BSF-specific hashing algorithm for demultiplexing
        @type hash_algorithm: bool | None
        @param library_path: Library annotation file path
        @type library_path: str | None
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str | None
        @param samples_directory: BSF samples directory
        @type samples_directory: str | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param java_classpath_illumina2bam: Illumina2Bam tools Java Class Path directory
        @type java_classpath_illumina2bam: str | None
        @param java_archive_picard: Picard tools Java Archive (JAR) file path
        @type java_archive_picard: str | None
        @param lanes: Number of lanes on the flow cell
        @type lanes: int | None
        @param force: Force de-multiplexing with a Library Annotation sheet failing validation
        @type force: bool | None
        """
        super(BamIndexDecoder, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            e_mail=e_mail,
            debug=debug,
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
    def get_experiment_directory(self):
        """Get the experiment directory.

        The experiment directory is a concatenation of the samples directory and the project name.
        @return: Experiment directory
        @rtype: str | None
        """
        if self.samples_directory and self.project_name:
            return os.path.join(self.samples_directory, self.project_name)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.illumina_to_bam_tools.BamIndexDecoder} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
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

    def run(self):
        """Run the C{bsf.analyses.illumina_to_bam_tools.BamIndexDecoder} analysis to
        decode an archive BAM file produced with Illumina2Bam tools into sample-specific BAM files.
        """

        def run_get_sample_file_name(sample_name):
            """Private function to format sample-specific BAM file names (i.e. project_lane#sample.bam).

            @param sample_name:
            @type sample_name: str
            @return: Sample-specific BAM file name
            @rtype: str
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        # a Picard tools *tab-separated* value (TSV) sample sheet.
        # lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        # barcode_sequence, barcode_name, library_name, sample_name, description

        super(BamIndexDecoder, self).run()

        # Load from the configuration file and override with the default if necessary.

        # Define the sequences and samples directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
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
            raise Exception(
                'The BamIndexDecoder samples_directory {!r} does not exist.'.format(self.samples_directory))

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the library annotation sheet.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(self.library_path):
            raise Exception('Library annotation file {!r} does not exist.'.format(self.library_path))

        # Load the LibraryAnnotationSheet and validate.

        library_annotation_sheet = LibraryAnnotationSheet.from_file_path(file_path=self.library_path)
        """ @type library_annotation_sheet: LibraryAnnotationSheet """

        validation_messages = library_annotation_sheet.validate(lanes=self.lanes)

        if validation_messages:
            if self.force:
                warnings.warn('Validation of library annotation sheet {!r}:\n{}'.
                              format(self.library_path, validation_messages))
            else:
                raise Exception('Validation of library annotation sheet {!r}:\n{}'.
                                format(self.library_path, validation_messages))

        # Get the Illumina2Bam tools Java Class Path directory.

        if not self.java_classpath_illumina2bam:
            self.java_classpath_illumina2bam = JavaClassPath.get_illumina2bam()
            if not self.java_classpath_illumina2bam:
                raise Exception("An 'BamIndexDecoder' analysis requires a "
                                "'java_classpath_illumina2bam' configuration option.")

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception("An 'BamIndexDecoder' analysis requires a "
                                "'java_archive_picard' configuration option.")

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

        cell_dependency_list = list()
        """ @type cell_dependency_list: list[str] """

        for lane_int in sorted(library_annotation_dict):
            lane_str = str(lane_int)
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
                # Add a row to the lane-specific tab-delimited IlluminaToBamTools BamIndexDecoder barcode file.

                bam_index_decoder_sheet.row_dicts.append({
                    'barcode_sequence': row_dict['barcode_sequence_1'] + row_dict['barcode_sequence_2'],
                    'barcode_name': row_dict['sample_name'],
                    'library_name': row_dict['library_name'],
                    'sample_name': row_dict['sample_name'],
                    'description': str(),
                })

                # Add a row to the flow cell-specific sample annotation sheet.

                sample_annotation_sheet.row_dicts.append({
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
                    # FIXME: For the moment there is no way to get the ReadStructure.
                    # bsf.illumina.RunInformation.get_picard_read_structure() would require a
                    # IlluminaRunFolder object.
                    'PairedReads Structure': '',
                    'Reads1 File': os.path.join(
                        os.path.basename(experiment_directory),
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    'Reads1 Name': '_'.join((self.project_name, lane_str, row_dict['sample_name'])),
                    'Reads2 File': '',
                    'Reads2 Name': '',
                })

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the Stage.
            runnable_lane = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))

            # TODO: It would be good to extend the Runnable so that it holds dependencies on other Runnable objects
            # and that it could be submitted to a Stage so that the Executable gets automatically created and submitted.
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
                # USE_HASH_DEMULTIPLEXING defaults to false.
                if self.hash_algorithm:
                    runnable_step.add_itb_option(key='USE_HASH_DEMULTIPLEXING', value='true')
                # BARCODE_STARTBASE_INDEX defaults to 0.
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

                # Move the samples directory into the experiment directory.

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
                # Unfortunately, the RGLB and RGSM fields cannot easily be overridden since all options need
                # to be set. However, the RGPU and other information is only available after the Illumina2bam stage
                # has completed, which may be well after submission of this analysis. The solution could be another
                # script to rad the BAM file and propagate RG information.

                # Create the experiment directory if it does not exist already.

                runnable_step = RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Move the samples directory into the experiment directory.

                runnable_step = RunnableStepMove(
                    name='move_samples_directory',
                    source_path=file_path_lane.samples_directory,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Add a symbolic link to the BSF Sequence Archive file within the samples directory.

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

        runnable_cell = self.add_runnable_consecutive(
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
