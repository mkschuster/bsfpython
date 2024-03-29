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
"""The :py:mod:`bsf.analyses.picard` module provides classes modelling Picard analyses data files and data directories.
"""
import logging
import os
import shutil
import stat
import warnings
import weakref
from argparse import ArgumentParser
from typing import Optional

from pysam import AlignmentFile

from bsf.analyses.illumina_to_bam_tools import BamIndexDecoder, LibraryAnnotationSheet
from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.executables.cloud import RunnableStepAzureBlockBlobUpload
from bsf.executables.collection import RunnableStepCollectionPruneFastq
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.ngs import Collection, ProcessedRunFolder, Project, Sample, PairedReads, Reads
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import RunnableStep, RunnableStepChangeMode, RunnableStepMakeDirectory, \
    RunnableStepMove, RunnableStepPicard, RunnableStepSleep
from bsf.standards import get_irf_path, SafeFileName, Configuration, StandardFilePath, JavaArchive, Operator, \
    VendorQualityFilter

module_logger = logging.getLogger(name=__name__)


class PicardIlluminaRunFolder(Analysis):
    """The :py:class:`bsf.analyses.picard.PicardIlluminaRunFolder` class models
    Picard Analyses acting on :emphasis:`Illumina Run Folder` entities.

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
    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    :ivar force: Request processing of incomplete :emphasis:`Illumina Run Folder` instances.
    :type force: bool | None
    """

    name = 'Picard PicardIlluminaRunFolder Analysis'
    prefix = 'picard_illumina_run_folder'

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
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_cell(), project_name))

    @classmethod
    def get_prefix_lane(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_lane(), project_name, lane))

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
            java_archive_picard: Optional[str] = None,
            force: Optional[bool] = None):
        """Initialise a :py:class:`bsf.analyses.picard.PicardIlluminaRunFolder` object.

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
        :param basecalls_directory: An :emphasis:`Illumina Run Folder` (IRF) :literal:`BaseCalls` directory path,
            defaults to :literal:`IRF/Data/Intensities/BaseCalls`.
        :type basecalls_directory: str | None
        :param experiment_name: Experiment name (i.e., flow cell identifier) normally automatically read from
            :emphasis:`Illumina Run Folder` parameters.
        :type experiment_name: str | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param force: Request processing of incomplete :emphasis:`Illumina Run Folder` instances.
        :type force: bool | None
        """
        super(PicardIlluminaRunFolder, self).__init__(
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

        self.run_directory = run_directory
        self.intensity_directory = intensity_directory
        self.basecalls_directory = basecalls_directory
        self.experiment_name = experiment_name
        self.java_archive_picard = java_archive_picard
        self.force = force

        self._irf: Optional[RunFolder] = None

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.picard.PicardIlluminaRunFolder` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(PicardIlluminaRunFolder, self).set_configuration(configuration=configuration, section=section)

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

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.PicardIlluminaRunFolder` object.
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(f"A {self.name!s} requires a 'run_directory' configuration option.")

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=StandardFilePath.get_illumina_run(absolute=True))

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(f"The 'run_directory' {self.run_directory!r} does not exist.")

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(f'The Illumina Run Folder {self.run_directory!r} is not complete.')

        # Define an 'Intensities' directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Illumina Run Folder path.

        if self.intensity_directory:
            self.intensity_directory = self.configuration.get_absolute_path(
                file_path=self.intensity_directory,
                default_path=self.run_directory)
        else:
            self.intensity_directory = os.path.join(self.run_directory, 'Data', 'Intensities')

        # Check that the "Intensities" directory exists.

        if not os.path.isdir(self.intensity_directory):
            raise Exception(f"The 'intensity_directory' {self.intensity_directory!r} does not exist.")

        # Define a 'BaseCalls' directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
        # Check if an absolute path has been provided, if not,
        # automatically prepend the "Intensities" directory path.

        if self.basecalls_directory:
            self.basecalls_directory = self.configuration.get_absolute_path(
                file_path=self.basecalls_directory,
                default_path=self.intensity_directory)
        else:
            self.basecalls_directory = os.path.join(self.intensity_directory, 'BaseCalls')

        # Check that the BaseCalls directory exists.

        if not os.path.isdir(self.basecalls_directory):
            raise Exception(f"The 'basecalls_directory' {self.basecalls_directory!r} does not exist.")

        self._irf = RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g., BSF_0000) is used as the prefix for archive BAM files.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = self._irf.run_parameters.get_experiment_name
            if not self.experiment_name:
                raise Exception("The 'experiment_name' was not provided and could not be read from the "
                                "Illumina Run Folder configuration.")

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, self._irf.run_information.flow_cell))

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        # Call the run method of the super class after the project_name has been defined.

        super(PicardIlluminaRunFolder, self).run()

        return


class ExtractIlluminaBarcodesSheet(AnnotationSheet):
    """The :py:class:`bsf.analyses.picard.ExtractIlluminaBarcodesSheet` class represents a
    :emphasis:`Tab-Separated Value` (TSV) table of library information for the
    :py:class:`bsf.analyses.picard.ExtractIlluminaBarcodes` class.
    """

    _file_type = 'excel-tab'

    _field_name_list = [
        'barcode_sequence_1',
        'barcode_sequence_2',
        'barcode_name',
        'library_name',
    ]


class IlluminaBasecallsToSamSheet(AnnotationSheet):
    """The :py:class:`bsf.analyses.picard.IlluminaBasecallsToSamSheet` class represents a
    :emphasis:`Tab-Separated Value` (TSV) table of library information for the
    :py:class:`bsf.analyses.picard.ExtractIlluminaBarcodes` class.
    """

    _file_type = 'excel-tab'

    _field_name_list = [
        'OUTPUT',
        'SAMPLE_ALIAS',
        'LIBRARY_NAME',
        'BARCODE_1',
        'BARCODE_2',
    ]

    def adjust(self, barcode_length_tuple: tuple[int, int], unassigned_file_path: str) -> None:
        """Adjust a :py:class:`bsf.analyses.picard.IlluminaBasecallsToSamSheet` object.

            - Insert the file path for unassigned reads.
            - Remove :literal:`BARCODE_N` index columns, which annotation is all empty.

        :param barcode_length_tuple: A Python :py:class:`tuple` object of
            Python :py:class:`int` (barcode length) objects.
        :type barcode_length_tuple: (int, int)
        :param unassigned_file_path: A file path for unassigned reads.
        :type unassigned_file_path: str
        """
        # The IlluminaBasecallsToSamSheet needs adjusting ...
        if len(self.row_dict_list) == 1 and len(self.row_dict_list[0]['BARCODE_1']) == 0 and len(
                self.row_dict_list[0]['BARCODE_2']) == 0:
            # ... if a single sample, but neither BARCODE_1 nor BARCODE_2 were defined,
            # BARCODE_1 needs setting to 'N'.
            self.row_dict_list[0]['BARCODE_1'] = 'N'
        else:
            # ... in all other cases, a last row for unmatched barcode sequences needs adding.
            self.row_dict_list.append({
                'BARCODE_1': 'N' * barcode_length_tuple[0],
                'BARCODE_2': 'N' * barcode_length_tuple[1],
                'OUTPUT': unassigned_file_path,
                'SAMPLE_ALIAS': 'Unmatched',
                'LIBRARY_NAME': self.row_dict_list[0]['LIBRARY_NAME'],
            })

        # Adjust the IlluminaBaseCallsToSamSheet and remove any BARCODE_N columns not represented
        # in the barcode length list.
        for index in range(0, 1 + 1):
            if barcode_length_tuple[index] == 0:
                # Remove the 'BARCODE_N' field from the list of field names.
                if 'BARCODE_' + str(index + 1) in self.field_name_list:
                    self.field_name_list.remove('BARCODE_' + str(index + 1))
                # Remove the 'BARCODE_N' entry from each row_dict object, since csv.DictWriter requires it.
                for row_dict in self.row_dict_list:
                    row_dict.pop('BARCODE_' + str(index + 1), None)

        return


class FilePathExtractIlluminaCell(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathExtractIlluminaCell` class models flow cell-specific file paths.

    See also classes:
        - :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell`
        - :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell`

    :ivar prefix_cell: A non-standard, flow cell-specific (i.e., project_name) prefix.
    :type prefix_cell: str
    :ivar sample_annotation_sheet_csv: A Sample Annotation Sheet CSV file path.
    :type sample_annotation_sheet_csv: str
    """

    def __init__(self, prefix: str, project_name: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathExtractIlluminaCell` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        """
        super(FilePathExtractIlluminaCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the flow cell-specific prefix accordingly.

        self.prefix_cell = project_name

        self.sample_annotation_sheet_csv = self.prefix_cell + '_samples.csv'

        return


class FilePathExtractIlluminaLane(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathExtractIlluminaLane` class models lane-specific file paths.

    See also classes:
        - :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane`
        - :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane`

    :ivar prefix_lane: A non-standard, lane-specific (i.e., project_name and lane) prefix.
    :type prefix_lane: str
    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar samples_directory: A directory storing sample-specific unaligned BAM files.
    :type samples_directory: str
    :ivar barcode_tsv: A barcode TSV file for the Picard :literal:`ExtractIlluminaBarcodes` tool.
    :type barcode_tsv: str
    :ivar library_tsv: A library TSV file for the Picard :literal:`IlluminaBasecallsToSam` tool.
    :type library_tsv: str
    :ivar metrics_tsv: A metrics TSV file path.
    :type metrics_tsv: str
    :ivar metrics_fraction_pdf: A lane-specific Picard :literal:`ExtractIlluminaBarcodes`
        fraction metrics PDF file path.
    :type metrics_fraction_pdf: str
    :ivar metrics_fraction_png: A lane-specific Picard :literal:`ExtractIlluminaBarcodes`
        fraction metrics PNG file path.
    :type metrics_fraction_png: str
    :ivar metrics_number_pdf: A lane-specific Picard :literal:`ExtractIlluminaBarcodes`
        number metrics PDF file path.
    :type metrics_number_pdf: str
    :ivar metrics_number_png: A lane-specific Picard :literal:`ExtractIlluminaBarcodes`
        number metrics PNG file path.
    :type metrics_number_png: str
    """

    def __init__(self, prefix: str, project_name: str, lane: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathExtractIlluminaLane` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        """
        super(FilePathExtractIlluminaLane, self).__init__(prefix=prefix)

        # Non-standard project_name_lane prefix
        self.prefix_lane = '_'.join((project_name, lane))

        self.output_directory = self.prefix_lane + '_output'
        self.samples_directory = self.prefix_lane + '_samples'
        self.barcode_tsv = self.prefix_lane + '_barcode.tsv'
        self.library_tsv = self.prefix_lane + '_library.tsv'
        self.metrics_tsv = self.prefix_lane + '_metrics.tsv'
        self.metrics_fraction_pdf = self.prefix_lane + '_metrics_fraction.pdf'
        self.metrics_fraction_png = self.prefix_lane + '_metrics_fraction.png'
        self.metrics_number_pdf = self.prefix_lane + '_metrics_number.pdf'
        self.metrics_number_png = self.prefix_lane + '_metrics_number.png'

        return


class ExtractIlluminaRunFolder(PicardIlluminaRunFolder):
    """The :py:class:`bsf.analyses.picard.ExtractIlluminaRunFolder` class extracts data from an
    :emphasis:`Illumina Run Folder`.

    The analysis is based on Picard :literal:`ExtractIlluminaBarcodes` and Picard :literal:`IlluminaBasecallsToSam`.

    :ivar samples_directory: A directory storing sample-specific unaligned BAM files.
    :type samples_directory: str | None
    :ivar library_path: A library annotation file path.
    :type library_path: str | None
    :ivar mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
    :type mode_directory: str | None
    :ivar mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
    :type mode_file: str | None
    :ivar max_mismatches: A maximum number of mismatches.
    :type max_mismatches: int | None
    :ivar min_base_quality: A minimum base quality score.
    :type min_base_quality: int | None
    :ivar sequencing_centre: A sequencing centre code.
    :type sequencing_centre: str | None
    :ivar lanes: A number of lanes on a flow cell.
    :type lanes: int | None
    :ivar vendor_quality_filter: Request vendor quality filtering.
    :type vendor_quality_filter: bool | None
    """

    name = 'Picard Extract Illumina Run Folder Analysis'
    prefix = 'extract_illumina_run_folder'

    @classmethod
    def get_file_path_cell(cls, project_name: str) -> FilePathExtractIlluminaCell:
        """Get a :py:class:`bsf.analyses.picard.FilePathExtractIlluminaCell` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A :py:class:`bsf.analyses.picard.FilePathExtractIlluminaCell` object.
        :rtype: FilePathExtractIlluminaCell
        """
        return FilePathExtractIlluminaCell(
            prefix=cls.get_prefix_cell(project_name=project_name),
            project_name=project_name)

    @classmethod
    def get_file_path_lane(cls, project_name: str, lane: str) -> FilePathExtractIlluminaLane:
        """Get a :py:class:`bsf.analyses.picard.FilePathExtractIlluminaLane` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A :py:class:`bsf.analyses.picard.FilePathExtractIlluminaLane` object.
        :rtype: FilePathExtractIlluminaLane
        """
        return FilePathExtractIlluminaLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane)

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
            java_archive_picard: Optional[str] = None,
            force: Optional[bool] = None,
            library_path: Optional[str] = None,
            samples_directory: Optional[str] = None,
            mode_directory: Optional[str] = None,
            mode_file: Optional[str] = None,
            max_mismatches: Optional[int] = None,
            min_base_quality: Optional[int] = None,
            sequencing_centre: Optional[str] = None,
            lanes: Optional[int] = None,
            vendor_quality_filter: Optional[bool] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.ExtractIlluminaRunFolder` object.

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
        :param basecalls_directory: An :emphasis:`Illumina Run Folder` (IRF) :literal:`BaseCalls` directory path,
            defaults to :literal:`IRF/Data/Intensities/BaseCalls`.
        :type basecalls_directory: str | None
        :param experiment_name: Experiment name (i.e., flow cell identifier) normally automatically read from
            :emphasis:`Illumina Run Folder` parameters.
        :type experiment_name: str | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param force: Request processing of incomplete :emphasis:`Illumina Run Folder` instances.
        :type force: bool | None
        :param samples_directory: A directory storing sample-specific unaligned BAM files.
        :type samples_directory: str | None
        :param library_path: A library annotation file path.
        :type library_path: str | None
        :param mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
        :type mode_directory: str | None
        :param mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
        :type mode_file: str | None
        :param max_mismatches: A maximum number of mismatches.
        :type max_mismatches: int | None
        :param min_base_quality: A minimum base quality score.
        :type min_base_quality: int | None
        :param sequencing_centre: A sequencing centre code.
        :type sequencing_centre: str | None
        :param lanes: A number of lanes on the flow cell.
        :type lanes: int | None
        :param vendor_quality_filter: Request vendor quality filtering.
        :type vendor_quality_filter: bool | None
        """
        super(ExtractIlluminaRunFolder, self).__init__(
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
            sample_list=sample_list,
            run_directory=run_directory,
            intensity_directory=intensity_directory,
            basecalls_directory=basecalls_directory,
            experiment_name=experiment_name,
            java_archive_picard=java_archive_picard,
            force=force)

        self.samples_directory = samples_directory
        self.library_path = library_path
        self.mode_directory = mode_directory
        self.mode_file = mode_file
        self.max_mismatches = max_mismatches
        self.min_base_quality = min_base_quality
        self.sequencing_centre = sequencing_centre
        self.lanes = lanes
        self.vendor_quality_filter = vendor_quality_filter

        return

    @property
    def get_experiment_directory(self) -> Optional[str]:
        """Get an experiment directory.

        Experiment directory names are a concatenation of the sequences directory and the project name.

        :return: An experiment directory.
        :rtype: str | None
        """
        if self.samples_directory and self.project_name:
            return os.path.join(self.samples_directory, self.project_name)

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.picard.ExtractIlluminaRunFolder` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(ExtractIlluminaRunFolder, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'library_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_path = configuration.config_parser.get(section=section, option=option)

        option = 'samples_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.samples_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        option = 'max_mismatches'
        if configuration.config_parser.has_option(section=section, option=option):
            self.max_mismatches = configuration.config_parser.getint(section=section, option=option)

        option = 'min_base_quality'
        if configuration.config_parser.has_option(section=section, option=option):
            self.min_base_quality = configuration.config_parser.getint(section=section, option=option)

        option = 'sequencing_centre'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequencing_centre = configuration.config_parser.get(section=section, option=option)

        option = 'lanes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.lanes = configuration.config_parser.getint(section=section, option=option)

        option = 'vendor_quality_filter'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vendor_quality_filter = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.ExtractIlluminaRunFolder` object.
        """

        def run_get_sample_file_name(sample_name: str) -> str:
            """Private function to format sample-specific BAM file names (i.e., project_lane#sample.bam).

            :param sample_name: A sample name.
            :type sample_name: str
            :return: A sample-specific BAM file name.
            :rtype: str
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # Start of the run() method body.

        super(ExtractIlluminaRunFolder, self).run()

        self.samples_directory = self.configuration.get_absolute_path(
            file_path=self.samples_directory,
            default_path=StandardFilePath.get_samples(absolute=True))

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(f"The 'samples_directory' {self.samples_directory!r} does not exist.")

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get sequencing centre information.

        if not self.sequencing_centre:
            self.sequencing_centre = Operator.get_sequencing_centre()
            if not self.sequencing_centre:
                raise Exception(f"A {self.name!s} requires a 'sequencing_centre' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self.vendor_quality_filter is None:
            self.vendor_quality_filter = VendorQualityFilter.get_vendor_quality_filter(
                flow_cell_type=self._irf.run_parameters.get_flow_cell_type)

        # Get the library annotation sheet.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(self.library_path):
            raise Exception(f"The 'library_path' {self.library_path!r} is not a valid file.")

        # Load the LibraryAnnotationSheet and validate.

        library_annotation_sheet: LibraryAnnotationSheet = LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)

        if self.lanes is None:
            message_list = library_annotation_sheet.validate(
                lanes=self._irf.run_information.flow_cell_layout.lane_count)
        else:
            message_list = library_annotation_sheet.validate(
                lanes=self.lanes)

        if message_list:
            if self.force:
                warnings.warn(f'Validation of library annotation sheet {self.library_path!r}:\n'
                              '\n'.join(message_list), UserWarning)
            else:
                raise Exception(f'Validation of library annotation sheet {self.library_path!r}:\n'
                                '\n'.join(message_list))

        library_annotation_dict = library_annotation_sheet.get_annotation_dict()
        library_barcode_dict = library_annotation_sheet.get_barcode_length_dict()

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        file_path_cell = self.get_file_path_cell(project_name=self.project_name)

        # NOTE: Use a SampleAnnotationSheet as long as the Collection does not have a method to
        # write a de-normalised table.
        # Create a SampleAnnotationSheet in the project directory and
        # eventually transfer it into the experiment_directory.

        sample_annotation_sheet = BamIndexDecoder.get_sample_annotation_sheet(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        # For each lane in the library_annotation_dict ...
        # TODO: For the moment this depends on the lanes (keys) defined in the LibraryAnnotationSheet.
        # Not all lanes may thus get extracted.
        # TODO: For NextSeq instruments, it would be sufficient to require annotation for only lane one and
        # copy information to lanes two to four internally.

        cell_dependency_list: list[str] = list()

        for lane_int in sorted(library_annotation_dict):
            lane_str = str(lane_int)

            # Determine the read structure from the first line of the lane-specific library annotation or,
            # if that fails from the run configuration via the bsf.illumina.RunInformation object.

            if 'read_structure' in library_annotation_dict[lane_int][0] and \
                    library_annotation_dict[lane_int][0]['read_structure']:
                read_structure = library_annotation_dict[lane_int][0]['read_structure']
            else:
                read_structure = self._irf.run_information.get_picard_read_structure

            file_path_lane = self.get_file_path_lane(project_name=self.project_name, lane=lane_str)

            # BARCODE_FILE
            eib_sheet = ExtractIlluminaBarcodesSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.barcode_tsv))

            # LIBRARY_PARAMS
            ibs_sheet = IlluminaBasecallsToSamSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.library_tsv))

            if library_barcode_dict[lane_int][1] > 0:
                barcode_number = 2
            elif library_barcode_dict[lane_int][0] > 0:
                barcode_number = 1
            else:
                barcode_number = 0

            # Sort each lane by sample name.
            flow_cell_dict_list = library_annotation_dict[lane_int]
            flow_cell_dict_list.sort(key=lambda item: item['sample_name'])

            for row_dict in flow_cell_dict_list:
                # Add a row to the lane-specific Picard ExtractIlluminaBarcodesSheet.

                eib_sheet.row_dict_list.append({
                    'barcode_sequence_1': row_dict['barcode_sequence_1'],
                    'barcode_sequence_2': row_dict['barcode_sequence_2'],
                    'barcode_name': row_dict['sample_name'],
                    'library_name': row_dict['library_name'],
                })

                # Add a row to the lane-specific Picard IlluminaBasecallsToSamSheet.

                ibs_sheet.row_dict_list.append({
                    'BARCODE_1': row_dict['barcode_sequence_1'],
                    'BARCODE_2': row_dict['barcode_sequence_2'],
                    'OUTPUT': os.path.join(
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    'SAMPLE_ALIAS': row_dict['sample_name'],
                    'LIBRARY_NAME': row_dict['library_name'],
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

            # Adjust the IlluminaBaseCallsToSamSheet by adding an entry for unassigned reads and
            # constraining the columns to the number of index reads.
            ibs_sheet.adjust(
                barcode_length_tuple=library_barcode_dict[lane_int],
                unassigned_file_path=os.path.join(
                    file_path_lane.samples_directory,
                    run_get_sample_file_name(sample_name='0')))

            # Write the lane-specific Picard ExtractIlluminaBarcodesSheet and Picard IlluminaBasecallsToSamSheet.

            if barcode_number > 0:
                eib_sheet.to_file_path()

            ibs_sheet.to_file_path()

            # Create a Runnable and Executable for the lane stage.

            runnable_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)

            cell_dependency_list.append(executable_lane.name)

            # Create an output_directory in the project_directory.

            runnable_step = RunnableStepMakeDirectory(
                name='make_output_directory',
                directory_path=file_path_lane.output_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Create a samples_directory in the project_directory.

            runnable_step = RunnableStepMakeDirectory(
                name='make_samples_directory',
                directory_path=file_path_lane.samples_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Create a RunnableStep for Picard ExtractIlluminaBarcodes, only if index (barcode) reads are present.

            if barcode_number > 0:
                runnable_step = RunnableStepPicard(
                    name='picard_extract_illumina_barcodes',
                    java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='ExtractIlluminaBarcodes')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # Required Arguments
                # BARCODE [null]
                # BARCODE_FILE [null]
                runnable_step.add_picard_option(key='BARCODE_FILE', value=file_path_lane.barcode_tsv)
                # BASECALLS_DIR [null]
                runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
                # LANE [null]
                runnable_step.add_picard_option(key='LANE', value=lane_str)
                # READ_STRUCTURE [null]
                runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)

                # Optional Tool Arguments
                # COMPRESS_OUTPUTS [false]
                runnable_step.add_picard_option(key='COMPRESS_OUTPUTS', value='true')
                # DISTANCE_MODE [HAMMING]
                # INPUT_PARAMS_FILE [null]
                # MAX_MISMATCHES [1]
                if self.max_mismatches is not None:
                    # Maximum mismatches for a barcode to be considered a match. Default value: '1'.
                    runnable_step.add_picard_option(key='MAX_MISMATCHES', value=str(self.max_mismatches))
                # MAX_NO_CALLS [2]
                # Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.
                # METRICS_FILE [null]
                runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
                # MIN_MISMATCH_DELTA [1]
                # Minimum difference between number of mismatches in the best and
                # second best barcodes for a barcode to be considered a match.
                # MINIMUM_BASE_QUALITY [0]
                if self.min_base_quality is not None:
                    # Minimum base quality. Any barcode bases falling below this quality will be considered
                    # a mismatch even in the bases match. Default value: '0'.
                    runnable_step.add_picard_option(key='MINIMUM_BASE_QUALITY', value=str(self.min_base_quality))
                # MINIMUM_QUALITY [2]
                # The minimum quality (after transforming 0s to 1s) expected from reads.
                # If qualities are lower than this value, an error is thrown.The default of 2 is what the
                # Illumina specification describes as the minimum, but in practice the value has been observed lower.
                # NUM_PROCESSORS [1]
                runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
                # OUTPUT_DIR [null] Default to BASECALLS_DIR
                runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.output_directory)

                # Optional Common Arguments
                # COMPRESSION_LEVEL [5]
                # CREATE_INDEX [false]
                # CREATE_MD5_FILE [false]
                # MAX_RECORDS_IN_RAM [500000]
                # QUIET [false]
                # REFERENCE_SEQUENCE [null]
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_lane.temporary_directory_path(absolute=False))
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]
                # VALIDATION_STRINGENCY [STRICT]
                # VERBOSITY [INFO]

                # Plot the metrics file.

                runnable_step = RunnableStep(
                    name='plot_metrics',
                    program='bsf_illumina_demultiplex_sam.R')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='file-path', value=file_path_lane.metrics_tsv)

            # Picard IlluminaBasecallsToSam

            # Create a RunnableStep for Picard IlluminaBasecallsToSam.

            runnable_step = RunnableStepPicard(
                name='picard_illumina_basecalls_to_sam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx2G',
                java_jar_path=self.java_archive_picard,
                picard_command='IlluminaBasecallsToSam')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Required Arguments
            # BARCODE_PARAMS [null] deprecated in favour of LIBRARY_PARAMS
            # BASECALLS_DIR [null]
            runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
            # LANE [null]
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # LIBRARY_PARAMS [null]
            runnable_step.add_picard_option(key='LIBRARY_PARAMS', value=file_path_lane.library_tsv)
            # OUTPUT [null] deprecated in favour of LIBRARY_PARAMS
            # READ_STRUCTURE [null]
            runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)
            # RUN_BARCODE [null]
            runnable_step.add_picard_option(
                key='RUN_BARCODE',
                value=self._irf.run_parameters.get_flow_cell_barcode)
            # SAMPLE_ALIAS [null] deprecated in favour of LIBRARY_PARAMS
            # SEQUENCING_CENTER [null]
            runnable_step.add_picard_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)

            # Optional Tool Arguments
            # ADAPTERS_TO_CHECK [INDEXED, DUAL_INDEXED, NEXTERA_V2, FLUIDIGM]
            # APPLY_EAMSS_FILTER [true]
            # BARCODE_POPULATION_STRATEGY [ORPHANS_ONLY]
            # BARCODES_DIR [null] Defaults to BASECALLS_DIR
            if barcode_number > 0:
                runnable_step.add_picard_option(key='BARCODES_DIR', value=file_path_lane.output_directory)
            # COMPRESS_OUTPUTS [false]
            # DISTANCE_MODE [HAMMING]
            # FIRST_TILE [null]
            # FIVE_PRIME_ADAPTER [null]
            # IGNORE_UNEXPECTED_BARCODES [false]
            # INCLUDE_BARCODE_QUALITY [false]
            # INCLUDE_BC_IN_RG_TAG [false]
            # INCLUDE_NON_PF_READS [true]
            if self.vendor_quality_filter:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
            else:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
            # INPUT_PARAMS_FILE [null] the LIBRARY_PARAMS file for IlluminaBasecallsToSam
            # LIBRARY_NAME [null] deprecated in favour of LIBRARY_PARAMS
            # MATCH_BARCODES_INLINE [false]
            # MAX_MISMATCHES [1]
            # MAX_NO_CALLS [1]
            # MAX_READS_IN_RAM_PER_TILE [1200000]
            # METRICS_FILE [null]
            # MIN_MISMATCH_DELTA [1]
            # MINIMUM_BASE_QUALITY [0]
            # MINIMUM_QUALITY [1]
            # MOLECULAR_INDEX_BASE_QUALITY_TAG [QX]
            # MOLECULAR_INDEX_TAG [RX]
            # NUM_PROCESSORS [0]
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # PLATFORM [illumina]
            # PLATFORM The name of the sequencing technology that produced the read.
            # NOTE: IlluminaToBam defaults to 'ILLUMINA'.
            # runnable_step.add_picard_option(key='PLATFORM', value='ILLUMINA')
            # PROCESS_SINGLE_TILE [null]
            # READ_GROUP_ID [null]
            runnable_step.add_picard_option(
                key='READ_GROUP_ID',
                value='_'.join((self._irf.run_parameters.get_flow_cell_barcode, lane_str)))
            # RUN_START_DATE [null]
            # NOTE: The ISO date format still does not work for Picard tools 2.6.1. Sigh.
            # runnable_step.add_picard_option(key='RUN_START_DATE', value=self._irf.run_information.get_iso_date)
            # NOTE: The only date format that seems to work is MM/DD/YYYY. Why?
            runnable_step.add_picard_option(
                key='RUN_START_DATE',
                value='/'.join((
                    self._irf.run_information.date[2:4],
                    self._irf.run_information.date[4:6],
                    '20' + self._irf.run_information.date[0:2])))
            # SORT [true]
            # TAG_PER_MOLECULAR_INDEX [null]
            # THREE_PRIME_ADAPTER [null]
            # TILE_LIMIT [null]

            # Optional Common Arguments
            # COMPRESSION_LEVEL [5]
            runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # MAX_RECORDS_IN_RAM [500000]
            # QUIET [false]
            # REFERENCE_SEQUENCE [null]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # VALIDATION_STRINGENCY [STRICT]
            # VERBOSITY [INFO]

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

            if barcode_number > 0:
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

        # Finally, write the flow cell-specific SampleAnnotationSheet to the internal file path,
        # but keep the order of field names defined above.

        sample_annotation_sheet.to_file_path(adjust_field_names=False)

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
            mode: Optional[str] = None,
            irf_path: Optional[str] = None,
            force: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a :py:class:`bsf.analyses.picard.ExtractIlluminaRunFolder` analysis.

        This analysis requires either a :literal:`configuration` argument or a :literal:`irf_path` argument.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param mode: A sequencer mode.
        :type mode: str | None
        :param irf_path: An :emphasis:`Illumina Run Folder` path.
        :type irf_path: str | None
        :param force: Request processing of imperfect library annotation sheets.
        :type force: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if configuration_path == Configuration.get_global_file_path():
            if not irf_path:
                raise Exception('The --irf argument is required if configuration is not set.')

        analysis: ExtractIlluminaRunFolder = cls.from_config_file_path(config_path=configuration_path)

        if irf_path:
            analysis.run_directory = irf_path

        if mode:
            if mode == 'high':
                analysis.lanes = 8
            elif mode == 'rapid':
                analysis.lanes = 2
            elif mode == 'miseq':
                analysis.lanes = 1
            elif mode == 'nextseq':
                analysis.lanes = 4
            elif mode == 'novaseq':
                analysis.lanes = 4
            else:
                raise Exception(f'The mode option {mode!r} is not supported.')

        if force:
            analysis.force = force

        analysis.run()
        analysis.check_state()
        analysis.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis.name)
        print('Project name:           ', analysis.project_name)
        print('Project directory:      ', analysis.project_directory)
        print('Illumina run directory: ', analysis.run_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a :py:class:`bsf.analyses.picard.ExtractIlluminaRunFolder` analysis.

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
            '--mode',
            help='HiSeq run mode i.e., high (high-output) or rapid (rapid run)')

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
            mode=name_space.mode,
            irf_path=name_space.irf_path,
            force=name_space.force)


class FilePathIlluminaMultiplexSamLane(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathIlluminaMultiplexSamLane` class models files in a directory.

    :ivar unsorted_bam: An :literal:`unsorted` BAM file path.
    :type unsorted_bam: str
    :ivar unsorted_md5: An :literal:`unsorted` BAM file MD5 checksum file path.
    :type unsorted_md5: str
    :ivar sorted_bam: A :literal:`sorted` BAM file path.
    :type sorted_bam: str
    :ivar sorted_md5: A :literal:`sorted` BAM file MD5 checksum file path.
    :type sorted_md5: str
    :ivar archive_bam: An :literal:`archive` BAM file path.
    :type archive_bam: str
    :ivar archive_md5: An :literal:`archive` BAM file MD5 checksum file path.
    :type archive_md5: str
    """

    def __init__(self, prefix: str, project_name: str, lane: str, experiment_directory: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathIlluminaMultiplexSamLane` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param experiment_directory: An experiment-specific directory.
        :type experiment_directory: str
        """
        super(FilePathIlluminaMultiplexSamLane, self).__init__(prefix=prefix)

        self.unsorted_bam = prefix + '_unsorted.bam'
        self.unsorted_md5 = prefix + '_unsorted.bam.md5'
        self.sorted_bam = prefix + '_sorted.bam'
        self.sorted_md5 = prefix + '_sorted.bam.md5'
        # The final BAM and MD5 checksum files are non-standard in that they do not contain a prefix and
        # reside in the experiment_directory.
        self.final_bam = '_'.join((project_name, lane)) + '.bam'
        self.final_md5 = '_'.join((project_name, lane)) + '.bam.md5'
        self.archive_bam = os.path.join(experiment_directory, self.final_bam)
        self.archive_md5 = os.path.join(experiment_directory, self.final_md5)
        # The Azure Storage Blob Service always uses URL-compliant slash characters as path separators.
        self.cloud_bam = '/'.join((project_name, self.final_bam))
        self.cloud_md5 = '/'.join((project_name, self.final_md5))

        return


class IlluminaMultiplexSam(PicardIlluminaRunFolder):
    """The :py:class:`bsf.analyses.picard.IlluminaMultiplexSam` class represents the
    Picard :literal:`IlluminaBasecallsToMultiplexSam` analysis.

    :ivar sequencing_centre: A sequencing centre code.
    :type sequencing_centre: str | None
    :ivar sequences_directory: A directory storing lane-specific unaligned BAM files.
    :type sequences_directory: str | None
    :ivar mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
    :type mode_directory: str | None
    :ivar mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
    :type mode_file: str | None
    :ivar eamss_filter: Request Illumina EAMSS or Read Segment Quality Control Metric filtering.
    :type eamss_filter: bool | None
    :ivar vendor_quality_filter: Request vendor quality filtering.
    :type vendor_quality_filter: bool | None
    :ivar compression_level: A Zlib compression level.
    :type compression_level: int | None
    :ivar cloud_account: A :literal:`Microsoft Azure Storage Account` name.
    :type cloud_account: str | None
    :ivar cloud_container: A :literal:`Microsoft Azure Blob Service` container name.
    :type cloud_container: str | None
    :ivar cloud_concurrency: A maximum number of concurrent network connections.
    :type cloud_concurrency: int | None
    """

    name = 'Picard IlluminaMultiplexSam Analysis'
    prefix = 'illumina_multiplex_sam'

    @classmethod
    def get_stage_name_cloud(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cloud'))

    @classmethod
    def get_prefix_cloud(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` (prefix) object representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix) object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_cloud(), project_name, lane))

    @classmethod
    def get_file_path_lane(
            cls,
            project_name: str,
            lane: str,
            experiment_directory: str) -> FilePathIlluminaMultiplexSamLane:
        """Get a :py:class:`bsf.analyses.picard.FilePathIlluminaMultiplexSamLane` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param experiment_directory: An experiment directory.
        :type experiment_directory: str
        :return: A :py:class:`bsf.analyses.picard.FilePathIlluminaMultiplexSamLane` object.
        :rtype: FilePathIlluminaMultiplexSamLane
        """
        return FilePathIlluminaMultiplexSamLane(
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
            java_archive_picard: Optional[str] = None,
            force: Optional[bool] = None,
            sequencing_centre: Optional[str] = None,
            sequences_directory: Optional[str] = None,
            mode_directory: Optional[str] = None,
            mode_file: Optional[str] = None,
            eamss_filter: Optional[bool] = None,
            vendor_quality_filter: Optional[bool] = None,
            compression_level: Optional[int] = None,
            cloud_account: Optional[str] = None,
            cloud_container: Optional[str] = None,
            cloud_concurrency: Optional[int] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.IlluminaMultiplexSam` object.

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
        :param basecalls_directory: An :emphasis:`Illumina Run Folder` (IRF) :literal:`BaseCalls` directory path,
            defaults to :literal:`IRF/Data/Intensities/BaseCalls`.
        :type basecalls_directory: str | None
        :param experiment_name: An experiment name (i.e., flow cell identifier) normally automatically read from
            :emphasis:`Illumina Run Folder` parameters.
        :type experiment_name: str | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param force: Request processing of incomplete :emphasis:`Illumina Run Folder` instances.
        :type force: bool | None
        :param sequencing_centre: A sequencing centre code.
        :type sequencing_centre: str | None
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str | None
        :param mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory.
        :type mode_directory: str | None
        :param mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file.
        :type mode_file: str | None
        :param eamss_filter: Request Illumina EAMSS or Read Segment Quality Control Metric filtering.
        :type eamss_filter: bool | None
        :param vendor_quality_filter: Request vendor quality filtering.
        :type vendor_quality_filter: bool | None
        :param compression_level: A Zlib compression level.
        :type compression_level: int | None
        :param cloud_account: A :literal:`Microsoft Azure Storage Account` name.
        :type cloud_account: str | None
        :param cloud_container: A :literal:`Microsoft Azure Blob Service` container name.
        :type cloud_container: str | None
        :param cloud_concurrency: A maximum number of concurrent network connections.
        :type cloud_concurrency: int | None
        """
        super(IlluminaMultiplexSam, self).__init__(
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
            sample_list=sample_list,
            run_directory=run_directory,
            intensity_directory=intensity_directory,
            basecalls_directory=basecalls_directory,
            experiment_name=experiment_name,
            java_archive_picard=java_archive_picard,
            force=force)

        self.sequencing_centre = sequencing_centre
        self.sequences_directory = sequences_directory
        self.mode_directory = mode_directory
        self.mode_file = mode_file
        self.eamss_filter = eamss_filter
        self.vendor_quality_filter = vendor_quality_filter
        self.compression_level = compression_level
        self.cloud_account = cloud_account
        self.cloud_container = cloud_container
        self.cloud_concurrency = cloud_concurrency

        return

    @property
    def get_experiment_directory(self) -> Optional[str]:
        """Get the experiment directory.

        The experiment directory is a concatenation of the sequences directory and the project name.

        :return: An experiment directory.
        :rtype: str | None
        """
        if self.sequences_directory and self.project_name:
            return os.path.join(self.sequences_directory, self.project_name)

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.picard.IlluminaMultiplexSam` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(IlluminaMultiplexSam, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'sequencing_centre'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequencing_centre = configuration.config_parser.get(section=section, option=option)

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        option = 'eamss_filter'
        if configuration.config_parser.has_option(section=section, option=option):
            self.eamss_filter = configuration.config_parser.getboolean(section=section, option=option)

        option = 'vendor_quality_filter'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vendor_quality_filter = configuration.config_parser.getboolean(section=section, option=option)

        option = 'compression_level'
        if configuration.config_parser.has_option(section=section, option=option):
            self.compression_level = configuration.config_parser.getint(section=section, option=option)

        option = 'cloud_account'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_account = configuration.config_parser.get(section=section, option=option)

        option = 'cloud_container'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_container = configuration.config_parser.get(section=section, option=option)

        option = 'cloud_concurrency'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_concurrency = configuration.config_parser.getint(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.IlluminaMultiplexSam` object.
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(f"A {self.name!s} requires a 'run_directory' configuration option.")

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=StandardFilePath.get_illumina_run(absolute=True))

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(f"The 'run_directory' {self.run_directory!r} does not exist.")

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(f'The Illumina Run Folder {self.run_directory!r} is not complete.')

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
                raise Exception(f"A {self.name!s} requires a 'sequencing_centre' configuration option.")

        # Define the sequence directory in which to create the experiment directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
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

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self.vendor_quality_filter is None:
            self.vendor_quality_filter = VendorQualityFilter.get_vendor_quality_filter(
                flow_cell_type=irf.run_parameters.get_flow_cell_type)

        # Check that the cloud concurrency value (maximum number of network connections) is an integer.

        if not self.cloud_concurrency:
            self.cloud_concurrency = 1

        # Call the run method of the super class after the project_name has been defined.

        super(IlluminaMultiplexSam, self).run()

        cell_dependency_list = list()

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())
        stage_cloud = self.get_stage(name=self.get_stage_name_cloud())

        for lane_int in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            file_path_lane = self.get_file_path_lane(
                project_name=self.project_name,
                lane=lane_str,
                experiment_directory=experiment_directory)

            runnable_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))

            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)
            # executable_lane.dependencies.append(
            #     IlluminaRunFolderRestore.get_prefix_compress_base_calls(
            #         project_name=self.project_name,
            #         lane=lane_str))
            # Add the dependency for the cell-specific process.
            cell_dependency_list.append(executable_lane.name)

            # Only submit this Executable if the final result file does not exist.
            # absolute_sorted_md5 = os.path.join(experiment_directory, file_path_lane.sorted_md5)
            # if os.path.exists(absolute_sorted_md5) and os.path.getsize(absolute_sorted_md5):
            #     executable_lane.submit = False

            # Run Picard IlluminaBasecallsToUndemuxSam.

            runnable_step = RunnableStepPicard(
                name='picard_illumina_basecalls_to_multiplex_sam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx32G',
                java_jar_path=self.java_archive_picard,
                picard_command='IlluminaBasecallsToMultiplexSam',
                option_pair=True)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Required Arguments
            # LANE [null]
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # OUTPUT [null]
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.sorted_bam)
            # RUN_DIR [null]
            runnable_step.add_picard_option(key='RUN_DIR', value=self.run_directory)

            # APPLY_EAMSS_FILTER [true]
            if self.eamss_filter:
                runnable_step.add_picard_option(key='APPLY_EAMSS_FILTER', value='true')
            else:
                runnable_step.add_picard_option(key='APPLY_EAMSS_FILTER', value='false')
            # DEFLATER_THREADS [0]
            # NOTE: Only available in version: 2.18.24-CeMM
            runnable_step.add_picard_option(key='DEFLATER_THREADS', value='16')
            # FIRST_TILE [null]
            # FORCE_GC [true]
            # INCLUDE_NON_PF_READS [true]
            if self.vendor_quality_filter:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
            else:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
            # NUM_PROCESSORS [0]
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # MAX_READS_IN_RAM_PER_TILE [1200000]
            # MINIMUM_QUALITY [2]
            # SEQUENCING_CENTER [BI]
            runnable_step.add_picard_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)
            # TILE_LIMIT [null]

            # Optional Common Arguments
            # COMPRESSION_LEVEL [5]
            if self.compression_level is not None:
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value=str(self.compression_level))
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # MAX_RECORDS_IN_RAM [500000]
            # QUIET [false]
            # REFERENCE_SEQUENCE [null]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # VALIDATION_STRINGENCY [STRICT]
            # VERBOSITY [INFO]

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

            # Sleep for 60.0 seconds to allow for file system synchronisation.

            runnable_step = RunnableStepSleep(
                name='sleep',
                sleep_time=60.0)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Upload the archive BAM file into the block blob storage.

            if self.cloud_account and self.cloud_container:
                runnable_cloud = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_cloud(project_name=self.project_name, lane=lane_str),
                        working_directory=self.project_directory))
                executable_cloud = self.set_stage_runnable(stage=stage_cloud, runnable=runnable_cloud)
                executable_cloud.dependencies.append(executable_lane.name)

                # Upload the unaligned BAM file from its archive location.
                runnable_step = RunnableStepAzureBlockBlobUpload(
                    name='blob_upload_bam',
                    account_name=self.cloud_account,
                    container_name=self.cloud_container,
                    source_path=file_path_lane.archive_bam,
                    target_path=file_path_lane.cloud_bam,
                    max_concurrency=self.cloud_concurrency)
                runnable_cloud.add_runnable_step(runnable_step=runnable_step)

                # Upload the MD5 checksum file from its archive location.
                runnable_step = RunnableStepAzureBlockBlobUpload(
                    name='blob_upload_md5',
                    account_name=self.cloud_account,
                    container_name=self.cloud_container,
                    source_path=file_path_lane.archive_md5,
                    target_path=file_path_lane.cloud_md5,
                    max_concurrency=self.cloud_concurrency)
                runnable_cloud.add_runnable_step(runnable_step=runnable_step)

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
        """Console function to submit a :py:class:`bsf.analyses.picard.IlluminaMultiplexSam` analysis.

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

        analysis: IlluminaMultiplexSam = cls.from_config_file_path(config_path=configuration_path)

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
        """Console entry point to submit a :py:class:`bsf.analyses.picard.IlluminaMultiplexSam` analysis.

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


class FilePathIlluminaDemultiplexSamCell(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell` class models
    flow cell-specific file paths.

    See also classes:
        - :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderCell`
        - :py:class:`bsf.analyses.picard.FilePathExtractIlluminaCell`

    :ivar prefix_cell: A non-standard, flow cell-specific (i.e., project_name) prefix.
    :type prefix_cell: str
    :ivar sample_annotation_sheet_csv: A Sample Annotation Sheet CSV file path.
    :type sample_annotation_sheet_csv: str
    """

    def __init__(self, prefix: str, project_name: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        """
        super(FilePathIlluminaDemultiplexSamCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the flow cell-specific prefix accordingly.

        self.prefix_cell = project_name

        self.sample_annotation_sheet_csv = self.prefix_cell + '_samples.csv'

        return


class FilePathIlluminaDemultiplexSamLane(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane` class models lane-specific file paths.

    See also classes:
        - :py:class:`bsf.analyses.illumina_to_bam_tools.FilePathBamIndexDecoderLane`
        - :py:class:`bsf.analyses.picard.FilePathExtractIlluminaLane`

    :ivar prefix_lane: A non-standard, lane-specific (i.e., project_name and lane) prefix.
    :type prefix_lane: str
    :ivar project_barcode: A project-specific barcode CSV file.
    :type project_barcode: str
    :ivar library_tsv: A library TSV file for the Picard :literal:`IlluminaBamDemux` tool.
    :type library_tsv: str
    :ivar metrics_tsv: A lane-specific metrics TSV file.
    :type metrics_tsv: str
    :ivar metrics_fraction_pdf: A lane-specific Picard :literal:`IlluminaBamDemux` fraction metrics PDF file path.
    :type metrics_fraction_pdf: str
    :ivar metrics_fraction_png: A lane-specific Picard :literal:`IlluminaBamDemux` fraction metrics PNG file path.
    :type metrics_fraction_png: str
    :ivar metrics_number_pdf: A lane-specific Picard :literal:`IlluminaBamDemux` number metrics PDF file path.
    :type metrics_number_pdf: str
    :ivar metrics_number_png: A lane-specific Picard :literal:`IlluminaBamDemux` number metrics PNG file path.
    :type metrics_number_png: str
    :ivar samples_directory: A directory storing sample-specific unaligned BAM files.
    :type samples_directory: str
    :ivar archive_bam: An archive BAM file path.
    :type archive_bam: str
    :ivar archive_md5: An archive BAM MD5 checksum file path.
    :type archive_md5: str
    """

    def __init__(self, prefix: str, project_name: str, lane: str, sequences_directory: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str
        """
        super(FilePathIlluminaDemultiplexSamLane, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the lane-specific prefix accordingly.

        self.prefix_lane = '_'.join((project_name, lane))

        self.project_barcode = self.prefix_lane + '_barcode.tsv'
        # self.barcode_tsv = self.prefix_lane + '_barcode.tsv'
        self.library_tsv = self.prefix_lane + '_library.tsv'
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


class IlluminaDemultiplexSamSheet(AnnotationSheet):
    """The :py:class:`bsf.analyses.picard.IlluminaDemultiplexSamSheet` class represents a
    :emphasis:`Tab-Separated Value` (TSV) table of library information for the
    :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` object.
    """

    _file_type = 'excel-tab'

    # FIXME: This class is only required, because Picard IlluminaBamDemux uses SAMPLE_NAME rather than SAMPLE_ALIAS.
    # Once fixed, the IlluminaBasecallsToSamSheet class could be used.
    _field_name_list = [
        'OUTPUT',
        'SAMPLE_NAME',
        'LIBRARY_NAME',
        'BARCODE_1',
        'BARCODE_2',
    ]

    def adjust(self, barcode_length_tuple: tuple[int, int]) -> None:
        """Adjust a :py:class:`bsf.analyses.picard.IlluminaDemultiplexSamSheet` object.

        Remove :literal:`BARCODE_N` index columns, which annotation is all empty.

        :param barcode_length_tuple: A Python :py:class:`tuple` object of
            Python :py:class:`int` (barcode 1 length) and
            Python :py:class:`int` (barcode 2 length) objects.
        :type barcode_length_tuple: (int, int)
        """
        # Adjust the IlluminaDemultiplexSamSheet and remove any BARCODE_N columns not represented
        # in the barcode length list.
        for index in range(0, 1 + 1):
            if barcode_length_tuple[index] == 0:
                # Remove the 'BARCODE_N' field from the list of field names.
                if 'BARCODE_' + str(index + 1) in self.field_name_list:
                    self.field_name_list.remove('BARCODE_' + str(index + 1))
                # Remove the 'BARCODE_N' entry from each row_dict object, since csv.DictWriter requires it.
                for row_dict in self.row_dict_list:
                    row_dict.pop('BARCODE_' + str(index + 1), None)

        return


class IlluminaDemultiplexSam(Analysis):
    """The :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` class represents the logic to
    decode sequence archive BAM files into sample-specific BAM files via the Picard :literal:`IlluminaSamDemux` tool.

    :ivar library_path: A library annotation file path.
    :type library_path: str | None
    :ivar run_directory: An :emphasis:`Illumina Run Folder` (IRF) directory path.
    :type run_directory: str | None
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
    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    :ivar compression_level: A Zlib compression level.
    :type compression_level: int | None
    :ivar deflater_threads: A number of deflater threads.
    :type deflater_threads: int | None
    :ivar matching_threads: A number of matching threads.
    :type matching_threads: int | None
    :ivar lanes: A number of lanes on the flow cell.
    :type lanes: int | None
    :ivar lane_list: A Python :py:class:`list` object of Python :py:class:`str` (lane number) objects to process.
    :type lane_list: list[str] | None
    :ivar force: Request de-multiplexing with a Library Annotation sheet failing validation.
    :type force: bool | None
    """

    name = 'Picard IlluminaDemultiplexSam Analysis'
    prefix = 'illumina_demultiplex_sam'

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
    def get_prefix_cell(cls, project_name) -> str:
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_cell(), project_name))

    @classmethod
    def get_prefix_lane(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_lane(), project_name, lane))

    @classmethod
    def get_file_path_cell(cls, project_name: str) -> FilePathIlluminaDemultiplexSamCell:
        """Get a :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell` object.
        :rtype: FilePathIlluminaDemultiplexSamCell
        """
        return FilePathIlluminaDemultiplexSamCell(
            prefix=cls.get_prefix_cell(project_name=project_name),
            project_name=project_name)

    @classmethod
    def get_file_path_lane(
            cls,
            project_name: str,
            lane: str,
            sequences_directory: str) -> FilePathIlluminaDemultiplexSamLane:
        """Get a :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str
        :return: A :py:class:`bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane` object.
        :rtype: FilePathIlluminaDemultiplexSamLane
        """
        return FilePathIlluminaDemultiplexSamLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane,
            sequences_directory=sequences_directory)

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
            library_path: Optional[str] = None,
            run_directory: Optional[str] = None,
            sequences_directory: Optional[str] = None,
            samples_directory: Optional[str] = None,
            mode_directory: Optional[str] = None,
            mode_file: Optional[str] = None,
            java_archive_picard: Optional[str] = None,
            compression_level: Optional[int] = None,
            deflater_threads: Optional[int] = None,
            matching_threads: Optional[int] = None,
            lanes: Optional[int] = None,
            lane_list: Optional[list[str]] = None,
            force: Optional[bool] = None):
        """Initialise a :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` object.

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
        :param library_path: A library annotation file path.
        :type library_path: str | None
        :param run_directory: An :emphasis:`Illumina Run Folder` (IRF) directory path.
        :type run_directory: str | None
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
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param compression_level: A Zlib compression level.
        :type compression_level: int | None
        :param deflater_threads: A number of deflater threads.
        :type deflater_threads: int | None
        :param matching_threads: A number of matching threads.
        :type matching_threads: int | None
        :param lanes: A number of lanes on the flow cell.
        :type lanes: int | None
        :param lane_list: A Python :py:class:`list` object of Python :py:class:`str` (lane number) objects to process.
        :type lane_list: list[str] | None
        :param force: Request de-multiplexing with a Library Annotation sheet failing validation.
        :type force: bool | None
        """
        super(IlluminaDemultiplexSam, self).__init__(
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

        self.library_path = library_path
        # FIXME: Phase out the run_directory?
        # The run_directory is only required to instantiate a bsf.illumina.RunFolder object.
        # If this could be done at run time via the @RG PU field in the archive BAM file,
        # the run_folder option would no longer be required.
        self.run_directory = run_directory
        self.sequences_directory = sequences_directory
        self.samples_directory = samples_directory
        self.mode_directory = mode_directory
        self.mode_file = mode_file
        self.java_archive_picard = java_archive_picard
        self.compression_level = compression_level
        self.deflater_threads = deflater_threads
        self.matching_threads = matching_threads
        self.lanes = lanes
        self.lane_list = lane_list
        self.force = force

        return

    @property
    def get_experiment_directory(self) -> Optional[str]:
        """Get an experiment directory.

        An experiment directory is a concatenation of the samples directory and the project name.

        :return: An experiment directory.
        :rtype: str | None
        """
        if self.samples_directory and self.project_name:
            return os.path.join(self.samples_directory, self.project_name)

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(IlluminaDemultiplexSam, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'library_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_path = configuration.config_parser.get(section=section, option=option)

        option = 'illumina_run_folder'
        if configuration.config_parser.has_option(section=section, option=option):
            self.run_directory = configuration.config_parser.get(section=section, option=option)

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

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        option = 'compression_level'
        if configuration.config_parser.has_option(section=section, option=option):
            self.compression_level = configuration.config_parser.getint(section=section, option=option)

        option = 'deflater_threads'
        if configuration.config_parser.has_option(section=section, option=option):
            self.deflater_threads = configuration.config_parser.getint(section=section, option=option)

        option = 'matching_threads'
        if configuration.config_parser.has_option(section=section, option=option):
            self.matching_threads = configuration.config_parser.getint(section=section, option=option)

        option = 'lanes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.lanes = configuration.config_parser.getint(section=section, option=option)

        option = 'lane_list'
        if configuration.config_parser.has_option(section=section, option=option):
            self.lane_list = configuration.get_list_from_csv(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` object to
        decode an archive BAM file produced with IlluminaMultiplexSam into sample-specific BAM files.

        The standard BSF Python :literal:`comma-separated` value sample sheet needs to be transformed into
        a Picard tools :literal:`tab-separated` value (TSV) sample sheet.

        Columns:

        - lane
        - barcode_sequence_1
        - barcode_sequence_2
        - sample_name
        - library_name
        - barcode_sequence
        - barcode_name
        - library_name
        - sample_name
        - description
        """

        def run_get_sample_file_name(sample_name: str) -> str:
            """Private function to format sample-specific BAM file names (i.e., project_lane#sample.bam).

            :param sample_name: A sample name.
            :type sample_name: str
            :return: A sample-specific BAM file name.
            :rtype: str
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # Start of the run() method body.

        # Configuration variable to use a bsf.ngs.Collection object rather than a
        # bsf.annotation.AnnotationSheet object to write the flow cell-specific
        # sample annotation sheet.

        use_collection = False

        # Define an Illumina Run Folder directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(f"A {self.name!s} requires a 'run_directory' configuration option.")

        # Illumina Run Folders are kept either in the 'active' or the 'archive' location.
        # Upon re-running of the de-multiplexing stage, the IRF may have been archived already.

        self.run_directory = get_irf_path(name=self.run_directory)

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(f"The 'run_directory' {self.run_directory!r} does not exist.")

        irf: RunFolder = RunFolder.from_file_path(file_path=self.run_directory)

        if not self.project_name:
            self.project_name = '_'.join((irf.run_parameters.get_experiment_name, irf.run_information.flow_cell))

        # Only call the super-class method once the project name is defined with certainty.

        super(IlluminaDemultiplexSam, self).run()

        # Define the sequences and samples directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
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
            raise Exception(f"The 'library_path' {self.library_path!r} does not exist.")

        # Load the LibraryAnnotationSheet and validate.

        library_annotation_sheet: LibraryAnnotationSheet = LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)

        if self.lanes is None:
            message_list = library_annotation_sheet.validate(
                lanes=irf.run_information.flow_cell_layout.lane_count)
        else:
            message_list = library_annotation_sheet.validate(
                lanes=self.lanes)

        if message_list:
            if self.force:
                warnings.warn(f'Validation of library annotation sheet {self.library_path!r}:\n'
                              '\n'.join(message_list), UserWarning)
            else:
                raise Exception(f'Validation of library annotation sheet {self.library_path!r}:\n'
                                '\n'.join(message_list))

        library_annotation_dict = library_annotation_sheet.get_annotation_dict()
        library_barcode_dict = library_annotation_sheet.get_barcode_length_dict()

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        file_path_cell = self.get_file_path_cell(project_name=self.project_name)

        # NOTE: Use a SampleAnnotationSheet as long as the Collection does not have a method to
        # write a de-normalised table.
        # Create a SampleAnnotationSheet in the project directory and
        # eventually transfer it into the experiment_directory.

        sample_annotation_sheet = BamIndexDecoder.get_sample_annotation_sheet(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        # Create a NGS Collection and (default) NGS ProcessedRunFolder object to write a SampleAnnotationSheet.

        collection = Collection(
            name='Samples',
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        prf = ProcessedRunFolder(name=ProcessedRunFolder.default_name)
        collection.add_processed_run_folder(prf=prf)

        cell_dependency_list: list[str] = list()

        for lane_int in sorted(library_annotation_dict):
            lane_str = str(lane_int)

            # Determine the read structure from the first line of the lane-specific library annotation or,
            # if that fails from the run configuration via the bsf.illumina.RunInformation object.

            if 'read_structure' in library_annotation_dict[lane_int][0] and \
                    library_annotation_dict[lane_int][0]['read_structure']:
                read_structure = library_annotation_dict[lane_int][0]['read_structure']
            else:
                read_structure = irf.run_information.get_picard_read_structure

            file_path_lane = self.get_file_path_lane(
                project_name=self.project_name,
                lane=lane_str,
                sequences_directory=self.sequences_directory)

            # LIBRARY_PARAMS
            ibs_sheet = IlluminaDemultiplexSamSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.library_tsv))

            for row_dict in library_annotation_dict[lane_int]:
                # Add a row to the lane-specific Picard IlluminaDemultiplexSamSheet.

                ibs_sheet.row_dict_list.append({
                    'OUTPUT': os.path.join(
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    'SAMPLE_NAME': row_dict['sample_name'],
                    'LIBRARY_NAME': row_dict['library_name'],
                    'BARCODE_1': row_dict['barcode_sequence_1'],
                    'BARCODE_2': row_dict['barcode_sequence_2'],
                })

                # Get or create a (default) NGS Project.

                if Project.default_name in prf.project_dict:
                    project = prf.project_dict[Project.default_name]
                else:
                    project = Project(name=Project.default_name)
                    prf.add_project(project=project)

                # Get or create a NGS Sample.

                if row_dict['sample_name'] in project.sample_dict:
                    sample = project.sample_dict[row_dict['sample_name']]
                else:
                    sample = Sample(name=row_dict['sample_name'])
                    project.add_sample(sample=sample)

                # Create one NGS Reads object per row.

                reads_1 = Reads(
                    name='_'.join((self.project_name, lane_str, row_dict['sample_name'])),
                    file_path=os.path.join(
                        os.path.basename(experiment_directory),
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])))

                # Create one NGS PairedReads object per NGS Reads object.

                paired_reads = PairedReads(
                    reads_1=reads_1,
                    index_1=row_dict['barcode_sequence_1'],
                    index_2=row_dict['barcode_sequence_2'])
                paired_reads.add_annotation(key='Flow Cell', value=self.project_name)
                paired_reads.add_annotation(key='Lane', value=lane_str)
                paired_reads.add_annotation(key='Flow Cell Lane', value='_'.join((self.project_name, lane_str)))
                paired_reads.add_annotation(key='Project Name', value=row_dict['library_name'])
                paired_reads.add_annotation(key='Project Size', value=row_dict['library_size'])
                paired_reads.add_annotation(key='Structure', value=read_structure)

                # Finally, add the NGS PairedReads object to the NGS Sample object.

                sample.add_paired_reads(paired_reads=paired_reads)

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

            # Adjust the IlluminaDemultiplexSamSheet by constraining the columns to the number of index reads.

            ibs_sheet.adjust(barcode_length_tuple=library_barcode_dict[lane_int])

            # Create a Runnable and Executable for the lane stage.

            runnable_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    working_directory=self.project_directory))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)
            executable_lane.dependencies.append(
                IlluminaMultiplexSam.get_prefix_lane(
                    project_name=self.project_name,
                    lane=lane_str))

            # If a list of lanes was defined, ignore those lanes not on the list.
            if self.lane_list and lane_str not in self.lane_list:
                executable_lane.submit = False

            if executable_lane.submit:
                # Only if this Executable actually gets submitted, because a status file indicating completion
                # is not available or a list of lanes was defined and this one i snot on it ...

                # Set a dependency for the cell stage.
                cell_dependency_list.append(executable_lane.name)

                # Write the lane-specific BamIndexDecoderSheet to the internal file path.
                ibs_sheet.to_file_path()

            # Create a samples_directory in the project_directory.

            runnable_step = RunnableStepMakeDirectory(
                name='make_project_lane_samples_directory',
                directory_path=file_path_lane.samples_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Run Picard IlluminaSamDemux.

            runnable_step = RunnableStepPicard(
                name='picard_illumina_demultiplex_sam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx5G',
                java_jar_path=self.java_archive_picard,
                picard_command='IlluminaSamDemux',
                option_pair=True)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Set htsjdk system properties at the level of the JavaVM process.
            # NOTE: Version 2.18.24-CeMM has a DEFLATER_THREADS option for parallel compression.
            # runnable_step.add_option_pair_short(key='Dsamjdk.use_async_io_read_samtools', value='TRUE')
            # runnable_step.add_option_pair_short(key='Dsamjdk.use_async_io_write_samtools', value='TRUE')

            # Required Arguments
            # INPUT [null]
            runnable_step.add_picard_option(key='INPUT', value=file_path_lane.archive_bam)
            # LIBRARY_PARAMS [null]
            runnable_step.add_picard_option(key='LIBRARY_PARAMS', value=file_path_lane.library_tsv)
            # METRICS_FILE [null]
            runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
            # OUTPUT_DIR [null] Defaults to the same directory as INPUT.
            runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.samples_directory)

            # Optional Tool Arguments
            # ADAPTERS_TO_CHECK [INDEXED, DUAL_INDEXED, NEXTERA_V2, FLUIDIGM]
            # BARCODE_TAG_NAME [BC]
            # BARCODE_QUALITY_TAG_NAME [QT]
            # CELL_INDEX_TAG [CB]
            # CELL_INDEX_BASE_QUALITY_TAG [CY]
            # DEFLATER_THREADS [0]
            # NOTE: Only available in version: 2.18.24-CeMM
            if self.deflater_threads:
                runnable_step.add_picard_option(key='DEFLATER_THREADS', value=str(self.deflater_threads))
            # FIVE_PRIME_ADAPTER [null]
            # MATCHING_THREADS [1]
            # NOTE: Only available in version: 2.18.24-CeMM
            if self.matching_threads:
                runnable_step.add_picard_option(key='MATCHING_THREADS', value=str(self.matching_threads))
            # MAX_MISMATCHES [1]
            if 'barcode_mismatches' in library_annotation_dict[lane_int][0] and \
                    library_annotation_dict[lane_int][0]['barcode_mismatches']:
                # If the barcode_mismatches field exists and has a meaningful value ...
                runnable_step.add_picard_option(
                    key='MAX_MISMATCHES',
                    value=library_annotation_dict[lane_int][0]['barcode_mismatches'])
            elif library_barcode_dict[lane_int][1] > 0:
                # If a second barcode is defined ...
                runnable_step.add_picard_option(key='MAX_MISMATCHES', value='2')
            # MAX_NO_CALLS [2]
            # MIN_MISMATCH_DELTA [1]
            # MINIMUM_BASE_QUALITY [0]
            # MOLECULAR_INDEX_BASE_QUALITY_TAG [QX]
            # MOLECULAR_INDEX_TAG [RX]
            # OUTPUT_FORMAT [bam]
            # OUTPUT_PREFIX [null] Defaults to the @RG ID attribute of the multiplexed BAM file
            runnable_step.add_picard_option(key='OUTPUT_PREFIX', value='_'.join((self.project_name, lane_str)))
            # READ_STRUCTURE [null]
            runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)
            # TAG_PER_MOLECULAR_INDEX [null]
            # THREE_PRIME_ADAPTER [null]

            # Optional Common Arguments
            # COMPRESSION_LEVEL [5]
            if self.compression_level is not None:
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value=str(self.compression_level))
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # MAX_RECORDS_IN_RAM [500000]
            # QUIET [false]
            # REFERENCE_SEQUENCE [null]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # VALIDATION_STRINGENCY [STRICT]
            # VERBOSITY [INFO]

            # Plot the metrics file.

            runnable_step = RunnableStep(
                name='plot_metrics',
                program='bsf_illumina_demultiplex_sam.R')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_long(key='file-path', value=file_path_lane.metrics_tsv)

            # Create the experiment directory if it does not exist already.

            runnable_step = RunnableStepMakeDirectory(
                name='make_experiment_directory',
                directory_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move the sample directory into the experiment directory.

            runnable_step = RunnableStepMove(
                name='move_project_lane_samples_directory',
                source_path=file_path_lane.samples_directory,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move the metrics file and the plots into the experiment directory.

            runnable_step = RunnableStepMove(
                name='move_metrics_tsv',
                source_path=file_path_lane.metrics_tsv,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = RunnableStepMove(
                name='move_metrics_fraction_pdf',
                source_path=file_path_lane.metrics_fraction_pdf,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = RunnableStepMove(
                name='move_metrics_fraction_png',
                source_path=file_path_lane.metrics_fraction_png,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = RunnableStepMove(
                name='move_metrics_number_pdf',
                source_path=file_path_lane.metrics_number_pdf,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = RunnableStepMove(
                name='move_metrics_number_png',
                source_path=file_path_lane.metrics_number_png,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

        # Finally, write the flow cell-specific Collection to the internal file path.

        if use_collection:
            collection.to_sas_path(
                file_path=os.path.join(
                    self.project_directory,
                    file_path_cell.sample_annotation_sheet_csv),
                name='Samples')
        else:
            # Write the SampleAnnotationSheet to its internal file path,
            # but keep the order of field names defined above.

            sample_annotation_sheet.to_file_path(adjust_field_names=False)

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

        if os.path.exists(collection.file_path):
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
    def console_reset(
            cls,
            experiment_name: str,
            lane_numbers: str,
            dry_run: bool = False) -> int:
        """Console function to reset the demultiplexing so that the
        :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam`
        analysis can be re-run with an altered library annotation sheet.

        - Find the project folder and remove:

          - :literal:`illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}_completed.txt`
          - :literal:`illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}.pkl`
          - If lane-specific files were scheduled for removal, also remove the cell-specific status files.

            - :literal:`illumina_demultiplex_sam_cell_{experiment_name}_completed.txt`
            - :literal:`illumina_demultiplex_sam_cell_{experiment_name}.pkl`

        - Find the project status folder and remove:

          - :literal:`illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}_{SLURM_JOB_ID}.err`
          - :literal:`illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}_{SLURM_JOB_ID}.out`
          - If lane-specific files were scheduled for removal, also remove the cell-specific status files.

            - :literal:`illumina_demultiplex_sam_cell_{experiment_name}_{SLURM_JOB_ID}.err`
            - :literal:`illumina_demultiplex_sam_cell_{experiment_name}_{SLURM_JOB_ID}.out`

        - Find the (de-multiplexed) sample folder and remove:

          - :literal:`{experiment_name}_{lane_number}_samples`
          - :literal:`{experiment_name}_{lane_number}_fraction.pdf`
          - :literal:`{experiment_name}_{lane_number}_fraction.png`
          - :literal:`{experiment_name}_{lane_number}_numbers.pdf`
          - :literal:`{experiment_name}_{lane_number}_numbers.png`
          - If any of the above get removed, also remove:

            - :literal:`{experiment_name}_samples.csv`

        :param experiment_name: An experiment-name.
        :type experiment_name: str
        :param lane_numbers: A comma-separated list of lane numbers.
        :type lane_numbers: str
        :param dry_run: Request a dry run.
        :type dry_run: bool
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """

        def _change_file_mode_user_writable(local_file_path):
            """Private function to change a file mode by adding :py:const:`stat.S_IWUSR`.

            :param local_file_path: A file path.
            :type local_file_path: str
            """
            try:
                file_path_result = os.stat(local_file_path)
            except FileExistsError:
                module_logger.debug('File path %r does not exist.', local_file_path)
                return

            if not file_path_result.st_mode & stat.S_IWUSR:
                module_logger.debug(
                    'Resetting mode %0o to %0o for %r',
                    file_path_result.st_mode,
                    file_path_result.st_mode | stat.S_IWUSR,
                    local_file_path)

                if not dry_run:
                    os.chmod(local_file_path, file_path_result.st_mode | stat.S_IWUSR)

            return

        lane_number_list = lane_numbers.split(',')

        # 1. Identify the project folder first. If it still exists, this is a recent project.
        project_path = os.path.join(StandardFilePath.get_projects(absolute=True), experiment_name)

        if os.path.isdir(project_path):
            module_logger.debug('Project folder %r exists.', project_path)

            # Find
            # '{experiment_name}_{lane_number}_samples' directories.
            directory_path_list = []
            for lane_number in lane_number_list:
                directory_path = os.path.join(project_path, f'{experiment_name}_{lane_number}_samples')
                if os.path.isdir(directory_path):
                    directory_path_list.append(directory_path)

            if directory_path_list:
                # List all directory paths to be removed.
                print(f'Removing {experiment_name} demultiplexing directories:')
                for directory_path in directory_path_list:
                    print('  ' + os.path.basename(directory_path))

                result = input('Proceed? [y/n] ')
                if result == 'y':
                    for directory_path in directory_path_list:
                        module_logger.debug('Removing %r', directory_path)
                        if not dry_run:
                            shutil.rmtree(directory_path)

            # Find and remove:
            # - illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}_completed.txt
            # - illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}.pkl
            file_path_list = []
            for lane_number in lane_number_list:
                prefix = IlluminaDemultiplexSam.get_prefix_lane(project_name=experiment_name, lane=str(lane_number))

                file_path = os.path.join(project_path, prefix + '_completed.txt')
                if os.path.isfile(file_path):
                    file_path_list.append(file_path)

                file_path = os.path.join(project_path, prefix + '.pkl')
                if os.path.isfile(file_path):
                    file_path_list.append(file_path)

            if file_path_list:
                # If lane-specific files were scheduled for removal, find and remove the cell-specific status files.
                # - illumina_demultiplex_sam_cell_{experiment_name}_completed.txt
                # - illumina_demultiplex_sam_cell_{experiment_name}.pkl
                prefix = IlluminaDemultiplexSam.get_prefix_cell(project_name=experiment_name)
                file_path = os.path.join(project_path, prefix + '_completed.txt')
                if os.path.isfile(file_path):
                    file_path_list.append(file_path)

                file_path = os.path.join(project_path, prefix + '.pkl')
                if os.path.isfile(file_path):
                    file_path_list.append(file_path)

                # List all file paths to be removed.
                print(f'Removing {experiment_name} project status files:')
                for file_path in file_path_list:
                    print('  ' + os.path.basename(file_path))

                result = input('Proceed? [y/n] ')
                if result == 'y':
                    for file_path in file_path_list:
                        module_logger.debug('Removing %r', file_path)
                        if not dry_run:
                            os.remove(file_path)

            output_path = os.path.join(project_path, 'bsfpython_slurm_output')
            if os.path.isdir(output_path):
                module_logger.debug('Project output folder %r exists.', output_path)

                # Find and remove
                # - illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}_{SLURM_JOB_ID}.err
                # - illumina_demultiplex_sam_lane_{experiment_name}_{lane_number}_{SLURM_JOB_ID}.out
                file_path_list = []
                file_name_list = os.listdir(output_path)
                for lane_number in lane_number_list:
                    prefix = IlluminaDemultiplexSam.get_prefix_lane(project_name=experiment_name, lane=str(lane_number))
                    for file_name in file_name_list:
                        if file_name.startswith(prefix):
                            file_path_list.append(os.path.join(output_path, file_name))

                if file_path_list:
                    # If lane-specific files were scheduled for removal, find and remove the cell-specific status files.
                    # - illumina_demultiplex_sam_cell_{experiment_name}_{SLURM_JOB_ID}.err
                    # - illumina_demultiplex_sam_cell_{experiment_name}_{SLURM_JOB_ID}.out
                    prefix = IlluminaDemultiplexSam.get_prefix_cell(project_name=experiment_name)
                    for file_name in file_name_list:
                        if file_name.startswith(prefix):
                            file_path_list.append(os.path.join(output_path, file_name))

                    print(f'Removing {experiment_name} project output files:')
                    for file_path in file_path_list:
                        print('  ' + os.path.basename(file_path))

                    result = input('Proceed? [y/n] ')
                    if result == 'y':
                        for file_path in file_path_list:
                            module_logger.debug('Removing %r', file_path)
                            if not dry_run:
                                os.remove(file_path)
        else:
            print(f'This is not a recent project since folder {experiment_name} no longer exists.')
            result = input('Proceed? [y/n] ')
            if result != 'y':
                return 0

        # 2. Identify the sample folder.
        sample_path = os.path.join(StandardFilePath.get_samples(absolute=True), experiment_name)

        if os.path.isdir(sample_path):
            module_logger.debug('Sample folder %r exists.', sample_path)
            file_path_list = []
            file_name_list = os.listdir(sample_path)
            for lane_number in lane_number_list:
                prefix = '_'.join((experiment_name, lane_number))
                for file_name in file_name_list:
                    if file_name.startswith(prefix):
                        file_path_list.append(os.path.join(sample_path, file_name))

            if file_path_list:
                # If lane-specific files were scheduled for removal, find and remove the sample annotation table.
                # - {experiment_name}_samples.csv
                file_path = os.path.join(sample_path, experiment_name + '_samples.csv')
                if os.path.isfile(file_path):
                    file_path_list.append(file_path)

                print(f'Removing {experiment_name} sample files and folders:')
                for file_path in file_path_list:
                    print('  ' + os.path.basename(file_path))

                result = input('Proceed? [y/n] ')
                if result == 'y':
                    for directory_path, directory_name_list, file_name_list in os.walk(sample_path):
                        for directory_name in directory_name_list:
                            module_logger.debug('Changing directory %r -> %r', directory_path, directory_name)
                            _change_file_mode_user_writable(
                                local_file_path=os.path.join(directory_path, directory_name))

                        for file_name in file_name_list:
                            module_logger.debug('Changing file %r -> %r', directory_path, file_name)
                            _change_file_mode_user_writable(local_file_path=os.path.join(directory_path, file_name))

                    # Reset also the top-level sample directory.
                    _change_file_mode_user_writable(local_file_path=sample_path)

                    for file_path in file_path_list:
                        module_logger.debug('Removing %r', file_path)
                        if not dry_run:
                            if os.path.isdir(file_path):
                                shutil.rmtree(file_path)
                            if os.path.isfile(file_path):
                                os.remove(file_path)

        return 0

    @classmethod
    def entry_point_reset(cls) -> int:
        """Console entry point to reset the demultiplexing so that the
        :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam`
        analysis can be re-run with an altered library annotation sheet.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description='IlluminaDemultiplexSam reset script.')

        argument_parser.add_argument(
            '--logging-level',
            default='WARNING',
            choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
            help='logging level [WARNING]')

        argument_parser.add_argument(
            '--experiment-name',
            required=True,
            help='experiment name')

        argument_parser.add_argument(
            '--lane-numbers',
            required=True,
            help='lane number(s)')

        argument_parser.add_argument(
            '--dry-run',
            action='store_true',
            help='dry run')

        name_space = argument_parser.parse_args()

        if name_space.logging_level:
            logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
            logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

            logging.basicConfig(level=name_space.logging_level)

        return cls.console_reset(
            experiment_name=name_space.experiment_name,
            lane_numbers=name_space.lane_numbers,
            dry_run=name_space.dry_run)

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            project_name: Optional[str] = None,
            irf_path: Optional[str] = None,
            mode: Optional[str] = None,
            sequences_directory: Optional[str] = None,
            lane_list: Optional[str] = None,
            library_path: Optional[str] = None,
            force: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` analysis.

        This analysis requires either a :literal:`configuration` argument or
        a :literal:`irf_path` argument. The :literal:`project_name` argument can be automatically determined from the
        :literal:`irf_path` argument, while the
        :literal:`library_path` argument can be automatically determined from the
        :literal:`project_name` argument.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param project_name: A project name.
        :type project_name: str | None
        :param irf_path: An :emphasis:`Illumina Run Folder` path.
        :type irf_path: str | None
        :param sequences_directory: A directory storing lane-specific unaligned BAM files.
        :type sequences_directory: str | None
        :param lane_list: A comma-separated list of lane numbers.
        :type lane_list: str
        :param mode: A sequencer mode.
        :type mode: str | None
        :param library_path: A library path.
        :type library_path: str | None
        :param force: Request processing of imperfect library annotation sheets.
        :type force: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if configuration_path == Configuration.get_global_file_path():
            # The --project-name value can be determined from the --irf option, while the
            # --library-path option can be determined from the --project-name option.
            if not (project_name or irf_path):
                raise Exception('The --project-name or --irf arguments are required if configuration is not set.')

        analysis: IlluminaDemultiplexSam = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (INI) file.')

            analysis.project_name = project_name

        if irf_path:
            analysis.run_directory = irf_path

        if sequences_directory:
            analysis.sequences_directory = sequences_directory

        if lane_list:
            analysis.lane_list = lane_list.split(',')

        # FIXME: Remove the mode option.
        # Since the script requires an IRF option, the lanes can be read from the bsf.illumina.RunFolder.
        # TODO: Synchronise to other de-multiplexing modules?
        if mode:
            if mode == 'high':
                analysis.lanes = 8
            elif mode == 'rapid':
                analysis.lanes = 2
            elif mode == 'miseq':
                analysis.lanes = 1
            elif mode == 'nextseq':
                analysis.lanes = 4
            elif mode == 'novaseq':
                analysis.lanes = 4
            else:
                raise Exception(f'The --mode option {mode!r} is not supported.')

        if force:
            analysis.force = force

        if library_path:
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
        """Console entry point to submit a :py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` analysis.

        This analysis requires either a positional :literal:`configuration` argument or
        a :literal:`--irf` argument. The :literal:`--project-name` argument can be automatically determined from the
        :literal:`--irf` argument, while the
        :literal:`--library-path` argument can be automatically determined from the
        :literal:`--project-name` argument.

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

        # FIXME: For the moment, the IRF needs passing in.
        argument_parser.add_argument(
            '--irf',
            help='Illumina Run Folder name or file path',
            dest='irf_path')

        argument_parser.add_argument(
            '--project-name',
            help='project name (i.e., flow cell identifier)')

        argument_parser.add_argument(
            '--sequences-path',
            help='sequence archive directory path')

        argument_parser.add_argument(
            '--library-path',
            help='library annotation sheet (CSV) file path')

        argument_parser.add_argument(
            '--lane-list',
            help='comma-separated list of lanes to process')

        argument_parser.add_argument(
            '--mode',
            help='HiSeq run mode (i.e., "high" (high-output) or "rapid" (rapid run))')

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
            irf_path=name_space.irf_path,
            sequences_directory=name_space.sequences_path,
            mode=name_space.mode,
            lane_list=name_space.lane_list,
            library_path=name_space.library_path,
            force=name_space.force)


class FilePathCollectHiSeqXPfFailMetricsLane(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane` class models files in a directory.

    :ivar summary_tsv: A summary metrics TSV file.
    :type summary_tsv: str
    :ivar detailed_tsv: A detailed metrics TSV file.
    :type detailed_tsv: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathCollectHiSeqXPfFailMetricsLane, self).__init__(prefix=prefix)

        self.summary_tsv = prefix + '.pffail_summary_metrics'
        self.detailed_tsv = prefix + '.pffail_detailed_metrics'

        return


class CollectHiSeqXPfFailMetrics(PicardIlluminaRunFolder):
    """The :py:class:`bsf.analyses.picard.CollectHiSeqXPfFailMetrics` class represents the logic to run the
    Picard :literal:`CollectHiSeqXPfFailMetrics` tool.
    """

    name = 'Picard CollectHiSeqXPfFailMetrics Analysis'
    prefix = 'picard_hiseq_x_pf_fail'

    @classmethod
    def get_file_path_lane(cls, project_name: str, lane: str) -> FilePathCollectHiSeqXPfFailMetricsLane:
        """Get a :py:class:`bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A :py:class:`bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane` object.
        :rtype: FilePathCollectHiSeqXPfFailMetricsLane
        """
        return FilePathCollectHiSeqXPfFailMetricsLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane))

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.CollectHiSeqXPfFailMetrics` object.
        """
        super(CollectHiSeqXPfFailMetrics, self).run()

        # Picard CollectHiSeqXPfFailMetrics

        stage_lane = self.get_stage(name=self.get_stage_name_lane())

        cell_dependency_list = list()

        for lane_int in range(0 + 1, self._irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            file_path_lane = self.get_file_path_lane(project_name=self.project_name, lane=lane_str)

            runnable_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_lane(
                        project_name=self.project_name,
                        lane=lane_str),
                    working_directory=self.project_directory))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)

            # Add the dependency for the cell-specific process.

            cell_dependency_list.append(executable_lane.name)

            runnable_step = RunnableStepPicard(
                name='collect_hiseq_x_fail_metrics',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx4G',
                java_jar_path=self.java_archive_picard,
                picard_command='CollectHiSeqXPfFailMetrics')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Required Arguments
            # BASECALLS_DIR [null]
            runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
            # OUTPUT [null]
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.prefix)
            # LANE [null]
            runnable_step.add_picard_option(key='LANE', value=lane_str)

            # Optional Tool Arguments
            # N_CYCLES [24] Should match Illumina RTA software.
            # NUM_PROCESSORS [1]
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # PROB_EXPLICIT_READS [0.0]

            # Optional Common Arguments
            # COMPRESSION_LEVEL [5]
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            # MAX_RECORDS_IN_RAM [500000]
            # QUIET [false]
            # REFERENCE_SEQUENCE [null]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # VALIDATION_STRINGENCY [STRICT]
            # VERBOSITY [INFO]

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
        """Console function to submit a :py:class:`bsf.analyses.picard.CollectHiSeqXPfFailMetrics` analysis.

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

        analysis: CollectHiSeqXPfFailMetrics = cls.from_config_file_path(config_path=configuration_path)

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

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a :py:class:`bsf.analyses.picard.CollectHiSeqXPfFailMetrics` analysis.

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


class FilePathDownsampleSamReadGroup(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathDownsampleSamReadGroup` models files in a directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar downsampled_bam: A down-sampled BAM file path.
    :type downsampled_bam: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathDownsampleSamReadGroup` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathDownsampleSamReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.downsampled_bam = prefix + '_downsampled.bam'

        return


class DownsampleSam(Analysis):
    """The :py:class:`bsf.analyses.picard.DownsampleSam` class represents the logic to run the
    Picard :literal:`DownsampleSam` tool.

    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    """

    name = 'Picard DownsampleSam Analysis'
    prefix = 'picard_downsample_sam'

    @classmethod
    def get_stage_name_read_group(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name: str) -> str:
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name: str) -> FilePathDownsampleSamReadGroup:
        """Get a :py:class:`bsf.analyses.picard.FilePathDownsampleSamReadGroup` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A :py:class:`bsf.analyses.picard.FilePathDownsampleSamReadGroup` object.
        :rtype: FilePathDownsampleSamReadGroup
        """
        return FilePathDownsampleSamReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

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
            java_archive_picard: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.DownsampleSam` object.

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
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        """
        super(DownsampleSam, self).__init__(
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

        self.java_archive_picard = java_archive_picard

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.picard.DownsampleSam` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(DownsampleSam, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.DownsampleSam` object.

        This method changes the :py:class:`bsf.analyses.picard.DownsampleSam.collection` attribute
        (:py:class:`bsf.ngs.Collection`) to update with FASTQ file paths.
        """

        def run_read_comparisons() -> None:
            """Private function to read a :py:class:`bsf.annotation.AnnotationSheet` specifying comparisons
            from a CSV file path.

            This implementation just adds all :py:class:`bsf.ngs.Sample` objects from the
            :py:attr:`bsf.analysis.Analysis.collection` instance variable
            (i.e., :py:class:`bsf.ngs.Collection` object) to the
            :py:attr:`bsf.analysis.Analysis.sample_list` instance variable.
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        super(DownsampleSam, self).run()

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        run_read_comparisons()

        # Picard DownsampleSam

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())

        for sample in self.sample_list:
            module_logger.debug('Sample name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    module_logger.debug('PairedReads name: %r', paired_reads.get_name())
                    module_logger.log(logging.DEBUG - 2, 'PairedReads: %r', paired_reads)

                    # Apply some sanity checks.

                    if paired_reads.reads_1 is None and paired_reads.reads_2 is not None:
                        raise Exception(f"The PairedReads object {paired_reads.get_name()!r} has a 'reads_2', "
                                        f"but no 'reads_1' object.")

                    reads = paired_reads.reads_1
                    if not (reads.file_path.endswith('.bam') or reads.file_path.endswith('.sam')):
                        raise Exception(f'The {self.name!s} can only work on BAM or SAM files. '
                                        f'Reads.file_path: {reads.file_path!r}')

                    prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    # Keep the original BAM file and modify the file_path in the Reads object.
                    bam_file_path = reads.file_path
                    reads.file_path = os.path.join(self.project_directory, file_path_read_group.downsampled_bam)
                    # (file_root, file_extension) = os.path.splitext(bam_file_path)
                    # reads.file_path = file_root + 'downsampled' + file_extension

                    # Create a Runnable for running the Picard DownsampleSam analysis.

                    runnable_picard_dss = self.add_runnable(
                        runnable=ConsecutiveRunnable(
                            name=prefix_read_group,
                            working_directory=self.project_directory))

                    # Create an Executable for running the Picard SamToFastq Runnable.

                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_picard_dss)

                    # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                    # runnable_step = RunnableStepMakeDirectory(
                    #    name='make_directory',
                    #    directory_path=file_path_read_group.output_directory)
                    # runnable_picard_dss.add_runnable_step(runnable_step=runnable_step)

                    # Create a new RunnableStep for the Picard DownsampleSam program.

                    runnable_step = RunnableStepPicard(
                        name='picard_downsample_sam',
                        java_temporary_path=runnable_picard_dss.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx2G',
                        java_jar_path=self.java_archive_picard,
                        picard_command='DownsampleSam')
                    runnable_picard_dss.add_runnable_step(runnable_step=runnable_step)

                    # Required Arguments
                    # INPUT [null]
                    runnable_step.add_picard_option(key='INPUT', value=bam_file_path)
                    # OUTPUT [null]
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.downsampled_bam)

                    # Optional Tool Arguments
                    # ACCURACY [1.0E-4]
                    # METRICS_FILE [null]
                    # PROBABILITY [1.0]
                    # FIXME: Add to the configuration file and documentation.
                    if 'DownsampleSam Probability' in paired_reads.annotation_dict:
                        runnable_step.add_picard_option(
                            key='PROBABILITY',
                            value=paired_reads.annotation_dict['DownsampleSam Probability'][0])
                    # RANDOM_SEED [1]
                    # REFERENCE_SEQUENCE [null]
                    # STRATEGY [ConstantMemory]

                    # Optional Common Arguments
                    # COMPRESSION_LEVEL [5]
                    # CREATE_INDEX [false]
                    # CREATE_MD5_FILE [false]
                    # MAX_RECORDS_IN_RAM [500000]
                    # QUIET [false]
                    # TMP_DIR [null]
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_picard_dss.temporary_directory_path(absolute=False))
                    # USE_JDK_DEFLATER [false]
                    # USE_JDK_INFLATER [false]
                    # VALIDATION_STRINGENCY [STRICT]
                    # VERBOSITY [INFO]

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'picard_downsample_sam_samples.csv'))),
            name='_'.join((self.project_name, 'picard_downsample_sam')))
        annotation_sheet.to_file_path()

        return

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            project_name: Optional[str] = None,
            sas_path: Optional[str] = None,
            *args, **kwargs) -> int:
        """Console function to submit a :py:class:`bsf.analyses.picard.DownsampleSam` analysis.

        This analysis requires either a :literal:`configuration` argument or
        :literal:`project_name` and :literal:`sas_path` arguments.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param project_name: A project name.
        :type project_name: str
        :param sas_path: A sample annotation sheet file path.
        :type sas_path: str
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if configuration_path == Configuration.get_global_file_path():
            if not project_name:
                raise Exception('The --project-name argument is required if configuration is not set.')
            if not sas_path:
                raise Exception('The --sas-path argument is required if configuration is not set.')

        analysis: DownsampleSam = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (INI) file.')

            analysis.project_name = project_name

        if sas_path:
            if not sas_path.endswith('.csv'):
                raise Exception('The --sas-path option requires a comma-separated value (CSV) file.')

            analysis.sas_file = sas_path

        analysis.run()
        analysis.check_state()
        analysis.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis.name)
        print('Project name:      ', analysis.project_name)
        print('Input directory:   ', analysis.input_directory)
        print('Project directory: ', analysis.project_directory)
        print('Genome directory:  ', analysis.genome_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a :py:class:`bsf.analyses.picard.DownsampleSam` analysis.

        This analysis requires either a positional :literal:`configuration` argument or
        :literal:`--project-name` and :literal:`--sas-path` arguments.

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
            help='project name')

        argument_parser.add_argument(
            '--sas-path',
            help='sample annotation sheet (*.csv) file path')

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
            sas_path=name_space.sas_path)


class FilePathSamToFastqReadGroup(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathSamToFastqReadGroup` class models read group-specific
    Picard :literal:`SamToFastq` files.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathSamToFastqReadGroup` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathSamToFastqReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class FilePathSamToFastqProject(FilePath):
    """The :py:class:`bsf.analyses.picard.FilePathSamToFastqProject` class models project-specific
    Picard :literal:`SamToFastq` files.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar sas_path_old: An old Sample Annotation Sheet file path.
    :type sas_path_old: str
    :ivar sas_path_new: A new Sample Annotation Sheet file path.
    :type sas_path_new: str
    """

    def __init__(self, prefix: str, prefix_analysis: str, project_name: str) -> None:
        """Initialise a :py:class:`bsf.analyses.picard.FilePathSamToFastqProject` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param prefix_analysis: A :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix_analysis: str
        :param project_name: A project name.
        :type project_name: str
        """
        super(FilePathSamToFastqProject, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sas_path_old = '_'.join((project_name, prefix_analysis, 'original.csv'))
        self.sas_path_new = '_'.join((project_name, prefix_analysis, 'samples.csv'))

        return


class SamToFastq(Analysis):
    """The :py:class:`bsf.analyses.picard.SamToFastq` class represents the logic to run the
    Picard :literal:`SamToFastq` tool.

    :cvar stage_name_read_group: A :py:attr:`bsf.analysis.Stage.name` for read group-specific
        :py:class:`bsf.procedure.Runnable` objects.
    :type stage_name_read_group: str
    :cvar stage_name_project: A :py:attr:`bsf.analysis.Stage.name` for project-specific
        :py:class:`bsf.procedure.Runnable` objects.
    :type stage_name_project: str
    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    :ivar include_non_pass_filter_reads: Request including non-pass filter reads.
    :type include_non_pass_filter_reads: bool | None
    :ivar drop_read_1: Request dropping of read 1.
    :type drop_read_1: bool | None
    :ivar drop_read_2: Request dropping of read 2.
    :type drop_read_2: bool | None
    """

    name = 'Picard SamToFastq Analysis'
    prefix = 'picard_sam_to_fastq'

    stage_name_read_group = '_'.join((prefix, 'read_group'))
    stage_name_project = '_'.join((prefix, 'project'))

    @classmethod
    def get_stage_name_read_group(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_stage_name_project(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'project'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name: str) -> str:
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_prefix_project(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` (prefix) object  representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_project(), project_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name: str) -> FilePathSamToFastqReadGroup:
        """Get a :py:class:`bsf.analyses.picard.FilePathSamToFastqReadGroup` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A :py:class:`bsf.analyses.picard.FilePathSamToFastqReadGroup` object.
        :rtype: FilePathSamToFastqReadGroup
        """
        return FilePathSamToFastqReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

    @classmethod
    def get_file_path_project(cls, project_name: str, prefix_analysis: str) -> FilePathSamToFastqProject:
        """Get a :py:class:`bsf.analyses.picard.FilePathSamToFastqProject` object.

        :param project_name: A project name.
        :type project_name: str
        :param prefix_analysis: A :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix_analysis: str
        :return: A :py:class:`bsf.analyses.picard.FilePathSamToFastqProject` object.
        :rtype: FilePathSamToFastqProject
        """
        return FilePathSamToFastqProject(
            prefix=cls.get_prefix_project(project_name=project_name),
            prefix_analysis=prefix_analysis,
            project_name=project_name)

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
            java_archive_picard: Optional[str] = None,
            include_non_pf_reads: Optional[bool] = None,
            drop_read_1: Optional[bool] = None,
            drop_read_2: Optional[bool] = None):
        """Initialise a :py:class:`bsf.analyses.picard.SamToFastq` object.

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
        :param report_style_path: A eport style :literal:`CSS` file path.
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
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param include_non_pf_reads: Request including non-pass filter reads.
        :type include_non_pf_reads: bool | None
        :param drop_read_1: Request dropping of read 1.
        :type drop_read_1: bool | None
        :param drop_read_2: Request dropping of read 2.
        :type drop_read_2: bool | None
        """
        super(SamToFastq, self).__init__(
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

        self.java_archive_picard = java_archive_picard
        self.include_non_pass_filter_reads = include_non_pf_reads
        self.drop_read_1 = drop_read_1
        self.drop_read_2 = drop_read_2

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.picard.SamToFastq` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """

        super(SamToFastq, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        option = 'include_non_pass_filter_reads'
        if configuration.config_parser.has_option(section=section, option=option):
            self.include_non_pass_filter_reads = configuration.config_parser.getboolean(section=section, option=option)

        option = 'drop_read_1'
        if configuration.config_parser.has_option(section=section, option=option):
            self.drop_read_1 = configuration.config_parser.getboolean(section=section, option=option)

        option = 'drop_read_2'
        if configuration.config_parser.has_option(section=section, option=option):
            self.drop_read_2 = configuration.config_parser.getboolean(section=section, option=option)

        return

    def _read_comparisons(self) -> None:
        """Private method to read a :py:class:`bsf.annotation.AnnotationSheet` object specifying comparisons from a
        CSV file path.

        This implementation just adds all :py:class:`bsf.ngs.Sample` objects from the
        :py:attr:`bsf.analysis.Analysis.collection` instance variable (i.e., :py:class:`bsf.ngs.Collection`) to the
        :py:attr:`bsf.analysis.Analysis.sample_list` instance variable.
        """

        self.sample_list.extend(self.collection.get_all_samples())

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.picard.SamToFastq` object.

        This method changes the :py:class:`bsf.analyses.picard.DownsampleSam.collection` attribute
        (:py:class:`bsf.ngs.Collection`) to update with FASTQ file paths.
        """

        # Start of the run() method body.

        super(SamToFastq, self).run()

        # Get the Picard tools Java Archive (JAR) file path.

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        self._read_comparisons()

        # Picard SamToFastq

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())
        stage_project = self.get_stage(name=self.get_stage_name_project())

        project_dependency_list: list[str] = list()

        for sample in self.sample_list:
            module_logger.debug('Sample name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    module_logger.debug('PairedReads name: %r', paired_reads.get_name())
                    module_logger.log(logging.DEBUG - 2, 'PairedReads: %r', paired_reads)

                    # Apply some sanity checks.

                    if paired_reads.reads_1 is None and paired_reads.reads_2 is not None:
                        raise Exception(f"PairedReads object {paired_reads.get_name()!r} has a 'reads_2', "
                                        f"but no 'reads_1' object.")

                    reads = paired_reads.reads_1
                    if reads.file_path.endswith('.bam'):
                        bam_file_path = reads.file_path
                        prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                        file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                        # Get the SAM header of a BAM file to extract the read group (@RG), amongst other things.

                        # Open the BAM file, while not checking sequence (@SQ) entries.
                        # De-multiplexed, unaligned BAM files have no reference sequence dictionary.

                        alignment_file = AlignmentFile(filename=reads.file_path, mode='rb', check_sq=False)

                        alignment_header_dict = alignment_file.header.to_dict()

                        rg_dict: dict[str, str]
                        for rg_index, rg_dict in enumerate(alignment_header_dict['RG']):
                            platform_unit = SafeFileName.get_safe_file_name(file_name=rg_dict['PU'])

                            read_group_list = ['@RG']
                            read_group_list.extend([':'.join((key, value)) for key, value in rg_dict.items()])

                            if rg_index == 0:
                                # For the first read group, modify the PairedReads object in place.
                                paired_reads.read_group = '\\t'.join(read_group_list)
                                paired_reads.reads_1.name = platform_unit + '_1'

                                paired_reads.reads_1.file_path = os.path.join(
                                    self.project_directory,
                                    file_path_read_group.output_directory,
                                    platform_unit + '_1.fastq.gz')

                                paired_reads.reads_2 = Reads(
                                    name=platform_unit + '_2',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_read_group.output_directory,
                                        platform_unit + '_2.fastq.gz'),
                                    file_type=paired_reads.reads_1.file_type,
                                    lane=paired_reads.reads_1.lane,
                                    read=paired_reads.reads_1.read,
                                    chunk=paired_reads.reads_1.chunk,
                                    weak_reference_paired_reads=weakref.ReferenceType(paired_reads))

                                # Retain the original BAM file path as annotation.
                                paired_reads.add_annotation(key='BAM File', value=bam_file_path)
                            else:
                                # For further read groups, create new PairedReads objects.
                                reads_1 = Reads(
                                    name=platform_unit + '_1',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_read_group.output_directory,
                                        platform_unit + '_1.fastq.gz'),
                                    file_type=paired_reads.reads_1.file_type,
                                    lane=paired_reads.reads_1.lane,
                                    read='R1',
                                    chunk=paired_reads.reads_1.chunk)

                                reads_2 = Reads(
                                    name=platform_unit + '_2',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_read_group.output_directory,
                                        platform_unit + '_2.fastq.gz'),
                                    file_type=paired_reads.reads_1.file_type,
                                    lane=paired_reads.reads_1.lane,
                                    read='R2',
                                    chunk=paired_reads.reads_1.chunk)

                                new_paired_reads = PairedReads(
                                    reads_1=reads_1,
                                    reads_2=reads_2,
                                    read_group='\\t'.join(read_group_list))

                                # Retain the original BAM file path as annotation.
                                new_paired_reads.add_annotation(key='BAM File', value=bam_file_path)

                                sample.add_paired_reads(paired_reads=new_paired_reads)

                        alignment_file.close()

                        # Create a Runnable for running the Picard SamToFastq analysis.

                        runnable_read_group = self.add_runnable(
                            runnable=ConsecutiveRunnable(
                                name=prefix_read_group,
                                working_directory=self.project_directory))
                        self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                        # Record the Executable.name for the project dependency.

                        project_dependency_list.append(runnable_read_group.name)

                        # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                        runnable_step = RunnableStepMakeDirectory(
                            name='mkdir',
                            directory_path=file_path_read_group.output_directory)
                        runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                        # Create a new RunnableStep for the Picard SamToFastq program.

                        runnable_step = RunnableStepPicard(
                            name='picard_sam_to_fastq',
                            java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx2G',
                            java_jar_path=self.java_archive_picard,
                            picard_command='SamToFastq')
                        runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                        # Required Arguments
                        # FASTQ [null]
                        # INPUT [null]
                        runnable_step.add_picard_option(key='INPUT', value=bam_file_path)

                        # Optional Tool Arguments
                        # CLIPPING_ACTION [null]
                        runnable_step.add_picard_option(key='CLIPPING_ACTION', value='N')
                        # CLIPPING_ATTRIBUTE [null]
                        runnable_step.add_picard_option(key='CLIPPING_ATTRIBUTE', value='XT')
                        # CLIPPING_MIN_LENGTH [0]
                        # COMPRESS_OUTPUTS_PER_RG [false]
                        runnable_step.add_picard_option(key='COMPRESS_OUTPUTS_PER_RG', value='true')
                        # INCLUDE_NON_PF_READS [false]
                        if self.include_non_pass_filter_reads:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
                        else:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
                        # INCLUDE_NON_PRIMARY_ALIGNMENTS [false]
                        # INTERLEAVE [false]
                        # OUTPUT_DIR [null]
                        runnable_step.add_picard_option(
                            key='OUTPUT_DIR',
                            value=file_path_read_group.output_directory)
                        # OUTPUT_PER_RG [false]
                        runnable_step.add_picard_option(key='OUTPUT_PER_RG', value='true')
                        # QUALITY [null]
                        # RE_REVERSE [true]
                        # READ1_MAX_BASES_TO_WRITE [null]
                        # READ1_TRIM [0]
                        # READ2_MAX_BASES_TO_WRITE [null]
                        # READ2_TRIM [0]
                        # RG_TAG [PU]
                        # SECOND_END_FASTQ [null]
                        # UNPAIRED_FASTQ [null]

                        # Optional Common Arguments
                        # COMPRESSION_LEVEL [5]
                        runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                        # CREATE_INDEX [false]
                        # CREATE_MD5_FILE [false]
                        runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                        # MAX_RECORDS_IN_RAM [500000]
                        # QUIET [false]
                        # REFERENCE_SEQUENCE [null]
                        # TMP_DIR [null]
                        runnable_step.add_picard_option(
                            key='TMP_DIR',
                            value=runnable_read_group.temporary_directory_path(absolute=False))
                        # USE_JDK_DEFLATER [false]
                        # USE_JDK_INFLATER [false]
                        # VALIDATION_STRINGENCY [STRICT]
                        # VERBOSITY [INFO]
                        runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

        # Create a Runnable for pruning the sample annotation sheet.

        prefix_project = self.get_prefix_project(project_name=self.project_name)

        file_path_project = self.get_file_path_project(project_name=self.project_name, prefix_analysis=self.prefix)

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(self.project_directory, file_path_project.sas_path_old),
            name=prefix_project)
        annotation_sheet.to_file_path()

        runnable_project = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=prefix_project,
                working_directory=self.project_directory))
        executable_project = self.set_stage_runnable(
            stage=stage_project,
            runnable=runnable_project)
        executable_project.dependencies.extend(project_dependency_list)

        # Create a new RunnableStep.

        runnable_step = RunnableStepCollectionPruneFastq(
            name='prune_sample_annotation_sheet',
            obsolete_file_path_list=[
                # file_path_project.sas_path_old,
            ],
            file_path_old=file_path_project.sas_path_old,
            file_path_new=file_path_project.sas_path_new,
            minimum_size=1024,
            drop_read_1=self.drop_read_1,
            drop_read_2=self.drop_read_2)
        runnable_project.add_runnable_step(runnable_step=runnable_step)

        return

    def prune(self) -> None:
        """Prune the project directory by replacing FASTQ files with (empty) status files (:literal:`*.truncated`).

        The status files are recognised by the
        :py:class:`bsf.executables.collection.RunnableStepCollectionPruneFastq` class
        so that :py:class:`bsf.ngs.Reads` and :py:class:`bsf.ngs.PairedReads` objects are kept in the
        sample annotation sheet.
        """
        super(SamToFastq, self).run()

        self._read_comparisons()

        module_logger.debug('Prune!')

        for sample in self.sample_list:
            module_logger.debug('Sample name: %r', sample.name)
            module_logger.log(logging.DEBUG - 1, 'Sample: %r', sample)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    module_logger.debug('PairedReads name: %r', paired_reads.get_name())
                    module_logger.log(logging.DEBUG - 2, 'PairedReads: %r', paired_reads)

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    def prune_file_path(reads: Reads) -> None:
                        """Private function to prune files of a Reads object.

                        :param reads: A :py:class:`bsf.ngs.Reads` object.
                        :type reads: Reads
                        """
                        if os.path.isabs(reads.file_path):
                            file_path = reads.file_path
                        else:
                            file_path = os.path.join(
                                self.genome_directory,
                                file_path_read_group.output_directory,
                                reads.file_path)

                        if os.path.exists(file_path):
                            open(file=file_path + '.truncated', mode='wt').close()
                            os.remove(file_path)

                        # Also remove an eventual MD5 checksum file.
                        file_path += '.md5'
                        if os.path.exists(file_path):
                            os.remove(file_path)

                        return

                    if paired_reads.reads_1 is not None:
                        prune_file_path(reads=paired_reads.reads_1)

                    if paired_reads.reads_2 is not None:
                        prune_file_path(reads=paired_reads.reads_2)

        return

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            project_name: Optional[str] = None,
            sas_path: Optional[str] = None,
            drop_read_1: Optional[bool] = None,
            drop_read_2: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a :py:class:`bsf.analyses.picard.SamToFastq` analysis.

        This analysis requires either a :literal:`configuration` argument or
        :literal:`project_name` and :literal:`sas_path` arguments.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param project_name: A project name.
        :type project_name: str
        :param sas_path: A sample annotation sheet file path.
        :type sas_path: str
        :param drop_read_1: Request dropping read 1.
        :type drop_read_1: bool | None
        :param drop_read_2: Request dropping read 2.
        :type drop_read_2: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        if configuration_path == Configuration.get_global_file_path():
            if not project_name:
                raise Exception('The --project-name argument is required if configuration is not set.')
            if not sas_path:
                raise Exception('The --sas-path argument is required if configuration is not set.')

        analysis: SamToFastq = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (INI) file.')

            analysis.project_name = project_name

        if sas_path:
            if not sas_path.endswith('.csv'):
                raise Exception('The --sas-path option requires a comma-separated value (CSV) file.')

            analysis.sas_file = sas_path

        if drop_read_1:
            analysis.drop_read_1 = drop_read_1

        if drop_read_2:
            analysis.drop_read_2 = drop_read_2

        if stage_name == 'prune':
            analysis.prune()
        else:
            analysis.run()
            analysis.check_state()
            analysis.submit(name=stage_name, drms_submit=drms_submit)

            print(analysis.name)
            print('Project name:      ', analysis.project_name)
            print('Input directory:   ', analysis.input_directory)
            print('Project directory: ', analysis.project_directory)
            print('Genome directory:  ', analysis.genome_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a :py:class:`bsf.analyses.picard.SamToFastq` analysis.

        This analysis requires either a positional :literal:`configuration` argument or
        :literal:`--project-name` and :literal:`--sas-path` arguments.

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
            help='project name')

        argument_parser.add_argument(
            '--sas-path',
            help='sample annotation sheet (*.csv) file path')

        argument_parser.add_argument(
            '--drop-read-1',
            action='store_true',
            help='drop read 1')

        argument_parser.add_argument(
            '--drop-read-2',
            action='store_true',
            help='drop read 2')

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
            sas_path=name_space.sas_path,
            drop_read_1=name_space.drop_read_1,
            drop_read_2=name_space.drop_read_2)
