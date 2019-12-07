# -*- coding: utf-8 -*-
"""Picard Analysis module.

A package of classes and methods modelling Picard analyses data files and data directories.
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
from __future__ import print_function

import os
import re
import sys
import warnings
import weakref

import pysam

import bsf.analyses.illumina_to_bam_tools
import bsf.analysis
import bsf.annotation
import bsf.executables.cloud
import bsf.illumina
import bsf.ngs
import bsf.procedure
import bsf.process
import bsf.standards


class PicardIlluminaRunFolder(bsf.analysis.Analysis):
    """The C{bsf.analyses.picard.PicardIlluminaRunFolder} class of Picard Analyses acting on Illumina Run Folders.

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
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | None
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool | None
    """

    name = 'Picard PicardIlluminaRunFolder Analysis'
    prefix = 'picard_illumina_run_folder'

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
            classpath_picard=None,
            force=False):
        """Initialise a C{bsf.analyses.picard.PicardIlluminaRunFolder} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
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
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
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
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
        @return:
        @rtype:
        """
        super(PicardIlluminaRunFolder, self).__init__(
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

        self.run_directory = run_directory
        self.intensity_directory = intensity_directory
        self.basecalls_directory = basecalls_directory
        self.experiment_name = experiment_name
        self.classpath_picard = classpath_picard
        self.force = force

        self._irf = None
        """ @type _irf: bsf.illumina.RunFolder | None """

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.picard.PicardIlluminaRunFolder} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
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

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.picard.PicardIlluminaRunFolder} C{bsf.analysis.Analysis}.

        @return:
        @rtype:
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(
                'The ' + self.name + " requires a 'run_directory' configuration option.")

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=bsf.standards.FilePath.get_illumina_run(absolute=True))

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The ' + self.name + " 'run_directory' " + repr(self.run_directory) + ' is not a valid directory.')

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise bsf.illumina.RunFolderNotComplete(
                'The Illumina Run Folder ' + repr(self.run_directory) + ' is not complete.')

        # Define an 'Intensities' directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Illumina Run Folder path.

        if self.intensity_directory:
            self.intensity_directory = self.configuration.get_absolute_path(
                file_path=self.intensity_directory,
                default_path=self.run_directory)
        else:
            self.intensity_directory = os.path.join(self.run_directory, 'Data', 'Intensities')

        # Check that the Intensities directory exists.

        if not os.path.isdir(self.intensity_directory):
            raise Exception(
                'The ' + self.name + " 'intensity_directory' " + repr(self.intensity_directory) +
                ' is not a valid directory.')

        # Define a 'BaseCalls' directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Intensities directory path.

        if self.basecalls_directory:
            self.basecalls_directory = self.configuration.get_absolute_path(
                file_path=self.basecalls_directory,
                default_path=self.intensity_directory)
        else:
            self.basecalls_directory = os.path.join(self.intensity_directory, 'BaseCalls')

        # Check that the BaseCalls directory exists.

        if not os.path.isdir(self.basecalls_directory):
            raise Exception(
                'The ' + self.name + " 'basecalls_directory' " + repr(self.basecalls_directory) +
                ' is not a valid directory.')

        self._irf = bsf.illumina.RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g. BSF_0000) is used as the prefix for archive BAM files.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = self._irf.run_parameters.get_experiment_name
            if not self.experiment_name:
                raise Exception(
                    'The ' + self.name + " 'experiment_name' was not provided " +
                    'and could not be read from the Illumina Run Folder configuration.')

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, self._irf.run_information.flow_cell))

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception(
                    'The ' + self.name + " requires a 'classpath_picard' configuration option.")

        # Call the run method of the super class after the project_name has been defined.

        super(PicardIlluminaRunFolder, self).run()

        return


class ExtractIlluminaBarcodesSheet(bsf.annotation.AnnotationSheet):
    """The C{bsf.analyses.picard.ExtractIlluminaBarcodesSheet} class represents a
    Tab-Separated Value (TSV) table of library information for the
    C{bsf.analyses.picard.ExtractIlluminaBarcodes} C{bsf.analysis.Analysis}.

    Attributes:
    """

    _file_type = 'excel-tab'

    _header_line = True

    _field_names = [
        'barcode_sequence_1',
        'barcode_sequence_2',
        'barcode_name',
        'library_name',
    ]

    _test_methods = dict()


class IlluminaBasecallsToSamSheet(bsf.annotation.AnnotationSheet):
    """The C{bsf.analyses.picard.IlluminaBasecallsToSamSheet} class represents a
    Tab-Separated Value (TSV) table of library information for the
    C{bsf.analyses.picard.ExtractIlluminaBarcodes} C{bsf.analysis.Analysis}.

    Attributes:
    """

    _file_type = 'excel-tab'

    _header_line = True

    _field_names = [
        'OUTPUT',
        'SAMPLE_ALIAS',
        'LIBRARY_NAME',
        'BARCODE_1',
        'BARCODE_2',
    ]

    _test_methods = dict()

    def adjust(self, barcode_length_tuple, unassigned_file_path):
        """Adjust the C{bsf.analyses.picard.IlluminaBasecallsToSamSheet}.

        Insert the file path for unassigned reads.
        Remove BARCODE_N index columns, which annotation is all empty.
        @param barcode_length_tuple: Python C{tuple} of Python C{int} barcode length objects
        @type barcode_length_tuple: (int, int)
        @param unassigned_file_path: File path for unassigned reads
        @type unassigned_file_path: str
        @return:
        @rtype:
        """
        # The IlluminaBasecallsToSamSheet needs adjusting ...
        if len(self.row_dicts) == 1 and len(self.row_dicts[0]['BARCODE_1']) == 0 and len(
                self.row_dicts[0]['BARCODE_2']) == 0:
            # ... if a single sample, but neither BARCODE_1 nor BARCODE_2 were defined,
            # BARCODE_1 needs setting to 'N'.
            self.row_dicts[0]['BARCODE_1'] = 'N'
        else:
            # ... in all other cases, a last row for unmatched barcode sequences needs adding.
            self.row_dicts.append({
                'BARCODE_1': 'N' * barcode_length_tuple[0],
                'BARCODE_2': 'N' * barcode_length_tuple[1],
                'OUTPUT': unassigned_file_path,
                'SAMPLE_ALIAS': 'Unmatched',
                'LIBRARY_NAME': self.row_dicts[0]['LIBRARY_NAME'],
            })

        # Further adjust the IlluminaBaseCallsToSamSheet and remove any BARCODE_N columns not represented
        # in the barcode length list.
        for index in range(0, 1 + 1):
            if barcode_length_tuple[index] == 0:
                # Remove the 'BARCODE_N' field from the list of field names.
                if 'BARCODE_' + str(index + 1) in self.field_names:
                    self.field_names.remove('BARCODE_' + str(index + 1))
                # Remove the 'BARCODE_N' entry from each row_dict object, since csv.DictWriter requires it.
                for row_dict in self.row_dicts:
                    row_dict.pop('BARCODE_' + str(index + 1), None)

        return


class FilePathExtractIlluminaCell(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathExtractIlluminaCell} models flow cell-specific file paths.

    Attributes:
    @ivar prefix_cell: Non-standard, flow cell-specific (i.e. project_name) prefix
    @type prefix_cell: str
    @ivar sample_annotation_sheet_csv: Sample Annotation Sheet CSV file
    @type sample_annotation_sheet_csv: str
    @see: FilePathBamIndexDecoderCell
    @see: FilePathIlluminaDemultiplexSamCell
    """

    def __init__(self, prefix, project_name):
        """Initialise a C{bsf.analyses.picard.FilePathExtractIlluminaCell} object.

        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @return:
        @rtype:
        """
        super(FilePathExtractIlluminaCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the flow cell-specific prefix accordingly.

        self.prefix_cell = project_name

        self.sample_annotation_sheet_csv = self.prefix_cell + '_samples.csv'

        return


class FilePathExtractIlluminaLane(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathExtractIlluminaLane} models lane-specific file paths.

    Attributes:
    @ivar prefix_lane: Non-standard, lane-specific (i.e. project_name and lane) prefix
    @type prefix_lane: str
    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar samples_directory: Samples directory
    @type samples_directory: str
    @ivar barcode_tsv: Barcode TSV file for Picard I{ExtractIlluminaBarcodes}
    @type barcode_tsv: str
    @ivar library_tsv: Library TSV file for Picard I{IlluminaBasecallsToSam}
    @type library_tsv: str
    @ivar metrics_tsv: Metrics TSV file
    @type metrics_tsv: str
    @ivar metrics_fraction_pdf: Lane-specific Picard I{ExtractIlluminaBarcodes} fraction metrics PDF file
    @type metrics_fraction_pdf: str
    @ivar metrics_fraction_png: Lane-specific Picard I{ExtractIlluminaBarcodes} fraction metrics PNG file
    @type metrics_fraction_png: str
    @ivar metrics_number_pdf: Lane-specific Picard I{ExtractIlluminaBarcodes} number metrics PDF file
    @type metrics_number_pdf: str
    @ivar metrics_number_png: Lane-specific Picard I{ExtractIlluminaBarcodes} number metrics PNG file
    @type metrics_number_png: str
    @see: FilePathBamIndexDecoderLane
    @see: FilePathIlluminaDemultiplexSamLane
    """

    def __init__(self, prefix, project_name, lane):
        """Initialise a C{bsf.analyses.picard.FilePathExtractIlluminaLane} object.

        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @return:
        @rtype:
        """
        super(FilePathExtractIlluminaLane, self).__init__(prefix=prefix)

        # Non-standard Non-standard project_name_lane prefix
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
    """The C{bsf.analyses.picard.ExtractIlluminaRunFolder} class extracts data from an Illumina Run Folder.

    The analysis is based on Picard I{ExtractIlluminaBarcodes} and Picard I{IlluminaBasecallsToSam}.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar samples_directory: BSF samples directory
    @type samples_directory: str | None
    @ivar library_path: Library annotation file path
    @type library_path: str | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar max_mismatches: Maximum number of mismatches
    @type max_mismatches: int | None
    @ivar min_base_quality: Minimum base quality
    @type min_base_quality: int | None
    @ivar sequencing_centre: Sequencing centre
    @type sequencing_centre: str | None
    @ivar lanes: Number of lanes on the flow cell
    @type lanes: int | None
    @ivar vendor_quality_filter: Vendor quality filter
    @type vendor_quality_filter: bool
    """

    name = 'Picard Extract Illumina Run Folder Analysis'
    prefix = 'extract_illumina_run_folder'

    @classmethod
    def get_file_path_cell(cls, project_name):
        """Get a C{FilePathExtractIlluminaCell} object.

        @param project_name: Project name
        @type project_name: str
        @return: C{FilePathExtractIlluminaCell} object
        @rtype: FilePathExtractIlluminaCell
        """
        return FilePathExtractIlluminaCell(
            prefix=cls.get_prefix_cell(project_name=project_name),
            project_name=project_name)

    @classmethod
    def get_file_path_lane(cls, project_name, lane):
        """Get a C{FilePathExtractIlluminaLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @return: C{FilePathExtractIlluminaLane} object
        @rtype: FilePathExtractIlluminaLane
        """
        return FilePathExtractIlluminaLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane)

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
            classpath_picard=None,
            force=False,
            library_path=None,
            samples_directory=None,
            mode_directory=None,
            mode_file=None,
            max_mismatches=None,
            min_base_quality=None,
            sequencing_centre=None,
            lanes=None,
            vendor_quality_filter=None):
        """Initialise a C{bsf.analyses.picard.ExtractIlluminaRunFolder} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
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
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
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
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
        @param samples_directory: BSF samples directory
        @type samples_directory: str | None
        @param library_path: Library annotation file path
        @type library_path: str | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param max_mismatches: Maximum number of mismatches
        @type max_mismatches: int | None
        @param min_base_quality: Minimum base quality
        @type min_base_quality: int | None
        @param sequencing_centre: Sequencing centre
        @type sequencing_centre: str | None
        @param lanes: Number of lanes on the flow cell
        @type lanes: int | None
        @param vendor_quality_filter: Vendor quality filter
        @type vendor_quality_filter: bool | None
        @return:
        @rtype:
        """
        super(ExtractIlluminaRunFolder, self).__init__(
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
            sample_list=sample_list,
            run_directory=run_directory,
            intensity_directory=intensity_directory,
            basecalls_directory=basecalls_directory,
            experiment_name=experiment_name,
            classpath_picard=classpath_picard,
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
    def get_experiment_directory(self):
        """Get the experiment directory.

        The experiment directory is a concatenation of the sequences directory and the project name.
        @return: Experiment directory
        @rtype: str | None
        """
        if self.samples_directory and self.project_name:
            return os.path.join(self.samples_directory, self.project_name)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.picard.ExtractIlluminaRunFolder} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
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

    def run(self):
        """Run a C{bsf.analyses.picard.ExtractIlluminaRunFolder} C{bsf.analysis.Analysis}.

        @return:
        @rtype:
        """

        def run_get_sample_file_name(sample_name):
            """Private function to format sample-specific BAM file names (i.e. project_lane#sample.bam).

            @param sample_name:
            @type sample_name: str
            @return:
            @rtype: str
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # Start of the run() method body.

        super(ExtractIlluminaRunFolder, self).run()

        self.samples_directory = self.configuration.get_absolute_path(
            file_path=self.samples_directory,
            default_path=bsf.standards.FilePath.get_samples(absolute=True))

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(
                'The ' + self.name + " 'samples_directory' " + repr(self.samples_directory) +
                ' is not a valid directory.')

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get sequencing centre information.

        if not self.sequencing_centre:
            self.sequencing_centre = bsf.standards.Operator.get_sequencing_centre()
            if not self.sequencing_centre:
                raise Exception(
                    'The ' + self.name + " requires a 'sequencing_centre' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self.vendor_quality_filter is None:
            self.vendor_quality_filter = bsf.standards.VendorQualityFilter.get_vendor_quality_filter(
                flow_cell_type=self._irf.run_parameters.get_flow_cell_type)

        # Get the library annotation sheet.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(self.library_path):
            raise Exception(
                'The ' + self.name + " 'library_path' " + repr(self.library_path) + ' is not a valid file.')

        # Load the LibraryAnnotationSheet and validate.

        library_annotation_sheet = bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)
        """ @type library_annotation_sheet: bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet """

        if self.lanes is None:
            validation_messages = library_annotation_sheet.validate(
                lanes=self._irf.run_information.flow_cell_layout.lane_count)
        else:
            validation_messages = library_annotation_sheet.validate(
                lanes=self.lanes)

        if validation_messages:
            if self.force:
                warnings.warn('Validation of library annotation sheet ' +
                              repr(self.library_path) + ':\n' + validation_messages)
            else:
                raise Exception('Validation of library annotation sheet ' +
                                repr(self.library_path) + ':\n' + validation_messages)

        library_annotation_dict = library_annotation_sheet.get_annotation_dict()
        library_barcode_dict = library_annotation_sheet.get_barcode_length_dict()

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        file_path_cell = self.get_file_path_cell(project_name=self.project_name)

        # NOTE: Use a SampleAnnotationSheet as long as the Collection does not have a method to
        # write a de-normalised table.
        # Create a SampleAnnotationSheet in the project directory and
        # eventually transfer it into the experiment_directory.

        sample_annotation_sheet = bsf.analyses.illumina_to_bam_tools.BamIndexDecoder.get_sample_annotation_sheet(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        # For each lane in the library_annotation_dict ...
        # TODO: For the moment this depends on the lanes (keys) defined in the LibraryAnnotationSheet.
        # Not all lanes may thus get extracted.
        # TODO: For NextSeq instruments, it would be sufficient to require annotation for only lane one and
        # copy information to lanes two to four internally.

        cell_dependency_list = list()
        """ @type cell_dependency_list: list[str] """

        for lane_int in sorted(library_annotation_dict):
            lane_str = str(lane_int)
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

                eib_sheet.row_dicts.append({
                    'barcode_sequence_1': row_dict['barcode_sequence_1'],
                    'barcode_sequence_2': row_dict['barcode_sequence_2'],
                    'barcode_name': row_dict['sample_name'],
                    'library_name': row_dict['library_name'],
                })

                # Add a row to the lane-specific Picard IlluminaBasecallsToSamSheet.

                ibs_sheet.row_dicts.append({
                    'BARCODE_1': row_dict['barcode_sequence_1'],
                    'BARCODE_2': row_dict['barcode_sequence_2'],
                    'OUTPUT': os.path.join(
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    'SAMPLE_ALIAS': row_dict['sample_name'],
                    'LIBRARY_NAME': row_dict['library_name'],
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
                    'PairedReads Structure': self._irf.run_information.get_picard_read_structure,
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

            runnable_lane = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)

            cell_dependency_list.append(executable_lane.name)

            # Create an output_directory in the project_directory.

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_output_directory',
                directory_path=file_path_lane.output_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Create a samples_directory in the project_directory.

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_samples_directory',
                directory_path=file_path_lane.samples_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Create a RunnableStep for Picard ExtractIlluminaBarcodes, only if index (barcode) reads are present.

            if barcode_number > 0:
                runnable_step = bsf.process.RunnableStepPicard(
                    name='picard_extract_illumina_barcodes',
                    java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    picard_classpath=self.classpath_picard,
                    picard_command='ExtractIlluminaBarcodes')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                # BASECALLS_DIR Required
                runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
                # OUTPUT_DIR [null] Default to BASECALLS_DIR
                runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.output_directory)
                # LANE Required
                runnable_step.add_picard_option(key='LANE', value=lane_str)
                # READ_STRUCTURE Required
                runnable_step.add_picard_option(
                    key='READ_STRUCTURE',
                    value=self._irf.run_information.get_picard_read_structure)
                # BARCODE [null]
                # BARCODE_FILE Required
                runnable_step.add_picard_option(key='BARCODE_FILE', value=file_path_lane.barcode_tsv)
                # METRICS_FILE Required
                runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
                # MAX_MISMATCHES [1]
                if self.max_mismatches is not None:
                    # Maximum mismatches for a barcode to be considered a match. Default value: '1'.
                    runnable_step.add_picard_option(key='MAX_MISMATCHES', value=str(self.max_mismatches))
                # MIN_MISMATCH_DELTA [1]
                # Minimum difference between number of mismatches in the best and
                # second best barcodes for a barcode to be considered a match.
                # MAX_NO_CALLS [2]
                # Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.
                # MINIMUM_BASE_QUALITY [0]
                if self.min_base_quality is not None:
                    # Minimum base quality. Any barcode bases falling below this quality will be considered
                    # a mismatch even in the bases match. Default value: '0'.
                    runnable_step.add_picard_option(key='MINIMUM_BASE_QUALITY', value=str(self.min_base_quality))
                # MINIMUM_QUALITY [2]
                # The minimum quality (after transforming 0s to 1s) expected from reads.
                # If qualities are lower than this value, an error is thrown.The default of 2 is what the
                # Illumina's spec describes as the minimum, but in practice the value has been observed lower.
                # COMPRESS_OUTPUTS [false]
                runnable_step.add_picard_option(key='COMPRESS_OUTPUTS', value='true')
                # NUM_PROCESSORS [1]
                runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_lane.temporary_directory_path(absolute=False))
                # VERBOSITY [INFO]
                # QUIET [false]
                # VALIDATION_STRINGENCY [STRICT]
                # COMPRESSION_LEVEL [5]
                # MAX_RECORDS_IN_RAM [500000]
                # CREATE_INDEX [false]
                # CREATE_MD5_FILE [false]
                # REFERENCE_SEQUENCE [null]
                # GA4GH_CLIENT_SECRETS [client_secrets.json]
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]
                # OPTIONS_FILE Required

                # Plot the metrics file.

                runnable_step = bsf.process.RunnableStep(
                    name='plot_metrics',
                    program='bsf_illumina_demultiplex_sam.R')
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='file-path', value=file_path_lane.metrics_tsv)

            # Picard IlluminaBasecallsToSam

            # Create a RunnableStep for Picard IlluminaBasecallsToSam.

            runnable_step = bsf.process.RunnableStepPicard(
                name='picard_illumina_basecalls_to_sam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx2G',
                picard_classpath=self.classpath_picard,
                picard_command='IlluminaBasecallsToSam')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # BASECALLS_DIR Required
            runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
            # BARCODES_DIR [null] Defaults to BASECALLS_DIR
            if barcode_number > 0:
                runnable_step.add_picard_option(key='BARCODES_DIR', value=file_path_lane.output_directory)
            # LANE Required
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # OUTPUT Deprecated
            # RUN_BARCODE Required
            runnable_step.add_picard_option(
                key='RUN_BARCODE',
                value=self._irf.run_parameters.get_flow_cell_barcode)
            # SAMPLE_ALIAS Deprecated
            # READ_GROUP_ID [null]
            runnable_step.add_picard_option(
                key='READ_GROUP_ID',
                value='_'.join((self._irf.run_parameters.get_flow_cell_barcode, lane_str)))
            # LIBRARY_NAME Deprecated
            # SEQUENCING_CENTER Required
            runnable_step.add_picard_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)
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
            # PLATFORM [illumina]
            # PLATFORM The name of the sequencing technology that produced the read.
            # NOTE: IlluminaToBam defaults to 'ILLUMINA'.
            # runnable_step.add_picard_option(key='PLATFORM', value='ILLUMINA')
            # INCLUDE_BC_IN_RG_TAG [false]
            # READ_STRUCTURE Required
            runnable_step.add_picard_option(
                key='READ_STRUCTURE',
                value=self._irf.run_information.get_picard_read_structure)
            # BARCODE_PARAMS Deprecated
            # LIBRARY_PARAMS Required
            runnable_step.add_picard_option(key='LIBRARY_PARAMS', value=file_path_lane.library_tsv)
            # ADAPTERS_TO_CHECK [INDEXED, DUAL_INDEXED, NEXTERA_V2, FLUIDIGM]
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # FIVE_PRIME_ADAPTER [null]
            # THREE_PRIME_ADAPTER [null]
            # NUM_PROCESSORS [0]
            # FIRST_TILE [null]
            # TILE_LIMIT [null]
            # PROCESS_SINGLE_TILE [null]
            # FORCE_GC [true]
            # APPLY_EAMSS_FILTER [true]
            # MAX_READS_IN_RAM_PER_TILE [1200000]
            # MINIMUM_QUALITY [2]
            # INCLUDE_NON_PF_READS [true]
            if self.vendor_quality_filter:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
            else:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
            # IGNORE_UNEXPECTED_BARCODES [false]
            # MOLECULAR_INDEX_TAG [RX]
            # MOLECULAR_INDEX_BASE_QUALITY_TAG [QX]
            # TAG_PER_MOLECULAR_INDEX [null]
            # BARCODE_POPULATION_STRATEGY [ORPHANS_ONLY]
            # INCLUDE_BARCODE_QUALITY [false]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # VERBOSITY [INFO]
            # QUIET [false]
            # VALIDATION_STRINGENCY [STRICT]
            # COMPRESSION_LEVEL [5]
            runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
            # MAX_RECORDS_IN_RAM [500000]
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # REFERENCE_SEQUENCE [null]
            # GA4GH_CLIENT_SECRETS [client_secrets.json]
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # OPTIONS_FILE Required

            # Create the experiment directory if it does not exist already.

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move the samples directory into the experiment directory.

            runnable_step = bsf.process.RunnableStepMove(
                name='move_samples_directory',
                source_path=file_path_lane.samples_directory,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move the metrics file into the experiment directory.

            if barcode_number > 0:
                runnable_step = bsf.process.RunnableStepMove(
                    name='move_metrics_tsv',
                    source_path=file_path_lane.metrics_tsv,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step = bsf.process.RunnableStepMove(
                    name='move_metrics_pdf',
                    source_path=file_path_lane.metrics_number_pdf,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step = bsf.process.RunnableStepMove(
                    name='move_metrics_png',
                    source_path=file_path_lane.metrics_number_png,
                    target_path=experiment_directory)
                runnable_lane.add_runnable_step(runnable_step=runnable_step)

        # Finally, write the flow cell-specific SampleAnnotationSheet to the internal file path,
        # but keep the order of field names defined above.

        sample_annotation_sheet.to_file_path(adjust_field_names=False)

        # Create a flow-cell specific Runnable.

        runnable_cell = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        # Move the SampleAnnotationSheet from the project_directory to the experiment_directory.

        if os.path.exists(sample_annotation_sheet.file_path):
            runnable_step = bsf.process.RunnableStepMove(
                name='move_sample_annotation',
                source_path=file_path_cell.sample_annotation_sheet_csv,
                target_path=experiment_directory)
            runnable_cell.add_runnable_step(runnable_step=runnable_step)

        # Change directory and file access permissions.

        runnable_step = bsf.process.RunnableStepChangeMode(
            name='chmod',
            file_path=experiment_directory,
            mode_directory=self.mode_directory,
            mode_file=self.mode_file)
        runnable_cell.add_runnable_step(runnable_step=runnable_step)

        return


class FilePathIlluminaMultiplexSamLane(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathIlluminaMultiplexSamLane} class models files in a directory.

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
        """Initialise a C{bsf.analyses.picard.FilePathIlluminaMultiplexSamLane} object.

        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param experiment_directory: Experiment-specific directory
        @type experiment_directory: str
        @return:
        @rtype:
        """
        super(FilePathIlluminaMultiplexSamLane, self).__init__(prefix=prefix)

        self.unsorted_bam = prefix + '_unsorted.bam'
        self.unsorted_md5 = prefix + '_unsorted.bam.md5'
        self.sorted_bam = prefix + '_sorted.bam'
        self.sorted_md5 = prefix + '_sorted.bam.md5'
        # The final BAM and MD5 files are non-standard in that they do not contain a prefix and
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
    """The C{bsf.analyses.picard.IlluminaMultiplexSam} class represents the
    Picard I{IlluminaBasecallsToMultiplexSam} analysis.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar sequencing_centre: Sequencing centre
    @type sequencing_centre: str | None
    @ivar sequences_directory: Sequences directory to store archive BAM files
    @type sequences_directory: str | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar eamss_filter: Apply the Illumina EAMSS or Read Segment Quality Control Metric filter
    @type eamss_filter: bool | None
    @ivar vendor_quality_filter: Vendor quality filter
    @type vendor_quality_filter: bool
    @ivar compression_level: (Zlib) Compression level
    @type compression_level: int | None
    @ivar cloud_account: I{Microsoft Azure Storage Account} name
    @type cloud_account: str | None
    @ivar cloud_container: I{Microsoft Azure Blob Service} container name
    @type cloud_container: str | None
    """

    name = 'Picard IlluminaMultiplexSam Analysis'
    prefix = 'illumina_multiplex_sam'

    @classmethod
    def get_stage_name_cloud(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cloud'))

    @classmethod
    def get_prefix_cloud(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_cloud(), project_name, lane))

    @classmethod
    def get_file_path_lane(cls, project_name, lane, experiment_directory):
        """Get a C{FilePathIlluminaMultiplexSamLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param experiment_directory: BSF experiment directory
        @type experiment_directory: str
        @return: C{FilePathIlluminaMultiplexSamLane} object
        @rtype: FilePathIlluminaMultiplexSamLane
        """
        return FilePathIlluminaMultiplexSamLane(
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
            classpath_picard=None,
            force=False,
            sequencing_centre=None,
            sequences_directory=None,
            mode_directory=None,
            mode_file=None,
            eamss_filter=None,
            vendor_quality_filter=None,
            compression_level=None,
            cloud_account=None,
            cloud_container=None):
        """Initialise a C{bsf.analyses.picard.IlluminaMultiplexSam} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
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
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
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
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
        @param sequencing_centre: Sequencing centre
        @type sequencing_centre: str | None
        @param sequences_directory: Sequences directory to store archive BAM files
        @type sequences_directory: str | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param eamss_filter: Apply the Illumina EAMSS or Read Segment Quality Control Metric filter
        @type eamss_filter: bool | None
        @param vendor_quality_filter: Vendor quality filter
        @type vendor_quality_filter: bool | None
        @param compression_level: (Zlib) Compression level
        @type compression_level: int | None
        @param cloud_account: I{Microsoft Azure Storage Account} name
        @type cloud_account: str | None
        @param cloud_container: I{Microsoft Azure Blob Service} container name
        @type cloud_container: str | None
        """
        super(IlluminaMultiplexSam, self).__init__(
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
            sample_list=sample_list,
            run_directory=run_directory,
            intensity_directory=intensity_directory,
            basecalls_directory=basecalls_directory,
            experiment_name=experiment_name,
            classpath_picard=classpath_picard,
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
        """Set instance variables of a C{bsf.analyses.picard.IlluminaMultiplexSam} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
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

        return

    def run(self):
        """Run this C{bsf.analyses.picard.IlluminaMultiplexSam} C{bsf.analysis.Analysis}.

        @return:
        @rtype:
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(
                'The ' + self.name + "requires a 'run_directory' configuration option.")

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=bsf.standards.FilePath.get_illumina_run(absolute=True))

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The ' + self.name + " 'run_directory' " + repr(self.run_directory) + ' is not a valid directory.')

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise bsf.illumina.RunFolderNotComplete(
                'The Illumina Run Folder ' + repr(self.run_directory) + ' is not complete.')

        irf = bsf.illumina.RunFolder.from_file_path(file_path=self.run_directory)

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
            self.sequencing_centre = bsf.standards.Operator.get_sequencing_centre()
            if not self.sequencing_centre:
                raise Exception('The ' + self.name + "requires a 'sequencing_centre' configuration option.")

        # Define the sequences directory in which to create the experiment directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # An absolute path cannot be prepended.

        if self.sequences_directory:
            self.sequences_directory = self.configuration.get_absolute_path(
                file_path=self.sequences_directory)
        else:
            self.sequences_directory = bsf.standards.FilePath.get_sequences(absolute=True)

        # As a safety measure, to prevent creation of rogue directory paths, the sequences_directory has to exist.

        if not os.path.isdir(self.sequences_directory):
            raise Exception('The ' + self.name + " 'sequences_directory' " + repr(self.sequences_directory) +
                            ' does not provide a valid directory.')

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the Picard tools Java Archive (JAR) class path directory

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('The ' + self.name + " requires a 'classpath_picard' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self.vendor_quality_filter is None:
            self.vendor_quality_filter = bsf.standards.VendorQualityFilter.get_vendor_quality_filter(
                flow_cell_type=irf.run_parameters.get_flow_cell_type)

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

            runnable_lane = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
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

            runnable_step = bsf.process.RunnableStepPicard(
                name='picard_illumina_basecalls_to_multiplex_sam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx24G',
                picard_classpath=self.classpath_picard,
                picard_command='IlluminaBasecallsToMultiplexSam')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # RUN_DIR Required
            runnable_step.add_picard_option(key='RUN_DIR', value=self.run_directory)
            # LANE Required
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # OUTPUT Required
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.sorted_bam)
            # SEQUENCING_CENTER [BI]
            runnable_step.add_picard_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)
            # NUM_PROCESSORS [0]
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # FIRST_TILE [null]
            # TILE_LIMIT [null]
            # FORCE_GC [true]
            # APPLY_EAMSS_FILTER [true]
            if self.eamss_filter:
                runnable_step.add_picard_option(key='APPLY_EAMSS_FILTER', value='true')
            else:
                runnable_step.add_picard_option(key='APPLY_EAMSS_FILTER', value='false')
            # MAX_READS_IN_RAM_PER_TILE [1200000]
            # MINIMUM_QUALITY [2]
            # INCLUDE_NON_PF_READS [true]
            if self.vendor_quality_filter:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
            else:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # VERBOSITY [INFO]
            # QUIET [false]
            # VALIDATION_STRINGENCY [STRICT]
            # COMPRESSION_LEVEL [5]
            if self.compression_level is not None:
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value=str(self.compression_level))
            # MAX_RECORDS_IN_RAM [500000]
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # REFERENCE_SEQUENCE [null]
            # GA4GH_CLIENT_SECRETS [client_secrets.json]
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # OPTIONS_FILE Required
            # DEFLATER_THREADS [0]
            # NOTE: Only available in version: 2.18.24-CeMM
            runnable_step.add_picard_option(key='DEFLATER_THREADS', value='16')

            # Create the experiment directory if it does not exist already.

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move and rename the final, sorted BAM file.

            runnable_step = bsf.process.RunnableStepMove(
                name='move_sorted_bam',
                source_path=file_path_lane.sorted_bam,
                target_path=file_path_lane.archive_bam)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move and rename the checksum file.

            runnable_step = bsf.process.RunnableStepMove(
                name='move_sorted_md5',
                source_path=file_path_lane.sorted_md5,
                target_path=file_path_lane.archive_md5)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Upload the archive BAM file into the block blob storage.

            if self.cloud_account and self.cloud_container:
                runnable_cloud = self.add_runnable_consecutive(
                    runnable=bsf.procedure.ConsecutiveRunnable(
                        name=self.get_prefix_cloud(project_name=self.project_name, lane=lane_str),
                        code_module='bsf.runnables.generic',
                        working_directory=self.project_directory))
                executable_cloud = self.set_stage_runnable(stage=stage_cloud, runnable=runnable_cloud)
                executable_cloud.dependencies.append(executable_lane.name)

                # Upload the unaligned BAM file from its archive location.
                runnable_step = bsf.executables.cloud.RunnableStepAzureBlockBlobUpload(
                    name='blob_upload_bam',
                    account_name=self.cloud_account,
                    container_name=self.cloud_container,
                    source_path=file_path_lane.archive_bam,
                    target_path=file_path_lane.cloud_bam)
                runnable_cloud.add_runnable_step(runnable_step=runnable_step)

                # Upload the MD5 checksum file from its archive location.
                runnable_step = bsf.executables.cloud.RunnableStepAzureBlockBlobUpload(
                    name='blob_upload_md5',
                    account_name=self.cloud_account,
                    container_name=self.cloud_container,
                    source_path=file_path_lane.archive_md5,
                    target_path=file_path_lane.cloud_md5)
                runnable_cloud.add_runnable_step(runnable_step=runnable_step)

        # Add another flow cell-specific Runnable to reset directory and file mode permissions if requested.

        runnable_cell = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))

        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        runnable_step = bsf.process.RunnableStepChangeMode(
            name='chmod',
            file_path=experiment_directory,
            mode_directory=self.mode_directory,
            mode_file=self.mode_file)
        runnable_cell.add_runnable_step(runnable_step=runnable_step)

        return


class FilePathIlluminaDemultiplexSamCell(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell} models flow cell-specific file paths.

    Attributes:
    @ivar prefix_cell: Non-standard, flow cell-specific (i.e. project_name) prefix
    @type prefix_cell: str
    @ivar sample_annotation_sheet_csv: Sample Annotation Sheet CSV file
    @type sample_annotation_sheet_csv: str
    @see: FilePathBamIndexDecoderCell
    @see: FilePathExtractIlluminaCell
    """

    def __init__(self, prefix, project_name):
        """Initialise a C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell} object.

        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @return:
        @rtype:
        """
        super(FilePathIlluminaDemultiplexSamCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the flow cell-specific prefix accordingly.

        self.prefix_cell = project_name

        self.sample_annotation_sheet_csv = self.prefix_cell + '_samples.csv'

        return


class FilePathIlluminaDemultiplexSamLane(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane} models lane-specific file paths.

    Attributes:
    @ivar prefix_lane: Non-standard, lane-specific (i.e. project_name and lane) prefix
    @type prefix_lane: str
    @ivar project_barcode: Project-specific barcode CSV file
    @type project_barcode: str
    @ivar library_tsv: Library TSV file for Picard I{IlluminaBamDemux}
    @type library_tsv: str
    @ivar metrics_tsv: Lane-specific metrics TSV file
    @type metrics_tsv: str
    @ivar metrics_fraction_pdf: Lane-specific Picard I{IlluminaBamDemux} fraction metrics PDF file
    @type metrics_fraction_pdf: str
    @ivar metrics_fraction_png: Lane-specific Picard I{IlluminaBamDemux} fraction metrics PNG file
    @type metrics_fraction_png: str
    @ivar metrics_number_pdf: Lane-specific Picard I{IlluminaBamDemux} number metrics PDF file
    @type metrics_number_pdf: str
    @ivar metrics_number_png: Lane-specific Picard I{IlluminaBamDemux} number metrics PNG file
    @type metrics_number_png: str
    @ivar samples_directory: Sample directory
    @type samples_directory: str
    @ivar archive_bam: Archive BAM file
    @type archive_bam: str
    @ivar archive_md5: Archive BAM MD5 check sum
    @type archive_md5: str
    @see: FilePathBamIndexDecoderLane
    @see: FilePathExtractIlluminaLane
    """

    def __init__(self, prefix, project_name, lane, sequences_directory):
        """Initialise a C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane} object.

        @param prefix: Prefix
        @type prefix: str
        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str
        @return:
        @rtype:
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


class IlluminaDemultiplexSamSheet(bsf.annotation.AnnotationSheet):
    """The C{bsf.analyses.picard.IlluminaDemultiplexSamSheet} class represents a
    Tab-Separated Value (TSV) table of library information for the
    C{bsf.analyses.picard.IlluminaDemultiplexSam} C{bsf.analysis.Analysis}.

    Attributes:
    """

    _file_type = 'excel-tab'

    _header_line = True

    # FIXME: This class is only required, because Picard IlluminaBamDemux uses SAMPLE_NAME rather than SAMPLE_ALIAS.
    # Once fixed, the IlluminaBasecallsToSamSheet class could be used.
    _field_names = [
        'OUTPUT',
        'SAMPLE_NAME',
        'LIBRARY_NAME',
        'BARCODE_1',
        'BARCODE_2',
    ]

    _test_methods = dict()

    def adjust(self, barcode_length_tuple):
        """Adjust a C{bsf.analyses.picard.IlluminaDemultiplexSamSheet}.

        Remove BARCODE_N index columns, which annotation is all empty.
        @param barcode_length_tuple: Python C{tuple} of Python C{int} (barcode 1 length) and
            Python C{int} (barcode 2 length) objects
        @type barcode_length_tuple: (int, int)
        @return:
        @rtype:
        """
        # Adjust the IlluminaDemultiplexSamSheet and remove any BARCODE_N columns not represented
        # in the barcode length list.
        for index in range(0, 1 + 1):
            if barcode_length_tuple[index] == 0:
                # Remove the 'BARCODE_N' field from the list of field names.
                if 'BARCODE_' + str(index + 1) in self.field_names:
                    self.field_names.remove('BARCODE_' + str(index + 1))
                # Remove the 'BARCODE_N' entry from each row_dict object, since csv.DictWriter requires it.
                for row_dict in self.row_dicts:
                    row_dict.pop('BARCODE_' + str(index + 1), None)

        return


class IlluminaDemultiplexSam(bsf.analysis.Analysis):
    """The C{bsf.analyses.picard.IlluminaDemultiplexSam} class represents the logic to
    decode sequence archive BAM files into sample-specific BAM files via the Picard I{IlluminaSamDemux} tool.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar library_path: Library annotation file path
    @type library_path: str | None
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | None
    @ivar sequences_directory: BSF sequences directory
    @type sequences_directory: str | None
    @ivar samples_directory: BSF samples directory
    @type samples_directory: str | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | None
    @ivar compression_level: (Zlib) Compression level
    @type compression_level: int | None
    @ivar lanes: Number of lanes on the flow cell
    @type lanes: int | None
    @ivar force: Force de-multiplexing with a Library Annotation sheet failing validation
    @type force: bool | None
    """

    name = 'Picard IlluminaDemultiplexSam Analysis'
    prefix = 'illumina_demultiplex_sam'

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
        """Get a C{FilePathIlluminaDemultiplexSamCell} object.

        @param project_name: Project name
        @type project_name: str
        @return: C{FilePathIlluminaDemultiplexSamCell} object
        @rtype: FilePathIlluminaDemultiplexSamCell
        """
        return FilePathIlluminaDemultiplexSamCell(
            prefix=cls.get_prefix_cell(project_name=project_name),
            project_name=project_name)

    @classmethod
    def get_file_path_lane(cls, project_name, lane, sequences_directory):
        """Get a C{FilePathIlluminaDemultiplexSamLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str
        @return: C{FilePathIlluminaDemultiplexSamLane} object
        @rtype: FilePathIlluminaDemultiplexSamLane
        """
        return FilePathIlluminaDemultiplexSamLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane),
            project_name=project_name,
            lane=lane,
            sequences_directory=sequences_directory)

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
            library_path=None,
            run_directory=None,
            sequences_directory=None,
            samples_directory=None,
            mode_directory=None,
            mode_file=None,
            classpath_picard=None,
            compression_level=None,
            lanes=None,
            lane_list=None,
            force=None):
        """Initialise a C{bsf.analyses.picard.IlluminaDemultiplexSam} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
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
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param library_path: Library annotation file path
        @type library_path: str | None
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | None
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str | None
        @param samples_directory: BSF samples directory
        @type samples_directory: str | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | None
        @param compression_level: (Zlib) Compression level
        @type compression_level: int | None
        @param lanes: Number of lanes on the flow cell
        @type lanes: int | None
        @param lane_list: Python C{list} of Python C{str} lane number objects to process.
        @type lane_list: list[str]
        @param force: Force de-multiplexing with a Library Annotation sheet failing validation
        @type force: bool | None
        @return:
        @rtype:
        """
        super(IlluminaDemultiplexSam, self).__init__(
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
        self.classpath_picard = classpath_picard
        self.compression_level = compression_level
        self.lanes = lanes
        self.lane_list = lane_list
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
        """Set instance variables of a C{bsf.analyses.picard.IlluminaDemultiplexSam} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
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

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'compression_level'
        if configuration.config_parser.has_option(section=section, option=option):
            self.compression_level = configuration.config_parser.getint(section=section, option=option)

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

    def run(self):
        """Run the C{bsf.analyses.picard.IlluminaDemultiplexSam} analysis to
        decode an archive BAM file produced with IlluminaMultiplexSam into sample-specific BAM files.

        The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        a Picard tools *tab-separated* value (TSV) sample sheet.
        lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        barcode_sequence, barcode_name, library_name, sample_name, description
        @return:
        @rtype:
        """

        def run_get_sample_file_name(sample_name):
            """Private function to format sample-specific BAM file names (i.e. project_lane#sample.bam).

            @param sample_name:
            @type sample_name: str
            @return:
            @rtype: str
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # Start of the run() method body.

        # Configuration variable to use a bsf.ngs.Collection object rather than a
        # bsf.annotation.AnnotationSheet object to write the flow cell-specific
        # sample annotation sheet.

        use_collection = False

        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception('The ' + self.name + " requires a 'run_directory' configuration option.")

        # Illumina Run Folders are kept either in the 'active' or the 'archive' location.
        # Upon re-running of the de-multiplexing stage, the IRF may have been archived already.

        self.run_directory = bsf.illumina.RunFolder.absolute_file_path(name=self.run_directory)

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception('The ' + self.name + " 'run_directory' " + repr(self.run_directory) +
                            ' is not a valid directory.')

        irf = bsf.illumina.RunFolder.from_file_path(file_path=self.run_directory)

        if not self.project_name:
            self.project_name = '_'.join((irf.run_parameters.get_experiment_name, irf.run_information.flow_cell))

        # Only call the super-class method once the project name is defined with certainty.

        super(IlluminaDemultiplexSam, self).run()

        # Define the sequences and samples directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.sequences_directory:
            self.sequences_directory = self.project_name

        self.sequences_directory = self.configuration.get_absolute_path(
            file_path=self.sequences_directory,
            default_path=bsf.standards.FilePath.get_sequences(absolute=True))

        self.samples_directory = self.configuration.get_absolute_path(
            file_path=self.samples_directory,
            default_path=bsf.standards.FilePath.get_samples(absolute=True))

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception('The ' + self.name + " 'samples_directory' " + repr(self.samples_directory) +
                            ' is not a valid directory.')

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the library annotation sheet.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(self.library_path):
            raise Exception('The ' + self.name + " 'library_path' " + repr(self.library_path) +
                            ' is not a valid file.')

        # Load the LibraryAnnotationSheet and validate.

        library_annotation_sheet = bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)

        if self.lanes is None:
            validation_messages = library_annotation_sheet.validate(
                lanes=irf.run_information.flow_cell_layout.lane_count)
        else:
            validation_messages = library_annotation_sheet.validate(
                lanes=self.lanes)

        if validation_messages:
            if self.force:
                warnings.warn('Validation of library annotation sheet ' + repr(self.library_path) + ':\n' +
                              validation_messages)
            else:
                raise Exception('Validation of library annotation sheet ' + repr(self.library_path) + ':\n' +
                                validation_messages)

        library_annotation_dict = library_annotation_sheet.get_annotation_dict()
        library_barcode_dict = library_annotation_sheet.get_barcode_length_dict()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('The ' + self.name + " requires a 'classpath_picard' configuration option.")

        stage_lane = self.get_stage(name=self.get_stage_name_lane())
        stage_cell = self.get_stage(name=self.get_stage_name_cell())

        file_path_cell = self.get_file_path_cell(project_name=self.project_name)

        # NOTE: Use a SampleAnnotationSheet as long as the Collection does not have a method to
        # write a de-normalised table.
        # Create a SampleAnnotationSheet in the project directory and
        # eventually transfer it into the experiment_directory.

        sample_annotation_sheet = bsf.analyses.illumina_to_bam_tools.BamIndexDecoder.get_sample_annotation_sheet(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        # Create a Collection and ProcessedRunFolder object to write a SampleAnnotationSheet.

        collection = bsf.ngs.Collection(
            name='Samples',
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        prf = bsf.ngs.ProcessedRunFolder(name=bsf.ngs.ProcessedRunFolder.default_name)

        collection.add_processed_run_folder(prf=prf)

        cell_dependency_list = list()
        """ @type cell_dependency_list: list[str] """

        for lane_int in sorted(library_annotation_dict):
            lane_str = str(lane_int)

            # If a list of lanes was defined, ignore lanes not on the list.
            if self.lane_list and lane_str not in self.lane_list:
                continue

            file_path_lane = self.get_file_path_lane(
                project_name=self.project_name,
                lane=lane_str,
                sequences_directory=self.sequences_directory)

            # experiment_samples_directory = os.path.join(experiment_directory, file_path_lane.samples_directory)

            # LIBRARY_PARAMS
            ibs_sheet = IlluminaDemultiplexSamSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.library_tsv))

            for row_dict in library_annotation_dict[lane_int]:
                # Add a row to the lane-specific Picard IlluminaDemultiplexSamSheet.

                ibs_sheet.row_dicts.append({
                    'BARCODE_1': row_dict['barcode_sequence_1'],
                    'BARCODE_2': row_dict['barcode_sequence_2'],
                    'OUTPUT': os.path.join(
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    # FIXME: IlluminaBamDemux does not seem standard and requires SAMPLE_NAME rather than SAMPLE_ALIAS?
                    # 'SAMPLE_ALIAS': row_dict['sample_name'],
                    'SAMPLE_NAME': row_dict['sample_name'],
                    'LIBRARY_NAME': row_dict['library_name'],
                })

                # Get or create a (default) NGS Project.

                if bsf.ngs.Sample.default_name in prf.project_dict:
                    project = prf.project_dict[bsf.ngs.Sample.default_name]
                else:
                    project = bsf.ngs.Project(name=bsf.ngs.Sample.default_name)
                    prf.add_project(project=project)

                # Get or create a NGS Sample.

                if row_dict['sample_name'] in project.sample_dict:
                    sample = project.sample_dict[row_dict['sample_name']]
                else:
                    sample = bsf.ngs.Sample(name=row_dict['sample_name'])
                    project.add_sample(sample=sample)

                sample.add_annotation(key='Project Name', value=row_dict['library_name'])
                sample.add_annotation(key='Project Size', value=row_dict['library_size'])

                # Create one NGS Reads object per row.

                reads_1 = bsf.ngs.Reads(
                    name='_'.join((self.project_name, lane_str, row_dict['sample_name'])),
                    file_path=os.path.join(
                        os.path.basename(experiment_directory),
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])))

                # Create one NGS PairedReads object per NGS Reads object.

                paired_reads = bsf.ngs.PairedReads(
                    reads_1=reads_1,
                    index_1=row_dict['barcode_sequence_1'],
                    index_2=row_dict['barcode_sequence_2'])
                paired_reads.add_annotation(key='Flow Cell', value=self.project_name)
                paired_reads.add_annotation(key='Lane', value=lane_str)
                paired_reads.add_annotation(key='Flow Cell Lane', value='_'.join((self.project_name, lane_str)))
                paired_reads.add_annotation(key='Structure', value=irf.run_information.get_picard_read_structure)

                # Finally, add the NGS PairedReads object to the NGS Sample object.

                sample.add_paired_reads(paired_reads=paired_reads)

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
                    'PairedReads ReadGroup': '',
                    'PairedReads Structure': irf.run_information.get_picard_read_structure,
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

            runnable_lane = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)
            executable_lane.dependencies.append(
                IlluminaMultiplexSam.get_prefix_lane(
                    project_name=self.project_name,
                    lane=lane_str))

            cell_dependency_list.append(executable_lane.name)

            if executable_lane.submit:
                # Only if this Executable actually gets submitted ...
                # Write the lane-specific BamIndexDecoderSheet to the internal file path.

                ibs_sheet.to_file_path()

            # Create a samples_directory in the project_directory.

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_project_lane_samples_directory',
                directory_path=file_path_lane.samples_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Run Picard IlluminaSamDemux.

            runnable_step = bsf.process.RunnableStepPicard(
                name='picard_illumina_demultiplex_sam',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx5G',
                picard_classpath=self.classpath_picard,
                picard_command='IlluminaSamDemux')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # INPUT Required
            runnable_step.add_picard_option(key='INPUT', value=file_path_lane.archive_bam)
            # OUTPUT_DIR [null] Defaults to the same directory as INPUT.
            runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.samples_directory)
            # OUTPUT_PREFIX [null] Defaults to the @RG ID attribute of the multiplexed BAM file
            runnable_step.add_picard_option(key='OUTPUT_PREFIX', value='_'.join((self.project_name, lane_str)))
            # OUTPUT_FORMAT [bam]
            # BARCODE_TAG_NAME [BC]
            # BARCODE_QUALITY_TAG_NAME [QT]
            # LIBRARY_PARAMS Required
            runnable_step.add_picard_option(key='LIBRARY_PARAMS', value=file_path_lane.library_tsv)
            # METRICS_FILE Required
            runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
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
            # MIN_MISMATCH_DELTA [1]
            # MAX_NO_CALLS [2]
            # MINIMUM_BASE_QUALITY [0]
            # READ_STRUCTURE [null]
            if 'read_structure' in library_annotation_dict[lane_int][0] and \
                    library_annotation_dict[lane_int][0]['read_structure']:
                runnable_step.add_picard_option(
                    key='READ_STRUCTURE',
                    value=library_annotation_dict[lane_int][0]['read_structure'])
            else:
                runnable_step.add_picard_option(
                    key='READ_STRUCTURE',
                    value=irf.run_information.get_picard_read_structure)
            # ADAPTERS_TO_CHECK [INDEXED, DUAL_INDEXED, NEXTERA_V2, FLUIDIGM]
            # FIVE_PRIME_ADAPTER [null]
            # THREE_PRIME_ADAPTER [null]
            # MOLECULAR_INDEX_TAG [RX]
            # MOLECULAR_INDEX_BASE_QUALITY_TAG [QX]
            # TAG_PER_MOLECULAR_INDEX [null]
            # CELL_INDEX_TAG [CB]
            # CELL_INDEX_BASE_QUALITY_TAG [CY]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))
            # VERBOSITY [INFO]
            # QUIET [false]
            # VALIDATION_STRINGENCY [STRICT]
            # COMPRESSION_LEVEL [5]
            if self.compression_level is not None:
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value=str(self.compression_level))
            # MAX_RECORDS_IN_RAM [500000]
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # REFERENCE_SEQUENCE [null]
            # GA4GH_CLIENT_SECRETS [client_secrets.json]
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # OPTIONS_FILE Required
            # DEFLATER_THREADS [0]
            # NOTE: Only available in version: 2.18.24-CeMM
            runnable_step.add_picard_option(key='DEFLATER_THREADS', value='16')

            # Plot the metrics file.

            runnable_step = bsf.process.RunnableStep(
                name='plot_metrics',
                program='bsf_illumina_demultiplex_sam.R')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_long(key='file-path', value=file_path_lane.metrics_tsv)

            # Create the experiment directory if it does not exist already.

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_experiment_directory',
                directory_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move the samples directory into the experiment directory.

            runnable_step = bsf.process.RunnableStepMove(
                name='move_project_lane_samples_directory',
                source_path=file_path_lane.samples_directory,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # Move the metrics file into the experiment directory.

            runnable_step = bsf.process.RunnableStepMove(
                name='move_metrics_tsv',
                source_path=file_path_lane.metrics_tsv,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = bsf.process.RunnableStepMove(
                name='move_metrics_pdf',
                source_path=file_path_lane.metrics_number_pdf,
                target_path=experiment_directory)
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = bsf.process.RunnableStepMove(
                name='move_metrics_png',
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

        runnable_cell = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        # Move the SampleAnnotationSheet from the project_directory to the experiment_directory.

        if os.path.exists(collection.file_path):
            runnable_step = bsf.process.RunnableStepMove(
                name='move_sample_annotation',
                source_path=file_path_cell.sample_annotation_sheet_csv,
                target_path=experiment_directory)
            runnable_cell.add_runnable_step(runnable_step=runnable_step)

        # Change directory and file access permissions.

        runnable_step = bsf.process.RunnableStepChangeMode(
            name='chmod',
            file_path=experiment_directory,
            mode_directory=self.mode_directory,
            mode_file=self.mode_file)
        runnable_cell.add_runnable_step(runnable_step=runnable_step)

        return


class FilePathCollectHiSeqXPfFailMetricsLane(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane} models files in a directory.

    Attributes:
    @ivar summary_tsv: Summary metrics TSV file
    @type summary_tsv: str
    @ivar detailed_tsv: Detailed metrics TSV file
    @type detailed_tsv: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane} object.

        @param prefix: Prefix
        @type prefix: str
        @return:
        @rtype:
        """
        super(FilePathCollectHiSeqXPfFailMetricsLane, self).__init__(prefix=prefix)

        self.summary_tsv = prefix + '.pffail_summary_metrics'
        self.detailed_tsv = prefix + '.pffail_detailed_metrics'

        return


class CollectHiSeqXPfFailMetrics(PicardIlluminaRunFolder):
    """The C{bsf.analyses.picard.CollectHiSeqXPfFailMetrics} class represents Picard CollectHiSeqXPfFailMetrics.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    """

    name = 'Picard CollectHiSeqXPfFailMetrics Analysis'
    prefix = 'picard_hiseq_x_pf_fail'

    @classmethod
    def get_file_path_lane(cls, project_name, lane):
        """Get a C{FilePathCollectHiSeqXPfFailMetricsLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @return: C{FilePathCollectHiSeqXPfFailMetricsLane} object
        @rtype: FilePathCollectHiSeqXPfFailMetricsLane
        """
        return FilePathCollectHiSeqXPfFailMetricsLane(
            prefix=cls.get_prefix_lane(project_name=project_name, lane=lane))

    def run(self):
        """Run the C{bsf.analyses.picard.CollectHiSeqXPfFailMetrics} C{bsf.analysis.Analysis}.

        @return:
        @rtype:
        """
        super(CollectHiSeqXPfFailMetrics, self).run()

        # Picard CollectHiSeqXPfFailMetrics

        stage_lane = self.get_stage(name=self.get_stage_name_lane())

        cell_dependency_list = list()

        for lane_int in range(0 + 1, self._irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            file_path_lane = self.get_file_path_lane(project_name=self.project_name, lane=lane_str)

            runnable_lane = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_lane(
                        project_name=self.project_name,
                        lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)

            # Add the dependency for the cell-specific process.

            cell_dependency_list.append(executable_lane.name)

            runnable_step = bsf.process.RunnableStepPicard(
                name='collect_hiseq_x_fail_metrics',
                java_temporary_path=runnable_lane.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx4G',
                picard_classpath=self.classpath_picard,
                picard_command='CollectHiSeqXPfFailMetrics')
            runnable_lane.add_runnable_step(runnable_step=runnable_step)

            # BASECALLS_DIR Required
            runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
            # OUTPUT Required
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.prefix)
            # PROB_EXPLICIT_READS [0.0]
            # LANE Required
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # NUM_PROCESSORS [1]
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # N_CYCLES [24] Should match Illumina RTA software.
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_lane.temporary_directory_path(absolute=False))

        return


class FilePathDownsampleSamReadGroup(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathDownsampleSamReadGroup} models files in a directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar downsampled_bam: Down-sampled BAM file
    @type downsampled_bam: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathDownsampleSamReadGroup} object.

        @param prefix: Prefix
        @type prefix: str
        @return:
        @rtype:
        """
        super(FilePathDownsampleSamReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.downsampled_bam = prefix + '_downsampled.bam'

        return


class DownsampleSam(bsf.analysis.Analysis):
    """The C{bsf.analyses.picard.DownsampleSam} class represents the logic to run the Picard DownsampleSam analysis.

    Attributes:

    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | None
    """

    name = 'Picard DownsampleSam Analysis'
    prefix = 'picard_downsample_sam'

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name):
        """Get a C{FilePathDownsampleSamReadGroup} object.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: C{FilePathDownsampleSamReadGroup} object
        @rtype: FilePathDownsampleSamReadGroup
        """
        return FilePathDownsampleSamReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

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
            classpath_picard=None):
        """Initialise a C{bsf.analyses.picard.DownsampleSam} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
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
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | None
        @return:
        @rtype:
        """
        super(DownsampleSam, self).__init__(
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

        self.classpath_picard = classpath_picard

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.picard.DownsampleSam} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(DownsampleSam, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.picard.DownsampleSam} C{bsf.analysis.Analysis}.

        This method changes the C{bsf.ngs.Collection} object of this C{bsf.analysis.Analysis}
        to update with FASTQ file paths.
        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.analysis.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.analysis.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        super(DownsampleSam, self).run()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception("A 'DownsampleSam' analysis requires a "
                                "'classpath_picard' configuration option.")

        run_read_comparisons()

        # Picard DownsampleSam

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    # Apply some sanity checks.

                    if paired_reads.reads_1 is None and paired_reads.reads_2 is not None:
                        raise Exception('PairedReads object with a reads_2, but no reads_1 object.')

                    reads = paired_reads.reads_1
                    if not (reads.file_path.endswith('.bam') or reads.file_path.endswith('.sam')):
                        raise Exception(
                            'Picard DownsampleSam can only work on BAM or SAM files. ' + reads.file_path)

                    prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    # Keep the original BAM file and modify the file_path in the Reads object.
                    bam_file_path = reads.file_path
                    reads.file_path = os.path.join(self.project_directory, file_path_read_group.downsampled_bam)
                    # (file_root, file_extension) = os.path.splitext(bam_file_path)
                    # reads.file_path = file_root + 'downsampled' + file_extension

                    # Create a Runnable for running the Picard DownsampleSam analysis.

                    runnable_picard_dss = self.add_runnable_consecutive(
                        runnable=bsf.procedure.ConsecutiveRunnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory))

                    # Create an Executable for running the Picard SamToFastq Runnable.

                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_picard_dss)

                    # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                    # runnable_step = bsf.process.RunnableStepMakeDirectory(
                    #    name='make_directory',
                    #    directory_path=file_path_read_group.output_directory)
                    # runnable_picard_dss.add_runnable_step(runnable_step=runnable_step)

                    # Create a new RunnableStep for the Picard DownsampleSam program.

                    runnable_step = bsf.process.RunnableStepPicard(
                        name='picard_downsample_sam',
                        java_temporary_path=runnable_picard_dss.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx2G',
                        picard_classpath=self.classpath_picard,
                        picard_command='DownsampleSam')
                    runnable_picard_dss.add_runnable_step(runnable_step=runnable_step)

                    # INPUT Required
                    runnable_step.add_picard_option(key='INPUT', value=bam_file_path)
                    # OUTPUT Required
                    runnable_step.add_picard_option(
                        key='OUTPUT',
                        value=file_path_read_group.downsampled_bam)
                    # FIXME: Add to the configuration file and documentation.
                    # STRATEGY [ConstantMemory]
                    # RANDOM_SEED [1]
                    # PROBABILITY [1.0]
                    if 'DownsampleSam Probability' in paired_reads.annotation_dict:
                        runnable_step.add_picard_option(
                            key='PROBABILITY',
                            value=paired_reads.annotation_dict['DownsampleSam Probability'][0])
                    # ACCURACY [1.0E-4]
                    # METRICS_FILE [null]
                    # TMP_DIR [null]
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_picard_dss.temporary_directory_path(absolute=False))
                    # VERBOSITY [INFO]
                    # QUIET [false]
                    # VALIDATION_STRINGENCY [STRICT]
                    # COMPRESSION_LEVEL [5]
                    # MAX_RECORDS_IN_RAM [500000]
                    # CREATE_INDEX [false]
                    # CREATE_MD5_FILE [false]
                    # REFERENCE_SEQUENCE [null]
                    # GA4GH_CLIENT_SECRETS [client_secrets.json]
                    # USE_JDK_DEFLATER [false]
                    # USE_JDK_INFLATER [false]
                    # OPTIONS_FILE Required

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'picard_downsample_sam_samples.csv'))),
            name='_'.join((self.project_name, 'picard_downsample_sam')))
        annotation_sheet.to_file_path()

        return


class FilePathSamToFastqReadGroup(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathSamToFastqReadGroup} class models read group-specific Picard SamToFastq files.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathSamToFastqReadGroup} object.

        @param prefix: Prefix
        @type prefix: str
        @return:
        @rtype:
        """
        super(FilePathSamToFastqReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class FilePathSamToFastqProject(bsf.procedure.FilePath):
    """The C{bsf.analyses.picard.FilePathSamToFastqProject} class models project-specific Picard SamToFastq files.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar sas_path_old: Old Sample Annotation Sheet file path
    @type sas_path_old: str
    @ivar sas_path_new: New Sample Annotation Sheet file path
    @type sas_path_new: str
    """

    def __init__(self, prefix, prefix_analysis, project_name):
        """Initialise a C{bsf.analyses.picard.FilePathSamToFastqProject} object.

        @param prefix: Prefix
        @type prefix: str
        @param prefix_analysis: C{bsf.analysis.Analysis.prefix}
        @type prefix_analysis: str
        @param project_name: Project name
        @type project_name: str
        @return:
        @rtype:
        """
        super(FilePathSamToFastqProject, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sas_path_old = '_'.join((project_name, prefix_analysis, 'original.csv'))
        self.sas_path_new = '_'.join((project_name, prefix_analysis, 'samples.csv'))

        return


class SamToFastq(bsf.analysis.Analysis):
    """The C{bsf.analyses.picard.SamToFastq} class represents the logic to run the Picard SamToFastq analysis.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_read_group: C{bsf.analysis.Stage.name} for read group-specific C{bsf.procedure.Runnable} objects
    @type stage_name_read_group: str
    @cvar stage_name_project: C{bsf.analysis.Stage.name} for project-specific C{bsf.procedure.Runnable} objects
    @type stage_name_project: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | None
    @ivar include_non_pf_reads: Include non-pass filer reads
    @type include_non_pf_reads: bool | None
    @ivar drop_read_1: Drop read 1
    @type drop_read_1: bool | None
    @ivar drop_read_2: Drop read 2
    @type drop_read_2: bool | None
    """

    name = 'Picard SamToFastq Analysis'
    prefix = 'picard_sam_to_fastq'

    stage_name_read_group = '_'.join((prefix, 'read_group'))
    stage_name_project = '_'.join((prefix, 'project'))

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_stage_name_project(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'project'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_prefix_project(cls, project_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: Project name
        @type project_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_project(), project_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name):
        """Get a C{FilePathSamToFastqReadGroup} object.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: C{FilePathSamToFastqReadGroup} object
        @rtype: FilePathSamToFastqReadGroup
        """
        return FilePathSamToFastqReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

    @classmethod
    def get_file_path_project(cls, project_name, prefix_analysis):
        """Get a C{FilePathSamToFastqProject} object.

        @param project_name: Project name
        @type project_name: str
        @param prefix_analysis: C{bsf.analysis.Analysis.prefix}
        @type prefix_analysis: str
        @return: C{FilePathSamToFastqProject} object
        @rtype: FilePathSamToFastqProject
        """
        return FilePathSamToFastqProject(
            prefix=cls.get_prefix_project(project_name=project_name),
            prefix_analysis=prefix_analysis,
            project_name=project_name)

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
            classpath_picard=None,
            include_non_pf_reads=None,
            drop_read_1=None,
            drop_read_2=None):
        """Initialise a C{bsf.analyses.picard.SamToFastq} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
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
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | None
        @param include_non_pf_reads: Include non-pass filer reads
        @type include_non_pf_reads: bool | None
        @param drop_read_1: Drop read 1
        @type drop_read_1: bool | None
        @param drop_read_2: Drop read 2
        @type drop_read_2: bool | None
        @return:
        @rtype:
        """
        super(SamToFastq, self).__init__(
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

        self.classpath_picard = classpath_picard
        self.include_non_pass_filter_reads = include_non_pf_reads
        self.drop_read_1 = drop_read_1
        self.drop_read_2 = drop_read_2

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.picard.SamToFastq} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(SamToFastq, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

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

    def _read_comparisons(self):
        """Private method to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

        This implementation just adds all C{bsf.ngs.Sample} objects from the
        C{bsf.analysis.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
        C{bsf.analysis.Analysis.sample_list} instance variable.
        @return:
        @rtype:
        """

        self.sample_list.extend(self.collection.get_all_samples())

        return

    def run(self):
        """Run the C{bsf.analyses.picard.SamToFastq} C{bsf.analysis.Analysis}.

        This method changes the C{bsf.ngs.Collection} object of this C{bsf.analysis.Analysis}
        to update with FASTQ file paths.
        @return:
        @rtype:
        """

        # Start of the run() method body.

        super(SamToFastq, self).run()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('The ' + self.name + " requires a 'classpath_picard' configuration option.")

        self._read_comparisons()

        # Picard SamToFastq

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())
        stage_project = self.get_stage(name=self.get_stage_name_project())

        project_dependency_list = list()
        """ @type project_dependency_list: list[str] """

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    # Apply some sanity checks.

                    if paired_reads.reads_1 is None and paired_reads.reads_2 is not None:
                        raise Exception('PairedReads object with a reads_2, but no reads_2 object.')

                    reads = paired_reads.reads_1
                    if reads.file_path.endswith('.bam'):
                        bam_file_path = reads.file_path
                        prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                        file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                        # Get the SAM header of a BAM file to extract the read group (@RG), amongst other things.

                        # Open the BAM file, while not checking sequence (@SQ) entries.
                        # De-multiplexed, unaligned BAM files have no reference sequence dictionary.

                        alignment_file = pysam.AlignmentFile(reads.file_path, 'rb', check_sq=False)

                        for read_group_dict in alignment_file.header['RG']:
                            """ @type read_group_dict: dict[str, str] """
                            # The makeFileNameSafe() method of htsjdk.samtools.util.IOUtil uses the following pattern:
                            # [\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]
                            platform_unit = re.sub(
                                pattern='[\\s!"#$%&\'()*/:;<=>?@\\[\\]\\\\^`{|}~]',
                                repl='_',
                                string=read_group_dict['PU'])
                            read_group_list = ['@RG']
                            read_group_list.extend([':'.join((key, value)) for key, value in read_group_dict.items()])
                            if read_group_dict == alignment_file.header['RG'][0]:
                                # Use the '==' rather than the 'is' operator, since dictionaries do not seem to be
                                # at the same memory address.
                                # For the first read group, modify the PairedReads object in place.
                                paired_reads.read_group = '\\t'.join(read_group_list)
                                paired_reads.reads_1.name = platform_unit + '_1'
                                paired_reads.reads_1.file_path = os.path.join(
                                    self.project_directory,
                                    file_path_read_group.output_directory,
                                    platform_unit + '_1.fastq.gz')
                                paired_reads.reads_2 = bsf.ngs.Reads(
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
                                reads_1 = bsf.ngs.Reads(
                                    name=platform_unit + '_1',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_read_group.output_directory,
                                        platform_unit + '_1.fastq.gz'),
                                    file_type=paired_reads.reads_1.file_type,
                                    lane=paired_reads.reads_1.lane,
                                    read='R1',
                                    chunk=paired_reads.reads_1.chunk)
                                reads_2 = bsf.ngs.Reads(
                                    name=platform_unit + '_2',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_read_group.output_directory,
                                        platform_unit + '_2.fastq.gz'),
                                    file_type=paired_reads.reads_1.file_type,
                                    lane=paired_reads.reads_1.lane,
                                    read='R2',
                                    chunk=paired_reads.reads_1.chunk)
                                new_paired_reads = bsf.ngs.PairedReads(
                                    reads_1=reads_1,
                                    reads_2=reads_2,
                                    read_group='\\t'.join(read_group_list))
                                # Retain the original BAM file path as annotation.
                                new_paired_reads.add_annotation(key='BAM File', value=bam_file_path)

                                sample.add_paired_reads(paired_reads=new_paired_reads)

                        alignment_file.close()

                        # Create a Runnable for running the Picard SamToFastq analysis.

                        runnable_read_group = self.add_runnable_consecutive(
                            runnable=bsf.procedure.ConsecutiveRunnable(
                                name=prefix_read_group,
                                code_module='bsf.runnables.generic',
                                working_directory=self.project_directory))
                        self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                        # Record the Executable.name for the project dependency.

                        project_dependency_list.append(runnable_read_group.name)

                        # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                        runnable_step = bsf.process.RunnableStepMakeDirectory(
                            name='mkdir',
                            directory_path=file_path_read_group.output_directory)
                        runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                        # Create a new RunnableStep for the Picard SamToFastq program.

                        runnable_step = bsf.process.RunnableStepPicard(
                            name='picard_sam_to_fastq',
                            java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx2G',
                            picard_classpath=self.classpath_picard,
                            picard_command='SamToFastq')
                        runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                        # INPUT Required
                        runnable_step.add_picard_option(key='INPUT', value=bam_file_path)
                        # FASTQ Required
                        # SECOND_END_FASTQ [null]
                        # UNPAIRED_FASTQ [null]
                        # OUTPUT_PER_RG [false]
                        runnable_step.add_picard_option(key='OUTPUT_PER_RG', value='true')
                        # COMPRESS_OUTPUTS_PER_RG [false]
                        runnable_step.add_picard_option(key='COMPRESS_OUTPUTS_PER_RG', value='true')
                        # RG_TAG [PU]
                        # OUTPUT_DIR [null]
                        runnable_step.add_picard_option(
                            key='OUTPUT_DIR',
                            value=file_path_read_group.output_directory)
                        # RE_REVERSE [true]
                        # INTERLEAVE [false]
                        # INCLUDE_NON_PF_READS [false]
                        if self.include_non_pass_filter_reads:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
                        else:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
                        # CLIPPING_ATTRIBUTE [null]
                        runnable_step.add_picard_option(key='CLIPPING_ATTRIBUTE', value='XT')
                        # CLIPPING_ACTION [null]
                        runnable_step.add_picard_option(key='CLIPPING_ACTION', value='N')
                        # CLIPPING_MIN_LENGTH [0]
                        # READ1_TRIM [0]
                        # READ1_MAX_BASES_TO_WRITE [null]
                        # READ2_TRIM [0]
                        # READ2_MAX_BASES_TO_WRITE [null]
                        # QUALITY [null]
                        # INCLUDE_NON_PRIMARY_ALIGNMENTS [false]
                        # TMP_DIR [null]
                        runnable_step.add_picard_option(
                            key='TMP_DIR',
                            value=runnable_read_group.temporary_directory_path(absolute=False))
                        # VERBOSITY [INFO]
                        runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                        # QUIET [false]
                        # VALIDATION_STRINGENCY [STRICT]
                        # COMPRESSION_LEVEL [5]
                        runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                        # MAX_RECORDS_IN_RAM [500000]
                        # CREATE_INDEX [false]
                        # CREATE_MD5_FILE [false]
                        runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                        # REFERENCE_SEQUENCE [null]
                        # GA4GH_CLIENT_SECRETS [client_secrets.json|
                        # USE_JDK_DEFLATER [false]
                        # USE_JDK_INFLATER [false]
                        # OPTIONS_FILE Required

        # Create a Runnable for pruning the sample annotation sheet.

        prefix_project = self.get_prefix_project(project_name=self.project_name)

        file_path_project = self.get_file_path_project(project_name=self.project_name, prefix_analysis=self.prefix)

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(self.project_directory, file_path_project.sas_path_old),
            name=prefix_project)
        annotation_sheet.to_file_path()

        runnable_project = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=prefix_project,
                code_module='bsf.runnables.picard_sam_to_fastq_sample_sheet',
                working_directory=self.project_directory))
        executable_project = self.set_stage_runnable(
            stage=stage_project,
            runnable=runnable_project)
        executable_project.dependencies.extend(project_dependency_list)

        # Create a new RunnableStep.

        runnable_step = bsf.process.RunnableStep(
            name='prune_sample_annotation_sheet')
        runnable_project.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_option_long(key='sas_path_old', value=file_path_project.sas_path_old)
        runnable_step.add_option_long(key='sas_path_new', value=file_path_project.sas_path_new)
        runnable_step.add_option_long(key='minimum_size', value='1024')
        if self.drop_read_1:
            runnable_step.add_option_long(key='drop_read_1', value='True')
        if self.drop_read_2:
            runnable_step.add_option_long(key='drop_read_2', value='True')

        return

    def prune(self):
        """Prune the Analysis by replacing FASTQ files with (empty) status files (*.truncated).

        The status files are recognised by the bsf.runnables.picard_sam_to_fastq_sample_sheet module
        so that Reads and PairedReads objects are kept in the sample annotation sheet.
        @return:
        @rtype:
        """
        super(SamToFastq, self).run()

        self._read_comparisons()

        if self.debug > 0:
            print('Prune!')

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    def prune_file_path(reads):
                        """Private function to prune files of a Reads object.

                        @param reads: Reads
                        @type reads: bsf.ngs.Reads
                        @return:
                        @rtype:
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

                        # Also remove an eventual MD5 sum file.
                        file_path += '.md5'
                        if os.path.exists(file_path):
                            os.remove(file_path)

                        return

                    if paired_reads.reads_1 is not None:
                        prune_file_path(reads=paired_reads.reads_1)

                    if paired_reads.reads_2 is not None:
                        prune_file_path(reads=paired_reads.reads_2)

        return
