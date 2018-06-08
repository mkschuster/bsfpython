"""bsf.analyses.picard

A package of classes and methods modelling Picard analyses data files and data directories.
"""

#
# Copyright 2013 - 2018 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF Python.
#
# BSF Python is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF Python is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function

import os
import re
import sys
import warnings
import weakref

import pysam

import bsf
import bsf.analyses.illumina_to_bam_tools
import bsf.annotation
import bsf.illumina
import bsf.ngs
import bsf.process
import bsf.standards


class PicardIlluminaRunFolder(bsf.Analysis):
    """The C{bsf.analyses.picard.PicardIlluminaRunFolder} class of Picard Analyses acting on Illumina Run Folders.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_lane: C{bsf.Stage.name} for the lane-specific stage
    @type stage_name_lane: str
    @cvar stage_name_cell: C{bsf.Stage.name} for the flow cell-specific stage
    @type stage_name_cell: str
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | unicode | None
    @ivar intensity_directory: File path to the I{Intensities} directory,
        defaults to I{illumina_run_folder/Data/Intensities}
    @type intensity_directory: str | unicode | None
    @ivar basecalls_directory: File path to the I{BaseCalls} directory,
        defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
    @type basecalls_directory: str | unicode | None
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode | None
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool | None
    """

    name = 'Picard PicardIlluminaRunFolder Analysis'
    prefix = 'picard_illumina_run_folder'

    stage_name_lane = '_'.join((prefix, 'lane'))
    stage_name_cell = '_'.join((prefix, 'cell'))

    @classmethod
    def get_prefix_cell(cls, project_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects.

        @param project_name: A project name
        @type project_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_cell, project_name))

    @classmethod
    def get_prefix_lane(cls, project_name, lane):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_lane, project_name, lane))

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
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode | None
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | unicode | None
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | unicode | None
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
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

        # Get Illumina Run Folder information.

        option = 'illumina_run_folder'
        if configuration.config_parser.has_option(section=section, option=option):
            self.run_directory = configuration.config_parser.get(section=section, option=option)

        option = 'intensity_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.intensity_directory = configuration.config_parser.get(section=section, option=option)

        option = 'basecalls_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.basecalls_directory = configuration.config_parser.get(section=section, option=option)

        # Get the experiment name.

        option = 'experiment_name'
        if configuration.config_parser.has_option(section=section, option=option):
            self.experiment_name = configuration.config_parser.get(section=section, option=option)

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.picard.PicardIlluminaRunFolder} C{bsf.Analysis}.

        @return:
        @rtype:
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
            default_path=bsf.standards.Default.absolute_illumina_run())

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception('The Illumina run directory {!r} does not exist.'.
                            format(self.run_directory))

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise bsf.illumina.RunFolderNotComplete('The Illumina run directory {!r} is not complete.'.
                                                    format(self.run_directory))

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
            raise Exception('The Intensity directory {!r} does not exist.'.
                            format(self.intensity_directory))

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
            raise Exception('The BaseCalls directory {!r} does not exist.'.
                            format(self.basecalls_directory))

        self._irf = bsf.illumina.RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g. BSF_0000) is used as the prefix for archive BAM files.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = self._irf.run_parameters.get_experiment_name
            if not self.experiment_name:
                raise Exception('Could not get an experiment_name from the Illumina Run Folder.')

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, self._irf.run_information.flow_cell))

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception("A 'PicardIlluminaRunFolder' analysis requires a "
                                "'classpath_picard' configuration option.")

        # Call the run method of the super class after the project_name has been defined.

        super(PicardIlluminaRunFolder, self).run()

        return


class ExtractIlluminaBarcodesSheet(bsf.annotation.AnnotationSheet):
    """The C{bsf.analyses.picard.ExtractIlluminaBarcodesSheet} class represents a
    Tab-Separated Value (TSV) table of library information for the
    C{bsf.analyses.picard.ExtractIlluminaBarcodes} C{bsf.Analysis}.

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
    C{bsf.analyses.picard.ExtractIlluminaBarcodes} C{bsf.Analysis}.

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

    def adjust(self, unassigned_file_path, index_read_number):
        """Adjust the C{bsf.analyses.picard.IlluminaBasecallsToSamSheet} by inserting the
        file path for unassigned reads and by adjusting the index columns to the number of
        index reads.

        @param unassigned_file_path: File path for unassigned reads
        @type unassigned_file_path: str | unicode
        @param index_read_number: Number of index reads
        @type index_read_number: int
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
            # ... in all other cases, since a last row for unmatched barcode sequences needs adding.
            self.row_dicts.append({
                'BARCODE_1': 'N',
                'BARCODE_2': '',
                'OUTPUT': unassigned_file_path,
                'SAMPLE_ALIAS': 'Unmatched',
                'LIBRARY_NAME': self.row_dicts[0]['LIBRARY_NAME'],
            })

        # Further adjust the IlluminaBaseCallsToSamSheet and remove any BARCODE_N columns not represented
        # in the read structure.
        for index in range(0, 1 + 1):
            if index + 1 > index_read_number:
                # Remove the 'BARCODE_N' field from the list of field names.
                if 'BARCODE_' + str(index + 1) in self.field_names:
                    self.field_names.remove('BARCODE_' + str(index + 1))
                # Remove the 'BARCODE_N' entry from each row_dict object, since csv.DictWriter requires it.
                for row_dict in self.row_dicts:
                    row_dict.pop('BARCODE_' + str(index + 1), None)

        return


class FilePathExtractIlluminaCell(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathExtractIlluminaCell} models files in a directory.

    Attributes:
    @ivar sample_annotation_sheet_csv: Sample Annotation Sheet CSV file
    @type sample_annotation_sheet_csv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathExtractIlluminaCell} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathExtractIlluminaCell, self).__init__(prefix=prefix)

        self.sample_annotation_sheet_csv = prefix + '_samples.csv'

        return


class FilePathExtractIlluminaLane(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathExtractIlluminaLane} models files in a directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar samples_directory: Samples directory
    @type samples_directory: str | unicode
    @ivar barcode_tsv: Barcode TSV file for Picard I{ExtractIlluminaBarcodes}
    @type barcode_tsv: str | unicode
    @ivar metrics_tsv: Metrics TSV file
    @type metrics_tsv: str | unicode
    @ivar library_tsv: Library TSV file for Picard I{IlluminaBasecallsToSam}
    @type library_tsv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathExtractIlluminaLane} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathExtractIlluminaLane, self).__init__(prefix=prefix)

        self.output_directory = prefix + '_output'
        self.samples_directory = prefix + '_samples'
        self.barcode_tsv = prefix + '_barcode.tsv'
        self.metrics_tsv = prefix + '_metrics.tsv'
        self.library_tsv = prefix + '_library.tsv'

        return


class ExtractIlluminaRunFolder(PicardIlluminaRunFolder):
    """The C{bsf.analyses.picard.ExtractIlluminaRunFolder} class extracts data from an Illumina Run Folder.

    The analysis is based on Picard I{ExtractIlluminaBarcodes} and Picard I{IlluminaBasecallsToSam}.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_lane: C{bsf.Stage.name} for the lane-specific stage
    @type stage_name_lane: str
    @cvar stage_name_cell: C{bsf.Stage.name} for the flow cell-specific stage
    @type stage_name_cell: str
    @ivar samples_directory: BSF samples directory
    @type samples_directory: str | unicode | None
    @ivar library_path: Library annotation file path
    @type library_path: str | unicode | None
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
    @ivar vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
    @type vendor_quality_filter: dict[str, bool]
    """

    name = 'Picard Extract Illumina Run Folder Analysis'
    prefix = 'extract_illumina_run_folder'

    stage_name_lane = '_'.join((prefix, 'lane'))
    stage_name_cell = '_'.join((prefix, 'cell'))

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
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode | None
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | unicode | None
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | unicode | None
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
        @param samples_directory: BSF samples directory
        @type samples_directory: str | unicode | None
        @param library_path: Library annotation file path
        @type library_path: str | unicode | None
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
        @param vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
        @type vendor_quality_filter: dict[str, bool] | None
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

        if vendor_quality_filter is None:
            self.vendor_quality_filter = dict()
        else:
            self.vendor_quality_filter = vendor_quality_filter

        return

    @property
    def get_experiment_directory(self):
        """Get the experiment directory.

        The experiment directory is a concatenation of the sequences directory and the project name.
        @return: Experiment directory
        @rtype: str | unicode | None
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

        # Get the library annotation file.

        option = 'library_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_path = configuration.config_parser.get(section=section, option=option)

        # Get the general samples directory.

        option = 'samples_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.samples_directory = configuration.config_parser.get(section=section, option=option)

        # Get directory access mode permission bits.

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        # Get file access mode permission bits.

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        # Get the maximum number of mismatches.

        option = 'max_mismatches'
        if configuration.config_parser.has_option(section=section, option=option):
            self.max_mismatches = configuration.config_parser.getint(section=section, option=option)

        # Get the minimum base quality.

        option = 'min_base_quality'
        if configuration.config_parser.has_option(section=section, option=option):
            self.min_base_quality = configuration.config_parser.getint(section=section, option=option)

        # Get sequencing centre information.

        option = 'sequencing_centre'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequencing_centre = configuration.config_parser.get(section=section, option=option)

        option = 'lanes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.lanes = configuration.config_parser.getint(section=section, option=option)

        # Read the VendorQualityFilter section, which consists of flow cell chemistry type keys and boolean values
        # to set filtering.

        vqf_section = '.'.join((section, 'VendorQualityFilter'))

        if configuration.config_parser.has_section(section=vqf_section):
            for option_name in configuration.config_parser.options(section=vqf_section):
                self.vendor_quality_filter[option_name] = configuration.config_parser.getboolean(
                    section=vqf_section,
                    option=option_name)

        return

    def run(self):
        """Run a C{bsf.analyses.picard.ExtractIlluminaRunFolder} C{bsf.Analysis}.

        @return:
        @rtype:
        """

        def run_get_sample_file_name(sample_name):
            """Private function to format sample-specific BAM file names (i.e. project_lane#sample.bam).

            @param sample_name:
            @type sample_name: str | unicode
            @return:
            @rtype: str | unicode
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # Start of the run() method body.

        super(ExtractIlluminaRunFolder, self).run()

        self.samples_directory = self.configuration.get_absolute_path(
            file_path=self.samples_directory,
            default_path=bsf.standards.Default.absolute_samples())

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(
                'The ExtractIlluminaRunFolder samples_directory {!r} does not exist.'.format(self.samples_directory))

        experiment_directory = self.get_experiment_directory

        # Get sequencing centre information.

        if not self.sequencing_centre:
            self.sequencing_centre = bsf.standards.Operator.get_sequencing_centre()
            if not self.sequencing_centre:
                raise Exception("An 'ExtractIlluminaRunFolder' analysis requires a "
                                "'sequencing_centre' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self._irf.run_parameters.get_flow_cell_type not in self.vendor_quality_filter:
            raise Exception('Flow cell chemistry type {!r} not defined.'.
                            format(self._irf.run_parameters.get_flow_cell_type))

        # Get the library annotation sheet.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(path=self.library_path):
            raise Exception('Library annotation file {!r} does not exist.'.format(self.library_path))

        stage_lane = self.get_stage(name=self.stage_name_lane)
        stage_cell = self.get_stage(name=self.stage_name_cell)

        # Load the library annotation sheet file and validate.

        library_annotation_sheet = bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)
        """ @type library_annotation_sheet: bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet """

        validation_messages = library_annotation_sheet.validate(lanes=self.lanes)

        if validation_messages:
            if self.force:
                warnings.warn('Validation of library annotation sheet {!r}:\n{}'.
                              format(self.library_path, validation_messages))
            else:
                raise Exception('Validation of library annotation sheet {!r}:\n{}'.
                                format(self.library_path, validation_messages))

        flow_cell_dict = library_annotation_sheet.index_by_lane()

        file_path_cell = FilePathExtractIlluminaCell(prefix=self.project_name)

        # Create a Sample Annotation Sheet in the project directory and
        # eventually transfer it into the experiment_directory.
        sample_annotation_sheet = bsf.ngs.SampleAnnotationSheet(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))

        # For each lane in the flow_cell_dict ...
        # TODO: For the moment this depends on the lanes (keys) defined in the LibraryAnnotationSheet.
        # Not all lanes may thus get extracted.
        # TODO: For NextSeq instruments, it would be sufficient to require annotation for only lane one and
        # copy information to lanes two to four internally.

        cell_dependency_list = list()
        """ @type cell_dependency_list: list[str] """

        for lane_str in sorted(flow_cell_dict):
            file_path_lane = FilePathExtractIlluminaLane(
                prefix=self.get_prefix_lane(
                    project_name=self.project_name,
                    lane=lane_str))

            # BARCODE_FILE
            eib_sheet = ExtractIlluminaBarcodesSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.barcode_tsv))

            # LIBRARY_PARAMS
            ibs_sheet = IlluminaBasecallsToSamSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.library_tsv))

            # Initialise a list of barcode sequence lengths required for the read structure calculation below.
            bc_length_list = list()
            """ @type bc_length_list: list[int] """

            # Sort each lane by sample name.
            flow_cell_dict_list = flow_cell_dict[lane_str]
            flow_cell_dict_list.sort(cmp=lambda x, y: cmp(x['sample_name'], y['sample_name']))

            for row_dict in flow_cell_dict_list:
                # Determine and check the length of the barcode sequences.
                for index in range(0, 1 + 1):
                    if len(bc_length_list) == index:
                        # If this is the first barcode, assign it.
                        bc_length_list.append(len(row_dict['barcode_sequence_' + str(index + 1)]))
                    else:
                        # If this a subsequent barcode, check it.
                        bc_length = len(row_dict['barcode_sequence_' + str(index + 1)])
                        if bc_length != bc_length_list[index]:
                            # Barcode lengths do not match ...
                            warnings.warn(
                                'The length ({}) of barcode {} {!r} does not match ' +
                                'the length ({}) of previous barcodes.'.format(
                                    bc_length,
                                    index + 1,
                                    row_dict['barcode_sequence_' + str(index + 1)],
                                    bc_length_list[index]),
                                UserWarning)

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
                    'File Type': 'Automatic',
                    'ProcessedRunFolder Name': self.project_name,
                    'Project Name': row_dict['library_name'],
                    'Project Size': row_dict['library_size'],
                    'Sample Name': row_dict['sample_name'],
                    'PairedReads Index 1': row_dict['barcode_sequence_1'],
                    'PairedReads Index 2': row_dict['barcode_sequence_2'],
                    # TODO: It would be good to add a RunnableStep to populate the ReadGroup.
                    'PairedReads ReadGroup': '',
                    'Reads1 Name': '_'.join((self.project_name, lane_str, row_dict['sample_name'])),
                    'Reads1 File': os.path.join(
                        os.path.basename(experiment_directory),
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    'Reads2 Name': '',
                    'Reads2 File': '',
                })

            # Calculate the read structure string from the IRF and the bc_length_list above ...

            read_structure = str()
            index_read_index = 0  # Number of index reads.
            # Instantiate and sort a new list of RunInformationRead objects.
            run_information_read_list = list(self._irf.run_information.run_information_read_list)
            run_information_read_list.sort(cmp=lambda x, y: cmp(x.number, y.number))
            for run_information_read in run_information_read_list:
                if run_information_read.index:
                    # For an index read ...
                    read_structure += str(bc_length_list[index_read_index]) + 'B'
                    if run_information_read.cycles < bc_length_list[index_read_index]:
                        read_structure += str(run_information_read.cycles - bc_length_list[index_read_index]) + 'S'
                    index_read_index += 1  # Increment to the next barcode read
                else:
                    # For a template read ...
                    read_structure += str(run_information_read.cycles) + 'T'

            # Adjust the IlluminaBaseCallsToSamSheet by adding an entry for unassigned reads and constrain
            # the columns to the number of index reads.
            ibs_sheet.adjust(
                unassigned_file_path=os.path.join(
                    file_path_lane.samples_directory,
                    run_get_sample_file_name(sample_name='0')),
                index_read_number=index_read_index)

            # Write the lane-specific Picard ExtractIlluminaBarcodesSheet and Picard IlluminaBasecallsToSamSheet.

            if index_read_index > 0:
                eib_sheet.to_file_path()

            ibs_sheet.to_file_path()

            # Create a Runnable and Executable for the lane stage.

            runnable_lane = self.add_runnable(
                runnable=bsf.Runnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory,
                    file_path_object=file_path_lane))
            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)

            cell_dependency_list.append(executable_lane.name)

            # Create an output_directory in the project_directory.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_output_directory',
                    directory_path=file_path_lane.output_directory))

            # Create a samples_directory in the project_directory.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_samples_directory',
                    directory_path=file_path_lane.samples_directory))

            # Create a RunnableStep for Picard ExtractIlluminaBarcodes, only if index (barcode) reads are present.

            if index_read_index > 0:
                runnable_step = runnable_lane.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_extract_illumina_barcodes',
                        java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx2G',
                        picard_classpath=self.classpath_picard,
                        picard_command='ExtractIlluminaBarcodes'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
                runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.output_directory)
                runnable_step.add_picard_option(key='LANE', value=lane_str)
                runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)
                runnable_step.add_picard_option(key='BARCODE_FILE', value=file_path_lane.barcode_tsv)
                runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
                if self.max_mismatches is not None:
                    # Maximum mismatches for a barcode to be considered a match. Default value: '1'.
                    runnable_step.add_picard_option(key='MAX_MISMATCHES', value=str(self.max_mismatches))
                # MIN_MISMATCH_DELTA: Minimum difference between number of mismatches in the best and
                # second best barcodes for a barcode to be considered a match. Default value: '1'.
                # MAX_NO_CALLS Maximum allowable number of no-calls in a barcode read before it is
                # considered unmatchable. Default value: '2'.
                if self.min_base_quality is not None:
                    # Minimum base quality. Any barcode bases falling below this quality will be considered
                    # a mismatch even in the bases match. Default value: '0'.
                    runnable_step.add_picard_option(key='MINIMUM_BASE_QUALITY', value=str(self.min_base_quality))
                # MINIMUM_QUALITY The minimum quality (after transforming 0s to 1s) expected from reads.
                # If qualities are lower than this value, an error is thrown.The default of 2 is what the
                # Illumina's spec describes as the minimum, but in practice the value has been observed lower.
                # Default value: '2'.
                runnable_step.add_picard_option(key='COMPRESS_OUTPUTS', value='true')
                runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_lane.get_relative_temporary_directory_path)
                # VERBOSITY defaults to 'INFO'.
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                # MAX_RECORDS_IN_RAM  defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                # REFERENCE_SEQUENCE
                # GA4GH_CLIENT_SECRETS defaults to 'client_secrets.json'.
                # USE_JDK_DEFLATER
                # USE_JDK_INFLATER
                # OPTIONS_FILE

            # Picard IlluminaBasecallsToSam

            # Create a RunnableStep for Picard IlluminaBasecallsToSam.

            runnable_step = runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepPicard(
                    name='picard_illumina_basecalls_to_sam',
                    java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    picard_classpath=self.classpath_picard,
                    picard_command='IlluminaBasecallsToSam'))
            """ @type runnable_step: bsf.process.RunnableStepPicard """
            runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
            if index_read_index > 0:
                runnable_step.add_picard_option(key='BARCODES_DIR', value=file_path_lane.output_directory)
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # OUTPUT is deprecated.
            runnable_step.add_picard_option(
                key='RUN_BARCODE',
                value=self._irf.run_parameters.get_flow_cell_barcode)
            # SAMPLE_ALIAS is deprecated.
            runnable_step.add_picard_option(
                key='READ_GROUP_ID',
                value='_'.join((self._irf.run_parameters.get_flow_cell_barcode, lane_str)))
            # LIBRARY_NAME is deprecated.
            runnable_step.add_picard_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)
            # NOTE: The ISO date format still does not work for Picard tools 2.6.1. Sigh.
            # runnable_step.add_picard_option(key='RUN_START_DATE', value=self._irf.run_information.get_iso_date)
            # NOTE: The only date format that seems to work is MM/DD/YYYY. Why?
            runnable_step.add_picard_option(
                key='RUN_START_DATE',
                value='/'.join((
                    self._irf.run_information.date[2:4],
                    self._irf.run_information.date[4:6],
                    '20' + self._irf.run_information.date[0:2])))
            # PLATFORM The name of the sequencing technology that produced the read. Default value: illumina.
            # NOTE: IlluminaToBam defaults to 'ILLUMINA'.
            # runnable_step.add_picard_option(key='PLATFORM', value='ILLUMINA')
            runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)
            # BARCODE_PARAMS is deprecated.
            runnable_step.add_picard_option(key='LIBRARY_PARAMS', value=file_path_lane.library_tsv)
            # ADAPTERS_TO_CHECK defaults to [INDEXED, DUAL_INDEXED, NEXTERA_V2, FLUIDIGM].
            # runnable_step.add_picard_option(key='ADAPTERS_TO_CHECK', value='')
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # FIRST_TILE defaults to 'null'.
            # TILE_LIMIT defaults to 'null'.
            # FORCE_GC defaults to 'true'.
            # APPLY_EAMSS_FILTER
            # MAX_READS_IN_RAM_PER_TILE
            # MINIMUM_QUALITY defaults to '2'.
            # INCLUDE_NON_PF_READS
            if not self.vendor_quality_filter[self._irf.run_parameters.get_flow_cell_type]:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
            # IGNORE_UNEXPECTED_BARCODES defaults to 'false'.
            # MOLECULAR_INDEX_TAG defaults to 'RX'.
            # MOLECULAR_INDEX_BASE_QUALITY_TAG defaults to 'QX'.
            # TAG_PER_MOLECULAR_INDEX The list of tags to store each molecular index.
            runnable_step.add_picard_option(key='TMP_DIR', value=runnable_lane.get_relative_temporary_directory_path)
            # VERBOSITY defaults to 'INFO'.
            # QUIET defaults to 'false'.
            # VALIDATION_STRINGENCY defaults to 'STRICT'.
            # COMPRESSION_LEVEL defaults to '5'.
            runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
            # MAX_RECORDS_IN_RAM  defaults to '500000'.
            # CREATE_INDEX defaults to 'false'.
            # CREATE_MD5_FILE defaults to 'false'.
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # REFERENCE_SEQUENCE
            # GA4GH_CLIENT_SECRETS defaults to 'client_secrets.json'.
            # USE_JDK_DEFLATER
            # USE_JDK_INFLATER
            # OPTIONS_FILE

            # Create the experiment directory if it does not exist already.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=experiment_directory))

            # Move the samples directory into the experiment directory.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_samples_directory',
                    source_path=file_path_lane.samples_directory,
                    target_path=experiment_directory))

            # Move the metrics file into the experiment directory.

            if index_read_index > 0:
                runnable_lane.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_metrics_tsv',
                        source_path=file_path_lane.metrics_tsv,
                        target_path=experiment_directory))

        # Finally, write the flow cell-specific SampleAnnotationSheet to the internal file path.

        sample_annotation_sheet.to_file_path()

        # Create a flow-cell specific Runnable.

        runnable_cell = self.add_runnable(
            runnable=bsf.Runnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_object=file_path_cell))
        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        # Move the Sample Annotation Sheet from the project_directory to the experiment_directory.

        if os.path.exists(sample_annotation_sheet.file_path):
            runnable_cell.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_sample_annotation',
                    source_path=file_path_cell.sample_annotation_sheet_csv,
                    target_path=experiment_directory))

        # Change directory and file access permissions.

        runnable_cell.add_runnable_step(
            runnable_step=bsf.process.RunnableStepChangeMode(
                name='chmod',
                file_path=experiment_directory,
                mode_directory=self.mode_directory,
                mode_file=self.mode_file))

        return


class FilePathIlluminaMultiplexSam(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathIlluminaMultiplexSam} class models files in a directory.

    Attributes:
    @ivar unsorted_bam: Unsorted BAM file
    @type unsorted_bam: str | unicode
    @ivar unsorted_md5: Unsorted BAM file MD5 check sum
    @type unsorted_md5: str | unicode
    @ivar sorted_bam: Sorted BAM file
    @type sorted_bam: str | unicode
    @ivar sorted_md5: Sorted BAM file MD5 check sum
    @type sorted_md5: str | unicode
    @ivar archive_bam: Archive BAM file
    @type archive_bam: str | unicode
    @ivar archive_md5: Archive BAM file MD5 check sum
    @type archive_md5: str | unicode
    """

    def __init__(self, project_name, lane, experiment_directory):
        """Initialise a C{bsf.analyses.picard.FilePathIlluminaMultiplexSam} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param experiment_directory: Experiment-specific directory
        @type experiment_directory: str | unicode
        @return:
        @rtype:
        """
        prefix = IlluminaMultiplexSam.get_prefix_lane(project_name=project_name, lane=lane)

        super(FilePathIlluminaMultiplexSam, self).__init__(prefix=prefix)

        self.unsorted_bam = prefix + '_unsorted.bam'
        self.unsorted_md5 = prefix + '_unsorted.bam.md5'
        self.sorted_bam = prefix + '_sorted.bam'
        self.sorted_md5 = prefix + '_sorted.bam.md5'
        # The final BAM and MD5 files are non-standard in that they do not contain a prefix and
        # reside in the experiment_directory.
        self.archive_bam = os.path.join(experiment_directory, '_'.join((project_name, lane)) + '.bam')
        self.archive_md5 = os.path.join(experiment_directory, '_'.join((project_name, lane)) + '.bam.md5')

        return


class IlluminaMultiplexSam(PicardIlluminaRunFolder):
    """The C{bsf.analyses.picard.IlluminaMultiplexSam} class represents the
    Picard IlluminaMultiplexSam analysis.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_lane: C{bsf.Stage.name} for the lane-specific stage
    @type stage_name_lane: str
    @cvar stage_name_cell: C{bsf.Stage.name} for the flow cell-specific stage
    @type stage_name_cell: str
    @ivar sequencing_centre: Sequencing centre
    @type sequencing_centre: str | None
    @ivar sequences_directory: Sequences directory to store archive BAM files
    @type sequences_directory: str | unicode | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar eamss_filter: Apply the Illumina EAMSS or Read Segment Quality Control Metric filter
    @type eamss_filter: bool | None
    @ivar vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
    @type vendor_quality_filter: dict[str, bool] | None
    """

    name = 'Picard IlluminaMultiplexSam Analysis'
    prefix = 'illumina_multiplex_sam'

    stage_name_lane = '_'.join((prefix, 'lane'))
    stage_name_cell = '_'.join((prefix, 'cell'))

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
            vendor_quality_filter=None):
        """Initialise a C{bsf.analyses.picard.IlluminaMultiplexSam} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode | None
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | unicode | None
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | unicode | None
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
        @param sequencing_centre: Sequencing centre
        @type sequencing_centre: str | None
        @param sequences_directory: Sequences directory to store archive BAM files
        @type sequences_directory: str | unicode | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param eamss_filter: Apply the Illumina EAMSS or Read Segment Quality Control Metric filter
        @type eamss_filter: bool | None
        @param vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
        @type vendor_quality_filter: dict[str, bool] | None
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

        if vendor_quality_filter is None:
            self.vendor_quality_filter = dict()
        else:
            self.vendor_quality_filter = vendor_quality_filter

        return

    @property
    def get_experiment_directory(self):
        """Get the experiment directory.

        The experiment directory is a concatenation of the sequences directory and the project name.
        @return: Experiment directory
        @rtype: str | unicode | None
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

        # Get sequencing centre information.

        option = 'sequencing_centre'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequencing_centre = configuration.config_parser.get(section=section, option=option)

        # Get the sequences directory information.

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        # Get directory access mode permission bits.

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        # Get file access mode permission bits.

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        # Get the EAMSS filter.

        option = 'eamss_filter'
        if configuration.config_parser.has_option(section=section, option=option):
            self.eamss_filter = configuration.config_parser.getboolean(section=section, option=option)

        # Read the VendorQualityFilter section, which consists of flow cell chemistry type keys and boolean values
        # to set filtering.

        vqf_section = '.'.join((section, 'VendorQualityFilter'))

        if configuration.config_parser.has_section(section=vqf_section):
            for option_name in configuration.config_parser.options(section=vqf_section):
                self.vendor_quality_filter[option_name] = configuration.config_parser.getboolean(
                    section=vqf_section,
                    option=option_name)

        return

    def run(self):
        """Run this C{bsf.analyses.picard.IlluminaMultiplexSam} C{bsf.Analysis}.

        @return:
        @rtype:
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
            default_path=bsf.standards.Default.absolute_illumina_run())

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The Illumina run directory {!r} does not exist.'.format(self.run_directory))

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise bsf.illumina.RunFolderNotComplete(
                'The Illumina run directory {!r} is not complete.'.format(self.run_directory))

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
            self.sequences_directory = bsf.standards.Default.absolute_sequences()

        # As a safety measure, to prevent creation of rogue directory paths, the sequences_directory has to exist.

        if not os.path.isdir(self.sequences_directory):
            raise Exception(
                'The IlluminaToBam sequences_directory {!r} does not exist.'.format(self.sequences_directory))

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the Picard tools Java Archive (JAR) class path directory

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception("An 'IlluminaToBam' analysis requires a "
                                "'classpath_picard' configuration option.")

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if irf.run_parameters.get_flow_cell_type not in self.vendor_quality_filter:
            raise Exception('Flow cell chemistry type {!r} not defined.'.format(irf.run_parameters.get_flow_cell_type))

        # Call the run method of the super class after the project_name has been defined.

        super(IlluminaMultiplexSam, self).run()

        cell_dependency_list = list()

        stage_lane = self.get_stage(name=self.stage_name_lane)
        stage_cell = self.get_stage(name=self.stage_name_cell)

        for lane_int in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            file_path_lane = FilePathIlluminaMultiplexSam(
                project_name=self.project_name,
                lane=lane_str,
                experiment_directory=experiment_directory)

            runnable_lane = self.add_runnable(
                runnable=bsf.Runnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory,
                    file_path_object=file_path_lane))

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

            runnable_step = runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepPicard(
                    name='picard_illumina_basecalls_to_undemux_sam',
                    java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx46G',
                    picard_classpath=self.classpath_picard,
                    picard_command='IlluminaBasecallsToUndemuxSam'))
            """ @type runnable_step: bsf.process.RunnableStepPicard """
            # RUN_DIR is required.
            runnable_step.add_picard_option(key='RUN_DIR', value=self.run_directory)
            # LANE is required.
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # OUTPUT is required.
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_lane.sorted_bam)
            # SEQUENCING_CENTER defaults to 'SC' for Sanger Center.
            runnable_step.add_picard_option(key='SEQUENCING_CENTER', value=self.sequencing_centre)
            # READ_STRUCTURE
            # ADAPTERS_TO_CHECK
            # FIVE_PRIME_ADAPTER
            # THREE_PRIME_ADAPTER
            # NUM_PROCESSORS
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # FIRST_TILE
            # TILE_LIMIT
            # FORCE_GC
            # APPLY_EAMSS_FILTER
            if self.eamss_filter:
                runnable_step.add_picard_option(key='APPLY_EAMSS_FILTER', value='true')
            else:
                runnable_step.add_picard_option(key='APPLY_EAMSS_FILTER', value='false')
            # MAX_READS_IN_RAM_PER_TILE
            # MINIMUM_QUALITY
            # INCLUDE_NON_PF_READS
            if self.vendor_quality_filter[irf.run_parameters.get_flow_cell_type]:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
            else:
                runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
            # MOLECULAR_INDEX_TAG
            # MOLECULAR_INDEX_BASE_QUALITY_TAG
            # TAG_PER_MOLECULAR_INDEX
            # TMP_DIR
            runnable_step.add_picard_option(key='TMP_DIR', value=runnable_lane.get_relative_temporary_directory_path)
            # VERBOSITY defaults to 'INFO'.
            # QUIET defaults to 'false'.
            # VALIDATION_STRINGENCY defaults to 'STRICT'.
            # COMPRESSION_LEVEL defaults to '5'.
            runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
            # MAX_RECORDS_IN_RAM  defaults to '500000'.
            # CREATE_INDEX defaults to 'false'.
            # CREATE_MD5_FILE defaults to 'false'.
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # REFERENCE_SEQUENCE
            # GA4GH_CLIENT_SECRETS defaults to 'client_secrets.json'.
            # USE_JDK_DEFLATER
            # USE_JDK_INFLATER
            # OPTIONS_FILE

            # Create the experiment directory if it does not exist already.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=experiment_directory))

            # Move and rename the final, sorted BAM file.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_sorted_bam',
                    source_path=file_path_lane.sorted_bam,
                    target_path=file_path_lane.archive_bam))

            # Move and rename the checksum file.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_sorted_md5',
                    source_path=file_path_lane.sorted_md5,
                    target_path=file_path_lane.archive_md5))

        # Add another flow cell-specific Runnable to reset directory and file mode permissions if requested.

        runnable_cell = self.add_runnable(
            runnable=bsf.Runnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))

        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        runnable_cell.add_runnable_step(
            runnable_step=bsf.process.RunnableStepChangeMode(
                name='chmod',
                file_path=experiment_directory,
                mode_directory=self.mode_directory,
                mode_file=self.mode_file))

        return


class FilePathIlluminaDemultiplexSamCell(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell} models flow cell-specific files.

    Attributes:
    @ivar sample_annotation_sheet_csv: Sample Annotation Sheet CSV file
    @type sample_annotation_sheet_csv: str | unicode
    """

    def __init__(self, project_name):
        """Initialise a C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamCell} object.

        @param project_name: Project name
        @type project_name: str
        @return:
        @rtype:
        """
        prefix = IlluminaDemultiplexSam.get_prefix_cell(project_name=project_name)

        super(FilePathIlluminaDemultiplexSamCell, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the prefix accordingly.
        prefix = project_name
        self.sample_annotation_sheet_csv = prefix + '_samples.csv'

        return


class FilePathIlluminaDemultiplexSamLane(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane} models lane-specific files.

    Attributes:
    @ivar project_barcode: Project-specific barcode CSV file
    @type project_barcode: str | unicode
    @ivar library_tsv: Library TSV file for Picard I{IlluminaBamDemux}
    @type library_tsv: str | unicode
    @ivar metrics_tsv: Lane-specific metrics TSV file
    @type metrics_tsv: str | unicode
    @ivar samples_directory: Sample directory
    @type samples_directory: str | unicode
    @ivar archive_bam: Archive BAM file
    @type archive_bam: str | unicode
    @ivar archive_md5: Archive BAM MD5 check sum
    @type archive_md5: str | unicode
    """

    def __init__(self, project_name, lane, sequences_directory):
        """Initialise a C{bsf.analyses.picard.FilePathIlluminaDemultiplexSamLane} object.

        @param project_name: Project name
        @type project_name: str
        @param lane: Lane
        @type lane: str
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str | unicode
        @return:
        @rtype:
        """
        prefix = IlluminaDemultiplexSam.get_prefix_lane(project_name=project_name, lane=lane)

        super(FilePathIlluminaDemultiplexSamLane, self).__init__(prefix=prefix)

        # All paths are non-standard in that they do not contain the regular Analysis Stage prefix.
        # Re-define the prefix accordingly.
        prefix = '_'.join((project_name, lane))
        self.project_barcode = prefix + '_barcode.tsv'
        # self.barcode_tsv = prefix + '_barcode.tsv'
        self.library_tsv = prefix + '_library.tsv'
        self.metrics_tsv = prefix + '_metrics.tsv'
        self.samples_directory = prefix + '_samples'
        # The input (sequence archive) files are non-standard in that they do not contain a prefix and
        # reside in the sequences_directory.
        self.archive_bam = os.path.join(sequences_directory, prefix + '.bam')
        self.archive_md5 = os.path.join(sequences_directory, prefix + '.bam.md5')

        return


class IlluminaDemultiplexSamSheet(bsf.annotation.AnnotationSheet):
    """The C{bsf.analyses.picard.IlluminaDemultiplexSamSheet} class represents a
    Tab-Separated Value (TSV) table of library information for the
    C{bsf.analyses.picard.IlluminaDemultiplexSam} C{bsf.Analysis}.

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

    def adjust(self, index_read_number):
        """Adjust the C{bsf.analyses.picard.IlluminaBasecallsToSamSheet} by inserting the
        file path for unassigned reads and by adjusting the index columns to the number of
        index reads.

        #@param unassigned_file_path: File path for unassigned reads
        #@type unassigned_file_path: str | unicode
        @param index_read_number: Number of index reads
        @type index_read_number: int
        @return:
        @rtype:
        """
        # The IlluminaBasecallsToSamSheet needs adjusting ...
        # if len(self.row_dicts) == 1 and len(self.row_dicts[0]['BARCODE_1']) == 0 and len(
        #         self.row_dicts[0]['BARCODE_2']) == 0:
        #     # ... if a single sample, but neither BARCODE_1 nor BARCODE_2 were defined,
        #     # BARCODE_1 needs setting to 'N'.
        #     self.row_dicts[0]['BARCODE_1'] = 'N'
        # else:
        #     # ... in all other cases, since a last row for unmatched barcode sequences needs adding.
        #     self.row_dicts.append({
        #         'BARCODE_1': 'N',
        #         'BARCODE_2': '',
        #         'OUTPUT': unassigned_file_path,
        #         'SAMPLE_NAME': 'Unmatched',
        #         'LIBRARY_NAME': self.row_dicts[0]['LIBRARY_NAME'],
        #     })

        # Further adjust the IlluminaBaseCallsToSamSheet and remove any BARCODE_N columns not represented
        # in the read structure.
        for index in range(0, 1 + 1):
            if index + 1 > index_read_number:
                # Remove the 'BARCODE_N' field from the list of field names.
                if 'BARCODE_' + str(index + 1) in self.field_names:
                    self.field_names.remove('BARCODE_' + str(index + 1))
                # Remove the 'BARCODE_N' entry from each row_dict object, since csv.DictWriter requires it.
                for row_dict in self.row_dicts:
                    row_dict.pop('BARCODE_' + str(index + 1), None)

        return


class IlluminaDemultiplexSam(bsf.Analysis):
    """The C{bsf.analyses.picard.IlluminaDemultiplexSam} class represents the logic to
    decode sequence archive BAM files into sample-specific BAM files.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_lane: C{bsf.Stage.name} for the lane-specific stage
    @type stage_name_lane: str
    @cvar stage_name_cell: C{bsf.Stage.name} for the flow cell-specific stage
    @type stage_name_cell: str
    @ivar library_path: Library annotation file path
    @type library_path: str | unicode | None
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | unicode | None
    @ivar sequences_directory: BSF sequences directory
    @type sequences_directory: str | unicode | None
    @ivar samples_directory: BSF samples directory
    @type samples_directory: str | unicode | None
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode | None
    @ivar lanes: Number of lanes on the flow cell
    @type lanes: int | None
    @ivar force: Force de-multiplexing with a Library Annotation sheet failing validation
    @type force: bool | None
    """

    name = 'Picard IlluminaDemultiplexSam Analysis'
    prefix = 'illumina_demultiplex_sam'

    stage_name_lane = '_'.join((prefix, 'lane'))
    stage_name_cell = '_'.join((prefix, 'cell'))

    @classmethod
    def get_prefix_cell(cls, project_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects.

        @param project_name: A project name
        @type project_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_cell, project_name))

    @classmethod
    def get_prefix_lane(cls, project_name, lane):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_lane, project_name, lane))

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
            lanes=8,
            force=False):
        """Initialise a C{bsf.analyses.picard.IlluminaDemultiplexSam} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param library_path: Library annotation file path
        @type library_path: str | unicode | None
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode | None
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str | unicode | None
        @param samples_directory: BSF samples directory
        @type samples_directory: str | unicode | None
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
        @param lanes: Number of lanes on the flow cell
        @type lanes: int | None
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
        self.lanes = lanes
        self.force = force

        return

    @property
    def get_experiment_directory(self):
        """Get the experiment directory.

        The experiment directory is a concatenation of the samples directory and the project name.
        @return: Experiment directory
        @rtype: str | unicode | None
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

        # Get the library annotation file.

        option = 'library_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_path = configuration.config_parser.get(section=section, option=option)

        # Get Illumina Run Folder information.

        option = 'illumina_run_folder'
        if configuration.config_parser.has_option(section=section, option=option):
            self.run_directory = configuration.config_parser.get(section=section, option=option)

        # Get the general sequences directory.

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        # Get the general samples directory.

        option = 'samples_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.samples_directory = configuration.config_parser.get(section=section, option=option)

        # Get directory access mode permission bits.

        option = 'mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_directory = configuration.config_parser.get(section=section, option=option)

        # Get file access mode permission bits.

        option = 'mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mode_file = configuration.config_parser.get(section=section, option=option)

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'lanes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.lanes = configuration.config_parser.getint(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.picard.IlluminaDemultiplexSam} analysis to
        decode an archive BAM file produced with IlluminaMultiplexSam into sample-specific BAM files.
        @return:
        @rtype:
        """

        def run_get_sample_file_name(sample_name):
            """Private function to format sample-specific BAM file names (i.e. project_lane#sample.bam).

            @param sample_name:
            @type sample_name: str | unicode
            @return:
            @rtype: str | unicode
            """
            # FIXME: Replace with a @property method in the FilePathLane object?
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        # a Picard tools *tab-separated* value (TSV) sample sheet.
        # lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        # barcode_sequence, barcode_name, library_name, sample_name, description

        super(IlluminaDemultiplexSam, self).run()

        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception('An Illumina run directory or file path has not been defined.')

        # self.run_directory = self.configuration.get_absolute_path(
        #     file_path=self.run_directory,
        #     default_path=bsf.standards.Default.absolute_illumina_run())
        # FIXME: Lookup the Illumina Run Folder.
        # Illumina Run Folders are kept either in the 'active' or the 'archive' location.
        # Upon re-running of the de-multiplexing stage, the IRF may have been archived already.

        # Lookup the Illumina Run Folder in the 'active' or 'archive' locations.
        self.run_directory = bsf.illumina.RunFolder.absolute_file_path(name=self.run_directory)

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception('The Illumina run directory {!r} does not exist.'.
                            format(self.run_directory))

        irf = bsf.illumina.RunFolder.from_file_path(file_path=self.run_directory)

        # Define the sequences and samples directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.sequences_directory:
            self.sequences_directory = self.project_name

        self.sequences_directory = self.configuration.get_absolute_path(
            file_path=self.sequences_directory,
            default_path=bsf.standards.Default.absolute_sequences())

        self.samples_directory = self.configuration.get_absolute_path(
            file_path=self.samples_directory,
            default_path=bsf.standards.Default.absolute_samples())

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(
                'The ' + self.name + ' samples_directory ' + repr(self.samples_directory) + ' does not exist.')

        # Get the experiment_directory once.

        experiment_directory = self.get_experiment_directory

        # Get the library annotation sheet.
        # The library annotation sheet is deliberately not passed in via sas_file,
        # as the bsf.Analysis.run() method reads that option into a BSF Collection object.

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = self.configuration.get_absolute_path(file_path=self.library_path)

        if not os.path.exists(path=self.library_path):
            raise Exception('Library annotation file ' + repr(self.library_path) + ' does not exist.')

        # Load the library annotation sheet file and validate.

        library_annotation_sheet = bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet.from_file_path(
            file_path=self.library_path)
        """ @type library_annotation_sheet: bsf.analyses.illumina_to_bam_tools.LibraryAnnotationSheet """

        validation_messages = library_annotation_sheet.validate(lanes=self.lanes)

        if validation_messages:
            if self.force:
                warnings.warn('Validation of library annotation sheet ' + repr(self.library_path) + ':\n' +
                              validation_messages)
            else:
                raise Exception('Validation of library annotation sheet ' + self.library_path + ':\n' +
                                validation_messages)

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('The ' + self.name + " requires a 'classpath_picard' configuration option.")

        stage_lane = self.get_stage(name=self.stage_name_lane)
        stage_cell = self.get_stage(name=self.stage_name_cell)

        flow_cell_dict = library_annotation_sheet.index_by_lane()

        file_path_cell = FilePathIlluminaDemultiplexSamCell(project_name=self.project_name)

        # Create a Collection and ProcessedRunFolder object to write a SampleAnnotationSheet.

        prf = bsf.ngs.ProcessedRunFolder(name=bsf.ngs.ProcessedRunFolder.default_name)

        collection = bsf.ngs.Collection(
            name='Samples',
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv))
        collection.add_processed_run_folder(prf=prf)

        cell_dependency_list = list()
        """ @type cell_dependency_list: list[str] """

        for lane_str in sorted(flow_cell_dict):
            file_path_lane = FilePathIlluminaDemultiplexSamLane(
                project_name=self.project_name,
                lane=lane_str,
                sequences_directory=self.sequences_directory)

            # experiment_samples_directory = os.path.join(experiment_directory, file_path_lane.samples_directory)

            # LIBRARY_PARAMS
            ibs_sheet = IlluminaDemultiplexSamSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.library_tsv))

            # Initialise a list of barcode sequence lengths required for the read structure calculation below.
            bc_length_list = list()
            """ @type bc_length_list: list[int] """

            for row_dict in flow_cell_dict[lane_str]:
                # Determine and check the length of the barcode sequences.
                for index in range(0, 1 + 1):
                    if len(bc_length_list) == index:
                        # If this is the first barcode, assign it.
                        bc_length_list.append(len(row_dict['barcode_sequence_' + str(index + 1)]))
                    else:
                        # If this a subsequent barcode, check it.
                        bc_length = len(row_dict['barcode_sequence_' + str(index + 1)])
                        if bc_length != bc_length_list[index]:
                            # Barcode lengths do not match ...
                            warnings.warn(
                                'The length ({}) of barcode {} {!r} does not match ' +
                                'the length ({}) of previous barcodes.'.format(
                                    bc_length,
                                    index + 1,
                                    row_dict['barcode_sequence_' + str(index + 1)],
                                    bc_length_list[index]),
                                UserWarning)

                # Add a row to the lane-specific Picard IlluminaBasecallsToSamSheet.

                ibs_sheet.row_dicts.append({
                    'BARCODE_1': row_dict['barcode_sequence_1'],
                    'BARCODE_2': row_dict['barcode_sequence_2'],
                    'OUTPUT': os.path.join(
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])),
                    # FIXME: IlluminaBamDemux does no seem standard and requires a SAMPLE_NAME?
                    # 'SAMPLE_ALIAS': row_dict['sample_name'],
                    'SAMPLE_NAME': row_dict['sample_name'],
                    'LIBRARY_NAME': row_dict['library_name'],
                })

                # Insert the ReadGroup into the Collection object.

                reads_1 = bsf.ngs.Reads(
                    name='_'.join((self.project_name, lane_str, row_dict['sample_name'])),
                    file_path=os.path.join(
                        os.path.basename(experiment_directory),
                        file_path_lane.samples_directory,
                        run_get_sample_file_name(sample_name=row_dict['sample_name'])))

                paired_reads = bsf.ngs.PairedReads(
                    reads_1=reads_1,
                    index_1=row_dict['barcode_sequence_1'],
                    index_2=row_dict['barcode_sequence_2'])
                # TODO: Annotate the ReadGroup with Run information?

                sample = bsf.ngs.Sample(name=row_dict['sample_name'])
                sample.add_paired_reads(paired_reads=paired_reads)

                project = bsf.ngs.Project(name=row_dict['library_name'])
                project.add_annotation(key='Size', value=row_dict['library_size'])
                project.add_sample(sample=sample)

                prf.add_project(project=project)

            # Calculate the read structure string from the IRF and the bc_length_list above ...

            read_structure = str()
            index_read_index = 0  # Number of index reads.
            # Instantiate and sort a new list of RunInformationRead objects.
            # FIXME: How to make the Illumina Run Folder available?
            # The IRF could be read form the @RG PU field, but that is not available at the time the run is configured.
            # Alternatively, the IRF has to be passed in.
            run_information_read_list = list(irf.run_information.run_information_read_list)
            run_information_read_list.sort(cmp=lambda x, y: cmp(x.number, y.number))
            for run_information_read in run_information_read_list:
                if run_information_read.index:
                    # For an index read ...
                    read_structure += str(bc_length_list[index_read_index]) + 'B'
                    if run_information_read.cycles < bc_length_list[index_read_index]:
                        read_structure += str(run_information_read.cycles - bc_length_list[index_read_index]) + 'S'
                    index_read_index += 1  # Increment to the next barcode read
                else:
                    # For a template read ...
                    read_structure += str(run_information_read.cycles) + 'T'

            # Adjust the IlluminaBaseCallsToSamSheet by adding an entry for unassigned reads and constrain
            # the columns to the number of index reads.
            ibs_sheet.adjust(
                # unassigned_file_path=os.path.join(
                #     file_path_lane.samples_directory,
                #     run_get_sample_file_name(sample_name='0')),
                index_read_number=index_read_index)

            # Create a Runnable and Executable for the lane stage.

            runnable_lane = self.add_runnable(
                runnable=bsf.Runnable(
                    name=self.get_prefix_lane(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory,
                    file_path_object=file_path_lane))
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

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_project_lane_samples_directory',
                    directory_path=file_path_lane.samples_directory))

            runnable_step = runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepPicard(
                    name='picard_extract_illumina_barcodes',
                    java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    picard_classpath=self.classpath_picard,
                    picard_command='IlluminaBamDemux'))
            """ @type runnable_step: bsf.process.RunnableStepPicard """
            # INPUT
            runnable_step.add_picard_option(key='INPUT', value=file_path_lane.archive_bam)
            # OUTPUT_DIR
            runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.samples_directory)
            # OUTPUT_PREFIX defaults to the @RG id field of the multiplexed BAM file
            runnable_step.add_picard_option(key='OUTPUT_PREFIX', value='_'.join((self.project_name, lane_str)))
            # OUTPUT_FORMAT
            # BARCODE_TAG_NAME defaults to 'BC'
            # BARCODE_QUALITY_TAG_NAME defaults to 'QT'
            # LIBRARY_PARAMS
            runnable_step.add_picard_option(key='LIBRARY_PARAMS', value=file_path_lane.library_tsv)
            # METRICS_FILE
            runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
            # MAX_MISMATCHES defaults to '1'
            # MIN_MISMATCH_DELTA defaults to '1'
            # MAX_NO_CALLS defaults to '2'
            # MINIMUM_BASE_QUALITY defaults to '2'
            # READ_STRUCTURE
            runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)
            runnable_step.add_picard_option(key='TMP_DIR', value=runnable_lane.get_relative_temporary_directory_path)
            # VERBOSITY defaults to 'INFO'.
            # QUIET defaults to 'false'.
            # VALIDATION_STRINGENCY defaults to 'STRICT'.
            # COMPRESSION_LEVEL defaults to '5'.
            runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
            # MAX_RECORDS_IN_RAM  defaults to '500000'.
            # CREATE_INDEX defaults to 'false'.
            # CREATE_MD5_FILE defaults to 'false'.
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # REFERENCE_SEQUENCE
            # GA4GH_CLIENT_SECRETS defaults to 'client_secrets.json'.
            # USE_JDK_DEFLATER
            # USE_JDK_INFLATER
            # OPTIONS_FILE

            # Create the experiment directory if it does not exist already.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_experiment_directory',
                    directory_path=experiment_directory))

            # Move the samples directory into the experiment directory.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_project_lane_samples_directory',
                    source_path=file_path_lane.samples_directory,
                    target_path=experiment_directory))

            # Move the metrics file into the experiment directory.

            runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_metrics_tsv',
                    source_path=file_path_lane.metrics_tsv,
                    target_path=experiment_directory))

        # Finally, write the flow cell-specific Collection to the internal file path.

        collection.to_sas_path(
            file_path=os.path.join(
                self.project_directory,
                file_path_cell.sample_annotation_sheet_csv),
            name='Samples')

        # Create a flow-cell specific Runnable.

        runnable_cell = self.add_runnable(
            runnable=bsf.Runnable(
                name=self.get_prefix_cell(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_object=file_path_cell))
        executable_cell = self.set_stage_runnable(
            stage=stage_cell,
            runnable=runnable_cell)
        executable_cell.dependencies.extend(cell_dependency_list)

        # Move the Sample Annotation Sheet from the project_directory to the experiment_directory.

        if os.path.exists(collection.file_path):
            runnable_cell.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMove(
                    name='move_sample_annotation',
                    source_path=file_path_cell.sample_annotation_sheet_csv,
                    target_path=experiment_directory))

        # Change directory and file access permissions.

        runnable_cell.add_runnable_step(
            runnable_step=bsf.process.RunnableStepChangeMode(
                name='chmod',
                file_path=experiment_directory,
                mode_directory=self.mode_directory,
                mode_file=self.mode_file))

        return


class FilePathCollectHiSeqXPfFailMetricsLane(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane} models files in a directory.

    Attributes:
    @ivar summary_tsv: Summary metrics TSV file
    @type summary_tsv: str | unicode
    @ivar detailed_tsv: Detailed metrics TSV file
    @type detailed_tsv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathCollectHiSeqXPfFailMetricsLane} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCollectHiSeqXPfFailMetricsLane, self).__init__(prefix=prefix)

        self.summary_tsv = prefix + '.pffail_summary_metrics'
        self.detailed_tsv = prefix + '.pffail_detailed_metrics'

        return


class CollectHiSeqXPfFailMetrics(PicardIlluminaRunFolder):
    """The C{bsf.analyses.picard.CollectHiSeqXPfFailMetrics} class represents Picard CollectHiSeqXPfFailMetrics.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_lane: C{bsf.Stage.name} for the lane-specific stage
    @type stage_name_lane: str
    """

    name = 'Picard CollectHiSeqXPfFailMetrics Analysis'
    prefix = 'picard_hiseq_x_pf_fail'

    stage_name_lane = '_'.join((prefix, 'lane'))

    def run(self):
        """Run the C{bsf.analyses.picard.CollectHiSeqXPfFailMetrics} C{bsf.Analysis}.

        @return:
        @rtype:
        """
        super(CollectHiSeqXPfFailMetrics, self).run()

        # Picard CollectHiSeqXPfFailMetrics

        stage_lane = self.get_stage(name=self.stage_name_lane)

        cell_dependency_list = list()

        for lane_int in range(0 + 1, self._irf.run_information.flow_cell_layout.lane_count + 1):
            lane_str = str(lane_int)

            prefix_lane = self.get_prefix_lane(
                project_name=self.project_name,
                lane=lane_str)

            file_path_lane = FilePathCollectHiSeqXPfFailMetricsLane(prefix=prefix_lane)

            # NOTE: The bsf.Runnable.name has to match the Executable.name that gets submitted via the Stage.
            runnable_lane = self.add_runnable(
                runnable=bsf.Runnable(
                    name=self.get_prefix_lane(
                        project_name=self.project_name,
                        lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory,
                    file_path_object=file_path_lane))

            executable_lane = self.set_stage_runnable(
                stage=stage_lane,
                runnable=runnable_lane)

            # Add the dependency for the cell-specific process.

            cell_dependency_list.append(executable_lane.name)

            runnable_step = runnable_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepPicard(
                    name='collect_hiseq_x_fail_metrics',
                    java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    picard_classpath=self.classpath_picard,
                    picard_command='CollectHiSeqXPfFailMetrics'))
            """ @type runnable_step: bsf.process.RunnableStepPicard """
            # BASECALLS_DIR is required.
            runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
            # OUTPUT is required.
            runnable_step.add_picard_option(key='OUTPUT', value=prefix_lane)
            # LANE is required.
            runnable_step.add_picard_option(key='LANE', value=lane_str)
            # NUM_PROCESSORS defaults to '1'.
            runnable_step.add_picard_option(key='NUM_PROCESSORS', value=str(stage_lane.threads))
            # N_CYCLES defaults to '24'. Should match Illumina RTA software.

        return


class FilePathDownsampleSam(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathDownsampleSam} models files in a directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar downsampled_bam: Down-sampled BAM file
    @type downsampled_bam: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathDownsampleSam} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathDownsampleSam, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.downsampled_bam = prefix + '_downsampled.bam'

        return


class DownsampleSam(bsf.Analysis):
    """The C{bsf.analyses.picard.DownsampleSam} class represents the logic to run the Picard DownsampleSam analysis.

    Attributes:

    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode | None
    """

    name = 'Picard DownsampleSam Analysis'
    prefix = 'picard_downsample_sam'

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
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
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

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.picard.DownsampleSam} C{bsf.Analysis}.

        This method changes the C{bsf.ngs.Collection} object of this C{bsf.Analysis} to update with FASTQ file paths.
        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.Analysis.sample_list} instance variable.
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

        stage_picard_dss = self.get_stage(name='picard_downsample_sam')

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

                    prefix_read_group = '_'.join((stage_picard_dss.name, paired_reads_name))

                    file_path_read_group = FilePathDownsampleSam(prefix=prefix_read_group)

                    # Keep the original BAM file and modify the file_path in the Reads object.
                    bam_file_path = reads.file_path
                    reads.file_path = os.path.join(self.project_directory, file_path_read_group.downsampled_bam)
                    # (file_root, file_extension) = os.path.splitext(bam_file_path)
                    # reads.file_path = file_root + 'downsampled' + file_extension

                    # Create a Runnable for running the Picard DownsampleSam analysis.

                    runnable_picard_dss = self.add_runnable(
                        runnable=bsf.Runnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory,
                            file_path_object=file_path_read_group))

                    # Create an Executable for running the Picard SamToFastq Runnable.

                    self.set_stage_runnable(stage=stage_picard_dss, runnable=runnable_picard_dss)

                    # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                    # runnable_picard_dss.add_runnable_step(
                    #     runnable_step=bsf.process.RunnableStepMakeDirectory(
                    #         name='make_directory',
                    #         directory_path=file_path_read_group.output_directory))

                    # Create a new RunnableStep for the Picard DownsampleSam program.

                    runnable_step = runnable_picard_dss.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepPicard(
                            name='picard_downsample_sam',
                            java_temporary_path=runnable_picard_dss.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx2G',
                            picard_classpath=self.classpath_picard,
                            picard_command='DownsampleSam'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    runnable_step.add_picard_option(key='INPUT', value=bam_file_path)
                    runnable_step.add_picard_option(
                        key='OUTPUT',
                        value=file_path_read_group.downsampled_bam)
                    # FIXME: Add to the configuration file and documentation.
                    if 'DownsampleSam Probability' in paired_reads.annotation_dict:
                        runnable_step.add_picard_option(
                            key='PROBABILITY',
                            value=paired_reads.annotation_dict['DownsampleSam Probability'][0])
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_picard_dss.get_relative_temporary_directory_path)
                    # VERBOSITY defaults to 'INFO'.
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    # QUIET defaults to 'false'.
                    runnable_step.add_picard_option(key='QUIET', value='false')
                    # VALIDATION_STRINGENCY defaults to 'STRICT'.
                    runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')
                    # COMPRESSION_LEVEL defaults to '5'.
                    # MAX_RECORDS_IN_RAM defaults to '500000'.
                    # CREATE_INDEX defaults to 'false'.
                    # CREATE_MD5_FILE defaults to 'false'.
                    # USE_JDK_DEFLATER
                    # USE_JDK_INFLATER
                    # OPTIONS_FILE

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'picard_downsample_sam_samples.csv'))),
            name='_'.join((self.project_name, 'picard_downsample_sam')))

        annotation_sheet.to_file_path()

        return


class FilePathSamToFastqReadGroup(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathSamToFastqReadGroup} class models read group-specific Picard SamToFastq files.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.picard.FilePathSamToFastqReadGroup} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathSamToFastqReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class FilePathSamToFastqProject(bsf.FilePath):
    """The C{bsf.analyses.picard.FilePathSamToFastqProject} class models project-specific Picard SamToFastq files.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar sas_path_old: Old Sample Annotation Sheet file path
    @type sas_path_old: str | unicode
    @ivar sas_path_new: New Sample Annotation Sheet file path
    @type sas_path_new: str | unicode
    """

    def __init__(self, prefix, prefix_analysis, project_name):
        """Initialise a C{bsf.analyses.picard.FilePathSamToFastqProject} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @param prefix_analysis: C{bsf.Analysis.prefix}
        @type prefix_analysis: str
        @param project_name: Project name
        @type project_name: str
        @return:
        @rtype
        """
        super(FilePathSamToFastqProject, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sas_path_old = '_'.join((project_name, prefix_analysis, 'original.csv'))
        self.sas_path_new = '_'.join((project_name, prefix_analysis, 'samples.csv'))

        return


class SamToFastq(bsf.Analysis):
    """The C{bsf.analyses.picard.SamToFastq} class represents the logic to run the Picard SamToFastq analysis.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_read_group: C{bsf.Stage.name} for read group-specific C{bsf.Runnable} objects
    @type stage_name_read_group: str
    @cvar stage_name_project: C{bsf.Stage.name} for project-specific C{bsf.Runnable} objects
    @type stage_name_project: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: None | str | unicode
    @ivar include_non_pf_reads: Include non-pass filer reads
    @type include_non_pf_reads: bool | None
    """

    name = 'Picard SamToFastq Analysis'
    prefix = 'picard_sam_to_fastq'

    stage_name_read_group = '_'.join((prefix, 'read_group'))
    stage_name_project = '_'.join((prefix, 'project'))

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
            include_non_pf_reads=None):
        """Initialise a C{bsf.analyses.picard.SamToFastq} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
        @param include_non_pf_reads: Include non-pass filer reads
        @type include_non_pf_reads: bool | None
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

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'include_non_pass_filter_reads'
        if configuration.config_parser.has_option(section=section, option=option):
            self.include_non_pass_filter_reads = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.picard.SamToFastq} C{bsf.Analysis}.

        This method changes the C{bsf.ngs.Collection} object of this C{bsf.Analysis} to update with FASTQ file paths.
        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        super(SamToFastq, self).run()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception("A 'SamToFastq' analysis requires a "
                                "'classpath_picard' configuration option.")

        run_read_comparisons()

        # Picard SamToFastq

        stage_read_group = self.get_stage(name=self.stage_name_read_group)
        stage_project = self.get_stage(name=self.stage_name_project)

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
                        prefix_read_group = '_'.join((stage_read_group.name, paired_reads_name))

                        file_path_read_group = FilePathSamToFastqReadGroup(prefix=prefix_read_group)

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
                            read_group_list.extend(map(lambda x: ':'.join((x, read_group_dict[x])),
                                                       read_group_dict.iterkeys()))
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

                        runnable_read_group = self.add_runnable(
                            runnable=bsf.Runnable(
                                name=prefix_read_group,
                                code_module='bsf.runnables.generic',
                                working_directory=self.project_directory,
                                file_path_object=file_path_read_group))
                        self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                        # Record the Executable.name for the project dependency.

                        project_dependency_list.append(runnable_read_group.name)

                        # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                        runnable_read_group.add_runnable_step(
                            runnable_step=bsf.process.RunnableStepMakeDirectory(
                                name='mkdir',
                                directory_path=file_path_read_group.output_directory))

                        # Create a new RunnableStep for the Picard SamToFastq program.

                        runnable_step = runnable_read_group.add_runnable_step(
                            runnable_step=bsf.process.RunnableStepPicard(
                                name='picard_sam_to_fastq',
                                java_temporary_path=runnable_read_group.get_relative_temporary_directory_path,
                                java_heap_maximum='Xmx2G',
                                picard_classpath=self.classpath_picard,
                                picard_command='SamToFastq'))
                        """ @type runnable_step: bsf.process.RunnableStepPicard """
                        runnable_step.add_picard_option(key='INPUT', value=bam_file_path)
                        # FASTQ
                        # SECOND_END_FASTQ
                        # UNPAIRED_FASTQ
                        runnable_step.add_picard_option(key='OUTPUT_PER_RG', value='true')
                        runnable_step.add_picard_option(key='COMPRESS_OUTPUTS_PER_RG', value='true')
                        # RG_TAG
                        runnable_step.add_picard_option(
                            key='OUTPUT_DIR',
                            value=file_path_read_group.output_directory)
                        # RE_REVERSE
                        # INTERLEAVE
                        if self.include_non_pass_filter_reads:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
                        else:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
                        # CLIPPING_ATTRIBUTE
                        # CLIPPING_ACTION
                        # CLIPPING_MIN_LENGTH
                        # READ1_TRIM
                        # READ1_MAX_BASES_TO_WRITE
                        # READ2_TRIM
                        # READ2_MAX_BASES_TO_WRITE
                        # QUALITY
                        # INCLUDE_NON_PRIMARY_ALIGNMENTS
                        runnable_step.add_picard_option(
                            key='TMP_DIR',
                            value=runnable_read_group.get_relative_temporary_directory_path)
                        # VERBOSITY defaults to 'INFO'.
                        runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                        # QUIET defaults to 'false'.
                        runnable_step.add_picard_option(key='QUIET', value='false')
                        # VALIDATION_STRINGENCY defaults to 'STRICT'.
                        runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')
                        # COMPRESSION_LEVEL defaults to '5'.
                        runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                        # MAX_RECORDS_IN_RAM defaults to '500000'.
                        # CREATE_INDEX defaults to 'false'.
                        # CREATE_MD5_FILE defaults to 'false'.
                        runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                        # REFERENCE_SEQUENCE
                        # GA4GH_CLIENT_SECRETS defaults to 'client_secrets.json'.
                        # USE_JDK_DEFLATER
                        # USE_JDK_INFLATER
                        # OPTIONS_FILE

        # Create a Runnable for pruning the sample annotation sheet.

        prefix_project = '_'.join((stage_project.name, self.project_name))

        file_path_project = FilePathSamToFastqProject(
            prefix=prefix_project,
            prefix_analysis=self.prefix,
            project_name=self.project_name)

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(self.project_directory, file_path_project.sas_path_old),
            name=prefix_project)
        annotation_sheet.to_file_path()

        runnable_project = self.add_runnable(
            runnable=bsf.Runnable(
                name=prefix_project,
                code_module='bsf.runnables.picard_sam_to_fastq_sample_sheet',
                working_directory=self.project_directory,
                file_path_object=file_path_project))
        executable_project = self.set_stage_runnable(
            stage=stage_project,
            runnable=runnable_project)
        executable_project.dependencies.extend(project_dependency_list)

        # Create a new RunnableStep.

        runnable_step_project = runnable_project.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='prune_sample_annotation_sheet'))

        runnable_step_project.add_option_long(key='sas_path_old', value=file_path_project.sas_path_old)
        runnable_step_project.add_option_long(key='sas_path_new', value=file_path_project.sas_path_new)
        runnable_step_project.add_option_long(key='minimum_size', value='1024')

        return
