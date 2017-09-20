"""bsf.analyses.picard

A package of classes and methods modelling Picard analyses data files and data directories.
"""

#
# Copyright 2013 - 2016 Michael K. Schuster
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


import os
import re
import warnings
import weakref

from bsf import Analysis, FilePath, Runnable
from bsf.analyses.illumina_to_bam_tools import LibraryAnnotationSheet
from bsf.annotation import AnnotationSheet
from bsf.ngs import Reads, PairedReads, SampleAnnotationSheet
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.process import RunnableStep, RunnableStepChangeMode, RunnableStepPicard, RunnableStepMakeDirectory,\
    RunnableStepMove
from bsf.standards import Default

import pysam


class PicardIlluminaRunFolder(Analysis):
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
    @type run_directory: str | unicode
    @ivar intensity_directory: File path to the I{Intensities} directory,
        defaults to I{illumina_run_folder/Data/Intensities}
    @type intensity_directory: str | unicode
    @ivar basecalls_directory: File path to the I{BaseCalls} directory,
        defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
    @type basecalls_directory: str | unicode
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool
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
        @type run_directory: str | unicode
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | unicode
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | unicode
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
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

        if run_directory is None:
            self.run_directory = str()
        else:
            self.run_directory = run_directory

        if intensity_directory is None:
            self.intensity_directory = str()
        else:
            self.intensity_directory = intensity_directory

        if basecalls_directory is None:
            self.basecalls_directory = str()
        else:
            self.basecalls_directory = basecalls_directory

        if experiment_name is None:
            self.experiment_name = str()
        else:
            self.experiment_name = experiment_name

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if force is None:
            self.force = False
        else:
            assert isinstance(force, bool)
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

        default = Default.get_global_default()

        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception('An Illumina run directory or file path has not been defined.')

        self.run_directory = Default.get_absolute_path(
            file_path=self.run_directory,
            default_path=Default.absolute_runs_illumina())

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception('The Illumina run directory {!r} does not exist.'.
                            format(self.run_directory))

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete('The Illumina run directory {!r} is not complete.'.
                                       format(self.run_directory))

        # Define an 'Intensities' directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Illumina Run Folder path.

        if self.intensity_directory:
            self.intensity_directory = Default.get_absolute_path(
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
            self.basecalls_directory = Default.get_absolute_path(
                file_path=self.basecalls_directory,
                default_path=self.intensity_directory)
        else:
            self.basecalls_directory = os.path.join(self.intensity_directory, 'BaseCalls')

        # Check that the BaseCalls directory exists.

        if not os.path.isdir(self.basecalls_directory):
            raise Exception('The BaseCalls directory {!r} does not exist.'.
                            format(self.basecalls_directory))

        self._irf = RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g. BSF_0000) is used as the prefix for archive BAM files.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = self._irf.run_parameters.get_experiment_name

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, self._irf.run_information.flow_cell))

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        # Call the run method of the super class after the project_name has been defined.

        super(PicardIlluminaRunFolder, self).run()

        return


class ExtractIlluminaBarcodesSheet(AnnotationSheet):
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


class IlluminaBasecallsToSamSheet(AnnotationSheet):
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


class FilePathExtractIlluminaCell(FilePath):

    def __init__(self, prefix):

        super(FilePathExtractIlluminaCell, self).__init__(prefix=prefix)

        self.sample_annotation_sheet_csv = prefix + '_samples.csv'

        return


class FilePathExtractIlluminaLane(FilePath):

    def __init__(self, prefix):

        super(FilePathExtractIlluminaLane, self).__init__(prefix=prefix)

        self.output_directory = prefix + '_output'
        self.samples_directory = prefix + '_samples'
        self.barcode_tsv = prefix + '_barcode.tsv'
        self.metrics_tsv = prefix + '_metrics.tsv'
        self.library_tsv = prefix + '_library.tsv'

        return


class ExtractIlluminaRunFolder(PicardIlluminaRunFolder):
    """The C{bsf.analyses.picard.ExtractIlluminaRunFolder} class represents to extract data from an Illumina Run Folder.

    The analysis is based on Picard ExtractIlluminaBarcodes and Picard IlluminaBasecallsToSam.

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
    @type samples_directory: str | unicode
    @ivar library_path: Library annotation file path
    @type library_path: str | unicode
    @ivar experiment_directory: Experiment directory
    @type experiment_directory: str | unicode
    @ivar mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_directory: str
    @ivar mode_file: Comma-separated list of file permission bit names according to the C{stat} module
    @type mode_file: str
    @ivar max_mismatches: Maximum number of mismatches
    @type max_mismatches: str
    @ivar min_base_quality: Minimum base quality
    @type min_base_quality: str
    @ivar sequencing_centre: Sequencing centre
    @type sequencing_centre: str
    @ivar lanes: Number of lanes on the flow cell
    @type lanes: int
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
            experiment_directory=None,
            mode_directory=None,
            mode_file=None,
            max_mismatches=None,
            min_base_quality=None,
            sequencing_centre=None,
            lanes=8,
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
        @type run_directory: str | unicode
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | unicode
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | unicode
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
        @param library_path: Library annotation file path
        @type library_path: str | unicode
        @param samples_directory: BSF samples directory
        @type samples_directory: str | unicode
        @param experiment_directory: Experiment directory
        @type experiment_directory: str | unicode
        @param mode_directory: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_directory: str
        @param mode_file: Comma-separated list of file permission bit names according to the C{stat} module
        @type mode_file: str
        @param max_mismatches: Maximum number of mismatches
        @type max_mismatches: str
        @param min_base_quality: Minimum base quality
        @type min_base_quality: str
        @param sequencing_centre: Sequencing centre
        @type sequencing_centre: str
        @param lanes: Number of lanes on the flow cell
        @type lanes: int
        @param vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
        @type vendor_quality_filter: dict[str, bool]
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

        if samples_directory is None:
            self.samples_directory = str()
        else:
            self.samples_directory = samples_directory

        if library_path is None:
            self.library_path = str()
        else:
            self.library_path = library_path

        if experiment_directory is None:
            self.experiment_directory = str()
        else:
            self.experiment_directory = experiment_directory

        # Can be None.
        self.mode_directory = mode_directory

        # Can be None.
        self.mode_file = mode_file

        if max_mismatches is None:
            self.max_mismatches = str()
        else:
            self.max_mismatches = max_mismatches

        if min_base_quality is None:
            self.min_base_quality = str()
        else:
            self.min_base_quality = min_base_quality

        if sequencing_centre is None:
            self.sequencing_centre = str()
        else:
            self.sequencing_centre = sequencing_centre

        if lanes is None:
            self.lanes = int()
        else:
            assert isinstance(lanes, int)
            self.lanes = lanes

        if vendor_quality_filter is None:
            self.vendor_quality_filter = dict()
        else:
            self.vendor_quality_filter = vendor_quality_filter

        return

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

        # For the moment, the experiment_directory cannot be configured.
        # It is automatically assembled from self.samples_directory and self.project_name by the run() method.

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
            self.max_mismatches = configuration.config_parser.get(section=section, option=option)

        # Get the minimum base quality.

        option = 'min_base_quality'
        if configuration.config_parser.has_option(section=section, option=option):
            self.min_base_quality = configuration.config_parser.get(section=section, option=option)

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

        def run_sample_file_name(sample_name):
            """Private function to format sample-specific BAM file names (i.e. project_lane#sample.bam).

            @param sample_name:
            @type sample_name: str | unicode
            @return:
            @rtype: str | unicode
            """
            return self.project_name + '_' + lane_str + '#' + sample_name + '.bam'

        # Start of the run() method body.

        super(ExtractIlluminaRunFolder, self).run()

        default = Default.get_global_default()

        self.samples_directory = Default.get_absolute_path(
            file_path=self.samples_directory,
            default_path=Default.absolute_samples())

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception('The ExtractIlluminaRunFolder samples_directory {!r} does not exist.'.
                            format(self.samples_directory))

        self.experiment_directory = os.path.join(self.samples_directory, self.project_name)

        # Get sequencing centre information.

        if not self.sequencing_centre:
            self.sequencing_centre = default.operator_sequencing_centre

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if self._irf.run_parameters.get_flow_cell_type not in self.vendor_quality_filter:
            raise Exception('Flow cell chemistry type {!r} not defined.'.
                            format(self._irf.run_parameters.get_flow_cell_type))

        # Get the library annotation sheet.
        # The library annotation sheet is deliberately not passed in via sas_file,
        # as the Analysis.run() method reads that option into a BSF Collection object.

        self.library_path = Default.get_absolute_path(file_path=self.library_path)

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = os.path.normpath(path=self.library_path)

        if not os.path.exists(path=self.library_path):
            raise Exception('Library annotation file {!r} does not exist.'.
                            format(self.library_path))

        stage_lane = self.get_stage(name=self.stage_name_lane)
        stage_cell = self.get_stage(name=self.stage_name_cell)

        # Read the LibraryAnnotationSheet and populate a flow cell dict indexed on the lane number ...

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

        flow_cell_dict = dict()
        """ @type flow_cell_dict: dict[str, list[dict[str, str | unicode]]] """

        for row_dict in library_annotation_sheet.row_dicts:
            if row_dict['lane'] not in flow_cell_dict:
                flow_cell_dict[row_dict['lane']] = list()
            flow_cell_dict[row_dict['lane']].append(row_dict)

        file_path_cell = FilePathExtractIlluminaCell(prefix=self.project_name)

        # Create a Sample Annotation Sheet in the project directory and
        # eventually transfer it into the experiment_directory.
        sample_annotation_sheet = SampleAnnotationSheet(
            file_path=os.path.join(self.project_directory, file_path_cell.sample_annotation_sheet_csv))

        # For each lane in the flow_cell_dict ...
        # TODO: For the moment this depends on the lanes (keys) defined in the LibraryAnnotationSheet.
        # Not all lanes may thus get extracted.
        # TODO: For NextSeq instruments, it would be sufficient to require annotation for only lane one and
        # copy information to lanes two to four internally.

        cell_dependency_list = list()

        lane_str_list = flow_cell_dict.keys()
        lane_str_list.sort(cmp=lambda x, y: cmp(x, y))

        for lane_str in lane_str_list:
            # The lane_str represents the lane number as a Python str.
            prefix_lane = '_'.join((self.project_name, lane_str))

            file_path_lane = FilePathExtractIlluminaLane(prefix=prefix_lane)

            # BARCODE_FILE
            eib_sheet = ExtractIlluminaBarcodesSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.barcode_tsv))

            # LIBRARY_PARAMS
            ibs_sheet = IlluminaBasecallsToSamSheet(
                file_path=os.path.join(self.project_directory, file_path_lane.library_tsv))

            # Initialise a list of barcode sequence lengths.
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
                                'The length ({}) of barcode {} {!r} does not match '
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
                        run_sample_file_name(sample_name=row_dict['sample_name'])),
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
                        os.path.basename(self.experiment_directory),
                        file_path_lane.samples_directory,
                        run_sample_file_name(sample_name=row_dict['sample_name'])),
                    'Reads2 Name': '',
                    'Reads2 File': '',
                })

            # The IlluminaBasecallsToSamSheet needs adjusting ...

            if len(ibs_sheet.row_dicts) == 1 and len(ibs_sheet.row_dicts[0]['BARCODE_1']) == 0 and len(
                    ibs_sheet.row_dicts[0]['BARCODE_2']) == 0:
                # ... if a single sample, but neither BARCODE_1 nor BARCODE_2 were defined,
                # BARCODE_1 needs setting to 'N'.
                ibs_sheet.row_dicts[0]['BARCODE_1'] = 'N'
            else:
                # ... in all other cases as a last row for unmatched barcode sequences needs adding.
                ibs_sheet.row_dicts.append({
                    'BARCODE_1': 'N',
                    'BARCODE_2': '',
                    'OUTPUT': os.path.join(
                        file_path_lane.samples_directory,
                        run_sample_file_name(sample_name='0')),
                    'SAMPLE_ALIAS': 'Unmatched',
                    'LIBRARY_NAME': ibs_sheet.row_dicts[0]['LIBRARY_NAME'],
                })

            # Calculate the read structure string from the IRF and the bc_length_list above ...

            read_structure = str()
            index_read_index = 0  # Number of index reads.
            # Instantiate and sort a new list of RunInformationRead objects.
            run_information_read_list = list(self._irf.run_information.reads)
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

            # Further adjust the IlluminaBaseCallsToSamSheet and remove any BARCODE_N columns not represented
            # in the read structure.

            for index in range(0, 1 + 1):
                if index + 1 > index_read_index:
                    # Remove the 'BARCODE_N' filed from the list of field names.
                    if 'BARCODE_' + str(index + 1) in ibs_sheet.field_names:
                        ibs_sheet.field_names.remove('BARCODE_' + str(index + 1))
                    # Remove the 'BARCODE_N' entry form each row dict object, since csv.DictWriter requires it.
                    for row_dict in ibs_sheet.row_dicts:
                        row_dict.pop('BARCODE_' + str(index + 1), None)

            # Write the lane-specific Picard ExtractIlluminaBarcodesSheet and Picard IlluminaBasecallsToSamSheet.

            if index_read_index > 0:
                eib_sheet.to_file_path()

            ibs_sheet.to_file_path()

            # Create a Runnable and Executable for the lane stage.

            runnable_lane = self.add_runnable(
                runnable=Runnable(
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
                runnable_step=RunnableStepMakeDirectory(
                    name='make_output_directory',
                    directory_path=file_path_lane.output_directory))

            # Create a samples_directory in the project_directory.

            runnable_lane.add_runnable_step(
                runnable_step=RunnableStepMakeDirectory(
                    name='make_samples_directory',
                    directory_path=file_path_lane.samples_directory))

            # Create a RunnableStep for Picard ExtractIlluminaBarcodes, only if index (barcode) reads are present.

            if index_read_index > 0:
                runnable_step = runnable_lane.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='picard_extract_illumina_barcodes',
                        java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx2G',
                        picard_classpath=self.classpath_picard,
                        picard_command='ExtractIlluminaBarcodes'))
                """ @type runnable_step: RunnableStepPicard """
                runnable_step.add_picard_option(key='BASECALLS_DIR', value=self.basecalls_directory)
                runnable_step.add_picard_option(key='OUTPUT_DIR', value=file_path_lane.output_directory)
                runnable_step.add_picard_option(key='LANE', value=lane_str)
                runnable_step.add_picard_option(key='READ_STRUCTURE', value=read_structure)
                runnable_step.add_picard_option(key='BARCODE_FILE', value=file_path_lane.barcode_tsv)
                runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_lane.metrics_tsv)
                if self.max_mismatches:
                    # Maximum mismatches for a barcode to be considered a match. Default value: '1'.
                    runnable_step.add_picard_option(key='MAX_MISMATCHES', value=self.max_mismatches)
                # MIN_MISMATCH_DELTA: Minimum difference between number of mismatches in the best and
                # second best barcodes for a barcode to be considered a match. Default value: '1'.
                # MAX_NO_CALLS Maximum allowable number of no-calls in a barcode read before it is
                # considered unmatchable. Default value: '2'.
                if self.min_base_quality:
                    # Minimum base quality. Any barcode bases falling below this quality will be considered
                    # a mismatch even in the bases match. Default value: '0'.
                    runnable_step.add_picard_option(key='MINIMUM_BASE_QUALITY', value=self.min_base_quality)
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

            # Picard IlluminaBasecallsToSam

            # Create a RunnableStep for Picard IlluminaBasecallsToSam.

            runnable_step = runnable_lane.add_runnable_step(
                runnable_step=RunnableStepPicard(
                    name='picard_illumina_basecalls_to_sam',
                    java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    picard_classpath=self.classpath_picard,
                    picard_command='IlluminaBasecallsToSam'))
            """ @type runnable_step: RunnableStepPicard """
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
            # NOTE: The only date format that seems to work is mm/dd/yyyy. Why?
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

            # Create the experiment directory if it does not exist already.

            runnable_lane.add_runnable_step(
                runnable_step=RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=self.experiment_directory))

            # Move the samples directory into the experiment directory.

            runnable_lane.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='move_samples_directory',
                    source_path=file_path_lane.samples_directory,
                    target_path=self.experiment_directory))

            # Move the metrics file into the experiment directory.

            if index_read_index > 0:
                runnable_lane.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_metrics_tsv',
                        source_path=file_path_lane.metrics_tsv,
                        target_path=self.experiment_directory))

        # Finally, write the flow cell-specific SampleAnnotationSheet to the internal file path.

        sample_annotation_sheet.to_file_path()

        # Create a flow-cell specific Runnable.

        runnable_cell = self.add_runnable(
            runnable=Runnable(
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
                runnable_step=RunnableStepMove(
                    name='move_sample_annotation',
                    source_path=file_path_cell.sample_annotation_sheet_csv,
                    target_path=self.experiment_directory))

        # Change directory and file access permissions.

        runnable_cell.add_runnable_step(
            runnable_step=RunnableStepChangeMode(
                name='chmod',
                file_path=self.experiment_directory,
                mode_directory=self.mode_directory,
                mode_file=self.mode_file))

        return


class FilePathCollectHiSeqXPfFailMetricsLane(FilePath):

    def __init__(self, prefix):

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

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the Stage.
            runnable_lane = self.add_runnable(
                runnable=Runnable(
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
                runnable_step=RunnableStepPicard(
                    name='illumina_to_bam',
                    java_temporary_path=runnable_lane.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    picard_classpath=self.classpath_picard,
                    picard_command='CollectHiSeqXPfFailMetrics'))
            """ @type runnable_step: RunnableStepPicard """
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


class FilePathDownsampleSam(FilePath):

    def __init__(self, prefix):

        super(FilePathDownsampleSam, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.downsampled_bam = prefix + '_downsampled.bam'

        return


class DownsampleSam(Analysis):
    """The C{bsf.analyses.picard.DownsampleSam} class represents the logic to run the Picard DownsampleSam analysis.

    Attributes:

    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
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
        @type classpath_picard: str | unicode
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

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
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

        default = Default.get_global_default()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        run_read_comparisons()

        # Picard DownsampleSam

        stage_picard_dss = self.get_stage(name='picard_downsample_sam')

        for sample in self.sample_list:
            if self.debug > 0:
                print self, 'Sample name:', sample.name
                print sample.trace(level=1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print self, 'PairedReads name:', paired_reads.get_name()

                    # Apply some sanity checks.

                    if paired_reads.reads_2 and not paired_reads.reads_1:
                        raise Exception('PairedReads object with reads_1 but no reads_2 object.', UserWarning)

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
                        runnable=Runnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory,
                            file_path_object=file_path_read_group))

                    # Create an Executable for running the Picard SamToFastq Runnable.

                    self.set_stage_runnable(stage=stage_picard_dss, runnable=runnable_picard_dss)

                    # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                    if False:
                        runnable_picard_dss.add_runnable_step(
                            runnable_step=RunnableStepMakeDirectory(
                                name='mkdir',
                                directory_path=file_path_read_group.output_directory))

                    # Create a new RunnableStep for the Picard DownsampleSam program.

                    runnable_step = runnable_picard_dss.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='picard_donwsample_sam',
                            java_temporary_path=runnable_picard_dss.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx2G',
                            picard_classpath=self.classpath_picard,
                            picard_command='DownsampleSam'))
                    """ @type runnable_step: RunnableStepPicard """
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
                    # OPTIONS_FILE

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'picard_downsample_sam_samples.csv'))),
            name='_'.join((self.project_name, 'picard_downsample_sam')))

        annotation_sheet.to_file_path()

        return


class FilePathSamToFastq(FilePath):

    def __init__(self, prefix):

        super(FilePathSamToFastq, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class SamToFastq(Analysis):
    """The C{bsf.analyses.picard.SamToFastq} class represents the logic to run the Picard SamToFastq analysis.

    Attributes:

    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar include_non_pf_reads: Include non-pass filer reads
    @type include_non_pf_reads: bool
    """

    name = 'Picard SamToFastq Analysis'
    prefix = 'picard_sam_to_fastq'

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
            include_non_pf_reads=False):
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
        @type classpath_picard: str | unicode
        @param include_non_pf_reads: Include non-pass filer reads
        @type include_non_pf_reads: bool
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

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if include_non_pf_reads is None:
            self.include_non_pass_filter_reads = False
        else:
            assert isinstance(include_non_pf_reads, bool)
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

        default = Default.get_global_default()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        run_read_comparisons()

        prune_sas_dependencies = list()

        # Picard SamToFastq

        stage_picard_stf = self.get_stage(name='picard_sam_to_fastq')

        for sample in self.sample_list:
            if self.debug > 0:
                print self, 'Sample name:', sample.name
                print sample.trace(level=1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print self, 'PairedReads name:', paired_reads.get_name()

                    # Apply some sanity checks.

                    if paired_reads.reads_2 and not paired_reads.reads_1:
                        raise Exception('PairedReads object with reads_1 but no reads_2 object.', UserWarning)

                    reads = paired_reads.reads_1
                    if reads.file_path.endswith('.bam'):
                        bam_file_path = reads.file_path
                        prefix_picard_stf = '_'.join((stage_picard_stf.name, paired_reads_name))

                        file_path_picard_stf = FilePathSamToFastq(prefix=prefix_picard_stf)

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
                                                       read_group_dict.keys()))
                            if read_group_dict == alignment_file.header['RG'][0]:
                                # Use the '==' rather than the 'is' operator, since dictionaries do not seem to be
                                # at the same memory address.
                                # For the first read group, modify the PairedReads object in place.
                                paired_reads.read_group = '\\t'.join(read_group_list)
                                paired_reads.reads_1.name = platform_unit + '_1'
                                paired_reads.reads_1.file_path = os.path.join(
                                    self.project_directory,
                                    file_path_picard_stf.output_directory,
                                    platform_unit + '_1.fastq')
                                paired_reads.reads_2 = Reads(
                                    name=platform_unit + '_2',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_picard_stf.output_directory,
                                        platform_unit + '_2.fastq'),
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
                                        file_path_picard_stf.output_directory,
                                        platform_unit + '_1.fastq'),
                                    file_type=paired_reads.reads_1.file_type,
                                    lane=paired_reads.reads_1.lane,
                                    read='R1',
                                    chunk=paired_reads.reads_1.chunk)
                                reads_2 = Reads(
                                    name=platform_unit + '_2',
                                    file_path=os.path.join(
                                        self.project_directory,
                                        file_path_picard_stf.output_directory,
                                        platform_unit + '_2.fastq'),
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

                        runnable_picard_stf = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_picard_stf,
                                code_module='bsf.runnables.generic',
                                working_directory=self.project_directory,
                                file_path_object=file_path_picard_stf))

                        # Create an Executable for running the Picard SamToFastq Runnable.

                        self.set_stage_runnable(stage=stage_picard_stf, runnable=runnable_picard_stf)

                        # Record the Executable.name for the prune_sas dependency.

                        prune_sas_dependencies.append(runnable_picard_stf.name)

                        # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                        runnable_picard_stf.add_runnable_step(
                            runnable_step=RunnableStepMakeDirectory(
                                name='mkdir',
                                directory_path=file_path_picard_stf.output_directory))

                        # Create a new RunnableStep for the Picard SamToFastq program.

                        runnable_step = runnable_picard_stf.add_runnable_step(
                            runnable_step=RunnableStepPicard(
                                name='picard_sam_to_fastq',
                                java_temporary_path=runnable_picard_stf.get_relative_temporary_directory_path,
                                java_heap_maximum='Xmx2G',
                                picard_classpath=self.classpath_picard,
                                picard_command='SamToFastq'))
                        """ @type runnable_step: RunnableStepPicard """
                        runnable_step.add_picard_option(key='INPUT', value=bam_file_path)
                        runnable_step.add_picard_option(key='OUTPUT_PER_RG', value='true')
                        runnable_step.add_picard_option(
                            key='OUTPUT_DIR',
                            value=file_path_picard_stf.output_directory)
                        # RE_REVERSE
                        # INTERLEAVE
                        if self.include_non_pass_filter_reads:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='true')
                        else:
                            runnable_step.add_picard_option(key='INCLUDE_NON_PF_READS', value='false')
                        # CLIPPING_ATTRIBUTE
                        # CLIPPING_ACTION
                        # READ1_TRIM
                        # READ1_MAX_BASES_TO_WRITE
                        # READ2_TRIM
                        # READ2_MAX_BASES_TO_WRITE
                        # INCLUDE_NON_PRIMARY_ALIGNMENTS
                        runnable_step.add_picard_option(
                            key='TMP_DIR',
                            value=runnable_picard_stf.get_relative_temporary_directory_path)
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
                        # OPTIONS_FILE

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'picard_sam_to_fastq_original.csv'))),
            name='_'.join((self.project_name, 'picard_sam_to_fastq')))

        annotation_sheet.to_file_path()

        # Create a Runnable for pruning the sample annotation sheet.

        prefix_prune_sas = '_'.join((stage_picard_stf.name, self.project_name))

        file_path_prune_sas = FilePathSamToFastq(prefix=prefix_prune_sas)

        runnable_prune_sas = self.add_runnable(
            runnable=Runnable(
                name=prefix_prune_sas,
                code_module='bsf.runnables.picard_sam_to_fastq_sample_sheet',
                working_directory=self.project_directory,
                file_path_object=file_path_prune_sas))

        # Create an Executable for running the Runnable for pruning the sample annotation sheet.

        executable_prune_sas = self.set_stage_runnable(
            stage=stage_picard_stf,
            runnable=runnable_prune_sas)
        executable_prune_sas.dependencies.extend(prune_sas_dependencies)

        # Create a new RunnableStep.

        prune_sas = runnable_prune_sas.add_runnable_step(
            runnable_step=RunnableStep(
                name='prune_sample_annotation_sheet'))

        prune_sas.add_option_long(key='sas_path', value=annotation_sheet.file_path)

        return annotation_sheet
