"""bsf.analyses.illumina_to_bam_tools

A package of classes and methods supporting analyses of the Illumina2Bam-Tools package.
"""

#
# Copyright 2013 - 2015 Michael K. Schuster
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


import errno
import os
import string
import warnings

from bsf import Analysis, Command, Configuration, Default, DRMS, Executable, Runnable, RunnableStep
from bsf.analyses.illumina_run_folder import IlluminaRunFolderRestore
from bsf.annotation import BamIndexDecoderSheet, LibraryAnnotationSheet, SampleAnnotationSheet
from bsf.illumina import RunFolder, RunFolderNotComplete


class IlluminaToBam(Analysis):
    """The C{IlluminaToBam} class represents the logic to convert Illumina BCL to a BAM or SAM files.

    Attributes:
    @cvar drms_name_illumina_to_bam: C{DRMS.name} for the C{IlluminaToBam} C{Analysis} stage
    @type drms_name_illumina_to_bam: str
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | unicode
    @ivar intensity_directory: File path to the I{Intensities} directory,
        defaults to I{illumina_run_folder/Data/Intensities}
    @type intensity_directory: str | unicode
    @ivar basecalls_directory: File path to the I{BaseCalls} directory,
        defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
    @type basecalls_directory: str | unicode
    @ivar experiment_name: Experiment name (i.e. flow-cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str
    @ivar sequencing_centre: Sequencing centre
    @type sequencing_centre: str
    @ivar sequences_directory: Sequences directory to store archive BAM files
    @type sequences_directory: str | unicode
    @ivar experiment_directory: Experiment-specific directory
    @type experiment_directory: str | unicode
    @ivar classpath_illumina2bam: Illumina2Bam tools Java Archive (JAR) class path directory
    @type classpath_illumina2bam: str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool
    """

    drms_name_illumina_to_bam = 'illumina_to_bam'

    @classmethod
    def get_prefix_illumina_to_bam(cls, project_name, lane):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_illumina_to_bam, project_name, lane), sep='_')

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{IlluminaToBam} object from a UNIX-style configuration file via the C{Configuration} class.

        @param config_path: UNIX-style configuration file
        @type config_path: str | unicode
        @return: IlluminaToBam
        @rtype: IlluminaToBam
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{IlluminaToBam} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{IlluminaToBam}
        @rtype: IlluminaToBam
        """

        assert isinstance(configuration, Configuration)

        itb = cls(configuration=configuration)

        # A "bsf.analyses.IlluminaToBamTools.IlluminaToBam" section specifies defaults
        # for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        itb.set_configuration(itb.configuration, section=section)

        return itb

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 run_directory=None, intensity_directory=None, basecalls_directory=None,
                 experiment_name=None, sequencing_centre=None,
                 sequences_directory=None, experiment_directory=None,
                 classpath_illumina2bam=None, classpath_picard=None,
                 force=False):
        """Initialise a C{IlluminaToBam} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{Analysis}-wide project directory,
            normally under the C{Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{Analysis}-wide genome directory,
            normally under the C{Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param drms_list: Python C{list} of C{DRMS} objects
        @type drms_list: list
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{tuple} objects of C{Sample} objects
        @type comparisons: dict
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode
        @param intensity_directory: File path to the I{Intensities} directory,
            defaults to I{illumina_run_folder/Data/Intensities}
        @type intensity_directory: str | unicode
        @param basecalls_directory: File path to the I{BaseCalls} directory,
            defaults to I{illumina_run_folder/Data/Intensities/BaseCalls}
        @type basecalls_directory: str | unicode
        @param experiment_name: Experiment name (i.e. flow-cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param sequencing_centre: Sequencing centre
        @type sequencing_centre: str
        @param sequences_directory: Sequences directory to store archive BAM files
        @type sequences_directory: str | unicode
        @param experiment_directory: Experiment-specific directory
        @type experiment_directory: str | unicode
        @param classpath_illumina2bam: Illumina2Bam tools Java Archive (JAR) class path directory
        @type classpath_illumina2bam: str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
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
            drms_list=drms_list,
            collection=collection,
            comparisons=comparisons,
            samples=samples)

        # Sub-class specific ...

        if run_directory:
            self.run_directory = run_directory
        else:
            self.run_directory = str()

        if intensity_directory:
            self.intensity_directory = intensity_directory
        else:
            self.intensity_directory = str()

        if basecalls_directory:
            self.basecalls_directory = basecalls_directory
        else:
            self.basecalls_directory = str()

        if experiment_name:
            self.experiment_name = experiment_name
        else:
            self.experiment_name = str()

        if sequencing_centre:
            self.sequencing_centre = sequencing_centre
        else:
            self.sequencing_centre = str()

        if sequences_directory:
            self.sequences_directory = sequences_directory
        else:
            self.sequences_directory = str()

        if experiment_directory:
            self.experiment_directory = experiment_directory
        else:
            self.experiment_directory = str()

        if classpath_illumina2bam:
            self.classpath_illumina2bam = classpath_illumina2bam
        else:
            self.classpath_illumina2bam = str()

        if classpath_picard:
            self.classpath_picard = classpath_picard
        else:
            self.classpath_picard = str()

        self.force = force

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{IlluminaToBam} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        super(IlluminaToBam, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get Illumina Run Folder information.

        if configuration.config_parser.has_option(section=section, option='illumina_run_folder'):
            self.run_directory = configuration.config_parser.get(
                section=section,
                option='illumina_run_folder')

        if configuration.config_parser.has_option(section=section, option='intensity_directory'):
            self.intensity_directory = configuration.config_parser.get(
                section=section,
                option='intensity_directory')

        if configuration.config_parser.has_option(section=section, option='basecalls_directory'):
            self.basecalls_directory = configuration.config_parser.get(
                section=section,
                option='basecalls_directory')

        # Get the experiment name.

        if configuration.config_parser.has_option(section=section, option='experiment_name'):
            self.experiment_name = configuration.config_parser.get(
                section=section,
                option='experiment_name')

        # Get sequencing centre information.

        if configuration.config_parser.has_option(section=section, option='sequencing_centre'):
            self.sequencing_centre = configuration.config_parser.get(
                section=section,
                option='sequencing_centre')

        # Get the sequences directory information.

        if configuration.config_parser.has_option(section=section, option='sequences_directory'):
            self.sequences_directory = configuration.config_parser.get(
                section=section,
                option='sequences_directory')

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_illumina2bam'):
            self.classpath_illumina2bam = configuration.config_parser.get(
                section=section,
                option='classpath_illumina2bam')

        # Get the Picard tools Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_picard'):
            self.classpath_picard = configuration.config_parser.get(
                section=section,
                option='classpath_picard')

        if configuration.config_parser.has_option(section=section, option='force'):
            self.force = configuration.config_parser.getboolean(
                section=section,
                option='force')

    def run(self):
        """Run this C{IlluminaToBam} C{Analysis}.

        Convert an Illumina flow-cell into lane-specific archive BAM files.

        To convert an Illumina flow cell, Illumina2bam is run first, setting the SAM Read Group (@RG)
        library name (LB) and sample name (SM) to 'flow-cell identifier.lane'.
        The resulting archive BAM file is then sorted by query name with Picard SortSam.
        """

        default = Default.get_global_default()

        # Define an Illumina Run Folder directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception('An Illumina run directory or file path has not been defined.')

        self.run_directory = os.path.expanduser(path=self.run_directory)
        self.run_directory = os.path.expandvars(path=self.run_directory)

        if not os.path.isabs(self.run_directory):
            self.run_directory = os.path.join(Default.absolute_runs_illumina(), self.run_directory)

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The Illumina run directory {!r} does not exist.'.format(self.run_directory))

        # Check that the Illumina Run Folder is complete.

        if not (os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(
                'The Illumina run directory {!r} is not complete.'.format(self.run_directory))

        # Define an 'Intensities' directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend the Illumina Run Folder path.

        if self.intensity_directory:
            intensity_directory = str(self.intensity_directory)
            intensity_directory = os.path.expandvars(intensity_directory)
            intensity_directory = os.path.expanduser(intensity_directory)
            if not os.path.isabs(intensity_directory):
                intensity_directory = os.path.join(self.run_directory, intensity_directory)
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
            basecalls_directory = str(self.basecalls_directory)
            basecalls_directory = os.path.expanduser(basecalls_directory)
            basecalls_directory = os.path.expandvars(basecalls_directory)
            if not os.path.isabs(basecalls_directory):
                basecalls_directory = os.path.join(intensity_directory, basecalls_directory)
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
            self.project_name = string.join(words=(self.experiment_name, irf.run_information.flow_cell), sep='_')

        # Get sequencing centre information.

        if not self.sequencing_centre:
            self.sequencing_centre = default.operator_sequencing_centre

        # Define the sequences directory in which to create the experiment directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,

        self.sequences_directory = os.path.expanduser(path=self.sequences_directory)
        self.sequences_directory = os.path.expandvars(path=self.sequences_directory)

        if not os.path.isabs(self.sequences_directory):
            self.sequences_directory = os.path.join(Default.absolute_sequences(), self.sequences_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the sequences_directory has to exist.

        if not os.path.isdir(self.sequences_directory):
            raise Exception(
                'The IlluminaToBam sequences_directory {!r} does not exist.'.format(self.sequences_directory))

        self.experiment_directory = os.path.join(self.sequences_directory, self.project_name)

        if not os.path.isdir(self.experiment_directory):
            try:
                os.makedirs(self.experiment_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        if not self.classpath_illumina2bam:
            self.classpath_illumina2bam = default.classpath_illumina2bam

        # Get the Picard tools Java Archive (JAR) class path directory

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        # Call the run method of the super class after the project_name has been defined.

        super(IlluminaToBam, self).run()

        drms_illumina_to_bam = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_illumina_to_bam,
            working_directory=self.project_directory,
            analysis=self))

        for lane in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):

            lane_str = str(lane)

            file_path_dict = dict(
                illumina_directory=self.run_directory,  # contains full path information
                sequences_directory=self.sequences_directory,  # contains full path information
                experiment_directory=self.experiment_directory,  # contains full path information
                sorted_bam=os.path.join(
                    self.experiment_directory,
                    string.join((self.project_name, lane_str, 'sorted.bam'), sep='_')),
                sorted_md5=os.path.join(
                    self.experiment_directory,
                    string.join((self.project_name, lane_str, 'sorted.bam.md5'), sep='_')),
                unsorted_bam=os.path.join(
                    self.experiment_directory,
                    string.join((self.project_name, lane_str, 'unsorted.bam'), sep='_')),
                unsorted_md5=os.path.join(
                    self.experiment_directory,
                    string.join((self.project_name, lane_str, 'unsorted.bam.md5'), sep='_'))
            )

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            runnable_illumina_to_bam = self.add_runnable(runnable=Runnable(
                name=self.get_prefix_illumina_to_bam(project_name=self.project_name, lane=lane_str),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict))

            # TODO: The Runnable class could have dependencies just like the Executable class so that they could be
            # passed on upon creation of the Executable from the Runnable via Executable.from_analysis_runnable().
            executable_illumina_to_bam = drms_illumina_to_bam.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_illumina_to_bam.name))

            executable_illumina_to_bam.dependencies.append(
                IlluminaRunFolderRestore.get_prefix_compress_base_calls(project_name=self.project_name, lane=lane_str))

            # Only submit this Executable if the final result file does not exist.
            if (os.path.exists(file_path_dict['sorted_md5'])
                    and os.path.getsize(file_path_dict['sorted_md5'])):
                executable_illumina_to_bam.submit = False

            # Run Illumina2Bam tools Illumina2bam.

            java_process = runnable_illumina_to_bam.add_runnable_step(runnable_step=RunnableStep(
                name='illumina_to_bam',
                program='java',
                sub_command=Command(command=str())))

            java_process.add_switch_short(key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_illumina2bam, 'Illumina2bam.jar'))
            java_process.add_switch_short(key='Xmx4G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=runnable_illumina_to_bam.get_relative_temporary_directory_path)

            sub_command = java_process.sub_command

            if self.intensity_directory:
                # Only set the RUN_FOLDER option, if a separate 'Intensities' directory has been configured.
                # The default is to use the directory two up from the INTENSITY_DIR.
                sub_command.add_option_pair(
                    key='RUN_FOLDER',
                    value=self.run_directory)
            sub_command.add_option_pair(
                key='INTENSITY_DIR',
                value=intensity_directory)
            if self.basecalls_directory:
                # Only set the BASECALLS_DIR option, if a separate 'BaseCalls' directory has been configured.
                # The default is to use the 'BaseCalls' directory under the INTENSITY_DIR.
                sub_command.add_option_pair(
                    key='BASECALLS_DIR',
                    value=basecalls_directory)
            sub_command.add_option_pair(
                key='LANE',
                value=lane_str)
            sub_command.add_option_pair(
                key='OUTPUT',
                value=file_path_dict['unsorted_bam'])
            sub_command.add_option_pair(
                key='GENERATE_SECONDARY_BASE_CALLS',
                value='false')
            sub_command.add_option_pair(
                key='PF_FILTER',
                value='false')
            sub_command.add_option_pair(
                key='READ_GROUP_ID',
                value=string.join((irf.run_information.flow_cell, lane_str), sep='_'))
            # SAMPLE_ALIAS
            sub_command.add_option_pair(
                key='LIBRARY_NAME',
                value=string.join((irf.run_information.flow_cell, lane_str), sep='_'))
            # STUDY_NAME
            # PLATFORM_UNIT
            # RUN_START_DATE
            sub_command.add_option_pair(
                key='SEQUENCING_CENTER',
                value=self.sequencing_centre)
            # PLATFORM
            # FIRST_TILE
            # TILE_LIMIT
            # BARCODE_SEQUENCE_TAG_NAME
            # BARCODE_QUALITY_TAG_NAME
            # SECOND_BARCODE_SEQUENCE_TAG_NAME
            # SECOND_BARCODE_QUALITY_TAG_NAME
            # FIRST_CYCLE
            # FINAL_CYCLE
            # FIRST_INDEX_CYCLE
            # FINAL_INDEX_CYCLE
            sub_command.add_option_pair(
                key='TMP_DIR',
                value=runnable_illumina_to_bam.get_relative_temporary_directory_path)
            sub_command.add_option_pair(
                key='VERBOSITY',
                value='WARNING')
            # QUIET
            # VALIDATION_STRINGENCY
            # COMPRESSION_LEVEL
            sub_command.add_option_pair(
                key='MAX_RECORDS_IN_RAM',
                value='2000000')
            sub_command.add_option_pair(
                key='CREATE_INDEX',
                value='false')
            sub_command.add_option_pair(
                key='CREATE_MD5_FILE',
                value='true')
            # OPTIONS_FILE

            # Run Picard SortSam

            java_process = runnable_illumina_to_bam.add_runnable_step(runnable_step=RunnableStep(
                name='picard_sort_sam',
                program='java',
                sub_command=Command(command=str()),
                obsolete_file_path_list=[
                    file_path_dict['unsorted_bam'],
                    file_path_dict['unsorted_md5']
                ]))

            java_process.add_switch_short(key='d64')
            java_process.add_option_short(key='jar', value=os.path.join(self.classpath_picard, 'SortSam.jar'))
            java_process.add_switch_short(key='Xmx4G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=runnable_illumina_to_bam.get_relative_temporary_directory_path)

            sub_command = java_process.sub_command

            sub_command.add_option_pair(
                key='INPUT',
                value=file_path_dict['unsorted_bam'])
            sub_command.add_option_pair(
                key='OUTPUT',
                value=file_path_dict['sorted_bam'])
            sub_command.add_option_pair(
                key='SORT_ORDER',
                value='queryname')
            sub_command.add_option_pair(
                key='TMP_DIR',
                value=runnable_illumina_to_bam.get_relative_temporary_directory_path)
            sub_command.add_option_pair(
                key='VERBOSITY',
                value='WARNING')
            # QUIET
            # VALIDATION_STRINGENCY
            # COMPRESSION_LEVEL
            sub_command.add_option_pair(
                key='MAX_RECORDS_IN_RAM',
                value='2000000')
            sub_command.add_option_pair(
                key='CREATE_INDEX',
                value='false')
            sub_command.add_option_pair(
                key='CREATE_MD5_FILE',
                value='true')
            # OPTIONS_FILE


class BamIndexDecoder(Analysis):
    """The C{BamIndexDecoder} class represents the logic to decode sequence archive BAM files into
    sample-specific BAM files.

    Attributes:
    @cvar drms_name_bam_index_decoder: C{DRMS.name} for the C{BamIndexDecoder} C{Analysis} stage
    @type drms_name_bam_index_decoder: str
    @ivar library_path: Library annotation file path
    @type library_path: str | unicode
    @ivar sequences_directory: BSF sequences directory
    @type sequences_directory: str | unicode
    @ivar samples_directory: BSF samples directory
    @type samples_directory: str | unicode
    @ivar experiment_directory: Experiment directory
    @type experiment_directory: str | unicode
    @ivar classpath_illumina2bam: Illumina2Bam tools Java Archive (JAR) class path directory
    @type classpath_illumina2bam: str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar lanes: Number of lanes on the flow cell
    @type lanes: int
    @ivar force: Force de-multiplexing with a Library Annotation sheet failing validation
    @type force: bool
    """

    drms_name_bam_index_decoder = 'bam_index_decoder'

    @classmethod
    def get_prefix_bam_index_decoder(cls, project_name, lane):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_bam_index_decoder, project_name, lane), sep='_')

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{BamIndexDecoder} object from a UNIX-style configuration file via the C{Configuration} class.

        @param config_path: UNIX-style configuration file
        @type config_path: str | unicode
        @return: C{BamIndexDecoder}
        @rtype: BamIndexDecoder
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{BamIndexDecoder} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{BamIndexDecoder}
        @rtype: BamIndexDecoder
        """

        assert isinstance(configuration, Configuration)

        itb = cls(configuration=configuration)

        # A "bsf.analyses.IlluminaToBamTools.BamIndexDecoder" section specifies defaults
        # for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        itb.set_configuration(itb.configuration, section=section)

        return itb

    def __init__(self, configuration=None, project_name=None, genome_version=None, input_directory=None,
                 output_directory=None, project_directory=None, genome_directory=None, e_mail=None, debug=0,
                 drms_list=None, collection=None, comparisons=None, samples=None, library_path=None,
                 sequences_directory=None, samples_directory=None, experiment_directory=None,
                 classpath_illumina2bam=None, classpath_picard=None, lanes=8, force=False):
        """Initialise a C{BamIndexDecoder} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{Analysis}-wide project directory,
            normally under the C{Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{Analysis}-wide genome directory,
            normally under the C{Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param drms_list: Python C{list} of C{DRMS} objects
        @type drms_list: list
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{tuple} objects of C{Sample} objects
        @type comparisons: dict
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        @param library_path: Library annotation file path
        @type library_path: str | unicode
        @param sequences_directory: BSF sequences directory
        @type sequences_directory: str | unicode
        @param samples_directory: BSF samples directory
        @type samples_directory: str | unicode
        @param experiment_directory: Experiment directory
        @type experiment_directory: str | unicode
        @param classpath_illumina2bam: Illumina2Bam tools Java Archive (JAR) class path directory
        @type classpath_illumina2bam: str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param lanes: Number of lanes on the flow cell
        @type lanes: int
        @param force: Force de-multiplexing with a Library Annotation sheet failing validation
        @type force: bool
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
            drms_list=drms_list,
            collection=collection,
            comparisons=comparisons,
            samples=samples)

        # Sub-class specific ...

        if library_path:
            self.library_path = library_path
        else:
            self.library_path = str()

        if sequences_directory:
            self.sequences_directory = sequences_directory
        else:
            self.sequences_directory = str()

        if samples_directory:
            self.samples_directory = samples_directory
        else:
            self.samples_directory = str()

        if experiment_directory:
            self.experiment_directory = experiment_directory
        else:
            self.experiment_directory = str()

        if classpath_illumina2bam:
            self.classpath_illumina2bam = classpath_illumina2bam
        else:
            self.classpath_illumina2bam = str()

        if classpath_picard:
            self.classpath_picard = classpath_picard
        else:
            self.classpath_picard = str()

        self.lanes = lanes

        self.force = force

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{BamIndexDecoder} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        super(BamIndexDecoder, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the library annotation file.

        if configuration.config_parser.has_option(section=section, option='library_path'):
            self.library_path = configuration.config_parser.get(section=section, option='library_path')

        # Get the BSF samples directory.

        if configuration.config_parser.has_option(section=section, option='sequences_directory'):
            self.sequences_directory = configuration.config_parser.get(
                section=section,
                option='sequences_directory')

        if configuration.config_parser.has_option(section=section, option='samples_directory'):
            self.samples_directory = configuration.config_parser.get(
                section=section,
                option='samples_directory')

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_illumina2bam'):
            self.classpath_illumina2bam = configuration.config_parser.get(
                section=section,
                option='classpath_illumina2bam')

        # Get the Picard tools Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_picard'):
            self.classpath_picard = configuration.config_parser.get(
                section=section,
                option='classpath_picard')

        if configuration.config_parser.has_option(section=section, option='lanes'):
            self.lanes = configuration.config_parser.getint(
                section=section,
                option='lanes')

        if configuration.config_parser.has_option(section=section, option='force'):
            self.force = configuration.config_parser.getboolean(
                section=section,
                option='force')

    def run(self):
        """Run the C{BamIndexDecoder} analysis to decode an archive BAM file produced with Illumina2Bam tools into
        sample-specific BAM files.
        """

        # The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        # a Picard tools *tab-separated* value (TSV) sample sheet.
        # lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        # barcode_sequence, barcode_name, library_name, sample_name, description

        super(BamIndexDecoder, self).run()

        # Read configuration options.

        # config_parser = self.configuration.config_parser
        # config_section = self.configuration.section_from_instance(self)

        default = Default.get_global_default()

        # Load from the configuration file and override with the default if necessary.

        # Define the sequences and samples directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.sequences_directory:
            self.sequences_directory = self.project_name

        self.sequences_directory = os.path.expanduser(path=self.sequences_directory)
        self.sequences_directory = os.path.expandvars(path=self.sequences_directory)

        if not os.path.isabs(self.sequences_directory):
            self.sequences_directory = os.path.join(Default.absolute_sequences(), self.sequences_directory)

        self.samples_directory = os.path.expanduser(path=self.samples_directory)
        self.samples_directory = os.path.expandvars(path=self.samples_directory)

        if not os.path.isabs(self.samples_directory):
            self.samples_directory = os.path.join(Default.absolute_samples(), self.samples_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(
                'The BamIndexDecoder samples_directory {!r} does not exist.'.format(self.samples_directory))

        self.experiment_directory = os.path.join(self.samples_directory, self.project_name)

        if not os.path.isdir(self.experiment_directory):
            try:
                os.makedirs(self.experiment_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        # Get the library annotation sheet.
        # The library annotation sheet is deliberately not passed in via sas_file,
        # as the Analysis.run() method reads that option into a BSF Collection object.

        self.library_path = os.path.expanduser(path=self.library_path)
        self.library_path = os.path.expandvars(path=self.library_path)

        if not self.library_path:
            self.library_path = string.join(words=(self.project_name, 'libraries.csv'), sep='_')

        if not os.path.exists(path=self.library_path):
            raise Exception('Library annotation file {!r} does not exist.'.format(self.library_path))

        # Load the library annotation sheet file and validate.

        library_annotation_sheet = LibraryAnnotationSheet.from_file_path(file_path=self.library_path)

        validation_messages = library_annotation_sheet.validate(lanes=self.lanes)

        if validation_messages:
            if self.force:
                warnings.warn('Validation of library annotation sheet {!r}:\n{}'.
                              format(self.library_path, validation_messages))
            else:
                raise Exception('Validation of library annotation sheet {!r}:\n{}'.
                                format(self.library_path, validation_messages))

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        if not self.classpath_illumina2bam:
            self.classpath_illumina2bam = default.classpath_illumina2bam

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        drms_bam_index_decoder = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_bam_index_decoder,
            working_directory=self.project_directory,
            analysis=self))

        index_by_lane = dict()

        for row_dict in library_annotation_sheet.row_dicts:
            if row_dict['lane'] in index_by_lane:
                lane_list = index_by_lane[row_dict['lane']]
            else:
                lane_list = list()
                index_by_lane[row_dict['lane']] = lane_list
            lane_list.append(row_dict)

        sample_annotation_sheet = SampleAnnotationSheet(
            file_path=os.path.join(
                self.experiment_directory,
                string.join(words=(self.project_name, 'samples.csv'), sep='_')))

        keys = index_by_lane.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            # The key represents the lane number as a Python str.

            file_path_dict = dict(
                samples_directory=os.path.join(
                    self.experiment_directory,
                    string.join(words=(self.project_name, key, 'samples'), sep='_')),
                barcode=os.path.join(
                    self.experiment_directory,
                    string.join(words=(self.project_name, key, 'barcode.tsv'), sep='_')),
                metrics=os.path.join(
                    self.experiment_directory,
                    string.join(words=(self.project_name, key, 'metrics.tsv'), sep='_')),
                input=os.path.join(
                    self.sequences_directory,
                    string.join(words=(self.project_name, key, 'sorted.bam'), sep='_'))
            )

            # Create the samples directory if it does not exist.

            if not os.path.isdir(file_path_dict['samples_directory']):
                try:
                    os.makedirs(file_path_dict['samples_directory'])
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            # Do not check whether the sorted BAM file exists, because at the time of
            # BamIndexDecoder submission the IlluminaToBam analysis may not have finished.
            #
            # if not os.path.exists(file_path_dict['input']):
            # raise Exception('Sequence archive BAM file {!r} does not exist.'.format(file_path_dict['input']))

            require_decoding = 0
            bam_index_decoder_sheet = BamIndexDecoderSheet(file_path=file_path_dict['barcode'])

            for row_dict in index_by_lane[key]:

                if len(row_dict['barcode_sequence_1']) or len(row_dict['barcode_sequence_2']):
                    require_decoding = 1

                # Add a row to the lane-specific tab-delimited IlluminaToBamTools BamIndexDecoder barcode file.

                bam_index_decoder_sheet.row_dicts.append(dict(
                    barcode_sequence=row_dict['barcode_sequence_1'] + row_dict['barcode_sequence_2'],
                    barcode_name=row_dict['sample_name'],
                    library_name=row_dict['library_name'],
                    sample_name=row_dict['sample_name'],
                    description=str()
                ))

                # Write the flow-cell-specific sample annotation sheet.
                sample_dict = dict(
                    ProcessedRunFolder=self.project_name,
                    Project=row_dict['library_name'],
                    Sample=row_dict['sample_name'],
                    Reads1=string.join(words=(self.project_name, key, row_dict['sample_name']), sep='_'),
                    File1=os.path.join(
                        file_path_dict['samples_directory'],
                        '{}_{}#{}.bam'.format(self.project_name, key, row_dict['sample_name'])),
                    LibrarySize=row_dict['library_size'],
                    Barcode1=row_dict['barcode_sequence_1'],
                    Barcode2=row_dict['barcode_sequence_2'])

                sample_annotation_sheet.row_dicts.append(sample_dict)

            # Write the lane-specific BamIndexDecoderSheet to the internal file path.
            bam_index_decoder_sheet.write_to_file()

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            runnable_bam_index_decoder = self.add_runnable(runnable=Runnable(
                name=self.get_prefix_bam_index_decoder(project_name=self.project_name, lane=key),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict))

            # TODO: It would be good to extend the Runnable so that it holds dependencies on other Runnable objects
            # and that it could be submitted to a DRMS so that the Executable gets automatically created and submitted.
            executable_bam_index_decoder = drms_bam_index_decoder.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_bam_index_decoder.name))

            executable_bam_index_decoder.dependencies.append(
                IlluminaToBam.get_prefix_illumina_to_bam(project_name=self.project_name, lane=key))

            # Only submit this Executable if the final result file does not exist.
            if (os.path.exists(
                    os.path.join(self.project_directory, file_path_dict['metrics']))
                and os.path.getsize(
                    os.path.join(self.project_directory, file_path_dict['metrics']))):
                executable_bam_index_decoder.submit = False

            if require_decoding:

                # Run the BamIndexDecoder if there is at least one line containing a barcode sequence.

                java_process = runnable_bam_index_decoder.add_runnable_step(runnable_step=RunnableStep(
                    name='bam_index_decoder',
                    program='java',
                    sub_command=Command(command=str()),
                    obsolete_file_path_list=[
                        file_path_dict['barcode']
                    ]))

                java_process.add_switch_short(key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_illumina2bam, 'BamIndexDecoder.jar'))
                java_process.add_switch_short(key='Xmx4G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=runnable_bam_index_decoder.get_relative_temporary_directory_path)

                sub_command = java_process.sub_command

                sub_command.add_option_pair(
                    key='INPUT',
                    value=file_path_dict['input'])
                # OUTPUT
                sub_command.add_option_pair(
                    key='OUTPUT_DIR',
                    value=file_path_dict['samples_directory'])
                sub_command.add_option_pair(
                    key='OUTPUT_PREFIX',
                    value=string.join(words=(self.project_name, key), sep='_'))
                sub_command.add_option_pair(
                    key='OUTPUT_FORMAT',
                    value='bam')
                # BARCODE_TAG_NAME
                # BARCODE_QUALITY_TAG_NAME
                # BARCODE
                sub_command.add_option_pair(
                    key='BARCODE_FILE',
                    value=file_path_dict['barcode'])
                sub_command.add_option_pair(
                    key='METRICS_FILE',
                    value=file_path_dict['metrics'])
                # MAX_MISMATCHES
                # MIN_MISMATCH_DELTA
                # MAX_NO_CALLS
                # CONVERT_LOW_QUALITY_TO_NO_CALL
                # MAX_LOW_QUALITY_TO_CONVERT
                sub_command.add_option_pair(
                    key='TMP_DIR',
                    value=runnable_bam_index_decoder.get_relative_temporary_directory_path)
                sub_command.add_option_pair(
                    key='VERBOSITY',
                    value='WARNING')
                # QUIET
                # VALIDATION_STRINGENCY
                # COMPRESSION_LEVEL
                # MAX_RECORDS_IN_RAM
                sub_command.add_option_pair(
                    key='CREATE_INDEX',
                    value='false')
                sub_command.add_option_pair(
                    key='CREATE_MD5_FILE',
                    value='true')
                # OPTIONS_FILE

            else:

                # Run Picard CollectAlignmentSummaryMetrics if there is no line containing a barcode sequence.

                java_process = runnable_bam_index_decoder.add_runnable_step(runnable_step=RunnableStep(
                    name='picard_collect_alignment_summary_metrics',
                    program='java',
                    sub_command=Command(command=str())))

                java_process.add_switch_short(key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_picard, 'CollectAlignmentSummaryMetrics.jar'))
                java_process.add_switch_short(key='Xmx4G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=runnable_bam_index_decoder.get_relative_temporary_directory_path)

                sub_command = java_process.sub_command

                # MAX_INSERT_SIZE
                # ADAPTER_SEQUENCE
                sub_command.add_option_pair(
                    key='METRIC_ACCUMULATION_LEVEL',
                    value='READ_GROUP')
                # IS_BISULFITE_SEQUENCED
                sub_command.add_option_pair(
                    key='INPUT',
                    value=file_path_dict['input'])
                sub_command.add_option_pair(
                    key='OUTPUT',
                    value=file_path_dict['metrics'])
                # REFERENCE_SEQUENCE
                # ASSUME_SORTED
                # STOP_AFTER
                sub_command.add_option_pair(
                    key='TMP_DIR',
                    value=runnable_bam_index_decoder.get_relative_temporary_directory_path)
                sub_command.add_option_pair(
                    key='VERBOSITY',
                    value='WARNING')
                sub_command.add_option_pair(
                    key='QUIET',
                    value='false')
                sub_command.add_option_pair(
                    key='VALIDATION_STRINGENCY',
                    value='STRICT')
                sub_command.add_option_pair(
                    key='COMPRESSION_LEVEL',
                    value='5')
                sub_command.add_option_pair(
                    key='MAX_RECORDS_IN_RAM',
                    value='4000000')
                sub_command.add_option_pair(
                    key='CREATE_INDEX',
                    value='true')
                sub_command.add_option_pair(
                    key='CREATE_MD5_FILE',
                    value='true')
                # OPTIONS_FILE

                # TODO: It would be even better to run Picard AddOrReplaceReadGroups to get a correct SAM header.
                # Add a symbolic link to the BSF Sequence Archive file within the samples directory.
                file_path_dict['link'] = os.path.join(
                    file_path_dict['samples_directory'],
                    '{}_{}#{}.bam'.format(self.project_name, key, index_by_lane[key][0]['sample_name']))

        # Finally, write the flow-cell-specific SampleAnnotationSheet to the internal file path.

        sample_annotation_sheet.write_to_file()
