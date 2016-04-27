"""bsf.analyses.illumina_to_bam_tools

A package of classes and methods supporting analyses of the Illumina2Bam-Tools package.
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


import errno
import os
import warnings

from bsf import Analysis, DRMS, Runnable
from bsf.analyses.illumina_run_folder import IlluminaRunFolderRestore
from bsf.annotation import BamIndexDecoderSheet, LibraryAnnotationSheet, SampleAnnotationSheet
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.process import Command, Executable, RunnableStep, RunnableStepJava, RunnableStepPicard, RunnableStepLink, \
    RunnableStepMakeDirectory, RunnableStepMove
from bsf.standards import Configuration, Default


class RunnableStepIlluminaToBam(RunnableStepJava):
    """The C{RunnableStepIlluminaToBam} class represents a C{RunnableStepJava} specific to IlluminaToBam tools.
    IlluminaToBam tools use the old Picard tools interface where each algorithm is implemented as a separate
    Java Archive (JAR) file

    Attributes:
    @ivar java_jar_path: Java archive file path
    @type java_jar_path: str | unicode
    @ivar itb_classpath: IlluminaToBam class path
    @type itb_classpath: str | unicode
    @ivar itb_command: IlluminaToBam command
    @type itb_command: str
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_maximum=None,
            java_jar_path=None,
            itb_classpath=None,
            itb_command=None):
        """Create a C{RunnableStep} for an IlluminaToBam algorithm.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @param itb_classpath: IlluminaToBam class path
        @type itb_classpath: str | unicode
        @param itb_command: IlluminaToBam command
        @type itb_command: str
        @return: C{RunnableStep}
        @rtype: RunnableStep
        """

        super(RunnableStepIlluminaToBam, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path,
            java_heap_maximum=java_heap_maximum,
            java_jar_path=java_jar_path)

        # Set the IlluminaToBam classpath and the IlluminaToBam Java archive.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(itb_classpath, itb_command + '.jar'))

        # Introduce another empty sub-command to separate Illumina2Bam-tools options.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command()

        return

    def add_itb_option(self, key, value, override=False):
        """Add an IlluminaToBam option.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


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
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
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
    @ivar vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
    @type vendor_quality_filter: dict[str, bool]
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
        return '_'.join((cls.drms_name_illumina_to_bam, project_name, lane))

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
            drms_list=None,
            collection=None,
            comparisons=None,
            samples=None,
            run_directory=None,
            intensity_directory=None,
            basecalls_directory=None,
            experiment_name=None,
            sequencing_centre=None,
            sequences_directory=None,
            experiment_directory=None,
            classpath_illumina2bam=None,
            classpath_picard=None,
            vendor_quality_filter=None,
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
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
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
        @param vendor_quality_filter: Python C{dict} of flow cell chemistry type and Python bool value for filtering
        @type vendor_quality_filter: dict[str, bool]
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
        @return:
        @rtype:
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

        if sequencing_centre is None:
            self.sequencing_centre = str()
        else:
            self.sequencing_centre = sequencing_centre

        if sequences_directory is None:
            self.sequences_directory = str()
        else:
            self.sequences_directory = sequences_directory

        if experiment_directory is None:
            self.experiment_directory = str()
        else:
            self.experiment_directory = experiment_directory

        if classpath_illumina2bam is None:
            self.classpath_illumina2bam = str()
        else:
            self.classpath_illumina2bam = classpath_illumina2bam

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if vendor_quality_filter is None:
            self.vendor_quality_filter = dict()
        else:
            self.vendor_quality_filter = vendor_quality_filter

        if force is None:
            self.force = False
        else:
            assert isinstance(force, bool)
            self.force = force

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{IlluminaToBam} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        super(IlluminaToBam, self).set_configuration(configuration=configuration, section=section)

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

        # Get sequencing centre information.

        option = 'sequencing_centre'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequencing_centre = configuration.config_parser.get(section=section, option=option)

        # Get the sequences directory information.

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        option = 'classpath_illumina2bam'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_illumina2bam = configuration.config_parser.get(section=section, option=option)

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

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
        """Run this C{IlluminaToBam} C{Analysis}.

        Convert an Illumina flow cell into lane-specific archive BAM files.

        To convert an Illumina flow cell, Illumina2bam is run first, setting the SAM Read Group (@RG)
        library name (LB) and sample name (SM) to 'flow cell identifier.lane'.
        The resulting archive BAM file is then sorted by query name with Picard SortSam.
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
            intensity_directory = Default.get_absolute_path(
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
            basecalls_directory = Default.get_absolute_path(
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
            self.sequencing_centre = default.operator_sequencing_centre

        # Define the sequences directory in which to create the experiment directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # An absolute path cannot be prepended.

        if self.sequences_directory:
            self.sequences_directory = Default.get_absolute_path(
                file_path=self.sequences_directory)
        else:
            self.sequences_directory = Default.absolute_sequences()

        # As a safety measure, to prevent creation of rogue directory paths, the sequences_directory has to exist.

        if not os.path.isdir(self.sequences_directory):
            raise Exception(
                'The IlluminaToBam sequences_directory {!r} does not exist.'.format(self.sequences_directory))

        self.experiment_directory = os.path.join(self.sequences_directory, self.project_name)

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        if not self.classpath_illumina2bam:
            self.classpath_illumina2bam = default.classpath_illumina2bam

        # Get the Picard tools Java Archive (JAR) class path directory

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        # Check that the flow cell chemistry type is defined in the vendor quality filter.

        if irf.run_parameters.get_flow_cell_type not in self.vendor_quality_filter:
            raise Exception('Flow cell chemistry type {!r} not defined.'.format(irf.run_parameters.get_flow_cell_type))

        # Call the run method of the super class after the project_name has been defined.

        super(IlluminaToBam, self).run()

        drms_illumina_to_bam = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_illumina_to_bam,
                working_directory=self.project_directory,
                analysis=self))

        for lane in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):

            lane_str = str(lane)

            file_path_dict = dict(
                illumina_directory=self.run_directory,  # contains full path information
                sequences_directory=self.sequences_directory,  # contains full path information
                experiment_directory=self.experiment_directory,  # contains full path information
                unsorted_bam='_'.join((self.project_name, lane_str, 'unsorted.bam')),
                unsorted_md5='_'.join((self.project_name, lane_str, 'unsorted.bam.md5')),
                sorted_bam='_'.join((self.project_name, lane_str, 'sorted.bam')),
                sorted_md5='_'.join((self.project_name, lane_str, 'sorted.bam.md5')),
                lane_bam='{}_{:d}.bam'.format(self.project_name, lane),
                lane_md5='{}_{:d}.bam.md5'.format(self.project_name, lane),
            )

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            runnable_illumina_to_bam = self.add_runnable(
                runnable=Runnable(
                    name=self.get_prefix_illumina_to_bam(project_name=self.project_name, lane=lane_str),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory,
                    file_path_dict=file_path_dict))

            # TODO: The Runnable class could have dependencies just like the Executable class so that they could be
            # passed on upon creation of the Executable from the Runnable via Executable.from_analysis_runnable().
            executable_illumina_to_bam = self.set_drms_runnable(
                drms=drms_illumina_to_bam,
                runnable=runnable_illumina_to_bam)
            executable_illumina_to_bam.dependencies.append(
                IlluminaRunFolderRestore.get_prefix_compress_base_calls(
                    project_name=self.project_name,
                    lane=lane_str))

            # Only submit this Executable if the final result file does not exist.
            if os.path.exists(file_path_dict['sorted_md5']) and os.path.getsize(file_path_dict['sorted_md5']):
                executable_illumina_to_bam.submit = False

            # Run Illumina2Bam tools Illumina2bam.

            runnable_step = runnable_illumina_to_bam.add_runnable_step(
                runnable_step=RunnableStepIlluminaToBam(
                    name='illumina_to_bam',
                    java_temporary_path=runnable_illumina_to_bam.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    itb_classpath=self.classpath_illumina2bam,
                    itb_command='Illumina2bam'))
            assert isinstance(runnable_step, RunnableStepIlluminaToBam)
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
            runnable_step.add_itb_option(key='OUTPUT', value=file_path_dict['unsorted_bam'])
            # GENERATE_SECONDARY_BASE_CALLS defaults to 'false'.
            # PF_FILTER defaults to 'true'.
            if not self.vendor_quality_filter[irf.run_parameters.get_flow_cell_type]:
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
            # TMP_DIR
            runnable_step.add_itb_option(
                key='TMP_DIR',
                value=runnable_illumina_to_bam.get_relative_temporary_directory_path)
            # VERBOSITY defaults to 'INFO'.
            runnable_step.add_itb_option(key='VERBOSITY', value='WARNING')
            # QUIET defaults to 'false'.
            # VALIDATION_STRINGENCY defaults to 'STRICT'.
            # COMPRESSION_LEVEL defaults to '5'.
            runnable_step.add_itb_option(key='COMPRESSION_LEVEL', value='9')
            # MAX_RECORDS_IN_RAM defaults to '500000'.
            # CREATE_INDEX defaults to 'false'.
            # CREATE_MD5_FILE defaults to 'false'.
            # OPTIONS_FILE

            # Run Picard SortSam

            runnable_step = runnable_illumina_to_bam.add_runnable_step(
                runnable_step=RunnableStepPicard(
                    name='picard_sort_sam',
                    obsolete_file_path_list=[
                        file_path_dict['unsorted_bam'],
                        file_path_dict['unsorted_md5'],
                    ],
                    java_temporary_path=runnable_illumina_to_bam.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    picard_classpath=self.classpath_picard,
                    picard_command='SortSam'))
            assert isinstance(runnable_step, RunnableStepPicard)
            # INPUT is required.
            runnable_step.add_picard_option(key='INPUT', value=file_path_dict['unsorted_bam'])
            # OUTPUT is required.
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict['sorted_bam'])
            # SORT_ORDER is required.
            runnable_step.add_picard_option(key='SORT_ORDER', value='queryname')
            # TMP_DIR
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_illumina_to_bam.get_relative_temporary_directory_path)
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

            runnable_illumina_to_bam.add_runnable_step(
                runnable_step=RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=self.experiment_directory))

            # Move and rename the final, sorted BAM file.

            runnable_illumina_to_bam.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='move_sorted_bam',
                    source_path=file_path_dict['sorted_bam'],
                    target_path=os.path.join(self.experiment_directory, file_path_dict['lane_bam'])))

            # Move and rename the checksum file.

            runnable_illumina_to_bam.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='move_sorted_md5',
                    source_path=file_path_dict['sorted_md5'],
                    target_path=os.path.join(self.experiment_directory, file_path_dict['lane_md5'])))

        return


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
        return '_'.join((cls.drms_name_bam_index_decoder, project_name, lane))

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
            drms_list=None,
            collection=None,
            comparisons=None,
            samples=None,
            library_path=None,
            sequences_directory=None,
            samples_directory=None,
            experiment_directory=None,
            classpath_illumina2bam=None,
            classpath_picard=None,
            lanes=8,
            force=False):
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
        @return:
        @rtype:
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

        if library_path is None:
            self.library_path = str()
        else:
            self.library_path = library_path

        if sequences_directory is None:
            self.sequences_directory = str()
        else:
            self.sequences_directory = sequences_directory

        if samples_directory is None:
            self.samples_directory = str()
        else:
            self.samples_directory = samples_directory

        if experiment_directory is None:
            self.experiment_directory = str()
        else:
            self.experiment_directory = experiment_directory

        if classpath_illumina2bam is None:
            self.classpath_illumina2bam = str()
        else:
            self.classpath_illumina2bam = classpath_illumina2bam

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if lanes is None:
            self.lanes = int()
        else:
            assert isinstance(lanes, int)
            self.lanes = lanes

        if force is None:
            self.force = False
        else:
            assert isinstance(force, bool)
            self.force = force

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{BamIndexDecoder} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(BamIndexDecoder, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the library annotation file.

        option = 'library_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_path = configuration.config_parser.get(section=section, option=option)

        # Get the BSF samples directory.

        option = 'sequences_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sequences_directory = configuration.config_parser.get(section=section, option=option)

        option = 'samples_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.samples_directory = configuration.config_parser.get(section=section, option=option)

        # Get the Illumina2Bam tools Java Archive (JAR) class path directory.

        option = 'classpath_illumina2bam'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_illumina2bam = configuration.config_parser.get(section=section, option=option)

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
        """Run the C{BamIndexDecoder} analysis to decode an archive BAM file produced with Illumina2Bam tools into
        sample-specific BAM files.
        @return:
        @rtype:
        """

        # The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        # a Picard tools *tab-separated* value (TSV) sample sheet.
        # lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        # barcode_sequence, barcode_name, library_name, sample_name, description

        super(BamIndexDecoder, self).run()

        default = Default.get_global_default()

        # Load from the configuration file and override with the default if necessary.

        # Define the sequences and samples directory.
        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.sequences_directory:
            self.sequences_directory = self.project_name

        self.sequences_directory = Default.get_absolute_path(
            file_path=self.sequences_directory,
            default_path=Default.absolute_sequences())

        self.samples_directory = Default.get_absolute_path(
            file_path=self.samples_directory,
            default_path=Default.absolute_samples())

        # As a safety measure, to prevent creation of rogue directory paths, the samples_directory has to exist.

        if not os.path.isdir(self.samples_directory):
            raise Exception(
                'The BamIndexDecoder samples_directory {!r} does not exist.'.format(self.samples_directory))

        self.experiment_directory = os.path.join(self.samples_directory, self.project_name)

        # Get the library annotation sheet.
        # The library annotation sheet is deliberately not passed in via sas_file,
        # as the Analysis.run() method reads that option into a BSF Collection object.

        self.library_path = os.path.expanduser(path=self.library_path)
        self.library_path = os.path.expandvars(path=self.library_path)

        if not self.library_path:
            self.library_path = '_'.join((self.project_name, 'libraries.csv'))

        self.library_path = os.path.normpath(path=self.library_path)

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

        drms_bam_index_decoder = self.add_drms(
            drms=DRMS.from_analysis(
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

        # Create a Sample Annotation Sheet in the project directory and
        # eventually transfer it into the experiment_directory.
        sample_annotation_name = '_'.join((self.project_name, 'samples.csv'))
        sample_annotation_sheet = SampleAnnotationSheet(
            file_path=os.path.join(self.project_directory, sample_annotation_name))
        sample_annotation_transferred = 0

        keys = index_by_lane.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            # The key represents the lane number as a Python str.

            file_path_dict = dict(
                project_barcode='_'.join((self.project_name, key, 'barcode.tsv')),
                samples_directory='_'.join((self.project_name, key, 'samples')),
                barcode_tsv='_'.join((self.project_name, key, 'barcode.tsv')),
                metrics_tsv='_'.join((self.project_name, key, 'metrics.tsv')),
                input_bam=os.path.join(
                    self.sequences_directory,
                    '{}_{}.bam'.format(self.project_name, key)),
            )

            # Do not check whether the sorted input BAM file exists, because at the time of
            # BamIndexDecoder submission the IlluminaToBam analysis may not have finished.
            #
            # if not os.path.exists(file_path_dict['input_bam']):
            # raise Exception('Sequence archive BAM file {!r} does not exist.'.format(file_path_dict['input_bam']))

            barcode_number = 0
            bam_index_decoder_sheet = BamIndexDecoderSheet(
                file_path=os.path.join(self.project_directory, file_path_dict['barcode_tsv']))

            for row_dict in index_by_lane[key]:

                # Determine the number of barcodes on the basis of the first line.
                # All other lines of a lane have to use the same number and lengths of barcodes.

                if not barcode_number:
                    if len(row_dict['barcode_sequence_1']):
                        barcode_number += 1
                    if len(row_dict['barcode_sequence_2']):
                        barcode_number += 1

                # Add a row to the lane-specific tab-delimited IlluminaToBamTools BamIndexDecoder barcode file.

                bam_index_decoder_sheet.row_dicts.append(dict(
                    barcode_sequence=row_dict['barcode_sequence_1'] + row_dict['barcode_sequence_2'],
                    barcode_name=row_dict['sample_name'],
                    library_name=row_dict['library_name'],
                    sample_name=row_dict['sample_name'],
                    description=str()
                ))

                # Add a row to the flow cell-specific sample annotation sheet.

                sample_dict = dict({
                    'File Type': 'Automatic',
                    'ProcessedRunFolder Name': self.project_name,
                    'Project Name': row_dict['library_name'],
                    'Project Size': row_dict['library_size'],
                    'Sample Name': row_dict['sample_name'],
                    'PairedReads Barcode1': row_dict['barcode_sequence_1'],
                    'PairedReads Barcode2': row_dict['barcode_sequence_2'],
                    # TODO: It would be good to add a RunnableStep to populate the ReadGroup.
                    'PairedReads ReadGroup': '',
                    'Reads1 Name': '_'.join((self.project_name, key, row_dict['sample_name'])),
                    'Reads1 File': os.path.join(
                        os.path.basename(self.experiment_directory),
                        file_path_dict['samples_directory'],
                        '{}_{}#{}.bam'.format(self.project_name, key, row_dict['sample_name'])),
                    'Reads2 Name': '',
                    'Reads2 File': '',
                })

                sample_annotation_sheet.row_dicts.append(sample_dict)

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            runnable_bam_index_decoder = self.add_runnable(
                runnable=Runnable(
                    name=self.get_prefix_bam_index_decoder(project_name=self.project_name, lane=key),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory,
                    file_path_dict=file_path_dict))

            # TODO: It would be good to extend the Runnable so that it holds dependencies on other Runnable objects
            # and that it could be submitted to a DRMS so that the Executable gets automatically created and submitted.
            executable_bam_index_decoder = self.set_drms_runnable(
                drms=drms_bam_index_decoder,
                runnable=runnable_bam_index_decoder)
            executable_bam_index_decoder.dependencies.append(
                IlluminaToBam.get_prefix_illumina_to_bam(
                    project_name=self.project_name,
                    lane=key))

            if executable_bam_index_decoder.submit:
                # Only if this Executable actually gets submitted ...
                # Create the samples directory in the project_directory if it does not exist.

                project_samples_path = os.path.join(self.project_directory, file_path_dict['samples_directory'])
                if not os.path.isdir(project_samples_path):
                    try:
                        os.makedirs(project_samples_path)
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise

                # Write the lane-specific BamIndexDecoderSheet to the internal file path.

                bam_index_decoder_sheet.write_to_file()

            if barcode_number:

                # Run the BamIndexDecoder if there is at least one line containing a barcode sequence.

                runnable_step = runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepIlluminaToBam(
                        name='illumina_to_bam',
                        java_temporary_path=runnable_bam_index_decoder.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        itb_classpath=self.classpath_illumina2bam,
                        itb_command='BamIndexDecoder'))
                assert isinstance(runnable_step, RunnableStepIlluminaToBam)
                # INPUT is required
                runnable_step.add_itb_option(key='INPUT', value=file_path_dict['input_bam'])
                # OUTPUT is required, but cannot be used together with OUTPUT_FORMAT, OUTPUT_PREFIX and OUTPUT_DIR.
                # OUTPUT_DIR is required, but cannot be used with OUTPUT.
                runnable_step.add_itb_option(key='OUTPUT_DIR', value=file_path_dict['samples_directory'])
                # OUTPUT_PREFIX is required, but cannot be used with OUTPUT.
                runnable_step.add_itb_option(key='OUTPUT_PREFIX', value='_'.join((self.project_name, key)))
                # OUTPUT_FORMAT is required, but cannot be used with OUTPUT.
                runnable_step.add_itb_option(key='OUTPUT_FORMAT', value='bam')
                # BARCODE_TAG_NAME defaults to 'BC'.
                # BARCODE_QUALITY_TAG_NAME defaults to 'QT'.
                # BARCODE cannot be used with BARCODE_FILE.
                # BARCODE_FILE is required, but cannot be used with BARCODE.
                runnable_step.add_itb_option(key='BARCODE_FILE', value=file_path_dict['barcode_tsv'])
                # METRICS_FILE is required.
                runnable_step.add_itb_option(key='METRICS_FILE', value=file_path_dict['metrics_tsv'])
                # MAX_MISMATCHES defaults to '1'.
                # if len(index_by_lane[key][0]['barcode_mismatches']):
                #     runnable_step.add_itb_option(
                #         key='MAX_MISMATCHES',
                #         value=index_by_lane[key][0]['barcode_mismatches'])
                # elif barcode_number == 2:
                if barcode_number == 2:
                    runnable_step.add_itb_option(key='MAX_MISMATCHES', value='2')
                # MIN_MISMATCH_DELTA defaults to '1'.
                # MAX_NO_CALLS defaults to '2'.
                # CONVERT_LOW_QUALITY_TO_NO_CALL defaults to 'false'.
                # MAX_LOW_QUALITY_TO_CONVERT defaults to '15'.
                # TMP_DIR
                runnable_step.add_itb_option(
                    key='TMP_DIR',
                    value=runnable_bam_index_decoder.get_relative_temporary_directory_path)
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

                # Create the experiment directory if it does not exist already.

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=self.experiment_directory))

                # Move the samples directory into the experiment directory.

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_samples_directory',
                        source_path=file_path_dict['samples_directory'],
                        target_path=self.experiment_directory))

                # Move the metrics file into the experiment directory.

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_metrics_tsv',
                        source_path=file_path_dict['metrics_tsv'],
                        target_path=self.experiment_directory))

            else:

                # Run Picard CollectAlignmentSummaryMetrics if there is no line containing a barcode sequence.

                runnable_step = runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='picard_collect_alignment_summary_metrics',
                        obsolete_file_path_list=[
                            file_path_dict['barcode_tsv'],
                        ],
                        java_temporary_path=runnable_bam_index_decoder.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        picard_classpath=self.classpath_picard,
                        picard_command='CollectAlignmentSummaryMetrics'))
                assert isinstance(runnable_step, RunnableStepPicard)
                # MAX_INSERT_SIZE defaults to '100000'.
                # ADAPTER_SEQUENCE
                # METRIC_ACCUMULATION_LEVEL.
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP')
                # IS_BISULFITE_SEQUENCED defaults to 'false'.
                # INPUT is required.
                runnable_step.add_picard_option(key='INPUT', value=file_path_dict['input_bam'])
                # OUTPUT is required.
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict['metrics_tsv'])
                # REFERENCE_SEQUENCE defaults to 'null'.
                # ASSUME_SORTED defaults to 'true'.
                # STOP_AFTER defaults to '0'.
                # TMP_DIR
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_bam_index_decoder.get_relative_temporary_directory_path)
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

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=self.experiment_directory))

                # Move the samples directory into the experiment directory.

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_samples_directory',
                        source_path=file_path_dict['samples_directory'],
                        target_path=self.experiment_directory))

                # Add a symbolic link to the BSF Sequence Archive file within the samples directory.

                experiment_samples_directory = os.path.join(
                    self.experiment_directory,
                    file_path_dict['samples_directory'])

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='link',
                        source_path=os.path.relpath(file_path_dict['input_bam'], experiment_samples_directory),
                        target_path=os.path.join(
                            experiment_samples_directory,
                            '{}_{}#{}.bam'.format(self.project_name, key,
                                                  index_by_lane[key][0]['sample_name']))))

                # Move the metrics file into the experiment directory.

                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_metrics_tsv',
                        source_path=file_path_dict['metrics_tsv'],
                        target_path=self.experiment_directory))

            # Move the Sample Annotation Sheet once.

            if not sample_annotation_transferred:
                sample_annotation_transferred += 1
                runnable_bam_index_decoder.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_sample_annotation',
                        source_path=os.path.join(self.project_directory, sample_annotation_name),
                        target_path=self.experiment_directory))

        # Finally, write the flow cell-specific SampleAnnotationSheet to the internal file path.

        sample_annotation_sheet.write_to_file()

        return
