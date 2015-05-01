"""bsf.analyses.illumina_run_folder

A package of classes and methods supporting analyses to archive and restore Illumina Run Folders.
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


import os
import string

from bsf import Analysis, Command, Configuration, Default, DRMS, Executable, Runnable, RunnableStep
from bsf.illumina import RunFolder, RunFolderNotComplete


class IlluminaRunFolderArchive(Analysis):
    """The C{IlluminaRunFolderArchive} class represents the logic to archive an Illumina Run Folder in a format
    suitable for magnetic tape libraries.

    Attributes:
    @cvar drms_name_pre_process_folder: C{DRMS.name} for the C{IlluminaRunFolderArchive} C{Analysis} stage
    @type drms_name_pre_process_folder: str
    @cvar drms_name_post_process_folder: C{DRMS.name} for the C{IlluminaRunFolderArchive} C{Analysis} stage
    @type drms_name_post_process_folder: str
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | unicode
    @ivar experiment_name: Experiment name (i.e. flow-cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool
    """

    drms_name_pre_process_folder = 'irf_archive_pre_process_folder'
    drms_name_post_process_folder = 'irf_archive_post_process_folder'

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{IlluminaRunFolderArchive} object from a UNIX-style configuration file via the
        C{Configuration} class.

        @param config_path: UNIX-style configuration file
        @type config_path: str | unicode
        @return: C{IlluminaRunFolderArchive}
        @rtype: IlluminaRunFolderArchive
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{IlluminaRunFolderArchive} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{IlluminaRunFolderArchive}
        @rtype: IlluminaRunFolderArchive
        """

        assert isinstance(configuration, Configuration)

        irf_archive = cls(configuration=configuration)

        # A "bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive" section specifies defaults
        # for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        irf_archive.set_configuration(irf_archive.configuration, section=section)

        return irf_archive

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 run_directory=None, experiment_name=None,
                 force=False):
        """Initialise a C{IlluminaRunFolderArchive} object.

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
        @param experiment_name: Experiment name (i.e. flow-cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
        """

        super(IlluminaRunFolderArchive, self).__init__(
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

        if experiment_name:
            self.experiment_name = experiment_name
        else:
            self.experiment_name = str()

        self.force = force

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{IlluminaRunFolderArchive} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        super(IlluminaRunFolderArchive, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get Illumina Run Folder information.

        if configuration.config_parser.has_option(section=section, option='illumina_run_folder'):
            self.run_directory = configuration.config_parser.get(
                section=section,
                option='illumina_run_folder')

        # Get the experiment name.

        if configuration.config_parser.has_option(section=section, option='experiment_name'):
            self.experiment_name = configuration.config_parser.get(
                section=section,
                option='experiment_name')

        if configuration.config_parser.has_option(section=section, option='force'):
            self.force = configuration.config_parser.getboolean(
                section=section,
                option='force')

    def run(self):
        """Run this C{IlluminaRunFolderArchive} C{Analysis}.

        Archive an Illumina Run Folder in a format suitable for magnetic tape libraries.

        1. Check if the Illumina Run has finished by testing for an
            RTAComplete.txt file.
        2. Check if an archive process is already running by testing for an
            archive directory.
        3. Create an archive directory.
        4. Size the native Illumina Run Folder via the du utility.
        5. Reset the file permissions for all directories via the find utility.
        6. Reset the file permissions for all regular files via the find utility.
        X. Zip all Data/Intensities/BaseCalls/L00[1-8]/CX.1/*.bcl files.
            find CX.1 -name '*.bcl' -execdir gzip --best --verbose {} \+
        X. Zip all Data/RTALogs/*.txt files
            gzip --best --recursive Data/RTALogs/
        7. Run the GNU tar utility over each Data/Intensities/L00[1-8] directory,
            before deleting the directory.
        8. Run the GNU tar utility over the remaining Data/Intensities directory,
            before deleting the Intensities directory.
        9. Run the GNU tar utility over the remaining Illumina Run folder.
        10. Record the archive file sizes via the ls utility.
        """

        # default = Default.get_global_default()

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

        # Check whether the Illumina Run Folder is complete.
        # 1. Check whether the RTAComplete.txt file exists in the Illumina Run Folder
        # to prevent archiving and deleting of an active folder.
        # Otherwise, require force to start archiving.

        if not os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) and not self.force:
            raise RunFolderNotComplete(
                'The Illumina run directory {!r} is not complete.'.format(self.run_directory))

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

        super(IlluminaRunFolderArchive, self).run()

        drms_pre_process_folder = DRMS.from_analysis(
            name=self.drms_name_pre_process_folder,
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(drms_pre_process_folder)

        drms_post_process_folder = DRMS.from_analysis(
            name=self.drms_name_post_process_folder,
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(drms_post_process_folder)

        # Pre-process on folder level.

        pre_process_folder_prefix = string.join(words=(drms_pre_process_folder.name, self.project_name), sep='_')

        file_path_dict = dict(
            temporary_directory=string.join((pre_process_folder_prefix, 'temporary'), sep='_'),
        )

        # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
        runnable_pre_process_folder = Runnable(
            name=pre_process_folder_prefix,
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict)
        self.add_runnable(runnable=runnable_pre_process_folder)

        # Create an Executable for pre-processing the folder.

        executable_pre_process_folder = Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_pre_process_folder.name)
        drms_pre_process_folder.add_executable(executable=executable_pre_process_folder)

        # TODO: The pre-processing does not depend on any other process, does it?
        # executable_pre_process_folder.dependencies.extend(vc_process_cohort_dependencies)

        # TODO:
        # 0. Check whether Picard ExtractIlluminaBarcodes has written any
        # s_<lane>_<tile>_barcode.txt(.gz) files into the BaseCalls directory.
        # Keeping them is rather pointless and they should be removed.
        # http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes

        # 5. Reset all file permissions for directories.

        reset_directory_permissions = RunnableStep(
            name='reset_directory_permissions',
            program='find',
            sub_command=Command(command=self.run_directory))
        runnable_pre_process_folder.add_runnable_step(runnable_step=reset_directory_permissions)

        reset_directory_permissions.add_option_short(key='type', value='d')
        reset_directory_permissions.add_option_short(key='execdir', value='chmod u=rwx,go=rx {} +')

        # 6. Reset all file permissions for regular files.

        reset_file_permissions = RunnableStep(
            name='reset_file_permissions',
            program='find',
            sub_command=Command(command=self.run_directory))
        runnable_pre_process_folder.add_runnable_step(runnable_step=reset_file_permissions)

        reset_file_permissions.add_option_short(key='type', value='f')
        reset_file_permissions.add_option_short(key='execdir', value='chmod u=rw,go=r {} +')

        # 7. Compress all files in the Logs and Logs/IALogs directories.
        # TODO: Check, whether the Illumina SAV needs them.

        compress_logs = RunnableStep(
            name='compress_logs',
            program='gzip')
        runnable_pre_process_folder.add_runnable_step(runnable_step=compress_logs)

        compress_logs.add_switch_long(key='best')
        compress_logs.add_switch_long(key='recursive')
        compress_logs.arguments.append(os.path.join(self.run_directory, 'Logs'))

        # 8. Compress all files in the Data/RTALogs directory.
        # TODO: Check, whether the Illumina SAV needs them.

        compress_rta_logs = RunnableStep(
            name='compress_rta_logs',
            program='gzip')
        runnable_pre_process_folder.add_runnable_step(runnable_step=compress_rta_logs)

        compress_rta_logs.add_switch_long(key='best')
        compress_rta_logs.add_switch_long(key='recursive')
        compress_rta_logs.arguments.append(os.path.join(self.run_directory, 'Data', 'RTALogs'))

        # Process per lane.

        for lane in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):

            lane_str = str(lane)
            prefix_lane = string.join(words=(drms_pre_process_folder.name, self.project_name, lane_str), sep='_')

            base_call_lane_path = os.path.join(
                self.run_directory,
                'Data',
                'Intensities',
                'BaseCalls',
                'L{:03d}'.format(lane))

            gzip_bcls = RunnableStep(
                name='gzip_bcls',
                program='find',
                sub_command=Command(command=base_call_lane_path))

            gzip_bcls.add_option_short(key='name', value='*.bcl')
            gzip_bcls.add_option_short(key='execdir', value='gzip --best {} +')

            # Run GNU Tar over the intensities directory.
            # tar -c -f "${archive_prefix}_${lane_name}.tar" "${lane_directory}/"
            archive_intensities = RunnableStep(
                name='archive_intensities',
                program='tar')
            archive_intensities.add_switch_long(key='create')
            # TODO: Need the archive folder.
            # archive_intensities.add_option_long(key='file', value=)
            # Exclude the *.clocs files from the intensities archive, since they are required for
            # extracting base calls.
            archive_intensities.add_option_long(key='exclude', value='*.clocs')
            # TODO: Where are *.locs files stored?
            # Since *_pos.txt files seem to be in Data/Intensities, they do not need excluding.
            # archive_intensities.add_option_long(key='exclude', value='*_pos.txt')

        # Post-process on folder level.

        post_process_folder_prefix = string.join(words=(drms_post_process_folder.name, self.project_name), sep='_')

        file_path_dict_post_process_folder = dict(
            temporary_directory=string.join((post_process_folder_prefix, 'temporary'), sep='_'),
        )

        # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
        runnable_post_process_folder = Runnable(
            name=post_process_folder_prefix,
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict)
        self.add_runnable(runnable=runnable_post_process_folder)

        # Create an Executable for pre-processing the folder.

        executable_post_process_folder = Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_post_process_folder.name)
        drms_post_process_folder.add_executable(executable=executable_post_process_folder)

        # TODO: The post-processing depend on other processes.
        # executable_pre_process_folder.dependencies.extend(vc_process_cohort_dependencies)

        # Archive the Illumina Run Folder

        archive_run_folder = RunnableStep(
            name='archive_run_folder',
            program='tar')
        runnable_post_process_folder.add_runnable_step(runnable_step=archive_run_folder)

        archive_run_folder.add_switch_long(key='create')
        # TODO: Need the archive folder.
        # archive_run_folder.add_option_long(key='file', value=)

        # X. Set file permissions.
        # chmod -R a-w,o-rx

        set_permissions = RunnableStep(
            name='set_permissions',
            program='chmod')
        runnable_post_process_folder.add_runnable_step(runnable_step=set_permissions)

        set_permissions.add_switch_long(key='recursive')
        set_permissions.arguments.append('a-w,o-rx')
        set_permissions.arguments.append(self.run_directory)


class IlluminaRunFolderRestore(Analysis):
    """The C{IlluminaRunFolderRestore} class represents the logic to restore an Illumina Run Folder from a format
    suitable for magnetic tape libraries.

    Attributes:
    @cvar drms_name_extract_archive: C{DRMS.name} for the C{IlluminaRunFolderRestore} C{Analysis} stage
    @type drms_name_extract_archive: str
    @cvar maximum_lane_number: Maximum number of lanes
    @type maximum_lane_number: int
    @ivar archive_directory: File path to an archive directory
    @type archive_directory: str | unicode
    @ivar illumina_directory: File path to the directory of I{Illumina Run Folder} directories
    @type illumina_directory: str | unicode
    @ivar experiment_name: Experiment name (i.e. flow-cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str
    @ivar extract_intensities: Extract cluster intensity file (*.cif) directories
    @type extract_intensities: bool
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool
    @ivar _run_directory_name: Illumina Run Folder directory name
    @type _run_directory_name: str
    @ivar _run_directory_path: Illumina Run Folder directory path
    @type _run_directory_path: str
    """

    drms_name_extract_archive = 'irf_restore_extract_archive'
    drms_name_compress_base_calls = 'irf_restore_compress_base_calls'
    drms_name_compress_logs = 'irf_restore_compress_logs'

    maximum_lane_number = 8

    @classmethod
    def get_prefix_extract_archive(cls, project_name, lane):
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_extract_archive, project_name, lane), sep='_')

    @classmethod
    def get_prefix_compress_base_calls(cls, project_name, lane):
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_compress_base_calls, project_name, lane), sep='_')

    @classmethod
    def get_prefix_compress_logs(cls, project_name):
        return string.join(words=(cls.drms_name_compress_logs, project_name), sep='_')

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{IlluminaRunFolderRestore} object from a UNIX-style configuration file via the
        C{Configuration} class.

        @param config_path: UNIX-style configuration file
        @type config_path: str | unicode
        @return: C{IlluminaRunFolderRestore}
        @rtype: IlluminaRunFolderRestore
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{IlluminaRunFolderRestore} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{IlluminaRunFolderRestore}
        @rtype: IlluminaRunFolderRestore
        """

        assert isinstance(configuration, Configuration)

        irf_restore = cls(configuration=configuration)

        # A "bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore" section specifies defaults
        # for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        irf_restore.set_configuration(irf_restore.configuration, section=section)

        return irf_restore

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 archive_directory=None, illumina_directory=None, experiment_name=None,
                 extract_intensities=None, force=False):
        """Initialise a C{IlluminaRunFolderRestore} object.

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
        @param archive_directory: File path to an archive directory
        @type archive_directory: str | unicode
        @param illumina_directory: File path to the directory of I{Illumina Run Folder} directories
        @type illumina_directory: str | unicode
        @param experiment_name: Experiment name (i.e. flow-cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param extract_intensities: Extract cluster intensity file (*.cif) directories
        @type extract_intensities: bool
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
        """

        super(IlluminaRunFolderRestore, self).__init__(
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

        if archive_directory:
            self.archive_directory = archive_directory
        else:
            self.archive_directory = str()

        if illumina_directory:
            self.illumina_directory = illumina_directory
        else:
            self.illumina_directory = str()

        if experiment_name:
            self.experiment_name = experiment_name
        else:
            self.experiment_name = str()

        if extract_intensities:
            self.extract_intensities = extract_intensities
        else:
            self.extract_intensities = False

        self.force = force

        self._run_directory_name = None
        self._run_directory_path = None

    @property
    def get_run_directory_name(self):

        if not self._run_directory_name:
            self.archive_directory = os.path.normpath(self.archive_directory)
            # Strip the '_archive' suffix i.e. the last 8 characters to get the run directory.
            self._run_directory_name = os.path.basename(self.archive_directory)[:-8]

        return self._run_directory_name

    @property
    def get_run_directory_path(self):

        if not self._run_directory_path:
            self._run_directory_path = os.path.join(self.illumina_directory, self.get_run_directory_name)

        return self._run_directory_path

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{IlluminaRunFolderRestore} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        super(IlluminaRunFolderRestore, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the archive directory.

        if configuration.config_parser.has_option(section=section, option='archive_directory'):
            self.archive_directory = configuration.config_parser.get(
                section=section,
                option='archive_directory')

        # Get Illumina Run Folder information.

        if configuration.config_parser.has_option(section=section, option='illumina_directory'):
            self.illumina_directory = configuration.config_parser.get(
                section=section,
                option='illumina_directory')

        # Get the experiment name.

        if configuration.config_parser.has_option(section=section, option='experiment_name'):
            self.experiment_name = configuration.config_parser.get(
                section=section,
                option='experiment_name')

        if configuration.config_parser.has_option(section=section, option='expand_intensities'):
            self.extract_intensities = configuration.config_parser.getboolean(
                section=section,
                option='expand_intensities')

        if configuration.config_parser.has_option(section=section, option='force'):
            self.force = configuration.config_parser.getboolean(
                section=section,
                option='force')

    def run(self):
        """Run this C{IlluminaRunFolderRunnable} C{Analysis}.

        Restore an Illumina Run Folder from a format suitable for magnetic tape libraries.
        """

        if not os.path.isdir(self.archive_directory):
            raise Exception('The Illumina run archive {!r} does not exist.'.format(self.archive_directory))

        if not self.illumina_directory:
            self.illumina_directory = Default.absolute_runs_illumina()

        if not os.path.isdir(self.illumina_directory):
            raise Exception('The directory of Illumina Run Folder directories {!r} does not exist.'.
                            format(self.illumina_directory))

        # The Illumina Run Folder must *not* already exist unless the force option has been set.

        if os.path.exists(self.get_run_directory_path) and not self.force:
            raise Exception('The Illumina Run Folder directory {!r} exists already.'.
                            format(self.get_run_directory_path))

        archive_name = dict(
            folder=string.join(words=(self.get_run_directory_name, 'Folder.tar'), sep='_'),
            intensities=string.join(words=(self.get_run_directory_name, 'Intensities.tar'), sep='_'))

        for lane in range(0 + 1, self.maximum_lane_number + 1):
            archive_name['L{:03d}'.format(lane)] = string.join(
                words=(self.get_run_directory_name, 'L{:03d}.tar'.format(lane)),
                sep='_')

        file_path_dict = dict()

        for key in archive_name.keys():
            file_path_dict[key] = os.path.join(self.archive_directory, archive_name[key])

        # At least the *_Folder.tar, *_Intensities.tar and *_L001.tar files have to be there.

        if not os.path.exists(file_path_dict['folder']):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(archive_name['folder']))

        if not os.path.exists(file_path_dict['intensities']):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(archive_name['intensities']))

        if not os.path.exists(file_path_dict['L001']):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(archive_name['L001']))

        super(IlluminaRunFolderRestore, self).run()

        drms_extract_archive = DRMS.from_analysis(
            name=self.drms_name_extract_archive,
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(drms_extract_archive)

        drms_compress_base_calls = DRMS.from_analysis(
            name=self.drms_name_compress_base_calls,
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(drms_compress_base_calls)

        drms_compress_logs = DRMS.from_analysis(
            name=self.drms_name_compress_logs,
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(drms_compress_logs)

        # Extract the *_Folder.tar file.

        runnable_extract_folder = Runnable(
            name=self.get_prefix_extract_archive(project_name=self.project_name, lane='folder'),
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict)
        self.add_runnable(runnable=runnable_extract_folder)

        executable_extract_folder = Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_extract_folder.name)
        drms_extract_archive.add_executable(executable=executable_extract_folder)

        extract_folder = RunnableStep(name='extract_folder', program='tar')
        runnable_extract_folder.add_runnable_step(runnable_step=extract_folder)

        extract_folder.add_switch_long(key='extract')
        extract_folder.add_option_long(key='directory', value=self.illumina_directory)
        extract_folder.add_option_long(key='file', value=file_path_dict['folder'])

        # Compress all files in the Log and Data/RTALogs directories.

        runnable_compress_logs = Runnable(
            name=self.get_prefix_compress_logs(project_name=self.project_name),
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict)
        self.add_runnable(runnable=runnable_compress_logs)

        executable_compress_logs = Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_compress_logs.name)
        drms_compress_logs.add_executable(executable=executable_compress_logs)
        # Wait for the restore of the *_Folder.tar file.
        executable_compress_logs.dependencies.append(executable_extract_folder.name)

        compress_logs = RunnableStep(
            name='compress_logs',
            program='gzip')
        runnable_compress_logs.add_runnable_step(runnable_step=compress_logs)

        compress_logs.add_switch_long(key='best')
        compress_logs.add_switch_long(key='recursive')
        compress_logs.arguments.append(os.path.join(self.get_run_directory_path, 'Logs'))

        compress_rta_logs = RunnableStep(
            name='compress_rta_logs',
            program='gzip')
        runnable_compress_logs.add_runnable_step(runnable_step=compress_rta_logs)

        compress_rta_logs.add_switch_long(key='best')
        compress_rta_logs.add_switch_long(key='recursive')
        compress_rta_logs.arguments.append(os.path.join(self.get_run_directory_path, 'Data', 'RTALogs'))

        # Extract the *_intensities.tar file.

        runnable_extract_intensities = Runnable(
            name=self.get_prefix_extract_archive(project_name=self.project_name, lane='intensities'),
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict)
        self.add_runnable(runnable=runnable_extract_intensities)

        executable_extract_intensities = Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_extract_intensities.name)
        drms_extract_archive.add_executable(executable=executable_extract_intensities)

        extract_intensities = RunnableStep(name='extract_intensities', program='tar')
        runnable_extract_intensities.add_runnable_step(runnable_step=extract_intensities)

        extract_intensities.add_switch_long(key='extract')
        extract_intensities.add_option_long(key='directory', value=self.illumina_directory)
        extract_intensities.add_option_long(key='file', value=file_path_dict['intensities'])

        # Unpack the *_L001.tar to *_L008.tar files if they exist.

        for lane in range(0 + 1, self.maximum_lane_number + 1):

            # HiSeq and MiSeq instruments offer various run modes with differing number of lanes.
            if not os.path.exists(file_path_dict['L{:03d}'.format(lane)]):
                continue

            runnable_extract_lane = Runnable(
                name=self.get_prefix_extract_archive(project_name=self.project_name, lane=str(lane)),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict)
            self.add_runnable(runnable=runnable_extract_lane)

            executable_extract_lane = Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_extract_lane.name)
            drms_extract_archive.add_executable(executable=executable_extract_lane)

            extract_lane = RunnableStep(name='extract_lane', program='tar')
            runnable_extract_lane.add_runnable_step(runnable_step=extract_lane)

            extract_lane.add_switch_long(key='extract')
            extract_lane.add_option_long(key='directory', value=self.illumina_directory)
            extract_lane.add_option_long(key='file', value=file_path_dict['L{:03d}'.format(lane)])

            if not self.extract_intensities:
                extract_lane.add_option_long(key='exclude', value='C*.1')

            # Create one process per lane to compress the base call (*.bcl) files.

            runnable_compress_base_calls = Runnable(
                name=self.get_prefix_compress_base_calls(project_name=self.project_name, lane=str(lane)),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict)
            self.add_runnable(runnable=runnable_compress_base_calls)

            executable_compress_base_calls = Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_compress_base_calls.name)
            drms_compress_base_calls.add_executable(executable=executable_compress_base_calls)
            # Wait for the restore of the *_Intensities.tar file.
            executable_compress_base_calls.dependencies.append(executable_extract_intensities.name)

            # Since find has a somewhat broken syntax, definition of the RunnableStep gets a bit complex.
            # The following command line needs to be run
            # find Data/Intensities/BaseCalls/L000 -name '*.bcl' -execdir gzip --best {} +
            compress_base_calls = RunnableStep(
                name='compress_base_calls',
                program='find',
                sub_command=Command(
                    command=os.path.join(
                        self.get_run_directory_path, 'Data', 'Intensities', 'BaseCalls', 'L{:03d}'.format(lane)),
                    sub_command=Command(
                        command='-execdir',
                        sub_command=Command(
                            command='gzip'))))
            runnable_compress_base_calls.add_runnable_step(runnable_step=compress_base_calls)

            find_command = compress_base_calls.sub_command  # directory option
            find_command.add_option_short(key='name', value='*.bcl')
            exec_dir_command = find_command.sub_command  # -execdir option
            gzip_command = exec_dir_command.sub_command  # gzip command
            gzip_command.add_switch_long(key='best')
            gzip_command.arguments.append('{}')
            gzip_command.arguments.append('+')
