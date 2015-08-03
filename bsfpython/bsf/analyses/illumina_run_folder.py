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


import errno
import os
import string

from bsf import Analysis, Command, Configuration, Default, DRMS, Executable, Runnable, RunnableStep, RunnableStepSleep
from bsf.illumina import RunFolder, RunFolderNotComplete


class IlluminaRunFolderArchive(Analysis):
    """The C{IlluminaRunFolderArchive} class represents the logic to archive an Illumina Run Folder in a format
    suitable for magnetic tape libraries.

    Attributes:
    @cvar drms_name_pre_process: C{DRMS.name} for the C{IlluminaRunFolderArchive} C{Analysis} stage
    @type drms_name_pre_process: str
    @type drms_name_base_calls: C{DRMS.name} for the C{IlluminaRunFolderArchive} C{Analysis} stage
    @cvar drms_name_base_calls: str
    @type drms_name_intensities: C{DRMS.name} for the C{IlluminaRunFolderArchive} C{Analysis} stage
    @cvar drms_name_intensities: str
    @cvar drms_name_archive_folder: C{DRMS.name} for the C{IlluminaRunFolderArchive} C{Analysis} stage
    @type drms_name_archive_folder: str
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | unicode
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool
    """

    drms_name_pre_process = 'irf_archive_pre_process'
    drms_name_base_calls = 'irf_archive_base_calls'
    drms_name_intensities = 'irf_archive_intensities'
    drms_name_archive_folder = 'irf_archive_folder'

    compress_archive_files = True

    @classmethod
    def get_prefix_pre_process(cls, project_name):
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_pre_process, project_name), sep='_')

    @classmethod
    def get_prefix_base_calls(cls, project_name, lane):
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_base_calls, project_name, lane), sep='_')

    @classmethod
    def get_prefix_intensities(cls, project_name, lane):
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_intensities, project_name, lane), sep='_')

    @classmethod
    def get_prefix_archive_folder(cls, project_name):
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
        return string.join(words=(cls.drms_name_archive_folder, project_name), sep='_')

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
                 archive_directory=None, run_directory=None, experiment_name=None,
                 force=False):
        """Initialise an C{IlluminaRunFolderArchive} C{Analysis}.

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
        @param archive_directory: Archive directory
        @type archive_directory: str | unicode
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
        @return:
        @rtype:
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

        if archive_directory:
            self.archive_directory = archive_directory
        else:
            self.archive_directory = str()

        if run_directory:
            self.run_directory = run_directory
        else:
            self.run_directory = str()

        if experiment_name:
            self.experiment_name = experiment_name
        else:
            self.experiment_name = str()

        self.force = force

        self._run_name = str()

        return

    @property
    def get_run_name(self):
        """Get the Illumina Run Folder name.

        @return: Illumina Run Folder name
        @rtype: str | unicode
        """

        if not self._run_name:
            self._run_name = os.path.basename(self.run_directory)

        return self._run_name

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{IlluminaRunFolderArchive} C{Analysis} via a section of a
        C{Configuration} object.

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

        super(IlluminaRunFolderArchive, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        if configuration.config_parser.has_option(section=section, option='archive_directory'):
            self.archive_directory = configuration.config_parser.get(
                section=section,
                option='archive_directory')

        # Get Illumina Run Folder information.

        if configuration.config_parser.has_option(section=section, option='run_directory'):
            self.run_directory = configuration.config_parser.get(
                section=section,
                option='run_directory')

        # Get the experiment name.

        if configuration.config_parser.has_option(section=section, option='experiment_name'):
            self.experiment_name = configuration.config_parser.get(
                section=section,
                option='experiment_name')

        # Get the force flag.

        if configuration.config_parser.has_option(section=section, option='force'):
            self.force = configuration.config_parser.getboolean(
                section=section,
                option='force')

        return

    def run(self):
        """Run this C{IlluminaRunFolderArchive} C{Analysis}.

        Archive an I{Illumina Run Folder} in a format suitable for magnetic tape libraries.

            1. Check if the Illumina Run has finished by testing for an
                RTAComplete.txt file.
            2. Check if an archive process is already running by testing for an
                archive directory.
            3. Create an archive directory.
            4. Reset the file permissions for all directories via the find utility.
                find . -type d -execdir chmod u=rwx,g=rx,o= {} \+
            5. Reset the file permissions for all regular files via the find utility.
                find . -type f -execdir chmod u=rw,g=r,o= {} \+
            6. Compress all files in the IRF/Logs/ directory.
                gzip --best --recursive Logs/
            7. Compress all files in the IRF/Data/RTALogs/ directory if it exists.
                gzip --best --recursive IRF/Data/RTALogs/
            8. Compress all IRF/Data/Intensities/BaseCalls/L00[1-8]/C1.1/*.bcl files.
                find IRF/Data/Intensities/BaseCalls/L00[1-8] -name '*.bcl' -execdir gzip --best --verbose {} \+
            9. Run the GNU Tar utility over each IRF/Data/Intensities/L00[1-8]/ directory,
               but exclude compressed cluster locations (*.clocs) files.
            10. Run the GNU Tar utility over the remaining Illumina Run folder,
                but exclude directories with cluster intensity (*.cif) files.
            11. Calculate an MD5 checksum.
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

        self.run_directory = Default.get_absolute_path(
            file_path=self.run_directory,
            default_path=Default.absolute_runs_illumina())

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The Illumina run directory {!r} does not exist.'.format(self.run_directory))

        # Check whether the Illumina Run Folder is complete.
        # 1. Check whether the IRF/RTAComplete.txt file exists in the Illumina Run Folder
        # to prevent archiving and deleting of an incomplete folder.
        # Alternatively, require force to start archiving.

        if not (os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(
                'The Illumina run directory {!r} is not complete.'.format(self.run_directory))

        # Define an Illumina Run Folder archive directory.

        if self.archive_directory:
            # If a relative path to an archive directory has been explicitly defined,
            # prepend it with the parent directory of the Illumina Run Folder (run_folder)
            # to have the run and archive directories in the same directory.
            self.archive_directory = Default.get_absolute_path(
                file_path=self.archive_directory,
                default_path=os.path.dirname(self.run_directory))

            # Raise an Exception if the archive directory does not exist at this stage.
            if not os.path.isdir(self.archive_directory):
                raise Exception('The archive directory {!r} does not exist.'.format(self.archive_directory))
        else:
            # If an archive directory has not been defined, simply append 'archive' to the run directory.
            self.archive_directory = string.join(words=(self.run_directory, 'archive'), sep='_')

        # Check that the directory above the archive directory exists to avoid creation of rogue paths.

        if not os.path.isdir(os.path.dirname(self.archive_directory)):
            raise Exception('The directory above the archive directory {!r} does not exist.'.
                            format(os.path.dirname(self.archive_directory)))

        # 2. Check if a process is already running by testing for an archive directory.
        # 3. Create the archive directory.

        if os.path.isdir(self.archive_directory) and not self.force:
            raise Exception('An archive directory {!r} exists already.'.format(self.archive_directory))
        else:
            try:
                os.makedirs(self.archive_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        irf = RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g. BSF_0000) is used as part of the project_name.
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

        drms_pre_process_folder = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_pre_process,
            working_directory=self.project_directory,
            analysis=self))

        drms_compress_base_calls = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_base_calls,
            working_directory=self.project_directory,
            analysis=self))

        drms_archive_intensities = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_intensities,
            working_directory=self.project_directory,
            analysis=self))

        drms_archive_folder = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_archive_folder,
            working_directory=self.project_directory,
            analysis=self))

        # Pre-process on folder level.

        runnable_pre_process_folder = self.add_runnable(
            runnable=Runnable(
                name=self.get_prefix_pre_process(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))

        executable_pre_process_folder = drms_pre_process_folder.add_executable(
            executable=Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_pre_process_folder.name))

        # executable_pre_process_folder.dependencies.extend()

        # TODO:
        # 0. Check whether Picard ExtractIlluminaBarcodes has written any
        # s_<lane>_<tile>_barcode.txt(.gz) files into the BaseCalls directory.
        # Keeping them is rather pointless and they should be removed.
        # http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes

        # 4. Reset all file permissions for directories.

        reset_directory_permissions = runnable_pre_process_folder.add_runnable_step(
            runnable_step=RunnableStep(
                name='reset_directory_permissions',
                program='find',
                sub_command=Command(
                    command=self.run_directory,
                    sub_command=Command(
                        command='-execdir',
                        sub_command=Command(
                            command='chmod')))))

        find_command = reset_directory_permissions.sub_command  # directory option
        find_command.add_option_short(key='type', value='d')
        exec_dir_command = find_command.sub_command  # -execdir option
        chmod_command = exec_dir_command.sub_command  # chmod command
        chmod_command.arguments.append('u=rwx,g=rx,o=')
        chmod_command.arguments.append('{}')
        chmod_command.arguments.append('+')

        # 5. Reset all file permissions for regular files.

        reset_file_permissions = runnable_pre_process_folder.add_runnable_step(
            runnable_step=RunnableStep(
                name='reset_file_permissions',
                program='find',
                sub_command=Command(
                    command=self.run_directory,
                    sub_command=Command(
                        command='-execdir',
                        sub_command=Command(
                            command='chmod')))))

        find_command = reset_file_permissions.sub_command  # directory option
        find_command.add_option_short(key='type', value='f')
        exec_dir_command = find_command.sub_command  # -execdir option
        chmod_command = exec_dir_command.sub_command  # chmod command
        chmod_command.arguments.append('u=rw,g=r,o=')
        chmod_command.arguments.append('{}')
        chmod_command.arguments.append('+')

        # 6. Compress all files in the IRF/Logs/ and IRF/Logs/IALogs/ directories.

        compress_logs = runnable_pre_process_folder.add_runnable_step(
            runnable_step=RunnableStep(
                name='compress_logs',
                program='gzip'))

        compress_logs.add_switch_long(key='best')
        compress_logs.add_switch_long(key='recursive')
        compress_logs.arguments.append(os.path.join(self.run_directory, 'Logs'))

        # 7. Compress all files in the IRF/Data/RTALogs/ directory if it exists.
        #    It does not on the HiSeq 3000/4000 platform.

        if os.path.isdir(os.path.join(self.run_directory, 'Data', 'RTALogs')):
            compress_rta_logs = runnable_pre_process_folder.add_runnable_step(runnable_step=RunnableStep(
                name='compress_rta_logs',
                program='gzip'))

            compress_rta_logs.add_switch_long(key='best')
            compress_rta_logs.add_switch_long(key='recursive')
            compress_rta_logs.arguments.append(os.path.join(self.run_directory, 'Data', 'RTALogs'))

        # 7. Compress all files in the IRF/RTALogs/ directory if it exists.
        #    It only exists on the HiSeq 3000/4000 platform.

        if os.path.isdir(os.path.join(self.run_directory, 'RTALogs')):
            compress_rta_logs = runnable_pre_process_folder.add_runnable_step(runnable_step=RunnableStep(
                name='compress_rta_logs',
                program='gzip'))

            compress_rta_logs.add_switch_long(key='best')
            compress_rta_logs.add_switch_long(key='recursive')
            compress_rta_logs.arguments.append(os.path.join(self.run_directory, 'RTALogs'))

        # Process per lane.

        # Cluster intensity file (*.cif) directories, if present, need excluding from archiving at a later stage.

        exclude_intensities_patterns = list()
        archive_folder_dependencies = list()

        for lane in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            # Process the IRF/Data/Intensities/BaseCalls/ directory.

            runnable_base_calls = self.add_runnable(
                runnable=Runnable(
                    name=self.get_prefix_base_calls(project_name=self.project_name, lane=str(lane)),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))

            executable_base_calls = drms_compress_base_calls.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_base_calls.name))

            # Set a dependency on the executable_pre_process_folder
            executable_base_calls.dependencies.append(executable_pre_process_folder.name)

            # 8. Compress all base call (*.bcl) files.

            compress_base_calls = runnable_base_calls.add_runnable_step(
                runnable_step=RunnableStep(
                    name='compress_base_calls',
                    program='find',
                    sub_command=Command(
                        command=os.path.join(
                            self.run_directory, 'Data', 'Intensities', 'BaseCalls', 'L{:03d}'.format(lane)),
                        sub_command=Command(
                            command='-execdir',
                            sub_command=Command(
                                command='gzip')))))

            find_command = compress_base_calls.sub_command  # directory option
            find_command.add_option_short(key='name', value='*.bcl')
            exec_dir_command = find_command.sub_command  # -execdir option
            gzip_command = exec_dir_command.sub_command  # gzip command
            gzip_command.add_switch_long(key='best')
            gzip_command.arguments.append('{}')
            gzip_command.arguments.append('+')

            # Record dependencies for the archive run folder analysis stage.
            archive_folder_dependencies.append(executable_base_calls.name)

            if os.path.exists(os.path.join(self.run_directory, 'Data', 'Intensities', 'L{:03d}'.format(lane))):

                # Process IRF/Data/Intensities/L00[1-8]/ directories if they exist.

                runnable_intensities = self.add_runnable(
                    runnable=Runnable(
                        name=self.get_prefix_intensities(project_name=self.project_name, lane=str(lane)),
                        code_module='bsf.runnables.generic',
                        working_directory=self.project_directory))

                executable_intensities = drms_archive_intensities.add_executable(
                    executable=Executable.from_analysis_runnable(
                        analysis=self,
                        runnable_name=runnable_intensities.name))

                # Set a dependency on the executable_pre_process_folder
                executable_intensities.dependencies.append(executable_pre_process_folder.name)

                # 9. Run GNU Tar over the IRF/Data/Intensities/L00[1-8]/ directories, but exclude
                # (compressed) cluster locations (*.clocs and *.locs) files
                # that are essential for extracting base calls and must be archived with the
                # IRF/Data/Intensities/BaseCalls/ directory.
                # Since *_pos.txt files are in IRF/Data/Intensities/, they do not need excluding.

                archive_file_path = string.join(words=(self.get_run_name, 'L{:03d}'.format(lane)), sep='_')
                if self.compress_archive_files:
                    archive_file_path = string.join(words=(archive_file_path, 'tar', 'gz'), sep='.')
                else:
                    archive_file_path = string.join(words=(archive_file_path, 'tar'), sep='.')
                archive_file_path = os.path.join(self.archive_directory, archive_file_path)

                archive_intensities = runnable_intensities.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='archive_intensities',
                        program='tar'))

                archive_intensities.add_switch_long(key='create')
                archive_intensities.add_option_long(key='directory', value=os.path.dirname(self.run_directory))
                archive_intensities.add_option_long(key='file', value=archive_file_path)
                archive_intensities.add_option_long(key='exclude', value='*.clocs')
                archive_intensities.add_option_long(key='exclude', value='*.locs')
                if self.compress_archive_files:
                    archive_intensities.add_switch_long(key='gzip')

                # Archiving needs the relative path.
                archive_intensities.arguments.append(
                    os.path.join(
                        os.path.basename(self.run_directory),
                        'Data',
                        'Intensities',
                        'L{:03d}'.format(lane)))

                # Record dependencies for the archive run folder analysis stage.
                # Since cluster intensity files are no longer automatically deleted,
                # but just excluded from archiving the folder, those dependencies are no longer required.
                # archive_folder_dependencies.append(executable_intensities.name)

                # Record an exclude pattern with the relative path to lane and cycle-specific
                # cluster intensity file (*.cif) directories. On the HiSeq 2000 platform, the lane directories
                # contain cluster locations (*.clocs) files that are essential in base call extraction and
                # need archiving with the IRF/Data/Intensities/BaseCalls folder.
                # By default GNU Tar treats exclusion members as globbing patterns.
                exclude_intensities_patterns.append(
                    os.path.join(
                        os.path.basename(self.run_directory),
                        'Data',
                        'Intensities',
                        'L{:03d}'.format(lane),
                        'C*'))

                # Calculate an MD5 checksum.

                md5_sum = runnable_intensities.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='md5sum',
                        program='md5sum',
                        stdout_path=archive_file_path + '.md5'))

                md5_sum.add_switch_long(key='binary')

                md5_sum.arguments.append(archive_file_path)

        # Process the whole run folder.

        runnable_archive_folder = self.add_runnable(
            runnable=Runnable(
                name=self.get_prefix_archive_folder(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))

        executable_archive_folder = drms_archive_folder.add_executable(
            executable=Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_archive_folder.name))

        executable_archive_folder.dependencies.extend(archive_folder_dependencies)

        # 10. Run the GNU Tar utility over the remaining Illumina Run folder,
        #     but exclude directories with cluster intensity (*.cif) files.

        archive_file_path = self.get_run_name
        if self.compress_archive_files:
            archive_file_path = string.join(words=(archive_file_path, 'tar', 'gz'), sep='.')
        else:
            archive_file_path = string.join(words=(archive_file_path, 'tar'), sep='.')
        archive_file_path = os.path.join(self.archive_directory, archive_file_path)

        archive_folder = runnable_archive_folder.add_runnable_step(
            runnable_step=RunnableStep(
                name='archive_folder',
                program='tar'))

        archive_folder.add_switch_long(key='create')
        archive_folder.add_option_long(key='directory', value=os.path.dirname(self.run_directory))
        archive_folder.add_option_long(key='file', value=archive_file_path)
        for pattern in exclude_intensities_patterns:
            archive_folder.add_option_long(key='exclude', value=pattern)
        if self.compress_archive_files:
            archive_folder.add_switch_long(key='gzip')

        archive_folder.arguments.append(os.path.basename(self.run_directory))

        # 11. Calculate an MD5 checksum.

        md5_sum = runnable_archive_folder.add_runnable_step(
            runnable_step=RunnableStep(
                name='md5_sum',
                program='md5sum',
                stdout_path=archive_file_path + '.md5'))

        md5_sum.add_switch_long(key='binary')

        md5_sum.arguments.append(archive_file_path)

        return


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
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
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
        """Get a process-specific prefix for a C{Runnable} or C{Executable} of this C{Analysis}.

        @param project_name: A project name
        @type project_name: str
        @return: The process-specific prefix for an C{Executable} or C{Runnable} of this C{Analysis}
        @rtype: str
        """
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
        """Initialise an C{IlluminaRunFolderRestore} C{Analysis}.

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
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str
        @param extract_intensities: Extract cluster intensity file (*.cif) directories
        @type extract_intensities: bool
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool
        @return:
        @rtype:
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

        return

    @property
    def get_run_directory_name(self):
        """Get the Illumina Run Folder name.

        @return: Illumina Run Folder name
        @rtype: str | unicode
        """

        if not self._run_directory_name:
            self.archive_directory = os.path.normpath(self.archive_directory)
            # Strip the '_archive' suffix i.e. the last 8 characters to get the run directory.
            self._run_directory_name = os.path.basename(self.archive_directory)[:-8]

        return self._run_directory_name

    @property
    def get_run_directory_path(self):
        """Get the Illumina Run Folder path.

        @return: Illumina Run Folder path
        @rtype: str | unicode
        """

        if not self._run_directory_path:
            self._run_directory_path = os.path.join(self.illumina_directory, self.get_run_directory_name)

        return self._run_directory_path

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{IlluminaRunFolderRestore} C{Analysis} via a section of a
        C{Configuration} object.

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

        return

    def run(self):
        """Run this C{IlluminaRunFolderRunnable} C{Analysis}.

        Restore an Illumina Run Folder from a format suitable for magnetic tape libraries.
            1. Extract the IRF_Folder.tar file.
            2. Extract the IRF_Intensities.tar file with a 60 seconds delay.
            3. Extract each IRF_L00[1-8].tar file with a 90 seconds delay.
        @return:
        @rtype:
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

        # At least the IRF_Folder.tar, IRF_Intensities.tar and IRF_L001.tar files have to be there.

        if not os.path.exists(file_path_dict['folder']):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(archive_name['folder']))

        if not os.path.exists(file_path_dict['intensities']):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(archive_name['intensities']))

        if not os.path.exists(file_path_dict['L001']):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(archive_name['L001']))

        super(IlluminaRunFolderRestore, self).run()

        drms_extract_archive = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_extract_archive,
            working_directory=self.project_directory, analysis=self))

        drms_compress_base_calls = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_compress_base_calls,
            working_directory=self.project_directory,
            analysis=self))

        drms_compress_logs = self.add_drms(drms=DRMS.from_analysis(
            name=self.drms_name_compress_logs,
            working_directory=self.project_directory,
            analysis=self))

        # Extract the IRF_Folder.tar file.

        runnable_extract_folder = self.add_runnable(runnable=Runnable(
            name=self.get_prefix_extract_archive(project_name=self.project_name, lane='folder'),
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict))

        executable_extract_folder = drms_extract_archive.add_executable(executable=Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_extract_folder.name))

        extract_folder = runnable_extract_folder.add_runnable_step(runnable_step=RunnableStep(
            name='extract_folder',
            program='tar'))

        extract_folder.add_switch_long(key='extract')
        extract_folder.add_option_long(key='directory', value=self.illumina_directory)
        extract_folder.add_option_long(key='file', value=file_path_dict['folder'])

        # Compress all files in the IRF/Logs and IRF/Data/RTALogs directories.

        runnable_compress_logs = self.add_runnable(runnable=Runnable(
            name=self.get_prefix_compress_logs(project_name=self.project_name),
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict))

        executable_compress_logs = drms_compress_logs.add_executable(executable=Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_compress_logs.name))

        # Wait for the extraction of the IRF_Folder.tar file.
        executable_compress_logs.dependencies.append(executable_extract_folder.name)

        compress_logs = runnable_compress_logs.add_runnable_step(runnable_step=RunnableStep(
            name='compress_logs',
            program='gzip'))

        compress_logs.add_switch_long(key='best')
        compress_logs.add_switch_long(key='recursive')
        compress_logs.arguments.append(os.path.join(self.get_run_directory_path, 'Logs'))

        compress_rta_logs = runnable_compress_logs.add_runnable_step(runnable_step=RunnableStep(
            name='compress_rta_logs',
            program='gzip'))

        compress_rta_logs.add_switch_long(key='best')
        compress_rta_logs.add_switch_long(key='recursive')
        compress_rta_logs.arguments.append(os.path.join(self.get_run_directory_path, 'Data', 'RTALogs'))

        # Extract the IRF_intensities.tar file.

        runnable_extract_intensities = self.add_runnable(runnable=Runnable(
            name=self.get_prefix_extract_archive(project_name=self.project_name, lane='intensities'),
            code_module='bsf.runnables.generic',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict))

        executable_extract_intensities = drms_extract_archive.add_executable(
            executable=Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_extract_intensities.name))

        # Sleep 60 seconds to allow the first process to create all directories.
        runnable_extract_intensities.add_runnable_step(runnable_step=RunnableStepSleep(name='sleep', sleep_time=60))

        extract_intensities = runnable_extract_intensities.add_runnable_step(runnable_step=RunnableStep(
            name='extract_intensities',
            program='tar'))

        extract_intensities.add_switch_long(key='extract')
        extract_intensities.add_option_long(key='directory', value=self.illumina_directory)
        extract_intensities.add_option_long(key='file', value=file_path_dict['intensities'])

        # Unpack the IRF_L001.tar to IRF_L008.tar files if they exist.

        for lane in range(0 + 1, self.maximum_lane_number + 1):

            # HiSeq and MiSeq instruments offer various run modes with differing number of lanes.
            if not os.path.exists(file_path_dict['L{:03d}'.format(lane)]):
                continue

            runnable_extract_lane = self.add_runnable(runnable=Runnable(
                name=self.get_prefix_extract_archive(project_name=self.project_name, lane=str(lane)),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict))

            drms_extract_archive.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_extract_lane.name))

            # Sleep 90 seconds to allow the first process to create all directories.
            runnable_extract_lane.add_runnable_step(runnable_step=RunnableStepSleep(name='sleep', sleep_time=90))

            extract_lane = runnable_extract_lane.add_runnable_step(runnable_step=RunnableStep(
                name='extract_lane',
                program='tar'))

            extract_lane.add_switch_long(key='extract')
            extract_lane.add_option_long(key='directory', value=self.illumina_directory)
            extract_lane.add_option_long(key='file', value=file_path_dict['L{:03d}'.format(lane)])

            if not self.extract_intensities:
                extract_lane.add_option_long(key='exclude', value='C*.1')

            # Create one process per lane to compress the base call (*.bcl) files.

            runnable_compress_base_calls = self.add_runnable(runnable=Runnable(
                name=self.get_prefix_compress_base_calls(project_name=self.project_name, lane=str(lane)),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict))

            executable_compress_base_calls = drms_compress_base_calls.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_compress_base_calls.name))

            # Wait for the extraction of the IRF_Intensities.tar file.
            executable_compress_base_calls.dependencies.append(executable_extract_intensities.name)

            # Since find has a somewhat broken syntax, definition of the RunnableStep gets a bit complex.
            # The following command line needs to be run
            # find IRF/Data/Intensities/BaseCalls/L000 -name '*.bcl' -execdir gzip --best {} +
            compress_base_calls = runnable_compress_base_calls.add_runnable_step(runnable_step=RunnableStep(
                name='compress_base_calls',
                program='find',
                sub_command=Command(
                    command=os.path.join(
                        self.get_run_directory_path, 'Data', 'Intensities', 'BaseCalls', 'L{:03d}'.format(lane)),
                    sub_command=Command(
                        command='-execdir',
                        sub_command=Command(
                            command='gzip')))))

            find_command = compress_base_calls.sub_command  # directory option
            find_command.add_option_short(key='name', value='*.bcl')
            exec_dir_command = find_command.sub_command  # -execdir option
            gzip_command = exec_dir_command.sub_command  # gzip command
            gzip_command.add_switch_long(key='best')
            gzip_command.arguments.append('{}')
            gzip_command.arguments.append('+')

        return
