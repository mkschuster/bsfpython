# -*- coding: utf-8 -*-
"""Illumina Run Folder Analysis module

A package of classes and methods supporting analyses to archive and restore Illumina Run Folders.
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

import bsf
import bsf.illumina
import bsf.procedure
import bsf.process
import bsf.standards


class IlluminaRunFolderArchive(bsf.Analysis):
    """The C{bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive} class represents the logic to archive
    an Illumina Run Folder in a format suitable for magnetic tape libraries.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar compress_archive_files: Compress archive files with GNU Zip
    @type compress_archive_files: bool
    @ivar archive_directory: Archive directory
    @type archive_directory: str | unicode | None
    @ivar run_directory: File path to an I{Illumina Run Folder}
    @type run_directory: str | unicode | None
    @ivar experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
        Illumina Run Folder parameters
    @type experiment_name: str | None
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool | None
    """

    name = 'Illumina Run Folder Archive Analysis'
    prefix = 'irf_archive'

    compress_archive_files = True

    @classmethod
    def get_stage_name_pre_process(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'pre_process'))

    @classmethod
    def get_stage_name_base_calls(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'base_calls'))

    @classmethod
    def get_stage_name_intensities(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'intensities'))

    @classmethod
    def get_stage_name_archive_folder(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'folder'))

    @classmethod
    def get_prefix_pre_process(cls, project_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_pre_process(), project_name))

    @classmethod
    def get_prefix_base_calls(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_base_calls(), project_name, lane))

    @classmethod
    def get_prefix_intensities(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_intensities(), project_name, lane))

    @classmethod
    def get_prefix_archive_folder(cls, project_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_archive_folder(), project_name))

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
            archive_directory=None,
            run_directory=None,
            experiment_name=None,
            force=None):
        """Initialise a C{bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive} C{bsf.Analysis}.

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
        @param archive_directory: Archive directory
        @type archive_directory: str | unicode | None
        @param run_directory: File path to an I{Illumina Run Folder}
        @type run_directory: str | unicode | None
        @param experiment_name: Experiment name (i.e. flow cell identifier) normally automatically read from
            Illumina Run Folder parameters
        @type experiment_name: str | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
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
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list)

        # Sub-class specific ...

        self.archive_directory = archive_directory
        self.run_directory = run_directory
        self.experiment_name = experiment_name
        self.force = force

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive} C{bsf.Analysis}
        via a section of a C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(IlluminaRunFolderArchive, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'archive_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.archive_directory = configuration.config_parser.get(section=section, option=option)

        # Get Illumina Run Folder information.

        option = 'run_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.run_directory = configuration.config_parser.get(section=section, option=option)

        # Get the experiment name.

        option = 'experiment_name'
        if configuration.config_parser.has_option(section=section, option=option):
            self.experiment_name = configuration.config_parser.get(section=section, option=option)

        # Get the force flag.

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run a C{bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive} C{bsf.Analysis}.

        Archive an I{Illumina Run Folder} in a format suitable for magnetic tape libraries.

        Pre-Process Illumina Run Folder:
            1. Check if the Illumina Run has finished by testing for an
                RTAComplete.txt file.
            2. Check if an archive process is already running by testing for an
                archive directory.
            3. Reset the file permissions for all directories via the find utility.
                find . -type d -execdir chmod u=rwx,g=rx,o= {} \+
            4. Reset the file permissions for all regular files via the find utility.
                find . -type f -execdir chmod u=rw,g=r,o= {} \+
            5. Compress all files in the IRF/Logs/ directory.
                gzip --best --recursive Logs/
            6. Compress all files in the IRF/Data/RTALogs/ directory if it exists.
                gzip --best --recursive IRF/Data/RTALogs/
            7. Compress all files in the IRF/RTALogs/ directory if it exists.
                gzip --best --recursive IRF/RTALogs/
            8. Create the archive directory.
        Lane specific:
            1. Compress all IRF/Data/Intensities/BaseCalls/L00[1-8]/C1.1/*.bcl files.
                find IRF/Data/Intensities/BaseCalls/L00[1-8] -name '*.bcl' -execdir gzip --best --verbose {} \+
            2. Run the GNU Tar utility over each IRF/Data/Intensities/L00[1-8]/ directory,
               but exclude compressed cluster locations (*.clocs) files.
            3. Calculate a MD5 checksum.
        Illumina Run Folder-specific:
            1. Run the GNU Tar utility over the remaining Illumina Run folder,
               but exclude directories with cluster intensity (*.cif) files.
            2. Calculate a MD5 checksum.
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
            default_path=bsf.standards.FilePath.get_illumina_run(absolute=True))

        # The Illumina Run Folder name would also be available from the bsf.illumina.RunFolder.get_name property below,
        # but using the file path provides more flexibility outside of the fixed runParameters.xml configuration file.

        run_name = os.path.basename(self.run_directory)

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(
                'The Illumina run directory {!r} does not exist.'.format(self.run_directory))

        # Check whether the Illumina Run Folder is complete.
        # 1. Check whether the IRF/RTAComplete.txt file exists in the Illumina Run Folder
        # to prevent archiving and deleting of an incomplete folder.
        # Alternatively, require force to start archiving.

        if not (os.path.exists(path=os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise bsf.illumina.RunFolderNotComplete(
                'The Illumina run directory {!r} is not complete.'.format(self.run_directory))

        # Define an Illumina Run Folder archive directory.

        if self.archive_directory:
            # If a relative path to an archive directory has been explicitly defined,
            # prepend it with the parent directory of the Illumina Run Folder (run_folder)
            # to have the run and archive directories in the same directory.
            self.archive_directory = self.configuration.get_absolute_path(
                file_path=self.archive_directory,
                default_path=os.path.dirname(self.run_directory))
        else:
            # If an archive directory has not been defined, simply append 'archive' to the run directory path.
            self.archive_directory = '_'.join((self.run_directory, 'archive'))

        # Check that the directory above the archive directory exists to avoid creation of rogue paths.

        if not os.path.isdir(os.path.dirname(self.archive_directory)):
            raise Exception('The directory above the archive directory {!r} does not exist.'.
                            format(os.path.dirname(self.archive_directory)))

        # 2. Check if a process is already running by testing for an archive directory.

        if os.path.isdir(self.archive_directory) and not self.force:
            raise Exception('An archive directory {!r} exists already.'.format(self.archive_directory))

        irf = bsf.illumina.RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g. BSF_0000) is used as part of the project_name.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = irf.run_parameters.get_experiment_name
            if not self.experiment_name:
                raise Exception('An experiment_name has not been defined.')

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, irf.run_information.flow_cell))

        super(IlluminaRunFolderArchive, self).run()

        stage_pre_process_folder = self.get_stage(name=self.get_stage_name_pre_process())
        stage_compress_base_calls = self.get_stage(name=self.get_stage_name_base_calls())
        stage_archive_intensities = self.get_stage(name=self.get_stage_name_intensities())
        stage_archive_folder = self.get_stage(name=self.get_stage_name_archive_folder())

        # Pre-process on folder level.

        runnable_pre_process_folder = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_pre_process(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))

        executable_pre_process_folder = self.set_stage_runnable(
            stage=stage_pre_process_folder,
            runnable=runnable_pre_process_folder)

        # 0. Check whether Picard ExtractIlluminaBarcodes has written any
        # s_<lane>_<tile>_barcode.txt(.gz) files into the BaseCalls directory.
        # Keeping them is rather pointless and they should be removed.
        # http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes

        # 3. Reset all file permissions for directories.

        reset_directory_permissions = runnable_pre_process_folder.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='reset_directory_permissions',
                program='find',
                sub_command=bsf.process.Command(
                    program=self.run_directory,
                    sub_command=bsf.process.Command(
                        program='-execdir',
                        sub_command=bsf.process.Command(
                            program='chmod')))))

        find_command = reset_directory_permissions.sub_command  # directory option
        find_command.add_option_short(key='type', value='d')
        exec_dir_command = find_command.sub_command  # -execdir option
        chmod_command = exec_dir_command.sub_command  # chmod command
        chmod_command.arguments.append('u=rwx,g=rx,o=')
        chmod_command.arguments.append('{}')
        chmod_command.arguments.append('+')

        # 4. Reset all file permissions for regular files.

        reset_file_permissions = runnable_pre_process_folder.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='reset_file_permissions',
                program='find',
                sub_command=bsf.process.Command(
                    program=self.run_directory,
                    sub_command=bsf.process.Command(
                        program='-execdir',
                        sub_command=bsf.process.Command(
                            program='chmod')))))

        find_command = reset_file_permissions.sub_command  # directory option
        find_command.add_option_short(key='type', value='f')
        exec_dir_command = find_command.sub_command  # -execdir option
        chmod_command = exec_dir_command.sub_command  # chmod command
        chmod_command.arguments.append('u=rw,g=r,o=')
        chmod_command.arguments.append('{}')
        chmod_command.arguments.append('+')

        # 5. Compress all files in the IRF/Logs/ and IRF/Logs/IALogs/ directories.
        #    The NextSeq compresses into a single Logs.zip file.

        if irf.run_parameters.get_instrument_type not in ('NextSeq',):
            compress_logs = runnable_pre_process_folder.add_runnable_step(
                runnable_step=bsf.process.RunnableStep(
                    name='compress_logs',
                    program='gzip'))

            compress_logs.add_switch_long(key='best')
            compress_logs.add_switch_long(key='recursive')
            compress_logs.arguments.append(os.path.join(self.run_directory, 'Logs'))

        # 6. Compress all files in the IRF/Data/RTALogs/ directory if it exists.
        #    It does not on the HiSeq 3000/4000 and NovaSeq platforms.

        if os.path.isdir(os.path.join(self.run_directory, 'Data', 'RTALogs')):
            compress_rta_logs = runnable_pre_process_folder.add_runnable_step(
                runnable_step=bsf.process.RunnableStep(
                    name='compress_rta_logs',
                    program='gzip'))

            compress_rta_logs.add_switch_long(key='best')
            compress_rta_logs.add_switch_long(key='recursive')
            compress_rta_logs.arguments.append(os.path.join(self.run_directory, 'Data', 'RTALogs'))

        # 7. Compress all files in the IRF/RTALogs/ directory if it exists.
        #    It only exists on the HiSeq 3000/4000 platform.

        if os.path.isdir(os.path.join(self.run_directory, 'RTALogs')):
            compress_rta_logs = runnable_pre_process_folder.add_runnable_step(
                runnable_step=bsf.process.RunnableStep(
                    name='compress_rta_logs',
                    program='gzip'))

            compress_rta_logs.add_switch_long(key='best')
            compress_rta_logs.add_switch_long(key='recursive')
            compress_rta_logs.arguments.append(os.path.join(self.run_directory, 'RTALogs'))

        # 8. Create the archive directory.

        runnable_pre_process_folder.add_runnable_step(
            runnable_step=bsf.process.RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=self.archive_directory))

        # Process per lane.

        # Cluster intensity file (*.cif) directories, if present, need excluding from archiving at a later stage.

        exclude_intensities_patterns = list()
        """ @type exclude_intensities_patterns: list[str | unicode] """
        archive_folder_dependencies = list()
        """ @type archive_folder_dependencies: list[str] """

        archive_folder_dependencies.append(runnable_pre_process_folder.name)

        for lane_int in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            # MiSeq, NextSeq and NovaSeq instruments do not need lane processing.
            if irf.run_parameters.get_instrument_type in ('MiSeq', 'NextSeq', 'NovaSeq'):
                continue

            # Process the IRF/Data/Intensities/BaseCalls/ directory.

            runnable_base_calls = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_base_calls(project_name=self.project_name, lane=str(lane_int)),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))
            executable_base_calls = self.set_stage_runnable(
                stage=stage_compress_base_calls,
                runnable=runnable_base_calls)

            # Set a dependency on the executable_pre_process_folder
            executable_base_calls.dependencies.append(executable_pre_process_folder.name)

            # 1. Compress all base call (*.bcl) files.

            compress_base_calls = runnable_base_calls.add_runnable_step(
                runnable_step=bsf.process.RunnableStep(
                    name='compress_base_calls',
                    program='find',
                    sub_command=bsf.process.Command(
                        program=os.path.join(
                            self.run_directory, 'Data', 'Intensities', 'BaseCalls',
                            'L{:03d}'.format(lane_int)),
                        sub_command=bsf.process.Command(
                            program='-execdir',
                            sub_command=bsf.process.Command(
                                program='gzip')))))

            find_command = compress_base_calls.sub_command  # directory option
            find_command.add_option_short(key='name', value='*.bcl')
            exec_dir_command = find_command.sub_command  # -execdir option
            gzip_command = exec_dir_command.sub_command  # gzip command
            gzip_command.add_switch_long(key='best')
            gzip_command.arguments.append('{}')
            gzip_command.arguments.append('+')

            # Record dependencies for the archive run folder analysis stage.
            archive_folder_dependencies.append(executable_base_calls.name)

            if os.path.exists(os.path.join(self.run_directory, 'Data', 'Intensities', 'L{:03d}'.format(lane_int))):

                # Process IRF/Data/Intensities/L00[1-8]/ directories if they exist.

                runnable_intensities = self.add_runnable_consecutive(
                    runnable=bsf.procedure.ConsecutiveRunnable(
                        name=self.get_prefix_intensities(project_name=self.project_name, lane=str(lane_int)),
                        code_module='bsf.runnables.generic',
                        working_directory=self.project_directory))
                executable_intensities = self.set_stage_runnable(
                    stage=stage_archive_intensities,
                    runnable=runnable_intensities)
                executable_intensities.dependencies.append(executable_pre_process_folder.name)

                # 2. Run GNU Tar over the IRF/Data/Intensities/L00[1-8]/ directories, but exclude
                # (compressed) cluster locations (*.clocs and *.locs) files
                # that are essential for extracting base calls and must be archived with the
                # IRF/Data/Intensities/BaseCalls/ directory.
                # Since *_pos.txt files are in IRF/Data/Intensities/, they do not need excluding.

                archive_file_path = '_'.join((run_name, 'L{:03d}'.format(lane_int)))
                if self.compress_archive_files:
                    archive_file_path = '.'.join((archive_file_path, 'tar', 'gz'))
                else:
                    archive_file_path = '.'.join((archive_file_path, 'tar'))
                archive_file_path = os.path.join(self.archive_directory, archive_file_path)

                archive_intensities = runnable_intensities.add_runnable_step(
                    runnable_step=bsf.process.RunnableStep(
                        name='archive_intensities',
                        program='tar'))

                archive_intensities.add_switch_long(key='create')
                archive_intensities.add_option_long(key='directory', value=os.path.dirname(self.run_directory))
                archive_intensities.add_option_long(key='file', value=archive_file_path)
                archive_intensities.add_option_long(key='exclude', value='*.clocs', override=True)
                archive_intensities.add_option_long(key='exclude', value='*.locs', override=True)
                if self.compress_archive_files:
                    archive_intensities.add_switch_long(key='gzip')

                # Archiving needs the relative path.
                archive_intensities.arguments.append(
                    os.path.join(
                        os.path.basename(self.run_directory),
                        'Data',
                        'Intensities',
                        'L{:03d}'.format(lane_int)))

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
                        'L{:03d}'.format(lane_int),
                        'C*'))

                # 3. Calculate a MD5 checksum.

                md5_sum = runnable_intensities.add_runnable_step(
                    runnable_step=bsf.process.RunnableStep(
                        name='md5sum',
                        program='md5sum',
                        stdout_path=archive_file_path + '.md5'))

                md5_sum.add_switch_long(key='binary')

                md5_sum.arguments.append(archive_file_path)

        # Process the whole run folder.

        runnable_archive_folder = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_archive_folder(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_archive_folder = self.set_stage_runnable(
            stage=stage_archive_folder,
            runnable=runnable_archive_folder)
        executable_archive_folder.dependencies.extend(archive_folder_dependencies)

        # 1. Run the GNU Tar utility over the remaining Illumina Run folder,
        #     but exclude directories with cluster intensity (*.cif) files.

        archive_file_path = run_name
        if self.compress_archive_files:
            archive_file_path = '.'.join((archive_file_path, 'tar', 'gz'))
        else:
            archive_file_path = '.'.join((archive_file_path, 'tar'))
        archive_file_path = os.path.join(self.archive_directory, archive_file_path)

        archive_folder = runnable_archive_folder.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='archive_folder',
                program='tar'))

        archive_folder.add_switch_long(key='create')
        archive_folder.add_option_long(key='directory', value=os.path.dirname(self.run_directory))
        archive_folder.add_option_long(key='file', value=archive_file_path)
        for pattern in exclude_intensities_patterns:
            archive_folder.add_option_long(key='exclude', value=pattern, override=True)
        if self.compress_archive_files:
            archive_folder.add_switch_long(key='gzip')

        archive_folder.arguments.append(os.path.basename(self.run_directory))

        # 2. Calculate a MD5 checksum.

        md5_sum = runnable_archive_folder.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='md5_sum',
                program='md5sum',
                stdout_path=archive_file_path + '.md5'))

        md5_sum.add_switch_long(key='binary')

        md5_sum.arguments.append(archive_file_path)

        return


class FilePathIlluminaRunFolderRestore(bsf.procedure.FilePath):
    """The C{bsf.analyses.illumina_run_folder.FilePathIlluminaRunFolderRestore} models files in an archive directory.

    Attributes:
    @ivar folder: Folder GNU Tar archive file
    @type folder: str | unicode
    @ivar intensities: Intensities GNU Tar archive file
    @type intensities: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.illumina_run_folder.FilePathIlluminaRunFolderRestore} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathIlluminaRunFolderRestore, self).__init__(prefix=prefix)

        self.folder = '_'.join((prefix, 'Folder.tar'))
        self.intensities = '_'.join((prefix, 'Intensities.tar'))
        # There are additional attributes L001 to L008 depending on the flow cell layout.

        return


class IlluminaRunFolderRestore(bsf.Analysis):
    """The C{bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore} class represents the logic to restore an
    Illumina Run Folder from a format suitable for magnetic tape libraries.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar maximum_lane_number: Maximum number of lanes
    @type maximum_lane_number: int
    @ivar archive_directory: File path to an archive directory
    @type archive_directory: str | unicode | None
    @ivar illumina_directory: File path to the directory of I{Illumina Run Folder} directories
    @type illumina_directory: str | unicode | None
    @ivar extract_intensities: Extract cluster intensity file (*.cif) directories
    @type extract_intensities: bool | None
    @ivar force: Force processing of incomplete Illumina Run Folders
    @type force: bool | None
    """

    name = 'Illumina Run Folder Restore Analysis'
    prefix = 'irf_restore'

    maximum_lane_number = 8

    @classmethod
    def get_stage_name_extract_archive(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'extract_archive'))

    @classmethod
    def get_stage_name_compress_base_calls(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'compress_base_calls'))

    @classmethod
    def get_stage_name_compress_logs(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'compress_logs'))

    @classmethod
    def get_prefix_extract_archive(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_extract_archive(), project_name, lane))

    @classmethod
    def get_prefix_compress_base_calls(cls, project_name, lane):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @param lane: A lane number
        @type lane: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_compress_base_calls(), project_name, lane))

    @classmethod
    def get_prefix_compress_logs(cls, project_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param project_name: A project name
        @type project_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_compress_logs(), project_name))

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
            archive_directory=None,
            illumina_directory=None,
            extract_intensities=False,
            force=False):
        """Initialise a C{bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore} C{bsf.Analysis}.

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
        @param archive_directory: File path to an archive directory
        @type archive_directory: str | unicode | None
        @param illumina_directory: File path to the directory of I{Illumina Run Folder} directories
        @type illumina_directory: str | unicode | None
        @param extract_intensities: Extract cluster intensity file (*.cif) directories
        @type extract_intensities: bool | None
        @param force: Force processing of incomplete Illumina Run Folders
        @type force: bool | None
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
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list)

        # Sub-class specific ...

        self.archive_directory = archive_directory
        self.illumina_directory = illumina_directory
        self.extract_intensities = extract_intensities
        self.force = force

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore} C{bsf.Analysis}
        via a section of a C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(IlluminaRunFolderRestore, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the archive directory.

        option = 'archive_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.archive_directory = configuration.config_parser.get(section=section, option=option)

        # Get Illumina Run Folder information.

        option = 'illumina_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.illumina_directory = configuration.config_parser.get(section=section, option=option)

        option = 'expand_intensities'
        if configuration.config_parser.has_option(section=section, option=option):
            self.extract_intensities = configuration.config_parser.getboolean(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore} C{bsf.Analysis}.

        Restore an Illumina Run Folder from a format suitable for magnetic tape libraries.
            1. Extract the IRF_Folder.tar file.
            2. Extract the IRF_Intensities.tar file with a 60 seconds delay.
            3. Extract each IRF_L00[1-8].tar file with a 90 seconds delay.
        @return:
        @rtype:
        """
        if not self.archive_directory:
            raise Exception('The archive_directory has not been defined.')

        self.archive_directory = self.configuration.get_absolute_path(
            file_path=self.archive_directory)

        if not os.path.isdir(self.archive_directory):
            raise Exception('The Illumina run archive {!r} does not exist.'.format(self.archive_directory))

        self.illumina_directory = self.configuration.get_absolute_path(
            file_path=self.illumina_directory,
            default_path=bsf.standards.FilePath.get_illumina_run(absolute=True))

        if not os.path.isdir(self.illumina_directory):
            raise Exception('The directory of Illumina Run Folder directories {!r} does not exist.'.
                            format(self.illumina_directory))

        # The Illumina Run Folder must *not* already exist unless the force option has been set.

        run_directory_name = os.path.basename(self.archive_directory)[:-8]
        run_directory_path = os.path.join(self.illumina_directory, run_directory_name)

        if os.path.exists(run_directory_path) and not self.force:
            raise Exception('The Illumina Run Folder directory {!r} exists already.'.
                            format(run_directory_path))

        # A FilePath object cannot be used, because instruments have a varable number of lanes.
        file_path_dict = {
            'folder': '_'.join((run_directory_name, 'Folder.tar')),
            'intensities': '_'.join((run_directory_name, 'Intensities.tar')),
        }

        # Add additional, lane-specific keys.
        for lane_int in range(0 + 1, self.maximum_lane_number + 1):
            file_path_dict['L{:03d}'.format(lane_int)] = \
                '_'.join((run_directory_name, 'L{:03d}.tar'.format(lane_int)))

        # At least the IRF_Folder.tar, IRF_Intensities.tar and IRF_L001.tar files have to be there.

        if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['folder'])):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(file_path_dict['folder']))

        if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['intensities'])):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(file_path_dict['intensities']))

        if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['L001'])):
            raise Exception('Illumina Run Archive file {!r} is missing.'.format(file_path_dict['L001']))

        super(IlluminaRunFolderRestore, self).run()

        stage_extract_archive = self.get_stage(name=self.get_stage_name_extract_archive())
        stage_compress_base_calls = self.get_stage(name=self.get_stage_name_compress_base_calls())
        stage_compress_logs = self.get_stage(name=self.get_stage_name_compress_logs())

        # Extract the IRF_Folder.tar file.

        runnable_extract_folder = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_extract_archive(project_name=self.project_name, lane='folder'),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_extract_folder = self.set_stage_runnable(
            stage=stage_extract_archive,
            runnable=runnable_extract_folder)

        extract_folder = runnable_extract_folder.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='extract_folder',
                program='tar'))

        extract_folder.add_switch_long(key='extract')
        extract_folder.add_option_long(key='directory', value=self.illumina_directory)
        extract_folder.add_option_long(
            key='file',
            value=os.path.join(self.archive_directory, file_path_dict['folder']))

        # Compress all files in the IRF/Logs and IRF/Data/RTALogs directories.

        runnable_compress_logs = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_compress_logs(project_name=self.project_name),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_compress_logs = self.set_stage_runnable(
            stage=stage_compress_logs,
            runnable=runnable_compress_logs)
        executable_compress_logs.dependencies.append(executable_extract_folder.name)

        compress_logs = runnable_compress_logs.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='compress_logs',
                program='gzip'))

        compress_logs.add_switch_long(key='best')
        compress_logs.add_switch_long(key='recursive')
        compress_logs.arguments.append(os.path.join(run_directory_path, 'Logs'))

        compress_rta_logs = runnable_compress_logs.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='compress_rta_logs',
                program='gzip'))

        compress_rta_logs.add_switch_long(key='best')
        compress_rta_logs.add_switch_long(key='recursive')
        compress_rta_logs.arguments.append(os.path.join(run_directory_path, 'Data', 'RTALogs'))

        # Extract the IRF_intensities.tar file.

        runnable_extract_intensities = self.add_runnable_consecutive(
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=self.get_prefix_extract_archive(project_name=self.project_name, lane='intensities'),
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_extract_intensities = self.set_stage_runnable(
            stage=stage_extract_archive,
            runnable=runnable_extract_intensities)

        # Sleep 60 seconds to allow the first process to create all directories.
        runnable_extract_intensities.add_runnable_step(
            runnable_step=bsf.process.RunnableStepSleep(
                name='sleep',
                sleep_time=60.0))

        extract_intensities = runnable_extract_intensities.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='extract_intensities',
                program='tar'))

        extract_intensities.add_switch_long(key='extract')
        extract_intensities.add_option_long(key='directory', value=self.illumina_directory)
        extract_intensities.add_option_long(
            key='file',
            value=os.path.join(self.archive_directory, file_path_dict['intensities']))

        # Unpack the IRF_L001.tar to IRF_L008.tar files if they exist.

        for lane_int in range(0 + 1, self.maximum_lane_number + 1):

            # HiSeq and MiSeq instruments offer various run modes with differing number of lanes.
            if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['L{:03d}'.format(lane_int)])):
                continue

            runnable_extract_lane = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_extract_archive(project_name=self.project_name, lane=str(lane_int)),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))
            self.set_stage_runnable(
                stage=stage_extract_archive,
                runnable=runnable_extract_lane)

            # Sleep 90 seconds to allow the first process to create all directories.
            runnable_extract_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStepSleep(
                    name='sleep',
                    sleep_time=90.0))

            extract_lane = runnable_extract_lane.add_runnable_step(
                runnable_step=bsf.process.RunnableStep(
                    name='extract_lane',
                    program='tar'))

            extract_lane.add_switch_long(key='extract')
            extract_lane.add_option_long(key='directory', value=self.illumina_directory)
            extract_lane.add_option_long(
                key='file',
                value=os.path.join(self.archive_directory, file_path_dict['L{:03d}'.format(lane_int)]))

            if not self.extract_intensities:
                extract_lane.add_option_long(key='exclude', value='C*.1', override=True)

            # Create one process per lane to compress the base call (*.bcl) files.

            runnable_compress_base_calls = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=self.get_prefix_compress_base_calls(project_name=self.project_name, lane=str(lane_int)),
                    code_module='bsf.runnables.generic',
                    working_directory=self.project_directory))
            executable_compress_base_calls = self.set_stage_runnable(
                stage=stage_compress_base_calls,
                runnable=runnable_compress_base_calls)
            executable_compress_base_calls.dependencies.append(executable_extract_intensities.name)

            # Since find has a somewhat broken syntax, definition of the bsf.process.RunnableStep gets a bit complex.
            # The following command line needs to be run
            # find IRF/Data/Intensities/BaseCalls/L000 -name '*.bcl' -execdir gzip --best {} +
            compress_base_calls = runnable_compress_base_calls.add_runnable_step(
                runnable_step=bsf.process.RunnableStep(
                    name='compress_base_calls',
                    program='find',
                    sub_command=bsf.process.Command(
                        program=os.path.join(
                            run_directory_path, 'Data', 'Intensities', 'BaseCalls',
                            'L{:03d}'.format(lane_int)),
                        sub_command=bsf.process.Command(
                            program='-execdir',
                            sub_command=bsf.process.Command(
                                program='gzip')))))

            find_command = compress_base_calls.sub_command  # directory option
            find_command.add_option_short(key='name', value='*.bcl')
            exec_dir_command = find_command.sub_command  # -execdir option
            gzip_command = exec_dir_command.sub_command  # gzip command
            gzip_command.add_switch_long(key='best')
            gzip_command.arguments.append('{}')
            gzip_command.arguments.append('+')

        return
