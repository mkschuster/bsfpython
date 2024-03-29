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
"""The :py:mod:`bsf.analyses.illumina_run_folder` module provides classes and methods supporting
the archiving and restoring of :emphasis:`Illumina Run Folder` (IRF) objects.
"""
import logging
import os
from argparse import ArgumentParser
from typing import Optional

from bsf.analysis import Analysis, Stage
from bsf.connector import ConnectorFile, ConnectorPipe, ConcurrentProcess
from bsf.executables.cloud import RunnableStepAzureBlockBlobUpload
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.ngs import Collection, Sample
from bsf.procedure import FilePath, ConsecutiveRunnable, ConcurrentRunnable
from bsf.process import Command, \
    RunnableStep, RunnableStepChangeMode, RunnableStepMakeDirectory, RunnableStepMove, RunnableStepSleep
from bsf.standards import Configuration, StandardFilePath


class IlluminaRunFolderArchive(Analysis):
    """The :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive` class represents the logic to convert
    an :emphasis:`Illumina Run Folder` (IRF) into one or more
    `GNU Tar <https://www.gnu.org/software/tar/>`_ archives.

    :cvar compress_archive_files: Request compressing archive files via
        `GNU Zip <https://www.gnu.org/software/gzip/>`_.
    :type compress_archive_files: bool
    :ivar archive_directory: An archive directory path.
    :type archive_directory: str | None
    :ivar run_directory: An :emphasis:`Illumina Run Folder` (IRF) directory path.
    :type run_directory: str | None
    :ivar experiment_name: An experiment name (i.e., flow cell identifier) normally automatically read from
        :emphasis:`Illumina Run Folder` parameters.
    :type experiment_name: str | None
    :ivar irf_mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory of the
            archived :emphasis:`Illumina Run Folder` (IRF).
    :type irf_mode_directory: str | None
    :ivar sav_mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory of the
            archived :emphasis:`Sequence Analysis Viewer` (SAV).
    :type sav_mode_directory: str | None
    :ivar irf_mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file of the
            archived :emphasis:`Illumina Run Folder` (IRF).
    :type irf_mode_file: str | None
    :ivar sav_mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file of the
            archived :emphasis:`Sequence Analysis Viewer` (SAV).
    :type sav_mode_file: str | None
    :ivar cloud_account: A :emphasis:`Microsoft Azure Storage Account` name.
    :type cloud_account: str | None
    :ivar cloud_container: A :emphasis:`Microsoft Azure Blob Service` container name.
    :type cloud_container: str | None
    :ivar cloud_path_prefix: A blob file path prefix.
    :type cloud_path_prefix: str | None
    :ivar cloud_concurrency: A maximum number of concurrent network connections.
    :type cloud_concurrency: int | None
    :ivar force: Request processing of incomplete :emphasis:`Illumina Run Folder` objects.
    :type force: bool | None
    """

    name = 'Illumina Run Folder Archive Analysis'
    prefix = 'irf_archive'

    compress_archive_files = True

    @classmethod
    def get_stage_name_pre_process(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'pre_process'))

    @classmethod
    def get_stage_name_base_calls(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'base_calls'))

    @classmethod
    def get_stage_name_intensities(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'intensities'))

    @classmethod
    def get_stage_name_archive_folder(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'folder'))

    @classmethod
    def get_stage_name_post_process(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'post_process'))

    @classmethod
    def get_stage_name_sav(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'sav'))

    @classmethod
    def get_stage_name_cloud(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cloud'))

    @classmethod
    def get_prefix_pre_process(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_pre_process(), project_name))

    @classmethod
    def get_prefix_base_calls(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_base_calls(), project_name, lane))

    @classmethod
    def get_prefix_intensities(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_intensities(), project_name, lane))

    @classmethod
    def get_prefix_archive_folder(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_archive_folder(), project_name))

    @classmethod
    def get_prefix_post_process(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_post_process(), project_name))

    @classmethod
    def get_prefix_sav(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_sav(), project_name))

    @classmethod
    def get_prefix_cloud(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_cloud(), project_name))

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
            archive_directory: Optional[str] = None,
            run_directory: Optional[str] = None,
            experiment_name: Optional[str] = None,
            irf_mode_directory: Optional[str] = None,
            sav_mode_directory: Optional[str] = None,
            irf_mode_file: Optional[str] = None,
            sav_mode_file: Optional[str] = None,
            cloud_account: Optional[str] = None,
            cloud_container: Optional[str] = None,
            cloud_path_prefix: Optional[str] = None,
            cloud_concurrency: Optional[int] = None,
            force: Optional[bool] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive` object.

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
        :param archive_directory: An archive directory path.
        :type archive_directory: str | None
        :param run_directory: An :emphasis:`Illumina Run Folder` (IRF) directory path.
        :type run_directory: str | None
        :param experiment_name: An experiment name (i.e., flow cell identifier) normally automatically read from
            :emphasis:`Illumina Run Folder` parameters.
        :type experiment_name: str | None
        :param irf_mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory of the
            archived :emphasis:`Illumina Run Folder` (IRF).
        :type irf_mode_directory: str | None
        :param sav_mode_directory: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each directory of the
            archived :emphasis:`Sequence Analysis Viewer` (SAV).
        :type sav_mode_directory: str | None
        :param irf_mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file of the
            archived :emphasis:`Illumina Run Folder` (IRF).
        :type irf_mode_file: str | None
        :param sav_mode_file: A comma-separated list of file permission bit names defined in the
            Python :py:mod:`stat` module applied to each file of the
            archived :emphasis:`Sequence Analysis Viewer` (SAV).
        :type sav_mode_file: str | None
        :param cloud_account: A :emphasis:`Microsoft Azure Storage Account` name.
        :type cloud_account: str | None
        :param cloud_container: A :emphasis:`Microsoft Azure Blob Service` container name.
        :type cloud_container: str | None
        :param cloud_path_prefix: A blob file path prefix.
        :type cloud_path_prefix: str | None
        :param cloud_concurrency: A maximum number of concurrent network connections.
        :type cloud_concurrency: int | None
        :param force: Request processing of incomplete :emphasis:`Illumina Run Folder` objects.
        :type force: bool | None
        """
        super(IlluminaRunFolderArchive, self).__init__(
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

        # Sub-class specific ...

        self.archive_directory = archive_directory
        self.run_directory = run_directory
        self.experiment_name = experiment_name
        self.irf_mode_directory = irf_mode_directory
        self.sav_mode_directory = sav_mode_directory
        self.irf_mode_file = irf_mode_file
        self.sav_mode_file = sav_mode_file
        self.cloud_account = cloud_account
        self.cloud_container = cloud_container
        self.cloud_path_prefix = cloud_path_prefix
        self.cloud_concurrency = cloud_concurrency
        self.force = force

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
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

        option = 'irf_mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.irf_mode_directory = configuration.config_parser.get(section=section, option=option)

        option = 'sav_mode_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sav_mode_directory = configuration.config_parser.get(section=section, option=option)

        option = 'irf_mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.irf_mode_file = configuration.config_parser.get(section=section, option=option)

        option = 'sav_mode_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sav_mode_file = configuration.config_parser.get(section=section, option=option)

        option = 'cloud_account'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_account = configuration.config_parser.get(section=section, option=option)

        option = 'cloud_container'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_container = configuration.config_parser.get(section=section, option=option)

        option = 'cloud_path_prefix'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_path_prefix = configuration.config_parser.get(section=section, option=option)

        option = 'cloud_concurrency'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cloud_concurrency = configuration.config_parser.getint(section=section, option=option)

        option = 'force'
        if configuration.config_parser.has_option(section=section, option=option):
            self.force = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive` object.

        Archive an :emphasis:`Illumina Run Folder` in a format suitable for magnetic tape libraries.

        - Pre-Process an :emphasis:`Illumina Run Folder`:

          1. Check if the :emphasis:`Illumina Run Folder` is complete by testing for an
             :literal:`RTAComplete.txt` file.
          2. Check if an archive process is already running by testing for an
             archive directory.
          3. Reset the file permissions for all directories via the find utility.
             :literal:`find . -type d -execdir chmod u=rwx,g=rx,o= {} +`
          4. Reset the file permissions for all regular files via the find utility.
             :literal:`find . -type f -execdir chmod u=rw,g=r,o= {} +`
          5. Compress all files in the :literal:`IRF/Logs/` directory.
             :literal:`gzip --best --recursive Logs/`
          6. Compress all files in the :literal:`IRF/Data/RTALogs/` directory if it exists.
             :literal:`gzip --best --recursive IRF/Data/RTALogs/`
          7. Compress all files in the :literal:`IRF/RTALogs/` directory if it exists.
             :literal:`gzip --best --recursive IRF/RTALogs/`
          8. Create the archive directory.

        - Lane specific:

          1. Compress all :literal:`IRF/Data/Intensities/BaseCalls/L00[1-8]/C1.1/*.bcl` files.
             :literal:`find IRF/Data/Intensities/BaseCalls/L00[1-8]
             -name '*.bcl' -execdir gzip --best --verbose {} +`
          2. Run the GNU Tar utility over each :literal:`IRF/Data/Intensities/L00[1-8]/` directory,
             but exclude compressed cluster locations (:literal:`*.clocs`) files.
          3. Calculate a MD5 checksum.

        - :emphasis:`Illumina Run Folder`-specific:

          1. Run the GNU Tar utility over the remaining :emphasis:`Illumina Run Folder`,
             but exclude directories with cluster intensity (:literal:`*.cif`) files.
          2. Calculate a MD5 checksum.
        """
        # Define an Illumina Run Folder directory.
        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        if not self.run_directory:
            raise Exception(f"A {self.name!s} requires a 'run_directory' configuration option.")

        self.run_directory = self.configuration.get_absolute_path(
            file_path=self.run_directory,
            default_path=StandardFilePath.get_illumina_run(absolute=True))

        # The Illumina Run Folder name would also be available from the bsf.illumina.RunFolder.get_name property below,
        # but using the file path provides more flexibility aside the fixed runParameters.xml configuration file.

        run_name = os.path.basename(self.run_directory)

        # Check that the Illumina Run Folder exists.

        if not os.path.isdir(self.run_directory):
            raise Exception(f"The Illumina 'run_directory' {self.run_directory!r} does not exist.")

        # Check whether the Illumina Run Folder is complete.
        # Check whether the IRF/RTAComplete.txt file exists in the Illumina Run Folder
        # to prevent archiving and deleting of an incomplete folder.
        # Alternatively, require the force instance variable to start archiving.

        if not (os.path.exists(os.path.join(self.run_directory, 'RTAComplete.txt')) or self.force):
            raise RunFolderNotComplete(f"The Illumina 'run_directory' {self.run_directory!r} is not complete.")

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
            raise Exception(f"The directory above the 'archive_directory' "
                            f"{os.path.dirname(self.archive_directory)!r} does not exist.")

        irf = RunFolder.from_file_path(file_path=self.run_directory)

        # The experiment name (e.g., BSF_0000) is used as part of the project_name.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if not self.experiment_name:
            self.experiment_name = irf.run_parameters.get_experiment_name
            if not self.experiment_name:
                raise Exception(f"A {self.name!s} requires an 'experiment_name' configuration option.")

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = '_'.join((self.experiment_name, irf.run_information.flow_cell))

        # Check that the cloud concurrency value (maximum number of network connections) is an integer.

        if not self.cloud_concurrency:
            self.cloud_concurrency = 1

        super(IlluminaRunFolderArchive, self).run()

        stage_pre_process_folder = self.get_stage(name=self.get_stage_name_pre_process())
        stage_compress_base_calls = self.get_stage(name=self.get_stage_name_base_calls())
        stage_archive_intensities = self.get_stage(name=self.get_stage_name_intensities())
        stage_archive_folder = self.get_stage(name=self.get_stage_name_archive_folder())
        stage_post_process_folder = self.get_stage(name=self.get_stage_name_post_process())
        stage_sav = self.get_stage(name=self.get_stage_name_sav())
        stage_cloud = self.get_stage(name=self.get_stage_name_cloud())

        # Pre-process on folder level.

        runnable_pre_process_folder = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_pre_process(project_name=self.project_name),
                working_directory=self.project_directory))
        executable_pre_process_folder = self.set_stage_runnable(
            stage=stage_pre_process_folder,
            runnable=runnable_pre_process_folder)

        # Check whether Picard ExtractIlluminaBarcodes has written any
        # s_<lane>_<tile>_barcode.txt(.gz) files into the BaseCalls directory.
        # Keeping them is rather pointless, and they should be removed.
        # http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes

        # Reset all file permissions on the Illumina Run Folder.

        runnable_step = RunnableStepChangeMode(
            name='reset_permissions',
            file_path=self.run_directory,
            mode_directory=self.irf_mode_directory,
            mode_file=self.irf_mode_file)
        runnable_pre_process_folder.add_runnable_step(runnable_step=runnable_step)

        # Compress all files in the IRF/Logs/ and IRF/Logs/IALogs/ directories.
        # The NextSeq compresses into a single Logs.zip file.

        if irf.run_parameters.get_instrument_type not in ('NextSeq',):
            runnable_step = RunnableStep(
                name='compress_logs',
                program='pigz')
            runnable_pre_process_folder.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_long(key='best')
            runnable_step.add_option_long(key='processes', value=str(stage_pre_process_folder.threads))
            runnable_step.add_switch_long(key='recursive')
            runnable_step.arguments.append(os.path.join(self.run_directory, 'Logs'))

        # Compress all files in the IRF/Data/RTALogs/ directory if it exists.
        # It does not on the HiSeq 3000/4000 and NovaSeq platforms.

        if os.path.isdir(os.path.join(self.run_directory, 'Data', 'RTALogs')):
            runnable_step = RunnableStep(
                name='compress_rta_logs',
                program='pigz')
            runnable_pre_process_folder.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_long(key='best')
            runnable_step.add_option_long(key='processes', value=str(stage_pre_process_folder.threads))
            runnable_step.add_switch_long(key='recursive')
            runnable_step.arguments.append(os.path.join(self.run_directory, 'Data', 'RTALogs'))

        # Compress all files in the IRF/RTALogs/ directory if it exists.
        # It only exists on the HiSeq 3000/4000 platform.

        if os.path.isdir(os.path.join(self.run_directory, 'RTALogs')):
            runnable_step = RunnableStep(
                name='compress_rta_logs',
                program='pigz')
            runnable_pre_process_folder.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_long(key='best')
            runnable_step.add_option_long(key='processes', value=str(stage_pre_process_folder.threads))
            runnable_step.add_switch_long(key='recursive')
            runnable_step.arguments.append(os.path.join(self.run_directory, 'RTALogs'))

        # Create the archive directory.

        runnable_step = RunnableStepMakeDirectory(
            name='make_directory',
            directory_path=self.archive_directory)
        runnable_pre_process_folder.add_runnable_step(runnable_step=runnable_step)

        # Process per lane.

        # Cluster intensity file (*.cif) directories, if present, need excluding from archiving at a later stage.

        exclude_intensities_patterns: list[str] = list()
        archive_folder_dependencies: list[str] = list()

        archive_folder_dependencies.append(runnable_pre_process_folder.name)

        for lane_int in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):
            if not irf.has_compressed_base_calls():
                # Process the IRF/Data/Intensities/BaseCalls/ directory.

                runnable_base_calls = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_base_calls(project_name=self.project_name, lane=str(lane_int)),
                        working_directory=self.project_directory))
                executable_base_calls = self.set_stage_runnable(
                    stage=stage_compress_base_calls,
                    runnable=runnable_base_calls)

                # Set a dependency on the executable_pre_process_folder
                executable_base_calls.dependencies.append(executable_pre_process_folder.name)

                # Compress all base call (*.bcl) files.

                runnable_step = RunnableStep(
                    name='compress_base_calls',
                    program='find',
                    sub_command=Command(
                        program=os.path.join(
                            self.run_directory, 'Data', 'Intensities', 'BaseCalls',
                            'L{:03d}'.format(lane_int)),
                        sub_command=Command(
                            program='-execdir',
                            sub_command=Command(
                                program='pigz'))))
                runnable_base_calls.add_runnable_step(runnable_step=runnable_step)

                find_command = runnable_step.sub_command  # directory option
                find_command.add_option_short(key='name', value='*.bcl')
                exec_dir_command = find_command.sub_command  # -execdir option
                pigz_command = exec_dir_command.sub_command  # pigz command
                pigz_command.add_switch_long(key='best')
                pigz_command.add_option_long(key='processes', value=str(stage_compress_base_calls.threads))
                pigz_command.arguments.append('{}')
                pigz_command.arguments.append('+')

                # Record dependencies for the archive run folder analysis stage.
                archive_folder_dependencies.append(executable_base_calls.name)

            if irf.has_intensities() and os.path.exists(
                    os.path.join(self.run_directory, 'Data', 'Intensities', 'L{:03d}'.format(lane_int))):

                # Process IRF/Data/Intensities/L00[1-8]/ directories if they exist.

                runnable_intensities = self.add_runnable(
                    runnable=ConcurrentRunnable(
                        name=self.get_prefix_intensities(project_name=self.project_name, lane=str(lane_int)),
                        working_directory=self.project_directory))
                executable_intensities = self.set_stage_runnable(
                    stage=stage_archive_intensities,
                    runnable=runnable_intensities)
                executable_intensities.dependencies.append(executable_pre_process_folder.name)

                # Run GNU Tar over the IRF/Data/Intensities/L00[1-8]/ directories, but exclude
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

                runnable_step = RunnableStep(
                    name='archive_intensities_tar',
                    program='tar',
                    stdout=ConnectorPipe())
                runnable_intensities.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_switch_long(key='create')
                runnable_step.add_option_long(key='directory', value=os.path.dirname(self.run_directory))
                runnable_step.add_option_long(key='exclude', value='*.clocs', override=True)
                runnable_step.add_option_long(key='exclude', value='*.locs', override=True)

                # Archiving needs the relative path.
                runnable_step.arguments.append(
                    os.path.join(
                        os.path.basename(self.run_directory),
                        'Data',
                        'Intensities',
                        'L{:03d}'.format(lane_int)))

                runnable_step = RunnableStep(
                    name='archive_intensities_pigz',
                    program='pigz',
                    stdin=ConcurrentProcess(name='archive_intensities_tar', connection='stdout'),
                    stdout=ConnectorFile(file_path=archive_file_path, file_mode='wb'))
                runnable_intensities.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_switch_long(key='best')
                runnable_step.add_switch_long(key='stdout')
                runnable_step.add_option_long(key='processes', value=str(stage_archive_intensities.threads))

                # Record dependencies for the archive run folder analysis stage.
                # Since cluster intensity files are no longer automatically deleted,
                # but just excluded from archiving the folder, those dependencies are no longer required.
                # archive_folder_dependencies.append(executable_intensities.name)

                # Record an exclude-pattern with the relative path to lane and cycle-specific
                # cluster intensity file (*.cif) directories. On the HiSeq 2000 platform, the lane directories
                # contain cluster locations (*.clocs) files that are essential in base call extraction and
                # need archiving with the IRF/Data/Intensities/BaseCalls folder.
                # By default, GNU Tar treats exclusion members as globbing patterns.
                exclude_intensities_patterns.append(
                    os.path.join(
                        os.path.basename(self.run_directory),
                        'Data',
                        'Intensities',
                        'L{:03d}'.format(lane_int),
                        'C*'))

                # Calculate a MD5 checksum.

                runnable_step = RunnableStep(
                    name='md5sum',
                    program='md5sum',
                    stdout=ConnectorFile(file_path=archive_file_path + '.md5', file_mode='wt'))
                runnable_intensities.add_runnable_step_epilogue(runnable_step=runnable_step)

                runnable_step.add_switch_long(key='binary')

                runnable_step.arguments.append(archive_file_path)

        # Process the whole run folder.

        runnable_archive_folder = self.add_runnable(
            runnable=ConcurrentRunnable(
                name=self.get_prefix_archive_folder(project_name=self.project_name),
                working_directory=self.project_directory))
        executable_archive_folder = self.set_stage_runnable(
            stage=stage_archive_folder,
            runnable=runnable_archive_folder)
        executable_archive_folder.dependencies.extend(archive_folder_dependencies)

        # Run the GNU Tar utility over the remaining Illumina Run folder,
        # but exclude directories with cluster intensity (*.cif) files.

        archive_file_path = run_name
        if self.compress_archive_files:
            archive_file_path = '.'.join((archive_file_path, 'tar', 'gz'))
        else:
            archive_file_path = '.'.join((archive_file_path, 'tar'))
        archive_file_path = os.path.join(self.archive_directory, archive_file_path)

        runnable_step = RunnableStep(
            name='archive_folder_tar',
            program='tar',
            stdout=ConnectorPipe())
        runnable_archive_folder.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='create')
        runnable_step.add_option_long(key='directory', value=os.path.dirname(self.run_directory))
        for pattern in exclude_intensities_patterns:
            runnable_step.add_option_long(key='exclude', value=pattern, override=True)

        runnable_step.arguments.append(os.path.basename(self.run_directory))

        runnable_step = RunnableStep(
            name='archive_folder_pigz',
            program='pigz',
            stdin=ConcurrentProcess(name='archive_folder_tar', connection='stdout'),
            stdout=ConnectorFile(file_path=archive_file_path, file_mode='wb'))
        runnable_archive_folder.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='best')
        runnable_step.add_switch_long(key='stdout')
        runnable_step.add_option_long(key='processes', value=str(stage_archive_folder.threads))

        # Post-process the archive folder.

        runnable_post_process_folder = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_post_process(project_name=self.project_name),
                working_directory=self.project_directory))
        executable_post_process_folder = self.set_stage_runnable(
            stage=stage_post_process_folder,
            runnable=runnable_post_process_folder)
        executable_post_process_folder.dependencies.append(runnable_archive_folder.name)

        # Calculate a MD5 checksum.

        runnable_step = RunnableStep(
            name='md5_sum',
            program='md5sum',
            stdout=ConnectorFile(file_path=archive_file_path + '.md5', file_mode='wt'))
        runnable_post_process_folder.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='binary')

        runnable_step.arguments.append(archive_file_path)

        # Extract files relevant for the Illumina Sequence Analysis Viewer (SAV).

        runnable_sav = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_sav(project_name=self.project_name),
                working_directory=self.project_directory))
        executable_sav = self.set_stage_runnable(
            stage=stage_sav,
            runnable=runnable_sav)
        executable_sav.dependencies.append(runnable_post_process_folder.name)

        runnable_step = RunnableStep(
            name='extract_folder_tar',
            program='tar')
        runnable_sav.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='extract')
        runnable_step.add_option_long(key='file', value=archive_file_path)
        runnable_step.add_switch_long(key='ignore-failed-read')
        runnable_step.add_switch_long(key='wildcards')

        for file_path in irf.get_sav_archive_paths():
            runnable_step.arguments.append(file_path)

        # Move the extracted Run Folder into the Illumina Sequence Analysis Viewer (SAV) directory.

        sav_path = Configuration.get_absolute_path(
            file_path=run_name + '_sav',
            default_path=StandardFilePath.get_illumina_sav(absolute=True))

        runnable_step = RunnableStepMove(
            name='move_sav',
            source_path=run_name,
            target_path=sav_path)
        runnable_sav.add_runnable_step(runnable_step=runnable_step)

        # Adjust the file permissions of the SAV folder.

        runnable_step = RunnableStepChangeMode(
            name='set_sav_permissions',
            file_path=sav_path,
            mode_directory=self.sav_mode_directory,
            mode_file=self.sav_mode_file)
        runnable_sav.add_runnable_step(runnable_step=runnable_step)

        # Upload the archive run folder into the block blob storage.

        if self.cloud_account and self.cloud_container:
            runnable_cloud = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_cloud(project_name=self.project_name),
                    working_directory=self.project_directory))
            executable_cloud = self.set_stage_runnable(
                stage=stage_cloud,
                runnable=runnable_cloud)
            executable_cloud.dependencies.append(runnable_post_process_folder.name)

            # The target (blob) path is the base name and optionally the cloud path prefix.
            target_path = os.path.basename(archive_file_path)
            if self.cloud_path_prefix:
                # The Azure Storage Blob Service always uses URL-compliant slash characters as path separators.
                target_path = '/'.join((self.cloud_path_prefix, target_path))

            # Upload the GNU Tar file.

            runnable_step = RunnableStepAzureBlockBlobUpload(
                name='blob_upload_tar',
                account_name=self.cloud_account,
                container_name=self.cloud_container,
                source_path=archive_file_path,
                target_path=target_path,
                max_concurrency=self.cloud_concurrency)
            runnable_cloud.add_runnable_step(runnable_step=runnable_step)

            # Upload the MD5 checksum file.

            runnable_step = RunnableStepAzureBlockBlobUpload(
                name='blob_upload_md5',
                account_name=self.cloud_account,
                container_name=self.cloud_container,
                source_path=archive_file_path + '.md5',
                target_path=target_path + '.md5',
                max_concurrency=self.cloud_concurrency)
            runnable_cloud.add_runnable_step(runnable_step=runnable_step)

        return

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            project_name: Optional[str] = None,
            archive_directory: Optional[str] = None,
            irf_path: Optional[str] = None,
            force: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a
        :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive` analysis.

        This analysis requires either a :literal:`configuration_path` argument or a :literal:`irf_path` argument.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param project_name: A project name.
        :type project_name: str | None
        :param archive_directory: An archive directory.
        :type archive_directory: str | None
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

        analysis: IlluminaRunFolderArchive = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (INI) file.')

            analysis.project_name = project_name

        if archive_directory:
            analysis.archive_directory = archive_directory

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
        print('Archive directory:      ', analysis.archive_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a
        :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderArchive` analysis.

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
            '--project-name',
            help='project name (i.e., instrument run identifier)')

        argument_parser.add_argument(
            '--archive-directory',
            help='archive directory')

        argument_parser.add_argument(
            '--irf',
            help='Illumina Run Folder name or file path',
            dest='irf_path')

        argument_parser.add_argument(
            '--force',
            action='store_true',
            help='force processing even if a run folder exists already')

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
            archive_directory=name_space.archive_directory,
            irf_path=name_space.irf_path,
            force=name_space.force)


class FilePathIlluminaRunFolderRestore(FilePath):
    """The :py:class:`bsf.analyses.illumina_run_folder.FilePathIlluminaRunFolderRestore` class models files in an
    archive directory.

    :ivar folder: GNU Tar archive file of the Illumina Run Folder.
    :type folder: str
    :ivar intensities: GNU Tar archive file of the intensities.
    :type intensities: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.illumina_run_folder.FilePathIlluminaRunFolderRestore` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathIlluminaRunFolderRestore, self).__init__(prefix=prefix)

        self.folder = '_'.join((prefix, 'Folder.tar'))
        self.intensities = '_'.join((prefix, 'Intensities.tar'))
        # There are additional attributes L001 to L008 depending on the flow cell layout.

        return


class IlluminaRunFolderRestore(Analysis):
    """The :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore` class represents the logic to
    restore an :emphasis:`Illumina Run Folder` from a format suitable for magnetic tape libraries.

    :cvar maximum_lane_number: A maximum number of lanes.
    :type maximum_lane_number: int
    :ivar archive_directory: An archive directory path.
    :type archive_directory: str | None
    :ivar illumina_directory: A directory of :emphasis:`Illumina Run Folder` directories path.
    :type illumina_directory: str | None
    :ivar extract_intensities: Request extracting cluster intensity file (:literal:`*.cif`) directories.
    :type extract_intensities: bool | None
    :ivar force: Request processing of incomplete :emphasis:`Illumina Run Folder` objects.
    :type force: bool | None
    """

    name = 'Illumina Run Folder Restore Analysis'
    prefix = 'irf_restore'

    maximum_lane_number = 8

    @classmethod
    def get_stage_name_extract_archive(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'extract_archive'))

    @classmethod
    def get_stage_name_compress_base_calls(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'compress_base_calls'))

    @classmethod
    def get_stage_name_compress_logs(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'compress_logs'))

    @classmethod
    def get_prefix_extract_archive(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_extract_archive(), project_name, lane))

    @classmethod
    def get_prefix_compress_base_calls(cls, project_name: str, lane: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :param lane: A lane number.
        :type lane: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_compress_base_calls(), project_name, lane))

    @classmethod
    def get_prefix_compress_logs(cls, project_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param project_name: A project name.
        :type project_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_compress_logs(), project_name))

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
            archive_directory: Optional[str] = None,
            illumina_directory: Optional[str] = None,
            extract_intensities: Optional[bool] = False,
            force: Optional[bool] = False):
        """Initialise a :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore` object.

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
        :param archive_directory: An archive directory path.
        :type archive_directory: str | None
        :param illumina_directory: A directory of :emphasis:`Illumina Run Folder` directories path.
        :type illumina_directory: str | None
        :param extract_intensities: Request extracting cluster intensity file (:literal:`*.cif`) directories.
        :type extract_intensities: bool | None
        :param force: Request processing of incomplete :emphasis:`Illumina Run Folder` objects.
        :type force: bool | None
        """
        super(IlluminaRunFolderRestore, self).__init__(
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

        # Sub-class specific ...

        self.archive_directory = archive_directory
        self.illumina_directory = illumina_directory
        self.extract_intensities = extract_intensities
        self.force = force

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
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

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore` object.

        Restore an Illumina Run Folder from a format suitable for magnetic tape libraries.
            1. Extract the :literal:`IRF_Folder.tar` file.
            2. Extract the :literal:`IRF_Intensities.tar` file with a 60 seconds delay.
            3. Extract each :literal:`IRF_L00[1-8].tar` file with a 90 seconds delay.
        """
        if not self.archive_directory:
            raise Exception(f"A {self.name!s} requires an 'archive_directory' configuration option.")

        self.archive_directory = self.configuration.get_absolute_path(
            file_path=self.archive_directory)

        if not os.path.isdir(self.archive_directory):
            raise Exception(f"The 'archive_directory' {self.archive_directory!r} does not exist.")

        self.illumina_directory = self.configuration.get_absolute_path(
            file_path=self.illumina_directory,
            default_path=StandardFilePath.get_illumina_run(absolute=True))

        if not os.path.isdir(self.illumina_directory):
            raise Exception(f"The 'illumina_directory' {self.illumina_directory!r} does not exist.")

        # The Illumina Run Folder must *not* already exist unless the force option has been set.

        run_directory_name = os.path.basename(self.archive_directory)[:-8]
        run_directory_path = os.path.join(self.illumina_directory, run_directory_name)

        if os.path.exists(run_directory_path) and not self.force:
            raise Exception(f"The Illumina Run Folder directory {run_directory_path!r} does exist already.")

        # A FilePath object cannot be used, because instruments have a variable number of lanes.
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
            raise Exception(f"The Illumina Run Archive file {file_path_dict['folder']!r} is missing.")

        if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['intensities'])):
            raise Exception(f"The Illumina Run Archive file {file_path_dict['intensities']!r} is missing.")

        if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['L001'])):
            raise Exception(f"The Illumina Run Archive file {file_path_dict['L001']!r} is missing.")

        super(IlluminaRunFolderRestore, self).run()

        stage_extract_archive = self.get_stage(name=self.get_stage_name_extract_archive())
        stage_compress_base_calls = self.get_stage(name=self.get_stage_name_compress_base_calls())
        stage_compress_logs = self.get_stage(name=self.get_stage_name_compress_logs())

        # Extract the IRF_Folder.tar file.

        runnable_extract_folder = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_extract_archive(project_name=self.project_name, lane='folder'),
                working_directory=self.project_directory))
        executable_extract_folder = self.set_stage_runnable(
            stage=stage_extract_archive,
            runnable=runnable_extract_folder)

        runnable_step = RunnableStep(
            name='extract_folder',
            program='tar')
        runnable_extract_folder.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='extract')
        runnable_step.add_option_long(key='directory', value=self.illumina_directory)
        runnable_step.add_option_long(
            key='file',
            value=os.path.join(self.archive_directory, file_path_dict['folder']))

        # Compress all files in the IRF/Logs and IRF/Data/RTALogs directories.

        runnable_compress_logs = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_compress_logs(project_name=self.project_name),
                working_directory=self.project_directory))
        executable_compress_logs = self.set_stage_runnable(
            stage=stage_compress_logs,
            runnable=runnable_compress_logs)
        executable_compress_logs.dependencies.append(executable_extract_folder.name)

        runnable_step = RunnableStep(
            name='compress_logs',
            program='gzip')
        runnable_compress_logs.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='best')
        runnable_step.add_switch_long(key='recursive')
        runnable_step.arguments.append(os.path.join(run_directory_path, 'Logs'))

        runnable_step = RunnableStep(
            name='compress_rta_logs',
            program='gzip')
        runnable_compress_logs.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='best')
        runnable_step.add_switch_long(key='recursive')
        runnable_step.arguments.append(os.path.join(run_directory_path, 'Data', 'RTALogs'))

        # Extract the IRF_intensities.tar file.

        runnable_extract_intensities = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_extract_archive(project_name=self.project_name, lane='intensities'),
                working_directory=self.project_directory))
        executable_extract_intensities = self.set_stage_runnable(
            stage=stage_extract_archive,
            runnable=runnable_extract_intensities)

        # Sleep 60 seconds to allow the first process to create all directories.
        runnable_step = RunnableStepSleep(
            name='sleep',
            sleep_time=60.0)
        runnable_extract_intensities.add_runnable_step(runnable_step=runnable_step)

        runnable_step = RunnableStep(
            name='extract_intensities',
            program='tar')
        runnable_extract_intensities.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='extract')
        runnable_step.add_option_long(key='directory', value=self.illumina_directory)
        runnable_step.add_option_long(
            key='file',
            value=os.path.join(self.archive_directory, file_path_dict['intensities']))

        # Unpack the IRF_L001.tar to IRF_L008.tar files if they exist.

        for lane_int in range(0 + 1, self.maximum_lane_number + 1):

            # HiSeq and MiSeq instruments offer various run modes with differing number of lanes.
            if not os.path.exists(os.path.join(self.archive_directory, file_path_dict['L{:03d}'.format(lane_int)])):
                continue

            runnable_extract_lane = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_extract_archive(project_name=self.project_name, lane=str(lane_int)),
                    working_directory=self.project_directory))
            self.set_stage_runnable(
                stage=stage_extract_archive,
                runnable=runnable_extract_lane)

            # Sleep 90 seconds to allow the first process to create all directories.
            runnable_step = RunnableStepSleep(
                name='sleep',
                sleep_time=90.0)
            runnable_extract_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step = RunnableStep(
                name='extract_lane',
                program='tar')
            runnable_extract_lane.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_long(key='extract')
            runnable_step.add_option_long(key='directory', value=self.illumina_directory)
            runnable_step.add_option_long(
                key='file',
                value=os.path.join(self.archive_directory, file_path_dict['L{:03d}'.format(lane_int)]))

            if not self.extract_intensities:
                runnable_step.add_option_long(key='exclude', value='C*.1', override=True)

            # Create one process per lane to compress the base call (*.bcl) files.

            runnable_compress_base_calls = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_compress_base_calls(project_name=self.project_name, lane=str(lane_int)),
                    working_directory=self.project_directory))
            executable_compress_base_calls = self.set_stage_runnable(
                stage=stage_compress_base_calls,
                runnable=runnable_compress_base_calls)
            executable_compress_base_calls.dependencies.append(executable_extract_intensities.name)

            # Since find has a somewhat broken syntax, definition of the bsf.process.RunnableStep gets a bit complex.
            # The following command line needs to be run
            # find IRF/Data/Intensities/BaseCalls/L000 -name '*.bcl' -execdir gzip --best {} +
            runnable_step = RunnableStep(
                name='compress_base_calls',
                program='find',
                sub_command=Command(
                    program=os.path.join(
                        run_directory_path, 'Data', 'Intensities', 'BaseCalls',
                        'L{:03d}'.format(lane_int)),
                    sub_command=Command(
                        program='-execdir',
                        sub_command=Command(
                            program='gzip'))))
            runnable_compress_base_calls.add_runnable_step(runnable_step=runnable_step)

            find_command = runnable_step.sub_command  # directory option
            find_command.add_option_short(key='name', value='*.bcl')
            exec_dir_command = find_command.sub_command  # -execdir option
            gzip_command = exec_dir_command.sub_command  # gzip command
            gzip_command.add_switch_long(key='best')
            gzip_command.arguments.append('{}')
            gzip_command.arguments.append('+')

        return

    @classmethod
    def console_submit(
            cls,
            configuration_path: str,
            stage_name: Optional[str] = None,
            drms_submit: Optional[bool] = None,
            project_name: Optional[str] = None,
            archive_directory: Optional[str] = None,
            extract_intensities: Optional[bool] = None,
            force: Optional[bool] = None,
            *args, **kwargs) -> int:
        """Console function to submit a
        :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore` analysis.

        :param configuration_path: A configuration `INI <https://en.wikipedia.org/wiki/INI_file>`_ file path.
        :type configuration_path: str
        :param stage_name: A :py:class:`bsf.analysis.Stage` name.
        :type stage_name: str | None
        :param drms_submit: Request submitting into the DRMS.
        :type drms_submit: bool | None
        :param project_name: A project name.
        :type project_name: str | None
        :param archive_directory: An archive directory.
        :type archive_directory: str | None
        :param extract_intensities: Request extracting intensities of an :emphasis:`Illumina Run Folder`.
        :type extract_intensities: bool | None
        :param force: Request processing incomplete :emphasis:`Illumina Run Folder` objects.
        :type force: bool | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        analysis: IlluminaRunFolderRestore = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (INI) file.')

            analysis.project_name = project_name

        if archive_directory:
            analysis.archive_directory = archive_directory

        if extract_intensities:
            analysis.extract_intensities = extract_intensities

        if force:
            analysis.force = force

        analysis.run()
        analysis.check_state()
        analysis.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis.name)
        print('Project name:           ', analysis.project_name)
        print('Project directory:      ', analysis.project_directory)
        print('Illumina run directory: ', analysis.illumina_directory)
        print('Archive directory:      ', analysis.archive_directory)

        return 0

    @classmethod
    def entry_point_submit(cls) -> int:
        """Console entry point to submit a
        :py:class:`bsf.analyses.illumina_run_folder.IlluminaRunFolderRestore` analysis.

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
            required=True,
            help='project name i.e. instrument run identifier')

        argument_parser.add_argument(
            '--archive-directory',
            required=True,
            help='archive directory')

        argument_parser.add_argument(
            '--extract-intensities',
            action='store_true',
            help='extract cluster intensity (CIF) files')

        argument_parser.add_argument(
            '--force',
            action='store_true',
            help='force processing even if a run folder exists already')

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
            archive_directory=name_space.archive_directory,
            extract_intensities=name_space.extract_intensities,
            force=name_space.force)
