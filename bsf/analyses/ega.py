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
"""The :py:mod:`bsf.analyses.ega` module supports the European Genome Phenome Archive (EGA) Cryptor tool.
"""
import logging
import os
from argparse import ArgumentParser
from typing import Optional

from bsf.analysis import Analysis, Stage
from bsf.ngs import Collection, Sample
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, RunnableStepJava, RunnableStepLink, RunnableStepMakeDirectory
from bsf.standards import Configuration

module_logger = logging.getLogger(name=__name__)


class FilePathEGACryptorReadGroup(FilePath):
    """The :py:class:`bsf.analyses.ega.FilePathEGACryptorReadGroup` models read group-specific EGA Cryptor file paths.
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.ega.FilePathEGACryptorReadGroup` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathEGACryptorReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class EGACryptor(Analysis):
    """The :py:class:`bsf.analyses.ega.EGACryptor` class models an EGA Cryptor.

    :ivar java_archive_ega_cryptor: An EGA Cryptor tool :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_ega_cryptor: str | None
    """

    name = 'EGA Cryptor Analysis'
    prefix = 'ega_cryptor'

    @classmethod
    def get_stage_name_read_group(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name: str) -> FilePathEGACryptorReadGroup:
        """Get a :py:class:`bsf.analyses.ega.FilePathEGACryptorReadGroup` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: :py:class:`bsf.analyses.ega.FilePathEGACryptorReadGroup` object.
        :rtype: FilePathEGACryptorReadGroup
        """
        return FilePathEGACryptorReadGroup(prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

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
            java_archive_ega_cryptor: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.ega.EGACryptor` object.

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
        :param report_style_path: A report style  :literal:`CSS` file path.
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
        :param java_archive_ega_cryptor: An EGA Cryptor tool :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_ega_cryptor: str | None
        """
        super(EGACryptor, self).__init__(
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

        self.java_archive_ega_cryptor = java_archive_ega_cryptor

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.ega.EGACryptor` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(EGACryptor, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the EGA Cryptor tool Java Archive (JAR) file path.

        option = 'java_archive_ega_cryptor'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_ega_cryptor = configuration.config_parser.get(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.ega.EGACryptor` object.
        """
        # Always check each BSF PairedReads object separately.
        replicate_grouping = False

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

        super(EGACryptor, self).run()

        run_read_comparisons()

        if not self.java_archive_ega_cryptor:
            raise Exception(f"A {self.name!s} requires a 'java_archive_ega_cryptor' configuration option.")

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            module_logger.debug('Sample name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            # The EGACryptor analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping, exclude=False)
            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    module_logger.debug('PairedReads name: %r', paired_reads.get_name())

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    # Create a Runnable and an Executable for running the EGA Cryptor analysis.

                    runnable_read_group = self.add_runnable(
                        runnable=ConsecutiveRunnable(
                            name=self.get_prefix_read_group(read_group_name=paired_reads_name),
                            working_directory=self.project_directory))
                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                    # Create a new RunnableStepMakeDirectory in preparation of the EGA Cryptor program.

                    runnable_step = RunnableStepMakeDirectory(
                        name='mkdir',
                        directory_path=file_path_read_group.output_directory)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    link_target_path = os.path.join(
                        file_path_read_group.output_directory,
                        os.path.basename(paired_reads.reads_1.file_path))

                    runnable_step = RunnableStepLink(
                        name='link',
                        source_path=paired_reads.reads_1.file_path,
                        target_path=link_target_path)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    # Create a RunnableStep to run the Java-based EGA Cryptor.

                    runnable_step = RunnableStepJava(
                        name='ega_cryptor',
                        java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx4G',
                        java_jar_path=self.java_archive_ega_cryptor)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    # Use a sequence of sub-Command objects to separate options that have to appear
                    # in a particular order. Sigh!

                    runnable_step.sub_command.sub_command = Command()
                    sub_command = runnable_step.sub_command.sub_command
                    sub_command.add_switch_short(key='p')

                    sub_command.sub_command = Command()
                    sub_command = sub_command.sub_command
                    sub_command.add_switch_short(key='fm')

                    sub_command.sub_command = Command()
                    sub_command = sub_command.sub_command
                    sub_command.add_switch_short(key='file')

                    sub_command.arguments.append(link_target_path)

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
        """Console function to submit a :py:class:`bsf.analyses.ega.EGACryptor` analysis.

        This analysis requires either a :literal:`configuration` argument or
        :literal:`project_name` and :literal:`sas_path` arguments.

        :param configuration_path: A UNIX configuration (*.ini) file path.
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

        analysis: EGACryptor = cls.from_config_file_path(config_path=configuration_path)

        if project_name:
            if project_name.endswith('.ini'):
                raise Exception('The --project-name option should not be a configuration (*.ini) file.')

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
        """Console entry point to submit a :py:class:`bsf.analyses.ega.EGACryptor` analysis.

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
            help=f'configuration (*.ini) file path [{Configuration.get_global_file_path()!s}]')

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
