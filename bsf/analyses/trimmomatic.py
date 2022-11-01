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
"""The :py:mod`bsf.analyses.trimmomatic` module classes supporting the
`Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ tool.
"""
import logging
import os
from typing import List

from bsf.analysis import Analysis, Stage
from bsf.executables.collection import RunnableStepCollectionPruneFastq
from bsf.ngs import Collection, Sample, PairedReads, Reads
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, RunnableStep, RunnableStepJava, RunnableStepMakeDirectory
from bsf.standards import Configuration, JavaArchive

module_logger = logging.getLogger(name=__name__)


class FilePathTrimmomaticReadGroup(FilePath):
    """The :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticReadGroup` class models
    read group-specific Trimmomatic files.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar trim_log_tsv: A Trimmomatic trim log Tab-Separated Value (TSV) file path.
    :type trim_log_tsv: str
    :ivar summary_tsv: A summary Tab-Separated Value (TSV) file path.
    :type summary_tsv: str
    :ivar coverage_png: A coverage Portable Network Graphics (PNG) file path.
    :type coverage_png: str
    :ivar frequency_png: A frequency Portable Network Graphics (PNG) file path.
    :type frequency_png: str
    :ivar surviving_png: A surviving Portable Network Graphics (PNG) file path.
    :type surviving_png: str
    :ivar reads_1p: A First Reads paired file path.
    :type reads_1p: str | None
    :ivar reads_1u: First Reads unpaired file path.
    :type reads_1u: str | None
    :ivar reads_2p: Second Reads paired file path.
    :type reads_2p: str | None
    :ivar reads_2u: Second Reads unpaired file path.
    :type reads_2u: str | None
    """

    def __init__(self, prefix):
        """Initialise a :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticReadGroup` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathTrimmomaticReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix
        # Automatic GNU Zip-compression of trim log files does not work.
        self.trim_log_tsv = prefix + '_trim_log.tsv'
        self.trim_log_tsv_gz = prefix + '_trim_log.tsv.gz'
        self.summary_tsv = prefix + '_summary.tsv'  # Defined by the R script.
        self.coverage_png = prefix + '_coverage.png'  # Defined by the R script.
        self.frequency_png = prefix + '_frequency.png'  # Defined by the R script.
        self.surviving_png = prefix + '_surviving.png'  # Defined by the R script.
        self.reads_1p = None
        self.reads_1u = None
        self.reads_2p = None
        self.reads_2u = None

        return


class FilePathTrimmomaticProject(FilePath):
    """The :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticProject` class models
    project-specific Trimmomatic files.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar sas_path_old: An old Sample Annotation Sheet file path.
    :type sas_path_old: str
    :ivar sas_path_new: A new Sample Annotation Sheet file path.
    :type sas_path_new: str
    """

    def __init__(self, prefix, prefix_analysis, project_name):
        """Initialise a :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticProject` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        :param prefix_analysis: A :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix_analysis: str
        :param project_name: A project name.
        :type project_name: str
        """
        super(FilePathTrimmomaticProject, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sas_path_old = '_'.join((project_name, prefix_analysis, 'original.csv'))
        self.sas_path_new = '_'.join((project_name, prefix_analysis, 'samples.csv'))

        return


class Trimmomatic(Analysis):
    """The :py:class:`bsf.analyses.trimmomatic.Trimmomatic` class represents the logic to run the Trimmomatic analysis.

    :ivar adapter_path: An adapter definition file path.
    :type adapter_path: str | None
    :ivar trimming_step_pe_list: A colon-separated list of Trimmomatic steps for paired-end data.
    :type trimming_step_pe_list: list[str] | None
    :ivar trimming_step_se_list: A colon-separated list of Trimmomatic steps for single-end data.
    :type trimming_step_se_list: list[str] | None
    :ivar java_archive_trimmomatic: A Trimmomatic tool Java Archive (JAR) file path.
    :type java_archive_trimmomatic: str | None
    """

    name = 'Trimmomatic Analysis'
    prefix = 'trimmomatic'

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_stage_name_summary(cls):
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'summary'))

    @classmethod
    def get_stage_name_project(cls):
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'project'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_prefix_summary(cls):
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return cls.get_stage_name_summary()

    @classmethod
    def get_prefix_project(cls):
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return cls.get_stage_name_project()

    @classmethod
    def get_file_path_read_group(cls, read_group_name):
        """Get a :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticReadGroup` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticReadGroup` object.
        :rtype: FilePathTrimmomaticReadGroup
        """
        return FilePathTrimmomaticReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

    @classmethod
    def get_file_path_project(cls, project_name, prefix_analysis):
        """Get a :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticProject` object.

        :param project_name: A project name.
        :type project_name: str
        :param prefix_analysis: A :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix_analysis: str
        :return: A :py:class:`bsf.analyses.trimmomatic.FilePathTrimmomaticProject` object.
        :rtype: FilePathTrimmomaticProject
        """
        return FilePathTrimmomaticProject(
            prefix=cls.get_prefix_project(),
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
            report_style_path=None,
            report_header_path=None,
            report_footer_path=None,
            e_mail=None,
            stage_list=None,
            collection=None,
            sample_list=None,
            adapter_path=None,
            trimming_step_pe_list=None,
            trimming_step_se_list=None,
            java_archive_trimmomatic=None):
        """Initialise a :py:class:`bsf.analyses.trimmomatic.Trimmomatic` object.

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
        :param report_style_path: Report :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: Report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: Report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a UCSC Genome Browser Track Hub.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param adapter_path: An adapter definition file path.
        :type adapter_path: str | None
        :param trimming_step_pe_list: A colon-separated list of Trimmomatic steps for paired-end data.
        :type trimming_step_pe_list: list[str] | None
        :param trimming_step_se_list: A colon-separated list of Trimmomatic steps for single-end data.
        :type trimming_step_se_list: list[str] | None
        :param java_archive_trimmomatic: A Trimmomatic tool Java Archive (JAR) file path.
        :type java_archive_trimmomatic: str | None
        """
        super(Trimmomatic, self).__init__(
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

        self.adapter_path = adapter_path
        self.trimming_step_pe_list = trimming_step_pe_list
        self.trimming_step_se_list = trimming_step_se_list
        self.java_archive_trimmomatic = java_archive_trimmomatic

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a :py:class:`bsf.analyses.trimmomatic.Trimmomatic` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(Trimmomatic, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'adapter_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.adapter_path = configuration.config_parser.get(section=section, option=option)

        # Get the list of default trimming steps for paired-end data mode.

        option = 'trimming_steps_pe'
        if configuration.config_parser.has_option(section=section, option=option):
            self.trimming_step_pe_list = configuration.get_list_from_csv(section=section, option=option)

        # Get the list of default trimming steps for single-end data mode.

        option = 'trimming_steps_se'
        if configuration.config_parser.has_option(section=section, option=option):
            self.trimming_step_se_list = configuration.get_list_from_csv(section=section, option=option)

        # Get the Trimmomatic tool Java Archive (JAR) file path.

        option = 'java_archive_trimmomatic'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_trimmomatic = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a :py:class:`bsf.analyses.trimmomatic.Trimmomatic` object.

        This method changes the :py:class:`bsf.analyses.trimmomatic.Trimmomatic.collection` attribute
        (:py:class:`bsf.ngs.Collection`) to update with FASTQ file paths.
        """

        def run_adjust_illumina_clip_path(trimming_step_list):
            """Private function to adjust the adapter FASTA file path of :literal:`ILLUMINACLIP` trimming steps.

            If the file path is not absolute, prepend the value of the adapter_path
            instance variable.

            :param trimming_step_list: A Python :py:class:`list` object of
                Python :py:class:`str` (trimming step) objects.
            :type trimming_step_list: list[str]
            """
            for i in range(0, len(trimming_step_list)):
                if trimming_step_list[i].startswith('ILLUMINACLIP'):
                    component_list = trimming_step_list[i].split(':')
                    if not os.path.isabs(component_list[1]):
                        component_list[1] = os.path.join(self.adapter_path, component_list[1])
                    trimming_step_list[i] = ':'.join(component_list)
            return

        def run_read_comparisons():
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

        super(Trimmomatic, self).run()

        # Get the Trimmomatic tool Java Archive (JAR) file path.

        if not self.java_archive_trimmomatic:
            self.java_archive_trimmomatic = JavaArchive.get_trimmomatic()
            if not self.java_archive_trimmomatic:
                raise Exception(f"A {self.name!s} requires a 'java_archive_trimmomatic' configuration option.")

        if not (self.adapter_path and os.path.isabs(self.adapter_path)):
            self.adapter_path = os.path.join(os.path.dirname(self.java_archive_trimmomatic), 'adapters')

        if self.trimming_step_pe_list is None:
            raise Exception(f"A {self.name!s} requires a 'trimming_steps_pe' configuration option.")

        if self.trimming_step_se_list is None:
            raise Exception(f"A {self.name!s} requires a 'trimming_steps_se' configuration option.")

        run_adjust_illumina_clip_path(trimming_step_list=self.trimming_step_pe_list)
        run_adjust_illumina_clip_path(trimming_step_list=self.trimming_step_se_list)

        run_read_comparisons()

        # Trimmomatic

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())
        stage_summary = self.get_stage(name=self.get_stage_name_summary())
        stage_project = self.get_stage(name=self.get_stage_name_project())

        project_dependency_list: List[str] = list()

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            module_logger.debug('Sample name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            sample_step_list: List[str] = list()
            if 'Trimmomatic Steps' in sample.annotation_dict:
                for trimming_step in sample.annotation_dict['Trimmomatic Steps']:
                    sample_step_list.extend(self.configuration.list_from_csv(csv_string=trimming_step))
                run_adjust_illumina_clip_path(trimming_step_list=sample_step_list)

            # The Trimmomatic analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    module_logger.debug('PairedReads name: %r', paired_reads.get_name())
                    module_logger.log(logging.DEBUG - 2, 'PairedReads: %r', paired_reads)

                    paired_reads_step_list = list()
                    if 'Trimmomatic Steps' in paired_reads.annotation_dict:
                        for trimming_step in paired_reads.annotation_dict['Trimmomatic Steps']:
                            paired_reads_step_list.extend(self.configuration.list_from_csv(csv_string=trimming_step))
                        run_adjust_illumina_clip_path(trimming_step_list=paired_reads_step_list)

                    # Apply some sanity checks.

                    # Maybe this case should be allowed after Trimmomatic trimming,
                    # where only the second Read survives.
                    if paired_reads.reads_1 is None and paired_reads.reads_2 is not None:
                        raise Exception(f"The PairedReads object {paired_reads.get_name()!r} has a 'reads_2', "
                                        f"but no reads_1 object.")

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    # Create a Runnable and an Executable for running the Trimmomatic analysis.

                    runnable_read_group = self.add_runnable_consecutive(
                        runnable=ConsecutiveRunnable(
                            name=self.get_prefix_read_group(read_group_name=paired_reads_name),
                            working_directory=self.project_directory))
                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                    # Record the Executable.name for the project dependency.

                    project_dependency_list.append(runnable_read_group.name)

                    # Create a new RunnableStepMakeDirectory in preparation of the Trimmomatic program.

                    runnable_step_read_group = RunnableStepMakeDirectory(
                        name='mkdir',
                        directory_path=file_path_read_group.output_directory)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step_read_group)

                    # Create a RunnableStep for the Trimmomatic program.

                    runnable_step_read_group = RunnableStepJava(
                        name='trimmomatic',
                        java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx4G',
                        java_jar_path=self.java_archive_trimmomatic)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step_read_group)

                    if paired_reads.reads_2 is None:
                        runnable_step_read_group.sub_command.sub_command = Command(program='SE')
                    else:
                        runnable_step_read_group.sub_command.sub_command = Command(program='PE')

                    # Add options to the sub command.
                    sub_command = runnable_step_read_group.sub_command.sub_command
                    sub_command.add_option_short(key='trimlog', value=file_path_read_group.trim_log_tsv)

                    if paired_reads.reads_2 is None:
                        file_path_read_group.reads_1u = os.path.join(
                            file_path_read_group.output_directory,
                            paired_reads.reads_1.name + 'U.fastq.gz')

                        sub_command.arguments.append(paired_reads.reads_1.file_path)
                        sub_command.arguments.append(file_path_read_group.reads_1u)

                        # Update unpaired Reads information.

                        paired_reads.reads_1.name += 'U'
                        paired_reads.reads_1.file_path = os.path.join(
                            self.genome_directory,
                            file_path_read_group.reads_1u)
                    else:
                        file_path_read_group.reads_1p = os.path.join(
                            file_path_read_group.output_directory,
                            paired_reads.reads_1.name + 'P.fastq.gz')
                        file_path_read_group.reads_1u = os.path.join(
                            file_path_read_group.output_directory,
                            paired_reads.reads_1.name + 'U.fastq.gz')
                        file_path_read_group.reads_2p = os.path.join(
                            file_path_read_group.output_directory,
                            paired_reads.reads_2.name + 'P.fastq.gz')
                        file_path_read_group.reads_2u = os.path.join(
                            file_path_read_group.output_directory,
                            paired_reads.reads_2.name + 'U.fastq.gz')

                        sub_command.arguments.append(paired_reads.reads_1.file_path)
                        sub_command.arguments.append(paired_reads.reads_2.file_path)
                        sub_command.arguments.append(file_path_read_group.reads_1p)
                        sub_command.arguments.append(file_path_read_group.reads_1u)
                        sub_command.arguments.append(file_path_read_group.reads_2p)
                        sub_command.arguments.append(file_path_read_group.reads_2u)

                        # Update paired Reads information.

                        paired_reads.reads_1.name += 'P'
                        paired_reads.reads_1.file_path = os.path.join(
                            self.genome_directory,
                            file_path_read_group.reads_1p)
                        paired_reads.reads_2.name += 'P'
                        paired_reads.reads_2.file_path = os.path.join(
                            self.genome_directory,
                            file_path_read_group.reads_2p)

                        # Add unpaired Reads 1 and 2 as separate PairedReads objects to this sample.

                        sample.add_paired_reads(
                            paired_reads=PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=Reads(
                                    name=paired_reads.reads_1.name[:-1] + 'U',
                                    file_path=os.path.join(
                                        self.genome_directory,
                                        file_path_read_group.reads_1u)),
                                exclude=paired_reads.exclude,
                                index_1=paired_reads.index_1,
                                index_2=paired_reads.index_2,
                                read_group=paired_reads.read_group))

                        sample.add_paired_reads(
                            paired_reads=PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=Reads(
                                    name=paired_reads.reads_2.name[:-1] + 'U',
                                    file_path=os.path.join(
                                        self.genome_directory,
                                        file_path_read_group.reads_2u)),
                                exclude=paired_reads.exclude,
                                index_1=paired_reads.index_1,
                                index_2=paired_reads.index_2,
                                read_group=paired_reads.read_group))

                    # Append trimming steps in order of specificity read group, sample and default.

                    if len(paired_reads_step_list):
                        sub_command.arguments.extend(paired_reads_step_list)
                    elif len(sample_step_list):
                        sub_command.arguments.extend(sample_step_list)
                    else:
                        if paired_reads.reads_2 is None:
                            sub_command.arguments.extend(self.trimming_step_se_list)
                        else:
                            sub_command.arguments.extend(self.trimming_step_pe_list)

                    # Create a new RunnableStep to compress the trim log file.

                    runnable_step_read_group = RunnableStep(
                        name='compress_logs',
                        program='pigz')
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step_read_group)

                    runnable_step_read_group.add_switch_long(key='best')
                    runnable_step_read_group.add_option_long(key='processes', value=str(stage_summary.threads))
                    runnable_step_read_group.arguments.append(file_path_read_group.trim_log_tsv)

                    # Create a Runnable for the bsf_trimmomatic_summary.R analysis.

                    prefix_summary = '_'.join((stage_summary.name, paired_reads_name))

                    runnable_summary = self.add_runnable_consecutive(
                        runnable=ConsecutiveRunnable(
                            name=prefix_summary,
                            working_directory=self.project_directory))
                    executable_summary = self.set_stage_runnable(stage=stage_summary, runnable=runnable_summary)
                    executable_summary.dependencies.append(runnable_read_group.name)

                    # Create a new RunnableStep to aggregate the trim log file.

                    runnable_step_summary = RunnableStep(
                        name='trimmomatic_summary',
                        program='bsf_trimmomatic_summary.R',
                        obsolete_file_path_list=[
                            # file_path_read_group.trim_log_tsv_gz,
                        ])
                    runnable_summary.add_runnable_step(runnable_step=runnable_step_summary)

                    runnable_step_summary.add_option_long(
                        key='file-path',
                        value=file_path_read_group.trim_log_tsv_gz)

        # Create a Runnable for pruning the sample annotation sheet.

        prefix_project = '_'.join((stage_project.name, self.project_name))

        file_path_project = self.get_file_path_project(project_name=self.project_name, prefix_analysis=self.prefix)

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(self.project_directory, file_path_project.sas_path_old),
            name=prefix_project)
        annotation_sheet.to_file_path()

        runnable_project = self.add_runnable_consecutive(
            runnable=ConsecutiveRunnable(
                name=prefix_project,
                working_directory=self.project_directory))
        executable_project = self.set_stage_runnable(
            stage=stage_project,
            runnable=runnable_project)
        executable_project.dependencies.extend(project_dependency_list)

        # Create a new RunnableStep.

        runnable_step_project = RunnableStepCollectionPruneFastq(
            name='prune_sample_annotation_sheet',
            obsolete_file_path_list=[
                # file_path_project.sas_path_old,
            ],
            file_path_old=file_path_project.sas_path_old,
            file_path_new=file_path_project.sas_path_new,
            minimum_size=1024)
        runnable_project.add_runnable_step(runnable_step=runnable_step_project)

        return

    def report(self):
        """Create a :literal:`XHTML 1.0` report and a :literal:`UCSC Genome Browser Track Hub`.
        """

        # Create a symbolic link containing the project name and a UUID.
        self.create_public_project_link()

        # Write a HTML document.

        report_list: List[str] = list()

        report_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
        report_list.append('\n')

        report_list.append('<h2 id="aliquot_and_sample_level">Aliquot and Sample Level</h2>\n')
        report_list.append('\n')
        report_list.append('<table id="aliquot_and_sample_table">\n')
        report_list.append('<thead>\n')
        report_list.append('<tr>\n')
        report_list.append('<th>Sample</th>\n')
        report_list.append('<th>Aliquot</th>\n')
        report_list.append('<th>Coverage</th>\n')
        report_list.append('<th>Frequency</th>\n')
        report_list.append('<th>Summary</th>\n')
        report_list.append('</tr>\n')
        report_list.append('</thead>\n')
        report_list.append('<tbody>\n')

        for sample in self.sample_list:
            # The Trimmomatic analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=False)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            report_list.append('<tr>\n')
            report_list.append('<td class="left">' + sample.name + '</td>\n')
            report_list.append('<td class="left"></td>\n')  # Aliquot
            report_list.append('<td class="center"></td>\n')  # Coverage PNG
            report_list.append('<td class="center"></td>\n')  # Frequency PNG
            report_list.append('<td class="center"></td>\n')  # Summary TSV
            report_list.append('</tr>\n')

            # This analysis is special in that read group names carry 'P' or 'U' suffices and samples carry additional
            # read groups after trimming that do no longer correspond to the initial Runnable objects. Sigh.
            # Transiently create a Python dict without the suffix and sort its keys (bsf.ngs.PairedReads.name).

            for paired_reads_name in sorted(dict(map(lambda item: (item[:-1], True), paired_reads_dict.keys()))):
                prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                # The second read may still not be there.
                if prefix_read_group not in self.runnable_dict:
                    continue

                runnable_read_group = self.runnable_dict[prefix_read_group]
                file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                report_list.append('<tr>\n')
                # Sample
                report_list.append('<td class="left"></td>\n')
                # Aliquot
                report_list.append('<td class="left">' + paired_reads_name + '</td>\n')
                # Coverage
                report_list.append('<td class="center">')
                report_list.append('<a href="' + file_path_read_group.coverage_png + '">')
                report_list.append('<img alt="Coverage ' + runnable_read_group.name + '"')
                report_list.append(' src="' + file_path_read_group.coverage_png + '"')
                report_list.append(' height="100" width="100" />')
                report_list.append('</a>')
                report_list.append('</td>\n')
                # Frequency
                report_list.append('<td class="center">')
                report_list.append('<a href="' + file_path_read_group.frequency_png + '">PNG</a>')
                report_list.append('</td>\n')
                # The frequency plots provide little information that does not necessarily justify
                # adding another set of images onto the HTML report.
                report_list.append('<td class="center">')
                report_list.append('<a href="' + file_path_read_group.summary_tsv + '">TSV</a>')
                report_list.append('</td>\n')
                report_list.append('</tr>\n')

        report_list.append('</tbody>\n')
        report_list.append('</table>\n')
        report_list.append('\n')

        self.report_to_file(content=report_list)

        return
