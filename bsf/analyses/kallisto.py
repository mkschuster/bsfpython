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
"""The :py:mod:`bsf.analyses.kallisto` module provides classes and methods supporting
`Kallisto <https://pachterlab.github.io/kallisto/manual>`_ analyses.
"""
import logging
import os

from bsf.analysis import Analysis, Stage
from bsf.ngs import Collection, Sample
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, RunnableStep, RunnableStepMakeDirectory
from bsf.standards import Configuration, StandardFilePath, Transcriptome

module_logger = logging.getLogger(name=__name__)


class FilePathSample(FilePath):
    """The :py:class:`bsf.analyses.kallisto.FilePathSample` class models files in a
    sample-specific Kallisto directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    """

    def __init__(self, prefix):
        """Initialise a :py:class:`bsf.analyses.kallisto.FilePathSample` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathSample, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundance_h5 = os.path.join(prefix, 'abundance.h5')
        self.abundance_tsv = os.path.join(prefix, 'abundance.tsv')
        self.run_info_json = os.path.join(prefix, 'run_info.json')

        return


class Kallisto(Analysis):
    """The :py:class:`bsf.analyses.kallisto.Kallisto` class provides support for the Kallisto pseudo-aligner.

    :ivar transcriptome_version: A transcriptome version.
    :type transcriptome_version: str | None
    :ivar transcriptome_index_path: A transcriptome index file path.
    :type transcriptome_index_path: str | None
    :ivar fragment_length_value: A fragment length.
    :type fragment_length_value: str | None
    :ivar fragment_length_standard_deviation: A fragment length standard deviation.
    :type fragment_length_standard_deviation: str | None
    :ivar bias_correction: Request sequence-specific bias correction.
    :type bias_correction: bool
    :ivar bootstrap_samples: A number of bootstrap samples.
    :type bootstrap_samples: str
    """

    name = 'Kallisto Analysis'
    prefix = 'kallisto'

    @classmethod
    def get_stage_name_sample(cls):
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'sample'))

    @classmethod
    def get_prefix_sample(cls, sample_name):
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_sample(), sample_name))

    @classmethod
    def get_file_path_sample(cls, sample_name):
        """Get a :py:class:`bsf.analyses.kallisto.FilePathSample` object from this or a subclass.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.kallisto.FilePathSample` object or subclass thereof.
        :rtype: FilePathSample
        """
        return FilePathSample(prefix=cls.get_prefix_sample(sample_name=sample_name))

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
            transcriptome_version=None,
            transcriptome_index_path=None,
            fragment_length_value=None,
            fragment_length_standard_deviation=None,
            bias_correction=None,
            bootstrap_samples=None):
        """Initialise a :py:class:`bsf.analyses.kallisto.Kallisto` object.

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
        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str | None
        :param transcriptome_index_path: A transcriptome index file path.
        :type transcriptome_index_path: str | None
        :param fragment_length_value: A fragment length.
        :type fragment_length_value: str | None
        :param fragment_length_standard_deviation: A fragment length standard deviation.
        :type fragment_length_standard_deviation: str | None
        :param bias_correction: Request sequence-specific bias correction.
        :type bias_correction: bool
        :param bootstrap_samples: A number of bootstrap samples.
        :type bootstrap_samples: str
        """
        super(Kallisto, self).__init__(
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

        self.transcriptome_version = transcriptome_version
        self.transcriptome_index_path = transcriptome_index_path
        self.fragment_length_value = fragment_length_value
        self.fragment_length_standard_deviation = fragment_length_standard_deviation
        self.bias_correction = bias_correction
        self.bootstrap_samples = bootstrap_samples

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a :py:class:`bsf.analyses.kallisto.Kallisto` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(Kallisto, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'fragment_length_value'
        if configuration.config_parser.has_option(section=section, option=option):
            self.fragment_length_value = configuration.config_parser.get(section=section, option=option)

        option = 'fragment_length_standard_deviation'
        if configuration.config_parser.has_option(section=section, option=option):
            self.fragment_length_standard_deviation = configuration.config_parser.get(section=section, option=option)

        option = 'bias_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.bias_correction = configuration.config_parser.getboolean(section=section, option=option)

        option = 'bootstrap_samples'
        if configuration.config_parser.has_option(section=section, option=option):
            self.bootstrap_samples = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a :py:class:`bsf.analyses.kallisto.Kallisto` object.
        """

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

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Kallisto requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception(f"A {self.name!s} requires a valid 'genome_version' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the bsf.analysis.Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception(f'Sample annotation sheet {self.sas_file!r} does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[self.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation sheet in the current working directory.')

        super(Kallisto, self).run()

        if not self.transcriptome_index_path:
            self.transcriptome_index_path = os.path.join(
                StandardFilePath.get_resource_transcriptome_index(
                    transcriptome_version=self.transcriptome_version,
                    transcriptome_index='kallisto'),
                self.transcriptome_version + '.idx')

        if not os.path.exists(self.transcriptome_index_path):
            raise Exception(f'The transcriptome index file path {self.transcriptome_index_path!r} does not exist.')

        run_read_comparisons()

        stage_sample = self.get_stage(name=self.get_stage_name_sample())

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            module_logger.debug('Sample name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=True, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            prefix_sample = self.get_prefix_sample(sample_name=sample.name)

            file_path_sample = FilePathSample(prefix=prefix_sample)

            runnable_sample = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=prefix_sample,
                    working_directory=self.genome_directory))
            # executable_sample =
            self.set_stage_runnable(
                stage=stage_sample,
                runnable=runnable_sample)
            # Set dependencies on previous Runnable or bsf.process.Executable objects.
            # executable_sample.dependencies.append(executable_sample.name)

            runnable_step = RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=file_path_sample.output_directory)
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step = RunnableStep(
                name='kallisto',
                program='kallisto',
                sub_command=Command(name='quant', program='quant'))
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            # Read configuration section [bsf.analyses.kallisto.Kallisto.kallisto]
            self.set_runnable_step_configuration(runnable_step=runnable_step)

            runnable_step.sub_command.add_option_long(key='index', value=self.transcriptome_index_path)
            runnable_step.sub_command.add_option_long(key='output-dir', value=file_path_sample.output_directory)
            runnable_step.sub_command.add_option_long(key='threads', value=str(stage_sample.threads))

            if self.bootstrap_samples:
                runnable_step.sub_command.add_option_long(key='bootstrap-samples', value=self.bootstrap_samples)

            if 'Library Type' in sample.annotation_dict:
                if len(sample.annotation_dict['Library Type']) > 1:
                    module_logger.warning(
                        "More than one annotation for sample %r and key 'Library Type'.",
                        sample.name)

                library_type = sample.annotation_dict['Library Type'][0]
                if library_type not in ('unstranded', 'first', 'second'):
                    raise Exception(f"The 'Library Type' annotation for sample {sample.name!r} "
                                    "has to be one of 'unstranded', 'first' or 'second'.")

                if library_type == 'second':
                    runnable_step.sub_command.add_switch_long(key='rf-stranded')
                elif library_type == 'first':
                    runnable_step.sub_command.add_switch_long(key='fr-stranded')

            single_end_mode = True
            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    module_logger.debug('PairedReads name:', paired_reads.get_name())
                    module_logger.log(logging.DEBUG - 2, 'PairedReads: %r', paired_reads)
                    if paired_reads.reads_1 is not None:
                        runnable_step.sub_command.arguments.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2 is not None:
                        runnable_step.sub_command.arguments.append(paired_reads.reads_2.file_path)
                        single_end_mode = False

            if single_end_mode:
                runnable_step.sub_command.add_switch_long(key='single')

                if 'Fragment Length Value' in sample.annotation_dict:
                    if len(sample.annotation_dict['Fragment Length Value']) > 1:
                        module_logger.warning(
                            "More than one annotation for sample %r and key 'Fragment Length Value'.",
                            sample.name)
                    runnable_step.sub_command.add_option_long(
                        key='fragment-length',
                        value=sample.annotation_dict['Fragment Length Value'][0])
                elif self.fragment_length_value:
                    runnable_step.sub_command.add_option_long(
                        key='fragment-length',
                        value=self.fragment_length_value)

                if 'Fragment Length Standard Deviation' in sample.annotation_dict:
                    if len(sample.annotation_dict['Fragment Length Standard Deviation']) > 1:
                        module_logger.warning(
                            "More than one annotation for sample %r and key 'Fragment Length Standard Deviation'.",
                            sample.name)
                    runnable_step.sub_command.add_option_long(
                        key='sd',
                        value=sample.annotation_dict['Fragment Length Standard Deviation'][0])
                elif self.fragment_length_standard_deviation:
                    runnable_step.sub_command.add_option_long(
                        key='sd',
                        value=self.fragment_length_standard_deviation)

        return
