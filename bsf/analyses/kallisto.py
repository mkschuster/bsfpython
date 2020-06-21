# -*- coding: utf-8 -*-
"""Kallisto Analysis module.

A package of classes and methods supporting Kallisto analyses.
See https://pachterlab.github.io/kallisto/manual
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
import sys
import warnings

import bsf.analysis
import bsf.procedure
import bsf.process
import bsf.standards


class FilePathSample(bsf.procedure.FilePath):
    """The C{bsf.analyses.kallisto.FilePathSample} models files in a sample-specific Kallisto directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.kallisto.FilePathSample} object

        @param prefix: Prefix
        @type prefix: str
        @return:
        @rtype:
        """
        super(FilePathSample, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundance_h5 = os.path.join(prefix, 'abundance.h5')
        self.abundance_tsv = os.path.join(prefix, 'abundance.tsv')
        self.run_info_json = os.path.join(prefix, 'run_info.json')

        return


class Kallisto(bsf.analysis.Analysis):
    """Kallisto C{bsf.analysis.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str | None
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | None
    @ivar fragment_length_value: Fragment length
    @type fragment_length_value: str | None
    @ivar fragment_length_standard_deviation: Fragment length standard deviation
    @type fragment_length_standard_deviation: str | None
    @ivar bias_correction: Enable sequence-specific bias correction
    @type bias_correction: bool
    @ivar bootstrap_samples: Number of bootstrap samples
    @type bootstrap_samples: str
    """

    name = 'Kallisto Analysis'
    prefix = 'kallisto'

    @classmethod
    def get_stage_name_sample(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'sample'))

    @classmethod
    def get_prefix_sample(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: Sample name
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_sample(), sample_name))

    @classmethod
    def get_file_path_sample(cls, sample_name):
        """Get a C{FilePathSample} object from this or a sub-class.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: C{FilePathSample} or sub-class object
        @rtype: FilePathSample
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
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            transcriptome_version=None,
            transcriptome_index_path=None,
            fragment_length_value=None,
            fragment_length_standard_deviation=None,
            bias_correction=None,
            bootstrap_samples=None):
        """Initialise a C{bsf.analyses.kallisto.Kallisto} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.analysis.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.analysis.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.analysis.Analysis}-wide project directory,
            normally under the C{bsf.analysis.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.analysis.Analysis}-wide genome directory,
            normally under the C{bsf.analysis.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.analysis.Stage} objects
        @type stage_list: list[bsf.analysis.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | None
        @param fragment_length_value: Fragment length
        @type fragment_length_value: str | None
        @param fragment_length_standard_deviation: Fragment length standard deviation
        @type fragment_length_standard_deviation: str | None
        @param bias_correction: Enable sequence-specific bias correction
        @type bias_correction: bool
        @param bootstrap_samples: Number of bootstrap samples
        @type bootstrap_samples: str
        @return:
        @rtype:
        """
        super(Kallisto, self).__init__(
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

        self.transcriptome_version = transcriptome_version
        self.transcriptome_index_path = transcriptome_index_path
        self.fragment_length_value = fragment_length_value
        self.fragment_length_standard_deviation = fragment_length_standard_deviation
        self.bias_correction = bias_correction
        self.bootstrap_samples = bootstrap_samples

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.rnaseq.Tuxedo} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
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
        """Run this C{bsf.analyses.rnaseq.Tuxedo} analysis.
        @return:
        @rtype:
        """
        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.analysis.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.analysis.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """
            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Kallisto requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A ' + self.name + " requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = bsf.standards.Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a valid 'transcriptome_version' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the bsf.analysis.Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[self.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        super(Kallisto, self).run()

        if not self.transcriptome_index_path:
            self.transcriptome_index_path = os.path.join(
                bsf.standards.FilePath.get_resource_transcriptome_index(
                    transcriptome_version=self.transcriptome_version,
                    transcriptome_index='kallisto'),
                self.transcriptome_version + '.idx')

        if not os.path.exists(self.transcriptome_index_path):
            raise Exception('The ' + self.name + ' transcriptome index file path does not exist.\n' +
                            self.transcriptome_index_path)

        run_read_comparisons()

        stage_sample = self.get_stage(name=self.get_stage_name_sample())

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=True, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            prefix_sample = self.get_prefix_sample(sample_name=sample.name)

            file_path_sample = FilePathSample(prefix=prefix_sample)

            runnable_sample = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=prefix_sample,
                    working_directory=self.genome_directory,
                    debug=self.debug))
            # executable_sample =
            self.set_stage_runnable(
                stage=stage_sample,
                runnable=runnable_sample)
            # Set dependencies on previous Runnable or bsf.process.Executable objects.
            # executable_sample.dependencies.append(executable_sample.name)

            runnable_step = bsf.process.RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=file_path_sample.output_directory)
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step = bsf.process.RunnableStep(
                name='kallisto',
                program='kallisto',
                sub_command=bsf.process.Command(name='quant', program='quant'))
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
                    warnings.warn('More than one annotation for sample ' + repr(sample.name) +
                                  " and key 'Library Type'.")

                library_type = sample.annotation_dict['Library Type'][0]
                if library_type not in ('unstranded', 'first', 'second'):
                    raise Exception("The 'Library Type' annotation for sample " + repr(sample.name) +
                                    " has to be one of 'unstranded', 'first' or 'second'.")

                if library_type == 'second':
                    runnable_step.sub_command.add_switch_long(key='rf-stranded')
                elif library_type == 'first':
                    runnable_step.sub_command.add_switch_long(key='fr-stranded')

            single_end_mode = True
            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())
                    if paired_reads.reads_1 is not None:
                        runnable_step.sub_command.arguments.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2 is not None:
                        runnable_step.sub_command.arguments.append(paired_reads.reads_2.file_path)
                        single_end_mode = False

            if single_end_mode:
                runnable_step.sub_command.add_switch_long(key='single')

                if 'Fragment Length Value' in sample.annotation_dict:
                    if len(sample.annotation_dict['Fragment Length Value']) > 1:
                        warnings.warn('More than one annotation for sample ' + repr(sample.name) +
                                      " and key 'Fragment Length Value'.")
                    runnable_step.sub_command.add_option_long(
                        key='fragment-length',
                        value=sample.annotation_dict['Fragment Length Value'][0])
                elif self.fragment_length_value:
                    runnable_step.sub_command.add_option_long(
                        key='fragment-length',
                        value=self.fragment_length_value)

                if 'Fragment Length Standard Deviation' in sample.annotation_dict:
                    if len(sample.annotation_dict['Fragment Length Standard Deviation']) > 1:
                        warnings.warn('More than one annotation for sample ' + repr(sample.name) +
                                      " and key 'Fragment Length Standard Deviation'.")
                    runnable_step.sub_command.add_option_long(
                        key='sd',
                        value=sample.annotation_dict['Fragment Length Standard Deviation'][0])
                elif self.fragment_length_standard_deviation:
                    runnable_step.sub_command.add_option_long(
                        key='sd',
                        value=self.fragment_length_standard_deviation)

        return
