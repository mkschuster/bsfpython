# -*- coding: utf-8 -*-
"""Trimmomatic Analysis module.

A package of classes and methods supporting the Trimmomatic tool.
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

import bsf.analysis
import bsf.executables.collection
import bsf.ngs
import bsf.procedure
import bsf.process
import bsf.standards


class FilePathTrimmomaticReadGroup(bsf.procedure.FilePath):
    """The C{bsf.analyses.trimmomatic.FilePathTrimmomaticReadGroup} models read group-specific Trimmomatic files.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar trim_log_tsv: Trimmomatic trim log Tab-Separated Value (TSV) file path
    @type trim_log_tsv: str
    @ivar summary_tsv: Summary Tab-Separated Value (TSV) file path
    @type summary_tsv: str
    @ivar coverage_png: Coverage Portable Network Graphics (PNG) file path
    @type coverage_png: str
    @ivar frequency_png: Frequency Portable Network Graphics (PNG) file path
    @type frequency_png: str
    @ivar surviving_png: Surviving Portable Network Graphics (PNG) file path
    @type surviving_png: str
    @ivar reads_1p: First Reads paired
    @type reads_1p: str | None
    @ivar reads_1u: First Reads unpaired
    @type reads_1u: str | None
    @ivar reads_2p: Second Reads paired
    @type reads_2p: str | None
    @ivar reads_2u: Second Reads unpaired
    @type reads_2u: str | None
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.trimmomatic.FilePathTrimmomaticReadGroup} object.

        @param prefix: Prefix
        @type prefix: str
        @return:
        @rtype:
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


class FilePathTrimmomaticProject(bsf.procedure.FilePath):
    """The C{bsf.analyses.trimmomatic.FilePathTrimmomaticProject} class models project-specific Trimmomatic files.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar sas_path_old: Old Sample Annotation Sheet file path
    @type sas_path_old: str
    @ivar sas_path_new: New Sample Annotation Sheet file path
    @type sas_path_new: str
    """

    def __init__(self, prefix, prefix_analysis, project_name):
        """Initialise a C{bsf.analyses.trimmomatic.FilePathTrimmomaticProject} object.

        @param prefix: Prefix
        @type prefix: str
        @param prefix_analysis: C{bsf.analysis.Analysis.prefix}
        @type prefix_analysis: str
        @param project_name: Project name
        @type project_name: str
        @return:
        @rtype:
        """
        super(FilePathTrimmomaticProject, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sas_path_old = '_'.join((project_name, prefix_analysis, 'original.csv'))
        self.sas_path_new = '_'.join((project_name, prefix_analysis, 'samples.csv'))

        return


class Trimmomatic(bsf.analysis.Analysis):
    """The C{bsf.analyses.trimmomatic.Trimmomatic} class represents the logic to run the Trimmomatic analysis.

    Attributes:

    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar adapter_path: Adapter file path
    @type adapter_path: str | None
    @ivar trimming_step_pe_list: Colon-separated Trimmomatic steps for paired-end data
    @type trimming_step_pe_list: list[str] | None
    @ivar trimming_step_se_list: Colon-separated Trimmomatic steps for single-end data
    @type trimming_step_se_list: list[str] | None
    @ivar classpath_trimmomatic: Trimmomatic tool Java Archive (JAR) class path directory
    @type classpath_trimmomatic: str | None
    """

    name = 'Trimmomatic Analysis'
    prefix = 'trimmomatic'

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_stage_name_summary(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'summary'))

    @classmethod
    def get_stage_name_project(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'project'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_prefix_summary(cls):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return cls.get_stage_name_summary()

    @classmethod
    def get_prefix_project(cls):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return cls.get_stage_name_project()

    @classmethod
    def get_file_path_read_group(cls, read_group_name):
        """Get a C{FilePathTrimmomaticReadGroup} object.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: C{FilePathTrimmomaticReadGroup} object
        @rtype: FilePathTrimmomaticReadGroup
        """
        return FilePathTrimmomaticReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

    @classmethod
    def get_file_path_project(cls, project_name, prefix_analysis):
        """Get a C{FilePathTrimmomaticProject} object.

        @param project_name: Project name
        @type project_name: str
        @param prefix_analysis: C{bsf.analysis.Analysis.prefix}
        @type prefix_analysis: str
        @return: C{FilePathTrimmomaticProject} object
        @rtype: FilePathTrimmomaticProject
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
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            adapter_path=None,
            trimming_step_pe_list=None,
            trimming_step_se_list=None,
            classpath_trimmomatic=None):
        """Initialise a C{bsf.analyses.trimmomatic.Trimmomatic} object.

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
        @param adapter_path: Adapter file path
        @type adapter_path: str | None
        @param trimming_step_pe_list: Colon-separated Trimmomatic steps for paired-end data
        @type trimming_step_pe_list: list[str] | None
        @param trimming_step_se_list: Colon-separated Trimmomatic steps for single-end data
        @type trimming_step_se_list: list[str] | None
        @param classpath_trimmomatic: Trimmomatic tool Java Archive (JAR) class path directory
        @type classpath_trimmomatic: str | None
        @return:
        @rtype:
        """
        super(Trimmomatic, self).__init__(
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

        self.adapter_path = adapter_path
        self.trimming_step_pe_list = trimming_step_pe_list
        self.trimming_step_se_list = trimming_step_se_list
        self.classpath_trimmomatic = classpath_trimmomatic

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.trimmomatic.Trimmomatic} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
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

        # Get the Trimmomatic tool Java Archive (JAR) class path directory.

        option = 'classpath_trimmomatic'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_trimmomatic = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.trimmomatic.Trimmomatic} C{bsf.analysis.Analysis}.

        This method changes the C{bsf.ngs.Collection} object of this C{bsf.analysis.Analysis}
        to update with FASTQ file paths.
        @return:
        @rtype:
        """

        def run_adjust_illumina_clip_path(trimming_step_list):
            """Private function to adjust the adapter FASTA file path of ILLUMINACLIP trimming steps.

            If the file path is not absolute, prepend the value of the adapter_path
            instance variable.
            @param trimming_step_list: Python C{list} of trimming steps.
            @type trimming_step_list: list[str]
            @return:
            @rtype:
            """
            for i in range(0, len(trimming_step_list)):
                if trimming_step_list[i].startswith('ILLUMINACLIP'):
                    component_list = trimming_step_list[i].split(':')
                    if not os.path.isabs(component_list[1]):
                        component_list[1] = os.path.join(self.adapter_path, component_list[1])
                    trimming_step_list[i] = ':'.join(component_list)
            return

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

        super(Trimmomatic, self).run()

        # Get the Trimmomatic tool Java Archive (JAR) class path directory.

        if not self.classpath_trimmomatic:
            self.classpath_trimmomatic = bsf.standards.JavaClassPath.get_trimmomatic()
            if not self.classpath_trimmomatic:
                raise Exception('A ' + self.name + "requires a 'classpath_trimmomatic' configuration option.")

        if not (self.adapter_path and os.path.isabs(self.adapter_path)):
            self.adapter_path = os.path.join(os.path.dirname(self.classpath_trimmomatic), 'adapters')

        if self.trimming_step_pe_list is None:
            raise Exception('A ' + self.name + " requires a 'trimming_steps_pe' configuration option.")

        if self.trimming_step_se_list is None:
            raise Exception('A ' + self.name + " requires a 'trimming_steps_se' configuration option.")

        run_adjust_illumina_clip_path(trimming_step_list=self.trimming_step_pe_list)
        run_adjust_illumina_clip_path(trimming_step_list=self.trimming_step_se_list)

        run_read_comparisons()

        # Trimmomatic

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())
        stage_summary = self.get_stage(name=self.get_stage_name_summary())
        stage_project = self.get_stage(name=self.get_stage_name_project())

        project_dependency_list = list()
        """ @type project_dependency_list: list[str] """

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            sample_step_list = list()
            """ @type sample_step_list: list[str] """
            if 'Trimmomatic Steps' in sample.annotation_dict:
                for trimming_step in sample.annotation_dict['Trimmomatic Steps']:
                    sample_step_list.extend(self.configuration.list_from_csv(csv_string=trimming_step))
                run_adjust_illumina_clip_path(trimming_step_list=sample_step_list)

            # The Trimmomatic analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=False)

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    paired_reads_step_list = list()
                    if 'Trimmomatic Steps' in paired_reads.annotation_dict:
                        for trimming_step in paired_reads.annotation_dict['Trimmomatic Steps']:
                            paired_reads_step_list.extend(self.configuration.list_from_csv(csv_string=trimming_step))
                        run_adjust_illumina_clip_path(trimming_step_list=paired_reads_step_list)

                    # Apply some sanity checks.

                    # Maybe this case should be allowed after Trimmomatic trimming,
                    # where only the second Read survives.
                    if paired_reads.reads_1 is None and paired_reads.reads_2 is not None:
                        raise Exception('PairedReads object with a reads_2, but no reads_1 object.')

                    file_path_read_group = self.get_file_path_read_group(read_group_name=paired_reads_name)

                    # Create a Runnable and an Executable for running the Trimmomatic analysis.

                    runnable_read_group = self.add_runnable_consecutive(
                        runnable=bsf.procedure.ConsecutiveRunnable(
                            name=self.get_prefix_read_group(read_group_name=paired_reads_name),
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory))
                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                    # Record the Executable.name for the project dependency.

                    project_dependency_list.append(runnable_read_group.name)

                    # Create a new RunnableStepMakeDirectory in preparation of the Trimmomatic program.

                    runnable_step_read_group = bsf.process.RunnableStepMakeDirectory(
                        name='mkdir',
                        directory_path=file_path_read_group.output_directory)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step_read_group)

                    # Create a RunnableStep for the Trimmomatic program.

                    runnable_step_read_group = bsf.process.RunnableStepJava(
                        name='trimmomatic',
                        java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx4G',
                        java_jar_path=self.classpath_trimmomatic)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step_read_group)

                    if paired_reads.reads_2 is None:
                        runnable_step_read_group.sub_command.sub_command = bsf.process.Command(program='SE')
                    else:
                        runnable_step_read_group.sub_command.sub_command = bsf.process.Command(program='PE')

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
                            paired_reads=bsf.ngs.PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=bsf.ngs.Reads(
                                    name=paired_reads.reads_1.name[:-1] + 'U',
                                    file_path=os.path.join(
                                        self.genome_directory,
                                        file_path_read_group.reads_1u)),
                                exclude=paired_reads.exclude,
                                index_1=paired_reads.index_1,
                                index_2=paired_reads.index_2,
                                read_group=paired_reads.read_group))

                        sample.add_paired_reads(
                            paired_reads=bsf.ngs.PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=bsf.ngs.Reads(
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

                    runnable_step_read_group = bsf.process.RunnableStep(
                        name='compress_logs',
                        program='pigz')
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step_read_group)

                    runnable_step_read_group.add_switch_long(key='best')
                    runnable_step_read_group.add_option_long(key='processes', value=str(stage_summary.threads))
                    runnable_step_read_group.arguments.append(file_path_read_group.trim_log_tsv)

                    # Create a Runnable for the bsf_trimmomatic_summary.R analysis.

                    prefix_summary = '_'.join((stage_summary.name, paired_reads_name))

                    runnable_summary = self.add_runnable_consecutive(
                        runnable=bsf.procedure.ConsecutiveRunnable(
                            name=prefix_summary,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory))
                    executable_summary = self.set_stage_runnable(stage=stage_summary, runnable=runnable_summary)
                    executable_summary.dependencies.append(runnable_read_group.name)

                    # Create a new RunnableStep to aggregate the trim log file.

                    runnable_step_summary = bsf.process.RunnableStep(
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
            runnable=bsf.procedure.ConsecutiveRunnable(
                name=prefix_project,
                code_module='bsf.runnables.generic',
                working_directory=self.project_directory))
        executable_project = self.set_stage_runnable(
            stage=stage_project,
            runnable=runnable_project)
        executable_project.dependencies.extend(project_dependency_list)

        # Create a new RunnableStep.

        runnable_step_project = bsf.executables.collection.RunnableStepCollectionPruneFastq(
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
        """Create a C{bsf.analyses.trimmomatic.Trimmomatic} report in HTML format and a
        UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # Create a symbolic link containing the project name and a UUID.
        self.create_public_project_link()

        # Write a HTML document.

        report_list = list()
        """ @type report_list: list[str] """

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
