"""bsf.analyses.bowtie

A package of classes and methods supporting Bowtie alignment analyses.
"""

#
# Copyright 2013 - 2018 Michael K. Schuster
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

from __future__ import print_function

import os
import re
import sys
import warnings

import bsf
import bsf.process


class Bowtie1(bsf.Analysis):
    """The C{bsf.analyses.bowtie.Bowtie1} class represents the logic to run the Bowtie1 short read aligner.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_align: C{bsf.Stage.name} for the Bowtie1 alignment C{bsf.Analysis} stage
    @type stage_name_align: str
    @ivar replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects for processing or run them separately
    @type replicate_grouping: bool | None
    @ivar index_basename: Bowtie index basename
    @type index_basename: str | unicode | None
    """

    name = 'Bowtie1 Analysis'
    prefix = 'bowtie1'

    stage_name_align = '_'.join((prefix, 'align'))

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
            replicate_grouping=None,
            index_basename=None):
        """Initialise a C{bsf.analyses.bowtie.Bowtie1}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str | None
        @param genome_version: Genome version
        @type genome_version: str | None
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str | None
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str | None
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str | None
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str | None
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str | None
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage] | None
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection | None
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample] | None
        @param replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects for processing or
            run them separately
        @type replicate_grouping: bool | None
        @param index_basename: Bowtie index basename
        @type index_basename: str | unicode | None
        @return:
        @rtype:
        """
        super(Bowtie1, self).__init__(
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

        self.replicate_grouping = replicate_grouping
        self.index_basename = index_basename

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.bowtie.Bowtie1} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(Bowtie1, self).set_configuration(configuration=configuration, section=section)

        config_parser = configuration.config_parser

        option = 'replicate_grouping'
        if config_parser.has_option(section=section, option=option):
            self.replicate_grouping = config_parser.getboolean(section=section, option=option)

        # Get the genome database.

        option = 'index_basename'
        if config_parser.has_option(section=section, option=option):
            self.index_basename = config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a C{bsf.analyses.bowtie.Bowtie1} analysis.

        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(path=self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = '_'.join((self.project_name, self.prefix, 'samples.csv'))
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # The Bowtie1 aligner requires a genome version.

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a 'genome_version' configuration option.")

        super(Bowtie1, self).run()

        if not self.index_basename:
            self.index_basename = os.path.join(
                bsf.standards.FilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie1'),
                self.genome_version)

        run_read_comparisons()

        stage_align = self.get_stage(name=self.stage_name_align)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', repr(sample.name))
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            for paired_reads_name in sorted(paired_reads_dict):
                if not paired_reads_dict[paired_reads_name]:
                    # Skip replicate keys, which PairedReads objects have all been excluded.
                    continue

                prefix_align = '_'.join((stage_align.name, paired_reads_name))

                file_path_align = FilePathAlign(prefix=prefix_align)

                # Create a Runnable and Executable for processing each Bowtie1 alignment.

                runnable_align = self.add_runnable(
                    runnable=bsf.Runnable(
                        name=prefix_align,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        file_path_object=file_path_align,
                        debug=self.debug))
                self.set_stage_runnable(
                    stage=stage_align,
                    runnable=runnable_align)

                runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_align.output_directory))

                runnable_step = runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStep(
                        name='bowtie1',
                        program='bowtie'))
                """ @type runnable_step: bsf.process.RunnableStep """
                runnable_step.arguments.append(self.index_basename)

                file_path_list_1 = list()
                file_path_list_2 = list()

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if paired_reads.reads_1 is not None:
                        file_path_list_1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2 is not None:
                        file_path_list_2.append(paired_reads.reads_2.file_path)

                if len(file_path_list_1) and not len(file_path_list_2):
                    # For Bowtie1 unpaired reads are an argument, paired come with options -1 <m1> and -2 <m2>.
                    runnable_step.arguments.append(','.join(file_path_list_1))
                elif len(file_path_list_1) and len(file_path_list_2):
                    runnable_step.add_option_short(key='1', value=','.join(file_path_list_1))
                if len(file_path_list_2):
                    runnable_step.add_option_short(key='2', value=','.join(file_path_list_2))

                runnable_step.arguments.append(file_path_align.aligned_sam)

        return


class FilePathAlign(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.bowtie.FilePathAlign} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'aligned.sam')))

        self.cleaned_sam = os.path.join(prefix, '_'.join((prefix, 'cleaned.sam')))

        self.aligned_bai = os.path.join(prefix, '_'.join((prefix, 'aligned.bai')))
        self.aligned_bam = os.path.join(prefix, '_'.join((prefix, 'aligned.bam')))
        self.aligned_lnk = os.path.join(prefix, '_'.join((prefix, 'aligned.bam.bai')))
        self.aligned_md5 = os.path.join(prefix, '_'.join((prefix, 'aligned.bam.md5')))

        return


class FilePathIndex(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.bowtie.FilePathIndex} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathIndex, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.merged_bai = os.path.join(prefix, '_'.join((prefix, 'merged.bai')))
        self.merged_bam = os.path.join(prefix, '_'.join((prefix, 'merged.bam')))
        self.merged_md5 = os.path.join(prefix, '_'.join((prefix, 'merged.bam.md5')))

        self.merged_bai_link_source = '_'.join((prefix, 'merged.bai'))
        self.merged_bai_link_target = os.path.join(prefix, '_'.join((prefix, 'merged.bam.bai')))

        self.sorted_bai = os.path.join(prefix, '_'.join((prefix, 'sorted.bai')))
        self.sorted_bam = os.path.join(prefix, '_'.join((prefix, 'sorted.bam')))
        self.sorted_md5 = os.path.join(prefix, '_'.join((prefix, 'sorted.bam.md5')))

        self.sorted_bai_link_source = '_'.join((prefix, 'sorted.bai'))
        self.sorted_bai_link_target = os.path.join(prefix, '_'.join((prefix, 'sorted.bam.bai')))

        return


class FilePathMerge(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.bowtie.FilePathMerge} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathMerge, self).__init__(prefix=prefix)

        self.output_directory = prefix

        # For the final files, just end in the suffix.
        self.merged_bai = os.path.join(prefix, prefix + '.bai')
        self.merged_bam = os.path.join(prefix, prefix + '.bam')
        self.merged_md5 = os.path.join(prefix, prefix + '.bam.md5')

        self.merged_bai_link_source = prefix + '.bai'
        self.merged_bai_link_target = os.path.join(prefix, prefix + '.bam.bai')

        self.alignment_summary_metrics = os.path.join(prefix, '_'.join((prefix, 'alignment_summary_metrics.tsv')))

        return


class FilePathSummary(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.bowtie.FilePathSummary} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathSummary, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class Bowtie2(bsf.Analysis):
    """The C{bsf.analyses.bowtie.Bowtie2} class represents the logic to run the Bowtie2 short read aligner.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_align: C{bsf.Stage.name} for the alignment C{bsf.Analysis} stage
    @type stage_name_align: str
    @cvar stage_name_index: C{bsf.Stage.name} for the indexing C{bsf.Analysis} stage
    @type stage_name_index: str
    @cvar stage_name_merge: C{bsf.Stage.name} for the merging C{bsf.Analysis} stage
    @type stage_name_merge: str
    @cvar stage_name_summary: C{bsf.Stage.name} for the summary C{bsf.Analysis} stage
    @type stage_name_summary: str
    @ivar genome_fasta: Genome FASTA file
    @type genome_fasta: str | unicode | None
    @ivar genome_index: Bowtie genome index basename (i.e. genome without *.fasta)
    @type genome_index: str | unicode | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: None | str | unicode
    """

    name = 'Bowtie2 Analysis'
    prefix = 'bowtie2'

    stage_name_align = '_'.join((prefix, 'align'))
    stage_name_index = '_'.join((prefix, 'index'))
    stage_name_merge = '_'.join((prefix, 'merge'))
    stage_name_summary = '_'.join((prefix, 'summary'))

    @classmethod
    def get_prefix_align(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_align, paired_reads_name))

    @classmethod
    def get_prefix_index(cls, bam_file_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param bam_file_name: Unmapped BAM file name
        @type bam_file_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_index, bam_file_name))

    @classmethod
    def get_prefix_merge(cls, sample_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_merge, sample_name))

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
            genome_fasta=None,
            genome_index=None,
            classpath_picard=None):
        """Initialise a C{bsf.analyses.bowtie.Bowtie1}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str | None
        @param genome_version: Genome version
        @type genome_version: str | None
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str | None
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str | None
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str | None
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str | None
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str | None
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage] | None
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection | None
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample] | None
        @param genome_fasta: Genome FASTA file
        @type genome_fasta: str | unicode | None
        @param genome_index: Bowtie genome index basename (i.e. genome without *.fasta)
        @type genome_index: str | unicode | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: None | str | unicode
        @return:
        @rtype:
        """
        super(Bowtie2, self).__init__(
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

        self.genome_fasta = genome_fasta
        self.genome_index = genome_index
        self.classpath_picard = classpath_picard

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.bowtie.Bowtie2} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(Bowtie2, self).set_configuration(configuration=configuration, section=section)

        config_parser = configuration.config_parser

        # Get the genome database.

        option = 'genome_fasta'
        if config_parser.has_option(section=section, option=option):
            self.genome_fasta = config_parser.get(section=section, option=option)

        option = 'genome_index'
        if config_parser.has_option(section=section, option=option):
            self.genome_index = config_parser.get(section=section, option=option)

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a C{bsf.analyses.bowtie.Bowtie2} analysis.

        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        def run_get_file_name(_paired_reads):
            """Get the BAM file base name annotation of a PairedReads object.

            @param _paired_reads: C{bsf.ngs.PairedReads} object
            @type _paired_reads: bsf.ngs.PairedReads
            @return: Platform unit
            @rtype: (str, str | unicode) | None
            """
            if 'BAM File' in _paired_reads.annotation_dict and paired_reads.annotation_dict['BAM File']:
                bam_name, bam_extension = os.path.splitext(
                    os.path.basename(paired_reads.annotation_dict['BAM File'][0L]))

                # TODO: This should be centralised.
                # The makeFileNameSafe() method of htsjdk.samtools.util.IOUtil uses the following pattern:
                # [\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]
                bam_name = re.sub(
                    pattern='[\\s!"#$%&\'()*/:;<=>?@\\[\\]\\\\^`{|}~]',
                    repl='_',
                    string=bam_name)

                return bam_name, paired_reads.annotation_dict['BAM File'][0L]

        # Start of the run() method body.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(path=self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = '_'.join((self.project_name, self.prefix, 'samples.csv'))
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # The Bowtie2 aligner requires a genome version.

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a 'genome_version' configuration option.")

        super(Bowtie2, self).run()

        if not self.genome_fasta:
            # Set the genome_index to None, to point to the genome directory that contains the
            # Picard sequence dictionary.
            self.genome_fasta = bsf.standards.FilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index=None)

        if not self.genome_index:
            self.genome_index = os.path.join(
                bsf.standards.FilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('A ' + self.name + " requires a 'classpath_picard' configuration option.")

        run_read_comparisons()

        stage_align = self.get_stage(name=self.stage_name_align)
        stage_index = self.get_stage(name=self.stage_name_index)
        stage_merge = self.get_stage(name=self.stage_name_merge)
        stage_summary = self.get_stage(name=self.stage_name_summary)

        prefix_summary = stage_summary.name

        file_path_summary = FilePathSummary(prefix=prefix_summary)

        runnable_merge_list = list()
        """ @type runnable_merge_list: list[bsf.Runnable] """

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', repr(sample.name))
                sys.stdout.writelines(sample.trace(level=1))

            # To run Picard MergeBamAlignment, all alignments from a BAM file need merging into one.

            bam_file_dict = dict()
            """ @type bam_file_dict: dict[str | unicode, (str | unicode, list[bsf.Runnable])] """

            runnable_index_list = list()
            """ @type runnable_index_list: list[bsf.Runnable] """

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            for paired_reads_name in sorted(paired_reads_dict):
                if not paired_reads_dict[paired_reads_name]:
                    # Skip replicate keys, which PairedReads objects have all been excluded.
                    continue

                if len(paired_reads_dict[paired_reads_name]) > 1:
                    raise Exception('Cannot align and process more than one PairedReads object at a time.')

                # Record the Runnable for the Picard MergeBamAlignment step.
                paired_reads = paired_reads_dict[paired_reads_name][0L]

                # Get the file paths for Reads1 and Reads2 anc check for FASTQ files.
                if paired_reads.reads_1 is None:
                    raise Exception('A ' + self.name + ' requires a Reads1 object.')
                else:
                    if paired_reads.reads_1.file_path is not None:
                        if not (paired_reads.reads_1.file_path.endswith('.fastq') or
                                paired_reads.reads_1.file_path.endswith('.fastq.gz')):
                            raise Exception('A ' + self.name + ' requires a (GNU Zip compressed) FASTQ file.')
                    file_path_1 = paired_reads.reads_1.file_path

                if paired_reads.reads_2 is None:
                    file_path_2 = None
                else:
                    if paired_reads.reads_2.file_path is not None:
                        if not (paired_reads.reads_2.file_path.endswith('.fastq') or
                                paired_reads.reads_2.file_path.endswith('.fastq.gz')):
                            raise Exception('A ' + self.name + ' requires a (GNU Zip compressed) FASTQ file.')
                    file_path_2 = paired_reads.reads_2.file_path

                ######################
                # 1. Alignment Stage #
                ######################

                prefix_align = self.get_prefix_align(paired_reads_name=paired_reads_name)

                file_path_align = FilePathAlign(prefix=prefix_align)

                # Create a Runnable and Executable for processing each Bowtie1 alignment.

                runnable_align = self.add_runnable(
                    runnable=bsf.Runnable(
                        name=prefix_align,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        file_path_object=file_path_align,
                        debug=self.debug))
                self.set_stage_runnable(
                    stage=stage_align,
                    runnable=runnable_align)

                # Record the alignment Runnable under the BAM file name
                # to merge the sub-alignments of all (trimmed) FASTQ files that
                # resulted from an unmapped BAM file.
                bam_file_name, bam_file_path = run_get_file_name(_paired_reads=paired_reads)
                if bam_file_name:
                    if bam_file_name in bam_file_dict:
                        bam_file_dict[bam_file_name][1L].append(runnable_align)
                    else:
                        bam_file_dict[bam_file_name] = (bam_file_path, [runnable_align])
                else:
                    warnings.warn("PairedReads object '" + paired_reads.get_name() + "without 'BAM File' annotation.")

                # Make a directory.

                runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_align.output_directory))

                # Run the Bowtie2 aligner.

                runnable_step = runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStep(
                        name='bowtie2',
                        program='bowtie2'))
                """ @type runnable_step: bsf.process.RunnableStep """
                runnable_step.add_option_short(key='S', value=file_path_align.aligned_sam)
                runnable_step.add_option_short(key='x', value=self.genome_index)
                runnable_step.add_option_long(key='threads', value=str(stage_align.threads))

                # NOTE: The following options are properties of the Sample,
                # PairedReads and Reads objects.
                # runnable_step.add_switch_short(key='q')
                # runnable_step.add_switch_short(key='phred33')

                if file_path_1 and not file_path_2:
                    runnable_step.add_option_short(key='U', value=file_path_1)
                else:
                    runnable_step.add_option_short(key='1', value=file_path_1)
                    runnable_step.add_option_short(key='2', value=file_path_2)

                # Run Picard CleanSam to convert the aligned SAM file into a cleaned SAM file.

                runnable_step = runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_clean_sam',
                        obsolete_file_path_list=[
                            file_path_align.aligned_sam,
                        ],
                        java_temporary_path=runnable_align.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='CleanSam'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                # INPUT []
                runnable_step.add_picard_option(key='INPUT', value=file_path_align.aligned_sam)
                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_align.cleaned_sam)
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_align.get_relative_temporary_directory_path)
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET [false]
                # VALIDATION_STRINGENCY [STRICT]
                # COMPRESSION_LEVEL [5]
                # MAX_RECORDS_IN_RAM [500000]
                # CREATE_INDEX [false]
                # CREATE_MD5_FILE [false]
                # REFERENCE_SEQUENCE [null]
                # GA4GH_CLIENT_SECRETS [client_secrets.json]
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]

                # Run Picard SortSam to convert the cleaned SAM file into a coordinate sorted BAM file.

                runnable_step = runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_sort_sam',
                        obsolete_file_path_list=[
                            file_path_align.cleaned_sam,
                        ],
                        java_temporary_path=runnable_align.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='SortSam'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                # INPUT []
                runnable_step.add_picard_option(key='INPUT', value=file_path_align.cleaned_sam)
                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_align.aligned_bam)
                # SORT_ORDER []
                runnable_step.add_picard_option(key='SORT_ORDER', value='coordinate')
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_align.get_relative_temporary_directory_path)
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET [false]
                # VALIDATION_STRINGENCY [STRICT]
                # COMPRESSION_LEVEL [5]
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                # MAX_RECORDS_IN_RAM [500000]
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                # CREATE_INDEX [false]
                runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                # CREATE_MD5_FILE [false]
                runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                # REFERENCE_SEQUENCE [null]
                # GA4GH_CLIENT_SECRETS [client_secrets.json]
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]

            #####################
            # 2. Indexing Stage #
            #####################

            # Merge all aligned BAM files of an original, unaligned BAM file.

            for bam_file_name, (bam_file_path, runnable_list) in bam_file_dict.iteritems():
                # Create a Runnable and Executable for merging each platform unit.
                prefix_index = self.get_prefix_index(bam_file_name=bam_file_name)

                file_path_index = FilePathIndex(prefix=prefix_index)

                runnable_index = self.add_runnable(
                    runnable=bsf.Runnable(
                        name=prefix_index,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        file_path_object=file_path_index,
                        debug=self.debug))
                executable_index = self.set_stage_runnable(
                    stage=stage_index,
                    runnable=runnable_index)

                # Record the dependency for all alignment Runnable objects.
                for runnable_align in runnable_list:
                    executable_index.dependencies.append(runnable_align.name)

                runnable_index_list.append(runnable_index)

                runnable_index.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_index.output_directory))

                if len(runnable_list) == 1L:
                    runnable_align = runnable_list[0L]
                    file_path_align = runnable_align.file_path_object
                    """ @type file_path_align: FilePathAlign """
                    # For a single ReadPair, just rename the files.
                    runnable_index.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_bam',
                            source_path=file_path_align.aligned_bam,
                            target_path=file_path_index.merged_bam))
                    runnable_index.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_bai',
                            source_path=file_path_align.aligned_bai,
                            target_path=file_path_index.merged_bai))
                    runnable_index.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_md5',
                            source_path=file_path_align.aligned_md5,
                            target_path=file_path_index.merged_md5))
                else:
                    runnable_step = runnable_index.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepPicard(
                            name='picard_merge_sam_files',
                            java_temporary_path=runnable_index.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx4G',
                            picard_classpath=self.classpath_picard,
                            picard_command='MergeSamFiles'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    # INPUT []
                    for runnable_align in runnable_list:
                        file_path_align = runnable_align.file_path_object
                        """ @type file_path_align: FilePathAlign """
                        runnable_step.add_picard_option(
                            key='INPUT',
                            value=file_path_align.aligned_bam,
                            override=True)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_bam)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_bai)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_md5)
                    # OUTPUT []
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_index.merged_bam)
                    # SORT_ORDER [coordinate]
                    # ASSUME_SORTED [false]
                    # MERGE_SEQUENCE_DICTIONARIES [false]
                    # USE_THREADING [false]
                    runnable_step.add_picard_option(key='USE_THREADING', value='true')
                    # COMMENT [null]
                    # INTERVALS [null]
                    # TMP_DIR [null]
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_index.get_relative_temporary_directory_path)
                    # VERBOSITY [INFO]
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    # QUIET [false]
                    # VALIDATION_STRINGENCY [STRICT]
                    # COMPRESSION_LEVEL [5]
                    runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                    # MAX_RECORDS_IN_RAM [500000]
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                    # CREATE_INDEX [false]
                    runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                    # CREATE_MD5_FILE [false]
                    runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                    # REFERENCE_SEQUENCE [null]
                    # GA4GH_CLIENT_SECRETS [client_secrets.json]
                    # USE_JDK_DEFLATER [false]
                    # USE_JDK_INFLATER [false]

                # Run Picard MergeBamAlignment to annotate the aligned with the unaligned BAM file.

                runnable_step = runnable_index.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_merge_bam_alignment',
                        obsolete_file_path_list=[
                            file_path_index.merged_bam,
                            file_path_index.merged_bai,
                            file_path_index.merged_md5,
                        ],
                        java_temporary_path=runnable_index.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx12G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='MergeBamAlignment'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                # UNMAPPED_BAM []
                runnable_step.add_picard_option(key='UNMAPPED_BAM', value=bam_file_path)
                # ALIGNED_BAM [null]
                runnable_step.add_picard_option(key='ALIGNED_BAM', value=file_path_index.merged_bam)
                # READ1_ALIGNED_BAM [null]
                # READ2_ALIGNED_BAM [null]
                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_index.sorted_bam)
                # PROGRAM_RECORD_ID [null]
                # PROGRAM_GROUP_VERSION [null]
                # PROGRAM_GROUP_COMMAND_LINE [null]
                # PROGRAM_GROUP_NAME [null]
                # CLIP_ADAPTERS [true]
                # IS_BISULFITE_SEQUENCE [false]
                # ALIGNED_READS_ONLY [false]
                # MAX_INSERTIONS_OR_DELETIONS [1]
                runnable_step.add_picard_option(key='MAX_INSERTIONS_OR_DELETIONS', value='-1')
                # ATTRIBUTES_TO_RETAIN [null]
                # ATTRIBUTES_TO_REMOVE [null]
                # ATTRIBUTES_TO_REVERSE [OQ, U2]
                # ATTRIBUTES_TO_REVERSE_COMPLEMENT [E2, SQ]
                # READ1_TRIM [0]
                # READ2_TRIM [0]
                # EXPECTED_ORIENTATIONS [null]
                # ALIGNER_PROPER_PAIR_FLAGS [false]
                # SORT_ORDER [coordinate]
                # PRIMARY_ALIGNMENT_STRATEGY [BestMapq]
                # CLIP_OVERLAPPING_READS [true]
                # INCLUDE_SECONDARY_ALIGNMENTS [true]
                # ADD_MATE_CIGAR [true]
                # UNMAP_CONTAMINANT_READS [false]
                # MIN_UNCLIPPED_BASES [32]
                # MATCHING_DICTIONARY_TAGS [M5, LN]
                # UNMAPPED_READ_STRATEGY [DO_NOT_CHANGE]
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_index.get_relative_temporary_directory_path)
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET [false]
                # VALIDATION_STRINGENCY [STRICT]
                # COMPRESSION_LEVEL [5]
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                # MAX_RECORDS_IN_RAM [500000]
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                # CREATE_INDEX [false]
                runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                # CREATE_MD5_FILE [false]
                runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                # REFERENCE_SEQUENCE [null]
                runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.genome_fasta)
                # GA4GH_CLIENT_SECRETS [client_secrets.json]
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]

            ####################
            # 3. Merging Stage #
            ####################

            # For more than one ReadPair object the aligned BAM files need merging into Sample-specific ones.

            prefix_merge = self.get_prefix_merge(sample_name=sample.name)

            file_path_merge = FilePathMerge(prefix=prefix_merge)

            runnable_merge = self.add_runnable(
                runnable=bsf.Runnable(
                    name=prefix_merge,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    file_path_object=file_path_merge,
                    debug=self.debug))
            executable_merge = self.set_stage_runnable(
                stage=stage_merge,
                runnable=runnable_merge)

            runnable_merge_list.append(runnable_merge)

            # Add dependencies on Runnable objects of the indexing stage.
            for runnable_index in runnable_index_list:
                executable_merge.dependencies.append(runnable_index.name)

            runnable_merge.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=file_path_merge.output_directory))

            if len(runnable_index_list) == 1L:
                runnable_index = runnable_index_list[0L]
                file_path_index = runnable_index.file_path_object
                """ @type file_path_index: FilePathIndex """
                # For a single ReadPair, just rename the files.
                runnable_merge.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_bam',
                        source_path=file_path_index.sorted_bam,
                        target_path=file_path_merge.merged_bam))
                runnable_merge.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_bai',
                        source_path=file_path_index.sorted_bai,
                        target_path=file_path_merge.merged_bai))
                runnable_merge.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_md5',
                        source_path=file_path_index.sorted_md5,
                        target_path=file_path_merge.merged_md5))
            else:
                # Run Picard MergeSamFiles on each BAM file.
                runnable_step = runnable_merge.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_merge_sam_files',
                        java_temporary_path=runnable_merge.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        picard_classpath=self.classpath_picard,
                        picard_command='MergeSamFiles'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                # INPUT []
                for runnable_index in runnable_index_list:
                    file_path_index = runnable_index.file_path_object
                    """ @type file_path_index: FilePathIndex """
                    runnable_step.obsolete_file_path_list.append(file_path_index.sorted_bam)
                    runnable_step.obsolete_file_path_list.append(file_path_index.sorted_bai)
                    runnable_step.obsolete_file_path_list.append(file_path_index.sorted_md5)
                    runnable_step.add_picard_option(key='INPUT', value=file_path_index.sorted_bam, override=True)

                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_merge.merged_bam)
                # SORT_ORDER [coordinate]
                # ASSUME_SORTED [false]
                # MERGE_SEQUENCE_DICTIONARIES [false]
                # USE_THREADING [false]
                runnable_step.add_picard_option(key='USE_THREADING', value='true')
                # COMMENT [null]
                # INTERVALS [null]
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_merge.get_relative_temporary_directory_path)
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET [false]
                # VALIDATION_STRINGENCY [STRICT]
                # COMPRESSION_LEVEL [5]
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                # MAX_RECORDS_IN_RAM [500000]
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                # CREATE_INDEX [false]
                runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                # CREATE_MD5_FILE [false]
                runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                # REFERENCE_SEQUENCE [null]
                # GA4GH_CLIENT_SECRETS [client_secrets.json]
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]

            # Create a symbolic link from the Picard-style *.bai file to a samtools-style *.bam.bai file.

            runnable_merge.add_runnable_step(
                runnable_step=bsf.process.RunnableStepLink(
                    name='link',
                    source_path=file_path_merge.merged_bai_link_source,
                    target_path=file_path_merge.merged_bai_link_target))

            # Run the Picard CollectAlignmentSummaryMetrics analysis.

            runnable_step = runnable_merge.add_runnable_step(
                runnable_step=bsf.process.RunnableStepPicard(
                    name='picard_collect_alignment_summary_metrics',
                    java_temporary_path=runnable_merge.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    picard_classpath=self.classpath_picard,
                    picard_command='CollectAlignmentSummaryMetrics'))
            """ @type runnable_step: bsf.process.RunnableStepPicard """
            # MAX_INSERT_SIZE [100000]
            # EXPECTED_PAIR_ORIENTATIONS [FR]
            # ADAPTER_SEQUENCE [
            #   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT,
            #   AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG,
            #   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT,
            #   AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG,
            #   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT,
            #   AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
            # ]
            # METRIC_ACCUMULATION_LEVEL [ALL_READS]
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
            # IS_BISULFITE_SEQUENCED [false]
            # INPUT []
            runnable_step.add_picard_option(key='INPUT', value=file_path_merge.merged_bam)
            # OUTPUT []
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_merge.alignment_summary_metrics)
            # ASSUME_SORTED [true]
            # STOP_AFTER [0]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_merge.get_relative_temporary_directory_path)
            # VERBOSITY [INFO]
            runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
            # QUIET [false]
            # VALIDATION_STRINGENCY [STRICT]
            # COMPRESSION_LEVEL [5]
            # MAX_RECORDS_IN_RAM [500000]
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            # REFERENCE_SEQUENCE [null]
            runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.genome_fasta)
            # GA4GH_CLIENT_SECRETS [client_secrets.json]
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]

        # Create a Runnable and Executable for the Bowtie2 summary.

        runnable_summary = self.add_runnable(
            runnable=bsf.Runnable(
                name=prefix_summary,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                file_path_object=file_path_summary,
                debug=self.debug))
        executable_summary = self.set_stage_runnable(
            stage=stage_summary,
            runnable=runnable_summary)

        # Add dependencies on Runnable objects of the merging stage.
        for runnable_merge in runnable_merge_list:
            executable_summary.dependencies.append(runnable_merge.name)

        # runnable_summary.add_runnable_step(
        #     runnable_step=bsf.process.RunnableStep(
        #         name='summary',
        #         program='bsf_bowtie2_summary.R'))
        # """ @type runnable_step: bsf.process.RunnableStep """

        return
