# -*- coding: utf-8 -*-
"""Aligner Analysis module

A package of classes and methods supporting Aligner analyses.
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

from __future__ import print_function

import os
import re
import sys
import warnings

import bsf
import bsf.annotation


class FilePathAlign(bsf.FilePath):
    """The C{bsf.analyses.aligner.FilePathAlign} class models file paths at the alignment stage.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar stderr_txt: Text file to capture STDERR of the aligner
    @type stderr_txt: str | unicode
    @ivar stdout_txt: Text file to capture STDOUT of the aligner
    @type stdout_txt: str | unicode
    @ivar aligned_sam: Aligned sequence alignment map (SAM) file path
    @type aligned_sam: str | unicode
    @ivar cleaned_sam: Cleaned sequence alignment map (SAM) file path
    @type cleaned_sam: str | unicode
    @ivar aligned_bai: Aligned binary alignment map index (BAI) file path
    @type aligned_bai: str | unicode
    @ivar aligned_bam: Aligned binary alignment map (BAM) file path
    @type aligned_bam: str | unicode
    @ivar aligned_md5: Aligned binary alignment map (BAM) file path
    @type aligned_md5: str | unicode
    @ivar aligned_bai_link_source: Symbolic link source of the aligned binary alignment map index (BAI) file path
    @type aligned_bai_link_source: str | unicode
    @ivar aligned_bai_link_target: Symbolic link target of the aligned binary alignment map index (BAI) file path
    @type aligned_bai_link_target: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.aligner.FilePathAlign} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.stderr_txt = os.path.join(prefix, '_'.join((prefix, 'stderr.txt')))
        self.stdout_txt = os.path.join(prefix, '_'.join((prefix, 'stdout.txt')))

        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'aligned.sam')))

        self.cleaned_sam = os.path.join(prefix, '_'.join((prefix, 'cleaned.sam')))

        self.aligned_bai = os.path.join(prefix, '_'.join((prefix, 'aligned.bai')))
        self.aligned_bam = os.path.join(prefix, '_'.join((prefix, 'aligned.bam')))
        self.aligned_md5 = os.path.join(prefix, '_'.join((prefix, 'aligned.bam.md5')))

        self.aligned_bai_link_source = '_'.join((prefix, 'aligned.bai'))
        self.aligned_bai_link_target = os.path.join(prefix, '_'.join((prefix, 'aligned.bam.bai')))

        return


class FilePathReadGroup(bsf.FilePath):
    """The C{bsf.analyses.aligner.FilePathReadGroup} class models file paths at the read group processing stage.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar merged_bai: Merged binary alignment map index (BAI) file path
    @type merged_bai: str | unicode
    @ivar merged_bam: Merged binary alignment map (BAM) file path
    @type merged_bam: str | unicode
    @ivar merged_md5: Merged binary alignment map (BAM) file path
    @type merged_md5: str | unicode
    @ivar merged_bai_link_source: Symbolic link source of the merged binary alignment map index (BAI) file path
    @type merged_bai_link_source: str | unicode
    @ivar merged_bai_link_target: Symbolic link target of the merged binary alignment map index (BAI) file path
    @type merged_bai_link_target: str | unicode
    @ivar sorted_bai: Sorted binary alignment map index (BAI) file path
    @type sorted_bai: str | unicode
    @ivar sorted_bam: Sorted binary alignment map (BAM) file path
    @type sorted_bam: str | unicode
    @ivar sorted_md5: Sorted binary alignment map (BAM) file path
    @type sorted_md5: str | unicode
    @ivar sorted_bai_link_source: Symbolic link source of the sorted binary alignment map index (BAI) file path
    @type sorted_bai_link_source: str | unicode
    @ivar sorted_bai_link_target: Symbolic link target of the sorted binary alignment map index (BAI) file path
    @type sorted_bai_link_target: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.aligner.FilePathReadGroup} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathReadGroup, self).__init__(prefix=prefix)

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


class FilePathSample(bsf.FilePath):
    """The C{bsf.analyses.aligner.FilePathSample} class models file paths at the sample processing stage.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar merged_bai: Merged binary alignment map index (BAI) file path
    @type merged_bai: str | unicode
    @ivar merged_bam: Merged binary alignment map (BAM) file path
    @type merged_bam: str | unicode
    @ivar merged_md5: Merged binary alignment map (BAM) file path
    @type merged_md5: str | unicode
    @ivar merged_bai_link_source: Symbolic link source of the merged binary alignment map index (BAI) file path
    @type merged_bai_link_source: str | unicode
    @ivar merged_bai_link_target: Symbolic link target of the merged binary alignment map index (BAI) file path
    @type merged_bai_link_target: str | unicode
    @ivar alignment_summary_metrics_tsv: Picard Alignment Summary Metrics file path
    @type alignment_summary_metrics_tsv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.aligner.FilePathSample} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathSample, self).__init__(prefix=prefix)

        self.output_directory = prefix

        # For the final files, just end in the suffix.
        self.merged_bai = os.path.join(prefix, prefix + '.bai')
        self.merged_bam = os.path.join(prefix, prefix + '.bam')
        self.merged_md5 = os.path.join(prefix, prefix + '.bam.md5')

        self.merged_bai_link_source = prefix + '.bai'
        self.merged_bai_link_target = os.path.join(prefix, prefix + '.bam.bai')

        self.alignment_summary_metrics_tsv = os.path.join(prefix, '_'.join((prefix, 'alignment_summary_metrics.tsv')))

        return


class FilePathSummary(bsf.FilePath):
    """The C{bsf.analyses.aligner.FilePathSummary} class models file paths at the summary stage.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.aligner.FilePathSummary} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathSummary, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.read_group_to_sample_tsv = prefix + '_read_group_to_sample.tsv'

        return


class Aligner(bsf.Analysis):
    """The C{bsf.analyses.aligner.Aligner} class represents the logic to run a (short read) aligner.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar genome_fasta: Genome FASTA file
    @type genome_fasta: str | unicode | None
    @ivar genome_index: Genome index
    @type genome_index: str | unicode | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: None | str | unicode
    """

    name = 'Aligner Analysis'
    prefix = 'aligner'

    @classmethod
    def get_stage_name_align(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'align'))

    @classmethod
    def get_stage_name_index(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'index'))

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_stage_name_sample(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'sample'))

    @classmethod
    def get_stage_name_summary(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'summary'))

    @classmethod
    def get_prefix_align(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_align(), paired_reads_name))

    @classmethod
    def get_prefix_index(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_index(), paired_reads_name))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_prefix_sample(cls, sample_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_sample(), sample_name))

    @classmethod
    def get_prefix_summary(cls):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return cls.get_stage_name_summary()

    @classmethod
    def get_file_path_align(cls, paired_reads_name):
        """Get a C{FilePathAlign} object from this or a sub-class.

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: C{FilePathAlign} or sub-class object
        @rtype: FilePathAlign
        """
        return FilePathAlign(prefix=cls.get_prefix_align(paired_reads_name=paired_reads_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name):
        """Get a C{FilePathReadGroup} object from this or a sub-class.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: C{FilePathReadGroup} or sub-class object
        @rtype: FilePathReadGroup
        """
        return FilePathReadGroup(prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

    @classmethod
    def get_file_path_sample(cls, sample_name):
        """Get a C{FilePathSample} object from this or a sub-class.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: C{FilePathSample} or sub-class object
        @rtype: FilePathSample
        """
        return FilePathSample(prefix=cls.get_prefix_sample(sample_name=sample_name))

    @classmethod
    def get_file_path_summary(cls):
        """Get a C{FilePathSummary} object from this or a sub-class.

        @return: C{FilePathSummary} or sub-class object
        @rtype: FilePathSummary
        """
        return FilePathSummary(prefix=cls.get_prefix_summary())

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
        """Initialise a C{bsf.analyses.aligner.Aligner}.

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
        @param genome_index: Genome index
        @type genome_index: str | unicode | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: None | str | unicode
        @return:
        @rtype:
        """
        super(Aligner, self).__init__(
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
        """Set instance variables of a C{bsf.analyses.aligner.Aligner} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(Aligner, self).set_configuration(configuration=configuration, section=section)

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

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add one or more Aligner-specific C{bsf.process.RunnableStep} objects to the C{bsf.Runnable}.

        @param runnable_align: C{bsf.Runnable}
        @type runnable_align: bsf.Runnable
        @param stage_align: C{bsf.Stage}
        @type stage_align: bsf.Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | unicode | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | unicode | None
        @return:
        @rtype:
        """
        return

    def add_runnable_step_summary(self, runnable_summary, stage_summary):
        """Add one or more Aligner-specific C{bsf.process.RunnableStep} objects to the C{bsf.Runnable}.

        @param runnable_summary: C{bsf.Runnable}
        @type runnable_summary: bsf.Runnable
        @param stage_summary: C{bsf.Stage}
        @type stage_summary: bsf.Stage
        @return:
        @rtype:
        """
        return

    def run(self):
        """Run a C{bsf.analyses.aligner.Aligner} analysis.

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

        def run_get_unmapped_bam(_paired_reads):
            """Get the unmapped BAM file annotation of a C{bsf.ngs.PairedReads} object.

            @param _paired_reads: C{bsf.ngs.PairedReads} object
            @type _paired_reads: bsf.ngs.PairedReads
            @return: Python tuple of Python str (base name) and Python str or unicode (unmapped BAM file path)
            @rtype: (str, str | unicode | None)
            """
            if 'BAM File' in _paired_reads.annotation_dict and _paired_reads.annotation_dict['BAM File']:
                # Consider only the first list component.
                bam_name, bam_extension = os.path.splitext(
                    os.path.basename(_paired_reads.annotation_dict['BAM File'][0]))

                # TODO: This should be centralised.
                # The makeFileNameSafe() method of htsjdk.samtools.util.IOUtil uses the following pattern:
                # [\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]
                bam_name = re.sub(
                    pattern='[\\s!"#$%&\'()*/:;<=>?@\\[\\]\\\\^`{|}~]',
                    repl='_',
                    string=bam_name)

                return bam_name, _paired_reads.annotation_dict['BAM File'][0]
            else:
                warnings.warn("PairedReads object '" + _paired_reads.get_name() + "without 'BAM File' annotation.")
                return _paired_reads.get_name(), None

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the bsf.Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(path=self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[self.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # The Aligner analysis requires a genome version.

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a 'genome_version' configuration option.")

        super(Aligner, self).run()

        if not self.genome_fasta:
            # Set the genome_index to None, to point to the genome directory that contains the
            # Picard sequence dictionary.
            self.genome_fasta = bsf.standards.FilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index=None)

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('A ' + self.name + " requires a 'classpath_picard' configuration option.")

        run_read_comparisons()

        stage_align = self.get_stage(name=self.get_stage_name_align())
        stage_index = self.get_stage(name=self.get_stage_name_index())
        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())
        stage_sample = self.get_stage(name=self.get_stage_name_sample())
        stage_summary = self.get_stage(name=self.get_stage_name_summary())

        file_path_summary = self.get_file_path_summary()

        # Create an annotation sheet linking sample name and read group name, which is required for the
        # summary script.

        annotation_sheet = bsf.annotation.AnnotationSheet(
            file_path=os.path.join(self.genome_directory, file_path_summary.read_group_to_sample_tsv),
            file_type='excel-tab',
            name='star_aligner_read_group',
            field_names=['sample', 'read_group'])

        runnable_sample_list = list()
        """ @type runnable_sample_list: list[bsf.Runnable] """

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', repr(sample.name))
                sys.stdout.writelines(sample.trace(level=1))

            # To run Picard MergeBamAlignment, all alignments from a BAM file need merging into one.

            unmapped_bam_file_dict = dict()
            """ @type unmapped_bam_file_dict: dict[str | unicode, (str | unicode, list[bsf.Runnable])] """

            runnable_read_group_list = list()
            """ @type runnable_read_group_list: list[bsf.Runnable] """

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
                paired_reads = paired_reads_dict[paired_reads_name][0]

                # Get the file paths for Reads1 and Reads2 and check for FASTQ files.
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

                ###################
                # Alignment Stage #
                ###################

                annotation_sheet.row_dicts.append({'sample': sample.name, 'read_group': paired_reads.get_name()})

                # Create a Runnable and Executable for alignment and processing.

                file_path_align = self.get_file_path_align(paired_reads_name=paired_reads_name)

                runnable_align = self.add_runnable(
                    runnable=bsf.Runnable(
                        name=self.get_prefix_align(paired_reads_name=paired_reads_name),
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        file_path_object=file_path_align,
                        debug=self.debug))
                self.set_stage_runnable(
                    stage=stage_align,
                    runnable=runnable_align)

                # Make a directory.

                runnable_align.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_align.output_directory))

                # Run the Aligner.

                self.add_runnable_step_aligner(
                    runnable_align=runnable_align,
                    stage_align=stage_align,
                    file_path_1=file_path_1,
                    file_path_2=file_path_2)

                ##################
                # Indexing Stage #
                ##################

                # Create a Runnable and Executable for indexing each read group.

                prefix_index = self.get_prefix_index(paired_reads_name=paired_reads_name)

                # file_path_index = FilePathIndex(prefix=prefix_index)

                runnable_index = self.add_runnable(
                    runnable=bsf.Runnable(
                        name=prefix_index,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        file_path_object=file_path_align,
                        debug=self.debug))
                executable_index = self.set_stage_runnable(
                    stage=stage_index,
                    runnable=runnable_index)
                executable_index.dependencies.append(runnable_align.name)

                # Record the indexing Runnable under the unmapped BAM file name
                # to merge the sub-alignments of all (trimmed) FASTQ files that
                # resulted from an unmapped BAM file.
                bam_file_name, bam_file_path = run_get_unmapped_bam(_paired_reads=paired_reads)

                if bam_file_name in unmapped_bam_file_dict:
                    unmapped_bam_file_dict[bam_file_name][1].append(runnable_index)
                else:
                    unmapped_bam_file_dict[bam_file_name] = (bam_file_path, [runnable_index])

                # Run Picard CleanSam to convert the aligned SAM file into a cleaned SAM file.

                runnable_step = runnable_index.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_clean_sam',
                        obsolete_file_path_list=[
                            file_path_align.aligned_sam,
                        ],
                        java_temporary_path=runnable_index.temporary_directory_path(absolute=False),
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
                    value=runnable_index.temporary_directory_path(absolute=False))
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

                runnable_step = runnable_index.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_sort_sam',
                        obsolete_file_path_list=[
                            file_path_align.cleaned_sam,
                        ],
                        java_temporary_path=runnable_index.temporary_directory_path(absolute=False),
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
                    value=runnable_index.temporary_directory_path(absolute=False))
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

            ###############################
            # Read Group Processing Stage #
            ###############################

            # Merge all aligned BAM files of each read group of an unaligned BAM file.

            for bam_file_name, (bam_file_path, runnable_index_list) in unmapped_bam_file_dict.iteritems():
                # Create a Runnable and Executable for merging each read group.

                file_path_read_group = self.get_file_path_read_group(read_group_name=bam_file_name)

                runnable_read_group = self.add_runnable(
                    runnable=bsf.Runnable(
                        name=self.get_prefix_read_group(read_group_name=bam_file_name),
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        file_path_object=file_path_read_group,
                        debug=self.debug))
                executable_read_group = self.set_stage_runnable(
                    stage=stage_read_group,
                    runnable=runnable_read_group)

                # Record the dependency for all alignment Runnable objects.
                for runnable_index in runnable_index_list:
                    executable_read_group.dependencies.append(runnable_index.name)

                runnable_read_group_list.append(runnable_read_group)

                runnable_read_group.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_read_group.output_directory))

                if len(runnable_index_list) == 1:
                    runnable_index = runnable_index_list[0]
                    file_path_align = runnable_index.file_path_object
                    """ @type file_path_align: FilePathAlign """
                    # For a single ReadPair, just rename the files.
                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_bai',
                            source_path=file_path_align.aligned_bai,
                            target_path=file_path_read_group.merged_bai))
                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_bam',
                            source_path=file_path_align.aligned_bam,
                            target_path=file_path_read_group.merged_bam))
                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_md5',
                            source_path=file_path_align.aligned_md5,
                            target_path=file_path_read_group.merged_md5))
                else:
                    runnable_step = runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepPicard(
                            name='picard_merge_sam_files',
                            java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx4G',
                            picard_classpath=self.classpath_picard,
                            picard_command='MergeSamFiles'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    # INPUT []
                    for runnable_index in runnable_index_list:
                        file_path_align = runnable_index.file_path_object
                        """ @type file_path_align: FilePathAlign """
                        runnable_step.add_picard_option(
                            key='INPUT',
                            value=file_path_align.aligned_bam,
                            override=True)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_bai)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_bam)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_md5)
                    # OUTPUT []
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.merged_bam)
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
                        value=runnable_read_group.temporary_directory_path(absolute=False))
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

                if bam_file_path is None:
                    # Without an unmapped BAM file, just rename the files.
                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_bai',
                            source_path=file_path_read_group.merged_bai,
                            target_path=file_path_read_group.sorted_bai))
                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_bam',
                            source_path=file_path_read_group.merged_bam,
                            target_path=file_path_read_group.sorted_bam))
                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMove(
                            name='move_md5',
                            source_path=file_path_read_group.merged_md5,
                            target_path=file_path_read_group.sorted_md5))
                else:
                    # Run Picard MergeBamAlignment to annotate the merged aligned with the unaligned BAM file.

                    runnable_step = runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepPicard(
                            name='picard_merge_bam_alignment',
                            obsolete_file_path_list=[
                                file_path_read_group.merged_bai,
                                file_path_read_group.merged_bam,
                                file_path_read_group.merged_md5,
                            ],
                            java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx12G',
                            java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                            picard_command='MergeBamAlignment'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    # UNMAPPED_BAM []
                    runnable_step.add_picard_option(key='UNMAPPED_BAM', value=bam_file_path)
                    # ALIGNED_BAM [null]
                    runnable_step.add_picard_option(key='ALIGNED_BAM', value=file_path_read_group.merged_bam)
                    # READ1_ALIGNED_BAM [null]
                    # READ2_ALIGNED_BAM [null]
                    # OUTPUT []
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.sorted_bam)
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
                        value=runnable_read_group.temporary_directory_path(absolute=False))
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

            ###########################
            # Sample Processing Stage #
            ###########################

            # For more than one ReadPair object, the aligned BAM files need merging into Sample-specific ones.

            file_path_sample = self.get_file_path_sample(sample_name=sample.name)

            runnable_sample = self.add_runnable(
                runnable=bsf.Runnable(
                    name=self.get_prefix_sample(sample_name=sample.name),
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    file_path_object=file_path_sample,
                    debug=self.debug))
            executable_sample = self.set_stage_runnable(
                stage=stage_sample,
                runnable=runnable_sample)

            runnable_sample_list.append(runnable_sample)

            # Add dependencies on Runnable objects of the read group processing stage.

            for runnable_read_group in runnable_read_group_list:
                executable_sample.dependencies.append(runnable_read_group.name)

            runnable_sample.add_runnable_step(
                runnable_step=bsf.process.RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=file_path_sample.output_directory))

            if len(runnable_read_group_list) == 1:
                runnable_read_group = runnable_read_group_list[0]
                file_path_read_group = runnable_read_group.file_path_object
                """ @type file_path_read_group: FilePathReadGroup """
                # For a single ReadPair, just rename the files.
                runnable_sample.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_bai',
                        source_path=file_path_read_group.sorted_bai,
                        target_path=file_path_sample.merged_bai))
                runnable_sample.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_bam',
                        source_path=file_path_read_group.sorted_bam,
                        target_path=file_path_sample.merged_bam))
                runnable_sample.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepMove(
                        name='move_md5',
                        source_path=file_path_read_group.sorted_md5,
                        target_path=file_path_sample.merged_md5))
            else:
                # Run Picard MergeSamFiles on each BAM file.
                runnable_step = runnable_sample.add_runnable_step(
                    runnable_step=bsf.process.RunnableStepPicard(
                        name='picard_merge_sam_files',
                        java_temporary_path=runnable_sample.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx4G',
                        picard_classpath=self.classpath_picard,
                        picard_command='MergeSamFiles'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                # INPUT []
                for runnable_read_group in runnable_read_group_list:
                    file_path_read_group = runnable_read_group.file_path_object
                    """ @type file_path_read_group: FilePathReadGroup """
                    runnable_step.obsolete_file_path_list.append(file_path_read_group.sorted_bai)
                    runnable_step.obsolete_file_path_list.append(file_path_read_group.sorted_bam)
                    runnable_step.obsolete_file_path_list.append(file_path_read_group.sorted_md5)
                    runnable_step.add_picard_option(key='INPUT', value=file_path_read_group.sorted_bam, override=True)

                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_sample.merged_bam)
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
                    value=runnable_sample.temporary_directory_path(absolute=False))
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

            runnable_sample.add_runnable_step(
                runnable_step=bsf.process.RunnableStepLink(
                    name='link',
                    source_path=file_path_sample.merged_bai_link_source,
                    target_path=file_path_sample.merged_bai_link_target))

            # Run the Picard CollectAlignmentSummaryMetrics analysis.

            runnable_step = runnable_sample.add_runnable_step(
                runnable_step=bsf.process.RunnableStepPicard(
                    name='picard_collect_alignment_summary_metrics',
                    java_temporary_path=runnable_sample.temporary_directory_path(absolute=False),
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
            runnable_step.add_picard_option(key='INPUT', value=file_path_sample.merged_bam)
            # OUTPUT []
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_sample.alignment_summary_metrics_tsv)
            # ASSUME_SORTED [true]
            # STOP_AFTER [0]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_sample.temporary_directory_path(absolute=False))
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

        #################
        # Summary Stage #
        #################

        # Write the AnnotationSheet to disk.

        annotation_sheet.to_file_path()

        # Create a Runnable and Executable for the summary.

        runnable_summary = self.add_runnable(
            runnable=bsf.Runnable(
                name=self.get_prefix_summary(),
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                file_path_object=file_path_summary,
                debug=self.debug))
        executable_summary = self.set_stage_runnable(
            stage=stage_summary,
            runnable=runnable_summary)

        # Add dependencies on Runnable objects of the merging stage.
        for runnable_sample in runnable_sample_list:
            executable_summary.dependencies.append(runnable_sample.name)

        self.add_runnable_step_summary(runnable_summary=runnable_summary, stage_summary=stage_summary)

        return
