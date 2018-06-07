"""bsf.analyses.star_aligner

A package of classes and methods supporting the spliced Transcripts Alignment
to a Reference (STAR) aligner by Alexander Dobin.

Project:  https://github.com/alexdobin/STAR
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

import pysam

from bsf import Analysis, FilePath, Runnable
from bsf.annotation import AnnotationSheet
from bsf.process import RunnableStep, RunnableStepLink, RunnableStepMove, RunnableStepPicard
from bsf.standards import JavaClassPath


class FilePathStarAlign(FilePath):
    def __init__(self, prefix):
        super(FilePathStarAlign, self).__init__(prefix=prefix)

        self.aligned_sam = prefix + '_Aligned.out.sam'
        self.splice_junctions_tsv = prefix + '_SJ.out.tab'

        return


class FilePathStarIndex(FilePath):
    def __init__(self, prefix):
        super(FilePathStarIndex, self).__init__(prefix=prefix)

        self.aligned_bam = prefix + '_Aligned.bam'
        self.aligned_bai = prefix + '_Aligned.bai'
        self.aligned_md5 = prefix + '_Aligned.bam.md5'
        self.cleaned_sam = prefix + '_Cleaned.sam'

        return


class FilePathStarMerge(FilePath):
    def __init__(self, prefix):
        super(FilePathStarMerge, self).__init__(prefix=prefix)

        self.merged_bam = prefix + '.bam'
        self.merged_bai = prefix + '.bai'
        self.merged_lnk = prefix + '.bam.bai'
        self.merged_md5 = prefix + '.bam.md5'
        self.merged_tsv = prefix + '.tsv'
        self.merged_pdf = prefix + '.pdf'

        return


class FilePathStarSummary(FilePath):
    def __init__(self, prefix):
        super(FilePathStarSummary, self).__init__(prefix=prefix)

        self.read_group_to_sample_tsv = prefix + '_read_group_to_sample.tsv'
        self.alignment_read_group_pdf = prefix + '_alignment_read_group.pdf'
        self.alignment_read_group_png = prefix + '_alignment_read_group.png'
        self.alignment_sample_png = prefix + '_alignment_sample.png'
        self.alignment_sample_pdf = prefix + '_alignment_sample.pdf'
        self.junction_fraction_read_group_pdf = prefix + '_junction_fraction_read_group.pdf'
        self.junction_fraction_read_group_png = prefix + '_junction_fraction_read_group.png'
        self.junction_fraction_sample_pdf = prefix + '_junction_fraction_sample.pdf'
        self.junction_fraction_sample_png = prefix + '_junction_fraction_sample.png'
        self.junction_number_read_group_pdf = prefix + '_junction_number_read_group.pdf'
        self.junction_number_read_group_png = prefix + '_junction_number_read_group.png'
        self.junction_number_sample_pdf = prefix + '_junction_number_sample.pdf'
        self.junction_number_sample_png = prefix + '_junction_number_sample.png'
        self.mapped_fraction_read_group_pdf = prefix + '_mapped_fraction_read_group.pdf'
        self.mapped_fraction_read_group_png = prefix + '_mapped_fraction_read_group.png'
        self.mapped_fraction_sample_png = prefix + '_mapped_fraction_sample.png'
        self.mapped_fraction_sample_pdf = prefix + '_mapped_fraction_sample.pdf'
        self.mapped_number_read_group_pdf = prefix + '_mapped_number_read_group.pdf'
        self.mapped_number_read_group_png = prefix + '_mapped_number_read_group.png'
        self.mapped_number_sample_png = prefix + '_mapped_number_sample.png'
        self.mapped_number_sample_pdf = prefix + '_mapped_number_sample.pdf'
        self.table_read_group = prefix + '_table_read_group.tsv'
        self.table_sample = prefix + '_table_sample.tsv'

        return


class StarAligner(Analysis):
    """STAR Aligner C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_align: C{bsf.Stage.name} for the alignment stage
    @type stage_name_align: str
    @cvar stage_name_index: C{bsf.Stage.name} for the indexing stage
    @type stage_name_index: str
    @cvar stage_name_merge: C{bsf.Stage.name} for the merging stage
    @type stage_name_merge: str
    @cvar stage_name_summary: C{bsf.Stage.name} for the summary stage
    @type stage_name_summary: str
    @ivar replicate_grouping: Group all replicates into a single STAR process
    @type replicate_grouping: bool
    @ivar index_directory: Genome directory with STAR indices
    @type index_directory: str | unicode
    @ivar stranded: Stranded RNA-seq protocol 'yes', 'no' 'reverse'
    @type stranded: str
    @ivar transcriptome_gtf: GTF file path of transcriptome annotation
    @type transcriptome_gtf: str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    """

    name = 'STAR Aligner Analysis'
    prefix = 'star_aligner'

    stage_name_align = '_'.join((prefix, 'align'))
    stage_name_index = '_'.join((prefix, 'index'))
    stage_name_merge = '_'.join((prefix, 'merge'))
    stage_name_summary = '_'.join((prefix, 'summary'))

    @classmethod
    def get_prefix_star_aligner_align(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: Replicate key
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_align, paired_reads_name))

    @classmethod
    def get_prefix_star_aligner_index(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: Replicate key
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_index, paired_reads_name))

    @classmethod
    def get_prefix_star_aligner_merge(cls, sample_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param sample_name: Sample name
        @type sample_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_merge, sample_name))

    @classmethod
    def get_prefix_star_aligner_summary(cls, sample_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param sample_name: Sample name
        @type sample_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_summary, sample_name))

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
            replicate_grouping=False,
            index_directory=None,
            stranded=None,
            transcriptome_gtf=None,
            classpath_picard=None):
        """Initialise a C{bsf.analyses.rna_seq.StarAligner} object.

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
        @param replicate_grouping: Group all replicates into a single STAR process
        @type replicate_grouping: bool
        @param index_directory: Genome directory with STAR indices
        @type index_directory: str | unicode
        @param stranded: Stranded RNA-seq protocol 'yes', 'no' 'reverse'
        @type stranded: str
        @param transcriptome_gtf: GTF file path of transcriptome annotation
        @type transcriptome_gtf: str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        """

        super(StarAligner, self).__init__(
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

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if index_directory is None:
            self.index_directory = str()
        else:
            self.index_directory = index_directory

        if stranded is None:
            self.stranded = str()
        else:
            self.stranded = stranded

        if transcriptome_gtf is None:
            self.transcriptome_gtf = str()
        else:
            self.transcriptome_gtf = transcriptome_gtf

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.rna_seq.StarAligner} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(StarAligner, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'index_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.index_directory = configuration.config_parser.get(section=section, option=option)

        option = 'stranded'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf = configuration.config_parser.get(section=section, option=option)

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.StarAligner} analysis.

        Although the STAR aligner can directly count reads according to its splice junction database,
        more than one read group may need aligning so that the count tables had to be combined.
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

        super(StarAligner, self).run()

        # Get global defaults.

        # The STAR Aligner requires a genome version.

        if not self.genome_version:
            raise Exception('A STAR Aligner analysis requires a genome_version configuration option.')

        if not self.index_directory:
            raise Exception('A STAR Aligner analysis requires an index_directory configuration option.')

        if self.stranded:
            if self.stranded not in ('yes', 'no', 'reverse'):
                raise Exception(
                    'The STAR Aligner configuration option "stranded" can only be '
                    '"yes", "no" or "reverse", not "' + self.stranded + '".')
        else:
            self.stranded = 'yes'

        if not self.transcriptome_gtf:
            raise Exception('A STAR Aligner analysis requires a transcriptome_gtf configuration option.')

        if not self.classpath_picard:
            self.classpath_picard = JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception("An 'StarAligner' analysis requires a "
                                "'classpath_picard' configuration option.")

        run_read_comparisons()

        stage_align = self.get_stage(name=self.stage_name_align)
        stage_index = self.get_stage(name=self.stage_name_index)
        stage_merge = self.get_stage(name=self.stage_name_merge)
        stage_summary = self.get_stage(name=self.stage_name_summary)

        prefix_summary = stage_summary.name

        file_path_summary = FilePathStarSummary(prefix=prefix_summary)

        # Create an annotation sheet linking sample name and read group name, which is required for the
        # summary script.

        annotation_sheet = AnnotationSheet(
            file_path=os.path.join(self.genome_directory, file_path_summary.read_group_to_sample_tsv),
            file_type='excel-tab',
            name='star_aligner_read_group',
            field_names=['sample', 'read_group'])

        runnable_merge_list = list()
        """ @type runnable_merge_list: list[bsf.Runnable] """

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                print(sample.trace(1))

            runnable_index_list = list()
            """ @type runnable_index_list: list[bsf.Runnable] """

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    ######################
                    # 1. Alignment Stage #
                    ######################

                    prefix_align = '_'.join((stage_align.name, paired_reads.get_name()))

                    # STAR-specific file paths

                    file_path_align = FilePathStarAlign(prefix=prefix_align)

                    # Create a Runnable and Executable for the STAR aligner.

                    runnable_align = self.add_runnable(
                        runnable=Runnable(
                            name=prefix_align,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            cache_directory=self.cache_directory,
                            file_path_object=file_path_align,
                            debug=self.debug))
                    executable_align = self.set_stage_runnable(
                        stage=stage_align,
                        runnable=runnable_align)

                    runnable_step = runnable_align.add_runnable_step(
                        runnable_step=RunnableStep(
                            name='STAR',
                            program='STAR'))
                    """ @type runnable_step: bsf.process.RunnableStep """
                    self.set_runnable_step_configuration(runnable_step=runnable_step)
                    runnable_step.add_option_long(key='runThreadN', value=str(stage_align.threads))
                    runnable_step.add_option_long(key='genomeDir', value=self.index_directory)
                    runnable_step.add_option_long(key='outFileNamePrefix', value=prefix_align + '_')
                    # NOTE: The STAR aligner command line interface is seriously broken,
                    # as the readFilesIn option requires two values.
                    # Hence, use class bsf.argument.OptionMultiLong via wrapper Command.add_option_multi_long().
                    if paired_reads.reads_2 is None:
                        runnable_step.add_option_long(
                            key='readFilesIn',
                            value=paired_reads.reads_1.file_path)
                    else:
                        runnable_step.add_option_multi_long(
                            key='readFilesIn',
                            value=' '.join((paired_reads.reads_1.file_path, paired_reads.reads_2.file_path)))
                    if paired_reads.reads_1.file_path.endswith('fastq.gz'):
                        runnable_step.add_option_long(key='readFilesCommand', value='zcat')

                    # If the original BAM file was annotated in the PairedReads object, the read group can be set.

                    if 'BAM File' in paired_reads.annotation_dict:
                        # There should be only one BAM file linked to the PairedReads object.
                        bam_file_path = paired_reads.annotation_dict['BAM File'][0]
                        if os.path.exists(bam_file_path):
                            alignment_file = pysam.AlignmentFile(bam_file_path, 'rb', check_sq=False)

                            # Add the @RG line.
                            for read_group_dict in alignment_file.header['RG']:
                                """ @type read_group_dict: dict[str, str] """
                                # There should also be only one read group.
                                key_list = read_group_dict.keys()
                                # Remove the 'ID' from the key_list, as it has to go first.
                                key_list.remove('ID')
                                read_group_str = ':'.join(('ID', read_group_dict['ID']))
                                for key in key_list:
                                    read_group_str += ' "' + ':'.join((key, read_group_dict[key])) + '"'
                                runnable_step.add_option_long(key='outSAMattrRGline', value=read_group_str)

                            # Add @PG lines.
                            # The STAR aligner allows only a single @PG line, which is not terribly helpful.
                            # for program_dict in alignment_file.header['PG']:
                            #     """ @type program_dict: dict[str, str] """
                            #     key_list = program_dict.keys()
                            #     program_str = ''
                            #     for key in key_list:
                            #         program_str += ' "' + ':'.join((key, program_dict[key])) + '"'
                            #     runnable_step.add_option_long(key='outSAMheaderPG', value=program_str)

                    annotation_sheet.row_dicts.append({'sample': sample.name, 'read_group': paired_reads.get_name()})

                    #####################
                    # 2. Indexing Stage #
                    #####################

                    prefix_index = '_'.join((stage_index.name, paired_reads.get_name()))

                    file_path_index = FilePathStarIndex(prefix=prefix_index)

                    runnable_index = self.add_runnable(
                        runnable=Runnable(
                            name=prefix_index,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            cache_directory=self.cache_directory,
                            file_path_object=file_path_index,
                            debug=self.debug))
                    executable_index = self.set_stage_runnable(
                        stage=stage_index,
                        runnable=runnable_index)
                    executable_index.dependencies.append(executable_align.name)

                    runnable_index_list.append(runnable_index)

                    runnable_step = runnable_index.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='picard_clean_sam',
                            obsolete_file_path_list=[file_path_align.aligned_sam],
                            java_temporary_path=runnable_index.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx2G',
                            picard_classpath=self.classpath_picard,
                            picard_command='CleanSam'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    runnable_step.add_picard_option(key='INPUT', value=file_path_align.aligned_sam)
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_index.cleaned_sam)
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_index.get_relative_temporary_directory_path)
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    runnable_step.add_picard_option(key='QUIET', value='false')
                    runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')

                    runnable_step = runnable_index.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='picard_sort_sam',
                            obsolete_file_path_list=[file_path_index.cleaned_sam],
                            java_temporary_path=runnable_index.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            picard_classpath=self.classpath_picard,
                            picard_command='SortSam'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    runnable_step.add_picard_option(key='INPUT', value=file_path_index.cleaned_sam)
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_index.aligned_bam)
                    runnable_step.add_picard_option(key='SORT_ORDER', value='coordinate')
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_index.get_relative_temporary_directory_path)
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    runnable_step.add_picard_option(key='QUIET', value='false')
                    runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')
                    runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                    runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                    runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')

                    # Run GNU Zip over the rather large splice junction table.

                    runnable_step = runnable_index.add_runnable_step(
                        runnable_step=RunnableStep(
                            name='gzip',
                            program='gzip'))
                    """ @type runnable_step: bsf.process.RunnableStep """
                    runnable_step.add_switch_long(key='best')
                    runnable_step.arguments.append(file_path_align.splice_junctions_tsv)

            ####################
            # 3. Merging Stage #
            ####################

            # For more than one ReadPair object the aligned BAM files need merging into Sample-specific ones.

            prefix_merge = '_'.join((stage_merge.name, sample.name))

            file_path_merge = FilePathStarMerge(prefix=prefix_merge)

            runnable_merge = self.add_runnable(
                runnable=Runnable(
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

            if len(runnable_index_list) == 1:
                runnable_index = runnable_index_list[0]
                file_path_index = runnable_index.file_path_object
                """ @type file_path_index: FilePathStarIndex """
                # For a single ReadPair, just rename the files.
                runnable_merge.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_bam',
                        source_path=file_path_index.aligned_bam,
                        target_path=file_path_merge.merged_bam))
                runnable_merge.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_bai',
                        source_path=file_path_index.aligned_bai,
                        target_path=file_path_merge.merged_bai))
                runnable_merge.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='move_md5',
                        source_path=file_path_index.aligned_md5,
                        target_path=file_path_merge.merged_md5))
            else:
                # Run Picard MergeSamFiles on each BAM file.
                runnable_step = runnable_merge.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='picard_merge_sam_files',
                        java_temporary_path=runnable_merge.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx2G',
                        picard_classpath=self.classpath_picard,
                        picard_command='MergeSamFiles'))
                """ @type runnable_step: bsf.process.RunnableStepPicard """
                for runnable_index in runnable_index_list:
                    file_path_index = runnable_index.file_path_object
                    """ @type file_path_index: FilePathStarIndex """
                    runnable_step.obsolete_file_path_list.append(file_path_index.aligned_bam)
                    runnable_step.obsolete_file_path_list.append(file_path_index.aligned_bai)
                    runnable_step.obsolete_file_path_list.append(file_path_index.aligned_md5)
                    runnable_step.add_picard_option(key='INPUT', value=file_path_index.aligned_bam, override=True)

                runnable_step.add_picard_option(key='OUTPUT', value=file_path_merge.merged_bam)
                runnable_step.add_picard_option(key='SORT_ORDER', value='coordinate')
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_merge.get_relative_temporary_directory_path)
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                runnable_step.add_picard_option(key='QUIET', value='false')
                runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')

            # Create a symbolic link from the Picard-style *.bai file to a samtools-style *.bam.bai file.

            runnable_merge.add_runnable_step(
                runnable_step=RunnableStepLink(
                    name='link',
                    source_path=file_path_merge.merged_bai,
                    target_path=file_path_merge.merged_lnk))

        # Write the AnnotationSheet to disk.

        annotation_sheet.to_file_path()

        # Create a Runnable and Executable for the STAR summary.

        runnable_summary = self.add_runnable(
            runnable=Runnable(
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

        runnable_summary.add_runnable_step(
            runnable_step=RunnableStep(
                name='summary',
                program='bsf_star_aligner_summary.R'))
        """ @type runnable_step: bsf.process.RunnableStep """

        return

    def report(self):
        """Create a report.

        @return:
        @rtype:
        """

        def report_html():
            """Private function to create a HTML report.

            @return:
            @rtype:
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str | unicode] """

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.append('<p id="ucsc_track_hub">')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<h2 id="sample_section">Sample Table</h2>\n')
            str_list.append('\n')
            str_list.append('<table id="sample_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Sample</th>\n')
            str_list.append('<th>BAM</th>\n')
            str_list.append('<th>BAI</th>\n')
            str_list.append('<th>MD5</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping,
                                                                exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                runnable_merge = self.runnable_dict[
                    '_'.join((self.stage_name_merge, sample.name))]
                file_path_merge = runnable_merge.file_path_object
                """ @type file_path_merge: FilePathStarMerge """

                str_list.append('<tr>\n')
                # Sample
                str_list.append('<td class="left">' + sample.name + '</td>\n')
                # BAM
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_merge.merged_bam + '">')
                str_list.append('<abbr title="Binary Alignment/Map">BAM</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # BAI
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_merge.merged_bai + '">')
                str_list.append('<abbr title="Binary Alignment/Map Index">BAI</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # MD5
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_merge.merged_md5 + '">')
                str_list.append('<abbr title="Message Digest 5 Checksum">MD5</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            str_list.append('<h2 id="qc_plots">QC Plots</h2>\n')
            str_list.append('\n')
            str_list.append('<table id="qc_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Sample</th>\n')
            str_list.append('<th>Read Group</th>\n')
            str_list.append('<th>Metrics</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            runnable_summary = self.runnable_dict[self.stage_name_summary]
            file_path_summary = runnable_summary.file_path_object
            """ @type file_path_summary: FilePathStarSummary """

            # Alignment Summary Plots
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.alignment_sample_pdf + '">')
            str_list.append('<img alt="Mapped - Sample"')
            str_list.append(' src="' + file_path_summary.alignment_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.alignment_read_group_pdf + '">')
            str_list.append('<img alt="Mapped - Read Group"')
            str_list.append(' src="' + file_path_summary.alignment_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Mapped</td>\n')
            str_list.append('</tr>\n')

            # Mapped Fraction Plots
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.mapped_fraction_sample_pdf + '">')
            str_list.append('<img alt="Mapped Fraction - Sample"')
            str_list.append(' src="' + file_path_summary.mapped_fraction_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.mapped_fraction_read_group_pdf + '">')
            str_list.append('<img alt="Mapped Fraction - Read Group"')
            str_list.append(' src="' + file_path_summary.mapped_fraction_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Mapped Fraction</td>\n')
            str_list.append('</tr>\n')

            # Mapped Number Plots
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.mapped_number_sample_pdf + '">')
            str_list.append('<img alt="Mapped Number - Sample"')
            str_list.append(' src="' + file_path_summary.mapped_number_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.mapped_number_read_group_pdf + '">')
            str_list.append('<img alt="Mapped Number - Read Group"')
            str_list.append(' src="' + file_path_summary.mapped_number_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Mapped Number</td>\n')
            str_list.append('</tr>\n')

            # Junction Fraction Plots
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.junction_fraction_sample_pdf + '">')
            str_list.append('<img alt="Junction Fraction - Sample"')
            str_list.append(' src="' + file_path_summary.junction_fraction_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.junction_fraction_read_group_pdf + '">')
            str_list.append('<img alt="Junction Fraction - Read Group"')
            str_list.append(' src="' + file_path_summary.junction_fraction_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Junction Fraction</td>\n')
            str_list.append('</tr>\n')

            # Junction Number Plots
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.junction_number_sample_pdf + '">')
            str_list.append('<img alt="Junction Number - Sample"')
            str_list.append(' src="' + file_path_summary.junction_number_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.junction_number_read_group_pdf + '">')
            str_list.append('<img alt="Junction Number - Read Group"')
            str_list.append(' src="' + file_path_summary.junction_number_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Junction Number</td>\n')
            str_list.append('</tr>\n')

            # Summary Tables
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.table_sample + '">')
            str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.table_read_group + '">')
            str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Summary</td>\n')
            str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

            return

        def report_hub():
            """Private function to create a UCSC Track Hub.

            @return:
            @rtype:
            """

            str_list = list()
            """ @type str_list: list[str | unicode] """

            # Group via UCSC super tracks.

            str_list.append('track Alignments\n')
            str_list.append('shortLabel Alignments\n')
            str_list.append('longLabel STAR alignments\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignments\n')
            str_list.append('\n')

            # Sample-specific tracks

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                runnable_merge = self.runnable_dict[
                    '_'.join((self.stage_name_merge, sample.name))]
                file_path_merge = runnable_merge.file_path_object
                """ @type file_path_merge: FilePathStarMerge """

                #
                # Add a trackDB entry for each Tophat accepted_hits.bam file.
                #
                # Common settings
                str_list.append('track ' + sample.name + '_alignments\n')
                str_list.append('type bam\n')
                str_list.append('shortLabel ' + sample.name + '_alignments\n')
                str_list.append('longLabel ' + sample.name + ' STAR alignments\n')
                str_list.append('bigDataUrl ' + file_path_merge.merged_bam + '\n')
                # str_list.append('html ...\n')
                str_list.append('visibility dense\n')

                # Common optional settings
                str_list.append('color 0,0,0\n')

                # bam/cram - Compressed Sequence Alignment track settings
                # None

                # Composite track settings
                str_list.append('parent Alignments\n')
                str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return
