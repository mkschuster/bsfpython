"""bsf.analyses.star_aligner

A package of classes and methods supporting the spliced Transcripts Alignment
to a Reference (STAR) aligner by Alexander Dobin.

Project:  https://github.com/alexdobin/STAR
"""

#
# Copyright 2013 - 2016 Michael K. Schuster
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


import os

import pysam

from bsf import Analysis, FilePath, Runnable
from bsf.process import RunnableStep, RunnableStepLink, RunnableStepMove, RunnableStepPicard
from bsf.standards import Default


class FilePathStarAlign(FilePath):

    def __init__(self, prefix):

        super(FilePathStarAlign, self).__init__(prefix=prefix)

        self.aligned_sam = prefix + '_Aligned.out.sam'

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
            comparisons=None,
            sample_list=None,
            replicate_grouping=False,
            index_directory=None,
            stranded=None,
            transcriptome_gtf=None,
            classpath_picard=None):
        """Initialise a C{bsf.analyses.rna_seq.Tuxedo} object.

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
        @param comparisons: Python C{dict} of Python C{str} (comparison name) key objects and
            Python C{tuple} value objects of C{bsf.ngs.Sample.name} and Python C{list} of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, (bsf.ngs.Sample.name, list[bsf.ngs.Sample])]
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
            comparisons=comparisons,
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
        """Set instance variables of a C{bsf.analyses.rna_seq.Tuxedo} object via a section of a
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

    def _read_comparisons(self):

        self.sample_list.extend(self.collection.get_all_samples())

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.Tuxedo} analysis.

        Although the STAR aligner can directly count reads according to its splice junction database,
        more than one read group may need aligning so that the count tables had to be combined.
        @return:
        @rtype:
        """

        super(StarAligner, self).run()

        # Get global defaults.

        default = Default.get_global_default()

        # The STAR Aligner requires a genome version.

        if not self.genome_version:
            raise Exception('A STAR Aligner analysis requires a genome_version configuration option.')

        if not self.index_directory:
            raise Exception('A STAR Aligner analysis requires an index_directory configuration option.')

        if self.stranded:
            if self.stranded not in ('yes', 'no', 'reverse'):
                raise Exception(
                    'The STAR Aligner configuration option "stranded" can only be '
                    '"yes", "no" or "reverse", not {}.'.format(self.stranded))
        else:
            self.stranded = 'yes'

        if not self.transcriptome_gtf:
            raise Exception('A STAR Aligner analysis requires a transcriptome_gtf configuration option.')

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        self._read_comparisons()

        stage_align = self.get_stage(name=self.stage_name_align)
        stage_index = self.get_stage(name=self.stage_name_index)
        stage_merge = self.get_stage(name=self.stage_name_merge)

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            runnable_index_list = list()
            """ @type runnable_index_list: list[bsf.Runnable] """

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

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
                    if paired_reads.reads_2 is None:
                        runnable_step.add_option_long(
                            key='readFilesIn',
                            value=paired_reads.reads_1.file_path)
                    else:
                        runnable_step.add_option_long(
                            key='readFilesIn',
                            value=' '.join((paired_reads.reads_1.file_path, paired_reads.reads_2.file_path)))

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
                                read_group_str = "ID:{}".format(read_group_dict['ID'])
                                for key in key_list:
                                    read_group_str += ' "{}:{}"'.format(key, read_group_dict[key])
                                runnable_step.add_option_long(key='outSAMattrRGline', value=read_group_str)

                            # Add @PG lines.
                            # The STAR aligner allows only a single @PG line, which is not terribly helpful.
                            # for program_dict in alignment_file.header['PG']:
                            #     """ @type program_dict: dict[str, str] """
                            #     key_list = program_dict.keys()
                            #     program_str = ''
                            #     for key in key_list:
                            #         program_str += ' "{}:{}"'.format(key, program_dict[key])
                            #     runnable_step.add_option_long(key='outSAMheaderPG', value=program_str)

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
                        RunnableStepPicard(
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
                        RunnableStepPicard(
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

            # Add dependencies on Runnable objects of the indexing stage.
            for runnable_index in runnable_index_list:
                executable_merge.dependencies.append(runnable_index.name)

            if len(runnable_index_list) == 1:
                runnable_index = runnable_index_list[0]
                file_path_index = runnable_index.file_path_object
                """ @type file_path_index: FilePathStarIndex """
                # For a single ReadPair, just rename the files.
                runnable_merge.add_runnable_step(
                    RunnableStepMove(
                        name='move_bam',
                        source_path=file_path_index.aligned_bam,
                        target_path=file_path_merge.merged_bam))
                runnable_merge.add_runnable_step(
                    RunnableStepMove(
                        name='move_bai',
                        source_path=file_path_index.aligned_bai,
                        target_path=file_path_merge.merged_bai))
                runnable_merge.add_runnable_step(
                    RunnableStepMove(
                        name='move_md5',
                        source_path=file_path_index.aligned_md5,
                        target_path=file_path_merge.merged_md5))
            else:
                # Run Picard MergeSamFiles on each BAM file.
                runnable_step = runnable_merge.add_runnable_step(
                    RunnableStepPicard(
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

            if False:
                # Create a RunnableStep for htseq-count.
                # For the moment, the counting infrastructure in Bioconductor GenomeAlignment is used
                # that directly yields a RangedSummarizedExperiment object suitable for downstream analysis.
                runnable_step = runnable_merge.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='count',
                        program='htseq-count',
                        stdout_path=file_path_merge.merged_tsv))
                runnable_step.add_option_long(key='format', value='bam')
                runnable_step.add_option_long(key='order', value='pos')
                runnable_step.add_option_long(key='stranded', value=self.stranded)
                runnable_step.add_switch_long(key='quiet')

                runnable_step.arguments.append(file_path_merge.merged_bam)
                runnable_step.arguments.append(self.transcriptome_gtf)

            if False:
                # The htseq-qa script no longer seems compatible with the Python matplotlib 1.5.3 module
                runnable_step = runnable_merge.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='qa',
                        program='htseq-qa'))
                runnable_step.add_option_long(key='type', value='bam')
                runnable_step.add_option_long(key='outfile', value=file_path_merge.merged_pdf)
                # readlength can be guessed from the file
                # gamma defaults to 0.3
                # nosplit defaults ot false
                # NOTE: The 'maxqual' option defaults to 41, but the HiSeq 3000/4000 platform may set 42 as maximum.
                runnable_step.add_option_long(key='maxqual', value='42')

                runnable_step.arguments.append(file_path_merge.merged_bam)

        return
