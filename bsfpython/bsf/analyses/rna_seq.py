"""bsf.analyses.rna_seq

A package of classes and methods supporting RNA-Seq analyses.
"""

#
# Copyright 2013 - 2017 Michael K. Schuster
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


import errno
import os
from pickle import Pickler, HIGHEST_PROTOCOL
import re

from bsf import Analysis, FilePath, Runnable
from bsf.analyses.star_aligner import StarAligner
from bsf.annotation import AnnotationSheet, TuxedoSamplePairSheet
from bsf.executables import TopHat
from bsf.process import Executable, RunnableStep, RunnableStepLink
from bsf.standards import Default


class FilePathTophat(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathTophat} models files in a sample-specific TopHat directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar accepted_hits_bam: TopHat accepted hits BAM file
    @type accepted_hits_bam: str | unicode
    @ivar accepted_hits_bam_link_source: TopHat accepted hits BAM file symbolic link source
    @type accepted_hits_bam_link_source: str | unicode
    @ivar accepted_hits_bam_link_target: TopHat accepted hits BAM file symbolic link target
    @type accepted_hits_bam_link_target: str | unicode
    @ivar accepted_hits_bai: TopHat accepted hits BAI file
    @type accepted_hits_bai: str | unicode
    @ivar accepted_hits_bai_link_source: TopHat accepted hits BAI file symbolic link source
    @type accepted_hits_bai_link_source: str | unicode
    @ivar accepted_hits_bai_link_target: TopHat accepted hits BAI file symbolic link target
    @type accepted_hits_bai_link_target: str | unicode
    @ivar accepted_hits_bw: TopHat accepted hits bigWig file
    @type accepted_hits_bw: str | unicode
    @ivar align_summary: TopHat align summary file
    @type align_summary: str | unicode
    @ivar deletions_bb: TopHat deletions bigBed file
    @type deletions_bb: str | unicode
    @ivar deletions_bed: TopHat deletions BED file
    @type deletions_bed: str | unicode
    @ivar insertions_bb: TopHat insertions bigBed file
    @type insertions_bb: str | unicode
    @ivar insertions_bed: TopHat insertions BED file
    @type insertions_bed: str | unicode
    @ivar junctions_bb: TopHat junctions bigBed file
    @type junctions_bb: str | unicode
    @ivar junctions_bed: TopHat junctions BED file
    @type junctions_bed: str | unicode
    @ivar prep_reads_info: TopHat prepare reads information file
    @type prep_reads_info: str | unicode
    @ivar unmapped_bam: TopHat unmapped BAM file
    @type unmapped_bam: str | unicode
    @ivar unmapped_bam_link_source: TopHat unmapped BAM file symbolic link source
    @type unmapped_bam_link_source: str | unicode
    @ivar unmapped_bam_link_target: TopHat unmapped BAM file symbolic link target
    @type unmapped_bam_link_target: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathTophat} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathTophat, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.accepted_hits_bam = os.path.join(prefix, 'accepted_hits.bam')
        self.accepted_hits_bam_link_source = 'accepted_hits.bam'
        self.accepted_hits_bam_link_target = os.path.join(prefix, prefix + '_accepted_hits.bam')
        self.accepted_hits_bai = os.path.join(prefix, 'accepted_hits.bam.bai')
        self.accepted_hits_bai_link_source = 'accepted_hits.bam.bai'
        self.accepted_hits_bai_link_target = os.path.join(prefix, prefix + '_accepted_hits.bam.bai')
        self.accepted_hits_bw = os.path.join(prefix, 'accepted_hits.bw')
        self.align_summary = os.path.join(prefix, 'align_summary.txt')
        self.deletions_bb = os.path.join(prefix, 'deletions.bb')
        self.deletions_bed = os.path.join(prefix, 'deletions.bed')
        self.insertions_bb = os.path.join(prefix, 'insertions.bb')
        self.insertions_bed = os.path.join(prefix, 'insertions.bed')
        self.junctions_bb = os.path.join(prefix, 'junctions.bb')
        self.junctions_bed = os.path.join(prefix, 'junctions.bed')
        self.prep_reads_info = os.path.join(prefix, 'prep_reads.info')
        self.unmapped_bam = os.path.join(prefix, 'unmapped.bam')
        self.unmapped_bam_link_source = 'unmapped.bam'
        self.unmapped_bam_link_target = os.path.join(prefix, prefix + '_unmapped.bam')

        return


class FilePathCufflinks(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCufflinks} models files in a sample-specific Cufflinks directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar fpkm_tracking_genes_tsv: Cufflinks FPKM tracking genes tab-separated value (TSV) file
    @type fpkm_tracking_genes_tsv: str | unicode
    @ivar fpkm_tracking_isoforms_tsv: Cufflinks FPKM tracking isoforms tab-separated value (TSV) file
    @type fpkm_tracking_isoforms_tsv: str | unicode
    @ivar skipped_gtf: Cufflinks skipped regions GTF file
    @type skipped_gtf: str | unicode
    @ivar skipped_gtf_link_source: Cufflinks skipped regions GTF symbolic links source
    @type skipped_gtf_link_source: str | unicode
    @ivar skipped_gtf_link_target: Cufflinks skipped regions GTF symbolic links target
    @type skipped_gtf_link_target: str | unicode
    @ivar temporary_big_gene_prediction: Temporary UCSC big gene prediction (bigGenePred) file
    @type temporary_big_gene_prediction: str | unicode
    @ivar temporary_gene_prediction: Temporary UCSC gene prediction (genePred) file
    @type temporary_gene_prediction: str | unicode
    @ivar temporary_sorted_tsv: Temporary sorted tab-separated value (TSV) file
    @type temporary_sorted_tsv: str | unicode
    @ivar transcripts_bb: Cufflinks transcript assembly bigBed file
    @type transcripts_bb: str | unicode
    @ivar transcripts_bb_link_source: Cufflinks transcript assembly bigBed symbolic link source
    @type transcripts_bb_link_source: str | unicode
    @ivar transcripts_bb_link_target: Cufflinks transcript assembly bigBed symbolic link target
    @type transcripts_bb_link_target: str | unicode
    @ivar transcripts_gtf: Cufflinks transcript assembly GTF file
    @type transcripts_gtf: str | unicode
    @ivar transcripts_gtf_link_source: Cufflinks transcript assembly GTF symbolic links source
    @type transcripts_gtf_link_source: str | unicode
    @ivar transcripts_gtf_link_target: Cufflinks transcript assembly GTF symbolic links target
    @type transcripts_gtf_link_target: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCufflinks} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCufflinks, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.fpkm_tracking_genes_tsv = os.path.join(prefix, prefix + '_genes_fpkm_tracking.tsv')
        self.fpkm_tracking_isoforms_tsv = os.path.join(prefix, prefix + '_isoforms_fpkm_tracking.tsv')
        self.skipped_gtf = os.path.join(prefix, 'skipped.gtf')
        self.skipped_gtf_link_source = 'skipped.gtf'
        self.skipped_gtf_link_target = os.path.join(prefix, prefix + '_skipped.gtf')
        self.temporary_big_gene_prediction = os.path.join(prefix, 'transcripts_big_gene_prediction.tsv')
        self.temporary_gene_prediction = os.path.join(prefix, 'transcripts_gene_prediction.tsv')
        self.temporary_sorted_tsv = os.path.join(prefix, 'transcripts_sorted.tsv')
        self.transcripts_bb = os.path.join(prefix, 'transcripts.bb')
        self.transcripts_bb_link_source = 'transcripts.bb'
        self.transcripts_bb_link_target = os.path.join(prefix, prefix + '_transcripts.bb')
        self.transcripts_gtf = os.path.join(prefix, 'transcripts.gtf')
        self.transcripts_gtf_link_source = 'transcripts.gtf'
        self.transcripts_gtf_link_target = os.path.join(prefix, prefix + '_transcripts.gtf')

        return


class FilePathCuffmerge(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffmerge} models files in a comparison-specific Cuffmerge directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar assembly_txt: Assembly text file
    @type assembly_txt: str | unicode
    @ivar merged_bb: Cuffmerge transcript assembly bigBed file
    @type merged_bb: str | unicode
    @ivar merged_bb_link_source: Cuffmerge transcript assembly bigBed symbolic link source
    @type merged_bb_link_source: str | unicode
    @ivar merged_bb_link_target: Cuffmerge transcript assembly bigBed symbolic link target
    @type merged_bb_link_target: str | unicode
    @ivar merged_gtf: Cuffmerge merged GTF file
    @type merged_gtf: str | unicode
    @ivar merged_gtf_link_source: Cuffmerge transcript assembly GTF symbolic links source
    @type merged_gtf_link_source: str | unicode
    @ivar merged_gtf_link_target: Cuffmerge transcript assembly GTF symbolic links target
    @type merged_gtf_link_target: str | unicode
    @ivar temporary_gene_prediction: Temporary UCSC gene prediction (genePred) file
    @type temporary_gene_prediction: str | unicode
    @ivar temporary_big_gene_prediction: Temporary UCSC big gene prediction (bigGenePred) file
    @type temporary_big_gene_prediction: str | unicode
    @ivar temporary_sorted_tsv: Temporary sorted tab-separated value (TSV) file
    @type temporary_sorted_tsv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCufflinks} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffmerge, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.assembly_txt = '_'.join((prefix, 'assembly.txt'))
        self.merged_bb = os.path.join(prefix, 'merged.bb')
        self.merged_bb_link_source = 'merged.bb'
        self.merged_bb_link_target = os.path.join(prefix, prefix + '_merged.bb')
        self.merged_gtf = os.path.join(prefix, 'merged.gtf')
        self.merged_gtf_link_source = 'merged.gtf'
        self.merged_gtf_link_target = os.path.join(prefix, prefix + '_merged.gtf')
        self.temporary_gene_prediction = os.path.join(prefix, 'merged_gene_prediction.tsv')
        self.temporary_big_gene_prediction = os.path.join(prefix, 'merged_big_gene_prediction.tsv')
        self.temporary_sorted_tsv = os.path.join(prefix, 'merged_sorted.tsv')

        return


class FilePathCuffquant(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffquant} models files in a sample-specific Cuffquant directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar abundances: Cuffquant abundances file
    @type abundances: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCuffquant} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffquant, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances = os.path.join(prefix, 'abundances.cxb')
        # TODO: Remove after rearranging _create_tophat_...() and _create_cuffmerge_...() into run().
        # self.tophat_accepted_hits = os.path.join('_'.join(('rnaseq_tophat', paired_reads_name)), 'accepted_hits.bam')
        self.tophat_accepted_hits = ''

        return


class FilePathCuffnorm(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffnorm} models files in a comparison-specific Cuffnorm directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCuffnorm} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffnorm, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class FilePathCuffdiff(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffdiff} models files in a comparison-specific Cuffdiff directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCuffdiff} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffdiff, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class Tuxedo(Analysis):
    """Tuxedo RNASeq C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_run_tophat: C{bsf.Stage.name} for the run Tophat stage
    @type stage_name_run_tophat: str
    @cvar stage_name_process_tophat: C{bsf.Stage.name} for the process Tophat stage
    @type stage_name_process_tophat: str
    @cvar stage_name_run_cufflinks: C{bsf.Stage.name} for the run Cufflinks stage
    @type stage_name_run_cufflinks: str
    @cvar stage_name_process_cufflinks: C{bsf.Stage.name} for the process Cufflinks stage
    @type stage_name_process_cufflinks: str
    @cvar stage_name_run_cuffmerge: C{bsf.Stage.name} for the Cuffmerge stage
    @type stage_name_run_cuffmerge: str
    @cvar stage_name_run_cuffquant: C{bsf.Stage.name} for the Cuffquant stage
    @type stage_name_run_cuffquant: str
    @cvar stage_name_run_cuffnorm: C{bsf.Stage.name} for the Cuffnorm stage
    @type stage_name_run_cuffnorm: str
    @cvar stage_name_run_cuffdiff: C{bsf.Stage.name} for the Cuffdiff stage
    @type stage_name_run_cuffdiff: str
    @cvar stage_name_process_cuffdiff: C{bsf.Stage.name} for the process Cuffdiff stage
    @type stage_name_process_cuffdiff: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar cmp_file: Comparison file
    @type cmp_file: str | unicode
    @ivar genome_index_path: Bowtie genome index path
    @type genome_index_path: str | unicode
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | unicode
    @ivar transcriptome_gtf_path: Reference transcriptome GTF file path
    @type transcriptome_gtf_path: str | unicode
    @ivar mask_gtf_path: GTF file path to mask transcripts
    @type mask_gtf_path: str | unicode
    @ivar multi_read_correction: Apply multi-read correction
    @type multi_read_correction: bool
    @ivar library_type: Library type
        Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
    @type library_type: str
    @ivar novel_transcripts: Assemble novel transcripts
    @type novel_transcripts: bool
    @ivar no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
    @type no_length_correction: bool
    """

    name = 'RNA-seq Analysis'
    prefix = 'rnaseq'

    # Replicate stage
    stage_name_run_tophat = '_'.join((prefix, 'run_tophat'))
    stage_name_process_tophat = '_'.join((prefix, 'process_tophat'))
    stage_name_run_cufflinks = '_'.join((prefix, 'run_cufflinks'))
    stage_name_process_cufflinks = '_'.join((prefix, 'process_cufflinks'))

    # Comparison stage
    stage_name_run_cuffmerge = '_'.join((prefix, 'cuffmerge'))
    stage_name_run_cuffquant = '_'.join((prefix, 'cuffquant'))
    stage_name_run_cuffnorm = '_'.join((prefix, 'cuffnorm'))
    stage_name_run_cuffdiff = '_'.join((prefix, 'cuffdiff'))
    stage_name_process_cuffdiff = '_'.join((prefix, 'process_cuffdiff'))

    @classmethod
    def get_prefix_rnaseq_run_tophat(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: Replicate key
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_tophat, paired_reads_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffmerge(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffmerge, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffquant(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffquant, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffnorm(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffnorm, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffdiff(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffdiff, comparison_name))

    @classmethod
    def get_prefix_rnaseq_process_cuffdiff(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_process_cuffdiff, comparison_name))

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
            cmp_file=None,
            genome_index_path=None,
            genome_fasta_path=None,
            transcriptome_version=None,
            transcriptome_index_path=None,
            transcriptome_gtf_path=None,
            mask_gtf_path=None,
            multi_read_correction=None,
            library_type=None,
            novel_transcripts=None,
            no_length_correction=False):
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
        @type comparisons: dict[str, (str, list[bsf.ngs.Sample])]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @param genome_index_path: Bowtie genome index path
        @type genome_index_path: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode
        @param mask_gtf_path: GTF file path to mask transcripts
        @type mask_gtf_path: str | unicode
        @param multi_read_correction: Apply multi-read correction
        @type multi_read_correction: bool
        @param library_type: Library type
            Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
        @type library_type: str
        @param novel_transcripts: Assemble novel transcripts
        @type novel_transcripts: bool
        @param no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
        @type no_length_correction: bool
        @return:
        @rtype:
        """

        super(Tuxedo, self).__init__(
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

        if cmp_file is None:
            self.cmp_file = str()
        else:
            self.cmp_file = cmp_file

        if genome_index_path is None:
            self.genome_index_path = str()
        else:
            self.genome_index_path = genome_index_path

        if genome_fasta_path is None:
            self.genome_fasta_path = str()
        else:
            self.genome_fasta_path = genome_fasta_path

        if transcriptome_version is None:
            self.transcriptome_version = str()
        else:
            self.transcriptome_version = transcriptome_version

        if transcriptome_index_path is None:
            self.transcriptome_index_path = str()
        else:
            self.transcriptome_index_path = transcriptome_index_path

        if transcriptome_gtf_path is None:
            self.transcriptome_gtf_path = str()
        else:
            self.transcriptome_gtf_path = transcriptome_gtf_path

        if mask_gtf_path is None:
            self.mask_gtf_path = str()
        else:
            self.mask_gtf_path = mask_gtf_path

        if multi_read_correction is None:
            self.multi_read_correction = False
        else:
            assert isinstance(multi_read_correction, bool)
            self.multi_read_correction = multi_read_correction

        if library_type is None:
            self.library_type = str()
        else:
            self.library_type = library_type

        if novel_transcripts is None:
            self.novel_transcripts = True
        else:
            self.novel_transcripts = novel_transcripts

        if no_length_correction is None:
            self.no_length_correction = False
        else:
            assert isinstance(no_length_correction, bool)
            self.no_length_correction = no_length_correction

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

        super(Tuxedo, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cmp_file = configuration.config_parser.get(section=section, option=option)

        option = 'genome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'mask_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mask_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'multi_read_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.multi_read_correction = configuration.config_parser.getboolean(section=section, option=option)

        option = 'library_type'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_type = configuration.config_parser.get(section=section, option=option)

        option = 'novel_transcripts'
        if configuration.config_parser.has_option(section=section, option=option):
            self.novel_transcripts = configuration.config_parser.getboolean(section=section, option=option)

        option = 'no_length_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.no_length_correction = configuration.config_parser.getboolean(section=section, option=option)

        return

    def _read_comparisons(self, cmp_file):
        """Read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

        All C{bsf.ngs.Sample} objects referenced in a comparison are added from the C{bsf.ngs.Collection} to the
        C{bsf.Analysis} object.

            - Column headers for CASAVA folders:
                - Treatment/Control/Point N ProcessedRunFolder:
                    - CASAVA processed run folder name or
                    - C{bsf.Analysis.input_directory} by default
                - Treatment/Control/Point N Project:
                    - CASAVA Project name or
                    - C{bsf.Analysis.project_name} by default
                - Treatment/Control/Point N Sample:
                    - CASAVA Sample name, no default
            - Column headers for independent samples:
                - Treatment/Control/Point N Sample:
                - Treatment/Control/Point N Reads:
                - Treatment/Control/Point N File:
        @param cmp_file: Comparisons file path
        @type cmp_file: str | unicode
        @return:
        @rtype:
        """

        if self.debug > 1:
            print '{!r} method _read_comparisons:'.format(self)

        regular_expression = re.compile(pattern='\W')

        annotation_sheet = AnnotationSheet.from_file_path(file_path=cmp_file)

        # TODO: Adjust by introducing a new class RNASeqComparisonSheet(AnnotationSheet) in this module?
        for row_dict in annotation_sheet.row_dicts:
            # In addition to defining samples, allow also the definition of groups in comparison files.
            # If the row dictionary has a 'Group' key, then the Sample in the same row gets added to the group.
            # So, 'ProcessedRunFolder', 'Project', 'Sample', 'Group' defines the groups, while ...
            # 'Control Group','Treatment Group' defines a comparison, as does ...
            # 'Control Group','Treatment ProcessedRunFolder','Treatment Project','Treatment Sample'

            # Get Sample objects for classical 'Control' and 'Treatment' keys,
            # before looking up a series of 'Point N' keys.
            i = -2
            key = str()
            comparison_groups = list()
            """ @type comparison_groups: list[(str, list[bsf.ngs.Sample])] """
            while True:
                i += 1
                if i == -1:
                    prefix = 'Control'
                elif i == 0:
                    prefix = 'Treatment'
                else:
                    prefix = 'Point {}'.format(i)
                # Get Sample objects for 'Point N' keys for as long as they are defined.
                # The Collection.get_sample_from_row_dict method can return one or more Sample objects,
                # depending on 'Group' or 'Sample' column entries.
                # In RNA-Seq experiments, entire pools of Sample objects (replicates) are compared with each other.
                group_name, group_samples = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                if group_name and len(group_samples):
                    key += group_name
                    key += '__'
                    comparison_groups.append((group_name, group_samples))
                    # Also expand each Python list of bsf.ngs.Sample objects to get all those bsf.ngs.Sample objects
                    # that this bsf.Analysis needs considering.
                    for sample in group_samples:
                        if self.debug > 1:
                            print '  {} Sample name: {!r} file_path: {!r}'.format(prefix, sample.name, sample.file_path)
                            # print sample.trace(1)
                        self.add_sample(sample=sample)
                elif i < 1:
                    # A Control and Treatment prefix is not required.
                    continue
                else:
                    # Only break if there is no further 'Point N' prefix.
                    break

            # For a successful comparison, more than one Sample (pool) has to be defined.

            if len(comparison_groups) < 2:
                if self.debug > 1:
                    print 'Comparison line with less than two Sample or Group keys. {!r}'. \
                        format(row_dict)
                continue

            if 'Comparison Name' in row_dict and row_dict['Comparison Name']:
                # For ridiculously large comparisons involving loads of groups or samples a comparison name can be
                # explicitly specified. Any non-word characters get replaced by underscore characters.
                key = re.sub(pattern=regular_expression, repl='_', string=row_dict['Comparison Name'])
            else:
                # Truncate the last '__' separator off the key string.
                key = key.rstrip('_')

            self.comparisons[key] = comparison_groups
            """ @type comparisons: dic[str, (str, list[bsf.ngs.Sample])] """

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.Tuxedo} analysis.
        @return:
        @rtype:
        """

        super(Tuxedo, self).run()

        # Get global defaults.

        default = Default.get_global_default()

        # Method configuration with regards to Cuffquant and Cuffdiff.
        run_cuffquant_before_cuffdiff = False

        # Tuxedo requires a genome version.

        if not self.genome_version:
            raise Exception('A Tuxedo analysis requires a genome_version configuration option.')

        # Tuxedo requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A Tuxedo analysis requires a transcriptome_version configuration option.')

        if self.cmp_file:
            # A comparison file path has been provided.
            if self.cmp_file == '*groups*':
                # The special file name *groups* creates the comparisons on the basis of an
                # all-against-all group comparison.
                # Without a comparison file path, simply add all Sample objects from the Collection.
                self.sample_list.extend(self.collection.get_all_samples())
                # Create a global comparison by adding all sample groups.
                comparison_groups = list()
                """ @type comparison_groups: list[(str, list[bsf.ngs.Sample])] """
                for key in self.collection.sample_group_dict.keys():
                    comparison_groups.append((key, self.collection.sample_group_dict[key]))
                # Sort the list of comparison groups by group name.
                comparison_groups.sort(cmp=lambda x, y: cmp(x[0], y[0]))
                self.comparisons['global'] = comparison_groups
            elif self.cmp_file == '*samples*':
                # The special file name *samples* creates the comparisons on the basis of an
                # all-against-all sample comparison.
                # Without a comparison file path, simply add all Sample objects from the Collection.
                self.sample_list.extend(self.collection.get_all_samples())
                # Create a global comparison by adding all samples under their sample name as group name.
                comparison_groups = list()
                """ @type comparison_groups: list[(str, list[bsf.ngs.Sample])] """
                for sample in self.sample_list:
                    # Add a tuple of group (i.e. sample) name and a Python list of the Sample object.
                    comparison_groups.append((sample.name, [sample]))
                # Sort the list of comparison groups by Sample.name.
                comparison_groups.sort(cmp=lambda x, y: cmp(x[0], y[0]))
                self.comparisons['global'] = comparison_groups
            else:
                # A comparison file path has been provided.
                # Expand an eventual user part i.e. on UNIX ~ or ~user and
                # expand any environment variables i.e. on UNIX ${NAME} or $NAME
                # Check if an absolute path has been provided, if not,
                # automatically prepend standard directory paths.
                self.cmp_file = os.path.expanduser(path=self.cmp_file)
                self.cmp_file = os.path.expandvars(path=self.cmp_file)
                if not os.path.isabs(self.cmp_file) and not os.path.exists(path=self.cmp_file):
                    self.cmp_file = os.path.join(self.project_directory, self.cmp_file)
                # Read and process the comparison file, which includes adding only those sample objects,
                # which are referenced in a comparison.
                self._read_comparisons(cmp_file=self.cmp_file)
        else:
            # Without a comparison file path, simply add all Sample objects from the Collection.
            # This means that only the initial pipeline stages, but not the comparison stage gets run.
            self.sample_list.extend(self.collection.get_all_samples())

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Sample objects are only put into the
        # bsf.Analysis.samples list by the _read_comparisons method.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        # Get the genome index

        if not self.genome_index_path:
            self.genome_index_path = os.path.join(
                Default.absolute_genome_resource(genome_version=self.genome_version),
                default.indices['bowtie2'],
                self.genome_version)

        if not self.genome_fasta_path:
            self.genome_fasta_path = self.genome_index_path + '.fa'

        # FIXME: For the moment use the FAI file instead of the UCSC chromosome sizes file.
        genome_sizes_path = self.genome_fasta_path + '.fai'

        if not os.path.exists(self.genome_fasta_path):
            raise Exception("Genome FASTA file path {!r} does not exists.".format(self.genome_fasta_path))

        # Define a reference transcriptome index directory or a GTF file path.

        if self.transcriptome_index_path:
            # Check if the transcriptome_index_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_index_path = Default.get_absolute_path(
                file_path=self.transcriptome_index_path,
                default_path=Default.absolute_transcriptome_resource(transcriptome_version=''))

            if not os.path.isdir(self.transcriptome_index_path):
                raise Exception("Reference transcriptome index directory {!r} does not exist.".
                                format(self.transcriptome_index_path))

            transcriptome_index = os.path.basename(self.transcriptome_index_path)

            # Does an indices_for_TopHat directory exist?
            transcriptome_index_path = os.path.join(self.transcriptome_index_path, default.indices['tophat'])
            if os.path.isdir(transcriptome_index_path):
                self.transcriptome_index_path = transcriptome_index_path

            # Finally, set the transcriptome GTF file path.
            # The tophat --transcript-index process puts a GFF file into the index directory
            # that really is a GTF file. A symbolic link to a GTF file is needed to make the
            # process cuffdiff script work.
            # For the moment, use the symbolic link in the indices_for_TopHat directory.

            self.transcriptome_gtf_path = os.path.join(
                self.transcriptome_index_path,
                '.'.join((transcriptome_index, 'gtf')))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        elif self.transcriptome_gtf_path:
            # Check if transcriptome_gtf_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_gtf_path = Default.get_absolute_path(
                file_path=self.transcriptome_gtf_path,
                default_path=Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        else:
            # Neither was provided, automatically discover on the basis of the transcriptome version.
            self.transcriptome_index_path = os.path.join(
                Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version),
                default.indices['tophat'],
                self.transcriptome_version,  # TopHat puts the transcriptome index into a sub directory.
                self.transcriptome_version)

            self.transcriptome_gtf_path = os.path.join(
                Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version),
                default.indices['tophat'],
                self.transcriptome_version) + '.gtf'
            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file path {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        if not self.transcriptome_gtf_path:
            raise Exception("Reference transcriptome GTF file not defined.")

        # Read configuration options.

        # TODO: Move the ConfigParser code.
        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        stage_run_tophat = self.get_stage(name=self.stage_name_run_tophat)
        stage_process_tophat = self.get_stage(name=self.stage_name_process_tophat)
        stage_run_cufflinks = self.get_stage(name=self.stage_name_run_cufflinks)
        stage_process_cufflinks = self.get_stage(name=self.stage_name_process_cufflinks)
        stage_run_cuffmerge = self.get_stage(name=self.stage_name_run_cuffmerge)
        stage_run_cuffquant = self.get_stage(name=self.stage_name_run_cuffquant)
        stage_run_cuffnorm = self.get_stage(name=self.stage_name_run_cuffnorm)
        stage_run_cuffdiff = self.get_stage(name=self.stage_name_run_cuffdiff)
        stage_process_cuffdiff = self.get_stage(name=self.stage_name_process_cuffdiff)

        runnable_run_cufflinks_list = list()
        """ @type runnable_run_cufflinks_list: list[Runnable] """

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                # Create a Tophat Runnable per paired_reads_name.

                # TODO: Activate the new code once the bsf_run_rnaseq_tophat.py script has been retired.

                # prefix_run_tophat = '_'.join((stage_run_tophat.name, paired_reads_name))

                # The output directory deviates from the prefix_run_tophat that itself is based on
                # stage_run_tophat.name. Both rnaseq_run_tophat and rnaseq_process_tophat processes should
                # use the same rnaseq_tophat directory.
                prefix_tophat = '_'.join(('rnaseq_tophat', paired_reads_name))

                file_path_tophat = FilePathTophat(prefix=prefix_tophat)

                # runnable_run_tophat = self.add_runnable(
                #         runnable=Runnable(
                #                 name=self.get_prefix_rnaseq_run_tophat(paired_reads_name=paired_reads_name),
                #                 code_module='bsf.runnables.generic',
                #                 working_directory=self.genome_directory,
                #                 file_path_object=file_path_tophat,
                #                 debug=self.debug))
                # executable_run_tophat = self.set_stage_runnable(
                #         stage=stage_run_tophat,
                #         runnable=runnable_run_tophat)
                #
                # Create a new Tophat RunnableStep.

                # tophat = runnable_run_tophat.add_runnable_step(
                #     runnable_step=RunnableStep(
                #     name='tophat2',
                #     program='tophat2'))

                tophat = TopHat(
                    name='_'.join(('rnaseq_tophat', paired_reads_name)),
                    analysis=self)

                # Set tophat options.

                tophat.add_option_long(
                    key='GTF',
                    value=self.transcriptome_gtf_path)
                if self.transcriptome_index_path:
                    tophat.add_option_long(
                        key='transcriptome-index',
                        value=self.transcriptome_index_path)
                tophat.add_option_long(
                    key='output-dir',
                    value=file_path_tophat.output_directory)
                tophat.add_option_long(
                    key='num-threads',
                    value=str(stage_run_tophat.threads))
                # TODO: These really are properties of the Reads, PairedReads or Sample objects rather than an Analysis.
                # TODO: Move the ConfigParser code.
                if config_parser.has_option(section=config_section, option='insert_size'):
                    insert_size = config_parser.getint(section=config_section, option='insert_size')
                    read_length = config_parser.getint(section=config_section, option='read_length')
                    mate_inner_dist = insert_size - 2 * read_length
                    tophat.add_option_long(
                        key='mate-inner-dist',
                        value=str(mate_inner_dist))
                # TODO: Move the ConfigParser code.
                if config_parser.has_option(section=config_section, option='mate-std-dev'):
                    tophat.add_option_long(
                        key='mate-std-dev',
                        value=config_parser.getint(section=config_section, option='mate-std-dev'))
                if self.library_type:
                    tophat.add_option_long(
                        key='library-type',
                        value=self.library_type)
                # The TopHat coverage search finds additional "GT-AG" introns, but is only recommended for
                # short reads (< 45 bp) and small read numbers (<= 10 M).
                # TODO: This option should possibly become configurable per sample.
                tophat.add_switch_long(key='no-coverage-search')
                # TODO: Set -rg-* options to back fill the read group from Illumina2bam.

                # Set rnaseq_tophat arguments.

                tophat.arguments.append(self.genome_index_path)

                # Set rnaseq_tophat arguments for reads1 and reads2.

                reads_1_file_path_list = list()
                """ @type reads_1_file_path_list: list[str | unicode] """
                reads_2_file_path_list = list()
                """ @type reads_2_file_path_list: list[str | unicode] """

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads_1:
                        reads_1_file_path_list.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads_2_file_path_list.append(paired_reads.reads_2.file_path)

                tophat.arguments.append(','.join(reads_1_file_path_list))
                tophat.arguments.append(','.join(reads_2_file_path_list))

                # Create a new rnaseq_run_tophat Executable.
                # TODO: The following code block is required as long as the bsf_run_rnaseq_tophat.py script
                # has not been retired.

                pickler_dict_run_tophat = {
                    'prefix': stage_run_tophat.name,
                    'replicate_key': paired_reads_name,
                    'tophat_executable': tophat,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(stage_run_tophat.name, paired_reads_name))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_run_tophat)
                pickler_file.close()

                executable_run_tophat = stage_run_tophat.add_executable(
                    executable=Executable(
                        name='_'.join((stage_run_tophat.name, paired_reads_name)),
                        program='bsf_run_rnaseq_tophat.py'))
                # Set dependencies on previous Runnable or Executable objects.
                # None.

                # Set rnaseq_run_tophat options.
                self.set_command_configuration(command=executable_run_tophat)
                executable_run_tophat.add_option_long(key='pickler_path', value=pickler_path)
                executable_run_tophat.add_option_long(key='debug', value=str(self.debug))

                # Only submit this Executable if the "align_summary.txt" file does not exist.
                file_path_temporary = os.path.join(self.genome_directory, file_path_tophat.align_summary)
                if os.path.exists(file_path_temporary) and os.path.getsize(file_path_temporary) > 0:
                    executable_run_tophat.submit = False

                # TODO: End of code block.

                # Create a new rnaseq_process_tophat Executable.

                executable_process_tophat = stage_process_tophat.add_executable(
                    executable=Executable(
                        name='_'.join((stage_process_tophat.name, paired_reads_name)),
                        program='bsf_rnaseq_process_tophat2.sh'))
                # Set dependencies on previous Runnable or Executable objects.
                executable_process_tophat.dependencies.append(executable_run_tophat.name)

                # Set rnaseq_process_tophat options.

                # Set rnaseq_process_tophat arguments.
                self.set_command_configuration(command=executable_process_tophat)
                executable_process_tophat.arguments.append(file_path_tophat.output_directory)
                executable_process_tophat.arguments.append(genome_sizes_path)

                # Only submit this Executable if the "accepted_hits.bam.bai" file does not exist.
                file_path_temporary = os.path.join(
                    self.genome_directory,
                    file_path_tophat.accepted_hits_bai_link_source)
                if os.path.exists(file_path_temporary) and os.path.getsize(file_path_temporary):
                    executable_process_tophat.submit = False

                # TODO: Switch from an external Bash script to a set of Runnable and RunnableStep objects.
                # Since the Bash script includes Perl code to reset the BED score field to 0, rather than
                # re-scale it properly, it would be good to write a new bsf.runnables.process_tophat module
                # to implement this in Python code.

                # Create a process_tophat Runnable.

                # file_path_process_tophat = FilePathProcessTophat(prefix=prefix_tophat)
                #
                # runnable_process_tophat = self.add_runnable(runnable=Runnable(
                #     name='rnaseq_process_tophat',
                #     code_module='bsf.runnables.generic',
                #     working_directory=self.genome_directory,
                #     file_path_object=file_path_process_tophat,
                #     debug=self.debug))

                # Create a Cufflinks Runnable per paired_reads_name.

                prefix_run_cufflinks = '_'.join((stage_run_cufflinks.name, paired_reads_name))

                prefix_cufflinks = '_'.join(('rnaseq_cufflinks', paired_reads_name))
                # The output directory deviates from the prefix_run_cufflinks that itself is based on
                # stage_run_cufflinks.name. Both rnaseq_run_cufflinks and rnaseq_process_cufflinks processes
                # should use the same rnaseq_cufflinks directory.

                file_path_cufflinks = FilePathCufflinks(prefix=prefix_cufflinks)

                runnable_run_cufflinks = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_run_cufflinks,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        file_path_object=file_path_cufflinks,
                        debug=self.debug))
                executable_run_cufflinks = self.set_stage_runnable(
                    stage=stage_run_cufflinks,
                    runnable=runnable_run_cufflinks)
                # Set dependencies for subsequent Runnable or Executable objects.
                runnable_run_cufflinks_list.append(runnable_run_cufflinks)
                # Set dependencies on previous Runnable or Executable objects.
                executable_run_cufflinks.dependencies.append(executable_run_tophat.name)

                # Create a new Cufflinks RunnableStep.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='cufflinks',
                        program='cufflinks'))
                """ @type runnable_step: RunnableStep """

                # Cufflinks has a GTF option, in which case it will not assemble
                # novel transcripts and a GTF-guide option in which case it will
                # assemble novel transcripts.

                if self.novel_transcripts:
                    runnable_step.add_option_long(
                        key='GTF-guide',
                        value=self.transcriptome_gtf_path)
                else:
                    runnable_step.add_option_long(
                        key='GTF',
                        value=self.transcriptome_gtf_path)

                if self.mask_gtf_path:
                    runnable_step.add_option_long(
                        key='mask-file',
                        value=self.mask_gtf_path)

                runnable_step.add_option_long(
                    key='frag-bias-correct',
                    value=self.genome_fasta_path)

                if self.multi_read_correction:
                    runnable_step.add_switch_long(
                        key='multi-read-correct')

                if self.library_type:
                    runnable_step.add_option_long(
                        key='library-type',
                        value=self.library_type)

                # Cufflinks has a --library-norm-method option, but only one option (classic-fpkm) seems supported.

                if self.no_length_correction:
                    runnable_step.add_switch_long(
                        key='no-length-correction')

                runnable_step.add_option_long(
                    key='output-dir',
                    value=file_path_cufflinks.output_directory)

                runnable_step.add_option_long(
                    key='num-threads',
                    value=str(stage_run_cufflinks.threads))

                runnable_step.add_switch_long(
                    key='quiet')

                runnable_step.add_switch_long(
                    key='no-update-check')

                # Set Cufflinks arguments.

                runnable_step.arguments.append(file_path_tophat.accepted_hits_bam)

                # Convert the resulting transcripts GTF file into a UCSC genePred file.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='gtf_to_gp',
                        program='gtfToGenePred'))
                """ @type runnable_step: RunnableStep """
                runnable_step.add_switch_short(key='genePredExt')
                runnable_step.arguments.append(file_path_cufflinks.transcripts_gtf)
                runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)

                # Convert the UCSC genePred into a UCSC bigGenePred file.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='gp_to_bgp',
                        program='genePredToBigGenePred',
                        obsolete_file_path_list=[
                            file_path_cufflinks.temporary_gene_prediction,
                        ]))
                """ @type runnable_step: RunnableStep """
                runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)
                runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)

                # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bed_sort',
                        program='bedSort',
                        obsolete_file_path_list=[
                            file_path_cufflinks.temporary_big_gene_prediction,
                        ]))
                """ @type runnable_step: RunnableStep """
                runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)
                runnable_step.arguments.append(file_path_cufflinks.temporary_sorted_tsv)

                # Convert the sorted UCSC bigGenePred into a bigBed file.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bgp_to_bb',
                        program='bedToBigBed',
                        obsolete_file_path_list=[
                            file_path_cufflinks.temporary_sorted_tsv,
                        ]))
                """ @type runnable_step: RunnableStep """
                # FIXME: It would be good to allow options with and without an equal sign.
                # i.e. --type=bed12+8 versus --type bed12+8, which does not work.
                # How could this work with configuration.ini files?
                # FIXME: This needs to be configurable.
                # Either another instance variable is required or via a generic configuration option
                # for the RunnableStep.
                # TODO: The location of the autoSQL file needs to be configurable.
                runnable_step.add_switch_short(key='as=/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
                runnable_step.add_switch_short(key='tab')
                runnable_step.add_switch_short(key='type=bed12+8')
                runnable_step.arguments.append(file_path_cufflinks.temporary_sorted_tsv)
                runnable_step.arguments.append(genome_sizes_path)
                runnable_step.arguments.append(file_path_cufflinks.transcripts_bb)

                # Add a symbolic link for the transcripts bigBed file, that includes a sample name prefix.
                runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='link_transcripts_bb',
                        source_path=file_path_cufflinks.transcripts_bb_link_source,
                        target_path=file_path_cufflinks.transcripts_bb_link_target))

                # Add a symbolic link for the transcripts GTF file, that includes a sample name prefix.
                runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='link_transcripts_gtf',
                        source_path=file_path_cufflinks.transcripts_gtf_link_source,
                        target_path=file_path_cufflinks.transcripts_gtf_link_target))

        # Create one process_cufflinks Executable to process all sub-directories.

        if len(runnable_run_cufflinks_list):
            executable_process_cufflinks = stage_process_cufflinks.add_executable(
                executable=Executable(
                    name=stage_process_cufflinks.name,
                    program='bsf_rnaseq_process_cufflinks.R'))
            # Set dependencies on previous Runnable or Executable objects.
            for runnable_run_cufflinks in runnable_run_cufflinks_list:
                executable_process_cufflinks.dependencies.append(runnable_run_cufflinks.name)

            # Set process_cufflinks options.
            self.set_command_configuration(command=executable_process_cufflinks)
            executable_process_cufflinks.add_option_long(
                key='gtf-reference',
                value=self.transcriptome_gtf_path)
            executable_process_cufflinks.add_option_long(
                key='genome-version',
                value=self.genome_version)

        # The Cuffmerge process generates temporary files in the working directory that
        # have the same name for each comparison. If more than one Cuffmerge process runs at the same time,
        # contention occurs where more than one process writes to the same file and one process deletes files
        # expected to exist by another process.
        # Circumvent such a situation by introducing dependencies on previous Cuffmerge processes. Sigh.
        # TODO: Report this to the Cufflinks author.
        previous_cuffmerge_name = str()

        comparison_keys = self.comparisons.keys()
        comparison_keys.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_key in comparison_keys:
            # TODO: Should the comparison prefix also include the project name or number?
            prefix_cuffmerge = self.get_prefix_rnaseq_run_cuffmerge(comparison_name=comparison_key)

            file_path_cuffmerge = FilePathCuffmerge(prefix=prefix_cuffmerge)

            runnable_run_cuffmerge = self.add_runnable(
                runnable=Runnable(
                    name=prefix_cuffmerge,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_object=file_path_cuffmerge,
                    debug=self.debug))
            executable_run_cuffmerge = self.set_stage_runnable(
                stage=stage_run_cuffmerge,
                runnable=runnable_run_cuffmerge)
            # Set a dependency on the previous Cuffmerge process to avoid file contention.
            if previous_cuffmerge_name:
                executable_run_cuffmerge.dependencies.append(previous_cuffmerge_name)
            previous_cuffmerge_name = executable_run_cuffmerge.name

            # Create a new Cuffmerge RunnableStep.

            runnable_step_cuffmerge = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='cuffmerge',
                    program='cuffmerge'))
            """ @type runnable_step_cuffmerge: RunnableStep """

            # Set rnaseq_cuffmerge options.

            runnable_step_cuffmerge.add_option_long(
                key='output-dir',
                value=file_path_cuffmerge.output_directory)
            runnable_step_cuffmerge.add_option_long(
                key='num-threads',
                value=str(stage_run_cuffmerge.threads))
            runnable_step_cuffmerge.add_option_long(
                key='ref-gtf',
                value=self.transcriptome_gtf_path)
            runnable_step_cuffmerge.add_option_long(
                key='ref-sequence',
                value=self.genome_fasta_path)

            # Set rnaseq_cuffmerge arguments.

            # Create an assembly manifest file to merge all replicates of each Sample object.
            # This file requires an absolute path, because the working directory is not set at the stage of
            # job submission.

            assembly_path = os.path.join(self.genome_directory, file_path_cuffmerge.assembly_txt)
            assembly_file = open(assembly_path, 'w')

            # Process rnaseq_cuffmerge and rnaseq_cuffdiff arguments in parallel.

            # Create a Python list of Python list objects of Cuffquant abundances per comparison group.
            cuffdiff_cuffnorm_abundances = list()
            """ @type cuffdiff_cuffnorm_abundances: list[list[str | unicode]] """
            cuffdiff_cuffnorm_alignments = list()
            """ @type cuffdiff_cuffnorm_alignments: list[list[str | unicode]] """
            cuffdiff_cuffnorm_dependencies = list()
            """ @type cuffdiff_cuffnorm_dependencies: list[str] """
            cuffdiff_cuffnorm_labels = list()
            """ @type cuffdiff_cuffnorm_labels: list[str] """

            for group_name, group_samples in self.comparisons[comparison_key]:
                """ @type group_name: str """
                """ @type group_samples: list[Sample] """
                per_group_abundances_list = list()
                """ @type per_group_abundances_list: list[str | unicode] """
                per_group_alignments_list = list()
                """ @type per_group_alignments_list: list[str | unicode] """
                # Count samples that remain after removing excluded PairedReads objects.
                sample_count = 0

                for sample in group_samples:
                    """ @type sample: bsf.ngs.Sample """
                    paired_reads_dict = sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping,
                        exclude=True)

                    paired_reads_name_list = paired_reads_dict.keys()
                    if len(paired_reads_name_list):
                        sample_count += 1
                    else:
                        # Skip Sample objects, which PairedReads objects have all been excluded.
                        continue
                    paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                    for paired_reads_name in paired_reads_name_list:
                        # Add the Cufflinks assembled transcripts to the Cuffmerge manifest.

                        transcripts_path = os.path.join(
                            self.genome_directory,
                            '_'.join(('rnaseq_cufflinks', paired_reads_name)),
                            'transcripts.gtf')
                        assembly_file.write(transcripts_path + '\n')

                        # Wait for each Cufflinks replicate to finish, before Cuffmerge can run.

                        executable_run_cuffmerge.dependencies.append(
                            '_'.join((self.stage_name_run_cufflinks, paired_reads_name)))

                        # Create a Cuffquant Runnable per comparison (comparison_key) and replicate (paired_reads_name)
                        # on the basis of the Cuffmerge GTF file.

                        prefix_cuffquant = '_'.join((stage_run_cuffquant.name, comparison_key, paired_reads_name))

                        file_path_cuffquant = FilePathCuffquant(prefix=prefix_cuffquant)
                        # TODO: Remove after rearranging this method into the run method.
                        file_path_cuffquant.tophat_accepted_hits = \
                            os.path.join('_'.join(('rnaseq_tophat', paired_reads_name)), 'accepted_hits.bam')

                        runnable_run_cuffquant = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_cuffquant,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                file_path_object=file_path_cuffquant,
                                debug=self.debug))
                        executable_run_cuffquant = self.set_stage_runnable(
                            stage=stage_run_cuffquant,
                            runnable=runnable_run_cuffquant)

                        # Each Cuffquant process depends on Cuffmerge.

                        executable_run_cuffquant.dependencies.append(executable_run_cuffmerge.name)

                        # Create a new cuffquant RunnableStep.

                        runnable_step_cuffquant = runnable_run_cuffquant.add_runnable_step(
                            runnable_step=RunnableStep(
                                name='cuffquant',
                                program='cuffquant'))
                        """ @type runnable_step_cuffquant: RunnableStep """

                        # Set Cuffquant options.

                        runnable_step_cuffquant.add_option_long(
                            key='output-dir',
                            value=file_path_cuffquant.output_directory)
                        runnable_step_cuffquant.add_option_long(
                            key='num-threads',
                            value=str(stage_run_cuffquant.threads))
                        if self.mask_gtf_path:
                            runnable_step_cuffquant.add_option_long(
                                key='mask-file',
                                value=self.mask_gtf_path)
                        runnable_step_cuffquant.add_option_long(
                            key='frag-bias-correct',
                            value=self.genome_fasta_path)
                        if self.multi_read_correction:
                            runnable_step_cuffquant.add_switch_long(
                                key='multi-read-correct')
                        if self.library_type:
                            runnable_step_cuffquant.add_option_long(
                                key='library-type',
                                value=self.library_type)
                        if self.no_length_correction:
                            runnable_step_cuffquant.add_switch_long(
                                key='no-length-correction')
                        runnable_step_cuffquant.add_switch_long(
                            key='quiet')
                        runnable_step_cuffquant.add_switch_long(
                            key='no-update-check')

                        # Set Cuffquant arguments.
                        # Add the Cuffmerge GTF file and the TopHat BAM file as Cuffquant arguments.

                        runnable_step_cuffquant.arguments.append(file_path_cuffmerge.merged_gtf)
                        runnable_step_cuffquant.arguments.append(file_path_cuffquant.tophat_accepted_hits)

                        # Add the Cuffquant abundances file to the Cuffdiff list.

                        per_group_abundances_list.append(file_path_cuffquant.abundances)

                        # Add the TopHat BAM file to the Cuffdiff alignments list.

                        per_group_alignments_list.append(file_path_cuffquant.tophat_accepted_hits)

                        # Add the Cuffquant Runnable process name to the Cuffdiff and Cuffnorm dependencies list.

                        cuffdiff_cuffnorm_dependencies.append(executable_run_cuffquant.name)

                if sample_count:
                    cuffdiff_cuffnorm_labels.append(group_name)
                    cuffdiff_cuffnorm_abundances.append(per_group_abundances_list)
                    cuffdiff_cuffnorm_alignments.append(per_group_alignments_list)

            assembly_file.close()

            # Add the assembly manifest file as Cuffmerge argument.
            runnable_step_cuffmerge.arguments.append(assembly_path)

            # Convert the resulting merged GTF file into a UCSC genePred file.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='gtf_to_gp',
                    program='gtfToGenePred'))
            """ @type runnable_step: RunnableStep """
            runnable_step.add_switch_short(key='genePredExt')
            runnable_step.arguments.append(file_path_cuffmerge.merged_gtf)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)

            # Convert the UCSC genePred into a UCSC bigGenePred file.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='gp_to_bgp',
                    program='genePredToBigGenePred',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.temporary_gene_prediction,
                    ]))
            """ @type runnable_step: RunnableStep """
            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)

            # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='bed_sort',
                    program='bedSort',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.temporary_big_gene_prediction,
                    ]))
            """ @type runnable_step: RunnableStep """
            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)

            # Convert the sorted UCSC bigGenePred into a bigBed file.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='bgp_to_bb',
                    program='bedToBigBed',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.temporary_sorted_tsv,
                    ]))
            """ @type runnable_step: RunnableStep """
            # FIXME: It would be good to allow options with and without an equal sign.
            # i.e. --type=bed12+8 versus --type bed12+8, which does not work.
            # How could this work with configuration.ini files?
            # FIXME: This needs to be configurable.
            # Either another instance variable is required or via a generic configuration option
            # for the RunnableStep.
            # TODO: The location of the autoSQL file needs to be configurable.
            runnable_step.add_switch_short(key='as=/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
            runnable_step.add_switch_short(key='tab')
            runnable_step.add_switch_short(key='type=bed12+8')
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)
            runnable_step.arguments.append(genome_sizes_path)
            runnable_step.arguments.append(file_path_cuffmerge.merged_bb)

            # Add a symbolic link for the merged bigBed file, that includes a comparison name prefix.
            runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStepLink(
                    name='link_merged_bb',
                    source_path=file_path_cuffmerge.merged_bb_link_source,
                    target_path=file_path_cuffmerge.merged_bb_link_target))

            # Add a symbolic link for the merged GTF file, that includes a comparison name prefix.
            runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStepLink(
                    name='link_merged_gtf',
                    source_path=file_path_cuffmerge.merged_gtf_link_source,
                    target_path=file_path_cuffmerge.merged_gtf_link_target))

            # Create a Cuffnorm Runnable per comparison.

            prefix_run_cuffnorm = '_'.join((stage_run_cuffnorm.name, comparison_key))

            file_path_run_cuffnorm = FilePathCuffnorm(prefix=prefix_run_cuffnorm)

            runnable_run_cuffnorm = self.add_runnable(
                runnable=Runnable(
                    name=prefix_run_cuffnorm,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_object=file_path_run_cuffnorm,
                    debug=self.debug))
            executable_run_cuffnorm = self.set_stage_runnable(
                stage=stage_run_cuffnorm,
                runnable=runnable_run_cuffnorm)

            executable_run_cuffnorm.dependencies.extend(cuffdiff_cuffnorm_dependencies)

            # Create a new Cuffnorm RunnableStep.

            runnable_step_cuffnorm = runnable_run_cuffnorm.add_runnable_step(
                runnable_step=RunnableStep(
                    name='cuffnorm',
                    program='cuffnorm'))
            """ @type runnable_step: RunnableStep """

            # Set Cuffnorm options.

            runnable_step_cuffnorm.add_option_long(
                key='output-dir',
                value=file_path_run_cuffnorm.output_directory)
            runnable_step_cuffnorm.add_option_long(
                key='labels',
                value=','.join(cuffdiff_cuffnorm_labels))
            runnable_step_cuffnorm.add_option_long(
                key='num-threads',
                value=str(stage_run_cuffnorm.threads))
            if self.library_type:
                runnable_step_cuffnorm.add_option_long(
                    key='library-type',
                    value=self.library_type)
            runnable_step_cuffnorm.add_switch_long(
                key='quiet')
            runnable_step_cuffnorm.add_switch_long(
                key='no-update-check')

            # Add the Cuffmerge GTF file as first Cuffnorm argument.
            runnable_step_cuffnorm.arguments.append(file_path_cuffmerge.merged_gtf)

            # Add the Cuffquant abundances files per point as Cuffnorm arguments.
            for per_group_abundances_list in cuffdiff_cuffnorm_abundances:
                runnable_step_cuffnorm.arguments.append(','.join(per_group_abundances_list))

            # Create a Cuffdiff Runnable per comparison.

            prefix_run_cuffdiff = '_'.join((stage_run_cuffdiff.name, comparison_key))

            file_path_run_cuffdiff = FilePathCuffdiff(prefix=prefix_run_cuffdiff)

            runnable_run_cuffdiff = self.add_runnable(
                runnable=Runnable(
                    name=prefix_run_cuffdiff,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_object=file_path_run_cuffdiff,
                    debug=self.debug))
            executable_run_cuffdiff = self.set_stage_runnable(
                stage=stage_run_cuffdiff,
                runnable=runnable_run_cuffdiff)

            if run_cuffquant_before_cuffdiff:
                # Add all executable_run_cuffquant dependencies to the executable_run_cuffdiff process.
                executable_run_cuffdiff.dependencies.extend(cuffdiff_cuffnorm_dependencies)
            else:
                # Add the executable_run_cuffmerge dependency to the executable_run_cuffdiff process.
                executable_run_cuffdiff.dependencies.append(executable_run_cuffmerge.name)

            # Create a new Cuffdiff RunnableStep.

            runnable_step_cuffdiff = runnable_run_cuffdiff.add_runnable_step(
                runnable_step=RunnableStep(
                    name='cuffdiff',
                    program='cuffdiff'))
            """ @type runnable_step: RunnableStep """

            # Set Cuffdiff options.

            runnable_step_cuffdiff.add_option_long(
                key='output-dir',
                value=file_path_run_cuffdiff.output_directory)
            runnable_step_cuffdiff.add_option_long(
                key='labels',
                value=','.join(cuffdiff_cuffnorm_labels))
            runnable_step_cuffdiff.add_option_long(
                key='num-threads',
                value=str(stage_run_cuffdiff.threads))
            if self.mask_gtf_path:
                runnable_step_cuffdiff.add_option_long(
                    key='mask-file',
                    value=self.mask_gtf_path)
            runnable_step_cuffdiff.add_option_long(
                key='frag-bias-correct',
                value=self.genome_fasta_path)
            if self.multi_read_correction:
                runnable_step_cuffdiff.add_switch_long(
                    key='multi-read-correct')
            if self.library_type:
                runnable_step_cuffdiff.add_option_long(
                    key='library-type',
                    value=self.library_type)
            if self.no_length_correction:
                runnable_step_cuffdiff.add_switch_long(
                    key='no-length-correction')
            runnable_step_cuffdiff.add_switch_long(
                key='quiet')
            runnable_step_cuffdiff.add_switch_long(
                key='no-update-check')

            # Add the Cuffmerge GTF file as first Cuffdiff argument.
            runnable_step_cuffdiff.arguments.append(file_path_cuffmerge.merged_gtf)

            # Cuffdiff seems to have a problem with Cuffquant abundances files in that the isoforms.count_tracking
            # files show ridiculously low numbers such as 1e-319 for some splice variants. Usually, other splice
            # variants in the same cluster seem fine.
            if run_cuffquant_before_cuffdiff:
                # Add the Cuffquant abundances files per comparison group as Cuffdiff arguments.
                for per_group_abundances_list in cuffdiff_cuffnorm_abundances:
                    runnable_step_cuffdiff.arguments.append(','.join(per_group_abundances_list))
            else:
                # Add the TopHat BAM files per comparison group as Cuffdiff arguments.
                for per_group_alignments_list in cuffdiff_cuffnorm_alignments:
                    runnable_step_cuffdiff.arguments.append(','.join(per_group_alignments_list))

            # Create a new rnaseq_process_cuffdiff Executable.

            executable_process_cuffdiff = stage_process_cuffdiff.add_executable(
                executable=Executable(
                    name='_'.join((stage_process_cuffdiff.name, comparison_key)),
                    program='bsf_rnaseq_process_cuffdiff.R'))
            executable_process_cuffdiff.dependencies.append(executable_run_cuffdiff.name)

            # Set rnaseq_process_cuffdiff options.
            self.set_command_configuration(command=executable_process_cuffdiff)
            executable_process_cuffdiff.add_option_long(
                key='comparison-name',
                value=comparison_key)
            executable_process_cuffdiff.add_option_long(
                key='gtf-assembly',
                value=file_path_cuffmerge.merged_gtf)
            executable_process_cuffdiff.add_option_long(
                key='gtf-reference',
                value=self.transcriptome_gtf_path)
            executable_process_cuffdiff.add_option_long(
                key='genome-version',
                value=self.genome_version)

            # Set rnaseq_process_cuffdiff arguments.

            # None so far.

        return

    def report(self):
        """Create a C{bsf.analyses.rna_seq.Tuxedo} report in HTML format and a UCSC Genome Browser Track Hub.
        @return:
        @rtype:
        """

        def relative_anchor(prefix, suffix, text):
            """Create a relative HTML anchor element.

            <a href="prefix/prefix_suffix">text</a>
            @param prefix: Prefix
            @type prefix: str
            @param suffix: Suffix
            @type suffix: str
            @param text: Link text
            @type text: str
            @return: Relative URL
            @rtype: str
            """
            return '<a href="' + prefix + '/' + prefix + '_' + suffix + '">' + text + '</a>'

        def relative_image(prefix, suffix, text):
            """Create a relative HTML image element

            <img alt="text" src="prefix/prefix_suffix" height="80" width="80" />
            @param prefix: Prefix
            @type prefix: str
            @param suffix: Suffix
            @type suffix: str
            @param text: Alternative text
            @type text: str
            @return:
            """
            return '<img' + \
                   ' alt="' + text + '"' + \
                   ' src="' + prefix + '/' + prefix + '_' + suffix + '"' + \
                   ' height="80"' + \
                   ' width="80"' + \
                   ' />'

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        # This code only needs the public URL.

        hub_list = list()
        """ @type hub_list: list[str | unicode] """

        # Write a HTML document.

        report_list = list()
        """ @type report_list: list[str | unicode] """

        report_list += '<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n'
        report_list += '\n'

        # TopHat and Cufflinks table.

        report_list += '<h2 id="transcriptome_browsing">Transcriptome Browsing</h2>\n'
        report_list += '\n'

        report_list += '<h3 id="read_alignments">Read Alignments</h3>\n'
        report_list += '\n'

        report_list += '<p id ="tophat">\n'
        # http://tophat.cbcb.umd.edu/manual.html
        report_list += '<strong><a href="http://ccb.jhu.edu/software/tophat/index.shtml">TopHat</a></strong> '
        report_list += 'aligns RNA-Seq reads to a genome in order to identify '
        report_list += 'exon-exon splice junctions. It is built on the ultra fast\n'
        report_list += 'short read mapping program\n'
        report_list += '<strong>'
        report_list += '<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie 2</a>'
        report_list += '</strong>.\n'
        report_list += '</p>\n'
        report_list += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = {
            'db': self.genome_version,
            'hubUrl': '/'.join((Default.url_absolute_projects(), link_name, 'rnaseq_hub.txt')),
        }

        report_list += '<p id="track_hub">\n'
        report_list += 'View TopHat <strong>read alignments</strong> tracks for each sample\n'
        report_list += 'in their genomic context via the project-specific\n'
        report_list += 'UCSC Genome Browser Track Hub '
        report_list += '<a href="' + self.ucsc_track_url(options_dict=options_dict) + '">' \
                       + self.project_name + '</a>.\n'
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<p>\n'
        report_list += '<a href="rnaseq_tophat_alignment_summary.pdf">'
        report_list += '<img ' \
                       'alt="TopHat Alignment Summary" ' \
                       'id="tophat_alignment_summary_img" ' \
                       'src="rnaseq_tophat_alignment_summary.png" ' \
                       'height="80" ' \
                       'width="80" ' \
                       '/>'
        report_list += '</a>\n'
        report_list += 'Alignment summary statistics <a href="rnaseq_tophat_alignment_summary.tsv">TSV</a>\n'
        report_list += '</p>\n'

        report_list += '<h3 id="alignment_events">Splice Junctions, Insertions and Deletions</h3>\n'
        report_list += '\n'

        report_list += '<p>\n'
        report_list += 'TopHat reports <strong>splice junctions</strong> on the basis of RNA-Seq\n'
        report_list += 'read alignments in UCSC BED track format.\n'
        report_list += 'Each junction consists of two connected BED blocks,\n'
        report_list += 'where each block is as long as the maximal overhang\n'
        report_list += 'of any read spanning the junction. The score is\n'
        report_list += 'the number of alignments spanning the junction.\n'
        report_list += 'UCSC BED tracks of <strong>insertions</strong> and\n'
        report_list += '<strong>deletions</strong> are also reported by TopHat.\n'
        report_list += '</p>\n'

        report_list += '<p>\n'
        report_list += 'View the corresponding TopHat tracks for junctions, deletions and insertions\n'
        report_list += 'for each sample in their genomic context via the project-specific\n'
        report_list += 'UCSC Genome Browser Track Hub '
        report_list += '<a href="' + self.ucsc_track_url(options_dict=options_dict) + '">' + \
                       self.project_name + '</a>.\n'
        report_list += '</p>\n'
        report_list += '\n'

        # report_list += '<p>\n'
        # report_list += 'Follow the links below to attach\n'
        # report_list += 'Tophat junction, deletion and insertion annotation to the\n'
        # report_list += 'UCSC Genome Browser. Since each file needs transferring to\n'
        # report_list += 'the UCSC site, subsequent pages will take some time to load.\n'
        # report_list += '</p>\n'

        report_list += '<h2 id="gene_expression_profiles">Gene Expression Profiles</h2>\n'
        report_list += '\n'

        report_list += '<p id="cufflinks">\n'
        # http://cufflinks.cbcb.umd.edu/howitworks.html
        report_list += '<strong><a href="http://cole-trapnell-lab.github.io/cufflinks/">Cufflinks</a></strong>\n'
        report_list += 'assembles aligned RNA-Seq reads into transcripts,\n'
        report_list += 'estimates their abundances, and tests for differential\n'
        report_list += 'expression and regulation transcriptome-wide.\n'
        report_list += 'It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of\n'
        report_list += 'transcripts. Cufflinks then estimates the relative abundances of these transcripts based\n'
        report_list += 'on how many reads support each one, taking into account biases in library preparation '
        report_list += 'protocols.\n'
        report_list += '</p>\n'

        report_list += '<p>\n'
        report_list += 'The Cufflinks <strong>assembled transcripts</strong> can be attached to the \n'
        report_list += 'UCSC Genome Browser, by following the "Transcript Assembly" links\n'
        report_list += 'below.\n'
        report_list += 'The isoforms.fpkm_tracking and genes.fpkm_tracking files\n'
        report_list += 'contain the estimated isoform or gene expression values in the generic\n'
        # http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
        report_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                       'fpkm-tracking-format">FPKM Tracking format</a>.\n'
        report_list += 'The isoforms.count_tracking and genes.count_tracking files\n'
        report_list += 'contain the scaled isoform or gene count values in the generic\n'
        report_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                       'count-tracking-format">Count Tracking format</a>.\n'
        report_list += '</p>\n'

        report_list += '<p>\n'
        report_list += 'Please see a more detailed description of\n'
        # http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output
        report_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                       'output-formats-used-in-the-cufflinks-suite">Cufflinks output</a>.\n'
        report_list += '</p>\n'

        report_list += '<table id="gene_expression_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Sample</th>\n'
        report_list += '<th>Assembled Transcripts</th>\n'
        report_list += '<th>Gene FPKM</th>\n'
        report_list += '<th>Transcript FPKM</th>\n'
        report_list += '<th>Genes (Symbols)</th>\n'
        report_list += '<th>Isoforms (Symbols)</th>\n'
        report_list += '<th>Aligned BAM file</th>\n'
        report_list += '<th>Aligned BAI file</th>\n'
        report_list += '<th>Unaligned BAM file</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        # Group via UCSC super tracks.

        hub_list += 'track Alignments\n'
        hub_list += 'shortLabel Alignments\n'
        hub_list += 'longLabel TopHat RNA-Seq read alignments\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group alignments\n'
        hub_list += '\n'

        hub_list += 'track Assemblies\n'
        hub_list += 'shortLabel Assemblies\n'
        hub_list += 'longLabel Cuffmerge transcript structures\n'
        hub_list += 'visibility full\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group assemblies\n'
        hub_list += '\n'

        hub_list += 'track Coverage\n'
        hub_list += 'shortLabel Coverage\n'
        hub_list += 'longLabel TopHat RNA-Seq alignment coverage\n'
        hub_list += 'visibility full\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group coverage\n'
        hub_list += '\n'

        hub_list += 'track Deletions\n'
        hub_list += 'shortLabel Deletions\n'
        hub_list += 'longLabel TopHat RNA-Seq deletions\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group alignments\n'
        hub_list += '\n'

        hub_list += 'track Insertions\n'
        hub_list += 'shortLabel Insertions\n'
        hub_list += 'longLabel TopHat RNA-Seq insertions\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group alignments\n'
        hub_list += '\n'

        hub_list += 'track Junctions\n'
        hub_list += 'shortLabel Junctions\n'
        hub_list += 'longLabel TopHat RNA-Seq splice junctions\n'
        hub_list += 'visibility show\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group alignments\n'
        hub_list += '\n'

        hub_list += 'track Transcripts\n'
        hub_list += 'shortLabel Transcripts\n'
        hub_list += 'longLabel Cufflinks transcript structures\n'
        hub_list += 'visibility show\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group transcripts\n'
        hub_list += '\n'

        for sample in self.sample_list:
            if self.debug > 0:
                print repr(self) + ' Sample name: ' + sample.name
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                # TopHat produces accepted_hits.bam, deletions.bb,
                # insertions.bb and junctions.bb files.

                #
                # Add a trackDB entry for each accepted_hits.bam file.
                #

                # Common trackDb settings.

                hub_list += 'track ' + paired_reads_name + '_alignments\n'
                hub_list += 'type bam\n'
                hub_list += 'shortLabel ' + paired_reads_name + '_alignments\n'
                hub_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq read alignments\n'
                hub_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/accepted_hits.bam\n'
                hub_list += 'visibility dense\n'
                # hub_list += 'html ...\n'

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # Compressed Sequence Alignment track settings.

                # None so far.

                # Composite track settings.

                hub_list += 'parent Alignments\n'
                hub_list += '\n'

                #
                # Add a trackDB entry for each accepted_hits.bw file.
                #

                # Common trackDB settings.

                hub_list += 'track ' + paired_reads_name + '_coverage\n'
                # TODO: The bigWig type must declare the expected signal range.
                # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                hub_list += 'type bigWig\n'
                hub_list += 'shortLabel ' + paired_reads_name + '_coverage\n'
                hub_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq alignment coverage\n'
                hub_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/accepted_hits.bw\n'
                hub_list += 'visibility full\n'
                # hub_list += 'html ...\n'

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # bigWig - Signal graphing track settings.

                hub_list += 'alwaysZero on\n'
                hub_list += 'autoScale on\n'
                hub_list += 'graphTypeDefault bar\n'
                hub_list += 'maxHeightPixels 100:60:20\n'
                # hub_list += 'maxWindowToQuery 10000000\n'
                # hub_list += 'smoothingWindow 5\n'
                # hub_list += 'transformFunc NONE\n'
                # hub_list += 'viewLimits 0:45\n'
                # hub_list += 'viewLimitsMax 0:50\n'
                # hub_list += 'windowingFunction maximum\n'
                # hub_list += 'yLineMark <#>\n'
                # hub_list += 'yLineOnOff on \n'
                # hub_list += 'gridDefault on\n'

                # Composite track settings.

                hub_list += 'parent Coverage\n'
                hub_list += 'centerLabelsDense off\n'
                hub_list += '\n'

                #
                # Add a trackDB entry for each deletions.bb file.
                #

                hub_list += 'track ' + paired_reads_name + '_deletions\n'
                hub_list += 'type bigBed\n'
                hub_list += 'shortLabel ' + paired_reads_name + '_deletions\n'
                hub_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq deletions\n'
                hub_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/deletions.bb\n'
                hub_list += 'visibility hide\n'
                # hub_list += 'html ...\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # Composite track settings.

                hub_list += 'parent Deletions\n'
                hub_list += '\n'

                # Insertions

                hub_list += 'track insertions_' + paired_reads_name + '\n'
                hub_list += 'type bigBed\n'
                hub_list += 'shortLabel ' + paired_reads_name + '_insertions\n'
                hub_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq insertions\n'
                hub_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/insertions.bb\n'
                hub_list += 'visibility hide\n'
                # hub_list += 'html ...\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # Composite track settings.

                hub_list += 'parent Insertions\n'
                hub_list += '\n'

                # Junctions

                hub_list += 'track ' + paired_reads_name + '_junctions\n'
                hub_list += 'type bigBed\n'
                hub_list += 'shortLabel ' + paired_reads_name + '_junctions\n'
                hub_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq splice junctions\n'
                hub_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/junctions.bb\n'
                hub_list += 'visibility pack\n'
                # hub_list += 'html ...\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # Composite track settings.

                hub_list += 'parent Junctions\n'
                hub_list += '\n'

                # Transcripts

                hub_list += 'track ' + paired_reads_name + '_transcripts\n'
                hub_list += 'type bigGenePred\n'
                hub_list += 'shortLabel ' + paired_reads_name + '_transcripts\n'
                hub_list += 'longLabel ' + paired_reads_name + ' Cufflinks transcript assembly\n'
                hub_list += 'bigDataUrl rnaseq_cufflinks_' + paired_reads_name + '/transcripts.bb\n'
                hub_list += 'visibility hide\n'
                # hub_list += 'html ...\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # Composite track settings.

                hub_list += 'parent Transcripts\n'
                hub_list += '\n'

                # Cufflinks produces genes.fpkm_tracking, isoforms.fpkm_tracking,
                # skipped.gtf and transcripts.gtf.

                path_prefix = 'rnaseq_cufflinks_' + paired_reads_name

                report_list += '<tr>\n'
                report_list += '<td class="left">' + paired_reads_name + '</td>\n'
                report_list += '<td class="center">'
                report_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='transcripts.gtf',
                    text='Transcript Assembly')
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + path_prefix + '/genes.fpkm_tracking">Genes FPKM</a>'
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + path_prefix + '/isoforms.fpkm_tracking">Isoforms FPKM</a>'
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_fpkm_tracking.tsv',
                    text='Genes (Symbols)')
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_fpkm_tracking.tsv',
                    text='Isoforms (Symbols)')

                # TODO: The aligned BAM and BAI files and the unaligned BAM file are currently non standard.
                # The files have a "rnasq_tophat_" prefix, but are in the "rnaseq_cufflinks_" directory.
                # This will be resolved when the process topHta step gets re-engineered.
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + path_prefix + '/rnaseq_tophat_' + paired_reads_name + \
                               '_accepted_hits.bam">Aligned BAM</a>'
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + path_prefix + '/rnaseq_tophat_' + paired_reads_name + \
                               '_accepted_hits.bam.bai">Aligned BAI</a>'
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + path_prefix + '/rnaseq_tophat_' + paired_reads_name + \
                               '_unmapped.bam">Unaligned BAM</a>'
                report_list += '</td>\n'
                report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        # Cuffdiff produces cds_exp.diff, gene_exp.diff, isoform_exp.diff
        # promoters.diff, splicing.diff and tss_group_exp.diff amongst many others.

        report_list += '<h2 id="differential_expression">Differential Expression</h2>\n'
        report_list += '\n'

        report_list += '<p id="cuffdiff">\n'
        # http://cufflinks.cbcb.umd.edu/howitworks.html#diff
        report_list += '<strong><a href="http://cole-trapnell-lab.github.io/cufflinks/">Cufflinks</a></strong>\n'
        report_list += 'finds significant changes in transcript\n'
        report_list += 'expression, splicing, and promoter use.'
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<h3 id="all_genes">All Genes</h3>\n'

        report_list += '<table id="differential_expression_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Comparison</th>\n'
        report_list += '<th>Samples</th>\n'
        report_list += '<th>Replicates</th>\n'
        report_list += '<th>Coding Sequences</th>\n'
        report_list += '<th>Genes</th>\n'
        report_list += '<th>Isoforms</th>\n'
        report_list += '<th>Promoters</th>\n'
        report_list += '<th>Splicing</th>\n'
        report_list += '<th>Transcription Start Sites</th>\n'
        report_list += '<th>Gene FPKM Replicates</th>\n'
        report_list += '<th>Gene Count Replicates</th>\n'
        report_list += '<th>Isoform FPKM Replicates</th>\n'
        report_list += '<th>Isoform Count Replicates</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        comparison_keys = self.comparisons.keys()
        comparison_keys.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_key in comparison_keys:
            # Assemblies Super Track

            hub_list += 'track ' + comparison_key + '_assembly\n'
            hub_list += 'type bigGenePred\n'
            hub_list += 'shortLabel ' + comparison_key + '_assembly\n'
            hub_list += 'longLabel ' + comparison_key + ' Cufflinks transcript assembly\n'
            hub_list += 'bigDataUrl rnaseq_cuffmerge_' + comparison_key + '/merged.bb\n'
            hub_list += 'visibility pack\n'
            # hub_list += 'html ...\n'
            # 'html' is missing from the common settings.

            # Common optional settings.

            hub_list += 'color 0,0,0\n'

            # Composite track settings.

            hub_list += 'parent Assemblies\n'
            hub_list += '\n'

            path_prefix = 'rnaseq_process_cuffdiff_' + comparison_key

            # Link to comparison-specific symbolic links in the directory after cummeRbund processing.

            report_list += '<tr>\n'
            report_list += '<td class="left">' + comparison_key + '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='samples.tsv',
                text='Samples')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='replicates.tsv',
                text='Replicates')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='cds_exp_diff.tsv',
                text='Coding Sequences')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_exp_diff.tsv',
                text='<strong>Genes</strong>')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_exp_diff.tsv',
                text='Isoforms')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='promoters_diff.tsv',
                text='Promoters')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='splicing_diff.tsv',
                text='Splicing')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='tss_group_exp_diff.tsv',
                text='Transcription Start Sites')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_fpkm_replicates.tsv',
                text='Gene FPKM Replicates')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_counts_replicates.tsv',
                text='Gene Count Replicates')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_fpkm_replicates.tsv',
                text='Isoform FPKM Replicates')
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_counts_replicates.tsv',
                text='Isoform Count Replicates')
            report_list += '</td>\n'
            report_list += '</tr>\n'

            # Read sample pair information if available.

            sample_pair_path = os.path.join(
                self.genome_directory,
                path_prefix,
                '_'.join((path_prefix, 'sample_pairs.tsv')))

            if os.path.exists(sample_pair_path):

                sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                for row_dict in sample_pair_sheet.row_dicts:
                    report_list += '<tr>\n'
                    report_list += '<td></td>\n'  # Comparison
                    report_list += '<td class="left" colspan="3">'
                    report_list += '<strong>' + \
                                   row_dict['V1'] + \
                                   '</strong> versus <strong>' + \
                                   row_dict['V2'] + \
                                   '</strong>'
                    report_list += '</td>\n'  # Sample
                    report_list += '<td class="center">'
                    report_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_diff.tsv')),
                        text='<strong>Genes</strong>')
                    report_list += '</td>\n'
                    report_list += '<td class="center">'
                    report_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='_'.join((row_dict['V1'], row_dict['V2'], 'isoforms_diff.tsv')),
                        text='Isoforms')
                    report_list += '</td>\n'
                    report_list += '<td class="left" colspan="5"></td>\n'
                    report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        report_list += '<h3 id="significant_genes">Significant Genes</h3>\n'

        report_list += '<table id="significant_genes_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Comparison</th>\n'
        report_list += '<th>Genes</th>\n'
        report_list += '<th>Isoforms</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        for comparison_key in comparison_keys:
            path_prefix = 'rnaseq_process_cuffdiff_' + comparison_key

            report_list += '<tr>\n'
            report_list += '<td class="left">' + comparison_key + '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_significance_matrix.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_significance_matrix.png',
                    text='Significance Matrix Plot - Genes - ' + comparison_key))
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_significance_matrix.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_significance_matrix.png',
                    text='Significance Matrix Plot - Isoforms - ' + comparison_key))
            report_list += '</td>\n'
            report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'

        # Show cummeRbund quality plots.

        report_list += '<h2 id="quality_plots">Quality Plots</h2>\n'
        report_list += '\n'

        report_list += '<p>\n'
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<table id="quality_plots_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Comparison</th>\n'
        report_list += '<th>Dispersion Plot - Genes</th>\n'
        report_list += '<th>Dispersion Plot - Isoforms</th>\n'
        report_list += '<th>Squared Coefficient of Variation - Genes</th>\n'
        report_list += '<th>Squared Coefficient of Variation - Isoforms</th>\n'
        report_list += '<th>Density Plot without Replicates - Genes</th>\n'
        report_list += '<th>Density Plot with Replicates - Genes</th>\n'
        report_list += '<th>Density Plot without Replicates - Isoforms</th>\n'
        report_list += '<th>Density Plot with Replicates - Isoforms</th>\n'
        report_list += '<th>Box Plot without Replicates - Genes</th>\n'
        report_list += '<th>Box Plot with Replicates - Genes</th>\n'
        report_list += '<th>Box Plot without Replicates - Isoforms</th>\n'
        report_list += '<th>Box Plot with Replicates - Isoforms</th>\n'
        report_list += '<th>Scatter Matrix Plot - Genes</th>\n'
        report_list += '<th>Scatter Matrix Plot - Isoforms</th>\n'
        report_list += '<th>Dendrogram Plot</th>\n'
        report_list += '<th>Volcano Matrix Plot - Genes</th>\n'
        report_list += '<th>Multidimensional Scaling Plot - Genes</th>\n'
        report_list += '<th>Principal Component Analysis Plot - Genes</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        for comparison_key in comparison_keys:
            path_prefix = 'rnaseq_process_cuffdiff_' + comparison_key

            report_list += '<tr>\n'
            report_list += '<td class="left">' + comparison_key + '</td>\n'

            # Dispersion Plots for Genes and Isoforms

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_dispersion.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_dispersion.png',
                    text='Dispersion Plot - Genes - ' + comparison_key))
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_dispersion.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_dispersion.png',
                    text='Dispersion Plot - Isoforms - ' + comparison_key))
            report_list += '</td>\n'

            # Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms

            if os.path.exists(
                    path=os.path.join(self.genome_directory, path_prefix, path_prefix + '_genes_scv.png')):
                report_list += '<td class="center">'
                report_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_scv.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_scv.png',
                        text='Squared Coefficient of Variation (SCV) - Genes - ' + comparison_key))
                report_list += '</td>\n'
            else:
                report_list += '<td class="center"></td>\n'

            if os.path.exists(
                    path=os.path.join(self.genome_directory, path_prefix, path_prefix + '_isoforms_scv.png')):
                report_list += '<td class="center">'
                report_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_scv.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_scv.png',
                        text='Squared Coefficient of Variation (SCV) - Isoforms - ' + comparison_key))
                report_list += '</td>\n'
            else:
                report_list += '<td class="center"></td>\n'

            # Density Plots for Genes without and with Replicates

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_density_wo_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_density_wo_replicates.png',
                    text='Density Plot without Replicates - Genes - ' + comparison_key))
            report_list += '</td>\n'

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_density_w_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_density_w_replicates.png',
                    text='Density Plot with Replicates - Genes - ' + comparison_key))
            report_list += '</td>\n'

            # Density Plots for Isoforms without and with Replicates

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_density_wo_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_density_wo_replicates.png',
                    text='Density Plot without Replicates - Isoforms - ' + comparison_key))
            report_list += '</td>\n'

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_density_w_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_density_w_replicates.png',
                    text='Density Plot with Replicates - Isoforms - ' + comparison_key))
            report_list += '</td>\n'

            # Box Plots for Genes without and with Replicates

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_box_wo_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_box_wo_replicates.png',
                    text='Box Plot without Replicates - Genes - ' + comparison_key))
            report_list += '</td>\n'

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_box_w_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_box_w_replicates.png',
                    text='Box Plot with Replicates - Genes - ' + comparison_key))
            report_list += '</td>\n'

            # Box Plots for Isoforms with and without Replicates

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_box_wo_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_box_wo_replicates.png',
                    text='Box Plot without Replicates - Isoforms - ' + comparison_key))
            report_list += '</td>\n'

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_box_w_replicates.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_box_w_replicates.png',
                    text='Box Plot with Replicates - Isoforms - ' + comparison_key))
            report_list += '</td>\n'

            # Scatter Matrix Plot for Genes and Isoforms

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_scatter_matrix.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_scatter_matrix.png',
                    text='Scatter Matrix Plot - Genes - ' + comparison_key))
            report_list += '</td>\n'

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='isoforms_scatter_matrix.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='isoforms_scatter_matrix.png',
                    text='Scatter Matrix Plot - Isoforms - ' + comparison_key))
            report_list += '</td>\n'

            # Dendrogram Plot for Genes

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_dendrogram.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_dendrogram.png',
                    text='Dendrogram Plot - Genes - ' + comparison_key))
            report_list += '</td>\n'

            # Volcano Matrix Plot for Genes

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_volcano_matrix.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_volcano_matrix.png',
                    text='Volcano Matrix Plot - Genes - ' + comparison_key))
            report_list += '</td>\n'

            # Multidimensional Scaling Plot for Genes

            if os.path.exists(
                    path=os.path.join(self.genome_directory, path_prefix, path_prefix + '_genes_mds.png')):
                report_list += '<td class="center">'
                report_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_mds.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_mds.png',
                        text='Multidimensional Scaling Plot - Genes - ' + comparison_key))
                report_list += '</td>\n'
            else:
                report_list += '<td></td>\n'

            # Principal Component Analysis Plot for Genes

            report_list += '<td class="center">'
            report_list += relative_anchor(
                prefix=path_prefix,
                suffix='genes_pca.pdf',
                text=relative_image(
                    prefix=path_prefix,
                    suffix='genes_pca.png',
                    text='Principal Component Analysis Plot - Genes - ' + comparison_key))
            report_list += '</td>\n'

            report_list += '</tr>\n'

            # Read sample pair information if available.

            sample_pair_path = os.path.join(
                self.genome_directory,
                path_prefix,
                '_'.join((path_prefix, 'sample_pairs.tsv')))

            if os.path.exists(sample_pair_path):
                sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                for row_dict in sample_pair_sheet.row_dicts:
                    report_list += '<tr>\n'

                    report_list += '<td class="left"></td>\n'
                    report_list += '<td  class="left" colspan="10">'
                    report_list += '<strong>' + \
                                   row_dict['V1'] + \
                                   '</strong> versus <strong>' + \
                                   row_dict['V2'] + \
                                   '</strong>'
                    report_list += '</td>\n'

                    report_list += '<td class="center">'
                    report_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_scatter.pdf')),
                        text=relative_image(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_scatter.png')),
                            text='Scatter Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2']))
                    report_list += '</td>\n'

                    report_list += '<td class="center"></td>\n'

                    report_list += '<td class="center">'
                    report_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='_'.join((row_dict['V1'], row_dict['V2'], 'maplot.pdf')),
                        text=relative_image(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'maplot.png')),
                            text='M vs A Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2']))
                    report_list += '</td>\n'

                    report_list += '<td class="center">'
                    report_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_volcano.pdf')),
                        text=relative_image(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_volcano.png')),
                            text='Volcano Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2']))
                    report_list += '</td>\n'

                    report_list += '<td  class="center" colspan="4"></td>\n'

                    report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        self.report_to_file(content=report_list)
        self.ucsc_hub_to_file(content=hub_list)

        return


class FilePathDESeq(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathDESeq} models files in a comparison-specific DESeq directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathDESeq} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathDESeq, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class DESeq(Analysis):
    """DESeq RNASeq C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_deseq: C{bsf.Stage.name} for the run DESeq stage
    @type stage_name_deseq: str
    """

    name = 'RNA-seq Analysis'
    prefix = 'rnaseq'

    # Replicate stage
    stage_name_deseq = '_'.join((prefix, 'deseq'))

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
            cmp_file=None,
            genome_fasta_path=None,
            transcriptome_gtf_path=None,
            transcriptome_index_path=None,
            design_formula=None,
            mask_gtf_path=None,
            multi_read_correction=None,
            library_type=None,
            no_length_correction=False):
        """Initialise a C{bsf.analyses.rna_seq.DESeq} object.

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
        @type comparisons: dict[str, (str, list[bsf.ngs.Sample])]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode
        @param design_formula: Design formula
        @type design_formula: str | unicode
        @param mask_gtf_path: GTF file path to mask transcripts
        @type mask_gtf_path: str | unicode
        @param multi_read_correction: Apply multi-read correction
        @type multi_read_correction: bool
        @param library_type: Library type
            Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
        @type library_type: str
        @param no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
        @type no_length_correction: bool
        @return:
        @rtype:
        """

        super(DESeq, self).__init__(
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

        if cmp_file is None:
            self.cmp_file = str()
        else:
            self.cmp_file = cmp_file

        if genome_fasta_path is None:
            self.genome_fasta_path = str()
        else:
            self.genome_fasta_path = genome_fasta_path

        if transcriptome_gtf_path is None:
            self.transcriptome_gtf_path = str()
        else:
            self.transcriptome_gtf_path = transcriptome_gtf_path

        if transcriptome_index_path is None:
            self.transcriptome_index_path = str()
        else:
            self.transcriptome_index_path = transcriptome_index_path

        if design_formula is None:
            self.design_formula = str()
        else:
            self.design_formula = design_formula

        if mask_gtf_path is None:
            self.mask_gtf_path = str()
        else:
            self.mask_gtf_path = mask_gtf_path

        if multi_read_correction is None:
            self.multi_read_correction = False
        else:
            assert isinstance(multi_read_correction, bool)
            self.multi_read_correction = multi_read_correction

        if library_type is None:
            self.library_type = str()
        else:
            self.library_type = library_type

        if no_length_correction is None:
            self.no_length_correction = False
        else:
            assert isinstance(no_length_correction, bool)
            self.no_length_correction = no_length_correction

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

        super(DESeq, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cmp_file = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'design_formula'
        if configuration.config_parser.has_option(section=section, option=option):
            self.design_formula = configuration.config_parser.get(section=section, option=option)

        option = 'mask_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mask_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'multi_read_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.multi_read_correction = configuration.config_parser.getboolean(section=section, option=option)

        option = 'library_type'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_type = configuration.config_parser.get(section=section, option=option)

        option = 'no_length_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.no_length_correction = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.Tuxedo} analysis.
        @return:
        @rtype:
        """

        super(DESeq, self).run()

        # Tuxedo requires a genome version.

        if not self.genome_version:
            raise Exception('A Tuxedo analysis requires a genome_version configuration option.')

        # For DESeq, all samples need adding to the Analysis regardless.
        # TODO: Sort out the comparison issue.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        prefix = 'rnaseq_deseq_global'

        comparison_directory = os.path.join(self.genome_directory, prefix)

        if not os.path.isdir(comparison_directory):
            try:
                os.makedirs(comparison_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        annotation_sheet = AnnotationSheet(
            file_path=os.path.join(comparison_directory, prefix + '_comparisons.tsv'),
            file_type='excel-tab',
            name='DESeq Annotation',
            header=True)

        # Re-index the sample group dict by sample name.
        sample_dict = dict()
        """ @type sample_dict: dict[str, list[str]] """
        for group_name, sample_list in self.collection.sample_group_dict.iteritems():
            for sample in sample_list:
                if sample.name not in sample_dict:
                    sample_dict[sample.name] = list()
                sample_dict[sample.name].append(group_name)

        for sample in self.sample_list:
            if self.debug > 0:
                print repr(self) + ' Sample name: ' + sample.name
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            row_dict = dict()
            """ @type row_dict: dict[str, str | unicode] """
            # Set the group from the original Sample Annotation Sheet.
            row_dict['sample'] = sample.name
            row_dict['group'] = sample_dict[sample.name][0]
            row_dict['bam_path'] = StarAligner.get_prefix_star_aligner_merge(sample_name=sample.name) + '.bam'
            row_dict['bai_path'] = StarAligner.get_prefix_star_aligner_merge(sample_name=sample.name) + '.bai'
            # Set additional columns from the Sample Annotation Sheet prefixed with "Sample DESeq ".
            for annotation_key in filter(lambda x: x.startswith('DESeq '), sample.annotation_dict.keys()):
                row_dict[annotation_key[6:]] = sample.annotation_dict[annotation_key][0]

                # for paired_reads_name in paired_reads_name_list:
                # Create a ??? Runnable per paired_reads_name.

                # prefix_run_tophat = '_'.join((stage_run_tophat.name, paired_reads_name))

                # file_path_deseq = FilePathDESeq(prefix=prefix)

            annotation_sheet.row_dicts.append(row_dict)
        # TODO: Create a directory for the comparison.
        # TODO: Write a data frame defining the experiment for DESeq2.

        # Process all row_dict objects to get the superset of field names.
        for row_dict in annotation_sheet.row_dicts:
            for key in row_dict.iterkeys():
                if key not in annotation_sheet.field_names:
                    annotation_sheet.field_names.append(key)

        # Write the Annotation Sheet to disk.
        annotation_sheet.to_file_path()

        return
