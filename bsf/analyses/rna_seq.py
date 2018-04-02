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
import pickle
import re
import warnings

from bsf import Analysis, FilePath, Runnable
from bsf.analyses.star_aligner import StarAligner
from bsf.annotation import AnnotationSheet
from bsf.executables import TopHat
from bsf.ngs import SampleGroup
from bsf.process import Executable, RunnableStep, RunnableStepCopy, RunnableStepLink, RunnableStepMakeDirectory
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
        self.abundances_tsv = prefix + '_abundances.tsv'

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
        self.abundances_tsv = prefix + '_abundances.tsv'
        self.alignments_tsv = prefix + '_alignments.tsv'

        return


class TuxedoSamplePairSheet(AnnotationSheet):
    """The C{bsf.analyses.rna_seq.TuxedoSamplePairSheet} class represents C{bsf.ngs.Sample} pairs.

    The C{bsf.ngs.Sample} pairs are defined by the C{bsf_rnaseq_process_cuffdiff.R} script.

    Attributes:
    """

    _file_type = 'excel-tab'

    _field_names = [
        'V1',
        'V2',
    ]

    _test_methods = dict()


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
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode
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

        @param paired_reads_name: PairedReads name
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_tophat, paired_reads_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffmerge(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffmerge, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffquant(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffquant, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffnorm(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffnorm, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffdiff(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffdiff, comparison_name))

    @classmethod
    def get_prefix_rnaseq_process_cuffdiff(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison name
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
            sample_list=None,
            replicate_grouping=False,
            comparison_path=None,
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
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode
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
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if comparison_path is None:
            self.comparison_path = str()
        else:
            self.comparison_path = comparison_path

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

        self._comparison_dict = dict()
        """ @type _comparison_dict: dict[str, list[bsf.ngs.SampleGroup]] """

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
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

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

    def run(self):
        """Run this C{bsf.analyses.rna_seq.Tuxedo} analysis.
        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

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
            @return:
            @rtype:
            """
            if self.comparison_path:
                # A comparison file path was provided.
                if self.comparison_path == '*groups*':
                    # The special file name *groups* creates TuxedoComparison objects on the basis of an
                    # all-against-all group comparison.
                    # Without a comparison file path, simply add all Sample objects from the Collection.
                    self.sample_list.extend(self.collection.get_all_samples(exclude=True))

                    # Create a global comparison by adding all sample groups.
                    _sample_group_list = list()
                    """ @type _sample_group_list: list[bsf.ngs.SampleGroup] """

                    for _group_name, _sample_list in self.collection.sample_group_dict.iteritems():
                        _sample_group = SampleGroup(name=_group_name, sample_list=_sample_list)
                        # SampleGroup objects are only useful, if at least one Sample object is not excluded.
                        if not _sample_group.is_excluded():
                            _sample_group_list.append(_sample_group)

                    # Sort the list of comparison groups by group name.
                    _sample_group_list.sort(cmp=lambda x, y: cmp(x.name, y.name))
                    # Set the comparison name to 'global'.
                    self._comparison_dict['global'] = _sample_group_list
                elif self.comparison_path == '*samples*':
                    # The special file name *samples* creates TuxedoComparison objects on the basis of an
                    # all-against-all sample comparison.
                    # Without a comparison file path, simply add all Sample objects from the Collection.
                    self.sample_list.extend(self.collection.get_all_samples(exclude=True))

                    # Create a global comparison by adding all samples under their sample name as group name.
                    _sample_group_list = list()
                    """ @type _sample_group_list: list[bsf.ngs.SampleGroup] """

                    for _sample in self.sample_list:
                        # Sample objects are only useful, if at least one PairedReads object is not excluded.
                        if not _sample.is_excluded():
                            _sample_group_list.append(SampleGroup(name=_sample.name, sample_list=[_sample]))

                    # Sort the list of comparison groups by group name.
                    _sample_group_list.sort(cmp=lambda x, y: cmp(x.name, y.name))
                    # Set the comparison name to 'global'.
                    self._comparison_dict['global'] = _sample_group_list
                else:
                    # A comparison file path was provided.
                    self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)
                    # Read and process the comparison file, which includes adding only those Sample objects,
                    # which are referenced in a comparison.
                    annotation_sheet = AnnotationSheet.from_file_path(file_path=self.comparison_path)
                    regular_expression = re.compile(pattern='\W')

                    for row_dict in annotation_sheet.row_dicts:
                        _sample_group_list = list()
                        """ @type _sample_group_list: list[bsf.ngs.SampleGroup] """
                        _comparison_name_list = list()
                        """ @type _comparison_name_list: list[str] """
                        # In addition to defining samples, allow also the definition of groups in comparison files.
                        # If the row dictionary has a 'Group' key, then the Sample in the same row gets added to
                        # the group. So,
                        # 'ProcessedRunFolder', 'Project', 'Sample', 'Group' defines the groups, while ...
                        # 'Control Group','Treatment Group' defines a comparison, as does ...
                        # 'Control Group','Treatment ProcessedRunFolder','Treatment Project','Treatment Sample'

                        # Get Sample objects for classical 'Control' and 'Treatment' keys,
                        # before looking up a series of 'Point N' keys.
                        i = -2
                        while True:
                            i += 1
                            if i == -1:
                                prefix = 'Control'
                            elif i == 0:
                                prefix = 'Treatment'
                            else:
                                prefix = 'Point ' + str(i)
                            # Get Sample objects for 'Point N' keys for as long as they are defined.
                            # The Collection.get_sample_from_row_dict method can return one or more Sample objects,
                            # depending on 'Group' or 'Sample' column entries.
                            # In RNA-Seq experiments, entire pools of Sample objects (replicates) are compared
                            # with each other.
                            _group_name, _sample_list_old = self.collection.get_samples_from_row_dict(
                                row_dict=row_dict,
                                prefix=prefix)
                            _sample_list_new = list()
                            """ @type _sample_list_new: list[bsf.ngs.Sample] """
                            if _group_name and len(_sample_list_old):
                                # Sample objects are unly useful, if at least one PairedReads object is not excluded.
                                for _sample in _sample_list_old:
                                    if not _sample.is_excluded():
                                        _sample_list_new.append(_sample)

                                if len(_sample_list_new):
                                    _comparison_name_list.append(_group_name)
                                    _sample_group_list.append(SampleGroup(
                                        name=_group_name,
                                        sample_list=_sample_list_new))
                                    # Also expand each Python list of bsf.ngs.Sample objects to get all those
                                    # bsf.ngs.Sample objects that this bsf.Analysis needs considering.
                                    for _sample in _sample_list_new:
                                        self.add_sample(sample=_sample)
                                        if self.debug > 1:
                                            print '  ', prefix, 'Sample name:', _sample.name, \
                                                'file_path:', _sample.file_path
                                        if self.debug > 2:
                                            print sample.trace(1)
                            elif i < 1:
                                # A Control and Treatment prefix is not required.
                                continue
                            else:
                                # Only break if there is no further 'Point N' prefix.
                                break

                        if 'Comparison Name' in row_dict and row_dict['Comparison Name']:
                            # For ridiculously large comparisons involving loads of groups or samples a
                            # comparison name can be explicitly specified.
                            # Any non-word characters get replaced by underscore characters.
                            _comparison_name = re.sub(pattern=regular_expression, repl='_',
                                                      string=row_dict['Comparison Name'])
                        else:
                            _comparison_name = '__'.join(_comparison_name_list)

                        # Sort the list of comparison groups by group name.
                        _sample_group_list.sort(cmp=lambda x, y: cmp(x.name, y.name))
                        # Set the comparison name.
                        self._comparison_dict[_comparison_name] = _sample_group_list
            else:
                # Without a comparison file path, simply add a comparison 'global' with a SampleGroup 'global'
                # with all Sample objects from the Collection.
                # This means that most pipeline stages with the exception of Cuffdiff can run.
                self.sample_list.extend(self.collection.get_all_samples(exclude=True))
                self._comparison_dict['global'] = [SampleGroup(
                    name='global',
                    sample_list=self.collection.get_all_samples(exclude=True))]

            # Check for comparisons without SampleGroup objects or SampleGroup objects without Sample.
            # FIXME: The complication is that Sample or ReadGroup objects could be excluded from the Analysis.
            for _comparison_name, _sample_group_list in self._comparison_dict.iteritems():
                if self.debug > 0:
                    print 'Comparison name:', _comparison_name
                    print 'SampleGroup list:'
                if len(_sample_group_list) < 1:
                    warnings.warn('Comparison ' + _comparison_name + ' without SampleGroup objects', UserWarning)
                for _sample_group in _sample_group_list:
                    if self.debug > 0:
                        print '  SampleGroup name:', _sample_group.name
                        print '  SampleGroup Sample list:'
                    if len(_sample_group.sample_list) < 1:
                        warnings.warn('SampleGroup ' + _sample_group.name + ' without Sample objects')
                    for _sample in _sample_group.sample_list:
                        if self.debug > 0:
                            print '    Sample name:', _sample.name

            return

        def run_write_annotation(annotation_path, annotation_dict):
            """Private function to write a sample annotation file for Cuffdiff or Cuffnorm to disk.

            @param annotation_path: Annotation file path
            @type annotation_path: str | unicode
            @param annotation_dict: Annotation dict
            @type annotation_dict: dict[str, list[str | unicode]]
            """
            _annotation_file = open(annotation_path, 'w')
            _annotation_file.write('sample_id\tgroup_label\n')
            for _group_name, per_group_list in annotation_dict.iteritems():
                for _file_path in per_group_list:
                    _annotation_file.write(_file_path + '\t' + _group_name + '\n')
            _annotation_file.close()

            return

        # Start of the run() method body.

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

        run_read_comparisons()

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

        genome_sizes_path = self.genome_fasta_path + '.fai'

        if not os.path.exists(self.genome_fasta_path):
            raise Exception('Genome FASTA file path {!r} does not exists.'.format(self.genome_fasta_path))

        # Define a reference transcriptome index directory or a GTF file path.

        if self.transcriptome_index_path:
            # Check if the transcriptome_index_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_index_path = self.configuration.get_absolute_path(
                file_path=self.transcriptome_index_path,
                default_path=Default.absolute_transcriptome_resource(transcriptome_version=''))

            if not os.path.isdir(self.transcriptome_index_path):
                raise Exception('Reference transcriptome index directory {!r} does not exist.'.
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
                raise Exception('Reference transcriptome GTF file {!r} does not exist.'.
                                format(self.transcriptome_gtf_path))
        elif self.transcriptome_gtf_path:
            # Check if transcriptome_gtf_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_gtf_path = self.configuration.get_absolute_path(
                file_path=self.transcriptome_gtf_path,
                default_path=Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception('Reference transcriptome GTF file {!r} does not exist.'.
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
                raise Exception('Reference transcriptome GTF file path {!r} does not exist.'.
                                format(self.transcriptome_gtf_path))

        if not self.transcriptome_gtf_path:
            raise Exception('Reference transcriptome GTF file not defined.')

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

        # Sort the Python list of Sample objects by the Sample.name.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        for sample in self.sample_list:
            if self.debug > 0:
                print self, 'Sample name:', sample.name
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
                # The TopHat coverage search finds additional 'GT-AG' introns, but is only recommended for
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
                        print self, 'PairedReads name:', paired_reads.get_name()

                    if paired_reads.reads_1:
                        reads_1_file_path_list.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads_2_file_path_list.append(paired_reads.reads_2.file_path)

                # Pass lists of files into Tophat, regardless of whether read 2 exists.
                # The bsf_run_rnaseq_tophat.py script eventually removes an empty read 2 argument.
                tophat.arguments.append(','.join(reads_1_file_path_list))
                tophat.arguments.append(','.join(reads_2_file_path_list))

                # Create a new rnaseq_run_tophat Executable.
                # TODO: The following code block is required as long as the bsf_run_rnaseq_tophat.py script
                # has not been retired.

                if self.debug > 0:
                    print 'Tophat Executable'
                    print tophat.trace(level=1)

                pickler_dict_run_tophat = {
                    'prefix': stage_run_tophat.name,
                    'replicate_key': paired_reads_name,
                    'tophat_executable': tophat,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    stage_run_tophat.name + '_' + paired_reads_name + '.pkl')
                pickler_file = open(pickler_path, 'wb')
                pickler = pickle.Pickler(file=pickler_file, protocol=pickle.HIGHEST_PROTOCOL)
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

                # Only submit this Executable if the 'align_summary.txt' file does not exist.
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

                # Only submit this Executable if the 'accepted_hits.bam.bai' file does not exist.
                file_path_temporary = os.path.join(
                    self.genome_directory,
                    file_path_tophat.accepted_hits_bai)
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
                # TODO: The location of the autoSQL file needs to be configurable.
                runnable_step.add_option_pair_short(key='as', value='/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
                runnable_step.add_switch_short(key='tab')
                runnable_step.add_option_pair_short(key='type', value='bed12+8')
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
            # TODO: Convert this into a Runnable and RunnableStep.
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
        executable_cuffmerge_dict = dict()
        """ @type executable_cuffmerge_dict: dict[str, bsf.process.Executable] """

        comparison_name_list = self._comparison_dict.keys()
        comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

        if self.debug > 0:
            print 'Tuxedo Comparison dict key list', comparison_name_list

        for comparison_name in comparison_name_list:
            if self.debug > 0:
                print '  Comparison name:', comparison_name

            sample_group_list = self._comparison_dict[comparison_name]

            if len(sample_group_list) == 0:
                continue

            # Process rnaseq_cuffmerge and rnaseq_cuffdiff arguments in parallel.
            # Check that the comparison contains at least one sample group.

            cuffdiff_cuffnorm_abundances_dict = dict()
            """ @type cuffdiff_cuffnorm_abundances_dict: dict[str, list[str | unicode]] """
            cuffdiff_cuffnorm_alignments_dict = dict()
            """ @type cuffdiff_cuffnorm_alignments_dict: dict[str, list[str | unicode]] """
            cuffdiff_cuffnorm_dependencies = list()
            """ @type cuffdiff_cuffnorm_dependencies: list[str] """
            cuffmerge_cuffnorm_submit = len(sample_group_list) >= 1
            """ @type cuffmerge_cuffnorm_submit: bool """
            cuffmerge_transcript_gtf_list = list()
            """ @type cuffmerge_transcript_gtf_list: list[str | unicode] """

            # TODO: Should the comparison prefix also include the project name or number?
            prefix_cuffmerge = self.get_prefix_rnaseq_run_cuffmerge(comparison_name=comparison_name)

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
            # Submit the Executable if the status file AND the sample group list above supports it.
            executable_run_cuffmerge.submit &= cuffmerge_cuffnorm_submit
            # Set a dependency on all other Cuffmerge process to avoid file contention.
            executable_cuffmerge_dict[prefix_cuffmerge] = executable_run_cuffmerge

            if self.novel_transcripts:
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

                # Add the assembly manifest file as Cuffmerge argument.
                # The file will be written below.
                runnable_step_cuffmerge.arguments.append(file_path_cuffmerge.assembly_txt)
            else:
                # If novel transcripts are not assembled, create RunnableStep objects to create the output directory
                # and copy the reference transcriptome GTF file.

                runnable_run_cuffmerge.add_runnable_step(
                    runnable_step=RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_cuffmerge.output_directory))

                runnable_run_cuffmerge.add_runnable_step(
                    runnable_step=RunnableStepCopy(
                        name='copy',
                        source_path=self.transcriptome_gtf_path,
                        target_path=file_path_cuffmerge.merged_gtf))

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
            # TODO: The location of the autoSQL file needs to be configurable.
            runnable_step.add_option_pair_short(key='as', value='/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
            runnable_step.add_switch_short(key='tab')
            runnable_step.add_option_pair_short(key='type', value='bed12+8')
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

            for sample_group in sample_group_list:
                if self.debug > 0:
                    print '    SampleGroup name:', sample_group.name

                per_group_abundances_list = list()
                """ @type per_group_abundances_list: list[str | unicode] """
                per_group_alignments_list = list()
                """ @type per_group_alignments_list: list[str | unicode] """

                for sample in sample_group.sample_list:
                    if self.debug > 0:
                        print '      Sample name:', sample.name

                    paired_reads_dict = sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping,
                        exclude=True)

                    paired_reads_name_list = paired_reads_dict.keys()
                    paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                    for paired_reads_name in paired_reads_name_list:
                        if self.debug > 0:
                            print '        PairedReads name:', paired_reads_name
                        # Add the Cufflinks assembled transcripts GTF to the Cuffmerge manifest.
                        cuffmerge_transcript_gtf_list.append(
                            os.path.join('_'.join(('rnaseq_cufflinks', paired_reads_name)), 'transcripts.gtf') + '\n')

                        # Wait for each Cufflinks PairedReads to finish, before Cuffmerge can run.

                        executable_run_cuffmerge.dependencies.append(
                            '_'.join((self.stage_name_run_cufflinks, paired_reads_name)))

                        # Create a Cuffquant Runnable per comparison (comparison_name) and
                        # PairedReads (paired_reads_name) on the basis of the Cuffmerge GTF file.

                        prefix_cuffquant = '_'.join((stage_run_cuffquant.name, comparison_name, paired_reads_name))

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

                cuffdiff_cuffnorm_abundances_dict[sample_group.name] = per_group_abundances_list
                cuffdiff_cuffnorm_alignments_dict[sample_group.name] = per_group_alignments_list

            if self.novel_transcripts:
                # Write a Cuffmerge assembly manifest file to merge all transcriptome GTF files of each Sample object.
                # This requires an absolute path, because the working directory is not set at the stage of
                # job submission.
                assembly_file = open(os.path.join(self.genome_directory, file_path_cuffmerge.assembly_txt), 'w')
                assembly_file.writelines(cuffmerge_transcript_gtf_list)
                assembly_file.close()

            # Create a Cuffnorm Runnable per comparison.

            prefix_run_cuffnorm = '_'.join((stage_run_cuffnorm.name, comparison_name))

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
            # Submit the Executable if the status file AND the sample group list above supports it.
            executable_run_cuffnorm.submit &= cuffmerge_cuffnorm_submit
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
            runnable_step_cuffnorm.add_switch_long(
                key='use-sample-sheet')

            # Add the Cuffmerge GTF file as first Cuffnorm argument.
            runnable_step_cuffnorm.arguments.append(file_path_cuffmerge.merged_gtf)

            # Add an abundances annotation file as second Cuffnorm argument.
            # Writing a Cuffnorm abundances TSV file requires an absolute path,
            # because the working directory is not set at the current stage of job submission.
            run_write_annotation(
                annotation_path=os.path.join(self.genome_directory, file_path_run_cuffnorm.abundances_tsv),
                annotation_dict=cuffdiff_cuffnorm_abundances_dict)
            runnable_step_cuffnorm.arguments.append(file_path_run_cuffnorm.abundances_tsv)

            if len(self._comparison_dict[comparison_name]) >= 2:
                # Create a Cuffdiff Runnable per comparison if there are at least two SampleGroup objects.

                prefix_run_cuffdiff = '_'.join((stage_run_cuffdiff.name, comparison_name))

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
                runnable_step_cuffdiff.add_switch_long(
                    key='use-sample-sheet')

                # Add the Cuffmerge GTF file as first Cuffdiff argument.
                runnable_step_cuffdiff.arguments.append(file_path_cuffmerge.merged_gtf)

                # Add an abundances or alignment annotation file as second Cuffdiff argument.
                # Cuffdiff seems to have a problem with Cuffquant abundances files in that the isoforms.count_tracking
                # files show ridiculously low numbers such as 1e-319 for some splice variants. Usually, other splice
                # variants in the same cluster seem fine.
                if run_cuffquant_before_cuffdiff:
                    # Writing a Cuffdiff abundances TSV file requires an absolute path,
                    # because the working directory is not set at the current stage of job submission.
                    run_write_annotation(
                        annotation_path=os.path.join(self.genome_directory, file_path_run_cuffdiff.abundances_tsv),
                        annotation_dict=cuffdiff_cuffnorm_abundances_dict)
                    runnable_step_cuffdiff.arguments.append(file_path_run_cuffdiff.abundances_tsv)
                else:
                    # Writing a Cuffdiff alignments TSV file requires an absolute path,
                    # because the working directory is not set at the current stage of job submission.
                    run_write_annotation(
                        annotation_path=os.path.join(self.genome_directory, file_path_run_cuffdiff.alignments_tsv),
                        annotation_dict=cuffdiff_cuffnorm_alignments_dict)
                    runnable_step_cuffdiff.arguments.append(file_path_run_cuffdiff.alignments_tsv)

                # Create a new rnaseq_process_cuffdiff Executable.

                executable_process_cuffdiff = stage_process_cuffdiff.add_executable(
                    executable=Executable(
                        name='_'.join((stage_process_cuffdiff.name, comparison_name)),
                        program='bsf_rnaseq_process_cuffdiff.R'))
                executable_process_cuffdiff.dependencies.append(executable_run_cuffdiff.name)

                # Set rnaseq_process_cuffdiff options.
                self.set_command_configuration(command=executable_process_cuffdiff)
                executable_process_cuffdiff.add_option_long(
                    key='comparison-name',
                    value=comparison_name)
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

        # Finally, set dependencies on all other Cuffmerge Executable objects to avoid file contention.
        for prefix_cuffmerge in executable_cuffmerge_dict.keys():
            for executable_cuffmerge in executable_cuffmerge_dict.values():
                if prefix_cuffmerge != executable_cuffmerge.name:
                    executable_cuffmerge.dependencies.append(prefix_cuffmerge)

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

            str_list += '<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n'
            str_list += '\n'

            # TopHat and Cufflinks table.

            str_list += '<h2 id="transcriptome_browsing">Transcriptome Browsing</h2>\n'
            str_list += '\n'

            str_list += '<h3 id="read_alignments">Read Alignments</h3>\n'
            str_list += '\n'

            str_list += '<p id ="tophat">'
            # http://tophat.cbcb.umd.edu/manual.html
            str_list += '<a href="http://ccb.jhu.edu/software/tophat/index.shtml"><strong>TopHat</strong></a> '
            str_list += 'aligns RNA-Seq reads to a genome in order to identify '
            str_list += 'exon-exon splice junctions. It is built on the ultra fast '
            str_list += 'short read mapping program '
            str_list += '<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">'
            str_list += '<strong>Bowtie 2</strong>'
            str_list += '</a>.'
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<p id="track_hub">'
            str_list += 'View TopHat <strong>read alignments</strong> tracks for each sample\n'
            str_list += 'in their genomic context via the project-specific '
            str_list += self.ucsc_hub_html_anchor(link_path=link_path)
            str_list += '.'
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<p>'
            str_list += '<a href="rnaseq_tophat_alignment_summary.pdf">'
            str_list += '<img '
            str_list += 'alt="TopHat Alignment Summary" '
            str_list += 'id="tophat_alignment_summary_img" '
            str_list += 'src="rnaseq_tophat_alignment_summary.png" '
            str_list += 'height="80" '
            str_list += 'width="80" '
            str_list += '/>'
            str_list += '</a> '
            str_list += 'Alignment summary statistics <a href="rnaseq_tophat_alignment_summary.tsv">TSV</a>'
            str_list += '</p>\n'

            str_list += '<h3 id="alignment_events">Splice Junctions, Insertions and Deletions</h3>\n'
            str_list += '\n'

            str_list += '<p>'
            str_list += 'TopHat reports <strong>splice junctions</strong> on the basis of RNA-Seq\n'
            str_list += 'read alignments in UCSC BED track format.\n'
            str_list += 'Each junction consists of two connected BED blocks,\n'
            str_list += 'where each block is as long as the maximal overhang\n'
            str_list += 'of any read spanning the junction. The score is\n'
            str_list += 'the number of alignments spanning the junction.\n'
            str_list += 'UCSC BED tracks of <strong>insertions</strong> and\n'
            str_list += '<strong>deletions</strong> are also reported by TopHat.'
            str_list += '</p>\n'

            str_list += '<p>'
            str_list += 'View the corresponding TopHat tracks for junctions, deletions and insertions\n'
            str_list += 'for each sample in their genomic context via the project-specific\n'
            str_list += self.ucsc_hub_html_anchor(link_path=link_path)
            str_list += '.'
            str_list += '</p>\n'
            str_list += '\n'

            # str_list += '<p>'
            # str_list += 'Follow the links below to attach\n'
            # str_list += 'Tophat junction, deletion and insertion annotation to the\n'
            # str_list += 'UCSC Genome Browser. Since each file needs transferring to\n'
            # str_list += 'the UCSC site, subsequent pages will take some time to load.'
            # str_list += '</p>\n'

            str_list += '<h2 id="gene_expression_profiles">Gene Expression Profiles</h2>\n'
            str_list += '\n'

            str_list += '<p id="cufflinks">'
            # http://cufflinks.cbcb.umd.edu/howitworks.html
            str_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/"><strong>Cufflinks</strong></a>\n'
            str_list += 'assembles aligned RNA-Seq reads into transcripts,\n'
            str_list += 'estimates their abundances, and tests for differential\n'
            str_list += 'expression and regulation transcriptome-wide.\n'
            str_list += 'It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of\n'
            str_list += 'transcripts. Cufflinks then estimates the relative abundances of these transcripts based\n'
            str_list += 'on how many reads support each one, taking into account biases in library preparation '
            str_list += 'protocols.'
            str_list += '</p>\n'

            str_list += '<p>'
            str_list += 'The Cufflinks <strong>assembled transcripts</strong> can be attached to the \n'
            str_list += 'UCSC Genome Browser, by following the "Transcript Assembly" links\n'
            str_list += 'below.\n'
            str_list += 'The isoforms.fpkm_tracking and genes.fpkm_tracking files\n'
            str_list += 'contain the estimated isoform or gene expression values in the generic\n'
            # http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
            str_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                        'fpkm-tracking-format">FPKM Tracking format</a>.\n'
            str_list += 'The isoforms.count_tracking and genes.count_tracking files\n'
            str_list += 'contain the scaled isoform or gene count values in the generic\n'
            str_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                        'count-tracking-format">Count Tracking format</a>.'
            str_list += '</p>\n'

            str_list += '<p>'
            str_list += 'Please see a more detailed description of\n'
            # http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output
            str_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                        'output-formats-used-in-the-cufflinks-suite">Cufflinks output</a>.'
            str_list += '</p>\n'

            str_list += '<table id="gene_expression_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Sample</th>\n'
            str_list += '<th>Assembled Transcripts</th>\n'
            str_list += '<th>Gene FPKM</th>\n'
            str_list += '<th>Transcript FPKM</th>\n'
            str_list += '<th>Genes (Symbols)</th>\n'
            str_list += '<th>Isoforms (Symbols)</th>\n'
            str_list += '<th>Aligned BAM file</th>\n'
            str_list += '<th>Aligned BAI file</th>\n'
            str_list += '<th>Unaligned BAM file</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            for sample in self.sample_list:
                if self.debug > 0:
                    print self, 'Sample name:', sample.name
                    print sample.trace(1)

                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                paired_reads_name_list = paired_reads_dict.keys()
                if not len(paired_reads_name_list):
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue
                paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for paired_reads_name in paired_reads_name_list:
                    # Cufflinks produces genes.fpkm_tracking, isoforms.fpkm_tracking,
                    # skipped.gtf and transcripts.gtf.

                    path_prefix = 'rnaseq_cufflinks_' + paired_reads_name

                    str_list += '<tr>\n'
                    # Sample
                    str_list += '<td class="left">' + paired_reads_name + '</td>\n'
                    # Assembled Transcripts
                    str_list += '<td class="center">'
                    str_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='transcripts.gtf',
                        text='Transcript Assembly')
                    str_list += '</td>\n'
                    # Gene FPKM
                    str_list += '<td class="center">'
                    str_list += '<a href="' + path_prefix + '/genes.fpkm_tracking">Genes FPKM</a>'
                    str_list += '</td>\n'
                    # Transcript FPKM
                    str_list += '<td class="center">'
                    str_list += '<a href="' + path_prefix + '/isoforms.fpkm_tracking">Isoforms FPKM</a>'
                    str_list += '</td>\n'
                    # Genes (Symbols)
                    str_list += '<td class="center">'
                    str_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='genes_fpkm_tracking.tsv',
                        text='Genes (Symbols)')
                    str_list += '</td>\n'
                    # Isoforms (Symbols)
                    str_list += '<td class="center">'
                    str_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='isoforms_fpkm_tracking.tsv',
                        text='Isoforms (Symbols)')

                    # TODO: The aligned BAM and BAI files and the unaligned BAM file are currently non standard.
                    # The files have a 'rnaseq_tophat_' prefix, but are in the 'rnaseq_cufflinks_' directory.
                    # This will be resolved when the process_tophat step gets re-engineered.
                    # Aligned BAM file
                    str_list += '</td>\n'
                    str_list += '<td class="center">'
                    str_list += '<a href="' + path_prefix + '/rnaseq_tophat_' + paired_reads_name + \
                                '_accepted_hits.bam">Aligned BAM</a>'
                    str_list += '</td>\n'
                    # Aligned BAI file
                    str_list += '<td class="center">'
                    str_list += '<a href="' + path_prefix + '/rnaseq_tophat_' + paired_reads_name + \
                                '_accepted_hits.bam.bai">Aligned BAI</a>'
                    str_list += '</td>\n'
                    # Unaligned BAM file
                    str_list += '<td class="center">'
                    str_list += '<a href="' + path_prefix + '/rnaseq_tophat_' + paired_reads_name + \
                                '_unmapped.bam">Unaligned BAM</a>'
                    str_list += '</td>\n'
                    str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            # Cuffdiff produces cds_exp.diff, gene_exp.diff, isoform_exp.diff
            # promoters.diff, splicing.diff and tss_group_exp.diff amongst many others.

            str_list += '<h2 id="differential_expression">Differential Expression</h2>\n'
            str_list += '\n'

            str_list += '<p id="cuffdiff">'
            # http://cufflinks.cbcb.umd.edu/howitworks.html#diff
            str_list += '<a href="http://cole-trapnell-lab.github.io/cufflinks/"><strong>Cufflinks</strong></a>\n'
            str_list += 'finds significant changes in transcript\n'
            str_list += 'expression, splicing, and promoter use.'
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<h3 id="all_genes">All Genes</h3>\n'

            str_list += '<table id="differential_expression_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Comparison</th>\n'
            str_list += '<th>Samples</th>\n'
            str_list += '<th>Replicates</th>\n'
            str_list += '<th>Coding Sequences</th>\n'
            str_list += '<th>Genes</th>\n'
            str_list += '<th>Isoforms</th>\n'
            str_list += '<th>Promoters</th>\n'
            str_list += '<th>Splicing</th>\n'
            str_list += '<th>Transcription Start Sites</th>\n'
            str_list += '<th>Gene FPKM Replicates</th>\n'
            str_list += '<th>Gene Count Replicates</th>\n'
            str_list += '<th>Isoform FPKM Replicates</th>\n'
            str_list += '<th>Isoform Count Replicates</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for comparison_name in comparison_name_list:
                path_prefix = 'rnaseq_process_cuffdiff_' + comparison_name

                # Link to comparison-specific symbolic links in the directory after cummeRbund processing.

                str_list += '<tr>\n'
                # Comparison
                str_list += '<td class="left">' + comparison_name + '</td>\n'
                # Samples
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='samples.tsv',
                    text='Samples')
                str_list += '</td>\n'
                # Replicates
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='replicates.tsv',
                    text='Replicates')
                str_list += '</td>\n'
                # Coding Sequences
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='cds_exp_diff.tsv',
                    text='Coding Sequences')
                str_list += '</td>\n'
                # Genes
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_exp_diff.tsv',
                    text='<strong>Genes</strong>')
                str_list += '</td>\n'
                # Isoforms
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_exp_diff.tsv',
                    text='Isoforms')
                str_list += '</td>\n'
                # Promoters
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='promoters_diff.tsv',
                    text='Promoters')
                str_list += '</td>\n'
                # Splicing
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='splicing_diff.tsv',
                    text='Splicing')
                str_list += '</td>\n'
                # Transcription Start Sites
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='tss_group_exp_diff.tsv',
                    text='Transcription Start Sites')
                str_list += '</td>\n'
                # Gene FPKM Replicates
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_fpkm_replicates.tsv',
                    text='Gene FPKM Replicates')
                str_list += '</td>\n'
                # Gene Count Replicates
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_counts_replicates.tsv',
                    text='Gene Count Replicates')
                str_list += '</td>\n'
                # Isoform FPKM Replicates
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_fpkm_replicates.tsv',
                    text='Isoform FPKM Replicates')
                str_list += '</td>\n'
                # Isoform Count Replicates
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_counts_replicates.tsv',
                    text='Isoform Count Replicates')
                str_list += '</td>\n'
                str_list += '</tr>\n'

                # Read sample pair information if available.

                sample_pair_path = os.path.join(
                    self.genome_directory,
                    path_prefix,
                    '_'.join((path_prefix, 'sample_pairs.tsv')))

                if os.path.exists(sample_pair_path):
                    sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                    for row_dict in sample_pair_sheet.row_dicts:
                        str_list += '<tr>\n'
                        # Comparison
                        str_list += '<td></td>\n'
                        # Samples
                        # Replicates
                        # Coding Sequences
                        str_list += '<td class="left" colspan="3">'
                        str_list += '<strong>' + row_dict['V1'] + '</strong>'
                        str_list += ' versus '
                        str_list += '<strong>' + row_dict['V2'] + '</strong>'
                        str_list += '</td>\n'
                        # Genes
                        str_list += '<td class="center">'
                        str_list += relative_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_diff.tsv')),
                            text='<strong>Genes</strong>')
                        str_list += '</td>\n'
                        # Isoforms
                        str_list += '<td class="center">'
                        str_list += relative_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'isoforms_diff.tsv')),
                            text='Isoforms')
                        str_list += '</td>\n'
                        # Promoters
                        # Splicing
                        # Transcription Start Sites
                        # Gene FPKM Replicates
                        # Gene Count Replicates
                        # Isoform FPKM Replicates
                        # Isoform Count Replicates
                        str_list += '<td class="left" colspan="7"></td>\n'
                        str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            str_list += '<h3 id="significant_genes">Significant Genes</h3>\n'

            str_list += '<table id="significant_genes_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Comparison</th>\n'
            str_list += '<th>Genes</th>\n'
            str_list += '<th>Isoforms</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            for comparison_name in comparison_name_list:
                path_prefix = 'rnaseq_process_cuffdiff_' + comparison_name

                str_list += '<tr>\n'
                # Comparison
                str_list += '<td class="left">' + comparison_name + '</td>\n'
                # Genes
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_significance_matrix.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_significance_matrix.png',
                        text='Significance Matrix Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'
                # Isoforms
                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_significance_matrix.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_significance_matrix.png',
                        text='Significance Matrix Plot - Isoforms - ' + comparison_name))
                str_list += '</td>\n'
                str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'

            # Show cummeRbund quality plots.

            str_list += '<h2 id="quality_plots">Quality Plots</h2>\n'
            str_list += '\n'

            str_list += '<p>'
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<table id="quality_plots_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Comparison</th>\n'
            str_list += '<th>Dispersion Plot - Genes</th>\n'
            str_list += '<th>Dispersion Plot - Isoforms</th>\n'
            str_list += '<th>Squared Coefficient of Variation - Genes</th>\n'
            str_list += '<th>Squared Coefficient of Variation - Isoforms</th>\n'
            str_list += '<th>Density Plot without Replicates - Genes</th>\n'
            str_list += '<th>Density Plot with Replicates - Genes</th>\n'
            str_list += '<th>Density Plot without Replicates - Isoforms</th>\n'
            str_list += '<th>Density Plot with Replicates - Isoforms</th>\n'
            str_list += '<th>Box Plot without Replicates - Genes</th>\n'
            str_list += '<th>Box Plot with Replicates - Genes</th>\n'
            str_list += '<th>Box Plot without Replicates - Isoforms</th>\n'
            str_list += '<th>Box Plot with Replicates - Isoforms</th>\n'
            str_list += '<th>Scatter Matrix Plot - Genes</th>\n'
            str_list += '<th>Scatter Matrix Plot - Isoforms</th>\n'
            str_list += '<th>Dendrogram Plot</th>\n'
            str_list += '<th>Volcano Matrix Plot - Genes</th>\n'
            str_list += '<th>Multidimensional Scaling Plot - Genes</th>\n'
            str_list += '<th>Principal Component Analysis Plot - Genes</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            for comparison_name in comparison_name_list:
                path_prefix = 'rnaseq_process_cuffdiff_' + comparison_name

                str_list += '<tr>\n'
                # Comparison
                str_list += '<td class="left">' + comparison_name + '</td>\n'

                # Dispersion Plots for Genes and Isoforms

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_dispersion.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_dispersion.png',
                        text='Dispersion Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_dispersion.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_dispersion.png',
                        text='Dispersion Plot - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                # Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms

                str_list += '<td class="center">'
                if os.path.exists(
                        path=os.path.join(self.genome_directory, path_prefix, path_prefix + '_genes_scv.png')):
                    str_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='genes_scv.pdf',
                        text=relative_image(
                            prefix=path_prefix,
                            suffix='genes_scv.png',
                            text='Squared Coefficient of Variation (SCV) - Genes - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                if os.path.exists(
                        path=os.path.join(self.genome_directory, path_prefix, path_prefix + '_isoforms_scv.png')):
                    str_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='isoforms_scv.pdf',
                        text=relative_image(
                            prefix=path_prefix,
                            suffix='isoforms_scv.png',
                            text='Squared Coefficient of Variation (SCV) - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                # Density Plots for Genes without and with Replicates

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_density_wo_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_density_wo_replicates.png',
                        text='Density Plot without Replicates - Genes - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_density_w_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_density_w_replicates.png',
                        text='Density Plot with Replicates - Genes - ' + comparison_name))
                str_list += '</td>\n'

                # Density Plots for Isoforms without and with Replicates

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_density_wo_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_density_wo_replicates.png',
                        text='Density Plot without Replicates - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_density_w_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_density_w_replicates.png',
                        text='Density Plot with Replicates - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                # Box Plots for Genes without and with Replicates

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_box_wo_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_box_wo_replicates.png',
                        text='Box Plot without Replicates - Genes - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_box_w_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_box_w_replicates.png',
                        text='Box Plot with Replicates - Genes - ' + comparison_name))
                str_list += '</td>\n'

                # Box Plots for Isoforms with and without Replicates

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_box_wo_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_box_wo_replicates.png',
                        text='Box Plot without Replicates - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_box_w_replicates.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_box_w_replicates.png',
                        text='Box Plot with Replicates - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                # Scatter Matrix Plot for Genes and Isoforms

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_scatter_matrix.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_scatter_matrix.png',
                        text='Scatter Matrix Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_scatter_matrix.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='isoforms_scatter_matrix.png',
                        text='Scatter Matrix Plot - Isoforms - ' + comparison_name))
                str_list += '</td>\n'

                # Dendrogram Plot for Genes

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_dendrogram.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_dendrogram.png',
                        text='Dendrogram Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'

                # Volcano Matrix Plot for Genes

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_volcano_matrix.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_volcano_matrix.png',
                        text='Volcano Matrix Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'

                # Multidimensional Scaling Plot for Genes

                str_list += '<td class="center">'
                if os.path.exists(
                        path=os.path.join(self.genome_directory, path_prefix, path_prefix + '_genes_mds.png')):
                    str_list += relative_anchor(
                        prefix=path_prefix,
                        suffix='genes_mds.pdf',
                        text=relative_image(
                            prefix=path_prefix,
                            suffix='genes_mds.png',
                            text='Multidimensional Scaling Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'

                # Principal Component Analysis Plot for Genes

                str_list += '<td class="center">'
                str_list += relative_anchor(
                    prefix=path_prefix,
                    suffix='genes_pca.pdf',
                    text=relative_image(
                        prefix=path_prefix,
                        suffix='genes_pca.png',
                        text='Principal Component Analysis Plot - Genes - ' + comparison_name))
                str_list += '</td>\n'

                str_list += '</tr>\n'

                # Read sample pair information if available.

                sample_pair_path = os.path.join(
                    self.genome_directory,
                    path_prefix,
                    '_'.join((path_prefix, 'sample_pairs.tsv')))

                if os.path.exists(sample_pair_path):
                    sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                    for row_dict in sample_pair_sheet.row_dicts:
                        str_list += '<tr>\n'

                        # Comparison
                        str_list += '<td class="left"></td>\n'
                        str_list += '<td class="left" colspan="10">'
                        str_list += '<strong>' + row_dict['V1'] + '</strong>'
                        str_list += ' versus '
                        str_list += '<strong>' + row_dict['V2'] + '</strong>'
                        str_list += '</td>\n'

                        str_list += '<td class="center">'
                        str_list += relative_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_scatter.pdf')),
                            text=relative_image(
                                prefix=path_prefix,
                                suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_scatter.png')),
                                text='Scatter Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2']))
                        str_list += '</td>\n'

                        str_list += '<td class="center"></td>\n'

                        str_list += '<td class="center">'
                        str_list += relative_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'maplot.pdf')),
                            text=relative_image(
                                prefix=path_prefix,
                                suffix='_'.join((row_dict['V1'], row_dict['V2'], 'maplot.png')),
                                text='M vs A Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2']))
                        str_list += '</td>\n'

                        str_list += '<td class="center">'
                        str_list += relative_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_volcano.pdf')),
                            text=relative_image(
                                prefix=path_prefix,
                                suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_volcano.png')),
                                text='Volcano Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2']))
                        str_list += '</td>\n'

                        str_list += '<td class="center" colspan="4"></td>\n'

                        str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

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

            str_list += 'track Alignments\n'
            str_list += 'shortLabel Alignments\n'
            str_list += 'longLabel TopHat RNA-Seq read alignments\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group alignments\n'
            str_list += '\n'

            str_list += 'track Assemblies\n'
            str_list += 'shortLabel Assemblies\n'
            str_list += 'longLabel Cuffmerge transcript structures\n'
            str_list += 'visibility full\n'
            str_list += 'superTrack on\n'
            str_list += 'group assemblies\n'
            str_list += '\n'

            str_list += 'track Coverage\n'
            str_list += 'shortLabel Coverage\n'
            str_list += 'longLabel TopHat RNA-Seq alignment coverage\n'
            str_list += 'visibility full\n'
            str_list += 'superTrack on\n'
            str_list += 'group coverage\n'
            str_list += '\n'

            str_list += 'track Deletions\n'
            str_list += 'shortLabel Deletions\n'
            str_list += 'longLabel TopHat RNA-Seq deletions\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group alignments\n'
            str_list += '\n'

            str_list += 'track Insertions\n'
            str_list += 'shortLabel Insertions\n'
            str_list += 'longLabel TopHat RNA-Seq insertions\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group alignments\n'
            str_list += '\n'

            str_list += 'track Junctions\n'
            str_list += 'shortLabel Junctions\n'
            str_list += 'longLabel TopHat RNA-Seq splice junctions\n'
            str_list += 'visibility show\n'
            str_list += 'superTrack on\n'
            str_list += 'group alignments\n'
            str_list += '\n'

            str_list += 'track Transcripts\n'
            str_list += 'shortLabel Transcripts\n'
            str_list += 'longLabel Cufflinks transcript structures\n'
            str_list += 'visibility show\n'
            str_list += 'superTrack on\n'
            str_list += 'group transcripts\n'
            str_list += '\n'

            # Sample-specific tracks

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                paired_reads_name_list = paired_reads_dict.keys()
                if not len(paired_reads_name_list):
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue
                paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for paired_reads_name in paired_reads_name_list:
                    #
                    # Add a trackDB entry for each Tophat accepted_hits.bam file.
                    #
                    # Common settings
                    str_list += 'track ' + paired_reads_name + '_alignments\n'
                    str_list += 'type bam\n'
                    str_list += 'shortLabel ' + paired_reads_name + '_alignments\n'
                    str_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq read alignments\n'
                    str_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/accepted_hits.bam\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility dense\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bam/cram - Compressed Sequence Alignment track settings
                    # None

                    # Composite track settings
                    str_list += 'parent Alignments\n'
                    str_list += '\n'

                    #
                    # Add a trackDB entry for each Tophat accepted_hits.bw file.
                    #
                    # Common settings
                    str_list += 'track ' + paired_reads_name + '_coverage\n'
                    # TODO: The bigWig type must declare the expected signal range.
                    # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                    str_list += 'type bigWig\n'
                    str_list += 'shortLabel ' + paired_reads_name + '_coverage\n'
                    str_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq alignment coverage\n'
                    str_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/accepted_hits.bw\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility full\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bigWig - Signal graphing track settings
                    str_list += 'alwaysZero on\n'
                    str_list += 'autoScale on\n'
                    str_list += 'graphTypeDefault bar\n'
                    str_list += 'maxHeightPixels 100:60:20\n'
                    # str_list += 'maxWindowToQuery 10000000\n'
                    # str_list += 'smoothingWindow 5\n'
                    # str_list += 'transformFunc NONE\n'
                    # str_list += 'viewLimits 0:45\n'
                    # str_list += 'viewLimitsMax 0:50\n'
                    # str_list += 'windowingFunction maximum\n'
                    # str_list += 'yLineMark <#>\n'
                    # str_list += 'yLineOnOff on \n'
                    # str_list += 'gridDefault on\n'

                    # Composite track settings
                    str_list += 'parent Coverage\n'
                    str_list += 'centerLabelsDense off\n'
                    str_list += '\n'

                    #
                    # Add a trackDB entry for each Tophat deletions.bb file.
                    #
                    # Common settings
                    str_list += 'track ' + paired_reads_name + '_deletions\n'
                    str_list += 'type bigBed\n'
                    str_list += 'shortLabel ' + paired_reads_name + '_deletions\n'
                    str_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq deletions\n'
                    str_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/deletions.bb\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility hide\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list += 'parent Deletions\n'
                    str_list += '\n'

                    #
                    # Add a trackDB entry for each Tophat insertions.bb file.
                    #
                    # Common settings
                    str_list += 'track insertions_' + paired_reads_name + '\n'
                    str_list += 'type bigBed\n'
                    str_list += 'shortLabel ' + paired_reads_name + '_insertions\n'
                    str_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq insertions\n'
                    str_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/insertions.bb\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility hide\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list += 'parent Insertions\n'
                    str_list += '\n'

                    #
                    # Add a trackDB entry for each Tophat junctions.bb file.
                    #
                    # Common settings
                    str_list += 'track ' + paired_reads_name + '_junctions\n'
                    str_list += 'type bigBed\n'
                    str_list += 'shortLabel ' + paired_reads_name + '_junctions\n'
                    str_list += 'longLabel ' + paired_reads_name + ' TopHat RNA-Seq splice junctions\n'
                    str_list += 'bigDataUrl rnaseq_tophat_' + paired_reads_name + '/junctions.bb\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility pack\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list += 'parent Junctions\n'
                    str_list += '\n'

                    #
                    # Add a trackDB entry for each Tophat transcripts.bb file.
                    #
                    # Common settings
                    str_list += 'track ' + paired_reads_name + '_transcripts\n'
                    str_list += 'type bigGenePred\n'
                    str_list += 'shortLabel ' + paired_reads_name + '_transcripts\n'
                    str_list += 'longLabel ' + paired_reads_name + ' Cufflinks transcript assembly\n'
                    str_list += 'bigDataUrl rnaseq_cufflinks_' + paired_reads_name + '/transcripts.bb\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility hide\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bigGenePred - Gene Annotations settings
                    # None

                    # Composite track settings
                    str_list += 'parent Transcripts\n'
                    str_list += '\n'

            # Comparison-specific tracks

            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for comparison_name in comparison_name_list:
                #
                # Add a trackDB entry for each Cuffmerge merged.bb file.
                #
                # Common settings
                str_list += 'track ' + comparison_name + '_assembly\n'
                str_list += 'type bigGenePred\n'
                str_list += 'shortLabel ' + comparison_name + '_assembly\n'
                str_list += 'longLabel ' + comparison_name + ' Cufflinks transcript assembly\n'
                str_list += 'bigDataUrl rnaseq_cuffmerge_' + comparison_name + '/merged.bb\n'
                # str_list += 'html ...\n'
                str_list += 'visibility pack\n'

                # Common optional settings
                str_list += 'color 0,0,0\n'

                # bigGenePred - Gene Annotations settings
                # None

                # Composite track settings
                str_list += 'parent Assemblies\n'
                str_list += '\n'

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

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
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode
    @ivar contrast_path: Contrast file path
    @type contrast_path: str | unicode
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode
    @ivar transcriptome_gtf_path: Reference transcriptome GTF file path
    @type transcriptome_gtf_path: str | unicode
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | unicode
    @ivar design_formula: Design formula
    @type design_formula: str | unicode
    @ivar library_type: Library type
        Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
    @type library_type: str
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
            sample_list=None,
            replicate_grouping=False,
            comparison_path=None,
            contrast_path=None,
            genome_fasta_path=None,
            transcriptome_gtf_path=None,
            transcriptome_index_path=None,
            design_formula=None,
            library_type=None):
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
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode
        @param contrast_path: Contrast file path
        @type contrast_path: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode
        @param design_formula: Design formula
        @type design_formula: str | unicode
        @param library_type: Library type
            Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
        @type library_type: str
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
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if comparison_path is None:
            self.comparison_path = str()
        else:
            self.comparison_path = comparison_path

        if contrast_path is None:
            self.contrast_path = str()
        else:
            self.contrast_path = contrast_path

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

        if library_type is None:
            self.library_type = str()
        else:
            self.library_type = library_type

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
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

        option = 'ctr_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.contrast_path = configuration.config_parser.get(section=section, option=option)

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

        option = 'library_type'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_type = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.DESeq} analysis.
        @return:
        @rtype:
        """

        super(DESeq, self).run()

        # DESeq requires a genome version.

        if not self.genome_version:
            raise Exception('A DESeq analysis requires a genome_version configuration option.')

        # For DESeq, all samples need adding to the Analysis regardless.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        # Read the designs (comparison) file.

        self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)

        design_sheet = AnnotationSheet.from_file_path(
            file_path=self.comparison_path,
            file_type='excel',
            name='DESeq Design Table')
        # Convert the CSV configuration table into a TSV analysis table.
        design_sheet.file_type = 'excel-tab'

        # Read the contrasts file.

        self.contrast_path = self.configuration.get_absolute_path(file_path=self.contrast_path)

        contrast_sheet = AnnotationSheet.from_file_path(
            file_path=self.contrast_path,
            file_type='excel',
            name='DESeq Contrast Table')
        # Convert the CSV configuration table into a TSV analysis table.
        contrast_sheet.file_type = 'excel-tab'

        # TODO: Adjust by introducing a new class RNASeqComparisonSheet(AnnotationSheet) in this module?
        for design_row_dict in design_sheet.row_dicts:
            design_name = design_row_dict['design']

            prefix = '_'.join((self.stage_name_deseq, design_name))

            comparison_directory = os.path.join(self.genome_directory, prefix)

            if not os.path.isdir(comparison_directory):
                try:
                    os.makedirs(comparison_directory)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            annotation_sheet = AnnotationSheet(
                file_path=os.path.join(comparison_directory, prefix + '_samples.tsv'),
                file_type='excel-tab',
                # file_type='excel',
                name='DESeq Sample Annotation',
                header=True)

            for sample in self.sample_list:
                if self.debug > 0:
                    print self, 'Sample name:', sample.name
                    print sample.trace(1)

                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                paired_reads_name_list = paired_reads_dict.keys()
                if not len(paired_reads_name_list):
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue
                paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                row_dict = {
                    'sample': sample.name,
                    'bam_path': StarAligner.get_prefix_star_aligner_merge(sample_name=sample.name) + '.bam',
                    'bai_path': StarAligner.get_prefix_star_aligner_merge(sample_name=sample.name) + '.bai',
                }
                """ @type row_dict: dict[str, str | unicode] """
                # Set additional columns from the Sample Annotation Sheet prefixed with 'Sample DESeq *'.
                for annotation_key in filter(lambda x: x.startswith('DESeq '), sample.annotation_dict.keys()):
                    row_dict[annotation_key[6:]] = sample.annotation_dict[annotation_key][0]

                annotation_sheet.row_dicts.append(row_dict)

            # Process all row_dict objects to get the superset of field names.
            for row_dict in annotation_sheet.row_dicts:
                for key in row_dict.iterkeys():
                    if key not in annotation_sheet.field_names:
                        annotation_sheet.field_names.append(key)

            # Write the DESeq Sample Annotation Sheet to disk.
            annotation_sheet.to_file_path()

            # Write the DESeq Design Annotation Sheet to disk.
            design_sheet.file_path = os.path.join(comparison_directory, prefix + '_designs.tsv')
            design_sheet.to_file_path()

            # Write the DESeq Contrast Annotation Sheet to disk.
            contrast_sheet.file_path = os.path.join(comparison_directory, prefix + '_contrasts.tsv')
            contrast_sheet.to_file_path()

        return
