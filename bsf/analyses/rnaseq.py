# -*- coding: utf-8 -*-
"""RNA-seq Analysis module

A package of classes and methods supporting RNA-seq analyses.
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

import errno
import os
import pickle
import re
import sys
import warnings

import bsf.analyses.hisat
import bsf.analyses.star_aligner
import bsf.analysis
import bsf.annotation
import bsf.connector
import bsf.executables
import bsf.ngs
import bsf.procedure
import bsf.process
import bsf.standards


class FilePathTophat(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathTophat} models files in a sample-specific TopHat directory.

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
        """Initialise a C{bsf.analyses.rnaseq.FilePathTophat} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
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


class FilePathCufflinks(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathCufflinks} models files in a sample-specific Cufflinks directory.

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
    @ivar temporary_slopped_tsv: Temporary slopped (bedtools slop) tab-separated value (TSV) file
    @type temporary_slopped_tsv: str | unicode
    @ivar temporary_fixed_tsv: Temporary splopped and fixed (tab at end) tab-separated value (TSV) file
    @type temporary_fixed_tsv: str | unicode
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
        """Initialise a C{bsf.analyses.rnaseq.FilePathCufflinks} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
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
        self.temporary_slopped_tsv = os.path.join(prefix, 'transcripts_slopped.tsv')
        self.temporary_fixed_tsv = os.path.join(prefix, 'transcripts_fixed.tsv')
        self.transcripts_bb = os.path.join(prefix, 'transcripts.bb')
        self.transcripts_bb_link_source = 'transcripts.bb'
        self.transcripts_bb_link_target = os.path.join(prefix, prefix + '_transcripts.bb')
        self.transcripts_gtf = os.path.join(prefix, 'transcripts.gtf')
        self.transcripts_gtf_link_source = 'transcripts.gtf'
        self.transcripts_gtf_link_target = os.path.join(prefix, prefix + '_transcripts.gtf')

        return


class FilePathCuffmerge(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathCuffmerge} models files in a comparison-specific Cuffmerge directory.

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
    @ivar cuffcompare_prefix: Cuffcompare output prefix, including the cuffmerge directory path
    @type cuffcompare_prefix: str | unicode
    @ivar cuffcompare_combined_gtf: Cuffcompare merged GTF file
    @type cuffcompare_combined_gtf: str | unicode
    @ivar cuffcompare_loci: Cuffcomapre loci file
    @type cuffcompare_loci: str | unicode
    @ivar cuffcompare_stats: Cuffcompare stats file
    @type cuffcompare_stats: str | unicode
    @ivar cuffcompare_tracking: Cuffcompare tracking file
    @type cuffcompare_tracking: str | unicode
    @ivar cuffcompare_refmap: Cuffcompare merged.gtf.refmap
    @type cuffcompare_refmap: str | unicode
    @ivar cuffcompare_tmap: Cuffcompare merged.gtf.tmap
    @type cuffcompare_tmap: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rnaseq.FilePathCuffmerge} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
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
        self.cuffcompare_prefix = os.path.join(prefix, 'cuffcmp')
        self.cuffcompare_combined_gtf = self.cuffcompare_prefix + '.combined.gtf'
        self.cuffcompare_loci = self.cuffcompare_prefix + '.loci'
        self.cuffcompare_stats = self.cuffcompare_prefix + '.stats'
        self.cuffcompare_tracking = self.cuffcompare_prefix + '.tracking'
        self.cuffcompare_refmap = self.cuffcompare_prefix + '.merged.gtf.refmap'
        self.cuffcompare_tmap = self.cuffcompare_prefix + '.merged.gtf.tmap'

        return


class FilePathCuffquant(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathCuffquant} models files in a sample-specific Cuffquant directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar abundances: Cuffquant abundances file
    @type abundances: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rnaseq.FilePathCuffquant} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathCuffquant, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances = os.path.join(prefix, 'abundances.cxb')

        return


class FilePathCuffnorm(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathCuffnorm} models files in a comparison-specific Cuffnorm directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar abundances_tsv: Abundances TSV file
    @type abundances_tsv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rnaseq.FilePathCuffnorm} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathCuffnorm, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances_tsv = prefix + '_abundances.tsv'

        return


class FilePathCuffdiff(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathCuffdiff} models files in a comparison-specific Cuffdiff directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rnaseq.FilePathCuffdiff} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathCuffdiff, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances_tsv = prefix + '_abundances.tsv'
        self.alignments_tsv = prefix + '_alignments.tsv'

        return


class FilePathProcessCuffdiff(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathProcessCuffdiff} models files in a comparison-specific Cuffdiff directory.

    Attributes:
    """

    pass


class FilePathMonocle(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathMonocle} models files in a comparison-specific Monocle directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar annotation_tsv: Monocle annotation TSV
    @type annotation_tsv: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rnaseq.FilePathMonocle} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathMonocle, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.annotation_tsv = prefix + '_annotation.tsv'

        return


class TuxedoSamplePairSheet(bsf.annotation.AnnotationSheet):
    """The C{bsf.analyses.rnaseq.TuxedoSamplePairSheet} class represents C{bsf.ngs.Sample} pairs.

    The C{bsf.ngs.Sample} pairs are defined by the C{bsf_rnaseq_process_cuffdiff.R} script.

    Attributes:
    """

    _file_type = 'excel-tab'

    _field_names = [
        'V1',
        'V2',
    ]

    _test_methods = dict()


class Tuxedo(bsf.analysis.Analysis):
    """Tuxedo RNASeq C{bsf.analysis.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool | None
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode | None
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode | None
    @ivar genome_index_path: Bowtie genome index path
    @type genome_index_path: str | unicode | None
    @ivar genome_sizes_path: Reference genome sizes file path
    @type genome_sizes_path: str | unicode | None
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str | None
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | unicode | None
    @ivar transcriptome_gtf_path: Reference transcriptome GTF file path
    @type transcriptome_gtf_path: str | unicode | None
    @ivar mask_gtf_path: GTF file path to mask transcripts
    @type mask_gtf_path: str | unicode | None
    @ivar multi_read_correction: Apply multi-read correction
    @type multi_read_correction: bool | None
    @ivar library_type: Library type
        Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
    @type library_type: str | None
    @ivar novel_transcripts: Assemble novel transcripts
    @type novel_transcripts: bool | None
    @ivar false_discovery_rate: False discovery rate (FDR) threshold
    @type false_discovery_rate: float | None
    @ivar no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
    @type no_length_correction: bool | None
    @ivar aligner: Alignment program
    @type aligner: str | None
    """

    name = 'RNA-seq Analysis'
    prefix = 'rnaseq'

    @classmethod
    def get_stage_name_run_tophat(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'run_tophat'))

    @classmethod
    def get_stage_name_process_tophat(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'process_tophat'))

    @classmethod
    def get_stage_name_run_cufflinks(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'run_cufflinks'))

    @classmethod
    def get_stage_name_process_cufflinks(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'process_cufflinks'))

    # Comparison stage

    @classmethod
    def get_stage_name_run_cuffmerge(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cuffmerge'))

    @classmethod
    def get_stage_name_run_cuffquant(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cuffquant'))

    @classmethod
    def get_stage_name_run_cuffnorm(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cuffnorm'))

    @classmethod
    def get_stage_name_run_cuffdiff(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'cuffdiff'))

    @classmethod
    def get_stage_name_process_cuffdiff(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'process_cuffdiff'))

    @classmethod
    def get_stage_name_monocle(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'monocle'))

    @classmethod
    def get_prefix_run_tophat(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: Sample name
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_run_tophat(), sample_name))

    @classmethod
    def get_prefix_process_tophat(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: Sample name
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_process_tophat(), sample_name))

    @classmethod
    def get_prefix_run_cufflinks(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: Sample name
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cufflinks(), sample_name))

    @classmethod
    def get_prefix_process_cufflinks(cls):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return cls.get_stage_name_process_cufflinks()

    @classmethod
    def get_prefix_run_cuffmerge(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffmerge(), comparison_name))

    @classmethod
    def get_prefix_run_cuffquant(cls, comparison_name, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @param sample_name: Sample name
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffquant(), comparison_name, sample_name))

    @classmethod
    def get_prefix_run_cuffnorm(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffnorm(), comparison_name))

    @classmethod
    def get_prefix_run_cuffdiff(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffdiff(), comparison_name))

    @classmethod
    def get_prefix_process_cuffdiff(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_process_cuffdiff(), comparison_name))

    @classmethod
    def get_prefix_monocle(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_monocle(), comparison_name))

    @classmethod
    def get_file_path_run_tophat(cls, sample_name):
        """Get a C{FilePathTophat} object.

        The prefix is non-standard, as I{rnaseq_run_tophat} and I{rnaseq_process_tophat}
        use the same I{rnaseq_tophat} prefix.
        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathTophat} object
        @rtype FilePathTophat
        """
        return FilePathTophat(
            prefix='_'.join(('rnaseq_tophat', sample_name)))

    @classmethod
    def get_file_path_process_tophat(cls, sample_name):
        """Get a C{FilePathTophat} object.

        The prefix is non-standard, as I{rnaseq_run_tophat} and I{rnaseq_process_tophat}
        use the same I{rnaseq_tophat} prefix.
        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathTophat} object
        @rtype FilePathTophat
        """
        return FilePathTophat(
            prefix='_'.join(('rnaseq_tophat', sample_name)))

    @classmethod
    def get_file_path_run_cufflinks(cls, sample_name):
        """Get a C{FilePathCufflinks} object.

        The prefix is non-standard, as I{rnaseq_run_cufflinks} and I{rnaseq_process_cufflinks}
        use the same I{rnaseq_cufflinks} prefix.
        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathCufflinks}
        @rtype: FilePathCufflinks
        """
        return FilePathCufflinks(
            prefix='_'.join(('rnaseq_cufflinks', sample_name)))

    @classmethod
    def get_file_path_process_cufflinks(cls, sample_name):
        """Get a C{FilePathCufflinks} object.

        The prefix is non-standard, as I{rnaseq_run_cufflinks} and I{rnaseq_process_cufflinks}
        use the same I{rnaseq_cufflinks} prefix.
        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathCufflinks}
        @rtype: FilePathCufflinks
        """
        return FilePathCufflinks(
            prefix='_'.join(('rnaseq_cufflinks', sample_name)))

    @classmethod
    def get_file_path_cuffmerge(cls, comparison_name):
        """Get a C{FilePathCuffmerge} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathCuffmerge}
        @rtype: FilePathCuffmerge
        """
        return FilePathCuffmerge(
            prefix=cls.get_prefix_run_cuffmerge(comparison_name=comparison_name))

    @classmethod
    def get_file_path_cuffquant(cls, comparison_name, sample_name):
        """Get a C{FilePathCuffquant} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathCuffquant}
        @rtype: FilePathCuffquant
        """
        return FilePathCuffquant(
            prefix=cls.get_prefix_run_cuffquant(comparison_name=comparison_name, sample_name=sample_name))

    @classmethod
    def get_file_path_cuffnorm(cls, comparison_name):
        """Get a C{FilePathCuffnorm} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathCuffnorm}
        @rtype: FilePathCuffnorm
        """
        return FilePathCuffnorm(
            prefix=cls.get_prefix_run_cuffnorm(comparison_name=comparison_name))

    @classmethod
    def get_file_path_run_cuffdiff(cls, comparison_name):
        """Get a C{FilePathCuffdiff} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathCuffdiff}
        @rtype: FilePathCuffdiff
        """
        return FilePathCuffdiff(
            prefix=cls.get_prefix_run_cuffdiff(comparison_name=comparison_name))

    @classmethod
    def get_file_path_process_cuffdiff(cls, comparison_name):
        """Get a C{FilePathProcessCuffdiff} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathProcessCuffdiff}
        @rtype: FilePathProcessCuffdiff
        """
        return FilePathProcessCuffdiff(
            prefix=cls.get_prefix_process_cuffdiff(comparison_name=comparison_name))

    @classmethod
    def get_file_path_monocle(cls, comparison_name):
        """Get a C{FilePathMonocle} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathMonocle}
        @rtype: FilePathMonocle
        """
        return FilePathMonocle(
            prefix=cls.get_prefix_monocle(comparison_name=comparison_name))

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
            comparison_path=None,
            genome_fasta_path=None,
            genome_index_path=None,
            genome_sizes_path=None,
            transcriptome_version=None,
            transcriptome_index_path=None,
            transcriptome_gtf_path=None,
            mask_gtf_path=None,
            multi_read_correction=None,
            library_type=None,
            novel_transcripts=None,
            false_discovery_rate=None,
            no_length_correction=None,
            aligner=None):
        """Initialise a C{bsf.analyses.rnaseq.Tuxedo} object.

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
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool | None
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode | None
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode | None
        @param genome_index_path: Bowtie genome index path
        @type genome_index_path: str | unicode | None
        @param genome_sizes_path: Reference genome sizes file path
        @type genome_sizes_path: str | unicode | None
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode | None
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode | None
        @param mask_gtf_path: GTF file path to mask transcripts
        @type mask_gtf_path: str | unicode | None
        @param multi_read_correction: Apply multi-read correction
        @type multi_read_correction: bool | None
        @param library_type: Library type
            Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
        @type library_type: str | None
        @param novel_transcripts: Assemble novel transcripts
        @type novel_transcripts: bool | None
        @param false_discovery_rate: False discovery rate (FDR) threshold
        @type false_discovery_rate: float | None
        @param no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
        @type no_length_correction: bool | None
        @param aligner: Alignment program
        @type aligner: str | None
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

        self.replicate_grouping = replicate_grouping
        self.comparison_path = comparison_path
        self.genome_fasta_path = genome_fasta_path
        self.genome_index_path = genome_index_path
        self.genome_sizes_path = genome_sizes_path
        self.transcriptome_version = transcriptome_version
        self.transcriptome_index_path = transcriptome_index_path
        self.transcriptome_gtf_path = transcriptome_gtf_path
        self.mask_gtf_path = mask_gtf_path
        self.multi_read_correction = multi_read_correction
        self.library_type = library_type
        self.novel_transcripts = novel_transcripts
        self.false_discovery_rate = false_discovery_rate
        self.no_length_correction = no_length_correction
        self.aligner = aligner

        self._comparison_dict = dict()
        """ @type _comparison_dict: dict[str, list[bsf.ngs.SampleGroup]] """

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
        super(Tuxedo, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_sizes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_sizes_path = configuration.config_parser.get(section=section, option=option)

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

        option = 'false_discovery_rate'
        if configuration.config_parser.has_option(section=section, option=option):
            self.false_discovery_rate = configuration.config_parser.getfloat(section=section, option=option)

        option = 'no_length_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.no_length_correction = configuration.config_parser.getboolean(section=section, option=option)

        option = 'aligner'
        if configuration.config_parser.has_option(section=section, option=option):
            self.aligner = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.rnaseq.Tuxedo} analysis.
        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            All C{bsf.ngs.Sample} objects referenced in a comparison are added from the C{bsf.ngs.Collection} to the
            C{bsf.analysis.Analysis} object.

                - Column headers for CASAVA folders:
                    - Treatment/Control/Point N ProcessedRunFolder:
                        - CASAVA processed run folder name or
                        - C{bsf.analysis.Analysis.input_directory} by default
                    - Treatment/Control/Point N Project:
                        - CASAVA Project name or
                        - C{bsf.analysis.Analysis.project_name} by default
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

                    for _group_name, _sample_list in self.collection.sample_group_dict.items():
                        _sample_group = bsf.ngs.SampleGroup(name=_group_name, sample_list=_sample_list)
                        # SampleGroup objects are only useful, if at least one Sample object is not excluded.
                        if not _sample_group.is_excluded():
                            _sample_group_list.append(_sample_group)

                    # Sort the list of comparison groups by SampleGroup.name.
                    _sample_group_list.sort(key=lambda item: item.name)
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
                            _sample_group_list.append(bsf.ngs.SampleGroup(name=_sample.name, sample_list=[_sample]))

                    # Sort the list of comparison groups by SampleGroup.name.
                    _sample_group_list.sort(key=lambda item: item.name)
                    # Set the comparison name to 'global'.
                    self._comparison_dict['global'] = _sample_group_list
                else:
                    # A comparison file path was provided.
                    self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)
                    # Read and process the comparison file, which includes adding only those Sample objects,
                    # which are referenced in a comparison.
                    annotation_sheet = bsf.annotation.AnnotationSheet.from_file_path(file_path=self.comparison_path)
                    regular_expression = re.compile(pattern='\\W')

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
                            # In RNA-seq experiments, entire pools of Sample objects (replicates) are compared
                            # with each other.
                            _group_name, _sample_list_old = self.collection.get_samples_from_row_dict(
                                row_dict=row_dict,
                                prefix=prefix)
                            _sample_list_new = list()
                            """ @type _sample_list_new: list[bsf.ngs.Sample] """
                            if _group_name and len(_sample_list_old):
                                # Sample objects are only useful, if at least one PairedReads object is not excluded.
                                for _sample in _sample_list_old:
                                    if not _sample.is_excluded():
                                        _sample_list_new.append(_sample)

                                if len(_sample_list_new):
                                    _comparison_name_list.append(_group_name)
                                    _sample_group_list.append(bsf.ngs.SampleGroup(
                                        name=_group_name,
                                        sample_list=_sample_list_new))
                                    # Also expand each Python list of bsf.ngs.Sample objects to get all those
                                    # bsf.ngs.Sample objects that this bsf.analysis.Analysis needs considering.
                                    for _sample in _sample_list_new:
                                        self.add_sample(sample=_sample)
                                        if self.debug > 1:
                                            print('  ', prefix, 'Sample name:', _sample.name,
                                                  'file_path:', _sample.file_path)
                                        if self.debug > 2:
                                            sys.stdout.writelines(sample.trace(level=1))
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

                        # Sort the list of comparison groups by SampleGroup.name.
                        _sample_group_list.sort(key=lambda item: item.name)
                        # Set the comparison name.
                        self._comparison_dict[_comparison_name] = _sample_group_list
            else:
                # Without a comparison file path, simply add a comparison 'global' with a SampleGroup 'global'
                # with all Sample objects from the Collection.
                # This means that most pipeline stages with the exception of Cuffdiff can run.
                self.sample_list.extend(self.collection.get_all_samples(exclude=True))
                self._comparison_dict['global'] = [bsf.ngs.SampleGroup(
                    name='global',
                    sample_list=self.collection.get_all_samples(exclude=True))]

            # Check for comparisons without SampleGroup objects or SampleGroup objects without Sample.
            # FIXME: The complication is that Sample or ReadGroup objects could be excluded from the Analysis.
            for _comparison_name, _sample_group_list in self._comparison_dict.items():
                if self.debug > 0:
                    print('Comparison name:', _comparison_name)
                    print('SampleGroup list:')
                if len(_sample_group_list) < 1:
                    warnings.warn('Comparison ' + _comparison_name + ' without SampleGroup objects', UserWarning)
                for _sample_group in _sample_group_list:
                    if self.debug > 0:
                        print('  SampleGroup name:', _sample_group.name)
                        print('  SampleGroup Sample list:')
                    if len(_sample_group.sample_list) < 1:
                        warnings.warn('SampleGroup ' + _sample_group.name + ' without Sample objects')
                    for _sample in _sample_group.sample_list:
                        if self.debug > 0:
                            print('    Sample name:', _sample.name)

            return

        def run_write_annotation(annotation_path, annotation_dict):
            """Private function to write a sample annotation file for Cuffdiff or Cuffnorm to disk.

            @param annotation_path: Annotation file path
            @type annotation_path: str | unicode
            @param annotation_dict: Annotation dict
            @type annotation_dict: dict[str, list[str | unicode]]
            """
            with open(file=annotation_path, mode='wt') as _annotation_file:
                _annotation_file.write('sample_id\tgroup_label\n')
                for _group_name in sorted(annotation_dict):
                    for _file_path in annotation_dict[_group_name]:
                        _annotation_file.write(_file_path + '\t' + _group_name + '\n')

            return

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Tuxedo requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A ' + self.name + " requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = bsf.standards.Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a valid 'transcriptome_version' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[self.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # Get the comparison annotation sheet before calling the run() method of the Analysis super-class.

        if not self.comparison_path:
            # A comparison path was not provided, check if a standard file exists in this directory.
            self.comparison_path = self.get_annotation_file(prefix_list=[self.prefix], suffix='comparisons.csv')
            if self.comparison_path:
                if not os.path.exists(self.comparison_path):
                    self.comparison_path = None
                    if self.debug > 0:
                        print('Standard comparison file not in current working directory:', self.comparison_path)
                else:
                    if self.debug > 0:
                        print('Standard comparison file in current working directory:', self.comparison_path)

        super(Tuxedo, self).run()

        # Method configuration with regards to Cuffquant and Cuffdiff.
        run_cuffquant_before_cuffdiff = False

        run_read_comparisons()

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        # Get the genome, FASTA, index and sizes.

        if not self.genome_index_path:
            self.genome_index_path = os.path.join(
                bsf.standards.FilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        if not self.genome_fasta_path:
            self.genome_fasta_path = bsf.standards.FilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if not os.path.exists(self.genome_fasta_path):
            raise Exception('Genome FASTA file path {!r} does not exists.'.format(self.genome_fasta_path))

        if not self.genome_sizes_path:
            self.genome_sizes_path = bsf.standards.FilePath.get_resource_genome_fasta_index(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if not os.path.exists(self.genome_sizes_path):
            raise Exception('Genome sizes file path {!r} does not exists.'.format(self.genome_sizes_path))

        # Define a reference transcriptome index directory or a GTF file path.

        if self.transcriptome_index_path:
            # Check if the transcriptome_index_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_index_path = self.configuration.get_absolute_path(
                file_path=self.transcriptome_index_path,
                default_path=bsf.standards.FilePath.get_resource_transcriptome(
                    transcriptome_version=None,
                    absolute=True))

            if not os.path.isdir(self.transcriptome_index_path):
                raise Exception('Reference transcriptome index directory {!r} does not exist.'.
                                format(self.transcriptome_index_path))

            transcriptome_index = os.path.basename(self.transcriptome_index_path)

            # Does an indices_for_TopHat directory exist?
            transcriptome_index_path = os.path.join(
                self.transcriptome_index_path,
                bsf.standards.Index.get(option='tophat2'))
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
            # prepend the default transcriptome directory.
            self.transcriptome_gtf_path = self.configuration.get_absolute_path(
                file_path=self.transcriptome_gtf_path,
                default_path=bsf.standards.FilePath.get_resource_transcriptome(
                    transcriptome_version=self.transcriptome_version,
                    absolute=True))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception('Reference transcriptome GTF file {!r} does not exist.'.
                                format(self.transcriptome_gtf_path))
        else:
            # Neither was provided, automatically discover on the basis of the transcriptome version.
            self.transcriptome_index_path = os.path.join(
                bsf.standards.FilePath.get_resource_transcriptome_index(
                    transcriptome_version=self.transcriptome_version,
                    transcriptome_index='tophat2'),
                self.transcriptome_version,  # TopHat puts the transcriptome index into a sub directory.
                self.transcriptome_version)

            self.transcriptome_gtf_path = bsf.standards.FilePath.get_resource_transcriptome_gtf(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='tophat2',
                basic=False)

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception('Reference transcriptome GTF file path {!r} does not exist.'.
                                format(self.transcriptome_gtf_path))

        if not self.transcriptome_gtf_path:
            raise Exception('Reference transcriptome GTF file not defined.\n' +
                            'A ' + self.name + " requires a 'transcriptome_index' or 'transcriptome_gtf' " +
                            "configuration option.")

        if not self.library_type:
            raise Exception('A ' + self.name + " requires a 'library_type' configuration option.")

        library_type_tuple = ('fr-unstranded', 'fr-firststrand', 'fr-secondstrand')
        if self.library_type not in library_type_tuple:
            raise Exception("Invalid 'library_type' configuration option: " + self.library_type + '\n' +
                            'Supported types: ' + repr(library_type_tuple))

        if self.aligner:
            # Defaults to '',
            aligner_tuple = ('hisat2', 'star', 'tophat2')
            if self.aligner not in aligner_tuple:
                raise Exception("Invalid 'aligner' configuration option: " + self.aligner + '\n' +
                                'Supported aligners: ' + repr(aligner_tuple))

        # Read configuration options.

        # TODO: Move the ConfigParser code.
        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        stage_run_tophat = self.get_stage(name=self.get_stage_name_run_tophat())
        stage_process_tophat = self.get_stage(name=self.get_stage_name_process_tophat())
        stage_run_cufflinks = self.get_stage(name=self.get_stage_name_run_cufflinks())
        stage_process_cufflinks = self.get_stage(name=self.get_stage_name_process_cufflinks())
        stage_run_cuffmerge = self.get_stage(name=self.get_stage_name_run_cuffmerge())
        stage_run_cuffquant = self.get_stage(name=self.get_stage_name_run_cuffquant())
        stage_run_cuffnorm = self.get_stage(name=self.get_stage_name_run_cuffnorm())
        stage_run_cuffdiff = self.get_stage(name=self.get_stage_name_run_cuffdiff())
        stage_process_cuffdiff = self.get_stage(name=self.get_stage_name_process_cuffdiff())
        # stage_monocle = self.get_stage(name=self.get_stage_name_monocle())

        runnable_run_cufflinks_list = list()
        """ @type runnable_run_cufflinks_list: list[Runnable] """

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            if self.aligner == 'hisat2':
                run_cufflinks_dependency = bsf.analyses.hisat.Hisat2.get_prefix_sample(
                    sample_name=sample.name)
            elif self.aligner == 'star':
                run_cufflinks_dependency = bsf.analyses.star_aligner.StarAligner.get_prefix_sample(
                    sample_name=sample.name)
            else:  # tophat2 is the default case.
                # Create a Tophat Runnable per Sample.name.

                # TODO: Activate the new code once the bsf_run_rnaseq_tophat.py script has been retired.

                prefix_run_tophat = self.get_prefix_run_tophat(sample_name=sample.name)

                file_path_run_tophat = self.get_file_path_run_tophat(sample_name=sample.name)

                # runnable_run_tophat = self.add_runnable(
                #         runnable=bsf.procedure.ConsecutiveRunnable(
                #                 name=self.get_prefix_rnaseq_run_tophat(sample_name=sample.name),
                #                 code_module='bsf.runnables.generic',
                #                 working_directory=self.genome_directory,
                #                 debug=self.debug))
                # executable_run_tophat = self.set_stage_runnable(
                #         stage=stage_run_tophat,
                #         runnable=runnable_run_tophat)
                #
                # Create a new Tophat bsf.process.RunnableStep.

                runnable_step = bsf.process.RunnableStep(
                    name='tophat2',
                    program='tophat2')
                # runnable_run_tophat.add_runnable_step(runnable_step=runnable_step)

                # Read configuration section [bsf.analyses.rnaseq.Tuxedo.tophat2]
                self.set_runnable_step_configuration(runnable_step=runnable_step)

                # Set tophat options.

                runnable_step.add_option_long(
                    key='GTF',
                    value=self.transcriptome_gtf_path)
                if self.transcriptome_index_path:
                    runnable_step.add_option_long(
                        key='transcriptome-index',
                        value=self.transcriptome_index_path)
                runnable_step.add_option_long(
                    key='output-dir',
                    value=file_path_run_tophat.output_directory)
                runnable_step.add_option_long(
                    key='num-threads',
                    value=str(stage_run_tophat.threads))
                # TODO: These really are properties of the Reads, PairedReads or Sample objects
                # rather than an Analysis.
                # TODO: Move the ConfigParser code.
                if config_parser.has_option(section=config_section, option='insert_size'):
                    insert_size = config_parser.getint(section=config_section, option='insert_size')
                    read_length = config_parser.getint(section=config_section, option='read_length')
                    mate_inner_dist = insert_size - 2 * read_length
                    runnable_step.add_option_long(
                        key='mate-inner-dist',
                        value=str(mate_inner_dist))
                # TODO: Move the ConfigParser code.
                if config_parser.has_option(section=config_section, option='mate-std-dev'):
                    runnable_step.add_option_long(
                        key='mate-std-dev',
                        value=config_parser.getint(section=config_section, option='mate-std-dev'))
                if self.library_type:
                    runnable_step.add_option_long(
                        key='library-type',
                        value=self.library_type)
                # The TopHat coverage search finds additional 'GT-AG' introns, but is only recommended for
                # short reads (< 45 bp) and small read numbers (<= 10 M).
                # TODO: This option should possibly become configurable per sample.
                runnable_step.add_switch_long(key='no-coverage-search')
                # TODO: Set -rg-* options to back fill the read group from Illumina2bam.

                # Set rnaseq_tophat arguments.

                runnable_step.arguments.append(self.genome_index_path)

                # Set rnaseq_tophat arguments for reads1 and reads2.

                reads_1_file_path_list = list()
                """ @type reads_1_file_path_list: list[str | unicode] """
                reads_2_file_path_list = list()
                """ @type reads_2_file_path_list: list[str | unicode] """

                for paired_reads_name in sorted(paired_reads_dict):
                    for paired_reads in paired_reads_dict[paired_reads_name]:
                        if self.debug > 0:
                            print(self, 'PairedReads name:', paired_reads.get_name())

                        if paired_reads.reads_1 is not None:
                            reads_1_file_path_list.append(paired_reads.reads_1.file_path)
                        if paired_reads.reads_2 is not None:
                            reads_2_file_path_list.append(paired_reads.reads_2.file_path)

                # Pass lists of files into Tophat, regardless of whether read 2 exists.
                # The bsf_run_rnaseq_tophat.py script eventually removes an empty read 2 argument.
                runnable_step.arguments.append(','.join(reads_1_file_path_list))
                runnable_step.arguments.append(','.join(reads_2_file_path_list))

                # Create a new rnaseq_run_tophat bsf.process.Executable.
                # TODO: The following code block is required as long as the bsf_run_rnaseq_tophat.py script
                # has not been retired.

                if self.debug > 0:
                    print('Tophat bsf.process.Executable')
                    sys.stdout.writelines(runnable_step.trace(level=1))

                pickler_dict_run_tophat = {
                    'prefix': stage_run_tophat.name,
                    'replicate_key': sample.name,
                    'runnable_step': runnable_step,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    stage_run_tophat.name + '_' + sample.name + '.pkl')
                with open(file=pickler_path, mode='wb') as pickler_file:
                    pickler = pickle.Pickler(file=pickler_file, protocol=pickle.HIGHEST_PROTOCOL)
                    pickler.dump(pickler_dict_run_tophat)

                executable_run_tophat = stage_run_tophat.add_executable(
                    executable=bsf.process.Executable(
                        name=prefix_run_tophat,
                        program='bsf_run_rnaseq_tophat.py'))
                # Set dependencies on previous Runnable or bsf.process.Executable objects.
                # None.

                # Set rnaseq_run_tophat options.
                executable_run_tophat.add_option_long(key='pickler_path', value=pickler_path)
                executable_run_tophat.add_option_long(key='debug', value=str(self.debug))

                # Only submit this bsf.process.Executable if the 'align_summary.txt' file does not exist.
                file_path_temporary = os.path.join(self.genome_directory, file_path_run_tophat.align_summary)
                if os.path.exists(file_path_temporary) and os.path.getsize(file_path_temporary) > 0:
                    executable_run_tophat.submit = False

                # TODO: End of code block.

                # Create a process_tophat Runnable per sample.name.

                prefix_process_tophat = self.get_prefix_process_tophat(sample_name=sample.name)

                runnable_process_tophat = self.add_runnable_consecutive(runnable=bsf.procedure.ConsecutiveRunnable(
                    name=prefix_process_tophat,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    debug=self.debug))
                executable_process_tophat = self.set_stage_runnable(
                    stage=stage_process_tophat,
                    runnable=runnable_process_tophat)
                # Set dependencies on previous Runnable or bsf.process.Executable objects.
                executable_process_tophat.dependencies.append(executable_run_tophat.name)

                # TODO: Switch from an external Bash script to a set of Runnable and RunnableStep objects.
                # Since the Bash script includes Perl code to reset the BED score field to 0, rather than
                # re-scale it properly, it would be good to write a new bsf.runnables.process_tophat module
                # to implement this in Python code.

                runnable_step = bsf.process.RunnableStep(
                    name='process_tophat',
                    program='bsf_rnaseq_process_tophat2.bash')
                runnable_process_tophat.add_runnable_step(runnable_step=runnable_step)

                # Read configuration section [bsf.analyses.rnaseq.Tuxedo.process_tophat]
                self.set_runnable_step_configuration(runnable_step=runnable_step)

                runnable_step.arguments.append(file_path_run_tophat.output_directory)
                runnable_step.arguments.append(self.genome_fasta_path)
                runnable_step.arguments.append(self.genome_sizes_path)

                # Only submit this bsf.process.Executable if the 'accepted_hits.bam.bai' file does not exist.
                # NOTE: This test is only required for older projects.
                # Newer ones make use of the Runnable and RunnableStep infrastructure.
                file_path_temporary = os.path.join(
                    self.genome_directory,
                    file_path_run_tophat.accepted_hits_bai)
                if os.path.exists(file_path_temporary) and os.path.getsize(file_path_temporary):
                    executable_process_tophat.submit = False

                run_cufflinks_dependency = executable_run_tophat.name

            # Create a run_cufflinks Runnable per Sample name.

            prefix_run_cufflinks = self.get_prefix_run_cufflinks(sample_name=sample.name)

            file_path_cufflinks = self.get_file_path_run_cufflinks(sample_name=sample.name)

            runnable_run_cufflinks = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=prefix_run_cufflinks,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    debug=self.debug))
            # Set dependencies for subsequent Runnable or bsf.process.Executable objects.
            runnable_run_cufflinks_list.append(runnable_run_cufflinks)
            executable_run_cufflinks = self.set_stage_runnable(
                stage=stage_run_cufflinks,
                runnable=runnable_run_cufflinks)
            # Set dependencies on previous Runnable or bsf.process.Executable objects.
            executable_run_cufflinks.dependencies.append(run_cufflinks_dependency)

            # Create a new Cufflinks bsf.process.RunnableStep.

            runnable_step = bsf.process.RunnableStep(
                name='cufflinks',
                program='cufflinks')
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # General Options:
            # --output-dir write all output files to this directory [.]
            runnable_step.add_option_long(
                key='output-dir',
                value=file_path_cufflinks.output_directory)
            # --num-threads number of threads used during analysis [1]
            runnable_step.add_option_long(
                key='num-threads',
                value=str(stage_run_cufflinks.threads))
            # --seed value of random number generator seed [0]
            # Cufflinks has a GTF option, in which case it will not assemble
            # novel transcripts and a GTF-guide option in which case it will
            # assemble novel transcripts.
            if self.novel_transcripts:
                # --GTF-guide use reference transcript annotation to guide assembly [NULL]
                runnable_step.add_option_long(
                    key='GTF-guide',
                    value=self.transcriptome_gtf_path)
            else:
                # --GTF quantify against reference transcript annotations [NULL]
                runnable_step.add_option_long(
                    key='GTF',
                    value=self.transcriptome_gtf_path)
            # --mask-file ignore all alignments within transcripts in this file
            if self.mask_gtf_path:
                runnable_step.add_option_long(
                    key='mask-file',
                    value=self.mask_gtf_path)
            # --frag-bias-correct use bias correction - reference fasta required [NULL]
            runnable_step.add_option_long(
                key='frag-bias-correct',
                value=self.genome_fasta_path)
            # --multi-read-correct use 'rescue method' for multi-reads (more accurate) [FALSE]
            if self.multi_read_correction:
                runnable_step.add_switch_long(
                    key='multi-read-correct')
            # --library-type library prep used for input reads [fr-unstranded]
            if self.library_type:
                runnable_step.add_option_long(
                    key='library-type',
                    value=self.library_type)
            # --library-norm-method Method used to normalize library sizes [classic-fpkm]

            # Advanced Abundance Estimation Options:
            # --frag-len-mean average fragment length (unpaired reads only) [200]
            # --frag-len-std-dev fragment length std deviation (unpaired reads only) [80]
            # --max-mle-iterations maximum iterations allowed for MLE calculation [5000]
            # --compatible-hits-norm count hits compatible with reference RNAs only [FALSE]
            # --total-hits-norm count all hits for normalization [TRUE]
            # --num-frag-count-draws Number of fragment generation samples [100]
            # --num-frag-assign-draws Number of fragment assignment samples per generation [50]
            # --max-frag-multihits Maximum number of alignments allowed per fragment [unlim]
            # --no-effective-length-correction No effective length correction [FALSE]
            # --no-length-correction No length correction [FALSE]
            if self.no_length_correction:
                runnable_step.add_switch_long(
                    key='no-length-correction')

            # Advanced Assembly Options:
            # ...

            # Advanced Reference Annotation Guided Assembly Options:
            # ...

            # Advanced Program Behavior Options:
            # --verbose log-friendly verbose processing (no progress bar) [FALSE]
            # --quiet log-friendly quiet processing (no progress bar) [FALSE]
            runnable_step.add_switch_long(
                key='quiet')
            # --no-update-check do not contact server to check for update availability [FALSE]
            runnable_step.add_switch_long(
                key='no-update-check')

            # Set Cufflinks arguments.

            if self.aligner == "hisat2":
                runnable_step.arguments.append(
                    bsf.analyses.hisat.Hisat2.get_file_path_sample(sample_name=sample.name).sample_bam)
            elif self.aligner == 'star':
                runnable_step.arguments.append(
                    bsf.analyses.star_aligner.StarAligner.get_file_path_sample(sample_name=sample.name).sample_bam)
            else:
                runnable_step.arguments.append(
                    self.get_file_path_run_tophat(sample_name=sample.name).accepted_hits_bam)

            # Convert the resulting transcripts GTF file into a UCSC genePred file.

            runnable_step = bsf.process.RunnableStep(
                name='gtf_to_gp',
                program='gtfToGenePred')
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_short(key='genePredExt')
            runnable_step.arguments.append(file_path_cufflinks.transcripts_gtf)
            runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)

            # Convert the UCSC genePred into a UCSC bigGenePred file.

            runnable_step = bsf.process.RunnableStep(
                name='gp_to_bgp',
                program='genePredToBigGenePred',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_gene_prediction,
                ])
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)
            runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)

            # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

            runnable_step = bsf.process.RunnableStep(
                name='bed_sort',
                program='bedSort',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_big_gene_prediction,
                ])
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)
            runnable_step.arguments.append(file_path_cufflinks.temporary_sorted_tsv)

            # Run bedtools slop on the sorted UCSC bigGenePred file to constrain to chromosome coordinates.

            runnable_step = bsf.process.RunnableStep(
                name='bedtools_slop',
                program='bedtools',
                sub_command=bsf.process.Command(program='slop'),
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_sorted_tsv,
                ],
                stdout=bsf.connector.ConnectorFile(file_path=file_path_cufflinks.temporary_slopped_tsv, file_mode='wt'))
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.sub_command.add_option_short(key='b', value='0')
            runnable_step.sub_command.add_option_short(key='i', value=file_path_cufflinks.temporary_sorted_tsv)
            runnable_step.sub_command.add_option_short(key='g', value=self.genome_sizes_path)

            # Since bedtools slop, at least in version v2.27.1, looses the last tab character,
            # it has to be put back with sed, before UCSC bedToBigBed can run.

            runnable_step = bsf.process.RunnableStep(
                name='sed',
                program='sed',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_slopped_tsv,
                ],
                stdout=bsf.connector.ConnectorFile(file_path=file_path_cufflinks.temporary_fixed_tsv, file_mode='wt'))
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_short(key='e', value='s/$/\\t/')

            runnable_step.arguments.append(file_path_cufflinks.temporary_slopped_tsv)

            # Convert the sorted UCSC bigGenePred into a bigBed file.

            runnable_step = bsf.process.RunnableStep(
                name='bgp_to_bb',
                program='bedToBigBed',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_fixed_tsv,
                ])
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # TODO: The location of the autoSQL file needs to be configurable.
            runnable_step.add_option_pair_short(key='as', value='/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
            runnable_step.add_switch_short(key='tab')
            runnable_step.add_option_pair_short(key='type', value='bed12+8')
            runnable_step.arguments.append(file_path_cufflinks.temporary_fixed_tsv)
            runnable_step.arguments.append(self.genome_sizes_path)
            runnable_step.arguments.append(file_path_cufflinks.transcripts_bb)

            # Add a symbolic link for the transcripts bigBed file, that includes a sample name prefix.
            runnable_step = bsf.process.RunnableStepLink(
                name='link_transcripts_bb',
                source_path=file_path_cufflinks.transcripts_bb_link_source,
                target_path=file_path_cufflinks.transcripts_bb_link_target)
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # Add a symbolic link for the transcripts GTF file, that includes a sample name prefix.
            runnable_step = bsf.process.RunnableStepLink(
                name='link_transcripts_gtf',
                source_path=file_path_cufflinks.transcripts_gtf_link_source,
                target_path=file_path_cufflinks.transcripts_gtf_link_target)
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

        # Create one process_cufflinks bsf.process.Executable to process all sub-directories.

        if len(runnable_run_cufflinks_list):
            prefix_process_cufflinks = self.get_prefix_process_cufflinks()

            runnable_process_cufflinks = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=prefix_process_cufflinks,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    debug=self.debug))
            executable_process_cufflinks = self.set_stage_runnable(
                stage=stage_process_cufflinks,
                runnable=runnable_process_cufflinks)
            # Set dependencies on previous Runnable or bsf.process.Executable objects.
            for runnable_run_cufflinks in runnable_run_cufflinks_list:
                executable_process_cufflinks.dependencies.append(runnable_run_cufflinks.name)

            runnable_step = bsf.process.RunnableStep(
                name='process_cufflinks',
                program='bsf_rnaseq_process_cufflinks.R')
            runnable_process_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # Read configuration section [bsf.analyses.rnaseq.Tuxedo.process_cufflinks]
            self.set_runnable_step_configuration(runnable_step=runnable_step)

            runnable_step.add_option_long(
                key='gtf-reference',
                value=self.transcriptome_gtf_path)
            runnable_step.add_option_long(
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

        for comparison_name in sorted(self._comparison_dict):
            if self.debug > 0:
                print('  Comparison name:', comparison_name)

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
            prefix_run_cuffmerge = self.get_prefix_run_cuffmerge(comparison_name=comparison_name)

            file_path_cuffmerge = self.get_file_path_cuffmerge(comparison_name=comparison_name)

            runnable_run_cuffmerge = self.add_runnable_consecutive(
                runnable=bsf.procedure.ConsecutiveRunnable(
                    name=prefix_run_cuffmerge,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    debug=self.debug))
            executable_run_cuffmerge = self.set_stage_runnable(
                stage=stage_run_cuffmerge,
                runnable=runnable_run_cuffmerge)
            # Submit the bsf.process.Executable if the status file AND the sample group list above supports it.
            executable_run_cuffmerge.submit &= cuffmerge_cuffnorm_submit
            # Set a dependency on all other Cuffmerge process to avoid file contention.
            executable_cuffmerge_dict[prefix_run_cuffmerge] = executable_run_cuffmerge

            if self.novel_transcripts:
                # Create a new Cuffmerge bsf.process.RunnableStep.

                runnable_step_cuffmerge = bsf.process.RunnableStep(
                    name='cuffmerge',
                    program='cuffmerge')
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step_cuffmerge)

                # Set rnaseq_cuffmerge options.

                # --output-dir Directory where merged assembly will be written [./merged_asm]
                runnable_step_cuffmerge.add_option_long(
                    key='output-dir',
                    value=file_path_cuffmerge.output_directory)
                # --ref-gtf An optional "reference" annotation GTF [NULL]
                runnable_step_cuffmerge.add_option_long(
                    key='ref-gtf',
                    value=self.transcriptome_gtf_path)
                # --ref-sequence <seq_dir>/<seq_fasta> Genomic DNA sequences for the reference
                runnable_step_cuffmerge.add_option_long(
                    key='ref-sequence',
                    value=self.genome_fasta_path)
                # --min-isoform-fraction <0-1.0> Discard isoforms with abundance below this [0.05]
                # --num-threads Use this many threads to merge assemblies [1]
                runnable_step_cuffmerge.add_option_long(
                    key='num-threads',
                    value=str(stage_run_cuffmerge.threads))
                # --keep-tmp Keep all intermediate files during merge [FALSE]

                # Set rnaseq_cuffmerge arguments.

                # Add the assembly manifest file as Cuffmerge argument.
                # The file will be written below.
                runnable_step_cuffmerge.arguments.append(file_path_cuffmerge.assembly_txt)
            else:
                # If novel transcripts are not assembled, create RunnableStep objects to create the output directory
                # and copy the reference transcriptome GTF file.

                runnable_step = bsf.process.RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=file_path_cuffmerge.output_directory)
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

                runnable_step = bsf.process.RunnableStepCopy(
                    name='copy',
                    source_path=self.transcriptome_gtf_path,
                    target_path=file_path_cuffmerge.merged_gtf)
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

                # Run cuffcompare in a self-comparison mode to get 'tss_id' and 'p_id' attributes populated.

                runnable_step = bsf.process.RunnableStep(
                    name='cuffcompare',
                    program='cuffcompare',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.merged_gtf,
                    ])
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_switch_short(key='C')  # include 'contained' transcripts
                runnable_step.add_switch_short(key='G')  # generic GFF input fields, i.e. not a Cufflinks GTF
                runnable_step.add_option_short(key='o', value=file_path_cuffmerge.cuffcompare_prefix)
                runnable_step.add_option_short(key='r', value=self.transcriptome_gtf_path)  # reference GTF
                runnable_step.add_option_short(key='s', value=self.genome_fasta_path)  # reference sequence
                runnable_step.arguments.append(file_path_cuffmerge.merged_gtf)

                # Move 'cuffcmp.combined.gtf' to 'merged.gtf' and delete obsolete files.
                runnable_step = bsf.process.RunnableStepMove(
                    name='move_combined_gtf',
                    source_path=file_path_cuffmerge.cuffcompare_combined_gtf,
                    target_path=file_path_cuffmerge.merged_gtf,
                    obsolete_file_path_list=[
                        file_path_cuffmerge.cuffcompare_loci,
                        file_path_cuffmerge.cuffcompare_stats,
                        file_path_cuffmerge.cuffcompare_tracking,
                        file_path_cuffmerge.cuffcompare_refmap,
                        file_path_cuffmerge.cuffcompare_tmap,
                    ])
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            # Convert the resulting merged GTF file into a UCSC genePred file.

            runnable_step = bsf.process.RunnableStep(
                name='gtf_to_gp',
                program='gtfToGenePred')
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_short(key='genePredExt')
            runnable_step.arguments.append(file_path_cuffmerge.merged_gtf)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)

            # Convert the UCSC genePred into a UCSC bigGenePred file.

            runnable_step = bsf.process.RunnableStep(
                name='gp_to_bgp',
                program='genePredToBigGenePred',
                obsolete_file_path_list=[
                    file_path_cuffmerge.temporary_gene_prediction,
                ])
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)

            # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

            runnable_step = bsf.process.RunnableStep(
                name='bed_sort',
                program='bedSort',
                obsolete_file_path_list=[
                    file_path_cuffmerge.temporary_big_gene_prediction,
                ])
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)

            # Convert the sorted UCSC bigGenePred into a bigBed file.

            runnable_step = bsf.process.RunnableStep(
                name='bgp_to_bb',
                program='bedToBigBed',
                obsolete_file_path_list=[
                    file_path_cuffmerge.temporary_sorted_tsv,
                ])
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            # TODO: The location of the autoSQL file needs to be configurable.
            runnable_step.add_option_pair_short(key='as', value='/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
            runnable_step.add_switch_short(key='tab')
            runnable_step.add_option_pair_short(key='type', value='bed12+8')
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)
            runnable_step.arguments.append(self.genome_sizes_path)
            runnable_step.arguments.append(file_path_cuffmerge.merged_bb)

            # Add a symbolic link for the merged bigBed file, that includes a comparison name prefix.
            runnable_step = bsf.process.RunnableStepLink(
                name='link_merged_bb',
                source_path=file_path_cuffmerge.merged_bb_link_source,
                target_path=file_path_cuffmerge.merged_bb_link_target)
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            # Add a symbolic link for the merged GTF file, that includes a comparison name prefix.
            runnable_step = bsf.process.RunnableStepLink(
                name='link_merged_gtf',
                source_path=file_path_cuffmerge.merged_gtf_link_source,
                target_path=file_path_cuffmerge.merged_gtf_link_target)
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            file_path_monocle = self.get_file_path_monocle(comparison_name=comparison_name)

            monocle_annotation_sheet = bsf.annotation.AnnotationSheet(
                file_path=os.path.join(self.genome_directory, file_path_monocle.annotation_tsv),
                file_type='excel-tab',
                header=True)

            for sample_group in sample_group_list:
                if self.debug > 0:
                    print('    SampleGroup name:', sample_group.name)

                per_group_abundances_list = list()
                """ @type per_group_abundances_list: list[str | unicode] """
                per_group_alignments_list = list()
                """ @type per_group_alignments_list: list[str | unicode] """

                for sample in sample_group.sample_list:
                    if self.debug > 0:
                        print('      Sample name:', sample.name)

                    paired_reads_dict = sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping,
                        exclude=True)

                    if not paired_reads_dict:
                        # Skip Sample objects, which PairedReads objects have all been excluded.
                        continue

                    # Add the Cufflinks assembled transcripts GTF to the Cuffmerge manifest.
                    cuffmerge_transcript_gtf_list.append(
                        os.path.join('_'.join(('rnaseq_cufflinks', sample.name)), 'transcripts.gtf') + '\n')

                    # Wait for each Cufflinks PairedReads to finish, before Cuffmerge can run.

                    executable_run_cuffmerge.dependencies.append(
                        self.get_prefix_run_cufflinks(sample_name=sample.name))

                    # Create a Cuffquant Runnable per comparison (comparison_name) and
                    # Sample.name on the basis of the Cuffmerge GTF file.

                    file_path_run_cuffquant = self.get_file_path_cuffquant(
                        comparison_name=comparison_name,
                        sample_name=sample.name)

                    runnable_run_cuffquant = self.add_runnable_consecutive(
                        runnable=bsf.procedure.ConsecutiveRunnable(
                            name=self.get_prefix_run_cuffquant(
                                comparison_name=comparison_name,
                                sample_name=sample.name),
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            debug=self.debug))
                    executable_run_cuffquant = self.set_stage_runnable(
                        stage=stage_run_cuffquant,
                        runnable=runnable_run_cuffquant)
                    # Each Cuffquant process depends on Cuffmerge.
                    executable_run_cuffquant.dependencies.append(executable_run_cuffmerge.name)

                    # Create a new cuffquant bsf.process.RunnableStep.

                    runnable_step_cuffquant = bsf.process.RunnableStep(
                        name='cuffquant',
                        program='cuffquant')
                    runnable_run_cuffquant.add_runnable_step(runnable_step=runnable_step_cuffquant)

                    # Set Cuffquant options.

                    # General Options:
                    # --output-dir write all output files to this directory [.]
                    runnable_step_cuffquant.add_option_long(
                        key='output-dir',
                        value=file_path_run_cuffquant.output_directory)
                    # --mask-file ignore all alignment within transcripts in this file [NULL]
                    if self.mask_gtf_path:
                        runnable_step_cuffquant.add_option_long(
                            key='mask-file',
                            value=self.mask_gtf_path)
                    # --frag-bias-correct use bias correction - reference fasta required [NULL]
                    runnable_step_cuffquant.add_option_long(
                        key='frag-bias-correct',
                        value=self.genome_fasta_path)
                    # --multi-read-correct use 'rescue method' for multi-reads [FALSE]
                    if self.multi_read_correction:
                        runnable_step_cuffquant.add_switch_long(
                            key='multi-read-correct')
                    # --num-threads number of threads used during quantification [1]
                    runnable_step_cuffquant.add_option_long(
                        key='num-threads',
                        value=str(stage_run_cuffquant.threads))
                    # --library-type Library prep used for input reads [fr-unstranded]
                    if self.library_type:
                        runnable_step_cuffquant.add_option_long(
                            key='library-type',
                            value=self.library_type)

                    # Advanced Options:
                    # --frag-len-mean average fragment length (unpaired reads only) [200]
                    # --frag-len-std-dev fragment length std deviation (unpaired reads only) [80]
                    # --min-alignment-count minimum number of alignments in a locus for testing [10]
                    # --max-mle-iterations maximum iterations allowed for MLE calculation [5000]
                    # --verbose log-friendly verbose processing (no progress bar) [FALSE]
                    # --quiet log-friendly quiet processing (no progress bar) [FALSE]
                    runnable_step_cuffquant.add_switch_long(
                        key='quiet')
                    # --seed value of random number generator seed [0]
                    # --no-update-check do not contact server to check for update availability [FALSE]
                    runnable_step_cuffquant.add_switch_long(
                        key='no-update-check')
                    # --max-bundle-frags maximum fragments allowed in a bundle before skipping [500000]
                    # --max-frag-multihits Maximum number of alignments allowed per fragment [unlim]
                    # --no-effective-length-correction No effective length correction [FALSE]
                    # --no-length-correction No length correction [FALSE]
                    if self.no_length_correction:
                        runnable_step_cuffquant.add_switch_long(
                            key='no-length-correction')

                    # Set Cuffquant arguments.
                    # Add the Cuffmerge GTF file and the TopHat BAM file as Cuffquant arguments.
                    # Add the TopHat BAM file to the Cuffdiff alignments list.

                    runnable_step_cuffquant.arguments.append(file_path_cuffmerge.merged_gtf)
                    if self.aligner == 'hisat2':
                        runnable_step_cuffquant.arguments.append(
                            bsf.analyses.hisat.Hisat2.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                        per_group_alignments_list.append(
                            bsf.analyses.hisat.Hisat2.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                    elif self.aligner == 'star':
                        runnable_step_cuffquant.arguments.append(
                            bsf.analyses.star_aligner.StarAligner.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                        per_group_alignments_list.append(
                            bsf.analyses.star_aligner.StarAligner.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                    else:
                        runnable_step_cuffquant.arguments.append(
                            self.get_file_path_run_tophat(sample_name=sample.name).accepted_hits_bam)
                        per_group_alignments_list.append(
                            self.get_file_path_run_tophat(sample_name=sample.name).accepted_hits_bam)

                    # Add the Cuffquant abundances file to the Cuffdiff list.

                    per_group_abundances_list.append(file_path_run_cuffquant.abundances)

                    # Add the Cuffquant Runnable process name to the Cuffdiff and Cuffnorm dependencies list.

                    cuffdiff_cuffnorm_dependencies.append(executable_run_cuffquant.name)

                    # Write Monocle annotation.
                    # FIXME: ReadGroup versus Sample level
                    # Depending on the replicate_grouping instance variable, the abundance file can be on
                    # the read_group or sample level, while Monocle annotation will always be on the sample level.
                    monocle_row_dict = {
                        'file': file_path_run_cuffquant.abundances,
                        # 'sample_name' is used by Monocle in the plot_cell_clusters() function internally.
                        'original_name': sample.name
                    }
                    """ @type monocle_row_dict: dict[str, str | unicode] """
                    # Set additional columns from the Sample Annotation Sheet prefixed with 'Sample Monocle *'.
                    for annotation_key in filter(
                            lambda x: x.startswith('Monocle '), sample.annotation_dict.keys()):
                        monocle_row_dict[annotation_key[8:]] = sample.annotation_dict[annotation_key][0]

                    monocle_annotation_sheet.row_dicts.append(monocle_row_dict)

                cuffdiff_cuffnorm_abundances_dict[sample_group.name] = per_group_abundances_list
                cuffdiff_cuffnorm_alignments_dict[sample_group.name] = per_group_alignments_list

            if self.novel_transcripts:
                # Write a Cuffmerge assembly manifest file to merge all transcriptome GTF files of each Sample object.
                # This requires an absolute path, because the working directory is not set at the stage of
                # job submission.
                with open(
                        file=os.path.join(self.genome_directory, file_path_cuffmerge.assembly_txt),
                        mode='wt') as assembly_file:
                    assembly_file.writelines(cuffmerge_transcript_gtf_list)

            if len(self._comparison_dict[comparison_name]) >= 2:
                # Create a Cuffnorm Runnable per comparison, if there are at least two SampleGroup objects.

                file_path_run_cuffnorm = self.get_file_path_cuffnorm(comparison_name=comparison_name)

                runnable_run_cuffnorm = self.add_runnable_consecutive(
                    runnable=bsf.procedure.ConsecutiveRunnable(
                        name=self.get_prefix_run_cuffnorm(comparison_name=comparison_name),
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        debug=self.debug))
                executable_run_cuffnorm = self.set_stage_runnable(
                    stage=stage_run_cuffnorm,
                    runnable=runnable_run_cuffnorm)
                # Submit the bsf.process.Executable if the status file AND the sample group list above supports it.
                executable_run_cuffnorm.submit &= cuffmerge_cuffnorm_submit
                executable_run_cuffnorm.dependencies.extend(cuffdiff_cuffnorm_dependencies)

                # Create a new Cuffnorm bsf.process.RunnableStep.

                runnable_step_cuffnorm = bsf.process.RunnableStep(
                    name='cuffnorm',
                    program='cuffnorm')
                runnable_run_cuffnorm.add_runnable_step(runnable_step=runnable_step_cuffnorm)

                # Set Cuffnorm options.

                # General Options:
                # --output-dir write all output files to this directory [.]
                runnable_step_cuffnorm.add_option_long(
                    key='output-dir',
                    value=file_path_run_cuffnorm.output_directory)
                # --labels comma-separated list of condition labels []
                # --norm-standards-file Housekeeping/spike genes to normalize libraries [NULL]
                # --num-threads number of threads used during quantification [1]
                runnable_step_cuffnorm.add_option_long(
                    key='num-threads',
                    value=str(stage_run_cuffnorm.threads))
                # --library-type Library prep used for input reads [fr-unstranded]
                if self.library_type:
                    runnable_step_cuffnorm.add_option_long(
                        key='library-type',
                        value=self.library_type)
                # --library-norm-method Method used to normalize library sizes [geometric]
                # --output-format Format for output tables [simple-table]

                # Advanced Options:
                # --compatible-hits-norm count hits compatible with reference RNAs only [TRUE]
                # --total-hits-norm count all hits for normalization [FALSE]
                # --quiet log-friendly quiet processing (no progress bar) [FALSE]
                runnable_step_cuffnorm.add_switch_long(key='quiet')
                # --seed value of random number generator seed [0]
                # --no-update-check do not contact server to check for update availability [FALSE]
                runnable_step_cuffnorm.add_switch_long(key='no-update-check')

                # Undocumented Options:
                # --use-sample-sheet
                runnable_step_cuffnorm.add_switch_long(key='use-sample-sheet')

                # Add the Cuffmerge GTF file as first Cuffnorm argument.
                runnable_step_cuffnorm.arguments.append(file_path_cuffmerge.merged_gtf)

                # Add an abundances annotation file as second Cuffnorm argument.
                # Writing a Cuffnorm abundances TSV file requires an absolute path,
                # because the working directory is not set at the current stage of job submission.
                run_write_annotation(
                    annotation_path=os.path.join(self.genome_directory, file_path_run_cuffnorm.abundances_tsv),
                    annotation_dict=cuffdiff_cuffnorm_abundances_dict)
                runnable_step_cuffnorm.arguments.append(file_path_run_cuffnorm.abundances_tsv)

                # Adjust field names to the annotation before writing to disk.
                monocle_annotation_sheet.to_file_path(adjust_field_names=True)

            if len(self._comparison_dict[comparison_name]) >= 2:
                # Create a Cuffdiff Runnable per comparison, if there are at least two SampleGroup objects.

                file_path_run_cuffdiff = self.get_file_path_run_cuffdiff(comparison_name=comparison_name)

                runnable_run_cuffdiff = self.add_runnable_consecutive(
                    runnable=bsf.procedure.ConsecutiveRunnable(
                        name=self.get_prefix_run_cuffdiff(comparison_name=comparison_name),
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
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

                # Set the environment variable 'LANG' to 'C'.
                runnable_step = bsf.process.RunnableStepSetEnvironment(
                    name='set_environment',
                    key='LANG',
                    value='C')
                runnable_run_cuffdiff.add_runnable_step(runnable_step=runnable_step)

                # Create a new Cuffdiff bsf.process.RunnableStep.

                runnable_step_run_cuffdiff = bsf.process.RunnableStep(
                    name='cuffdiff',
                    program='cuffdiff')
                runnable_run_cuffdiff.add_runnable_step(runnable_step=runnable_step_run_cuffdiff)

                # Set Cuffdiff options.

                # General Options:

                # --output-dir write all output files to this directory [.]
                runnable_step_run_cuffdiff.add_option_long(
                    key='output-dir',
                    value=file_path_run_cuffdiff.output_directory)
                # --labels comma-separated list of condition labels []
                # --FDR False discovery rate used in testing [0.05]
                if self.false_discovery_rate is not None:
                    runnable_step_run_cuffdiff.add_option_long(
                        key='FDR',
                        value=str(self.false_discovery_rate))
                # --mask-file ignore all alignment within transcripts in this file [NULL]
                if self.mask_gtf_path:
                    runnable_step_run_cuffdiff.add_option_long(
                        key='mask-file',
                        value=self.mask_gtf_path)
                # --contrast-file Perform the contrasts specified in this file [NULL]
                # --frag-bias-correct use bias correction - reference fasta required [NULL]
                runnable_step_run_cuffdiff.add_option_long(
                    key='frag-bias-correct',
                    value=self.genome_fasta_path)
                # --multi-read-correct use 'rescue method' for multi-reads [FALSE]
                if self.multi_read_correction:
                    runnable_step_run_cuffdiff.add_switch_long(
                        key='multi-read-correct')
                # --num-threads number of threads used during quantification [1]
                runnable_step_run_cuffdiff.add_option_long(
                    key='num-threads',
                    value=str(stage_run_cuffdiff.threads))
                # --no-diff Don't generate differential analysis files [FALSE]
                # --no-js-tests Don't perform isoform switching tests [FALSE]
                # --time-series treat samples as a time-series [FALSE]
                # --library-type Library prep used for input reads [fr-unstranded]
                if self.library_type:
                    runnable_step_run_cuffdiff.add_option_long(
                        key='library-type',
                        value=self.library_type)
                # --dispersion-method Method used to estimate dispersion models [pooled]
                # --library-norm-method Method used to normalize library sizes [geometric]

                # Advanced Options:
                # --frag-len-mean average fragment length (unpaired reads only) [200]
                # --frag-len-std-dev fragment length std deviation (unpaired reads only) [80]
                # --min-alignment-count minimum number of alignments in a locus for testing [10]
                # --max-mle-iterations maximum iterations allowed for MLE calculation [5000]
                # --compatible-hits-norm count hits compatible with reference RNAs only [TRUE]
                # --total-hits-norm count all hits for normalization [FALSE]
                # --verbose log-friendly verbose processing (no progress bar) [FALSE]
                # --quiet log-friendly quiet processing (no progress bar) [FALSE]
                runnable_step_run_cuffdiff.add_switch_long(key='quiet')
                # --seed value of random number generator seed [0]
                # --no-update-check do not contact server to check for update availability [FALSE]
                # --max-bundle-frags maximum fragments allowed in a bundle before skipping [500000]
                # --num-frag-count-draws Number of fragment generation samples [100]
                # --num-frag-assign-draws Number of fragment assignment samples per generation [50]
                # --max-frag-multihits Maximum number of alignments allowed per fragment [unlim]
                # --min-reps-for-js-test Replicates needed for relative isoform shift testing [3]
                # --no-effective-length-correction No effective length correction [FALSE]
                # --no-length-correction No length correction [FALSE]
                if self.no_length_correction:
                    runnable_step_run_cuffdiff.add_switch_long(key='no-length-correction')
                # --no-update-check do not contact server to check for update availability [FALSE]
                runnable_step_run_cuffdiff.add_switch_long(key='no-update-check')

                # Undocumented Options:
                # --use-sample-sheet
                runnable_step_run_cuffdiff.add_switch_long(key='use-sample-sheet')

                # Add the Cuffmerge GTF file as first Cuffdiff argument.
                runnable_step_run_cuffdiff.arguments.append(file_path_cuffmerge.merged_gtf)

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
                    runnable_step_run_cuffdiff.arguments.append(file_path_run_cuffdiff.abundances_tsv)
                else:
                    # Writing a Cuffdiff alignments TSV file requires an absolute path,
                    # because the working directory is not set at the current stage of job submission.
                    run_write_annotation(
                        annotation_path=os.path.join(self.genome_directory, file_path_run_cuffdiff.alignments_tsv),
                        annotation_dict=cuffdiff_cuffnorm_alignments_dict)
                    runnable_step_run_cuffdiff.arguments.append(file_path_run_cuffdiff.alignments_tsv)

                prefix_process_cuffdiff = self.get_prefix_process_cuffdiff(comparison_name=comparison_name)

                runnable_process_cuffdiff = self.add_runnable_consecutive(
                    runnable=bsf.procedure.ConsecutiveRunnable(
                        name=prefix_process_cuffdiff,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        debug=self.debug))
                executable_process_cuffdiff = self.set_stage_runnable(
                    stage=stage_process_cuffdiff,
                    runnable=runnable_process_cuffdiff)
                executable_process_cuffdiff.dependencies.append(executable_run_cuffdiff.name)

                runnable_step_process_cuffdiff = bsf.process.RunnableStep(
                    name='process_cuffdiff',
                    program='bsf_rnaseq_process_cuffdiff.R')
                runnable_process_cuffdiff.add_runnable_step(runnable_step=runnable_step_process_cuffdiff)

                # Set rnaseq_process_cuffdiff options.
                # Read configuration section [bsf.analyses.rnaseq.Tuxedo.process_cuffdiff]
                self.set_runnable_step_configuration(runnable_step=runnable_step_process_cuffdiff)

                runnable_step_process_cuffdiff.add_option_long(
                    key='comparison-name',
                    value=comparison_name)
                runnable_step_process_cuffdiff.add_option_long(
                    key='gtf-assembly',
                    value=file_path_cuffmerge.merged_gtf)
                runnable_step_process_cuffdiff.add_option_long(
                    key='gtf-reference',
                    value=self.transcriptome_gtf_path)
                runnable_step_process_cuffdiff.add_option_long(
                    key='genome-version',
                    value=self.genome_version)

        # Finally, set dependencies on all other Cuffmerge bsf.process.Executable objects to avoid file contention.
        for prefix_run_cuffmerge in executable_cuffmerge_dict:
            for executable_cuffmerge in executable_cuffmerge_dict.values():
                if prefix_run_cuffmerge != executable_cuffmerge.name:
                    executable_cuffmerge.dependencies.append(prefix_run_cuffmerge)

        return

    def report(self):
        """Create a C{bsf.analyses.rnaseq.Tuxedo} report in HTML format and a UCSC Genome Browser Track Hub.
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

            str_list.append('<h1 id="' + self.prefix + '_analysis">')
            str_list.append(self.project_name + ' ' + self.name)
            str_list.append('</h1>\n')
            str_list.append('\n')

            str_list.extend(self.get_html_genome(genome_version=self.genome_version))
            str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
            str_list.append('\n')

            # TopHat and Cufflinks table.

            str_list.append('<h2 id="transcriptome_browsing">Transcriptome Browsing</h2>\n')
            str_list.append('\n')

            str_list.append('<h3 id="read_alignments">Read Alignments</h3>\n')
            str_list.append('\n')

            if self.aligner == 'hisat2':
                str_list.append('<p id="hisat2">')
                str_list.append('<a href="https://ccb.jhu.edu/software/hisat2/index.shtml">')
                str_list.append('<strong>HISAT2</strong>')
                str_list.append('</a> ')
                str_list.append('aligns RNA-seq reads to a reference genome in order to identify ')
                str_list.append('exon-exon splice junctions. ')
                # str_list.append('<br />\n')
                str_list.append('Please see the ')
                str_list.append('<a href="' + bsf.analyses.hisat.Hisat2.prefix + '_report.html">')
                str_list.append(self.project_name + ' ' + bsf.analyses.hisat.Hisat2.name)
                str_list.append('</a> report for quality plots and ')
                str_list.append('a link to alignment visualisation in the UCSC Genome Browser.\n')
                str_list.append('</p>\n')
                str_list.append('\n')
            elif self.aligner == 'star':
                str_list.append('<p id="star_aligner">')
                str_list.append('<a href="https://github.com/alexdobin/STAR">STAR</a> ')
                str_list.append('aligns RNA-seq reads to a reference genome in order to identify ')
                str_list.append('exon-exon splice junctions. ')
                # str_list.append('<br />\n')
                str_list.append('Please see the ')
                str_list.append('<a href="' + bsf.analyses.star_aligner.StarAligner.prefix + '_report.html">')
                str_list.append(self.project_name + ' ' + bsf.analyses.star_aligner.StarAligner.name)
                str_list.append('</a> report for quality plots and ')
                str_list.append('a link to alignment visualisation in the UCSC Genome Browser.\n')
                str_list.append('</p>\n')
                str_list.append('\n')
            else:
                str_list.append('<p id ="tophat">')
                # http://tophat.cbcb.umd.edu/manual.html
                str_list.append('<a href="http://ccb.jhu.edu/software/tophat/index.shtml"><strong>TopHat</strong></a> ')
                str_list.append('aligns RNA-seq reads to a genome in order to identify ')
                str_list.append('exon-exon splice junctions. It is built on the ultra fast ')
                str_list.append('short read mapping program ')
                str_list.append('<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">')
                str_list.append('<strong>Bowtie 2</strong>')
                str_list.append('</a>.')
                str_list.append('</p>\n')
                str_list.append('\n')

                str_list.append('<p id="track_hub">')
                str_list.append('View TopHat <strong>read alignments</strong> tracks for each sample\n')
                str_list.append('in their genomic context via the project-specific ')
                str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
                str_list.append('.')
                str_list.append('</p>\n')
                str_list.append('\n')

                str_list.append('<p>')
                str_list.append('<a href="rnaseq_tophat_alignment_summary.pdf">')
                str_list.append('<img ')
                str_list.append('alt="TopHat Alignment Summary" ')
                str_list.append('id="tophat_alignment_summary_img" ')
                str_list.append('src="rnaseq_tophat_alignment_summary.png" ')
                str_list.append('height="80" ')
                str_list.append('width="80" ')
                str_list.append('/>')
                str_list.append('</a> ')
                str_list.append('Alignment summary statistics <a href="rnaseq_tophat_alignment_summary.tsv">TSV</a>')
                str_list.append('</p>\n')

                str_list.append('<h3 id="alignment_events">Splice Junctions, Insertions and Deletions</h3>\n')
                str_list.append('\n')

                str_list.append('<p>')
                str_list.append('TopHat reports <strong>splice junctions</strong> on the basis of RNA-seq\n')
                str_list.append('read alignments in UCSC BED track format.\n')
                str_list.append('Each junction consists of two connected BED blocks,\n')
                str_list.append('where each block is as long as the maximal overhang\n')
                str_list.append('of any read spanning the junction. The score is\n')
                str_list.append('the number of alignments spanning the junction.\n')
                str_list.append('UCSC BED tracks of <strong>insertions</strong> and\n')
                str_list.append('<strong>deletions</strong> are also reported by TopHat.')
                str_list.append('</p>\n')

                str_list.append('<p>')
                str_list.append('View the corresponding TopHat tracks for junctions, deletions and insertions\n')
                str_list.append('for each sample in their genomic context via the project-specific\n')
                str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
                str_list.append('.')
                str_list.append('</p>\n')
                str_list.append('\n')

                # str_list.append('<p>')
                # str_list.append('Follow the links below to attach\n')
                # str_list.append('Tophat junction, deletion and insertion annotation to the\n')
                # str_list.append('UCSC Genome Browser. Since each file needs transferring to\n')
                # str_list.append('the UCSC site, subsequent pages will take some time to load.')
                # str_list.append('</p>\n')

            str_list.append('<h2 id="raw_gene_expression_profiles">Raw Gene Expression Profiles</h2>\n')
            str_list.append('\n')

            str_list.append('<p id="cufflinks">')
            # http://cufflinks.cbcb.umd.edu/howitworks.html
            str_list.append('<a href="http://cole-trapnell-lab.github.io/cufflinks/"><strong>Cufflinks</strong></a>\n')
            str_list.append('assembles aligned RNA-seq reads into transcripts,\n')
            str_list.append('estimates their abundances, and tests for differential\n')
            str_list.append('expression and regulation transcriptome-wide.\n')
            str_list.append('It accepts aligned RNA-seq reads and assembles the alignments into a parsimonious set\n')
            str_list.append('of transcripts. Cufflinks then estimates the relative abundances of these transcripts\n')
            str_list.append('based on how many reads support each one, taking into account biases in library\n')
            str_list.append('preparation protocols.')
            str_list.append('</p>\n')

            str_list.append('<p>')
            str_list.append('<strong>Please note:</strong>\n')
            str_list.append('Since Cufflinks estimates transcriptome expression profiles on a per sample ')
            str_list.append('(or replicate) basis, the resulting raw FPKM values cannot directly be compared ')
            str_list.append('between samples.\n')
            str_list.append('Please see the <a href="#differential_expression_table">Differential Expression</a> ')
            str_list.append('section below for tables of normalised FPKM values that allow for direct ')
            str_list.append('sample comparisons.')
            str_list.append('</p>\n')

            str_list.append('<p>')
            str_list.append('The Cufflinks <strong>assembled transcripts</strong> can be attached to the \n')
            str_list.append('UCSC Genome Browser, by following the &laquo;Transcript Assembly&raquo; links\n')
            str_list.append('below.\n')
            str_list.append('The isoforms.fpkm_tracking and genes.fpkm_tracking files\n')
            str_list.append('contain the estimated isoform or gene expression values in the generic\n')
            # http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
            str_list.append('<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' +
                            'fpkm-tracking-format">FPKM Tracking format</a>.\n')
            str_list.append('The isoforms.count_tracking and genes.count_tracking files\n')
            str_list.append('contain the scaled isoform or gene count values in the generic\n')
            str_list.append('<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' +
                            'count-tracking-format">Count Tracking format</a>. ')
            str_list.append('Please see a more detailed description of\n')
            # http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output
            str_list.append('<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' +
                            'output-formats-used-in-the-cufflinks-suite">Cufflinks output</a> format.')
            str_list.append('</p>\n')

            str_list.append('<table id="gene_expression_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Sample</th>\n')
            str_list.append('<th>Assembled Transcripts</th>\n')
            str_list.append('<th>Gene FPKM</th>\n')
            str_list.append('<th>Transcript FPKM</th>\n')
            str_list.append('<th>Genes (Symbols)</th>\n')
            str_list.append('<th>Isoforms (Symbols)</th>\n')
            str_list.append('<th>Aligned BAM</th>\n')
            str_list.append('<th>Aligned BAI</th>\n')
            if self.aligner == 'hisat2':
                str_list.append('<th>Aligned MD5</th>\n')
            elif self.aligner == 'star':
                str_list.append('<th>Aligned MD5</th>\n')
            else:
                str_list.append('<th>Unaligned BAM</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for sample in self.sample_list:
                if self.debug > 0:
                    print(self, 'Sample name:', sample.name)
                    sys.stdout.writelines(sample.trace(level=1))

                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                if self.aligner == 'hisat2' or self.aligner == 'star':
                    if self.aligner == 'hisat2':
                        file_path_aligner_sample = bsf.analyses.hisat.Hisat2.get_file_path_sample(
                            sample_name=sample.name)
                    elif self.aligner == 'star':
                        file_path_aligner_sample = bsf.analyses.star_aligner.StarAligner.get_file_path_sample(
                            sample_name=sample.name)
                    else:
                        raise Exception()

                    path_prefix = 'rnaseq_cufflinks_' + sample.name

                    str_list.append('<tr>\n')

                    # Sample
                    str_list.append('<td class="left">' + sample.name + '</td>\n')
                    # Assembled Transcripts
                    str_list.append('<td class="center">')
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='transcripts.gtf',
                        text='Transcript Assembly'))
                    str_list.append('</td>\n')

                    # Gene FPKM
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/genes.fpkm_tracking">Genes FPKM</a>')
                    str_list.append('</td>\n')

                    # Transcript FPKM
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/isoforms.fpkm_tracking">Isoforms FPKM</a>')
                    str_list.append('</td>\n')

                    # Genes (Symbols)
                    str_list.append('<td class="center">')
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='genes_fpkm_tracking.tsv',
                        text='Genes (Symbols)'))
                    str_list.append('</td>\n')

                    # Isoforms (Symbols)
                    str_list.append('<td class="center">')
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='isoforms_fpkm_tracking.tsv',
                        text='Isoforms (Symbols)'))
                    str_list.append('</td>\n')

                    # Aligned BAM file
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + file_path_aligner_sample.sample_bam + '">')
                    str_list.append('<abbr title="Binary Alignment/Map">BAM</abbr>')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Aligned BAI file
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + file_path_aligner_sample.sample_bai + '">')
                    str_list.append('<abbr title="Binary Alignment/Map Index">BAI</abbr>')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Aligned BAI file
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + file_path_aligner_sample.sample_md5 + '">')
                    str_list.append('<abbr title="Message Digest 5 Checksum">MD5</abbr>')
                    str_list.append('</a>')
                    str_list.append('</td>\n')
                    str_list.append('</tr>\n')
                else:
                    # Cufflinks produces genes.fpkm_tracking, isoforms.fpkm_tracking,
                    # skipped.gtf and transcripts.gtf.

                    path_prefix = 'rnaseq_cufflinks_' + sample.name

                    str_list.append('<tr>\n')

                    # Sample
                    str_list.append('<td class="left">' + sample.name + '</td>\n')
                    # Assembled Transcripts
                    str_list.append('<td class="center">')
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='transcripts.gtf',
                        text='Transcript Assembly'))
                    str_list.append('</td>\n')

                    # Gene FPKM
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/genes.fpkm_tracking">Genes FPKM</a>')
                    str_list.append('</td>\n')

                    # Transcript FPKM
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/isoforms.fpkm_tracking">Isoforms FPKM</a>')
                    str_list.append('</td>\n')

                    # Genes (Symbols)
                    str_list.append('<td class="center">')
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='genes_fpkm_tracking.tsv',
                        text='Genes (Symbols)'))
                    str_list.append('</td>\n')

                    # Isoforms (Symbols)
                    str_list.append('<td class="center">')
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='isoforms_fpkm_tracking.tsv',
                        text='Isoforms (Symbols)'))
                    str_list.append('</td>\n')

                    # TODO: The aligned BAM and BAI files and the unaligned BAM file are currently non standard.
                    # The files have a 'rnaseq_tophat_' prefix, but are in the 'rnaseq_cufflinks_' directory.
                    # This will be resolved when the process_tophat step gets re-engineered.

                    # Aligned BAM file
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/rnaseq_tophat_' + sample.name +
                                    '_accepted_hits.bam">Aligned BAM</a>')
                    str_list.append('</td>\n')

                    # Aligned BAI file
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/rnaseq_tophat_' + sample.name +
                                    '_accepted_hits.bam.bai">Aligned BAI</a>')
                    str_list.append('</td>\n')

                    # Unaligned BAM file
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + path_prefix + '/rnaseq_tophat_' + sample.name +
                                    '_unmapped.bam">Unaligned BAM</a>')
                    str_list.append('</td>\n')

                    str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Cuffdiff produces cds_exp.diff, gene_exp.diff, isoform_exp.diff
            # promoters.diff, splicing.diff and tss_group_exp.diff amongst many others.

            str_list.append('<h2 id="differential_expression">Differential Expression</h2>\n')
            str_list.append('\n')

            str_list.append('<p id="cuffdiff">')
            # http://cufflinks.cbcb.umd.edu/howitworks.html#diff
            str_list.append('<a href="http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html">')
            str_list.append('<strong>Cuffdiff</strong>')
            str_list.append('</a>\n')
            str_list.append('finds significant changes in transcript\n')
            str_list.append('expression, splicing, and promoter use.')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<p>')
            str_list.append('A comparison involves one or more sample groups with one or more replicates each.\n')
            str_list.append('Cuffdiff normalises all replicates and performs an all-against-all sample groups ')
            str_list.append('comparison.\nTypically, a single comparison including all sample groups is set up, ')
            str_list.append('but sometimes, technical constraints (i.e. memory requirements) require setting up ')
            str_list.append('several comparisons, each with fewer sample groups, in parallel.\n')
            str_list.append('Since Cuffdiff normalises all samples within a comparison, FPKM values of all samples ')
            str_list.append('in a particular comparison can be directly compared.')
            str_list.append('</p>\n')

            str_list.append('<h3 id="all_genes">All Genes</h3>\n')

            str_list.append('<p>')
            str_list.append('The links in the &laquo;Genes&raquo; column provide pairwise slices of the ')
            str_list.append('all-against-all comparisons in tab-separated value (TSV) format, ')
            str_list.append('which could be imported into spreadsheet software for post-filtering.\n')
            str_list.append('Each table contains information about a pair of sample groups annotated in ')
            str_list.append('&laquo;sample_1&raquo; and &laquo;sample_2&raquo;, their corresponding FPKM values in ')
            str_list.append('&laquo;value_1&raquo; and &laquo;value_2&raquo;, log2-fold changes of ')
            str_list.append('&laquo;sample_2&raquo; over &laquo;sample_1&raquo;, as well as p-value and ')
            str_list.append('multiple testing-corrected (Benjamini-Hochberg) q-values.\n')
            str_list.append('By default, the false discovery rate (FDR) threshold is set at 0.05, ')
            str_list.append('which is also the basis for the assignment of &laquo;yes&raquo; or &laquo;no&raquo; ')
            str_list.append('in the &laquo;significant&raquo; column.\n')
            str_list.append('Depending on the biological system, the threshold could be dynamically raised or ')
            str_list.append('lowered in a post-filtering approach.\n')
            str_list.append('Please see a more detailed description of the ')
            str_list.append('<a href="http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html' +
                            '#differential-expression-tests">Cuffdiff differential expression tests</a> format.</p>\n')
            str_list.append('\n')
            str_list.append('<p>Since some spreadsheet programmes have no concept of positive or negative infinity ')
            str_list.append('values (+Inf, -Inf), the &laquo;Genes&raquo; and &laquo;Isoforms&raquo; tables contain ')
            str_list.append('three additional ranking columns for ')
            str_list.append('(1) the effect size (rank_log2_fold_change), ')
            str_list.append('(2) the absolute level (rank_value) and ')
            str_list.append('(3) the statistical significance (rank_q_value).\n')
            str_list.append('Ranking in these columns is correct and could be combined to down-vote huge ')
            str_list.append('log2-fold changes that are the result of tiny absolute values and therefore ')
            str_list.append('biologically meaningless.\n')
            str_list.append('A maximum rank of combinations of columns could work.')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="differential_expression_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Samples</th>\n')
            str_list.append('<th>Replicates</th>\n')
            str_list.append('<th>Coding Sequences</th>\n')
            str_list.append('<th>Genes</th>\n')
            str_list.append('<th>Isoforms</th>\n')
            str_list.append('<th>Promoters</th>\n')
            str_list.append('<th>Splicing</th>\n')
            str_list.append('<th>Transcription Start Sites</th>\n')
            str_list.append('<th>Gene FPKM Replicates</th>\n')
            str_list.append('<th>Gene Count Replicates</th>\n')
            str_list.append('<th>Isoform FPKM Replicates</th>\n')
            str_list.append('<th>Isoform Count Replicates</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            # Since the sorted list of comparison names is used several times below, sort it only once.
            comparison_name_list = sorted(self._comparison_dict)

            for comparison_name in comparison_name_list:
                path_prefix = 'rnaseq_process_cuffdiff_' + comparison_name

                # Link to comparison-specific symbolic links in the directory after cummeRbund processing.

                str_list.append('<tr>\n')
                # Comparison
                str_list.append('<td class="left">' + comparison_name + '</td>\n')
                # Samples
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='samples.tsv',
                    text='Samples'))
                str_list.append('</td>\n')
                # Replicates
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='replicates.tsv',
                    text='Replicates'))
                str_list.append('</td>\n')
                # Coding Sequences
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='cds_exp_diff.tsv',
                    text='Coding Sequences'))
                str_list.append('</td>\n')
                # Genes
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_exp_diff.tsv',
                    text='<strong>Genes</strong>'))
                str_list.append('</td>\n')
                # Isoforms
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_exp_diff.tsv',
                    text='Isoforms'))
                str_list.append('</td>\n')
                # Promoters
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='promoters_diff.tsv',
                    text='Promoters'))
                str_list.append('</td>\n')
                # Splicing
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='splicing_diff.tsv',
                    text='Splicing'))
                str_list.append('</td>\n')
                # Transcription Start Sites
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='tss_group_exp_diff.tsv',
                    text='Transcription Start Sites'))
                str_list.append('</td>\n')
                # Gene FPKM Replicates
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_fpkm_replicates.tsv',
                    text='Gene FPKM Replicates'))
                str_list.append('</td>\n')
                # Gene Count Replicates
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_counts_replicates.tsv',
                    text='Gene Count Replicates'))
                str_list.append('</td>\n')
                # Isoform FPKM Replicates
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_fpkm_replicates.tsv',
                    text='Isoform FPKM Replicates'))
                str_list.append('</td>\n')
                # Isoform Count Replicates
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_counts_replicates.tsv',
                    text='Isoform Count Replicates'))
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

                # Read sample pair information if available.

                sample_pair_path = os.path.join(
                    self.genome_directory,
                    path_prefix,
                    '_'.join((path_prefix, 'sample_pairs.tsv')))

                if os.path.exists(sample_pair_path):
                    sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                    for row_dict in sample_pair_sheet.row_dicts:
                        str_list.append('<tr>\n')
                        # Comparison
                        str_list.append('<td></td>\n')
                        # Samples
                        # Replicates
                        # Coding Sequences
                        str_list.append('<td class="left" colspan="3">')
                        str_list.append('<strong>' + row_dict['V1'] + '</strong>')
                        str_list.append(' versus ')
                        str_list.append('<strong>' + row_dict['V2'] + '</strong>')
                        str_list.append('</td>\n')
                        # Genes
                        str_list.append('<td class="center">')
                        str_list.append(self.get_html_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_diff.tsv')),
                            text='<strong>Genes</strong>'))
                        str_list.append('</td>\n')
                        # Isoforms
                        str_list.append('<td class="center">')
                        str_list.append(self.get_html_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'isoforms_diff.tsv')),
                            text='Isoforms'))
                        str_list.append('</td>\n')
                        # Promoters
                        # Splicing
                        # Transcription Start Sites
                        # Gene FPKM Replicates
                        # Gene Count Replicates
                        # Isoform FPKM Replicates
                        # Isoform Count Replicates
                        str_list.append('<td class="left" colspan="7"></td>\n')
                        str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            str_list.append('<h3 id="significant_genes">Significant Genes</h3>\n')

            str_list.append('<table id="significant_genes_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Genes</th>\n')
            str_list.append('<th>Isoforms</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in comparison_name_list:
                path_prefix = 'rnaseq_process_cuffdiff_' + comparison_name

                str_list.append('<tr>\n')
                # Comparison
                str_list.append('<td class="left">' + comparison_name + '</td>\n')
                # Genes
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_significance_matrix.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_significance_matrix.png',
                        text='Significance Matrix Plot - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')
                # Isoforms
                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_significance_matrix.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_significance_matrix.png',
                        text='Significance Matrix Plot - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')

            # Show cummeRbund quality plots.

            str_list.append('<h2 id="quality_plots">Quality Plots</h2>\n')
            str_list.append('\n')

            str_list.append('<p>')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="quality_plots_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Dispersion Plot - Genes</th>\n')
            str_list.append('<th>Dispersion Plot - Isoforms</th>\n')
            str_list.append('<th>Squared Coefficient of Variation - Genes</th>\n')
            str_list.append('<th>Squared Coefficient of Variation - Isoforms</th>\n')
            str_list.append('<th>Density Plot without Replicates - Genes</th>\n')
            str_list.append('<th>Density Plot with Replicates - Genes</th>\n')
            str_list.append('<th>Density Plot without Replicates - Isoforms</th>\n')
            str_list.append('<th>Density Plot with Replicates - Isoforms</th>\n')
            str_list.append('<th>Box Plot without Replicates - Genes</th>\n')
            str_list.append('<th>Box Plot with Replicates - Genes</th>\n')
            str_list.append('<th>Box Plot without Replicates - Isoforms</th>\n')
            str_list.append('<th>Box Plot with Replicates - Isoforms</th>\n')
            str_list.append('<th>Scatter Matrix Plot - Genes</th>\n')
            str_list.append('<th>Scatter Matrix Plot - Isoforms</th>\n')
            str_list.append('<th>Dendrogram Plot</th>\n')
            str_list.append('<th>Volcano Matrix Plot - Genes</th>\n')
            str_list.append('<th>Multidimensional Scaling Plot - Genes</th>\n')
            str_list.append('<th>Principal Component Analysis Plot - Genes</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in comparison_name_list:
                path_prefix = 'rnaseq_process_cuffdiff_' + comparison_name

                str_list.append('<tr>\n')
                # Comparison
                str_list.append('<td class="left">' + comparison_name + '</td>\n')

                # Dispersion Plots for Genes and Isoforms

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_dispersion.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_dispersion.png',
                        text='Dispersion Plot - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_dispersion.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_dispersion.png',
                        text='Dispersion Plot - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms

                str_list.append('<td class="center">')
                if os.path.exists(os.path.join(self.genome_directory, path_prefix, path_prefix + '_genes_scv.png')):
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='genes_scv.pdf',
                        text=self.get_html_image(
                            prefix=path_prefix,
                            suffix='genes_scv.png',
                            text='Squared Coefficient of Variation (SCV) - Genes - ' + comparison_name,
                            height='80',
                            width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                if os.path.exists(os.path.join(self.genome_directory, path_prefix, path_prefix + '_isoforms_scv.png')):
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='isoforms_scv.pdf',
                        text=self.get_html_image(
                            prefix=path_prefix,
                            suffix='isoforms_scv.png',
                            text='Squared Coefficient of Variation (SCV) - Isoforms - ' + comparison_name,
                            height='80',
                            width='80')))
                str_list.append('</td>\n')

                # Density Plots for Genes without and with Replicates

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_density_wo_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_density_wo_replicates.png',
                        text='Density Plot without Replicates - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_density_w_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_density_w_replicates.png',
                        text='Density Plot with Replicates - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Density Plots for Isoforms without and with Replicates

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_density_wo_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_density_wo_replicates.png',
                        text='Density Plot without Replicates - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_density_w_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_density_w_replicates.png',
                        text='Density Plot with Replicates - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Box Plots for Genes without and with Replicates

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_box_wo_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_box_wo_replicates.png',
                        text='Box Plot without Replicates - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_box_w_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_box_w_replicates.png',
                        text='Box Plot with Replicates - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Box Plots for Isoforms with and without Replicates

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_box_wo_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_box_wo_replicates.png',
                        text='Box Plot without Replicates - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_box_w_replicates.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_box_w_replicates.png',
                        text='Box Plot with Replicates - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Scatter Matrix Plot for Genes and Isoforms

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_scatter_matrix.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_scatter_matrix.png',
                        text='Scatter Matrix Plot - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='isoforms_scatter_matrix.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='isoforms_scatter_matrix.png',
                        text='Scatter Matrix Plot - Isoforms - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Dendrogram Plot for Genes

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_dendrogram.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_dendrogram.png',
                        text='Dendrogram Plot - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Volcano Matrix Plot for Genes

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_volcano_matrix.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_volcano_matrix.png',
                        text='Volcano Matrix Plot - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # Multidimensional Scaling Plot for Genes

                str_list.append('<td class="center">')
                if os.path.exists(os.path.join(self.genome_directory, path_prefix, path_prefix + '_genes_mds.png')):
                    str_list.append(self.get_html_anchor(
                        prefix=path_prefix,
                        suffix='genes_mds.pdf',
                        text=self.get_html_image(
                            prefix=path_prefix,
                            suffix='genes_mds.png',
                            text='Multidimensional Scaling Plot - Genes - ' + comparison_name,
                            height='80',
                            width='80')))
                str_list.append('</td>\n')

                # Principal Component Analysis Plot for Genes

                str_list.append('<td class="center">')
                str_list.append(self.get_html_anchor(
                    prefix=path_prefix,
                    suffix='genes_pca.pdf',
                    text=self.get_html_image(
                        prefix=path_prefix,
                        suffix='genes_pca.png',
                        text='Principal Component Analysis Plot - Genes - ' + comparison_name,
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                str_list.append('</tr>\n')

                # Read sample pair information if available.

                sample_pair_path = os.path.join(
                    self.genome_directory,
                    path_prefix,
                    '_'.join((path_prefix, 'sample_pairs.tsv')))

                if os.path.exists(sample_pair_path):
                    sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                    for row_dict in sample_pair_sheet.row_dicts:
                        str_list.append('<tr>\n')

                        # Comparison
                        str_list.append('<td class="left"></td>\n')
                        str_list.append('<td class="left" colspan="10">')
                        str_list.append('<strong>' + row_dict['V1'] + '</strong>')
                        str_list.append(' versus ')
                        str_list.append('<strong>' + row_dict['V2'] + '</strong>')
                        str_list.append('</td>\n')

                        str_list.append('<td class="center">')
                        str_list.append(self.get_html_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_scatter.pdf')),
                            text=self.get_html_image(
                                prefix=path_prefix,
                                suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_scatter.png')),
                                text='Scatter Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2'],
                                height='80',
                                width='80')))
                        str_list.append('</td>\n')

                        str_list.append('<td class="center"></td>\n')

                        str_list.append('<td class="center">')
                        str_list.append(self.get_html_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_ma.pdf')),
                            text=self.get_html_image(
                                prefix=path_prefix,
                                suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_ma.png')),
                                text='M vs A Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2'],
                                height='80',
                                width='80')))
                        str_list.append('</td>\n')

                        str_list.append('<td class="center">')
                        str_list.append(self.get_html_anchor(
                            prefix=path_prefix,
                            suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_volcano.pdf')),
                            text=self.get_html_image(
                                prefix=path_prefix,
                                suffix='_'.join((row_dict['V1'], row_dict['V2'], 'genes_volcano.png')),
                                text='Volcano Plot on genes ' + row_dict['V1'] + ' versus ' + row_dict['V2'],
                                height='80',
                                width='80')))
                        str_list.append('</td>\n')

                        str_list.append('<td class="center" colspan="4"></td>\n')

                        str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

            if self.debug > 0:
                print('Report list:')
                print(repr(str_list))

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
            str_list.append('longLabel TopHat RNA-seq read alignments\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignments\n')
            str_list.append('\n')

            str_list.append('track Assemblies\n')
            str_list.append('shortLabel Assemblies\n')
            str_list.append('longLabel Cuffmerge transcript structures\n')
            str_list.append('visibility full\n')
            str_list.append('superTrack on\n')
            str_list.append('group assemblies\n')
            str_list.append('\n')

            str_list.append('track Coverage\n')
            str_list.append('shortLabel Coverage\n')
            str_list.append('longLabel TopHat RNA-seq alignment coverage\n')
            str_list.append('visibility full\n')
            str_list.append('superTrack on\n')
            str_list.append('group coverage\n')
            str_list.append('\n')

            str_list.append('track Deletions\n')
            str_list.append('shortLabel Deletions\n')
            str_list.append('longLabel TopHat RNA-seq deletions\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignments\n')
            str_list.append('\n')

            str_list.append('track Insertions\n')
            str_list.append('shortLabel Insertions\n')
            str_list.append('longLabel TopHat RNA-seq insertions\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignments\n')
            str_list.append('\n')

            str_list.append('track Junctions\n')
            str_list.append('shortLabel Junctions\n')
            str_list.append('longLabel TopHat RNA-seq splice junctions\n')
            str_list.append('visibility show\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignments\n')
            str_list.append('\n')

            str_list.append('track Transcripts\n')
            str_list.append('shortLabel Transcripts\n')
            str_list.append('longLabel Cufflinks transcript structures\n')
            str_list.append('visibility show\n')
            str_list.append('superTrack on\n')
            str_list.append('group transcripts\n')
            str_list.append('\n')

            # Sample-specific tracks

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                if self.aligner == 'hisat2':
                    pass
                elif self.aligner == 'star':
                    pass
                else:
                    #
                    # Add a trackDB entry for each Tophat accepted_hits.bam file.
                    #
                    # Common settings
                    str_list.append('track ' + sample.name + '_alignments\n')
                    str_list.append('type bam\n')
                    str_list.append('shortLabel ' + sample.name + '_alignments\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq read alignments\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/accepted_hits.bam\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility dense\n')

                    # Common optional settings
                    str_list.append('color 0,0,0\n')

                    # bam/cram - Compressed Sequence Alignment track settings
                    # None

                    # Composite track settings
                    str_list.append('parent Alignments\n')
                    str_list.append('\n')

                    #
                    # Add a trackDB entry for each Tophat accepted_hits.bw file.
                    #
                    # Common settings
                    str_list.append('track ' + sample.name + '_coverage\n')
                    # TODO: The bigWig type must declare the expected signal range.
                    # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                    str_list.append('type bigWig\n')
                    str_list.append('shortLabel ' + sample.name + '_coverage\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq alignment coverage\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/accepted_hits.bw\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility full\n')

                    # Common optional settings
                    str_list.append('color 0,0,0\n')

                    # bigWig - Signal graphing track settings
                    str_list.append('alwaysZero on\n')
                    str_list.append('autoScale on\n')
                    str_list.append('graphTypeDefault bar\n')
                    str_list.append('maxHeightPixels 100:60:20\n')
                    # str_list.append('maxWindowToQuery 10000000\n')
                    # str_list.append('smoothingWindow 5\n')
                    # str_list.append('transformFunc NONE\n')
                    # str_list.append('viewLimits 0:45\n')
                    # str_list.append('viewLimitsMax 0:50\n')
                    # str_list.append('windowingFunction maximum\n')
                    # str_list.append('yLineMark <#>\n')
                    # str_list.append('yLineOnOff on \n')
                    # str_list.append('gridDefault on\n')

                    # Composite track settings
                    str_list.append('parent Coverage\n')
                    str_list.append('centerLabelsDense off\n')
                    str_list.append('\n')

                    #
                    # Add a trackDB entry for each Tophat deletions.bb file.
                    #
                    # Common settings
                    str_list.append('track ' + sample.name + '_deletions\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_deletions\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq deletions\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/deletions.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')

                    # Common optional settings
                    str_list.append('color 0,0,0\n')

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list.append('parent Deletions\n')
                    str_list.append('\n')

                    #
                    # Add a trackDB entry for each Tophat insertions.bb file.
                    #
                    # Common settings
                    str_list.append('track insertions_' + sample.name + '\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_insertions\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq insertions\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/insertions.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')

                    # Common optional settings
                    str_list.append('color 0,0,0\n')

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list.append('parent Insertions\n')
                    str_list.append('\n')

                    #
                    # Add a trackDB entry for each Tophat junctions.bb file.
                    #
                    # Common settings
                    str_list.append('track ' + sample.name + '_junctions\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_junctions\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq splice junctions\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/junctions.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility pack\n')

                    # Common optional settings
                    str_list.append('color 0,0,0\n')

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list.append('parent Junctions\n')
                    str_list.append('\n')

                    #
                    # Add a trackDB entry for each Tophat transcripts.bb file.
                    #
                    # Common settings
                    str_list.append('track ' + sample.name + '_transcripts\n')
                    str_list.append('type bigGenePred\n')
                    str_list.append('shortLabel ' + sample.name + '_transcripts\n')
                    str_list.append('longLabel ' + sample.name + ' Cufflinks transcript assembly\n')
                    str_list.append('bigDataUrl rnaseq_cufflinks_' + sample.name + '/transcripts.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')

                    # Common optional settings
                    str_list.append('color 0,0,0\n')

                    # bigGenePred - Gene Annotations settings
                    # None

                    # Composite track settings
                    str_list.append('parent Transcripts\n')
                    str_list.append('\n')

            # Comparison-specific tracks

            for comparison_name in sorted(self._comparison_dict):
                #
                # Add a trackDB entry for each Cuffmerge merged.bb file.
                #
                # Common settings
                str_list.append('track ' + comparison_name + '_assembly\n')
                str_list.append('type bigGenePred\n')
                str_list.append('shortLabel ' + comparison_name + '_assembly\n')
                str_list.append('longLabel ' + comparison_name + ' Cufflinks transcript assembly\n')
                str_list.append('bigDataUrl rnaseq_cuffmerge_' + comparison_name + '/merged.bb\n')
                # str_list.append('html ...\n')
                str_list.append('visibility pack\n')

                # Common optional settings
                str_list.append('color 0,0,0\n')

                # bigGenePred - Gene Annotations settings
                # None

                # Composite track settings
                str_list.append('parent Assemblies\n')
                str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return


class FilePathDESeq(bsf.procedure.FilePath):
    """The C{bsf.analyses.rnaseq.FilePathDESeq} models files in a comparison-specific DESeq directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rnaseq.FilePathDESeq} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathDESeq, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class DESeq(bsf.analysis.Analysis):
    """DESeq RNASeq C{bsf.analysis.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode | None
    @ivar contrast_path: Contrast file path
    @type contrast_path: str | unicode | None
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode | None
    @ivar transcriptome_gtf_path: Reference transcriptome GTF file path
    @type transcriptome_gtf_path: str | unicode | None
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | unicode | None
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str | None
    """

    name = 'DESeq RNA-seq Analysis'
    prefix = '_'.join(('rnaseq', 'deseq'))

    @classmethod
    def get_stage_name_analysis(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'analysis'))

    @classmethod
    def get_stage_name_results(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'results'))

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
            transcriptome_version=None):
        """Initialise a C{bsf.analyses.rnaseq.DESeq} object.

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
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode | None
        @param contrast_path: Contrast file path
        @type contrast_path: str | unicode | None
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode | None
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode | None
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode | None
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

        self.replicate_grouping = replicate_grouping
        self.comparison_path = comparison_path
        self.contrast_path = contrast_path
        self.genome_fasta_path = genome_fasta_path
        self.transcriptome_gtf_path = transcriptome_gtf_path
        self.transcriptome_index_path = transcriptome_index_path
        self.transcriptome_version = transcriptome_version

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

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.rnaseq.DESeq} analysis.

        @return:
        @rtype:
        """

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # DESeq requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A ' + self.name + " requires a 'transcriptome_version' configuration option.")

        if not self.genome_version:
            self.genome_version = bsf.standards.Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a valid 'transcriptome_version' configuration option.")

        # Get the annotation sheets before calling the run() method of the Analysis super-class.
        # If file paths were not provided, try to find them in the current directory.
        # The complication is that either the Tuxedo.prefix or the DESeq.prefix could be used.

        # Get the sample annotation sheet.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(
                prefix_list=[DESeq.prefix, Tuxedo.prefix],
                suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # Get the design annotation sheet.

        if self.comparison_path:
            self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)
            if not os.path.exists(self.comparison_path):
                raise Exception(
                    'Comparison (design) annotation file ' + repr(self.comparison_path) + ' does not exist.')
        else:
            self.comparison_path = self.get_annotation_file(
                prefix_list=[DESeq.prefix, Tuxedo.prefix],
                suffix='designs.csv')
            if not self.comparison_path:
                raise Exception('No suitable comparison (design) annotation file in the current working directory.')

        # Get the contrast annotation sheet.

        if self.contrast_path:
            self.contrast_path = self.configuration.get_absolute_path(file_path=self.contrast_path)
            if not os.path.exists(self.contrast_path):
                raise Exception('Contrast annotation file ' + repr(self.contrast_path) + ' does not exist.')
        else:
            self.contrast_path = self.get_annotation_file(
                prefix_list=[DESeq.prefix, Tuxedo.prefix],
                suffix='contrasts.csv')
            if not self.contrast_path:
                raise Exception('No suitable contrast annotation file in the current working directory.')

        super(DESeq, self).run()

        # For DESeq, all samples need adding to the Analysis regardless.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        # Read the designs (comparison) file.

        design_sheet = bsf.annotation.AnnotationSheet.from_file_path(
            file_path=self.comparison_path,
            file_type='excel',
            name='DESeq Design Table')

        # Read the contrasts file.

        contrast_sheet = bsf.annotation.AnnotationSheet.from_file_path(
            file_path=self.contrast_path,
            file_type='excel',
            name='DESeq Contrast Table')

        # TODO: Adjust by introducing a new class RNASeqComparisonSheet(AnnotationSheet) in this module?
        for design_row_dict in design_sheet.row_dicts:
            design_name = design_row_dict['design']

            prefix = '_'.join((self.prefix, design_name))

            comparison_directory = os.path.join(self.genome_directory, prefix)

            if not os.path.isdir(comparison_directory):
                try:
                    os.makedirs(comparison_directory)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            annotation_sheet = bsf.annotation.AnnotationSheet(
                file_path=os.path.join(comparison_directory, prefix + '_samples.tsv'),
                file_type='excel-tab',
                # file_type='excel',
                name='DESeq Sample Annotation',
                header=True)

            # Sort the Python list of Sample objects by Sample.name.

            self.sample_list.sort(key=lambda item: item.name)

            for sample in self.sample_list:
                if self.debug > 0:
                    print(self, 'Sample name:', sample.name)
                    sys.stdout.writelines(sample.trace(level=1))

                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                file_path_sample = bsf.analyses.star_aligner.StarAligner.get_file_path_sample(sample_name=sample.name)

                row_dict = {
                    'bam_path': file_path_sample.sample_bam,
                    'bai_path': file_path_sample.sample_bai,
                }
                # If the sample annotation sheet contains a "DESeq sample" variable,
                # technical replicates need collapsing. Set the original bsf.ngs.Sample.name as
                # 'run', the DESeq sample name will be filled in from the annotation dict below.
                if "DESeq sample" in sample.annotation_dict:
                    row_dict['run'] = sample.name
                else:
                    row_dict['sample'] = sample.name
                """ @type row_dict: dict[str, str | unicode] """
                # Set additional columns from the Sample Annotation Sheet prefixed with 'Sample DESeq *'.
                for annotation_key in filter(lambda x: x.startswith('DESeq '), sample.annotation_dict.keys()):
                    row_dict[annotation_key[6:]] = sample.annotation_dict[annotation_key][0]

                annotation_sheet.row_dicts.append(row_dict)

            # Write the DESeq Sample Annotation Sheet to disk.
            annotation_sheet.to_file_path(adjust_field_names=True)

            # Re-write the Annotation Sheet objects for DESeq designs and contrasts to new file paths in the
            # genome directory. Convert the CSV configuration tables into TSV analysis tables.
            design_sheet.file_path = os.path.join(comparison_directory, prefix + '_designs.tsv')
            design_sheet.file_type = 'excel-tab'
            design_sheet.to_file_path()

            contrast_sheet.file_path = os.path.join(comparison_directory, prefix + '_contrasts.tsv')
            contrast_sheet.file_type = 'excel-tab'
            contrast_sheet.to_file_path()

        return

    def report(self):
        """Create a C{bsf.analyses.rnaseq.DESeq} report in HTML format.
        @return:
        @rtype:
        """

        def relative_image_source(prefix, suffix):
            """Get a relative HTML image source path.

            @param prefix: Prefix
            @type prefix: str
            @param suffix: Suffix
            @type suffix: str
            @return: Relative HTML image source path
            @rtype: str | unicode
            """
            return prefix + '/' + prefix + '_' + suffix

        def report_html():
            """Private function to create a HTML report.

            @return:
            @rtype:
            """

            # Create a symbolic link containing the project name and a UUID.
            self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str | unicode] """

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.extend(self.get_html_genome(genome_version=self.genome_version))
            str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
            str_list.append('\n')

            str_list.append('<p id="star_aligner">')
            str_list.append('<a href="https://github.com/alexdobin/STAR">STAR</a> ')
            str_list.append('aligns RNA-seq reads to a reference genome in order to identify ')
            str_list.append('exon-exon splice junctions.\n')
            # str_list.append('<br />\n')
            str_list.append('Please see the ')
            str_list.append('<a href="' + bsf.analyses.star_aligner.StarAligner.prefix + '_report.html">')
            str_list.append(self.project_name + ' ' + bsf.analyses.star_aligner.StarAligner.name)
            str_list.append('</a> report for quality plots and ')
            str_list.append('a link to alignment visualisation in the UCSC Genome Browser.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<p id="deseq2">')
            str_list.append('The <a href="https://bioconductor.org/">Bioconductor</a> ')
            str_list.append('<a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">DESeq2</a> ')
            str_list.append('package estimates variance-mean dependence in count data from ')
            str_list.append('high-throughput sequencing assays and tests for differential expression ')
            str_list.append('based on a model using the negative binomial distribution.')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Likelihood Ratio Testing (LRT) Table

            # Read the designs table as a backup, in case the design-specific LRT summary table is not available.

            design_sheet = bsf.annotation.AnnotationSheet.from_file_path(
                file_path=self.comparison_path,
                file_type='excel',
                name='DESeq Design Table')

            str_list.append('<h2 id="lrt">Likelihood Ratio Testing (LRT)</h2>\n')
            str_list.append('\n')

            str_list.append('<p>')
            str_list.append('A ')
            str_list.append('<a href="https://en.wikipedia.org/wiki/Likelihood-ratio_test">Likelihood-ratio test</a> ')
            str_list.append('compares the goodness of fit of two models, a null (i.e. full) model and an alternative ')
            str_list.append('(i.e. reduced) model. The likelihood ratio expresses how many times more likely ')
            str_list.append('the data fits one model than the other. ')
            str_list.append('Since each gene is modelled, the LRT allows identifying those genes that benefit ')
            str_list.append('from a variable or factor dropped in the reduced model.')
            str_list.append('</p>\n')
            str_list.append('<p>')
            str_list.append('Please see also the DESeq2 manual section on ')
            str_list.append('<a href="https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/')
            str_list.append('inst/doc/DESeq2.html#likelihood-ratio-test">Likelihood-ratio testing</a>.')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="lrt_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Model Formula</th>\n')
            str_list.append('<th class="left">Reduced Name</th>\n')
            str_list.append('<th class="left">Reduced Formula</th>\n')
            str_list.append('<th class="left">Differential Genes</th>\n')
            str_list.append('<th class="left">Significant Genes</th>\n')
            str_list.append('<th class="left">Significant Number</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_row_dict in design_sheet.row_dicts:
                design_prefix = '_'.join((self.prefix, design_row_dict['design']))
                lrt_summary_path = os.path.join(
                    self.genome_directory,
                    design_prefix,
                    design_prefix + '_lrt_summary.tsv')
                if os.path.exists(lrt_summary_path):
                    lrt_summary_sheet = bsf.annotation.AnnotationSheet.from_file_path(
                        file_path=lrt_summary_path,
                        file_type='excel-tab',
                        name='DESeq LRT Summary Table')
                    lrt_row_dict_list = lrt_summary_sheet.row_dicts
                else:
                    lrt_row_dict_list = list()
                    for reduced_tuple in design_row_dict['reduced_formulas'].split(';'):
                        reduced_name, reduced_formula = reduced_tuple.split(':')
                        lrt_row_dict_list.append({
                            'design': design_row_dict['design'],
                            'full_formula': design_row_dict['full_formula'],
                            'reduced_name': reduced_name,
                            'reduced_formula': reduced_formula,
                            'significant': '',
                        })

                # for reduced_tuple in design_row_dict['reduced_formulas'].split(';'):
                #     reduced_name, reduced_formula = reduced_tuple.split(':')
                for lrt_row_dict in lrt_row_dict_list:
                    str_list.append('<tr>\n')
                    str_list.append('<td>' + lrt_row_dict['design'] + '</td>\n')
                    str_list.append('<td>' + lrt_row_dict['full_formula'] + '</td>\n')
                    str_list.append('<td>' + lrt_row_dict['reduced_name'] + '</td>\n')
                    str_list.append('<td>' + lrt_row_dict['reduced_formula'] + '</td>\n')
                    # Differential genes
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('lrt', lrt_row_dict['reduced_name'] + '.tsv')),
                                        text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                                    '</td>\n')
                    # Significant genes
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('lrt', lrt_row_dict['reduced_name'], 'significant.tsv')),
                                        text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                                    '</td>\n')
                    if lrt_row_dict['significant']:
                        str_list.append('<td class="right">{:,}</td>\n'.format(int(lrt_row_dict['significant'])))
                    else:
                        str_list.append('<td class="right"></td>\n')
                    str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Contrast table.

            str_list.append('<h2 id="contrasts">Contrasts</h2>\n')
            str_list.append('\n')

            str_list.append('<table id="contrasts_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Contrast</th>\n')
            str_list.append('<th class="left">Differential Genes</th>\n')
            str_list.append('<th class="left">Significant Genes</th>\n')
            str_list.append('<th class="left">Significant Number</th>\n')
            str_list.append('<th class="left">MA Plot</th>\n')
            str_list.append('<th class="left">Volcano Plot</th>\n')
            str_list.append('<th class="left">Numerator</th>\n')
            str_list.append('<th class="left">Denominator</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            # Read the contrast file as a backup, in case a design-specific summary table was not available.

            contrast_sheet = bsf.annotation.AnnotationSheet.from_file_path(
                file_path=self.contrast_path,
                file_type='excel',
                name='DESeq Contrasts Table')

            # Re-index the contrast sheet on design names.
            contrast_dict = dict()
            """ @type contrast_dict: dict[str, list] """

            for contrast_row_dict in contrast_sheet.row_dicts:
                # Exclude contrasts if requested.
                if contrast_sheet.get_boolean(row_dict=contrast_row_dict, key='Exclude'):
                    continue

                if contrast_row_dict['Design'] not in contrast_dict:
                    contrast_dict[contrast_row_dict['Design']] = list()
                contrast_dict[contrast_row_dict['Design']].append(contrast_row_dict)

            for design_name in sorted(contrast_dict):
                design_prefix = '_'.join((self.prefix, design_name))
                contrasts_summary_path = os.path.join(
                    self.genome_directory,
                    design_prefix,
                    design_prefix + '_contrasts_summary.tsv')
                if os.path.exists(contrasts_summary_path):
                    summary_sheet = bsf.annotation.AnnotationSheet.from_file_path(
                        file_path=contrasts_summary_path,
                        file_type='excel-tab',
                        name='DESeq Contrasts Summary Table')
                    summary_row_dicts = summary_sheet.row_dicts
                else:
                    summary_row_dicts = contrast_dict[design_name]

                for row_dict in summary_row_dicts:
                    str_list.append('<tr>\n')
                    # Design
                    str_list.append('<td class="left">' + row_dict['Design'] + '</td>\n')
                    # Contrast
                    str_list.append('<td class="left">' + row_dict['Label'] + '</td>\n')
                    # TSV
                    numerator = row_dict['Numerator'].replace(',', '_')
                    denominator = row_dict['Denominator'].replace(',', '_')
                    if not denominator:
                        denominator = 'intercept'
                    # Differential Genes
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('contrast', numerator, 'against', denominator,
                                                         'genes.tsv')),
                                        text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                                    '</td>\n')
                    # Significant Genes
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('contrast', numerator, 'against', denominator,
                                                         'significant.tsv')),
                                        text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                                    '</td>\n')
                    # Significant Number
                    if 'Significant' in row_dict:
                        str_list.append('<td class="right">{:,}</td>\n'.format(int(row_dict['Significant'])))
                    else:
                        str_list.append('<td></td>\n')
                    # MA Plot
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('contrast', numerator, 'against', denominator, 'ma.pdf')),
                                        text=self.get_html_image(
                                            prefix=design_prefix,
                                            suffix='_'.join(('contrast', numerator, 'against', denominator, 'ma.png')),
                                            text='MA plot',
                                            height='80',
                                            width='80')) +
                                    '</td>\n')
                    # Volcano Plot
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join((
                                            'contrast', numerator, 'against', denominator, 'volcano.pdf')),
                                        text=self.get_html_image(
                                            prefix=design_prefix,
                                            suffix='_'.join((
                                                'contrast', numerator, 'against', denominator, 'volcano.png')),
                                            text='Volcano plot',
                                            height='80',
                                            width='80')) +
                                    '</td>\n')
                    # Numerator
                    str_list.append('<td class="left">' + row_dict['Numerator'] + '</td>\n')
                    # Denominator
                    str_list.append('<td class="left">' + row_dict['Denominator'] + '</td>\n')
                    str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Link Enrichr and Heatmap reports.

            str_list.append('<h2 id="accessory_reports">Accessory Reports</h2>\n')
            str_list.append('\n')

            str_list.append('<p>Accessory reports provide annotated expression heat map plots and\n')
            str_list.append('<a href="https://amp.pharm.mssm.edu/Enrichr/">Enrichr</a> results of\n')
            str_list.append('up- and down-regulated genes for selected libraries.</p>\n')

            str_list.append('<table id="contrasts_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Enrichr Report</th>\n')
            str_list.append('<th class="left">Heatmap Report</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_name in sorted(contrast_dict):
                str_list.append('<tr>\n')
                str_list.append('<td>' + design_name + '</td>\n')

                enrichr_prefix = '_'.join((self.prefix, design_name, 'enrichr'))
                if os.path.exists(os.path.join(
                        self.genome_directory,
                        enrichr_prefix,
                        enrichr_prefix + '_report.html')):
                    str_list.append('<td><a href="' +
                                    '/'.join((enrichr_prefix, enrichr_prefix + '_report.html')) +
                                    '">HTML</a></td>\n')
                else:
                    str_list.append('<td></td>\n')

                heatmap_prefix = '_'.join((self.prefix, design_name, 'heatmap'))
                if os.path.exists(os.path.join(
                        self.genome_directory,
                        heatmap_prefix,
                        heatmap_prefix + '_report_model.html')):
                    str_list.append('<td><a href="' +
                                    '/'.join((heatmap_prefix, heatmap_prefix + '_report_model.html')) +
                                    '">HTML</a></td>\n')
                else:
                    str_list.append('<td></td>\n')

                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Link MDS and PCA plots.

            str_list.append('<h2 id="plots">Plots</h2>\n')
            str_list.append('\n')

            str_list.append('<p>')
            str_list.append('For each design, multi-dimensional scaling (MDS), principal component analysis (PCA) ')
            str_list.append('and heatmap plots are provided for combinations of variables or factors. ')
            str_list.append('Variance by principal component plots show the distribution of the variance for ')
            str_list.append('a maximum of the first hundred components. ')
            str_list.append('Two sets of plots are available, based on data in model-aware or blind mode. ')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="plot_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">MDS Blind</th>\n')
            str_list.append('<th class="left">MDS Model</th>\n')
            str_list.append('<th class="left">PCA Blind</th>\n')
            str_list.append('<th class="left">PCA Model</th>\n')
            str_list.append('<th class="left">Heatmap Blind</th>\n')
            str_list.append('<th class="left">Heatmap Model</th>\n')
            # str_list.append('<th>Plot Aesthetics</th>\n')
            str_list.append('<th class="left">Cook\'s Distances</th>\n')
            str_list.append('<th class="left">FPKM Density</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_row_dict in design_sheet.row_dicts:
                design_prefix = '_'.join((self.prefix, design_row_dict['design']))
                for plot_instance in design_row_dict['plot_aes'].split('|'):
                    ggplot_geom, ggplot_aes_list = plot_instance.split(':')
                    plot_path = '__'.join((ggplot_geom[5:], ggplot_aes_list.replace(',', '_').replace('=', '_')))
                    str_list.append('<tr>\n')
                    # Design
                    str_list.append('<td>' + design_row_dict['design'] + '</td>\n')
                    for plot_type in ('mds', 'pca'):
                        for model_type in ('blind', 'model'):
                            plot_path_pdf = '_'.join((plot_type, plot_path, model_type + '.pdf'))
                            plot_path_png = '_'.join((plot_type, plot_path, model_type + '.png'))
                            str_list.append('<td>')
                            if os.path.exists(
                                    os.path.join(
                                        self.genome_directory,
                                        relative_image_source(prefix=design_prefix, suffix=plot_path_png))):
                                str_list.append(
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix=plot_path_pdf,
                                        text=self.get_html_image(
                                            prefix=design_prefix,
                                            suffix=plot_path_png,
                                            text=plot_type + ' plot',
                                            height='80',
                                            width='80')))
                            str_list.append('</td>\n')
                    # Heatmap plots
                    plot_path = '_'.join(map(lambda x: x.split('=')[1], ggplot_aes_list.split(',')))
                    plot_type = 'heatmap'
                    for model_type in ('blind', 'model'):
                        plot_path_pdf = '_'.join((plot_type, plot_path, model_type + '.pdf'))
                        plot_path_png = '_'.join((plot_type, plot_path, model_type + '.png'))
                        str_list.append('<td>')
                        if os.path.exists(
                                os.path.join(
                                    self.genome_directory,
                                    relative_image_source(prefix=design_prefix, suffix=plot_path_png))):
                            str_list.append(
                                self.get_html_anchor(
                                    prefix=design_prefix,
                                    suffix=plot_path_pdf,
                                    text=self.get_html_image(
                                        prefix=design_prefix,
                                        suffix=plot_path_png,
                                        text=plot_type + ' plot',
                                        height='80',
                                        width='80')))
                        str_list.append('</td>\n')

                    # Plot Aesthetics
                    # str_list.append('<td>' + plot_instance + '</td>\n')
                    # Cook's Distances
                    str_list.append('<td></td>\n')
                    # FPKM Density
                    str_list.append('<td></td>\n')
                    str_list.append('</tr>\n')

                # Add a line with the variance per principal component plots.
                str_list.append('<tr>\n')
                str_list.append('<td>' + design_row_dict['design'] + '</td>\n')
                str_list.append('<td></td>\n')  # MDS Blind
                str_list.append('<td></td>\n')  # MDS Model

                # PCA Variance per component for Blind and Model.
                plot_type = 'pca'
                plot_path = 'variance'
                for model_type in ('blind', 'model'):
                    plot_path_pdf = '_'.join((plot_type, plot_path, model_type + '.pdf'))
                    plot_path_png = '_'.join((plot_type, plot_path, model_type + '.png'))
                    str_list.append('<td>')
                    if os.path.exists(
                            os.path.join(
                                self.genome_directory,
                                relative_image_source(prefix=design_prefix, suffix=plot_path_png))):
                        str_list.append(
                            self.get_html_anchor(
                                prefix=design_prefix,
                                suffix=plot_path_pdf,
                                text=self.get_html_image(
                                    prefix=design_prefix,
                                    suffix=plot_path_png,
                                    text=plot_type + ' plot',
                                    height='80',
                                    width='80')))
                    str_list.append('</td>\n')

                str_list.append('<td></td>\n')  # Heatmap Blind
                str_list.append('<td></td>\n')  # Heatmap Model

                # Cook's Distance Box Plot
                str_list.append('<td>\n')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='cooks_distances.pdf',
                    text=self.get_html_image(
                        prefix=design_prefix,
                        suffix='cooks_distances.png',
                        text='Cook\'s distance box plot',
                        height='80',
                        width='80')))
                str_list.append('</td>\n')

                # FPKM Density
                str_list.append('<td>\n')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='fpkm_density.pdf',
                    text=self.get_html_image(
                        prefix=design_prefix,
                        suffix='fpkm_density.png',
                        text='FPKM density plot',
                        height='80',
                        width='80')))
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Additional tables

            str_list.append('<h2 id="tables">Additional Tables</h2>\n')
            str_list.append('\n')

            str_list.append('<table id="table_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Raw Counts</th>\n')
            str_list.append('<th class="left">VST Blind Counts</th>\n')
            str_list.append('<th class="left">VST Model Counts</th>\n')
            str_list.append('<th class="left">FPKMs</th>\n')
            str_list.append('<th class="left">Samples</th>\n')
            str_list.append('<th class="left">Annotation</th>\n')
            str_list.append('<th class="left">Contrasts</th>\n')
            str_list.append('<th class="left">LRT</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_row_dict in design_sheet.row_dicts:
                design_prefix = '_'.join((self.prefix, design_row_dict['design']))

                str_list.append('<tr>\n')

                str_list.append('<td>' + design_row_dict['design'] + '</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='counts_raw.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='counts_vst_blind.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='counts_vst_model.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='fpkms.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='samples.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='annotation.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='contrasts_summary.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('<td>')
                str_list.append(self.get_html_anchor(
                    prefix=design_prefix,
                    suffix='lrt_summary.tsv',
                    text='<abbr title="Tab-Separated Value">TSV</abbr>'))
                str_list.append('</td>\n')

                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

            if self.debug > 0:
                print('Report list:')
                print(repr(str_list))

            return

        report_html()

        return
