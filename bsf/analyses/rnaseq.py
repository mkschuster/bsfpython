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
"""The :py:mod:`bsf.analyses.rnaseq` module provides two main classes supporting RNA-seq analyses.

    - The :py:class:`bsf.analyses.rnaseq.DESeq` class models an RNA-seq analysis based on the
        `Bioconductor <https://bioconductor.org/>`_
        `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ package.

    - The :py:class:`bsf.analyses.rnaseq.Tuxedo` class models an RNA-seq analysis based on the
        `Cufflinks <http://cole-trapnell-lab.github.io/cufflinks/>`_ package.
"""
import errno
import logging
import os
import pickle
import re
from typing import Optional

from bsf.analyses.hisat import Hisat2
from bsf.analyses.kallisto import Kallisto
from bsf.analyses.star import Star
from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.connector import ConnectorFile
from bsf.ngs import Collection, Sample, SampleGroup
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, Executable, \
    RunnableStep, RunnableStepCopy, RunnableStepLink, RunnableStepMakeDirectory, RunnableStepMove, \
    RunnableStepSetEnvironment
from bsf.standards import Configuration, StandardFilePath, Index, Transcriptome

module_logger = logging.getLogger(name=__name__)


class FilePathTophat(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathTophat` class models files in a sample-specific TopHat directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar accepted_hits_bam: TopHat accepted hits BAM file
    :type accepted_hits_bam: str
    :ivar accepted_hits_bam_link_source: TopHat accepted hits BAM file symbolic link source
    :type accepted_hits_bam_link_source: str
    :ivar accepted_hits_bam_link_target: TopHat accepted hits BAM file symbolic link target
    :type accepted_hits_bam_link_target: str
    :ivar accepted_hits_bai: TopHat accepted hits BAI file
    :type accepted_hits_bai: str
    :ivar accepted_hits_bai_link_source: TopHat accepted hits BAI file symbolic link source
    :type accepted_hits_bai_link_source: str
    :ivar accepted_hits_bai_link_target: TopHat accepted hits BAI file symbolic link target
    :type accepted_hits_bai_link_target: str
    :ivar accepted_hits_bw: TopHat accepted hits bigWig file
    :type accepted_hits_bw: str
    :ivar align_summary: TopHat align summary file
    :type align_summary: str
    :ivar deletions_bb: TopHat deletions bigBed file
    :type deletions_bb: str
    :ivar deletions_bed: TopHat deletions BED file
    :type deletions_bed: str
    :ivar insertions_bb: TopHat insertions bigBed file
    :type insertions_bb: str
    :ivar insertions_bed: TopHat insertions BED file
    :type insertions_bed: str
    :ivar junctions_bb: TopHat junctions bigBed file
    :type junctions_bb: str
    :ivar junctions_bed: TopHat junctions BED file
    :type junctions_bed: str
    :ivar prep_reads_info: TopHat prepare reads information file
    :type prep_reads_info: str
    :ivar unmapped_bam: TopHat unmapped BAM file
    :type unmapped_bam: str
    :ivar unmapped_bam_link_source: TopHat unmapped BAM file symbolic link source
    :type unmapped_bam_link_source: str
    :ivar unmapped_bam_link_target: TopHat unmapped BAM file symbolic link target
    :type unmapped_bam_link_target: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathTophat` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
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
    """The :py:class:`bsf.analyses.rnaseq.FilePathCufflinks` class models files in a
    sample-specific Cufflinks directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar fpkm_tracking_genes_tsv: Cufflinks FPKM tracking genes :emphasis:`Tab-Separated Value` (TSV) file
    :type fpkm_tracking_genes_tsv: str
    :ivar fpkm_tracking_isoforms_tsv: Cufflinks FPKM tracking isoforms :emphasis:`Tab-Separated Value` (TSV) file
    :type fpkm_tracking_isoforms_tsv: str
    :ivar skipped_gtf: Cufflinks skipped regions GTF file
    :type skipped_gtf: str
    :ivar skipped_gtf_link_source: Cufflinks skipped regions GTF symbolic link source
    :type skipped_gtf_link_source: str
    :ivar skipped_gtf_link_target: Cufflinks skipped regions GTF symbolic link target
    :type skipped_gtf_link_target: str
    :ivar temporary_big_gene_prediction: Temporary UCSC big gene prediction (bigGenePred) file
    :type temporary_big_gene_prediction: str
    :ivar temporary_gene_prediction: Temporary UCSC gene prediction (genePred) file
    :type temporary_gene_prediction: str
    :ivar temporary_sorted_tsv: Temporary sorted :emphasis:`Tab-Separated Value` (TSV) file
    :type temporary_sorted_tsv: str
    :ivar temporary_slopped_tsv: Temporary slopped (bedtools slop) :emphasis:`Tab-Separated Value` (TSV) file
    :type temporary_slopped_tsv: str
    :ivar temporary_fixed_tsv: Temporary slopped and fixed (tab at end) :emphasis:`Tab-Separated Value` (TSV) file
    :type temporary_fixed_tsv: str
    :ivar transcripts_bb: Cufflinks transcript assembly bigBed file
    :type transcripts_bb: str
    :ivar transcripts_bb_link_source: Cufflinks transcript assembly bigBed symbolic link source
    :type transcripts_bb_link_source: str
    :ivar transcripts_bb_link_target: Cufflinks transcript assembly bigBed symbolic link target
    :type transcripts_bb_link_target: str
    :ivar transcripts_gtf: Cufflinks transcript assembly GTF file
    :type transcripts_gtf: str
    :ivar transcripts_gtf_link_source: Cufflinks transcript assembly GTF symbolic link source
    :type transcripts_gtf_link_source: str
    :ivar transcripts_gtf_link_target: Cufflinks transcript assembly GTF symbolic link target
    :type transcripts_gtf_link_target: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathCufflinks` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
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


class FilePathCuffmerge(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathCuffmerge` class models files in a
    comparison-specific Cuffmerge directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar assembly_txt: Assembly text file
    :type assembly_txt: str
    :ivar merged_bb: Cuffmerge transcript assembly bigBed file
    :type merged_bb: str
    :ivar merged_bb_link_source: Cuffmerge transcript assembly bigBed symbolic link source
    :type merged_bb_link_source: str
    :ivar merged_bb_link_target: Cuffmerge transcript assembly bigBed symbolic link target
    :type merged_bb_link_target: str
    :ivar merged_gtf: Cuffmerge merged GTF file
    :type merged_gtf: str
    :ivar merged_gtf_link_source: Cuffmerge transcript assembly GTF symbolic link source
    :type merged_gtf_link_source: str
    :ivar merged_gtf_link_target: Cuffmerge transcript assembly GTF symbolic link target
    :type merged_gtf_link_target: str
    :ivar temporary_gene_prediction: Temporary UCSC gene prediction (genePred) file
    :type temporary_gene_prediction: str
    :ivar temporary_big_gene_prediction: Temporary UCSC big gene prediction (bigGenePred) file
    :type temporary_big_gene_prediction: str
    :ivar temporary_sorted_tsv: Temporary sorted :emphasis:`Tab-Separated Value` (TSV) file
    :type temporary_sorted_tsv: str
    :ivar cuffcompare_prefix: Cuffcompare output prefix, including the cuffmerge directory path
    :type cuffcompare_prefix: str
    :ivar cuffcompare_combined_gtf: Cuffcompare merged GTF file
    :type cuffcompare_combined_gtf: str
    :ivar cuffcompare_loci: Cuffcompare loci file
    :type cuffcompare_loci: str
    :ivar cuffcompare_stats: Cuffcompare stats file
    :type cuffcompare_stats: str
    :ivar cuffcompare_tracking: Cuffcompare tracking file
    :type cuffcompare_tracking: str
    :ivar cuffcompare_refmap: Cuffcompare merged.gtf.refmap
    :type cuffcompare_refmap: str
    :ivar cuffcompare_tmap: Cuffcompare merged.gtf.tmap
    :type cuffcompare_tmap: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathCuffmerge` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
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


class FilePathCuffquant(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathCuffquant` class models files in a
    sample-specific Cuffquant directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar abundances: Cuffquant abundances file
    :type abundances: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathCuffquant` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathCuffquant, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances = os.path.join(prefix, 'abundances.cxb')

        return


class FilePathCuffnorm(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathCuffnorm` class models files in a
    comparison-specific Cuffnorm directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar abundances_tsv: Abundances TSV file
    :type abundances_tsv: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathCuffnorm` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathCuffnorm, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances_tsv = prefix + '_abundances.tsv'

        return


class FilePathCuffdiff(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathCuffdiff` class models files in a
    comparison-specific Cuffdiff directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathCuffdiff` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathCuffdiff, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances_tsv = prefix + '_abundances.tsv'
        self.alignments_tsv = prefix + '_alignments.tsv'

        return


class FilePathProcessCuffdiff(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathProcessCuffdiff` class models files in a
    comparison-specific Cuffdiff directory.
    """

    pass


class FilePathMonocle(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathMonocle` class models files in a
    comparison-specific Monocle directory.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar annotation_tsv: Monocle annotation TSV
    :type annotation_tsv: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathMonocle` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathMonocle, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.annotation_tsv = prefix + '_annotation.tsv'

        return


class TuxedoSamplePairSheet(AnnotationSheet):
    """The :py:class:`bsf.analyses.rnaseq.TuxedoSamplePairSheet` class represents
    :py:class:`bsf.ngs.Sample` object pairs.

    The :py:class:`bsf.ngs.Sample` object pairs are defined by the :literal:`bsf_rnaseq_process_cuffdiff.R` script.
    """

    _file_type = 'excel-tab'

    _field_name_list = [
        'V1',
        'V2',
    ]


class Tuxedo(Analysis):
    """The :py:class:`bsf.analyses.rnaseq.Tuxedo` class models an RNA-seq analysis based on the
    `Cufflinks <http://cole-trapnell-lab.github.io/cufflinks/>`_ package as part of the Tuxedo suite.

    :ivar replicate_grouping: Request grouping all replicates into a single Tophat and Cufflinks process
    :type replicate_grouping: bool | None
    :ivar comparison_path: Comparison file path
    :type comparison_path: str | None
    :ivar genome_fasta_path: Reference genome sequence :emphasis:`FASTA` file path
    :type genome_fasta_path: str | None
    :ivar genome_index_path: Bowtie genome index path
    :type genome_index_path: str | None
    :ivar genome_sizes_path: Reference genome sizes file path
    :type genome_sizes_path: str | None
    :ivar transcriptome_version: A transcriptome version.
    :type transcriptome_version: str | None
    :ivar transcriptome_gtf: A transcriptome annotation GTF file path.
    :type transcriptome_gtf: str | None
    :ivar transcriptome_index: A transcriptome index directory path.
    :type transcriptome_index: str | None
    :ivar insert_size: An insert size.
    :type insert_size: int | None
    :ivar insert_size_sd: An insert size standard deviation.
    :type insert_size_sd: int | None
    :ivar read_length: A read length.
    :type read_length: int | None
    :ivar mask_gtf_path: A GTF file path to mask transcripts.
    :type mask_gtf_path: str | None
    :ivar multi_read_correction: Apply multi-read correction
    :type multi_read_correction: bool | None
    :ivar library_type: Library type Cuffquant and Cuffdiff
        :literal:`fr-unstranded` (default), :literal:`fr-firststrand` or :literal:`fr-secondstrand`
    :type library_type: str | None
    :ivar novel_transcripts: Assemble novel transcripts
    :type novel_transcripts: bool | None
    :ivar false_discovery_rate: False discovery rate (FDR) threshold
    :type false_discovery_rate: float | None
    :ivar no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
    :type no_length_correction: bool | None
    :ivar aligner: Alignment program
    :type aligner: str | None
    :ivar ucsc_autosql_big_gene_prediction: A UCSC bigGenePred AutoSQL file path.
    :type ucsc_autosql_big_gene_prediction: str | None
    """

    name = 'RNA-seq Analysis'
    prefix = 'rnaseq'

    @classmethod
    def get_stage_name_run_tophat(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'run_tophat'))

    @classmethod
    def get_stage_name_process_tophat(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'process_tophat'))

    @classmethod
    def get_stage_name_run_cufflinks(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'run_cufflinks'))

    @classmethod
    def get_stage_name_process_cufflinks(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'process_cufflinks'))

    # Comparison stage

    @classmethod
    def get_stage_name_run_cuffmerge(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cuffmerge'))

    @classmethod
    def get_stage_name_run_cuffquant(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cuffquant'))

    @classmethod
    def get_stage_name_run_cuffnorm(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cuffnorm'))

    @classmethod
    def get_stage_name_run_cuffdiff(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'cuffdiff'))

    @classmethod
    def get_stage_name_process_cuffdiff(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'process_cuffdiff'))

    @classmethod
    def get_stage_name_monocle(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'monocle'))

    @classmethod
    def get_prefix_run_tophat(cls, sample_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_run_tophat(), sample_name))

    @classmethod
    def get_prefix_process_tophat(cls, sample_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_process_tophat(), sample_name))

    @classmethod
    def get_prefix_run_cufflinks(cls, sample_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cufflinks(), sample_name))

    @classmethod
    def get_prefix_process_cufflinks(cls) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return cls.get_stage_name_process_cufflinks()

    @classmethod
    def get_prefix_run_cuffmerge(cls, comparison_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffmerge(), comparison_name))

    @classmethod
    def get_prefix_run_cuffquant(cls, comparison_name: str, sample_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffquant(), comparison_name, sample_name))

    @classmethod
    def get_prefix_run_cuffnorm(cls, comparison_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffnorm(), comparison_name))

    @classmethod
    def get_prefix_run_cuffdiff(cls, comparison_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_run_cuffdiff(), comparison_name))

    @classmethod
    def get_prefix_process_cuffdiff(cls, comparison_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_process_cuffdiff(), comparison_name))

    @classmethod
    def get_prefix_monocle(cls, comparison_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_monocle(), comparison_name))

    @classmethod
    def get_file_path_run_tophat(cls, sample_name: str) -> FilePathTophat:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathTophat` object.

        The prefix is non-standard, as :literal:`rnaseq_run_tophat` and :literal:`rnaseq_process_tophat`
        use the same :literal:`rnaseq_tophat` prefix.

        :param sample_name: A sample name
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathTophat` object.
        :rtype: FilePathTophat
        """
        return FilePathTophat(
            prefix='_'.join(('rnaseq_tophat', sample_name)))

    @classmethod
    def get_file_path_process_tophat(cls, sample_name: str) -> FilePathTophat:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathTophat` object.

        The prefix is non-standard, as :literal:`rnaseq_run_tophat` and :literal:`rnaseq_process_tophat`
        use the same :literal:`rnaseq_tophat` prefix.

        :param sample_name: A sample name.
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathTophat` object.
        :rtype: FilePathTophat
        """
        return FilePathTophat(
            prefix='_'.join(('rnaseq_tophat', sample_name)))

    @classmethod
    def get_file_path_run_cufflinks(cls, sample_name: str) -> FilePathCufflinks:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathCufflinks` object.

        The prefix is non-standard, as :literal:`rnaseq_run_cufflinks` and :literal:`rnaseq_process_cufflinks`
        use the same :literal:`rnaseq_cufflinks` prefix.

        :param sample_name: A sample name.
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathCufflinks` object.
        :rtype: FilePathCufflinks
        """
        return FilePathCufflinks(
            prefix='_'.join(('rnaseq_cufflinks', sample_name)))

    @classmethod
    def get_file_path_process_cufflinks(cls, sample_name: str) -> FilePathCufflinks:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathCufflinks` object.

        The prefix is non-standard, as :literal:`rnaseq_run_cufflinks` and :literal:`rnaseq_process_cufflinks`
        use the same :literal:`rnaseq_cufflinks` prefix.

        :param sample_name: A sample name.
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathCufflinks` object.
        :rtype: FilePathCufflinks
        """
        return FilePathCufflinks(
            prefix='_'.join(('rnaseq_cufflinks', sample_name)))

    @classmethod
    def get_file_path_cuffmerge(cls, comparison_name: str) -> FilePathCuffmerge:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathCuffmerge` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathCuffmerge` object.
        :rtype: FilePathCuffmerge
        """
        return FilePathCuffmerge(
            prefix=cls.get_prefix_run_cuffmerge(comparison_name=comparison_name))

    @classmethod
    def get_file_path_cuffquant(cls, comparison_name: str, sample_name: str) -> FilePathCuffquant:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathCuffquant` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :param sample_name: A sample name
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathCuffquant` object.
        :rtype: FilePathCuffquant
        """
        return FilePathCuffquant(
            prefix=cls.get_prefix_run_cuffquant(comparison_name=comparison_name, sample_name=sample_name))

    @classmethod
    def get_file_path_cuffnorm(cls, comparison_name: str) -> FilePathCuffnorm:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathCuffnorm` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathCuffnorm` object.
        :rtype: FilePathCuffnorm
        """
        return FilePathCuffnorm(
            prefix=cls.get_prefix_run_cuffnorm(comparison_name=comparison_name))

    @classmethod
    def get_file_path_run_cuffdiff(cls, comparison_name: str) -> FilePathCuffdiff:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathCuffdiff` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathCuffdiff` object.
        :rtype: FilePathCuffdiff
        """
        return FilePathCuffdiff(
            prefix=cls.get_prefix_run_cuffdiff(comparison_name=comparison_name))

    @classmethod
    def get_file_path_process_cuffdiff(cls, comparison_name: str) -> FilePathProcessCuffdiff:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathProcessCuffdiff` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathProcessCuffdiff` object.
        :rtype: FilePathProcessCuffdiff
        """
        return FilePathProcessCuffdiff(
            prefix=cls.get_prefix_process_cuffdiff(comparison_name=comparison_name))

    @classmethod
    def get_file_path_monocle(cls, comparison_name: str) -> FilePathMonocle:
        """Get a :py:class:`bsf.analyses.rnaseq.FilePathMonocle` object.

        :param comparison_name: A comparison name.
        :type comparison_name: str
        :return: A :py:class:`bsf.analyses.rnaseq.FilePathMonocle` object.
        :rtype: FilePathMonocle
        """
        return FilePathMonocle(
            prefix=cls.get_prefix_monocle(comparison_name=comparison_name))

    def __init__(
            self,
            configuration: Optional[Configuration] = None,
            project_name: Optional[str] = None,
            genome_version: Optional[str] = None,
            input_directory: Optional[str] = None,
            output_directory: Optional[str] = None,
            project_directory: Optional[str] = None,
            genome_directory: Optional[str] = None,
            report_style_path: Optional[str] = None,
            report_header_path: Optional[str] = None,
            report_footer_path: Optional[str] = None,
            e_mail: Optional[str] = None,
            stage_list: Optional[list[Stage]] = None,
            collection: Optional[Collection] = None,
            sample_list: Optional[list[Sample]] = None,
            replicate_grouping: Optional[bool] = None,
            comparison_path: Optional[str] = None,
            genome_fasta_path: Optional[str] = None,
            genome_index_path: Optional[str] = None,
            genome_sizes_path: Optional[str] = None,
            transcriptome_version: Optional[str] = None,
            transcriptome_gtf: Optional[str] = None,
            transcriptome_index: Optional[str] = None,
            insert_size: Optional[int] = None,
            insert_size_sd: Optional[int] = None,
            read_length: Optional[int] = None,
            mask_gtf_path: Optional[str] = None,
            multi_read_correction: Optional[bool] = None,
            library_type: Optional[str] = None,
            novel_transcripts: Optional[bool] = None,
            false_discovery_rate: Optional[float] = None,
            no_length_correction: Optional[bool] = None,
            aligner: Optional[str] = None,
            ucsc_autosql_big_gene_prediction: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.Tuxedo` object.

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
        :param report_style_path: A report style :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: A report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: A report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a :emphasis:`UCSC Genome Browser Track Hub`.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param replicate_grouping: Request grouping all replicates into a single Tophat and Cufflinks process.
        :type replicate_grouping: bool | None
        :param comparison_path: Comparison file path
        :type comparison_path: str | None
        :param genome_fasta_path: Reference genome sequence :emphasis:`FASTA` file path
        :type genome_fasta_path: str | None
        :param genome_index_path: Bowtie genome index path
        :type genome_index_path: str | None
        :param genome_sizes_path: Reference genome sizes file path
        :type genome_sizes_path: str | None
        :param transcriptome_version: Transcriptome version
        :type transcriptome_version: str | None
        :param transcriptome_gtf: Transcriptome annotation GTF file path
        :type transcriptome_gtf: str | None
        :param transcriptome_index: Transcriptome index directory path
        :type transcriptome_index: str | None
        :param insert_size: Insert size
        :type insert_size: int | None
        :param insert_size_sd: Insert size standard deviation
        :type insert_size_sd: int | None
        :param read_length: Read length
        :type read_length: int | None
        :param mask_gtf_path: GTF file path to mask transcripts
        :type mask_gtf_path: str | None
        :param multi_read_correction: Apply multi-read correction
        :type multi_read_correction: bool | None
        :param library_type: Library type Cuffquant and Cuffdiff
            :literal:`fr-unstranded` (default), :literal:`fr-firststrand` or :literal:`fr-secondstrand`.
        :type library_type: str | None
        :param novel_transcripts: Assemble novel transcripts
        :type novel_transcripts: bool | None
        :param false_discovery_rate: False discovery rate (FDR) threshold
        :type false_discovery_rate: float | None
        :param no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
        :type no_length_correction: bool | None
        :param aligner: Alignment program
        :type aligner: str | None
        :param ucsc_autosql_big_gene_prediction: A UCSC bigGenePred AutoSQL file path.
        :type ucsc_autosql_big_gene_prediction: str | None
        """
        super(Tuxedo, self).__init__(
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

        self.replicate_grouping = replicate_grouping
        self.comparison_path = comparison_path
        self.genome_fasta_path = genome_fasta_path
        self.genome_index_path = genome_index_path
        self.genome_sizes_path = genome_sizes_path
        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf = transcriptome_gtf
        self.transcriptome_index = transcriptome_index
        self.insert_size = insert_size
        self.insert_size_sd = insert_size_sd
        self.read_length = read_length
        self.mask_gtf_path = mask_gtf_path
        self.multi_read_correction = multi_read_correction
        self.library_type = library_type
        self.novel_transcripts = novel_transcripts
        self.false_discovery_rate = false_discovery_rate
        self.no_length_correction = no_length_correction
        self.aligner = aligner
        self.ucsc_autosql_big_gene_prediction = ucsc_autosql_big_gene_prediction

        self._comparison_dict: dict[str, list[SampleGroup]] = dict()

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.rnaseq.Tuxedo` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
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

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index = configuration.config_parser.get(section=section, option=option)

        option = 'insert_size'
        if configuration.config_parser.has_option(section=section, option=option):
            self.insert_size = configuration.config_parser.getint(section=section, option=option)

        option = 'insert_std_dev'
        if configuration.config_parser.has_option(section=section, option=option):
            self.insert_size_sd = configuration.config_parser.getint(section=section, option=option)

        option = 'read_length'
        if configuration.config_parser.has_option(section=section, option=option):
            self.read_length = configuration.config_parser.getint(section=section, option=option)

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

        option = 'ucsc_autosql_big_gene_prediction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.ucsc_autosql_big_gene_prediction = configuration.get_list_from_csv(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.rnaseq.Tuxedo` object.
        """

        def run_read_comparisons() -> None:
            """Private function to read a :py:class:`bsf.annotation.AnnotationSheet` specifying comparisons
            from a CSV file path.

            All :py:class:`bsf.ngs.Sample` objects referenced in a comparison are added from the
            :py:attr:`bsf.analysis.Analysis.collection` instance variable
            (i.e., :py:class:`bsf.ngs.Collection` object) to the
            :py:attr:`bsf.analysis.Analysis.sample_list` instance variable.

            - Column headers for CASAVA folders

              - Treatment/Control/Point N ProcessedRunFolder

                - CASAVA processed run folder name or
                - :py:attr:`bsf.analysis.Analysis.input_directory` attribute by default

              - Treatment/Control/Point N Project

                - CASAVA Project name or
                - :py:attr:`bsf.analysis.Analysis.project_name` attribute by default

              - Treatment/Control/Point N Sample

                - CASAVA Sample name, no default

            - Column headers for independent samples

              - Treatment/Control/Point N Sample
              - Treatment/Control/Point N Reads
              - Treatment/Control/Point N File
            """
            if self.comparison_path:
                # A comparison file path was provided.
                if self.comparison_path == '*groups*':
                    # The special file name *groups* creates TuxedoComparison objects on the basis of an
                    # all-against-all group comparison.
                    # Without a comparison file path, simply add all Sample objects from the Collection.
                    self.sample_list.extend(self.collection.get_all_samples(exclude=True))

                    # Create a global comparison by adding all sample groups.
                    _sample_group_list: list[SampleGroup] = list()

                    for _group_name, _sample_list in self.collection.sample_group_dict.items():
                        _sample_group = SampleGroup(name=_group_name, sample_list=_sample_list)
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
                    _sample_group_list: list[SampleGroup] = list()

                    for _sample in self.sample_list:
                        # Sample objects are only useful, if at least one PairedReads object is not excluded.
                        if not _sample.is_excluded():
                            _sample_group_list.append(SampleGroup(name=_sample.name, sample_list=[_sample]))

                    # Sort the list of comparison groups by SampleGroup.name.
                    _sample_group_list.sort(key=lambda item: item.name)
                    # Set the comparison name to 'global'.
                    self._comparison_dict['global'] = _sample_group_list
                else:
                    # A comparison file path was provided.
                    self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)
                    # Read and process the comparison file, which includes adding only those Sample objects,
                    # which are referenced in a comparison.
                    annotation_sheet: AnnotationSheet = AnnotationSheet.from_file_path(file_path=self.comparison_path)
                    re_pattern = re.compile(pattern=r'\W')

                    for row_dict in annotation_sheet.row_dict_list:
                        _sample_group_list: list[SampleGroup] = list()
                        _comparison_name_list: list[str] = list()
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
                            _sample_list_new: list[Sample] = list()
                            if _group_name and len(_sample_list_old):
                                # Sample objects are only useful, if at least one PairedReads object is not excluded.
                                for _sample in _sample_list_old:
                                    if not _sample.is_excluded():
                                        _sample_list_new.append(_sample)

                                if len(_sample_list_new):
                                    _comparison_name_list.append(_group_name)
                                    _sample_group_list.append(SampleGroup(
                                        name=_group_name,
                                        sample_list=_sample_list_new))
                                    # Also expand each Python list of bsf.ngs.Sample objects to get all those
                                    # bsf.ngs.Sample objects that this bsf.analysis.Analysis needs considering.
                                    for _sample in _sample_list_new:
                                        self.add_sample(sample=_sample)
                                        module_logger.log(
                                            logging.DEBUG - 1,
                                            'prefix: %r Sample.name: %r Sample.file_path: %r',
                                            prefix, _sample.name, _sample.file_path)
                                        module_logger.log(logging.DEBUG - 2, 'Sample: %r', _sample)
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
                            _comparison_name = re.sub(pattern=re_pattern, repl='_', string=row_dict['Comparison Name'])
                        else:
                            _comparison_name = '__'.join(_comparison_name_list)

                        # Sort the list of comparison groups by SampleGroup.name.
                        _sample_group_list.sort(key=lambda item: item.name)
                        # Set the comparison name.
                        self._comparison_dict[_comparison_name] = _sample_group_list
            else:
                # Without a comparison file path, simply add a comparison 'global' with a SampleGroup 'global'
                # with all Sample objects from the Collection.
                # This means that most pipeline stages except Cuffdiff can run.
                self.sample_list.extend(self.collection.get_all_samples(exclude=True))
                self._comparison_dict['global'] = [SampleGroup(
                    name='global',
                    sample_list=self.collection.get_all_samples(exclude=True))]

            # Check for comparisons without SampleGroup objects or SampleGroup objects without Sample.
            # FIXME: The complication is that Sample or ReadGroup objects could be excluded from the Analysis.
            for _comparison_name, _sample_group_list in self._comparison_dict.items():
                module_logger.debug('Comparison name: %r', _comparison_name)
                module_logger.debug('SampleGroup list:')

                if len(_sample_group_list) < 1:
                    module_logger.warning('Comparison name %r without SampleGroup objects.', _comparison_name)

                for _sample_group in _sample_group_list:
                    module_logger.debug('SampleGroup.name: %r', _sample_group.name)
                    module_logger.debug('SampleGroup Sample list:')

                    if len(_sample_group.sample_list) < 1:
                        module_logger.warning('SampleGroup.name %r without Sample objects.', _sample_group.name)

                    for _sample in _sample_group.sample_list:
                        module_logger.debug('Sample.name: %r', _sample.name)

            return

        def run_write_annotation(annotation_path: str, annotation_dict: dict[str, list[str]]) -> None:
            """Private function to write a sample annotation file for Cuffdiff or Cuffnorm to disk.

            :param annotation_path: Annotation file path.
            :type annotation_path: str
            :param annotation_dict: Annotation :py:class:`dict` object.
            :type annotation_dict: dict[str, list[str]]
            """
            with open(file=annotation_path, mode='wt') as _output_text_io:
                _output_text_io.write('sample_id\tgroup_label\n')
                for _group_name in sorted(annotation_dict):
                    for _file_path in annotation_dict[_group_name]:
                        _output_text_io.write(_file_path + '\t' + _group_name + '\n')

            return

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Tuxedo requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception(f"A {self.name!s} requires a valid 'transcriptome_version' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception(f'The sample annotation sheet {self.sas_file!r} does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[self.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation sheet in the current working directory.')

        # Get the comparison annotation sheet before calling the run() method of the Analysis super-class.

        if not self.comparison_path:
            # A comparison path was not provided, check if a standard file exists in this directory.
            self.comparison_path = self.get_annotation_file(prefix_list=[self.prefix], suffix='comparisons.csv')
            if self.comparison_path:
                if not os.path.exists(self.comparison_path):
                    self.comparison_path = None
                    module_logger.debug(
                        'Standard comparison file not in current working directory: %r',
                        self.comparison_path)
                else:
                    module_logger.debug(
                        'Standard comparison file in current working directory: %r',
                        self.comparison_path)

        super(Tuxedo, self).run()

        # Method configuration regarding Cuffquant and Cuffdiff.
        run_cuffquant_before_cuffdiff = False

        run_read_comparisons()

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        # Get the genome, FASTA, index and sizes.

        if not self.genome_index_path:
            self.genome_index_path = os.path.join(
                StandardFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        if not self.genome_fasta_path:
            self.genome_fasta_path = StandardFilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if not os.path.exists(self.genome_fasta_path):
            raise Exception(f'The genome FASTA file path {self.genome_fasta_path!r} does not exist.')

        if not self.genome_sizes_path:
            self.genome_sizes_path = StandardFilePath.get_resource_genome_fasta_index(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if not os.path.exists(self.genome_sizes_path):
            raise Exception(f'The genome sizes file path {self.genome_sizes_path!r} does not exist.')

        # Define a reference transcriptome index directory or a GTF file path.

        if self.transcriptome_index:
            # Check if the transcriptome_index is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_index = self.configuration.get_absolute_path(
                file_path=self.transcriptome_index,
                default_path=StandardFilePath.get_resource_transcriptome(
                    transcriptome_version=None,
                    absolute=True))

            if not os.path.isdir(self.transcriptome_index):
                raise Exception(f'The reference transcriptome index directory {self.transcriptome_index!r} '
                                f'does not exist.')

            transcriptome_prefix = os.path.basename(self.transcriptome_index)

            # Does an indices_for_TopHat directory exist?
            transcriptome_index = os.path.join(
                self.transcriptome_index,
                Index.get(option='tophat2'))
            if os.path.isdir(transcriptome_index):
                self.transcriptome_index = transcriptome_index

            # Finally, set the transcriptome GTF file path.
            # The tophat --transcript-index process puts a GFF file into the index directory
            # that really is a GTF file. A symbolic link to a GTF file is needed to make the
            # process cuffdiff script work.
            # For the moment, use the symbolic link in the indices_for_TopHat directory.

            self.transcriptome_gtf = os.path.join(
                self.transcriptome_index,
                '.'.join((transcriptome_prefix, 'gtf')))

            if not os.path.exists(self.transcriptome_gtf):
                raise Exception(f'The reference transcriptome GTF file {self.transcriptome_gtf!r} does not exist.')
        elif self.transcriptome_gtf:
            # Check if transcriptome_gtf is absolute and if not,
            # prepend the default transcriptome directory.
            self.transcriptome_gtf = self.configuration.get_absolute_path(
                file_path=self.transcriptome_gtf,
                default_path=StandardFilePath.get_resource_transcriptome(
                    transcriptome_version=self.transcriptome_version,
                    absolute=True))

            if not os.path.exists(self.transcriptome_gtf):
                raise Exception(f'The reference transcriptome GTF file {self.transcriptome_gtf!r} does not exist.')
        else:
            # Neither was provided, automatically discover on the basis of the transcriptome version.
            self.transcriptome_index = os.path.join(
                StandardFilePath.get_resource_transcriptome_index(
                    transcriptome_version=self.transcriptome_version,
                    transcriptome_index='tophat2'),
                self.transcriptome_version,  # Tophat puts the transcriptome index into a subdirectory.
                self.transcriptome_version)  # Tophat uses a transcriptome prefix to lookup index files.

            self.transcriptome_gtf = StandardFilePath.get_resource_transcriptome_gtf(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='tophat2',
                basic=False)

            if not os.path.exists(self.transcriptome_gtf):
                raise Exception(f'The reference transcriptome GTF file {self.transcriptome_gtf!r} does not exist.')

        if not self.transcriptome_gtf:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_index' or 'transcriptome_gtf' "
                            f"configuration option.")

        if not self.library_type:
            raise Exception(f"A {self.name!s} requires a 'library_type' configuration option.")

        library_type_tuple = ('fr-unstranded', 'fr-firststrand', 'fr-secondstrand')
        if self.library_type not in library_type_tuple:
            raise Exception(f"The 'library_type' configuration option {self.library_type!r} "
                            f"is not a member of {library_type_tuple!r}.")
        if self.aligner:
            # Defaults to '',
            aligner_tuple = ('hisat2', 'star', 'tophat2')
            if self.aligner not in aligner_tuple:
                raise Exception(f"The 'aligner' configuration option {self.aligner!r} "
                                f"is not a member of {aligner_tuple!r}.")

        # Check for a UCSC bigGenePred.as AutoSQL file.

        if not self.ucsc_autosql_big_gene_prediction:
            self.ucsc_autosql_big_gene_prediction = os.path.join(StandardFilePath.get_resource_ucsc(), 'bigGenePred.as')

        if not (self.ucsc_autosql_big_gene_prediction and os.path.exists(self.ucsc_autosql_big_gene_prediction)):
            raise Exception(f"The {self.name!r} requires a valid 'ucsc_autosql_big_gene_prediction' option.")

        # Read configuration options.

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

        runnable_run_cufflinks_list: list[ConsecutiveRunnable] = list()

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            module_logger.debug('Sample.name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            if self.aligner == 'hisat2':
                run_cufflinks_dependency = Hisat2.get_prefix_sample(
                    sample_name=sample.name)
            elif self.aligner == 'star':
                run_cufflinks_dependency = Star.get_prefix_sample(
                    sample_name=sample.name)
            else:  # tophat2 is the default case.
                # Create a Tophat Runnable per Sample.name.

                file_path_run_tophat = self.get_file_path_run_tophat(sample_name=sample.name)

                runnable_run_tophat = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_run_tophat(sample_name=sample.name),
                        working_directory=self.genome_directory))
                executable_run_tophat = self.set_stage_runnable(
                    stage=stage_run_tophat,
                    runnable=runnable_run_tophat)

                # NOTE: The rnaseq_run_tophat stage does not follow the standard.
                # Instead of adding the bsf.process.RunnableStep to the bsf.procedure.ConsecutiveRunnable,
                # it gets serialised into a separate pickler file.
                # Create a new Tophat bsf.process.RunnableStep.

                runnable_step = RunnableStep(
                    name='tophat2',
                    program='tophat2')

                # Read configuration section [bsf.analyses.rnaseq.Tuxedo.tophat2]
                self.set_runnable_step_configuration(runnable_step=runnable_step)

                # Set tophat options.

                runnable_step.add_option_long(key='GTF', value=self.transcriptome_gtf)

                if self.transcriptome_index:
                    runnable_step.add_option_long(key='transcriptome-index', value=self.transcriptome_index)

                runnable_step.add_option_long(key='output-dir', value=file_path_run_tophat.output_directory)

                runnable_step.add_option_long(key='num-threads', value=str(stage_run_tophat.threads))

                # TODO: These really are properties of the Reads, PairedReads or Sample objects rather than an Analysis.
                if self.insert_size and self.read_length:
                    runnable_step.add_option_long(
                        key='mate-inner-dist',
                        value=str(self.insert_size - 2 * self.read_length))

                if self.insert_size_sd:
                    runnable_step.add_option_long(key='mate-std-dev', value=str(self.insert_size_sd))

                if self.library_type:
                    runnable_step.add_option_long(key='library-type', value=self.library_type)

                # The TopHat coverage search finds additional 'GT-AG' introns, but is only recommended for
                # short reads (< 45 bp) and small read numbers (<= 10 M).
                # TODO: This option should possibly become configurable per sample.
                runnable_step.add_switch_long(key='no-coverage-search')

                # Set rnaseq_tophat arguments.

                runnable_step.arguments.append(self.genome_index_path)

                # Set rnaseq_tophat arguments for reads1 and reads2.

                reads_1_file_path_list: list[str] = list()
                reads_2_file_path_list: list[str] = list()

                for paired_reads_name in sorted(paired_reads_dict):
                    for paired_reads in paired_reads_dict[paired_reads_name]:
                        module_logger.debug('PairedReads name: %r', paired_reads.get_name())

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

                module_logger.log(logging.DEBUG - 1, 'Tophat Executable: %r', runnable_step)

                pickler_dict_run_tophat = {
                    'prefix': stage_run_tophat.name,
                    'replicate_key': sample.name,
                    'runnable_step': runnable_step,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    stage_run_tophat.name + '_' + sample.name + '_run_tophat.pkl')

                with open(file=pickler_path, mode='wb') as output_binary_io:
                    pickler = pickle.Pickler(file=output_binary_io, protocol=pickle.HIGHEST_PROTOCOL)
                    pickler.dump(pickler_dict_run_tophat)

                runnable_step = RunnableStep(name='run_tophat2', program='bsf_run_rnaseq_tophat.py')
                runnable_run_tophat.add_runnable_step(runnable_step=runnable_step)
                runnable_step.add_option_long(key='pickler_path', value=pickler_path)

                # Create a process_tophat Runnable per sample.name.

                runnable_process_tophat = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_process_tophat(sample_name=sample.name),
                        working_directory=self.genome_directory))
                executable_process_tophat = self.set_stage_runnable(
                    stage=stage_process_tophat,
                    runnable=runnable_process_tophat)
                # Set dependencies on previous Runnable or bsf.process.Executable objects.
                executable_process_tophat.dependencies.append(executable_run_tophat.name)

                # TODO: Switch from an external Bash script to a set of Runnable and RunnableStep objects.
                #  Since the Bash script includes Perl code to reset the BED score field to 0, rather than
                #  re-scale it properly, it would be good to write a new bsf.runnables.process_tophat module
                #  to implement this in Python code.

                runnable_step = RunnableStep(
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

            file_path_cufflinks = self.get_file_path_run_cufflinks(sample_name=sample.name)

            runnable_run_cufflinks = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_run_cufflinks(sample_name=sample.name),
                    working_directory=self.genome_directory))
            # Set dependencies for subsequent Runnable or bsf.process.Executable objects.
            runnable_run_cufflinks_list.append(runnable_run_cufflinks)
            executable_run_cufflinks = self.set_stage_runnable(
                stage=stage_run_cufflinks,
                runnable=runnable_run_cufflinks)
            # Set dependencies on previous Runnable or bsf.process.Executable objects.
            executable_run_cufflinks.dependencies.append(run_cufflinks_dependency)

            # Create a new Cufflinks bsf.process.RunnableStep.

            runnable_step = RunnableStep(
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
                    value=self.transcriptome_gtf)
            else:
                # --GTF quantify against reference transcript annotations [NULL]
                runnable_step.add_option_long(
                    key='GTF',
                    value=self.transcriptome_gtf)
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
                    Hisat2.get_file_path_sample(sample_name=sample.name).sample_bam)
            elif self.aligner == 'star':
                runnable_step.arguments.append(
                    Star.get_file_path_sample(sample_name=sample.name).sample_bam)
            else:
                runnable_step.arguments.append(
                    self.get_file_path_run_tophat(sample_name=sample.name).accepted_hits_bam)

            # Convert the resulting transcripts GTF file into a UCSC genePred file.

            runnable_step = RunnableStep(
                name='gtf_to_gp',
                program='gtfToGenePred')
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_short(key='genePredExt')
            runnable_step.arguments.append(file_path_cufflinks.transcripts_gtf)
            runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)

            # Convert the UCSC genePred into a UCSC bigGenePred file.

            runnable_step = RunnableStep(
                name='gp_to_bgp',
                program='genePredToBigGenePred',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_gene_prediction,
                ])
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)
            runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)

            # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

            runnable_step = RunnableStep(
                name='bed_sort',
                program='bedSort',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_big_gene_prediction,
                ])
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)
            runnable_step.arguments.append(file_path_cufflinks.temporary_sorted_tsv)

            # Run bedtools slop on the sorted UCSC bigGenePred file to constrain to chromosome coordinates.

            runnable_step = RunnableStep(
                name='bedtools_slop',
                program='bedtools',
                sub_command=Command(program='slop'),
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_sorted_tsv,
                ],
                stdout=ConnectorFile(file_path=file_path_cufflinks.temporary_slopped_tsv, file_mode='wt'))
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.sub_command.add_option_short(key='b', value='0')
            runnable_step.sub_command.add_option_short(key='i', value=file_path_cufflinks.temporary_sorted_tsv)
            runnable_step.sub_command.add_option_short(key='g', value=self.genome_sizes_path)

            # Since bedtools slop, at least in version v2.27.1, looses the last tab character,
            # it has to be put back with sed, before UCSC bedToBigBed can run.

            runnable_step = RunnableStep(
                name='sed',
                program='sed',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_slopped_tsv,
                ],
                stdout=ConnectorFile(file_path=file_path_cufflinks.temporary_fixed_tsv, file_mode='wt'))
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_short(key='e', value='s/$/\\t/')

            runnable_step.arguments.append(file_path_cufflinks.temporary_slopped_tsv)

            # Convert the sorted UCSC bigGenePred into a bigBed file.

            runnable_step = RunnableStep(
                name='bgp_to_bb',
                program='bedToBigBed',
                obsolete_file_path_list=[
                    file_path_cufflinks.temporary_fixed_tsv,
                ])
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_pair_short(key='as', value=self.ucsc_autosql_big_gene_prediction)
            runnable_step.add_switch_short(key='tab')
            runnable_step.add_option_pair_short(key='type', value='bed12+8')
            runnable_step.arguments.append(file_path_cufflinks.temporary_fixed_tsv)
            runnable_step.arguments.append(self.genome_sizes_path)
            runnable_step.arguments.append(file_path_cufflinks.transcripts_bb)

            # Add a symbolic link for the skipped.gtf file, which includes a sample name prefix.
            runnable_step = RunnableStepLink(
                name='link_skipped_gtf',
                source_path=file_path_cufflinks.skipped_gtf_link_source,
                target_path=file_path_cufflinks.skipped_gtf_link_target)
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # Add a symbolic link for the transcripts bigBed file, which includes a sample name prefix.
            runnable_step = RunnableStepLink(
                name='link_transcripts_bb',
                source_path=file_path_cufflinks.transcripts_bb_link_source,
                target_path=file_path_cufflinks.transcripts_bb_link_target)
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # Add a symbolic link for the transcripts GTF file, which includes a sample name prefix.
            runnable_step = RunnableStepLink(
                name='link_transcripts_gtf',
                source_path=file_path_cufflinks.transcripts_gtf_link_source,
                target_path=file_path_cufflinks.transcripts_gtf_link_target)
            runnable_run_cufflinks.add_runnable_step(runnable_step=runnable_step)

        # Create one process_cufflinks bsf.process.Executable to process all subdirectories.

        if len(runnable_run_cufflinks_list):
            runnable_process_cufflinks = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_process_cufflinks(),
                    working_directory=self.genome_directory))
            executable_process_cufflinks = self.set_stage_runnable(
                stage=stage_process_cufflinks,
                runnable=runnable_process_cufflinks)
            # Set dependencies on previous Runnable or bsf.process.Executable objects.
            for runnable_run_cufflinks in runnable_run_cufflinks_list:
                executable_process_cufflinks.dependencies.append(runnable_run_cufflinks.name)

            runnable_step = RunnableStep(
                name='process_cufflinks',
                program='bsf_rnaseq_process_cufflinks.R')
            runnable_process_cufflinks.add_runnable_step(runnable_step=runnable_step)

            # Read configuration section [bsf.analyses.rnaseq.Tuxedo.process_cufflinks]
            self.set_runnable_step_configuration(runnable_step=runnable_step)

            runnable_step.add_option_long(
                key='gtf-reference',
                value=self.transcriptome_gtf)
            runnable_step.add_option_long(
                key='genome-version',
                value=self.genome_version)

            runnable_step = RunnableStep(
                name='tophat_summary',
                program='bsf_rnaseq_tophat_summary.R')
            runnable_process_cufflinks.add_runnable_step(runnable_step=runnable_step)

        # The Cuffmerge process generates temporary files in the working directory that
        # have the same name for each comparison. If more than one Cuffmerge process runs at the same time,
        # contention occurs where more than one process writes to the same file and one process deletes files
        # expected to exist by another process.
        # Circumvent such a situation by introducing dependencies on previous Cuffmerge processes. Sigh.
        # TODO: Report this to the Cufflinks author.
        executable_cuffmerge_dict: dict[str, Executable] = dict()

        for comparison_name in sorted(self._comparison_dict):
            module_logger.debug('Comparison name: %r', comparison_name)

            sample_group_list = self._comparison_dict[comparison_name]

            if len(sample_group_list) == 0:
                continue

            # Process rnaseq_cuffmerge and rnaseq_cuffdiff arguments in parallel.
            # Check that the comparison contains at least one sample group.

            cuffdiff_cuffnorm_abundances_dict: dict[str, list[str]] = dict()
            cuffdiff_cuffnorm_alignments_dict: dict[str, list[str]] = dict()
            cuffdiff_cuffnorm_dependencies: list[str] = list()
            cuffmerge_cuffnorm_submit: bool = len(sample_group_list) >= 1
            cuffmerge_transcript_gtf_list: list[str] = list()

            # TODO: Should the comparison prefix also include the project name or number?
            prefix_run_cuffmerge = self.get_prefix_run_cuffmerge(comparison_name=comparison_name)

            file_path_cuffmerge = self.get_file_path_cuffmerge(comparison_name=comparison_name)

            runnable_run_cuffmerge = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=prefix_run_cuffmerge,
                    working_directory=self.genome_directory))
            executable_run_cuffmerge = self.set_stage_runnable(
                stage=stage_run_cuffmerge,
                runnable=runnable_run_cuffmerge)
            # Submit the bsf.process.Executable if the status file AND the sample group list above supports it.
            executable_run_cuffmerge.submit &= cuffmerge_cuffnorm_submit
            # Set a dependency on all other Cuffmerge process to avoid file contention.
            executable_cuffmerge_dict[prefix_run_cuffmerge] = executable_run_cuffmerge

            if self.novel_transcripts:
                # Create a new Cuffmerge bsf.process.RunnableStep.

                runnable_step_cuffmerge = RunnableStep(
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
                    value=self.transcriptome_gtf)
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

                runnable_step = RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=file_path_cuffmerge.output_directory)
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

                runnable_step = RunnableStepCopy(
                    name='copy',
                    source_path=self.transcriptome_gtf,
                    target_path=file_path_cuffmerge.merged_gtf)
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

                # Run cuffcompare in a self-comparison mode to get 'tss_id' and 'p_id' attributes populated.

                runnable_step = RunnableStep(
                    name='cuffcompare',
                    program='cuffcompare',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.merged_gtf,
                    ])
                runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_switch_short(key='C')  # include 'contained' transcripts
                runnable_step.add_switch_short(key='G')  # generic GFF input fields, (i.e., not a Cufflinks GTF)
                runnable_step.add_option_short(key='o', value=file_path_cuffmerge.cuffcompare_prefix)
                runnable_step.add_option_short(key='r', value=self.transcriptome_gtf)  # reference GTF
                runnable_step.add_option_short(key='s', value=self.genome_fasta_path)  # reference sequence
                runnable_step.arguments.append(file_path_cuffmerge.merged_gtf)

                # Move 'cuffcmp.combined.gtf' to 'merged.gtf' and delete obsolete files.
                runnable_step = RunnableStepMove(
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

            runnable_step = RunnableStep(
                name='gtf_to_gp',
                program='gtfToGenePred')
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_switch_short(key='genePredExt')
            runnable_step.arguments.append(file_path_cuffmerge.merged_gtf)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)

            # Convert the UCSC genePred into a UCSC bigGenePred file.

            runnable_step = RunnableStep(
                name='gp_to_bgp',
                program='genePredToBigGenePred',
                obsolete_file_path_list=[
                    file_path_cuffmerge.temporary_gene_prediction,
                ])
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)

            # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

            runnable_step = RunnableStep(
                name='bed_sort',
                program='bedSort',
                obsolete_file_path_list=[
                    file_path_cuffmerge.temporary_big_gene_prediction,
                ])
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)

            # Convert the sorted UCSC bigGenePred into a bigBed file.

            runnable_step = RunnableStep(
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
            runnable_step = RunnableStepLink(
                name='link_merged_bb',
                source_path=file_path_cuffmerge.merged_bb_link_source,
                target_path=file_path_cuffmerge.merged_bb_link_target)
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            # Add a symbolic link for the merged GTF file, that includes a comparison name prefix.
            runnable_step = RunnableStepLink(
                name='link_merged_gtf',
                source_path=file_path_cuffmerge.merged_gtf_link_source,
                target_path=file_path_cuffmerge.merged_gtf_link_target)
            runnable_run_cuffmerge.add_runnable_step(runnable_step=runnable_step)

            file_path_monocle = self.get_file_path_monocle(comparison_name=comparison_name)

            monocle_annotation_sheet = AnnotationSheet(
                file_path=os.path.join(self.genome_directory, file_path_monocle.annotation_tsv),
                file_type='excel-tab',
                header=True)

            for sample_group in sample_group_list:
                module_logger.debug('SampleGroup.name: %r', sample_group.name)

                per_group_abundances_list: list[str] = list()
                per_group_alignments_list: list[str] = list()

                for sample in sample_group.sample_list:
                    module_logger.debug('Sample.name: %r', sample.name)

                    paired_reads_dict = sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping,
                        exclude=True)

                    if not paired_reads_dict:
                        # Skip Sample objects, which PairedReads objects have all been excluded.
                        continue

                    # Add the Cufflinks assembled transcripts GTF to the Cuffmerge manifest.
                    cuffmerge_transcript_gtf_list.append(
                        os.path.join('_'.join(('rnaseq_cufflinks', sample.name)), 'transcripts.gtf') + '\n')

                    # Wait for all Cufflinks instances to finish, before Cuffmerge can run.

                    executable_run_cuffmerge.dependencies.append(
                        self.get_prefix_run_cufflinks(sample_name=sample.name))

                    # Create a Cuffquant Runnable per comparison (comparison_name) and
                    # Sample.name on the basis of the Cuffmerge GTF file.

                    file_path_run_cuffquant = self.get_file_path_cuffquant(
                        comparison_name=comparison_name,
                        sample_name=sample.name)

                    runnable_run_cuffquant = self.add_runnable(
                        runnable=ConsecutiveRunnable(
                            name=self.get_prefix_run_cuffquant(
                                comparison_name=comparison_name,
                                sample_name=sample.name),
                            working_directory=self.genome_directory))
                    executable_run_cuffquant = self.set_stage_runnable(
                        stage=stage_run_cuffquant,
                        runnable=runnable_run_cuffquant)
                    # Each Cuffquant process depends on Cuffmerge.
                    executable_run_cuffquant.dependencies.append(executable_run_cuffmerge.name)

                    # Create a new cuffquant bsf.process.RunnableStep.

                    runnable_step_cuffquant = RunnableStep(
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
                            Hisat2.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                        per_group_alignments_list.append(
                            Hisat2.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                    elif self.aligner == 'star':
                        runnable_step_cuffquant.arguments.append(
                            Star.get_file_path_sample(
                                sample_name=sample.name).sample_bam)
                        per_group_alignments_list.append(
                            Star.get_file_path_sample(
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
                    #  Depending on the replicate_grouping instance variable, the abundance file can be on
                    #  the read_group or sample level, while Monocle annotation will always be on the sample level.
                    monocle_row_dict: dict[str, str] = {
                        'file': file_path_run_cuffquant.abundances,
                        # 'sample_name' is used by Monocle in the plot_cell_clusters() function internally.
                        'original_name': sample.name
                    }
                    # Set additional columns from the Sample Annotation Sheet prefixed with 'Sample Monocle *'.
                    for annotation_key in filter(
                            lambda x: x.startswith('Monocle '), sample.annotation_dict.keys()):
                        monocle_row_dict[annotation_key[8:]] = sample.annotation_dict[annotation_key][0]

                    monocle_annotation_sheet.row_dict_list.append(monocle_row_dict)

                cuffdiff_cuffnorm_abundances_dict[sample_group.name] = per_group_abundances_list
                cuffdiff_cuffnorm_alignments_dict[sample_group.name] = per_group_alignments_list

            if self.novel_transcripts:
                # Write a Cuffmerge assembly manifest file to merge all transcriptome GTF files of each Sample object.
                # This requires an absolute path, because the working directory is not set at the stage of
                # job submission.
                with open(
                        file=os.path.join(self.genome_directory, file_path_cuffmerge.assembly_txt),
                        mode='wt') as output_text_io:
                    output_text_io.writelines(cuffmerge_transcript_gtf_list)

            if len(self._comparison_dict[comparison_name]) >= 2:
                # Create a Cuffnorm Runnable per comparison, if there are at least two SampleGroup objects.

                file_path_run_cuffnorm = self.get_file_path_cuffnorm(comparison_name=comparison_name)

                runnable_run_cuffnorm = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_run_cuffnorm(comparison_name=comparison_name),
                        working_directory=self.genome_directory))
                executable_run_cuffnorm = self.set_stage_runnable(
                    stage=stage_run_cuffnorm,
                    runnable=runnable_run_cuffnorm)
                # Submit the bsf.process.Executable if the status file AND the sample group list above supports it.
                executable_run_cuffnorm.submit &= cuffmerge_cuffnorm_submit
                executable_run_cuffnorm.dependencies.extend(cuffdiff_cuffnorm_dependencies)

                # Create a new Cuffnorm bsf.process.RunnableStep.

                runnable_step_cuffnorm = RunnableStep(
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

                # Add an abundance annotation file as second Cuffnorm argument.
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

                runnable_run_cuffdiff = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_run_cuffdiff(comparison_name=comparison_name),
                        working_directory=self.genome_directory))
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
                runnable_step = RunnableStepSetEnvironment(
                    name='set_environment',
                    key='LANG',
                    value='C')
                runnable_run_cuffdiff.add_runnable_step(runnable_step=runnable_step)

                # Create a new Cuffdiff bsf.process.RunnableStep.

                runnable_step_run_cuffdiff = RunnableStep(
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

                runnable_process_cuffdiff = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_process_cuffdiff(comparison_name=comparison_name),
                        working_directory=self.genome_directory))
                executable_process_cuffdiff = self.set_stage_runnable(
                    stage=stage_process_cuffdiff,
                    runnable=runnable_process_cuffdiff)
                executable_process_cuffdiff.dependencies.append(executable_run_cuffdiff.name)

                runnable_step_process_cuffdiff = RunnableStep(
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
                    value=self.transcriptome_gtf)
                runnable_step_process_cuffdiff.add_option_long(
                    key='genome-version',
                    value=self.genome_version)

        # Finally, set dependencies on all other Cuffmerge bsf.process.Executable objects to avoid file contention.
        for prefix_run_cuffmerge in executable_cuffmerge_dict:
            for executable_cuffmerge in executable_cuffmerge_dict.values():
                if prefix_run_cuffmerge != executable_cuffmerge.name:
                    executable_cuffmerge.dependencies.append(prefix_run_cuffmerge)

        return

    def report(self) -> None:
        """Create a :literal:`XHTML 1.0` report and a :emphasis:`UCSC Genome Browser Track Hub`.
        """

        def report_html() -> None:
            """Private function to create an :literal:`XHTML 1.0` report.
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list: list[str] = list()

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
                str_list.append('<a href="' + Hisat2.prefix + '_report.html">')
                str_list.append(self.project_name + ' ' + Hisat2.name)
                str_list.append('</a> report for quality plots and ')
                str_list.append('a link to alignment visualisation in the UCSC Genome Browser.\n')
                str_list.append('</p>\n')
                str_list.append('\n')
            elif self.aligner == 'star':
                str_list.append('<p id="star">')
                str_list.append('<a href="https://github.com/alexdobin/STAR">STAR</a> ')
                str_list.append('aligns RNA-seq reads to a reference genome in order to identify ')
                str_list.append('exon-exon splice junctions. ')
                # str_list.append('<br />\n')
                str_list.append('Please see the ')
                str_list.append('<a href="' + Star.prefix + '_report.html">')
                str_list.append(self.project_name + ' ' + Star.name)
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
                module_logger.debug('Sample.name: %r', sample.name)
                module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                if self.aligner == 'hisat2' or self.aligner == 'star':
                    if self.aligner == 'hisat2':
                        file_path_aligner_sample = Hisat2.get_file_path_sample(
                            sample_name=sample.name)
                    elif self.aligner == 'star':
                        file_path_aligner_sample = Star.get_file_path_sample(
                            sample_name=sample.name)
                    else:
                        raise Exception('Program error.')

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
                    # The Cufflinks tool creates genes.fpkm_tracking, isoforms.fpkm_tracking,
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

                    # TODO: The aligned BAM and BAI files and the unaligned BAM file are currently non-standard.
                    #  The files have a 'rnaseq_tophat_' prefix, but are in the 'rnaseq_cufflinks_' directory.
                    #  This will be resolved when the process_tophat step gets re-engineered.

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
            str_list.append('but sometimes, technical constraints (i.e., memory requirements) require setting up ')
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
                    sample_pair_sheet: TuxedoSamplePairSheet = TuxedoSamplePairSheet.from_file_path(
                        file_path=sample_pair_path)

                    for row_dict in sample_pair_sheet.row_dict_list:
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
                    sample_pair_sheet: TuxedoSamplePairSheet = TuxedoSamplePairSheet.from_file_path(
                        file_path=sample_pair_path)

                    for row_dict in sample_pair_sheet.row_dict_list:
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

            return

        def report_hub() -> None:
            """Private function to create a :emphasis:`UCSC Genome Browser Track Hub`.
            """

            str_list: list[str] = list()

            # Group via UCSC super tracks.

            str_list.append('track Alignments\n')
            str_list.append('shortLabel Alignments\n')
            str_list.append('longLabel TopHat RNA-seq read alignments\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Assemblies\n')
            str_list.append('shortLabel Assemblies\n')
            str_list.append('longLabel Cuffmerge transcript structures\n')
            str_list.append('visibility full\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Coverage\n')
            str_list.append('shortLabel Coverage\n')
            str_list.append('longLabel TopHat RNA-seq alignment coverage\n')
            str_list.append('visibility full\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Deletions\n')
            str_list.append('shortLabel Deletions\n')
            str_list.append('longLabel TopHat RNA-seq deletions\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Insertions\n')
            str_list.append('shortLabel Insertions\n')
            str_list.append('longLabel TopHat RNA-seq insertions\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Junctions\n')
            str_list.append('shortLabel Junctions\n')
            str_list.append('longLabel TopHat RNA-seq splice junctions\n')
            str_list.append('visibility show\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Transcripts\n')
            str_list.append('shortLabel Transcripts\n')
            str_list.append('longLabel Cufflinks transcript structures\n')
            str_list.append('visibility show\n')
            str_list.append('superTrack on\n')
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
                    # Alignments track
                    #
                    # Common track settings
                    str_list.append('track ' + sample.name + '_alignments\n')
                    str_list.append('type bam\n')
                    str_list.append('shortLabel ' + sample.name + '_alignments\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq read alignments\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/accepted_hits.bam\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility dense\n')
                    # Common optional track settings
                    str_list.append('color 0,0,0\n')
                    # bam/cram - Compressed Sequence Alignment track settings
                    # Supertrack settings
                    str_list.append('parent Alignments\n')
                    str_list.append('\n')

                    # Coverage track
                    #
                    # Common track settings
                    str_list.append('track ' + sample.name + '_coverage\n')
                    # TODO: The bigWig type must declare the expected signal range.
                    #  The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                    str_list.append('type bigWig\n')
                    str_list.append('shortLabel ' + sample.name + '_coverage\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq alignment coverage\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/accepted_hits.bw\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility full\n')
                    # Common optional track settings
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
                    # Supertrack settings
                    str_list.append('parent Coverage\n')
                    str_list.append('\n')

                    # Deletions track
                    #
                    # Common track settings
                    str_list.append('track ' + sample.name + '_deletions\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_deletions\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq deletions\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/deletions.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')
                    # Common optional track settings
                    str_list.append('color 0,0,0\n')
                    # bigBed - Item or region track settings
                    # Supertrack settings
                    str_list.append('parent Deletions\n')
                    str_list.append('\n')

                    # Insertions track
                    #
                    # Common track settings
                    str_list.append('track insertions_' + sample.name + '\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_insertions\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq insertions\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/insertions.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')
                    # Common optional track settings
                    str_list.append('color 0,0,0\n')
                    # bigBed - Item or region track settings
                    # Supertrack settings
                    str_list.append('parent Insertions\n')
                    str_list.append('\n')

                    # Junctions track
                    #
                    # Common track settings
                    str_list.append('track ' + sample.name + '_junctions\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_junctions\n')
                    str_list.append('longLabel ' + sample.name + ' TopHat RNA-seq splice junctions\n')
                    str_list.append('bigDataUrl rnaseq_tophat_' + sample.name + '/junctions.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility pack\n')
                    # Common optional track settings
                    str_list.append('color 0,0,0\n')
                    # bigBed - Item or region track settings
                    # Supertrack settings
                    str_list.append('parent Junctions\n')
                    str_list.append('\n')

                    # Transcripts track
                    #
                    # Common track settings
                    str_list.append('track ' + sample.name + '_transcripts\n')
                    str_list.append('type bigGenePred\n')
                    str_list.append('shortLabel ' + sample.name + '_transcripts\n')
                    str_list.append('longLabel ' + sample.name + ' Cufflinks transcript assembly\n')
                    str_list.append('bigDataUrl rnaseq_cufflinks_' + sample.name + '/transcripts.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')
                    # Common optional track settings
                    str_list.append('color 0,0,0\n')
                    # bigGenePred - Gene Annotations settings
                    # Supertrack settings
                    str_list.append('parent Transcripts\n')
                    str_list.append('\n')

            # Comparison-specific tracks

            for comparison_name in sorted(self._comparison_dict):
                # Assemblies track
                #
                # Common track settings
                str_list.append('track ' + comparison_name + '_assembly\n')
                str_list.append('type bigGenePred\n')
                str_list.append('shortLabel ' + comparison_name + '_assembly\n')
                str_list.append('longLabel ' + comparison_name + ' Cufflinks transcript assembly\n')
                str_list.append('bigDataUrl rnaseq_cuffmerge_' + comparison_name + '/merged.bb\n')
                # str_list.append('html ...\n')
                str_list.append('visibility pack\n')
                # Common optional track settings
                str_list.append('color 0,0,0\n')
                # bigGenePred - Gene Annotations settings
                # Supertrack settings
                str_list.append('parent Assemblies\n')
                str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return


class FilePathDESeq(FilePath):
    """The :py:class:`bsf.analyses.rnaseq.FilePathDESeq` class models files in a comparison-specific DESeq directory.
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.FilePathDESeq` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathDESeq, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class DESeq(Analysis):
    """The :py:class:`bsf.analyses.rnaseq.DESeq` class models an RNA-seq analysis based on the
    `Bioconductor <https://bioconductor.org/>`_
    `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ package.

    :ivar transcriptome_version: A transcriptome version.
    :type transcriptome_version: str | None
    :ivar transcriptome_granges: A Bioconductor :literal:`GenomicRanges::GRanges` transcriptome resource file path.
    :type transcriptome_granges: str | None
    :ivar comparison_path: A comparison file path.
    :type comparison_path: str | None
    :ivar contrast_path: A contrast file path.
    :type contrast_path: str | None
    :ivar only_counting: Request creating count matrices rather than running a full analysis.
    :type only_counting: bool | None
    """

    name = 'DESeq RNA-seq Analysis'
    prefix = '_'.join(('rnaseq', 'deseq'))

    @classmethod
    def get_stage_name_analysis(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'analysis'))

    @classmethod
    def get_stage_name_results(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'results'))

    @classmethod
    def get_prefix_analysis(cls, design_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param design_name: A design name.
        :type design_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_analysis(), design_name))

    @classmethod
    def get_prefix_results(cls, design_name: str) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param design_name: A design name.
        :type design_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_results(), design_name))

    def __init__(
            self,
            configuration: Optional[Configuration] = None,
            project_name: Optional[str] = None,
            genome_version: Optional[str] = None,
            input_directory: Optional[str] = None,
            output_directory: Optional[str] = None,
            project_directory: Optional[str] = None,
            genome_directory: Optional[str] = None,
            report_style_path: Optional[str] = None,
            report_header_path: Optional[str] = None,
            report_footer_path: Optional[str] = None,
            e_mail: Optional[str] = None,
            stage_list: Optional[list[Stage]] = None,
            collection: Optional[Collection] = None,
            sample_list: Optional[list[Sample]] = None,
            transcriptome_version: Optional[str] = None,
            transcriptome_granges: Optional[str] = None,
            comparison_path: Optional[str] = None,
            contrast_path: Optional[str] = None,
            only_counting: Optional[bool] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.rnaseq.DESeq` object.

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
        :param report_style_path: A report style :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: A report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: A report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a :emphasis:`UCSC Genome Browser Track Hub`.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str | None
        :param transcriptome_granges: A Bioconductor :literal:`GenomicRanges::GRanges` transcriptome resource file path.
        :type transcriptome_granges: str | None
        :param comparison_path: A comparison file path.
        :type comparison_path: str | None
        :param contrast_path: A contrast file path.
        :type contrast_path: str | None
        :param only_counting: Request creating count matrices rather than running a full analysis.
        :type only_counting: bool | None
        """

        super(DESeq, self).__init__(
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
        self.transcriptome_granges = transcriptome_granges
        self.comparison_path = comparison_path
        self.contrast_path = contrast_path
        self.only_counting = only_counting

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.rnaseq.DESeq` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """

        super(DESeq, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_granges'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_granges = configuration.config_parser.get(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

        option = 'ctr_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.contrast_path = configuration.config_parser.get(section=section, option=option)

        option = 'only_counting'
        if configuration.config_parser.has_option(section=section, option=option):
            self.only_counting = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.rnaseq.DESeq` object.
        """

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # DESeq requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_version' configuration option.")

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception(f"A {self.name!s} requires a valid 'transcriptome_version' configuration option.")

        # Get the annotation sheets before calling the run() method of the Analysis super-class.
        # If file paths were not provided, try to find them in the current directory.
        # The complication is that either the Tuxedo.prefix or the DESeq.prefix could be used.

        # Get the sample annotation sheet.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception(f'The sample annotation sheet {self.sas_file!r} does not exist.')
        else:
            self.sas_file = self.get_annotation_file(
                prefix_list=[DESeq.prefix, Tuxedo.prefix],
                suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation sheet in the current working directory.')

        # Get the design annotation sheet.

        if self.comparison_path:
            self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)
            if not os.path.exists(self.comparison_path):
                raise Exception(f'The comparison (design) annotation sheet {self.comparison_path!r} does not exist.')
        else:
            self.comparison_path = self.get_annotation_file(
                prefix_list=[DESeq.prefix, Tuxedo.prefix],
                suffix='designs.csv')
            if not self.comparison_path:
                raise Exception('No suitable comparison (design) annotation sheet in the current working directory.')

        # Get the contrast annotation sheet.

        if self.contrast_path:
            self.contrast_path = self.configuration.get_absolute_path(file_path=self.contrast_path)
        else:
            self.contrast_path = self.get_annotation_file(
                prefix_list=[DESeq.prefix, Tuxedo.prefix],
                suffix='contrasts.csv')

        super(DESeq, self).run()

        # Get the transcriptome GTF file.

        if not self.transcriptome_granges:
            self.transcriptome_granges = StandardFilePath.get_resource_transcriptome_granges(
                transcriptome_version=self.transcriptome_version)

        stage_analysis = self.get_stage(name=self.get_stage_name_analysis())
        stage_results = self.get_stage(name=self.get_stage_name_results())

        # For DESeq, all samples need adding to the Analysis regardless.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        # Read the designs (comparison) file.

        design_sheet: AnnotationSheet = AnnotationSheet.from_file_path(
            file_path=self.comparison_path,
            file_type='excel',
            name='DESeq Design Table')

        # Read the contrasts file.

        contrast_sheet: Optional[AnnotationSheet] = None

        if self.contrast_path:
            contrast_sheet = AnnotationSheet.from_file_path(
                file_path=self.contrast_path,
                file_type='excel',
                name='DESeq Contrast Table')

        # TODO: Adjust by introducing a new class RNASeqComparisonSheet(AnnotationSheet) in this module?
        for design_row_dict in design_sheet.row_dict_list:
            design_name = design_row_dict['design']

            prefix = '_'.join((self.prefix, design_name))

            comparison_directory = os.path.join(self.genome_directory, prefix)

            if not os.path.isdir(comparison_directory):
                try:
                    os.makedirs(comparison_directory)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise exception

            annotation_sheet = AnnotationSheet(
                file_path=os.path.join(comparison_directory, prefix + '_samples.tsv'),
                file_type='excel-tab',
                # file_type='excel',
                name='DESeq Sample Annotation',
                header=True)

            # Sort the Python list of Sample objects by Sample.name.

            self.sample_list.sort(key=lambda item: item.name)

            for sample in self.sample_list:
                module_logger.debug('Sample.name: %r', sample.name)
                module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                file_path_star_sample = Star.get_file_path_sample(
                    sample_name=sample.name)
                file_path_kallisto_sample = Kallisto.get_file_path_sample(
                    sample_name=sample.name)
                row_dict: dict[str, str] = {
                    'bam_path': file_path_star_sample.sample_bam,
                    'bai_path': file_path_star_sample.sample_bai,
                    'ah5_path': file_path_kallisto_sample.abundance_h5,
                }

                # If the sample annotation sheet contains a "DESeq sample" variable,
                # technical replicates need collapsing. Set the original bsf.ngs.Sample.name as
                # 'run', the DESeq sample name will be filled in from the annotation dict below.
                if "DESeq sample" in sample.annotation_dict:
                    row_dict['run'] = sample.name
                else:
                    row_dict['sample'] = sample.name

                # Set additional columns from the Sample Annotation Sheet prefixed with 'Sample DESeq *'.
                for annotation_key in filter(lambda x: x.startswith('DESeq '), sample.annotation_dict.keys()):
                    row_dict[annotation_key[6:]] = sample.annotation_dict[annotation_key][0]

                # Post process the UMIs variable to use 'FALSE' or 'TRUE' which can be parsed by R.
                if 'UMIs' in row_dict:
                    if annotation_sheet.get_boolean(row_dict=row_dict, key='UMIs'):
                        row_dict['UMIs'] = 'TRUE'
                    else:
                        row_dict['UMIs'] = 'FALSE'
                else:
                    row_dict['UMIs'] = 'FALSE'

                annotation_sheet.row_dict_list.append(row_dict)

            # Write the DESeq Sample Annotation Sheet to disk.
            annotation_sheet.to_file_path(adjust_field_names=True)

            # Re-write the Annotation Sheet objects for DESeq designs and contrasts to new file paths in the
            # genome directory. Convert the CSV configuration tables into TSV analysis tables.
            design_sheet.file_path = os.path.join(comparison_directory, prefix + '_designs.tsv')
            design_sheet.file_type = 'excel-tab'
            design_sheet.to_file_path()

            has_contrasts = False
            if self.contrast_path:
                # Write the full contrast sheet into the design-specific directory.
                contrast_sheet.file_path = os.path.join(comparison_directory, prefix + '_contrasts.tsv')
                contrast_sheet.file_type = 'excel-tab'
                contrast_sheet.to_file_path()

                # Search for design-specific contrasts upon which the results stage can be submitted.
                for contrast_row_dict in contrast_sheet.row_dict_list:
                    if 'Design' in contrast_row_dict \
                            and contrast_row_dict['Design'] == design_name \
                            and not contrast_sheet.get_boolean(row_dict=contrast_row_dict, key='Exclude'):
                        has_contrasts = True
                        break

            runnable_analysis = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_analysis(design_name=design_name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory))
            executable_analysis = self.set_stage_runnable(
                stage=stage_analysis,
                runnable=runnable_analysis)
            # Do not submit the Executable, if the design was excluded.
            if design_sheet.get_boolean(row_dict=design_row_dict, key='exclude'):
                executable_analysis.submit = False

            if self.only_counting:
                runnable_step = RunnableStep(
                    name='analysis',
                    program='bsf_rnaseq_deseq_rse.R')
            else:
                runnable_step = RunnableStep(
                    name='analysis',
                    program='bsf_rnaseq_deseq_analysis.R')
            runnable_analysis.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_long(key='design-name', value=design_name)
            runnable_step.add_option_long(key='transcriptome-path', value=self.transcriptome_granges)
            runnable_step.add_option_long(key='threads', value=str(stage_analysis.threads))

            # Run the results stage, if a contrast annotation sheet is already available and has design-specific rows.
            if has_contrasts:
                runnable_results = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_results(design_name=design_name),
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory))
                executable_results = self.set_stage_runnable(
                    stage=stage_results,
                    runnable=runnable_results)
                executable_results.dependencies.append(runnable_analysis.name)
                # Do not submit the Executable, if the design was excluded.
                if design_sheet.get_boolean(row_dict=design_row_dict, key='exclude'):
                    executable_results.submit = False

                runnable_step = RunnableStep(
                    name='results',
                    program='bsf_rnaseq_deseq_results.R')
                runnable_results.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='design-name', value=design_name)
                runnable_step.add_option_long(key='threads', value=str(stage_results.threads))

                runnable_step = RunnableStep(
                    name='enrichr',
                    program='bsf_rnaseq_deseq_enrichr.R')
                runnable_results.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='design-name', value=design_name)

                runnable_step = RunnableStep(
                    name='heatmap',
                    program='bsf_rnaseq_deseq_heatmap.R')
                runnable_results.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='design-name', value=design_name)

        return

    def report(self) -> None:
        """Create a :literal:`XHTML 1.0` report.
        """

        def image_source_exists(prefix: str, suffix: str) -> bool:
            """Test whether an image source (i.e., path) exists.

            :param prefix: A prefix.
            :type prefix: str
            :param suffix: A suffix.
            :type suffix: str
            :return: :py:const:`True` for existing image sources, :py:const`False` otherwise.
            :rtype: bool
            """
            return os.path.exists(
                os.path.join(
                    self.genome_directory,
                    prefix,
                    prefix + '_' + suffix))

        def report_html() -> None:
            """Private function to create a :literal:`XHTML 1.0` report.
            """
            # Read the design table as a backup, in case the design-specific LRT summary table is not available.
            # Exclude designs if requested.

            design_sheet: AnnotationSheet = AnnotationSheet.from_file_path(
                file_path=self.comparison_path,
                file_type='excel',
                name='DESeq Design Table')

            design_dict: dict[str, dict[str, str]] = {
                value['design']: value for value in design_sheet.row_dict_list if
                not design_sheet.get_boolean(row_dict=value, key='exclude')
            }

            # Create a symbolic link containing the project name and a UUID.
            self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list: list[str] = list()

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.extend(self.get_html_genome(genome_version=self.genome_version))
            str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
            str_list.append('\n')

            str_list.append('<p id="star">')
            str_list.append('<a href="https://github.com/alexdobin/STAR">STAR</a> ')
            str_list.append('aligns RNA-seq reads to a reference genome in order to identify ')
            str_list.append('exon-exon splice junctions.\n')
            # str_list.append('<br />\n')
            str_list.append('Please see the ')
            str_list.append('<a href="' + Star.prefix + '_report.html">')
            str_list.append(self.project_name + ' ' + Star.name)
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

            # Exploratory Analysis Table.

            str_list.append('<h2 id="exploratory_analysis">Exploratory Analysis</h2>\n')
            str_list.append('\n')

            str_list.append('<p>')
            str_list.append('For each design, multi-dimensional scaling (MDS), principal component analysis (PCA) ')
            str_list.append('and heatmap plots are provided for combinations of variables or factors. ')
            str_list.append('Variance by principal component plots show the distribution of the variance for ')
            str_list.append('a maximum of the first hundred components. ')
            str_list.append('Two sets of plots are available, based on data in model-aware or blind mode. ')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="exploratory_analysis_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">MDS Blind</th>\n')
            str_list.append('<th class="left">MDS Model</th>\n')
            str_list.append('<th class="left">PCA Blind</th>\n')
            str_list.append('<th class="left">PCA Model</th>\n')
            str_list.append('<th class="left">Heatmap Blind</th>\n')
            str_list.append('<th class="left">Heatmap Model</th>\n')
            str_list.append('<th class="left">Cook\'s Distances</th>\n')
            str_list.append('<th class="left">FPKM Density</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_name in sorted(design_dict):
                design_row_dict = design_dict[design_name]
                design_prefix = '_'.join((self.prefix, design_row_dict['design']))
                for plot_instance in design_row_dict['plot_aes'].split('|'):
                    plot_path = plot_instance.replace('geom_', '').replace(';', '__').translate(
                        str.maketrans(':=,', '___'))

                    str_list.append('<tr>\n')

                    # Design
                    str_list.append('<td>' + design_row_dict['design'] + '</td>\n')

                    # MDS and PCA plots for Blind and Model.
                    for plot_type in ('mds', 'pca'):
                        for model_type in ('blind', 'model'):
                            plot_path_pdf = '_'.join((plot_type, plot_path, model_type + '.pdf'))
                            plot_path_png = '_'.join((plot_type, plot_path, model_type + '.png'))
                            str_list.append('<td>')
                            if image_source_exists(prefix=design_prefix, suffix=plot_path_png):
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
                    plot_type = 'heatmap'
                    for model_type in ('blind', 'model'):
                        plot_path_pdf = '_'.join((plot_type, plot_path, model_type + '.pdf'))
                        plot_path_png = '_'.join((plot_type, plot_path, model_type + '.png'))
                        str_list.append('<td>')
                        if image_source_exists(prefix=design_prefix, suffix=plot_path_png):
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

                    # Cook's Distances
                    str_list.append('<td></td>\n')

                    # FPKM Density
                    str_list.append('<td></td>\n')

                    str_list.append('</tr>\n')

                # Add a line with the variance-per-principal-component and RIN score density plots.
                str_list.append('<tr>\n')

                # Design
                str_list.append('<td>' + design_row_dict['design'] + '</td>\n')

                # MDS Blind
                str_list.append('<td></td>\n')

                # MDS Model (use for RIN density plot if available)
                str_list.append('<td>')
                if image_source_exists(prefix=design_prefix, suffix='rin_density.png'):
                    str_list.append(self.get_html_anchor(
                        prefix=design_prefix,
                        suffix='rin_density.pdf',
                        text=self.get_html_image(
                            prefix=design_prefix,
                            suffix='rin_density.png',
                            text='RIN density plot',
                            height='80',
                            width='80')))
                str_list.append('</td>\n')

                # PCA Variance per Principal Component for Blind and Model.
                plot_type = 'pca'
                plot_path = 'variance'
                for model_type in ('blind', 'model'):
                    plot_path_pdf = '_'.join((plot_type, plot_path, model_type + '.pdf'))
                    plot_path_png = '_'.join((plot_type, plot_path, model_type + '.png'))
                    str_list.append('<td>')
                    if image_source_exists(prefix=design_prefix, suffix=plot_path_png):
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

                # Heatmap Blind
                str_list.append('<td></td>\n')

                # Heatmap Model
                str_list.append('<td></td>\n')

                # Cook's Distance Box Plot
                str_list.append('<td>\n')
                if image_source_exists(prefix=design_prefix, suffix='cooks_distances.png'):
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
                if image_source_exists(prefix=design_prefix, suffix='fpkm_density.png'):
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

            # Likelihood Ratio Testing (LRT) Table

            str_list.append('<h2 id="lrt">Likelihood Ratio Testing (LRT)</h2>\n')
            str_list.append('\n')

            str_list.append('<p>')
            str_list.append('A ')
            str_list.append('<a href="https://en.wikipedia.org/wiki/Likelihood-ratio_test">Likelihood-ratio test</a> ')
            str_list.append('compares the goodness of fit of two models, a null (i.e., full) model and an alternative ')
            str_list.append('(i.e., reduced) model. The likelihood ratio expresses how many times more likely ')
            str_list.append('the data fits one model than the other. ')
            str_list.append('Since each gene is modelled, the LRT allows identifying those genes that benefit ')
            str_list.append('from a variable or factor dropped in the reduced model. ')
            str_list.append('The intention behind LRT is to show that the terms included in the model are both, ')
            str_list.append('specific and relevant, as spurious terms unnecessarily reduce statistical power. ')
            str_list.append('For the differential expression calling in the ')
            str_list.append('<a href="#contrasts_table">Differential Expression Testing Table</a> below, ')
            str_list.append('the Wald test is applied.')
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
            # str_list.append('<th class="left">Effective Number</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_name in sorted(design_dict):
                design_row_dict = design_dict[design_name]
                design_prefix = '_'.join((self.prefix, design_row_dict['design']))
                lrt_summary_path = os.path.join(
                    self.genome_directory,
                    design_prefix,
                    design_prefix + '_lrt_summary.tsv')

                if os.path.exists(lrt_summary_path):
                    lrt_summary_sheet: AnnotationSheet = AnnotationSheet.from_file_path(
                        file_path=lrt_summary_path,
                        file_type='excel-tab',
                        name='DESeq LRT Summary Table')

                    lrt_row_dict_list = lrt_summary_sheet.row_dict_list
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
                                        suffix='_'.join(('lrt', lrt_row_dict['reduced_name'] + 'differential.tsv')),
                                        text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                                    '</td>\n')

                    # Significant genes
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('lrt', lrt_row_dict['reduced_name'], 'significant.tsv')),
                                        text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                                    '</td>\n')

                    # Significant Number
                    if lrt_row_dict['significant']:
                        str_list.append('<td class="right">{:,}</td>\n'.format(int(lrt_row_dict['significant'])))
                    else:
                        str_list.append('<td class="right"></td>\n')

                    # Effective Number
                    # if lrt_row_dict['effective']:
                    #     str_list.append('<td class="right">{:,}</td>\n'.format(int(lrt_row_dict['effective'])))
                    # else:
                    #     str_list.append('<td class="right"></td>\n')

                    str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Contrast table.

            str_list.append('<h2 id="contrasts">Differential Expression Testing</h2>\n')
            str_list.append('<p>')
            str_list.append('Differential expression calls are based on the full model formula ')
            str_list.append('shown in the <a href="#lrt_table">LRT Table</a> above, ')
            str_list.append('using the negative binomial distribution. ')
            str_list.append('Biologically meaningful contrasts (i.e., comparisons) were extracted from the model, ')
            str_list.append('log2-fold values where shrunk with the CRAN ')
            str_list.append('<a href="https://cran.r-project.org/package=ashr">ashr</a> package, ')
            str_list.append('while two-tailed p-values obtained from Wald testing were adjusted with the ')
            str_list.append('<a href="https://bioconductor.org/">Bioconductor</a> Independent Hypothesis Weighting ')
            str_list.append('(<a href="https://www.bioconductor.org/packages/IHW/">IHW</a>) package.')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="contrasts_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Contrast</th>\n')
            str_list.append('<th class="left">Differential Genes</th>\n')
            str_list.append('<th class="left">Significant Genes</th>\n')
            str_list.append('<th class="left">Significant Number</th>\n')
            str_list.append('<th class="left">Significant Up</th>\n')
            str_list.append('<th class="left">Significant Down</th>\n')
            # str_list.append('<th class="left">Effective Genes</th>\n')
            # str_list.append('<th class="left">Effective Number</th>\n')
            # str_list.append('<th class="left">Effective Up</th>\n')
            # str_list.append('<th class="left">Effective Down</th>\n')
            str_list.append('<th class="left">MA Plot</th>\n')
            str_list.append('<th class="left">Volcano Plot</th>\n')
            str_list.append('<th class="left">Numerator</th>\n')
            str_list.append('<th class="left">Denominator</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            # Read the contrasts from a design-specific summary file including numbers of significantly differentially
            # expressed genes or the project-specific contrasts file as a backup.

            contrast_sheet: AnnotationSheet = AnnotationSheet.from_file_path(
                file_path=self.contrast_path,
                file_type='excel',
                name='DESeq Contrasts Table')

            # For compatibility with design-specific summary files, re-index the contrast sheet on design names.
            contrast_dict: dict[str, list[dict[str, str]]] = dict()

            for contrast_row_dict in contrast_sheet.row_dict_list:
                if contrast_row_dict['Design'] not in contrast_dict:
                    contrast_dict[contrast_row_dict['Design']] = list()
                contrast_dict[contrast_row_dict['Design']].append(contrast_row_dict)

            # Print contrasts ordered by design names.
            for design_name in sorted(contrast_dict):
                # Exclude designs that are either not or no longer in the design dict since they were excluded.
                if design_name not in design_dict:
                    continue

                design_prefix = '_'.join((self.prefix, design_name))
                contrasts_summary_path = os.path.join(
                    self.genome_directory,
                    design_prefix,
                    design_prefix + '_contrasts_summary.tsv')
                if os.path.exists(contrasts_summary_path):
                    summary_sheet: AnnotationSheet = AnnotationSheet.from_file_path(
                        file_path=contrasts_summary_path,
                        file_type='excel-tab',
                        name='DESeq Contrasts Summary Table')
                    summary_row_dict_list = summary_sheet.row_dict_list
                else:
                    summary_row_dict_list = contrast_dict[design_name]

                for row_dict in summary_row_dict_list:
                    # Exclude contrasts if requested.
                    if contrast_sheet.get_boolean(row_dict=row_dict, key='Exclude'):
                        continue

                    str_list.append('<tr>\n')

                    # Design
                    str_list.append('<td class="left">' + row_dict['Design'] + '</td>\n')

                    # Contrast
                    str_list.append('<td class="left">' + row_dict['Label'] + '</td>\n')

                    # TSV
                    numerator = row_dict['Numerator'].replace(',', '_')
                    denominator = row_dict['Denominator'].replace(',', '_')
                    if not denominator or denominator == 'NA':
                        denominator = 'intercept'

                    # Differential Genes
                    str_list.append('<td>' +
                                    self.get_html_anchor(
                                        prefix=design_prefix,
                                        suffix='_'.join(('contrast', numerator, 'against', denominator,
                                                         'differential.tsv')),
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

                    # Significant Up genes
                    if 'SignificantUp' in row_dict:
                        str_list.append('<td class="right">{:,}</td>\n'.format(int(row_dict['SignificantUp'])))
                    else:
                        str_list.append('<td></td>\n')

                    # Significant Down genes
                    if 'SignificantDown' in row_dict:
                        str_list.append('<td class="right">{:,}</td>\n'.format(int(row_dict['SignificantDown'])))
                    else:
                        str_list.append('<td></td>\n')

                    # Effective Genes
                    # str_list.append('<td>' +
                    #                 self.get_html_anchor(
                    #                     prefix=design_prefix,
                    #                     suffix='_'.join(('contrast', numerator, 'against', denominator,
                    #                                      'effective.tsv')),
                    #                     text='<abbr title="Tab-Separated Value">TSV</abbr>') +
                    #                 '</td>\n')

                    # Effective Number
                    # if 'Effective' in row_dict:
                    #     str_list.append('<td class="right">{:,}</td>\n'.format(int(row_dict['Effective'])))
                    # else:
                    #     str_list.append('<td></td>\n')

                    # Effective Up genes
                    # if 'EffectiveUp' in row_dict:
                    #     str_list.append('<td class="right">{:,}</td>\n'.format(int(row_dict['EffectiveUp'])))
                    # else:
                    #     str_list.append('<td></td>\n')

                    # Effective Down genes
                    # if 'EffectiveDown' in row_dict:
                    #     str_list.append('<td class="right">{:,}</td>\n'.format(int(row_dict['EffectiveDown'])))
                    # else:
                    #     str_list.append('<td></td>\n')

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
                    str_list.append('<td class="left">')
                    if not row_dict['Denominator'] or row_dict['Denominator'] == 'NA':
                        str_list.append('intercept')
                    else:
                        str_list.append(row_dict['Denominator'])
                    str_list.append('</td>\n')

                    str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Link Enrichr, Heatmap and Volcano reports.

            str_list.append('<h2 id="accessory_reports">Accessory Reports</h2>\n')
            str_list.append('\n')

            str_list.append('<p>Accessory reports provide annotated expression heat map plots and\n')
            str_list.append('<a href="https://amp.pharm.mssm.edu/Enrichr/">Enrichr</a> results of\n')
            str_list.append('up- and down-regulated genes for selected libraries.</p>\n')

            str_list.append('<table id="accessory_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Enrichr Report</th>\n')
            str_list.append('<th class="left">Heatmap Report</th>\n')
            str_list.append('<th class="left">Volcano Report</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_name in sorted(design_dict):
                str_list.append('<tr>\n')

                # Design
                str_list.append('<td>' + design_name + '</td>\n')

                # Enrichr Report
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

                # Heatmap Report
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

                # Volcano Report
                volcano_prefix = '_'.join((self.prefix, design_name, 'volcano'))
                if os.path.exists(os.path.join(
                        self.genome_directory,
                        volcano_prefix,
                        volcano_prefix + '_report.html')):
                    str_list.append('<td><a href="' +
                                    '/'.join((volcano_prefix, volcano_prefix + '_report.html')) +
                                    '">HTML</a></td>\n')
                else:
                    str_list.append('<td></td>\n')

                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Additional tables

            str_list.append('<h2 id="tables">Additional Tables</h2>\n')
            str_list.append('\n')
            str_list.append('<p>')
            str_list.append('For each design, the following tables are provided in tab-separated value (TSV) format.')
            str_list.append('</p>\n')

            str_list.append('<ul>\n')
            str_list.append('<li>')
            str_list.append('<strong>Raw Counts:</strong> ')
            str_list.append('Raw (integer) counts per gene and sample.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>Normalised Counts:</strong> ')
            str_list.append('Normalised (real) counts per gene and sample ')
            str_list.append('scaled to the sequencing depth by DESeq2.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>VST Blind Counts:</strong> ')
            str_list.append('Counts on a log2 scale per gene and sample after ')
            str_list.append('variance stabilising transformation (VST) disregarding the ')
            str_list.append('generalised linear model design.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>VST Model Counts:</strong> ')
            str_list.append('Counts on a log2 scale per gene and sample after ')
            str_list.append('variance stabilising transformation (VST) considering the ')
            str_list.append('generalised linear model design.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>FPKMs:</strong> ')
            str_list.append('Fragments Per Kilobase Per Million Reads Mapped per gene and sample, ')
            str_list.append('calculated by DESeq2 based on the union of all exons of a gene and ')
            str_list.append('thus an underestimate.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>Samples:</strong> ')
            str_list.append('The sample annotation table for a particular design.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>Features:</strong> ')
            str_list.append('The feature (i.e., gene) annotation imported from the reference transcriptome.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>Contrasts:</strong> ')
            str_list.append('The contrast summary table.')
            str_list.append('</li>\n')
            str_list.append('<li>')
            str_list.append('<strong>LRT:</strong> ')
            str_list.append('The Likelihood Ratio Testing summary table.')
            str_list.append('</li>\n')
            str_list.append('</ul>\n')

            str_list.append('<table id="table_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th class="left">Design</th>\n')
            str_list.append('<th class="left">Raw Counts</th>\n')
            str_list.append('<th class="left">Normalised Counts</th>\n')
            str_list.append('<th class="left">VST Blind Counts</th>\n')
            str_list.append('<th class="left">VST Model Counts</th>\n')
            str_list.append('<th class="left">FPKMs</th>\n')
            str_list.append('<th class="left">Samples</th>\n')
            str_list.append('<th class="left">Features</th>\n')
            str_list.append('<th class="left">Contrasts</th>\n')
            str_list.append('<th class="left">LRT</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for design_name in sorted(design_dict):
                design_row_dict = design_dict[design_name]
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
                    suffix='counts_normalised.tsv',
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
                    suffix='features_gene.tsv',
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

            return

        report_html()

        return
