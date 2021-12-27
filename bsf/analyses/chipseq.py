# -*- coding: utf-8 -*-
"""ChIP-seq Analysis module.

A package of classes and methods supporting ChIP-Seq analyses.
"""
#  Copyright 2013 - 2021 Michael K. Schuster
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
from typing import Callable, Dict, List

from bsf.analyses.bowtie import Bowtie2
from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.connector import ConnectorFile
from bsf.ngs import Collection, Sample
from bsf.procedure import FilePath, ConsecutiveRunnable
from bsf.process import Command, RunnableStep, RunnableStepMakeDirectory, RunnableStepLink
from bsf.standards import Configuration, StandardFilePath, Genome, Transcriptome


class ChIPSeqComparison(object):
    """ChIP-Seq comparison annotation sheet.

    @ivar name: Comparison name
    @type name: str
    @ivar group: Group
    @type group: str
    @ivar c_name: Control name
    @type c_name: str | None
    @ivar t_name: Treatment name
    @type t_name: str | None
    @ivar c_samples: Python C{list} of control C{bsf.ngs.Sample} objects
    @type c_samples: list[Sample]
    @ivar t_samples: Python C{list} of treatment C{bsf.ngs.Sample} objects
    @type t_samples: list[Sample]
    @ivar factor: ChIP factor
    @type factor: str
    @ivar tissue: Tissue
    @type tissue: str
    @ivar condition: Condition
    @type condition: str
    @ivar treatment: Treatment
    @type treatment: str
    @ivar replicate: replicate number
    @type replicate: int
    @ivar diff_bind: Run the DiffBind analysis
    @type diff_bind: bool
    """

    def __init__(
            self,
            name,
            group,
            c_name,
            t_name,
            c_samples,
            t_samples,
            factor,
            tissue=None,
            condition=None,
            treatment=None,
            replicate=None,
            diff_bind=None):
        """Initialise a C{bsf.analyses.chipseq.ChIPSeqComparison}.

        @param name: Comparison name
        @type name: str
        @param group: Group
        @type group: str
        @param c_name: Control name
        @type c_name: str | None
        @param t_name: Treatment name
        @type t_name: str | None
        @param c_samples: Python C{list} of control C{bsf.ngs.Sample} objects
        @type c_samples: list[Sample] | None
        @param t_samples: Python C{list} of treatment C{bsf.ngs.Sample} objects
        @type t_samples: list[Sample] | None
        @param factor: ChIP factor
        @type factor: str | None
        @param tissue: Tissue
        @type tissue: str | None
        @param condition: Condition
        @type condition: str | None
        @param treatment: Treatment
        @type treatment: str | None
        @param replicate: replicate number
        @type replicate: int | None
        @param diff_bind: Run the DiffBind analysis
        @type diff_bind: bool | None
        """
        super(ChIPSeqComparison, self).__init__()

        # Condition', 'Treatment', 'Replicate',
        # 'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat'

        self.name = name
        self.group = group
        self.c_name = c_name
        self.t_name = t_name

        if c_samples is None:
            self.c_samples = list()
        else:
            self.c_samples = c_samples

        if t_samples is None:
            self.t_samples = list()
        else:
            self.t_samples = t_samples

        if factor is None:
            self.factor = str()
        else:
            self.factor = factor

        if tissue is None:
            self.tissue = str()
        else:
            self.tissue = tissue

        if condition is None:
            self.condition = str()
        else:
            self.condition = condition

        if treatment is None:
            self.treatment = str()
        else:
            self.treatment = treatment

        if replicate is None:
            self.replicate = 0
        else:
            self.replicate = replicate

        if diff_bind is None:
            self.diff_bind = True
        else:
            self.diff_bind = diff_bind

        return

    def get_key(self):
        """Get a C{ChIPSeqComparison} key based on the control and treatment pair name, or just the treatment name.

        ChIP-Seq experiments use the order treatment versus control in comparisons.
        @return: C{ChIPSeqComparison} key
        @rtype: str
        """
        if self.c_name and self.t_name:
            return '__'.join((self.t_name, self.c_name))
        else:
            return self.t_name

    def get_prefix_peak_calling(self):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return ChIPSeq.get_prefix_peak_calling(t_name=self.t_name, c_name=self.c_name)

    def get_file_path_peak_calling(self):
        """Get a C{FilePathPeakCalling} object from this or a sub-class.

        @return: C{FilePathPeakCalling} or sub-class object
        @rtype: FilePathPeakCalling
        """
        return ChIPSeq.get_file_path_peak_calling(t_name=self.t_name, c_name=self.c_name)


class FilePathAlignment(FilePath):
    """The C{bsf.analyses.chipseq.FilePathAlignment} models alignment file paths.

    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar sample_bam: Sample-specific BAM file path
    @type sample_bam: str
    @ivar sample_bai: Sample-specific BAI index file path
    @type sample_bai: str
    @ivar sample_md5: Sample-specific MD5 checksum file path
    @type sample_md5: str
    @ivar filter_metrics_tsv: deepTools alignmentSieve metrics
    @type filter_metrics_tsv: str
    @ivar coverage_bw: Coverage bigWig file path
    @type coverage_bw: str
    @ivar coverage_bwi_txt: Coverage bigWig information file path
    @type coverage_bwi_txt: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chipseq.FilePathAlignment} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathAlignment, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.sample_bam = prefix + '.bam'
        self.sample_bai = prefix + '.bam.bai'
        self.sample_md5 = prefix + '.bam.md5'

        self.filter_metrics_tsv = prefix + '_metrics.tsv'

        self.coverage_bw = prefix + '.bw'
        self.coverage_bwi_txt = prefix + '_bwi.txt'

        return


class FilePathPeakCalling(FilePath):
    """The C{bsf.analyses.chipseq.FilePathPeakCalling} models files in a comparison-specific MACS directory.

    @ivar output_directory: Output directory
    @type output_directory: str
    @ivar name_prefix: MACS --name option
    @type name_prefix: str
    @ivar control_bdg: Control (lambda) signal bedGraph
    @type control_bdg: str
    @ivar control_bw: Control (lambda) signal bigWig
    @type control_bw: str
    @ivar control_bwi: Control (lambda) signal bigWigInfo
    @type control_bwi: str
    @ivar treatment_bdg: Treatment signal bedGraph
    @type treatment_bdg: str
    @ivar treatment_bw: Treatment signal bigWig
    @type treatment_bw: str
    @ivar treatment_bwi: Treatment signal bigWigInfo
    @type treatment_bwi: str
    @ivar comparison_log_fe_bdg: Comparison signal as log10(fold-enrichment) bedGraph
    @type comparison_log_fe_bdg: str
    @ivar comparison_log_fe_bw: Comparison signal as log10(fold-enrichment) bigWig
    @type comparison_log_fe_bw: str
    @ivar comparison_log_fe_bwi: Comparison signal as log10(fold-enrichment) bigWigInfo
    @type comparison_log_fe_bwi: str
    @ivar comparison_ppois_bdg: Comparison signal as Poisson p-value bedGraph
    @type comparison_ppois_bdg: str
    @ivar comparison_ppois_bw: Comparison signal as Poisson p-value bigWig
    @type comparison_ppois_bw: str
    @ivar comparison_ppois_bwi: Comparison signal as Poisson p-value bigWigInfo
    @type comparison_ppois_bwi: str
    @ivar comparison_subtract_bdg: Comparison signal as subtraction bedGraph
    @type comparison_subtract_bdg: str
    @ivar comparison_subtract_bw: Comparison signal as subtraction bigWig
    @type comparison_subtract_bw: str
    @ivar comparison_subtract_bwi: Comparison signal as subtraction bigWigInfo
    @type comparison_subtract_bwi: str
    @ivar summits_bed: Peak summits BED
    @type summits_bed: str
    @ivar summits_bb: Peak summits bigBed
    @type summits_bb: str
    @ivar summits_bbi: Peak summits bigBedInfo
    @type summits_bbi: str
    @ivar broad_peaks_bed: Broad peaks BED
    @type broad_peaks_bed: str
    @ivar broad_peaks_bb: Broad peaks bigBed
    @type broad_peaks_bb: str
    @ivar broad_peaks_bbi: Broad peaks bigBedInfo
    @type broad_peaks_bbi: str
    @ivar gapped_peaks_bed: Gapped peaks BED
    @type gapped_peaks_bed: str
    @ivar gapped_peaks_bb: Gapped peaks bigBed
    @type gapped_peaks_bb: str
    @ivar gapped_peaks_bbi: Gapped peaks bigBedInfo
    @type gapped_peaks_bbi: str
    @ivar narrow_peaks_bed: Narrow peaks BED
    @type narrow_peaks_bed: str
    @ivar narrow_peaks_bb: Narrow peaks bigBed
    @type narrow_peaks_bb: str
    @ivar narrow_peaks_bbi: Narrow peaks bigBedInfo
    @type narrow_peaks_bbi: str
    @ivar peaks_broad: MACS peaks broadPeak
    @type peaks_broad: str
    @ivar peaks_gapped: MACS peaks gappedPeak
    @type peaks_gapped: str
    @ivar peaks_narrow: MACS peaks narrowPeak
    @type peaks_narrow: str
    @ivar peaks_tsv: MACS2 peaks tab-separated value (TSV)
    @type peaks_tsv: str
    @ivar peaks_xls: MACS2 peaks tab-separated value (TSV)
    @type peaks_xls: str
    @ivar model_r: MACS peak model as Rscript
    @type model_r: str
    @ivar model_pdf: MACS "Peak Model" and "Cross-Correlation" plots PDF
    @type model_pdf: str
    @ivar model_0_png: MACS "Peak Model" plot PNG
    @type model_0_png: str
    @ivar model_1_png: MACS "Cross-Correlation" plot PNG
    @type model_1_png: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chipseq.FilePathPeakCalling} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathPeakCalling, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.name_prefix = os.path.join(prefix, prefix)  # MACS2 --name option

        self.control_bdg = os.path.join(prefix, '_'.join((prefix, 'control_lambda.bdg')))
        self.control_bw = os.path.join(prefix, '_'.join((prefix, 'control_lambda.bw')))
        self.control_bwi = os.path.join(prefix, '_'.join((prefix, 'control_lambda_bwi.txt')))

        self.treatment_bdg = os.path.join(prefix, '_'.join((prefix, 'treat_pileup.bdg')))
        self.treatment_bw = os.path.join(prefix, '_'.join((prefix, 'treat_pileup.bw')))
        self.treatment_bwi = os.path.join(prefix, '_'.join((prefix, 'treat_pileup_bwi.txt')))

        self.comparison_log_fe_bdg = os.path.join(prefix, '_'.join((prefix, 'logFE.bdg')))
        self.comparison_log_fe_bw = os.path.join(prefix, '_'.join((prefix, 'logFE.bw')))
        self.comparison_log_fe_bwi = os.path.join(prefix, '_'.join((prefix, 'logFE_bwi.txt')))

        self.comparison_ppois_bdg = os.path.join(prefix, '_'.join((prefix, 'ppois.bdg')))
        self.comparison_ppois_bw = os.path.join(prefix, '_'.join((prefix, 'ppois.bw')))
        self.comparison_ppois_bwi = os.path.join(prefix, '_'.join((prefix, 'ppois_bwi.txt')))

        self.comparison_subtract_bdg = os.path.join(prefix, '_'.join((prefix, 'subtract.bdg')))
        self.comparison_subtract_bw = os.path.join(prefix, '_'.join((prefix, 'subtract.bw')))
        self.comparison_subtract_bwi = os.path.join(prefix, '_'.join((prefix, 'subtract_bwi.txt')))

        self.peaks_broad = os.path.join(prefix, '_'.join((prefix, 'peaks.broadPeak')))

        self.broad_peaks_bed = os.path.join(prefix, '_'.join((prefix, 'broad_peaks.bed')))
        self.broad_peaks_bb = os.path.join(prefix, '_'.join((prefix, 'broad_peaks.bb')))
        self.broad_peaks_bbi = os.path.join(prefix, '_'.join((prefix, 'broad_peaks_bbi.txt')))

        self.peaks_gapped = os.path.join(prefix, '_'.join((prefix, 'peaks.gappedPeak')))

        self.gapped_peaks_bed = os.path.join(prefix, '_'.join((prefix, 'gapped_peaks.bed')))
        self.gapped_peaks_bb = os.path.join(prefix, '_'.join((prefix, 'gapped_peaks.bb')))
        self.gapped_peaks_bbi = os.path.join(prefix, '_'.join((prefix, 'gapped_peaks_bbi.txt')))

        self.peaks_narrow = os.path.join(prefix, '_'.join((prefix, 'peaks.narrowPeak')))

        self.narrow_peaks_bed = os.path.join(prefix, '_'.join((prefix, 'narrow_peaks.bed')))
        self.narrow_peaks_bb = os.path.join(prefix, '_'.join((prefix, 'narrow_peaks.bb')))
        self.narrow_peaks_bbi = os.path.join(prefix, '_'.join((prefix, 'narrow_peaks_bbi.txt')))

        self.summits_bed = os.path.join(prefix, '_'.join((prefix, 'summits.bed')))
        self.summits_bb = os.path.join(prefix, '_'.join((prefix, 'summits.bb')))
        self.summits_bbi = os.path.join(prefix, '_'.join((prefix, 'summits_bbi.txt')))

        self.peaks_tsv = os.path.join(prefix, '_'.join((prefix, 'peaks.tsv')))
        self.peaks_xls = os.path.join(prefix, '_'.join((prefix, 'peaks.xls')))

        self.model_r = os.path.join(prefix, '_'.join((prefix, 'model.r')))
        self.model_pdf = os.path.join(prefix, '_'.join((prefix, 'model.pdf')))
        self.model_0_png = os.path.join(prefix, '_'.join((prefix, 'model-0.png')))
        self.model_1_png = os.path.join(prefix, '_'.join((prefix, 'model-1.png')))

        return


class FilePathChIPQC(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chipseq.FilePathChIPQC} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathChIPQC, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.report_html = os.path.join(prefix, 'ChIPQC.html')

        return


class FilePathDiffBind(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chipseq.FilePathDiffBind} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathDiffBind, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sample_annotation_sheet = prefix + '_samples.csv'
        self.correlation_read_counts_pdf = os.path.join(prefix, prefix + '_correlation_read_counts.pdf')
        self.correlation_read_counts_png = os.path.join(prefix, prefix + '_correlation_read_counts.png')
        self.correlation_peak_caller_score_pdf = os.path.join(prefix, prefix + '_correlation_peak_caller_score.pdf')
        self.correlation_peak_caller_score_png = os.path.join(prefix, prefix + '_correlation_peak_caller_score.png')
        self.correlation_analysis_pdf = os.path.join(prefix, prefix + '_correlation_analysis.pdf')
        self.correlation_analysis_png = os.path.join(prefix, prefix + '_correlation_analysis.png')
        self.pca_pdf = os.path.join(prefix, prefix + '_pca_plot.pdf')
        self.pca_png = os.path.join(prefix, prefix + '_pca_plot.png')
        self.contrasts_csv = os.path.join(prefix, prefix + '_contrasts.csv')
        self.genes_complete_tsv = os.path.join(prefix, prefix + '_peak_set_genes_complete.tsv')
        self.regions_pdf = os.path.join(prefix, prefix + '_peak_set_regions.pdf')
        self.regions_png = os.path.join(prefix, prefix + '_peak_set_regions.png')
        self.regions_tsv = os.path.join(prefix, prefix + '_peak_set_regions.tsv')
        self.peak_set_bed = os.path.join(prefix, prefix + '_peak_set.bed')
        self.peak_set_bb = os.path.join(prefix, prefix + '_peak_set.bb')

        return


class FilePathDiffBindContrast(FilePath):
    def __init__(self, prefix, group_1, group_2):
        """Initialise a C{bsf.analyses.chipseq.FilePathDiffBind} object

        @param prefix: Prefix
        @type prefix: str
        @param group_1: Group 1
        @type group_1: str
        @param group_2: Group 2
        @type group_2: str
        """
        super(FilePathDiffBindContrast, self).__init__(prefix=prefix)

        suffix = group_1 + '__' + group_2

        self.ma_plot_pdf = os.path.join(prefix, prefix + '_ma_plot_' + suffix + '.pdf')
        self.ma_plot_png = os.path.join(prefix, prefix + '_ma_plot_' + suffix + '.png')
        self.scatter_plot_pdf = os.path.join(prefix, prefix + '_scatter_plot_' + suffix + '.pdf')
        self.scatter_plot_png = os.path.join(prefix, prefix + '_scatter_plot_' + suffix + '.png')
        self.pca_plot_pdf = os.path.join(prefix, prefix + '_pca_plot_' + suffix + '.pdf')
        self.pca_plot_png = os.path.join(prefix, prefix + '_pca_plot_' + suffix + '.png')
        self.box_plot_pdf = os.path.join(prefix, prefix + '_box_plot_' + suffix + '.pdf')
        self.box_plot_png = os.path.join(prefix, prefix + '_box_plot_' + suffix + '.png')
        self.peaks_tsv = os.path.join(prefix, prefix + '_peaks_' + suffix + '.tsv')
        self.genes_complete_tsv = os.path.join(prefix, prefix + '_peaks_' + suffix + '_genes_complete.tsv')
        self.genes_significant_tsv = os.path.join(prefix, prefix + '_peaks_' + suffix + '_genes_significant.tsv')
        self.regions_pdf = os.path.join(prefix, prefix + '_peaks_' + suffix + '_regions.pdf')
        self.regions_png = os.path.join(prefix, prefix + '_peaks_' + suffix + '_regions.png')
        self.regions_tsv = os.path.join(prefix, prefix + '_peaks_' + suffix + '_regions.tsv')

        return


class ChIPSeqDiffBindSheet(AnnotationSheet):
    """ChIP-Seq Bioconductor DiffBind annotation sheet class.
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'SampleID',
        'Tissue',
        'Factor',
        'Condition',
        'Treatment',
        'Replicate',
        'bamReads',
        'bamControl',
        'ControlID',
        'Peaks',
        'PeakCaller',
        'PeakFormat',
        # 'ScoreCol',
        # 'LowerBetter',
        # 'Counts',
    ]

    _test_methods: Dict[str, List[Callable[[int, Dict[str, str], str], str]]] = {
        'SampleID': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Tissue': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Factor': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Condition': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Treatment': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Replicate': [
            AnnotationSheet.check_numeric,
        ],
        'ControlID': [
            AnnotationSheet.check_alphanumeric,
        ],
        'PeakCaller': [
            AnnotationSheet.check_alphanumeric,
        ],
        'PeakFormat': [
            AnnotationSheet.check_alphanumeric,
        ],
        'ScoreCol': [
            AnnotationSheet.check_alphanumeric,
        ],
        'LowerBetter': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Counts': [
            AnnotationSheet.check_alphanumeric
        ],
    }

    def sort(self):
        """Sort by columns I{Tissue}, I{Factor}, I{Condition}, I{Treatment} and I{Replicate}.
        """
        self.row_dicts.sort(
            key=lambda item: '_'.join((
                item['Tissue'],
                item['Factor'],
                item['Condition'],
                item['Treatment'],
                '{:06d}'.format(int(item['Replicate'])))))

        return

    def to_file_path(self, adjust_field_names=None):
        """Write a C{bsf.analyses.ChIPSeqDiffBindSheet} to a file.

        @param adjust_field_names: Clear and adjust the Python C{list} of Python C{str} field name objects
        @type adjust_field_names: bool
        """
        # Override the method from the super-class to automatically sort before writing to a file.

        self.sort()
        super(ChIPSeqDiffBindSheet, self).to_file_path(adjust_field_names=adjust_field_names)

        return


class ChIPSeq(Analysis):
    """The C{bsf.analyses.chipseq.ChIPSeq} class represents the logic to run a ChIP-Seq-specific Analysis.

    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar skip_alignment_sieve: Skip the DeepTools alignment sieve
    @type skip_alignment_sieve: bool
    @ivar replicate_grouping: Group all replicates into a single process
    @type replicate_grouping: bool | None
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | None
    @ivar genome_black_list: Genome black list file path
    @type genome_black_list: str | None
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | None
    @ivar genome_sizes_path: Reference genome (chromosome) sizes file path
    @type genome_sizes_path: str | None
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str | None
    @ivar transcriptome_gtf_path: Transcriptome GTF file path
    @type transcriptome_gtf_path: str | None
    @ivar transcriptome_txdb_path: Transcriptome TxDb file path
    @type transcriptome_txdb_path: str | None
    @ivar colour_default: Default UCSC Genome Browser Track Hub RGB colour
    @type colour_default: str | None
    @ivar colour_dict: Python C{dict} of
        Python C{str} factor name key and
        Python C{str} RGB colour value data
    @type colour_dict: dict[str, str] | None
    @ivar factor_default: Default factor
    @type factor_default: str
    """

    name = 'ChIP-seq Analysis'
    prefix = 'chipseq'
    skip_alignment_sieve = True

    @classmethod
    def get_stage_name_alignment(cls):
        """Get a Python C{str} for a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'alignment'))

    @classmethod
    def get_stage_name_peak_calling(cls):
        """Get a Python C{str} for a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'peak_calling'))

    @classmethod
    def get_stage_name_chipqc(cls):
        """Get a Python C{str} for a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'chipqc'))

    @classmethod
    def get_stage_name_diff_bind(cls):
        """Get a Python C{str} for a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'diff_bind'))

    @classmethod
    def get_prefix_alignment(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_alignment(), sample_name))

    @classmethod
    def get_prefix_peak_calling(cls, t_name, c_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param t_name: Treatment C{bsf.ngs.Sample.name}
        @type t_name: str
        @param c_name: Control C{bsf.ngs.Sample.name}
        @type c_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        if c_name:
            return cls.get_stage_name_peak_calling() + '_' + t_name + '__' + c_name
        else:
            return cls.get_stage_name_peak_calling() + '_' + t_name

    @classmethod
    def get_prefix_chipqc(cls, comparison_name, factor_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @param factor_name: Factor name
        @type factor_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_chipqc(), comparison_name, factor_name))

    @classmethod
    def get_prefix_diff_bind(cls, comparison_name, factor_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @param factor_name: Factor name
        @type factor_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_diff_bind(), comparison_name, factor_name))

    @classmethod
    def get_file_path_alignment(cls, sample_name):
        """Get a C{FilePathAlignment} object.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: C{FilePathAlignment}
        @rtype: FilePathAlignment
        """
        return FilePathAlignment(prefix=cls.get_prefix_alignment(sample_name=sample_name))

    @classmethod
    def get_file_path_peak_calling(cls, t_name, c_name):
        """Get a C{FilePathPeakCalling} object from this or a sub-class.

        @param t_name: C{bsf.ngs.Sample.name}
        @type t_name: str
        @param c_name: C{bsf.ngs.Sample.name}
        @type c_name: str
        @return: C{FilePathPeakCalling} or sub-class object
        @rtype: FilePathPeakCalling
        """
        return FilePathPeakCalling(
            prefix=cls.get_prefix_peak_calling(
                t_name=t_name,
                c_name=c_name))

    @classmethod
    def get_file_path_chipqc(cls, comparison_name, factor_name):
        """Get a C{FilePathChIPQC} object from this or a sub-class.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @param factor_name: Factor name
        @type factor_name: str
        @return: C{FilePathChIPQC} or sub-class object
        @rtype: FilePathChIPQC
        """
        return FilePathChIPQC(
            prefix=cls.get_prefix_chipqc(
                comparison_name=comparison_name,
                factor_name=factor_name))

    @classmethod
    def get_file_path_diff_bind(cls, comparison_name, factor_name):
        """Get a C{FilePathDiffBind} object from this or a sub-class.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @param factor_name: Factor name
        @type factor_name: str
        @return: C{FilePathDiffBind} or sub-class object
        @rtype: FilePathDiffBind
        """
        return FilePathDiffBind(
            prefix=cls.get_prefix_diff_bind(
                comparison_name=comparison_name,
                factor_name=factor_name))

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
            genome_black_list=None,
            genome_fasta_path=None,
            genome_sizes_path=None,
            genome_effective_size=None,
            transcriptome_version=None,
            transcriptome_gtf_path=None,
            transcriptome_txdb_path=None,
            colour_default=None,
            colour_dict=None,
            factor_default=None):
        """Initialise a C{bsf.analyses.chipseq.ChIPSeq}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
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
        @type stage_list: list[Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[Sample]
        @param replicate_grouping: Group all replicates into a single process
        @type replicate_grouping: bool
        @param comparison_path: Comparison file path
        @type comparison_path: str | None
        @param genome_black_list: Genome black list file path
        @type genome_black_list: str | None
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | None
        @param genome_sizes_path: Reference genome (chromosome) sizes file path
        @type genome_sizes_path: str | None
        @param genome_effective_size: Effective genome size
        @type genome_effective_size: str | None
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @param transcriptome_gtf_path: Transcriptome GTF file path
        @type transcriptome_gtf_path: str | None
        @param transcriptome_txdb_path: Transcriptome TxDb file path
        @type transcriptome_txdb_path: str | None
        @param colour_default: Default UCSC Genome Browser Track Hub RGB colour
        @type colour_default: str | None
        @param colour_dict: Python C{dict} of
            Python C{str} factor name key and
            Python C{str} RGB colour value data
        @type colour_dict: dict[str, str] | None
        @param factor_default: Default factor
        @type factor_default: str
        """
        super(ChIPSeq, self).__init__(
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
            self.replicate_grouping = True
        else:
            self.replicate_grouping = replicate_grouping

        self.comparison_path = comparison_path
        self.colour_default = colour_default
        self.colour_dict = colour_dict
        self.factor_default = factor_default
        self.genome_black_list = genome_black_list
        self.genome_fasta_path = genome_fasta_path
        self.genome_sizes_path = genome_sizes_path
        self.genome_effective_size = genome_effective_size
        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf_path = transcriptome_gtf_path
        self.transcriptome_txdb_path = transcriptome_txdb_path

        self._comparison_dict: Dict[str, Dict[str, ChIPSeqComparison]] = dict()

        self._factor_dict: Dict[str, Dict[str, List[ChIPSeqComparison]]] = dict()

        self._macs_version: int = 2

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.chipseq.ChIPSeq} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """
        super(ChIPSeq, self).set_configuration(configuration=configuration, section=section)

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_black_list'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_black_list = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_sizes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_sizes_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_effective_size'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_effective_size = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_txdb'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_txdb_path = configuration.config_parser.get(section=section, option=option)

        option = 'colour_default'
        if configuration.config_parser.has_option(section=section, option=option):
            self.colour_default = configuration.config_parser.get(section=section, option=option)

        track_hub_section = '.'.join((section, 'TrackHubColours'))
        if configuration.config_parser.has_section(section=track_hub_section):
            if self.colour_dict is None:
                self.colour_dict = dict()

            for option in configuration.config_parser.options(section=track_hub_section):
                self.colour_dict[option] = configuration.config_parser.get(
                    section=track_hub_section,
                    option=option)

        option = 'factor_default'
        if configuration.config_parser.has_option(section=section, option=option):
            self.factor_default = configuration.config_parser.get(section=section, option=option)

        return

    def get_colour(self, factor):
        """Get a UCSC Genome Browser Track Hub RGB colour by a factor name.

        @param factor: ChIP-seq factor name
        @type factor: str
        @return: Track Hub RGB colour
        @rtype: str
        """
        if factor in self.colour_dict:
            return self.colour_dict[factor]
        else:
            return self.colour_default

    def run(self):
        """Run a C{bsf.analyses.chipseq.ChIPSeq} C{bsf.analysis.Analysis}.
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file from disk.

                - Column headers for CASAVA folders:
                    - Treatment/Control ProcessedRunFolder:
                        - CASAVA processed run folder name or
                        - C{bsf.analysis.Analysis.input_directory} by default
                    - Treatment/Control Project:
                        - CASAVA Project name or
                        - C{bsf.analysis.Analysis.project_name} by default
                    - Treatment/Control Sample:
                        - CASAVA Sample name, no default
                - Column headers for independent samples:
                    - Treatment/Control Group:
                    - Treatment/Control Sample:
                    - Treatment/Control Reads:
                    - Treatment/Control File:
            """
            annotation_sheet = AnnotationSheet.from_file_path(file_path=self.comparison_path)

            # Unfortunately, two passes through the comparison sheet are required.
            # In the first one merge all Sample objects that share the name.
            # Merging Sample objects is currently the only way to pool PairedReads objects,
            # which is required for ChIP-Seq experiments.

            # Second pass, add all Sample objects mentioned in a comparison.

            for row_dict in annotation_sheet.row_dicts:
                if 'Name' in row_dict and row_dict['Name']:
                    comparison_name = row_dict['Name']
                else:
                    comparison_name = 'global'

                if 'Group' in row_dict and row_dict['Group']:
                    comparison_group = row_dict['Group']
                else:
                    raise Exception('The ChIPSeq comparison sheet needs a group name for each treatment sample.')

                # If undefined, bsf.ngs.Collection.get_samples_from_row_dict() returns an empty str and list object.
                c_name, c_sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Control')
                t_name, t_sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Treatment')

                # For a successful comparison, at least a Python list of Sample objects has to be defined.

                if not len(t_sample_list):
                    if self.debug > 1:
                        print(
                            'Redundant comparison line with Treatment:', repr(t_name),
                            'samples:', len(t_sample_list),
                            'and Control:', repr(c_name),
                            'samples:', len(c_sample_list))
                    continue

                # Add all control Sample or SampleGroup objects to the Sample list.

                for c_sample in c_sample_list:
                    if self.debug > 1:
                        print('  Control Sample name:', c_sample.name, 'file_path:', c_sample.file_path)
                    if self.debug > 2:
                        sys.stdout.writelines(c_sample.trace(level=1))
                        # Find the Sample in the unified sample dictionary.
                    self.add_sample(sample=c_sample)

                # Add all treatment Sample or SampleGroup objects to the Sample list.

                for t_sample in t_sample_list:
                    if self.debug > 1:
                        print('  Treatment Sample name:', t_sample.name, 'file_path:', t_sample.file_path)
                    if self.debug > 2:
                        sys.stdout.writelines(t_sample.trace(level=1))
                    self.add_sample(sample=t_sample)

                if 'Tissue' in row_dict:
                    tissue = row_dict['Tissue']
                else:
                    tissue = str()

                if 'Factor' in row_dict:
                    factor = row_dict['Factor']
                else:
                    factor = self.factor_default

                if 'Condition' in row_dict:
                    condition = row_dict['Condition']
                else:
                    condition = str()

                if 'Treatment' in row_dict:
                    treatment = row_dict['Treatment']
                else:
                    treatment = str()

                if 'DiffBind' in row_dict:
                    diff_bind = annotation_sheet.get_boolean(row_dict=row_dict, key='DiffBind')
                else:
                    diff_bind = True

                # Automatically create replicate numbers for the sample annotation sheet required by the
                # DiffBind Bioconductor package.
                # Use a first-level dict with replicate key data and second-level dict value data.
                # The second-level dict stores Treatment Sample key data and int value data.

                comparison = ChIPSeqComparison(
                    name=comparison_name,
                    group=comparison_group,
                    c_name=c_name,
                    t_name=t_name,
                    c_samples=c_sample_list,
                    t_samples=t_sample_list,
                    factor=factor,
                    tissue=tissue,
                    condition=condition,
                    treatment=treatment,
                    replicate=0,
                    diff_bind=diff_bind)

                comparison_pair = comparison.get_key()

                if comparison_name not in self._comparison_dict:
                    self._comparison_dict[comparison_name] = dict()

                self._comparison_dict[comparison_name][comparison_pair] = comparison

            # Sort the comparison keys alphabetically and assign DiffBind replicate numbers
            # into ChIPSeqComparison objects.

            for comparison_dict in self._comparison_dict.values():
                # Rearrange by comparison group and comparison key.
                comparison_group_dict: Dict[str, List[ChIPSeqComparison]] = dict()
                for comparison in comparison_dict.values():
                    if not comparison.diff_bind:
                        continue
                    if comparison.group not in comparison_group_dict:
                        comparison_group_dict[comparison.group] = list()
                    comparison_group_dict[comparison.group].append(comparison)

                # For each group, order by ChIPSeqComparison.get_key() and assign DiffBind replicate numbers.
                for comparison_list in comparison_group_dict.values():
                    for replicate, comparison in enumerate(sorted(comparison_list, key=lambda item: item.get_key())):
                        comparison.replicate = replicate + 1

            return

        def run_create_alignment_jobs():
            """Create alignment jobs.
            """
            for sample in self.sample_list:
                file_path_alignment = self.get_file_path_alignment(sample_name=sample.name)

                runnable_alignment = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_alignment(sample_name=sample.name),
                        working_directory=self.genome_directory,
                        debug=self.debug))
                self.set_stage_runnable(stage=stage_alignment, runnable=runnable_alignment)
                # No dependencies at this stage.

                if self.skip_alignment_sieve:
                    # Add RunnableStepLink objects for symbolic links for BAI and BAM files.

                    runnable_alignment.add_runnable_step(
                        runnable_step=RunnableStepLink(
                            name='link_bam',
                            source_path=Bowtie2.get_file_path_sample(sample_name=sample.name).sample_bam,
                            target_path=file_path_alignment.sample_bam))

                    runnable_alignment.add_runnable_step(
                        runnable_step=RunnableStepLink(
                            name='link_bai',
                            source_path=Bowtie2.get_file_path_sample(sample_name=sample.name).sample_bai,
                            target_path=file_path_alignment.sample_bai))
                else:
                    # Add a RunnableStep for the deepTools alignmentSieve tool.

                    runnable_step = RunnableStep(
                        name='alignment_sieve',
                        program='alignmentSieve')
                    runnable_alignment.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_option_long(
                        key='bam',
                        value=Bowtie2.get_file_path_sample(sample_name=sample.name).sample_bam)
                    runnable_step.add_option_long(key='outFile', value=file_path_alignment.sample_bam)
                    runnable_step.add_option_long(key='numberOfProcessors', value=str(stage_alignment.threads))
                    runnable_step.add_option_long(key='filterMetrics', value=file_path_alignment.filter_metrics_tsv)
                    runnable_step.add_option_long(key='blackListFileName', value=self.genome_black_list)

                    # Index the resulting filtered BAM file.

                    runnable_step = RunnableStep(
                        name='samtools_index',
                        program='samtools',
                        sub_command=Command(program='index'))
                    runnable_alignment.add_runnable_step(runnable_step=runnable_step)

                    sub_command = runnable_step.sub_command

                    sub_command.add_switch_short(key='b')
                    sub_command.add_option_short(key='@', value=str(stage_alignment.threads))

                    sub_command.arguments.append(file_path_alignment.sample_bam)

                # Add a RunnableStep for the deepTools bamCoverage tool.

                runnable_step = RunnableStep(
                    name='bam_coverage',
                    program='bamCoverage')
                runnable_alignment.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_option_long(key='bam', value=file_path_alignment.sample_bam)
                runnable_step.add_option_long(key='outFileName', value=file_path_alignment.coverage_bw)
                runnable_step.add_option_long(key='outFileFormat', value='bigwig')
                runnable_step.add_option_long(key='binSize', value='10')
                # deepTools checks for overlapping black list intervals.
                # The ENCODE hg38 black list (ENCFF419RSJ.bed) has such overlapping intervals,
                # but unfortunately, the Bioconductor GRanges::disjoin() method does not help.
                # runnable_step.add_option_long(key='blackListFileName', value=self.genome_black_list)
                runnable_step.add_option_long(key='numberOfProcessors', value=str(stage_alignment.threads))
                runnable_step.add_option_long(key='effectiveGenomeSize', value=self.genome_effective_size)
                # Normalise to reads per genomic content (1x normalization)
                runnable_step.add_option_long(key='normalizeUsing', value='RPGC')
                runnable_step.add_option_long(key='extendReads', value='175')

                # Capture the bigWig information for Track Hub generation.

                runnable_step = RunnableStep(
                    name='bigwiginfo',
                    program='bigWigInfo',
                    stdout=ConnectorFile(
                        file_path=file_path_alignment.coverage_bwi_txt,
                        file_mode='wt'))
                runnable_alignment.add_runnable_step(runnable_step=runnable_step)

                runnable_step.arguments.append(file_path_alignment.coverage_bw)

            return

        def run_create_macs1_jobs():
            """Create MACS1 peak caller jobs.
            """
            chipseq_comparison_key_list: List[str] = list()

            for comparison_name in sorted(self._comparison_dict):
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    # Since a particular peak calling pair may exists in more than one comparison group,
                    # make sure only unique pairs get submitted.
                    if chipseq_comparison.get_key() in chipseq_comparison_key_list:
                        continue
                    else:
                        chipseq_comparison_key_list.append(chipseq_comparison.get_key())

                    t_file_path_list = list()
                    for t_sample in chipseq_comparison.t_samples:
                        t_file_path_list.append(self.get_file_path_alignment(sample_name=t_sample.name).sample_bam)

                    c_file_path_list = list()
                    for c_sample in chipseq_comparison.c_samples:
                        c_file_path_list.append(self.get_file_path_alignment(sample_name=c_sample.name).sample_bam)

                    prefix_peak_calling = chipseq_comparison.get_prefix_peak_calling()

                    file_path_peak_calling = FilePathPeakCalling(prefix=prefix_peak_calling)

                    runnable_peak_calling = self.add_runnable_consecutive(
                        runnable=ConsecutiveRunnable(
                            name=prefix_peak_calling,
                            working_directory=self.genome_directory,
                            debug=self.debug))
                    executable_peak_calling = self.set_stage_runnable(
                        stage=stage_peak_calling,
                        runnable=runnable_peak_calling)

                    # Add a RunnableStep to create the output directory.

                    runnable_step = RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_peak_calling.output_directory)
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    # Add a RunnableStep for MACS14 call peak.

                    runnable_step = RunnableStep(
                        name='macs14_call_peak',
                        program='macs14',
                        sub_command=Command(program='callpeak'))
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    # Read RunnableStep options from configuration sections:
                    # [bsf.analyses.chipseq.ChIPSeq.macs14_call_peak]
                    # [bsf.analyses.chipseq.ChIPSeq.macs14_call_peak.callpeak]
                    # [bsf.analyses.chipseq.ChIPSeq.macs14_call_peak.callpeak.{factor}]
                    self.set_runnable_step_configuration(runnable_step=runnable_step, tag=chipseq_comparison.factor)

                    sub_command = runnable_step.sub_command

                    sub_command.add_option_multi_long(
                        key='treatment',
                        value=' '.join(t_file_path_list))

                    if c_file_path_list:
                        # Control (input) samples are optional.
                        sub_command.add_option_multi_long(
                            key='control',
                            value=' '.join(c_file_path_list))

                    # MACS14 can hopefully also cope with directories specified in the --name option, but
                    # the resulting R script has them set too. Hence the R script has to be started
                    # from the genome_directory. However, the R script needs re-writing anyway, because
                    # it would be better to use the PNG rather than the PDF device for plotting.
                    sub_command.sub_command.add_option_long(
                        key='name',
                        value=file_path_peak_calling.name_prefix)

                    # The 'gsize' option has to be specified via the configuration.ini file in section
                    # [bsf.analyses.chipseq.ChIPSeq.macs14_call_peak.callpeak].
                    sub_command.add_switch_long(key='single-profile')
                    sub_command.add_switch_long(key='call-subpeaks')
                    sub_command.add_switch_long(key='wig')

                    # Add a RunnableStep to process MACS14 output.

                    runnable_step = RunnableStep(
                        name='process_macs14',
                        program='bsf_chipseq_process_macs14.bash')
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    # Specify the output path as in the macs14 --name option.
                    runnable_step.arguments.append(prefix_peak_calling)
                    runnable_step.arguments.append(self.genome_sizes_path)

                    if os.path.exists(os.path.join(self.genome_directory, file_path_peak_calling.summits_bb)):
                        executable_peak_calling.submit = False

            return

        def run_create_macs2_jobs():
            """Create MACS2 peak caller jobs.
            """
            chipseq_comparison_key_list: List[str] = list()

            for comparison_name in sorted(self._comparison_dict):
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    # Since a particular peak calling pair may exists in more than one comparison group,
                    # make sure only unique pairs get submitted.
                    if chipseq_comparison.get_key() in chipseq_comparison_key_list:
                        continue
                    else:
                        chipseq_comparison_key_list.append(chipseq_comparison.get_key())

                    prefix_peak_calling = chipseq_comparison.get_prefix_peak_calling()

                    file_path_peak_calling = chipseq_comparison.get_file_path_peak_calling()

                    runnable_peak_calling = self.add_runnable_consecutive(
                        runnable=ConsecutiveRunnable(
                            name=prefix_peak_calling,
                            working_directory=self.genome_directory,
                            debug=self.debug))
                    executable_peak_calling = self.set_stage_runnable(
                        stage=stage_peak_calling,
                        runnable=runnable_peak_calling)

                    t_file_path_list = list()
                    for t_sample in chipseq_comparison.t_samples:
                        t_file_path_list.append(self.get_file_path_alignment(sample_name=t_sample.name).sample_bam)
                        executable_peak_calling.dependencies.append(
                            self.get_prefix_alignment(sample_name=t_sample.name))

                    c_file_path_list = list()
                    for c_sample in chipseq_comparison.c_samples:
                        c_file_path_list.append(self.get_file_path_alignment(sample_name=c_sample.name).sample_bam)
                        executable_peak_calling.dependencies.append(
                            self.get_prefix_alignment(sample_name=c_sample.name))

                    # Add a RunnableStep to create the output directory.

                    runnable_step = RunnableStepMakeDirectory(
                        name='make_directory',
                        directory_path=file_path_peak_calling.output_directory)
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    # Add a RunnableStep for MACS2 call peak.

                    runnable_step = RunnableStep(
                        name='macs2_call_peak',
                        program='macs2',
                        sub_command=Command(program='callpeak'))
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    # Read RunnableStep options from configuration sections:
                    # [bsf.analyses.chipseq.ChIPSeq.macs2_call_peak]
                    # [bsf.analyses.chipseq.ChIPSeq.macs2_call_peak.callpeak]
                    # [bsf.analyses.chipseq.ChIPSeq.macs2_call_peak.callpeak.{factor}]
                    self.set_runnable_step_configuration(runnable_step=runnable_step, tag=chipseq_comparison.factor)

                    sub_command = runnable_step.sub_command

                    sub_command.add_option_multi_long(
                        key='treatment',
                        value=' '.join(t_file_path_list))

                    if c_file_path_list:
                        # The control (input) samples are optional.
                        sub_command.add_option_multi_long(
                            key='control',
                            value=' '.join(c_file_path_list))
                    # --format ["AUTO"]
                    # --gsize ["hs"] Genome size
                    # The 'gsize' option has to be specified via the configuration.ini file in section
                    # [bsf.analyses.chipseq.ChIPSeq.macs2_call_peak.callpeak].
                    # --tsize [null] Tag size or read length
                    # --keep-dup [1]

                    # Output arguments
                    # --outdir [.]
                    sub_command.add_option_long(
                        key='outdir',
                        value=file_path_peak_calling.output_directory)
                    # --name ["NA"]
                    # MACS2 can cope with directories specified in the --name option, but
                    # the resulting R script has them set too. Hence the R script has to be started
                    # from the genome_directory. However, the R script needs re-writing anyway, because
                    # it would be better to use the PNG rather than the PDF device for plotting.
                    sub_command.add_option_long(key='name', value=prefix_peak_calling)

                    # --bdg [False]
                    sub_command.add_switch_long(key='bdg')
                    # --verbose [2]
                    # --trackline [False]
                    # --SPMR [False]
                    sub_command.add_switch_long(key='SPMR')

                    # Shifting model arguments
                    # --nomodel [False]
                    # --shift [0]
                    # --extsize [200]
                    # --bw [300]
                    # --mfold [5 50]
                    # --fix-bimodal [False]

                    # Peak calling arguments
                    # --qvalue [0.05]
                    # --pvalue [null]
                    # --scale-to ["small"]
                    # --ratio [ignore]
                    # --down-sample [False]
                    # --seed [null]
                    sub_command.add_option_long(
                        key='tempdir',
                        value=runnable_peak_calling.temporary_directory_path(absolute=False))
                    # --nolambda [null]
                    # --slocal [1000]
                    # --llocal [10000]
                    # --max-gap [null]
                    # --min-length [null]
                    # --broad [False]
                    # --broad-cutoff [0.1]
                    # --cutoff-analysis [False]

                    # Post-processing options
                    # --call-summits [False]
                    # --fe-cutoff [1.0]

                    # Other options
                    # --buffer-size [100000]

                    # Add a RunnableStep to compare bedGraph files.

                    runnable_step = RunnableStep(
                        name='macs2_bdg_cmp',
                        program='macs2',
                        sub_command=Command(program='bdgcmp'))
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    # Read RunnableStep options from configuration sections:
                    # [bsf.analyses.chipseq.ChIPSeq.macs2_bdg_cmp]
                    # [bsf.analyses.chipseq.ChIPSeq.macs2_bdg_cmp.bdgcmp]
                    self.set_runnable_step_configuration(runnable_step=runnable_step)

                    sub_command = runnable_step.sub_command

                    # --tfile
                    sub_command.add_option_long(
                        key='tfile',
                        value=file_path_peak_calling.treatment_bdg)

                    # --cfile
                    sub_command.add_option_long(
                        key='cfile',
                        value=file_path_peak_calling.control_bdg)

                    # --scaling-factor [1.0]
                    # --pseudocount [0.0]
                    sub_command.add_option_long(key='pseudocount', value='0.00001')

                    # --method [ppois] i.e. Poisson Pvalue -log10(pvalue), which yields data on a logarithmic scale
                    sub_command.add_option_multi_long(key='method', value='ppois subtract logFE')

                    # --outdir [.]
                    sub_command.add_option_long(key='outdir', value=file_path_peak_calling.output_directory)
                    # --o-prefix [null]
                    sub_command.add_option_long(key='o-prefix', value=prefix_peak_calling)
                    # --ofile [] Mutually exclusive with --o-prefix, must correspond to --method.

                    # Add a RunnableStep to process MACS2 output.

                    runnable_step = RunnableStep(
                        name='process_macs2',
                        program='bsf_chipseq_process_macs2.bash')
                    runnable_peak_calling.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.arguments.append(prefix_peak_calling)
                    runnable_step.arguments.append(self.genome_sizes_path)

                    if os.path.exists(os.path.join(self.genome_directory, file_path_peak_calling.narrow_peaks_bb)):
                        executable_peak_calling.submit = False

            return

        def run_create_diff_bind_jobs():
            """Create Bioconductor DiffBind jobs.
            """
            # Firstly, organise the ChIPSeqComparison objects by comparison name.
            for comparison_dict in self._comparison_dict.values():
                # Secondly, organise the ChIPSeqComparison objects by factor.
                for chipseq_comparison in comparison_dict.values():
                    if not chipseq_comparison.diff_bind:
                        continue
                    if chipseq_comparison.name not in self._factor_dict:
                        self._factor_dict[chipseq_comparison.name] = dict()
                    if chipseq_comparison.factor not in self._factor_dict[chipseq_comparison.name]:
                        self._factor_dict[chipseq_comparison.name][chipseq_comparison.factor] = list()
                    self._factor_dict[chipseq_comparison.name][chipseq_comparison.factor].append(chipseq_comparison)

            for comparison_name in sorted(self._factor_dict):
                for factor_name in sorted(self._factor_dict[comparison_name]):
                    if self.debug > 0:
                        print('chipseq factor_name:', factor_name)

                    # No comparison for less than two items.
                    if len(self._factor_dict[comparison_name][factor_name]) < 2:
                        continue

                    # Create a directory per factor.

                    prefix_diff_bind = self.get_prefix_diff_bind(
                        comparison_name=comparison_name,
                        factor_name=factor_name)

                    file_path_diff_bind = self.get_file_path_diff_bind(
                        comparison_name=comparison_name,
                        factor_name=factor_name)

                    job_dependencies = list()

                    # Create a new ChIPSeq DiffBind Sheet per factor.

                    dbs = ChIPSeqDiffBindSheet(
                        file_path=os.path.join(
                            self.genome_directory,
                            file_path_diff_bind.sample_annotation_sheet))

                    if self.debug > 0:
                        print('ChIPSeqDiffBindSheet file_path:', dbs.file_path)

                    for chipseq_comparison in sorted(
                            self._factor_dict[comparison_name][factor_name],
                            key=lambda item: item.c_name):
                        if not chipseq_comparison.diff_bind:
                            continue

                        t_file_path_list = list()
                        for t_sample in chipseq_comparison.t_samples:
                            t_file_path_list.append(self.get_file_path_alignment(sample_name=t_sample.name).sample_bam)

                        c_file_path_list = list()
                        for c_sample in chipseq_comparison.c_samples:
                            c_file_path_list.append(self.get_file_path_alignment(sample_name=c_sample.name).sample_bam)
                        else:
                            c_file_path_list.append('')

                        file_path_peak_calling = chipseq_comparison.get_file_path_peak_calling()

                        job_dependencies.append(chipseq_comparison.get_prefix_peak_calling())

                        dbs.row_dicts.append({
                            'SampleID': chipseq_comparison.t_name,
                            'Tissue': chipseq_comparison.tissue,
                            'Factor': chipseq_comparison.factor,
                            'Condition': chipseq_comparison.condition,
                            'Treatment': chipseq_comparison.treatment,
                            'Replicate': str(chipseq_comparison.replicate),
                            'bamReads': t_file_path_list[0],
                            'bamControl': c_file_path_list[0],
                            'ControlID': chipseq_comparison.c_name,
                            'Peaks': file_path_peak_calling.peaks_xls,
                            'PeakCaller': 'macs',
                            'PeakFormat': 'macs',
                            # 'ScoreCol': '',
                            # 'LowerBetter': '',
                            # 'Counts': '',
                        })

                    dbs.to_file_path()

                    runnable_diff_bind = self.add_runnable_consecutive(
                        runnable=ConsecutiveRunnable(
                            name=prefix_diff_bind,
                            working_directory=self.genome_directory,
                            debug=self.debug))
                    executable_diff_bind = self.set_stage_runnable(
                        stage=stage_diff_bind,
                        runnable=runnable_diff_bind)
                    executable_diff_bind.dependencies.extend(job_dependencies)

                    # Add a RunnableStep for Bioconductor DiffBind.

                    runnable_step = RunnableStep(
                        name='diff_bind',
                        program='bsf_chipseq_diffbind.R')
                    runnable_diff_bind.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_option_long(
                        key='comparison',
                        value=comparison_name)
                    runnable_step.add_option_long(
                        key='factor',
                        value=factor_name)
                    runnable_step.add_option_long(
                        key='sample-annotation',
                        value=file_path_diff_bind.sample_annotation_sheet)
                    runnable_step.add_option_long(
                        key='black-list',
                        value=self.genome_black_list)
                    runnable_step.add_option_long(
                        key='genome-version',
                        value=self.genome_version)

                    # Add a RunnableStep for Bioconductor ChIPpeakAnno

                    if self.transcriptome_gtf_path and self.transcriptome_txdb_path:
                        runnable_step = RunnableStep(
                            name='annotation',
                            program='bsf_chipseq_annotation.R')
                        runnable_diff_bind.add_runnable_step(runnable_step=runnable_step)

                        runnable_step.add_option_long(
                            key='comparison',
                            value=comparison_name)
                        runnable_step.add_option_long(
                            key='factor',
                            value=factor_name)
                        runnable_step.add_option_long(
                            key='gtf-reference',
                            value=self.transcriptome_gtf_path)
                        runnable_step.add_option_long(
                            key='txdb-path',
                            value=self.transcriptome_txdb_path)

                    # Add a RunnableStep for Bioconductor ChIPQC.

                    runnable_step = RunnableStep(
                        name='chipqc',
                        program='bsf_chipseq_chipqc.R')
                    runnable_diff_bind.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_option_long(
                        key='comparison',
                        value=comparison_name)
                    runnable_step.add_option_long(
                        key='factor',
                        value=factor_name)
                    runnable_step.add_option_long(
                        key='threads',
                        value=str(stage_diff_bind.threads))

                    # Add a RunnableStep to run UCSC bedSort on the DiffBind consensus peak set.

                    runnable_step = RunnableStep(name='bed_sort', program='bedSort')
                    runnable_diff_bind.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.arguments.append(file_path_diff_bind.peak_set_bed)
                    runnable_step.arguments.append(file_path_diff_bind.peak_set_bed)

                    # Add a RunnableStep to run UCSC BedToBigBed on the DiffBind consensus peak set.

                    runnable_step = RunnableStep(name='bed_to_big_bed', program='bedToBigBed')
                    runnable_diff_bind.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_option_pair_short(key='type', value='bed6')

                    runnable_step.arguments.append(file_path_diff_bind.peak_set_bed)
                    runnable_step.arguments.append(self.genome_sizes_path)
                    runnable_step.arguments.append(file_path_diff_bind.peak_set_bb)

            return

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # ChIPSeq requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A ' + self.name + " requires a 'transcriptome_version' configuration option.")

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a valid 'transcriptome_version' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[ChIPSeq.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # Get the comparison annotation sheet before calling the run() method of the Analysis super-class.

        if self.comparison_path:
            self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)
            if not os.path.exists(self.comparison_path):
                raise Exception('Comparison annotation file ' + repr(self.comparison_path) + ' does not exist.')
        else:
            self.comparison_path = self.get_annotation_file(prefix_list=[ChIPSeq.prefix], suffix='comparisons.csv')
            if not self.comparison_path:
                raise Exception('No suitable comparison annotation file in the current working directory.')

        if not self.transcriptome_gtf_path:
            self.transcriptome_gtf_path = StandardFilePath.get_resource_transcriptome_gtf(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='none',
                basic=True,
                absolute=True)

        if not self.transcriptome_txdb_path:
            self.transcriptome_txdb_path = StandardFilePath.get_resource_transcriptome_txdb(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='none',
                basic=True,
                absolute=True)

        super(ChIPSeq, self).run()

        if not self.colour_default:
            self.colour_default = '0,0,0'

        if self.colour_dict is None:
            raise Exception('A ' + self.name + " requires a TrackHubColours configuration section.")

        if not self.factor_default:
            self.factor_default = 'OTHER'

        if not self.genome_black_list:
            self.genome_black_list = StandardFilePath.get_resource_genome_black_list(
                genome_version=self.genome_version)

        if self.genome_sizes_path:
            self.genome_sizes_path = self.configuration.get_absolute_path(file_path=self.genome_sizes_path)
        else:
            self.genome_sizes_path = StandardFilePath.get_resource_genome_fasta_index(
                genome_version=self.genome_version)

        if not self.genome_effective_size:
            self.genome_effective_size = Genome.get_effective_size(genome_version=self.genome_version)

        stage_alignment = self.get_stage(name=self.get_stage_name_alignment())
        stage_peak_calling = self.get_stage(name=self.get_stage_name_peak_calling())
        stage_diff_bind = self.get_stage(name=self.get_stage_name_diff_bind())

        run_read_comparisons()

        run_create_alignment_jobs()
        if self._macs_version == 1:
            run_create_macs1_jobs()
        elif self._macs_version == 2:
            run_create_macs2_jobs()

        run_create_diff_bind_jobs()

        return

    def report(self):
        """Create a C{bsf.analyses.chipseq.ChIPSeq} report in HTML format and a UCSC Genome Browser Track Hub.
        """

        # contrast_field_names = ['', 'Group1', 'Members1', 'Group2', 'Members2', 'DB.edgeR']

        def report_html_1():
            """Private function to create a HTML report for MACS1.
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # Write a HTML document.

            str_list: List[str] = list()

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.extend(self.get_html_genome(genome_version=self.genome_version))
            str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
            str_list.append('\n')

            str_list.append('<p>\n')
            str_list.append('Next-Generation Sequencing reads are aligned with the short read aligner\n')
            str_list.append('<strong>')
            str_list.append('<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>')
            str_list.append('</strong>,\n')
            str_list.append('before peaks are called with ')
            str_list.append('<a href="https://macs3-project.github.io/MACS/">MACS 1.4</a>\n')
            str_list.append('on an enrichment and background sample pair.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<p id="bowtie2_report">\n')
            str_list.append('Please see the ')
            str_list.append('<a href="' + Bowtie2.prefix + '_report.html">')
            str_list.append(self.project_name + ' ' + Bowtie2.name)
            str_list.append('</a> report for quality plots and ')
            str_list.append('a link to alignment visualisation in the UCSC Genome Browser.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Construct an automatic UCSC Track Hub link.

            str_list.append('<p id="ucsc_track_hub">\n')
            str_list.append('View Bowtie2 <strong>read alignment</strong> tracks for each sample\n')
            str_list.append('in their genomic context via the project-specific\n')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('.')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Peak Calling section.

            str_list.append('<h2 id="peak_calling">Peak Calling</h2>\n')
            str_list.append('\n')

            str_list.append('<table id="peak_calling_table">\n')
            str_list.append('<thead>\n')

            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Peaks</th>\n')
            str_list.append('<th>Negative Peaks</th>\n')
            str_list.append('<th>R Model</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in sorted(self._comparison_dict):
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    t_file_path_list = list()
                    for t_sample in chipseq_comparison.t_samples:
                        t_file_path_list.append(self.get_file_path_alignment(sample_name=t_sample.name).sample_bam)

                    c_file_path_list = list()
                    for c_sample in chipseq_comparison.c_samples:
                        c_file_path_list.append(self.get_file_path_alignment(sample_name=c_sample.name).sample_bam)
                    else:
                        c_file_path_list.append('')

                    for t_sample in chipseq_comparison.t_samples:
                        t_paired_reads_dict = t_sample.get_all_paired_reads(
                            replicate_grouping=self.replicate_grouping)
                        for t_paired_reads_name in sorted(t_paired_reads_dict):
                            for c_sample in chipseq_comparison.c_samples:
                                c_paired_reads_dict = c_sample.get_all_paired_reads(
                                    replicate_grouping=self.replicate_grouping)
                                for c_paired_reads_name in sorted(c_paired_reads_dict):
                                    str_list.append('<tr>\n')
                                    # Comparison name
                                    str_list.append('<td>' + comparison_name + '</td>\n')
                                    # Peaks
                                    str_list.append('<td>')
                                    str_list.append('<a href="' + t_paired_reads_name + '__' + c_paired_reads_name +
                                                    '_peaks.xls">')
                                    str_list.append('Peaks ' + t_paired_reads_name + ' versus ' + c_paired_reads_name)
                                    str_list.append('</a>')
                                    str_list.append('</td>\n')
                                    # Negative Peaks
                                    str_list.append('<td>')
                                    str_list.append('<a href="' + t_paired_reads_name + '__' + c_paired_reads_name +
                                                    '_negative_peaks.xls">')
                                    str_list.append('Negative peaks ' + t_paired_reads_name + ' versus ' +
                                                    c_paired_reads_name)
                                    str_list.append('</a>')
                                    str_list.append('</td>\n')
                                    # R Model
                                    str_list.append('<td>')
                                    str_list.append('<a href="' + t_paired_reads_name + '__' + c_paired_reads_name +
                                                    '_model.r">')
                                    str_list.append('R model')
                                    str_list.append('</a>')
                                    str_list.append('</td>\n')

                                    str_list.append('</tr>\n')
                                    str_list.append('\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

            return

        def report_html_2():
            """Private function to create a HTML report for MACS2.
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # Write a HTML document.

            str_list: List[str] = list()

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.extend(self.get_html_genome(genome_version=self.genome_version))
            str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
            str_list.append('\n')

            str_list.append('<p>\n')
            str_list.append('Next-Generation Sequencing reads are aligned with the short read aligner\n')
            str_list.append('<strong>')
            str_list.append('<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>')
            str_list.append('</strong>,\n')
            str_list.append('before peaks are called with ')
            str_list.append('<a href="https://macs3-project.github.io/MACS/">MACS 2</a>\n')
            str_list.append('on an enrichment and background sample pair.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<p id="bowtie2_report">\n')
            str_list.append('Please see the ')
            str_list.append('<a href="' + Bowtie2.prefix + '_report.html">')
            str_list.append(self.project_name + ' ' + Bowtie2.name)
            str_list.append('</a> report for quality plots and ')
            str_list.append('a link to alignment visualisation in the UCSC Genome Browser.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Construct an automatic UCSC Track Hub link.

            str_list.append('<p id="ucsc_track_hub">\n')
            str_list.append('View Bowtie2 <strong>read alignment</strong> tracks for each sample\n')
            str_list.append('and <strong>ChIP-seq signal</strong> and background tracks\n')
            str_list.append('in their genomic context via the project-specific\n')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('.')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Alignment section

            str_list.append('<h2 id="alignment">Alignment</h2>\n')
            str_list.append('<p>')
            str_list.append('Bowtie2 alignments were post-processed by removing reads from the (ENCODE)\n')
            str_list.append('blacklist regions, before coverage tracks were normalised to one-fold genome coverage.')
            str_list.append('</p>\n')

            str_list.append('<table id="alignment_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Sample</th>\n')
            str_list.append('<th>BAM</th>\n')
            str_list.append('<th>BAI</th>\n')
            # str_list.append('<th>MD5</th>\n')
            str_list.append('<th>bigWig</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for sample in self.sample_list:
                file_path_alignment = self.get_file_path_alignment(sample_name=sample.name)
                str_list.append('<tr>\n')

                # Sample
                str_list.append('<td class="left">')
                str_list.append(sample.name)
                str_list.append('</td>\n')
                # BAM
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_alignment.sample_bam + '">')
                str_list.append('<abbr title="Binary Alignment/Map">BAM</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # BAI
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_alignment.sample_bai + '">')
                str_list.append('<abbr title="Binary Alignment/Map Index">BAI</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # MD5
                # str_list.append('<td class="center">')
                # str_list.append('<a href="' + file_path_alignment.sample_md5 + '">')
                # str_list.append('<abbr title="Message Digest 5 Checksum">MD5</abbr>')
                # str_list.append('</a>')
                # str_list.append('</td>\n')
                # bigWig
                str_list.append('<td>')
                str_list.append('<a href="' + file_path_alignment.coverage_bw + '">')
                str_list.append('<abbr title="UCSC Big Wiggle signal track">bigWig</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')

                str_list.append('</tr>\n')
                str_list.append('\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Peak Calling section.

            str_list.append('<h2 id="peak_calling">Peak Calling</h2>\n')
            str_list.append('\n')

            str_list.append('<table id="peak_calling_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Pair</th>\n')
            str_list.append('<th>Peaks</th>\n')
            str_list.append('<th>Peak Model</th>\n')
            str_list.append('<th>Cross Correlation</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in sorted(self._comparison_dict):
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    file_path_peak_calling = chipseq_comparison.get_file_path_peak_calling()

                    str_list.append('<tr>\n')
                    # Comparison name
                    str_list.append('<td>' + comparison_name + '</td>\n')
                    # Pair
                    str_list.append('<td>')
                    str_list.append('<strong>' + chipseq_comparison.t_name + '</strong>')
                    if chipseq_comparison.c_name:
                        str_list.append(' versus ')
                        str_list.append('<strong>' + chipseq_comparison.c_name + '</strong>')
                    str_list.append('</td>\n')
                    # Peaks
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_peak_calling.peaks_tsv + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                    str_list.append('</td>\n')
                    # Peak Model
                    str_list.append('<td>')
                    if os.path.exists(path=os.path.join(self.genome_directory, file_path_peak_calling.model_pdf)):
                        str_list.append('<a href="' + file_path_peak_calling.model_pdf + '">')
                        str_list.append('<img')
                        str_list.append(' alt="MACS2 peak model for treatment ' + chipseq_comparison.t_name + '"')
                        str_list.append(' src="' + file_path_peak_calling.model_0_png + '"')
                        str_list.append(' height="80" width="80" />')
                        str_list.append('</a>')
                    str_list.append('</td>')
                    # Cross-Correlation
                    str_list.append('<td>')
                    if os.path.exists(path=os.path.join(self.genome_directory, file_path_peak_calling.model_pdf)):
                        str_list.append('<a href="' + file_path_peak_calling.model_pdf + '">')
                        str_list.append('<img')
                        str_list.append(' alt="MACS2 cross-correlation for treatment ' +
                                        chipseq_comparison.t_name + '"')
                        str_list.append(' src="' + file_path_peak_calling.model_1_png + '"')
                        str_list.append(' height="80" width="80" />')
                        str_list.append('</a>')
                    str_list.append('</td>')

                    # R Model
                    # str_list.append('<td><a href="' + file_path_peak_calling.model_r + '">')
                    # str_list.append('R model')
                    # str_list.append('</a></td>\n')

                    str_list.append('</tr>\n')
                    str_list.append('\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Differential Binding Analysis section.

            str_list.append('<h2 id="differential_binding">Differential Binding Analysis</h2>\n')
            str_list.append('\n')

            str_list.append('<table id="differential_binding_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Factor and Contrast</th>\n')
            str_list.append('<th>Correlation Peak Caller</th>\n')
            str_list.append('<th>Correlation Peak Counts</th>\n')
            str_list.append('<th>Correlation Analysis</th>\n')
            str_list.append('<th>MA Plot</th>\n')
            str_list.append('<th>Scatter Plot</th>\n')
            str_list.append('<th>PCA Plot</th>\n')
            str_list.append('<th>Box Plot</th>\n')
            str_list.append('<th>Chromosome Regions</th>\n')
            str_list.append('<th>Differentially Bound Sites</th>\n')
            str_list.append('<th>ChIPQC Report</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in sorted(self._factor_dict):
                for factor_name in sorted(self._factor_dict[comparison_name]):
                    file_path_diff_bind = self.get_file_path_diff_bind(
                        comparison_name=comparison_name,
                        factor_name=factor_name)

                    file_path_chipqc = self.get_file_path_chipqc(
                        comparison_name=comparison_name,
                        factor_name=factor_name)

                    str_list.append('<tr>\n')

                    # Comparison name
                    str_list.append('<td>' + comparison_name + '</td>\n')
                    # ChIP Factor
                    str_list.append('<td><strong>' + factor_name + '</strong></td>\n')
                    # Correlation heat map of peak caller scores
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_diff_bind.correlation_peak_caller_score_pdf + '">')
                    str_list.append('<img')
                    str_list.append(' alt="DiffBind correlation analysis for factor ' + factor_name + '"')
                    str_list.append(' src="' + file_path_diff_bind.correlation_peak_caller_score_png + '"')
                    str_list.append(' height="80" width="80" />')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Correlation heat map of counts
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_diff_bind.correlation_read_counts_pdf + '">')
                    str_list.append('<img')
                    str_list.append(' alt="DiffBind correlation analysis for factor ' + factor_name + '"')
                    str_list.append(' src="' + file_path_diff_bind.correlation_read_counts_png + '"')
                    str_list.append(' height="80" width="80" />')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Correlation heat map of differential binding analysis
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_diff_bind.correlation_analysis_pdf + '">')
                    str_list.append('<img')
                    str_list.append(' alt="DiffBind correlation analysis for factor ' + factor_name + '"')
                    str_list.append(' src="' + file_path_diff_bind.correlation_analysis_png + '"')
                    str_list.append(' height="80" width="80" />')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    str_list.append('<td></td>\n')  # MA Plot
                    str_list.append('<td></td>\n')  # Scatter Plot

                    # Score-based PCA Plot
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_diff_bind.pca_pdf + '">')
                    str_list.append('<img')
                    str_list.append(' alt="Score-based PCA plot for factor ' + factor_name + '"')
                    str_list.append(' src="' + file_path_diff_bind.pca_png + '"')
                    str_list.append(' height="80" width="80" />')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    str_list.append('<td></td>\n')  # Box Plot

                    # Chromosome Regions
                    str_list.append('<td>')
                    if os.path.exists(os.path.join(self.genome_directory, file_path_diff_bind.regions_png)):
                        str_list.append('<a href="' + file_path_diff_bind.regions_pdf + '">')
                        str_list.append('<img')
                        str_list.append(' alt="Chromosome region plot for factor ' + factor_name + '"')
                        str_list.append(' src="' + file_path_diff_bind.regions_png + '"')
                        str_list.append(' height="80" width="80" />')
                        str_list.append('</a>')
                        str_list.append('<br />')
                        str_list.append('<a href="' + file_path_diff_bind.regions_tsv + '">')
                        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                        str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Differentially Bound Sites
                    str_list.append('<td>')
                    if os.path.exists(os.path.join(self.genome_directory, file_path_diff_bind.genes_complete_tsv)):
                        str_list.append('<a href="' + file_path_diff_bind.genes_complete_tsv + '">')
                        str_list.append('<abbr title="Tab-Separated Value">complete</abbr>')
                        str_list.append('</a>')
                    str_list.append('</td>\n')

                    # ChIPQC HTML Report
                    str_list.append('<td>')
                    if os.path.exists(os.path.join(self.genome_directory, file_path_chipqc.report_html)):
                        str_list.append('<a href="' + file_path_chipqc.report_html + '">')
                        str_list.append('<abbr title="Hypertext Markup Language">HTML</abbr>')
                        str_list.append('</a>')
                    str_list.append('</td>\n')

                    str_list.append('</tr>\n')

                    # Read the file of contrasts ...

                    file_path = os.path.join(self.genome_directory, file_path_diff_bind.contrasts_csv)

                    if not os.path.exists(file_path):
                        warnings.warn(
                            'Contrasts table does not exist: ' + file_path,
                            UserWarning)
                        continue

                    annotation_sheet = AnnotationSheet.from_file_path(
                        file_path=file_path,
                        file_type='excel')

                    for row_dict in annotation_sheet.row_dicts:
                        # DiffBind2 uses a 'Group1' variable, DiffBind3 just 'Group'.
                        group1_name = 'Group1'
                        if 'Group' in row_dict:
                            group1_name = 'Group'

                        file_path_diff_bind_contrast = FilePathDiffBindContrast(
                            prefix=self.get_prefix_diff_bind(
                                comparison_name=comparison_name,
                                factor_name=factor_name),
                            group_1=row_dict[group1_name],
                            group_2=row_dict['Group2'])

                        str_list.append('<tr>\n')

                        str_list.append('<td></td>\n')  # Comparison name
                        str_list.append('<td>' + row_dict[group1_name] + ' vs ' + row_dict['Group2'] + '</td>\n')
                        str_list.append('<td></td>\n')  # Correlation heat map of peak caller scores
                        str_list.append('<td></td>\n')  # Correlation heat map of counts
                        str_list.append('<td></td>\n')  # Correlation heat map of differential binding analysis

                        # MA Plot
                        str_list.append('<td>')
                        str_list.append('<a href="' + file_path_diff_bind_contrast.ma_plot_pdf + '">')
                        str_list.append('<img')
                        str_list.append(' alt="DiffBind MA plot for factor {} contrast {} vs {}"'.format(
                            factor_name, row_dict[group1_name], row_dict['Group2']))
                        str_list.append(' src="' + file_path_diff_bind_contrast.ma_plot_png + '"')
                        str_list.append(' height="80" width="80" />')
                        str_list.append('</a>')
                        str_list.append('</td>\n')

                        # Scatter Plot
                        str_list.append('<td>')
                        str_list.append('<a href="' + file_path_diff_bind_contrast.scatter_plot_pdf + '">')
                        str_list.append('<img')
                        str_list.append(' alt="DiffBind scatter for factor {} contrast {} vs {}"'.format(
                            factor_name, row_dict[group1_name], row_dict['Group2']))
                        str_list.append(' src="' + file_path_diff_bind_contrast.scatter_plot_png + '"')
                        str_list.append(' height="80" width="80" />')
                        str_list.append('</a>')
                        str_list.append('</td>\n')

                        # Principal Component Analysis Plot (optional)
                        str_list.append('<td>')
                        if os.path.exists(
                                os.path.join(
                                    self.genome_directory,
                                    file_path_diff_bind_contrast.pca_plot_png)):
                            str_list.append('<a href="' + file_path_diff_bind_contrast.pca_plot_pdf + '">')
                            str_list.append('<img')
                            str_list.append(' alt="DiffBind PCA plot for factor {} contrast {} vs {}"'.format(
                                factor_name, row_dict[group1_name], row_dict['Group2']))
                            str_list.append(' src="' + file_path_diff_bind_contrast.pca_plot_png + '"')
                            str_list.append(' height="80" width="80" />')
                            str_list.append('</a>')
                        str_list.append('</td>\n')

                        # Box Plot (optional)
                        str_list.append('<td>')
                        if os.path.exists(
                                os.path.join(self.genome_directory, file_path_diff_bind_contrast.box_plot_png)):
                            str_list.append('<a href="' + file_path_diff_bind_contrast.box_plot_pdf + '">')
                            str_list.append('<img')
                            str_list.append(' alt="DiffBind Box plot for factor {} contrast {} vs {}"'.format(
                                factor_name, row_dict[group1_name], row_dict['Group2']))
                            str_list.append(' src="' + file_path_diff_bind_contrast.box_plot_png + '"')
                            str_list.append(' height="80" width="80" />')
                            str_list.append('</a>')
                        str_list.append('</td>\n')

                        # Chromosome Regions
                        str_list.append('<td>')
                        # Insert the chromosome regions plot
                        if os.path.exists(os.path.join(
                                self.genome_directory,
                                file_path_diff_bind_contrast.regions_png)):
                            str_list.append('<a href="' + file_path_diff_bind_contrast.regions_pdf + '">')
                            str_list.append('<img')
                            str_list.append(' alt="Chromosome region plot for factor {} contrast {} vs {}"'.format(
                                factor_name, row_dict[group1_name], row_dict['Group2']))
                            str_list.append(' src="' + file_path_diff_bind_contrast.regions_png + '"')
                            str_list.append(' height="80" width="80" />')
                            str_list.append('</a>')
                            str_list.append('<br />')
                            str_list.append('<a href="' + file_path_diff_bind_contrast.regions_tsv + '">')
                            str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                            str_list.append('</a>')
                        str_list.append('</td>\n')

                        # Differentially Bound Sites
                        str_list.append('<td>')
                        # Try the annotated report (*_suffix_genes_complete.tsv) first,
                        # then the unannotated one (_suffix.tsv).
                        if os.path.exists(
                                os.path.join(self.genome_directory, file_path_diff_bind_contrast.genes_complete_tsv)):
                            str_list.append('<a href="' + file_path_diff_bind_contrast.genes_complete_tsv + '">')
                            str_list.append('<abbr title="Tab-Separated Value">complete</abbr>')
                            str_list.append('</a>')
                            str_list.append('<br />\n')
                            str_list.append('<a href="' + file_path_diff_bind_contrast.genes_significant_tsv + '">')
                            str_list.append('<abbr title="Tab-Separated Value">significant</abbr>')
                            str_list.append('</a>')
                        elif os.path.exists(
                                os.path.join(self.genome_directory, file_path_diff_bind_contrast.peaks_tsv)):
                            str_list.append('<a href="' + file_path_diff_bind_contrast.peaks_tsv + '">')
                            str_list.append('<abbr title="Tab-Separated Value">raw</abbr>')
                            str_list.append('</a>')
                        str_list.append('</td>\n')

                        # ChIPQC
                        str_list.append('<td></td>\n')

                        str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

            return

        def report_hub_1():
            """Private function to create a UCSC Track Hub for MACS1.
            """

            str_list: List[str] = list()

            chipseq_comparison_key_list: List[str] = list()

            for comparison_name in sorted(self._comparison_dict):
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    # Since a particular peak calling pair may exists in more than one comparison group,
                    # make sure only unique pairs get listed in the track database definition.
                    if chipseq_comparison.get_key() in chipseq_comparison_key_list:
                        continue
                    else:
                        chipseq_comparison_key_list.append(chipseq_comparison.get_key())

                    prefix = chipseq_comparison.get_prefix_peak_calling()

                    # Add UCSC trackDB entries for each treatment/control and absolute/normalised pair.

                    for treatment in (True, False):
                        if treatment:
                            state = 'treat'
                        else:
                            state = 'control'

                        for absolute in (True, False):
                            if absolute:
                                scaling = 'absolute'
                            else:
                                scaling = 'normalised'

                            #
                            # Add a UCSC trackDB entry for each bigWig file
                            #
                            # Common track settings
                            str_list.append('track ChIP_' + '_'.join((prefix, state, scaling)) + '\n')
                            str_list.append('type bigWig\n')
                            str_list.append('shortLabel ChIP_' + '_'.join((prefix, state, scaling)) + '\n')
                            str_list.append('longLabel ' + scaling.capitalize() +
                                            ' ChIP-Seq read counts for ' + state + ' of ' +
                                            chipseq_comparison.t_name + ' versus ' + chipseq_comparison.c_name + '\n')
                            str_list.append('bigDataUrl ')
                            str_list.append('/'.join((
                                prefix + '_MACS_wiggle', state,
                                '_'.join((prefix, state, 'afterfiting', 'all.bw')))) + '\n')
                            # str_list.append('html ...\n'
                            if treatment and not absolute:
                                str_list.append('visibility full\n')
                            else:
                                str_list.append('visibility hide\n')
                            # Common optional track settings
                            str_list.append('color ' + self.get_colour(factor=chipseq_comparison.factor) + '\n')
                            # bigWig - Signal graphing track settings
                            str_list.append('graphTypeDefault bar\n')
                            str_list.append('maxHeightPixels 100:60:20\n')
                            str_list.append('smoothingWindow off\n')
                            if absolute:
                                # Track with absolute scaling.
                                str_list.append('autoScale on\n')
                            else:
                                # Track with relative scaling.
                                str_list.append('autoScale off\n')
                                str_list.append('viewLimits 0:40\n')
                            str_list.append('\n')

                    #
                    # Add a UCSC trackDB entry for each bigBed peaks file
                    #
                    # Common track settings
                    str_list.append('track Peaks_' + prefix + '\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel Peaks_' + prefix + '\n')
                    str_list.append('longLabel ChIP-Seq peaks for ' +
                                    chipseq_comparison.t_name + ' versus ' + chipseq_comparison.c_name + '\n')
                    str_list.append('bigDataUrl ' + prefix + '_peaks.bb\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility pack\n')
                    # Common optional track settings
                    str_list.append('color ' + self.get_colour(factor=chipseq_comparison.factor) + '\n')
                    str_list.append('\n')

            # Add UCSC trackDB entries for each Bowtie2 BAM file.

            for sample in self.sample_list:
                file_path_alignment = self.get_file_path_alignment(sample_name=sample.name)
                # Common track settings
                str_list.append('track Alignment_' + sample.name + '\n')
                str_list.append('type bam\n')
                str_list.append('shortLabel Alignment_' + sample.name + '\n')
                str_list.append('longLabel Bowtie2 alignment of ' + sample.name + '\n')
                str_list.append('bigDataUrl ' + file_path_alignment.sample_bam + '\n')
                # str_list.append('html ...\n')
                str_list.append('visibility hide\n')
                # Common optional track settings
                # str_list.append('color ' + self.get_colour(factor=factor_name) + '\n')
                # bam/cram - Compressed Sequence Alignment track settings.
                str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        def report_hub_2():
            """Private function to create a UCSC Track Hub for MACS2.
            """
            # Composite tracks need tags and labels for all subGroupN entries that are defined.

            comparison_str = ''
            factor_dict = dict()

            for comparison_name in sorted(self._comparison_dict):
                comparison_str += ' ' + comparison_name + '=' + comparison_name
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    factor_dict[chipseq_comparison.factor] = True

            factor_str = ''
            for factor in factor_dict:
                factor_str += ' ' + factor + '=' + factor

            # 1. Composite track "alignment" (BAM)

            str_list_1: List[str] = list()

            # Common track settings
            str_list_1.append('track alignment\n')
            str_list_1.append('type bam\n')
            str_list_1.append('shortLabel QC Alignment\n')
            str_list_1.append('longLabel ChIP Read Alignment\n')
            # str_list_1.append('html ...\n')
            str_list_1.append('visibility hide\n')
            # Common optional track settings
            # bam/cram - Compressed Sequence Alignment track settings
            # Composite track settings
            str_list_1.append('compositeTrack on\n')
            str_list_1.append('allButtonPair on\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_1.append('centerLabelsDense on\n')
            # str_list_1.append('dragAndDrop subTracks\n')
            str_list_1.append('\n')

            # 2. Composite track "coverage" (bigWig)

            str_list_2: List[str] = list()

            # Common track settings
            str_list_2.append('track coverage\n')
            str_list_2.append('type bigWig\n')
            str_list_2.append('shortLabel QC Coverage\n')
            str_list_2.append('longLabel Normalised ChIP Read Alignment Coverage\n')
            # str_list_2.append('html ...\n')
            str_list_2.append('visibility hide\n')
            # Common optional track settings
            # bigWig - Signal graphing track settings
            str_list_2.append('alwaysZero off\n')
            str_list_2.append('autoScale off\n')
            str_list_2.append('graphTypeDefault bar\n')
            str_list_2.append('maxHeightPixels 100:60:20\n')
            # str_list_2.append('maxWindowToQuery 10000000\n')
            str_list_2.append('smoothingWindow 5\n')
            # str_list_2.append('transformFunc NONE\n')
            str_list_2.append('viewLimits 0:15\n')
            str_list_2.append('viewLimitsMax 0:40\n')
            str_list_2.append('windowingFunction maximum\n')
            # str_list_2.append('yLineMark <#>\n')
            # str_list_2.append('yLineOnOff on \n')
            # str_list_2.append('gridDefault on\n')
            # Composite track settings
            str_list_2.append('compositeTrack on\n')
            str_list_2.append('allButtonPair on\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_2.append('centerLabelsDense on\n')
            # str_list_2.append('dragAndDrop subTracks\n')
            str_list_2.append('\n')

            # 3. Composite track "background" (bigWig)

            str_list_3: List[str] = list()

            # Common track settings
            str_list_3.append('track background\n')
            str_list_3.append('type bigWig\n')
            str_list_3.append('shortLabel QC Background\n')
            str_list_3.append('longLabel ChIP Background Signal\n')
            # str_list_3.append('html ...\n')
            str_list_3.append('visibility hide\n')
            # Common optional track settings
            # bigWig - Signal graphing track settings
            str_list_3.append('alwaysZero off\n')
            str_list_3.append('autoScale off\n')
            str_list_3.append('graphTypeDefault bar\n')
            str_list_3.append('maxHeightPixels 100:60:20\n')
            # str_list_3.append('maxWindowToQuery 10000000\n')
            str_list_3.append('smoothingWindow 5\n')
            # str_list_3.append('transformFunc NONE\n')
            str_list_3.append('viewLimits 0:15\n')
            str_list_3.append('viewLimitsMax 0:40\n')
            str_list_3.append('windowingFunction maximum\n')
            # str_list_3.append('yLineMark <#>\n')
            # str_list_3.append('yLineOnOff on \n')
            # str_list_3.append('gridDefault on\n')
            # Composite track settings
            str_list_3.append('compositeTrack on\n')
            str_list_3.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_3.append('centerLabelsDense on\n')
            # str_list_3.append('dragAndDrop subTracks\n')
            str_list_3.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_3.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_3.append('dimensions dimX=comparison dimY=factor\n')
            str_list_3.append('sortOrder comparison=+ factor=+\n')
            str_list_3.append('\n')

            # 4. Composite track "enrichment" (bigWig)

            str_list_4: List[str] = list()

            # Common track settings
            str_list_4.append('track enrichment\n')
            str_list_4.append('type bigWig\n')
            str_list_4.append('shortLabel QC Enrichment\n')
            str_list_4.append('longLabel ChIP Enrichment Signal\n')
            # str_list_4.append('html ...\n')
            str_list_4.append('visibility hide\n')
            # Common optional track settings
            # bigWig - Signal graphing track settings
            str_list_4.append('alwaysZero off\n')
            str_list_4.append('autoScale off\n')
            str_list_4.append('graphTypeDefault bar\n')
            str_list_4.append('maxHeightPixels 100:60:20\n')
            # str_list_4.append('maxWindowToQuery 10000000\n')
            str_list_4.append('smoothingWindow 5\n')
            # str_list_4.append('transformFunc NONE\n')
            str_list_4.append('viewLimits 0:15\n')
            str_list_4.append('viewLimitsMax 0:40\n')
            str_list_4.append('windowingFunction maximum\n')
            # str_list_4.append('yLineMark <#>\n')
            # str_list_4.append('yLineOnOff on \n')
            # str_list_4.append('gridDefault on\n')
            # Composite track settings
            str_list_4.append('compositeTrack on\n')
            str_list_4.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_4.append('centerLabelsDense on\n')
            # str_list_4.append('dragAndDrop subTracks\n')
            str_list_4.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_4.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_4.append('dimensions dimX=comparison dimY=factor\n')
            str_list_4.append('sortOrder comparison=+ factor=+\n')
            str_list_4.append('\n')

            # 5. Composite track "intensity" (bigWig)

            str_list_5: List[str] = list()

            # Common track settings
            str_list_5.append('track intensity\n')
            str_list_5.append('type bigWig\n')
            str_list_5.append('shortLabel ChIP Intensity\n')
            str_list_5.append('longLabel ChIP Intensity\n')
            # str_list_5.append('html ...\n')
            str_list_5.append('visibility hide\n')
            # Common optional track settings
            # bigWig - Signal graphing track settings
            str_list_5.append('alwaysZero off\n')
            str_list_5.append('autoScale off\n')
            str_list_5.append('graphTypeDefault bar\n')
            str_list_5.append('maxHeightPixels 100:60:20\n')
            # str_list_5.append('maxWindowToQuery 10000000\n')
            str_list_5.append('smoothingWindow 5\n')
            # str_list_5.append('transformFunc NONE\n')
            str_list_5.append('viewLimits 0:15\n')
            str_list_5.append('viewLimitsMax 0:40\n')
            str_list_5.append('windowingFunction maximum\n')
            # str_list_5.append('yLineMark <#>\n')
            # str_list_5.append('yLineOnOff on \n')
            # str_list_5.append('gridDefault on\n')
            # Composite track settings
            str_list_5.append('compositeTrack on\n')
            str_list_5.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_5.append('centerLabelsDense on\n')
            # str_list_5.append('dragAndDrop subTracks\n')
            str_list_5.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_5.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_5.append('subGroup3 scale Scale' +
                              ' ppois=Poisson_pvalue logfe=Log10_fold_enrichment subtract=Subtraction\n')
            str_list_5.append('dimensions dimX=comparison dimY=factor dimA=scale\n')
            str_list_5.append('filterComposite dimA\n')
            str_list_5.append('sortOrder comparison=+ factor=+\n')
            str_list_5.append('\n')

            # 6. Composite track "narrow_peaks" (bigBed)

            str_list_6: List[str] = list()

            # Common track settings
            str_list_6.append('track narrow_peaks\n')
            str_list_6.append('type bigBed 6+4\n')
            str_list_6.append('shortLabel ChIP Narrow Peaks\n')
            str_list_6.append('longLabel ChIP Narrow Peaks\n')
            # str_list_6.append('html ...\n')
            str_list_6.append('visibility hide\n')
            # Common optional track settings
            # bigBed - Item or region track settings
            # Composite track settings
            str_list_6.append('compositeTrack on\n')
            str_list_6.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_6.append('centerLabelsDense on\n')
            # str_list_6.append('dragAndDrop subTracks\n')
            str_list_6.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_6.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_6.append('dimensions dimX=comparison dimY=factor\n')
            str_list_6.append('sortOrder comparison=+ factor=+\n')
            str_list_6.append('\n')

            # 7. Composite track "summits" (bigBed)

            str_list_7: List[str] = list()

            # Common track settings
            str_list_7.append('track summits\n')
            str_list_7.append('type bigBed 4+1\n')
            str_list_7.append('shortLabel ChIP Summits\n')
            str_list_7.append('longLabel ChIP Summits\n')
            # str_list_7.append('html ...\n')
            str_list_7.append('visibility hide\n')
            # Common optional track settings
            # bigBed - Item or region track settings
            # Composite track settings
            str_list_7.append('compositeTrack on\n')
            str_list_7.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_7.append('centerLabelsDense on\n')
            # str_list_7.append('dragAndDrop subTracks\n')
            str_list_7.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_7.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_7.append('dimensions dimX=comparison dimY=factor\n')
            str_list_7.append('sortOrder comparison=+ factor=+\n')
            str_list_7.append('\n')

            # 8. Composite track "broad_peaks" (bigBed)

            str_list_8: List[str] = list()

            # Common track settings
            str_list_8.append('track broad_peaks\n')
            str_list_8.append('type bigBed 6+3\n')
            str_list_8.append('shortLabel ChIP Broad Peaks\n')
            str_list_8.append('longLabel ChIP Broad Peaks\n')
            # str_list_8.append('html ...\n')
            str_list_8.append('visibility hide\n')
            # Common optional track settings
            # bigBed - Item or region track settings
            # Composite track settings
            str_list_8.append('compositeTrack on\n')
            str_list_8.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_8.append('centerLabelsDense on\n')
            # str_list_8.append('dragAndDrop subTracks\n')
            str_list_8.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_8.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_8.append('dimensions dimX=comparison dimY=factor\n')
            str_list_8.append('sortOrder comparison=+ factor=+\n')
            str_list_8.append('\n')

            # 9. Composite track "broad_peaks" (bigBed)

            str_list_9: List[str] = list()

            # Common track settings
            str_list_9.append('track gapped_peaks\n')
            str_list_9.append('type bigBed 12+3\n')
            str_list_9.append('shortLabel ChIP Gapped Peaks\n')
            str_list_9.append('longLabel ChIP Gapped Peaks\n')
            # str_list_9.append('html ...\n')
            str_list_9.append('visibility hide\n')
            # Common optional track settings
            # bigBed - Item or region track settings
            # Composite track settings
            str_list_9.append('compositeTrack on\n')
            str_list_9.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_9.append('centerLabelsDense on\n')
            # str_list_9.append('dragAndDrop subTracks\n')
            str_list_9.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_9.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_9.append('dimensions dimX=comparison dimY=factor\n')
            str_list_9.append('sortOrder comparison=+ factor=+\n')
            str_list_9.append('\n')

            # 10. Composite track "consensus_peaks" (bigBed)

            str_list_10: List[str] = list()

            # Common track settings
            str_list_10.append('track consensus_peaks\n')
            str_list_10.append('type bigBed 6\n')
            str_list_10.append('shortLabel ChIP Consensus Peaks\n')
            str_list_10.append('longLabel ChIP Consensus Peaks\n')
            # str_list_10.append('html ...\n')
            str_list_10.append('visibility hide\n')
            # Common optional track settings
            # bigBed - Item or region track settings
            # Composite track settings
            str_list_10.append('compositeTrack on\n')
            str_list_10.append('allButtonPair off\n')  # Has to be "off" to allow for configuration via a matrix.
            str_list_10.append('centerLabelsDense on\n')
            # str_list_10.append('dragAndDrop subTracks\n')
            str_list_10.append('subGroup1 comparison Comparison' + comparison_str + '\n')
            str_list_10.append('subGroup2 factor Factor' + factor_str + '\n')
            str_list_10.append('dimensions dimX=comparison dimY=factor\n')
            str_list_10.append('sortOrder comparison=+ factor=+\n')
            str_list_10.append('\n')

            # Add UCSC trackDB entries for each comparison.

            for comparison_name in sorted(self._comparison_dict):
                for comparison_pair in sorted(self._comparison_dict[comparison_name]):
                    chipseq_comparison = self._comparison_dict[comparison_name][comparison_pair]
                    factor_name = chipseq_comparison.factor
                    file_path_peak_calling = chipseq_comparison.get_file_path_peak_calling()

                    if chipseq_comparison.c_name:
                        prefix_short = '__'.join((chipseq_comparison.t_name, chipseq_comparison.c_name))
                        prefix_long = chipseq_comparison.t_name + ' versus ' + chipseq_comparison.c_name
                    else:
                        prefix_short = chipseq_comparison.t_name
                        prefix_long = chipseq_comparison.t_name

                    # Make track names unique by including the comparison name.

                    prefix_short = '_'.join((comparison_name, prefix_short))
                    prefix_long = ' '.join((comparison_name, prefix_long))

                    # Add a "background" sub-track for each NAME_control_lambda.bw file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.control_bw)):
                        # Common track settings
                        str_list_3.append('  track ' + prefix_short + '_background\n')
                        str_list_3.append('  ' + self.ucsc_hub_bigwig_info_signal_range(
                            file_path=file_path_peak_calling.control_bwi))
                        str_list_3.append('  shortLabel ' + prefix_short + '_background\n')
                        str_list_3.append('  longLabel ChIP background signal ' + prefix_long + '\n')
                        str_list_3.append('  bigDataUrl ' + file_path_peak_calling.control_bw + '\n')
                        # str_list_3.append('  html ...\n')
                        str_list_3.append('  visibility full\n')
                        # Common optional track settings
                        str_list_3.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigWig - Signal graphing track settings
                        # Composite track settings
                        str_list_3.append('  parent background on\n')
                        str_list_3.append('  subGroups comparison=' + comparison_name + ' factor=' + factor_name + '\n')
                        str_list_3.append('  \n')

                    # Add an "enrichment" sub-track for each NAME_treat_pileup.bw file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.treatment_bw)):
                        # Common track settings
                        str_list_4.append('  track ' + prefix_short + '_enrichment\n')
                        str_list_4.append('  ' + self.ucsc_hub_bigwig_info_signal_range(
                            file_path=file_path_peak_calling.treatment_bwi))
                        str_list_4.append('  shortLabel ' + prefix_short + '_enrichment\n')
                        str_list_4.append('  longLabel ChIP enrichment signal ' + prefix_long + '\n')
                        str_list_4.append('  bigDataUrl ' + file_path_peak_calling.treatment_bw + '\n')
                        # str_list_4.append('  html ...\n')
                        str_list_4.append('  visibility full\n')
                        # Common optional track settings
                        str_list_4.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigWig - Signal graphing track settings
                        # Composite track settings
                        str_list_4.append('  parent enrichment on\n')
                        str_list_4.append('  subGroups comparison=' + comparison_name + ' factor=' + factor_name + '\n')
                        str_list_4.append('  \n')

                    # Add an "intensity" sub-track for each NAME_ppois.bw file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.comparison_ppois_bw)):
                        # Common track settings
                        str_list_5.append('  track ' + prefix_short + '_ppois\n')
                        str_list_5.append('  ' + self.ucsc_hub_bigwig_info_signal_range(
                            file_path=file_path_peak_calling.comparison_ppois_bwi))
                        str_list_5.append('  shortLabel ' + prefix_short + '_ppois\n')
                        str_list_5.append('  longLabel ChIP intensity as Poisson p-value ' + prefix_long + '\n')
                        str_list_5.append('  bigDataUrl ' + file_path_peak_calling.comparison_ppois_bw + '\n')
                        # str_list_5.append('  html ...\n')
                        str_list_5.append('  visibility full\n')
                        # Common optional track settings
                        str_list_5.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigWig - Signal graphing track settings
                        str_list_5.append('  alwaysZero off\n')
                        str_list_5.append('  autoScale off\n')
                        str_list_5.append('  graphTypeDefault bar\n')
                        str_list_5.append('  maxHeightPixels 100:60:20\n')
                        # str_list_5.append('  maxWindowToQuery 10000000\n')
                        str_list_5.append('  smoothingWindow 5\n')
                        # str_list_5.append('  transformFunc NONE\n')
                        str_list_5.append('  viewLimits 0:15\n')
                        str_list_5.append('  viewLimitsMax 0:40\n')
                        str_list_5.append('  windowingFunction maximum\n')
                        # str_list_5.append('  yLineMark <#>\n')
                        # str_list_5.append('  yLineOnOff on \n')
                        # str_list_5.append('  gridDefault on\n')
                        # Composite track settings
                        str_list_5.append('  parent intensity on\n')
                        str_list_5.append('  subGroups' +
                                          ' comparison=' + comparison_name +
                                          ' factor=' + factor_name +
                                          ' scale=ppois\n')
                        str_list_5.append('  \n')

                    # Add an "intensity" sub-track for each NAME_subtract.bw file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.comparison_subtract_bw)):
                        # Common track settings
                        str_list_5.append('  track ' + prefix_short + '_subtract\n')
                        str_list_5.append('  ' + self.ucsc_hub_bigwig_info_signal_range(
                            file_path=file_path_peak_calling.comparison_subtract_bwi))
                        str_list_5.append('  shortLabel ' + prefix_short + '_subtract\n')
                        str_list_5.append('  longLabel ChIP intensity as subtraction ' + prefix_long + '\n')
                        str_list_5.append('  bigDataUrl ' + file_path_peak_calling.comparison_subtract_bw + '\n')
                        # str_list_5.append('  html ...\n')
                        str_list_5.append('  visibility full\n')
                        # Common optional track settings
                        str_list_5.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigWig - Signal graphing track settings
                        str_list_5.append('  alwaysZero off\n')
                        str_list_5.append('  autoScale off\n')
                        str_list_5.append('  graphTypeDefault bar\n')
                        str_list_5.append('  maxHeightPixels 100:60:20\n')
                        # str_list_5.append('  maxWindowToQuery 10000000\n')
                        str_list_5.append('  smoothingWindow 5\n')
                        # str_list_5.append('  transformFunc NONE\n')
                        str_list_5.append('  viewLimits 0:15\n')
                        str_list_5.append('  viewLimitsMax 0:40\n')
                        str_list_5.append('  windowingFunction maximum\n')
                        # str_list_5.append('  yLineMark <#>\n')
                        # str_list_5.append('  yLineOnOff on \n')
                        # str_list_5.append('  gridDefault on\n')
                        # Composite track settings
                        str_list_5.append('  parent intensity on\n')
                        str_list_5.append('  subGroups' +
                                          ' comparison=' + comparison_name +
                                          ' factor=' + factor_name +
                                          ' scale=subtract\n')
                        str_list_5.append('  \n')

                    # Add an "intensity" sub-track for each NAME_logFE.bw file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.comparison_log_fe_bw)):
                        # Common track settings
                        str_list_5.append('  track ' + prefix_short + '_log_fe\n')
                        str_list_5.append('  ' + self.ucsc_hub_bigwig_info_signal_range(
                            file_path=file_path_peak_calling.comparison_log_fe_bwi))
                        str_list_5.append('  shortLabel ' + prefix_short + '_log_fe\n')
                        str_list_5.append('  longLabel ChIP intensity as log10(fold-enrichment) ' + prefix_long + '\n')
                        str_list_5.append('  bigDataUrl ' + file_path_peak_calling.comparison_log_fe_bw + '\n')
                        # str_list_5.append('  html ...\n')
                        str_list_5.append('  visibility full\n')
                        # Common optional track settings
                        str_list_5.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigWig - Signal graphing track settings
                        str_list_5.append('  alwaysZero off\n')
                        str_list_5.append('  autoScale off\n')
                        str_list_5.append('  graphTypeDefault bar\n')
                        str_list_5.append('  maxHeightPixels 100:60:20\n')
                        # str_list_5.append('  maxWindowToQuery 10000000\n')
                        str_list_5.append('  smoothingWindow 5\n')
                        # str_list_5.append('  transformFunc NONE\n')
                        str_list_5.append('  viewLimits 0:15\n')
                        str_list_5.append('  viewLimitsMax 0:40\n')
                        str_list_5.append('  windowingFunction maximum\n')
                        # str_list_5.append('  yLineMark <#>\n')
                        # str_list_5.append('  yLineOnOff on \n')
                        # str_list_5.append('  gridDefault on\n')
                        # Composite track settings
                        str_list_5.append('  parent intensity on\n')
                        str_list_5.append('  subGroups' +
                                          ' comparison=' + comparison_name +
                                          ' factor=' + factor_name +
                                          ' scale=logfe\n')
                        str_list_5.append('  \n')

                    # Add a "peaks" sub-track for each NAME_peaks.bb file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.narrow_peaks_bb)):
                        # Common track settings
                        str_list_6.append('  track ' + prefix_short + '_narrow_peaks\n')
                        str_list_6.append('  type bigBed 6+4\n')
                        str_list_6.append('  shortLabel ' + prefix_short + '_narrow_peaks\n')
                        str_list_6.append('  longLabel ChIP narrow peaks ' + prefix_long + '\n')
                        str_list_6.append('  bigDataUrl ' + file_path_peak_calling.narrow_peaks_bb + '\n')
                        # str_list_6.append('  html ...\n')
                        str_list_6.append('  visibility squish\n')
                        # Common optional track settings
                        str_list_6.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigBed - Item or region track settings
                        # Composite track settings
                        str_list_6.append('  parent narrow_peaks on\n')
                        str_list_6.append('  subGroups comparison=' + comparison_name + ' factor=' + factor_name + '\n')
                        str_list_6.append('  \n')

                    # Add a "summits" sub-track for each NAME_summits.bb file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.summits_bb)):
                        # Common track settings
                        str_list_7.append('  track ' + prefix_short + '_summits\n')
                        str_list_7.append('  type bigBed 4+1\n')
                        str_list_7.append('  shortLabel ' + prefix_short + '_summits\n')
                        str_list_7.append('  longLabel ChIP summits ' + prefix_long + '\n')
                        str_list_7.append('  bigDataUrl ' + file_path_peak_calling.summits_bb + '\n')
                        # str_list_7.append('  html ...\n')
                        str_list_7.append('  visibility squish\n')
                        # Common optional track settings
                        str_list_7.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigBed - Item or region track settings
                        # Composite track settings
                        str_list_7.append('  parent summits on\n')
                        str_list_7.append('  subGroups comparison=' + comparison_name + ' factor=' + factor_name + '\n')
                        str_list_7.append('  \n')

                    # Add a "broad_peaks" sub-track for each NAME_broad_peaks.bd file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.broad_peaks_bb)):
                        # Common track settings
                        str_list_8.append('  track ' + prefix_short + '_broad_peaks\n')
                        str_list_8.append('  type bigBed 6+3\n')
                        str_list_8.append('  shortLabel ' + prefix_short + '_broad_peaks\n')
                        str_list_8.append('  longLabel ChIP broad peaks ' + prefix_long + '\n')
                        str_list_8.append('  bigDataUrl ' + file_path_peak_calling.broad_peaks_bb + '\n')
                        # str_list_8.append('  html ...\n')
                        str_list_8.append('  visibility squish\n')
                        # Common optional track settings
                        str_list_8.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigBed - Item or region track settings
                        # Composite track settings
                        str_list_8.append('  parent broad_peaks on\n')
                        str_list_8.append('  subGroups comparison=' + comparison_name + ' factor=' + factor_name + '\n')
                        str_list_8.append('  \n')

                    # Add a "gapped_peaks" sub-track for each NAME_gapped_peaks.bd file.

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_peak_calling.gapped_peaks_bb)):
                        # Common track settings
                        str_list_8.append('  track ' + prefix_short + '_gapped_peaks\n')
                        str_list_8.append('  type bigBed 12+3\n')
                        str_list_8.append('  shortLabel ' + prefix_short + '_gapped_peaks\n')
                        str_list_8.append('  longLabel ChIP gapped peaks ' + prefix_long + '\n')
                        str_list_8.append('  bigDataUrl ' + file_path_peak_calling.gapped_peaks_bb + '\n')
                        # str_list_8.append('  html ...\n')
                        str_list_8.append('  visibility squish\n')
                        # Common optional track settings
                        str_list_8.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigBed - Item or region track settings
                        # Composite track settings
                        str_list_8.append('  parent gapped_peaks on\n')
                        str_list_8.append('  subGroups comparison=' + comparison_name + ' factor=' + factor_name + '\n')
                        str_list_8.append('  \n')

            # Add UCSC trackDB entries for each alignment.

            for sample in self.sample_list:
                file_path_alignment = self.get_file_path_alignment(sample_name=sample.name)

                # Add an "alignment" sub-track for each NAME.bam file.

                # Common track settings
                str_list_1.append('  track ' + sample.name + '_alignment\n')
                str_list_1.append('  type bam\n')
                str_list_1.append('  shortLabel ' + sample.name + '_alignment\n')
                str_list_1.append('  longLabel ChIP read alignment ' + sample.name + '\n')
                str_list_1.append('  bigDataUrl ' + file_path_alignment.sample_bam + '\n')
                # str_list_1.append('  html ...\n')
                str_list_1.append('  visibility hide\n')
                # Common optional track settings
                # str_list_1.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                # bam - Compressed Sequence Alignment track settings
                # Composite track settings
                str_list_1.append('  parent alignment on\n')
                str_list_1.append('  \n')

                # Add a "coverage" sub-track for each NAME.bw file.

                # Common track settings
                str_list_2.append('  track ' + sample.name + '_coverage\n')
                str_list_2.append('  ' + self.ucsc_hub_bigwig_info_signal_range(
                    file_path=file_path_alignment.coverage_bwi_txt))
                str_list_2.append('  shortLabel ' + sample.name + '_coverage\n')
                str_list_2.append('  longLabel Normalised alignment coverage ' + sample.name + '\n')
                str_list_2.append('  bigDataUrl ' + file_path_alignment.coverage_bw + '\n')
                # str_list_2.append('  html ...\n')
                str_list_2.append('  visibility full\n')
                # Common optional track settings
                # str_list_2.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                # bigWig - Signal graphing track settings
                str_list_2.append('  alwaysZero off\n')
                str_list_2.append('  autoScale off\n')
                str_list_2.append('  graphTypeDefault bar\n')
                str_list_2.append('  maxHeightPixels 100:60:20\n')
                # str_list_2.append('  maxWindowToQuery 10000000\n')
                str_list_2.append('  smoothingWindow 5\n')
                # str_list_2.append('  transformFunc NONE\n')
                str_list_2.append('  viewLimits 0:15\n')
                str_list_2.append('  viewLimitsMax 0:40\n')
                str_list_2.append('  windowingFunction maximum\n')
                # str_list_2.append('  yLineMark <#>\n')
                # str_list_2.append('  yLineOnOff on \n')
                # str_list_2.append('  gridDefault on\n')
                # Composite track settings
                str_list_2.append('  parent coverage on\n')
                str_list_2.append('  \n')

            for comparison_name in sorted(self._factor_dict):
                for factor_name in sorted(self._factor_dict[comparison_name]):
                    if self.debug > 0:
                        print('chipseq factor_name:', factor_name)

                    # No comparison for less than two items.
                    if len(self._factor_dict[comparison_name][factor_name]) < 2:
                        continue

                    prefix_diff_bind = self.get_prefix_diff_bind(
                        comparison_name=comparison_name,
                        factor_name=factor_name)

                    file_path_diff_bind = FilePathDiffBind(prefix=prefix_diff_bind)

                    prefix_short = '_'.join((comparison_name, factor_name))

                    if os.path.exists(
                            os.path.join(self.genome_directory, file_path_diff_bind.peak_set_bb)):
                        # Common track settings
                        str_list_10.append('  track ' + prefix_diff_bind + '_consensus_peaks\n')
                        str_list_10.append('  type bigBed 6\n')
                        str_list_10.append('  shortLabel ' + prefix_short + '_consensus_peaks\n')
                        str_list_10.append('  longLabel ChIP consensus peaks ' + prefix_short + '\n')
                        str_list_10.append('  bigDataUrl ' + file_path_diff_bind.peak_set_bb + '\n')
                        # str_list_10.append('  html ...\n')
                        str_list_10.append('  visibility squish\n')
                        # Common optional track settings
                        str_list_10.append('  color ' + self.get_colour(factor=factor_name) + '\n')
                        # bigBed - Item or region track settings
                        # Composite track settings
                        str_list_10.append('  parent consensus_peaks on\n')
                        str_list_10.append('  subGroups comparison=' + comparison_name +
                                           ' factor=' + factor_name + '\n')
                        str_list_10.append('  \n')

            self.ucsc_hub_to_file(
                content=(
                        str_list_1 +
                        str_list_2 +
                        str_list_3 +
                        str_list_4 +
                        str_list_5 +
                        str_list_6 +
                        str_list_7 +
                        str_list_8 +
                        str_list_9 +
                        str_list_10
                )
            )

            return

        if self._macs_version == 1:
            report_html_1()
            report_hub_1()
        elif self._macs_version == 2:
            report_html_2()
            report_hub_2()

        return
