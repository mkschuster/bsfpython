# -*- coding: utf-8 -*-
"""Variant Calling Analysis module.

A package of classes and methods supporting variant calling analyses.
"""
#
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
import errno
import os
import pickle
import sys
import warnings
from subprocess import Popen

from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.argument import Argument
from bsf.connector import Connector, ConnectorFile
from bsf.executables.vcf import RunnableStepCsqToVep
from bsf.intervals import Container, get_genome_tiles, get_interval_tiles
from bsf.ngs import Collection, Sample
from bsf.procedure import FilePath, Runnable, ConsecutiveRunnable
from bsf.process import Command, Executable, \
    RunnableStep, RunnableStepJava, RunnableStepMove, RunnableStepLink, RunnableStepPicard
from bsf.standards import Configuration, StandardFilePath, EnsemblVEP, JavaArchive


class RunnableStepGATK(RunnableStepJava):
    """C{bsf.analyses.variant_calling.RunnableStepGATK} class representing a Genome Analysis Toolkit (GATK) program.

    Attributes:
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_minimum=None,
            java_heap_maximum=None,
            java_jar_path=None):
        """Initialise a C{bsf.analyses.variant_calling.RunnableStepGATK}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} into the C{bsf.analysis.Stage}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | None
        @param java_heap_minimum: Java heap minimum size (-Xms option)
        @type java_heap_minimum: str | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java Archive (JAR) file path
        @type java_jar_path: str | None
        """
        super(RunnableStepGATK, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path,
            java_heap_minimum=java_heap_minimum,
            java_heap_maximum=java_heap_maximum,
            java_jar_path=java_jar_path)

        # The GATK algorithm is then another empty sub-command.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command()

        return

    def add_gatk_option(self, key, value, override=False):
        """Add a C{bsf.argument.OptionPair} to a C{bsf.analyses.variant_calling.RunnableStepGATK}.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """

        return self.sub_command.sub_command.add_option_long(key=key, value=value, override=override)

    def add_gatk_switch(self, key):
        """Add a C{bsf.argument.SwitchLong} to a C{bsf.analyses.variant_calling.RunnableStepGATK}.

        @param key: Option key
        @type key: str
        """

        return self.sub_command.sub_command.add_switch_long(key=key)


class FilePathAlignment(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathAlignment} models files in a sample-specific directory.

    Attributes:
    @ivar aligned_bam: Alignment BAM file path
    @type aligned_bam: str
    @ivar aligned_bai: Alignment BAI file path
    @type aligned_bai: str
    @ivar aligned_md5: Alignment MD5 check sum file path
    @type aligned_md5: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAlignment} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathAlignment, self).__init__(prefix=prefix)

        self.aligned_bam = prefix + '.bam'
        self.aligned_bai = prefix + '.bai'
        self.aligned_md5 = prefix + '.bam.md5'

        return


class FilePathProcessReadGroup(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathProcessReadGroup} models files in a sample-specific directory.

    Attributes:
    @ivar trimmed_bam: Fulcrum Genomics (fgbio) TrimPrimer BAM file path
    @type trimmed_bam: str
    @ivar trimmed_bai: Fulcrum Genomics (fgbio) TrimPrimer BAI file path
    @type trimmed_bai: str
    @ivar trimmed_md5: Fulcrum Genomics (fgbio) TrimPrimer MD5 check sum file path
    @type trimmed_md5: str
    @ivar duplicates_marked_bam: Picard Mark Duplicates BAM file path
    @type duplicates_marked_bam: str
    @ivar duplicates_marked_bai: Picard Mark Duplicates BAI file path
    @type duplicates_marked_bai: str
    @ivar duplicates_marked_md5: Picard Mark Duplicates MD5 check sum file path
    @type duplicates_marked_md5: str
    @ivar duplicate_metrics: Picard Mark Duplicates Metrics TSV file path
    @type duplicate_metrics: str
    @ivar realigner_targets: GATK RealignerTargetCreator interval list file path
    @type realigner_targets: str
    @ivar realigned_bam: GATK re-aligned BAM file path
    @type realigned_bam: str
    @ivar realigned_bai: GATK re-aligned BAI file path
    @type realigned_bai: str
    @ivar realigned_md5: GATK re-aligned MD5 check sum file path
    @type realigned_md5: str
    @ivar recalibration_table_pre: GATK pre-recalibration table file path
    @type recalibration_table_pre: str
    @ivar recalibration_table_post: GATK post-recalibration table file path
    @type recalibration_table_post: str
    @ivar recalibration_plot: GATK Recalibration plot PDF file path
    @type recalibration_plot: str
    @ivar recalibrated_bam: Recalibrated BAM file path
    @type recalibrated_bam: str
    @ivar recalibrated_bai: Recalibrated BAI file path
    @type recalibrated_bai: str
    @ivar recalibrated_md5: Recalibrated BAM MD5 check sum file path
    @type recalibrated_md5: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathProcessReadGroup} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathProcessReadGroup, self).__init__(prefix=prefix)

        self.trimmed_bam = prefix + '_trimmed.bam'
        self.trimmed_bai = prefix + '_trimmed.bai'
        self.trimmed_md5 = prefix + '_trimmed.bam.md5'
        self.duplicates_marked_bam = prefix + '_duplicates_marked.bam'
        self.duplicates_marked_bai = prefix + '_duplicates_marked.bai'
        self.duplicates_marked_md5 = prefix + '_duplicates_marked.bam.md5'
        self.duplicate_metrics = prefix + '_duplicate_metrics.tsv'
        self.realigner_targets = prefix + '_realigner.interval_list'
        self.realigned_bam = prefix + '_realigned.bam'
        self.realigned_bai = prefix + '_realigned.bai'
        self.realigned_md5 = prefix + '_realigned.bam.md5'
        self.recalibration_table_pre = prefix + '_recalibration_pre.table'
        self.recalibration_table_post = prefix + '_recalibration_post.table'
        self.recalibration_plot = prefix + '_recalibration_report.pdf'
        self.recalibrated_bam = prefix + '_recalibrated.bam'
        self.recalibrated_bai = prefix + '_recalibrated.bai'
        self.recalibrated_md5 = prefix + '_recalibrated.bam.md5'

        return


class FilePathProcessSample(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathProcessSample} models files in a sample-specific directory.

    Attributes:
    @ivar merged_bam: Merged BAM file path
    @type merged_bam: str
    @ivar merged_bai: Merged BAI file path
    @type merged_bai: str
    @ivar merged_md5: Merged MD5 check sum file path
    @type merged_md5: str
    @ivar duplicates_marked_bam: Duplicates-marked BAM file path
    @type duplicates_marked_bam: str
    @ivar duplicates_marked_bai: Duplicates-marked BAI file path
    @type duplicates_marked_bai: str
    @ivar duplicates_marked_md5: Duplicates-marked MD5 check sum file path
    @type duplicates_marked_md5: str
    @ivar duplicate_metrics: Duplicate metrics TSV file path
    @type duplicate_metrics: str
    @ivar realigner_targets: Re-aligner target file path
    @type realigner_targets: str
    @ivar realigned_bam: Re-aligned BAM file path
    @type realigned_bam: str
    @ivar realigned_bai: Re-aligned BAI file path
    @type realigned_bai: str
    @ivar realigned_md5: Re-aligned MD5 check sum file path
    @type realigned_md5: str
    @ivar realigned_bam_bai: Re-aligned BAM BAI file path
    @type realigned_bam_bai: str
    @ivar alignment_summary_metrics: Alignment summary metrics TSV file path
    @type alignment_summary_metrics: str
    @ivar raw_variants_gvcf_vcf: Raw variants gVCF file path
    @type raw_variants_gvcf_vcf: str
    @ivar raw_variants_gvcf_tbi: Raw variants TBI file path
    @type raw_variants_gvcf_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathProcessSample} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathProcessSample, self).__init__(prefix=prefix)

        self.merged_bam = prefix + '_merged.bam'
        self.merged_bai = prefix + '_merged.bai'
        self.merged_md5 = prefix + '_merged.bam.md5'
        self.duplicates_marked_bam = prefix + '_duplicates_marked.bam'
        self.duplicates_marked_bai = prefix + '_duplicates_marked.bai'
        self.duplicates_marked_md5 = prefix + '_duplicates_marked.bam.md5'
        self.duplicate_metrics = prefix + '_duplicate_metrics.tsv'
        self.realigner_targets = prefix + '_realigner.intervals'
        self.realigned_bam = prefix + '_realigned.bam'
        self.realigned_bai = prefix + '_realigned.bai'
        self.realigned_md5 = prefix + '_realigned.bam.md5'
        self.realigned_bam_bai = prefix + '_realigned.bam.bai'
        self.alignment_summary_metrics = prefix + '_alignment_summary_metrics.tsv'
        self.raw_variants_gvcf_vcf = prefix + '_raw_variants.g.vcf.gz'
        self.raw_variants_gvcf_tbi = prefix + '_raw_variants.g.vcf.gz.tbi'

        return


class FilePathDiagnoseSample(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathDiagnoseSample} models files in a sample-specific directory.

    Attributes:
    @ivar diagnose_targets_vcf: Diagnose Targets VCF file path
    @type diagnose_targets_vcf: str
    @ivar diagnose_targets_tbi: Diagnose Targets TBI file path
    @type diagnose_targets_tbi: str
    @ivar missing_intervals: Missing intervals file path
    @type missing_intervals: str
    @ivar missing_report: Missing GATK Report file path
    @type missing_report: str
    @ivar callable_bed: Callable BED file path
    @type callable_bed: str
    @ivar callable_txt: Callable text file path
    @type callable_txt: str
    @ivar callable_bb: Callable BigBed file path
    @type callable_bb: str
    @ivar sorted_bed: Sorted BED file path
    @type sorted_bed: str
    @ivar non_callable_loci_tsv: Non-callable Loci TSV file path
    @type non_callable_loci_tsv: str
    @ivar non_callable_regions_tsv: Non-callable Regions TSV file path
    @type non_callable_regions_tsv: str
    @ivar non_callable_summary_tsv: Non-callable Summary TSV file path
    @type non_callable_summary_tsv: str
    @ivar hybrid_selection_metrics: Hybrid Selection Metrics TSV file path
    @type hybrid_selection_metrics: str
    @ivar insert_size_pdf: Insert size plot PDF file path
    @type insert_size_pdf: str
    @ivar insert_size_png: Insert size plot PNG file path
    @type insert_size_png: str
    @ivar insert_size_tsv: Insert size TSV file path
    @type insert_size_tsv: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathDiagnoseSample} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathDiagnoseSample, self).__init__(prefix=prefix)

        self.diagnose_targets_vcf = prefix + '_diagnose_targets.vcf.gz'
        self.diagnose_targets_tbi = prefix + '_diagnose_targets.vcf.gz.tbi'
        self.missing_intervals = prefix + '_missing.intervals'
        self.missing_report = prefix + '_missing.gatkreport'
        self.callable_bed = prefix + '_callable_loci.bed'
        self.callable_txt = prefix + '_callable_loci.txt'
        self.callable_bb = prefix + '_callable_loci.bb'
        self.sorted_bed = prefix + '_callable_sorted.bed'
        # Defined in bsf_variant_calling_coverage.R
        self.non_callable_loci_tsv = prefix + '_non_callable_loci.tsv'
        self.non_callable_regions_tsv = prefix + '_non_callable_regions.tsv'
        self.non_callable_summary_tsv = prefix + '_non_callable_summary.tsv'
        self.hybrid_selection_metrics = prefix + '_hybrid_selection_metrics.tsv'
        # Defined in bsf_variant_calling_insert_size.R
        self.insert_size_pdf = prefix + '_insert_size.pdf'
        self.insert_size_png = prefix + '_insert_size.png'
        self.insert_size_tsv = prefix + '_insert_size.tsv'

        return


class FilePathMergeCohort(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathMergeCohort} models files in a cohort-specific directory.

    Attributes:
    @ivar combined_gvcf_vcf: Combined gVCF file path
    @type combined_gvcf_vcf: str
    @ivar combined_gvcf_tbi: Combined TBI file path
    @type combined_gvcf_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathMergeCohort} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathMergeCohort, self).__init__(prefix=prefix)

        self.combined_gvcf_vcf = prefix + '_combined.g.vcf.gz'
        self.combined_gvcf_tbi = prefix + '_combined.g.vcf.gz.tbi'

        return


class FilePathGenotypeCohort(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathGenotypeCohort} models files in a cohort-specific directory.

    Attributes:
    @ivar genotyped_raw_vcf: Genotyped raw VCF file path
    @type genotyped_raw_vcf: str
    @ivar genotyped_raw_tbi: Genotyped raw TBI file path
    @type genotyped_raw_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathGenotypeCohort} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathGenotypeCohort, self).__init__(prefix=prefix)

        self.genotyped_raw_vcf = prefix + '_genotyped_raw_snp_raw_indel.vcf.gz'
        self.genotyped_raw_tbi = prefix + '_genotyped_raw_snp_raw_indel.vcf.gz.tbi'

        return


class FilePathProcessCohort(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathGenotypeCohort} models files in a cohort-specific directory.

    Attributes:
    @ivar recalibrated_snp_raw_indel_vcf: Recalibrated SNP raw InDel VCF file path
    @type recalibrated_snp_raw_indel_vcf: str
    @ivar recalibrated_snp_raw_indel_tbi: Recalibrated SNP raw InDel TBI file path
    @type recalibrated_snp_raw_indel_tbi: str
    @ivar recalibrated_snp_recalibrated_indel_vcf: Recalibrated SNP recalibrated InDel VCF file path
    @type recalibrated_snp_recalibrated_indel_vcf: str
    @ivar recalibrated_snp_recalibrated_indel_tbi: Recalibrated SNP recalibrated InDel TBI file path
    @type recalibrated_snp_recalibrated_indel_tbi: str
    @ivar multi_sample_vcf: Multi-sample VCF file path
    @type multi_sample_vcf: str
    @ivar multi_sample_tbi: Multi-sample TBI file path
    @type multi_sample_tbi: str
    @ivar recalibration_indel: Recalibration InDel file path
    @type recalibration_indel: str
    @ivar recalibration_snp: Recalibration SNP file path
    @type recalibration_snp: str
    @ivar tranches_indel: Tranches InDel file path
    @type tranches_indel: str
    @ivar tranches_snp: Tranches SNP file path
    @type tranches_snp: str
    @ivar plots_indel: Plots InDel R script file path
    @type plots_indel: str
    @ivar plots_snp: Plots SNP R script file path
    @type plots_snp: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathProcessCohort} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathProcessCohort, self).__init__(prefix=prefix)

        self.recalibrated_snp_raw_indel_vcf = prefix + '_recalibrated_snp_raw_indel.vcf.gz'
        self.recalibrated_snp_raw_indel_tbi = prefix + '_recalibrated_snp_raw_indel.vcf.gz.tbi'
        self.recalibrated_snp_recalibrated_indel_vcf = prefix + '_recalibrated_snp_recalibrated_indel.vcf.gz'
        self.recalibrated_snp_recalibrated_indel_tbi = prefix + '_recalibrated_snp_recalibrated_indel.vcf.gz.tbi'
        self.multi_sample_vcf = prefix + '_multi_sample.vcf.gz'
        self.multi_sample_tbi = prefix + '_multi_sample.vcf.gz.tbi'
        self.recalibration_indel = prefix + '_recalibration_indel.recal'
        self.recalibration_snp = prefix + '_recalibration_snp.recal'
        self.tranches_indel = prefix + '_recalibration_indel.tranches'
        self.tranches_snp = prefix + '_recalibration_snp.tranches'
        self.plots_indel = prefix + '_recalibration_indel.R'
        self.plots_snp = prefix + '_recalibration_snp.R'

        return


class FilePathAnnotateSnpEff(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathAnnotateSnpEff} models snpEff-annotated, cohort-specific files.

    Attributes:
    @ivar complete_vcf: Complete snpEff VCF file path
    @type complete_vcf: str
    @ivar complete_vcf_bgz: Bgzip-compressed complete snpEff VCF file path
    @type complete_vcf_bgz: str
    @ivar complete_vcf_tbi: Tabix-indexed complete snpEff VCF file path
    @type complete_vcf_tbi: str
    @ivar complete_genes: Complete snpEff genes annotation TXT file path
    @type complete_genes: str
    @ivar complete_stats: Complete snpEff statistics HTML file path
    @type complete_stats: str
    @ivar gatk_vcf: GATK-style snpEff VCF file path
    @type gatk_vcf: str
    @ivar gatk_vcf_bgz: Bgzip-compressed GATK-style snpEff VCF file path
    @type gatk_vcf_bgz: str
    @ivar gatk_vcf_tbi: Tabix-indexed GATK-style snpEff VCF file path
    @type gatk_vcf_tbi: str
    @ivar gatk_genes: GATK-style snpEff genes annotation TXT file path
    @type gatk_genes: str
    @ivar gatk_stats: GATK-style snpEff statistics HTML file path
    @type gatk_stats: str
    @ivar annotated_vcf: Annotated VCF file path
    @type annotated_vcf: str
    @ivar annotated_tbi: Annotated TBI file path
    @type annotated_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAnnotateSnpEff} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathAnnotateSnpEff, self).__init__(prefix=prefix)

        self.complete_vcf = prefix + '_complete.vcf'
        self.complete_vcf_bgz = prefix + '_complete.vcf.gz'
        self.complete_vcf_tbi = prefix + '_complete.vcf.gz.tbi'
        self.complete_genes = prefix + '_complete_summary.genes.txt'
        self.complete_stats = prefix + '_complete_summary.html'
        self.gatk_vcf = prefix + '_gatk.vcf'
        self.gatk_vcf_bgz = prefix + '_gatk.vcf.gz'
        self.gatk_vcf_tbi = prefix + '_gatk.vcf.gz.tbi'
        self.gatk_genes = prefix + '_gatk_summary.genes.txt'
        self.gatk_stats = prefix + '_gatk_summary.html'
        self.annotated_vcf = prefix + '_annotated.vcf.gz'
        self.annotated_tbi = prefix + '_annotated.vcf.gz.tbi'

        return


class FilePathAnnotateVEP(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathAnnotateVEP} models Ensembl VEP-annotated, cohort-specific files.

    Attributes:
    @ivar statistics: Ensembl VEP statistics HTML file path
    @type statistics: str
    @ivar complete_raw_vcf: Ensembl VEP complete raw VCF file path
    @type complete_raw_vcf: str
    @ivar complete_raw_vcf_bgz: Ensembl VEP complete, raw Bgzip-compressed VCF file path
    @type complete_raw_vcf_bgz: str
    @ivar complete_raw_vcf_tbi: Ensembl VEP complete, raw Tabix-indexed TBI file path
    @type complete_raw_vcf_tbi: str
    @ivar filtered_raw_vcf: Ensembl VEP filtered, raw VCF file path
    @type filtered_raw_vcf: str
    @ivar filtered_raw_vcf_bgz: Ensembl VEP filtered, raw Bgzip-compressed VCF file path
    @type filtered_raw_vcf_bgz: str
    @ivar filtered_raw_vcf_tbi: Ensembl VEP filtered, raw Tabix-indexed TBI file path
    @type filtered_raw_vcf_tbi: str
    @ivar complete_vcf_bgz: Ensembl VEP complete VCF file path
    @type complete_vcf_bgz: str
    @ivar complete_vcf_tbi: Ensembl VEP complete TBI file path
    @type complete_vcf_tbi: str
    @ivar filtered_vcf_bgz: Ensembl VEP filtered VCF file path
    @type filtered_vcf_bgz: str
    @ivar filtered_vcf_tbi: Ensembl VEP filtered TBI file path
    @type filtered_vcf_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAnnotateVEP} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathAnnotateVEP, self).__init__(prefix=prefix)

        self.statistics = prefix + '_statistics.html'
        # Complete VEP set raw
        self.complete_raw_vcf = prefix + '_complete_raw.vcf'
        self.complete_raw_vcf_bgz = prefix + '_complete_raw.vcf.gz'
        self.complete_raw_vcf_tbi = prefix + '_complete_raw.vcf.gz.tbi'
        # Filtered VEP set raw
        self.filtered_raw_vcf = prefix + '_filtered_raw.vcf'
        self.filtered_raw_vcf_bgz = prefix + '_filtered_raw.vcf.gz'
        self.filtered_raw_vcf_tbi = prefix + '_filtered_raw.vcf.gz.tbi'
        # Complete VEP set VCF.Filter-converted
        self.complete_vcf_bgz = prefix + '_complete.vcf.gz'
        self.complete_vcf_tbi = prefix + '_complete.vcf.gz.tbi'
        # Filtered VEP set VCF.Filter-converted
        self.filtered_vcf_bgz = prefix + '_filtered.vcf.gz'
        self.filtered_vcf_tbi = prefix + '_filtered.vcf.gz.tbi'

        return


class FilePathSplitCohort(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathSplitCohort} models sample-specific files.

    Attributes:
    @ivar sample_vcf: Sample VCF file path
    @type sample_vcf: str
    @ivar sample_tbi: Sample TBI file path
    @type sample_tbi: str
    @ivar sample_tsv: Sample TSV file path
    @type sample_tsv: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSplitCohort} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathSplitCohort, self).__init__(prefix=prefix)

        self.sample_vcf = prefix + '.vcf.gz'
        self.sample_tbi = prefix + '.vcf.gz.tbi'
        self.sample_tsv = prefix + '.tsv'

        return


class FilePathSummary(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSummary} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathSummary, self).__init__(prefix=prefix)

        self.alignment_metrics_sample_tsv = prefix + '_alignment_metrics_sample.tsv'
        self.alignment_metrics_read_group_tsv = prefix + '_alignment_metrics_read_group.tsv'
        self.alignment_absolute_sample_pdf = prefix + '_alignment_absolute_sample.pdf'
        self.alignment_absolute_sample_png = prefix + '_alignment_absolute_sample.png'
        self.alignment_absolute_read_group_pdf = prefix + '_alignment_absolute_read_group.pdf'
        self.alignment_absolute_read_group_png = prefix + '_alignment_absolute_read_group.png'
        self.alignment_percentage_sample_pdf = prefix + '_alignment_percentage_sample.pdf'
        self.alignment_percentage_sample_png = prefix + '_alignment_percentage_sample.png'
        self.alignment_percentage_read_group_pdf = prefix + '_alignment_percentage_read_group.pdf'
        self.alignment_percentage_read_group_png = prefix + '_alignment_percentage_read_group.png'
        self.duplication_metrics_sample_tsv = prefix + '_duplication_metrics_sample.tsv'
        self.duplication_percentage_sample_png = prefix + '_duplication_percentage_sample.png'
        self.duplication_percentage_sample_pdf = prefix + '_duplication_percentage_sample.pdf'
        self.duplication_levels_sample_pdf = prefix + '_duplication_levels_sample.pdf'
        self.duplication_levels_sample_png = prefix + '_duplication_levels_sample.png'
        self.hybrid_metrics_sample_tsv = prefix + '_hybrid_metrics_sample.tsv'
        self.hybrid_metrics_read_group_tsv = prefix + '_hybrid_metrics_read_group.tsv'
        self.hybrid_unique_percentage_sample_pdf = prefix + '_hybrid_unique_percentage_sample.pdf'
        self.hybrid_unique_percentage_sample_png = prefix + '_hybrid_unique_percentage_sample.png'
        self.hybrid_unique_percentage_read_group_pdf = prefix + '_hybrid_unique_percentage_read_group.pdf'
        self.hybrid_unique_percentage_read_group_png = prefix + '_hybrid_unique_percentage_read_group.png'
        self.hybrid_coverage_sample_pdf = prefix + '_hybrid_target_coverage_sample.pdf'
        self.hybrid_coverage_sample_png = prefix + '_hybrid_target_coverage_sample.png'
        self.hybrid_coverage_read_group_pdf = prefix + '_hybrid_target_coverage_read_group.pdf'
        self.hybrid_coverage_read_group_png = prefix + '_hybrid_target_coverage_read_group.png'
        self.hybrid_coverage_levels_sample_pdf = prefix + '_hybrid_target_coverage_levels_sample.pdf'
        self.hybrid_coverage_levels_sample_png = prefix + '_hybrid_target_coverage_levels_sample.png'
        self.hybrid_coverage_levels_read_group_pdf = prefix + '_hybrid_target_coverage_levels_read_group.pdf'
        self.hybrid_coverage_levels_read_group_png = prefix + '_hybrid_target_coverage_levels_read_group.png'
        self.hybrid_excluded_bases_read_group_pdf = prefix + '_hybrid_excluded_bases_read_group.pdf'
        self.hybrid_excluded_bases_read_group_png = prefix + '_hybrid_excluded_bases_read_group.png'
        self.hybrid_excluded_bases_sample_pdf = prefix + '_hybrid_excluded_bases_sample.pdf'
        self.hybrid_excluded_bases_sample_png = prefix + '_hybrid_excluded_bases_sample.png'
        self.non_callable_metrics_sample_tsv = prefix + '_non_callable_metrics_sample.tsv'
        self.non_callable_absolute_sample_pdf = prefix + '_non_callable_absolute_sample.pdf'
        self.non_callable_absolute_sample_png = prefix + '_non_callable_absolute_sample.png'
        self.non_callable_percentage_sample_pdf = prefix + '_non_callable_percentage_sample.pdf'
        self.non_callable_percentage_sample_png = prefix + '_non_callable_percentage_sample.png'

        return


class FilePathSomatic(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathSomatic} models files in somatic scatter and gathering.

    Attributes:
    @ivar somatic_vcf: Somatic VCF file path
    @type somatic_vcf: str
    @ivar somatic_tbi: Somatic TBI file path
    @type somatic_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSomatic} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathSomatic, self).__init__(prefix=prefix)

        self.somatic_vcf = prefix + '_somatic.vcf.gz'
        self.somatic_tbi = prefix + '_somatic.vcf.gz.tbi'

        return


class FilePathSomaticScatterGather(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathSomaticScatterGather} models files in somatic scatter and gathering.

    Attributes:
    @ivar somatic_vcf: Somatic VSF file path
    @type somatic_vcf: str
    @ivar somatic_tbi: Somatic TBI file path
    @type somatic_tbi: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSomaticScatterGather} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathSomaticScatterGather, self).__init__(prefix=prefix)

        self.somatic_vcf = prefix + '_somatic.vcf.gz'
        self.somatic_tbi = prefix + '_somatic.vcf.gz.tbi'

        return


class FilePathSplitSomatic(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathSplitSomatic} models sample-specific somatic variant calling files.

    Attributes:
    @ivar comparison_tsv: Comparison-specific TSV file path
    @type comparison_tsv: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSplitSomatic} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathSplitSomatic, self).__init__(prefix=prefix)

        self.comparison_tsv = prefix + '.tsv'

        return


class VariantCallingGATKComparison(object):
    """C{bsf.analyses.variant_calling.VariantCallingGATKComparison} class representing a somatic comparison.

    Attributes:
    @ivar normal_sample: Normal sample
    @type normal_sample: Sample | None
    @ivar tumor_sample: Tumour sample
    @type tumor_sample: Sample | None
    @ivar panel_of_normal_path: File path to a Panel-Of-Normal (PON) VCF file.
    @type panel_of_normal_path: str | None
    """

    def __init__(
            self,
            normal_sample=None,
            tumor_sample=None,
            panel_of_normal_path=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATKComparison} object.

        @param normal_sample: Normal sample
        @type normal_sample: Sample | None
        @param tumor_sample: Tumour sample
        @type tumor_sample: Sample | None
        @param panel_of_normal_path: File path to a Panel-Of-Normal (PON) VCF file.
        @type panel_of_normal_path: str | None
        """

        self.normal_sample = normal_sample
        self.tumor_sample = tumor_sample
        self.panel_of_normal_path = panel_of_normal_path

        return

    @property
    def get_name(self):
        """Get the name of a C{bsf.analyses.variant_calling.VariantCallingGATKComparison}.

        @return: Comparison name
        @rtype: str | None
        """
        if self.normal_sample and self.tumor_sample:
            return '__'.join((self.normal_sample.name, self.tumor_sample.name))
        # elif self.normal_sample:
        #     return self.normal_sample.name
        elif self.tumor_sample:
            return self.tumor_sample.name
        else:
            return


class VariantCallingGATKAmplicons(object):
    """C{bsf.analyses.variant_calling.VariantCallingGATKAmplicons} class representing amplicons.

    Attributes:
    @ivar name: Name
    @type name: str | None
    @ivar amplicons_path: Amplicons file path
    @type amplicons_path: str | None
    """

    @classmethod
    def from_sample(cls, sample):
        """Create a C{VariantCallingGATKAmplicons} object from a C{bsf.ngs.Sample} object.

        @param sample: C{bsf.ngs.Sample}
        @type sample: Sample
        @return: C{bsf.analyses.variant_calling.VariantCallingGATKAmplicons}
        @rtype: VariantCallingGATKAmplicons
        """
        amplicons = cls()

        if 'Amplicons Name' in sample.annotation_dict:
            amplicons_name_list = sample.annotation_dict['Amplicons Name']
            if len(amplicons_name_list) > 1:
                raise Exception('More than one Amplicons Name annotation per sample is not allowed.\n'
                                'Sample: {!r} Amplicons Name list: {!r}'.
                                format(sample.name, amplicons_name_list))
            amplicons.name = amplicons_name_list[0]

        if 'Amplicons' in sample.annotation_dict:
            amplicons_path_list = sample.annotation_dict['Amplicons']
            if len(amplicons_path_list) > 1:
                raise Exception('More than one Amplicons annotation per sample is not allowed.\n'
                                'Sample: {!r} Amplicons list: {!r}'.
                                format(sample.name, amplicons_path_list))
            amplicons.amplicons_path = amplicons_path_list[0]
            if amplicons.amplicons_path and not os.path.isabs(amplicons.amplicons_path):
                amplicons.calling_path = Configuration.get_absolute_path(
                    file_path=amplicons.calling_path,
                    default_path=StandardFilePath.get_resource_intervals(absolute=True))

        return amplicons

    def __init__(self, name=None, amplicons_path=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATKAmplicons} object.

        @param name: Name
        @type name: str | None
        @param amplicons_path: Calling intervals file path
        @type amplicons_path: str | None
        """
        self.name = name
        self.amplicons_path = amplicons_path

        return


class VariantCallingGATKCallingIntervals(object):
    """C{bsf.analyses.variant_calling.VariantCallingGATKCallingIntervals} class representing calling intervals.

    Attributes:
    @ivar name: Name
    @type name: str | None
    @ivar calling_path: Calling intervals file path
    @type calling_path: str | None
    """

    @classmethod
    def from_sample(cls, sample):
        """Create a C{VariantCallingGATKCallingIntervals} object from a C{bsf.ngs.Sample} object.

        @param sample: C{bsf.ngs.Sample}
        @type sample: Sample
        @return: C{bsf.analyses.variant_calling.VariantCallingGATKCallingIntervals}
        @rtype: VariantCallingGATKCallingIntervals
        """
        calling_intervals = cls()

        if 'Calling Name' in sample.annotation_dict:
            calling_name_list = sample.annotation_dict['Calling Name']
            if len(calling_name_list) > 1:
                raise Exception('More than one Calling Name annotation per sample is not allowed.\n'
                                'Sample: {!r} Calling Name list: {!r}'.
                                format(sample.name, calling_name_list))
            calling_intervals.name = calling_name_list[0]

        if 'Calling Intervals' in sample.annotation_dict:
            calling_interval_list = sample.annotation_dict['Calling Intervals']
            if len(calling_interval_list) > 1:
                raise Exception('More than one Calling Intervals annotation per sample is not allowed.\n'
                                'Sample: {!r} Calling Intervals list: {!r}'.
                                format(sample.name, calling_interval_list))
            calling_intervals.calling_path = calling_interval_list[0]
            if calling_intervals.calling_path and not os.path.isabs(calling_intervals.calling_path):
                calling_intervals.calling_path = Configuration.get_absolute_path(
                    file_path=calling_intervals.calling_path,
                    default_path=StandardFilePath.get_resource_intervals(absolute=True))

        return calling_intervals

    def __init__(self, name=None, calling_path=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATKCallingIntervals} object.

        @param name: Name
        @type name: str | None
        @param calling_path: Calling intervals file path
        @type calling_path: str | None
        """
        self.name = name
        self.calling_path = calling_path

        return


class VariantCallingGATKTargetIntervals(object):
    """C{bsf.analyses.variant_calling.VariantCallingGATKTargetIntervals} class representing target intervals.

    Attributes:
    @ivar name: Name
    @type name: str | None
    @ivar probes_path: Probe (bait) intervals file path
    @type probes_path: str | None
    @ivar targets_path: Target intervals file path
    @type targets_path: str | None
    """

    @classmethod
    def from_sample(cls, sample):
        """Create a C{VariantCallingGATKTargetIntervals} object from a C{bsf.ngs.Sample} object.

        @param sample: C{bsf.ngs.Sample}
        @type sample: Sample
        @return: C{bsf.analyses.variant_calling.VariantCallingGATKTargetIntervals}
        @rtype: VariantCallingGATKTargetIntervals
        """
        target_intervals = cls()

        if 'Target Name' in sample.annotation_dict:
            target_name_list = sample.annotation_dict['Target Name']
            if len(target_name_list) > 1:
                raise Exception('More than one Target Name annotation per sample is not allowed.\n'
                                'Sample: {!r} Target Name list: {!r}'.
                                format(sample.name, target_name_list))
            target_intervals.name = target_name_list[0]

        if 'Target Intervals' in sample.annotation_dict:
            target_interval_list = sample.annotation_dict['Target Intervals']
            if len(target_interval_list) > 1:
                raise Exception('More than one Target Intervals annotation per sample is not allowed.\n'
                                'Sample: {!r} Target Intervals list: {!r}'.
                                format(sample.name, target_interval_list))
            target_intervals.targets_path = target_interval_list[0]
            if target_intervals.targets_path and not os.path.isabs(target_intervals.targets_path):
                target_intervals.targets_path = Configuration.get_absolute_path(
                    file_path=target_intervals.targets_path,
                    default_path=StandardFilePath.get_resource_intervals(absolute=True))

        if 'Probe Intervals' in sample.annotation_dict:
            probe_interval_list = sample.annotation_dict['Probe Intervals']
            if len(probe_interval_list) > 1:
                raise Exception('More than one Probe Intervals annotation per sample is not allowed.\n'
                                'Sample: {!r} Probe Intervals list: {!r}'.
                                format(sample.name, probe_interval_list))
            target_intervals.probes_path = probe_interval_list[0]
            if target_intervals.probes_path and not os.path.isabs(target_intervals.probes_path):
                target_intervals.probes_path = Configuration.get_absolute_path(
                    file_path=target_intervals.probes_path,
                    default_path=StandardFilePath.get_resource_intervals(absolute=True))

        return target_intervals

    def __init__(self, name=None, probes_path=None, targets_path=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATKTargetIntervals} object.

        @param name: Name
        @type name: str | None
        @param probes_path: Probes (baits) interval file path
        @type probes_path: str | None
        @param targets_path: Targets interval file path
        @type targets_path: str | None
        """
        self.name = name
        self.probes_path = probes_path
        self.targets_path = targets_path

        return


class VariantCallingGATK(Analysis):
    """C{bsf.analyses.variant_calling.VariantCallingGATK} class representing the logic to run the
    Genome Analysis Toolkit (GATK).

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar variants_to_table: Run I{GATK VariantsToTable} to convert VCF to TSV files
    @type variants_to_table: bool
    @ivar replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects for processing or run them separately
    @type replicate_grouping: bool | None
    @ivar bwa_genome_db: Genome sequence file path with BWA index
    @type bwa_genome_db: str | None
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | None
    @ivar cohort_name: Cohort name
    @type cohort_name: str | None
    @ivar accessory_cohort_gvcfs: Python C{list} of Python C{str} (GVCF file path) objects
    @type accessory_cohort_gvcfs: list[str]
    @ivar skip_mark_duplicates: Skip the Picard MarkDuplicates step
    @type skip_mark_duplicates: bool | None
    @ivar skip_indel_realignment: Skip the GATK RealignerTargetCreator and GATK IndelRealigner steps
    @type skip_indel_realignment: bool | None
    @ivar known_sites_discovery: VCF file path for variant discovery via Haplotype Caller or Unified Genotyper
    @type known_sites_discovery: str | None
    @ivar known_sites_realignment: Python C{list} of Python C{str} VCF file paths
        for InDel realignment
    @type known_sites_realignment: list[str] | None
    @ivar known_sites_recalibration: Python C{list} of Python C{str} VCF file paths
        for Base Quality Score Recalibration (BQSR)
    @type known_sites_recalibration: list[str] | None
    @ivar known_somatic_discovery: Catalogue Of Somatic Mutations In Cancer (COSMIC) VCF file path
        for somatic variant discovery via MuTect2
    @type known_somatic_discovery: list[str] | None
    @ivar annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
        Python C{tuple} of
        Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
    @type annotation_resources_dict: dict[str, (str, list[str])] | None
    @ivar truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
    @type truth_sensitivity_filter_level_indel: str | None
    @ivar truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
    @type truth_sensitivity_filter_level_snp: str | None
    @ivar vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
    @type vqsr_skip_indel: bool | None
    @ivar vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
    @type vqsr_skip_snp: bool | None
    @ivar vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
    @type vqsr_resources_indel_dict: dict[str, dict[str, str]] | None
    @ivar vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
    @type vqsr_resources_snp_dict: dict[str, dict[str, str]] | None
    @ivar vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_indel_list: list[str] | None
    @ivar vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_snp_list: list[str] | None
    @ivar vqsr_bad_lod_cutoff_indel: LOD score cutoff for negative training set for INDELs
    @type vqsr_bad_lod_cutoff_indel: float | None
    @ivar vqsr_bad_lod_cutoff_snp: LOD score cutoff for negative training set for SNPs
    @type vqsr_bad_lod_cutoff_snp: float | None
    @ivar vqsr_max_gaussians_pos_indel: Maximum number of Gaussians in the positive training for INDELs
    @type vqsr_max_gaussians_pos_indel: int | None
    @ivar vqsr_max_gaussians_pos_snp: Maximum number of Gaussians in the positive training for SNPs
    @type vqsr_max_gaussians_pos_snp: int | None
    @ivar exclude_intervals_list: Python C{list} of Python C{str} (intervals) to exclude from the analysis
    @type exclude_intervals_list: list[str] | None
    @ivar include_intervals_list: Python C{list} of Python C{str} (intervals) to include in the analysis
    @type include_intervals_list: list[str] | None
    @ivar interval_padding: Interval padding
    @type interval_padding: int | None
    @ivar scatter_intervals_path: Picard ScatterIntervalsByNs interval list file path
    @type scatter_intervals_path: str | None
    @ivar number_of_tiles_cohort: Number of genomic tiles for scattering in stage variant_calling_process_cohort
    @type number_of_tiles_cohort: int | None
    @ivar number_of_chunks_cohort: Number of chunks for gathering in stage variant_calling_process_cohort
    @type number_of_chunks_cohort: int | None
    @ivar number_of_tiles_somatic: Number of genomic tiles for scattering in stage variant_calling_somatic
    @type number_of_tiles_somatic: int | None
    @ivar number_of_chunks_somatic: Number of chunks for gathering in stage variant_calling_somatic
    @type number_of_chunks_somatic: int | None
    @ivar gatk_bundle_version: GATK resource bundle version
    @type gatk_bundle_version: str | None
    @ivar snpeff_genome_version: snpEff genome version
    @type snpeff_genome_version: str | None
    @ivar genome_annotation_gtf: Genome annotation Gene Transfer Format (GTF) file path
    @type genome_annotation_gtf: str | None
    @ivar vep_annotation: Ensembl Variant Effect Predictor (VEP) annotation type (i.e. ensembl, refseq or merged)
    @type vep_annotation: str
    @ivar vep_assembly: Ensembl Variant Effect Predictor (VEP) assembly
    @type vep_assembly: str | None
    @ivar vep_cache: Ensembl Variant Effect Predictor (VEP) cache directory
    @type vep_cache: str | None
    @ivar vep_fasta: Ensembl Variant Effect Predictor (VEP) FASTA directory
    @type vep_fasta: str | None
    @ivar vep_plugin: Ensembl Variant Effect Predictor (VEP) plug-in directory
    @type vep_plugin: str | None
    @ivar vep_species: Ensembl Variant Effect Predictor (VEP) species
    @type vep_species: str | None
    @ivar vep_source: Ensembl Variant Effect Predictor (VEP) source directory
    @type vep_source: str | None
    @ivar vep_sql_user: Ensembl Variant Effect Predictor (VEP) SQL database user name
    @type vep_sql_user: str | None
    @ivar vep_sql_pass: Ensembl Variant Effect Predictor (VEP) SQL database password
    @type vep_sql_pass: str | None
    @ivar vep_sql_host: Ensembl Variant Effect Predictor (VEP) SQL host
    @type vep_sql_host: str | None
    @ivar vep_sql_port: Ensembl Variant Effect Predictor (VEP) SQL TCP/IP port
    @type vep_sql_port: str | None
    @ivar vep_ofc_path: Ensembl Variant Effect Predictor (VEP) output fields configuration (TSV) file path
    @type vep_ofc_path: str | None
    @ivar vep_soc_path: Ensembl Variant Effect Predictor (VEP) Sequence Ontology term (TSV) configuration file path
    @type vep_soc_path: str | None
    @ivar vep_refseq_alignments_path: Ensembl Variant Effect Predictor (VEP) RefSeq alignments (BAM) file path
    @type vep_refseq_alignments_path: str | None
    @ivar vep_plugin_cadd_path: Ensembl Variant Effect Predictor (VEP) CADD file path
    @type vep_plugin_cadd_path: str | None
    @ivar java_archive_fgbio: Fulcrum Genomics (fgbio) Java archive path
    @type java_archive_fgbio: str | None
    @ivar java_archive_gatk: Genome Analysis Tool Kit Java Archive (JAR) file path
    @type java_archive_gatk: str | None
    @ivar java_archive_picard: Picard tools Java Archive (JAR) file path
    @type java_archive_picard: str | None
    @ivar java_archive_snpeff: snpEff tool Java Archive (JAR) file path
    @type java_archive_snpeff: str | None
    @ivar java_archive_vcf_filter: VCF.Filter tool Java Archive (JAR) file path
    @type java_archive_vcf_filter: str | None
    @ivar _tile_region_cohort_list: Python C{list} of C{bsf.interval.Container} objects
    @type _tile_region_cohort_list: list[Container] | None
    @ivar _tile_region_somatic_list: Python C{list} of C{bsf.interval.Container} objects
    @type _tile_region_somatic_list: list[Container] | None
    """

    name = 'Variant Calling Analysis'
    prefix = 'variant_calling'
    variants_to_table = False

    @classmethod
    def get_stage_name_align_lane(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'align_lane'))

    @classmethod
    def get_stage_name_process_lane(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'process_lane'))

    @classmethod
    def get_stage_name_process_sample(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'process_sample'))

    @classmethod
    def get_stage_name_diagnose_sample(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'diagnose_sample'))

    @classmethod
    def get_stage_name_merge_cohort(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'merge_cohort'))

    @classmethod
    def get_stage_name_process_cohort(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'process_cohort'))

    @classmethod
    def get_stage_name_annotate_cohort_snpeff(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'annotate_cohort_snpeff'))

    @classmethod
    def get_stage_name_annotate_cohort_vep(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'annotate_cohort_vep'))

    @classmethod
    def get_stage_name_split_cohort_snpeff(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'split_cohort_snpeff'))

    @classmethod
    def get_stage_name_split_cohort_vep(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'split_cohort_vep'))

    @classmethod
    def get_stage_name_summary(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'summary'))

    @classmethod
    def get_stage_name_somatic(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'somatic'))

    @classmethod
    def get_stage_name_annotate_somatic_snpeff(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'annotate_somatic_snpeff'))

    @classmethod
    def get_stage_name_annotate_somatic_vep(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'annotate_somatic_vep'))

    @classmethod
    def get_stage_name_split_somatic_snpeff(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'split_somatic_snpeff'))

    @classmethod
    def get_stage_name_split_somatic_vep(cls):
        """Get a particular C{bsf.analysis.Stage.name}.

        @return: C{bsf.analysis.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'split_somatic_vep'))

    @classmethod
    def get_prefix_align_lane(cls, paired_reads_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_align_lane(), paired_reads_name))

    @classmethod
    def get_prefix_process_lane(cls, paired_reads_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_process_lane(), paired_reads_name))

    @classmethod
    def get_prefix_process_sample(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_process_sample(), sample_name))

    @classmethod
    def get_prefix_diagnose_sample(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_diagnose_sample(), sample_name))

    @classmethod
    def get_prefix_merge_cohort(cls, cohort_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_merge_cohort(), cohort_name))

    @classmethod
    def get_prefix_process_cohort(cls, cohort_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_process_cohort(), cohort_name))

    @classmethod
    def get_prefix_annotate_cohort_snpeff(cls, cohort_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_annotate_cohort_snpeff(), cohort_name))

    @classmethod
    def get_prefix_annotate_cohort_vep(cls, cohort_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_annotate_cohort_vep(), cohort_name))

    @classmethod
    def get_prefix_split_cohort_snpeff(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_split_cohort_snpeff(), sample_name))

    @classmethod
    def get_prefix_split_cohort_vep(cls, sample_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param sample_name: C{bsf.ngs.Sample.name}
        @type sample_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_split_cohort_vep(), sample_name))

    @classmethod
    def get_prefix_summary(cls, cohort_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_summary(), cohort_name))

    @classmethod
    def get_prefix_somatic(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_somatic(), comparison_name))

    @classmethod
    def get_prefix_annotate_somatic_snpeff(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_annotate_somatic_snpeff(), comparison_name))

    @classmethod
    def get_prefix_annotate_somatic_vep(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_annotate_somatic_vep(), comparison_name))

    @classmethod
    def get_prefix_split_somatic_snpeff(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_split_somatic_snpeff(), comparison_name))

    @classmethod
    def get_prefix_split_somatic_vep(cls, comparison_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_split_somatic_vep(), comparison_name))

    @classmethod
    def get_file_path_align_lane(cls, paired_reads_name):
        """Get a C{FilePathAlignment} object.

        @param paired_reads_name: PairedReads name
        @type paired_reads_name: str
        @return: C{FilePathAlignment} object
        @rtype: FilePathAlignment
        """
        return FilePathAlignment(
            prefix=cls.get_prefix_align_lane(paired_reads_name=paired_reads_name))

    @classmethod
    def get_file_path_process_read_group(cls, paired_reads_name):
        """Get a C{FilePathProcessReadGroup} object.

        @param paired_reads_name: PairedReads name
        @type paired_reads_name: str
        @return: C{FilePathProcessReadGroup} object
        @rtype: FilePathProcessReadGroup
        """
        return FilePathProcessReadGroup(
            prefix=cls.get_prefix_process_lane(paired_reads_name=paired_reads_name))

    @classmethod
    def get_file_path_process_sample(cls, sample_name):
        """Get a C{FilePathProcessSample} object.

        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathProcessSample} object
        @rtype: FilePathProcessSample
        """
        return FilePathProcessSample(
            prefix=cls.get_prefix_process_sample(sample_name=sample_name))

    @classmethod
    def get_file_path_diagnose_sample(cls, sample_name):
        """Get a C{FilePathDiagnoseSample} object.

        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathDiagnoseSample} object
        @rtype: FilePathDiagnoseSample
        """
        return FilePathDiagnoseSample(
            prefix=cls.get_prefix_diagnose_sample(sample_name=sample_name))

    @classmethod
    def get_file_path_merge_cohort(cls, cohort_name):
        """Get a C{FilePathMergeCohort} object.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: C{FilePathMergeCohort} object
        @rtype: FilePathMergeCohort
        """
        return FilePathMergeCohort(
            prefix=cls.get_prefix_merge_cohort(cohort_name=cohort_name))

    @classmethod
    def get_file_path_process_cohort(cls, cohort_name):
        """Get a C{FilePathProcessCohort} object.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: C{FilePathProcessCohort} object
        @rtype: FilePathProcessCohort
        """
        return FilePathProcessCohort(
            prefix=cls.get_prefix_process_cohort(cohort_name=cohort_name))

    @classmethod
    def get_file_path_annotate_cohort_snpeff(cls, cohort_name):
        """Get a C{FilePathAnnotateSnpEff} object.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: C{FilePathAnnotateSnpEff} object
        @rtype: FilePathAnnotateSnpEff
        """
        return FilePathAnnotateSnpEff(
            prefix=cls.get_prefix_annotate_cohort_snpeff(cohort_name=cohort_name))

    @classmethod
    def get_file_path_annotate_cohort_vep(cls, cohort_name):
        """Get a C{FilePathAnnotateVEP} object.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: C{FilePathAnnotateVEP} object
        @rtype: FilePathAnnotateVEP
        """
        return FilePathAnnotateVEP(
            prefix=cls.get_prefix_annotate_cohort_vep(cohort_name=cohort_name))

    @classmethod
    def get_file_path_annotate_somatic_snpeff(cls, comparison_name):
        """Get a C{FilePathAnnotateSnpEff} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathAnnotateSnpEff} object
        @rtype: FilePathAnnotateSnpEff
        """
        return FilePathAnnotateSnpEff(
            prefix=cls.get_prefix_annotate_somatic_snpeff(comparison_name=comparison_name))

    @classmethod
    def get_file_path_annotate_somatic_vep(cls, comparison_name):
        """Get a C{FilePathAnnotateVEP} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathAnnotateVEP} object
        @rtype: FilePathAnnotateVEP
        """
        return FilePathAnnotateVEP(
            prefix=cls.get_prefix_annotate_somatic_vep(comparison_name=comparison_name))

    @classmethod
    def get_file_path_split_cohort_snpeff(cls, sample_name):
        """Get a C{FilePathSplitCohort} object.

        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathSplitCohort} object
        @rtype: FilePathSplitCohort
        """
        return FilePathSplitCohort(
            prefix=cls.get_prefix_split_cohort_snpeff(sample_name=sample_name))

    @classmethod
    def get_file_path_split_cohort_vep(cls, sample_name):
        """Get a C{FilePathSplitCohort} object.

        @param sample_name: Sample name
        @type sample_name: str
        @return: C{FilePathSplitCohort} object
        @rtype: FilePathSplitCohort
        """
        return FilePathSplitCohort(
            prefix=cls.get_prefix_split_cohort_vep(sample_name=sample_name))

    @classmethod
    def get_file_path_summary(cls, cohort_name):
        """Get a C{FilePathSummary} object.

        @param cohort_name: Cohort name
        @type cohort_name: str
        @return: C{FilePathSummary} object
        @rtype: FilePathSummary
        """
        return FilePathSummary(
            prefix=cls.get_prefix_summary(cohort_name=cohort_name))

    @classmethod
    def get_file_path_somatic(cls, comparison_name):
        """Get a C{FilePathSomatic} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathSomatic} object
        @rtype: FilePathSomatic
        """
        return FilePathSomatic(
            prefix=cls.get_prefix_somatic(comparison_name=comparison_name))

    @classmethod
    def get_file_path_split_somatic_snpeff(cls, comparison_name):
        """Get a C{FilePathSplitSomatic} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathSplitSomatic} object
        @rtype: FilePathSplitSomatic
        """
        return FilePathSplitSomatic(
            prefix=cls.get_prefix_split_somatic_snpeff(comparison_name=comparison_name))

    @classmethod
    def get_file_path_split_somatic_vep(cls, comparison_name):
        """Get a C{FilePathSplitSomatic} object.

        @param comparison_name: Comparison name
        @type comparison_name: str
        @return: C{FilePathSplitSomatic} object
        @rtype: FilePathSplitSomatic
        """
        return FilePathSplitSomatic(
            prefix=cls.get_prefix_split_somatic_vep(comparison_name=comparison_name))

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
            bwa_genome_db=None,
            comparison_path=None,
            cohort_name=None,
            accessory_cohort_gvcfs=None,
            skip_mark_duplicates=None,
            skip_indel_realignment=None,
            known_sites_discovery=None,
            known_sites_realignment=None,
            known_sites_recalibration=None,
            known_somatic_discovery=None,
            annotation_resources_dict=None,
            truth_sensitivity_filter_level_indel=None,
            truth_sensitivity_filter_level_snp=None,
            vqsr_skip_indel=None,
            vqsr_skip_snp=None,
            vqsr_resources_indel_dict=None,
            vqsr_resources_snp_dict=None,
            vqsr_annotations_indel_list=None,
            vqsr_annotations_snp_list=None,
            vqsr_bad_lod_cutoff_indel=None,
            vqsr_bad_lod_cutoff_snp=None,
            vqsr_max_gaussians_pos_indel=4,
            vqsr_max_gaussians_pos_snp=None,
            exclude_intervals_list=None,
            include_intervals_list=None,
            interval_padding=None,
            scatter_intervals_path=None,
            number_of_tiles_cohort=None,
            number_of_chunks_cohort=None,
            number_of_tiles_somatic=None,
            number_of_chunks_somatic=None,
            gatk_bundle_version=None,
            snpeff_genome_version=None,
            genome_annotation_gtf=None,
            vep_annotation=None,
            vep_assembly=None,
            vep_cache=None,
            vep_fasta=None,
            vep_plugin=None,
            vep_species=None,
            vep_source=None,
            vep_sql_user=None,
            vep_sql_pass=None,
            vep_sql_host=None,
            vep_sql_port=None,
            vep_ofc_path=None,
            vep_soc_path=None,
            vep_refseq_alignments_path=None,
            vep_plugin_cadd_path=None,
            java_archive_fgbio=None,
            java_archive_gatk=None,
            java_archive_picard=None,
            java_archive_snpeff=None,
            java_archive_vcf_filter=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATK}.

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
        @param replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects
            for processing or run them separately
        @type replicate_grouping: bool | None
        @param bwa_genome_db: Genome sequence file path with BWA index
        @type bwa_genome_db: str | None
        @param comparison_path: Comparison file path
        @type comparison_path: str | None
        @param cohort_name: Cohort name
        @type cohort_name: str | None
        @param accessory_cohort_gvcfs: Python C{list} of Python C{str} (GVCF file path) objects
        @type accessory_cohort_gvcfs: list[str] | None
        @param skip_mark_duplicates: Skip the Picard MarkDuplicates step
        @type skip_mark_duplicates: bool | None
        @param skip_indel_realignment: Skip the GATK RealignerTargetCreator and GATK IndelRealigner steps
        @type skip_indel_realignment: bool | None
        @param known_sites_discovery: VCF file path for variant discovery via Haplotype Caller or Unified Genotyper
        @type known_sites_discovery: str | None
        @param known_sites_realignment: Python C{list} of Python C{str} VCF file paths
            for InDel realignment
        @type known_sites_realignment: list[str] | None
        @param known_sites_recalibration: Python C{list} of Python C{str} VCF file paths
            for Base Quality Score Recalibration (BQSR)
        @type known_sites_recalibration: list[str] | None
        @param known_somatic_discovery: Catalogue Of Somatic Mutations In Cancer (COSMIC) VCF file path
            for somatic variant discovery via MuTect2
        @type known_somatic_discovery: list[str] | None
        @param annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
            Python C{tuple} of
            Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
        @type annotation_resources_dict: dict[str, (str, list[str])] | None
        @param truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
        @type truth_sensitivity_filter_level_indel: str | None
        @param truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
        @type truth_sensitivity_filter_level_snp: str | None
        @param vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
        @type vqsr_skip_indel: bool | None
        @param vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
        @type vqsr_skip_snp: bool | None
        @param vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_indel_dict: dict[str, dict[str, str]] | None
        @param vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_snp_dict: dict[str, dict[str, str]] | None
        @param vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_indel_list: list[str] | None
        @param vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_snp_list: list[str] | None
        @param vqsr_bad_lod_cutoff_indel: LOD score cutoff for negative training set for INDELs
        @type vqsr_bad_lod_cutoff_indel: float | None
        @param vqsr_bad_lod_cutoff_snp: LOD score cutoff for negative training set for SNPs
        @type vqsr_bad_lod_cutoff_snp: float | None
        @param vqsr_max_gaussians_pos_indel: Maximum number of Gaussians in the positive training for INDELs
        @type vqsr_max_gaussians_pos_indel: int | None
        @param vqsr_max_gaussians_pos_snp: Maximum number of Gaussians in the positive training for SNPs
        @type vqsr_max_gaussians_pos_snp: int | None
        @param exclude_intervals_list: Python C{list} of Python C{str} (intervals) to exclude from the analysis
        @type exclude_intervals_list: list[str] | None
        @param include_intervals_list: Python C{list} of Python C{str} (intervals) to include in the analysis
        @type include_intervals_list: list[str] | None
        @param interval_padding: Interval padding
        @type interval_padding: int | None
        @param scatter_intervals_path: Picard ScatterIntervalsByNs interval list file path
        @type scatter_intervals_path: str | None
        @param number_of_tiles_cohort: Number of genomic tiles for scattering in stage variant_calling_process_cohort
        @type number_of_tiles_cohort: int | None
        @param number_of_chunks_cohort: Number of chunks for gathering in stage variant_calling_process_cohort
        @type number_of_chunks_cohort: int | None
        @param number_of_tiles_somatic: Number of genomic tiles for scattering in stage variant_calling_somatic
        @type number_of_tiles_somatic: int | None
        @param number_of_chunks_somatic: Number of chunks for gathering in stage variant_calling_somatic
        @type number_of_chunks_somatic: int | None
        @param gatk_bundle_version: GATK resource bundle version
        @type gatk_bundle_version: str | None
        @param snpeff_genome_version: snpEff genome version
        @type snpeff_genome_version: str | None
        @param genome_annotation_gtf: Genome annotation Gene Transfer Format (GTF) file path
        @type genome_annotation_gtf: str | None
        @param vep_annotation: Ensembl Variant Effect Predictor (VEP) annotation type (i.e. ensembl, refseq or merged)
        @type vep_annotation: str
        @param vep_assembly: Ensembl Variant Effect Predictor (VEP) assembly
        @type vep_assembly: str | None
        @param vep_fasta: Ensembl Variant Effect Predictor (VEP) FASTA directory
        @type vep_fasta: str | None
        @param vep_cache: Ensembl Variant Effect Predictor (VEP) cache directory
        @type vep_cache: str | None
        @param vep_plugin: Ensembl Variant Effect Predictor (VEP) plug-in directory
        @type vep_plugin: str | None
        @param vep_species: Ensembl Variant Effect Predictor (VEP) species
        @type vep_species: str | None
        @param vep_source: Ensembl Variant Effect Predictor (VEP) source directory
        @type vep_source: str | None
        @param vep_sql_user: Ensembl Variant Effect Predictor (VEP) SQL database user name
        @type vep_sql_user: str | None
        @param vep_sql_pass: Ensembl Variant Effect Predictor (VEP) SQL database password
        @type vep_sql_pass: str | None
        @param vep_sql_host: Ensembl Variant Effect Predictor (VEP) SQL host
        @type vep_sql_host: str | None
        @param vep_sql_port: Ensembl Variant Effect Predictor (VEP) SQL TCP/IP port
        @type vep_sql_port: str | None
        @param vep_ofc_path: Ensembl Variant Effect Predictor (VEP) output fields configuration (TSV) file path
        @type vep_ofc_path: str | None
        @param vep_soc_path: Ensembl Variant Effect Predictor (VEP) Sequence Ontology term (TSV) configuration file path
        @type vep_soc_path: str | None
        @param vep_refseq_alignments_path: Ensembl Variant Effect Predictor (VEP) RefSeq alignments (BAM) file path
        @type vep_refseq_alignments_path: str | None
        @param vep_plugin_cadd_path: Ensembl Variant Effect Predictor (VEP) CADD file path
        @type vep_plugin_cadd_path: str | None
        @param java_archive_fgbio: Fulcrum Genomics (fgbio) Java archive path
        @type java_archive_fgbio: str | None
        @param java_archive_gatk: Genome Analysis Tool Kit Java Archive (JAR) file path
        @type java_archive_gatk: str | None
        @param java_archive_picard: Picard tools Java Archive (JAR) file path
        @type java_archive_picard: str | None
        @param java_archive_snpeff: snpEff tool Java Archive (JAR) file path
        @type java_archive_snpeff: str | None
        @param java_archive_vcf_filter: VCF.Filter tool Java Archive (JAR) file path
        @type java_archive_vcf_filter: str | None
        """

        super(VariantCallingGATK, self).__init__(
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
        self.bwa_genome_db = bwa_genome_db
        self.comparison_path = comparison_path
        self.cohort_name = cohort_name
        self.accessory_cohort_gvcfs = accessory_cohort_gvcfs

        self.skip_mark_duplicates = skip_mark_duplicates
        self.skip_indel_realignment = skip_indel_realignment

        self.known_sites_discovery = known_sites_discovery
        self.known_sites_realignment = known_sites_realignment
        self.known_sites_recalibration = known_sites_recalibration
        self.known_somatic_discovery = known_somatic_discovery

        self.annotation_resources_dict = annotation_resources_dict

        self.truth_sensitivity_filter_level_indel = truth_sensitivity_filter_level_indel
        self.truth_sensitivity_filter_level_snp = truth_sensitivity_filter_level_snp

        self.vqsr_skip_indel = vqsr_skip_indel
        self.vqsr_skip_snp = vqsr_skip_snp
        self.vqsr_resources_indel_dict = vqsr_resources_indel_dict
        self.vqsr_resources_snp_dict = vqsr_resources_snp_dict
        self.vqsr_annotations_indel_list = vqsr_annotations_indel_list
        self.vqsr_annotations_snp_list = vqsr_annotations_snp_list
        self.vqsr_bad_lod_cutoff_indel = vqsr_bad_lod_cutoff_indel
        self.vqsr_bad_lod_cutoff_snp = vqsr_bad_lod_cutoff_snp
        self.vqsr_max_gaussians_pos_indel = vqsr_max_gaussians_pos_indel
        self.vqsr_max_gaussians_pos_snp = vqsr_max_gaussians_pos_snp

        self.exclude_intervals_list = exclude_intervals_list
        self.include_intervals_list = include_intervals_list
        self.interval_padding = interval_padding

        self.scatter_intervals_path = scatter_intervals_path
        self.number_of_tiles_cohort = number_of_tiles_cohort
        self.number_of_chunks_cohort = number_of_chunks_cohort
        self.number_of_tiles_somatic = number_of_tiles_somatic
        self.number_of_chunks_somatic = number_of_chunks_somatic
        self.gatk_bundle_version = gatk_bundle_version
        self.snpeff_genome_version = snpeff_genome_version
        self.genome_annotation_gtf = genome_annotation_gtf

        self.vep_annotation = vep_annotation
        self.vep_assembly = vep_assembly
        self.vep_cache = vep_cache
        self.vep_fasta = vep_fasta
        self.vep_plugin = vep_plugin
        self.vep_source = vep_source
        self.vep_species = vep_species
        self.vep_sql_user = vep_sql_user
        self.vep_sql_pass = vep_sql_pass
        self.vep_sql_host = vep_sql_host
        self.vep_sql_port = vep_sql_port
        self.vep_ofc_path = vep_ofc_path
        self.vep_soc_path = vep_soc_path
        self.vep_refseq_alignments_path = vep_refseq_alignments_path
        self.vep_plugin_cadd_path = vep_plugin_cadd_path

        self.java_archive_fgbio = java_archive_fgbio
        self.java_archive_gatk = java_archive_gatk
        self.java_archive_picard = java_archive_picard
        self.java_archive_snpeff = java_archive_snpeff
        self.java_archive_vcf_filter = java_archive_vcf_filter

        self._comparison_dict = dict()
        """ @type _comparison_dict: dict[str, VariantCallingGATKComparison] """

        self._cache_path_dict = None
        """ @type _cache_path_dict: dict[str, str] | None """

        # Initialise the Python list of genome tile regions with an empty region to run a single process by default.

        self._tile_region_cohort_list = None
        """ @type _tile_region_cohort_list: list[Container] | None """

        self._tile_region_somatic_list = None
        """ @type _tile_region_somatic_list: list[Container] | None """

        return

    @property
    def get_gatk_bundle_path(self):
        """Get the absolute GATK bundle directory
        C{bsf.standards.StandardFilePath.get_resource_gatk_bundle()} for the set
        C{bsf.analyses.variant_calling.VariantCallingGATK.gatk_bundle_version} and
        C{bsf.analyses.variant_calling.VariantCallingGATK.genome_version}.

        @return: Absolute GATK bundle directory
        @rtype: str
        """
        return StandardFilePath.get_resource_gatk_bundle(
            gatk_bundle_version=self.gatk_bundle_version,
            genome_version=self.genome_version,
            absolute=True)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.variant_calling.VariantCallingGATK} via a
        C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        def set_vqsr_configuration(vqsr_resources_dict, variation_type):
            """Private function to read variant quality score recalibration (VQSR) configuration information.

            Configuration options I{vqsr_resources_indel} and I{vqsr_resources_snp} provide a comma-separated list of
            resources that are to be used in the VQSR procedure. Each option needs to correspond to a sub-section of
            the C{configparser.ConfigParser} in C{bsf.standards.Configuration.config_parser}.
            Each sub-section needs options 'known', 'training', 'truth', 'prior' and 'file_path'.
            @param vqsr_resources_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
            @type vqsr_resources_dict: dict[str, dict[str, str]] | None
            @param variation_type: Variation type I{indel} or I{snp}
            @type variation_type: str
            @return: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
            @rtype: dict[str, dict[str, str]] | None
            """
            if variation_type not in ('indel', 'snp'):
                raise Exception("Variation type has to be 'indel' or 'snp', not " + repr(variation_type) + '.')

            # The vqsr_resources_indel|snp options of the current configuration section hold a comma-separated list
            # of resources that should correspond to a sub-section in the configuration file.
            vqsr_option = '_'.join(('vqsr_resources', variation_type))
            if config_parser.has_option(section=section, option=vqsr_option):
                if vqsr_resources_dict is None:
                    vqsr_resources_dict = dict()
                # Split the resource list on a comma, split white space characters and remove remaining empty strings.
                for resource_key in configuration.get_list_from_csv(section=section, option=vqsr_option):
                    # The VQSR resource section consists of section.vqsr_(indel|snp)_resource.
                    resource_section = '.'.join((section, '_'.join(('vqsr', variation_type, resource_key))))
                    if config_parser.has_section(section=resource_section):
                        if resource_key not in vqsr_resources_dict:
                            vqsr_resources_dict[resource_key] = dict()
                        resource_dict = vqsr_resources_dict[resource_key]

                        for resource_option in ('known', 'training', 'truth', 'prior', 'file_path'):
                            if config_parser.has_option(section=resource_section, option=resource_option):
                                resource_dict[resource_option] = config_parser.get(
                                    section=resource_section,
                                    option=resource_option)
                            else:
                                raise Exception(
                                    'Missing configuration option ' + repr(resource_option) +
                                    ' in section ' + repr(resource_section) + '.')
                    else:
                        raise Exception(
                            'Missing configuration section ' + repr(resource_section) +
                            ' declared in option ' + repr(vqsr_option) + ' ' +
                            repr(config_parser.get(section=section, option=vqsr_option)) + '.')

            return vqsr_resources_dict

        def set_annotation_configuration(annotation_resources_dict):
            """Private function to read variant annotation configuration information.

            @param annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
                Python C{tuple} of Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
            @type annotation_resources_dict: dict[str, (str, list[str])]
            @return: Python C{dict} of Python C{str} (annotation resource name) key and
                Python C{tuple} of Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
            @rtype: dict[str, (str, list[str])]
            """
            annotation_option = 'annotation_resources'
            if config_parser.has_option(section=section, option=annotation_option):
                if annotation_resources_dict is None:
                    annotation_resources_dict = dict()
                # Split the resource list on a comma, split white space characters and remove remaining empty strings.
                for resource_name in configuration.get_list_from_csv(section=section, option=annotation_option):
                    # The annotation resource section consists of section.annotation_resource.
                    resource_section = '.'.join((section, '_'.join(('annotation', resource_name))))
                    if config_parser.has_section(section=resource_section):
                        resource_option = 'file_path'
                        if config_parser.has_option(section=resource_section, option=resource_option):
                            file_path = config_parser.get(section=resource_section, option=resource_option)
                        else:
                            raise Exception(
                                'Missing configuration option ' + repr(resource_option) +
                                ' in configuration section ' + repr(resource_section) + '.')
                        resource_option = 'annotations'
                        if config_parser.has_option(section=resource_section, option=resource_option):
                            # Split the annotation list on a comma, split white space characters and
                            # remove remaining empty strings.
                            annotation_list = configuration.get_list_from_csv(
                                section=resource_section,
                                option=resource_option)
                        else:
                            raise Exception(
                                'Missing configuration option ' + repr(resource_option) +
                                ' in configuration section ' + repr(resource_section) + '.')
                        # Create a dict key and a tuple of a Python str and Python list.
                        annotation_resources_dict[resource_name] = file_path, annotation_list
                    else:
                        raise Exception(
                            'Missing configuration section ' + repr(resource_section) +
                            ' declared in option ' + repr(annotation_option) + ' ' +
                            repr(config_parser.get(section=section, option=annotation_option)) + '.')

            return annotation_resources_dict

        # Start of set_configuration() method body.

        super(VariantCallingGATK, self).set_configuration(configuration=configuration, section=section)

        config_parser = configuration.config_parser

        option = 'replicate_grouping'
        if config_parser.has_option(section=section, option=option):
            self.replicate_grouping = config_parser.getboolean(section=section, option=option)

        option = 'bwa_genome_db'
        if config_parser.has_option(section=section, option=option):
            self.bwa_genome_db = config_parser.get(section=section, option=option)

        option = 'cmp_file'
        if config_parser.has_option(section=section, option=option):
            self.comparison_path = config_parser.get(section=section, option=option)

        option = 'cohort_name'
        if config_parser.has_option(section=section, option=option):
            self.cohort_name = config_parser.get(section=section, option=option)

        option = 'accessory_cohort_gvcfs'
        if config_parser.has_option(section=section, option=option):
            self.accessory_cohort_gvcfs = configuration.get_list_from_csv(section=section, option=option)

        option = 'skip_mark_duplicates'
        if config_parser.has_option(section=section, option=option):
            self.skip_mark_duplicates = config_parser.getboolean(section=section, option=option)

        option = 'skip_indel_realignment'
        if config_parser.has_option(section=section, option=option):
            self.skip_indel_realignment = config_parser.getboolean(section=section, option=option)

        option = 'truth_sensitivity_filter_level_indel'
        if config_parser.has_option(section=section, option=option):
            self.truth_sensitivity_filter_level_indel = config_parser.get(section=section, option=option)

        option = 'truth_sensitivity_filter_level_snp'
        if config_parser.has_option(section=section, option=option):
            self.truth_sensitivity_filter_level_snp = config_parser.get(section=section, option=option)

        option = 'vqsr_skip_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_skip_indel = config_parser.getboolean(section=section, option=option)

        option = 'vqsr_skip_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_skip_snp = config_parser.getboolean(section=section, option=option)

        option = 'vqsr_annotations_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_annotations_indel_list = configuration.get_list_from_csv(section=section, option=option)

        option = 'vqsr_annotations_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_annotations_snp_list = configuration.get_list_from_csv(section=section, option=option)

        option = 'vqsr_bad_lod_cutoff_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_bad_lod_cutoff_indel = config_parser.getfloat(section=section, option=option)

        option = 'vqsr_bad_lod_cutoff_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_bad_lod_cutoff_snp = config_parser.getfloat(section=section, option=option)

        option = 'vqsr_max_gaussians_pos_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_max_gaussians_pos_indel = config_parser.getint(section=section, option=option)

        option = 'vqsr_max_gaussians_pos_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_max_gaussians_pos_snp = config_parser.getint(section=section, option=option)

        # Set VQSR resources and corresponding configuration sections for INDELs.

        self.vqsr_resources_indel_dict = set_vqsr_configuration(
            vqsr_resources_dict=self.vqsr_resources_indel_dict,
            variation_type='indel')

        # Set VQSR resources and corresponding configuration sections for SNPs.

        self.vqsr_resources_snp_dict = set_vqsr_configuration(
            vqsr_resources_dict=self.vqsr_resources_snp_dict,
            variation_type='snp')

        # Set additionally requested annotation resources for the GATK AnnotateVariants step.

        self.annotation_resources_dict = set_annotation_configuration(
            annotation_resources_dict=self.annotation_resources_dict)

        option = 'known_sites_discovery'
        if config_parser.has_option(section=section, option=option):
            self.known_sites_discovery = config_parser.get(section=section, option=option)

        option = 'known_sites_realignment'
        if config_parser.has_option(section=section, option=option):
            self.known_sites_realignment = configuration.get_list_from_csv(section=section, option=option)

        option = 'known_sites_recalibration'
        if config_parser.has_option(section=section, option=option):
            self.known_sites_recalibration = configuration.get_list_from_csv(section=section, option=option)

        option = 'known_somatic_discovery'
        if config_parser.has_option(section=section, option=option):
            self.known_somatic_discovery = configuration.get_list_from_csv(section=section, option=option)

        option = 'exclude_intervals'
        if config_parser.has_option(section=section, option=option):
            self.exclude_intervals_list = configuration.get_list_from_csv(section=section, option=option)

        option = 'include_intervals'
        if config_parser.has_option(section=section, option=option):
            self.include_intervals_list = configuration.get_list_from_csv(section=section, option=option)

        option = 'interval_padding'
        if config_parser.has_option(section=section, option=option):
            self.interval_padding = config_parser.getint(section=section, option=option)

        option = 'scatter_intervals_path'
        if config_parser.has_option(section=section, option=option):
            self.scatter_intervals_path = config_parser.get(section=section, option=option)

        option = 'number_of_tiles_cohort'
        if config_parser.has_option(section=section, option=option):
            self.number_of_tiles_cohort = config_parser.getint(section=section, option=option)

        option = 'number_of_chunks_cohort'
        if config_parser.has_option(section=section, option=option):
            self.number_of_chunks_cohort = config_parser.getint(section=section, option=option)

        option = 'number_of_tiles_somatic'
        if config_parser.has_option(section=section, option=option):
            self.number_of_tiles_somatic = config_parser.getint(section=section, option=option)

        option = 'number_of_chunks_somatic'
        if config_parser.has_option(section=section, option=option):
            self.number_of_chunks_somatic = config_parser.getint(section=section, option=option)

        option = 'gatk_bundle_version'
        if config_parser.has_option(section=section, option=option):
            self.gatk_bundle_version = config_parser.get(section=section, option=option)

        option = 'snpeff_genome_version'
        if config_parser.has_option(section=section, option=option):
            self.snpeff_genome_version = config_parser.get(section=section, option=option)

        option = 'genome_annotation_gtf'
        if config_parser.has_option(section=section, option=option):
            self.genome_annotation_gtf = config_parser.get(section=section, option=option)

        option = 'vep_annotation'
        if config_parser.has_option(section=section, option=option):
            self.vep_annotation = config_parser.get(section=section, option=option)

        option = 'vep_assembly'
        if config_parser.has_option(section=section, option=option):
            self.vep_assembly = config_parser.get(section=section, option=option)

        option = 'vep_species'
        if config_parser.has_option(section=section, option=option):
            self.vep_species = config_parser.get(section=section, option=option)

        option = 'vep_cache'
        if config_parser.has_option(section=section, option=option):
            self.vep_cache = config_parser.get(section=section, option=option)

        option = 'vep_fasta'
        if config_parser.has_option(section=section, option=option):
            self.vep_fasta = config_parser.get(section=section, option=option)

        option = 'vep_plugin'
        if config_parser.has_option(section=section, option=option):
            self.vep_plugin = config_parser.get(section=section, option=option)

        option = 'vep_source'
        if config_parser.has_option(section=section, option=option):
            self.vep_source = config_parser.get(section=section, option=option)

        option = 'vep_sql_user'
        if config_parser.has_option(section=section, option=option):
            self.vep_sql_user = config_parser.get(section=section, option=option)

        option = 'vep_sql_pass'
        if config_parser.has_option(section=section, option=option):
            self.vep_sql_pass = config_parser.get(section=section, option=option)

        option = 'vep_sql_host'
        if config_parser.has_option(section=section, option=option):
            self.vep_sql_host = config_parser.get(section=section, option=option)

        option = 'vep_sql_port'
        if config_parser.has_option(section=section, option=option):
            self.vep_sql_port = config_parser.get(section=section, option=option)

        option = 'vep_ofc_path'
        if config_parser.has_option(section=section, option=option):
            self.vep_ofc_path = config_parser.get(section=section, option=option)

        option = 'vep_soc_path'
        if config_parser.has_option(section=section, option=option):
            self.vep_soc_path = config_parser.get(section=section, option=option)

        option = 'vep_refseq_alignments_path'
        if config_parser.has_option(section=section, option=option):
            self.vep_refseq_alignments_path = config_parser.get(section=section, option=option)

        option = 'vep_plugin_cadd_path'
        if config_parser.has_option(section=section, option=option):
            self.vep_plugin_cadd_path = config_parser.get(section=section, option=option)

        option = 'java_archive_fgbio'
        if config_parser.has_option(section=section, option=option):
            self.java_archive_fgbio = config_parser.get(section=section, option=option)

        option = 'java_archive_gatk'
        if config_parser.has_option(section=section, option=option):
            self.java_archive_gatk = config_parser.get(section=section, option=option)

        option = 'java_archive_picard'
        if config_parser.has_option(section=section, option=option):
            self.java_archive_picard = config_parser.get(section=section, option=option)

        option = 'java_archive_snpeff'
        if config_parser.has_option(section=section, option=option):
            self.java_archive_snpeff = config_parser.get(section=section, option=option)

        option = 'java_archive_vcf_filter'
        if config_parser.has_option(section=section, option=option):
            self.java_archive_vcf_filter = config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a C{bsf.analyses.variant_calling.VariantCallingGATK} analysis.
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file from disk.

                - Column headers for CASAVA folders:
                    - Normal/Tumor ProcessedRunFolder:
                        - CASAVA processed run folder name or
                        - C{bsf.analysis.Analysis.input_directory} by default
                    - Normal/Tumor Project:
                        - CASAVA Project name or
                        - C{bsf.analysis.Analysis.project_name} by default
                    - Normal/Tumor Sample:
                        - CASAVA Sample name, no default
                - Column headers for independent samples:
                    - Normal/Tumor Sample:
                    - Normal/Tumor Reads:
                    - Normal/Tumor File:
                - PON Path:
                    - File path to a Panel-Of-Normal (PON) VCF file
            """
            # For variant calling, all samples need adding to the bsf.analysis.Analysis regardless.
            for _sample in self.collection.get_all_samples():
                self.add_sample(sample=_sample)

            if self.comparison_path:
                # A comparison file path has been provided.
                self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)

                annotation_sheet = AnnotationSheet.from_file_path(
                    file_path=self.comparison_path,
                    name='Somatic Comparisons')

                for row_dict in annotation_sheet.row_dicts:
                    if self.debug > 0:
                        print('Comparison sheet row_dict:', row_dict)

                    comparison = VariantCallingGATKComparison()

                    for prefix in ('Normal', 'Tumor'):
                        group_name, group_samples = self.collection.get_samples_from_row_dict(
                            row_dict=row_dict,
                            prefix=prefix)
                        if group_name and len(group_samples):
                            if len(group_samples) != 1:
                                raise Exception(
                                    'Got more than one Sample for class {!r} in comparison {!r}'.format(
                                        prefix,
                                        row_dict))

                            if prefix == 'Normal':
                                comparison.normal_sample = group_samples[0]
                            if prefix == 'Tumor':
                                comparison.tumor_sample = group_samples[0]

                    prefix = 'PON Path'
                    if prefix in row_dict and row_dict[prefix]:
                        comparison.panel_of_normal_path = row_dict[prefix]

                    # At least a tumor Sample has to be defined for the comparison to make sense.
                    if comparison.tumor_sample:
                        self._comparison_dict[comparison.get_name] = comparison

            return

        use_cache = False

        variants_to_table_fields = {
            'fixed': ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'),
            'haplotype_caller': {
                'info': ('AC', 'AF', 'AN', 'DB', 'DP', 'QD', 'VQSLOD', 'culprit'),
                'format': (
                    'AD', 'DP', 'GQ', 'GT', 'MIN_DP', 'PGT', 'PID', 'PL', 'RGQ', 'SB',
                ),
            },
            'mutect2': {
                'info': ('DB', 'MAX_ED', 'MIN_ED', 'NLOD', 'PON', 'RPA', 'RU', 'STR', 'TLOD'),
                'format': (
                    'AD', 'AF', 'ALT_F1R2', 'ALT_F2R1', 'DP', 'FOXOG', 'GQ', 'GT', 'PGT', 'PID', 'PL',
                    'QSS', 'REF_F1R2', 'REF_F2R1',
                ),
            },
            'snpeff': (
                'SNPEFF_EFFECT',
                'SNPEFF_IMPACT',
                'SNPEFF_FUNCTIONAL_CLASS',
                'SNPEFF_CODON_CHANGE',
                'SNPEFF_AMINO_ACID_CHANGE',
                'SNPEFF_GENE_NAME',
                'SNPEFF_GENE_BIOTYPE',
                'SNPEFF_TRANSCRIPT_ID',
                'SNPEFF_EXON_ID',
            ),
            'vep': (
                'VEP_AA_AF',
                'VEP_AF',
                'VEP_AFR_AF',
                'VEP_ALLELE_NUM',
                'VEP_AMR_AF',
                'VEP_APPRIS',
                'VEP_Allele',
                'VEP_Amino_acids',
                'VEP_BIOTYPE',
                'VEP_CADD_PHRED',
                'VEP_CADD_RAW',
                'VEP_CANONICAL',
                'VEP_CCDS',
                'VEP_CDS_position',
                'VEP_CLIN_SIG',
                'VEP_Codons',
                'VEP_Consequence',
                'VEP_DISTANCE',
                'VEP_DOMAINS',
                'VEP_EAS_AF',
                'VEP_EA_AF',
                'VEP_ENSP',
                'VEP_EUR_AF',
                'VEP_EXON_AFFECTED',
                'VEP_EXON_TOTAL',
                # 'VEP_ExAC_AF',
                # 'VEP_ExAC_AFR_AF',
                # 'VEP_ExAC_AMR_AF',
                # 'VEP_ExAC_Adj_AF',
                # 'VEP_ExAC_EAS_AF',
                # 'VEP_ExAC_FIN_AF',
                # 'VEP_ExAC_NFE_AF',
                # 'VEP_ExAC_OTH_AF',
                # 'VEP_ExAC_SAS_AF',
                'VEP_gnomAD_AF',
                'VEP_gnomAD_AFR_AF',
                'VEP_gnomAD_AMR_AF',
                'VEP_gnomAD_ASJ_AF',
                'VEP_gnomAD_EAS_AF',
                'VEP_gnomAD_FIN_AF',
                'VEP_gnomAD_NFE_AF',
                'VEP_gnomAD_OTH_AF',
                'VEP_gnomAD_SAS_AF',
                'VEP_Existing_variation',
                'VEP_FLAGS',
                'VEP_Feature',
                'VEP_Feature_type',
                'VEP_GENE_PHENO',
                'VEP_Gene',
                'VEP_HGNC_ID',
                'VEP_HGVS_OFFSET',
                'VEP_HGVSc',
                'VEP_HGVSp',
                'VEP_HIGH_INF_POS',
                'VEP_IMPACT',
                'VEP_INTRON_AFFECTED',
                'VEP_INTRON_TOTAL',
                'VEP_MAX_AF',
                'VEP_MAX_AF_POPS',
                'VEP_MOTIF_NAME',
                'VEP_MOTIF_POS',
                'VEP_MOTIF_SCORE_CHANGE',
                'VEP_PHENO',
                'VEP_PICK',
                'VEP_PUBMED',
                'VEP_PolyPhen',
                'VEP_Protein_position',
                'VEP_SAS_AF',
                'VEP_SIFT',
                'VEP_SOMATIC',
                'VEP_STRAND',
                'VEP_SWISSPROT',
                'VEP_SYMBOL',
                'VEP_SYMBOL_SOURCE',
                'VEP_TREMBL',
                'VEP_TSL',
                'VEP_UNIPARC',
                'VEP_VARIANT_CLASS',
                'VEP_cDNA_position_end',
                'VEP_cDNA_position_start',
            ),
        }

        def run_merge_cohort_scatter_gather(analysis_stage, cohort_runnable_dict, cohort_name):
            """Private method to hierarchically merge samples into cohorts using a scatter gather approach.

            This method merges GVCF file paths from individual process_sample bsf.procedure.Runnable objects or
            an accessory cohort.
            @param analysis_stage: C{bsf.analysis.Analysis} C{bsf.analysis.Stage}
            @type analysis_stage: Stage
            @param cohort_runnable_dict: Python C{dict} of Python C{str} key and Python C{list} of
                C{bsf.procedure.Runnable}, Python C{str} object value data
            @type cohort_runnable_dict: dict[str, list[(Runnable | str, str)]]
            @param cohort_name: Cohort name to select a Python list of C{bsf.procedure.Runnable} objects from the
                I{cohort_runnable_dict} Python C{dict}
            @type cohort_name: str
            @return: A Python C{tuple} of the final C{bsf.procedure.Runnable} of the gather stage and
                the gVCF file path
            @rtype: (Runnable, str)
            """

            # Private variables are prefixed with an underscore to avoid clashes with variables in the run() method.

            prefix_merge_cohort_final = '_'.join((analysis_stage.name, cohort_name))

            file_path_merge_cohort_final = FilePathMergeCohort(prefix=prefix_merge_cohort_final)

            # If the final TBI index file already exists, create the bsf.procedure.Runnable objects,
            # but do not submit their corresponding bsf.process.Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_merge_cohort_final.combined_gvcf_tbi)):
                final_index_exists = False
            else:
                final_index_exists = True

            # The cohort_object_list contains either bsf.procedure.Runnable objects from the process_sample stage or
            # Python str (GVCF file path) objects for accessory cohorts to be merged.
            cohort_object_list = cohort_runnable_dict[cohort_name]

            # Scatter
            runnable_scatter = None
            """ @type runnable_scatter: Runnable | None """
            runnable_scatter_list = list()
            """ @type runnable_scatter_list: list[Runnable] """

            tile_index_list = range(0, len(self._tile_region_cohort_list))

            for tile_index in tile_index_list:
                prefix_merge_cohort_scatter = '_'.join((analysis_stage.name, cohort_name, 'scatter', str(tile_index)))

                file_path_merge_cohort_scatter = FilePathMergeCohort(prefix=prefix_merge_cohort_scatter)

                runnable_scatter = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=prefix_merge_cohort_scatter,
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        debug=self.debug))
                executable_scatter = self.set_stage_runnable(
                    stage=analysis_stage,
                    runnable=runnable_scatter)
                # Submit the bsf.process.Executable only, if the final TBI index file does not exist,
                # but do not override the state set by the bsf.analysis.Analysis.set_stage_runnable() method.
                if not final_index_exists:
                    executable_scatter.submit = False
                for cohort_component, cohort_component_prefix in cohort_object_list:
                    # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
                    # Set them only for bsf.procedure.Runnable objects,
                    # but not for Python str (file path) objects.
                    if isinstance(cohort_component, Runnable):
                        executable_scatter.dependencies.append(cohort_component.name)

                runnable_scatter_list.append(runnable_scatter)

                if use_cache:
                    reference_scatter = runnable_scatter.get_cache_file_path(
                        file_path=self.bwa_genome_db,
                        absolute=True)
                else:
                    reference_scatter = self.bwa_genome_db

                # Run GATK CombineGVCFs

                _runnable_step = RunnableStepGATK(
                    name='merge_cohort_gatk_combine_gvcfs',
                    java_temporary_path=runnable_scatter.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx8G',
                    java_jar_path=self.java_archive_gatk)
                runnable_scatter.add_runnable_step(runnable_step=_runnable_step)

                _runnable_step.add_gatk_option(key='analysis_type', value='CombineGVCFs')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_scatter)
                if self.exclude_intervals_list:  # not None and not empty
                    for _interval in self.exclude_intervals_list:
                        _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                for _interval in self._tile_region_cohort_list[tile_index].interval_list:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if _interval.name:
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value=_interval.to_gatk_interval(),
                            override=True)
                for cohort_component, _file_path_gvcf in cohort_object_list:
                    # Regardless of the cohort_component class, set the file path.
                    _runnable_step.add_gatk_option(key='variant', value=_file_path_gvcf, override=True)
                _runnable_step.add_gatk_option(
                    key='out',
                    value=file_path_merge_cohort_scatter.combined_gvcf_vcf)

            # Gather
            runnable_gather = None
            """ @type runnable_gather: Runnable | None """

            if len(self._tile_region_cohort_list) == 1:
                # If there is only one tile, no need to gather.
                # Assign the sole scatter bsf.procedure.Runnable to the sole gather bsf.procedure.Runnable.
                runnable_gather = runnable_scatter
            else:
                # Gather by hierarchically merging by the number of chunks on the partitioned genome tile index list.
                # Initialise a list of bsf.procedure.Runnable objects and indices for the hierarchical merge.
                runnable_gather_list = runnable_scatter_list
                gather_level = 0
                while len(tile_index_list) > 1:
                    temporary_runnable_gather_list = list()
                    """ @type temporary_runnable_gather_list: list[Runnable] """
                    temporary_tile_index_list = list()
                    """ @type temporary_tile_index_list: list[int] """
                    # Partition the index list into chunks of given size.
                    partition_list = [tile_index_list[offset:offset + self.number_of_chunks_cohort]
                                      for offset in range(0, len(tile_index_list), self.number_of_chunks_cohort)]

                    for partition_index, chunk_index_list in enumerate(partition_list):
                        # The file prefix includes the level and partition index.
                        prefix_merge_cohort_gather = '_'.join((
                            analysis_stage.name,
                            cohort_name,
                            'gather',
                            str(gather_level),
                            str(partition_index)))

                        file_path_merge_cohort_gather = FilePathMergeCohort(prefix=prefix_merge_cohort_gather)

                        runnable_gather = self.add_runnable_consecutive(
                            runnable=ConsecutiveRunnable(
                                name=prefix_merge_cohort_gather,
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                debug=self.debug))
                        executable_gather = self.set_stage_runnable(
                            stage=analysis_stage,
                            runnable=runnable_gather)
                        # Submit the bsf.process.Executable only, if the final TBI index file does not exist,
                        # but do not override the state set by the bsf.analysis.Analysis.set_stage_runnable() method.
                        if not final_index_exists:
                            executable_gather.submit = False
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_runnable_gather_list.append(runnable_gather)
                        temporary_tile_index_list.append(partition_index)

                        if use_cache:
                            reference_gather = runnable_gather.get_cache_file_path(
                                file_path=self.bwa_genome_db,
                                absolute=True)
                        else:
                            reference_gather = self.bwa_genome_db

                        # GATK CatVariants bypasses the GATK engine and thus requires a completely different
                        # command line.
                        _runnable_step = RunnableStepJava(
                            name='merge_cohort_gatk_cat_variants',
                            sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                            java_temporary_path=runnable_gather.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx4G')
                        runnable_gather.add_runnable_step(runnable_step=_runnable_step)

                        _runnable_step.add_option_short(key='classpath', value=self.java_archive_gatk)
                        _sub_command = _runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        _sub_command.add_option_long(key='reference', value=reference_gather)
                        _sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_merge_cohort_gather.combined_gvcf_vcf)
                        _sub_command.add_switch_long(key='assumeSorted')
                        # Finally, process per chunk index.
                        for chunk_index in chunk_index_list:
                            _runnable_scatter_or_gather = runnable_gather_list[chunk_index]
                            _file_path_scatter_or_gather = FilePathMergeCohort(prefix=_runnable_scatter_or_gather.name)
                            # Set GATK option variant
                            _sub_command.add_option_long(
                                key='variant',
                                value=_file_path_scatter_or_gather.combined_gvcf_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            _runnable_step.obsolete_file_path_list.append(
                                _file_path_scatter_or_gather.combined_gvcf_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            _runnable_step.obsolete_file_path_list.append(
                                _file_path_scatter_or_gather.combined_gvcf_tbi)
                            # Set dependencies on preceding bsf.procedure.Runnable.name or
                            # bsf.process.Executable.name objects.
                            # Depend on the bsf.procedure.Runnable.name of the corresponding
                            # bsf.procedure.Runnable of the scattering above.
                            executable_gather.dependencies.append(_runnable_scatter_or_gather.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    runnable_gather_list = temporary_runnable_gather_list
                    tile_index_list = temporary_tile_index_list
                    gather_level += 1

            # For the last gather bsf.procedure.Runnable, move (rename) file paths to the top-level prefix and
            # adjust the FilePath object accordingly.
            # NOTE: The renaming of merged cohort files from *_scatter_0_* or *_gather_1_0_* to
            # *_combined.g.vcf causes a lot of extra work, because the file names can no longer be
            # deduced form the Runnable.name.
            file_path_merge_cohort_gather = FilePathMergeCohort(prefix=runnable_gather.name)

            _runnable_step = RunnableStepMove(
                name='merge_cohort_gather_move_vcf',
                source_path=file_path_merge_cohort_gather.combined_gvcf_vcf,
                target_path=file_path_merge_cohort_final.combined_gvcf_vcf)
            runnable_gather.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStepMove(
                name='merge_cohort_gather_move_tbi',
                source_path=file_path_merge_cohort_gather.combined_gvcf_tbi,
                target_path=file_path_merge_cohort_final.combined_gvcf_tbi)
            runnable_gather.add_runnable_step(runnable_step=_runnable_step)

            file_path_merge_cohort_gather.combined_gvcf_vcf = file_path_merge_cohort_final.combined_gvcf_vcf
            file_path_merge_cohort_gather.combined_gvcf_tbi = file_path_merge_cohort_final.combined_gvcf_tbi

            return runnable_gather, file_path_merge_cohort_final.combined_gvcf_vcf

        def run_genotype_cohort_scatter_gather(file_path_cohort_gvcf):
            """Private function to genotype a cohort in a scatter and gather approach.

            @param file_path_cohort_gvcf: Cohort-level gVCF file path
            @type file_path_cohort_gvcf: str
            @return: Final C{bsf.procedure.Runnable} of the gather stage
            @rtype: Runnable
            """

            file_path_genotype_cohort_final = FilePathGenotypeCohort(
                prefix=self.get_prefix_process_cohort(
                    cohort_name=self.cohort_name))

            # If the final TBI index file already exists, create the bsf.procedure.Runnable objects,
            # but do not submit their corresponding bsf.process.Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_genotype_cohort_final.genotyped_raw_tbi)):
                final_index_exists = False
            else:
                final_index_exists = True

            runnable_scatter = None
            """ @type runnable_scatter: Runnable | None """
            runnable_scatter_list = list()
            """ @type runnable_scatter_list: list[Runnable] """

            tile_index_list = range(0, len(self._tile_region_cohort_list))

            for tile_index in tile_index_list:
                prefix_process_cohort_scatter = '_'.join((
                    stage_process_cohort.name, self.cohort_name, 'scatter', str(tile_index)))

                file_path_genotype_cohort_scatter = FilePathGenotypeCohort(prefix=prefix_process_cohort_scatter)

                runnable_scatter = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=prefix_process_cohort_scatter,
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        debug=self.debug))
                executable_scatter = self.set_stage_runnable(
                    stage=stage_process_cohort,
                    runnable=runnable_scatter)
                # Submit the bsf.process.Executable only, if the final TBI index file does not exist,
                # but do not override the state set by the bsf.analysis.Analysis.set_stage_runnable() method.
                if not final_index_exists:
                    executable_scatter.submit = False
                # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
                executable_scatter.dependencies.append(runnable_merge_cohort.name)

                runnable_scatter_list.append(runnable_scatter)

                if use_cache:
                    reference_scatter = runnable_scatter.get_cache_file_path(
                        file_path=self.bwa_genome_db,
                        absolute=True)
                else:
                    reference_scatter = self.bwa_genome_db

                # Run the GATK GenotypeGVCFs analysis.

                _runnable_step = RunnableStepGATK(
                    name='process_cohort_gatk_genotype_gvcfs_scatter',
                    java_temporary_path=runnable_scatter.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx12G',
                    java_jar_path=self.java_archive_gatk)
                runnable_scatter.add_runnable_step(runnable_step=_runnable_step)

                _runnable_step.add_gatk_option(key='analysis_type', value='GenotypeGVCFs')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_scatter)
                if self.exclude_intervals_list:  # not None and not empty
                    for _interval in self.exclude_intervals_list:
                        _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                for _interval in self._tile_region_cohort_list[tile_index].interval_list:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if _interval.name:
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value=_interval.to_gatk_interval(),
                            override=True)
                if self.known_sites_discovery:
                    _runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
                _runnable_step.add_gatk_option(key='variant', value=file_path_cohort_gvcf)
                _runnable_step.add_gatk_option(key='out', value=file_path_genotype_cohort_scatter.genotyped_raw_vcf)

            # Gather
            runnable_gather = None
            """ @type runnable_gather: Runnable | None """

            if len(self._tile_region_cohort_list) == 1:
                # If there is only one tile, no need to gather.
                # Assign the sole scatter bsf.procedure.Runnable to the sole gather bsf.procedure.Runnable.
                runnable_gather = runnable_scatter
            else:
                # Gather by hierarchically merging by the number of chunks on the partitioned genome tile index list.
                # Initialise a list of bsf.procedure.Runnable objects and indices for the hierarchical merge.
                runnable_gather_list = runnable_scatter_list
                gather_level = 0
                while len(tile_index_list) > 1:
                    temporary_runnable_gather_list = list()
                    """ @type temporary_runnable_gather_list: list[Runnable] """
                    temporary_tile_index_list = list()
                    """ @type temporary_tile_index_list: list[int] """
                    # Partition the index list into chunks of given size.
                    partition_list = [tile_index_list[offset:offset + self.number_of_chunks_cohort]
                                      for offset in range(0, len(tile_index_list), self.number_of_chunks_cohort)]

                    for partition_index, chunk_index_list in enumerate(partition_list):
                        # The file prefix includes the level and partition index.
                        prefix_process_cohort_gather = '_'.join(
                            (stage_process_cohort.name,
                             self.cohort_name,
                             'gather',
                             str(gather_level),
                             str(partition_index)))

                        file_path_genotype_cohort_gather = FilePathGenotypeCohort(prefix=prefix_process_cohort_gather)

                        runnable_gather = self.add_runnable_consecutive(
                            runnable=ConsecutiveRunnable(
                                name=prefix_process_cohort_gather,
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                debug=self.debug))
                        executable_gather = self.set_stage_runnable(
                            stage=stage_process_cohort,
                            runnable=runnable_gather)
                        # Submit the bsf.process.Executable only, if the final TBI index file does not exist,
                        # but do not override the state set by the bsf.analysis.Analysis.set_stage_runnable() method.
                        if not final_index_exists:
                            executable_gather.submit = False
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_runnable_gather_list.append(runnable_gather)
                        temporary_tile_index_list.append(partition_index)

                        if use_cache:
                            reference_gather = runnable_gather.get_cache_file_path(
                                file_path=self.bwa_genome_db,
                                absolute=True)
                        else:
                            reference_gather = self.bwa_genome_db

                        # GATK CatVariants by-passes the GATK engine and thus requires a completely different
                        # command line.
                        _runnable_step = RunnableStepJava(
                            name='merge_cohort_gatk_cat_variants',
                            sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                            java_temporary_path=runnable_gather.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx4G')
                        runnable_gather.add_runnable_step(runnable_step=_runnable_step)

                        _runnable_step.add_option_short(key='classpath', value=self.java_archive_gatk)
                        _sub_command = _runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        _sub_command.add_option_long(key='reference', value=reference_gather)
                        _sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_genotype_cohort_gather.genotyped_raw_vcf)
                        _sub_command.add_switch_long(key='assumeSorted')
                        # Finally, add bsf.process.RunnableStep options, obsolete files and
                        # bsf.process.Executable dependencies per chunk index.
                        for chunk_index in chunk_index_list:
                            runnable_scatter_or_gather = runnable_gather_list[chunk_index]
                            _file_path_scatter_or_gather = FilePathGenotypeCohort(
                                prefix=runnable_scatter_or_gather.name)
                            # Set GATK option variant
                            _sub_command.add_option_long(
                                key='variant',
                                value=_file_path_scatter_or_gather.genotyped_raw_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            _runnable_step.obsolete_file_path_list.append(
                                _file_path_scatter_or_gather.genotyped_raw_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            _runnable_step.obsolete_file_path_list.append(
                                _file_path_scatter_or_gather.genotyped_raw_tbi)
                            # Set dependencies on preceding bsf.procedure.Runnable.name or
                            # bsf.process.Executable.name objects.
                            executable_gather.dependencies.append(runnable_scatter_or_gather.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    runnable_gather_list = temporary_runnable_gather_list
                    tile_index_list = temporary_tile_index_list
                    gather_level += 1

            # For the last gather bsf.procedure.Runnable, move (rename) file paths to the top-level prefix and
            # adjust the FilePath object accordingly.
            file_path_genotype_cohort_gather = FilePathGenotypeCohort(prefix=runnable_gather.name)

            _runnable_step = RunnableStepMove(
                name='process_cohort_gather_move_vcf',
                source_path=file_path_genotype_cohort_gather.genotyped_raw_vcf,
                target_path=file_path_genotype_cohort_final.genotyped_raw_vcf)
            runnable_gather.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStepMove(
                name='process_cohort_gather_move_tbi',
                source_path=file_path_genotype_cohort_gather.genotyped_raw_tbi,
                target_path=file_path_genotype_cohort_final.genotyped_raw_tbi)
            runnable_gather.add_runnable_step(runnable_step=_runnable_step)

            file_path_genotype_cohort_gather.genotyped_raw_vcf = file_path_genotype_cohort_final.genotyped_raw_vcf
            file_path_genotype_cohort_gather.genotyped_raw_tbi = file_path_genotype_cohort_final.genotyped_raw_tbi

            return runnable_gather

        def run_somatic_scatter_gather(comparison_key):
            """Private method to run somatic variant calling in a scatter and gather approach.

            @param comparison_key: C{VariantCallingGATKComparison.name}
            @type comparison_key: str
            @return: Final C{bsf.procedure.Runnable} of the gather stage
            @rtype: Runnable
            """

            # Private variables are prefixed with an underscore to avoid clashes with variables in the run() method.

            file_path_somatic_final = self.get_file_path_somatic(comparison_name=comparison_key)

            # If the final TBI index file already exists, create the bsf.procedure.Runnable objects,
            # but do not submit their corresponding bsf.process.Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_somatic_final.somatic_tbi)):
                final_index_exists = False
            else:
                final_index_exists = True

            comparison = self._comparison_dict[comparison_key]

            _calling_intervals = VariantCallingGATKCallingIntervals.from_sample(sample=comparison.tumor_sample)

            _target_intervals = VariantCallingGATKTargetIntervals.from_sample(sample=comparison.tumor_sample)

            # FIXME: New sample-specific code, which should replace the code below.
            # For the moment, the HaplotypeCaller is run on the entire interval file without additional
            # scatter gather. Somatic calling with MuTect2 has scatter gather implemented.
            #
            # if _calling_intervals.calling_path:
            #     # Sample-specific calling intervals are available ...
            #     _tile_region_somatic_list = get_interval_tiles(
            #         interval_path=_calling_intervals.calling_path,
            #         tile_number=self.number_of_tiles_somatic)
            # elif _target_intervals.targets_path:
            #     # Sample-specific target intervals are available ...
            #     _tile_region_somatic_list = get_interval_tiles(
            #         interval_path=_target_intervals.targets_path,
            #         tile_number=self.number_of_tiles_somatic)
            # elif self.scatter_intervals_path:
            #     # ACGTmer intervals are available ...
            #     _tile_region_somatic_list = get_interval_tiles(
            #         interval_path=self.scatter_intervals_path,
            #         tile_number=self.number_of_tiles_somatic)
            # else:
            #     # Only genome intervals are available ...
            #     _tile_region_somatic_list = get_genome_tiles(
            #         dictionary_path=dictionary_path,  # FIXME: implement globally?
            #         tile_number=self.number_of_tiles_somatic)

            # Scatter
            runnable_scatter = None
            """ @type runnable_scatter: Runnable | None """
            runnable_scatter_list = list()
            """ @type runnable_scatter_list: list[Runnable] """

            tile_index_list = range(0, len(self._tile_region_somatic_list))

            for tile_index in tile_index_list:
                prefix_somatic_scatter = '_'.join(
                    (self.get_stage_name_somatic(), comparison_key, 'scatter', str(tile_index)))

                file_path_somatic_scatter = FilePathSomaticScatterGather(prefix=prefix_somatic_scatter)

                runnable_scatter = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=prefix_somatic_scatter,
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        debug=self.debug))
                executable_scatter = self.set_stage_runnable(
                    stage=stage_somatic,
                    runnable=runnable_scatter)
                # Submit the bsf.process.Executable only, if the final TBI index file does not exist,
                # but do not override the state set by the bsf.analysis.Analysis.set_stage_runnable() method.
                if not final_index_exists:
                    executable_scatter.submit = False
                # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
                if comparison.normal_sample:
                    executable_scatter.dependencies.append(
                        self.get_prefix_process_sample(sample_name=comparison.normal_sample.name))
                if comparison.tumor_sample:
                    executable_scatter.dependencies.append(
                        self.get_prefix_process_sample(sample_name=comparison.tumor_sample.name))

                runnable_scatter_list.append(runnable_scatter)

                if use_cache:
                    reference_scatter = runnable_scatter.get_cache_file_path(
                        file_path=self.bwa_genome_db,
                        absolute=True)
                else:
                    reference_scatter = self.bwa_genome_db

                # Run GATK MuTect2

                _runnable_step = RunnableStepGATK(
                    name='somatic_gatk_mutect2_scatter',
                    java_temporary_path=runnable_scatter.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx4G',
                    java_jar_path=self.java_archive_gatk)
                runnable_scatter.add_runnable_step(runnable_step=_runnable_step)

                _runnable_step.add_gatk_option(key='analysis_type', value='MuTect2')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_scatter)
                if self.exclude_intervals_list:  # not None and not empty
                    for _interval in self.exclude_intervals_list:
                        _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                for _interval in self._tile_region_somatic_list[tile_index].interval_list:
                    # The list of tiles is initialised to at least one empty Interval to trigger at least one process.
                    # Test for a non-empty interval first before testing for more specialised cases.
                    if _interval:
                        # For a real genome interval on the tile list ...
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value=_interval.to_gatk_interval(),
                            override=True)
                    elif _calling_intervals.calling_path:
                        # If not running on genome tiles, the MuTect2 analysis could be run on calling intervals,
                        # preferentially.
                        _runnable_step.add_gatk_option(key='intervals', value=_calling_intervals.calling_path)
                        if self.interval_padding is not None:
                            _runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                    elif _target_intervals.targets_path:
                        # If not running on genome tiles, the MuTect2 analysis could be run on target intervals, only.
                        _runnable_step.add_gatk_option(key='intervals', value=_target_intervals.targets_path)
                        if self.interval_padding is not None:
                            _runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                if self.known_sites_discovery:
                    _runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
                if self.known_somatic_discovery:  # not None and not empty
                    for _file_path in self.known_somatic_discovery:
                        _runnable_step.add_gatk_option(key='cosmic', value=_file_path, override=True)

                # Find and add the FilePathProcessSample object for the 'normal' Sample object.
                if comparison.normal_sample:
                    _file_path_process_sample = self.get_file_path_process_sample(
                        sample_name=comparison.normal_sample.name)
                    _runnable_step.add_gatk_option(
                        key='input_file:normal',
                        value=_file_path_process_sample.realigned_bam)
                elif comparison.panel_of_normal_path:
                    _runnable_step.add_gatk_option(
                        key='normal_panel',
                        value=comparison.panel_of_normal_path)

                # Find and add the FilePathProcessSample object for the 'tumor' Sample object.
                if comparison.tumor_sample:
                    _file_path_process_sample = self.get_file_path_process_sample(
                        sample_name=comparison.tumor_sample.name)
                    _runnable_step.add_gatk_option(
                        key='input_file:tumor',
                        value=_file_path_process_sample.realigned_bam)

                _runnable_step.add_gatk_option(key='out', value=file_path_somatic_scatter.somatic_vcf)

            # Gather
            runnable_gather = None
            """ @type runnable_gather: Runnable | None """

            if len(self._tile_region_somatic_list) == 1:
                # If there is only one tile, no need to gather.
                # Assign the sole scatter bsf.procedure.Runnable to the sole gather bsf.procedure.Runnable.
                runnable_gather = runnable_scatter
            else:
                # Gather by hierarchically merging by the number of chunks on the partitioned genome tile index list.
                # Initialise a list of bsf.procedure.Runnable objects and indices for the hierarchical merge.
                runnable_gather_list = runnable_scatter_list
                gather_level = 0
                while len(tile_index_list) > 1:
                    temporary_runnable_gather_list = list()
                    """ @type temporary_runnable_gather_list: list[Runnable] """
                    temporary_tile_index_list = list()
                    """ @type temporary_tile_index_list: list[int] """
                    # Partition the index list into chunks of given size.
                    partition_list = [tile_index_list[offset:offset + self.number_of_chunks_somatic]
                                      for offset in range(0, len(tile_index_list), self.number_of_chunks_somatic)]

                    for partition_index, chunk_index_list in enumerate(partition_list):
                        # The file prefix includes the level and partition index.
                        prefix_somatic_gather = '_'.join((self.get_stage_name_somatic(), comparison_key, 'gather',
                                                          str(gather_level), str(partition_index)))

                        file_path_somatic_gather = FilePathSomaticScatterGather(prefix=prefix_somatic_gather)

                        runnable_gather = self.add_runnable_consecutive(
                            runnable=ConsecutiveRunnable(
                                name=prefix_somatic_gather,
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                debug=self.debug))
                        executable_gather = self.set_stage_runnable(
                            stage=stage_somatic,
                            runnable=runnable_gather)
                        # Submit the bsf.process.Executable only, if the final TBI index file does not exist,
                        # but do not override the state set by the bsf.analysis.Analysis.set_stage_runnable() method.
                        if not final_index_exists:
                            executable_gather.submit = False
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_runnable_gather_list.append(runnable_gather)
                        temporary_tile_index_list.append(partition_index)

                        if use_cache:
                            reference_gather = runnable_gather.get_cache_file_path(
                                file_path=self.bwa_genome_db,
                                absolute=True)
                        else:
                            reference_gather = self.bwa_genome_db

                        # GATK CatVariants by-passes the GATK engine and thus requires a completely different
                        # command line.
                        _runnable_step = RunnableStepJava(
                            name='somatic_gatk_cat_variants',
                            sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                            java_temporary_path=runnable_gather.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx4G')
                        runnable_gather.add_runnable_step(runnable_step=_runnable_step)

                        _runnable_step.add_option_short(key='classpath', value=self.java_archive_gatk)
                        _sub_command = _runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        _sub_command.add_option_long(key='reference', value=reference_gather)
                        _sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_somatic_gather.somatic_vcf)
                        _sub_command.add_switch_long(key='assumeSorted')
                        # Finally, add bsf.process.RunnableStep options, obsolete files and
                        # bsf.process.Executable dependencies per chunk index.
                        for chunk_index in chunk_index_list:
                            runnable_object = runnable_gather_list[chunk_index]
                            file_path_somatic_gather = FilePathSomaticScatterGather(prefix=runnable_object.name)
                            # Set GATK option variant
                            _sub_command.add_option_long(
                                key='variant',
                                value=file_path_somatic_gather.somatic_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            _runnable_step.obsolete_file_path_list.append(file_path_somatic_gather.somatic_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            _runnable_step.obsolete_file_path_list.append(file_path_somatic_gather.somatic_tbi)
                            # Set dependencies on preceding bsf.procedure.Runnable.name or
                            # bsf.process.Executable.name objects.
                            executable_gather.dependencies.append(runnable_object.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    runnable_gather_list = temporary_runnable_gather_list
                    tile_index_list = temporary_tile_index_list
                    gather_level += 1

            # For the last gather bsf.procedure.Runnable, move (rename) file paths to the top-level prefix and
            # adjust the FilePath object accordingly.
            file_path_somatic_gather = FilePathSomaticScatterGather(prefix=runnable_gather.name)

            _runnable_step = RunnableStepMove(
                name='somatic_gather_move_vcf',
                source_path=file_path_somatic_gather.somatic_vcf,
                target_path=file_path_somatic_final.somatic_vcf)
            runnable_gather.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStepMove(
                name='somatic_gather_move_tbi',
                source_path=file_path_somatic_gather.somatic_tbi,
                target_path=file_path_somatic_final.somatic_tbi)
            runnable_gather.add_runnable_step(runnable_step=_runnable_step)

            file_path_somatic_gather.somatic_vcf = file_path_somatic_final.somatic_vcf
            file_path_somatic_gather.somatic_tbi = file_path_somatic_final.somatic_tbi

            return runnable_gather

        def run_annotate_snpeff(prefix, vcf_file_path):
            """Private function to annotate a VCF file via the snpEff tool.

            This function is used in both, cohort and somatic variant calling annotation.
            @param prefix: Prefix
            @type prefix: str
            @param vcf_file_path: VCF file path
            @type vcf_file_path: str
            @return: C{bsf.procedure.Runnable}
            @rtype: Runnable
            """
            # snpEff                  (snpeff)
            # Bgzip                   (snpeff_bgzip)
            # Tabix                   (snpeff_tabix)
            # GATK VariantAnnotator   (gatk_variant_annotator)

            prefix_annotate = prefix

            file_path_annotate = FilePathAnnotateSnpEff(prefix=prefix_annotate)

            runnable_annotate = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=prefix_annotate,
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))

            if use_cache:
                reference_annotate = runnable_annotate.get_cache_file_path(
                    file_path=self.bwa_genome_db,
                    absolute=True)
            else:
                reference_annotate = self.bwa_genome_db

            # Run the snpEff tool for functional variant annotation in VCF output mode.

            _runnable_step = RunnableStepJava(
                name='snpeff_complete',
                stdout=ConnectorFile(file_path=file_path_annotate.complete_vcf, file_mode='wt'),
                java_temporary_path=runnable_annotate.temporary_directory_path(absolute=False),
                java_jar_path=self.java_archive_snpeff,
                java_heap_maximum='Xmx6G')
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            # Use a sequence of sub-Command objects to separate options that have to appear
            # in a particular order. Sigh!
            _runnable_step.sub_command.sub_command = Command(program='ann')
            _sub_command = _runnable_step.sub_command.sub_command
            _sub_command.add_switch_short(key='download')
            _sub_command.add_option_short(key='o', value='vcf')
            _sub_command.add_option_short(key='stats', value=file_path_annotate.complete_stats)
            _sub_command.add_option_short(
                key='config',
                value=os.path.join(os.path.dirname(self.java_archive_snpeff), 'snpEff.config'))

            _sub_command.arguments.append(self.snpeff_genome_version)
            _sub_command.arguments.append(vcf_file_path)

            # Compress and index the snpEff VCF file with bgzip and tabix, respectively.

            _runnable_step = RunnableStep(
                name='bgzip_complete',
                program='bgzip',
                arguments=[file_path_annotate.complete_vcf])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStep(
                name='tabix_complete',
                program='tabix',
                arguments=[file_path_annotate.complete_vcf_bgz])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            # Run the snpEff tool for functional variant annotation in GATK output mode.

            _runnable_step = RunnableStepJava(
                name='snpeff_gatk',
                stdout=ConnectorFile(file_path=file_path_annotate.gatk_vcf, file_mode='wt'),
                java_temporary_path=runnable_annotate.temporary_directory_path(absolute=False),
                java_jar_path=self.java_archive_snpeff,
                java_heap_maximum='Xmx6G')
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            # Use a sequence of sub-Command objects to separate options that have to appear
            # in a particular order. Sigh!
            _runnable_step.sub_command.sub_command = Command(program='ann')
            _sub_command = _runnable_step.sub_command.sub_command
            _sub_command.add_switch_short(key='download')
            _sub_command.add_option_short(key='o', value='gatk')
            _sub_command.add_option_short(key='stats', value=file_path_annotate.gatk_stats)
            _sub_command.add_option_short(
                key='config',
                value=os.path.join(os.path.dirname(self.java_archive_snpeff), 'snpEff.config'))

            _sub_command.arguments.append(self.snpeff_genome_version)
            _sub_command.arguments.append(vcf_file_path)

            # Compress and index the snpEff VCF file with bgzip and tabix, respectively.

            _runnable_step = RunnableStep(
                name='bgzip_gatk',
                program='bgzip',
                arguments=[file_path_annotate.gatk_vcf])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStep(
                name='tabix_gatk',
                program='tabix',
                arguments=[file_path_annotate.gatk_vcf_bgz])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step.add_option_long(key='preset', value='vcf')

            # Run the GATK VariantAnnotator analysis.

            _runnable_step = RunnableStepGATK(
                name='gatk_variant_annotator',
                java_temporary_path=runnable_annotate.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx4G',
                java_jar_path=self.java_archive_gatk)
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step.add_gatk_option(key='analysis_type', value='VariantAnnotator')
            _runnable_step.add_gatk_option(key='reference_sequence', value=reference_annotate)
            if self.known_sites_discovery:
                _runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)

            # Add annotation resources and their corresponding expression options.
            if self.annotation_resources_dict:  # not None and not empty
                for _resource_name, (_file_path, _annotation_list) in self.annotation_resources_dict.items():
                    if _file_path and _annotation_list:
                        _runnable_step.add_gatk_option(
                            key=':'.join(('resource', _resource_name)),
                            value=_file_path)
                        for _annotation in _annotation_list:
                            _runnable_step.add_gatk_option(
                                key='expression',
                                value='.'.join((_resource_name, _annotation)),
                                override=True)

            _runnable_step.add_gatk_option(key='variant', value=vcf_file_path)
            # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
            # _runnable_step.add_gatk_option(key='annotation', value='AlleleBalanceBySample')
            _runnable_step.add_gatk_option(key='annotation', value='SnpEff')
            _runnable_step.add_gatk_option(key='snpEffFile', value=file_path_annotate.gatk_vcf_bgz)
            _runnable_step.add_gatk_option(key='out', value=file_path_annotate.annotated_vcf)

            return runnable_annotate

        def run_annotate_vep(prefix, vcf_file_path):
            """Private function to annotate a VCF file via the Ensembl Variant Effect Predictor (VEP).

            This function is used in both, cohort and somatic variant calling annotation.
            @param prefix: Prefix
            @type prefix: str
            @param vcf_file_path: VCF file path
            @type vcf_file_path: str
            @return: C{bsf.procedure.Runnable}
            @rtype: Runnable
            """
            # Ensembl Variant Effect Predictor (ensembl_vep)
            # Bgzip                            (ensembl_vep_bgzip)
            # Tabix                            (ensembl_vep_tabix)
            # Ensembl Variant Effect Filter    (ensembl_filter)
            # Bgzip                            (ensembl_filter_bgzip)
            # Tabix                            (ensembl_filter_tabix)
            # VCF.Filter                       (vcf_filter_complete)
            # VCF.Filter                       (vcf_filter_filtered)

            prefix_annotate = prefix

            file_path_annotate = FilePathAnnotateVEP(prefix=prefix_annotate)

            runnable_annotate = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=prefix_annotate,
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))

            # reference_annotate = runnable_annotate.get_absolute_cache_file_path(
            #     file_path=self.bwa_genome_db)

            # if not os.path.exists(os.path.join(self.genome_directory, file_path_annotate.complete_vcf_tbi)):
            # Run the Ensembl Variant Effect Predictor script.

            _runnable_step = RunnableStep(
                name='ensembl_vep',
                program='perl',
                sub_command=Command())
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            # FIXME: Distinguish between cohort and somatic annotation for RunnableStep configuration.
            # (section annotate_cohort_ensembl_vep vs annotate_somatic_ensembl_vep)
            # Read RunnableStep options from configuration sections:
            # [bsf.analyses.variant_calling.VariantCallingGATK.ensembl_vep]
            # self.set_runnable_step_configuration(runnable_step=_runnable_step)
            _runnable_step.arguments.append(os.path.join(self.vep_source, 'vep'))
            _sub_command = _runnable_step.sub_command
            _sub_command.add_option_long(key='fork', value=str(stage_annotate_cohort_vep.threads))
            # Basic options
            _sub_command.add_switch_long(key='everything')
            # Input options
            _sub_command.add_option_long(key='species', value=self.vep_species)
            _sub_command.add_option_long(key='assembly', value=self.vep_assembly)
            _sub_command.add_option_long(key='input_file', value=vcf_file_path)
            _sub_command.add_option_long(key='format', value='vcf')  # Input file format
            _sub_command.add_option_long(key='output_file', value=file_path_annotate.complete_raw_vcf)
            _sub_command.add_switch_long(key='force_overwrite')
            _sub_command.add_option_long(key='stats_file', value=file_path_annotate.statistics)
            # Cache options
            _sub_command.add_switch_long(key='cache')
            _sub_command.add_switch_long(key='offline')  # VEP e91 option
            _sub_command.add_option_long(key='dir_cache', value=self.vep_cache)
            _sub_command.add_option_long(key='dir_plugins', value=self.vep_plugin)
            _sub_command.add_option_long(key='fasta_dir', value=self.vep_fasta)  # VEP e91 option
            if self.vep_annotation in ('refseq', 'merged'):
                # Options --refseq  or --merged.
                _sub_command.add_switch_long(key=self.vep_annotation)
            _sub_command.add_option_long(key='bam', value=self.vep_refseq_alignments_path)
            # Other annotation sources
            if self.vep_plugin_cadd_path:
                _sub_command.add_option_long(key='plugin', value=','.join(('CADD', self.vep_plugin_cadd_path)))
            # Output options
            _sub_command.add_switch_long(key='allele_number')
            _sub_command.add_switch_long(key='no_escape')  # Do not percent escape HGVS strings
            # Identifiers
            _sub_command.add_switch_long(key='hgvsg')
            # Co-located variants
            _sub_command.add_option_long(key='failed', value='1')
            # Data format options
            _sub_command.add_switch_long(key='vcf')
            # Filtering and QC options
            if self.vep_annotation == 'ensembl':
                # Since the --gencode_basic option suppresses all RefSeq annotation,
                # it can only be set for Ensembl annotation.
                _sub_command.add_switch_long(key='gencode_basic')
            _sub_command.add_switch_long(key='dont_skip')
            _sub_command.add_switch_long(key='allow_non_variant')
            _sub_command.add_switch_long(key='flag_pick_allele_gene')
            # Database options
            if self.vep_sql_user:
                _sub_command.add_option_long(key='user', value=self.vep_sql_user)
            if self.vep_sql_pass:
                _sub_command.add_option_long(key='password', value=self.vep_sql_pass)
            if self.vep_sql_host:
                _sub_command.add_option_long(key='host', value=self.vep_sql_host)
            if self.vep_sql_port:
                _sub_command.add_option_long(key='port', value=self.vep_sql_port)
            # Undocumented options
            _sub_command.add_switch_long(key='no_progress')
            _sub_command.add_option_long(
                key='tmpdir',
                value=runnable_annotate.temporary_directory_path(absolute=False))

            _runnable_step = RunnableStep(
                name='ensembl_vep_bgzip',
                program='bgzip',
                arguments=[file_path_annotate.complete_raw_vcf])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStep(
                name='ensembl_vep_tabix',
                program='tabix',
                arguments=[file_path_annotate.complete_raw_vcf_bgz])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step.add_option_long(key='preset', value='vcf')

            # if not os.path.exists(os.path.join(self.genome_directory, file_path_annotate.vep_filtered_vcf_tbi)):
            # Run the Ensembl Variant Effect Filter script.

            _runnable_step = RunnableStep(
                name='ensembl_filter',
                program='perl',
                sub_command=Command())
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            # FIXME: Distinguish between cohort and somatic annotation for RunnableStep configuration.
            # (section annotate_cohort_ensembl_filter vs annotate_somatic_ensembl_filter)
            # Read RunnableStep options from configuration sections:
            # [bsf.analyses.variant_calling.VariantCallingGATK.ensembl_filter]
            # self.set_runnable_step_configuration(runnable_step=_runnable_step)
            _runnable_step.arguments.append(os.path.join(self.vep_source, 'filter_vep'))
            _sub_command = _runnable_step.sub_command
            _sub_command.add_option_long(key='input_file', value=file_path_annotate.complete_raw_vcf_bgz)
            _sub_command.add_option_long(key='format', value='vcf')
            _sub_command.add_option_long(key='output_file', value=file_path_annotate.filtered_raw_vcf)
            _sub_command.add_switch_long(key='only_matched')
            _sub_command.add_option_long(key='filter', value='Consequence ne upstream_gene_variant', override=True)
            _sub_command.add_option_long(key='filter', value='Consequence ne downstream_gene_variant', override=True)
            _sub_command.add_option_long(key='filter', value='Consequence ne intron_variant', override=True)
            _sub_command.add_option_long(key='filter', value='BIOTYPE ne processed_transcript', override=True)
            # _sub_command.add_option_long(key='filter', value='CANONICAL eq YES', override=True)
            _sub_command.add_switch_long(key='force_overwrite')

            _runnable_step = RunnableStep(
                name='ensembl_filter_bgzip',
                program='bgzip',
                arguments=[file_path_annotate.filtered_raw_vcf])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStep(
                name='ensembl_filter_tabix',
                program='tabix',
                arguments=[file_path_annotate.filtered_raw_vcf_bgz])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step.add_option_long(key='preset', value='vcf')

            # Run the VCF Filter on the complete VEP set to convert (re-model) the CSQ field into
            # a set of independent INFO fields.

            # Convert the CSQ INFO annotation into a set of INFO VEP_* annotation
            # that is more accessible to down-stream tools.

            _runnable_step = RunnableStepCsqToVep(
                name='ensembl_complete_csq_to_vep',
                soc_path=self.vep_soc_path,
                ofc_path=self.vep_ofc_path,
                vcf_path_old=file_path_annotate.complete_raw_vcf_bgz,
                vcf_path_new=file_path_annotate.complete_vcf_bgz)
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStepCsqToVep(
                name='ensembl_filtered_csq_to_vep',
                soc_path=self.vep_soc_path,
                ofc_path=self.vep_ofc_path,
                vcf_path_old=file_path_annotate.filtered_raw_vcf_bgz,
                vcf_path_new=file_path_annotate.filtered_vcf_bgz)
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step = RunnableStep(
                name='ensembl_complete_csq_to_vep_tabix',
                program='tabix',
                arguments=[file_path_annotate.complete_vcf_bgz])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step.add_option_long(key='preset', value='vcf')

            _runnable_step = RunnableStep(
                name='ensembl_filtered_csq_to_vep_tabix',
                program='tabix',
                arguments=[file_path_annotate.filtered_vcf_bgz])
            runnable_annotate.add_runnable_step(runnable_step=_runnable_step)

            _runnable_step.add_option_long(key='preset', value='vcf')

            return runnable_annotate

        # Start of the run() method body.

        super(VariantCallingGATK, self).run()

        # Get global defaults.

        # VariantCallingGATK requires a genome version, which gets configured by the super-class.

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a 'genome_version' configuration option.")

        if not self.bwa_genome_db:
            raise Exception('A ' + self.name + " requires a 'bwa_genome_db' configuration option.")

        if not self.cohort_name:
            self.cohort_name = self.project_name  # The cohort_name used to default to just 'default'.

        if not self.gatk_bundle_version:
            raise Exception('A ' + self.name + " requires a 'gatk_bundle_version' configuration option.")

        if not self.snpeff_genome_version:
            raise Exception('A ' + self.name + " requires a 'snpeff_genome_version' configuration option.")

        if not self.vep_annotation:
            self.vep_annotation = 'ensembl'

        if self.vep_annotation not in ('ensembl', 'refseq', 'merged'):
            raise Exception(
                'The ' + self.name + " option 'vep_annotation' has to be 'ensembl', 'refseq' or 'merged', " +
                'not ' + repr(self.vep_annotation) + '.')

        if not self.vep_assembly:
            self.vep_assembly = EnsemblVEP.get_name_assembly(genome_version=self.genome_version)
            if not self.vep_assembly:
                raise Exception('A ' + self.name + " requires a 'vep_assembly' configuration option.")

        if not self.vep_cache:
            self.vep_cache = EnsemblVEP.get_directory_cache(genome_version=self.genome_version)
            if not self.vep_cache:
                raise Exception('A ' + self.name + " requires a 'vep_cache' configuration option.")

        if not self.vep_fasta:
            self.vep_fasta = EnsemblVEP.get_directory_fasta(genome_version=self.genome_version)
            if not self.vep_fasta:
                raise Exception('A ' + self.name + " requires a 'vep_fasta' configuration option.")

        if not self.vep_plugin:
            self.vep_plugin = EnsemblVEP.get_directory_plugin(genome_version=self.genome_version)
            if not self.vep_plugin:
                raise Exception('A ' + self.name + " requires a 'vep_plugin' configuration option.")

        if not self.vep_source:
            self.vep_source = EnsemblVEP.get_directory_source(genome_version=self.genome_version)
            if not self.vep_source:
                raise Exception('A ' + self.name + " requires a 'vep_source' configuration option.")

        if not self.vep_species:
            self.vep_species = EnsemblVEP.get_name_species(genome_version=self.genome_version)
            if not self.vep_species:
                raise Exception('A ' + self.name + " requires a 'vep_species' configuration option.")

        if not self.vep_sql_user:
            self.vep_sql_user = EnsemblVEP.get_sql_user(genome_version=self.genome_version)

        if not self.vep_sql_pass:
            self.vep_sql_pass = EnsemblVEP.get_sql_pass(genome_version=self.genome_version)

        if not self.vep_sql_host:
            self.vep_sql_host = EnsemblVEP.get_sql_host(genome_version=self.genome_version)

        if not self.vep_sql_port:
            self.vep_sql_port = EnsemblVEP.get_sql_port(genome_version=self.genome_version)

        if not self.vep_ofc_path:
            self.vep_ofc_path = EnsemblVEP.get_ofc_path(genome_version=self.genome_version)
            if not self.vep_ofc_path:
                raise Exception('A ' + self.name + " requires a 'vep_ofc_path' configuration option.")

        if not self.vep_soc_path:
            self.vep_soc_path = EnsemblVEP.get_soc_path(genome_version=self.genome_version)
            if not self.vep_soc_path:
                raise Exception('A ' + self.name + " requires a 'vep_soc_path' configuration option.")

        if not self.vep_refseq_alignments_path:
            self.vep_refseq_alignments_path = EnsemblVEP.get_refseq_alignments_path(
                genome_version=self.genome_version)
            if not self.vep_refseq_alignments_path:
                raise Exception('A ' + self.name + " requires a 'vep_refseq_alignments_path' configuration option.")

        if not self.vep_plugin_cadd_path:
            self.vep_plugin_cadd_path = EnsemblVEP.get_cadd_path(genome_version=self.genome_version)
        if self.vep_plugin_cadd_path and not os.path.isabs(self.vep_plugin_cadd_path):
            cadd_resource_path = StandardFilePath.get_resource_cadd(absolute=True)
            if cadd_resource_path:
                self.vep_plugin_cadd_path = os.path.join(cadd_resource_path, self.vep_plugin_cadd_path)

        if not self.java_archive_fgbio:
            self.java_archive_fgbio = JavaArchive.get_fgbio()
            if not self.java_archive_fgbio:
                raise Exception('A ' + self.name + " requires a 'java_archive_fgbio' configuration option.")

        if not self.java_archive_gatk:
            self.java_archive_gatk = JavaArchive.get_gatk()
            if not self.java_archive_gatk:
                raise Exception('A ' + self.name + " requires a 'java_archive_gatk' configuration option.")

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception('A ' + self.name + " requires a 'java_archive_picard' configuration option.")

        if not self.java_archive_snpeff:
            self.java_archive_snpeff = JavaArchive.get_snpeff()
            if not self.java_archive_snpeff:
                raise Exception('A ' + self.name + " requires a 'java_archive_snpeff' configuration option.")

        if not self.java_archive_vcf_filter:
            self.java_archive_vcf_filter = JavaArchive.get_vcf_filter()
            if not self.java_archive_vcf_filter:
                raise Exception('A ' + self.name + " requires a 'java_archive_vcf_filter' configuration option.")

        # Check for absolute paths and adjust if required before checking for existence.

        self.bwa_genome_db = self.configuration.get_absolute_path(
            file_path=self.bwa_genome_db,
            default_path=self.get_gatk_bundle_path)
        if not os.path.exists(self.bwa_genome_db):
            raise Exception('The file path ' + repr(self.bwa_genome_db) +
                            " in option 'bwa_genome_db' does not exist.")

        # GATK does a lot of read requests from the reference FASTA file.
        # Place it and the accompanying *.fasta.fai and *.dict files in the cache directory.
        if use_cache:
            self._cache_path_dict = {
                'reference_fasta': self.bwa_genome_db,
                'reference_fai': self.bwa_genome_db + '.fai',
                'reference_dict': os.path.splitext(self.bwa_genome_db)[0] + '.dict'
            }

        # List of accessory cohort GVCF file paths

        if self.accessory_cohort_gvcfs:  # not None and not empty
            for i, file_path in enumerate(self.accessory_cohort_gvcfs):
                file_path = self.configuration.get_absolute_path(
                    file_path=file_path,
                    default_path=StandardFilePath.get_projects(absolute=True))
                if os.path.exists(file_path):
                    self.accessory_cohort_gvcfs[i] = file_path
                else:
                    raise Exception('The file path ' + repr(file_path) +
                                    " in option 'accessory_cohort_gvcf' does not exist.")
                    # TODO: Check the cohorts so that their sample names do not clash.

        # Dict of annotation resources

        if self.annotation_resources_dict:  # not None and not empty
            for resource_name, (file_path, annotation_list) in self.annotation_resources_dict.items():
                file_path = self.configuration.get_absolute_path(
                    file_path=file_path,
                    default_path=self.get_gatk_bundle_path)
                if os.path.exists(file_path):
                    self.annotation_resources_dict[resource_name] = file_path, annotation_list
                else:
                    raise Exception('The file path ' + repr(file_path) +
                                    ' in annotation resource ' + repr(resource_name) + ' does not exist.')

        # Known sites for cohort variant discovery

        if self.known_sites_discovery:
            self.known_sites_discovery = self.configuration.get_absolute_path(
                file_path=self.known_sites_discovery,
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(self.known_sites_discovery):
                raise Exception('The file path ' + repr(self.known_sites_discovery) +
                                " in option 'known_sites_discovery' does not exist.")

        # List of known sites for InDel re-alignment

        if self.known_sites_realignment:  # not None and not empty
            for i, file_path in enumerate(self.known_sites_realignment):
                file_path = self.configuration.get_absolute_path(
                    file_path=file_path,
                    default_path=self.get_gatk_bundle_path)
                if os.path.exists(file_path):
                    self.known_sites_realignment[i] = file_path
                else:
                    raise Exception('The file path ' + repr(file_path) +
                                    " in option 'known_sites_realignment' does not exist.")

        # List of known sites for BQSR

        if self.known_sites_recalibration:  # not None and not empty
            for i, file_path in enumerate(self.known_sites_recalibration):
                file_path = self.configuration.get_absolute_path(
                    file_path=file_path,
                    default_path=self.get_gatk_bundle_path)
                if os.path.exists(file_path):
                    self.known_sites_recalibration[i] = file_path
                else:
                    raise Exception('The file path ' + repr(file_path) +
                                    " in option 'known_sites_recalibration' does not exist.")

        # List of known Catalogue Of Somatic Mutations In Cancer (COSMIC) sites for somatic variant calling

        if self.known_somatic_discovery:
            for i, file_path in enumerate(self.known_somatic_discovery):
                file_path = self.configuration.get_absolute_path(
                    file_path=file_path,
                    default_path=StandardFilePath.get_resource_cosmic(absolute=True))
                if os.path.exists(file_path):
                    self.known_somatic_discovery[i] = file_path
                else:
                    raise Exception('The file path ' + repr(file_path) +
                                    " in option 'known_somatic_discovery' does not exist.")

        # Dict of VQSR InDel resources

        if self.vqsr_resources_indel_dict:  # not None and not empty
            for resource_name, resource_dict in self.vqsr_resources_indel_dict.items():
                resource_dict['file_path'] = self.configuration.get_absolute_path(
                    file_path=resource_dict['file_path'],
                    default_path=self.get_gatk_bundle_path)
                if not os.path.exists(resource_dict['file_path']):
                    raise Exception('The file path ' + repr(resource_dict['file_path']) +
                                    " in option 'vqsr_resources_indel' " + repr(resource_name) + ' does not exist.')

        # Dict of VQSR SNP resources

        if self.vqsr_resources_snp_dict:  # not None and not empty
            for resource_name, resource_dict in self.vqsr_resources_snp_dict.items():
                resource_dict['file_path'] = self.configuration.get_absolute_path(
                    file_path=resource_dict['file_path'],
                    default_path=self.get_gatk_bundle_path)
                if not os.path.exists(resource_dict['file_path']):
                    raise Exception('The file path ' + repr(resource_dict['file_path']) +
                                    " in option 'vqsr_resources_snp' " + repr(resource_name) + ' does not exist.')

        # List of excluded intervals

        if self.exclude_intervals_list:  # not None and not empty
            for i, interval in enumerate(self.exclude_intervals_list):
                if interval.endswith('.intervals') or interval.endswith('.interval_list'):
                    # For Picard-style interval lists prepend the current directory if necessary.
                    if not os.path.isabs(interval):
                        interval = os.path.join(os.path.realpath(os.path.curdir), interval)
                    if os.path.exists(interval):
                        self.exclude_intervals_list[i] = interval
                    else:
                        raise Exception('The file path ' + repr(interval) +
                                        " in option 'exclude_intervals' does not exist.")

        # List of included intervals

        if self.include_intervals_list:  # not None and not empty
            for i, interval in enumerate(self.include_intervals_list):
                if interval.endswith('.intervals') or interval.endswith('.interval_list'):
                    # For Picard-style interval lists prepend the current directory if necessary.
                    if not os.path.isabs(interval):
                        interval = os.path.join(os.path.realpath(os.path.curdir), interval)
                    if os.path.exists(interval):
                        self.include_intervals_list[i] = interval
                    else:
                        raise Exception('The file path ' + repr(interval) +
                                        " in option 'include_intervals' does not exist.")

        # Genome Annotation GTF file path, defaults to the interval files directory.

        if self.genome_annotation_gtf and not os.path.isabs(self.genome_annotation_gtf):
            self.genome_annotation_gtf = self.configuration.get_absolute_path(
                file_path=self.genome_annotation_gtf,
                default_path=StandardFilePath.get_resource_intervals(absolute=True))
            # TODO: Use the transcriptome directory as the default location.
            # Create a new, genome-specific configuration object in the bsf.standards module.

        # Read comparisons for somatic mutation calling.
        run_read_comparisons()

        # Create genomic tiles for scatter gather approaches.
        # In case the number of tiles is not defined or 0, a Python list with a single Container with an empty Interval
        # gets returned by both functions.

        if self.scatter_intervals_path:
            if not os.path.exists(self.scatter_intervals_path):
                raise Exception('Picard ScatterIntervalsByNs interval file ' + repr(self.scatter_intervals_path) +
                                ' does not exist.')

            self._tile_region_cohort_list = get_interval_tiles(
                interval_path=self.scatter_intervals_path,
                tile_number=self.number_of_tiles_cohort)

            self._tile_region_somatic_list = get_interval_tiles(
                interval_path=self.scatter_intervals_path,
                tile_number=self.number_of_tiles_somatic)
        else:
            dictionary_path = os.path.splitext(self.bwa_genome_db)[0] + '.dict'

            if not os.path.exists(dictionary_path):
                raise Exception('Picard sequence dictionary ' + repr(dictionary_path) + ' does not exist.')

            self._tile_region_cohort_list = get_genome_tiles(
                dictionary_path=dictionary_path,
                tile_number=self.number_of_tiles_cohort)

            self._tile_region_somatic_list = get_genome_tiles(
                dictionary_path=dictionary_path,
                tile_number=self.number_of_tiles_somatic)

        stage_align_lane = self.get_stage(name=self.get_stage_name_align_lane())
        stage_process_lane = self.get_stage(name=self.get_stage_name_process_lane())
        stage_process_sample = self.get_stage(name=self.get_stage_name_process_sample())
        stage_diagnose_sample = self.get_stage(name=self.get_stage_name_diagnose_sample())
        stage_merge_cohort = self.get_stage(name=self.get_stage_name_merge_cohort())
        stage_process_cohort = self.get_stage(name=self.get_stage_name_process_cohort())
        stage_annotate_cohort_snpeff = self.get_stage(name=self.get_stage_name_annotate_cohort_snpeff())
        stage_annotate_cohort_vep = self.get_stage(name=self.get_stage_name_annotate_cohort_vep())
        stage_split_cohort_snpeff = self.get_stage(name=self.get_stage_name_split_cohort_snpeff())
        stage_split_cohort_vep = self.get_stage(name=self.get_stage_name_split_cohort_vep())
        stage_summary = self.get_stage(name=self.get_stage_name_summary())
        stage_somatic = self.get_stage(name=self.get_stage_name_somatic())
        stage_annotate_somatic_snpeff = self.get_stage(name=self.get_stage_name_annotate_somatic_snpeff())
        stage_annotate_somatic_vep = self.get_stage(name=self.get_stage_name_annotate_somatic_vep())
        stage_split_somatic_snpeff = self.get_stage(name=self.get_stage_name_split_somatic_snpeff())
        stage_split_somatic_vep = self.get_stage(name=self.get_stage_name_split_somatic_vep())

        # Create a Python dict of Python str (cohort name) key and Python list of
        # process_sample bsf.procedure.Runnable object value data.
        # This dictionary is required by the merge_cohort stage to hierarchically merge cohorts.

        runnable_process_sample_dict = dict()
        """ @type runnable_process_sample_dict: dict[str, list[(Runnable, str)]] """

        # Create a Python list of diagnose_sample bsf.procedure.Runnable objects.

        runnable_diagnose_sample_list = list()
        """ @type runnable_diagnose_sample_list: list[Runnable] """

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

            runnable_process_read_group_list = list()
            """ @type runnable_process_read_group_list: list[Runnable] """

            for paired_reads_name in sorted(paired_reads_dict):
                if not paired_reads_dict[paired_reads_name]:
                    # Skip names, which PairedReads objects have all been excluded.
                    continue

                #################################
                # Step 1: Align per read group. #
                #################################
                #
                # bsf_run_bwa.py
                # - Picard SamToFastq
                # - BWA MEM

                runnable_step = RunnableStep(
                    name='align_lane_bwa',
                    program='bwa',
                    sub_command=Command(program='mem'))
                # Instead of adding the bsf.process.RunnableStep to the bsf.procedure.Runnable,
                # it gets serialised into the pickler_file.

                # Read configuration sections
                # [bsf.analyses.variant_calling.VariantCallingGATK.align_lane_bwa]
                # [bsf.analyses.variant_calling.VariantCallingGATK.align_lane_bwa.mem]
                self.set_runnable_step_configuration(runnable_step=runnable_step)

                bwa_mem = runnable_step.sub_command

                # Set BWA mem options.

                # Allow as many threads as defined in the corresponding Stage.
                bwa_mem.add_option_short(key='t', value=str(stage_align_lane.threads))
                # Append FASTA/Q comment to SAM output.
                # Illumina-style FASTQ headers obey the following schema, which is not SAM compliant.
                # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> \
                # <read>:<is filtered>:<control number>:<sample number>
                # bwa_mem.add_switch_short(key='C')
                # Mark shorter split hits as secondary (for Picard compatibility).
                bwa_mem.add_switch_short(key='M')
                # Output errors only.
                bwa_mem.add_option_short(key='v', value='1')

                # Set BWA arguments.

                bwa_mem.arguments.append(self.bwa_genome_db)

                reads1 = list()
                reads2 = list()

                # Propagate the SAM read group information around FASTQ files if required.
                # Please note that only the first read group can be propagated per
                # PairedReads object.

                read_group = str()

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if paired_reads.reads_1 is not None:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2 is not None:
                        reads2.append(paired_reads.reads_2.file_path)
                    if not read_group and paired_reads.read_group:
                        read_group = paired_reads.read_group

                if read_group:
                    bwa_mem.add_option_short(key='R', value=read_group)

                if len(reads1) and not len(reads2):
                    bwa_mem.arguments.append(','.join(reads1))
                elif len(reads1) and len(reads2):
                    bwa_mem.arguments.append(','.join(reads1))
                    bwa_mem.arguments.append(','.join(reads2))
                elif not len(reads1) and len(reads2):
                    warnings.warn('Only second reads, but no first reads have been defined.')
                else:
                    warnings.warn('No reads have been defined.')

                # TODO: The name for the aligned BAM is constructed by the bsf_run_bwa.py script.
                # It is currently based on the stage_align_lane.name and paired_reads_name.
                # The script should also be changed to pre-set all file names beforehand.
                file_path_alignment = self.get_file_path_align_lane(paired_reads_name=paired_reads_name)

                # Normally, the bwa object would be pushed onto the Stage list.
                # Experimentally, use Pickler to serialize the bsf.process.Executable into a file.

                pickler_dict_align_lane = {
                    'prefix': stage_align_lane.name,
                    'replicate_key': paired_reads_name,
                    'java_archive_gatk': self.java_archive_gatk,
                    'java_archive_picard': self.java_archive_picard,
                    'runnable_step': runnable_step,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    stage_align_lane.name + '_' + paired_reads_name + '.pkl')
                with open(file=pickler_path, mode='wb') as pickler_file:
                    pickler = pickle.Pickler(file=pickler_file, protocol=pickle.HIGHEST_PROTOCOL)
                    pickler.dump(pickler_dict_align_lane)

                # Create a bsf_run_bwa.py job to run the pickled object.

                executable_align_lane = stage_align_lane.add_executable(
                    executable=Executable(
                        name='_'.join((stage_align_lane.name, paired_reads_name)),
                        program='bsf_run_bwa.py'))

                # Only submit this bsf.process.Executable if the final result file does not exist.
                if (os.path.exists(os.path.join(self.genome_directory, file_path_alignment.aligned_md5)) and
                        os.path.getsize(os.path.join(self.genome_directory, file_path_alignment.aligned_md5))):
                    executable_align_lane.submit = False
                # Check also for existence of a new-style bsf.procedure.Runnable status file.
                if os.path.exists(os.path.join(
                        stage_align_lane.working_directory,
                        '_'.join((stage_align_lane.name, paired_reads_name, 'completed.txt')))):
                    executable_align_lane.submit = False

                # Set executable_align_lane options.
                executable_align_lane.add_option_long(key='pickler_path', value=pickler_path)
                executable_align_lane.add_option_long(key='debug', value=str(self.debug))

                ###################################
                # Step 2: Process per read group. #
                ###################################
                #
                # Picard MarkDuplicates                 (process_lane_picard_mark_duplicates)
                # GATK RealignerTargetCreator           (process_lane_gatk_realigner_target_creator)
                # GATK IndelRealigner                   (process_lane_gatk_indel_realigner)
                # GATK BaseRecalibrator                 (process_lane_gatk_base_recalibrator_pre)
                # GATK BaseRecalibrator on-the-fly      (process_lane_gatk_base_recalibrator_post)
                # GATK AnalyzeCovariates                (process_lane_gatk_analyze_covariates)
                # GATK PrintReads                       (process_lane_gatk_print_reads)
                # Picard CollectAlignmentSummaryMetrics (process_lane_picard_collect_alignment_summary_metrics)

                # Lane-specific file paths

                file_path_process_read_group = self.get_file_path_process_read_group(
                    paired_reads_name=paired_reads_name)

                # Create a bsf.procedure.Runnable and bsf.process.Executable for processing each read group.

                runnable_process_lane = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_process_lane(paired_reads_name=paired_reads_name),
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        debug=self.debug))
                executable_process_lane = self.set_stage_runnable(
                    stage=stage_process_lane,
                    runnable=runnable_process_lane)
                # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
                executable_process_lane.dependencies.append(executable_align_lane.name)
                # Set dependencies for succeeding bsf.procedure.Runnable or bsf.process.Executable objects.
                runnable_process_read_group_list.append(runnable_process_lane)

                if use_cache:
                    reference_process_lane = runnable_process_lane.get_cache_file_path(
                        file_path=self.bwa_genome_db,
                        absolute=True)
                else:
                    reference_process_lane = self.bwa_genome_db

                # Run fgbio TrimPrimer if requested.

                amplicons = VariantCallingGATKAmplicons.from_sample(sample=sample)
                if amplicons.amplicons_path:
                    runnable_step = RunnableStepJava(
                        name='trim_primers',
                        java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx4G',
                        java_jar_path=self.java_archive_fgbio)
                    runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                    # Use a sequence of sub-Command objects to separate options that have to appear
                    # in a particular order. Sigh!

                    # Sub-command to separate the fgbio Java archive from general options
                    runnable_step.sub_command.sub_command = Command()
                    # async-io [false]
                    runnable_step.sub_command.sub_command.add_option_pair_long(key='async-io', value='true')
                    # compression [5]
                    # runnable_step.sub_command.sub_command.add_option_long(key='compression', value='9')
                    runnable_step.sub_command.sub_command.add_option_pair_long(
                        key='tmp-dir',
                        value=runnable_process_lane.temporary_directory_path(absolute=False))
                    # log-level [Info]
                    # sam-validation-stringency []

                    # fgbio 'TrimPrimer' command
                    runnable_step.sub_command.sub_command.sub_command = Command(program='TrimPrimers')
                    sub_command = runnable_step.sub_command.sub_command.sub_command
                    sub_command.add_option_pair_long(
                        key='input',
                        value=file_path_alignment.aligned_bam)
                    sub_command.add_option_pair_long(
                        key='output',
                        value=file_path_process_read_group.trimmed_bam)
                    sub_command.add_option_pair_long(key='primers', value=amplicons.amplicons_path)
                    sub_command.add_option_pair_long(key='hard-clip', value='true')
                    sub_command.add_option_pair_long(key='ref', value=reference_process_lane)
                    # slop [5]
                    # sort-order [input sort order]
                    # auto-trim-attributes [false]
                    # NOTE: This option trims all attributes, even the RG (read group)
                else:
                    file_path_process_read_group.trimmed_bam = file_path_alignment.aligned_bam
                    file_path_process_read_group.trimmed_bai = file_path_alignment.aligned_bai
                    file_path_process_read_group.trimmed_md5 = file_path_alignment.aligned_md5

                # Run the Picard MarkDuplicates analysis, unless configured to skip it.

                if self.skip_mark_duplicates:
                    file_path_process_read_group.duplicates_marked_bam = file_path_process_read_group.trimmed_bam
                    file_path_process_read_group.duplicates_marked_bai = file_path_process_read_group.trimmed_bai
                    file_path_process_read_group.duplicates_marked_md5 = file_path_process_read_group.trimmed_md5
                else:
                    # Run the Picard MarkDuplicates analysis.

                    runnable_step = RunnableStepPicard(
                        name='process_lane_picard_mark_duplicates',
                        obsolete_file_path_list=[
                            file_path_process_read_group.trimmed_bam,
                            file_path_process_read_group.trimmed_bai,
                            file_path_process_read_group.trimmed_md5,
                        ],
                        java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx4G',
                        java_jar_path=self.java_archive_picard,
                        picard_command='MarkDuplicates')
                    runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_picard_option(key='INPUT', value=file_path_process_read_group.trimmed_bam)
                    runnable_step.add_picard_option(
                        key='OUTPUT',
                        value=file_path_process_read_group.duplicates_marked_bam)
                    runnable_step.add_picard_option(
                        key='METRICS_FILE',
                        value=file_path_process_read_group.duplicate_metrics)
                    # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
                    # as otherwise, optical duplicates could not be detected. This is a consequence of using
                    # Illumina2bam rather than Picard ExtractIlluminaBarcodes, IlluminaBasecallsToFastq and
                    # IlluminaBasecallsToSam.
                    # See BioStar post: http://www.biostars.org/p/12538/
                    # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    runnable_step.add_picard_option(
                        key='READ_NAME_REGEX',
                        value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*')
                    runnable_step.add_picard_option(
                        key='OPTICAL_DUPLICATE_PIXEL_DISTANCE',
                        value='5000')
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_process_lane.temporary_directory_path(absolute=False))
                    # VERBOSITY defaults to 'INFO'.
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    # QUIET defaults to 'false'.
                    # VALIDATION_STRINGENCY defaults to 'STRICT'.
                    # COMPRESSION_LEVEL defaults to '5'.
                    runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                    # MAX_RECORDS_IN_RAM defaults to '500000'.
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='4000000')
                    # CREATE_INDEX defaults to 'false'.
                    runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                    # CREATE_MD5_FILE defaults to 'false'.
                    runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                    # OPTIONS_FILE

                # Run the GATK RealignerTargetCreator and GATK IndelRealigner analyses, unless configured to skip them.

                if self.skip_indel_realignment:
                    file_path_process_read_group.realigned_bam = file_path_process_read_group.duplicates_marked_bam
                    file_path_process_read_group.realigned_bai = file_path_process_read_group.duplicates_marked_bai
                    file_path_process_read_group.realigned_md5 = file_path_process_read_group.duplicates_marked_md5
                else:
                    # Run the GATK RealignerTargetCreator analysis as the first-pass walker
                    # for the GATK IndelRealigner analysis.

                    runnable_step = RunnableStepGATK(
                        name='process_lane_gatk_realigner_target_creator',
                        java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_gatk)
                    runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_gatk_option(key='analysis_type', value='RealignerTargetCreator')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                    if self.known_sites_realignment:  # not None and not empty
                        for file_path in self.known_sites_realignment:
                            runnable_step.add_gatk_option(key='known', value=file_path, override=True)
                    runnable_step.add_gatk_option(
                        key='input_file',
                        value=file_path_process_read_group.duplicates_marked_bam)
                    runnable_step.add_gatk_option(
                        key='out',
                        value=file_path_process_read_group.realigner_targets)

                    # Run the GATK IndelRealigner analysis as a second-pass walker after the
                    # GATK RealignerTargetCreator analysis.

                    runnable_step = RunnableStepGATK(
                        name='process_lane_gatk_indel_realigner',
                        obsolete_file_path_list=[
                            file_path_process_read_group.duplicates_marked_bam,
                            file_path_process_read_group.duplicates_marked_bai,
                            file_path_process_read_group.duplicates_marked_md5,
                        ],
                        java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_gatk)
                    runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_gatk_option(key='analysis_type', value='IndelRealigner')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                    runnable_step.add_gatk_switch(key='keep_program_records')
                    runnable_step.add_gatk_switch(key='generate_md5')
                    runnable_step.add_gatk_option(key='bam_compression', value='9')
                    if self.known_sites_realignment:  # not None and not empty
                        for file_path in self.known_sites_realignment:
                            runnable_step.add_gatk_option(key='knownAlleles', value=file_path, override=True)
                    runnable_step.add_gatk_option(
                        key='input_file',
                        value=file_path_process_read_group.duplicates_marked_bam)
                    runnable_step.add_gatk_option(
                        key='targetIntervals',
                        value=file_path_process_read_group.realigner_targets)
                    runnable_step.add_gatk_option(
                        key='out',
                        value=file_path_process_read_group.realigned_bam)

                # Run the GATK BaseRecalibrator analysis as a first-pass walker
                # for the GATK PrintReads analysis.

                runnable_step = RunnableStepGATK(
                    name='process_lane_gatk_base_recalibrator_pre',
                    java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_gatk)
                runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='BaseRecalibrator')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                if self.known_sites_recalibration:  # not None and not empty
                    for file_path in self.known_sites_recalibration:
                        runnable_step.add_gatk_option(key='knownSites', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_read_group.realigned_bam)
                runnable_step.add_gatk_option(key='out', value=file_path_process_read_group.recalibration_table_pre)

                # Run the GATK BaseRecalibrator on-the-fly recalibration analysis to generate plots.

                runnable_step = RunnableStepGATK(
                    name='process_lane_gatk_base_recalibrator_post',
                    java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_gatk)
                runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='BaseRecalibrator')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                if self.known_sites_recalibration:  # not None and not empty
                    for file_path in self.known_sites_recalibration:
                        runnable_step.add_gatk_option(key='knownSites', value=file_path, override=True)
                runnable_step.add_gatk_option(key='BQSR', value=file_path_process_read_group.recalibration_table_pre)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_read_group.realigned_bam)
                runnable_step.add_gatk_option(key='out', value=file_path_process_read_group.recalibration_table_post)

                # Run the GATK AnalyzeCovariates analysis to create a recalibration plot.

                runnable_step = RunnableStepGATK(
                    name='process_lane_gatk_analyze_covariates',
                    java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_gatk)
                runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='AnalyzeCovariates')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                runnable_step.add_gatk_option(
                    key='afterReportFile',
                    value=file_path_process_read_group.recalibration_table_post)
                runnable_step.add_gatk_option(
                    key='beforeReportFile',
                    value=file_path_process_read_group.recalibration_table_pre)
                runnable_step.add_gatk_option(
                    key='plotsReportFile',
                    value=file_path_process_read_group.recalibration_plot)
                # runnable_step.add_gatk_option(key='logging_level', value='DEBUG')

                # Run the GATK PrintReads analysis as second-pass walker after the GATK BaseRecalibrator analysis.

                runnable_step = RunnableStepGATK(
                    name='process_lane_gatk_print_reads',
                    obsolete_file_path_list=[
                        file_path_process_read_group.realigned_bam,
                        file_path_process_read_group.realigned_bai,
                        file_path_process_read_group.realigned_md5,
                    ],
                    java_temporary_path=runnable_process_lane.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_gatk)
                runnable_process_lane.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='PrintReads')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                runnable_step.add_gatk_switch(key='keep_program_records')
                runnable_step.add_gatk_switch(key='generate_md5')
                runnable_step.add_gatk_option(key='bam_compression', value='9')
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_read_group.realigned_bam)
                runnable_step.add_gatk_option(key='BQSR', value=file_path_process_read_group.recalibration_table_pre)
                runnable_step.add_gatk_option(key='out', value=file_path_process_read_group.recalibrated_bam)

            ###############################
            # Step 3: Process per sample. #
            ###############################
            #
            # Picard MergeSamFiles                      (process_sample_picard_merge_sam_files)
            # Picard MarkDuplicates                     (process_sample_picard_mark_duplicates)
            # GATK RealignerTargetCreator               (process_sample_gatk_realigner_target_creator)
            # GATK IndelRealigner                       (process_sample_gatk_indel_realigner)
            # Symbolic link                             (process_sample_link_bam_bai)
            # Picard CollectAlignmentSummaryMetrics     (process_sample_picard_collect_alignment_summary_metrics)
            # GATK HaplotypeCaller                      (process_sample_gatk_haplotype_caller)

            calling_intervals = VariantCallingGATKCallingIntervals.from_sample(sample=sample)

            target_intervals = VariantCallingGATKTargetIntervals.from_sample(sample=sample)

            file_path_process_sample = self.get_file_path_process_sample(sample_name=sample.name)

            # Create a bsf.procedure.Runnable and bsf.process.Executable for processing each Sample.

            runnable_process_sample = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_process_sample(sample_name=sample.name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))
            executable_process_sample = self.set_stage_runnable(
                stage=stage_process_sample,
                runnable=runnable_process_sample)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            for runnable_process_lane in runnable_process_read_group_list:
                executable_process_sample.dependencies.append(runnable_process_lane.name)

            if use_cache:
                reference_process_sample = runnable_process_sample.get_cache_file_path(
                    file_path=self.bwa_genome_db,
                    absolute=True)
            else:
                reference_process_sample = self.bwa_genome_db

            if len(runnable_process_read_group_list) == 1:
                # If there is only one read group, sample-level read processing can be skipped.
                # Rename files on the basis of the first and only list component.
                runnable_process_lane = runnable_process_read_group_list[0]
                file_path_process_read_group = FilePathProcessReadGroup(prefix=runnable_process_lane.name)

                # Rename the BAM file.
                runnable_step = RunnableStepMove(
                    name='process_sample_move_recalibrated_bam',
                    source_path=file_path_process_read_group.recalibrated_bam,
                    target_path=file_path_process_sample.realigned_bam)
                runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                # Rename the BAI file.
                runnable_step = RunnableStepMove(
                    name='process_sample_move_recalibrated_bai',
                    source_path=file_path_process_read_group.recalibrated_bai,
                    target_path=file_path_process_sample.realigned_bai)
                runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                # Rename the MD5 file.
                runnable_step = RunnableStepMove(
                    name='process_sample_move_recalibrated_md5',
                    source_path=file_path_process_read_group.recalibrated_md5,
                    target_path=file_path_process_sample.realigned_md5)
                runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                # Link the Picard Duplicate Metrics if it has been created.
                if not self.skip_mark_duplicates:
                    runnable_step = RunnableStepLink(
                        name='process_sample_link_duplicate_metrics',
                        source_path=file_path_process_read_group.duplicate_metrics,
                        target_path=file_path_process_sample.duplicate_metrics)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)
            else:
                # Run the Picard MergeSamFiles analysis.

                runnable_step = RunnableStepPicard(
                    name='process_sample_picard_merge_sam_files',
                    java_temporary_path=runnable_process_sample.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='MergeSamFiles')
                runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_picard_option(key='COMMENT', value='Merged from the following files:')
                for runnable_process_lane in runnable_process_read_group_list:
                    file_path_process_read_group = FilePathProcessReadGroup(prefix=runnable_process_lane.name)
                    runnable_step.add_picard_option(
                        key='COMMENT',
                        value=file_path_process_read_group.recalibrated_bam,
                        override=True)
                    runnable_step.add_picard_option(
                        key='INPUT',
                        value=file_path_process_read_group.recalibrated_bam,
                        override=True)
                    runnable_step.obsolete_file_path_list.append(file_path_process_read_group.recalibrated_bam)
                    runnable_step.obsolete_file_path_list.append(file_path_process_read_group.recalibrated_bai)
                    runnable_step.obsolete_file_path_list.append(file_path_process_read_group.recalibrated_md5)
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_process_sample.merged_bam)
                runnable_step.add_picard_option(key='USE_THREADING', value='true')
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_process_sample.temporary_directory_path(absolute=False))
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='4000000')
                # CREATE_INDEX defaults to 'false'.
                runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                # CREATE_MD5_FILE defaults to 'false'.
                runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                # OPTIONS_FILE

                # Run the Picard MarkDuplicates analysis, unless configured to skip it.

                if self.skip_mark_duplicates:
                    # Rename the files after Picard MergeSamFiles to files after Picard MarkDuplicates.
                    runnable_step = RunnableStepMove(
                        name='process_sample_move_merged_bam',
                        source_path=file_path_process_sample.merged_bam,
                        target_path=file_path_process_sample.duplicates_marked_bam)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step = RunnableStepMove(
                        name='process_sample_move_merged_bai',
                        source_path=file_path_process_sample.merged_bai,
                        target_path=file_path_process_sample.duplicates_marked_bai)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step = RunnableStepMove(
                        name='process_sample_move_merged_md5',
                        source_path=file_path_process_sample.merged_md5,
                        target_path=file_path_process_sample.duplicates_marked_md5)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)
                else:
                    # Run the Picard MarkDuplicates analysis.
                    # Optical duplicates should already have been flagged in the lane-specific processing step.

                    runnable_step = RunnableStepPicard(
                        name='process_sample_picard_mark_duplicates',
                        obsolete_file_path_list=[
                            file_path_process_sample.merged_bam,
                            file_path_process_sample.merged_bai,
                            file_path_process_sample.merged_md5,
                        ],
                        java_temporary_path=runnable_process_sample.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_picard,
                        picard_command='MarkDuplicates')
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_picard_option(
                        key='INPUT',
                        value=file_path_process_sample.merged_bam)
                    runnable_step.add_picard_option(
                        key='OUTPUT',
                        value=file_path_process_sample.duplicates_marked_bam)
                    runnable_step.add_picard_option(
                        key='METRICS_FILE',
                        value=file_path_process_sample.duplicate_metrics)
                    # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
                    # as otherwise optical duplicates could not be detected. This is a consequence of using
                    # Illumina2bam rather than Picard ExtractIlluminaBarcodes, IlluminaBasecallsToFastq and
                    # IlluminaBasecallsToSam.
                    # See BioStar post: http://www.biostars.org/p/12538/
                    # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    runnable_step.add_picard_option(
                        key='READ_NAME_REGEX',
                        value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*')
                    runnable_step.add_picard_option(
                        key='OPTICAL_DUPLICATE_PIXEL_DISTANCE',
                        value='5000')
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_process_sample.temporary_directory_path(absolute=False))
                    # VERBOSITY defaults to 'INFO'.
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    # QUIET defaults to 'false'.
                    # VALIDATION_STRINGENCY defaults to 'STRICT'.
                    # COMPRESSION_LEVEL defaults to '5'.
                    runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
                    # MAX_RECORDS_IN_RAM defaults to '500000'.
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='4000000')
                    # CREATE_INDEX defaults to 'false'.
                    runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                    # CREATE_MD5_FILE defaults to 'false'.
                    runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
                    # OPTIONS_FILE

                # Run the GATK RealignerTargetCreator and GATK IndelRealigner analyses, unless configured to skip them.

                if self.skip_indel_realignment:
                    # Rename files after Picard MarkDuplicates to files after GATK IndelRealigner.
                    runnable_step = RunnableStepMove(
                        name='process_sample_move_duplicates_marked_bam',
                        source_path=file_path_process_sample.duplicates_marked_bam,
                        target_path=file_path_process_sample.realigned_bam)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step = RunnableStepMove(
                        name='process_sample_move_duplicates_marked_bai',
                        source_path=file_path_process_sample.duplicates_marked_bai,
                        target_path=file_path_process_sample.realigned_bai)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step = RunnableStepMove(
                        name='process_sample_move_duplicates_marked_md5',
                        source_path=file_path_process_sample.duplicates_marked_md5,
                        target_path=file_path_process_sample.realigned_md5)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)
                else:
                    # Run the GATK RealignerTargetCreator analysis as the first-pass walker
                    # for the GATK IndelRealigner analysis.

                    runnable_step = RunnableStepGATK(
                        name='process_sample_gatk_realigner_target_creator',
                        java_temporary_path=runnable_process_sample.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_gatk)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_gatk_option(key='analysis_type', value='RealignerTargetCreator')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_sample)
                    if self.known_sites_realignment:  # not None and not empty
                        for file_path in self.known_sites_realignment:
                            runnable_step.add_gatk_option(key='known', value=file_path, override=True)
                    runnable_step.add_gatk_option(
                        key='input_file',
                        value=file_path_process_sample.duplicates_marked_bam)
                    runnable_step.add_gatk_option(key='out', value=file_path_process_sample.realigner_targets)

                    # Run the GATK IndelRealigner analysis as a second-pass walker
                    # after the GATK RealignerTargetCreator analysis.

                    runnable_step = RunnableStepGATK(
                        name='process_sample_gatk_indel_realigner',
                        obsolete_file_path_list=[
                            file_path_process_sample.duplicates_marked_bam,
                            file_path_process_sample.duplicates_marked_bai,
                            file_path_process_sample.duplicates_marked_md5,
                        ],
                        java_temporary_path=runnable_process_sample.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_gatk)
                    runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_gatk_option(key='analysis_type', value='IndelRealigner')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_sample)
                    runnable_step.add_gatk_switch(key='keep_program_records')
                    runnable_step.add_gatk_switch(key='generate_md5')
                    runnable_step.add_gatk_option(key='bam_compression', value='9')
                    if self.known_sites_realignment:  # not None and not empty
                        for file_path in self.known_sites_realignment:
                            runnable_step.add_gatk_option(key='knownAlleles', value=file_path, override=True)
                    runnable_step.add_gatk_option(
                        key='input_file',
                        value=file_path_process_sample.duplicates_marked_bam)
                    runnable_step.add_gatk_option(
                        key='targetIntervals',
                        value=file_path_process_sample.realigner_targets)
                    runnable_step.add_gatk_option(
                        key='out',
                        value=file_path_process_sample.realigned_bam)
                    # For debugging only.
                    # runnable_step.add_gatk_option(key='logging_level', value='DEBUG')

            # Set a symbolic link to a samtools-style (*.bam.bai) index file.
            # The UCSC Genome Browser suddenly requires it for its track hubs.

            runnable_step = RunnableStepLink(
                name='process_sample_link_bam_bai',
                source_path=file_path_process_sample.realigned_bai,
                target_path=file_path_process_sample.realigned_bam_bai)
            runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

            # Run the Picard CollectAlignmentSummaryMetrics analysis.

            runnable_step = RunnableStepPicard(
                name='process_sample_picard_collect_alignment_summary_metrics',
                java_temporary_path=runnable_process_sample.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx6G',
                java_jar_path=self.java_archive_picard,
                picard_command='CollectAlignmentSummaryMetrics')
            runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_picard_option(key='INPUT', value=file_path_process_sample.realigned_bam)
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_process_sample.alignment_summary_metrics)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
            runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_process_sample.temporary_directory_path(absolute=False))
            # VERBOSITY defaults to 'INFO'.
            runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
            # QUIET defaults to 'false'.
            # VALIDATION_STRINGENCY defaults to 'STRICT'.
            # COMPRESSION_LEVEL defaults to '5'.
            # MAX_RECORDS_IN_RAM defaults to '500000'.
            # CREATE_INDEX defaults to 'false'.
            # CREATE_MD5_FILE defaults to 'false'.
            # OPTIONS_FILE

            # Run the GATK HaplotypeCaller analysis per sample.

            runnable_step = RunnableStepGATK(
                name='process_sample_gatk_haplotype_caller',
                java_temporary_path=runnable_process_sample.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx8G',
                java_jar_path=self.java_archive_gatk)
            runnable_process_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='HaplotypeCaller')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_sample)
            if calling_intervals.calling_path:
                # If calling intervals are available, the Haplotype Caller analysis is preferentially run only on them.
                runnable_step.add_gatk_option(key='intervals', value=calling_intervals.calling_path)
                if self.interval_padding is not None:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            elif target_intervals.targets_path:
                # If target intervals are available, the Haplotype Caller analysis is preferentially run only on them.
                runnable_step.add_gatk_option(key='intervals', value=target_intervals.targets_path)
                if self.interval_padding is not None:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            elif self.scatter_intervals_path:
                # If a scatter interval file path is available assign it without padding.
                runnable_step.add_gatk_option(key='intervals', value=self.scatter_intervals_path)
            elif self.exclude_intervals_list:  # not None and not empty
                # If target intervals are not available, decoy sequences etc. can be excluded
                # from the genome-wide Haplotype Caller analysis.
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            runnable_step.add_gatk_option(key='genotyping_mode', value='DISCOVERY')
            runnable_step.add_gatk_option(key='emitRefConfidence', value='GVCF')
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
            runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
            runnable_step.add_gatk_option(key='out', value=file_path_process_sample.raw_variants_gvcf_vcf)

            # Finally, record the process_sample bsf.procedure.Runnable for the merge_cohort stage under the Sample's
            # 'Cohort Name' annotation, or if it does not exist, under the cohort name defined in the
            # bsf.analysis.Analysis in the configuration file.

            if 'Cohort Name' in sample.annotation_dict:
                cohort_key = sample.annotation_dict['Cohort Name'][0]
            else:
                cohort_key = self.cohort_name

            if cohort_key not in runnable_process_sample_dict:
                runnable_process_sample_dict[cohort_key] = list()
            _runnable_process_sample_list = runnable_process_sample_dict[cohort_key]
            _runnable_process_sample_list.append((
                runnable_process_sample,
                file_path_process_sample.raw_variants_gvcf_vcf))

            ################################
            # Step 4: Diagnose the sample. #
            ################################
            #
            # GATK CallableLoci                      (diagnose_sample_gatk_callable_loci)
            # UCSC bedSort                           (diagnose_sample_bed_sort)
            # UCSC bedToBigBed                       (diagnose_sample_bed_sort)
            # bsfR bsf_variant_calling_coverage.R    (diagnose_sample_coverage)
            # bsfR bsf_variant_calling_insert_size.R (diagnose_sample_insert_size)
            # For target intervals only:
            #   GATK DiagnoseTarget                  (diagnose_sample_gatk_diagnose_target)
            #   GATK QualifyMissingIntervals         (diagnose_sample_gatk_qualify_missing_intervals)
            #   Picard CollectHsMetrics              (diagnose_sample_picard_collect_hybrid_selection_metrics)

            file_path_diagnose_sample = self.get_file_path_diagnose_sample(sample_name=sample.name)

            # Create a bsf.procedure.Runnable and bsf.process.Executable for diagnosing each sample.

            runnable_diagnose_sample = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_diagnose_sample(sample_name=sample.name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))
            executable_diagnose_sample = self.set_stage_runnable(
                stage=stage_diagnose_sample,
                runnable=runnable_diagnose_sample)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            executable_diagnose_sample.dependencies.append(runnable_process_sample.name)
            # Set dependencies for succeeding bsf.procedure.Runnable or bsf.process.Executable objects.
            runnable_diagnose_sample_list.append(runnable_diagnose_sample)

            if use_cache:
                reference_diagnose_sample = runnable_diagnose_sample.get_cache_file_path(
                    file_path=self.bwa_genome_db,
                    absolute=True)
            else:
                reference_diagnose_sample = self.bwa_genome_db

            # Run the GATK CallableLoci analysis per sample.

            runnable_step = RunnableStepGATK(
                name='diagnose_sample_gatk_callable_loci',
                java_temporary_path=runnable_diagnose_sample.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx6G',
                java_jar_path=self.java_archive_gatk)
            runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

            # Read RunnableStep options from configuration sections:
            # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_callable_loci]
            # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_callable_loci.]
            # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_callable_loci..]
            self.set_runnable_step_configuration(runnable_step=runnable_step)
            # CommandLineGATK
            # Required Parameters
            runnable_step.add_gatk_option(key='analysis_type', value='CallableLoci')
            # Optional Inputs
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_diagnose_sample)
            if calling_intervals.calling_path:
                # If calling intervals are available, the Callable Loci analysis is preferentially run only on them.
                runnable_step.add_gatk_option(key='intervals', value=calling_intervals.calling_path)
                if self.interval_padding is not None:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            elif target_intervals.targets_path:
                # If target intervals are available, the Callable Loci analysis is run only on them.
                runnable_step.add_gatk_option(key='intervals', value=target_intervals.targets_path)
                if self.interval_padding is not None:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            elif self.exclude_intervals_list:  # not None and not empty
                # If target intervals are not available, decoy sequences etc. can be excluded
                # from the genome-wide Callable Loci analysis.
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
            # Required Outputs
            runnable_step.add_gatk_option(key='summary', value=file_path_diagnose_sample.callable_txt)
            # Optional Outputs
            runnable_step.add_gatk_option(key='out', value=file_path_diagnose_sample.callable_bed)

            # Run the UCSC bedSort tool.

            runnable_step = RunnableStep(
                name='diagnose_sample_bed_sort',
                program='bedSort')
            runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step.arguments.append(file_path_diagnose_sample.callable_bed)
            runnable_step.arguments.append(file_path_diagnose_sample.sorted_bed)

            # Run the UCSC bedToBigBed tool.

            runnable_step = RunnableStep(
                name='diagnose_sample_bed_to_big_bed',
                program='bedToBigBed',
                obsolete_file_path_list=[
                    file_path_diagnose_sample.sorted_bed,
                ])
            runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_pair_short(key='type', value='bed4')
            runnable_step.arguments.append(file_path_diagnose_sample.sorted_bed)
            runnable_step.arguments.append(reference_diagnose_sample + '.fai')
            runnable_step.arguments.append(file_path_diagnose_sample.callable_bb)

            # Run the bsfR bsf_variant_calling_coverage.R script.

            runnable_step = RunnableStep(
                name='diagnose_sample_coverage',
                program='bsf_variant_calling_coverage.R')
            runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

            # Read RunnableStep options from configuration sections:
            # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_coverage]
            self.set_runnable_step_configuration(runnable_step=runnable_step)
            runnable_step.add_option_long(key='exons', value=self.genome_annotation_gtf)
            runnable_step.add_option_long(key='callable-loci', value=file_path_diagnose_sample.callable_bed)
            if target_intervals.targets_path is None:
                # If a target interval path has not been defined, run without it.
                pass
            elif target_intervals.targets_path.endswith('.bed'):
                runnable_step.add_option_long(key='targets', value=target_intervals.targets_path)
            elif target_intervals.targets_path.endswith('.interval_list'):
                runnable_step.add_option_long(key='targets', value=target_intervals.targets_path[:-13] + 'bed')
            elif target_intervals.targets_path:
                runnable_step.add_option_long(key='targets', value=target_intervals.targets_path)

            # Run the bsfR bsf_variant_calling_insert_size.R script.

            runnable_step = RunnableStep(
                name='diagnose_sample_insert_size',
                program='bsf_variant_calling_insert_size.R')
            runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_option_long(key='file-path', value=file_path_process_sample.realigned_bam)

            if target_intervals.targets_path:
                # Run the GATK DiagnoseTarget analysis per sample, only if targets have been defined.

                runnable_step = RunnableStepGATK(
                    name='diagnose_sample_gatk_diagnose_target',
                    java_temporary_path=runnable_diagnose_sample.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx8G',
                    java_jar_path=self.java_archive_gatk)
                runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

                # Read RunnableStep options from configuration sections:
                # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_diagnose_target]
                # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_diagnose_target.]
                # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_diagnose_target..]
                self.set_runnable_step_configuration(runnable_step=runnable_step)
                # CommandLineGATK
                # Required Parameters
                runnable_step.add_gatk_option(key='analysis_type', value='DiagnoseTargets')
                # Optional Inputs
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_diagnose_sample)
                # The Diagnose Targets analysis is run on the target intervals, only.
                # This will not include excluded intervals, such as decoy sequences.
                runnable_step.add_gatk_option(key='intervals', value=target_intervals.targets_path)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
                # DiagnoseTargets
                # Optional Outputs
                runnable_step.add_gatk_option(
                    key='missing_intervals',
                    value=file_path_diagnose_sample.missing_intervals)
                runnable_step.add_gatk_option(key='out', value=file_path_diagnose_sample.diagnose_targets_vcf)
                # Optional Parameters

                # Run the GATK QualifyMissingIntervals analysis per sample, only if targets have been defined.

                runnable_step = RunnableStepGATK(
                    name='diagnose_sample_gatk_qualify_missing_intervals',
                    java_temporary_path=runnable_diagnose_sample.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx8G',
                    java_jar_path=self.java_archive_gatk)
                runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

                # Read RunnableStep options from configuration sections:
                # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_qualify_missing_intervals]
                # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_qualify_missing_intervals.]
                # [bsf.analyses.variant_calling.VariantCallingGATK.diagnose_sample_gatk_qualify_missing_intervals..]
                self.set_runnable_step_configuration(runnable_step=runnable_step)
                # CommandLineGATK
                # Required Parameters
                runnable_step.add_gatk_option(key='analysis_type', value='QualifyMissingIntervals')
                # Optional Inputs
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_diagnose_sample)
                # The Qualify Missing Intervals analysis is run only on the missing intervals
                # of the Diagnose Targets analysis above, regardless.
                # This will not include excluded intervals, such as decoy sequences.
                runnable_step.add_gatk_option(key='intervals', value=file_path_diagnose_sample.missing_intervals)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
                # QualifyMissingIntervals
                # Required Parameters
                runnable_step.add_gatk_option(key='targetsfile', value=target_intervals.targets_path)
                # Optional Outputs
                runnable_step.add_gatk_option(key='out', value=file_path_diagnose_sample.missing_report)
                # Optional Parameters
                if target_intervals.probes_path:
                    runnable_step.add_gatk_option(key='baitsfile', value=target_intervals.probes_path)

                # Run the Picard CollectHsMetrics analysis per sample, only if targets have been defined.

                runnable_step = RunnableStepPicard(
                    name='diagnose_sample_picard_collect_hybrid_selection_metrics',
                    java_temporary_path=runnable_diagnose_sample.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx12G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='CollectHsMetrics')
                runnable_diagnose_sample.add_runnable_step(runnable_step=runnable_step)

                if target_intervals.probes_path:
                    runnable_step.add_picard_option(key='BAIT_INTERVALS', value=target_intervals.probes_path)
                else:
                    runnable_step.add_picard_option(key='BAIT_INTERVALS', value=target_intervals.targets_path)
                if target_intervals.name:
                    runnable_step.add_picard_option(key='BAIT_SET_NAME', value=target_intervals.name)
                runnable_step.add_picard_option(key='TARGET_INTERVALS', value=target_intervals.targets_path)
                runnable_step.add_picard_option(key='INPUT', value=file_path_process_sample.realigned_bam)
                runnable_step.add_picard_option(
                    key='OUTPUT',
                    value=file_path_diagnose_sample.hybrid_selection_metrics)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
                runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
                # PER_TARGET_COVERAGE
                # TMP_DIR
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_diagnose_sample.temporary_directory_path(absolute=False))
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                # OPTIONS_FILE

        #########################################
        # Step 5: Hierarchically merge cohorts. #
        #########################################
        #
        # GATK CombineGVCFs for each cohort and each sample (merge_cohort_gatk_combine_gvcfs)
        # GATK CombineGVCFs for each cohort                 (merge_cohort_gatk_combine_gvcfs)
        #
        # The cohorts to merge in need to be configurable and it would be essential,
        # How should sample names be checked before the merge to avoid clashes?
        # Should sample annotation sheets be read or can the combined GVCF file be read in
        # to extract the actual sample names?

        # Create a Python dict of Python str (cohort name) key and Python list of bsf.procedure.Runnable value data.
        # Initialise a single key with the final cohort name and an empty list for merging the final cohort.

        runnable_merge_cohort_dict = {self.cohort_name: []}
        """ @type runnable_merge_cohort_dict: 
            dict[str, list[(Runnable | str, str)]] """

        runnable_merge_cohort_list = runnable_merge_cohort_dict[self.cohort_name]

        # Run the GATK CombineGVCFs analysis for each cohort and Sample defined in this project to build up
        # cohort-specific GVCF files.

        for cohort_key in runnable_process_sample_dict:
            runnable_merge_cohort_list.append(
                run_merge_cohort_scatter_gather(
                    analysis_stage=stage_merge_cohort,
                    cohort_runnable_dict=runnable_process_sample_dict,
                    cohort_name=cohort_key))

        # Run the GATK CombineGVCF analysis once more to merge all cohort-specific GVCF files defined in this project.

        if len(runnable_merge_cohort_list) == 1:
            # If the cohort-specific bsf.procedure.Runnable list has only one component,
            # the merge has already been completed.
            runnable_merge_cohort, file_path_merge_cohort_gvcf = runnable_merge_cohort_list[-1]
        elif len(runnable_merge_cohort_list) > 1:
            runnable_merge_cohort, file_path_merge_cohort_gvcf = run_merge_cohort_scatter_gather(
                analysis_stage=stage_merge_cohort,
                cohort_runnable_dict=runnable_merge_cohort_dict,
                cohort_name=self.cohort_name)
        else:
            raise Exception('Unexpected number of bsf.procedure.Runnable objects on the merge_cohort list.')

        # Run an additional GATK CombineGVCFs analysis to merge into a super-cohort, if defined.

        if self.accessory_cohort_gvcfs:  # not None and not empty
            # If accessory cohorts are defined, initialise a new Python dict of Python str cohort name key and
            # Python list of bsf.procedure.Runnable value data. Initialise the list with the last
            # bsf.procedure.Runnable object and extend with the the list of accessory cohort file names.
            # The run_merge_cohort_scatter_gather() method can cope with bsf.procedure.Runnable or
            # str objects.
            cohort_key = '_'.join((self.cohort_name, 'accessory'))
            runnable_merge_cohort_dict = {cohort_key: [(runnable_merge_cohort, file_path_merge_cohort_gvcf)]}
            runnable_merge_cohort_list = runnable_merge_cohort_dict[cohort_key]
            for accessory_cohort_gvcf in self.accessory_cohort_gvcfs:
                runnable_merge_cohort_list.append((accessory_cohort_gvcf, accessory_cohort_gvcf))

            runnable_merge_cohort, file_path_merge_cohort_gvcf = run_merge_cohort_scatter_gather(
                analysis_stage=stage_merge_cohort,
                cohort_runnable_dict=runnable_merge_cohort_dict,
                cohort_name=cohort_key)

            # prefix_merge_cohort = self.get_prefix_merge_cohort(cohort_name=cohort_key)
        else:
            # prefix_merge_cohort = self.get_prefix_merge_cohort(cohort_name=self.cohort_name)
            pass

        # Specify the final FilePathMergeCohort object from the final bsf.procedure.Runnable object.
        # file_path_merge_cohort = FilePathMergeCohort(prefix=prefix_merge_cohort)

        ###############################
        # Step 6: Process per cohort. #
        ###############################
        #
        # GATK CombineGVCFs                     (process_cohort_gatk_combine_gvcfs_accessory)
        # GATK GenotypeGVCFs                    (process_cohort_gatk_genotype_gvcfs)
        # GATK VariantRecalibrator for SNPs     (process_cohort_gatk_variant_recalibrator_snp)
        # GATK ApplyRecalibration for SNPs      (process_cohort_gatk_apply_recalibration_snp)
        # GATK VariantRecalibrator for INDELs   (process_cohort_gatk_variant_recalibrator_indel)
        # GATK ApplyRecalibration for INDELs    (process_cohort_gatk_apply_recalibration_indel)
        # GATK SelectVariants                   (process_cohort_gatk_select_variants_cohort)
        # snpEff                                (process_cohort_snpeff)
        # GATK VariantAnnotator                 (process_cohort_gatk_variant_annotator)

        file_path_genotype_cohort = FilePathGenotypeCohort(
            prefix=self.get_prefix_process_cohort(
                cohort_name=self.cohort_name))

        runnable_process_cohort_gather = run_genotype_cohort_scatter_gather(
            file_path_cohort_gvcf=file_path_merge_cohort_gvcf)

        file_path_process_cohort = self.get_file_path_process_cohort(cohort_name=self.cohort_name)

        # Create a bsf.procedure.Runnable and bsf.process.Executable for processing the cohort.

        runnable_process_cohort = self.add_runnable_consecutive(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_process_cohort(cohort_name=self.cohort_name),
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                cache_path_dict=self._cache_path_dict,
                debug=self.debug))
        executable_process_cohort = self.set_stage_runnable(
            stage=stage_process_cohort,
            runnable=runnable_process_cohort)
        # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
        executable_process_cohort.dependencies.append(runnable_process_cohort_gather.name)

        if use_cache:
            reference_process_cohort = runnable_process_cohort.get_cache_file_path(
                file_path=self.bwa_genome_db,
                absolute=True)
        else:
            reference_process_cohort = self.bwa_genome_db

        # Run the VQSR procedure on SNPs.

        if self.vqsr_skip_snp:
            file_path_process_cohort.recalibrated_snp_raw_indel_vcf = file_path_genotype_cohort.genotyped_raw_vcf
            file_path_process_cohort.recalibrated_snp_raw_indel_tbi = file_path_genotype_cohort.genotyped_raw_tbi
        else:

            # Run the GATK VariantRecalibrator analysis on SNPs.

            runnable_step = RunnableStepGATK(
                name='process_cohort_gatk_variant_recalibrator_snp',
                java_temporary_path=runnable_process_cohort.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx8G',
                java_jar_path=self.java_archive_gatk)
            runnable_process_cohort.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            runnable_step.add_gatk_option(key='mode', value='SNP')
            if self.vqsr_resources_snp_dict:  # not None and not empty
                for resource_name, resource_dict in self.vqsr_resources_snp_dict.items():
                    resource_option = 'resource:{},known={},training={},truth={},prior={}'.format(
                        resource_name,
                        resource_dict['known'],
                        resource_dict['training'],
                        resource_dict['truth'],
                        resource_dict['prior'])
                    runnable_step.add_gatk_option(
                        key=resource_option,
                        value=resource_dict['file_path'])
            if self.vqsr_annotations_snp_list:  # not None and not empty
                for annotation in self.vqsr_annotations_snp_list:
                    runnable_step.add_gatk_option(key='use_annotation', value=annotation, override=True)
            if self.vqsr_max_gaussians_pos_snp is not None:
                runnable_step.add_gatk_option(key='maxGaussians', value=str(self.vqsr_max_gaussians_pos_snp))
            runnable_step.add_gatk_option(key='input', value=file_path_genotype_cohort.genotyped_raw_vcf)
            runnable_step.add_gatk_option(key='recal_file', value=file_path_process_cohort.recalibration_snp)
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_process_cohort.tranches_snp)
            runnable_step.add_gatk_option(key='rscript_file', value=file_path_process_cohort.plots_snp)
            if self.vqsr_bad_lod_cutoff_snp is not None:
                runnable_step.add_gatk_option(key='badLodCutoff', value=str(self.vqsr_bad_lod_cutoff_snp))

            # Run the GATK ApplyRecalibration analysis on SNPs.

            runnable_step = RunnableStepGATK(
                name='process_cohort_gatk_apply_recalibration_snp',
                java_temporary_path=runnable_process_cohort.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx4G',
                java_jar_path=self.java_archive_gatk)
            runnable_process_cohort.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='ApplyRecalibration')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            runnable_step.add_gatk_option(key='mode', value='SNP')
            runnable_step.add_gatk_option(key='input', value=file_path_genotype_cohort.genotyped_raw_vcf)
            runnable_step.add_gatk_option(key='recal_file', value=file_path_process_cohort.recalibration_snp)
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_process_cohort.tranches_snp)
            runnable_step.add_gatk_option(key='out', value=file_path_process_cohort.recalibrated_snp_raw_indel_vcf)
            # The lodCutoff (VQSLOD score) filter is not applied for the moment.
            if self.truth_sensitivity_filter_level_snp:
                runnable_step.add_gatk_option(key='ts_filter_level', value=self.truth_sensitivity_filter_level_snp)

        # Run the VQSR procedure on INDELs.

        if self.vqsr_skip_indel:
            file_path_process_cohort.recalibrated_snp_recalibrated_indel_vcf = \
                file_path_process_cohort.recalibrated_snp_raw_indel_vcf
            file_path_process_cohort.recalibrated_snp_recalibrated_indel_tbi = \
                file_path_process_cohort.recalibrated_snp_raw_indel_tbi
        else:

            # Run the GATK VariantRecalibrator analysis on INDELs.

            runnable_step = RunnableStepGATK(
                name='process_cohort_gatk_variant_recalibrator_indel',
                java_temporary_path=runnable_process_cohort.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx8G',
                java_jar_path=self.java_archive_gatk)
            runnable_process_cohort.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            runnable_step.add_gatk_option(key='mode', value='INDEL')
            if self.vqsr_resources_indel_dict:  # not None and not empty
                for resource_name, resource_dict in self.vqsr_resources_indel_dict.items():
                    resource_option = 'resource:{},known={},training={},truth={},prior={}'.format(
                        resource_name,
                        resource_dict['known'],
                        resource_dict['training'],
                        resource_dict['truth'],
                        resource_dict['prior'])
                    runnable_step.add_gatk_option(
                        key=resource_option,
                        value=resource_dict['file_path'])
            if self.vqsr_annotations_indel_list:  # not None and not empty
                for annotation in self.vqsr_annotations_indel_list:
                    runnable_step.add_gatk_option(key='use_annotation', value=annotation, override=True)
            if self.vqsr_max_gaussians_pos_indel is not None:
                runnable_step.add_gatk_option(key='maxGaussians', value=str(self.vqsr_max_gaussians_pos_indel))
            runnable_step.add_gatk_option(key='input', value=file_path_process_cohort.recalibrated_snp_raw_indel_vcf)
            runnable_step.add_gatk_option(key='recal_file', value=file_path_process_cohort.recalibration_indel)
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_process_cohort.tranches_indel)
            runnable_step.add_gatk_option(key='rscript_file', value=file_path_process_cohort.plots_indel)
            if self.vqsr_bad_lod_cutoff_indel is not None:
                runnable_step.add_gatk_option(key='badLodCutoff', value=str(self.vqsr_bad_lod_cutoff_indel))

            # Run the GATK ApplyRecalibration analysis on INDELs.

            runnable_step = RunnableStepGATK(
                name='process_cohort_gatk_apply_recalibration_indel',
                java_temporary_path=runnable_process_cohort.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx4G',
                java_jar_path=self.java_archive_gatk)
            runnable_process_cohort.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='ApplyRecalibration')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            runnable_step.add_gatk_option(key='mode', value='INDEL')
            runnable_step.add_gatk_option(key='input', value=file_path_process_cohort.recalibrated_snp_raw_indel_vcf)
            runnable_step.add_gatk_option(key='recal_file', value=file_path_process_cohort.recalibration_indel)
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_process_cohort.tranches_indel)
            runnable_step.add_gatk_option(
                key='out',
                value=file_path_process_cohort.recalibrated_snp_recalibrated_indel_vcf)
            # The lodCutoff (VQSLOD score) filter is not applied for the moment.
            if self.truth_sensitivity_filter_level_indel:
                runnable_step.add_gatk_option(key='ts_filter_level', value=self.truth_sensitivity_filter_level_indel)

        # Run GATK SelectVariants, in case accessory GVCF files have been used, re-create a multi-sample VCF file
        # with just the samples in this cohort.

        if self.accessory_cohort_gvcfs:  # not None and not empty
            runnable_step = RunnableStepGATK(
                name='process_cohort_gatk_select_variants_cohort',
                java_temporary_path=runnable_process_cohort.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx4G',
                java_jar_path=self.java_archive_gatk)
            runnable_process_cohort.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            runnable_step.add_gatk_option(
                key='variant',
                value=file_path_process_cohort.recalibrated_snp_recalibrated_indel_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_process_cohort.multi_sample_vcf)
            for sample in self.sample_list:
                # Get all PairedReads objects solely to exclude samples without any.
                paired_reads_dict = sample.get_all_paired_reads(
                    replicate_grouping=self.replicate_grouping,
                    exclude=True)
                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue
                runnable_step.add_gatk_option(key='sample_name', value=sample.name, override=True)
            runnable_step.add_gatk_switch(key='excludeNonVariants')
        else:
            file_path_process_cohort.multi_sample_vcf = \
                file_path_process_cohort.recalibrated_snp_recalibrated_indel_vcf
            file_path_process_cohort.multi_sample_tbi = \
                file_path_process_cohort.recalibrated_snp_recalibrated_indel_tbi

        ################################
        # Step 7: Annotate the cohort. #
        ################################

        # Run snpEff annotation.

        runnable_annotate_cohort_snpeff = run_annotate_snpeff(
            prefix=self.get_prefix_annotate_cohort_snpeff(cohort_name=self.cohort_name),
            vcf_file_path=file_path_process_cohort.multi_sample_vcf)
        executable_annotate_cohort_snpeff = self.set_stage_runnable(
            stage=stage_annotate_cohort_snpeff,
            runnable=runnable_annotate_cohort_snpeff)
        # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
        executable_annotate_cohort_snpeff.dependencies.append(runnable_process_cohort.name)

        # Run Ensembl Variant Effect Predictor (VEP) annotation.

        runnable_annotate_cohort_vep = run_annotate_vep(
            prefix=self.get_prefix_annotate_cohort_vep(cohort_name=self.cohort_name),
            vcf_file_path=file_path_process_cohort.multi_sample_vcf)
        executable_annotate_cohort_vep = self.set_stage_runnable(
            stage=stage_annotate_cohort_vep,
            runnable=runnable_annotate_cohort_vep)
        # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
        executable_annotate_cohort_vep.dependencies.append(runnable_process_cohort.name)

        ######################################################
        # Step 8: Re-process and split the cohort by sample. #
        ######################################################
        #
        # GATK SelectVariants   (split_cohort_snpeff_gatk_select_variants_snpeff)
        # GATK VariantsToTable  (split_cohort_snpeff_gatk_variants_to_table_snpeff)
        # GATK SelectVariants   (split_cohort_gatk_select_variants_vep)
        # GATK VariantsToTable  (split_cohort_gatk_variants_to_table_vep)

        file_path_annotate_cohort_snpeff = self.get_file_path_annotate_cohort_snpeff(cohort_name=self.cohort_name)

        file_path_annotate_cohort_vep = self.get_file_path_annotate_cohort_vep(cohort_name=self.cohort_name)

        for sample in self.sample_list:
            # Get all PairedReads objects solely to exclude samples without any.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)
            if not len(paired_reads_dict):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            # Split the snpEff-annotated multi-sample VCF file.

            file_path_split_cohort_snpeff = self.get_file_path_split_cohort_snpeff(sample_name=sample.name)

            runnable_split_cohort_snpeff = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_split_cohort_snpeff(sample_name=sample.name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))
            executable_split_cohort_snpeff = self.set_stage_runnable(
                stage=stage_split_cohort_snpeff,
                runnable=runnable_split_cohort_snpeff)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            executable_split_cohort_snpeff.dependencies.append(runnable_annotate_cohort_snpeff.name)

            if use_cache:
                reference_split_cohort_snpeff = runnable_split_cohort_snpeff.get_cache_file_path(
                    file_path=self.bwa_genome_db,
                    absolute=True)
            else:
                reference_split_cohort_snpeff = self.bwa_genome_db

            # Run the GATK SelectVariants analysis to split multi-sample VCF files into one per sample.

            runnable_step = RunnableStepGATK(
                name='split_cohort_snpeff_gatk_select_variants_snpeff',
                java_temporary_path=runnable_split_cohort_snpeff.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx2G',
                java_jar_path=self.java_archive_gatk)
            runnable_split_cohort_snpeff.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_snpeff)
            runnable_step.add_gatk_option(key='variant', value=file_path_annotate_cohort_snpeff.annotated_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_snpeff.sample_vcf)
            runnable_step.add_gatk_option(key='sample_name', value=sample.name)
            runnable_step.add_gatk_switch(key='excludeNonVariants')

            # Run the GATK VariantsToTable analysis.

            if self.variants_to_table:
                runnable_step = RunnableStepGATK(
                    name='split_cohort_snpeff_gatk_variants_to_table_snpeff',
                    java_temporary_path=runnable_split_cohort_snpeff.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    java_jar_path=self.java_archive_gatk)
                runnable_split_cohort_snpeff.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_snpeff)
                runnable_step.add_gatk_option(key='variant', value=file_path_split_cohort_snpeff.sample_vcf)
                runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_snpeff.sample_tsv)
                runnable_step.add_gatk_switch(key='showFiltered')
                # Set of fixed VCF fields.
                for field_name in variants_to_table_fields['fixed']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK Haplotype Caller-specific INFO fields.
                for field_name in variants_to_table_fields['haplotype_caller']['info']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK Haplotype Caller-specific genotype fields.
                for field_name in variants_to_table_fields['haplotype_caller']['format']:
                    runnable_step.add_gatk_option(key='genotypeFields', value=field_name, override=True)
                # Set of snpEff-specific INFO fields.
                for field_name in variants_to_table_fields['snpeff']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Automatically add all fields defined for the Variant Annotator resources, above.
                if self.annotation_resources_dict:  # not None and not empty
                    for resource_name, (file_path, annotation_list) in self.annotation_resources_dict.items():
                        if file_path and annotation_list:
                            for annotation in annotation_list:
                                runnable_step.add_gatk_option(
                                    key='fields',
                                    value='.'.join((resource_name, annotation)),
                                    override=True)

            # Split the Ensembl Variant Effect Predictor-annotated multi-sample VCF file.

            file_path_split_cohort_vep = self.get_file_path_split_cohort_vep(sample_name=sample.name)

            runnable_split_cohort_vep = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_split_cohort_vep(sample_name=sample.name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))
            executable_split_cohort_vep = self.set_stage_runnable(
                stage=stage_split_cohort_vep,
                runnable=runnable_split_cohort_vep)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            executable_split_cohort_vep.dependencies.append(runnable_annotate_cohort_vep.name)

            if use_cache:
                reference_split_cohort_vep = runnable_split_cohort_vep.get_cache_file_path(
                    file_path=self.bwa_genome_db,
                    absolute=True)
            else:
                reference_split_cohort_vep = self.bwa_genome_db

            # Run the GATK SelectVariants analysis to split multi-sample VCF files into one per sample.

            runnable_step = RunnableStepGATK(
                name='split_cohort_vep_gatk_select_variants_vep',
                java_temporary_path=runnable_split_cohort_vep.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx2G',
                java_jar_path=self.java_archive_gatk)
            runnable_split_cohort_vep.add_runnable_step(runnable_step=runnable_step)

            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_vep)
            runnable_step.add_gatk_option(key='variant', value=file_path_annotate_cohort_vep.complete_vcf_bgz)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_vep.sample_vcf)
            runnable_step.add_gatk_option(key='sample_name', value=sample.name)
            runnable_step.add_gatk_switch(key='excludeNonVariants')

            # Run the GATK VariantsToTable analysis.

            if self.variants_to_table:
                runnable_step = RunnableStepGATK(
                    name='split_cohort_vep_gatk_variants_to_table_vep',
                    java_temporary_path=runnable_split_cohort_vep.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    java_jar_path=self.java_archive_gatk)
                runnable_split_cohort_vep.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_vep)
                runnable_step.add_gatk_option(key='variant', value=file_path_split_cohort_vep.sample_vcf)
                runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_vep.sample_tsv)
                runnable_step.add_gatk_switch(key='showFiltered')
                # Set of fixed VCF fields.
                for field_name in variants_to_table_fields['fixed']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK Haplotype Caller-specific INFO fields.
                for field_name in variants_to_table_fields['haplotype_caller']['info']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK Haplotype Caller-specific genotype fields.
                for field_name in variants_to_table_fields['haplotype_caller']['format']:
                    runnable_step.add_gatk_option(key='genotypeFields', value=field_name, override=True)
                # Set of Ensembl VEP-specific INFO fields.
                for field_name in variants_to_table_fields['vep']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)

        ###################################################
        # Step 9: Summarise the variant calling analysis. #
        ###################################################
        #
        # bsfR bsf_variant_calling_summary.R    (summary)

        # file_path_summary = FilePathSummary(prefix=self.get_prefix_summary(cohort_name=self.cohort_name))

        # Create a bsf.procedure.Runnable and bsf.process.Executable for summarising the cohort.

        runnable_summary = self.add_runnable_consecutive(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_summary(cohort_name=self.cohort_name),
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                cache_path_dict=self._cache_path_dict,
                debug=self.debug))
        executable_summary = self.set_stage_runnable(
            stage=stage_summary,
            runnable=runnable_summary)
        # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
        for runnable_diagnose_sample in runnable_diagnose_sample_list:
            executable_summary.dependencies.append(runnable_diagnose_sample.name)

        # Run the bsfR script to summarise the variant calling procedure.

        runnable_step = RunnableStep(
            name='summary',
            program='bsf_variant_calling_summary.R')
        runnable_summary.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_option_long(key='prefix', value=self.get_prefix_summary(cohort_name=self.cohort_name))

        #####################################
        # Step 10: Somatic variant calling. #
        #####################################
        #
        # GATK MuTect2          (somatic_gatk_mutect2)
        # snpEff                (somatic_snpeff)
        # GATK VariantAnnotator (somatic_gatk_variant_annotator)
        # GATK VariantsToTable  (somatic_gatk_variants_to_table)

        for comparison_name in sorted(self._comparison_dict):
            # Somatic variant calling-specific file paths

            file_path_somatic = self.get_file_path_somatic(comparison_name=comparison_name)

            runnable_somatic_gather = run_somatic_scatter_gather(comparison_key=comparison_name)

            # Create a bsf.procedure.Runnable for processing the somatic calls.

            runnable_somatic = self.add_runnable_consecutive(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_somatic(comparison_name=comparison_name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    debug=self.debug))
            executable_somatic = self.set_stage_runnable(
                stage=stage_somatic,
                runnable=runnable_somatic)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            executable_somatic.dependencies.append(runnable_somatic_gather.name)
            # For the moment, this bsf.procedure.Runnable has no bsf.process.RunnableStep objects assigned.
            # It is purely needed for accessing FilePath objects in the report() method below.
            # Never submit the corresponding executable.
            executable_somatic.submit = False

            # reference_somatic = runnable_somatic.get_absolute_cache_file_path(
            #     file_path=self.bwa_genome_db)

            ############################################
            # Step 11: Annotate somatic variant calls. #
            ############################################

            # Run snpEff annotation.

            runnable_annotate_somatic_snpeff = run_annotate_snpeff(
                prefix=self.get_prefix_annotate_somatic_snpeff(comparison_name=comparison_name),
                vcf_file_path=file_path_somatic.somatic_vcf)
            executable_annotate_somatic_snpeff = self.set_stage_runnable(
                stage=stage_annotate_somatic_snpeff,
                runnable=runnable_annotate_somatic_snpeff)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            executable_annotate_somatic_snpeff.dependencies.append(runnable_somatic_gather.name)

            file_path_annotate_somatic_snpeff = self.get_file_path_annotate_somatic_snpeff(
                comparison_name=comparison_name)

            # Run Ensembl Variant Effect Predictor (VEP) annotation.

            runnable_annotate_somatic_vep = run_annotate_vep(
                prefix=self.get_prefix_annotate_somatic_vep(comparison_name=comparison_name),
                vcf_file_path=file_path_somatic.somatic_vcf)
            executable_annotate_somatic_vep = self.set_stage_runnable(
                stage=stage_annotate_somatic_vep,
                runnable=runnable_annotate_somatic_vep)
            # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
            executable_annotate_somatic_vep.dependencies.append(runnable_somatic_gather.name)

            file_path_annotate_somatic_vep = self.get_file_path_annotate_somatic_vep(
                comparison_name=comparison_name)

            #########################################
            # Step 12: Split somatic variant calls. #
            #########################################

            # Split the somatic snpEff-annotated VCF file into a TSV file.

            if self.variants_to_table:
                file_path_split_somatic_snpeff = self.get_file_path_split_somatic_snpeff(
                    comparison_name=comparison_name)

                runnable_split_somatic_snpeff = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_split_somatic_snpeff(comparison_name=comparison_name),
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        debug=self.debug))
                executable_split_somatic_snpeff = self.set_stage_runnable(
                    stage=stage_split_somatic_snpeff,
                    runnable=runnable_split_somatic_snpeff)
                # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
                executable_split_somatic_snpeff.dependencies.append(runnable_annotate_somatic_snpeff.name)

                if use_cache:
                    reference_split_somatic_snpeff = runnable_split_somatic_snpeff.get_cache_file_path(
                        file_path=self.bwa_genome_db,
                        absolute=True)
                else:
                    reference_split_somatic_snpeff = self.bwa_genome_db

                # Run the GATK VariantsToTable analysis.

                runnable_step = RunnableStepGATK(
                    name='somatic_gatk_variants_to_table_snpeff',
                    java_temporary_path=runnable_split_somatic_snpeff.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    java_jar_path=self.java_archive_gatk)
                runnable_split_somatic_snpeff.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_somatic_snpeff)
                runnable_step.add_gatk_option(key='variant', value=file_path_annotate_somatic_snpeff.annotated_vcf)
                runnable_step.add_gatk_option(key='out', value=file_path_split_somatic_snpeff.comparison_tsv)
                runnable_step.add_gatk_switch(key='showFiltered')
                # Set of fixed VCF fields.
                for field_name in variants_to_table_fields['fixed']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK MuTect2 INFO fields.
                for field_name in variants_to_table_fields['mutect2']['info']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK MuTect2 FORMAT fields.
                for field_name in variants_to_table_fields['mutect2']['format']:
                    runnable_step.add_gatk_option(key='genotypeFields', value=field_name, override=True)
                # Set of snpEff INFO fields.
                for field_name in variants_to_table_fields['snpeff']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Automatically add all fields defined for the Variant Annotator resources, above.
                if self.annotation_resources_dict:  # not None and not empty
                    for resource_name, (file_path, annotation_list) in self.annotation_resources_dict.items():
                        if file_path and annotation_list:
                            for annotation in annotation_list:
                                runnable_step.add_gatk_option(
                                    key='fields',
                                    value='.'.join((resource_name, annotation)),
                                    override=True)

            # Split the somatic Ensembl VEP-annotated VCF file into a TSV file.

            if self.variants_to_table:
                file_path_split_somatic_vep = self.get_file_path_split_somatic_vep(comparison_name=comparison_name)

                runnable_split_somatic_vep = self.add_runnable_consecutive(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_split_somatic_vep(comparison_name=comparison_name),
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        debug=self.debug))
                executable_split_somatic_vep = self.set_stage_runnable(
                    stage=stage_split_somatic_vep,
                    runnable=runnable_split_somatic_vep)
                # Set dependencies on preceding bsf.procedure.Runnable.name or bsf.process.Executable.name objects.
                executable_split_somatic_vep.dependencies.append(runnable_annotate_somatic_vep.name)

                if use_cache:
                    reference_split_somatic_vep = runnable_split_somatic_vep.get_cache_file_path(
                        file_path=self.bwa_genome_db,
                        absolute=True)
                else:
                    reference_split_somatic_vep = self.bwa_genome_db

                # Run the GATK VariantsToTable analysis.

                runnable_step = RunnableStepGATK(
                    name='somatic_gatk_variants_to_table_vep',
                    java_temporary_path=runnable_split_somatic_vep.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    java_jar_path=self.java_archive_gatk)
                runnable_split_somatic_vep.add_runnable_step(runnable_step=runnable_step)

                runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_somatic_vep)
                runnable_step.add_gatk_option(key='variant', value=file_path_annotate_somatic_vep.complete_vcf_bgz)
                runnable_step.add_gatk_option(key='out', value=file_path_split_somatic_vep.comparison_tsv)
                runnable_step.add_gatk_switch(key='showFiltered')
                # Set of fixed VCF fields.
                for field_name in variants_to_table_fields['fixed']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK MuTect2 INFO fields.
                for field_name in variants_to_table_fields['mutect2']['info']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)
                # Set of GATK MuTect2 FORMAT fields.
                for field_name in variants_to_table_fields['mutect2']['format']:
                    runnable_step.add_gatk_option(key='genotypeFields', value=field_name, override=True)
                # Set of Ensembl VEP-specific INFO fields.
                for field_name in variants_to_table_fields['vep']:
                    runnable_step.add_gatk_option(key='fields', value=field_name, override=True)

        return

    def report(self):
        """Create a C{bsf.analyses.variant_calling.VariantCallingGATK} report in HTML format and a
        UCSC Genome Browser Track Hub.
        """

        def report_create_directory(path):
            """Private function to create a directory avoiding race conditions.

            @param path: Path
            @type path: str
            @return: Path
            @rtype: str
            """
            if not os.path.isdir(path):
                try:
                    os.makedirs(path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            return path

        def report_create_symbolic_link(source_path, target_path):
            """Private function to set symbolic links.

            @param source_path: Source path
            @type source_path: str
            @param target_path: Target path
            @type target_path: str
            """
            if not os.path.islink(target_path):
                # try:
                #     os.remove(target_path)
                # except OSError as exception:
                #     if exception.errno != errno.ENOENT:
                #         raise
                try:
                    os.symlink(source_path, target_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise

            return

        def report_link():
            """Private function to create simpler symbolic links in structured sub-directories.

            The simplified symbolic links facilitate file system browsing.
            """
            # Create a result directory as concatenation of bsf.analysis.Analysis.genome_version and
            # bsf.analysis.Analysis.prefix.

            directory_results = report_create_directory(
                path=os.path.join(self.project_directory, '_'.join((self.genome_version, self.prefix))))

            directory_results_by_cohort = report_create_directory(
                path=os.path.join(directory_results, 'by_cohort'))

            directory_results_by_pair = report_create_directory(
                path=os.path.join(directory_results, 'by_pair'))

            directory_results_by_sample = report_create_directory(
                path=os.path.join(directory_results, 'by_sample'))

            directory_results_by_type = report_create_directory(
                path=os.path.join(directory_results, 'by_type'))

            # Process per sample.

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

                file_path_process_sample = self.get_file_path_process_sample(sample_name=sample.name)
                # file_path_diagnose_sample = self.get_file_path_diagnose_sample(sample_name=sample.name)
                file_path_annotate_cohort_snpeff = self.get_file_path_annotate_cohort_snpeff(
                    cohort_name=self.cohort_name)
                file_path_annotate_cohort_vep = self.get_file_path_annotate_cohort_vep(cohort_name=self.cohort_name)
                file_path_split_cohort_snpeff = self.get_file_path_split_cohort_snpeff(sample_name=sample.name)
                file_path_split_cohort_vep = self.get_file_path_split_cohort_vep(sample_name=sample.name)

                # Create a sample-specific directory and add symbolic links to it.

                directory_sample = report_create_directory(
                    path=os.path.join(directory_results_by_sample, sample.name))

                for attribute, extension in (
                        ('realigned_bam', '.bam'),
                        ('realigned_bai', '.bam.bai'),
                        ('realigned_md5', '.bam.md5')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_process_sample, attribute)),
                            directory_sample),
                        target_path=os.path.join(directory_sample, sample.name + extension))

                for attribute, extension in (
                        ('sample_vcf', '_snpeff.vcf.gz'),
                        ('sample_tbi', '_snpeff.vcf.gz.tbi'),
                        ('sample_tsv', '_snpeff.tsv')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_cohort_snpeff, attribute)),
                            directory_sample),
                        target_path=os.path.join(directory_sample, sample.name + extension))

                for attribute, extension in (
                        ('sample_vcf', '_vep.vcf.gz'),
                        ('sample_tbi', '_vep.vcf.gz.tbi'),
                        ('sample_tsv', '_vep.tsv')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_cohort_vep, attribute)),
                            directory_sample),
                        target_path=os.path.join(directory_sample, sample.name + extension))

                # Create an alignment-specific directory and add symbolic links to it.

                directory_alignments = report_create_directory(
                    path=os.path.join(directory_results_by_type, 'alignments'))

                for attribute, extension in (
                        ('realigned_bam', '.bam'),
                        ('realigned_bai', '.bam.bai'),
                        ('realigned_md5', '.bam.md5')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_process_sample, attribute)),
                            directory_alignments),
                        target_path=os.path.join(directory_alignments, sample.name + extension))

                # Create a variant-specific directory and add symbolic links to it.

                directory_variants = report_create_directory(
                    path=os.path.join(directory_results_by_type, 'variants'))

                directory_variants_snpeff = report_create_directory(
                    path=os.path.join(directory_variants, 'snpeff'))

                for attribute, extension in (
                        ('sample_vcf', '_snpeff.vcf.gz'),
                        ('sample_tbi', '_snpeff.vcf.gz.tbi'),
                        ('sample_tsv', '_snpeff.tsv')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_cohort_snpeff, attribute)),
                            directory_variants_snpeff),
                        target_path=os.path.join(directory_variants_snpeff, sample.name + extension))

                directory_variants_vep = report_create_directory(
                    path=os.path.join(directory_variants, 'vep'))

                for attribute, extension in (
                        ('sample_vcf', '_vep.vcf.gz'),
                        ('sample_tbi', '_vep.vcf.gz.tbi'),
                        ('sample_tsv', '_vep.tsv')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_cohort_vep, attribute)),
                            directory_variants_vep),
                        target_path=os.path.join(directory_variants_vep, sample.name + extension))

                for attribute, extension in (
                        ('annotated_vcf', '_annotated.vcf.gz'),
                        ('annotated_tbi', '_annotated.vcf.gz.tbi'),
                        ('complete_vcf_bgz', '_snpeff.vcf.gz'),
                        ('complete_vcf_tbi', '_snpeff.vcf.gz.tbi'),
                        ('complete_genes', '_snpeff_summary.genes.txt'),
                        ('complete_stats', '_snpeff_summary.html')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_cohort_snpeff, attribute)),
                            directory_results_by_cohort),
                        target_path=os.path.join(directory_results_by_cohort, self.cohort_name + extension))

                for attribute, extension in (
                        ('complete_vcf_bgz', '_vep.vcf.gz'),
                        ('complete_vcf_tbi', '_vep.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_cohort_vep, attribute)),
                            directory_results_by_cohort),
                        target_path=os.path.join(directory_results_by_cohort, self.cohort_name + extension))

            # Process per (somatic) comparison

            for comparison_name in self._comparison_dict:
                file_path_somatic = self.get_file_path_somatic(comparison_name=comparison_name)

                for attribute, extension in (
                        ('somatic_vcf', '.vcf.gz'),
                        ('somatic_tbi', '.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_somatic, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                file_path_annotate_somatic_snpeff = self.get_file_path_annotate_somatic_snpeff(
                    comparison_name=comparison_name)

                for attribute, extension in (
                        ('complete_vcf_bgz', '_snpeff.vcf.gz'),
                        ('complete_vcf_tbi', '_snpeff.vcf.gz.tbi'),
                        ('complete_genes', '_snpeff_summary.genes.txt'),
                        ('complete_stats', '_snpeff_summary.html'),
                        ('annotated_vcf', '_annotated.vcf.gz'),
                        ('annotated_tbi', '_annotated.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_somatic_snpeff, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                file_path_split_somatic_snpeff = self.get_file_path_split_somatic_snpeff(
                    comparison_name=comparison_name)

                for attribute, extension in (
                        ('comparison_tsv', '_snpeff.tsv'),):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_somatic_snpeff, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                file_path_annotate_somatic_vep = self.get_file_path_annotate_somatic_vep(
                    comparison_name=comparison_name)

                for attribute, extension in (
                        ('complete_vcf_bgz', '_vep.vcf.gz'),
                        ('complete_vcf_tbi', '_vep.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_somatic_vep, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                file_path_split_somatic_vep = self.get_file_path_split_somatic_vep(
                    comparison_name=comparison_name)

                for attribute, extension in (
                        ('comparison_tsv', '_vep.tsv'),):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_somatic_vep, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

            return

        def report_html():
            """Private function to create a HTML report.
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str] """

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.append('<h2 id="genome_browsing">Genome Browsing</h2>\n')
            str_list.append('\n')

            str_list.append('<p id="ucsc_track_hub">')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<h2 id="read_group_and_sample_level">Read Group and Sample Level</h2>\n')
            str_list.append('\n')
            str_list.append('<table id="read_group_and_sample_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Sample</th>\n')
            str_list.append('<th>Variants<br />snpEff<br />Ensembl&nbsp;VEP</th>\n')
            str_list.append('<th>Alignments</th>\n')
            str_list.append('<th>Read Group</th>\n')
            str_list.append('<th>Duplicate Metrics</th>\n')
            str_list.append('<th>Alignment Summary Metrics</th>\n')
            str_list.append('<th>Hybrid Selection Metrics</th>\n')
            str_list.append('<th>Non-Callable Loci</th>\n')
            str_list.append('<th>Non-Callable Summary</th>\n')
            str_list.append('<th>Insert Size</th>\n')
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

                file_path_process_sample = self.get_file_path_process_sample(sample_name=sample.name)
                file_path_diagnose_sample = self.get_file_path_diagnose_sample(sample_name=sample.name)
                file_path_split_cohort_snpeff = self.get_file_path_split_cohort_snpeff(sample_name=sample.name)
                file_path_split_cohort_vep = self.get_file_path_split_cohort_vep(sample_name=sample.name)

                str_list.append('<tr>\n')
                # Sample
                str_list.append('<td class="left">' + sample.name + '</td>\n')
                # Variants
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_split_cohort_snpeff.sample_vcf + '">')
                str_list.append('<abbr title="Variant Calling Format">VCF</abbr>')
                str_list.append('</a>')
                str_list.append('&nbsp;')
                str_list.append('<a href="' + file_path_split_cohort_snpeff.sample_tbi + '">')
                str_list.append('<abbr title="Tabix Index">TBI</abbr>')
                str_list.append('</a>')
                if self.variants_to_table:
                    str_list.append('&nbsp;')
                    str_list.append('<a href="' + file_path_split_cohort_snpeff.sample_tsv + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                str_list.append('<br />')
                str_list.append('<a href="' + file_path_split_cohort_vep.sample_vcf + '">')
                str_list.append('<strong><abbr title="Variant Calling Format">VCF</abbr></strong>')
                str_list.append('</a>')
                str_list.append('&nbsp;')
                str_list.append('<a href="' + file_path_split_cohort_vep.sample_tbi + '">')
                str_list.append('<abbr title="Tabix Index">TBI</abbr>')
                str_list.append('</a>')
                if self.variants_to_table:
                    str_list.append('&nbsp;')
                    str_list.append('<a href="' + file_path_split_cohort_vep.sample_tsv + '">')
                    str_list.append('<strong><abbr title="Tab-Separated Value">TSV</abbr></strong>')
                    str_list.append('</a>')
                str_list.append('</td>\n')
                # Alignments
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_process_sample.realigned_bam + '">')
                str_list.append('<abbr title="Binary Alignment/Map">BAM</abbr>')
                str_list.append('</a>&nbsp;')
                str_list.append('<a href="' + file_path_process_sample.realigned_bai + '">')
                str_list.append('<abbr title="Binary Alignment/Map Index">BAI</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # Read Group
                str_list.append('<td class="left"></td>\n')
                # Duplicate Metrics
                str_list.append('<td class="center">')
                # This can be a sample-specific file or a symbolic link to the read group-specific file.
                if os.path.exists(
                        os.path.join(
                            self.genome_directory,
                            file_path_process_sample.duplicate_metrics)):
                    str_list.append('<a href="' + file_path_process_sample.duplicate_metrics + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                str_list.append('</td>\n')
                # Alignment Summary Metrics
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_process_sample.alignment_summary_metrics + '">')
                str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # Hybrid Selection Metrics
                str_list.append('<td class="center">')
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnose_sample.hybrid_selection_metrics)):
                    str_list.append('<a href="' + file_path_diagnose_sample.hybrid_selection_metrics + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                str_list.append('</td>\n')
                # Non-Callable Loci
                str_list.append('<td class="center">')
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnose_sample.non_callable_loci_tsv)):
                    str_list.append('<a href="' + file_path_diagnose_sample.callable_bed + '">')
                    str_list.append('<abbr title="Browser Extensible Data">BED</abbr>')
                    str_list.append('</a>&nbsp;')
                    str_list.append('<a href="' + file_path_diagnose_sample.callable_bb + '">')
                    str_list.append('<abbr title="Big Browser Extensible Data">BigBED</abbr>')
                    str_list.append('</a>&nbsp;')
                    # Do not link the more complex file_path_diagnose_sample.non_callable_regions_tsv
                    # file for the moment.
                    str_list.append('<a href="' + file_path_diagnose_sample.non_callable_regions_tsv + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                str_list.append('</td>\n')
                # Non-Callable Summary
                str_list.append('<td class="center">')
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnose_sample.non_callable_summary_tsv)):
                    str_list.append('<a href="' + file_path_diagnose_sample.non_callable_summary_tsv + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                str_list.append('</td>\n')
                # Insert size
                str_list.append('<td class="center">')
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnose_sample.insert_size_tsv)):
                    str_list.append('<a href="' + file_path_diagnose_sample.insert_size_pdf + '">')
                    str_list.append('<img alt="Insert Size per Sample ' + sample.name + '"')
                    str_list.append(' src="' + file_path_diagnose_sample.insert_size_png + '"')
                    str_list.append(' height="80" width="80" />')
                    str_list.append('</a>')
                    str_list.append('<a href="' + file_path_diagnose_sample.insert_size_tsv + '">')
                    str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                    str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

                for paired_reads_name in sorted(paired_reads_dict):
                    file_path_process_read_group = self.get_file_path_process_read_group(
                        paired_reads_name=paired_reads_name)

                    str_list.append('<tr>\n')
                    # Sample
                    str_list.append('<td class="left"></td>\n')
                    # Variants
                    str_list.append('<td class="center"></td>\n')
                    # Alignments
                    str_list.append('<td class="center"></td>\n')
                    # Read Group
                    str_list.append('<td class="left">' + paired_reads_name + '</td>\n')
                    # Duplicate Metrics
                    str_list.append('<td class="center">')
                    if os.path.isfile(
                            os.path.join(
                                self.genome_directory,
                                file_path_process_read_group.duplicate_metrics)):
                        str_list.append('<a href="' + file_path_process_read_group.duplicate_metrics + '">')
                        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                        str_list.append('</a>')
                    str_list.append('</td>\n')
                    # Alignment Summary Metrics
                    str_list.append('<td class="center">')
                    str_list.append('</td>\n')
                    # Hybrid Selection Metrics
                    str_list.append('<td class="center"></td>\n')
                    # Non-Callable Loci
                    str_list.append('<td class="center"></td>\n')
                    # Non-Callable Summary
                    str_list.append('<td class="center"></td>\n')
                    # Insert Size
                    str_list.append('<td class="center"></td>\n')
                    str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            str_list.append('<h2 id="cohort_level">Cohort Level</h2>\n')
            str_list.append('\n')
            str_list.append('<table id="cohort_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Cohort</th>\n')
            str_list.append('<th>Information</th>\n')
            str_list.append('<th>Comment</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            file_path_annotate_cohort_snpeff = self.get_file_path_annotate_cohort_snpeff(cohort_name=self.cohort_name)

            str_list.append('<tr>\n')
            str_list.append('<td class="left">' + self.cohort_name + '</td>\n')
            str_list.append('<td class="left">')
            str_list.append('<a href="' + file_path_annotate_cohort_snpeff.gatk_stats + '">')
            str_list.append('snpEff Summary Statistics')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">')
            str_list.append('<a href="' + file_path_annotate_cohort_snpeff.gatk_genes + '">')
            str_list.append('snpEff Summary Genes')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('</tr>\n')

            str_list.append('<tr>\n')
            str_list.append('<td class="left">' + self.cohort_name + '</td>\n')
            str_list.append('<td class="left">')
            str_list.append('snpEff-annotated multi-sample ')
            str_list.append('<a href="' + file_path_annotate_cohort_snpeff.gatk_vcf_bgz + '">')
            str_list.append('<abbr title="Variant Calling Format">VCF</abbr>')
            str_list.append('</a> and ')
            str_list.append('<a href="' + file_path_annotate_cohort_snpeff.gatk_vcf_tbi + '">')
            str_list.append('<abbr title="Tabix Index">TBI</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Functional annotation of all splice variants</td>\n')
            str_list.append('</tr>\n')

            str_list.append('<tr>\n')
            str_list.append('<td class="left">' + self.cohort_name + '</td>\n')
            str_list.append('<td class="left">')
            str_list.append('GATK-annotated multi-sample ')
            str_list.append('<a href="' + file_path_annotate_cohort_snpeff.annotated_vcf + '">')
            str_list.append('<abbr title="Variant Calling Format">VCF</abbr>')
            str_list.append('</a> and ')
            str_list.append('<a href="' + file_path_annotate_cohort_snpeff.annotated_tbi + '">')
            str_list.append('<abbr title="Tabix Index">TBI</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">')
            str_list.append('Functional annotation of only the most severely affected splice variant')
            str_list.append('</td>\n')
            str_list.append('</tr>\n')

            file_path_annotate_cohort_vep = self.get_file_path_annotate_cohort_vep(cohort_name=self.cohort_name)

            str_list.append('<tr>\n')
            str_list.append('<td class="left">' + self.cohort_name + '</td>\n')
            str_list.append('<td class="left">')
            str_list.append('<a href="' + file_path_annotate_cohort_vep.statistics + '">')
            str_list.append('Ensembl Variant Effect Predictor Summary Statistics')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">')
            str_list.append('</td>\n')
            str_list.append('</tr>\n')

            str_list.append('<tr>\n')
            str_list.append('<td class="left">' + self.cohort_name + '</td>\n')
            str_list.append('<td class="left">')
            str_list.append('Ensembl VEP-annotated multi-sample ')
            str_list.append('<a href="' + file_path_annotate_cohort_vep.complete_vcf_bgz + '">')
            str_list.append('<abbr title="Variant Calling Format">VCF</abbr>')
            str_list.append('</a> and ')
            str_list.append('<a href="' + file_path_annotate_cohort_vep.complete_vcf_tbi + '">')
            str_list.append('<abbr title="Tabix Index">TBI</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Functional annotation of all Ensembl splice variants</td>\n')
            str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Somatic variant calling.

            if self._comparison_dict:
                str_list.append('<h2 id="somatic_variants">Somatic Variants</h2>\n')
                str_list.append('\n')
                str_list.append('<table id="somatic_variants_table">\n')
                str_list.append('<thead>\n')
                str_list.append('<tr>\n')
                str_list.append('<th>Comparison</th>\n')
                str_list.append('<th>Variants<br />snpEff<br />VEP</th>\n')
                str_list.append('<th>Summary&nbsp;Statistics<br />snpEff<br />VEP</th>\n')
                str_list.append('</tr>\n')
                str_list.append('</thead>\n')
                str_list.append('<tbody>\n')

                for comparison_name in sorted(self._comparison_dict):
                    # Add snpEff annotation

                    file_path_annotate_somatic_snpeff = self.get_file_path_annotate_somatic_snpeff(
                        comparison_name=comparison_name)

                    file_path_split_somatic_snpeff = self.get_file_path_split_somatic_snpeff(
                        comparison_name=comparison_name)

                    # Add VEP annotation

                    file_path_annotate_somatic_vep = self.get_file_path_annotate_somatic_vep(
                        comparison_name=comparison_name)

                    file_path_split_somatic_vep = self.get_file_path_split_somatic_vep(
                        comparison_name=comparison_name)

                    str_list.append('<tr>\n')

                    # Comparison
                    str_list.append('<td class="left">' + comparison_name + '</td>\n')

                    # Variants
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + file_path_annotate_somatic_snpeff.annotated_vcf + '">')
                    str_list.append('<abbr title="Variant Calling Format">VCF</abbr>')
                    str_list.append('</a>')
                    str_list.append('&nbsp;')
                    str_list.append('<a href="' + file_path_annotate_somatic_snpeff.annotated_tbi + '">')
                    str_list.append('<abbr title="Tabix Index">TBI</abbr>')
                    str_list.append('</a>')
                    if self.variants_to_table:
                        str_list.append('&nbsp;')
                        str_list.append('<a href="' + file_path_split_somatic_snpeff.comparison_tsv + '">')
                        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
                        str_list.append('</a>')
                    str_list.append('<br />')
                    str_list.append('<a href="' + file_path_annotate_somatic_vep.complete_vcf_bgz + '">')
                    str_list.append('<strong><abbr title="Variant Calling Format">VCF</abbr></strong>')
                    str_list.append('</a>')
                    str_list.append('&nbsp;')
                    str_list.append('<a href="' + file_path_annotate_somatic_vep.complete_vcf_tbi + '">')
                    str_list.append('<abbr title="Tabix Index">TBI</abbr>')
                    str_list.append('</a>')
                    if self.variants_to_table:
                        str_list.append('&nbsp;')
                        str_list.append('<a href="' + file_path_split_somatic_vep.comparison_tsv + '">')
                        str_list.append('<strong><abbr title="Tab-Separated Value">TSV</abbr></strong>')
                        str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Summary Statistics
                    str_list.append('<td class="center">')
                    str_list.append('<a href="' + file_path_annotate_somatic_snpeff.gatk_stats + '">')
                    str_list.append('<abbr title="Hyper Text Markup Language">HTML</abbr>')
                    str_list.append('</a>&nbsp;')
                    str_list.append('<a href="' + file_path_annotate_somatic_snpeff.gatk_genes + '">')
                    str_list.append('<abbr title="Text">TXT</abbr>')
                    str_list.append('</a><br />')
                    str_list.append('<a href="' + file_path_annotate_somatic_vep.statistics + '">')
                    str_list.append('<strong><abbr title="Hyper Text Markup Language">HTML</abbr></strong>')
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

            file_path_summary = self.get_file_path_summary(cohort_name=self.cohort_name)

            # Alignment Summary - TSV
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.alignment_metrics_sample_tsv)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.alignment_metrics_sample_tsv + '">TSV</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.alignment_metrics_read_group_tsv + '">TSV</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">')
                str_list.append('<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' +
                                '#AlignmentSummaryMetrics">Alignment Summary</a>')
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            # Alignment Summary - Percent Aligned
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.alignment_percentage_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.alignment_percentage_sample_pdf + '">')
                str_list.append('<img alt="Alignment Summary - Percent Aligned per Sample"')
                str_list.append(' src="' + file_path_summary.alignment_percentage_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.alignment_percentage_read_group_pdf + '">')
                str_list.append('<img alt="Alignment Summary - Percent Aligned per Read Group"')
                str_list.append(' src="' + file_path_summary.alignment_percentage_read_group_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Alignment Summary - Percent Aligned</td>\n')
                str_list.append('</tr>\n')

            # Alignment Summary - Reads Aligned
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.alignment_absolute_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.alignment_absolute_sample_pdf + '">')
                str_list.append('<img alt="Alignment Summary - Reads Aligned per Sample"')
                str_list.append(' src="' + file_path_summary.alignment_absolute_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.alignment_absolute_read_group_pdf + '">')
                str_list.append('<img alt="Alignment Summary - Reads Aligned per Read Group"')
                str_list.append(' src="' + file_path_summary.alignment_absolute_read_group_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Alignment Summary - Reads Aligned</td>\n')
                str_list.append('</tr>\n')

            # Duplication - TSV
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.duplication_metrics_sample_tsv)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.duplication_metrics_sample_tsv + '">TSV</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center"></td>\n')
                str_list.append('<td class="left">')
                str_list.append('<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' +
                                '#DuplicationMetrics">Duplication</a>')
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            # Duplication - Fraction
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.duplication_percentage_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.duplication_percentage_sample_pdf + '">')
                str_list.append('<img alt="Duplication - Duplicated Reads per Sample"')
                str_list.append(' src="' + file_path_summary.duplication_percentage_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center"></td>\n')
                str_list.append('<td class="left">Duplication - Fraction</td>\n')
                str_list.append('</tr>\n')

            # Duplication - Levels
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.duplication_levels_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.duplication_levels_sample_pdf + '">')
                str_list.append('<img alt="Duplication - Duplication Levels per Sample"')
                str_list.append(' src="' + file_path_summary.duplication_levels_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center"></td>\n')
                str_list.append('<td class="left">Duplication - Levels</td>\n')
                str_list.append('</tr>\n')

            # Hybrid Selection - TSV
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.hybrid_metrics_sample_tsv)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_metrics_sample_tsv + '">TSV</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_metrics_read_group_tsv + '">TSV</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">')
                str_list.append('<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' +
                                '#HsMetrics">Hybrid Selection</a>')
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            # Hybrid Selection - Target Coverage Levels
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.hybrid_coverage_levels_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_coverage_levels_sample_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Mean Target Coverage Levels per Sample"')
                str_list.append(' src="' + file_path_summary.hybrid_coverage_levels_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_coverage_levels_read_group_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Mean Target Coverage Levels per Read Group"')
                str_list.append(' src="' + file_path_summary.hybrid_coverage_levels_read_group_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Hybrid Selection - Mean Target Coverage Levels</td>\n')
                str_list.append('</tr>\n')

            # Hybrid Selection - Mean Target Coverage
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.hybrid_coverage_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_coverage_sample_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Mean Target Coverage per Sample"')
                str_list.append(' src="' + file_path_summary.hybrid_coverage_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_coverage_read_group_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Mean Target Coverage per Read Group"')
                str_list.append(' src="' + file_path_summary.hybrid_coverage_read_group_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Hybrid Selection - Mean Target Coverage</td>\n')
                str_list.append('</tr>\n')

            # Hybrid Selection - Excluded Bases
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.hybrid_excluded_bases_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_excluded_bases_sample_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Percent Excluded Bases per Sample"')
                str_list.append(' src="' + file_path_summary.hybrid_excluded_bases_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_excluded_bases_read_group_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Percent Excluded Bases per Read Group"')
                str_list.append(' src="' + file_path_summary.hybrid_excluded_bases_read_group_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Hybrid Selection - Percent Excluded Bases</td>\n')
                str_list.append('</tr>\n')

            # Hybrid Selection - Percent Unique Reads
            # The plot is meaningless if Picard MarkDuplicates has not run.
            if not self.skip_mark_duplicates and os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.hybrid_unique_percentage_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_unique_percentage_sample_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Percent Unique Reads per Sample"')
                str_list.append(' src="' + file_path_summary.hybrid_unique_percentage_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.hybrid_unique_percentage_read_group_pdf + '">')
                str_list.append('<img alt="Hybrid Selection - Percent Unique Reads per Read Group"')
                str_list.append(' src="' + file_path_summary.hybrid_unique_percentage_read_group_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Hybrid Selection - Percent Unique Reads</td>\n')
                str_list.append('</tr>\n')

            # Non-Callable Loci - TSV
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.non_callable_metrics_sample_tsv)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.non_callable_metrics_sample_tsv + '">TSV</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center"></td>\n')
                str_list.append('<td class="left">Non-Callable Loci</td>\n')
                str_list.append('</tr>\n')

            # Non-Callable Loci - Fraction
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.non_callable_percentage_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.non_callable_percentage_sample_pdf + '">')
                str_list.append('<img alt="Non-Callable Loci - Fraction"')
                str_list.append(' src="' + file_path_summary.non_callable_percentage_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center"></td>\n')
                str_list.append('<td class="left">Non-Callable Loci - Fraction</td>\n')
                str_list.append('</tr>\n')

            # Non-Callable Loci - Number
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.non_callable_absolute_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.non_callable_absolute_sample_pdf + '">')
                str_list.append('<img alt="Non-Callable Loci - Number"')
                str_list.append(' src="' + file_path_summary.non_callable_absolute_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center"></td>\n')
                str_list.append('<td class="left">Non-Callable Loci - Number</td>\n')
                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

        def report_hub():
            """Private function to create a UCSC Track Hub.
            """

            str_list = list()
            """ @type str_list: list[str] """

            # Group via UCSC super tracks.

            str_list.append('track Alignments\n')
            str_list.append('shortLabel Alignments\n')
            str_list.append('longLabel BWA NGS read alignments\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Callable\n')
            str_list.append('shortLabel Callable\n')
            str_list.append('longLabel Callable\n')
            str_list.append('superTrack on\n')
            str_list.append('\n')

            str_list.append('track Variants\n')
            str_list.append('shortLabel Variants\n')
            str_list.append('longLabel Variant calls\n')
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

                file_path_process_sample = self.get_file_path_process_sample(sample_name=sample.name)
                file_path_diagnose_sample = self.get_file_path_diagnose_sample(sample_name=sample.name)
                file_path_split_cohort_snpeff = self.get_file_path_split_cohort_snpeff(sample_name=sample.name)
                file_path_split_cohort_vep = self.get_file_path_split_cohort_vep(sample_name=sample.name)

                #  Alignments track
                #
                # Common track settings
                str_list.append('track ' + sample.name + '_alignments\n')
                str_list.append('type bam\n')
                str_list.append('shortLabel ' + sample.name + '_alignments\n')
                str_list.append('longLabel ' + sample.name + 'BWA NGS read alignments\n')
                str_list.append('bigDataUrl ' + file_path_process_sample.realigned_bam + '\n')
                # str_list.append('html ...\n')
                str_list.append('visibility squish\n')
                # Common optional track settings
                str_list.append('color 0,0,0\n')
                # bam/cram - Compressed Sequence Alignment track settings
                # Supertrack settings
                str_list.append('parent Alignments\n')
                str_list.append('\n')

                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnose_sample.callable_bb)):
                    #  Callable track
                    #
                    # Common track settings
                    str_list.append('track ' + sample.name + '_callable\n')
                    str_list.append('type bigBed\n')
                    str_list.append('shortLabel ' + sample.name + '_callable\n')
                    str_list.append('longLabel ' + sample.name + 'callable\n')
                    str_list.append('bigDataUrl ' + file_path_diagnose_sample.callable_bb + '\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility squish\n')
                    # Common optional track settings
                    str_list.append('color 0,0,0\n')
                    # bigBed - Item or region track settings
                    # Supertrack settings
                    str_list.append('parent Callable\n')
                    str_list.append('\n')

                # snpEff Variants track
                #
                # Common track settings
                str_list.append('track ' + sample.name + '_snpeff\n')
                str_list.append('type vcfTabix\n')
                str_list.append('shortLabel ' + sample.name + '_snpeff\n')
                str_list.append('longLabel ' + sample.name + ' snpEff-annotated variant calls\n')
                str_list.append('bigDataUrl ' + file_path_split_cohort_snpeff.sample_vcf + '\n')
                # str_list.append('html ...\n')
                str_list.append('visibility dense\n')
                # Common optional track settings
                # vcfTabix - Variant Call Format (indexed by tabix) track settings
                # Supertrack settings
                str_list.append('parent Variants\n')
                str_list.append('\n')

                # Ensembl VEP Variants track
                #
                # Common track settings
                str_list.append('track ' + sample.name + '_vep\n')
                str_list.append('type vcfTabix\n')
                str_list.append('shortLabel ' + sample.name + '_vep\n')
                str_list.append('longLabel ' + sample.name + ' Ensembl VEP-annotated variant calls\n')
                str_list.append('bigDataUrl ' + file_path_split_cohort_vep.sample_vcf + '\n')
                # str_list.append('html ...\n')
                str_list.append('visibility dense\n')
                # Common optional track settings
                # vcfTabix - Variant Call Format (indexed by tabix) track settings
                # Supertrack settings
                str_list.append('parent Variants\n')
                str_list.append('\n')

            # Comparison-specific tracks

            self.ucsc_hub_to_file(content=str_list)

            return

        report_link()
        report_html()
        report_hub()

        return
