"""bsf.analyses.variant_calling

A package of classes and methods supporting variant calling analyses.
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
import math
import os
from pickle import Pickler, HIGHEST_PROTOCOL
import warnings

from bsf import Analysis, FilePath, Runnable
from bsf.annotation import AnnotationSheet
from bsf.executables import BWA
from bsf.process import Command, Executable, RunnableStep, RunnableStepJava, RunnableStepPicard, RunnableStepLink, \
    RunnableStepMove
from bsf.standards import Configuration, Default

import pysam


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
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_maximum=None,
            java_jar_path=None,
            gatk_classpath=None):
        """Initialise a C{bsf.analyses.variant_calling.RunnableStepGATK}.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{bsf.process.Executable} into the C{bsf.Stage}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepGATK, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path,
            java_heap_maximum=java_heap_maximum,
            java_jar_path=java_jar_path)

        # Set the GATK classpath and the GATK Java archive.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(gatk_classpath, 'GenomeAnalysisTK.jar'))

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
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_option_long(key=key, value=value, override=override)

    def add_gatk_switch(self, key):
        """Add a C{bsf.argument.SwitchLong} to a C{bsf.analyses.variant_calling.RunnableStepGATK}.

        @param key: Option key
        @type key: str
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_switch_long(key=key)


class FilePathAlignment(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathAlignment} models files in a sample-specific directory.

    Attributes:
    @ivar aligned_bam: Alignment BAM file path
    @type aligned_bam: str | unicode
    @ivar aligned_bai: Alignment BAI file path
    @type aligned_bai: str | unicode
    @ivar aligned_md5: Alignment MD5 check sum file path
    @type aligned_md5: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAlignment} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathAlignment, self).__init__(prefix=prefix)

        self.aligned_bam = prefix + '.bam'
        self.aligned_bai = prefix + '.bai'
        self.aligned_md5 = prefix + '.bam.md5'

        return


class FilePathProcessReadGroup(FilePath):
    """The C{bsf.analyses.variant_calling.FilePathAlignment} models files in a sample-specific directory.

    Attributes:
    @ivar duplicates_marked_bam: Picard Mark Duplicates BAM file
    @type duplicates_marked_bam: str | unicode
    @ivar duplicates_marked_bai: Picard Mark Duplicates BAI file
    @type duplicates_marked_bai: str | unicode
    @ivar duplicates_marked_md5: Picard Mark Duplicates MD5 check sum file
    @type duplicates_marked_md5: str | unicode
    @ivar duplicate_metrics: Picard Mark Duplicates metrics file
    @type duplicate_metrics: str | unicode
    @ivar realigner_targets: GATK
    @type realigner_targets: str | unicode
    @ivar realigned_bam: Re-aligned BAM file
    @type realigned_bam: str | unicode
    @ivar realigned_bai: Re-aligned BAI file
    @type realigned_bai: str | unicode
    @ivar realigned_md5: Re-aligned MD5 check sum file
    @type realigned_md5: str | unicode
    @ivar recalibration_table_pre:
    @type recalibration_table_pre: str | unicode
    @ivar recalibration_table_post:
    @type recalibration_table_post: str | unicode
    @ivar recalibration_plot:
    @type recalibration_plot: str | unicode
    @ivar recalibrated_bam: Recalibrated BAM file
    @type recalibrated_bam: str | unicode
    @ivar recalibrated_bai: Recalibrated BAI file
    @type recalibrated_bai: str | unicode
    @ivar recalibrated_md5: Recalibrated BAM MD5 check sum file
    @type recalibrated_md5: str | unicode
    @ivar alignment_summary_metrics: Picard Collect Alignment Summary metrics file
    @type alignment_summary_metrics: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathProcessReadGroup} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathProcessReadGroup, self).__init__(prefix=prefix)

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
        self.alignment_summary_metrics = prefix + '_alignment_summary_metrics.tsv'

        return


class FilePathProcessSample(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathProcessSample} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
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
        self.realigned_md5 = prefix + '_realigned_bam.md5'
        self.realigned_bam_bai = prefix + '_realigned.bam.bai'
        self.alignment_summary_metrics = prefix + '_alignment_summary_metrics.tsv'
        self.raw_variants_gvcf_vcf = prefix + '_raw_variants.g.vcf.gz'
        self.raw_variants_gvcf_idx = prefix + '_raw_variants.g.vcf.gz.tbi'

        return


class FilePathDiagnoseSample(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathDiagnoseSample} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathDiagnoseSample, self).__init__(prefix=prefix)

        self.diagnose_targets_vcf = prefix + '_diagnose_targets.vcf.gz'
        self.diagnose_targets_idx = prefix + '_diagnose_targets.vcf.gz.tbi'
        self.missing_intervals = prefix + '_missing.intervals'
        self.missing_report = prefix + '_missing.gatkreport'
        self.callable_bed = prefix + '_callable_loci.bed'
        self.callable_txt = prefix + '_callable_loci.txt'
        self.callable_bb = prefix + '_callable_loci.bb'
        self.sorted_bed = prefix + '_callable_sorted.bed'
        # Defined in bsf_variant_calling_coverage.R.
        self.non_callable_loci_tsv = prefix + '_non_callable_loci.tsv'
        self.non_callable_regions_tsv = prefix + '_non_callable_regions.tsv'
        self.non_callable_summary_tsv = prefix + '_non_callable_summary.tsv'
        self.hybrid_selection_metrics = prefix + '_hybrid_selection_metrics.tsv'

        return


class FilePathMergeCohort(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathMergeCohort} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathMergeCohort, self).__init__(prefix=prefix)

        self.combined_gvcf_vcf = prefix + '_combined.g.vcf.gz'
        self.combined_gvcf_tbi = prefix + '_combined.g.vcf.gz.tbi'
        self.partial_gvcf_vcf = prefix + '_partial.g.vcf.gz'
        self.partial_gvcf_tbi = prefix + '_partial.g.vcf.gz.tbi'

        return


class FilePathGenotypeCohort(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathGenotypeCohort} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathGenotypeCohort, self).__init__(prefix=prefix)

        self.genotyped_raw_vcf = prefix + '_genotyped_raw_snp_raw_indel.vcf.gz'
        self.genotyped_raw_tbi = prefix + '_genotyped_raw_snp_raw_indel.vcf.gz.tbi'
        # self.combined_vcf = prefix + '_combined.vcf.gz'
        # self.combined_tbi = prefix + '_combined.vcf.gz.tbi'
        self.partial_vcf = prefix + '_partial.vcf.gz'
        self.partial_tbi = prefix + '_partial.vcf.gz.tbi'

        return


class FilePathProcessCohort(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathProcessCohort} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathProcessCohort, self).__init__(prefix=prefix)

        # self.genotyped_raw_vcf = prefix + '_genotyped_raw_snp_raw_indel.vcf.gz'
        # self.genotyped_raw_tbi = prefix + '_genotyped_raw_snp_raw_indel.vcf.gz.tbi'
        self.recalibrated_snp_raw_indel_vcf = prefix + '_recalibrated_snp_raw_indel.vcf.gz'
        self.recalibrated_snp_raw_indel_idx = prefix + '_recalibrated_snp_raw_indel.vcf.gz.tbi'
        self.recalibrated_snp_recalibrated_indel_vcf = prefix + '_recalibrated_snp_recalibrated_indel.vcf.gz'
        self.recalibrated_snp_recalibrated_indel_idx = prefix + '_recalibrated_snp_recalibrated_indel.vcf.gz.tbi'
        self.multi_sample_vcf = prefix + '_multi_sample.vcf.gz'
        self.multi_sample_idx = prefix + '_multi_sample.vcf.gz.tbi'
        self.snpeff_vcf = prefix + '_snpeff.vcf'
        self.snpeff_idx = prefix + '_snpeff.vcf.idx'
        self.snpeff_vcf_bgz = prefix + '_snpeff.vcf.gz'
        self.snpeff_vcf_tbi = prefix + '_snpeff.vcf.gz.tbi'
        self.snpeff_stats = prefix + '_snpeff_summary.html'
        self.annotated_vcf = prefix + '_annotated.vcf.gz'
        self.annotated_tbi = prefix + '_annotated.vcf.gz.tbi'
        self.recalibration_indel = prefix + '_recalibration_indel.recal'
        self.recalibration_snp = prefix + '_recalibration_snp.recal'
        self.tranches_indel = prefix + '_recalibration_indel.tranches'
        self.tranches_snp = prefix + '_recalibration_snp.tranches'
        self.plots_indel = prefix + '_recalibration_indel.R'
        self.plots_snp = prefix + '_recalibration_snp.R'

        return


class FilePathAnnotateCohort(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAnnotateCohort} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathAnnotateCohort, self).__init__(prefix=prefix)

        self.vep_statistics = prefix + '_vep_statistics.html'
        # Complete VEP set raw
        self.vep_complete_raw_vcf = prefix + '_vep_complete_raw.vcf'
        self.vep_complete_raw_vcf_bgz = prefix + '_vep_complete_raw.vcf.gz'
        self.vep_complete_raw_vcf_tbi = prefix + '_vep_complete_raw.vcf.gz.tbi'
        # Filtered VEP set raw
        self.vep_filtered_raw_vcf = prefix + '_vep_filtered_raw.vcf'
        self.vep_filtered_raw_vcf_bgz = prefix + '_vep_filtered_raw.vcf.gz'
        self.vep_filtered_raw_vcf_tbi = prefix + '_vep_filtered_raw.vcf.gz.tbi'
        # VCF.Filter-converted set
        self.vep_complete_vcf_bgz = prefix + '_vep_complete.vcf.gz'
        self.vep_complete_vcf_tbi = prefix + '_vep_complete.vcf.gz.tbi'
        self.vep_filtered_vcf_bgz = prefix + '_vep_filtered.vcf.gz'
        self.vep_filtered_vcf_tbi = prefix + '_vep_filtered.vcf.gz.tbi'

        return


class FilePathSplitCohort(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSplitCohort} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathSplitCohort, self).__init__(prefix=prefix)

        self.sample_vcf = prefix + '.vcf.gz'
        self.sample_idx = prefix + '.vcf.gz.tbi'
        self.sample_tsv = prefix + '.tsv'

        return


class FilePathSummary(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSummary} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
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
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSomatic} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathSomatic, self).__init__(prefix=prefix)

        self.somatic_vcf = prefix + '_somatic.vcf.gz'
        self.somatic_idx = prefix + '_somatic.vcf.gz.tbi'
        self.snpeff_vcf = prefix + '_snpeff.vcf'
        self.snpeff_idx = prefix + '_snpeff.vcf.idx'
        self.snpeff_vcf_bgz = prefix + '_snpeff.vcf.gz'
        self.snpeff_vcf_tbi = prefix + '_snpeff.vcf.gz.tbi'
        self.snpeff_stats = prefix + '_snpeff_summary.html'
        self.snpeff_genes = prefix + '_snpeff_summary.genes.txt'
        self.annotated_vcf = prefix + '_annotated.vcf.gz'
        self.annotated_tbi = prefix + '_annotated.vcf.gz.tbi'
        self.annotated_tsv = prefix + '_annotated.tsv'

        return


class FilePathSomaticScatterGather(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSomaticScatterGather} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathSomaticScatterGather, self).__init__(prefix=prefix)

        self.partial_vcf = prefix + '_somatic.vcf.gz'
        self.partial_tbi = prefix + '_somatic.vcf.gz.tbi'

        return


class VariantCallingGATKComparison(object):
    """C{bsf.analyses.variant_calling.VariantCallingGATKComparison} class representing a somatic comparison.

    Attributes:
    @ivar normal_sample: Normal sample
    @type normal_sample: bsf.ngs.Sample
    @ivar tumor_sample: Tumour sample
    @type tumor_sample: bsf.ngs.Sample
    @ivar panel_of_normal_path: File path to a Panel-Of-Normal (PON) VCF file.
    @type panel_of_normal_path: str | unicode
    """

    def __init__(
            self,
            normal_sample=None,
            tumor_sample=None,
            panel_of_normal_path=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATKComparison} object.

        @param normal_sample: Normal sample
        @type normal_sample: bsf.ngs.Sample
        @param tumor_sample: Tumour sample
        @type tumor_sample: bsf.ngs.Sample
        @param panel_of_normal_path: File path to a Panel-Of-Normal (PON) VCF file.
        @type panel_of_normal_path: str | unicode
        """

        self.normal_sample = normal_sample  # Can be None.
        self.tumor_sample = tumor_sample  # Can be None.
        self.panel_of_normal_path = panel_of_normal_path  # Cen be None.

        return

    @property
    def get_name(self):
        """Get the name of a C{bsf.analyses.variant_calling.VariantCallingGATKComparison}.

        @return: Comparison name
        @rtype: str
        """
        name = str()

        if self.normal_sample is not None:
            name += self.normal_sample.name

        name += '__'

        if self.tumor_sample is not None:
            name += self.tumor_sample.name

        return name.strip('_')


class VariantCallingGATK(Analysis):
    """C{bsf.analyses.variant_calling.VariantCallingGATK} class representing the logic to run the
    Genome Analysis Toolkit (GATK).

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_align_lane: C{bsf.Stage.name} for the lane alignment stage
    @type stage_name_align_lane: str
    @cvar stage_name_process_lane: C{bsf.Stage.name} for the lane processing stage
    @type stage_name_process_lane: str
    @cvar stage_name_process_sample: C{bsf.Stage.name} for the sample processing stage
    @type stage_name_process_sample: str
    @cvar stage_name_diagnose_sample: C{bsf.Stage.name} for the sample diagnosis stage
    @type stage_name_diagnose_sample: str
    @cvar stage_name_process_cohort: C{bsf.Stage.name} for the cohort processing stage
    @type stage_name_process_cohort: str
    @cvar stage_name_annotate_cohort: C{bsf.Stage.name} for the cohort annotation stage
    @type stage_name_annotate_cohort: str
    @cvar stage_name_split_cohort: C{bsf.Stage.name} for the cohort splitting stage
    @type stage_name_split_cohort: str
    @cvar stage_name_summary: C{bsf.Stage.name} for the summary stage
    @type stage_name_summary: str
    @cvar stage_name_somatic: C{bsf.Stage.name} for the somatic stage
    @type stage_name_somatic: str
    @ivar replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects for processing or run them separately
    @type replicate_grouping: bool
    @ivar bwa_genome_db: Genome sequence file path with BWA index
    @type bwa_genome_db: str | unicode
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode
    @ivar cohort_name: Cohort name
    @type cohort_name: str
    @ivar accessory_cohort_gvcfs: Python C{list} of Python C{str} or C{unicode} (GVCF file path) objects
    @type accessory_cohort_gvcfs: list[str | unicode]
    @ivar skip_mark_duplicates: Skip the Picard MarkDuplicates step
    @type skip_mark_duplicates: bool
    @ivar skip_indel_realignment: Skip the GATK RealignerTargetCreator and GATK IndelRealigner steps
    @type skip_indel_realignment: bool
    @ivar known_sites_discovery: VCF file path for variant discovery via The Haplotype Caller or Unified Genotyper
    @type known_sites_discovery: str | unicode
    @ivar known_sites_realignment: Python C{list} of Python C{str} or C{unicode} (VCF file paths)
        for realignment
    @type known_sites_realignment: list[str | unicode]
    @ivar known_sites_recalibration: Python C{list} of Python C{str} or C{unicode} (VCF file paths)
        for recalibration
    @type known_sites_recalibration: list[str | unicode]
    @ivar known_somatic_discovery: Cosmic VCF file path for somatic variant discovery via MuTect2
    @type known_somatic_discovery: list[str | unicode]
    @ivar annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
        Python C{tuple} of
        Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
    @type annotation_resources_dict: dict[str, (str | unicode, list[str])]
    @ivar truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
    @type truth_sensitivity_filter_level_indel: str
    @ivar truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
    @type truth_sensitivity_filter_level_snp: str
    @ivar vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
    @type vqsr_skip_indel: bool, None
    @ivar vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
    @type vqsr_skip_snp: bool, None
    @ivar vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
    @type vqsr_resources_indel_dict: dict[str, dict[str, str | unicode]]
    @ivar vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
    @type vqsr_resources_snp_dict: dict[str, dict[str, str | unicode]]
    @ivar vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_indel_list: list[str]
    @ivar vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_snp_list: list[str]
    @ivar vqsr_bad_lod_cutoff_indel: LOD score cutoff for negative training set for INDELs
    @type vqsr_bad_lod_cutoff_indel: float
    @ivar vqsr_bad_lod_cutoff_snp: LOD score cutoff for negative training set for SNPs
    @type vqsr_bad_lod_cutoff_snp: float
    @ivar vqsr_max_gaussians_pos_indel: Maximum number of Gaussians in the positive training for INDELs
    @type vqsr_max_gaussians_pos_indel: int
    @ivar vqsr_max_gaussians_pos_snp: Maximum number of Gaussians in the positive training for SNPs
    @type vqsr_max_gaussians_pos_snp: int
    @ivar exclude_intervals_list: Python C{list} of Python C{str} (intervals) to exclude from the analysis
    @type exclude_intervals_list: list[str]
    @ivar include_intervals_list: Python C{list} of Python C{str} (intervals) to include in the analysis
    @type include_intervals_list: list[str]
    @ivar interval_padding: Interval padding
    @type interval_padding: int
    @ivar number_of_tiles_cohort: Number of genomic tiles for scattering in stage variant_calling_process_cohort
    @type number_of_tiles_cohort: int
    @ivar number_of_chunks_cohort: Number of chunks for gathering in stage variant_calling_process_cohort
    @type number_of_chunks_cohort: int
    @ivar number_of_tiles_somatic: Number of genomic tiles for scattering in stage variant_calling_somatic
    @type number_of_tiles_somatic: int
    @ivar number_of_chunks_somatic: Number of chunks for gathering in stage variant_calling_somatic
    @type number_of_chunks_somatic: int
    @ivar downsample_to_fraction: Down-sample to fraction
    @type downsample_to_fraction: str
    @ivar gatk_bundle_version: GATK resource bundle version
    @type gatk_bundle_version: str
    @ivar snpeff_genome_version: snpEff genome version
    @type snpeff_genome_version: str
    @ivar genome_annotation_gtf: Genome annotation Gene Transfer Format (GTF) file path
    @type genome_annotation_gtf: str | unicode
    @ivar classpath_gatk: Genome Analysis Tool Kit Java Archive (JAR) class path directory
    @type classpath_gatk: str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar classpath_snpeff: snpEff tool Java Archive (JAR) class path directory
    @type classpath_snpeff: str | unicode
    @ivar classpath_vcf_filter: VCF.Filter tool Java Archive (JAR) class path directory
    @type classpath_vcf_filter: str | unicode
    """

    name = 'Variant Calling Analysis'
    prefix = 'variant_calling'

    stage_name_align_lane = '_'.join((prefix, 'align_lane'))
    stage_name_process_lane = '_'.join((prefix, 'process_lane'))
    stage_name_process_sample = '_'.join((prefix, 'process_sample'))
    stage_name_diagnose_sample = '_'.join((prefix, 'diagnose_sample'))
    stage_name_merge_cohort = '_'.join((prefix, 'merge_cohort'))
    stage_name_process_cohort = '_'.join((prefix, 'process_cohort'))
    stage_name_annotate_cohort = '_'.join((prefix, 'annotate_cohort'))
    stage_name_split_cohort = '_'.join((prefix, 'split_cohort'))
    stage_name_summary = '_'.join((prefix, 'summary'))
    stage_name_somatic = '_'.join((prefix, 'somatic'))

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
            bwa_genome_db=None,
            comparison_path=None,
            cohort_name=None,
            accessory_cohort_gvcfs=None,
            skip_mark_duplicates=False,
            skip_indel_realignment=False,
            known_sites_discovery=None,
            known_sites_realignment=None,
            known_sites_recalibration=None,
            known_somatic_discovery=None,
            annotation_resources_dict=None,
            truth_sensitivity_filter_level_indel=None,
            truth_sensitivity_filter_level_snp=None,
            vqsr_skip_indel=False,
            vqsr_skip_snp=False,
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
            number_of_tiles_cohort=None,
            number_of_chunks_cohort=None,
            number_of_tiles_somatic=None,
            number_of_chunks_somatic=None,
            downsample_to_fraction=None,
            gatk_bundle_version=None,
            snpeff_genome_version=None,
            genome_annotation_gtf=None,
            classpath_gatk=None,
            classpath_picard=None,
            classpath_snpeff=None,
            classpath_vcf_filter=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATK}.

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
        @param comparisons: Python C{dict} of Python C{str} key and Python C{list} value objects of
            C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, list[bsf.ngs.Sample]]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects
            for processing or run them separately
        @type replicate_grouping: bool
        @param bwa_genome_db: Genome sequence file path with BWA index
        @type bwa_genome_db: str | unicode
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode
        @param cohort_name: Cohort name
        @type cohort_name: str
        @param accessory_cohort_gvcfs: Python C{list} of Python C{str} or C{unicode} (GVCF file path) objects
        @type accessory_cohort_gvcfs: list[str | unicode]
        @param skip_mark_duplicates: Skip the Picard MarkDuplicates step
        @type skip_mark_duplicates: bool
        @param skip_indel_realignment: Skip the GATK RealignerTargetCreator and GATK IndelRealigner steps
        @type skip_indel_realignment: bool
        @param known_sites_discovery: VCF file path for variant discovery via The Haplotype Caller or Unified Genotyper
        @type known_sites_discovery: str | unicode
        @param known_sites_realignment: Python C{list} of Python C{str} or C{unicode} (VCF file paths)
            for realignment
        @type known_sites_realignment: list[str | unicode]
        @param known_sites_recalibration: Python C{list} of Python C{str} or C{unicode} (VCF file paths)
            for recalibration
        @type known_sites_recalibration: list[str | unicode]
        @param known_somatic_discovery: Cosmic VCF file path for somatic variant discovery via MuTect2
        @type known_somatic_discovery: list[str | unicode]
        @param annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
            Python C{tuple} of
            Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
        @type annotation_resources_dict: dict[str, (str | unicode, list[str])]
        @param truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
        @type truth_sensitivity_filter_level_indel: str
        @param truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
        @type truth_sensitivity_filter_level_snp: str
        @param vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
        @type vqsr_skip_indel: bool, None
        @param vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
        @type vqsr_skip_snp: bool, None
        @param vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_indel_dict: dict[str, dict[str, str | unicode]]
        @param vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_snp_dict: dict[str, dict[str, str | unicode]]
        @param vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_indel_list: list[str]
        @param vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_snp_list: list[str]
        @param vqsr_bad_lod_cutoff_indel: LOD score cutoff for negative training set for INDELs
        @type vqsr_bad_lod_cutoff_indel: float
        @param vqsr_bad_lod_cutoff_snp: LOD score cutoff for negative training set for SNPs
        @type vqsr_bad_lod_cutoff_snp: float
        @param vqsr_max_gaussians_pos_indel: Maximum number of Gaussians in the positive training for INDELs
        @type vqsr_max_gaussians_pos_indel: int
        @param vqsr_max_gaussians_pos_snp: Maximum number of Gaussians in the positive training for SNPs
        @type vqsr_max_gaussians_pos_snp: int
        @param exclude_intervals_list: Python C{list} of Python C{str} (intervals) to exclude from the analysis
        @type exclude_intervals_list: list[str]
        @param include_intervals_list: Python C{list} of Python C{str} (intervals) to include in the analysis
        @type include_intervals_list: list[str]
        @param interval_padding: Interval padding
        @type interval_padding: int
        @param number_of_tiles_cohort: Number of genomic tiles for scattering in stage variant_calling_process_cohort
        @type number_of_tiles_cohort: int
        @param number_of_chunks_cohort: Number of chunks for gathering in stage variant_calling_process_cohort
        @type number_of_chunks_cohort: int
        @param number_of_tiles_somatic: Number of genomic tiles for scattering in stage variant_calling_somatic
        @type number_of_tiles_somatic: int
        @param number_of_chunks_somatic: Number of chunks for gathering in stage variant_calling_somatic
        @type number_of_chunks_somatic: int
        @param downsample_to_fraction: Down-sample to fraction
        @type downsample_to_fraction: str
        @param gatk_bundle_version: GATK resource bundle version
        @type gatk_bundle_version: str
        @param snpeff_genome_version: snpEff genome version
        @type snpeff_genome_version: str
        @param genome_annotation_gtf: Genome annotation Gene Transfer Format (GTF) file path
        @type genome_annotation_gtf: str | unicode
        @param classpath_gatk: Genome Analysis Tool Kit Java Archive (JAR) class path directory
        @type classpath_gatk: str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param classpath_snpeff: snpEff tool Java Archive (JAR) class path directory
        @type classpath_snpeff: str | unicode
        @param classpath_vcf_filter: VCF.Filter tool Java Archive (JAR) class path directory
        @type classpath_vcf_filter: str | unicode
        @return:
        @rtype:
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
            comparisons=comparisons,
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if bwa_genome_db is None:
            self.bwa_genome_db = str()
        else:
            self.bwa_genome_db = bwa_genome_db

        if comparison_path is None:
            self.comparison_path = str()
        else:
            self.comparison_path = comparison_path

        if cohort_name is None:
            self.cohort_name = str()
        else:
            self.cohort_name = cohort_name

        if accessory_cohort_gvcfs is None:
            self.accessory_cohort_gvcfs = list()
        else:
            self.accessory_cohort_gvcfs = accessory_cohort_gvcfs

        if skip_mark_duplicates is None:
            self.skip_mark_duplicates = False
        else:
            assert isinstance(skip_mark_duplicates, bool)
            self.skip_mark_duplicates = skip_mark_duplicates

        if skip_indel_realignment is None:
            self.skip_indel_realignment = False
        else:
            assert isinstance(skip_indel_realignment, bool)
            self.skip_indel_realignment = skip_indel_realignment

        if known_sites_discovery is None:
            self.known_sites_discovery = str()
        else:
            self.known_sites_discovery = known_sites_discovery

        if known_sites_realignment is None:
            self.known_sites_realignment = list()
        else:
            self.known_sites_realignment = known_sites_realignment

        if known_sites_recalibration is None:
            self.known_sites_recalibration = list()
        else:
            self.known_sites_recalibration = known_sites_recalibration

        if known_somatic_discovery is None:
            self.known_somatic_discovery = list()
        else:
            self.known_somatic_discovery = known_somatic_discovery

        if annotation_resources_dict is None:
            self.annotation_resources_dict = dict()
        else:
            self.annotation_resources_dict = annotation_resources_dict

        if truth_sensitivity_filter_level_indel is None:
            self.truth_sensitivity_filter_level_indel = str()
        else:
            self.truth_sensitivity_filter_level_indel = truth_sensitivity_filter_level_indel

        if truth_sensitivity_filter_level_snp is None:
            self.truth_sensitivity_filter_level_snp = str()
        else:
            self.truth_sensitivity_filter_level_snp = truth_sensitivity_filter_level_snp

        if vqsr_skip_indel is None:
            self.vqsr_skip_indel = False
        else:
            assert isinstance(vqsr_skip_indel, bool)
            self.vqsr_skip_indel = vqsr_skip_indel

        if vqsr_skip_snp is None:
            self.vqsr_skip_snp = False
        else:
            assert isinstance(vqsr_skip_snp, bool)
            self.vqsr_skip_snp = vqsr_skip_snp

        if vqsr_resources_indel_dict is None:
            self.vqsr_resources_indel_dict = dict()
        else:
            self.vqsr_resources_indel_dict = vqsr_resources_indel_dict

        if vqsr_resources_snp_dict is None:
            self.vqsr_resources_snp_dict = dict()
        else:
            self.vqsr_resources_snp_dict = vqsr_resources_snp_dict

        if vqsr_annotations_indel_list is None:
            self.vqsr_annotations_indel_list = list()
        else:
            self.vqsr_annotations_indel_list = vqsr_annotations_indel_list

        if vqsr_annotations_snp_list is None:
            self.vqsr_annotations_snp_list = list()
        else:
            self.vqsr_annotations_snp_list = vqsr_annotations_snp_list

        self.vqsr_bad_lod_cutoff_indel = vqsr_bad_lod_cutoff_indel  # Can be None.
        self.vqsr_bad_lod_cutoff_snp = vqsr_bad_lod_cutoff_snp  # Can be None.
        self.vqsr_max_gaussians_pos_indel = vqsr_max_gaussians_pos_indel  # Can be None.
        self.vqsr_max_gaussians_pos_snp = vqsr_max_gaussians_pos_snp  # Can be None.

        if exclude_intervals_list is None:
            self.exclude_intervals_list = list()
        else:
            self.exclude_intervals_list = exclude_intervals_list

        if include_intervals_list is None:
            self.include_intervals_list = list()
        else:
            self.include_intervals_list = include_intervals_list

        if interval_padding is None:
            self.interval_padding = int(x=0)
        else:
            assert isinstance(interval_padding, int)
            self.interval_padding = interval_padding

        if number_of_tiles_cohort is None:
            self.number_of_tiles_cohort = 0
        else:
            assert isinstance(number_of_tiles_cohort, int)
            self.number_of_tiles_cohort = number_of_tiles_cohort

        if number_of_chunks_cohort is None:
            self.number_of_chunks_cohort = 0
        else:
            assert isinstance(number_of_chunks_cohort, int)
            self.number_of_chunks_cohort = number_of_chunks_cohort

        if number_of_tiles_somatic is None:
            self.number_of_tiles_somatic = 0
        else:
            assert isinstance(number_of_tiles_somatic, int)
            self.number_of_tiles_somatic = number_of_tiles_somatic

        if number_of_chunks_somatic is None:
            self.number_of_chunks_somatic = 0
        else:
            assert isinstance(number_of_chunks_somatic, int)
            self.number_of_chunks_somatic = number_of_chunks_somatic

        if downsample_to_fraction is None:
            self.downsample_to_fraction = str()
        else:
            self.downsample_to_fraction = downsample_to_fraction

        if gatk_bundle_version is None:
            self.gatk_bundle_version = str()
        else:
            self.gatk_bundle_version = gatk_bundle_version

        if snpeff_genome_version is None:
            self.snpeff_genome_version = str()
        else:
            self.snpeff_genome_version = snpeff_genome_version

        if genome_annotation_gtf is None:
            self.genome_annotation_gtf = str()
        else:
            self.genome_annotation_gtf = genome_annotation_gtf

        if classpath_gatk is None:
            self.classpath_gatk = str()
        else:
            self.classpath_gatk = classpath_gatk

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if classpath_snpeff is None:
            self.classpath_snpeff = str()
        else:
            self.classpath_snpeff = classpath_snpeff

        if classpath_vcf_filter is None:
            self.classpath_vcf_filter = str()
        else:
            self.classpath_vcf_filter = classpath_vcf_filter

        self._cache_path_dict = None
        """ @type _cache_path_dict: dict[str, str | unicode] """

        # Initialise the Python list of genome tile regions with an empty region to run a single process by default.

        self._tile_region_cohort_list = [[('', 0, 0)]]
        """ @type _tile_region_cohort_list: list[list[(str, int, int)]] """

        self._tile_region_somatic_list = [[('', 0, 0)]]
        """ @type _tile_region_somatic_list: list[list[(str, int, int)]] """

        return

    @property
    def get_gatk_bundle_path(self):
        """Get the absolute GATK bundle directory C{bsf.standards.Default.absolute_gatk_bundle} for the set
        C{bsf.analyses.variant_calling.VariantCallingGATK.gatk_bundle_version} and
        C{bsf.analyses.variant_calling.VariantCallingGATK.genome_version}.

        @return: Absolute GATK bundle directory
        @rtype: str | unicode
        """
        return Default.absolute_gatk_bundle(
            gatk_bundle_version=self.gatk_bundle_version,
            genome_version=self.genome_version)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.variant_calling.VariantCallingGATK} via a
        C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        def set_vqsr_configuration(vqsr_resources_dict, variation_type):
            """Private function to read variant quality score recalibration (VQSR) configuration information.

            Configuration options I{vqsr_resources_indel} and I{vqsr_resources_snp} provide a comma-separated list of
            resources that are to be used in the VQSR procedure. Each option needs to correspond to a sub-section of
            the C{ConfigParser.SafeConfigParser} in C{bsf.standards.Configuration.config_parser}.
            Each sub-section needs options 'known', 'training', 'truth', 'prior' and 'file_path'.
            @param vqsr_resources_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
            @type vqsr_resources_dict: dict[str, dict[str, str | unicode]]
            @param variation_type: Variation type I{indel} or I{snp}
            @type variation_type: str
            @return:
            @rtype:
            """

            assert isinstance(vqsr_resources_dict, dict)
            assert isinstance(variation_type, str)

            if variation_type not in ('indel', 'snp'):
                raise Exception("Variation type has to be 'indel' or 'snp', not {!r}.".format(variation_type))

            # The vqsr_resources_indel|snp options of the current configuration section hold a comma-separated list
            # of resources that should correspond to a sub-section in the configuration file.
            vqsr_option = '_'.join(('vqsr_resources', variation_type))
            if config_parser.has_option(section=section, option=vqsr_option):
                # Split the resource list on a comma, split white space characters and remove remaining empty strings.
                for resource_key in filter(
                        lambda x: x != '',
                        map(
                            lambda x: x.strip(),
                            config_parser.get(section=section, option=vqsr_option).split(','))):
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
                                    'Missing configuration option {!r} in section {!r}'.format(
                                        resource_option,
                                        resource_section))
                    else:
                        raise Exception(
                            'Missing configuration section {!r} declared in option {!r} {!r}.'.format(
                                resource_section,
                                vqsr_option,
                                config_parser.get(section=section, option=vqsr_option)))

            return

        def set_annotation_configuration(annotation_resources_dict):
            """Private function to read variant annotation configuration information.

            @param annotation_resources_dict: Python dict of Python str (annotation resource name) key and
                Python tuple of Python str (file path) and Python list of Python str (annotation) value data
            @type annotation_resources_dict: dict
            @return:
            @rtype:
            """

            annotation_option = 'annotation_resources'
            if config_parser.has_option(section=section, option=annotation_option):
                # Split the resource list on a comma, split white space characters and remove remaining empty strings.
                for resource_key in filter(
                        lambda x: x != '',
                        map(
                            lambda x: x.strip(),
                            config_parser.get(section=section, option=annotation_option).split(','))):
                    # The annotation resource section consists of section.annotation_resource.
                    resource_section = '.'.join((section, '_'.join(('annotation', resource_key))))
                    if config_parser.has_section(section=resource_section):
                        resource_option = 'file_path'
                        if config_parser.has_option(section=resource_section, option=resource_option):
                            file_path = config_parser.get(section=resource_section, option=resource_option)
                        else:
                            raise Exception(
                                "Missing configuration option {!r} in configuration section {!r}.".format(
                                    resource_option,
                                    resource_section))
                        resource_option = 'annotations'
                        if config_parser.has_option(section=resource_section, option=resource_option):
                            # Split the annotation list on a comma, split white space characters and
                            # remove remaining empty strings.
                            annotation_list = filter(
                                lambda x: x != '',
                                map(
                                    lambda x: x.strip(), config_parser.get(
                                        section=resource_section,
                                        option=resource_option).split(',')))
                        else:
                            raise Exception(
                                "Missing configuration option {!r} in configuration section {!r}.".format(
                                    resource_option,
                                    resource_section))
                        # Create a dict key and a tuple of a Python str and Python list.
                        annotation_resources_dict[resource_key] = file_path, annotation_list
                    else:
                        raise Exception(
                            'Missing configuration section {!r} declared in option {!r} {!r}.'.format(
                                resource_section,
                                annotation_option,
                                config_parser.get(section=section, option=annotation_option)))

            return

        # Start of set_configuration() method body.

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        super(VariantCallingGATK, self).set_configuration(configuration=configuration, section=section)

        config_parser = configuration.config_parser

        option = 'replicate_grouping'
        if config_parser.has_option(section=section, option=option):
            self.replicate_grouping = config_parser.getboolean(section=section, option=option)

        # Get the genome database.

        option = 'bwa_genome_db'
        if config_parser.has_option(section=section, option=option):
            self.bwa_genome_db = config_parser.get(section=section, option=option)

        # Read a comparison file.

        option = 'cmp_file'
        if config_parser.has_option(section=section, option=option):
            self.comparison_path = config_parser.get(section=section, option=option)

        # Get the cohort name.

        option = 'cohort_name'
        if config_parser.has_option(section=section, option=option):
            self.cohort_name = config_parser.get(section=section, option=option)

        # Comma-separated list of GVCF files from accessory cohorts
        # that should be used in the recalibration procedure.

        option = 'accessory_cohort_gvcfs'
        if config_parser.has_option(section=section, option=option):
            self.accessory_cohort_gvcfs = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Get the skip mark duplicates option.

        option = 'skip_mark_duplicates'
        if config_parser.has_option(section=section, option=option):
            self.skip_mark_duplicates = config_parser.getboolean(section=section, option=option)

        # Get the skip INDEL realignment option.

        option = 'skip_indel_realignment'
        if config_parser.has_option(section=section, option=option):
            self.skip_indel_realignment = config_parser.getboolean(section=section, option=option)

        # Get the truth sensitivity filter level for INDELs.

        option = 'truth_sensitivity_filter_level_indel'
        if config_parser.has_option(section=section, option=option):
            self.truth_sensitivity_filter_level_indel = config_parser.get(section=section, option=option)

        # Get the truth sensitivity filter level for SNPs.

        option = 'truth_sensitivity_filter_level_snp'
        if config_parser.has_option(section=section, option=option):
            self.truth_sensitivity_filter_level_snp = config_parser.get(section=section, option=option)

        # Get the flag for skipping the Variant Quality Score Recalibration (VQSR) for INDELs.

        option = 'vqsr_skip_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_skip_indel = config_parser.getboolean(section=section, option=option)

        # Get the flag for skipping the Variant Quality Score Recalibration (VQSR) for SNPs.

        option = 'vqsr_skip_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_skip_snp = config_parser.getboolean(section=section, option=option)

        # Get the list of annotations for the Variant Quality Score Recalibration (VQSR) for INDELs.

        option = 'vqsr_annotations_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_annotations_indel_list = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Get the list of annotations for the Variant Quality Score Recalibration (VQSR) for SNPs.

        option = 'vqsr_annotations_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_annotations_snp_list = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Get the bad LOD cutoff for the negative training set for INDELs.

        option = 'vqsr_bad_lod_cutoff_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_bad_lod_cutoff_indel = config_parser.getfloat(section=section, option=option)

        # Get the bad LOD cutoff for the negative training set for SNPs.

        option = 'vqsr_bad_lod_cutoff_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_bad_lod_cutoff_snp = config_parser.getfloat(section=section, option=option)

        # Get the maximum number of Gaussians in the positive training for INDELs.

        option = 'vqsr_max_gaussians_pos_indel'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_max_gaussians_pos_indel = config_parser.getint(section=section, option=option)

        # Get the maximum number of Gaussians in the positive training for SNPs.

        option = 'vqsr_max_gaussians_pos_snp'
        if config_parser.has_option(section=section, option=option):
            self.vqsr_max_gaussians_pos_snp = config_parser.getint(section=section, option=option)

        # Set VQSR resources and corresponding configuration sections for INDELs.

        set_vqsr_configuration(vqsr_resources_dict=self.vqsr_resources_indel_dict, variation_type='indel')

        # Set VQSR resources and corresponding configuration sections for SNPs.

        set_vqsr_configuration(vqsr_resources_dict=self.vqsr_resources_snp_dict, variation_type='snp')

        # Set additionally requested annotation resources for the GATK AnnotateVariants step.

        set_annotation_configuration(annotation_resources_dict=self.annotation_resources_dict)

        # Single VCF file of known sites for the
        # GATK HaplotypeCaller and GenotypeGVCFs steps.

        option = 'known_sites_discovery'
        if config_parser.has_option(section=section, option=option):
            self.known_sites_discovery = config_parser.get(section=section, option=option)

        # Comma-separated list of VCF files with known variant sites for the
        # GATK RealignerTargetCreator and IndelRealigner steps.

        option = 'known_sites_realignment'
        if config_parser.has_option(section=section, option=option):
            self.known_sites_realignment = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Comma-separated list of VCF files with known variant sites for the
        # GATK BaseRecalibrator and PrintReads steps.

        option = 'known_sites_recalibration'
        if config_parser.has_option(section=section, option=option):
            self.known_sites_recalibration = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Comma-separated list of VCF files with known somatic variant sites for the
        # GATK MuTect2 steps.

        option = 'known_somatic_discovery'
        if config_parser.has_option(section=section, option=option):
            self.known_somatic_discovery = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Get the list of intervals to exclude.

        option = 'exclude_intervals'
        if config_parser.has_option(section=section, option=option):
            self.exclude_intervals_list = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Get the list of intervals to include.

        option = 'include_intervals'
        if config_parser.has_option(section=section, option=option):
            self.include_intervals_list = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    config_parser.get(section=section, option=option).split(',')))

        # Get the interval padding.

        option = 'interval_padding'
        if config_parser.has_option(section=section, option=option):
            self.interval_padding = config_parser.getint(section=section, option=option)

        # Get the number of tiles for variant_calling_process_cohort.

        option = 'number_of_tiles_cohort'
        if config_parser.has_option(section=section, option=option):
            self.number_of_tiles_cohort = config_parser.getint(section=section, option=option)

        # Get the number of chunks for variant_calling_process_cohort.

        option = 'number_of_chunks_cohort'
        if config_parser.has_option(section=section, option=option):
            self.number_of_chunks_cohort = config_parser.getint(section=section, option=option)

        # Get the number of tiles for variant_calling_somatic.

        option = 'number_of_tiles_cohort'
        if config_parser.has_option(section=section, option=option):
            self.number_of_tiles_cohort = config_parser.getint(section=section, option=option)

        # Get the number of chunks for variant_calling_somatic.

        option = 'number_of_chunks_cohort'
        if config_parser.has_option(section=section, option=option):
            self.number_of_chunks_cohort = config_parser.getint(section=section, option=option)

        # Get the down-sample to fraction information.

        option = 'downsample_to_fraction'
        if config_parser.has_option(section=section, option=option):
            self.downsample_to_fraction = config_parser.get(section=section, option=option)

        # Get the GATK bundle version.

        option = 'gatk_bundle_version'
        if config_parser.has_option(section=section, option=option):
            self.gatk_bundle_version = config_parser.get(section=section, option=option)

        # Get the snpEff genome version.

        option = 'snpeff_genome_version'
        if config_parser.has_option(section=section, option=option):
            self.snpeff_genome_version = config_parser.get(section=section, option=option)

        # Get the genome annotation Gene Transfer Format (GTF) file path.

        option = 'genome_annotation_gtf'
        if config_parser.has_option(section=section, option=option):
            self.genome_annotation_gtf = config_parser.get(section=section, option=option)

        # Get the Genome Analysis Tool Kit (GATK) Java Archive (JAR) class path directory.

        option = 'classpath_gatk'
        if config_parser.has_option(section=section, option=option):
            self.classpath_gatk = config_parser.get(section=section, option=option)

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if config_parser.has_option(section=section, option=option):
            self.classpath_picard = config_parser.get(section=section, option=option)

        # Get the snpEff tool Java Archive (JAR) class path directory.

        option = 'classpath_snpeff'
        if config_parser.has_option(section=section, option=option):
            self.classpath_snpeff = config_parser.get(section=section, option=option)

        # Get the VCF.Filter tool Java Archive (JAR) class path directory.

        option = 'classpath_vcf_filter'
        if config_parser.has_option(section=section, option=option):
            self.classpath_vcf_filter = config_parser.get(section=section, option=option)

        return

    def _read_comparisons(self, comparison_path):
        """Read a C{bsf.annotation.AnnotationSheet} CSV file from disk.

            - Column headers for CASAVA folders:
                - Normal/Tumor ProcessedRunFolder:
                    - CASAVA processed run folder name or
                    - C{bsf.Analysis.input_directory} by default
                - Normal/Tumor Project:
                    - CASAVA Project name or
                    - C{bsf.Analysis.project_name} by default
                - Normal/Tumor Sample:
                    - CASAVA Sample name, no default
            - Column headers for independent samples:
                - Normal/Tumor Sample:
                - Normal/Tumor Reads:
                - Normal/Tumor File:
            - PON Path:
                - File path to a Panel-Of-Normal (PON) VCF file
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode
        @return:
        @rtype:
        """

        assert isinstance(comparison_path, (str, unicode))

        # For variant calling, all samples need adding to the Analysis regardless.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        if not comparison_path:
            return

        annotation_sheet = AnnotationSheet.from_file_path(file_path=comparison_path, name='Somatic Comparisons')

        for row_dict in annotation_sheet.row_dicts:
            comparison = VariantCallingGATKComparison()

            if self.debug > 0:
                print "Comparison sheet row_dict {!r}".format(row_dict)

            for prefix in ('Normal', 'Tumor'):
                group_name, group_samples = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                if group_name and len(group_samples):
                    if len(group_samples) != 1:
                        raise Exception("Got more than one Sample for comparison {!r}".format(row_dict))

                    if prefix == 'Normal':
                        comparison.normal_sample = group_samples[0]
                    if prefix == 'Tumor':
                        comparison.tumor_sample = group_samples[0]

            prefix = 'PON Path'
            if prefix in row_dict and row_dict[prefix]:
                comparison.panel_of_normal_path = row_dict[prefix]

            # At least a tumor Sample has to be defined for the "comparison" to make sense.
            if comparison.tumor_sample is not None:
                self.comparisons[comparison.get_name] = comparison

        return

    def run(self):
        """Run a C{bsf.analyses.variant_calling.VariantCallingGATK} analysis.

        @return:
        @rtype:
        """

        def run_create_genome_tiles(tiles=1, width=0):
            """Private function to create genomic tiles for scattering and a partitioned list of indices for gathering.

            The tiles are created on the basis of a Picard sequence dictionary accompanying the genome FASTA file.
            If both tiles and width are 0, then the sequence regions serve as natural tiles.
            @param tiles: Number of tiles for scattering
            @type tiles: int
            @param width: Tile width for scattering
            @type width: int
            @return: Python C{list} of Python C{list} (tiles) of Python C{tuple} (tile region) of
                Python C{str} (sequence region), Python C{int} (start) and Python C{int} (end)
            @rtype: list[list[(str, int, int)]]
            """
            tile_region_list = list()
            """ @type tile_region_list: list[list[(str, int, int)]] """
            total_length = 0
            """ @type total_length: int """

            dict_path = os.path.splitext(self.bwa_genome_db)[0] + '.dict'
            if not os.path.exists(dict_path):
                raise Exception("Picard sequence dictionary {!r} does not exist.".format(dict_path))

            alignment_file = pysam.AlignmentFile(dict_path, 'r')
            # Summarise sequence lengths to get the total length.
            for sq_entry in alignment_file.header['SQ']:
                """ @type sq_entry: dict"""
                total_length += int(sq_entry['LN'])

            if tiles:
                tile_length = float(total_length) / float(tiles)
            elif width:
                tile_length = float(width)
            else:
                # The intervals are just the natural sequence regions.
                # Thus the start coordinate is always 1 and the end coordinate is the sequence length (@SQ SL).
                # 1 2 3 ... 7 8 9
                # start = 1
                # end = 9
                # length = end - start + 1 = 9 - 1 + 1 = 9
                for sq_entry in alignment_file.header['SQ']:
                    tile_region_list.append([(str(sq_entry['SN']), 1, int(sq_entry['LN']))])

                return tile_region_list

            current_length = 0.0
            """ @type current_length: float """
            sq_list = list()
            """ @type sq_list: list[(str, int, int)] """
            for sq_entry in alignment_file.header['SQ']:
                sq_start = 0.0
                """ @type sq_start: float """
                sq_length = float(sq_entry['LN'])
                """ @type sq_length: float """
                while sq_start < sq_length:  # float
                    # The sequence end is the minimum of the sequence start plus remaining tile length or
                    # the sequence length.
                    sq_end = min(sq_start + tile_length - current_length, sq_length)
                    """ @type seq_end: float """
                    sq_list.append((sq_entry['SN'], int(math.floor(sq_start + 1.0)), int(math.floor(sq_end))))
                    current_length += sq_end - sq_start
                    sq_start = sq_end

                    if math.floor(current_length) >= math.floor(tile_length):
                        # If a tile is complete, append the sequence list to the interval list and reset both
                        # list and length.
                        tile_region_list.append(sq_list)
                        sq_list = []
                        current_length = 0.0

            if len(sq_list):
                tile_region_list.append(sq_list)

            return tile_region_list

        def run_merge_cohort_scatter_gather(analysis_stage, cohort_runnable_dict, cohort_name):
            """Private method to hierarchically merge samples into cohorts using a scatter gather approach.

            This method merges GVCF file paths from individual process_sample Runnable objects or an accessory cohort.
            @param analysis_stage: C{Analysis} C{Stage}
            @type analysis_stage: bsf.Stage
            @param cohort_runnable_dict: Python C{dict} of Python C{str} key and Python C{list} of
                C{Runnable} object value data
            @type cohort_runnable_dict: dict[str, list[Runnable]]
            @param cohort_name: Cohort name to select a Python list of C{Runnable} objects from the
                I{cohort_runnable_dict} Python C{dict}
            @type cohort_name: str
            @return: Final C{Runnable} of the gather stage
            @rtype: Runnable
            """

            # Private variables are prefixed with an underscore to avoid clashes with variables in the run() method.

            prefix_merge_cohort = '_'.join((analysis_stage.name, cohort_name))

            _file_path_merge_cohort = FilePathMergeCohort(prefix=prefix_merge_cohort)

            # If the cohort index file already exists, create the Runnable objects, but do not submit their
            # corresponding Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, _file_path_merge_cohort.combined_gvcf_tbi)):
                _submit_runnable = False
            else:
                _submit_runnable = True

            # The cohort_object_list contains either Runnable objects from the process_sample stage or
            # Python str | unicode (GVCF file path) objects for accessory cohorts to be merged.
            cohort_object_list = cohort_runnable_dict[cohort_name]

            vc_runnable_merge_cohort_scatter_list = list()
            """ @type vc_runnable_merge_cohort_scatter_list: list[Runnable] """

            # Scatter by the number of genomic tiles.
            # Create a Runnable and Executable for each GATK CombineGVCFs from its Sample objects.
            runnable_merge_cohort_scatter = None
            for _tile_index in range(0, len(self._tile_region_cohort_list)):
                prefix_merge_cohort_scatter = '_'.join((analysis_stage.name, cohort_name, 'scatter', str(_tile_index)))

                file_path_merge_cohort_scatter = FilePathMergeCohort(prefix=prefix_merge_cohort_scatter)

                runnable_merge_cohort_scatter = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_merge_cohort_scatter,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        file_path_object=file_path_merge_cohort_scatter,
                        debug=self.debug))
                executable_merge_cohort_scatter = self.set_stage_runnable(
                    stage=analysis_stage,
                    runnable=runnable_merge_cohort_scatter)
                # Submit the Executable only, if the cohort index file does not exist.
                executable_merge_cohort_scatter.submit = _submit_runnable
                for cohort_component in cohort_object_list:
                    # Set dependencies on preceding Runnable.name or Executable.name objects.
                    # Set them only for Runnable objects, but not for Python str | unicode (file path) objects.
                    if isinstance(cohort_component, Runnable):
                        executable_merge_cohort_scatter.dependencies.append(cohort_component.name)

                vc_runnable_merge_cohort_scatter_list.append(runnable_merge_cohort_scatter)

                reference_merge_cohort_scatter = runnable_merge_cohort_scatter.get_absolute_cache_file_path(
                    file_path=self.bwa_genome_db)

                _runnable_step = runnable_merge_cohort_scatter.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='merge_cohort_gatk_combine_gvcfs',
                        java_temporary_path=runnable_merge_cohort_scatter.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        gatk_classpath=self.classpath_gatk))
                """ @type _runnable_step: RunnableStepGATK """
                _runnable_step.add_gatk_option(key='analysis_type', value='CombineGVCFs')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_merge_cohort_scatter)
                for _interval in self.exclude_intervals_list:
                    _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                # for _interval in self.include_intervals_list:
                #     _runnable_step.add_gatk_option(key='intervals', value=_interval, override=True)
                # if self.interval_padding:
                #     _runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for _region_tuple in self._tile_region_cohort_list[_tile_index]:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if _region_tuple[0]:
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value='{:s}:{:d}-{:d}'.format(_region_tuple[0], _region_tuple[1], _region_tuple[2]),
                            override=True)
                for cohort_component in cohort_object_list:
                    if isinstance(cohort_component, Runnable):
                        # Variant GVCF file paths are stored under
                        # the 'raw_variants_gvcf_vcf' attribute for the 'process_sample' stage or under
                        # the 'combined_gvcf_vcf' attribute for the 'merge_cohort' stage.
                        variant_path = ''
                        for variant_key in ('raw_variants_gvcf_vcf', 'combined_gvcf_vcf'):
                            try:
                                variant_path = getattr(cohort_component.file_path_object, variant_key)
                            except AttributeError:
                                pass
                        _runnable_step.add_gatk_option(key='variant', value=variant_path, override=True)
                    elif isinstance(cohort_component, (str, unicode)):
                        _runnable_step.add_gatk_option(key='variant', value=cohort_component, override=True)
                    else:
                        raise Exception('Unexpected object type on merge_cohort list.')
                _runnable_step.add_gatk_option(
                    key='out',
                    value=file_path_merge_cohort_scatter.partial_gvcf_vcf)

            # If there is only one tile, no need to gather, just rename file paths and return the (scatter) Runnable.

            if len(self._tile_region_cohort_list) == 1:
                file_path_merge_cohort_scatter = runnable_merge_cohort_scatter.file_path_object
                """ @type file_path_merge_cohort_scatter: FilePathMergeCohort """
                # Add the final, cohort-specific keys to the FilePathMergeCohort for the final (scatter) Runnable.
                file_path_merge_cohort_scatter.combined_gvcf_vcf = _file_path_merge_cohort.combined_gvcf_vcf
                file_path_merge_cohort_scatter.combined_gvcf_tbi = _file_path_merge_cohort.combined_gvcf_tbi

                runnable_merge_cohort_scatter.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='merge_cohort_gather_move_vcf',
                        source_path=file_path_merge_cohort_scatter.partial_gvcf_vcf,
                        target_path=file_path_merge_cohort_scatter.combined_gvcf_vcf))

                runnable_merge_cohort_scatter.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='merge_cohort_gather_move_tbi',
                        source_path=file_path_merge_cohort_scatter.partial_gvcf_tbi,
                        target_path=file_path_merge_cohort_scatter.combined_gvcf_tbi))

                return runnable_merge_cohort_scatter

            # Second, gather by the number of chunks on the partitioned genome tile index list.

            # Second, Merge chunks hierarchically.
            # Initialise a list of Runnable objects and indices for the hierarchical merge.
            vc_runnable_merge_cohort_gather_list = vc_runnable_merge_cohort_scatter_list
            vc_merge_cohort_gather_index_list = range(0, len(self._tile_region_cohort_list))
            # Global variable to keep and return the last Runnable.
            runnable_merge_cohort_gather = None
            """ @type runnable_merge_cohort_gather: Runnable """
            _level = 0
            while len(vc_merge_cohort_gather_index_list) > 1:
                _temporary_gather_runnable_list = list()
                """ @type _temporary_gather_runnable_list: list[Runnable] """
                _temporary_gather_index_list = list()
                """ @type _temporary_gather_index_list: list[int] """
                # Partition the index list into chunks of given size.
                _partition_list = [vc_merge_cohort_gather_index_list[_offset:_offset + self.number_of_chunks_cohort]
                                   for _offset in range(0,
                                                        len(vc_merge_cohort_gather_index_list),
                                                        self.number_of_chunks_cohort)]

                for _partition_index in range(0, len(_partition_list)):
                    _chunk_index_list = _partition_list[_partition_index]
                    # The file prefix includes the level and partition index.
                    prefix_merge_cohort_gather = '_'.join(
                        (analysis_stage.name, cohort_name, 'gather', str(_level), str(_partition_index)))

                    file_path_merge_cohort_gather = FilePathMergeCohort(prefix=prefix_merge_cohort_gather)

                    runnable_merge_cohort_gather = self.add_runnable(
                        runnable=Runnable(
                            name=prefix_merge_cohort_gather,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            cache_directory=self.cache_directory,
                            cache_path_dict=self._cache_path_dict,
                            file_path_object=file_path_merge_cohort_gather,
                            debug=self.debug))
                    executable_merge_cohort_gather = self.set_stage_runnable(
                        stage=analysis_stage,
                        runnable=runnable_merge_cohort_gather)
                    # Submit the Executable only, if the cohort index file does not exist.
                    executable_merge_cohort_gather.submit = _submit_runnable
                    # Dependencies on scatter processes are set based on genome tile indices below.
                    _temporary_gather_runnable_list.append(runnable_merge_cohort_gather)
                    _temporary_gather_index_list.append(_partition_index)

                    reference_merge_cohort_gather = runnable_merge_cohort_gather.get_absolute_cache_file_path(
                        file_path=self.bwa_genome_db)

                    # GATK CatVariants bypasses the GATK engine and thus requires a completely different command line.
                    _runnable_step = runnable_merge_cohort_gather.add_runnable_step(
                        runnable_step=RunnableStepJava(
                            name='merge_cohort_gatk_cat_variants',
                            sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                            java_temporary_path=runnable_merge_cohort_gather.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx4G'))
                    """ @type _runnable_step: RunnableStepJava """
                    _runnable_step.add_option_short(
                        key='classpath',
                        value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                    _sub_command = _runnable_step.sub_command
                    # Add the 'reference' not 'reference_sequence' option.
                    _sub_command.add_option_long(
                        key='reference',
                        value=reference_merge_cohort_gather)
                    _sub_command.add_option_long(
                        key='outputFile',
                        value=file_path_merge_cohort_gather.partial_gvcf_vcf)
                    _sub_command.add_switch_long(key='assumeSorted')
                    # Finally, process per chunk index.
                    for _chunk_index in _chunk_index_list:
                        _runnable_object = vc_runnable_merge_cohort_gather_list[_chunk_index]
                        _file_path_object = _runnable_object.file_path_object
                        """ @type _file_path_object: FilePathMergeCohort """
                        # Set GATK option variant
                        _sub_command.add_option_long(
                            key='variant',
                            value=_file_path_object.partial_gvcf_vcf,
                            override=True)
                        # Delete the *.g.vcf.gz file.
                        _runnable_step.obsolete_file_path_list.append(_file_path_object.partial_gvcf_vcf)
                        # Delete the *.g.vcf.gz.tbi file.
                        _runnable_step.obsolete_file_path_list.append(_file_path_object.partial_gvcf_tbi)
                        # Set dependencies on preceding Runnable.name or Executable.name objects.
                        # Depend on the Runnable.name of the corresponding Runnable of the scattering above.
                        executable_merge_cohort_gather.dependencies.append(_runnable_object.name)

                # Set the temporary index list as the new list and increment the merge level.
                vc_runnable_merge_cohort_gather_list = _temporary_gather_runnable_list
                vc_merge_cohort_gather_index_list = _temporary_gather_index_list
                _level += 1
            else:
                # For the last instance, additionally rename the final file.
                file_path_merge_cohort_gather = runnable_merge_cohort_gather.file_path_object
                """ @type file_path_merge_cohort_gather: FilePathMergeCohort """

                # Add cohort-specific keys to the file path dictionary.
                file_path_merge_cohort_gather.combined_gvcf_vcf = _file_path_merge_cohort.combined_gvcf_vcf
                file_path_merge_cohort_gather.combined_gvcf_tbi = _file_path_merge_cohort.combined_gvcf_tbi

                runnable_merge_cohort_gather.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='merge_cohort_gather_move_vcf',
                        source_path=file_path_merge_cohort_gather.partial_gvcf_vcf,
                        target_path=file_path_merge_cohort_gather.combined_gvcf_vcf))

                runnable_merge_cohort_gather.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='merge_cohort_gather_move_tbi',
                        source_path=file_path_merge_cohort_gather.partial_gvcf_tbi,
                        target_path=file_path_merge_cohort_gather.combined_gvcf_tbi))

            return runnable_merge_cohort_gather

        # Start of the run() method body.

        super(VariantCallingGATK, self).run()

        # Get global defaults.

        default = Default.get_global_default()

        # VariantCallingGATK requires a genome version, which gets configured by the super-class.

        if not self.genome_version:
            raise Exception("A 'VariantCallingGATK' analysis requires a 'genome_version' configuration option.")

        if not self.bwa_genome_db:
            raise Exception("A 'VariantCallingGATK' analysis requires a 'bwa_genome_db' configuration option.")

        if not self.cohort_name:
            self.cohort_name = self.project_name  # The cohort_name used to default to just 'default'.

        if not self.gatk_bundle_version:
            raise Exception("A 'VariantCallingGATK' analysis requires a 'gatk_bundle_version' configuration option.")

        if not self.snpeff_genome_version:
            raise Exception("A 'VariantCallingGATK' analysis requires a 'snpeff_genome_version' configuration option.")

        if not self.classpath_gatk:
            self.classpath_gatk = default.classpath_gatk

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        if not self.classpath_snpeff:
            self.classpath_snpeff = default.classpath_snpeff

        if not self.classpath_vcf_filter:
            self.classpath_vcf_filter = default.classpath_vcf_filter

        # Check for absolute paths and adjust if required before checking for existence.

        self.bwa_genome_db = Default.get_absolute_path(
            file_path=self.bwa_genome_db,
            default_path=self.get_gatk_bundle_path)
        if not os.path.exists(path=self.bwa_genome_db):
            raise Exception("The bwa_genome_db file {!r} does not exist.".format(self.bwa_genome_db))

        # GATK does a lot of read requests from the reference FASTA file.
        # Place it and the accompanying *.fasta.fai and *.dict files in the cache directory.
        self._cache_path_dict = {
            'reference_fasta': self.bwa_genome_db,
            'reference_fai': self.bwa_genome_db + '.fai',
            'reference_dict': os.path.splitext(self.bwa_genome_db)[0] + '.dict'
        }

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for file_path in self.accessory_cohort_gvcfs:
            file_path = Default.get_absolute_path(file_path=file_path, default_path=Default.absolute_projects())
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The accessory_cohort_gvcf file {!r} does not exist.'.format(file_path))
                # TODO: Check the cohorts so that their sample names do not clash.
        self.accessory_cohort_gvcfs = temporary_list

        for key in self.annotation_resources_dict.keys():
            file_path, annotation_list = self.annotation_resources_dict[key]
            file_path = Default.get_absolute_path(file_path=file_path, default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                self.annotation_resources_dict[key] = file_path, annotation_list
            else:
                raise Exception('The file path {!r} for annotation resource {!r} does not exist.'.
                                format(file_path, key))

        if self.known_sites_discovery:
            self.known_sites_discovery = Default.get_absolute_path(
                file_path=self.known_sites_discovery,
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(self.known_sites_discovery):
                raise Exception('The known_sites_discovery file {!r} does not exist.'.
                                format(self.known_sites_discovery))

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for file_path in self.known_sites_realignment:
            file_path = Default.get_absolute_path(file_path=file_path, default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The file path {!r} for known_sites_realignment does not exist.'.
                                format(file_path))
        self.known_sites_realignment = temporary_list

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for file_path in self.known_sites_recalibration:
            file_path = Default.get_absolute_path(file_path=file_path, default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The file path {!r} for known_sites_recalibration does not exist.'.
                                format(file_path))
        self.known_sites_recalibration = temporary_list

        for key in self.vqsr_resources_indel_dict:
            resource_dict = self.vqsr_resources_indel_dict[key]
            resource_dict['file_path'] = Default.get_absolute_path(
                file_path=resource_dict['file_path'],
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(resource_dict['file_path']):
                raise Exception('The file path {!r} for vqsr_resources_indel {!r} does not exist.'.
                                format(resource_dict['file_path'], key))

        for key in self.vqsr_resources_snp_dict:
            resource_dict = self.vqsr_resources_snp_dict[key]
            resource_dict['file_path'] = Default.get_absolute_path(
                file_path=resource_dict['file_path'],
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(resource_dict['file_path']):
                raise Exception('The file path {!r} for vqsr_resources_snp {!r} does not exist.'.
                                format(resource_dict['file_path'], key))

        # Exclude intervals

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for interval in self.exclude_intervals_list:
            if interval.endswith('.intervals') or interval.endswith('.interval_list'):
                # For Picard-style interval lists prepend the current directory if necessary.
                if not os.path.isabs(interval):
                    interval = os.path.join(os.path.realpath(os.path.curdir), interval)
                if os.path.exists(interval):
                    temporary_list.append(interval)
                else:
                    raise Exception('Exclude intervals file {!r} does not exist.'.format(interval))
            else:
                temporary_list.append(interval)
        self.exclude_intervals_list = temporary_list

        # Include intervals

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for interval in self.include_intervals_list:
            if interval.endswith('.intervals') or interval.endswith('.interval_list'):
                # For Picard-style interval lists prepend the current directory if necessary.
                if not os.path.isabs(interval):
                    interval = os.path.join(os.path.realpath(os.path.curdir), interval)
                if os.path.exists(interval):
                    temporary_list.append(interval)
                else:
                    raise Exception('Include intervals file {!r} does not exist.'.format(interval))
            else:
                temporary_list.append(interval)
        self.include_intervals_list = temporary_list

        # Genome Annotation GTF file path, defaults to the interval files directory.

        if self.genome_annotation_gtf and not os.path.isabs(self.genome_annotation_gtf):
            self.genome_annotation_gtf = Default.get_absolute_path(
                file_path=self.genome_annotation_gtf,
                default_path=Default.absolute_intervals())

        # Read comparisons for somatic mutation calling.
        self._read_comparisons(comparison_path=self.comparison_path)

        # Create genomic tiles for scatter gather approaches.
        if self.number_of_tiles_cohort:
            self._tile_region_cohort_list = run_create_genome_tiles(tiles=self.number_of_tiles_cohort, width=0)

        if self.number_of_tiles_somatic:
            self._tile_region_somatic_list = run_create_genome_tiles(tiles=self.number_of_tiles_somatic, width=0)

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        stage_align_lane = self.get_stage(name=self.stage_name_align_lane)
        stage_process_lane = self.get_stage(name=self.stage_name_process_lane)
        stage_process_sample = self.get_stage(name=self.stage_name_process_sample)
        stage_diagnose_sample = self.get_stage(name=self.stage_name_diagnose_sample)
        stage_merge_cohort = self.get_stage(name=self.stage_name_merge_cohort)
        stage_process_cohort = self.get_stage(name=self.stage_name_process_cohort)
        stage_annotate_cohort = self.get_stage(name=self.stage_name_annotate_cohort)
        stage_split_cohort = self.get_stage(name=self.stage_name_split_cohort)
        stage_summary = self.get_stage(name=self.stage_name_summary)
        stage_somatic = self.get_stage(name=self.stage_name_somatic)

        # Create a Python dict of Python str (cohort name) key and Python list of process_sample Runnable object
        # value data. This dictionary is required by the merge_cohort stage to hierarchically merge cohorts.

        vc_runnable_process_sample_dict = dict()
        """ @type vc_runnable_process_sample_dict: dict[str, list[Runnable]] """

        # Create a Python list of diagnose_sample Runnable objects.

        vc_runnable_diagnose_sample_list = list()
        """ @type vc_runnable_diagnose_sample_list: list[Runnable] """

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

            vc_runnable_process_read_group_list = list()
            """ @type vc_runnable_process_read_group_list: list[Runnable] """

            for paired_reads_name in paired_reads_name_list:
                if not len(paired_reads_dict[paired_reads_name]):
                    # Skip names, which PairedReads objects have all been excluded.
                    continue

                #################################
                # Step 1: Align per read group. #
                #################################
                #
                # bsf_run_bwa.py
                # - Picard SamToFastq
                # - BWA MEM

                bwa = BWA(name='variant_calling_bwa_' + paired_reads_name, analysis=self)
                # Instead of adding the BWA Executable to the Stage, it gets serialised into the pickler_file.
                # stage_bwa.add_executable(executable=bwa)

                bwa_mem = bwa.sub_command

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
                    if paired_reads.reads_1:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
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

                prefix_align_lane = '_'.join((stage_align_lane.name, paired_reads_name))
                # TODO: The name for the aligned BAM is constructed by the bsf_run_bwa.py script.
                # It is currently based on the stage_align_lane.name and paired_reads_name.
                # The script should also be changed to pre-set all file names beforehand.
                file_path_alignment = FilePathAlignment(prefix=prefix_align_lane)

                # Normally, the bwa object would be pushed onto the Stage list.
                # Experimentally, use Pickler to serialize the Executable into a file.

                pickler_dict_align_lane = {
                    'prefix': stage_align_lane.name,
                    'replicate_key': paired_reads_name,
                    'classpath_gatk': self.classpath_gatk,
                    'classpath_picard': self.classpath_picard,
                    'bwa_executable': bwa,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    stage_align_lane.name + '_' + paired_reads_name + '.pkl')
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_align_lane)
                pickler_file.close()

                # Create a bsf_run_bwa.py job to run the pickled object.

                run_bwa = stage_align_lane.add_executable(
                    executable=Executable(
                        name='_'.join((stage_align_lane.name, paired_reads_name)),
                        program='bsf_run_bwa.py'))

                # Only submit this Executable if the final result file does not exist.
                if (os.path.exists(os.path.join(self.genome_directory, file_path_alignment.aligned_md5)) and
                        os.path.getsize(os.path.join(self.genome_directory, file_path_alignment.aligned_md5))):
                    run_bwa.submit = False
                # Check also for existence of a new-style Runnable status file.
                if os.path.exists(os.path.join(
                        stage_align_lane.working_directory,
                        '_'.join((stage_align_lane.name, paired_reads_name, 'completed.txt')))):
                    run_bwa.submit = False

                    # Set run_bwa options.

                self.set_command_configuration(command=run_bwa)
                run_bwa.add_option_long(key='pickler_path', value=pickler_path)
                run_bwa.add_option_long(key='debug', value=str(self.debug))

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

                prefix_lane = '_'.join((stage_process_lane.name, paired_reads_name))

                # Lane-specific file paths

                file_path_process_read_group = FilePathProcessReadGroup(prefix=prefix_lane)

                # Create a Runnable and Executable for processing each read group.

                runnable_process_lane = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_lane,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        file_path_object=file_path_process_read_group,
                        debug=self.debug))
                executable_process_lane = self.set_stage_runnable(
                    stage=stage_process_lane,
                    runnable=runnable_process_lane)
                # Set dependencies on preceding Runnable.name or Executable.name objects.
                executable_process_lane.dependencies.append(run_bwa.name)
                # Set dependencies for succeeding Runnable or Executable objects.
                vc_runnable_process_read_group_list.append(runnable_process_lane)

                reference_process_lane = runnable_process_lane.get_absolute_cache_file_path(
                    file_path=self.bwa_genome_db)

                # Run the Picard MarkDuplicates analysis, unless configured to skip it.

                if self.skip_mark_duplicates:
                    file_path_process_read_group.duplicates_marked_bam = file_path_alignment.aligned_bam
                    file_path_process_read_group.duplicates_marked_bai = file_path_alignment.aligned_bai
                    file_path_process_read_group.duplicates_marked_md5 = file_path_alignment.aligned_md5
                else:
                    # Run the Picard MarkDuplicates analysis.

                    runnable_step = runnable_process_lane.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='process_lane_picard_mark_duplicates',
                            obsolete_file_path_list=[
                                file_path_alignment.aligned_bam,
                                file_path_alignment.aligned_bai,
                                file_path_alignment.aligned_md5,
                            ],
                            java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx4G',
                            picard_classpath=self.classpath_picard,
                            picard_command='MarkDuplicates'))
                    """ @type runnable_step: RunnableStepPicard """
                    runnable_step.add_picard_option(key='INPUT', value=file_path_alignment.aligned_bam)
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
                        key='TMP_DIR',
                        value=runnable_process_lane.get_relative_temporary_directory_path)
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

                    runnable_step = runnable_process_lane.add_runnable_step(
                        runnable_step=RunnableStepGATK(
                            name='process_lane_gatk_realigner_target_creator',
                            java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            gatk_classpath=self.classpath_gatk))
                    """ @type runnable_step: RunnableStepGATK """
                    runnable_step.add_gatk_option(key='analysis_type', value='RealignerTargetCreator')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                    if self.downsample_to_fraction:
                        runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                    # Run sequence processing steps always on the full genome.
                    # for interval in self.exclude_intervals_list:
                    #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                    # for interval in self.include_intervals_list:
                    #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                    # if self.interval_padding:
                    #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
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

                    runnable_step = runnable_process_lane.add_runnable_step(
                        runnable_step=RunnableStepGATK(
                            name='process_lane_gatk_indel_realigner',
                            obsolete_file_path_list=[
                                file_path_process_read_group.duplicates_marked_bam,
                                file_path_process_read_group.duplicates_marked_bai,
                                file_path_process_read_group.duplicates_marked_md5,
                            ],
                            java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            gatk_classpath=self.classpath_gatk))
                    """ @type runnable_step: RunnableStepGATK """
                    runnable_step.add_gatk_option(key='analysis_type', value='IndelRealigner')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                    if self.downsample_to_fraction:
                        runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                    # Run sequence processing steps always on the full genome.
                    # for interval in self.exclude_intervals_list:
                    #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                    # for interval in self.include_intervals_list:
                    #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                    # if self.interval_padding:
                    #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                    runnable_step.add_gatk_switch(key='keep_program_records')
                    runnable_step.add_gatk_switch(key='generate_md5')
                    runnable_step.add_gatk_option(key='bam_compression', value='9')
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

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_base_recalibrator_pre',
                        java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                runnable_step.add_gatk_option(key='analysis_type', value='BaseRecalibrator')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # Run sequence processing steps always on the full genome.
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_recalibration:
                    runnable_step.add_gatk_option(key='knownSites', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_read_group.realigned_bam)
                runnable_step.add_gatk_option(key='out', value=file_path_process_read_group.recalibration_table_pre)

                # Run the GATK BaseRecalibrator on-the-fly recalibration analysis to generate plots.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_base_recalibrator_post',
                        java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                runnable_step.add_gatk_option(key='analysis_type', value='BaseRecalibrator')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # Run sequence processing steps always on the full genome.
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_recalibration:
                    runnable_step.add_gatk_option(key='knownSites', value=file_path, override=True)
                runnable_step.add_gatk_option(key='BQSR', value=file_path_process_read_group.recalibration_table_pre)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_read_group.realigned_bam)
                runnable_step.add_gatk_option(key='out', value=file_path_process_read_group.recalibration_table_post)

                # Run the GATK AnalyzeCovariates analysis to create a recalibration plot.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_analyze_covariates',
                        java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                runnable_step.add_gatk_option(key='analysis_type', value='AnalyzeCovariates')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # Run sequence processing steps always on the full genome.
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
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

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_print_reads',
                        obsolete_file_path_list=[
                            file_path_process_read_group.realigned_bam,
                            file_path_process_read_group.realigned_bai,
                            file_path_process_read_group.realigned_md5,
                        ],
                        java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                runnable_step.add_gatk_option(key='analysis_type', value='PrintReads')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_lane)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # Run sequence processing steps always on the full genome.
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                runnable_step.add_gatk_switch(key='keep_program_records')
                runnable_step.add_gatk_switch(key='generate_md5')
                runnable_step.add_gatk_option(key='bam_compression', value='9')
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_read_group.realigned_bam)
                runnable_step.add_gatk_option(key='BQSR', value=file_path_process_read_group.recalibration_table_pre)
                runnable_step.add_gatk_option(key='out', value=file_path_process_read_group.recalibrated_bam)

                # Run the Picard CollectAlignmentSummaryMetrics analysis.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='process_lane_picard_collect_alignment_summary_metrics',
                        java_temporary_path=runnable_process_lane.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='CollectAlignmentSummaryMetrics'))
                """ @type runnable_step: RunnableStepPicard """
                runnable_step.add_picard_option(
                    key='INPUT',
                    value=file_path_process_read_group.recalibrated_bam)
                runnable_step.add_picard_option(
                    key='OUTPUT',
                    value=file_path_process_read_group.alignment_summary_metrics)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
                runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_process_lane.get_relative_temporary_directory_path)
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                # OPTIONS_FILE

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

            prefix_sample = '_'.join((stage_process_sample.name, sample.name))

            file_path_process_sample = FilePathProcessSample(prefix=prefix_sample)

            # Create a Runnable and Executable for processing each Sample.

            runnable_process_sample = self.add_runnable(
                runnable=Runnable(
                    name=prefix_sample,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_process_sample,
                    debug=self.debug))
            executable_process_sample = self.set_stage_runnable(
                stage=stage_process_sample,
                runnable=runnable_process_sample)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            for runnable_process_lane in vc_runnable_process_read_group_list:
                executable_process_sample.dependencies.append(runnable_process_lane.name)

            reference_process_sample = runnable_process_sample.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            if len(vc_runnable_process_read_group_list) == 1:
                # If there is only one read group, sample-level read processing can be skipped.
                # Rename files on the basis of the first and only list component.
                runnable_process_lane = vc_runnable_process_read_group_list[0]
                file_path_process_read_group = runnable_process_lane.file_path_object
                """ @type file_path_process_read_group: FilePathProcessReadGroup """

                # Rename the BAM file.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_sample_move_recalibrated_bam',
                        source_path=file_path_process_read_group.recalibrated_bam,
                        target_path=file_path_process_sample.realigned_bam))

                # Rename the BAI file.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_sample_move_recalibrated_bai',
                        source_path=file_path_process_read_group.recalibrated_bai,
                        target_path=file_path_process_sample.realigned_bai))

                # Rename the MD5 file.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_sample_move_recalibrated_md5',
                        source_path=file_path_process_read_group.recalibrated_md5,
                        target_path=file_path_process_sample.realigned_md5))

                # Link the Picard Alignment Summary Metrics.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='process_sample_link_alignment_metrics',
                        source_path=file_path_process_read_group.alignment_summary_metrics,
                        target_path=file_path_process_sample.alignment_summary_metrics))

                # Link the Picard Duplicate Metrics if it has been created.
                if not self.skip_mark_duplicates:
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepLink(
                            name='process_sample_link_duplicate_metrics',
                            source_path=file_path_process_read_group.duplicate_metrics,
                            target_path=file_path_process_sample.duplicate_metrics))
            else:
                # Run the Picard MergeSamFiles analysis.

                runnable_step = runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='process_sample_picard_merge_sam_files',
                        java_temporary_path=runnable_process_sample.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='MergeSamFiles'))
                """ @type runnable_step: RunnableStepPicard """
                runnable_step.add_picard_option(key='COMMENT', value='Merged from the following files:')
                for runnable_process_lane in vc_runnable_process_read_group_list:
                    file_path_process_read_group = runnable_process_lane.file_path_object
                    """ @type file_path_process_read_group: FilePathProcessReadGroup """
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
                    value=runnable_process_sample.get_relative_temporary_directory_path)
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
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='process_sample_move_merged_bam',
                            source_path=file_path_process_sample.merged_bam,
                            target_path=file_path_process_sample.duplicates_marked_bam))
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='process_sample_move_merged_bai',
                            source_path=file_path_process_sample.merged_bai,
                            target_path=file_path_process_sample.duplicates_marked_bai))
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='process_sample_move_merged_md5',
                            source_path=file_path_process_sample.merged_md5,
                            target_path=file_path_process_sample.duplicates_marked_md5))
                else:
                    # Run the Picard MarkDuplicates analysis.
                    # Optical duplicates should already have been flagged in the lane-specific processing step.

                    runnable_step = runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='process_sample_picard_mark_duplicates',
                            obsolete_file_path_list=[
                                file_path_process_sample.merged_bam,
                                file_path_process_sample.merged_bai,
                                file_path_process_sample.merged_md5,
                            ],
                            java_temporary_path=runnable_process_sample.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            picard_classpath=self.classpath_picard,
                            picard_command='MarkDuplicates'))
                    """ @type runnable_step: RunnableStepPicard """
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
                        key='TMP_DIR',
                        value=runnable_process_sample.get_relative_temporary_directory_path)
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
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='process_sample_move_duplicates_marked_bam',
                            source_path=file_path_process_sample.duplicates_marked_bam,
                            target_path=file_path_process_sample.realigned_bam))
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='process_sample_move_duplicates_marked_bai',
                            source_path=file_path_process_sample.duplicates_marked_bai,
                            target_path=file_path_process_sample.realigned_bai))
                    runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='process_sample_move_duplicates_marked_md5',
                            source_path=file_path_process_sample.duplicates_marked_md5,
                            target_path=file_path_process_sample.realigned_md5))
                else:
                    # Run the GATK RealignerTargetCreator analysis as the first-pass walker
                    # for the GATK IndelRealigner analysis.

                    runnable_step = runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepGATK(
                            name='process_sample_gatk_realigner_target_creator',
                            java_temporary_path=runnable_process_sample.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            gatk_classpath=self.classpath_gatk))
                    """ @type runnable_step: RunnableStepGATK """
                    runnable_step.add_gatk_option(key='analysis_type', value='RealignerTargetCreator')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_sample)
                    if self.downsample_to_fraction:
                        runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                    # Run sequence processing steps always on the full genome.
                    # for interval in self.exclude_intervals_list:
                    #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                    # for interval in self.include_intervals_list:
                    #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                    # if self.interval_padding:
                    #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                    for file_path in self.known_sites_realignment:
                        runnable_step.add_gatk_option(key='known', value=file_path, override=True)
                    runnable_step.add_gatk_option(
                        key='input_file',
                        value=file_path_process_sample.duplicates_marked_bam)
                    runnable_step.add_gatk_option(key='out', value=file_path_process_sample.realigner_targets)

                    # Run the GATK IndelRealigner analysis as a second-pass walker
                    # after the GATK RealignerTargetCreator analysis.

                    runnable_step = runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepGATK(
                            name='process_sample_gatk_indel_realigner',
                            obsolete_file_path_list=[
                                file_path_process_sample.duplicates_marked_bam,
                                file_path_process_sample.duplicates_marked_bai,
                                file_path_process_sample.duplicates_marked_md5,
                            ],
                            java_temporary_path=runnable_process_sample.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            gatk_classpath=self.classpath_gatk))
                    """ @type runnable_step: RunnableStepGATK """
                    runnable_step.add_gatk_option(key='analysis_type', value='IndelRealigner')
                    runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_sample)
                    if self.downsample_to_fraction:
                        runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                    # Run sequence processing steps always on the full genome.
                    # for interval in self.exclude_intervals_list:
                    #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                    # for interval in self.include_intervals_list:
                    #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                    # if self.interval_padding:
                    #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                    runnable_step.add_gatk_switch(key='keep_program_records')
                    runnable_step.add_gatk_switch(key='generate_md5')
                    runnable_step.add_gatk_option(key='bam_compression', value='9')
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

            runnable_process_sample.add_runnable_step(
                runnable_step=RunnableStepLink(
                    name='process_sample_link_bam_bai',
                    source_path=file_path_process_sample.realigned_bai,
                    target_path=file_path_process_sample.realigned_bam_bai))

            # Run the Picard CollectAlignmentSummaryMetrics analysis.

            runnable_step = runnable_process_sample.add_runnable_step(
                runnable_step=RunnableStepPicard(
                    name='process_sample_picard_collect_alignment_summary_metrics',
                    java_temporary_path=runnable_process_sample.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx6G',
                    picard_classpath=self.classpath_picard,
                    picard_command='CollectAlignmentSummaryMetrics'))
            """ @type runnable_step: RunnableStepPicard """
            runnable_step.add_picard_option(key='INPUT', value=file_path_process_sample.realigned_bam)
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_process_sample.alignment_summary_metrics)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
            runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_process_sample.get_relative_temporary_directory_path)
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

            runnable_step = runnable_process_sample.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_sample_gatk_haplotype_caller',
                    java_temporary_path=runnable_process_sample.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='HaplotypeCaller')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_sample)
            if self.downsample_to_fraction:
                runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
            # Run sequence processing steps always on the full genome.
            # for interval in self.exclude_intervals_list:
            #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            # for interval in self.include_intervals_list:
            #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            # if self.interval_padding:
            #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # The number of threads should be configurable, but multi-threading seems to cause the occasional problem.
            # runnable_step.add_gatk_option(key='num_cpu_threads_per_data_thread', value='1')
            # runnable_step.add_gatk_option(key='pair_hmm_implementation', value='VECTOR_LOGLESS_CACHING')
            runnable_step.add_gatk_option(key='genotyping_mode', value='DISCOVERY')
            # runnable_step.add_gatk_option(key='standard_min_confidence_threshold_for_emitting', value='10')
            # runnable_step.add_gatk_option(key='standard_min_confidence_threshold_for_calling', value='30')
            runnable_step.add_gatk_option(key='emitRefConfidence', value='GVCF')
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
            runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
            runnable_step.add_gatk_option(key='out', value=file_path_process_sample.raw_variants_gvcf_vcf)
            # Parameter to pass to the VCF/BCF IndexCreator
            # runnable_step.add_gatk_option(key='variant_index_type', value='LINEAR')
            # runnable_step.add_gatk_option(key='variant_index_parameter', value='128000')

            # Finally, record the process_sample Runnable for the merge_cohort stage under the Sample's
            # 'Cohort Name' annotation, or if it does not exist, under the cohort name defined in the
            # Analysis in the configuration file.

            if 'Cohort Name' in sample.annotation_dict:
                cohort_key = sample.annotation_dict['Cohort Name'][0]
            else:
                cohort_key = self.cohort_name

            if cohort_key not in vc_runnable_process_sample_dict:
                vc_runnable_process_sample_dict[cohort_key] = list()
            _runnable_process_sample_list = vc_runnable_process_sample_dict[cohort_key]
            _runnable_process_sample_list.append(runnable_process_sample)

            ################################
            # Step 4: Diagnose the sample. #
            ################################
            #
            # GATK DiagnoseTarget                   (diagnose_sample_gatk_diagnose_target)
            # GATK QualifyMissingIntervals          (diagnose_sample_gatk_qualify_missing_intervals)
            # GATK CallableLoci                     (diagnose_sample_gatk_callable_loci)
            # bsfR bsf_variant_calling_coverage.R   (diagnose_sample_coverage)
            # Picard CalculateHsMetrics             (diagnose_sample_picard_calculate_hybrid_selection_metrics)

            prefix_diagnose_sample = '_'.join((stage_diagnose_sample.name, sample.name))

            file_path_diagnose_sample = FilePathDiagnoseSample(prefix=prefix_diagnose_sample)

            target_interval_name = str()
            target_interval_path = str()
            probe_interval_path = str()

            if 'Target Name' in sample.annotation_dict:
                target_name_list = sample.annotation_dict['Target Name']
                if len(target_name_list) > 1:
                    warnings.warn('More than one set of Target Name annotations is currently not supported.\n'
                                  'Choosing the first one of {!r} for sample {!r}'.
                                  format(target_name_list, sample.name))
                target_interval_name = target_name_list[0]

            if 'Target Intervals' in sample.annotation_dict:
                target_interval_list = sample.annotation_dict['Target Intervals']
                if len(target_interval_list) > 1:
                    warnings.warn('More than one set of Target Interval annotations is currently not supported.\n'
                                  'Choosing the first one of {!r} for sample {!r}'.
                                  format(target_interval_list, sample.name))
                target_interval_path = target_interval_list[0]
                if target_interval_path and not os.path.isabs(target_interval_path):
                    target_interval_path = Default.get_absolute_path(
                        file_path=target_interval_path,
                        default_path=Default.absolute_intervals())

            if 'Probe Intervals' in sample.annotation_dict:
                probe_interval_list = sample.annotation_dict['Probe Intervals']
                if len(probe_interval_list) > 1:
                    warnings.warn('More than one set of Probe Interval annotations is currently not supported.\n'
                                  'Choosing the first one of ' +
                                  repr(probe_interval_list) +
                                  ' for sample ' +
                                  repr(sample.name))
                probe_interval_path = probe_interval_list[0]
                if probe_interval_path and not os.path.isabs(probe_interval_path):
                    probe_interval_path = Default.get_absolute_path(
                        file_path=probe_interval_path,
                        default_path=Default.absolute_intervals())

            # Create a Runnable and Executable for diagnosing each sample.

            runnable_diagnose_sample = self.add_runnable(
                runnable=Runnable(
                    name=prefix_diagnose_sample,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_diagnose_sample,
                    debug=self.debug))
            executable_diagnose_sample = self.set_stage_runnable(
                stage=stage_diagnose_sample,
                runnable=runnable_diagnose_sample)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_diagnose_sample.dependencies.append(runnable_process_sample.name)
            # Set dependencies for succeding Runnable or Executable objects.
            vc_runnable_diagnose_sample_list.append(runnable_diagnose_sample)

            reference_diagnose_sample = runnable_diagnose_sample.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            if target_interval_path:
                # Run the GATK DiagnoseTarget analysis per sample, only if targets have been defined.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='diagnose_sample_gatk_diagnose_target',
                        java_temporary_path=runnable_diagnose_sample.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                self.set_runnable_step_configuration(runnable_step=runnable_step)
                # CommandLineGATK
                # Required Parameters
                runnable_step.add_gatk_option(key='analysis_type', value='DiagnoseTargets')
                # Optional Inputs
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_diagnose_sample)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     sub_command.add_option_long(key='interval_padding', value=str(self.interval_padding))
                # The Diagnose Targets analysis is run on the target intervals, only.
                if target_interval_path:
                    runnable_step.add_gatk_option(key='intervals', value=target_interval_path)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
                # DiagnoseTargets
                # Optional Outputs
                runnable_step.add_gatk_option(
                    key='missing_intervals',
                    value=file_path_diagnose_sample.missing_intervals)
                runnable_step.add_gatk_option(key='out', value=file_path_diagnose_sample.diagnose_targets_vcf)
                # Optional Parameters

                # Run the GATK QualifyMissingIntervals analysis per sample, only if targets have been defined.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='diagnose_sample_gatk_qualify_missing_intervals',
                        java_temporary_path=runnable_diagnose_sample.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                self.set_runnable_step_configuration(runnable_step=runnable_step)
                # CommandLineGATK
                # Required Parameters
                runnable_step.add_gatk_option(key='analysis_type', value='QualifyMissingIntervals')
                # Optional Inputs
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_diagnose_sample)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     sub_command.add_option_long(key='interval_padding', value=str(self.interval_padding))
                # The Qualify Missing Intervals analysis is run on the missing intervals
                # of the Diagnose Targets analysis, regardless.
                runnable_step.add_gatk_option(key='intervals', value=file_path_diagnose_sample.missing_intervals)
                runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
                # QualifyMissingIntervals
                # Required Parameters
                if target_interval_path:
                    runnable_step.add_gatk_option(key='targetsfile', value=target_interval_path)
                # Optional Outputs
                runnable_step.add_gatk_option(key='out', value=file_path_diagnose_sample.missing_report)
                # Optional Parameters
                if probe_interval_path:
                    runnable_step.add_gatk_option(key='baitsfile', value=probe_interval_path)

            # Run the GATK CallableLoci analysis per sample.

            runnable_step = runnable_diagnose_sample.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='diagnose_sample_gatk_callable_loci',
                    java_temporary_path=runnable_diagnose_sample.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx6G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            self.set_runnable_step_configuration(runnable_step=runnable_step)
            # CommandLineGATK
            # Required Parameters
            runnable_step.add_gatk_option(key='analysis_type', value='CallableLoci')
            # Optional Inputs
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_diagnose_sample)
            if self.downsample_to_fraction:
                runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            # for interval in self.include_intervals_list:
            #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            # if self.interval_padding:
            #     sub_command.add_option_long(key='interval_padding', value=str(self.interval_padding))
            # The Callable Loci analysis is run on the target intervals, only.
            if target_interval_path:
                runnable_step.add_gatk_option(key='intervals', value=target_interval_path)
            runnable_step.add_gatk_option(key='input_file', value=file_path_process_sample.realigned_bam)
            # Required Outputs
            runnable_step.add_gatk_option(key='summary', value=file_path_diagnose_sample.callable_txt)
            # Optional Outputs
            runnable_step.add_gatk_option(key='out', value=file_path_diagnose_sample.callable_bed)

            # Run the bsfR bsf_variant_calling_coverage.R script.

            runnable_step = runnable_diagnose_sample.add_runnable_step(
                runnable_step=RunnableStep(
                    name='diagnose_sample_coverage',
                    program='bsf_variant_calling_coverage.R'))
            """ @type runnable_step: RunnableStep """
            self.set_runnable_step_configuration(runnable_step=runnable_step)
            runnable_step.add_option_long(key='exons', value=self.genome_annotation_gtf)
            runnable_step.add_option_long(key='callable-loci', value=file_path_diagnose_sample.callable_bed)
            if target_interval_path.endswith('.bed'):
                runnable_step.add_option_long(key='targets', value=target_interval_path)
            elif target_interval_path.endswith('.interval_list'):
                runnable_step.add_option_long(key='targets', value=target_interval_path[:-13] + 'bed')
            elif target_interval_path:
                runnable_step.add_option_long(key='targets', value=target_interval_path)
            # If a target interval path has not been defined, run without it.

            # Run the UCSC bedSort tool.

            runnable_step = runnable_diagnose_sample.add_runnable_step(
                runnable_step=RunnableStep(
                    name='diagnose_sample_bed_sort',
                    program='bedSort'))
            """ @type runnable_step: RunnableStep """
            runnable_step.arguments.append(file_path_diagnose_sample.callable_bed)
            runnable_step.arguments.append(file_path_diagnose_sample.sorted_bed)

            # Run the UCSC bedToBigBed tool.

            runnable_step = runnable_diagnose_sample.add_runnable_step(
                runnable_step=RunnableStep(
                    name='diagnose_sample_bed_to_big_bed',
                    program='bedToBigBed',
                    obsolete_file_path_list=[
                        file_path_diagnose_sample.sorted_bed,
                    ]))
            """ @type runnable_step: RunnableStep """
            # FIXME: It would be good to allow options with and without an equal sign.
            runnable_step.add_switch_short(key='type=bed4')
            runnable_step.arguments.append(file_path_diagnose_sample.sorted_bed)
            runnable_step.arguments.append(reference_diagnose_sample + '.fai')
            runnable_step.arguments.append(file_path_diagnose_sample.callable_bb)

            if target_interval_path:
                # Run the Picard CalculateHsMetrics analysis per sample, only if targets have been defined.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='diagnose_sample_picard_calculate_hybrid_selection_metrics',
                        java_temporary_path=runnable_diagnose_sample.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='CalculateHsMetrics'))
                """ @type runnable_step: RunnableStepPicard """
                if probe_interval_path:
                    runnable_step.add_picard_option(key='BAIT_INTERVALS', value=probe_interval_path)
                else:
                    runnable_step.add_picard_option(key='BAIT_INTERVALS', value=target_interval_path)
                if target_interval_name:
                    runnable_step.add_picard_option(key='BAIT_SET_NAME', value=target_interval_name)
                runnable_step.add_picard_option(key='TARGET_INTERVALS', value=target_interval_path)
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
                    value=runnable_diagnose_sample.get_relative_temporary_directory_path)
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

        # Create a Python dict of Python str (cohort name) key and Python list of Runnable value data.
        # Initialise a single key with the final cohort name and an empty list for merging the final cohort.

        vc_runnable_merge_cohort_dict = {self.cohort_name: []}
        vc_runnable_merge_cohort_list = vc_runnable_merge_cohort_dict[self.cohort_name]

        # Run the GATK CombineGVCFs analysis for each cohort and Sample defined in this project to build up
        # cohort-specific GVCF files.

        for cohort_key in vc_runnable_process_sample_dict.keys():
            vc_runnable_merge_cohort_list.append(
                run_merge_cohort_scatter_gather(
                    analysis_stage=stage_merge_cohort,
                    cohort_runnable_dict=vc_runnable_process_sample_dict,
                    cohort_name=cohort_key))

        # Run the GATK CombineGVCF analysis once more to merge all cohort-specific GVCF files defined in this project.

        if len(vc_runnable_merge_cohort_list) == 1:
            # If the cohort-specific Runnable list has only one component, the merge has already been completed.
            runnable_merge_cohort = vc_runnable_merge_cohort_list[-1]
        elif len(vc_runnable_merge_cohort_list) > 1:
            runnable_merge_cohort = run_merge_cohort_scatter_gather(
                analysis_stage=stage_merge_cohort,
                cohort_runnable_dict=vc_runnable_merge_cohort_dict,
                cohort_name=self.cohort_name)
        else:
            raise Exception("Unexpected number of Runnable objects on the merge_cohort list.")

        # Run an additional GATK CombineGVCFs analysis to merge into a super-cohort, if defined.

        if len(self.accessory_cohort_gvcfs):
            # If accessory cohorts are defined, initialise a new Python dict of Python str cohort name key and
            # Python list of Runnable value data. Initialise the list with teh last Runnable object and
            # extend with the the list of accessory cohort file names. The run_merge_cohort_scatter_gather() method
            # can cope with Runnable or str | unicode objects.
            cohort_key = '_'.join((self.cohort_name, 'accessory'))
            vc_runnable_merge_cohort_dict = {cohort_key: [runnable_merge_cohort]}
            vc_runnable_merge_cohort_list = vc_runnable_merge_cohort_dict[cohort_key]
            vc_runnable_merge_cohort_list.extend(self.accessory_cohort_gvcfs)

            runnable_merge_cohort = run_merge_cohort_scatter_gather(
                analysis_stage=stage_merge_cohort,
                cohort_runnable_dict=vc_runnable_merge_cohort_dict,
                cohort_name=cohort_key)

        # Specify the final FilePathMergeCohort object from the final Runnable object.
        file_path_merge_cohort = runnable_merge_cohort.file_path_object
        """ @type file_path_merge_cohort: FilePathMergeCohort """

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

        prefix_process_cohort = '_'.join((stage_process_cohort.name, self.cohort_name))

        file_path_genotype_cohort = FilePathGenotypeCohort(prefix=prefix_process_cohort)

        # Run GATK GenotypeGVCFs in a scatter and gather approach.

        if os.path.exists(os.path.join(self.genome_directory, file_path_genotype_cohort.genotyped_raw_tbi)):
            submit_runnable = False
        else:
            submit_runnable = True

        vc_runnable_process_cohort_scatter_list = list()
        """ @type vc_runnable_process_cohort_scatter_list: list[Runnable] """
        runnable_process_cohort_scatter = None
        """ @type runnable_process_cohort_scatter: Runnable """
        for tile_index in range(0, len(self._tile_region_cohort_list)):
            prefix_process_cohort_scatter = '_'.join((
                stage_process_cohort.name, self.cohort_name, 'scatter', str(tile_index)))

            file_path_genotype_cohort_scatter = FilePathGenotypeCohort(prefix=prefix_process_cohort_scatter)

            runnable_process_cohort_scatter = self.add_runnable(
                runnable=Runnable(
                    name=prefix_process_cohort_scatter,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_genotype_cohort_scatter,
                    debug=self.debug))
            executable_process_cohort_scatter = self.set_stage_runnable(
                stage=stage_process_cohort,
                runnable=runnable_process_cohort_scatter)
            executable_process_cohort_scatter.submit = submit_runnable
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_process_cohort_scatter.dependencies.append(runnable_merge_cohort.name)

            vc_runnable_process_cohort_scatter_list.append(runnable_process_cohort_scatter)

            reference_process_cohort_scatter = runnable_process_cohort_scatter.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the GATK GenotypeGVCFs analysis.

            runnable_step = runnable_process_cohort_scatter.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_genotype_gvcfs_scatter',
                    java_temporary_path=runnable_process_cohort_scatter.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx12G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='GenotypeGVCFs')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort_scatter)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            # for interval in self.include_intervals_list:
            #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            # if self.interval_padding:
            #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            for region_tuple in self._tile_region_cohort_list[tile_index]:
                # The list of tiles is initialised to an empty tile to trigger at least one process.
                # Do not assign an interval in such cases.
                if region_tuple[0]:
                    runnable_step.add_gatk_option(
                        key='intervals',
                        value='{:s}:{:d}-{:d}'.format(region_tuple[0], region_tuple[1], region_tuple[2]),
                        override=True)
            # Scatter gather is more robust than GATK multi-threading.
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
            runnable_step.add_gatk_option(key='variant', value=file_path_merge_cohort.combined_gvcf_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_genotype_cohort_scatter.partial_vcf)

        # Gather

        # If there is only one tile, no need to gather, just rename the file and return the Runnable.

        if len(self._tile_region_cohort_list) == 1:
            file_path_genotype_cohort_scatter = runnable_process_cohort_scatter.file_path_object
            """ @type file_path_genotype_cohort_scatter: FilePathGenotypeCohort """
            # Add cohort-specific keys to the file path dictionary.
            file_path_genotype_cohort_scatter.genotyped_raw_vcf = file_path_genotype_cohort.genotyped_raw_vcf
            file_path_genotype_cohort_scatter.genotyped_raw_tbi = file_path_genotype_cohort.genotyped_raw_tbi

            runnable_process_cohort_scatter.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='process_cohort_gather_move_vcf',
                    source_path=file_path_genotype_cohort_scatter.partial_vcf,
                    target_path=file_path_genotype_cohort_scatter.genotyped_raw_vcf))

            runnable_process_cohort_scatter.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='process_cohort_gather_move_tbi',
                    source_path=file_path_genotype_cohort_scatter.partial_tbi,
                    target_path=file_path_genotype_cohort_scatter.genotyped_raw_vcf))

            # return runnable_merge_cohort_scatter
            runnable_process_cohort_gather = runnable_process_cohort_scatter
        else:
            # Second, gather by the number of chunks on the partitioned genome tile index list.

            # Second, Merge chunks hierarchically.
            # Initialise a list of Runnable objects and indices for the hierarchical merge.
            vc_runnable_process_cohort_gather_list = vc_runnable_process_cohort_scatter_list
            vc_process_cohort_gather_index_list = range(0, len(self._tile_region_cohort_list))
            runnable_process_cohort_gather = None  # Global variable to keep and return the last Runnable.
            """ @type runnable_process_cohort_gather: Runnable """
            level = 0
            while len(vc_process_cohort_gather_index_list) > 1:
                temporary_gather_runnable_list = list()
                """ @type temporary_gather_runnable_list: list[Runnable] """
                temporary_gather_index_list = list()
                """ @type temporary_gather_index_list: list[int] """
                # Partition the index list into chunks of given size.
                partition_list = [vc_process_cohort_gather_index_list[offset:offset + self.number_of_chunks_cohort]
                                  for offset in range(0,
                                                      len(vc_process_cohort_gather_index_list),
                                                      self.number_of_chunks_cohort)]

                for partition_index in range(0, len(partition_list)):
                    chunk_index_list = partition_list[partition_index]
                    # The file prefix includes the level and partition index.
                    prefix_process_cohort_gather = '_'.join(
                        (stage_process_cohort.name, self.cohort_name, 'gather', str(level), str(partition_index)))

                    file_path_genotype_cohort_gather = FilePathGenotypeCohort(prefix=prefix_process_cohort_gather)

                    runnable_process_cohort_gather = self.add_runnable(
                        runnable=Runnable(
                            name=prefix_process_cohort_gather,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            cache_directory=self.cache_directory,
                            cache_path_dict=self._cache_path_dict,
                            file_path_object=file_path_genotype_cohort_gather,
                            debug=self.debug))
                    executable_process_cohort_gather = self.set_stage_runnable(
                        stage=stage_process_cohort,
                        runnable=runnable_process_cohort_gather)
                    executable_process_cohort_gather.submit = submit_runnable

                    # Dependencies on scatter processes are set based on genome tile indices below.
                    temporary_gather_runnable_list.append(runnable_process_cohort_gather)
                    temporary_gather_index_list.append(partition_index)

                    reference_process_cohort_gather = runnable_process_cohort_gather.get_absolute_cache_file_path(
                        file_path=self.bwa_genome_db)

                    # GATK CatVariants by-passes the GATK engine and thus requires a completely different command line.
                    runnable_step = runnable_process_cohort_gather.add_runnable_step(
                        runnable_step=RunnableStepJava(
                            name='merge_cohort_gatk_cat_variants',
                            sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                            java_temporary_path=runnable_process_cohort_gather.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx4G'))
                    """ @type runnable_step: RunnableStepJava """
                    runnable_step.add_option_short(
                        key='classpath',
                        value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                    sub_command = runnable_step.sub_command
                    # Add the 'reference' not 'reference_sequence' option.
                    sub_command.add_option_long(
                        key='reference',
                        value=reference_process_cohort_gather)
                    sub_command.add_option_long(
                        key='outputFile',
                        value=file_path_genotype_cohort_gather.partial_vcf)
                    sub_command.add_switch_long(key='assumeSorted')
                    # Finally, process per chunk index.
                    for chunk_index in chunk_index_list:
                        runnable_object = vc_runnable_process_cohort_gather_list[chunk_index]
                        file_path_object = runnable_object.file_path_object
                        """ @type file_path_object: FilePathGenotypeCohort """
                        # Set GATK option variant
                        sub_command.add_option_long(key='variant', value=file_path_object.partial_vcf, override=True)
                        # Delete the *.g.vcf.gz file.
                        runnable_step.obsolete_file_path_list.append(file_path_object.partial_vcf)
                        # Delete the *.g.vcf.gz.tbi file.
                        runnable_step.obsolete_file_path_list.append(file_path_object.partial_tbi)
                        # Set dependencies on preceding Runnable.name or Executable.name objects.
                        executable_process_cohort_gather.dependencies.append(runnable_object.name)

                # Set the temporary index list as the new list and increment the merge level.
                vc_runnable_process_cohort_gather_list = temporary_gather_runnable_list
                vc_process_cohort_gather_index_list = temporary_gather_index_list
                level += 1
            else:
                # For the last instance, additionally rename the final file.
                file_path_genotype_cohort_gather = runnable_process_cohort_gather.file_path_object
                """ @type file_path_genotype_cohort_gather: FilePathGenotypeCohort """

                # Add cohort-specific keys to the file path dictionary.
                file_path_genotype_cohort_gather.genotyped_raw_vcf = file_path_genotype_cohort.genotyped_raw_vcf
                file_path_genotype_cohort_gather.genotyped_raw_tbi = file_path_genotype_cohort.genotyped_raw_tbi

                runnable_process_cohort_gather.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_cohort_gather_move_vcf',
                        source_path=file_path_genotype_cohort_gather.partial_vcf,
                        target_path=file_path_genotype_cohort_gather.genotyped_raw_vcf))

                runnable_process_cohort_gather.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_cohort_gather_move_tbi',
                        source_path=file_path_genotype_cohort_gather.partial_tbi,
                        target_path=file_path_genotype_cohort_gather.genotyped_raw_tbi))

                # return runnable_process_cohort_gather

        file_path_process_cohort = FilePathProcessCohort(prefix=prefix_process_cohort)

        # Create a Runnable and Executable for processing the cohort.

        runnable_process_cohort = self.add_runnable(
            runnable=Runnable(
                name=prefix_process_cohort,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                cache_path_dict=self._cache_path_dict,
                file_path_object=file_path_process_cohort,
                debug=self.debug))
        executable_process_cohort = self.set_stage_runnable(
            stage=stage_process_cohort,
            runnable=runnable_process_cohort)
        # Set dependencies on preceding Runnable.name or Executable.name objects.
        executable_process_cohort.dependencies.append(runnable_process_cohort_gather.name)

        reference_process_cohort = runnable_process_cohort.get_absolute_cache_file_path(
            file_path=self.bwa_genome_db)

        # Run the VQSR procedure on SNPs.

        if self.vqsr_skip_snp:
            file_path_process_cohort.recalibrated_snp_raw_indel_vcf = file_path_genotype_cohort.genotyped_raw_vcf
            file_path_process_cohort.recalibrated_snp_raw_indel_idx = file_path_genotype_cohort.genotyped_raw_tbi
        else:

            # Run the GATK VariantRecalibrator analysis on SNPs.

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_variant_recalibrator_snp',
                    java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
            runnable_step.add_gatk_option(key='mode', value='SNP')
            for resource in self.vqsr_resources_snp_dict.keys():
                resource_option = 'resource:{},known={},training={},truth={},prior={}'.format(
                    resource,
                    self.vqsr_resources_snp_dict[resource]['known'],
                    self.vqsr_resources_snp_dict[resource]['training'],
                    self.vqsr_resources_snp_dict[resource]['truth'],
                    self.vqsr_resources_snp_dict[resource]['prior'])
                runnable_step.add_gatk_option(
                    key=resource_option,
                    value=self.vqsr_resources_snp_dict[resource]['file_path'])
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

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_apply_recalibration_snp',
                    java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='ApplyRecalibration')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
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
            file_path_process_cohort.recalibrated_snp_recalibrated_indel_idx = \
                file_path_process_cohort.recalibrated_snp_raw_indel_idx
        else:

            # Run the GATK VariantRecalibrator analysis on INDELs.

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_variant_recalibrator_indel',
                    java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
            runnable_step.add_gatk_option(key='mode', value='INDEL')
            for resource in self.vqsr_resources_indel_dict.keys():
                resource_option = 'resource:{},known={},training={},truth={},prior={}'.format(
                    resource,
                    self.vqsr_resources_indel_dict[resource]['known'],
                    self.vqsr_resources_indel_dict[resource]['training'],
                    self.vqsr_resources_indel_dict[resource]['truth'],
                    self.vqsr_resources_indel_dict[resource]['prior'])
                runnable_step.add_gatk_option(
                    key=resource_option,
                    value=self.vqsr_resources_indel_dict[resource]['file_path'])
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

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_apply_recalibration_indel',
                    java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='ApplyRecalibration')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
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

        if len(self.accessory_cohort_gvcfs):
            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_select_variants_cohort',
                    java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
            runnable_step.add_gatk_option(
                key='variant',
                value=file_path_process_cohort.recalibrated_snp_recalibrated_indel_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_process_cohort.multi_sample_vcf)
            for sample in self.sample_list:
                runnable_step.add_gatk_option(key='sample_name', value=sample.name, override=True)
            runnable_step.add_gatk_switch(key='excludeNonVariants')
        else:
            file_path_process_cohort.multi_sample_vcf = \
                file_path_process_cohort.recalibrated_snp_recalibrated_indel_vcf
            file_path_process_cohort.multi_sample_idx = \
                file_path_process_cohort.recalibrated_snp_recalibrated_indel_idx

        # Run the snpEff tool for functional variant annotation.

        runnable_step = runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='process_cohort_snpeff',
                program='java',
                sub_command=Command(program='eff')))
        """ @type runnable_step: RunnableStep """

        runnable_step.add_switch_short(
            key='d64')
        runnable_step.add_option_short(
            key='jar',
            value=os.path.join(self.classpath_snpeff, 'snpEff.jar'))
        runnable_step.add_switch_short(
            key='Xmx6G')
        runnable_step.add_option_pair(
            key='-Djava.io.tmpdir',
            value=runnable_process_cohort.get_relative_temporary_directory_path)
        runnable_step.stdout_path = file_path_process_cohort.snpeff_vcf

        sub_command = runnable_step.sub_command
        sub_command.add_switch_short(key='download')
        sub_command.add_option_short(key='o', value='gatk')
        sub_command.add_option_short(key='stats', value=file_path_process_cohort.snpeff_stats)
        sub_command.add_option_short(key='config', value=os.path.join(self.classpath_snpeff, 'snpEff.config'))

        sub_command.arguments.append(self.snpeff_genome_version)
        sub_command.arguments.append(file_path_process_cohort.multi_sample_vcf)

        # Automatically compress and index the snpEff VCF file with bgzip and tabix, respectively.
        # TODO: It would be better for the file system, if output could be directly piped into bgzip.

        runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='snpeff_bgzip',
                program='bgzip',
                arguments=[file_path_process_cohort.snpeff_vcf]))

        runnable_step = runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='snpeff_tabix',
                program='tabix',
                arguments=[file_path_process_cohort.snpeff_vcf_bgz]))
        """ @type runnable_step: RunnableStep """
        runnable_step.add_option_long(key='preset', value='vcf')

        # Run the GATK VariantAnnotator analysis.

        runnable_step = runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStepGATK(
                name='process_cohort_gatk_variant_annotator',
                java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                java_heap_maximum='Xmx4G',
                gatk_classpath=self.classpath_gatk))
        """ @type runnable_step: RunnableStepGATK """
        runnable_step.add_gatk_option(key='analysis_type', value='VariantAnnotator')
        runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
        for interval in self.exclude_intervals_list:
            runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
        for interval in self.include_intervals_list:
            runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
        if self.interval_padding:
            runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
        # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
        if self.known_sites_discovery:
            runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)

        # Add annotation resources and their corresponding expression options.
        for annotation_resource in self.annotation_resources_dict.keys():
            if len(self.annotation_resources_dict[annotation_resource][0]) \
                    and len(self.annotation_resources_dict[annotation_resource][1]):
                runnable_step.add_gatk_option(
                    key=':'.join(('resource', annotation_resource)),
                    value=self.annotation_resources_dict[annotation_resource][0])
                for annotation in self.annotation_resources_dict[annotation_resource][1]:
                    runnable_step.add_gatk_option(
                        key='expression',
                        value='.'.join((annotation_resource, annotation)),
                        override=True)

        runnable_step.add_gatk_option(key='variant', value=file_path_process_cohort.multi_sample_vcf)
        # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
        # runnable_step.add_gatk_option(key='annotation', value='AlleleBalanceBySample')
        runnable_step.add_gatk_option(key='annotation', value='SnpEff')
        runnable_step.add_gatk_option(key='snpEffFile', value=file_path_process_cohort.snpeff_vcf_bgz)
        runnable_step.add_gatk_option(key='out', value=file_path_process_cohort.annotated_vcf)

        ################################
        # Step 7: Annotate the cohort. #
        ################################
        #
        # Ensembl Variant Effect Predictor (annotate_cohort_ensembl_vep)
        # Bgzip                            (annotate_cohort_ensembl_vep_bgzip)
        # Tabix                            (annotate_cohort_ensembl_vep_tabix)
        # Ensembl Variant Effect Filter    (annotate_cohort_ensembl_filter)
        # Bgzip                            (annotate_cohort_ensembl_filter_bgzip)
        # Tabix                            (annotate_cohort_ensembl_filter_tabix)
        # VCF.Filter                       (annotate_cohort_vcf_filter_complete)
        # VCF.Filter                       (annotate_cohort_vcf_filter_filtered)

        prefix_annotate_cohort = '_'.join((stage_annotate_cohort.name, self.cohort_name))

        file_path_annotate_cohort = FilePathAnnotateCohort(prefix=prefix_annotate_cohort)

        runnable_annotate_cohort = self.add_runnable(
            runnable=Runnable(
                name=prefix_annotate_cohort,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                cache_path_dict=self._cache_path_dict,
                file_path_object=file_path_annotate_cohort,
                debug=self.debug))
        executable_annotate_cohort = self.set_stage_runnable(
            stage=stage_annotate_cohort,
            runnable=runnable_annotate_cohort)
        # Set dependencies on preceding Runnable.name or Executable.name objects.
        executable_annotate_cohort.dependencies.append(runnable_process_cohort.name)

        # reference_annotate_cohort = runnable_annotate_cohort.get_absolute_cache_file_path(
        #     file_path=self.bwa_genome_db)

        # if not os.path.exists(os.path.join(self.genome_directory, file_path_annotate_cohort.vep_complete_vcf_tbi)):
        # Run the Ensembl Variant Effect Predictor script.

        runnable_step = runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='annotate_cohort_ensembl_vep',
                program='perl',
                sub_command=Command()))
        """ @type runnable_step: RunnableStep """
        # self.set_runnable_step_configuration(runnable_step=runnable_step)
        runnable_step.arguments.append(os.path.join(default.directory_vep_src, 'vep'))
        sub_command = runnable_step.sub_command
        sub_command.add_switch_long(key='no_progress')
        sub_command.add_switch_long(key='everything')
        sub_command.add_option_long(key='species', value='homo_sapiens')  # TODO: Has to be configurable
        sub_command.add_option_long(key='assembly', value='GRCh37')  # TODO: Has to be configurable
        sub_command.add_option_long(key='input_file', value=file_path_process_cohort.multi_sample_vcf)
        sub_command.add_option_long(key='format', value='vcf')  # Input file format
        sub_command.add_option_long(key='output_file', value=file_path_annotate_cohort.vep_complete_raw_vcf)
        sub_command.add_option_long(key='stats_file', value=file_path_annotate_cohort.vep_statistics)
        sub_command.add_switch_long(key='force_overwrite')
        sub_command.add_switch_long(key='dont_skip')
        sub_command.add_switch_long(key='cache')
        sub_command.add_option_long(key='dir_cache', value=default.directory_vep_cache)
        sub_command.add_option_long(key='dir_plugins', value=default.directory_vep_plugins)
        sub_command.add_option_long(key='failed', value='1')
        sub_command.add_switch_long(key='vcf')
        sub_command.add_switch_long(key='allele_number')
        # sub_command.add_switch_long(key='no_escape')  # Do not percent escape HGVS strings
        sub_command.add_switch_long(key='allow_non_variant')
        sub_command.add_switch_long(key='flag_pick_allele_gene')
        sub_command.add_option_long(key='port', value='3337')
        sub_command.add_switch_long(key='gencode_basic')
        # sub_command.add_switch_long(key='write_cache')
        sub_command.add_option_long(  # TODO: Has to be configurable
            key='plugin',
            value='CADD,/scratch/lab_bsf/resources/CADD/b37/whole_genome_SNVs.tsv.gz')
        sub_command.add_option_long(
            key='tmpdir',
            value=runnable_annotate_cohort.get_relative_temporary_directory_path)

        runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='annotate_cohort_ensembl_vep_bgzip',
                program='bgzip',
                arguments=[file_path_annotate_cohort.vep_complete_raw_vcf]))

        runnable_step = runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='annotate_cohort_ensembl_vep_tabix',
                program='tabix',
                arguments=[file_path_annotate_cohort.vep_complete_raw_vcf_bgz]))
        """ @type runnable_step: RunnableStep """
        runnable_step.add_option_long(key='preset', value='vcf')

        # if not os.path.exists(os.path.join(self.genome_directory, file_path_annotate_cohort.vep_filtered_vcf_tbi)):
        # Run the Ensembl Variant Effect Filter script.

        runnable_step = runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='annotate_cohort_ensembl_filter',
                program='perl',
                sub_command=Command()))
        """ @type runnable_step: RunnableStep """
        # self.set_runnable_step_configuration(runnable_step=runnable_step)
        runnable_step.arguments.append(os.path.join(default.directory_vep_src, 'filter_vep'))
        sub_command = runnable_step.sub_command
        sub_command.add_option_long(key='input_file', value=file_path_annotate_cohort.vep_complete_raw_vcf_bgz)
        sub_command.add_option_long(key='format', value='vcf')
        sub_command.add_option_long(key='output_file', value=file_path_annotate_cohort.vep_filtered_raw_vcf)
        sub_command.add_switch_long(key='only_matched')
        sub_command.add_option_long(key='filter', value='Consequence ne upstream_gene_variant', override=True)
        sub_command.add_option_long(key='filter', value='Consequence ne downstream_gene_variant', override=True)
        sub_command.add_option_long(key='filter', value='Consequence ne intron_variant', override=True)
        sub_command.add_option_long(key='filter', value='BIOTYPE ne processed_transcript', override=True)
        # sub_command.add_option_long(key='filter', value='CANONICAL eq YES', override=True)
        sub_command.add_switch_long(key='force_overwrite')

        runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='annotate_cohort_ensembl_filter_bgzip',
                program='bgzip',
                arguments=[file_path_annotate_cohort.vep_filtered_raw_vcf]))

        runnable_step = runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='annotate_cohort_ensembl_filter_tabix',
                program='tabix',
                arguments=[file_path_annotate_cohort.vep_filtered_raw_vcf_bgz]))
        """ @type runnable_step: RunnableStep """
        runnable_step.add_option_long(key='preset', value='vcf')

        # Run the VCF Filter on the complete VEP set to convert (re-model) the CSQ field into
        # a set of independent INFO fields.

        runnable_step = runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStepJava(
                name='annotate_cohort_vcf_filter_complete',
                sub_command=Command(program='at.ac.oeaw.cemm.bsf.vcffilter.vep.vep2vcf'),
                java_temporary_path=runnable_annotate_cohort.get_relative_temporary_directory_path,
                java_heap_maximum='Xmx2G'))
        """ @type runnable_step: RunnableStepJava """
        runnable_step.add_option_short(
            key='classpath',
            value=os.path.join(self.classpath_vcf_filter, 'VCFFilter.jar'))
        sub_command = runnable_step.sub_command
        sub_command.add_option_pair(key='INPUT', value=file_path_annotate_cohort.vep_complete_raw_vcf_bgz)
        sub_command.add_option_pair(key='OUTPUT', value=file_path_annotate_cohort.vep_complete_vcf_bgz)

        # Run the VCF Filter on the filtered VEP set to convert (re-model) the CSQ field into
        # a set of independent INFO fields.

        runnable_step = runnable_annotate_cohort.add_runnable_step(
            runnable_step=RunnableStepJava(
                name='annotate_cohort_vcf_filter_filtered',
                sub_command=Command(program='at.ac.oeaw.cemm.bsf.vcffilter.vep.vep2vcf'),
                java_temporary_path=runnable_annotate_cohort.get_relative_temporary_directory_path,
                java_heap_maximum='Xmx2G'))
        """ @type runnable_step: RunnableStepJava """
        runnable_step.add_option_short(
            key='classpath',
            value=os.path.join(self.classpath_vcf_filter, 'VCFFilter.jar'))
        sub_command = runnable_step.sub_command
        sub_command.add_option_pair(key='INPUT', value=file_path_annotate_cohort.vep_filtered_raw_vcf_bgz)
        sub_command.add_option_pair(key='OUTPUT', value=file_path_annotate_cohort.vep_filtered_vcf_bgz)

        ######################################################
        # Step 8: Re-process and split the cohort by sample. #
        ######################################################
        #
        # GATK SelectVariants   (split_cohort_gatk_select_variants)
        # GATK VariantsToTable  (split_cohort_gatk_variants_to_table)

        for sample in self.sample_list:
            # Get all PairedReads objects solely to exclude samples without any.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)
            if not len(paired_reads_dict):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            prefix_split_cohort = '_'.join((stage_split_cohort.name, sample.name))

            file_path_split_cohort = FilePathSplitCohort(prefix=prefix_split_cohort)

            runnable_split_cohort = self.add_runnable(
                runnable=Runnable(
                    name=prefix_split_cohort,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_split_cohort,
                    debug=self.debug))
            executable_split_cohort = self.set_stage_runnable(
                stage=stage_split_cohort,
                runnable=runnable_split_cohort)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_split_cohort.dependencies.append(runnable_process_cohort.name)

            reference_split_cohort = runnable_split_cohort.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the GATK SelectVariants analysis to split multi-sample VCF files into one per sample.

            runnable_step = runnable_split_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_gatk_select_variants',
                    java_temporary_path=runnable_split_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_process_cohort.threads))
            runnable_step.add_gatk_option(key='variant', value=file_path_process_cohort.annotated_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort.sample_vcf)
            runnable_step.add_gatk_option(key='sample_name', value=sample.name)
            runnable_step.add_gatk_switch(key='excludeNonVariants')

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_split_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_gatk_variants_to_table',
                    java_temporary_path=runnable_split_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='variant', value=file_path_split_cohort.sample_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort.sample_tsv)
            runnable_step.add_gatk_switch(key='allowMissingData')
            runnable_step.add_gatk_switch(key='showFiltered')
            # Set of standard VCF fields.
            runnable_step.add_gatk_option(key='fields', value='CHROM', override=True)
            runnable_step.add_gatk_option(key='fields', value='POS', override=True)
            runnable_step.add_gatk_option(key='fields', value='ID', override=True)
            runnable_step.add_gatk_option(key='fields', value='REF', override=True)
            runnable_step.add_gatk_option(key='fields', value='ALT', override=True)
            runnable_step.add_gatk_option(key='fields', value='QUAL', override=True)
            runnable_step.add_gatk_option(key='fields', value='FILTER', override=True)
            #
            runnable_step.add_gatk_option(key='fields', value='AF', override=True)
            runnable_step.add_gatk_option(key='fields', value='VQSLOD', override=True)
            runnable_step.add_gatk_option(key='fields', value='culprit', override=True)
            # GATK Haplotype Caller genotype fields: GT:AD:DP:GQ:PL
            runnable_step.add_gatk_option(key='genotypeFields', value='AD', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='DP', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='GQ', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='GT', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='PGT', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='PID', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='PL', override=True)
            # Set of snpEff fields.
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_EFFECT', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_IMPACT', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_FUNCTIONAL_CLASS', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_CODON_CHANGE', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_AMINO_ACID_CHANGE', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_GENE_NAME', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_GENE_BIOTYPE', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_TRANSCRIPT_ID', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_EXON_ID', override=True)

            # Automatically add all fields defined for the Variant Annotator resources, above.
            for annotation_resource in self.annotation_resources_dict.keys():
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        runnable_step.add_gatk_option(
                            key='fields',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

        ###################################################
        # Step 9: Summarise the variant calling analysis. #
        ###################################################
        #
        # bsfR bsf_variant_calling_summary.R    (summary)

        prefix_summary = '_'.join((stage_summary.name, self.cohort_name))

        file_path_summary = FilePathSummary(prefix=prefix_summary)

        # Create a Runnable and Executable for summarising the cohort.

        runnable_summary = self.add_runnable(
            runnable=Runnable(
                name=prefix_summary,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory,
                cache_path_dict=self._cache_path_dict,
                file_path_object=file_path_summary,
                debug=self.debug))
        executable_summary = self.set_stage_runnable(
            stage=stage_summary,
            runnable=runnable_summary)
        # Set dependencies on preceding Runnable.name or Executable.name objects.
        for runnable_diagnose_sample in vc_runnable_diagnose_sample_list:
            executable_summary.dependencies.append(runnable_diagnose_sample.name)

        # Run the bsfR script to summarise the variant calling procedure.

        runnable_step = runnable_summary.add_runnable_step(
            runnable_step=RunnableStep(
                name='summary',
                program='bsf_variant_calling_summary.R'))
        """ @type runnable_step: RunnableStep """
        runnable_step.add_option_long(key='prefix', value=prefix_summary)

        #####################################
        # Step 10: Somatic variant calling. #
        #####################################
        #
        # GATK MuTect2          (somatic_gatk_mutect2)
        # snpEff                (somatic_snpeff)
        # GATK VariantAnnotator (somatic_gatk_variant_annotator)
        # GATK VariantsToTable  (somatic_gatk_variants_to_table)

        key_list = self.comparisons.keys()
        key_list.sort(cmp=lambda x, y: cmp(x, y))

        if self.debug > 0:
            print 'Somatic variant calling: ' + repr(key_list)

        for key in key_list:
            comparison = self.comparisons[key]
            """ @type comparison: VariantCallingGATKComparison """

            prefix_somatic = '_'.join((stage_somatic.name, key))

            # Somatic variant calling-specific file paths

            file_path_somatic = FilePathSomatic(prefix=prefix_somatic)

            # Run the GATK MuTect2 analysis to characterise somatic variants in a scatter gather approach.

            vc_runnable_somatic_scatter_list = list()
            """ @type vc_runnable_somatic_scatter_list: list[Runnable] """

            # If the somatic index file already exists, create the Runnable objects, but do not submit their
            # corresponding Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_somatic.annotated_tbi)):
                submit_runnable = False
            else:
                submit_runnable = True

            runnable_somatic_scatter = None
            for tile_index in range(0, len(self._tile_region_somatic_list)):
                prefix_somatic_scatter = '_'.join((
                    stage_somatic.name, key, 'scatter', str(tile_index)))

                file_path_somatic_scatter = FilePathSomaticScatterGather(prefix=prefix_somatic_scatter)

                runnable_somatic_scatter = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_somatic_scatter,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        file_path_object=file_path_somatic_scatter,
                        debug=self.debug))
                executable_somatic_scatter = self.set_stage_runnable(
                    stage=stage_somatic,
                    runnable=runnable_somatic_scatter)
                # Set dependencies on preceding Runnable.name or Executable.name objects.
                if comparison.normal_sample is not None:
                    executable_somatic_scatter.dependencies.append(
                        '_'.join((stage_process_sample.name, comparison.normal_sample.name)))
                if comparison.tumor_sample is not None:
                    executable_somatic_scatter.dependencies.append(
                        '_'.join((stage_process_sample.name, comparison.tumor_sample.name)))
                executable_somatic_scatter.submit = submit_runnable

                vc_runnable_somatic_scatter_list.append(runnable_somatic_scatter)

                reference_somatic_scatter = runnable_somatic_scatter.get_absolute_cache_file_path(
                    file_path=self.bwa_genome_db)

                runnable_step = runnable_somatic_scatter.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='somatic_gatk_mutect2_scatter',
                        java_temporary_path=runnable_somatic_scatter.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        gatk_classpath=self.classpath_gatk))
                """ @type runnable_step: RunnableStepGATK """
                runnable_step.add_gatk_option(key='analysis_type', value='MuTect2')
                runnable_step.add_gatk_option(key='reference_sequence', value=reference_somatic_scatter)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                # if self.interval_padding:
                #     runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for region_tuple in self._tile_region_somatic_list[tile_index]:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if region_tuple[0]:
                        runnable_step.add_gatk_option(
                            key='intervals',
                            value='{:s}:{:d}-{:d}'.format(region_tuple[0], region_tuple[1], region_tuple[2]),
                            override=True)
                if self.known_sites_discovery:
                    runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
                for file_path in self.known_somatic_discovery:
                    runnable_step.add_gatk_option(key='cosmic', value=file_path, override=True)

                # Find and add the FilePathProcessSample object for the 'normal' Sample object.
                if comparison.normal_sample is not None:
                    file_path_object = self.runnable_dict['_'.join((
                        stage_process_sample.name,
                        comparison.normal_sample.name))].file_path_object
                    """ @type file_path_object: FilePathProcessSample """
                    runnable_step.add_gatk_option(key='input_file:normal', value=file_path_object.realigned_bam)
                elif comparison.panel_of_normal_path is not None:
                    runnable_step.add_gatk_option(key='normal_panel', value=comparison.panel_of_normal_path)

                # Find and add the FilePathProcessSample object for the 'tumor' Sample object.
                if comparison.tumor_sample is not None:
                    file_path_object = self.runnable_dict['_'.join((
                        stage_process_sample.name,
                        comparison.tumor_sample.name))].file_path_object
                    """ @type file_path_object: FilePathProcessSample """
                    runnable_step.add_gatk_option(key='input_file:tumor', value=file_path_object.realigned_bam)

                runnable_step.add_gatk_option(key='out', value=file_path_somatic_scatter.partial_vcf)

            # Gather

            # If there is only one tile, no need to gather, just rename the file and return the Runnable.

            if len(self._tile_region_somatic_list) == 1:
                file_path_somatic_scatter = runnable_somatic_scatter.file_path_object
                """ @type file_path_somatic_scatter: FilePathSomaticScatterGather """
                # Add sample-specific keys to the file path dictionary.
                file_path_somatic_scatter.somatic_vcf = file_path_somatic.somatic_vcf
                file_path_somatic_scatter.somatic_idx = file_path_somatic.somatic_idx

                runnable_somatic_scatter.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='somatic_gather_move_vcf',
                        source_path=file_path_somatic_scatter.partial_vcf,
                        target_path=file_path_somatic_scatter.somatic_vcf))

                runnable_somatic_scatter.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='somatic_gather_move_tbi',
                        source_path=file_path_somatic_scatter.partial_tbi,
                        target_path=file_path_somatic_scatter.somatic_idx))

                # return runnable_somatic_scatter
                runnable_somatic_gather = runnable_somatic_scatter
            else:
                # Second, gather by the number of chunks on the partitioned genome tile index list.

                # Second, Merge chunks hierarchically.
                # Initialise a list of Runnable objects and indices for the hierarchical merge.
                vc_runnable_somatic_gather_list = vc_runnable_somatic_scatter_list
                vc_somatic_gather_index_list = range(0, len(self._tile_region_somatic_list))
                runnable_somatic_gather = None  # Global variable to keep and return the last Runnable.
                level = 0
                while len(vc_somatic_gather_index_list) > 1:
                    temporary_gather_runnable_list = list()
                    """ @type temporary_gather_runnable_list: list[Runnable] """
                    temporary_gather_index_list = list()
                    """ @type temporary_gather_index_list: list[int] """
                    # Partition the index list into chunks of given size.
                    partition_list = [vc_somatic_gather_index_list[offset:offset + self.number_of_chunks_cohort]
                                      for offset in range(0,
                                                          len(vc_somatic_gather_index_list),
                                                          self.number_of_chunks_cohort)]

                    for partition_index in range(0, len(partition_list)):
                        chunk_index_list = partition_list[partition_index]
                        # The file prefix includes the level and partition index.
                        prefix_somatic_gather = '_'.join(
                            (stage_somatic.name, key, 'gather', str(level), str(partition_index)))

                        file_path_somatic_gather = FilePathSomaticScatterGather(prefix=prefix_somatic_gather)

                        runnable_somatic_gather = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_somatic_gather,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                file_path_object=file_path_somatic_gather,
                                debug=self.debug))
                        executable_somatic_gather = self.set_stage_runnable(
                            stage=stage_somatic,
                            runnable=runnable_somatic_gather)
                        executable_somatic_gather.submit = submit_runnable
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_gather_runnable_list.append(runnable_somatic_gather)
                        temporary_gather_index_list.append(partition_index)

                        reference_somatic_gather = runnable_somatic_gather.get_absolute_cache_file_path(
                            file_path=self.bwa_genome_db)

                        # GATK CatVariants by-passes the GATK engine and thus requires a completely different
                        # command line.
                        runnable_step = runnable_somatic_gather.add_runnable_step(
                            runnable_step=RunnableStepJava(
                                name='somatic_gatk_cat_variants',
                                sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                                java_temporary_path=runnable_somatic_gather.get_relative_temporary_directory_path,
                                java_heap_maximum='Xmx4G'))
                        """ @type runnable_step: RunnableStepJava """
                        runnable_step.add_option_short(
                            key='classpath',
                            value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                        sub_command = runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        sub_command.add_option_long(
                            key='reference',
                            value=reference_somatic_gather)
                        sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_somatic_gather.partial_vcf)
                        sub_command.add_switch_long(key='assumeSorted')
                        # Finally, add RunnableStep options, obsolete files and Executable dependencies per chunk index.
                        for chunk_index in chunk_index_list:
                            runnable_object = vc_runnable_somatic_gather_list[chunk_index]
                            file_path_object = runnable_object.file_path_object
                            """ @type file_path_object: FilePathSomaticScatterGather """
                            # Set GATK option variant
                            sub_command.add_option_long(
                                key='variant',
                                value=file_path_object.partial_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            runnable_step.obsolete_file_path_list.append(file_path_object.partial_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            runnable_step.obsolete_file_path_list.append(file_path_object.partial_tbi)
                            # Set dependencies on preceding Runnable.name or Executable.name objects.
                            executable_somatic_gather.dependencies.append(runnable_object.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    vc_runnable_somatic_gather_list = temporary_gather_runnable_list
                    vc_somatic_gather_index_list = temporary_gather_index_list
                    level += 1
                else:
                    # For the last instance, additionally rename the final file.
                    file_path_somatic_gather = runnable_somatic_gather.file_path_object
                    """ @type file_path_somatic_gather: FilePathSomaticScatterGather """

                    # Add sample-specific keys to the file path dictionary.
                    file_path_somatic_gather.somatic_vcf = file_path_somatic.somatic_vcf
                    file_path_somatic_gather.somatic_idx = file_path_somatic.somatic_idx

                    runnable_somatic_gather.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='somatic_gather_move_vcf',
                            source_path=file_path_somatic_gather.partial_vcf,
                            target_path=file_path_somatic_gather.somatic_vcf))

                    runnable_somatic_gather.add_runnable_step(
                        runnable_step=RunnableStepMove(
                            name='somatic_gather_move_tbi',
                            source_path=file_path_somatic_gather.partial_tbi,
                            target_path=file_path_somatic_gather.somatic_idx))

                    # return runnable_somatic_gather

            runnable_somatic = self.add_runnable(
                runnable=Runnable(
                    name=prefix_somatic,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_somatic,
                    debug=self.debug))
            executable_somatic = self.set_stage_runnable(
                stage=stage_somatic,
                runnable=runnable_somatic)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_somatic.dependencies.append(runnable_somatic_gather.name)

            reference_somatic = runnable_somatic.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the snpEff tool for functional variant annotation.

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStep(
                    name='somatic_snpeff',
                    program='java',
                    sub_command=Command(program='eff')))
            """ @type runnable_step: RunnableStep """
            runnable_step.add_switch_short(
                key='d64')
            runnable_step.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_snpeff, 'snpEff.jar'))
            runnable_step.add_switch_short(
                key='Xmx4G')
            runnable_step.add_option_pair(
                key='-Djava.io.tmpdir',
                value=runnable_somatic.get_relative_temporary_directory_path)
            runnable_step.stdout_path = file_path_somatic.snpeff_vcf

            sub_command = runnable_step.sub_command
            sub_command.add_switch_short(key='download')
            sub_command.add_option_short(key='o', value='gatk')
            sub_command.add_option_short(key='stats', value=file_path_somatic.snpeff_stats)
            sub_command.add_option_short(key='config', value=os.path.join(self.classpath_snpeff, 'snpEff.config'))

            sub_command.arguments.append(self.snpeff_genome_version)
            sub_command.arguments.append(file_path_somatic.somatic_vcf)

            # Automatically compress and index the snpEff VCF file with bgzip and tabix, respectively.
            # TODO: It would be better for the file system, if output could be directly piped into bgzip.

            runnable_somatic.add_runnable_step(
                runnable_step=RunnableStep(
                    name='snpeff_bgzip',
                    program='bgzip',
                    arguments=[file_path_somatic.snpeff_vcf]))

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStep(
                    name='snpeff_tabix',
                    program='tabix',
                    arguments=[file_path_somatic.snpeff_vcf_bgz]))
            """ @type runnable_step: RunnableStep """
            runnable_step.add_option_long(key='preset', value='vcf')

            # Run the GATK VariantAnnotator analysis.

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_variant_annotator',
                    java_temporary_path=runnable_somatic.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantAnnotator')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_somatic)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # runnable_step.add_gatk_option(key='num_threads', value=str(stage_somatic.threads))
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)

            # Add annotation resources and their corresponding expression options.
            for annotation_resource in self.annotation_resources_dict.keys():
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    runnable_step.add_gatk_option(
                        key=':'.join(('resource', annotation_resource)),
                        value=self.annotation_resources_dict[annotation_resource][0])
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        runnable_step.add_gatk_option(
                            key='expression',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

            runnable_step.add_gatk_option(key='variant', value=file_path_somatic.somatic_vcf)
            # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
            # runnable_step.add_gatk_option(key='annotation', value='AlleleBalanceBySample')
            runnable_step.add_gatk_option(key='annotation', value='SnpEff')
            runnable_step.add_gatk_option(key='snpEffFile', value=file_path_somatic.snpeff_vcf_bgz)
            runnable_step.add_gatk_option(key='out', value=file_path_somatic.annotated_vcf)

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_variants_to_table',
                    java_temporary_path=runnable_somatic.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_somatic)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='variant', value=file_path_somatic.annotated_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_somatic.annotated_tsv)
            runnable_step.add_gatk_switch(key='allowMissingData')
            runnable_step.add_gatk_switch(key='showFiltered')
            # Set of standard VCF fields.
            runnable_step.add_gatk_option(key='fields', value='CHROM', override=True)
            runnable_step.add_gatk_option(key='fields', value='POS', override=True)
            runnable_step.add_gatk_option(key='fields', value='ID', override=True)
            runnable_step.add_gatk_option(key='fields', value='REF', override=True)
            runnable_step.add_gatk_option(key='fields', value='ALT', override=True)
            runnable_step.add_gatk_option(key='fields', value='QUAL', override=True)
            runnable_step.add_gatk_option(key='fields', value='FILTER', override=True)
            # GATK MuTect2 info fields.
            runnable_step.add_gatk_option(key='fields', value='NLOD', override=True)
            runnable_step.add_gatk_option(key='fields', value='TLOD', override=True)
            # GATK MuTect2 genotype fields: GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1
            runnable_step.add_gatk_option(key='genotypeFields', value='AD', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='AF', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='ALT_F1R2', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='ALT_F2R1', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='DP', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='FOXOG', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='GQ', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='GT', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='PGT', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='PID', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='PL', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='QSS', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='REF_F1R2', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='REF_F2R1', override=True)
            # Set of snpEff fields.
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_EFFECT', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_IMPACT', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_FUNCTIONAL_CLASS', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_CODON_CHANGE', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_AMINO_ACID_CHANGE', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_GENE_NAME', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_GENE_BIOTYPE', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_TRANSCRIPT_ID', override=True)
            runnable_step.add_gatk_option(key='fields', value='SNPEFF_EXON_ID', override=True)

            # Automatically add all fields defined for the Variant Annotator resources, above.
            for annotation_resource in self.annotation_resources_dict.keys():
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        runnable_step.add_gatk_option(
                            key='fields',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

        return

    def report(self):
        """Create a C{bsf.analyses.variant_calling.VariantCallingGATK} report in HTML format and a
        UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        def report_create_directory(path):
            """Private function to create a directory avoiding race conditions.

            @param path: Path
            @type path: str | unicode
            @return: Path
            @rtype: str | unicode
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
            @type source_path: str | unicode
            @param target_path: Target path
            @type target_path: str | unicode
            @return:
            @rtype:
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

        # Create simpler symbolic links in structured sub-directoreies to facilitate file system browsing.

        # Create a result directory as concatenation of bsf.Analysis.genome_version and bsf.Analysis.prefix.

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

            runnable_process_sample = self.runnable_dict[
                '_'.join((self.stage_name_process_sample, sample.name))]
            file_path_process_sample = runnable_process_sample.file_path_object
            """ @type file_path_process_sample: FilePathProcessSample """

            # runnable_diagnose_sample = self.runnable_dict[
            #     '_'.join((self.stage_name_diagnose_sample, sample.name))]
            # file_path_diagnose_sample = runnable_diagnose_sample.file_path_object
            """ @type file_path_diagnose_sample: FilePathDiagnoseSample """

            runnable_split_cohort = self.runnable_dict[
                '_'.join((self.stage_name_split_cohort, sample.name))]
            file_path_split_cohort = runnable_split_cohort.file_path_object
            """ @type file_path_split_cohort: FilePathSplitCohort """

            directory_sample = report_create_directory(
                path=os.path.join(directory_results_by_sample, sample.name))

            for key, extension in (
                    ('realigned_bam', '.bam'),
                    ('realigned_bai', '.bam.bai'),
                    ('realigned_md5', '.bam.md5')):
                report_create_symbolic_link(
                    source_path=os.path.relpath(
                        os.path.join(self.genome_directory, getattr(file_path_process_sample, key)),
                        directory_sample),
                    target_path=os.path.join(directory_sample, sample.name + extension))

            for key, extension in (
                    ('sample_vcf', '.vcf.gz'),
                    ('sample_idx', '.vcf.gz.tbi'),
                    ('sample_tsv', '.tsv')):
                report_create_symbolic_link(
                    source_path=os.path.relpath(
                        os.path.join(self.genome_directory, getattr(file_path_split_cohort, key)),
                        directory_sample),
                    target_path=os.path.join(directory_sample, sample.name + extension))

            directory_alignments = report_create_directory(
                path=os.path.join(directory_results_by_type, 'alignments'))

            for key, extension in (
                    ('realigned_bam', '.bam'),
                    ('realigned_bai', '.bam.bai'),
                    ('realigned_md5', '.bam.md5')):
                report_create_symbolic_link(
                    source_path=os.path.relpath(
                        os.path.join(self.genome_directory, getattr(file_path_process_sample, key)),
                        directory_alignments),
                    target_path=os.path.join(directory_alignments, sample.name + extension))

            directory_variants = report_create_directory(
                path=os.path.join(directory_results_by_type, 'variants'))

            for key, extension in (
                    ('sample_vcf', '.vcf.gz'),
                    ('sample_idx', '.vcf.gz.tbi'),
                    ('sample_tsv', '.tsv')):
                report_create_symbolic_link(
                    source_path=os.path.relpath(
                        os.path.join(self.genome_directory, getattr(file_path_split_cohort, key)),
                        directory_variants),
                    target_path=os.path.join(directory_variants, sample.name + extension))

        runnable_process_cohort = self.runnable_dict[
            '_'.join((self.stage_name_process_cohort, self.cohort_name))]
        file_path_process_cohort = runnable_process_cohort.file_path_object
        """ @type file_path_process_cohort: FilePathProcessCohort """

        for key, extension in (
                ('annotated_vcf', '_annotated.vcf.gz'),
                ('annotated_tbi', '_annotated.vcf.gz.tbi'),
                ('snpeff_vcf_bgz', '_snpeff.vcf.gz'),
                ('snpeff_vcf_tbi', '_snpeff.vcf.gz.tbi'),
                ('snpeff_stats', '_snpeff_summary.html')):
            report_create_symbolic_link(
                source_path=os.path.relpath(
                    os.path.join(self.genome_directory, getattr(file_path_process_cohort, key)),
                    directory_results_by_cohort),
                target_path=os.path.join(directory_results_by_cohort, self.cohort_name + extension))

        for comparison_name in self.comparisons.keys():
            runnable_somatic = self.runnable_dict['_'.join((self.stage_name_somatic, comparison_name))]
            file_path_somatic = runnable_somatic.file_path_object
            """ @type file_path_somatic: FilePathSomatic """

            for key, extension in (
                    ('somatic_vcf', '.vcf.gz'),
                    ('somatic_idx', '.vcf.gz.tbi'),
                    ('snpeff_vcf_bgz', '_snpeff.vcf.gz'),
                    ('snpeff_vcf_tbi', '_snpeff.vcf.gz.tbi'),
                    ('snpeff_genes', '_snpeff_summary.genes.txt'),
                    ('snpeff_stats', '_snpeff_summary.html'),
                    ('annotated_vcf', '_annotated.vcf.gz'),
                    ('annotated_tbi', '_annotated.vcf.gz.tbi'),
                    ('annotated_tsv', '_annotated.tsv')):
                report_create_symbolic_link(
                    source_path=os.path.relpath(
                        os.path.join(self.genome_directory, getattr(file_path_somatic, key)),
                        directory_results_by_pair),
                    target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

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

        report_list += '<h2 id="genome_browsing">Genome Browsing</h2>\n'
        report_list += '\n'

        # Resolve an eventual alias for the UCSC genome assembly name.

        if default.genome_aliases_ucsc_dict is not None and self.genome_version in default.genome_aliases_ucsc_dict:
            ucsc_genome_version = default.genome_aliases_ucsc_dict[self.genome_version]
        else:
            ucsc_genome_version = self.genome_version

        options_dict = {
            'db': ucsc_genome_version,
            'hubUrl': Default.url_absolute_projects() + '/' + link_name + '/variant_calling_hub.txt',
        }

        # The above code for resolving a UCSC Genome Browser genome assembly alias could be centralised in Analysis.
        report_list += '<p id="ucsc_track_hub">'
        report_list += 'UCSC Genome Browser Track Hub '
        report_list += '<a href="' + self.ucsc_track_url(options_dict=options_dict) + '">' + self.project_name + '</a>.'
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<h2 id="read_group_and_sample_level">Read Group and Sample Level</h2>\n'
        report_list += '\n'
        report_list += '<table id="read_group_and_sample_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Sample</th>\n'
        report_list += '<th>Variants</th>\n'
        report_list += '<th>Alignments</th>\n'
        report_list += '<th>Read Group</th>\n'
        report_list += '<th>Duplicate Metrics</th>\n'
        report_list += '<th>Alignment Summary Metrics</th>\n'
        report_list += '<th>Hybrid Selection Metrics</th>\n'
        report_list += '<th>Non-Callable Loci</th>\n'
        report_list += '<th>Non-Callable Summary</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        # Group via UCSC super tracks.

        hub_list += 'track Alignments\n'
        hub_list += 'shortLabel Alignments\n'
        hub_list += 'longLabel BWA NGS read alignments\n'
        hub_list += 'superTrack on show\n'
        hub_list += 'group alignments\n'
        hub_list += '\n'

        hub_list += 'track Callable\n'
        hub_list += 'shortLabel Callable\n'
        hub_list += 'longLabel Callable\n'
        hub_list += 'superTrack on show\n'
        hub_list += 'group callable\n'
        hub_list += '\n'

        hub_list += 'track Variants\n'
        hub_list += 'shortLabel Variants\n'
        hub_list += 'longLabel Variant calls\n'
        hub_list += 'superTrack on show\n'
        hub_list += 'group variants\n'
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

            runnable_process_sample = self.runnable_dict[
                '_'.join((self.stage_name_process_sample, sample.name))]
            file_path_process_sample = runnable_process_sample.file_path_object
            """ @type file_path_process_sample: FilePathProcessSample """

            runnable_diagnose_sample = self.runnable_dict[
                '_'.join((self.stage_name_diagnose_sample, sample.name))]
            file_path_diagnosis = runnable_diagnose_sample.file_path_object
            """ @type file_path_diagnosis: FilePathDiagnoseSample """

            runnable_split_cohort = self.runnable_dict[
                '_'.join((self.stage_name_split_cohort, sample.name))]
            file_path_split_cohort = runnable_split_cohort.file_path_object
            """ @type file_path_split_cohort: FilePathSplitCohort """

            # Alignment track
            # Common settings

            hub_list += 'track ' + sample.name + '_alignments\n'
            hub_list += 'type bam\n'
            hub_list += 'shortLabel ' + sample.name + '_alignments\n'
            hub_list += 'longLabel ' + sample.name + 'BWA NGS read alignments\n'
            hub_list += 'bigDataUrl ' + file_path_process_sample.realigned_bam + '\n'
            # track_output += 'html ...\n'
            hub_list += 'visibility squish\n'

            # Common optional settings.

            hub_list += 'color 0,0,0\n'

            # Compressed Sequence Alignment track settings.

            # None so far.

            # Composite track settings.

            hub_list += 'parent Alignments\n'
            hub_list += '\n'

            # Callable track

            if os.path.isfile(
                    os.path.join(
                        self.genome_directory,
                        file_path_diagnosis.callable_bb)):
                hub_list += 'track ' + sample.name + '_callable\n'
                hub_list += 'type bigBed\n'
                hub_list += 'shortLabel ' + sample.name + '_callable\n'
                hub_list += 'longLabel ' + sample.name + 'callable\n'
                hub_list += 'bigDataUrl ' + file_path_diagnosis.callable_bb + '\n'
                # track_output += 'html ...\n'
                hub_list += 'visibility squish\n'

                # Common optional settings.

                hub_list += 'color 0,0,0\n'

                # Composite track settings.

                hub_list += 'parent Callable\n'
                hub_list += '\n'

            # Variant track

            hub_list += 'track ' + sample.name + '_variants\n'
            hub_list += 'type vcfTabix\n'
            hub_list += 'shortLabel ' + sample.name + '_variants\n'
            hub_list += 'longLabel ' + sample.name + ' variant calls\n'
            hub_list += 'bigDataUrl ' + file_path_split_cohort.sample_vcf + '\n'
            # track_output += 'html ...\n'
            hub_list += 'visibility dense\n'

            # vcfTabix specific settings

            # None so far.

            # Composite track settings.

            hub_list += 'parent Variants\n'
            hub_list += '\n'

            report_list += '<tr>\n'
            report_list += '<td class="left">' + sample.name + '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_split_cohort.sample_vcf + '">VCF</a>&nbsp;'
            report_list += '<a href="' + file_path_split_cohort.sample_idx + '">TBI</a>&nbsp;'
            report_list += '<a href="' + file_path_split_cohort.sample_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_process_sample.realigned_bam + '">BAM</a>&nbsp;'
            report_list += '<a href="' + file_path_process_sample.realigned_bai + '">BAI</a>'
            report_list += '</td>\n'
            report_list += '<td class="left"></td>\n'  # Read Group
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_process_sample.duplicate_metrics + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_process_sample.alignment_summary_metrics + '">TSV</a>'
            report_list += '</td>\n'
            if os.path.isfile(
                    os.path.join(
                        self.genome_directory,
                        file_path_diagnosis.hybrid_selection_metrics)):
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_diagnosis.hybrid_selection_metrics + '">TSV</a>'
                report_list += '</td>\n'
            else:
                report_list += '<td class="center"></td>\n'
            report_list += '<td class="center">'
            if os.path.isfile(
                    os.path.join(
                        self.genome_directory,
                        file_path_diagnosis.non_callable_loci_tsv)):
                report_list += '<a href="' + file_path_diagnosis.callable_bed + '">BED</a>&nbsp;'
                report_list += '<a href="' + file_path_diagnosis.callable_bb + '">BigBED</a>&nbsp;'
                # Do not link the more complex file_path_diagnosis.non_callable_regions_tsv
                # file for the moment.
                report_list += '<a href="' + file_path_diagnosis.non_callable_regions_tsv + '">TSV</a>'
            report_list += '</td>\n'
            if os.path.isfile(
                    os.path.join(
                        self.genome_directory,
                        file_path_diagnosis.non_callable_summary_tsv)):
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_diagnosis.non_callable_summary_tsv + '">TSV</a>'
                report_list += '</td>\n'
            else:
                report_list += '<td class="center"></td>\n'
            report_list += '</tr>\n'

            for paired_reads_name in paired_reads_name_list:
                runnable_process_lane = self.runnable_dict[
                    '_'.join((self.stage_name_process_lane, paired_reads_name))]
                file_path_read_group = runnable_process_lane.file_path_object
                """ @type file_path_read_group: FilePathProcessReadGroup """

                report_list += '<tr>\n'
                report_list += '<td class="left"></td>\n'  # Sample
                report_list += '<td class="center"></td>\n'  # Variants VCF, TBI, TSV
                report_list += '<td class="center"></td>\n'  # Alignments BAM, BAI
                report_list += '<td class="left">' + paired_reads_name + '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_read_group.duplicate_metrics + '">TSV</a>'
                report_list += '</td>\n'
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_read_group.alignment_summary_metrics + '">TSV</a>'
                report_list += '</td>\n'
                report_list += '<td class="center"></td>\n'  # Hybrid Selection Metrics
                report_list += '<td class="center"></td>\n'  # Non-Callable Loci
                report_list += '<td class="center"></td>\n'  # Non-Callable Summary
                report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        report_list += '<h2 id="cohort_level">Cohort Level</h2>\n'
        report_list += '\n'
        report_list += '<table id="cohort_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Cohort</th>\n'
        report_list += '<th>Information</th>\n'
        report_list += '<th>Comment</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        runnable_process_cohort = self.runnable_dict[
            '_'.join((self.stage_name_process_cohort, self.cohort_name))]
        file_path_process_cohort = runnable_process_cohort.file_path_object
        """ @type file_path_process_cohort: FilePathProcessCohort """

        report_list += '<tr>\n'
        report_list += '<td class="left">' + self.cohort_name + '</td>\n'
        report_list += '<td class="left">'
        report_list += '<a href="' + file_path_process_cohort.snpeff_stats + '">snpEff Summary Statistics</a>'
        report_list += '</td>\n'
        report_list += '<td class="left"></td>\n'
        report_list += '</tr>\n'

        report_list += '<tr>\n'
        report_list += '<td class="left">' + self.cohort_name + '</td>\n'
        report_list += '<td class="left">'
        report_list += 'snpEff-annotated multi-sample '
        report_list += '<a href="' + file_path_process_cohort.snpeff_vcf_bgz + '">VCF</a> and '
        report_list += '<a href="' + file_path_process_cohort.snpeff_vcf_tbi + '">TBI</a>'
        report_list += '</td>\n'
        report_list += '<td class="left">Functional annotation of all splice variants</td>\n'
        report_list += '</tr>\n'

        report_list += '<tr>\n'
        report_list += '<td class="left">' + self.cohort_name + '</td>\n'
        report_list += '<td class="left">'
        report_list += 'GATK-annotated multi-sample '
        report_list += '<a href="' + file_path_process_cohort.annotated_vcf + '">VCF</a> and '
        report_list += '<a href="' + file_path_process_cohort.annotated_tbi + '">TBI</a>'
        report_list += '</td>\n'
        report_list += '<td class="left">Functional annotation of only the most severely affected splice variant</td>\n'
        report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        # Somatic variant calling.

        key_list = self.comparisons.keys()
        key_list.sort(cmp=lambda x, y: cmp(x, y))

        if len(key_list):
            report_list += '<h2 id="somatic_variants">Somatic Variants</h2>\n'
            report_list += '\n'
            report_list += '<table id="somatic_variants_table">\n'
            report_list += '<thead>\n'
            report_list += '<tr>\n'
            report_list += '<th>Comparison</th>\n'
            report_list += '<th>Annotated VCF</th>\n'
            report_list += '<th>Annotated TSV</th>\n'
            report_list += '<th>snpEff Summary Statistics</th>\n'
            report_list += '<th>snpEff Genes</th>\n'
            report_list += '</tr>\n'
            report_list += '</thead>\n'
            report_list += '<tbody>\n'

            for key in key_list:
                runnable_somatic = self.runnable_dict['_'.join((self.stage_name_somatic, key))]
                file_path_somatic = runnable_somatic.file_path_object
                """ @type file_path_somatic: FilePathSomatic """

                report_list += '<tr>\n'
                report_list += '<td>' + key + '</td>\n'
                report_list += '<td><a href="' + file_path_somatic.annotated_vcf + '">VCF</a></td>\n'
                report_list += '<td><a href="' + file_path_somatic.annotated_tsv + '">TSV</a></td>\n'
                report_list += '<td><a href="' + file_path_somatic.snpeff_stats + '">HTML</a></td>\n'
                report_list += '<td><a href="' + file_path_somatic.snpeff_genes + '">TXT</a></td>\n'
                report_list += '</tr>\n'

            report_list += '</tbody>\n'
            report_list += '</table>\n'
            report_list += '\n'

        report_list += '<h2 id="qc_plots">QC Plots</h2>\n'
        report_list += '\n'
        report_list += '<table id="qc_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Sample</th>\n'
        report_list += '<th>Read Group</th>\n'
        report_list += '<th>Metrics</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        runnable_summary = self.runnable_dict['_'.join((self.stage_name_summary, self.cohort_name))]
        file_path_summary = runnable_summary.file_path_object
        """ @type file_path_summary: FilePathSummary """

        # Alignment Summary - TSV
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.alignment_metrics_sample_tsv)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.alignment_metrics_sample_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.alignment_metrics_read_group_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">'
            report_list += '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' \
                           '#AlignmentSummaryMetrics">Alignment Summary</a>'
            report_list += '</td>\n'
            report_list += '</tr>\n'

        # Alignment Summary - Percent Aligned
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.alignment_percentage_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.alignment_percentage_sample_pdf + '">'
            report_list += '<img alt="Alignment Summary - Percent Aligned per Sample"'
            report_list += ' src="' + file_path_summary.alignment_percentage_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.alignment_percentage_read_group_pdf + '">'
            report_list += '<img alt="Alignment Summary - Percent Aligned per Read Group"'
            report_list += ' src="' + file_path_summary.alignment_percentage_read_group_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">Alignment Summary - Percent Aligned</td>\n'
            report_list += '</tr>\n'

        # Alignment Summary - Reads Aligned
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.alignment_absolute_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.alignment_absolute_sample_pdf + '">'
            report_list += '<img alt="Alignment Summary - Reads Aligned per Sample"'
            report_list += ' src="' + file_path_summary.alignment_absolute_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.alignment_absolute_read_group_pdf + '">'
            report_list += '<img alt="Alignment Summary - Reads Aligned per Read Group"'
            report_list += ' src="' + file_path_summary.alignment_absolute_read_group_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">Alignment Summary - Reads Aligned</td>\n'
            report_list += '</tr>\n'

        # Duplication - TSV
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.duplication_metrics_sample_tsv)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.duplication_metrics_sample_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="center"></td>\n'
            report_list += '<td class="left">'
            report_list += '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' \
                           '#DuplicationMetrics">Duplication</a>'
            report_list += '</td>\n'
            report_list += '</tr>\n'

        # Duplication - Fraction
        if os.path.exists(
                os.path.join(
                    self.genome_directory,
                    file_path_summary.duplication_percentage_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.duplication_percentage_sample_pdf + '">'
            report_list += '<img alt="Duplication - Duplicated Reads per Sample"'
            report_list += ' src="' + file_path_summary.duplication_percentage_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center"></td>\n'
            report_list += '<td class="left">Duplication - Fraction</td>\n'
            report_list += '</tr>\n'

        # Duplication - Levels
        if os.path.exists(
                os.path.join(
                    self.genome_directory,
                    file_path_summary.duplication_levels_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.duplication_levels_sample_pdf + '">'
            report_list += '<img alt="Duplication - Duplication Levels per Sample"'
            report_list += ' src="' + file_path_summary.duplication_levels_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center"></td>\n'
            report_list += '<td class="left">Duplication - Levels</td>\n'
            report_list += '</tr>\n'

        # Hybrid Selection - TSV
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.hybrid_metrics_sample_tsv)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_metrics_sample_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_metrics_read_group_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">'
            report_list += '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' \
                           '#HsMetrics">Hybrid Selection</a>'
            report_list += '</td>\n'
            report_list += '</tr>\n'

        # Hybrid Selection - Target Coverage Levels
        if os.path.exists(
                os.path.join(
                    self.genome_directory,
                    file_path_summary.hybrid_coverage_levels_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_coverage_levels_sample_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Mean Target Coverage Levels per Sample"'
            report_list += ' src="' + file_path_summary.hybrid_coverage_levels_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_coverage_levels_read_group_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Mean Target Coverage Levels per Read Group"'
            report_list += ' src="' + file_path_summary.hybrid_coverage_levels_read_group_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">Hybrid Selection - Mean Target Coverage Levels</td>\n'
            report_list += '</tr>\n'

        # Hybrid Selection - Mean Target Coverage
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.hybrid_coverage_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_coverage_sample_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Mean Target Coverage per Sample"'
            report_list += ' src="' + file_path_summary.hybrid_coverage_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_coverage_read_group_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Mean Target Coverage per Read Group"'
            report_list += ' src="' + file_path_summary.hybrid_coverage_read_group_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">Hybrid Selection - Mean Target Coverage</td>\n'
            report_list += '</tr>\n'

        # Hybrid Selection - Excluded Bases
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.hybrid_excluded_bases_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_excluded_bases_sample_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Percent Excluded Bases per Sample"'
            report_list += ' src="' + file_path_summary.hybrid_excluded_bases_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_excluded_bases_read_group_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Percent Excluded Bases per Read Group"'
            report_list += ' src="' + file_path_summary.hybrid_excluded_bases_read_group_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">Hybrid Selection - Percent Excluded Bases</td>\n'
            report_list += '</tr>\n'

        # Hybrid Selection - Percent Unique Reads
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.hybrid_unique_percentage_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_unique_percentage_sample_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Percent Unique Reads per Sample"'
            report_list += ' src="' + file_path_summary.hybrid_unique_percentage_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.hybrid_unique_percentage_read_group_pdf + '">'
            report_list += '<img alt="Hybrid Selection - Percent Unique Reads per Read Group"'
            report_list += ' src="' + file_path_summary.hybrid_unique_percentage_read_group_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="left">Hybrid Selection - Percent Unique Reads</td>\n'
            report_list += '</tr>\n'

        # Non-Callable Loci - TSV
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.non_callable_metrics_sample_tsv)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.non_callable_metrics_sample_tsv + '">TSV</a>'
            report_list += '</td>\n'
            report_list += '<td class="center"></td>\n'
            report_list += '<td class="left">Non-Callable Loci</td>\n'
            report_list += '</tr>\n'

        # Non-Callable Loci - Fraction
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.non_callable_percentage_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.non_callable_percentage_sample_pdf + '">'
            report_list += '<img alt="Non-Callable Loci - Fraction"'
            report_list += ' src="' + file_path_summary.non_callable_percentage_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center"></td>\n'
            report_list += '<td class="left">Non-Callable Loci - Fraction</td>\n'
            report_list += '</tr>\n'

        # Non-Callable Loci - Number
        if os.path.exists(os.path.join(
                self.genome_directory,
                file_path_summary.non_callable_absolute_sample_png)):
            report_list += '<tr>\n'
            report_list += '<td class="center">'
            report_list += '<a href="' + file_path_summary.non_callable_absolute_sample_pdf + '">'
            report_list += '<img alt="Non-Callable Loci - Number"'
            report_list += ' src="' + file_path_summary.non_callable_absolute_sample_png + '"'
            report_list += ' height="100" width="100" />'
            report_list += '</a>'
            report_list += '</td>\n'
            report_list += '<td class="center"></td>\n'
            report_list += '<td class="left">Non-Callable Loci - Number</td>\n'
            report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        self.report_to_file(content=report_list)
        self.ucsc_hub_to_file(content=hub_list)

        return
