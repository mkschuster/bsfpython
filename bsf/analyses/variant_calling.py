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
import pickle
import warnings

import pysam

from bsf import Analysis, FilePath, Runnable
from bsf.annotation import AnnotationSheet
from bsf.executables import BWA
from bsf.process import Command, Executable, RunnableStep, RunnableStepJava, RunnableStepPicard, RunnableStepLink, \
    RunnableStepMove
from bsf.standards import Configuration, Default, EnsemblVEP, JavaClassPath


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
        self.realigned_md5 = prefix + '_realigned.bam.md5'
        self.realigned_bam_bai = prefix + '_realigned.bam.bai'
        self.alignment_summary_metrics = prefix + '_alignment_summary_metrics.tsv'
        self.raw_variants_gvcf_vcf = prefix + '_raw_variants.g.vcf.gz'
        self.raw_variants_gvcf_tbi = prefix + '_raw_variants.g.vcf.gz.tbi'

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
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAnnotateSnpEff} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathAnnotateSnpEff, self).__init__(prefix=prefix)

        self.snpeff_vcf = prefix + '_snpeff.vcf'
        self.snpeff_vcf_bgz = prefix + '_snpeff.vcf.gz'
        self.snpeff_vcf_tbi = prefix + '_snpeff.vcf.gz.tbi'
        self.snpeff_genes = prefix + '_snpeff_summary.genes.txt'
        self.snpeff_stats = prefix + '_snpeff_summary.html'
        self.annotated_vcf = prefix + '_annotated.vcf.gz'
        self.annotated_tbi = prefix + '_annotated.vcf.gz.tbi'

        return


class FilePathAnnotateVEP(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathAnnotateVEP} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathAnnotateVEP, self).__init__(prefix=prefix)

        self.vep_statistics = prefix + '_vep_statistics.html'
        # Complete VEP set raw
        self.vep_complete_raw_vcf = prefix + '_vep_complete_raw.vcf'
        self.vep_complete_raw_vcf_bgz = prefix + '_vep_complete_raw.vcf.gz'
        self.vep_complete_raw_vcf_tbi = prefix + '_vep_complete_raw.vcf.gz.tbi'
        # Filtered VEP set raw
        self.vep_filtered_raw_vcf = prefix + '_vep_filtered_raw.vcf'
        self.vep_filtered_raw_vcf_bgz = prefix + '_vep_filtered_raw.vcf.gz'
        self.vep_filtered_raw_vcf_tbi = prefix + '_vep_filtered_raw.vcf.gz.tbi'
        # Complete VEP set VCF.Filter-converted
        self.vep_complete_vcf_bgz = prefix + '_vep_complete.vcf.gz'
        self.vep_complete_vcf_tbi = prefix + '_vep_complete.vcf.gz.tbi'
        # Filtered VEP set VCF.Filter-converted
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
        self.sample_tbi = prefix + '.vcf.gz.tbi'
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
        self.somatic_tbi = prefix + '_somatic.vcf.gz.tbi'

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

        self.somatic_vcf = prefix + '_somatic.vcf.gz'
        self.somatic_tbi = prefix + '_somatic.vcf.gz.tbi'

        return


class FilePathSplitSomatic(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.variant_calling.FilePathSplitSomatic} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathSplitSomatic, self).__init__(prefix=prefix)

        self.comparison_tsv = prefix + '.tsv'

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
        self.panel_of_normal_path = panel_of_normal_path  # Can be None.

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


class VariantCallingGATKTargetIntervals(object):
    """C{bsf.analyses.variant_calling.VariantCallingGATKTargetIntervals} class representing target intervals.

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar probes_path: Probes (baits) interval file path
    @type probes_path: str | unicode
    @ivar targets_path: Targets interval file path
    @type targets_path: str | unicode
    """

    @classmethod
    def from_sample(cls, sample):
        """Create a C{VariantCallingGATKTargetIntervals} object from a C{bsf.ngs.Sample} object.

        @param sample: C{bsf.ngs.Sample}
        @type sample: bsf.ngs.Sample
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
                    default_path=Default.absolute_intervals())

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
                    default_path=Default.absolute_intervals())

        return target_intervals

    def __init__(self, name=None, probes_path=None, targets_path=None):
        """Initialise a C{bsf.analyses.variant_calling.VariantCallingGATKTargetIntervals} object.

        @param name: Name
        @type name: str
        @param probes_path: Probes (baits) interval file path
        @type probes_path: str | unicode
        @param targets_path: Targets interval file path
        @type targets_path: str | unicode
        @return:
        @rtype:
        """
        self.name = name  # Can be None.
        self.probes_path = probes_path  # Can be None.
        self.targets_path = targets_path  # Can be None.

        return


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
    @cvar stage_name_annotate_cohort_snpeff: C{bsf.Stage.name} for the snpEff cohort annotation stage
    @type stage_name_annotate_cohort_snpeff: str
    @cvar stage_name_annotate_cohort_vep: C{bsf.Stage.name} for the Ensembl VEP cohort annotation stage
    @type stage_name_annotate_cohort_vep: str
    @cvar stage_name_split_cohort_snpeff: C{bsf.Stage.name} for the snpEff cohort splitting stage
    @type stage_name_split_cohort_snpeff: str
    @cvar stage_name_split_cohort_vep: C{bsf.Stage.name} for the Ensembl VEP cohort splitting stage
    @type stage_name_split_cohort_vep: str
    @cvar stage_name_summary: C{bsf.Stage.name} for the summary stage
    @type stage_name_summary: str
    @cvar stage_name_somatic: C{bsf.Stage.name} for the somatic stage
    @type stage_name_somatic: str
    @cvar stage_name_annotate_somatic_snpeff: C{bsf.Stage.name} for the snpEff somatic annotation stage
    @type stage_name_annotate_somatic_snpeff: str
    @cvar stage_name_annotate_somatic_vep: C{bsf.Stage.name} for the Ensembl VEP somatic annotation stage
    @type stage_name_annotate_somatic_vep: str
    @cvar stage_name_split_somatic_snpeff: C{bsf.Stage.name} for the snpEff somatic splitting stage
    @type stage_name_split_somatic_snpeff: str
    @cvar stage_name_split_somatic_vep: C{bsf.Stage.name} for the Ensembl VEP somatic splitting stage
    @type stage_name_split_somatic_vep: str
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
    @type vqsr_skip_indel: None | bool
    @ivar vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
    @type vqsr_skip_snp: None | bool
    @ivar vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
    @type vqsr_resources_indel_dict: dict[str, dict[str, str | unicode]]
    @ivar vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
    @type vqsr_resources_snp_dict: dict[str, dict[str, str | unicode]]
    @ivar vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_indel_list: list[str]
    @ivar vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_snp_list: list[str]
    @ivar vqsr_bad_lod_cutoff_indel: LOD score cutoff for negative training set for INDELs
    @type vqsr_bad_lod_cutoff_indel: None | float
    @ivar vqsr_bad_lod_cutoff_snp: LOD score cutoff for negative training set for SNPs
    @type vqsr_bad_lod_cutoff_snp: None | float
    @ivar vqsr_max_gaussians_pos_indel: Maximum number of Gaussians in the positive training for INDELs
    @type vqsr_max_gaussians_pos_indel: None | int
    @ivar vqsr_max_gaussians_pos_snp: Maximum number of Gaussians in the positive training for SNPs
    @type vqsr_max_gaussians_pos_snp: None | int
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
    @ivar gatk_bundle_version: GATK resource bundle version
    @type gatk_bundle_version: str
    @ivar snpeff_genome_version: snpEff genome version
    @type snpeff_genome_version: str
    @ivar genome_annotation_gtf: Genome annotation Gene Transfer Format (GTF) file path
    @type genome_annotation_gtf: str | unicode
    @ivar vep_assembly: Ensembl Variant Effect Predictor (VEP) assembly
    @type vep_assembly: None | str | unicode
    @ivar vep_cache: Ensembl Variant Effect Predictor (VEP) cache directory
    @type vep_cache: None | str | unicode
    @ivar vep_fasta: Ensembl Variant Effect Predictor (VEP) FASTA directory
    @type vep_fasta: None | str | unicode
    @ivar vep_plugin: Ensembl Variant Effect Predictor (VEP) plug-in directory
    @type vep_plugin: None | str | unicode
    @ivar vep_species: Ensembl Variant Effect Predictor (VEP) species
    @type vep_species: None | str | unicode
    @ivar vep_source: Ensembl Variant Effect Predictor (VEP) source directory
    @type vep_source: None | str | unicode
    @ivar vep_sql_user: Ensembl Variant Effect Predictor (VEP) SQL database user name
    @type vep_sql_user: None | str | unicode
    @ivar vep_sql_pass: Ensembl Variant Effect Predictor (VEP) SQL database password
    @type vep_sql_pass: None | str | unicode
    @ivar vep_sql_host: Ensembl Variant Effect Predictor (VEP) SQL host
    @type vep_sql_host: None | str | unicode
    @ivar vep_sql_port: Ensembl Variant Effect Predictor (VEP) SQL TCP/IP port
    @type vep_sql_port: None | str | unicode
    @ivar classpath_gatk: Genome Analysis Tool Kit Java Archive (JAR) class path directory
    @type classpath_gatk: None | str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: None | str | unicode
    @ivar classpath_snpeff: snpEff tool Java Archive (JAR) class path directory
    @type classpath_snpeff: None | str | unicode
    @ivar classpath_vcf_filter: VCF.Filter tool Java Archive (JAR) class path directory
    @type classpath_vcf_filter: None | str | unicode
    """

    name = 'Variant Calling Analysis'
    prefix = 'variant_calling'

    stage_name_align_lane = '_'.join((prefix, 'align_lane'))
    stage_name_process_lane = '_'.join((prefix, 'process_lane'))
    stage_name_process_sample = '_'.join((prefix, 'process_sample'))
    stage_name_diagnose_sample = '_'.join((prefix, 'diagnose_sample'))
    stage_name_merge_cohort = '_'.join((prefix, 'merge_cohort'))
    stage_name_process_cohort = '_'.join((prefix, 'process_cohort'))
    stage_name_annotate_cohort_snpeff = '_'.join((prefix, 'annotate_cohort_snpeff'))
    stage_name_annotate_cohort_vep = '_'.join((prefix, 'annotate_cohort_vep'))
    stage_name_split_cohort_snpeff = '_'.join((prefix, 'split_cohort_snpeff'))
    stage_name_split_cohort_vep = '_'.join((prefix, 'split_cohort_vep'))
    stage_name_summary = '_'.join((prefix, 'summary'))
    stage_name_somatic = '_'.join((prefix, 'somatic'))
    stage_name_annotate_somatic_snpeff = '_'.join((prefix, 'annotate_somatic_snpeff'))
    stage_name_annotate_somatic_vep = '_'.join((prefix, 'annotate_somatic_vep'))
    stage_name_split_somatic_snpeff = '_'.join((prefix, 'split_somatic_snpeff'))
    stage_name_split_somatic_vep = '_'.join((prefix, 'split_somatic_vep'))

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
            gatk_bundle_version=None,
            snpeff_genome_version=None,
            genome_annotation_gtf=None,
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
        @type vqsr_skip_indel: None | bool
        @param vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
        @type vqsr_skip_snp: None | bool
        @param vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_indel_dict: dict[str, dict[str, str | unicode]]
        @param vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_snp_dict: dict[str, dict[str, str | unicode]]
        @param vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_indel_list: list[str]
        @param vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_snp_list: list[str]
        @param vqsr_bad_lod_cutoff_indel: LOD score cutoff for negative training set for INDELs
        @type vqsr_bad_lod_cutoff_indel: None | float
        @param vqsr_bad_lod_cutoff_snp: LOD score cutoff for negative training set for SNPs
        @type vqsr_bad_lod_cutoff_snp: None | float
        @param vqsr_max_gaussians_pos_indel: Maximum number of Gaussians in the positive training for INDELs
        @type vqsr_max_gaussians_pos_indel: None | int
        @param vqsr_max_gaussians_pos_snp: Maximum number of Gaussians in the positive training for SNPs
        @type vqsr_max_gaussians_pos_snp: None | int
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
        @param gatk_bundle_version: GATK resource bundle version
        @type gatk_bundle_version: str
        @param snpeff_genome_version: snpEff genome version
        @type snpeff_genome_version: str
        @param genome_annotation_gtf: Genome annotation Gene Transfer Format (GTF) file path
        @type genome_annotation_gtf: str | unicode
        @param vep_assembly: Ensembl Variant Effect Predictor (VEP) assembly
        @type vep_assembly: None | str | unicode
        @param vep_fasta: Ensembl Variant Effect Predictor (VEP) FASTA directory
        @type vep_fasta: None | str | unicode
        @param vep_cache: Ensembl Variant Effect Predictor (VEP) cache directory
        @type vep_cache: None | str | unicode
        @param vep_plugin: Ensembl Variant Effect Predictor (VEP) plug-in directory
        @type vep_plugin: None | str | unicode
        @param vep_species: Ensembl Variant Effect Predictor (VEP) species
        @type vep_species: None | str | unicode
        @param vep_source: Ensembl Variant Effect Predictor (VEP) source directory
        @type vep_source: None | str | unicode
        @param vep_sql_user: Ensembl Variant Effect Predictor (VEP) SQL database user name
        @type vep_sql_user: None | str | unicode
        @param vep_sql_pass: Ensembl Variant Effect Predictor (VEP) SQL database password
        @type vep_sql_pass: None | str | unicode
        @param vep_sql_host: Ensembl Variant Effect Predictor (VEP) SQL host
        @type vep_sql_host: None | str | unicode
        @param vep_sql_port: Ensembl Variant Effect Predictor (VEP) SQL TCP/IP port
        @type vep_sql_port: None | str | unicode
        @param classpath_gatk: Genome Analysis Tool Kit Java Archive (JAR) class path directory
        @type classpath_gatk: None | str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: None | str | unicode
        @param classpath_snpeff: snpEff tool Java Archive (JAR) class path directory
        @type classpath_snpeff: None | str | unicode
        @param classpath_vcf_filter: VCF.Filter tool Java Archive (JAR) class path directory
        @type classpath_vcf_filter: None | str | unicode
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
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
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
            self.skip_mark_duplicates = skip_mark_duplicates

        if skip_indel_realignment is None:
            self.skip_indel_realignment = False
        else:
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
            self.vqsr_skip_indel = vqsr_skip_indel

        if vqsr_skip_snp is None:
            self.vqsr_skip_snp = False
        else:
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

        self.vqsr_bad_lod_cutoff_indel = vqsr_bad_lod_cutoff_indel
        self.vqsr_bad_lod_cutoff_snp = vqsr_bad_lod_cutoff_snp
        self.vqsr_max_gaussians_pos_indel = vqsr_max_gaussians_pos_indel
        self.vqsr_max_gaussians_pos_snp = vqsr_max_gaussians_pos_snp

        if exclude_intervals_list is None:
            self.exclude_intervals_list = list()
        else:
            self.exclude_intervals_list = exclude_intervals_list

        if include_intervals_list is None:
            self.include_intervals_list = list()
        else:
            self.include_intervals_list = include_intervals_list

        if interval_padding is None:
            self.interval_padding = 0
        else:
            self.interval_padding = interval_padding

        if number_of_tiles_cohort is None:
            self.number_of_tiles_cohort = 0
        else:
            self.number_of_tiles_cohort = number_of_tiles_cohort

        if number_of_chunks_cohort is None:
            self.number_of_chunks_cohort = 0
        else:
            self.number_of_chunks_cohort = number_of_chunks_cohort

        if number_of_tiles_somatic is None:
            self.number_of_tiles_somatic = 0
        else:
            self.number_of_tiles_somatic = number_of_tiles_somatic

        if number_of_chunks_somatic is None:
            self.number_of_chunks_somatic = 0
        else:
            self.number_of_chunks_somatic = number_of_chunks_somatic

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

        self.classpath_gatk = classpath_gatk
        self.classpath_picard = classpath_picard
        self.classpath_snpeff = classpath_snpeff
        self.classpath_vcf_filter = classpath_vcf_filter

        self._comparison_dict = dict()
        """ @type _comparison_dict: dict[str, bsf.analyses.variant_calling.VariantCallingGATKComparison] """

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
                                'Missing configuration option {!r} in configuration section {!r}.'.format(
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
                                'Missing configuration option {!r} in configuration section {!r}.'.format(
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

        option = 'number_of_tiles_somatic'
        if config_parser.has_option(section=section, option=option):
            self.number_of_tiles_somatic = config_parser.getint(section=section, option=option)

        # Get the number of chunks for variant_calling_somatic.

        option = 'number_of_chunks_somatic'
        if config_parser.has_option(section=section, option=option):
            self.number_of_chunks_somatic = config_parser.getint(section=section, option=option)

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

        # VEP names

        option = 'vep_assembly'
        if config_parser.has_option(section=section, option=option):
            self.vep_assembly = config_parser.get(section=section, option=option)

        option = 'vep_species'
        if config_parser.has_option(section=section, option=option):
            self.vep_species = config_parser.get(section=section, option=option)

        # VEP directories

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

        # VEP SQL database

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

    def run(self):
        """Run a C{bsf.analyses.variant_calling.VariantCallingGATK} analysis.

        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file from disk.

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
            @return:
            @rtype:
            """
            # For variant calling, all samples need adding to the Analysis regardless.
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
                        print 'Comparison sheet row_dict:', row_dict

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
                    if comparison.tumor_sample is not None:
                        self._comparison_dict[comparison.get_name] = comparison

            return

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
                'VEP_ExAC_AF',
                'VEP_ExAC_AFR_AF',
                'VEP_ExAC_AMR_AF',
                'VEP_ExAC_Adj_AF',
                'VEP_ExAC_EAS_AF',
                'VEP_ExAC_FIN_AF',
                'VEP_ExAC_NFE_AF',
                'VEP_ExAC_OTH_AF',
                'VEP_ExAC_SAS_AF',
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
                'VQSLOD',
                'culprit',
            ),
        }

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
                raise Exception('Picard sequence dictionary {!r} does not exist.'.format(dict_path))

            alignment_file = pysam.AlignmentFile(dict_path, 'r')
            # Summarise sequence lengths to get the total length.
            for sq_entry in alignment_file.header['SQ']:
                """ @type sq_entry: dict """
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
                C{Runnable}, Python C{str} of Python C{unicode} object value data
            @type cohort_runnable_dict: dict[str, list[Runnable | str | unicode]]
            @param cohort_name: Cohort name to select a Python list of C{Runnable} objects from the
                I{cohort_runnable_dict} Python C{dict}
            @type cohort_name: str
            @return: Final C{Runnable} of the gather stage
            @rtype: Runnable
            """

            # Private variables are prefixed with an underscore to avoid clashes with variables in the run() method.

            prefix_merge_cohort_final = '_'.join((analysis_stage.name, cohort_name))

            file_path_merge_cohort_final = FilePathMergeCohort(prefix=prefix_merge_cohort_final)

            # If the final TBI index file already exists, create the Runnable objects, but do not submit their
            # corresponding Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_merge_cohort_final.combined_gvcf_tbi)):
                final_index_exists = False
            else:
                final_index_exists = True

            # The cohort_object_list contains either Runnable objects from the process_sample stage or
            # Python str | unicode (GVCF file path) objects for accessory cohorts to be merged.
            cohort_object_list = cohort_runnable_dict[cohort_name]

            # Scatter
            runnable_scatter = None
            """ @type runnable_scatter: Runnable """
            runnable_scatter_list = list()
            """ @type runnable_scatter_list: list[Runnable] """

            tile_index_list = range(0, len(self._tile_region_cohort_list))

            for tile_index in tile_index_list:
                prefix_merge_cohort_scatter = '_'.join((analysis_stage.name, cohort_name, 'scatter', str(tile_index)))

                file_path_merge_cohort_scatter = FilePathMergeCohort(prefix=prefix_merge_cohort_scatter)

                runnable_scatter = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_merge_cohort_scatter,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        file_path_object=file_path_merge_cohort_scatter,
                        debug=self.debug))
                executable_scatter = self.set_stage_runnable(
                    stage=analysis_stage,
                    runnable=runnable_scatter)
                # Submit the Executable only, if the final TBI index file does not exist,
                # but do not override the state set by the Analysis.set_stage_runnable() method.
                if not final_index_exists:
                    executable_scatter.submit = False
                for cohort_component in cohort_object_list:
                    # Set dependencies on preceding Runnable.name or Executable.name objects.
                    # Set them only for Runnable objects, but not for Python str | unicode (file path) objects.
                    if isinstance(cohort_component, Runnable):
                        executable_scatter.dependencies.append(cohort_component.name)

                runnable_scatter_list.append(runnable_scatter)

                reference_scatter = runnable_scatter.get_absolute_cache_file_path(file_path=self.bwa_genome_db)

                # Run GATK CombineGVCFs

                _runnable_step = runnable_scatter.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='merge_cohort_gatk_combine_gvcfs',
                        java_temporary_path=runnable_scatter.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        gatk_classpath=self.classpath_gatk))
                """ @type _runnable_step: RunnableStepGATK """
                _runnable_step.add_gatk_option(key='analysis_type', value='CombineGVCFs')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_scatter)
                for _interval in self.exclude_intervals_list:
                    _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                for region_tuple in self._tile_region_cohort_list[tile_index]:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if region_tuple[0]:
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value='{:s}:{:d}-{:d}'.format(region_tuple[0], region_tuple[1], region_tuple[2]),
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
                    value=file_path_merge_cohort_scatter.combined_gvcf_vcf)

            # Gather
            runnable_gather = None
            """ @type runnable_gather: Runnable """

            if len(self._tile_region_cohort_list) == 1:
                # If there is only one tile, no need to gather.
                # Assign the sole scatter Runnable to the sole gather Runnable.
                runnable_gather = runnable_scatter
            else:
                # Gather by hierarchically merging by the number of chunks on the partitioned genome tile index list.
                # Initialise a list of Runnable objects and indices for the hierarchical merge.
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

                    for partition_index in range(0, len(partition_list)):
                        chunk_index_list = partition_list[partition_index]
                        # The file prefix includes the level and partition index.
                        prefix_merge_cohort_gather = '_'.join((
                            analysis_stage.name,
                            cohort_name,
                            'gather',
                            str(gather_level),
                            str(partition_index)))

                        file_path_merge_cohort_gather = FilePathMergeCohort(prefix=prefix_merge_cohort_gather)

                        runnable_gather = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_merge_cohort_gather,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                file_path_object=file_path_merge_cohort_gather,
                                debug=self.debug))
                        executable_gather = self.set_stage_runnable(
                            stage=analysis_stage,
                            runnable=runnable_gather)
                        # Submit the Executable only, if the final TBI index file does not exist,
                        # but do not override the state set by the Analysis.set_stage_runnable() method.
                        if not final_index_exists:
                            executable_gather.submit = False
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_runnable_gather_list.append(runnable_gather)
                        temporary_tile_index_list.append(partition_index)

                        reference_gather = runnable_gather.get_absolute_cache_file_path(file_path=self.bwa_genome_db)

                        # GATK CatVariants bypasses the GATK engine and thus requires a completely different
                        # command line.
                        _runnable_step = runnable_gather.add_runnable_step(
                            runnable_step=RunnableStepJava(
                                name='merge_cohort_gatk_cat_variants',
                                sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                                java_temporary_path=runnable_gather.get_relative_temporary_directory_path,
                                java_heap_maximum='Xmx4G'))
                        """ @type _runnable_step: RunnableStepJava """
                        _runnable_step.add_option_short(
                            key='classpath',
                            value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                        _sub_command = _runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        _sub_command.add_option_long(key='reference', value=reference_gather)
                        _sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_merge_cohort_gather.combined_gvcf_vcf)
                        _sub_command.add_switch_long(key='assumeSorted')
                        # Finally, process per chunk index.
                        for chunk_index in chunk_index_list:
                            _runnable_object = runnable_gather_list[chunk_index]
                            _file_path_object = _runnable_object.file_path_object
                            """ @type _file_path_object: FilePathMergeCohort """
                            # Set GATK option variant
                            _sub_command.add_option_long(
                                key='variant',
                                value=_file_path_object.combined_gvcf_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            _runnable_step.obsolete_file_path_list.append(_file_path_object.combined_gvcf_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            _runnable_step.obsolete_file_path_list.append(_file_path_object.combined_gvcf_tbi)
                            # Set dependencies on preceding Runnable.name or Executable.name objects.
                            # Depend on the Runnable.name of the corresponding Runnable of the scattering above.
                            executable_gather.dependencies.append(_runnable_object.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    runnable_gather_list = temporary_runnable_gather_list
                    tile_index_list = temporary_tile_index_list
                    gather_level += 1

            # For the last gather Runnable, move (rename) file paths to the top-level prefix and
            # adjust the FilePath object accordingly.
            file_path_merge_cohort_gather = runnable_gather.file_path_object
            """ @type file_path_merge_cohort_gather: FilePathMergeCohort """

            runnable_gather.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='merge_cohort_gather_move_vcf',
                    source_path=file_path_merge_cohort_gather.combined_gvcf_vcf,
                    target_path=file_path_merge_cohort_final.combined_gvcf_vcf))

            runnable_gather.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='merge_cohort_gather_move_tbi',
                    source_path=file_path_merge_cohort_gather.combined_gvcf_tbi,
                    target_path=file_path_merge_cohort_final.combined_gvcf_tbi))

            file_path_merge_cohort_gather.combined_gvcf_vcf = file_path_merge_cohort_final.combined_gvcf_vcf
            file_path_merge_cohort_gather.combined_gvcf_tbi = file_path_merge_cohort_final.combined_gvcf_tbi

            return runnable_gather

        def run_genotype_cohort_scatter_gather():
            """Private function to genotype a cohort in a scatter and gather approach.

            @return: Final C{Runnable} of the gather stage
            @rtype: Runnable
            """

            prefix_process_cohort_final = '_'.join((stage_process_cohort.name, self.cohort_name))

            file_path_genotype_cohort_final = FilePathGenotypeCohort(prefix=prefix_process_cohort_final)

            # If the final TBI index file already exists, create the Runnable objects, but do not submit their
            # corresponding Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_genotype_cohort_final.genotyped_raw_tbi)):
                final_index_exists = False
            else:
                final_index_exists = True

            runnable_scatter = None
            """ @type runnable_scatter: Runnable """
            runnable_scatter_list = list()
            """ @type runnable_scatter_list: list[Runnable] """

            tile_index_list = range(0, len(self._tile_region_cohort_list))

            for tile_index in tile_index_list:
                prefix_process_cohort_scatter = '_'.join((
                    stage_process_cohort.name, self.cohort_name, 'scatter', str(tile_index)))

                file_path_genotype_cohort_scatter = FilePathGenotypeCohort(prefix=prefix_process_cohort_scatter)

                runnable_scatter = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_process_cohort_scatter,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        file_path_object=file_path_genotype_cohort_scatter,
                        debug=self.debug))
                executable_scatter = self.set_stage_runnable(
                    stage=stage_process_cohort,
                    runnable=runnable_scatter)
                # Submit the Executable only, if the final TBI index file does not exist,
                # but do not override the state set by the Analysis.set_stage_runnable() method.
                if not final_index_exists:
                    executable_scatter.submit = False
                # Set dependencies on preceding Runnable.name or Executable.name objects.
                executable_scatter.dependencies.append(runnable_merge_cohort.name)

                runnable_scatter_list.append(runnable_scatter)

                reference_scatter = runnable_scatter.get_absolute_cache_file_path(file_path=self.bwa_genome_db)

                # Run the GATK GenotypeGVCFs analysis.

                _runnable_step = runnable_scatter.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_cohort_gatk_genotype_gvcfs_scatter',
                        java_temporary_path=runnable_scatter.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx12G',
                        gatk_classpath=self.classpath_gatk))
                """ @type _runnable_step: RunnableStepGATK """
                _runnable_step.add_gatk_option(key='analysis_type', value='GenotypeGVCFs')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_scatter)
                for _interval in self.exclude_intervals_list:
                    _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                for region_tuple in self._tile_region_cohort_list[tile_index]:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if region_tuple[0]:
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value='{:s}:{:d}-{:d}'.format(region_tuple[0], region_tuple[1], region_tuple[2]),
                            override=True)
                if self.known_sites_discovery:
                    _runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
                _runnable_step.add_gatk_option(key='variant', value=file_path_merge_cohort.combined_gvcf_vcf)
                _runnable_step.add_gatk_option(key='out', value=file_path_genotype_cohort_scatter.genotyped_raw_vcf)

            # Gather
            runnable_gather = None
            """ @type runnable_gather: Runnable """

            if len(self._tile_region_cohort_list) == 1:
                # If there is only one tile, no need to gather.
                # Assign the sole scatter Runnable to the sole gather Runnable.
                runnable_gather = runnable_scatter
            else:
                # Gather by hierarchically merging by the number of chunks on the partitioned genome tile index list.
                # Initialise a list of Runnable objects and indices for the hierarchical merge.
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

                    for partition_index in range(0, len(partition_list)):
                        chunk_index_list = partition_list[partition_index]
                        # The file prefix includes the level and partition index.
                        prefix_process_cohort_gather = '_'.join(
                            (stage_process_cohort.name,
                             self.cohort_name,
                             'gather',
                             str(gather_level),
                             str(partition_index)))

                        file_path_genotype_cohort_gather = FilePathGenotypeCohort(prefix=prefix_process_cohort_gather)

                        runnable_gather = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_process_cohort_gather,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                file_path_object=file_path_genotype_cohort_gather,
                                debug=self.debug))
                        executable_gather = self.set_stage_runnable(
                            stage=stage_process_cohort,
                            runnable=runnable_gather)
                        # Submit the Executable only, if the final TBI index file does not exist,
                        # but do not override the state set by the Analysis.set_stage_runnable() method.
                        if not final_index_exists:
                            executable_gather.submit = False
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_runnable_gather_list.append(runnable_gather)
                        temporary_tile_index_list.append(partition_index)

                        reference_gather = runnable_gather.get_absolute_cache_file_path(file_path=self.bwa_genome_db)

                        # GATK CatVariants by-passes the GATK engine and thus requires a completely different
                        # command line.
                        _runnable_step = runnable_gather.add_runnable_step(
                            runnable_step=RunnableStepJava(
                                name='merge_cohort_gatk_cat_variants',
                                sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                                java_temporary_path=runnable_gather.get_relative_temporary_directory_path,
                                java_heap_maximum='Xmx4G'))
                        """ @type _runnable_step: RunnableStepJava """
                        _runnable_step.add_option_short(
                            key='classpath',
                            value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                        _sub_command = _runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        _sub_command.add_option_long(key='reference', value=reference_gather)
                        _sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_genotype_cohort_gather.genotyped_raw_vcf)
                        _sub_command.add_switch_long(key='assumeSorted')
                        # Finally, add RunnableStep options, obsolete files and Executable dependencies per chunk index.
                        for chunk_index in chunk_index_list:
                            runnable_object = runnable_gather_list[chunk_index]
                            file_path_object = runnable_object.file_path_object
                            """ @type file_path_object: FilePathGenotypeCohort """
                            # Set GATK option variant
                            _sub_command.add_option_long(
                                key='variant',
                                value=file_path_object.genotyped_raw_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            _runnable_step.obsolete_file_path_list.append(file_path_object.genotyped_raw_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            _runnable_step.obsolete_file_path_list.append(file_path_object.genotyped_raw_tbi)
                            # Set dependencies on preceding Runnable.name or Executable.name objects.
                            executable_gather.dependencies.append(runnable_object.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    runnable_gather_list = temporary_runnable_gather_list
                    tile_index_list = temporary_tile_index_list
                    gather_level += 1

            # For the last gather Runnable, move (rename) file paths to the top-level prefix and
            # adjust the FilePath object accordingly.
            file_path_genotype_cohort_gather = runnable_gather.file_path_object
            """ @type file_path_genotype_cohort_gather: FilePathGenotypeCohort """

            runnable_gather.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='process_cohort_gather_move_vcf',
                    source_path=file_path_genotype_cohort_gather.genotyped_raw_vcf,
                    target_path=file_path_genotype_cohort_final.genotyped_raw_vcf))

            runnable_gather.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='process_cohort_gather_move_tbi',
                    source_path=file_path_genotype_cohort_gather.genotyped_raw_tbi,
                    target_path=file_path_genotype_cohort_final.genotyped_raw_tbi))

            file_path_genotype_cohort_gather.genotyped_raw_vcf = file_path_genotype_cohort_final.genotyped_raw_vcf
            file_path_genotype_cohort_gather.genotyped_raw_tbi = file_path_genotype_cohort_final.genotyped_raw_tbi

            return runnable_gather

        def run_somatic_scatter_gather(comparison_key):
            """Private method to run somatic variant calling in a scatter and gather approach.

            @param comparison_key: C{VariantCallingGATKComparison.name}
            @type comparison_key: str
            @return: Final C{Runnable} of the gather stage
            @rtype: Runnable
            """

            # Private variables are prefixed with an underscore to avoid clashes with variables in the run() method.

            prefix_somatic_final = '_'.join((stage_somatic.name, comparison_key))

            file_path_somatic_final = FilePathSomatic(prefix=prefix_somatic_final)

            # If the final TBI index file already exists, create the Runnable objects, but do not submit their
            # corresponding Executable objects.

            if os.path.exists(os.path.join(self.genome_directory, file_path_somatic_final.somatic_tbi)):
                final_index_exists = False
            else:
                final_index_exists = True

            comparison = self._comparison_dict[comparison_key]

            _target_intervals = VariantCallingGATKTargetIntervals.from_sample(sample=comparison.tumor_sample)

            # Scatter
            runnable_scatter = None
            """ @type runnable_scatter: Runnable """
            runnable_scatter_list = list()
            """ @type runnable_scatter_list: list[Runnable] """

            tile_index_list = range(0, len(self._tile_region_somatic_list))

            for tile_index in tile_index_list:
                prefix_somatic_scatter = '_'.join((stage_somatic.name, comparison_key, 'scatter', str(tile_index)))

                file_path_somatic_scatter = FilePathSomaticScatterGather(prefix=prefix_somatic_scatter)

                runnable_scatter = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_somatic_scatter,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        cache_path_dict=self._cache_path_dict,
                        file_path_object=file_path_somatic_scatter,
                        debug=self.debug))
                executable_scatter = self.set_stage_runnable(
                    stage=stage_somatic,
                    runnable=runnable_scatter)
                # Submit the Executable only, if the final TBI index file does not exist,
                # but do not override the state set by the Analysis.set_stage_runnable() method.
                if not final_index_exists:
                    executable_scatter.submit = False
                # Set dependencies on preceding Runnable.name or Executable.name objects.
                if comparison.normal_sample is not None:
                    executable_scatter.dependencies.append(
                        '_'.join((stage_process_sample.name, comparison.normal_sample.name)))
                if comparison.tumor_sample is not None:
                    executable_scatter.dependencies.append(
                        '_'.join((stage_process_sample.name, comparison.tumor_sample.name)))

                runnable_scatter_list.append(runnable_scatter)

                reference_scatter = runnable_scatter.get_absolute_cache_file_path(file_path=self.bwa_genome_db)

                # Run GATK MuTect2

                _runnable_step = runnable_scatter.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='somatic_gatk_mutect2_scatter',
                        java_temporary_path=runnable_scatter.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        gatk_classpath=self.classpath_gatk))
                """ @type _runnable_step: RunnableStepGATK """
                _runnable_step.add_gatk_option(key='analysis_type', value='MuTect2')
                _runnable_step.add_gatk_option(key='reference_sequence', value=reference_scatter)
                for _interval in self.exclude_intervals_list:
                    _runnable_step.add_gatk_option(key='excludeIntervals', value=_interval, override=True)
                for region_tuple in self._tile_region_somatic_list[tile_index]:
                    # The list of tiles is initialised to an empty tile to trigger at least one process.
                    # Do not assign an interval in such cases.
                    if region_tuple[0]:
                        _runnable_step.add_gatk_option(
                            key='intervals',
                            value='{:s}:{:d}-{:d}'.format(region_tuple[0], region_tuple[1], region_tuple[2]),
                            override=True)
                    elif _target_intervals.targets_path:
                        # If not running on genome tiles, the MuTect2 analysis is run on the target intervals, only.
                        _runnable_step.add_gatk_option(key='intervals', value=_target_intervals.targets_path)
                        _runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                if self.known_sites_discovery:
                    _runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
                for _file_path in self.known_somatic_discovery:
                    _runnable_step.add_gatk_option(key='cosmic', value=_file_path, override=True)

                # Find and add the FilePathProcessSample object for the 'normal' Sample object.
                if comparison.normal_sample is not None:
                    file_path_object = self.runnable_dict['_'.join((
                        stage_process_sample.name,
                        comparison.normal_sample.name))].file_path_object
                    """ @type file_path_object: FilePathProcessSample """
                    _runnable_step.add_gatk_option(key='input_file:normal', value=file_path_object.realigned_bam)
                elif comparison.panel_of_normal_path is not None:
                    _runnable_step.add_gatk_option(key='normal_panel', value=comparison.panel_of_normal_path)

                # Find and add the FilePathProcessSample object for the 'tumor' Sample object.
                if comparison.tumor_sample is not None:
                    file_path_object = self.runnable_dict['_'.join((
                        stage_process_sample.name,
                        comparison.tumor_sample.name))].file_path_object
                    """ @type file_path_object: FilePathProcessSample """
                    _runnable_step.add_gatk_option(key='input_file:tumor', value=file_path_object.realigned_bam)

                _runnable_step.add_gatk_option(key='out', value=file_path_somatic_scatter.somatic_vcf)

            # Gather
            runnable_gather = None
            """ @type runnable_gather: Runnable """

            if len(self._tile_region_somatic_list) == 1:
                # If there is only one tile, no need to gather.
                # Assign the sole scatter Runnable to the sole gather Runnable.
                runnable_gather = runnable_scatter
            else:
                # Gather by hierarchically merging by the number of chunks on the partitioned genome tile index list.
                # Initialise a list of Runnable objects and indices for the hierarchical merge.
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

                    for partition_index in range(0, len(partition_list)):
                        chunk_index_list = partition_list[partition_index]
                        # The file prefix includes the level and partition index.
                        prefix_somatic_gather = '_'.join((
                            stage_somatic.name,
                            comparison_key,
                            'gather',
                            str(gather_level),
                            str(partition_index)))

                        file_path_somatic_gather = FilePathSomaticScatterGather(prefix=prefix_somatic_gather)

                        runnable_gather = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_somatic_gather,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                cache_directory=self.cache_directory,
                                cache_path_dict=self._cache_path_dict,
                                file_path_object=file_path_somatic_gather,
                                debug=self.debug))
                        executable_gather = self.set_stage_runnable(
                            stage=stage_somatic,
                            runnable=runnable_gather)
                        # Submit the Executable only, if the final TBI index file does not exist,
                        # but do not override the state set by the Analysis.set_stage_runnable() method.
                        if not final_index_exists:
                            executable_gather.submit = False
                        # Dependencies on scatter processes are set based on genome tile indices below.
                        temporary_runnable_gather_list.append(runnable_gather)
                        temporary_tile_index_list.append(partition_index)

                        reference_gather = runnable_gather.get_absolute_cache_file_path(file_path=self.bwa_genome_db)

                        # GATK CatVariants by-passes the GATK engine and thus requires a completely different
                        # command line.
                        _runnable_step = runnable_gather.add_runnable_step(
                            runnable_step=RunnableStepJava(
                                name='somatic_gatk_cat_variants',
                                sub_command=Command(program='org.broadinstitute.gatk.tools.CatVariants'),
                                java_temporary_path=runnable_gather.get_relative_temporary_directory_path,
                                java_heap_maximum='Xmx4G'))
                        """ @type _runnable_step: RunnableStepJava """
                        _runnable_step.add_option_short(
                            key='classpath',
                            value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                        _sub_command = _runnable_step.sub_command
                        # Add the 'reference' not 'reference_sequence' option.
                        _sub_command.add_option_long(key='reference', value=reference_gather)
                        _sub_command.add_option_long(
                            key='outputFile',
                            value=file_path_somatic_gather.somatic_vcf)
                        _sub_command.add_switch_long(key='assumeSorted')
                        # Finally, add RunnableStep options, obsolete files and Executable dependencies per chunk index.
                        for chunk_index in chunk_index_list:
                            runnable_object = runnable_gather_list[chunk_index]
                            file_path_object = runnable_object.file_path_object
                            """ @type file_path_object: FilePathSomaticScatterGather """
                            # Set GATK option variant
                            _sub_command.add_option_long(
                                key='variant',
                                value=file_path_object.somatic_vcf,
                                override=True)
                            # Delete the *.g.vcf.gz file.
                            _runnable_step.obsolete_file_path_list.append(file_path_object.somatic_vcf)
                            # Delete the *.g.vcf.gz.tbi file.
                            _runnable_step.obsolete_file_path_list.append(file_path_object.somatic_tbi)
                            # Set dependencies on preceding Runnable.name or Executable.name objects.
                            executable_gather.dependencies.append(runnable_object.name)

                    # Set the temporary index list as the new list and increment the merge level.
                    runnable_gather_list = temporary_runnable_gather_list
                    tile_index_list = temporary_tile_index_list
                    gather_level += 1

            # For the last gather Runnable, move (rename) file paths to the top-level prefix and
            # adjust the FilePath object accordingly.
            file_path_somatic_gather = runnable_gather.file_path_object
            """ @type file_path_somatic_gather: FilePathSomaticScatterGather """

            runnable_gather.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='somatic_gather_move_vcf',
                    source_path=file_path_somatic_gather.somatic_vcf,
                    target_path=file_path_somatic_final.somatic_vcf))

            runnable_gather.add_runnable_step(
                runnable_step=RunnableStepMove(
                    name='somatic_gather_move_tbi',
                    source_path=file_path_somatic_gather.somatic_tbi,
                    target_path=file_path_somatic_final.somatic_tbi))

            file_path_somatic_gather.somatic_vcf = file_path_somatic_final.somatic_vcf
            file_path_somatic_gather.somatic_tbi = file_path_somatic_final.somatic_tbi

            return runnable_gather

        def run_annotate_snpeff(prefix, vcf_file_path):
            """Private function to annotate a VCF file via the snpEff tool.

            This function is used in both, cohort and somatic variant calling annotation.
            @param prefix: Prefix
            @type prefix: str
            @param vcf_file_path: VCF file path
            @type vcf_file_path: str | unicode
            @return: C{Runnable}
            @rtype: bsf.Runnable
            """
            # snpEff                  (snpeff)
            # Bgzip                   (snpeff_bgzip)
            # Tabix                   (snpeff_tabix)
            # GATK VariantAnnotator   (gatk_variant_annotator)

            prefix_annotate = prefix

            file_path_annotate = FilePathAnnotateSnpEff(prefix=prefix_annotate)

            runnable_annotate = self.add_runnable(
                runnable=Runnable(
                    name=prefix_annotate,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_annotate,
                    debug=self.debug))

            reference_annotate = runnable_annotate.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the snpEff tool for functional variant annotation.

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='snpeff',
                    program='java',
                    sub_command=Command(program='eff')))
            """ @type _runnable_step: RunnableStep """

            _runnable_step.add_switch_short(
                key='d64')
            _runnable_step.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_snpeff, 'snpEff.jar'))
            _runnable_step.add_switch_short(
                key='Xmx6G')
            _runnable_step.add_option_pair(
                key='-Djava.io.tmpdir',
                value=runnable_annotate.get_relative_temporary_directory_path)
            _runnable_step.stdout_path = file_path_annotate.snpeff_vcf

            _sub_command = _runnable_step.sub_command
            _sub_command.add_switch_short(key='download')
            _sub_command.add_option_short(key='o', value='gatk')
            _sub_command.add_option_short(key='stats', value=file_path_annotate.snpeff_stats)
            _sub_command.add_option_short(key='config', value=os.path.join(self.classpath_snpeff, 'snpEff.config'))

            _sub_command.arguments.append(self.snpeff_genome_version)
            _sub_command.arguments.append(vcf_file_path)

            # Automatically compress and index the snpEff VCF file with bgzip and tabix, respectively.

            runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='snpeff_bgzip',
                    program='bgzip',
                    arguments=[file_path_annotate.snpeff_vcf]))

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='snpeff_tabix',
                    program='tabix',
                    arguments=[file_path_annotate.snpeff_vcf_bgz]))
            """ @type _runnable_step: RunnableStep """
            _runnable_step.add_option_long(key='preset', value='vcf')

            # Run the GATK VariantAnnotator analysis.

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='gatk_variant_annotator',
                    java_temporary_path=runnable_annotate.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            """ @type _runnable_step: RunnableStepGATK """
            _runnable_step.add_gatk_option(key='analysis_type', value='VariantAnnotator')
            _runnable_step.add_gatk_option(key='reference_sequence', value=reference_annotate)
            if self.known_sites_discovery:
                _runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)

            # Add annotation resources and their corresponding expression options.
            for _annotation_resource in self.annotation_resources_dict.keys():
                if len(self.annotation_resources_dict[_annotation_resource][0]) \
                        and len(self.annotation_resources_dict[_annotation_resource][1]):
                    _runnable_step.add_gatk_option(
                        key=':'.join(('resource', _annotation_resource)),
                        value=self.annotation_resources_dict[_annotation_resource][0])
                    for _annotation in self.annotation_resources_dict[_annotation_resource][1]:
                        _runnable_step.add_gatk_option(
                            key='expression',
                            value='.'.join((_annotation_resource, _annotation)),
                            override=True)

            _runnable_step.add_gatk_option(key='variant', value=vcf_file_path)
            # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
            # _runnable_step.add_gatk_option(key='annotation', value='AlleleBalanceBySample')
            _runnable_step.add_gatk_option(key='annotation', value='SnpEff')
            _runnable_step.add_gatk_option(key='snpEffFile', value=file_path_annotate.snpeff_vcf_bgz)
            _runnable_step.add_gatk_option(key='out', value=file_path_annotate.annotated_vcf)

            return runnable_annotate

        def run_annotate_vep(prefix, vcf_file_path):
            """Private function to annotate a VCF file via the Ensembl Variant Effect Predictor (VEP).

            This function is used in both, cohort and somatic variant calling annotation.
            @param prefix: Prefix
            @type prefix: str
            @param vcf_file_path: VCF file path
            @type vcf_file_path: str | unicode
            @return: C{Runnable}
            @rtype: bsf.Runnable
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

            runnable_annotate = self.add_runnable(
                runnable=Runnable(
                    name=prefix_annotate,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_annotate,
                    debug=self.debug))

            # reference_annotate = runnable_annotate.get_absolute_cache_file_path(
            #     file_path=self.bwa_genome_db)

            # if not os.path.exists(os.path.join(self.genome_directory, file_path_annotate.vep_complete_vcf_tbi)):
            # Run the Ensembl Variant Effect Predictor script.

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='ensembl_vep',
                    program='perl',
                    sub_command=Command()))
            """ @type _runnable_step: RunnableStep """
            # self.set_runnable_step_configuration(runnable_step=_runnable_step)
            _runnable_step.arguments.append(os.path.join(self.vep_source, 'vep'))
            _sub_command = _runnable_step.sub_command
            # Basic options
            _sub_command.add_switch_long(key='everything')
            # Input options
            _sub_command.add_option_long(key='species', value=self.vep_species)
            _sub_command.add_option_long(key='assembly', value=self.vep_assembly)
            _sub_command.add_option_long(key='input_file', value=vcf_file_path)
            _sub_command.add_option_long(key='format', value='vcf')  # Input file format
            _sub_command.add_option_long(key='output_file', value=file_path_annotate.vep_complete_raw_vcf)
            _sub_command.add_switch_long(key='force_overwrite')
            _sub_command.add_option_long(key='stats_file', value=file_path_annotate.vep_statistics)
            # Cache options
            _sub_command.add_switch_long(key='cache')
            if False:
                # FIXME: For Ensembl VEP 91.
                _sub_command.add_switch_long(key='offline')
            _sub_command.add_option_long(key='dir_cache', value=self.vep_cache)
            _sub_command.add_option_long(key='dir_plugins', value=self.vep_plugin)
            if False:
                # FIXME: For Ensembl VEP 91.
                _sub_command.add_option_long(key='fasta_dir', value=self.vep_fasta)
            # Other annotation sources
            _sub_command.add_option_long(  # TODO: Has to be configurable
                key='plugin',
                value='CADD,/scratch/lab_bsf/resources/CADD/b37/whole_genome_SNVs.tsv.gz')
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
                value=runnable_annotate.get_relative_temporary_directory_path)

            runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='ensembl_vep_bgzip',
                    program='bgzip',
                    arguments=[file_path_annotate.vep_complete_raw_vcf]))

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='ensembl_vep_tabix',
                    program='tabix',
                    arguments=[file_path_annotate.vep_complete_raw_vcf_bgz]))
            """ @type _runnable_step: RunnableStep """
            _runnable_step.add_option_long(key='preset', value='vcf')

            # if not os.path.exists(os.path.join(self.genome_directory, file_path_annotate.vep_filtered_vcf_tbi)):
            # Run the Ensembl Variant Effect Filter script.

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='ensembl_filter',
                    program='perl',
                    sub_command=Command()))
            """ @type _runnable_step: RunnableStep """
            # self.set_runnable_step_configuration(runnable_step=_runnable_step)
            _runnable_step.arguments.append(os.path.join(self.vep_source, 'filter_vep'))
            _sub_command = _runnable_step.sub_command
            _sub_command.add_option_long(key='input_file', value=file_path_annotate.vep_complete_raw_vcf_bgz)
            _sub_command.add_option_long(key='format', value='vcf')
            _sub_command.add_option_long(key='output_file', value=file_path_annotate.vep_filtered_raw_vcf)
            _sub_command.add_switch_long(key='only_matched')
            _sub_command.add_option_long(key='filter', value='Consequence ne upstream_gene_variant', override=True)
            _sub_command.add_option_long(key='filter', value='Consequence ne downstream_gene_variant', override=True)
            _sub_command.add_option_long(key='filter', value='Consequence ne intron_variant', override=True)
            _sub_command.add_option_long(key='filter', value='BIOTYPE ne processed_transcript', override=True)
            # _sub_command.add_option_long(key='filter', value='CANONICAL eq YES', override=True)
            _sub_command.add_switch_long(key='force_overwrite')

            runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='ensembl_filter_bgzip',
                    program='bgzip',
                    arguments=[file_path_annotate.vep_filtered_raw_vcf]))

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStep(
                    name='ensembl_filter_tabix',
                    program='tabix',
                    arguments=[file_path_annotate.vep_filtered_raw_vcf_bgz]))
            """ @type _runnable_step: RunnableStep """
            _runnable_step.add_option_long(key='preset', value='vcf')

            # Run the VCF Filter on the complete VEP set to convert (re-model) the CSQ field into
            # a set of independent INFO fields.

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStepJava(
                    name='vcf_filter_complete',
                    sub_command=Command(program='at.ac.oeaw.cemm.bsf.vcffilter.vep.vep2vcf'),
                    java_temporary_path=runnable_annotate.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G'))
            """ @type _runnable_step: RunnableStepJava """
            _runnable_step.add_option_short(
                key='classpath',
                value=os.path.join(self.classpath_vcf_filter, 'VCFFilter.jar'))
            _sub_command = _runnable_step.sub_command
            _sub_command.add_option_pair(key='INPUT', value=file_path_annotate.vep_complete_raw_vcf_bgz)
            _sub_command.add_option_pair(key='OUTPUT', value=file_path_annotate.vep_complete_vcf_bgz)

            # Run the VCF Filter on the filtered VEP set to convert (re-model) the CSQ field into
            # a set of independent INFO fields.

            _runnable_step = runnable_annotate.add_runnable_step(
                runnable_step=RunnableStepJava(
                    name='vcf_filter_filtered',
                    sub_command=Command(program='at.ac.oeaw.cemm.bsf.vcffilter.vep.vep2vcf'),
                    java_temporary_path=runnable_annotate.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G'))
            """ @type _runnable_step: RunnableStepJava """
            _runnable_step.add_option_short(
                key='classpath',
                value=os.path.join(self.classpath_vcf_filter, 'VCFFilter.jar'))
            _sub_command = _runnable_step.sub_command
            _sub_command.add_option_pair(key='INPUT', value=file_path_annotate.vep_filtered_raw_vcf_bgz)
            _sub_command.add_option_pair(key='OUTPUT', value=file_path_annotate.vep_filtered_vcf_bgz)

            return runnable_annotate

        # Start of the run() method body.

        super(VariantCallingGATK, self).run()

        # Get global defaults.

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

        if not self.vep_assembly:
            self.vep_assembly = EnsemblVEP.get_name_assembly(genome_version=self.genome_version)
            if not self.vep_assembly:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'vep_assembly' configuration option.")

        if not self.vep_cache:
            self.vep_cache = EnsemblVEP.get_directory_cache(genome_version=self.genome_version)
            if not self.vep_cache:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'vep_cache' configuration option.")

        if not self.vep_fasta:
            self.vep_fasta = EnsemblVEP.get_directory_fasta(genome_version=self.genome_version)
            if not self.vep_fasta:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'vep_fasta' configuration option.")

        if not self.vep_plugin:
            self.vep_plugin = EnsemblVEP.get_directory_plugin(genome_version=self.genome_version)
            if not self.vep_plugin:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'vep_plugin' configuration option.")

        if not self.vep_source:
            self.vep_source = EnsemblVEP.get_directory_source(genome_version=self.genome_version)
            if not self.vep_source:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'vep_source' configuration option.")

        if not self.vep_species:
            self.vep_species = EnsemblVEP.get_name_species(genome_version=self.genome_version)
            if not self.vep_species:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'vep_species' configuration option.")

        if not self.vep_sql_user:
            self.vep_sql_user = EnsemblVEP.get_sql_user(genome_version=self.genome_version)

        if not self.vep_sql_pass:
            self.vep_sql_pass = EnsemblVEP.get_sql_pass(genome_version=self.genome_version)

        if not self.vep_sql_host:
            self.vep_sql_host = EnsemblVEP.get_sql_host(genome_version=self.genome_version)

        if not self.vep_sql_port:
            self.vep_sql_port = EnsemblVEP.get_sql_port(genome_version=self.genome_version)

        if not self.classpath_gatk:
            self.classpath_gatk = JavaClassPath.get_gatk()
            if not self.classpath_gatk:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'classpath_gatk' configuration option.")

        if not self.classpath_picard:
            self.classpath_picard = JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'classpath_picard' configuration option.")

        if not self.classpath_snpeff:
            self.classpath_snpeff = JavaClassPath.get_snpeff()
            if not self.classpath_snpeff:
                raise Exception("A 'VariantCallingGATK' analysis requires a 'classpath_snpeff' configuration option.")

        if not self.classpath_vcf_filter:
            self.classpath_vcf_filter = JavaClassPath.get_vcf_filter()
            if not self.classpath_vcf_filter:
                raise Exception("A 'VariantCallingGATK' analysis requires a "
                                "'classpath_vcf_filter' configuration option.")

        # Check for absolute paths and adjust if required before checking for existence.

        self.bwa_genome_db = self.configuration.get_absolute_path(
            file_path=self.bwa_genome_db,
            default_path=self.get_gatk_bundle_path)
        if not os.path.exists(path=self.bwa_genome_db):
            raise Exception('The bwa_genome_db file {!r} does not exist.'.format(self.bwa_genome_db))

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
            file_path = self.configuration.get_absolute_path(
                file_path=file_path,
                default_path=Default.absolute_projects())
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The accessory_cohort_gvcf file {!r} does not exist.'.format(file_path))
                # TODO: Check the cohorts so that their sample names do not clash.
        self.accessory_cohort_gvcfs = temporary_list

        for key in self.annotation_resources_dict.keys():
            file_path, annotation_list = self.annotation_resources_dict[key]
            file_path = self.configuration.get_absolute_path(
                file_path=file_path,
                default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                self.annotation_resources_dict[key] = file_path, annotation_list
            else:
                raise Exception('The file path {!r} for annotation resource {!r} does not exist.'.
                                format(file_path, key))

        if self.known_sites_discovery:
            self.known_sites_discovery = self.configuration.get_absolute_path(
                file_path=self.known_sites_discovery,
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(self.known_sites_discovery):
                raise Exception('The known_sites_discovery file {!r} does not exist.'.
                                format(self.known_sites_discovery))

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for file_path in self.known_sites_realignment:
            file_path = self.configuration.get_absolute_path(
                file_path=file_path,
                default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The file path {!r} for known_sites_realignment does not exist.'.
                                format(file_path))
        self.known_sites_realignment = temporary_list

        temporary_list = list()
        """ @type temporary_list: list[str | unicode] """
        for file_path in self.known_sites_recalibration:
            file_path = self.configuration.get_absolute_path(
                file_path=file_path,
                default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The file path {!r} for known_sites_recalibration does not exist.'.
                                format(file_path))
        self.known_sites_recalibration = temporary_list

        for key in self.vqsr_resources_indel_dict:
            resource_dict = self.vqsr_resources_indel_dict[key]
            resource_dict['file_path'] = self.configuration.get_absolute_path(
                file_path=resource_dict['file_path'],
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(resource_dict['file_path']):
                raise Exception('The file path {!r} for vqsr_resources_indel {!r} does not exist.'.
                                format(resource_dict['file_path'], key))

        for key in self.vqsr_resources_snp_dict:
            resource_dict = self.vqsr_resources_snp_dict[key]
            resource_dict['file_path'] = self.configuration.get_absolute_path(
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
            self.genome_annotation_gtf = self.configuration.get_absolute_path(
                file_path=self.genome_annotation_gtf,
                default_path=Default.absolute_intervals())

        # Read comparisons for somatic mutation calling.
        run_read_comparisons()

        # Create genomic tiles for scatter gather approaches.
        if self.number_of_tiles_cohort:
            self._tile_region_cohort_list = run_create_genome_tiles(tiles=self.number_of_tiles_cohort, width=0)

        if self.number_of_tiles_somatic:
            self._tile_region_somatic_list = run_create_genome_tiles(tiles=self.number_of_tiles_somatic, width=0)

        stage_align_lane = self.get_stage(name=self.stage_name_align_lane)
        stage_process_lane = self.get_stage(name=self.stage_name_process_lane)
        stage_process_sample = self.get_stage(name=self.stage_name_process_sample)
        stage_diagnose_sample = self.get_stage(name=self.stage_name_diagnose_sample)
        stage_merge_cohort = self.get_stage(name=self.stage_name_merge_cohort)
        stage_process_cohort = self.get_stage(name=self.stage_name_process_cohort)
        stage_annotate_cohort_snpeff = self.get_stage(name=self.stage_name_annotate_cohort_snpeff)
        stage_annotate_cohort_vep = self.get_stage(name=self.stage_name_annotate_cohort_vep)
        stage_split_cohort_snpeff = self.get_stage(name=self.stage_name_split_cohort_snpeff)
        stage_split_cohort_vep = self.get_stage(name=self.stage_name_split_cohort_vep)
        stage_summary = self.get_stage(name=self.stage_name_summary)
        stage_somatic = self.get_stage(name=self.stage_name_somatic)
        stage_annotate_somatic_snpeff = self.get_stage(name=self.stage_name_annotate_somatic_snpeff)
        stage_annotate_somatic_vep = self.get_stage(name=self.stage_name_annotate_somatic_vep)
        stage_split_somatic_snpeff = self.get_stage(name=self.stage_name_split_somatic_snpeff)
        stage_split_somatic_vep = self.get_stage(name=self.stage_name_split_somatic_vep)

        # Create a Python dict of Python str (cohort name) key and Python list of process_sample Runnable object
        # value data. This dictionary is required by the merge_cohort stage to hierarchically merge cohorts.

        runnable_process_sample_dict = dict()
        """ @type runnable_process_sample_dict: dict[str, list[Runnable]] """

        # Create a Python list of diagnose_sample Runnable objects.

        runnable_diagnose_sample_list = list()
        """ @type runnable_diagnose_sample_list: list[Runnable] """

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

            runnable_process_read_group_list = list()
            """ @type runnable_process_read_group_list: list[Runnable] """

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
                pickler = pickle.Pickler(file=pickler_file, protocol=pickle.HIGHEST_PROTOCOL)
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
                runnable_process_read_group_list.append(runnable_process_lane)

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
                        key='OPTICAL_DUPLICATE_PIXEL_DISTANCE',
                        value='5000')
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

            target_intervals = VariantCallingGATKTargetIntervals.from_sample(sample=sample)

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
            for runnable_process_lane in runnable_process_read_group_list:
                executable_process_sample.dependencies.append(runnable_process_lane.name)

            reference_process_sample = runnable_process_sample.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            if len(runnable_process_read_group_list) == 1:
                # If there is only one read group, sample-level read processing can be skipped.
                # Rename files on the basis of the first and only list component.
                runnable_process_lane = runnable_process_read_group_list[0]
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
                for runnable_process_lane in runnable_process_read_group_list:
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
                        key='OPTICAL_DUPLICATE_PIXEL_DISTANCE',
                        value='5000')
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
            if target_intervals.targets_path:
                # If target intervals are available, the Haplotype Caller analysis is run only on them.
                runnable_step.add_gatk_option(key='intervals', value=target_intervals.targets_path)
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            else:
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

            # Finally, record the process_sample Runnable for the merge_cohort stage under the Sample's
            # 'Cohort Name' annotation, or if it does not exist, under the cohort name defined in the
            # Analysis in the configuration file.

            if 'Cohort Name' in sample.annotation_dict:
                cohort_key = sample.annotation_dict['Cohort Name'][0]
            else:
                cohort_key = self.cohort_name

            if cohort_key not in runnable_process_sample_dict:
                runnable_process_sample_dict[cohort_key] = list()
            _runnable_process_sample_list = runnable_process_sample_dict[cohort_key]
            _runnable_process_sample_list.append(runnable_process_sample)

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

            prefix_diagnose_sample = '_'.join((stage_diagnose_sample.name, sample.name))

            file_path_diagnose_sample = FilePathDiagnoseSample(prefix=prefix_diagnose_sample)

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
            # Set dependencies for succeeding Runnable or Executable objects.
            runnable_diagnose_sample_list.append(runnable_diagnose_sample)

            reference_diagnose_sample = runnable_diagnose_sample.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

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
            if target_intervals.targets_path:
                # If target intervals are available, the Callable Loci analysis is run only on them.
                runnable_step.add_gatk_option(key='intervals', value=target_intervals.targets_path)
            else:
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
            runnable_step.add_option_pair_short(key='type', value='bed4')
            runnable_step.arguments.append(file_path_diagnose_sample.sorted_bed)
            runnable_step.arguments.append(reference_diagnose_sample + '.fai')
            runnable_step.arguments.append(file_path_diagnose_sample.callable_bb)

            # Run the bsfR bsf_variant_calling_coverage.R script.

            runnable_step = runnable_diagnose_sample.add_runnable_step(
                runnable_step=RunnableStep(
                    name='diagnose_sample_coverage',
                    program='bsf_variant_calling_coverage.R'))
            """ @type runnable_step: RunnableStep """
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

            runnable_step = runnable_diagnose_sample.add_runnable_step(
                runnable_step=RunnableStep(
                    name='diagnose_sample_insert_size',
                    program='bsf_variant_calling_insert_size.R'))
            """ @type runnable_step: RunnableStep """
            runnable_step.add_option_long(key='file-path', value=file_path_process_sample.realigned_bam)

            if target_intervals.targets_path:
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

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='diagnose_sample_picard_collect_hybrid_selection_metrics',
                        java_temporary_path=runnable_diagnose_sample.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='CollectHsMetrics'))
                """ @type runnable_step: RunnableStepPicard """
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

        runnable_merge_cohort_dict = {self.cohort_name: []}
        runnable_merge_cohort_list = runnable_merge_cohort_dict[self.cohort_name]

        # Run the GATK CombineGVCFs analysis for each cohort and Sample defined in this project to build up
        # cohort-specific GVCF files.

        for cohort_key in runnable_process_sample_dict.keys():
            runnable_merge_cohort_list.append(
                run_merge_cohort_scatter_gather(
                    analysis_stage=stage_merge_cohort,
                    cohort_runnable_dict=runnable_process_sample_dict,
                    cohort_name=cohort_key))

        # Run the GATK CombineGVCF analysis once more to merge all cohort-specific GVCF files defined in this project.

        if len(runnable_merge_cohort_list) == 1:
            # If the cohort-specific Runnable list has only one component, the merge has already been completed.
            runnable_merge_cohort = runnable_merge_cohort_list[-1]
        elif len(runnable_merge_cohort_list) > 1:
            runnable_merge_cohort = run_merge_cohort_scatter_gather(
                analysis_stage=stage_merge_cohort,
                cohort_runnable_dict=runnable_merge_cohort_dict,
                cohort_name=self.cohort_name)
        else:
            raise Exception('Unexpected number of Runnable objects on the merge_cohort list.')

        # Run an additional GATK CombineGVCFs analysis to merge into a super-cohort, if defined.

        if len(self.accessory_cohort_gvcfs):
            # If accessory cohorts are defined, initialise a new Python dict of Python str cohort name key and
            # Python list of Runnable value data. Initialise the list with teh last Runnable object and
            # extend with the the list of accessory cohort file names. The run_merge_cohort_scatter_gather() method
            # can cope with Runnable or str | unicode objects.
            cohort_key = '_'.join((self.cohort_name, 'accessory'))
            runnable_merge_cohort_dict = {cohort_key: [runnable_merge_cohort]}
            runnable_merge_cohort_list = runnable_merge_cohort_dict[cohort_key]
            runnable_merge_cohort_list.extend(self.accessory_cohort_gvcfs)

            runnable_merge_cohort = run_merge_cohort_scatter_gather(
                analysis_stage=stage_merge_cohort,
                cohort_runnable_dict=runnable_merge_cohort_dict,
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

        runnable_process_cohort_gather = run_genotype_cohort_scatter_gather()

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
            file_path_process_cohort.recalibrated_snp_raw_indel_tbi = file_path_genotype_cohort.genotyped_raw_tbi
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

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_variant_recalibrator_indel',
                    java_temporary_path=runnable_process_cohort.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_process_cohort)
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
            file_path_process_cohort.multi_sample_tbi = \
                file_path_process_cohort.recalibrated_snp_recalibrated_indel_tbi

        ################################
        # Step 7: Annotate the cohort. #
        ################################

        # Run snpEff annotation.

        prefix_annotate_cohort_snpeff = '_'.join((stage_annotate_cohort_snpeff.name, self.cohort_name))

        runnable_annotate_cohort_snpeff = run_annotate_snpeff(
            prefix=prefix_annotate_cohort_snpeff,
            vcf_file_path=file_path_process_cohort.multi_sample_vcf)
        executable_annotate_cohort_snpeff = self.set_stage_runnable(
            stage=stage_annotate_cohort_snpeff,
            runnable=runnable_annotate_cohort_snpeff)
        # Set dependencies on preceding Runnable.name or Executable.name objects.
        executable_annotate_cohort_snpeff.dependencies.append(runnable_process_cohort.name)

        # Run Ensembl Variant Effect Predictor (VEP) annotation.

        prefix_annotate_cohort_vep = '_'.join((stage_annotate_cohort_vep.name, self.cohort_name))

        runnable_annotate_cohort_vep = run_annotate_vep(
            prefix=prefix_annotate_cohort_vep,
            vcf_file_path=file_path_process_cohort.multi_sample_vcf)
        executable_annotate_cohort_vep = self.set_stage_runnable(
            stage=stage_annotate_cohort_vep,
            runnable=runnable_annotate_cohort_vep)
        # Set dependencies on preceding Runnable.name or Executable.name objects.
        executable_annotate_cohort_vep.dependencies.append(runnable_process_cohort.name)

        ######################################################
        # Step 8: Re-process and split the cohort by sample. #
        ######################################################
        #
        # GATK SelectVariants   (split_cohort_snpeff_gatk_select_variants_snpeff)
        # GATK VariantsToTable  (split_cohort_snpeff_gatk_variants_to_table_snpeff)
        # GATK SelectVariants   (split_cohort_gatk_select_variants_vep)
        # GATK VariantsToTable  (split_cohort_gatk_variants_to_table_vep)

        file_path_annotate_cohort_snpeff = runnable_annotate_cohort_snpeff.file_path_object
        """ @type file_path_annotate_cohort_snpeff: FilePathAnnotateSnpEff """

        file_path_annotate_cohort_vep = runnable_annotate_cohort_vep.file_path_object
        """ @type file_path_annotate_cohort_vep: FilePathAnnotateVEP """

        for sample in self.sample_list:
            # Get all PairedReads objects solely to exclude samples without any.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)
            if not len(paired_reads_dict):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            # Split the snpEff-annotated multi-sample VCF file.

            prefix_split_cohort_snpeff = '_'.join((stage_split_cohort_snpeff.name, sample.name))

            file_path_split_cohort_snpeff = FilePathSplitCohort(prefix=prefix_split_cohort_snpeff)

            runnable_split_cohort_snpeff = self.add_runnable(
                runnable=Runnable(
                    name=prefix_split_cohort_snpeff,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_split_cohort_snpeff,
                    debug=self.debug))
            executable_split_cohort_snpeff = self.set_stage_runnable(
                stage=stage_split_cohort_snpeff,
                runnable=runnable_split_cohort_snpeff)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_split_cohort_snpeff.dependencies.append(runnable_annotate_cohort_snpeff.name)

            reference_split_cohort_snpeff = runnable_split_cohort_snpeff.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the GATK SelectVariants analysis to split multi-sample VCF files into one per sample.

            runnable_step = runnable_split_cohort_snpeff.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_snpeff_gatk_select_variants_snpeff',
                    java_temporary_path=runnable_split_cohort_snpeff.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_snpeff)
            runnable_step.add_gatk_option(key='variant', value=file_path_annotate_cohort_snpeff.annotated_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_snpeff.sample_vcf)
            runnable_step.add_gatk_option(key='sample_name', value=sample.name)
            runnable_step.add_gatk_switch(key='excludeNonVariants')

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_split_cohort_snpeff.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_snpeff_gatk_variants_to_table_snpeff',
                    java_temporary_path=runnable_split_cohort_snpeff.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_snpeff)
            runnable_step.add_gatk_option(key='variant', value=file_path_split_cohort_snpeff.sample_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_snpeff.sample_tsv)
            runnable_step.add_gatk_switch(key='allowMissingData')
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
            for annotation_resource in self.annotation_resources_dict.keys():
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        runnable_step.add_gatk_option(
                            key='fields',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

            # Split the Ensembl Variant Effect Predictor-annotated multi-sample VCF file.

            prefix_split_cohort_vep = '_'.join((stage_split_cohort_vep.name, sample.name))

            file_path_split_cohort_vep = FilePathSplitCohort(prefix=prefix_split_cohort_vep)

            runnable_split_cohort_vep = self.add_runnable(
                runnable=Runnable(
                    name=prefix_split_cohort_vep,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_split_cohort_vep,
                    debug=self.debug))
            executable_split_cohort_vep = self.set_stage_runnable(
                stage=stage_split_cohort_vep,
                runnable=runnable_split_cohort_vep)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_split_cohort_vep.dependencies.append(runnable_annotate_cohort_vep.name)

            reference_split_cohort_vep = runnable_split_cohort_vep.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the GATK SelectVariants analysis to split multi-sample VCF files into one per sample.

            runnable_step = runnable_split_cohort_vep.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_vep_gatk_select_variants_vep',
                    java_temporary_path=runnable_split_cohort_vep.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_vep)
            runnable_step.add_gatk_option(key='variant', value=file_path_annotate_cohort_vep.vep_complete_vcf_bgz)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_vep.sample_vcf)
            runnable_step.add_gatk_option(key='sample_name', value=sample.name)
            runnable_step.add_gatk_switch(key='excludeNonVariants')

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_split_cohort_vep.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_vep_gatk_variants_to_table_vep',
                    java_temporary_path=runnable_split_cohort_vep.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_cohort_vep)
            runnable_step.add_gatk_option(key='variant', value=file_path_split_cohort_vep.sample_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_cohort_vep.sample_tsv)
            runnable_step.add_gatk_switch(key='allowMissingData')
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
        for runnable_diagnose_sample in runnable_diagnose_sample_list:
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

        comparison_name_list = self._comparison_dict.keys()
        comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

        if self.debug > 0:
            print 'Somatic variant calling comparison:', comparison_name_list

        for comparison_name in comparison_name_list:
            prefix_somatic = '_'.join((stage_somatic.name, comparison_name))

            # Somatic variant calling-specific file paths

            file_path_somatic = FilePathSomatic(prefix=prefix_somatic)

            runnable_somatic_gather = run_somatic_scatter_gather(comparison_key=comparison_name)

            # Create a Runnable for processing the somatic calls.

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
            # For the moment, this Runnable has no RunnableStep objects assigned.
            # It is purely needed for accessing FilePath objects in the report() method below.
            # Never submit the corresponding executable.
            executable_somatic.submit = False

            # reference_somatic = runnable_somatic.get_absolute_cache_file_path(
            #     file_path=self.bwa_genome_db)

            ############################################
            # Step 11: Annotate somatic variant calls. #
            ############################################

            # Run snpEff annotation.

            prefix_annotate_somatic_snpeff = '_'.join((stage_annotate_somatic_snpeff.name, comparison_name))

            runnable_annotate_somatic_snpeff = run_annotate_snpeff(
                prefix=prefix_annotate_somatic_snpeff,
                vcf_file_path=file_path_somatic.somatic_vcf)
            executable_annotate_somatic_snpeff = self.set_stage_runnable(
                stage=stage_annotate_somatic_snpeff,
                runnable=runnable_annotate_somatic_snpeff)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_annotate_somatic_snpeff.dependencies.append(runnable_somatic_gather.name)

            file_path_annotate_somatic_snpeff = runnable_annotate_somatic_snpeff.file_path_object
            """ @type file_path_annotate_somatic_snpeff: FilePathAnnotateSnpEff """

            # Run Ensembl Variant Effect Predictor (VEP) annotation.

            prefix_annotate_somatic_vep = '_'.join((stage_annotate_somatic_vep.name, comparison_name))

            runnable_annotate_somatic_vep = run_annotate_vep(
                prefix=prefix_annotate_somatic_vep,
                vcf_file_path=file_path_somatic.somatic_vcf)
            executable_annotate_somatic_vep = self.set_stage_runnable(
                stage=stage_annotate_somatic_vep,
                runnable=runnable_annotate_somatic_vep)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_annotate_somatic_vep.dependencies.append(runnable_somatic_gather.name)

            file_path_annotate_somatic_vep = runnable_annotate_somatic_vep.file_path_object
            """ @type file_path_annotate_somatic_vep: FilePathAnnotateVEP """

            #########################################
            # Step 12: Split somatic variant calls. #
            #########################################

            # Split the somatic snpEff-annotated VCF file into a TSV file.

            prefix_split_somatic_snpeff = '_'.join((stage_split_somatic_snpeff.name, comparison_name))

            file_path_split_somatic_snpeff = FilePathSplitSomatic(prefix=prefix_split_somatic_snpeff)

            runnable_split_somatic_snpeff = self.add_runnable(
                runnable=Runnable(
                    name=prefix_split_somatic_snpeff,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_split_somatic_snpeff,
                    debug=self.debug))
            executable_split_somatic_snpeff = self.set_stage_runnable(
                stage=stage_split_somatic_snpeff,
                runnable=runnable_split_somatic_snpeff)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_split_somatic_snpeff.dependencies.append(runnable_annotate_somatic_snpeff.name)

            reference_split_somatic_snpeff = runnable_split_somatic_snpeff.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_split_somatic_snpeff.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_variants_to_table',
                    java_temporary_path=runnable_split_somatic_snpeff.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_somatic_snpeff)
            runnable_step.add_gatk_option(key='variant', value=file_path_annotate_somatic_snpeff.annotated_vcf)
            runnable_step.add_gatk_option(key='out', value=file_path_split_somatic_snpeff.comparison_tsv)
            runnable_step.add_gatk_switch(key='allowMissingData')
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
            for annotation_resource in self.annotation_resources_dict.keys():
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        runnable_step.add_gatk_option(
                            key='fields',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

            # Split the somatic Ensembl VEP-annotated VCF file into a TSV file.

            prefix_split_somatic_vep = '_'.join((stage_split_somatic_vep.name, comparison_name))

            file_path_split_somatic_vep = FilePathSplitSomatic(prefix=prefix_split_somatic_vep)

            runnable_split_somatic_vep = self.add_runnable(
                runnable=Runnable(
                    name=prefix_split_somatic_vep,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory,
                    cache_path_dict=self._cache_path_dict,
                    file_path_object=file_path_split_somatic_vep,
                    debug=self.debug))
            executable_split_somatic_vep = self.set_stage_runnable(
                stage=stage_split_somatic_vep,
                runnable=runnable_split_somatic_vep)
            # Set dependencies on preceding Runnable.name or Executable.name objects.
            executable_split_somatic_vep.dependencies.append(runnable_annotate_somatic_vep.name)

            reference_split_somatic_vep = runnable_split_somatic_vep.get_absolute_cache_file_path(
                file_path=self.bwa_genome_db)

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_split_somatic_vep.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_variants_to_table',
                    java_temporary_path=runnable_split_somatic_vep.get_relative_temporary_directory_path,
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            """ @type runnable_step: RunnableStepGATK """
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=reference_split_somatic_vep)
            runnable_step.add_gatk_option(key='variant', value=file_path_annotate_somatic_vep.vep_complete_vcf_bgz)
            runnable_step.add_gatk_option(key='out', value=file_path_split_somatic_vep.comparison_tsv)
            runnable_step.add_gatk_switch(key='allowMissingData')
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

        def report_link():
            """Private function to create simpler symbolic links in structured sub-directories.

            The simplified symbolic links facilitate file system browsing.
            @return:
            @rtype:
            """
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

            # Process per sample.

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
                # paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                runnable_process_sample = self.runnable_dict[
                    '_'.join((self.stage_name_process_sample, sample.name))]
                file_path_process_sample = runnable_process_sample.file_path_object
                """ @type file_path_process_sample: FilePathProcessSample """

                # runnable_diagnose_sample = self.runnable_dict[
                #     '_'.join((self.stage_name_diagnose_sample, sample.name))]
                # file_path_diagnose_sample = runnable_diagnose_sample.file_path_object
                """ @type file_path_diagnose_sample: FilePathDiagnoseSample """

                runnable_annotate_cohort_snpeff = self.runnable_dict[
                    '_'.join((self.stage_name_annotate_cohort_snpeff, self.cohort_name))]
                file_path_annotate_cohort_snpeff = runnable_annotate_cohort_snpeff.file_path_object
                """ @type file_path_annotate_cohort_snpeff: FilePathAnnotateSnpEff """

                runnable_annotate_cohort_vep = self.runnable_dict[
                    '_'.join((self.stage_name_annotate_cohort_vep, self.cohort_name))]
                file_path_annotate_cohort_vep = runnable_annotate_cohort_vep.file_path_object
                """ @type file_path_annotate_cohort_vep: FilePathAnnotateVep """

                runnable_split_cohort_snpeff = self.runnable_dict[
                    '_'.join((self.stage_name_split_cohort_snpeff, sample.name))]
                file_path_split_cohort_snpeff = runnable_split_cohort_snpeff.file_path_object
                """ @type file_path_split_cohort_snpeff: FilePathSplitCohort """

                runnable_split_cohort_vep = self.runnable_dict[
                    '_'.join((self.stage_name_split_cohort_vep, sample.name))]
                file_path_split_cohort_vep = runnable_split_cohort_vep.file_path_object
                """ @type file_path_split_cohort_vep: FilePathSplitCohort """

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
                        ('snpeff_vcf_bgz', '_snpeff.vcf.gz'),
                        ('snpeff_vcf_tbi', '_snpeff.vcf.gz.tbi'),
                        ('snpeff_stats', '_snpeff_summary.html')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_cohort_snpeff, attribute)),
                            directory_results_by_cohort),
                        target_path=os.path.join(directory_results_by_cohort, self.cohort_name + extension))

                for attribute, extension in (
                        ('vep_complete_vcf_bgz', '_vep.vcf.gz'),
                        ('vep_complete_vcf_tbi', '_vep.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_cohort_vep, attribute)),
                            directory_results_by_cohort),
                        target_path=os.path.join(directory_results_by_cohort, self.cohort_name + extension))

            # Process per (somatic) comparison

            for comparison_name in self._comparison_dict.keys():
                runnable_somatic = self.runnable_dict[
                    '_'.join((self.stage_name_somatic, comparison_name))]
                file_path_somatic = runnable_somatic.file_path_object
                """ @type file_path_somatic: FilePathSomatic """

                for attribute, extension in (
                        ('somatic_vcf', '.vcf.gz'),
                        ('somatic_tbi', '.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_somatic, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                runnable_annotate_somatic_snpeff = self.runnable_dict[
                    '_'.join((self.stage_name_annotate_somatic_snpeff, comparison_name))]
                file_path_annotate_somatic_snpeff = runnable_annotate_somatic_snpeff.file_path_object
                """ @type file_path_annotate_somatic_snpeff: FilePathAnnotateSnpEff """

                for attribute, extension in (
                        ('snpeff_vcf_bgz', '_snpeff.vcf.gz'),
                        ('snpeff_vcf_tbi', '_snpeff.vcf.gz.tbi'),
                        ('snpeff_genes', '_snpeff_summary.genes.txt'),
                        ('snpeff_stats', '_snpeff_summary.html'),
                        ('annotated_vcf', '_annotated.vcf.gz'),
                        ('annotated_tbi', '_annotated.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_somatic_snpeff, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                runnable_split_somatic_snpeff = self.runnable_dict[
                    '_'.join((self.stage_name_split_somatic_snpeff, comparison_name))]
                file_path_split_somatic_snpeff = runnable_split_somatic_snpeff.file_path_object
                """ @type file_path_split_somatic_snpeff: FilePathSplitSomatic """

                for attribute, extension in (
                        ('comparison_tsv', '_snpeff.tsv'),):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_split_somatic_snpeff, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                runnable_annotate_somatic_vep = self.runnable_dict[
                    '_'.join((self.stage_name_annotate_somatic_vep, comparison_name))]
                file_path_annotate_somatic_vep = runnable_annotate_somatic_vep.file_path_object
                """ @type file_path_annotate_somatic_vep: FilePathAnnotateVEP """

                for attribute, extension in (
                        ('vep_complete_vcf_bgz', '_vep.vcf.gz'),
                        ('vep_complete_vcf_tbi', '_vep.vcf.gz.tbi')):
                    report_create_symbolic_link(
                        source_path=os.path.relpath(
                            os.path.join(self.genome_directory, getattr(file_path_annotate_somatic_vep, attribute)),
                            directory_results_by_pair),
                        target_path=os.path.join(directory_results_by_pair, comparison_name + extension))

                runnable_split_somatic_vep = self.runnable_dict[
                    '_'.join((self.stage_name_split_somatic_vep, comparison_name))]
                file_path_split_somatic_vep = runnable_split_somatic_vep.file_path_object
                """ @type file_path_split_somatic_vep: FilePathSplitSomatic """

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

            str_list += '<h2 id="genome_browsing">Genome Browsing</h2>\n'
            str_list += '\n'

            str_list += '<p id="ucsc_track_hub">'
            str_list += self.ucsc_hub_html_anchor(link_path=link_path)
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<h2 id="read_group_and_sample_level">Read Group and Sample Level</h2>\n'
            str_list += '\n'
            str_list += '<table id="read_group_and_sample_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Sample</th>\n'
            str_list += '<th>Variants<br />snpEff<br />Ensembl&nbsp;VEP</th>\n'
            str_list += '<th>Alignments</th>\n'
            str_list += '<th>Read Group</th>\n'
            str_list += '<th>Duplicate Metrics</th>\n'
            str_list += '<th>Alignment Summary Metrics</th>\n'
            str_list += '<th>Hybrid Selection Metrics</th>\n'
            str_list += '<th>Non-Callable Loci</th>\n'
            str_list += '<th>Non-Callable Summary</th>\n'
            str_list += '<th>Insert Size</th>\n'
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

                runnable_process_sample = self.runnable_dict[
                    '_'.join((self.stage_name_process_sample, sample.name))]
                file_path_process_sample = runnable_process_sample.file_path_object
                """ @type file_path_process_sample: FilePathProcessSample """

                runnable_diagnose_sample = self.runnable_dict[
                    '_'.join((self.stage_name_diagnose_sample, sample.name))]
                file_path_diagnosis = runnable_diagnose_sample.file_path_object
                """ @type file_path_diagnosis: FilePathDiagnoseSample """

                runnable_split_cohort_snpeff = self.runnable_dict[
                    '_'.join((self.stage_name_split_cohort_snpeff, sample.name))]
                file_path_split_cohort_snpeff = runnable_split_cohort_snpeff.file_path_object
                """ @type file_path_split_cohort_snpeff: FilePathSplitCohort """

                runnable_split_cohort_vep = self.runnable_dict[
                    '_'.join((self.stage_name_split_cohort_vep, sample.name))]
                file_path_split_cohort_vep = runnable_split_cohort_vep.file_path_object
                """ @type file_path_split_cohort_vep: FilePathSplitCohort """

                str_list += '<tr>\n'
                # Sample
                str_list += '<td class="left">' + sample.name + '</td>\n'
                # Variants
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_split_cohort_snpeff.sample_vcf + '">'
                str_list += '<abbr title="Variant Calling Format">VCF</abbr>'
                str_list += '</a>&nbsp;'
                str_list += '<a href="' + file_path_split_cohort_snpeff.sample_tbi + '">'
                str_list += '<abbr title="Tabix Index">TBI</abbr>'
                str_list += '</a>&nbsp;'
                str_list += '<a href="' + file_path_split_cohort_snpeff.sample_tsv + '">'
                str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                str_list += '</a><br />'
                str_list += '<a href="' + file_path_split_cohort_vep.sample_vcf + '">'
                str_list += '<strong><abbr title="Variant Calling Format">VCF</abbr></strong>'
                str_list += '</a>&nbsp;'
                str_list += '<a href="' + file_path_split_cohort_vep.sample_tbi + '">'
                str_list += '<abbr title="Tabix Index">TBI</abbr>'
                str_list += '</a>&nbsp;'
                str_list += '<a href="' + file_path_split_cohort_vep.sample_tsv + '">'
                str_list += '<strong><abbr title="Tab-Separated Value">TSV</abbr></strong>'
                str_list += '</a>'
                str_list += '</td>\n'
                # Alignments
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_process_sample.realigned_bam + '">'
                str_list += '<abbr title="Binary Alignment/Map">BAM</abbr>'
                str_list += '</a>&nbsp;'
                str_list += '<a href="' + file_path_process_sample.realigned_bai + '">'
                str_list += '<abbr title="Binary Alignment/Map Index">BAI</abbr>'
                str_list += '</a>'
                str_list += '</td>\n'
                # Read Group
                str_list += '<td class="left"></td>\n'
                # Duplicate Metrics
                str_list += '<td class="center">'
                # This can be a sample-specific file or a symbolic link to the read group-specific file.
                if os.path.exists(
                        os.path.join(
                            self.genome_directory,
                            file_path_process_sample.duplicate_metrics)):
                    str_list += '<a href="' + file_path_process_sample.duplicate_metrics + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a>'
                str_list += '</td>\n'
                # Alignment Summary Metrics
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_process_sample.alignment_summary_metrics + '">'
                str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                str_list += '</a>'
                str_list += '</td>\n'
                # Hybrid Selection Metrics
                str_list += '<td class="center">'
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnosis.hybrid_selection_metrics)):
                    str_list += '<a href="' + file_path_diagnosis.hybrid_selection_metrics + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a>'
                str_list += '</td>\n'
                # Non-Callable Loci
                str_list += '<td class="center">'
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnosis.non_callable_loci_tsv)):
                    str_list += '<a href="' + file_path_diagnosis.callable_bed + '">'
                    str_list += '<abbr title="Browser Extensible Data">BED</abbr>'
                    str_list += '</a>&nbsp;'
                    str_list += '<a href="' + file_path_diagnosis.callable_bb + '">'
                    str_list += '<abbr title="Big Browser Extensible Data">BigBED</abbr>'
                    str_list += '</a>&nbsp;'
                    # Do not link the more complex file_path_diagnosis.non_callable_regions_tsv
                    # file for the moment.
                    str_list += '<a href="' + file_path_diagnosis.non_callable_regions_tsv + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a>'
                str_list += '</td>\n'
                # Non-Callable Summary
                str_list += '<td class="center">'
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnosis.non_callable_summary_tsv)):
                    str_list += '<a href="' + file_path_diagnosis.non_callable_summary_tsv + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a>'
                str_list += '</td>\n'
                # Insert size
                str_list += '<td class="center">'
                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnosis.insert_size_tsv)):
                    str_list += '<a href="' + file_path_diagnosis.insert_size_pdf + '">'
                    str_list += '<img alt="Insert Size per Sample ' + sample.name + '"'
                    str_list += ' src="' + file_path_diagnosis.insert_size_png + '"'
                    str_list += ' height="80" width="80" />'
                    str_list += '</a>'
                    str_list += '<a href="' + file_path_diagnosis.insert_size_tsv + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a>'
                str_list += '</td>\n'
                str_list += '</tr>\n'

                for paired_reads_name in paired_reads_name_list:
                    runnable_process_lane = self.runnable_dict[
                        '_'.join((self.stage_name_process_lane, paired_reads_name))]
                    file_path_process_read_group = runnable_process_lane.file_path_object
                    """ @type file_path_process_read_group: FilePathProcessReadGroup """

                    str_list += '<tr>\n'
                    # Sample
                    str_list += '<td class="left"></td>\n'
                    # Variants
                    str_list += '<td class="center"></td>\n'
                    # Alignments
                    str_list += '<td class="center"></td>\n'
                    # Read Group
                    str_list += '<td class="left">' + paired_reads_name + '</td>\n'
                    # Duplicate Metrics
                    str_list += '<td class="center">'
                    if os.path.isfile(
                            os.path.join(
                                self.genome_directory,
                                file_path_process_read_group.duplicate_metrics)):
                        str_list += '<a href="' + file_path_process_read_group.duplicate_metrics + '">'
                        str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                        str_list += '</a>'
                    str_list += '</td>\n'
                    # Alignment Summary Metrics
                    str_list += '<td class="center">'
                    str_list += '<a href="' + file_path_process_read_group.alignment_summary_metrics + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a>'
                    str_list += '</td>\n'
                    # Hybrid Selection Metrics
                    str_list += '<td class="center"></td>\n'
                    # Non-Callable Loci
                    str_list += '<td class="center"></td>\n'
                    # Non-Callable Summary
                    str_list += '<td class="center"></td>\n'
                    # Insert Size
                    str_list += '<td class="center"></td>\n'
                    str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            str_list += '<h2 id="cohort_level">Cohort Level</h2>\n'
            str_list += '\n'
            str_list += '<table id="cohort_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Cohort</th>\n'
            str_list += '<th>Information</th>\n'
            str_list += '<th>Comment</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            runnable_annotate_cohort_snpeff = self.runnable_dict[
                '_'.join((self.stage_name_annotate_cohort_snpeff, self.cohort_name))]
            file_path_annotate_cohort_snpeff = runnable_annotate_cohort_snpeff.file_path_object
            """ @type file_path_annotate_cohort_snpeff: FilePathAnnotateSnpEff """

            str_list += '<tr>\n'
            str_list += '<td class="left">' + self.cohort_name + '</td>\n'
            str_list += '<td class="left">'
            str_list += '<a href="' + file_path_annotate_cohort_snpeff.snpeff_stats + '">'
            str_list += 'snpEff Summary Statistics'
            str_list += '</a>'
            str_list += '</td>\n'
            str_list += '<td class="left">'
            str_list += '<a href="' + file_path_annotate_cohort_snpeff.snpeff_genes + '">'
            str_list += 'snpEff Summary Genes'
            str_list += '</a>'
            str_list += '</td>\n'
            str_list += '</tr>\n'

            str_list += '<tr>\n'
            str_list += '<td class="left">' + self.cohort_name + '</td>\n'
            str_list += '<td class="left">'
            str_list += 'snpEff-annotated multi-sample '
            str_list += '<a href="' + file_path_annotate_cohort_snpeff.snpeff_vcf_bgz + '">'
            str_list += '<abbr title="Variant Calling Format">VCF</abbr>'
            str_list += '</a> and '
            str_list += '<a href="' + file_path_annotate_cohort_snpeff.snpeff_vcf_tbi + '">'
            str_list += '<abbr title="Tabix Index">TBI</abbr>'
            str_list += '</a>'
            str_list += '</td>\n'
            str_list += '<td class="left">Functional annotation of all splice variants</td>\n'
            str_list += '</tr>\n'

            str_list += '<tr>\n'
            str_list += '<td class="left">' + self.cohort_name + '</td>\n'
            str_list += '<td class="left">'
            str_list += 'GATK-annotated multi-sample '
            str_list += '<a href="' + file_path_annotate_cohort_snpeff.annotated_vcf + '">'
            str_list += '<abbr title="Variant Calling Format">VCF</abbr>'
            str_list += '</a> and '
            str_list += '<a href="' + file_path_annotate_cohort_snpeff.annotated_tbi + '">'
            str_list += '<abbr title="Tabix Index">TBI</abbr>'
            str_list += '</a>'
            str_list += '</td>\n'
            str_list += '<td class="left">'
            str_list += 'Functional annotation of only the most severely affected splice variant'
            str_list += '</td>\n'
            str_list += '</tr>\n'

            runnable_annotate_cohort_vep = self.runnable_dict[
                '_'.join((self.stage_name_annotate_cohort_vep, self.cohort_name))]
            file_path_annotate_cohort_vep = runnable_annotate_cohort_vep.file_path_object
            """ @type file_path_annotate_cohort_vep: FilePathAnnotateVEP """

            str_list += '<tr>\n'
            str_list += '<td class="left">' + self.cohort_name + '</td>\n'
            str_list += '<td class="left">'
            str_list += '<a href="' + file_path_annotate_cohort_vep.vep_statistics + '">'
            str_list += 'Ensembl Variant Effect Predictor Summary Statistics'
            str_list += '</a>'
            str_list += '</td>\n'
            str_list += '<td class="left">'
            str_list += '</td>\n'
            str_list += '</tr>\n'

            str_list += '<tr>\n'
            str_list += '<td class="left">' + self.cohort_name + '</td>\n'
            str_list += '<td class="left">'
            str_list += 'Ensembl VEP-annotated multi-sample '
            str_list += '<a href="' + file_path_annotate_cohort_vep.vep_complete_vcf_bgz + '">'
            str_list += '<abbr title="Variant Calling Format">VCF</abbr>'
            str_list += '</a> and '
            str_list += '<a href="' + file_path_annotate_cohort_vep.vep_complete_vcf_tbi + '">'
            str_list += '<abbr title="Tabix Index">TBI</abbr>'
            str_list += '</a>'
            str_list += '</td>\n'
            str_list += '<td class="left">Functional annotation of all Ensembl splice variants</td>\n'
            str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            # Somatic variant calling.

            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            if len(comparison_name_list):
                str_list += '<h2 id="somatic_variants">Somatic Variants</h2>\n'
                str_list += '\n'
                str_list += '<table id="somatic_variants_table">\n'
                str_list += '<thead>\n'
                str_list += '<tr>\n'
                str_list += '<th>Comparison</th>\n'
                str_list += '<th>Variants<br />snpEff<br />VEP</th>\n'
                str_list += '<th>Summary&nbsp;Statistics<br />snpEff<br />VEP</th>\n'
                str_list += '</tr>\n'
                str_list += '</thead>\n'
                str_list += '<tbody>\n'

                for comparison_name in comparison_name_list:
                    # runnable_somatic = self.runnable_dict[
                    #     '_'.join((self.stage_name_somatic, comparison_name))]
                    # file_path_somatic = runnable_somatic.file_path_object
                    # """ @type file_path_somatic: FilePathSomatic """

                    # Add snpEff annotation

                    runnable_annotate_somatic_snpeff = self.runnable_dict[
                        '_'.join((self.stage_name_annotate_somatic_snpeff, comparison_name))]
                    file_path_annotate_somatic_snpeff = runnable_annotate_somatic_snpeff.file_path_object
                    """ @type file_path_annotate_somatic_snpeff: FilePathAnnotateSnpEff """

                    runnable_split_somatic_snpeff = self.runnable_dict[
                        '_'.join((self.stage_name_split_somatic_snpeff, comparison_name))]
                    file_path_split_somatic_snpeff = runnable_split_somatic_snpeff.file_path_object
                    """ @type file_path_split_somatic_snpeff: FilePathSplitSomatic """

                    # Add VEP annotation

                    runnable_annotate_somatic_vep = self.runnable_dict[
                        '_'.join((self.stage_name_annotate_somatic_vep, comparison_name))]
                    file_path_annotate_somatic_vep = runnable_annotate_somatic_vep.file_path_object
                    """ @type file_path_annotate_somatic_vep: FilePathAnnotateVEP """

                    runnable_split_somatic_vep = self.runnable_dict[
                        '_'.join((self.stage_name_split_somatic_vep, comparison_name))]
                    file_path_split_somatic_vep = runnable_split_somatic_vep.file_path_object
                    """ @type file_path_split_somatic_vep: FilePathSplitSomatic """

                    str_list += '<tr>\n'

                    # Comparison
                    str_list += '<td class="left">' + comparison_name + '</td>\n'

                    # Variants
                    str_list += '<td class="center">'
                    str_list += '<a href="' + file_path_annotate_somatic_snpeff.annotated_vcf + '">'
                    str_list += '<abbr title="Variant Calling Format">VCF</abbr>'
                    str_list += '</a>&nbsp;'
                    str_list += '<a href="' + file_path_annotate_somatic_snpeff.annotated_tbi + '">'
                    str_list += '<abbr title="Tabix Index">TBI</abbr>'
                    str_list += '</a>&nbsp;'
                    str_list += '<a href="' + file_path_split_somatic_snpeff.comparison_tsv + '">'
                    str_list += '<abbr title="Tab-Separated Value">TSV</abbr>'
                    str_list += '</a><br />'
                    str_list += '<a href="' + file_path_annotate_somatic_vep.vep_complete_vcf_bgz + '">'
                    str_list += '<strong><abbr title="Variant Calling Format">VCF</abbr></strong>'
                    str_list += '</a>&nbsp;'
                    str_list += '<a href="' + file_path_annotate_somatic_vep.vep_complete_vcf_tbi + '">'
                    str_list += '<abbr title="Tabix Index">TBI</abbr>'
                    str_list += '</a>&nbsp;'
                    str_list += '<a href="' + file_path_split_somatic_vep.comparison_tsv + '">'
                    str_list += '<strong><abbr title="Tab-Separated Value">TSV</abbr></strong>'
                    str_list += '</a>'
                    str_list += '</td>\n'

                    # Summary Statistics
                    str_list += '<td class="center">'
                    str_list += '<a href="' + file_path_annotate_somatic_snpeff.snpeff_stats + '">'
                    str_list += '<abbr title="Hyper Text Markup Language">HTML</abbr>'
                    str_list += '</a>&nbsp;'
                    str_list += '<a href="' + file_path_annotate_somatic_snpeff.snpeff_genes + '">'
                    str_list += '<abbr title="Text">TXT</abbr>'
                    str_list += '</a><br />'
                    str_list += '<a href="' + file_path_annotate_somatic_vep.vep_statistics + '">'
                    str_list += '<strong><abbr title="Hyper Text Markup Language">HTML</abbr></strong>'
                    str_list += '</a>'
                    str_list += '</td>\n'

                    str_list += '</tr>\n'

                str_list += '</tbody>\n'
                str_list += '</table>\n'
                str_list += '\n'

            str_list += '<h2 id="qc_plots">QC Plots</h2>\n'
            str_list += '\n'
            str_list += '<table id="qc_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Sample</th>\n'
            str_list += '<th>Read Group</th>\n'
            str_list += '<th>Metrics</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            runnable_summary = self.runnable_dict['_'.join((self.stage_name_summary, self.cohort_name))]
            file_path_summary = runnable_summary.file_path_object
            """ @type file_path_summary: FilePathSummary """

            # Alignment Summary - TSV
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.alignment_metrics_sample_tsv)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.alignment_metrics_sample_tsv + '">TSV</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.alignment_metrics_read_group_tsv + '">TSV</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">'
                str_list += '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' \
                            '#AlignmentSummaryMetrics">Alignment Summary</a>'
                str_list += '</td>\n'
                str_list += '</tr>\n'

            # Alignment Summary - Percent Aligned
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.alignment_percentage_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.alignment_percentage_sample_pdf + '">'
                str_list += '<img alt="Alignment Summary - Percent Aligned per Sample"'
                str_list += ' src="' + file_path_summary.alignment_percentage_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.alignment_percentage_read_group_pdf + '">'
                str_list += '<img alt="Alignment Summary - Percent Aligned per Read Group"'
                str_list += ' src="' + file_path_summary.alignment_percentage_read_group_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">Alignment Summary - Percent Aligned</td>\n'
                str_list += '</tr>\n'

            # Alignment Summary - Reads Aligned
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.alignment_absolute_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.alignment_absolute_sample_pdf + '">'
                str_list += '<img alt="Alignment Summary - Reads Aligned per Sample"'
                str_list += ' src="' + file_path_summary.alignment_absolute_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.alignment_absolute_read_group_pdf + '">'
                str_list += '<img alt="Alignment Summary - Reads Aligned per Read Group"'
                str_list += ' src="' + file_path_summary.alignment_absolute_read_group_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">Alignment Summary - Reads Aligned</td>\n'
                str_list += '</tr>\n'

            # Duplication - TSV
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.duplication_metrics_sample_tsv)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.duplication_metrics_sample_tsv + '">TSV</a>'
                str_list += '</td>\n'
                str_list += '<td class="center"></td>\n'
                str_list += '<td class="left">'
                str_list += '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' \
                            '#DuplicationMetrics">Duplication</a>'
                str_list += '</td>\n'
                str_list += '</tr>\n'

            # Duplication - Fraction
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.duplication_percentage_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.duplication_percentage_sample_pdf + '">'
                str_list += '<img alt="Duplication - Duplicated Reads per Sample"'
                str_list += ' src="' + file_path_summary.duplication_percentage_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center"></td>\n'
                str_list += '<td class="left">Duplication - Fraction</td>\n'
                str_list += '</tr>\n'

            # Duplication - Levels
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.duplication_levels_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.duplication_levels_sample_pdf + '">'
                str_list += '<img alt="Duplication - Duplication Levels per Sample"'
                str_list += ' src="' + file_path_summary.duplication_levels_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center"></td>\n'
                str_list += '<td class="left">Duplication - Levels</td>\n'
                str_list += '</tr>\n'

            # Hybrid Selection - TSV
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.hybrid_metrics_sample_tsv)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_metrics_sample_tsv + '">TSV</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_metrics_read_group_tsv + '">TSV</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">'
                str_list += '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html' \
                            '#HsMetrics">Hybrid Selection</a>'
                str_list += '</td>\n'
                str_list += '</tr>\n'

            # Hybrid Selection - Target Coverage Levels
            if os.path.exists(
                    os.path.join(
                        self.genome_directory,
                        file_path_summary.hybrid_coverage_levels_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_coverage_levels_sample_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Mean Target Coverage Levels per Sample"'
                str_list += ' src="' + file_path_summary.hybrid_coverage_levels_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_coverage_levels_read_group_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Mean Target Coverage Levels per Read Group"'
                str_list += ' src="' + file_path_summary.hybrid_coverage_levels_read_group_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">Hybrid Selection - Mean Target Coverage Levels</td>\n'
                str_list += '</tr>\n'

            # Hybrid Selection - Mean Target Coverage
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.hybrid_coverage_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_coverage_sample_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Mean Target Coverage per Sample"'
                str_list += ' src="' + file_path_summary.hybrid_coverage_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_coverage_read_group_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Mean Target Coverage per Read Group"'
                str_list += ' src="' + file_path_summary.hybrid_coverage_read_group_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">Hybrid Selection - Mean Target Coverage</td>\n'
                str_list += '</tr>\n'

            # Hybrid Selection - Excluded Bases
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.hybrid_excluded_bases_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_excluded_bases_sample_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Percent Excluded Bases per Sample"'
                str_list += ' src="' + file_path_summary.hybrid_excluded_bases_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_excluded_bases_read_group_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Percent Excluded Bases per Read Group"'
                str_list += ' src="' + file_path_summary.hybrid_excluded_bases_read_group_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">Hybrid Selection - Percent Excluded Bases</td>\n'
                str_list += '</tr>\n'

            # Hybrid Selection - Percent Unique Reads
            # The plot is meaningless if Picard MarkDuplicates has not run.
            if not self.skip_mark_duplicates and os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.hybrid_unique_percentage_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_unique_percentage_sample_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Percent Unique Reads per Sample"'
                str_list += ' src="' + file_path_summary.hybrid_unique_percentage_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.hybrid_unique_percentage_read_group_pdf + '">'
                str_list += '<img alt="Hybrid Selection - Percent Unique Reads per Read Group"'
                str_list += ' src="' + file_path_summary.hybrid_unique_percentage_read_group_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="left">Hybrid Selection - Percent Unique Reads</td>\n'
                str_list += '</tr>\n'

            # Non-Callable Loci - TSV
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.non_callable_metrics_sample_tsv)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.non_callable_metrics_sample_tsv + '">TSV</a>'
                str_list += '</td>\n'
                str_list += '<td class="center"></td>\n'
                str_list += '<td class="left">Non-Callable Loci</td>\n'
                str_list += '</tr>\n'

            # Non-Callable Loci - Fraction
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.non_callable_percentage_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.non_callable_percentage_sample_pdf + '">'
                str_list += '<img alt="Non-Callable Loci - Fraction"'
                str_list += ' src="' + file_path_summary.non_callable_percentage_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center"></td>\n'
                str_list += '<td class="left">Non-Callable Loci - Fraction</td>\n'
                str_list += '</tr>\n'

            # Non-Callable Loci - Number
            if os.path.exists(os.path.join(
                    self.genome_directory,
                    file_path_summary.non_callable_absolute_sample_png)):
                str_list += '<tr>\n'
                str_list += '<td class="center">'
                str_list += '<a href="' + file_path_summary.non_callable_absolute_sample_pdf + '">'
                str_list += '<img alt="Non-Callable Loci - Number"'
                str_list += ' src="' + file_path_summary.non_callable_absolute_sample_png + '"'
                str_list += ' height="100" width="100" />'
                str_list += '</a>'
                str_list += '</td>\n'
                str_list += '<td class="center"></td>\n'
                str_list += '<td class="left">Non-Callable Loci - Number</td>\n'
                str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            self.report_to_file(content=str_list)

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
            str_list += 'longLabel BWA NGS read alignments\n'
            str_list += 'superTrack on show\n'
            str_list += 'group alignments\n'
            str_list += '\n'

            str_list += 'track Callable\n'
            str_list += 'shortLabel Callable\n'
            str_list += 'longLabel Callable\n'
            str_list += 'superTrack on show\n'
            str_list += 'group callable\n'
            str_list += '\n'

            str_list += 'track Variants\n'
            str_list += 'shortLabel Variants\n'
            str_list += 'longLabel Variant calls\n'
            str_list += 'superTrack on show\n'
            str_list += 'group variants\n'
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
                # paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                runnable_process_sample = self.runnable_dict[
                    '_'.join((self.stage_name_process_sample, sample.name))]
                file_path_process_sample = runnable_process_sample.file_path_object
                """ @type file_path_process_sample: FilePathProcessSample """

                runnable_diagnose_sample = self.runnable_dict[
                    '_'.join((self.stage_name_diagnose_sample, sample.name))]
                file_path_diagnose_sample = runnable_diagnose_sample.file_path_object
                """ @type file_path_diagnose_sample: FilePathDiagnoseSample """

                runnable_split_cohort_snpeff = self.runnable_dict[
                    '_'.join((self.stage_name_split_cohort_snpeff, sample.name))]
                file_path_split_cohort_snpeff = runnable_split_cohort_snpeff.file_path_object
                """ @type file_path_split_cohort_snpeff: FilePathSplitCohort """

                runnable_split_cohort_vep = self.runnable_dict[
                    '_'.join((self.stage_name_split_cohort_vep, sample.name))]
                file_path_split_cohort_vep = runnable_split_cohort_vep.file_path_object
                """ @type file_path_split_cohort_vep: FilePathSplitCohort """

                #
                #  Alignments track
                #
                # Common settings
                str_list += 'track ' + sample.name + '_alignments\n'
                str_list += 'type bam\n'
                str_list += 'shortLabel ' + sample.name + '_alignments\n'
                str_list += 'longLabel ' + sample.name + 'BWA NGS read alignments\n'
                str_list += 'bigDataUrl ' + file_path_process_sample.realigned_bam + '\n'
                # str_list += 'html ...\n'
                str_list += 'visibility squish\n'

                # Common optional settings
                str_list += 'color 0,0,0\n'

                # bam/cram - Compressed Sequence Alignment track settings
                # None

                # Composite track settings
                str_list += 'parent Alignments\n'
                str_list += '\n'

                if os.path.isfile(
                        os.path.join(
                            self.genome_directory,
                            file_path_diagnose_sample.callable_bb)):
                    #
                    #  Callable track
                    #
                    # Common settings
                    str_list += 'track ' + sample.name + '_callable\n'
                    str_list += 'type bigBed\n'
                    str_list += 'shortLabel ' + sample.name + '_callable\n'
                    str_list += 'longLabel ' + sample.name + 'callable\n'
                    str_list += 'bigDataUrl ' + file_path_diagnose_sample.callable_bb + '\n'
                    # str_list += 'html ...\n'
                    str_list += 'visibility squish\n'

                    # Common optional settings
                    str_list += 'color 0,0,0\n'

                    # bigBed - Item or region track settings
                    # None

                    # Composite track settings
                    str_list += 'parent Callable\n'
                    str_list += '\n'

                #
                # snpEff Variants track
                #
                # Common settings
                str_list += 'track ' + sample.name + '_snpeff\n'
                str_list += 'type vcfTabix\n'
                str_list += 'shortLabel ' + sample.name + '_snpeff\n'
                str_list += 'longLabel ' + sample.name + ' snpEff-annotated variant calls\n'
                str_list += 'bigDataUrl ' + file_path_split_cohort_snpeff.sample_vcf + '\n'
                # str_list += 'html ...\n'
                str_list += 'visibility dense\n'

                # Common optional settings

                # vcfTabix specific settings

                # Composite track settings
                str_list += 'parent Variants\n'
                str_list += '\n'

                #
                # Ensembl VEP Variants track
                #
                # Common settings
                str_list += 'track ' + sample.name + '_vep\n'
                str_list += 'type vcfTabix\n'
                str_list += 'shortLabel ' + sample.name + '_vep\n'
                str_list += 'longLabel ' + sample.name + ' Ensembl VEP-annotated variant calls\n'
                str_list += 'bigDataUrl ' + file_path_split_cohort_vep.sample_vcf + '\n'
                # str_list += 'html ...\n'
                str_list += 'visibility dense\n'

                # Common optional settings

                # vcfTabix specific settings

                # Composite track settings
                str_list += 'parent Variants\n'
                str_list += '\n'

            # Comparison-specific tracks

            self.ucsc_hub_to_file(content=str_list)

            return

        report_link()
        report_html()
        report_hub()

        return
