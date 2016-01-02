"""bsf.analyses.variant_calling

A package of classes and methods supporting variant calling analyses.
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

import os.path
from pickle import Pickler, HIGHEST_PROTOCOL
import warnings

from bsf import Analysis, Command, Configuration, Default, defaults, DRMS, Executable, Runnable, RunnableStep, \
    RunnableStepLink, RunnableStepMove, RunnableStepJava, RunnableStepPicard
from bsf.annotation import SampleAnnotationSheet
from bsf.data import PairedReads
from bsf.executables import BWA


class RunnableStepGATK(RunnableStepJava):
    """The C{RunnableStepGATK} class represents a C{RunnableStepJava} specific to the Genome Analysis Toolkit (GATK).

    Attributes:
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None,
                 java_temporary_path=None, java_heap_maximum=None, java_jar_path=None,
                 gatk_classpath=None):
        """Create a C{RunnableStep} for a GATK algorithm.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @return: C{RunnableStep}
        @rtype: RunnableStep
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
        """Add an option to the GATK command.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_option_long(key=key, value=value, override=override)

    def add_gatk_switch(self, key):
        """Add a switch to the GATK command.

        @param key: Option key
        @type key: str
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_switch_long(key=key)


class VariantCallingGATK(Analysis):
    """The C{VariantCallingGATK} class represents the logic to run the Genome Analysis Toolkit (GATK).

    Attributes:
    @cvar drms_name_align_lane: C{DRMS.name} for the lane alignment C{Analysis} stage
    @type drms_name_align_lane: str
    @cvar drms_name_process_lane: C{DRMS.name} for the lane processing C{Analysis} stage
    @type drms_name_process_lane: str
    @cvar drms_name_process_sample: C{DRMS.name} for the sample processing C{Analysis} stage
    @type drms_name_process_sample: str
    @cvar drms_name_diagnose_sample: C{DRMS.name} for the sample diagnosis C{Analysis} stage
    @type drms_name_diagnose_sample: str
    @cvar drms_name_process_cohort: C{DRMS.name} for the cohort processing C{Analysis} stage
    @type drms_name_process_cohort: str
    @cvar drms_name_split_cohort: C{DRMS.name} for the cohort splitting C{Analysis} stage
    @type drms_name_split_cohort: str
    @cvar drms_name_somatic: C{DRMS.name} for the somatic C{Analysis} stage
    @type drms_name_somatic: str
    @ivar replicate_grouping: Group individual C{PairedReads} objects for processing or run them separately
    @type replicate_grouping: bool
    @ivar comparison_path: Comparison file
    @type comparison_path: str | unicode
    @ivar cohort_name: Cohort name
    @type cohort_name: str
    @ivar accessory_cohort_gvcfs: Python C{list} of Python C{str} or C{unicode} (GVCF file path) objects
    @type accessory_cohort_gvcfs: list[str | unicode]
    @ivar skip_mark_duplicates: Skip the Picard MarkDuplicates step
    @type skip_mark_duplicates: bool
    @ivar known_sites_discovery: VCF file path for variant discovery via The Haplotype Caller or Unified Genotyper
    @type known_sites_discovery: str | unicode
    @ivar known_sites_realignment: Python C{list} of Python C{str} or C{unicode} (VCF file paths)
        for realignment
    @type known_sites_realignment: list[str | unicode]
    @ivar known_sites_recalibration: Python C{list} of Python C{str} or C{unicode} (VCF file paths)
        for recalibration
    @type known_sites_recalibration: list[str | unicode]
    @ivar annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
        Python C{tuple} of
        Python C{str} or C{unicode} (file path) and Python C{list} of Python C{str} (annotation) value data
    @type annotation_resources_dict: dict[str, (str | unicode, list[str])]
    @ivar truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
    @type truth_sensitivity_filter_level_indel: str
    @ivar truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
    @type truth_sensitivity_filter_level_snp: str
    @ivar vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
    @type vqsr_skip_indel: bool, None
    @ivar vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
    @type vqsr_skip_snp: bool, None
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
    @ivar downsample_to_fraction: Down-sample to fraction
    @type downsample_to_fraction: str
    @ivar gatk_bundle_version: GATK resource bundle version
    @type gatk_bundle_version: str
    @ivar snpeff_genome_version: snpEff genome version
    @type snpeff_genome_version: str
    @ivar classpath_gatk: Genome Analysis Tool Kit Java Archive (JAR) class path directory
    @type classpath_gatk: str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar classpath_snpeff: snpEff tool Java Archive (JAR) class path directory
    @type classpath_snpeff: str | unicode
    """

    drms_name_align_lane = 'variant_calling_align_lane'
    drms_name_process_lane = 'variant_calling_process_lane'
    drms_name_process_sample = 'variant_calling_process_sample'
    drms_name_diagnose_sample = 'variant_calling_diagnose_sample'
    drms_name_merge_cohort = 'variant_calling_merge_cohort'
    drms_name_process_cohort = 'variant_calling_process_cohort'
    drms_name_split_cohort = 'variant_calling_split_cohort'
    drms_name_somatic = 'variant_calling_somatic'

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 replicate_grouping=False, bwa_genome_db=None, comparison_path=None,
                 cohort_name=None, accessory_cohort_gvcfs=None, skip_mark_duplicates=False,
                 known_sites_discovery=None, known_sites_realignment=None, known_sites_recalibration=None,
                 known_somatic_discovery=None,
                 annotation_resources_dict=None,
                 truth_sensitivity_filter_level_indel=None,
                 truth_sensitivity_filter_level_snp=None,
                 vqsr_skip_indel=False, vqsr_skip_snp=False,
                 vqsr_resources_indel_dict=None, vqsr_resources_snp_dict=None,
                 vqsr_annotations_indel_list=None, vqsr_annotations_snp_list=None,
                 vqsr_bad_lod_cutoff_indel=None, vqsr_bad_lod_cutoff_snp=None,
                 vqsr_max_gaussians_pos_indel=4, vqsr_max_gaussians_pos_snp=None,
                 exclude_intervals_list=None,
                 include_intervals_list=None,
                 interval_padding=None,
                 downsample_to_fraction=None,
                 gatk_bundle_version=None, snpeff_genome_version=None,
                 classpath_gatk=None, classpath_picard=None, classpath_snpeff=None):
        """Initialise a C{VariantCallingGATK} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{Analysis}-wide project directory,
            normally under the C{Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{Analysis}-wide genome directory,
            normally under the C{Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param drms_list: Python C{list} of C{DRMS} objects
        @type drms_list: list[DRMS]
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{str} key and Python C{list} objects of C{Sample} objects
        @type comparisons: dict[str, list[Sample]]
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list[Sample]
        @param replicate_grouping: Group individual C{PairedReads} objects for processing or run them separately
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
        @param downsample_to_fraction: Down-sample to fraction
        @type downsample_to_fraction: str
        @param gatk_bundle_version: GATK resource bundle version
        @type gatk_bundle_version: str
        @param snpeff_genome_version: snpEff genome version
        @type snpeff_genome_version: str
        @param classpath_gatk: Genome Analysis Tool Kit Java Archive (JAR) class path directory
        @type classpath_gatk: str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param classpath_snpeff: snpEff tool Java Archive (JAR) class path directory
        @type classpath_snpeff: str | unicode
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
            drms_list=drms_list,
            collection=collection,
            comparisons=comparisons,
            samples=samples)

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
            self.interval_padding = int()
        else:
            assert isinstance(interval_padding, int)
            self.interval_padding = interval_padding

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

        return

    @property
    def get_gatk_bundle_path(self):
        """Get the absolute GATK bundle directory C{Default.absolute_gatk_bundle} for the set
        C{VariantCallingGATK.gatk_bundle_version} and C{VariantCallingGATK.genome_version}.

        @return: Absolute GATK bundle directory
        @rtype: str | unicode
        """
        return Default.absolute_gatk_bundle(
            gatk_bundle_version=self.gatk_bundle_version,
            genome_version=self.genome_version)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{VariantCallingGATK} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(VariantCallingGATK, self).set_configuration(configuration=configuration, section=section)

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        # Get the genome database.

        option = 'bwa_genome_db'
        if configuration.config_parser.has_option(section=section, option=option):
            self.bwa_genome_db = configuration.config_parser.get(section=section, option=option)

        # Read a comparison file.

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

        # Get the cohort name.

        option = 'cohort_name'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cohort_name = configuration.config_parser.get(section=section, option=option)

        # Comma-separated list of GVCF files from accessory cohorts
        # that should be used in the recalibration procedure.

        option = 'accessory_cohort_gvcfs'
        if configuration.config_parser.has_option(section=section, option=option):
            for file_path in configuration.config_parser.get(section=section, option=option).split(','):
                self.accessory_cohort_gvcfs.append(file_path.strip())  # Strip white space around commas.

        # Get the skip mark duplicates option.

        option = 'skip_mark_duplicates'
        if configuration.config_parser.has_option(section=section, option=option):
            self.skip_mark_duplicates = configuration.config_parser.getboolean(section=section, option=option)

        # Get the truth sensitivity filter level for INDELs.

        option = 'truth_sensitivity_filter_level_indel'
        if configuration.config_parser.has_option(section=section, option=option):
            self.truth_sensitivity_filter_level_indel = configuration.config_parser.get(section=section, option=option)

        # Get the truth sensitivity filter level for SNPs.

        option = 'truth_sensitivity_filter_level_snp'
        if configuration.config_parser.has_option(section=section, option=option):
            self.truth_sensitivity_filter_level_snp = configuration.config_parser.get(section=section, option=option)

        # Get the flag for skipping the Variant Quality Score Recalibration (VQSR) for INDELs.

        option = 'vqsr_skip_indel'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vqsr_skip_indel = configuration.config_parser.getboolean(section=section, option=option)

        # Get the flag for skipping the Variant Quality Score Recalibration (VQSR) for SNPs.

        option = 'vqsr_skip_snp'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vqsr_skip_snp = configuration.config_parser.getboolean(section=section, option=option)

        # Get the list of annotations for the Variant Quality Score Recalibration (VQSR) for INDELs.

        option = 'vqsr_annotations_indel'
        if configuration.config_parser.has_option(section=section, option=option):
            for annotation in configuration.config_parser.get(section=section, option=option).split(','):
                self.vqsr_annotations_indel_list.append(annotation.strip())  # Strip white space around commas.

        # Get the list of annotations for the Variant Quality Score Recalibration (VQSR) for SNPs.

        option = 'vqsr_annotations_snp'
        if configuration.config_parser.has_option(section=section, option=option):
            for annotation in configuration.config_parser.get(section=section, option=option).split(','):
                self.vqsr_annotations_snp_list.append(annotation.strip())  # Strip white space around commas.

        # Get the bad LOD cutoff for the negative training set for INDELs.

        option = 'vqsr_bad_lod_cutoff_indel'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vqsr_bad_lod_cutoff_indel = configuration.config_parser.getfloat(section=section, option=option)

        # Get the bad LOD cutoff for the negative training set for SNPs.

        option = 'vqsr_bad_lod_cutoff_snp'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vqsr_bad_lod_cutoff_snp = configuration.config_parser.getfloat(section=section, option=option)

        # Get the maximum number of Gaussians in the positive training for INDELs.

        option = 'vqsr_max_gaussians_pos_indel'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vqsr_max_gaussians_pos_indel = configuration.config_parser.getint(section=section, option=option)

        # Get the maximum number of Gaussians in the positive training for SNPs.

        option = 'vqsr_max_gaussians_pos_snp'
        if configuration.config_parser.has_option(section=section, option=option):
            self.vqsr_max_gaussians_pos_snp = configuration.config_parser.getint(section=section, option=option)

        self._read_vqsr_configuration(vqsr_resources_dict=self.vqsr_resources_indel_dict, variation_type='indel')
        self._read_vqsr_configuration(vqsr_resources_dict=self.vqsr_resources_snp_dict, variation_type='snp')

        # Read additionally requested annotation resources for the GATK AnnotateVariants step.

        # Python dict of Python str (annotation resource name) key and
        # Python tuple of
        # Python str (file path) and Python list of Python str (annotation) value data.

        option = 'annotation_resources'
        if configuration.config_parser.has_option(section=section, option=option):
            for annotation_resource in configuration.config_parser.get(section=section, option=option).split(','):
                annotation_resource = annotation_resource.strip()  # Strip white space around commas.
                resource_section = '_'.join((annotation_resource, 'resource'))
                if configuration.config_parser.has_section(section=resource_section):
                    annotation_list = list()
                    if configuration.config_parser.has_option(section=resource_section, option='file_path'):
                        file_path = configuration.config_parser.get(section=resource_section, option='file_path')
                    else:
                        raise Exception(
                            "Missing configuration option 'file_path' in configuration section {!r}.".
                            format(resource_section))
                    if configuration.config_parser.has_option(
                            section=resource_section,
                            option='annotations'):
                        for annotation in configuration.config_parser.get(
                                section=resource_section,
                                option='annotations').split(','):
                            annotation_list.append(annotation.strip())  # Strip white space around commas.
                    else:
                        raise Exception(
                            "Missing configuration option 'annotations' in configuration section {!r}.".
                            format(resource_section))
                    # Create a dict key and a tuple of a Python str and Python list.
                    if annotation_resource not in self.annotation_resources_dict:
                        self.annotation_resources_dict[annotation_resource] = file_path, annotation_list
                else:
                    raise Exception(
                        'Missing configuration section {!r} declared in option annotation_resources {!r}.'.
                        format(resource_section,
                               configuration.config_parser.get(section=section, option='annotation_resources')))

        # Single VCF file of known sites for the
        # GATK HaplotypeCaller and GenotypeGVCFs steps.

        option = 'known_sites_discovery'
        if configuration.config_parser.has_option(section=section, option=option):
            self.known_sites_discovery = configuration.config_parser.get(section=section, option=option)

        # Comma-separated list of VCF files with known variant sites for the
        # GATK RealignerTargetCreator and IndelRealigner steps.

        option = 'known_sites_realignment'
        if configuration.config_parser.has_option(section=section, option=option):
            for file_path in configuration.config_parser.get(section=section, option=option).split(','):
                self.known_sites_realignment.append(file_path.strip())  # Strip white space around commas.

        # Comma-separated list of VCF files with known variant sites for the
        # GATK BaseRecalibrator and PrintReads steps.

        option = 'known_sites_recalibration'
        if configuration.config_parser.has_option(section=section, option=option):
            for file_path in configuration.config_parser.get(section=section, option=option).split(','):
                self.known_sites_recalibration.append(file_path.strip())  # Strip white space around commas.

        # Comma-separated list of VCF files with known somatic variant sites for the
        # GATK MuTect2 steps.

        option = 'known_somatic_discovery'
        if configuration.config_parser.has_option(section=section, option=option):
            for file_path in configuration.config_parser.get(section=section, option=option).split(','):
                self.known_somatic_discovery.append(file_path.strip())  # Strip white space around commas.

        # Get the list of intervals to exclude.

        option = 'exclude_intervals'
        if configuration.config_parser.has_option(section=section, option=option):
            exclude_intervals = configuration.config_parser.get(section=section, option=option)
            # For comma-separated interval lists split into components on commas, strip white space
            # and push them onto the list individually.
            self.exclude_intervals_list.extend(map(lambda x: x.strip(), exclude_intervals.split(',')))

        # Get the list of intervals to include.

        option = 'include_intervals'
        if configuration.config_parser.has_option(section=section, option=option):
            include_intervals = configuration.config_parser.get(section=section, option=option)
            # For comma-separated interval lists split into components on commas, strip white space and
            # push them onto the list individually.
            self.include_intervals_list.extend(map(lambda x: x.strip(), include_intervals.split(',')))

        # Get the interval padding.

        option = 'interval_padding'
        if configuration.config_parser.has_option(section=section, option=option):
            self.interval_padding = configuration.config_parser.getint(section=section, option=option)

        # Get the down-sample to fraction information.

        option = 'downsample_to_fraction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.downsample_to_fraction = configuration.config_parser.get(section=section, option=option)

        # Get the GATK bundle version.

        option = 'gatk_bundle_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.gatk_bundle_version = configuration.config_parser.get(section=section, option=option)

        # Get the snpEff genome version.

        option = 'snpeff_genome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.snpeff_genome_version = configuration.config_parser.get(section=section, option=option)

        # Get the Genome Analysis Tool Kit (GATK) Java Archive (JAR) class path directory.

        option = 'classpath_gatk'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_gatk = configuration.config_parser.get(section=section, option=option)

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        # Get the snpEff tool Java Archive (JAR) class path directory.

        option = 'classpath_snpeff'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_snpeff = configuration.config_parser.get(section=section, option=option)

        return

    def _read_comparisons(self, comparison_path):
        """Read a C{SampleAnnotationSheet} CSV file from disk.

            - Column headers for CASAVA folders:
                - Treatment/Control ProcessedRunFolder:
                    - CASAVA processed run folder name or
                    - C{Analysis.input_directory} by default
                - Treatment/Control Project:
                    - CASAVA Project name or
                    - C{Analysis.project_name} by default
                - Treatment/Control Sample:
                    - CASAVA Sample name, no default
            - Column headers for independent samples:
                - Treatment/Control Sample:
                - Treatment/Control Reads:
                - Treatment/Control File:
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode
        @return:
        @rtype:
        """

        assert isinstance(comparison_path, (str, unicode, None))

        # For variant calling, all samples need adding to the Analysis object regardless.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        if comparison_path:
            sas = SampleAnnotationSheet.from_file_path(file_path=comparison_path, name='Somatic Comparisons')

            for row_dict in sas.row_dicts:
                key = str()
                comparison_groups = list()

                if self.debug > 0:
                    print "Comparison sheet row_dict {!r}".format(row_dict)

                for prefix in 'Normal', 'Tumor':
                    # TODO: Collection.get_samples_from_row_dict() may not work with normalised
                    # Sample Annotation Sheets.
                    group_name, group_samples = self.collection.get_samples_from_row_dict(
                        row_dict=row_dict, prefix=prefix)
                    if group_name and len(group_samples):
                        key += group_name
                        key += '__'
                        # key = '__'.join((key, group_name))
                        comparison_groups.append((group_name, group_samples))
                        # Also expand each Python list of Sample objects to get all those Sample objects
                        # that this Analysis needs considering.
                        for sample in group_samples:
                            if self.debug > 1:
                                print '  {} Sample name: {!r} file_path:{!r}'.\
                                    format(prefix, sample.name, sample.file_path)
                                # print sample.trace(1)
                            # self.add_sample(sample=sample)

                # Remove the last '__' from the key.
                self.comparisons[key[:-2]] = comparison_groups

        return

    def _read_vqsr_configuration(self, vqsr_resources_dict, variation_type=None):
        """Private method to read variant quality score recalibration (VQSR) configuration information.

        @param vqsr_resources_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_dict: dict[str, dict[str, str | unicode]]
        @param variation_type: Variation type I{indel} or I{snp}
        @type variation_type: str
        @return:
        @rtype:
        """

        if variation_type not in ('indel', 'snp'):
            raise Exception("Variation type has to be 'indel' or 'snp', not {!r}.".format(variation_type))

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        resource_option = '_'.join(('vqsr_resources', variation_type))
        if config_parser.has_option(section=config_section, option=resource_option):
            for resource in config_parser.get(section=config_section, option=resource_option).split(','):
                resource = resource.strip()
                resource_section = '_'.join(('vqsr', variation_type, resource))
                if config_parser.has_section(section=resource_section):
                    if resource in vqsr_resources_dict:
                        resource_dict = vqsr_resources_dict[resource]
                        assert isinstance(resource_dict, dict)
                    else:
                        resource_dict = dict()
                        vqsr_resources_dict[resource] = resource_dict
                    if config_parser.has_option(section=resource_section, option='known'):
                        resource_dict['known'] = config_parser.get(section=resource_section, option='known')
                    if config_parser.has_option(section=resource_section, option='training'):
                        resource_dict['training'] = config_parser.get(section=resource_section, option='training')
                    if config_parser.has_option(section=resource_section, option='truth'):
                        resource_dict['truth'] = config_parser.get(section=resource_section, option='truth')
                    if config_parser.has_option(section=resource_section, option='prior'):
                        resource_dict['prior'] = config_parser.get(section=resource_section, option='prior')
                    if config_parser.has_option(section=resource_section, option='file_path'):
                        resource_dict['file_path'] = config_parser.get(section=resource_section, option='file_path')
                else:
                    raise Exception(
                        'Missing configuration section {!r} declared in option {!r} {!r}.'.
                        format(resource_section, resource_option,
                               config_parser.get(section=config_section, option=resource_option)))

        return

    def run(self):
        """Run this C{VariantCallingGATK} analysis.
        @return:
        @rtype:
        """

        # Get global defaults.

        default = Default.get_global_default()

        super(VariantCallingGATK, self).run()

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

        # Check for absolute paths and adjust if required before checking for existence.

        self.bwa_genome_db = Default.get_absolute_path(
            file_path=self.bwa_genome_db,
            default_path=self.get_gatk_bundle_path)
        if not os.path.exists(path=self.bwa_genome_db):
            raise Exception("The bwa_genome_db file {!r} does not exist.".format(self.bwa_genome_db))

        temporary_list = list()
        for file_path in self.accessory_cohort_gvcfs:
            file_path = Default.get_absolute_path(file_path=file_path, default_path=Default.absolute_projects())
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The accessory_cohort_gvcf file {!r} does not exist.'.format(file_path))
            # TODO: Check the cohorts so that their sample names do not clash.
        self.accessory_cohort_gvcfs = temporary_list

        for key in self.annotation_resources_dict.keys():
            assert isinstance(key, str)
            file_path, annotation_list = self.annotation_resources_dict[key]
            assert isinstance(file_path, (str, unicode))
            assert isinstance(annotation_list, list)
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
        for file_path in self.known_sites_realignment:
            assert isinstance(file_path, (str, unicode))
            file_path = Default.get_absolute_path(file_path=file_path, default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The file path {!r} for known_sites_realignment does not exist.'.
                                format(file_path))
        self.known_sites_realignment = temporary_list

        temporary_list = list()
        for file_path in self.known_sites_recalibration:
            assert isinstance(file_path, (str, unicode))
            file_path = Default.get_absolute_path(file_path=file_path, default_path=self.get_gatk_bundle_path)
            if os.path.exists(file_path):
                temporary_list.append(file_path)
            else:
                raise Exception('The file path {!r} for known_sites_recalibration does not exist.'.
                                format(file_path))
        self.known_sites_recalibration = temporary_list

        for key in self.vqsr_resources_indel_dict:
            assert isinstance(key, str)
            resource_dict = self.vqsr_resources_indel_dict[key]
            assert isinstance(resource_dict, dict)
            resource_dict['file_path'] = Default.get_absolute_path(
                file_path=resource_dict['file_path'],
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(resource_dict['file_path']):
                raise Exception('The file path {!r} for vqsr_resources_indel {!r} does not exist.'.
                                format(resource_dict['file_path'], key))

        for key in self.vqsr_resources_snp_dict:
            assert isinstance(key, str)
            resource_dict = self.vqsr_resources_snp_dict[key]
            assert isinstance(resource_dict, dict)
            resource_dict['file_path'] = Default.get_absolute_path(
                file_path=resource_dict['file_path'],
                default_path=self.get_gatk_bundle_path)
            if not os.path.exists(resource_dict['file_path']):
                raise Exception('The file path {!r} for vqsr_resources_snp {!r} does not exist.'.
                                format(resource_dict['file_path'], key))

        # Exclude intervals

        temporary_list = list()
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

        # The comparison file can be absolute, in the same directory or in the project directory.
        if not os.path.isabs(self.comparison_path) and not os.path.exists(self.comparison_path):
            self.comparison_path = Default.get_absolute_path(
                file_path=self.comparison_path,
                default_path=self.project_directory)

        # Read comparisons for somatic mutation calling.
        self._read_comparisons(comparison_path=self.comparison_path)

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # bsf_run_bwa.py script.

        drms_align_lane = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_align_lane,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_process_lane Runnable.

        drms_process_lane = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_process_lane,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_process_sample Runnable.

        drms_process_sample = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_process_sample,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_diagnose_sample Runnable.

        drms_diagnose_sample = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_diagnose_sample,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_merge_cohort Runnable.

        drms_merge_cohort = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_merge_cohort,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_process_cohort Runnable.

        drms_process_cohort = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_process_cohort,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_split_cohort Runnable.

        drms_split_cohort = self.add_drms(
            drms=DRMS.from_analysis(
                name=self.drms_name_split_cohort,
                working_directory=self.genome_directory,
                analysis=self))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_somatic Runnable.

        drms_somatic = self.add_drms(
                drms=DRMS.from_analysis(
                    name=self.drms_name_somatic,
                    working_directory=self.genome_directory,
                    analysis=self))

        vc_merge_cohort_dependency_dict = dict()
        vc_merge_cohort_sample_dict = dict()

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            vc_process_sample_dependencies = list()
            vc_process_sample_replicates = list()

            # Sample.get_all_paired_reads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                # Step 1: Process per lane.

                bwa = BWA(name='variant_calling_bwa_{}'.format(replicate_key), analysis=self)
                # Instead of adding the BWA Executable to the DRMS, it gets serialised into the pickler_file.
                # bwa_drms.add_executable(executable=bwa)

                bwa_mem = bwa.sub_command

                # Set BWA mem options.

                # Allow as many threads as defined in the corresponding DRMS object.
                bwa_mem.add_option_short(key='t', value=str(drms_align_lane.threads))
                # Append FASTA/Q comment to SAM output.
                bwa_mem.add_switch_short(key='C')
                # Mark shorter split hits as secondary (for Picard compatibility).
                bwa_mem.add_switch_short(key='M')
                # Output warnings and errors only.
                bwa_mem.add_option_short(key='v', value='2')

                # Set BWA arguments.

                bwa_mem.arguments.append(self.bwa_genome_db)

                reads1 = list()
                reads2 = list()

                # Propagate the SAM read group information around FASTQ files if required.
                # Please note that only the first read group can be propagated per
                # PairedReads object.

                read_group = str()

                for paired_reads in replicate_dict[replicate_key]:
                    assert isinstance(paired_reads, PairedReads)
                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)
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

                file_path_align_lane = dict(
                    # TODO: The name for the aligned BAM is constructed by the bsf_run_bwa.py script.
                    # It is currently based on the drms_align_lane.name and replicate_key.
                    # The script should also be changed to pre-set all file names beforehand.
                    aligned_bam='{}_{}.bam'.format(drms_align_lane.name, replicate_key),
                    aligned_bai='{}_{}.bai'.format(drms_align_lane.name, replicate_key),
                    aligned_md5='{}_{}.bam.md5'.format(drms_align_lane.name, replicate_key))

                # Normally, the bwa object would be pushed onto the drms list.
                # Experimentally, use Pickler to serialize the Executable object into a file.

                pickler_dict_align_lane = dict()
                pickler_dict_align_lane['prefix'] = drms_align_lane.name
                pickler_dict_align_lane['replicate_key'] = replicate_key
                pickler_dict_align_lane['classpath_gatk'] = self.classpath_gatk
                pickler_dict_align_lane['classpath_picard'] = self.classpath_picard
                pickler_dict_align_lane['bwa_executable'] = bwa

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(drms_align_lane.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_align_lane)
                pickler_file.close()

                # Create a bsf_run_bwa.py job to run the pickled object.

                run_bwa = drms_align_lane.add_executable(executable=Executable.from_analysis(
                    name='_'.join((drms_align_lane.name, replicate_key)),
                    program='bsf_run_bwa.py',
                    analysis=self))

                # Only submit this Executable if the final result file does not exist.
                if (os.path.exists(os.path.join(self.genome_directory, file_path_align_lane['aligned_md5'])) and
                        os.path.getsize(os.path.join(self.genome_directory, file_path_align_lane['aligned_md5']))):
                    run_bwa.submit = False
                # Check also for existence of a new-style Runnable status file.
                if os.path.exists(os.path.join(
                        drms_align_lane.working_directory,
                        '_'.join((drms_align_lane.name, replicate_key, 'completed.txt')))):
                    run_bwa.submit = False

                    # Set run_bwa options.

                run_bwa.add_option_long(key='pickler_path', value=pickler_path)
                run_bwa.add_option_long(key='debug', value=str(self.debug))

                prefix_lane = '_'.join((drms_process_lane.name, replicate_key))

                # Lane-specific file paths

                file_path_dict_lane = dict(
                    temporary_directory=prefix_lane + '_temporary',
                    # TODO: The name for the aligned BAM is constructed by the bsf_run_bwa.py script.
                    # It is currently based on the drms_align_lane.name and replicate_key.
                    # The script should also be changed to pre-set all file names beforehand.
                    aligned_bam='{}_{}.bam'.format(drms_align_lane.name, replicate_key),
                    aligned_bai='{}_{}.bai'.format(drms_align_lane.name, replicate_key),
                    aligned_md5='{}_{}.bam.md5'.format(drms_align_lane.name, replicate_key),
                    duplicates_marked_bam=prefix_lane + '_duplicates_marked.bam',
                    duplicates_marked_bai=prefix_lane + '_duplicates_marked.bai',
                    duplicates_marked_md5=prefix_lane + '_duplicates_marked.bam.md5',
                    duplicate_metrics=prefix_lane + '_duplicate_metrics.tsv',
                    realigner_targets=prefix_lane + '_realigner.interval_list',
                    realigned_bam=prefix_lane + '_realigned.bam',
                    realigned_bai=prefix_lane + '_realigned.bai',
                    recalibration_table_pre=prefix_lane + '_recalibration_pre.table',
                    recalibration_table_post=prefix_lane + '_recalibration_post.table',
                    recalibration_plot=prefix_lane + '_recalibration_report.pdf',
                    recalibrated_bam=prefix_lane + '_recalibrated.bam',
                    recalibrated_bai=prefix_lane + '_recalibrated.bai',
                    alignment_summary_metrics=prefix_lane + '_alignment_summary_metrics.tsv')

                # Create a Runnable and Executable for processing each lane.

                runnable_process_lane = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_lane,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        file_path_dict=file_path_dict_lane,
                        debug=self.debug))
                executable_process_lane = drms_process_lane.add_executable(
                    executable=Executable.from_analysis_runnable(
                        analysis=self,
                        runnable_name=runnable_process_lane.name))
                executable_process_lane.dependencies.append(run_bwa.name)

                # Run the Picard MarkDuplicates analysis, unless configured to skip it.

                if self.skip_mark_duplicates:
                    file_path_dict_lane['duplicates_marked_bam'] = file_path_dict_lane['aligned_bam']
                    file_path_dict_lane['duplicates_marked_bai'] = file_path_dict_lane['aligned_bai']
                    file_path_dict_lane['duplicates_marked_md5'] = file_path_dict_lane['aligned_md5']
                else:
                    runnable_step = runnable_process_lane.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='process_lane_picard_mark_duplicates',
                            obsolete_file_path_list=[
                                # It may not be the best idea to remove the aligned BAM file from the previous
                                # lane-specific alignment step here. For the moment, keep pipeline steps independent
                                # from each other.
                                file_path_dict_lane['aligned_bam'],
                                file_path_dict_lane['aligned_bai'],
                                file_path_dict_lane['aligned_md5']
                            ],
                            java_temporary_path=file_path_dict_lane['temporary_directory'],
                            java_heap_maximum='Xmx6G',
                            picard_classpath=self.classpath_picard,
                            picard_command='MarkDuplicates'))
                    assert isinstance(runnable_step, RunnableStepPicard)
                    runnable_step.add_picard_option(key='INPUT', value=file_path_dict_lane['aligned_bam'])
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict_lane['duplicates_marked_bam'])
                    runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_dict_lane['duplicate_metrics'])
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
                    runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict_lane['temporary_directory'])
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

                # Run the GATK RealignerTargetCreator analysis as the first-pass walker
                # for the GATK IndelRealigner analysis.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_realigner_target_creator',
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='RealignerTargetCreator')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_realignment:
                    runnable_step.add_gatk_option(key='known', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_lane['duplicates_marked_bam'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_lane['realigner_targets'])

                # Run the GATK IndelRealigner analysis as a second-pass walker after the
                # GATK RealignerTargetCreator analysis.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_indel_realigner',
                        obsolete_file_path_list=[
                            file_path_dict_lane['duplicates_marked_bam'],
                            file_path_dict_lane['duplicates_marked_bai'],
                            file_path_dict_lane['duplicates_marked_md5'],
                        ],
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='IndelRealigner')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_realignment:
                    runnable_step.add_gatk_option(key='knownAlleles', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_lane['duplicates_marked_bam'])
                runnable_step.add_gatk_option(key='targetIntervals', value=file_path_dict_lane['realigner_targets'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_lane['realigned_bam'])

                # Run the GATK BaseRecalibrator analysis as a first-pass walker
                # for the GATK PrintReads analysis.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_base_recalibrator_pre',
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='BaseRecalibrator')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_recalibration:
                    runnable_step.add_gatk_option(key='knownSites', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_lane['realigned_bam'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_lane['recalibration_table_pre'])

                # Run the GATK BaseRecalibrator on-the-fly recalibration analysis to generate plots.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_base_recalibrator_post',
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='BaseRecalibrator')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_recalibration:
                    runnable_step.add_gatk_option(key='knownSites', value=file_path, override=True)
                runnable_step.add_gatk_option(key='BQSR', value=file_path_dict_lane['recalibration_table_pre'])
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_lane['realigned_bam'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_lane['recalibration_table_post'])

                # Run the GATK AnalyzeCovariates analysis to create a recalibration plot.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_analyze_covariates',
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='AnalyzeCovariates')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                runnable_step.add_gatk_option(
                    key='afterReportFile',
                    value=file_path_dict_lane['recalibration_table_post'])
                runnable_step.add_gatk_option(
                    key='beforeReportFile',
                    value=file_path_dict_lane['recalibration_table_pre'])
                runnable_step.add_gatk_option(
                    key='plotsReportFile',
                    value=file_path_dict_lane['recalibration_plot'])
                # runnable_step.add_gatk_option(key='logging_level', value='DEBUG')

                # Run the GATK PrintReads analysis as second-pass walker after the GATK BaseRecalibrator analysis.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_lane_gatk_print_reads',
                        obsolete_file_path_list=[
                            file_path_dict_lane['realigned_bam'],
                            file_path_dict_lane['realigned_bai'],
                        ],
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='PrintReads')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_lane['realigned_bam'])
                runnable_step.add_gatk_option(key='BQSR', value=file_path_dict_lane['recalibration_table_pre'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_lane['recalibrated_bam'])

                # Run the Picard CollectAlignmentSummaryMetrics analysis.

                runnable_step = runnable_process_lane.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='process_lane_picard_collect_alignment_summary_metrics',
                        java_temporary_path=file_path_dict_lane['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='CollectAlignmentSummaryMetrics'))
                assert isinstance(runnable_step, RunnableStepPicard)
                runnable_step.add_picard_option(key='INPUT', value=file_path_dict_lane['recalibrated_bam'])
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict_lane['alignment_summary_metrics'])
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
                runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
                runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict_lane['temporary_directory'])
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                # OPTIONS_FILE

                # Set dependencies for the next stage.
                vc_process_sample_dependencies.append(executable_process_lane.name)
                # Add the result of the variant_calling_process_lane Runnable.
                vc_process_sample_replicates.append((
                    file_path_dict_lane['recalibrated_bam'],
                    file_path_dict_lane['recalibrated_bai'],
                    file_path_dict_lane['alignment_summary_metrics'],
                    file_path_dict_lane['duplicate_metrics']))

            # Step 2: Process per sample.
            #
            #   Picard MergeSamFiles
            #   Picard MarkDuplicates
            #   GATK RealignerTargetCreator
            #   GATK IndelRealigner
            #   Picard CollectAlignmentSummaryMetrics
            #   GATK HaplotypeCaller

            prefix_sample = '_'.join((drms_process_sample.name, sample.name))

            file_path_dict_sample = dict(
                temporary_directory=prefix_sample + '_temporary',
                merged_bam=prefix_sample + '_merged.bam',
                merged_bai=prefix_sample + '_merged.bai',
                merged_md5=prefix_sample + '_merged.bam.md5',
                duplicates_marked_bam=prefix_sample + '_duplicates_marked.bam',
                duplicates_marked_bai=prefix_sample + '_duplicates_marked.bai',
                duplicates_marked_md5=prefix_sample + '_duplicates_marked.bam.md5',
                duplicate_metrics=prefix_sample + '_duplicate_metrics.tsv',
                realigner_targets=prefix_sample + '_realigner.intervals',
                realigned_bam=prefix_sample + '_realigned.bam',
                realigned_bai=prefix_sample + '_realigned.bai',
                alignment_summary_metrics=prefix_sample + '_alignment_summary_metrics.tsv',
                raw_variants_gvcf_vcf=prefix_sample + '_raw_variants.g.vcf.gz',
                raw_variants_gvcf_idx=prefix_sample + '_raw_variants.g.vcf.gz.tbi')

            # Create a Runnable and Executable for processing each Sample object.

            runnable_process_sample = self.add_runnable(
                runnable=Runnable(
                    name=prefix_sample,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_dict=file_path_dict_sample,
                    debug=self.debug))
            executable_process_sample = drms_process_sample.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_process_sample.name))
            executable_process_sample.dependencies.extend(vc_process_sample_dependencies)

            if len(vc_process_sample_replicates) == 1:
                # If there is only one replicate, just rename the BAM and BAI files.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_sample_move_recalibrated_bam',
                        source_path=vc_process_sample_replicates[0][0],
                        target_path=file_path_dict_sample['realigned_bam']))

                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepMove(
                        name='process_sample_move_recalibrated_bai',
                        source_path=vc_process_sample_replicates[0][1],
                        target_path=file_path_dict_sample['realigned_bai']))

                # Link the Picard Alignment Summary Metrics.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='process_sample_link_alignment_metrics',
                        source_path=vc_process_sample_replicates[0][2],
                        target_path=file_path_dict_sample['alignment_summary_metrics']))

                # Link the Picard Duplicate Metrics.
                runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='process_sample_link_duplicate_metrics',
                        source_path=vc_process_sample_replicates[0][3],
                        target_path=file_path_dict_sample['duplicate_metrics']))
            else:
                # Run the Picard MergeSamFiles analysis.

                runnable_step = runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='process_sample_picard_merge_sam_files',
                        java_temporary_path=file_path_dict_sample['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='MergeSamFiles'))
                assert isinstance(runnable_step, RunnableStepPicard)
                for file_path_tuple in vc_process_sample_replicates:
                    runnable_step.add_picard_option(key='INPUT', value=file_path_tuple[0], override=True)
                    runnable_step.obsolete_file_path_list.append(file_path_tuple[0])
                    runnable_step.obsolete_file_path_list.append(file_path_tuple[1])
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict_sample['merged_bam'])
                runnable_step.add_picard_option(key='USE_THREADING', value='true')
                runnable_step.add_picard_option(key='COMMENT', value='Merged from the following files:')
                for file_path_tuple in vc_process_sample_replicates:
                    runnable_step.add_picard_option(key='COMMENT', value=file_path_tuple[0], override=True)
                runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict_sample['temporary_directory'])
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
                # Optical duplicates should already have been flagged in the lane-specific processing step.

                if self.skip_mark_duplicates:
                    file_path_dict_sample['duplicates_marked_bam'] = file_path_dict_sample['merged_bam']
                    file_path_dict_sample['duplicates_marked_bai'] = file_path_dict_sample['merged_bai']
                    file_path_dict_sample['duplicates_marked_md5'] = file_path_dict_sample['merged_md5']
                else:
                    runnable_step = runnable_process_sample.add_runnable_step(
                        runnable_step=RunnableStepPicard(
                            name='process_sample_picard_mark_duplicates',
                            obsolete_file_path_list=[
                                file_path_dict_sample['merged_bam'],
                                file_path_dict_sample['merged_bai'],
                                file_path_dict_sample['merged_md5'],
                            ],
                            java_temporary_path=file_path_dict_sample['temporary_directory'],
                            java_heap_maximum='Xmx6G',
                            picard_classpath=self.classpath_picard,
                            picard_command='MarkDuplicates'))
                    assert isinstance(runnable_step, RunnableStepPicard)
                    runnable_step.add_picard_option(key='INPUT', value=file_path_dict_sample['merged_bam'])
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict_sample['duplicates_marked_bam'])
                    runnable_step.add_picard_option(
                        key='METRICS_FILE',
                        value=file_path_dict_sample['duplicate_metrics'])
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
                    runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict_sample['temporary_directory'])
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

                # Run the GATK RealignerTargetCreator analysis as the first-pass walker
                # for the GATK IndelRealigner analysis.

                runnable_step = runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_sample_gatk_realigner_target_creator',
                        java_temporary_path=file_path_dict_sample['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='RealignerTargetCreator')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_realignment:
                    runnable_step.add_gatk_option(key='known', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_sample['duplicates_marked_bam'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_sample['realigner_targets'])

                # Run the GATK IndelRealigner analysis as a second-pass walker
                # after the GATK RealignerTargetCreator analysis.

                runnable_step = runnable_process_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='process_sample_gatk_indel_realigner',
                        obsolete_file_path_list=[
                            file_path_dict_sample['duplicates_marked_bam'],
                            file_path_dict_sample['duplicates_marked_bai'],
                            file_path_dict_sample['duplicates_marked_md5']
                        ],
                        java_temporary_path=file_path_dict_sample['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='IndelRealigner')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                for interval in self.include_intervals_list:
                    runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                if self.interval_padding:
                    runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
                for file_path in self.known_sites_realignment:
                    runnable_step.add_gatk_option(key='knownAlleles', value=file_path, override=True)
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_sample['duplicates_marked_bam'])
                runnable_step.add_gatk_option(key='targetIntervals', value=file_path_dict_sample['realigner_targets'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_sample['realigned_bam'])
                # For debugging only.
                # runnable_step.add_gatk_option(key='logging_level', value='DEBUG')

            # Run the Picard CollectAlignmentSummaryMetrics analysis.

            runnable_step = runnable_process_sample.add_runnable_step(
                runnable_step=RunnableStepPicard(
                    name='process_sample_picard_collect_alignment_summary_metrics',
                    java_temporary_path=file_path_dict_sample['temporary_directory'],
                    java_heap_maximum='Xmx6G',
                    picard_classpath=self.classpath_picard,
                    picard_command='CollectAlignmentSummaryMetrics'))
            assert isinstance(runnable_step, RunnableStepPicard)
            runnable_step.add_picard_option(key='INPUT', value=file_path_dict_sample['realigned_bam'])
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict_sample['alignment_summary_metrics'])
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
            runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
            runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict_sample['temporary_directory'])
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
                    java_temporary_path=file_path_dict_sample['temporary_directory'],
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='HaplotypeCaller')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            if self.downsample_to_fraction:
                runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            # The number of threads should be configurable, but multi-threading seems to cause the occasional problem.
            # runnable_step.add_gatk_option(key='num_cpu_threads_per_data_thread', value='1')
            runnable_step.add_gatk_option(key='pair_hmm_implementation', value='VECTOR_LOGLESS_CACHING')
            runnable_step.add_gatk_option(key='genotyping_mode', value='DISCOVERY')
            runnable_step.add_gatk_option(key='standard_min_confidence_threshold_for_emitting', value='10')
            runnable_step.add_gatk_option(key='standard_min_confidence_threshold_for_calling', value='30')
            runnable_step.add_gatk_option(key='emitRefConfidence', value='GVCF')
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
            runnable_step.add_gatk_option(key='input_file', value=file_path_dict_sample['realigned_bam'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_sample['raw_variants_gvcf_vcf'])
            # Parameter to pass to the VCF/BCF IndexCreator
            runnable_step.add_gatk_option(key='variant_index_type', value='LINEAR')
            runnable_step.add_gatk_option(key='variant_index_parameter', value='128000')

            # Finally, record GVCF files and dependencies for the merge cohort stage.
            # Check if the Sample has annotation for a 'Cohort Name', if not use the cohort name defined for this
            # Analysis in the configuration file.

            if 'Cohort Name' in sample.annotation:
                cohort_key = sample.annotation['Cohort Name'][0]
            else:
                cohort_key = self.cohort_name

            # Record the Executable.name in the cohort dictionary.

            if cohort_key not in vc_merge_cohort_dependency_dict:
                vc_merge_cohort_dependency_dict[cohort_key] = list()
            sample_list = vc_merge_cohort_dependency_dict[cohort_key]
            sample_list.append(executable_process_sample.name)

            # Record the sample-specific GVCF file in the cohort dictionary.

            if cohort_key not in vc_merge_cohort_sample_dict:
                vc_merge_cohort_sample_dict[cohort_key] = list()
            sample_list = vc_merge_cohort_sample_dict[cohort_key]
            assert isinstance(sample_list, list)
            sample_list.append(file_path_dict_sample['raw_variants_gvcf_vcf'])

            # Diagnose the sample

            prefix_diagnosis = '_'.join((drms_diagnose_sample.name, sample.name))

            file_path_dict_diagnosis = dict(
                temporary_directory=prefix_diagnosis + '_temporary',
                realigned_bam=file_path_dict_sample['realigned_bam'],
                realigned_bai=file_path_dict_sample['realigned_bai'],
                diagnose_targets_vcf=prefix_diagnosis + '_diagnose_targets.vcf.gz',
                diagnose_targets_idx=prefix_diagnosis + '_diagnose_targets.vcf.gz.tbi',
                missing_intervals=prefix_diagnosis + '_missing.intervals',
                missing_report=prefix_diagnosis + '_missing.gatkreport',
                callable_bed=prefix_diagnosis + '_callable_loci.bed',
                callable_txt=prefix_diagnosis + '_callable_loci.txt',
                hybrid_selection_metrics=prefix_diagnosis + '_hybrid_selection_metrics.tsv',
            )

            target_interval_name = str()
            target_interval_path = str()
            probe_interval_path = str()

            if 'Target Name' in sample.annotation:
                target_name_list = sample.annotation['Target Name']
                if len(target_name_list) > 1:
                    warnings.warn('More than one set of Target Name annotations is currently not supported.\n'
                                  'Choosing the first one of {!r} for sample {!r}'.
                                  format(target_name_list, sample.name))
                target_interval_name = target_name_list[0]

            if 'Target Intervals' in sample.annotation:
                target_interval_list = sample.annotation['Target Intervals']
                if len(target_interval_list) > 1:
                    warnings.warn('More than one set of Target Interval annotations is currently not supported.\n'
                                  'Choosing the first one of {!r} for sample {!r}'.
                                  format(target_interval_list, sample.name))
                target_interval_path = target_interval_list[0]
                if target_interval_path and not os.path.isabs(target_interval_path):
                    target_interval_path = Default.get_absolute_path(
                        file_path=target_interval_path,
                        default_path=Default.absolute_intervals())

            if 'Probe Intervals' in sample.annotation:
                probe_interval_list = sample.annotation['Probe Intervals']
                if len(probe_interval_list) > 1:
                    warnings.warn('More than one set of Probe Interval annotations is currently not supported.\n'
                                  'Choosing the first one of {!r} for sample {!r}'.
                                  format(probe_interval_list, sample.name))
                probe_interval_path = probe_interval_list[0]
                if probe_interval_path and not os.path.isabs(probe_interval_path):
                    probe_interval_path = Default.get_absolute_path(
                        file_path=probe_interval_path,
                        default_path=Default.absolute_intervals())

            if target_interval_path:
                # Create a Runnable and Executable for diagnosing each sample.

                runnable_diagnose_sample = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_diagnosis,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        file_path_dict=file_path_dict_diagnosis,
                        debug=self.debug))
                executable_diagnose_sample = drms_diagnose_sample.add_executable(
                    executable=Executable.from_analysis_runnable(
                        analysis=self,
                        runnable_name=runnable_diagnose_sample.name))
                executable_diagnose_sample.dependencies.extend(vc_process_sample_dependencies)

                # Run the GATK DiagnoseTarget analysis per sample.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='diagnose_sample_gatk_diagnose_target',
                        java_temporary_path=file_path_dict_diagnosis['temporary_directory'],
                        java_heap_maximum='Xmx8G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='DiagnoseTargets')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                runnable_step.add_gatk_option(key='intervals', value=target_interval_path)
                # if self.interval_padding:
                #     sub_command.add_option_long(key='interval_padding', value=str(self.interval_padding))
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_diagnosis['realigned_bam'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_diagnosis['diagnose_targets_vcf'])
                runnable_step.add_gatk_option(
                    key='missing_intervals',
                    value=file_path_dict_diagnosis['missing_intervals'])

                # Run the GATK QualifyMissingIntervals analysis per sample.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='diagnose_sample_gatk_qualify_missing_intervals',
                        java_temporary_path=file_path_dict_diagnosis['temporary_directory'],
                        java_heap_maximum='Xmx8G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='QualifyMissingIntervals')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                runnable_step.add_gatk_option(key='intervals', value=file_path_dict_diagnosis['missing_intervals'])
                # if self.interval_padding:
                #     sub_command.add_option_long(key='interval_padding', value=str(self.interval_padding))
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_diagnosis['realigned_bam'])
                runnable_step.add_gatk_option(key='targetsfile', value=target_interval_path)
                # TODO: Add baitsfile option if available?
                runnable_step.add_gatk_option(key='out', value=file_path_dict_diagnosis['missing_report'])

                # Run the GATK CallableLoci analysis per sample.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepGATK(
                        name='diagnose_sample_gatk_callable_loci',
                        java_temporary_path=file_path_dict_diagnosis['temporary_directory'],
                        java_heap_maximum='Xmx8G',
                        gatk_classpath=self.classpath_gatk))
                assert isinstance(runnable_step, RunnableStepGATK)
                runnable_step.add_gatk_option(key='analysis_type', value='CallableLoci')
                runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    runnable_step.add_gatk_option(key='downsample_to_fraction', value=self.downsample_to_fraction)
                # for interval in self.exclude_intervals_list:
                #     runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
                # for interval in self.include_intervals_list:
                #     runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
                runnable_step.add_gatk_option(key='intervals', value=target_interval_path)
                # if self.interval_padding:
                #     sub_command.add_option_long(key='interval_padding', value=str(self.interval_padding))
                runnable_step.add_gatk_option(key='input_file', value=file_path_dict_diagnosis['realigned_bam'])
                runnable_step.add_gatk_option(key='out', value=file_path_dict_diagnosis['callable_bed'])
                runnable_step.add_gatk_option(key='summary', value=file_path_dict_diagnosis['callable_txt'])

                # Run the Picard CalculateHsMetrics analysis per sample.

                runnable_step = runnable_diagnose_sample.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='diagnose_sample_picard_calculate_hybrid_selection_metrics',
                        java_temporary_path=file_path_dict_diagnosis['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        picard_classpath=self.classpath_picard,
                        picard_command='CalculateHsMetrics'))
                assert isinstance(runnable_step, RunnableStepPicard)
                if probe_interval_path:
                    runnable_step.add_picard_option(key='BAIT_INTERVALS', value=probe_interval_path)
                else:
                    runnable_step.add_picard_option(key='BAIT_INTERVALS', value=target_interval_path)
                if target_interval_name:
                    runnable_step.add_picard_option(key='BAIT_SET_NAME', value=target_interval_name)
                runnable_step.add_picard_option(key='TARGET_INTERVALS', value=target_interval_path)
                runnable_step.add_picard_option(key='INPUT', value=file_path_dict_diagnosis['realigned_bam'])
                runnable_step.add_picard_option(
                    key='OUTPUT',
                    value=file_path_dict_diagnosis['hybrid_selection_metrics'])
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
                runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
                runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
                # PER_TARGET_COVERAGE
                # TMP_DIR
                runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict_diagnosis['temporary_directory'])
                # VERBOSITY defaults to 'INFO'.
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                # QUIET defaults to 'false'.
                # VALIDATION_STRINGENCY defaults to 'STRICT'.
                # COMPRESSION_LEVEL defaults to '5'.
                # MAX_RECORDS_IN_RAM defaults to '500000'.
                # CREATE_INDEX defaults to 'false'.
                # CREATE_MD5_FILE defaults to 'false'.
                # OPTIONS_FILE

        # Step 3: Hierarchically merge cohorts.
        #
        # GATK CombineGVCFs for each cohort and each sample
        # GATK CombineGVCFs for each cohort
        #
        # The cohorts to merge in need to be configurable and it would be essential,
        # How should sample names be checked before the merge to avoid clashes?
        # Should sample annotation sheets be read or can the combined GVCF file be read in
        # to extract the actual sample names?

        vc_merge_cohort_file_path_list = list()
        vc_merge_cohort_dependency_list = list()
        vc_merge_cohort_prefix_list = list()

        # Run the GATK CombineGVCFs analysis for each cohort and Sample defined in this project to build up
        # cohort-specific GVCF files.

        for cohort_key, sample_gvcf_list in vc_merge_cohort_sample_dict.iteritems():
            prefix_merge = '_'.join((drms_merge_cohort.name, cohort_key))

            vc_merge_cohort_prefix_list.append(prefix_merge)

            file_path_dict_merge = dict(
                temporary_directory=prefix_merge + '_temporary',
                combined_gvcf_vcf=prefix_merge + '_combined.g.vcf.gz',
                combined_gvcf_tbi=prefix_merge + '_combined.g.vcf.gz.tbi',
            )

            # Create a Runnable and Executable for merging each sub-cohort from its Sample objects.

            runnable_merge_cohort = self.add_runnable(
                runnable=Runnable(
                    name=prefix_merge,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_dict=file_path_dict_merge,
                    debug=self.debug))
            executable_merge_cohort = drms_merge_cohort.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_merge_cohort.name))
            executable_merge_cohort.dependencies.extend(vc_merge_cohort_dependency_dict[cohort_key])

            runnable_step = runnable_merge_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='merge_cohort_gatk_combine_gvcfs',
                    java_temporary_path=file_path_dict_merge['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='CombineGVCFs')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            for file_path in sample_gvcf_list:
                runnable_step.add_gatk_option(key='variant', value=file_path, override=True)
            runnable_step.add_gatk_option(key='out', value=file_path_dict_merge['combined_gvcf_vcf'])

            # Finally record dependencies for the process cohort stage.
            vc_merge_cohort_file_path_list.append(file_path_dict_merge['combined_gvcf_vcf'])
            vc_merge_cohort_dependency_list.append(executable_merge_cohort.name)

        # Run the GATK CombineGVCFs analysis once more to merge all cohort-specific GVCF files defined in this project.

        if len(vc_merge_cohort_sample_dict) > 1:
            prefix_merge = '_'.join((drms_merge_cohort.name, self.cohort_name))

            vc_merge_cohort_prefix_list.append(prefix_merge)

            file_path_dict_merge = dict(
                temporary_directory=prefix_merge + '_temporary',
                combined_gvcf_vcf=prefix_merge + '_combined.g.vcf.gz',
                combined_gvcf_tbi=prefix_merge + '_combined.g.vcf.gz.tbi',
            )

            # Create a Runnable and Executable for merging the final cohort from its sub-cohorts.

            runnable_merge_cohort = self.add_runnable(
                runnable=Runnable(
                    name=prefix_merge,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_dict=file_path_dict_merge,
                    debug=self.debug))
            executable_merge_cohort = drms_merge_cohort.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_merge_cohort.name))
            executable_merge_cohort.dependencies.extend(vc_merge_cohort_dependency_list)

            runnable_step = runnable_merge_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='merge_cohort_gatk_combine_gvcfs',
                    java_temporary_path=file_path_dict_merge['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='CombineGVCFs')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            for file_path in vc_merge_cohort_file_path_list:
                runnable_step.add_gatk_option(key='variant', value=file_path, override=True)
            runnable_step.add_gatk_option(key='out', value=file_path_dict_merge['combined_gvcf_vcf'])

            vc_merge_cohort_final_dependency = executable_merge_cohort.name
        else:
            vc_merge_cohort_final_dependency = vc_merge_cohort_dependency_list[0]

        # Step 4: Process per cohort.
        #
        #   GATK CombineGVCFs
        #   GATK GenotypeGVCFs
        #   GATK VariantRecalibrator for SNPs
        #   GATK VariantRecalibrator for INDELs
        #   GATK ApplyRecalibration for SNPs
        #   GATK ApplyRecalibration for INDELs

        prefix_cohort = '_'.join((drms_process_cohort.name, self.cohort_name))

        file_path_dict_cohort = dict(
            temporary_directory=prefix_cohort + '_temporary',
            # Combined GVCF file for the cohort defined in this project.
            combined_gvcf_vcf=vc_merge_cohort_prefix_list[-1] + '_combined.g.vcf.gz',
            combined_gvcf_idx=vc_merge_cohort_prefix_list[-1] + '_combined.g.vcf.gz.tbi',
            # Temporary GVCF file with other cohorts merged in to facilitate recalibration.
            temporary_gvcf_vcf=prefix_cohort + '_temporary.g.vcf.gz',
            temporary_gvcf_idx=prefix_cohort + '_temporary.g.vcf.gz.tbi',
            genotyped_raw_vcf=prefix_cohort + '_genotyped_raw_snp_raw_indel.vcf.gz',
            genotyped_raw_idx=prefix_cohort + '_genotyped_raw_snp_raw_indel.vcf.gz.tbi',
            recalibrated_snp_raw_indel_vcf=prefix_cohort + '_recalibrated_snp_raw_indel.vcf.gz',
            recalibrated_snp_raw_indel_idx=prefix_cohort + '_recalibrated_snp_raw_indel.vcf.gz.tbi',
            recalibrated_snp_recalibrated_indel_vcf=prefix_cohort + '_recalibrated_snp_recalibrated_indel.vcf.gz',
            recalibrated_snp_recalibrated_indel_idx=prefix_cohort + '_recalibrated_snp_recalibrated_indel.vcf.gz.tbi',
            multi_sample_vcf=prefix_cohort + '_multi_sample.vcf.gz',
            multi_sample_idx=prefix_cohort + '_multi_sample.vcf.gz.tbi',
            snpeff_vcf=prefix_cohort + '_snpeff.vcf',
            snpeff_idx=prefix_cohort + '_snpeff.vcf.idx',
            snpeff_stats=prefix_cohort + '_snpeff_summary.html',
            annotated_vcf=prefix_cohort + '_annotated.vcf.gz',
            annotated_idx=prefix_cohort + '_annotated.vcf.gz.tbi',
            recalibration_indel=prefix_cohort + '_recalibration_indel.recal',
            recalibration_snp=prefix_cohort + '_recalibration_snp.recal',
            tranches_indel=prefix_cohort + '_recalibration_indel.tranches',
            tranches_snp=prefix_cohort + '_recalibration_snp.tranches',
            plots_indel=prefix_cohort + '_recalibration_indel.R',
            plots_snp=prefix_cohort + '_recalibration_snp.R')

        # Create a Runnable and Executable for processing the cohort.

        runnable_process_cohort = self.add_runnable(
            runnable=Runnable(
                name=prefix_cohort,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                file_path_dict=file_path_dict_cohort,
                debug=self.debug))
        executable_process_cohort = drms_process_cohort.add_executable(
            executable=Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_process_cohort.name))
        executable_process_cohort.dependencies.append(vc_merge_cohort_final_dependency)

        # Run an additional GATK CombineGVCFs analysis to merge into a super-cohort.

        if len(self.accessory_cohort_gvcfs):
            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_combine_gvcfs_accessory',
                    java_temporary_path=file_path_dict_cohort['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='CombineGVCFs')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            for file_path in self.accessory_cohort_gvcfs:
                runnable_step.add_gatk_option(key='variant', value=file_path, override=True)
            runnable_step.add_gatk_option(
                key='variant',
                value=file_path_dict_cohort['combined_gvcf_vcf'],
                override=True)
            runnable_step.add_gatk_option(key='out', value=file_path_dict_cohort['temporary_gvcf_vcf'])

        # Run the GATK GenotypeGVCFs analysis.

        runnable_step = runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStepGATK(
                name='process_cohort_gatk_genotype_gvcfs',
                java_temporary_path=file_path_dict_cohort['temporary_directory'],
                java_heap_maximum='Xmx6G',
                gatk_classpath=self.classpath_gatk))
        assert isinstance(runnable_step, RunnableStepGATK)
        runnable_step.add_gatk_option(key='analysis_type', value='GenotypeGVCFs')
        runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
        for interval in self.exclude_intervals_list:
            runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
        for interval in self.include_intervals_list:
            runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
        if self.interval_padding:
            runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
        if self.known_sites_discovery:
            runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
        if len(self.accessory_cohort_gvcfs):
            runnable_step.add_gatk_option(key='variant', value=file_path_dict_cohort['temporary_gvcf_vcf'])
        else:
            runnable_step.add_gatk_option(key='variant', value=file_path_dict_cohort['combined_gvcf_vcf'])
        runnable_step.add_gatk_option(key='out', value=file_path_dict_cohort['genotyped_raw_vcf'])

        # Run the VQSR procedure on SNPs.

        if self.vqsr_skip_snp:
            file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'] = file_path_dict_cohort['genotyped_raw_vcf']
            file_path_dict_cohort['recalibrated_snp_raw_indel_idx'] = file_path_dict_cohort['genotyped_raw_idx']
        else:

            # Run the GATK VariantRecalibrator analysis on SNPs.

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_variant_recalibrator_snp',
                    java_temporary_path=file_path_dict_cohort['temporary_directory'],
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='mode', value='SNP')
            for resource in self.vqsr_resources_snp_dict.keys():
                resource_option = 'resource:{},known={},training={},truth={},prior={}'. \
                    format(resource,
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
            runnable_step.add_gatk_option(key='input', value=file_path_dict_cohort['genotyped_raw_vcf'])
            runnable_step.add_gatk_option(key='recal_file', value=file_path_dict_cohort['recalibration_snp'])
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_dict_cohort['tranches_snp'])
            runnable_step.add_gatk_option(key='rscript_file', value=file_path_dict_cohort['plots_snp'])
            if self.vqsr_bad_lod_cutoff_snp is not None:
                runnable_step.add_gatk_option(key='badLodCutoff', value=str(self.vqsr_bad_lod_cutoff_snp))

            # Run the GATK ApplyRecalibration analysis on SNPs.

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_apply_recalibration_snp',
                    java_temporary_path=file_path_dict_cohort['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='ApplyRecalibration')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='mode', value='SNP')
            runnable_step.add_gatk_option(key='input', value=file_path_dict_cohort['genotyped_raw_vcf'])
            runnable_step.add_gatk_option(key='recal_file', value=file_path_dict_cohort['recalibration_snp'])
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_dict_cohort['tranches_snp'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'])
            # The lodCutoff (VQSLOD score) filter is not applied for the moment.
            if self.truth_sensitivity_filter_level_snp:
                runnable_step.add_gatk_option(key='ts_filter_level', value=self.truth_sensitivity_filter_level_snp)

        # Run the VQSR procedure on INDELs.

        if self.vqsr_skip_indel:
            file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf'] = \
                file_path_dict_cohort['recalibrated_snp_raw_indel_vcf']
            file_path_dict_cohort['recalibrated_snp_recalibrated_indel_idx'] = \
                file_path_dict_cohort['recalibrated_snp_raw_indel_idx']
        else:

            # Run the GATK VariantRecalibrator analysis on INDELs.

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_variant_recalibrator_indel',
                    java_temporary_path=file_path_dict_cohort['temporary_directory'],
                    java_heap_maximum='Xmx8G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='VariantRecalibrator')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='mode', value='INDEL')
            for resource in self.vqsr_resources_indel_dict.keys():
                resource_option = 'resource:{},known={},training={},truth={},prior={}'. \
                    format(resource,
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
            runnable_step.add_gatk_option(key='input', value=file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'])
            runnable_step.add_gatk_option(key='recal_file', value=file_path_dict_cohort['recalibration_indel'])
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_dict_cohort['tranches_indel'])
            runnable_step.add_gatk_option(key='rscript_file', value=file_path_dict_cohort['plots_indel'])
            if self.vqsr_bad_lod_cutoff_indel is not None:
                runnable_step.add_gatk_option(key='badLodCutoff', value=str(self.vqsr_bad_lod_cutoff_indel))

            # Run the GATK ApplyRecalibration analysis on INDELs.

            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_apply_recalibration_indel',
                    java_temporary_path=file_path_dict_cohort['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='ApplyRecalibration')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='mode', value='INDEL')
            runnable_step.add_gatk_option(key='input', value=file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'])
            runnable_step.add_gatk_option(key='recal_file', value=file_path_dict_cohort['recalibration_indel'])
            runnable_step.add_gatk_option(key='tranches_file', value=file_path_dict_cohort['tranches_indel'])
            runnable_step.add_gatk_option(
                key='out',
                value=file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf'])
            # The lodCutoff (VQSLOD score) filter is not applied for the moment.
            if self.truth_sensitivity_filter_level_indel:
                runnable_step.add_gatk_option(key='ts_filter_level', value=self.truth_sensitivity_filter_level_indel)

        # In case accessory GVCF files have been used, re-create a multi-sample VCF file with just the samples
        # in this cohort.

        if len(self.accessory_cohort_gvcfs):
            runnable_step = runnable_process_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='process_cohort_gatk_select_variants_cohort',
                    java_temporary_path=file_path_dict_cohort['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))

            runnable_step.add_gatk_option(
                key='variant',
                value=file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_cohort['multi_sample_vcf'])
            for sample in self.samples:
                runnable_step.add_gatk_option(key='sample_name', value=sample.name, override=True)
            runnable_step.add_gatk_switch(key='excludeNonVariants')
        else:
            file_path_dict_cohort['multi_sample_vcf'] = \
                file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf']
            file_path_dict_cohort['multi_sample_idx'] = \
                file_path_dict_cohort['recalibrated_snp_recalibrated_indel_idx']

        # Run the snpEff tool for functional variant annotation.

        java_process = runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStep(
                name='process_cohort_snpeff',
                program='java',
                sub_command=Command(program='eff')))

        java_process.add_switch_short(
            key='d64')
        java_process.add_option_short(
            key='jar',
            value=os.path.join(self.classpath_snpeff, 'snpEff.jar'))
        java_process.add_switch_short(
            key='Xmx6G')
        java_process.add_option_pair(
            key='-Djava.io.tmpdir',
            value=file_path_dict_cohort['temporary_directory'])
        java_process.stdout_path = file_path_dict_cohort['snpeff_vcf']

        sub_command = java_process.sub_command
        sub_command.add_switch_short(key='download')
        sub_command.add_option_short(key='o', value='gatk')
        sub_command.add_option_short(key='stats', value=file_path_dict_cohort['snpeff_stats'])
        sub_command.add_option_short(key='config', value=os.path.join(self.classpath_snpeff, 'snpEff.config'))

        sub_command.arguments.append(self.snpeff_genome_version)
        sub_command.arguments.append(file_path_dict_cohort['multi_sample_vcf'])

        # Run the GATK VariantAnnotator analysis.

        runnable_step = runnable_process_cohort.add_runnable_step(
            runnable_step=RunnableStepGATK(
                name='process_cohort_gatk_variant_annotator',
                java_temporary_path=file_path_dict_cohort['temporary_directory'],
                java_heap_maximum='Xmx4G',
                gatk_classpath=self.classpath_gatk))
        assert isinstance(runnable_step, RunnableStepGATK)
        runnable_step.add_gatk_option(key='analysis_type', value='VariantAnnotator')
        runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
        for interval in self.exclude_intervals_list:
            runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
        for interval in self.include_intervals_list:
            runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
        if self.interval_padding:
            runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
        if self.known_sites_discovery:
            runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)

        # Add annotation resources and their corresponding expression options.
        for annotation_resource in self.annotation_resources_dict.keys():
            assert isinstance(annotation_resource, str)
            if len(self.annotation_resources_dict[annotation_resource][0]) \
                    and len(self.annotation_resources_dict[annotation_resource][1]):
                runnable_step.add_gatk_option(
                    key=':'.join(('resource', annotation_resource)),
                    value=self.annotation_resources_dict[annotation_resource][0])
                for annotation in self.annotation_resources_dict[annotation_resource][1]:
                    assert isinstance(annotation, str)
                    runnable_step.add_gatk_option(
                        key='expression',
                        value='.'.join((annotation_resource, annotation)),
                        override=True)

        runnable_step.add_gatk_option(key='variant', value=file_path_dict_cohort['multi_sample_vcf'])
        # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
        # runnable_step.add_gatk_option(key='annotation', value='AlleleBalanceBySample')
        runnable_step.add_gatk_option(key='annotation', value='SnpEff')
        runnable_step.add_gatk_option(key='snpEffFile', value=file_path_dict_cohort['snpeff_vcf'])
        runnable_step.add_gatk_option(key='out', value=file_path_dict_cohort['annotated_vcf'])

        # Re-process and split the cohort by sample.

        for sample in self.samples:

            prefix_split = '_'.join((drms_split_cohort.name, sample.name))

            file_path_dict_split = dict(
                temporary_directory=prefix_split + '_temporary',
                sample_vcf=prefix_split + '.vcf.gz',
                sample_idx=prefix_split + '.vcf.gz.tbi',
                sample_tsv=prefix_split + '.tsv')

            runnable_split_cohort = self.add_runnable(
                runnable=Runnable(
                    name=prefix_split,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_dict=file_path_dict_split,
                    debug=self.debug))
            executable_split_cohort = drms_split_cohort.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_split_cohort.name))
            executable_split_cohort.dependencies.append(executable_process_cohort.name)

            # Run the GATK SelectVariants analysis to split multi-sample VCF files into one per sample.

            runnable_step = runnable_split_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_gatk_select_variants',
                    java_temporary_path=file_path_dict_split['temporary_directory'],
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='SelectVariants')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))

            runnable_step.add_gatk_option(key='variant', value=file_path_dict_cohort['annotated_vcf'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_split['sample_vcf'])
            runnable_step.add_gatk_option(key='sample_name', value=sample.name)
            runnable_step.add_gatk_switch(key='excludeNonVariants')

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_split_cohort.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='split_cohort_gatk_variants_to_table',
                    java_temporary_path=file_path_dict_split['temporary_directory'],
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='variant', value=file_path_dict_split['sample_vcf'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_split['sample_tsv'])
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
            runnable_step.add_gatk_option(key='fields', value='VQSLOD', override=True)
            runnable_step.add_gatk_option(key='fields', value='AF', override=True)
            # GATK Haplotype Caller genotype fields: GT:AD:DP:GQ:PL
            runnable_step.add_gatk_option(key='genotypeFields', value='GT', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='AD', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='DP', override=True)
            runnable_step.add_gatk_option(key='genotypeFields', value='GQ', override=True)
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
                assert isinstance(annotation_resource, str)
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        assert isinstance(annotation, str)
                        runnable_step.add_gatk_option(
                            key='fields',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

        # Somatic variant calling.

        key_list = self.comparisons.keys()
        key_list.sort(cmp=lambda x, y: cmp(x, y))

        if self.debug > 0:
            print "Somatic variant calling: {!r}".format(key_list)

        for key in key_list:

            # The list of samples must contain exactly one normal and one tumor sample.
            if len(self.comparisons[key]) != 2:
                continue

            prefix_somatic = '_'.join((drms_somatic.name, key))

            # Somatic variant calling-specific file paths

            file_path_dict_somatic = dict(
                temporary_directory=prefix_somatic + '_temporary',
                somatic_vcf=prefix_somatic + '_somatic.vcf.gz',
                somatic_idx=prefix_somatic + '_somatic.vcf.gz.tbi',
                snpeff_vcf=prefix_somatic + '_snpeff.vcf',
                snpeff_idx=prefix_somatic + '_snpeff.vcf.idx',
                snpeff_stats=prefix_somatic + '_snpeff_summary.html',
                snpeff_genes=prefix_somatic + '_snpeff_summary.genes.txt',
                annotated_vcf=prefix_somatic + '_annotated.vcf.gz',
                annotated_idx=prefix_somatic + '_annotated.vcf.gz.tbi',
                annotated_tsv=prefix_somatic + '_annotated.tsv')

            runnable_somatic = self.add_runnable(
                runnable=Runnable(
                    name=prefix_somatic,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_dict=file_path_dict_somatic,
                    debug=self.debug))
            executable_somatic = drms_somatic.add_executable(
                executable=Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_somatic.name))
            executable_somatic.dependencies.append(
                    'variant_calling_process_sample_' + self.comparisons[key][0][1][0].name)
            executable_somatic.dependencies.append(
                    'variant_calling_process_sample_' + self.comparisons[key][-1][1][0].name)

            # Run the GATK MuTect2 analysis to characterise somatic variants.

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_mutect2',
                    java_temporary_path=file_path_dict_somatic['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='MuTect2')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)
            for file_path in self.known_somatic_discovery:
                runnable_step.add_gatk_option(key='cosmic', value=file_path, override=True)
            runnable_step.add_gatk_option(
                key='input_file:normal',
                value='variant_calling_process_sample_{}_realigned.bam'.format(self.comparisons[key][0][1][0].name))
            runnable_step.add_gatk_option(
                key='input_file:tumor',
                value='variant_calling_process_sample_{}_realigned.bam'.format(self.comparisons[key][-1][1][0].name))
            runnable_step.add_gatk_option(key='out', value=file_path_dict_somatic['somatic_vcf'])

            # Run the snpEff tool for functional variant annotation.

            java_process = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStep(
                    name='somatic_snpeff',
                    program='java',
                    sub_command=Command(program='eff')))

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_snpeff, 'snpEff.jar'))
            java_process.add_switch_short(
                key='Xmx6G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_somatic['temporary_directory'])
            java_process.stdout_path = file_path_dict_somatic['snpeff_vcf']

            sub_command = java_process.sub_command
            sub_command.add_switch_short(key='download')
            sub_command.add_option_short(key='o', value='gatk')
            sub_command.add_option_short(key='stats', value=file_path_dict_somatic['snpeff_stats'])
            sub_command.add_option_short(key='config', value=os.path.join(self.classpath_snpeff, 'snpEff.config'))

            sub_command.arguments.append(self.snpeff_genome_version)
            sub_command.arguments.append(file_path_dict_somatic['somatic_vcf'])

            # Run the GATK VariantAnnotator analysis.

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_variant_annotator',
                    java_temporary_path=file_path_dict_somatic['temporary_directory'],
                    java_heap_maximum='Xmx4G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='VariantAnnotator')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            if self.known_sites_discovery:
                runnable_step.add_gatk_option(key='dbsnp', value=self.known_sites_discovery)

            # Add annotation resources and their corresponding expression options.
            for annotation_resource in self.annotation_resources_dict.keys():
                assert isinstance(annotation_resource, str)
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    runnable_step.add_gatk_option(
                        key=':'.join(('resource', annotation_resource)),
                        value=self.annotation_resources_dict[annotation_resource][0])
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        assert isinstance(annotation, str)
                        runnable_step.add_gatk_option(
                            key='expression',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

            runnable_step.add_gatk_option(key='variant', value=file_path_dict_somatic['somatic_vcf'])
            # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
            # runnable_step.add_gatk_option(key='annotation', value='AlleleBalanceBySample')
            runnable_step.add_gatk_option(key='annotation', value='SnpEff')
            runnable_step.add_gatk_option(key='snpEffFile', value=file_path_dict_somatic['snpeff_vcf'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_somatic['annotated_vcf'])

            # Run the GATK VariantsToTable analysis.

            runnable_step = runnable_somatic.add_runnable_step(
                runnable_step=RunnableStepGATK(
                    name='somatic_gatk_variants_to_table',
                    java_temporary_path=file_path_dict_somatic['temporary_directory'],
                    java_heap_maximum='Xmx2G',
                    gatk_classpath=self.classpath_gatk))
            assert isinstance(runnable_step, RunnableStepGATK)
            runnable_step.add_gatk_option(key='analysis_type', value='VariantsToTable')
            runnable_step.add_gatk_option(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                runnable_step.add_gatk_option(key='excludeIntervals', value=interval, override=True)
            for interval in self.include_intervals_list:
                runnable_step.add_gatk_option(key='intervals', value=interval, override=True)
            if self.interval_padding:
                runnable_step.add_gatk_option(key='interval_padding', value=str(self.interval_padding))
            runnable_step.add_gatk_option(key='variant', value=file_path_dict_somatic['annotated_vcf'])
            runnable_step.add_gatk_option(key='out', value=file_path_dict_somatic['annotated_tsv'])
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
                assert isinstance(annotation_resource, str)
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        assert isinstance(annotation, str)
                        runnable_step.add_gatk_option(
                            key='fields',
                            value='.'.join((annotation_resource, annotation)),
                            override=True)

        return

    def report(self):
        """Create a C{VariantCallingGATK} report in HTML format and a UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        # This code only needs the public URL.

        track_output = str()

        # Write a HTML document.

        output = str()

        output += defaults.web.html_header(title='{} Variant Calling Analysis'.format(self.project_name))
        output += '<body>\n'
        output += '\n'

        output += '<h1>{} Variant Calling Analysis</h1>\n'.format(self.project_name)
        output += '\n'

        output += '<h2>UCSC Alignment Track Hub</h2>\n'
        output += '\n'

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/variant_calling_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output += '<p>UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.</p>\n'.format(
            defaults.web.ucsc_track_url(
                options_dict=options_dict,
                host_name=default.ucsc_host_name),
            self.project_name)
        output += '\n'

        output += '<h2>Sample and Aliquot Level</h2>\n'
        output += '\n'
        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr>\n'
        output += '<th>Sample</th>\n'
        output += '<th>Variants VCF file</th>\n'
        output += '<th>Variants TSV file</th>\n'
        output += '<th>Aligned BAM file</th>\n'
        output += '<th>Aligned BAI file</th>\n'
        output += '<th>Aliquot</th>\n'
        output += '<th>Duplicate Metrics</th>\n'
        output += '<th>Alignment Summary Metrics</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        # Group via UCSC super tracks.

        track_output += 'track Alignments\n'
        track_output += 'shortLabel Alignments\n'
        track_output += 'longLabel BWA NGS read alignments\n'
        track_output += 'superTrack on show\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Variants\n'
        track_output += 'shortLabel Variants\n'
        track_output += 'longLabel Variant calls\n'
        track_output += 'superTrack on show\n'
        track_output += 'group variants\n'
        track_output += '\n'

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            runnable_process_sample = self.runnable_dict[
                '_'.join((self.drms_name_process_sample, sample.name))]
            assert isinstance(runnable_process_sample, Runnable)
            runnable_split_cohort = self.runnable_dict[
                '_'.join((self.drms_name_split_cohort, sample.name))]
            assert isinstance(runnable_split_cohort, Runnable)

            # Alignment track
            # Common settings

            track_output += 'track {}_alignments\n'. \
                format(sample.name)
            track_output += 'type bam\n'
            track_output += 'shortLabel {}_alignments\n'. \
                format(sample.name)
            track_output += 'longLabel {} BWA NGS read alignments\n'. \
                format(sample.name)
            track_output += 'bigDataUrl ./{}\n'. \
                format(runnable_process_sample.file_path_dict['realigned_bam'])
            # track_output += 'html {}\n'.format()
            track_output += 'visibility squish\n'

            # Common optional settings.

            track_output += 'color {}\n'. \
                format('0,0,0')

            # Compressed Sequence Alignment track settings.

            # None so far.

            # Composite track settings.

            track_output += 'parent Alignments\n'
            track_output += '\n'

            # Variant track

            track_output += 'track {}_variants\n'.format(sample.name)
            track_output += 'type vcfTabix\n'
            track_output += 'shortLabel {}_variants\n'.format(sample.name)
            track_output += 'longLabel {} variant calls\n'.format(sample.name)
            track_output += 'bigDataUrl ./{}\n'.format(runnable_split_cohort.file_path_dict['sample_vcf'])
            # track_output += 'html {}\n'.format()
            track_output += 'visibility dense\n'

            # vcfTabix specific settings

            # None so far.

            # Composite track settings.

            track_output += 'parent Variants\n'
            track_output += '\n'

            output += '<tr>\n'
            output += '<td>{}</td>\n'. \
                format(sample.name)
            output += '<td><a href="{}">VCF</a></td>\n'. \
                format(runnable_split_cohort.file_path_dict['sample_vcf'])
            output += '<td><a href="{}">TSV</a></td>\n'. \
                format(runnable_split_cohort.file_path_dict['sample_tsv'])
            output += '<td><a href="{}">BAM</a></td>\n'. \
                format(runnable_process_sample.file_path_dict['realigned_bam'])
            output += '<td><a href="{}">BAI</a></td>\n'. \
                format(runnable_process_sample.file_path_dict['realigned_bai'])
            output += '<td></td>\n'  # Aliquot
            output += '<td><a href="{}">TSV</a></td>\n'. \
                format(runnable_process_sample.file_path_dict['duplicate_metrics'])
            output += '<td><a href="{}">TSV</a></td>\n'. \
                format(runnable_process_sample.file_path_dict['alignment_summary_metrics'])
            output += '</tr>\n'

            # bsf.data.Sample.get_all_paired_reads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of bsf.data.PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                runnable_process_lane = self.runnable_dict[
                    '_'.join((self.drms_name_process_lane, replicate_key))]
                assert isinstance(runnable_process_lane, Runnable)
                output += '<tr>\n'
                output += '<td></td>\n'  # Sample
                output += '<td></td>\n'  # Variants VCF
                output += '<td></td>\n'  # Variants TSV
                output += '<td></td>\n'  # Aligned BAM
                output += '<td></td>\n'  # Aligned BAI
                output += '<td>{}</td>\n'.format(replicate_key)
                output += '<td><a href="{}">TSV</a></td>\n'. \
                    format(runnable_process_lane.file_path_dict['duplicate_metrics'])
                output += '<td><a href="{}">TSV</a></td>\n'. \
                    format(runnable_process_lane.file_path_dict['alignment_summary_metrics'])
                output += '</tr>\n'

        output += '</tbody>\n'
        output += '</table>\n'
        output += '\n'

        output += '<h2>Cohort Level</h2>\n'
        output += '\n'
        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr>\n'
        output += '<th>Cohort</th>\n'
        output += '<th>Information</th>\n'
        output += '<th>Comment</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        runnable_process_cohort = self.runnable_dict[
            '_'.join((self.drms_name_process_cohort, self.cohort_name))]
        assert isinstance(runnable_process_cohort, Runnable)

        output += '<tr>\n'
        output += '<td>{}</td>\n'. \
            format(self.cohort_name)
        output += '<td><a href="{}">snpEff Summary Statistics</a></td>\n'. \
            format(runnable_process_cohort.file_path_dict['snpeff_stats'])
        output += '<td></td>\n'
        output += '</tr>\n'

        output += '<tr>\n'
        output += '<td>{}</td>\n'. \
            format(self.cohort_name)
        output += '<td>snpEff-annotated multi-sample <a href="{}">VCF</a></td>\n'. \
            format(runnable_process_cohort.file_path_dict['snpeff_vcf'])
        output += '<td>Functional annotation of all splice variants</td>\n'
        output += '</tr>\n'

        output += '<tr>\n'
        output += '<td>{}</td>\n'. \
            format(self.cohort_name)
        output += '<td>GATK-annotated multi-sample <a href="{}">VCF</a></td>\n'. \
            format(runnable_process_cohort.file_path_dict['annotated_vcf'])
        output += '<td>Functional annotation of only the most severely affected splice variant</td>\n'
        output += '</tr>\n'

        output += '</tbody>\n'
        output += '</table>\n'
        output += '\n'

        # Somatic variant calling.

        key_list = self.comparisons.keys()
        key_list.sort(cmp=lambda x, y: cmp(x, y))

        if len(key_list):
            output += '<h2>Somatic Variants</h2>\n'
            output += '\n'
            output += '<table>\n'
            output += '<thead>\n'
            output += '<tr>\n'
            output += '<th>Comparison</th>\n'
            output += '<th>Annotated VCF</th>\n'
            output += '<th>Annotated TSV</th>\n'
            output += '<th>snpEff Summary Statistics</th>\n'
            output += '<th>snpEff Genes</th>\n'
            output += '</tr>\n'
            output += '</thead>\n'
            output += '<tbody>\n'

            for key in key_list:
                # The list of samples must contain exactly one normal and one tumor sample.
                if len(self.comparisons[key]) != 2:
                    continue

                runnable_somatic = self.runnable_dict['_'.join((self.drms_name_somatic, key))]
                assert isinstance(runnable_somatic, Runnable)
                output += '<tr>\n'
                output += '<td>{}</td>\n'.format(key)
                output += '<td><a href="{}">VCF</a></td>\n'.format(runnable_somatic.file_path_dict['annotated_vcf'])
                output += '<td><a href="{}">TSV</a></td>\n'.format(runnable_somatic.file_path_dict['annotated_tsv'])
                output += '<td><a href="{}">HTML</a></td>\n'.format(runnable_somatic.file_path_dict['snpeff_stats'])
                output += '<td><a href="{}">TXT</a></td>\n'.format(runnable_somatic.file_path_dict['snpeff_genes'])
                output += '</tr>\n'

            output += '</tbody>\n'
            output += '</table>\n'
            output += '\n'

        output += '<h2>QC Plots</h2>\n'
        output += '\n'
        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr><th>Sample</th><th>Aliquot</th><th>Metrics</th></tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        # Alignment Summary - Percent Aligned
        if os.path.exists(os.path.join(
                self.genome_directory, 'variant_calling_summary_alignment_percentage_sample.png')):
            output += '<tr>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_alignment_percentage_sample.pdf">'
            output += '<img alt="Alignment Summary - Percent Aligned per Sample" ' \
                      'src="variant_calling_summary_alignment_percentage_sample.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_alignment_percentage_read_group.pdf">'
            output += '<img alt="Alignment Summary - Percent Aligned per Aliquot" ' \
                      'src="variant_calling_summary_alignment_percentage_read_group.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>Alignment Summary - Percent Aligned</td>\n'
            output += '</tr>\n'

        # Alignment Summary - Reads Aligned
        if os.path.exists(os.path.join(
                self.genome_directory, 'variant_calling_summary_alignment_reads_sample.png')):
            output += '<tr>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_alignment_reads_sample.pdf">'
            output += '<img alt="Alignment Summary - Reads Aligned per Sample" ' \
                      'src="variant_calling_summary_alignment_reads_sample.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_alignment_reads_read_group.pdf">'
            output += '<img alt="Alignment Summary - Reads Aligned per Aliquot" ' \
                      'src="variant_calling_summary_alignment_reads_read_group.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>Alignment Summary - Reads Aligned</td>\n'
            output += '</tr>\n'

        # Hybrid Selection - Mean Target Coverage
        if os.path.exists(os.path.join(
                self.genome_directory, 'variant_calling_summary_hybrid_target_coverage_sample.png')):
            output += '<tr>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_hybrid_target_coverage_sample.pdf">'
            output += '<img alt="Hybrid Selection - Mean Target Coverage per Sample" ' \
                      'src="variant_calling_summary_hybrid_target_coverage_sample.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_hybrid_target_coverage_read_group.pdf">'
            output += '<img alt="Hybrid Selection - Mean Target Coverage per Aliquot" ' \
                      'src="variant_calling_summary_hybrid_target_coverage_read_group.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>Hybrid Selection - Mean Target Coverage</td>\n'
            output += '</tr>\n'

        # Hybrid Selection - Percent Unique Reads
        if os.path.exists(os.path.join(
                self.genome_directory, 'variant_calling_summary_hybrid_unique_percentage_sample.png')):
            output += '<tr>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_hybrid_unique_percentage_sample.pdf">'
            output += '<img alt="Hybrid Selection - Percent Unique Reads per Sample" ' \
                      'src="variant_calling_summary_hybrid_unique_percentage_sample.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>'
            output += '<a href="variant_calling_summary_hybrid_unique_percentage_read_group.pdf">'
            output += '<img alt="Hybrid Selection - Percent Unique Reads per Aliquot" ' \
                      'src="variant_calling_summary_hybrid_unique_percentage_read_group.png" ' \
                      'height="100" width="100" />'
            output += '</a>'
            output += '</td>\n'
            output += '<td>Hybrid Selection - Percent Unique Reads</td>\n'
            output += '</tr>\n'

        output += '</tbody>\n'
        output += '</table>\n'
        output += '\n'
        output += '</body>\n'
        output += defaults.web.html_footer()

        file_path = os.path.join(self.genome_directory, 'variant_calling_report.html')

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        # Create the UCSC Genome Browser Track Hub.

        self.ucsc_hub_write_hub(prefix='variant_calling')
        self.ucsc_hub_write_genomes(prefix='variant_calling')
        self.ucsc_hub_write_tracks(output=track_output, prefix='variant_calling')

        return
