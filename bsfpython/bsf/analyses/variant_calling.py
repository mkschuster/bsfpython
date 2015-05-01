"""bsf.analyses.variant_calling

A package of classes and methods supporting variant calling analyses.
"""

#
# Copyright 2014 Michael K. Schuster
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
import string
import warnings

from bsf import Analysis, Command, Configuration, Default, defaults, DRMS, Executable, Runnable, RunnableStep
from bsf.annotation import SampleAnnotationSheet
from bsf.data import PairedReads
from bsf.executables import BWA


class VariantCallingGATK(Analysis):
    """The C{VariantCallingGATK} class represents the logic to run the Genome Analysis Toolkit (GATK).

    Attributes:
    @cvar drms_name_align_lane: C{DRMS.name} for the lane alignment C{Analysis} stage
    @type drms_name_align_lane: str
    @cvar drms_name_process_lane: C{DRMS.name} for the lane processing C{Analysis} stage
    @type drms_name_process_lane: str
    @cvar drms_name_process_sample: C{DRMS.name} for the sample processing C{Analysis} stage
    @type drms_name_process_sample: str
    @cvar drms_name_process_cohort: C{DRMS.name} for the cohort processing C{Analysis} stage
    @type drms_name_process_cohort: str
    @cvar drms_name_split_cohort: C{DRMS.name} for the cohort splitting C{Analysis} stage
    @type drms_name_split_cohort: str
    @ivar replicate_grouping: Group individual C{PairedReads} objects for processing or run them separately
    @type replicate_grouping: bool
    @ivar comparison_path: Comparison file
    @type comparison_path: str | unicode
    @ivar cohort_name: Cohort name
    @type cohort_name: str
    @ivar accessory_cohort_gvcfs: Python C{list} of Python C{str} | C{unicode} (GVCF file path) objects
    @type accessory_cohort_gvcfs: list
    @ivar skip_mark_duplicates: Skip the Picard MarkDuplicates step
    @type skip_mark_duplicates: bool
    @ivar known_sites_discovery: VCF file path for variant discovery via The Haplotype Caller or Unified Genotyper
    @type known_sites_discovery: str | unicode
    @ivar known_sites_realignment: Python C{list} of Python C{str} | C{unicode} (VCF file paths)
        for realignment
    @type known_sites_realignment:list
    @ivar known_sites_recalibration: Python C{list} of Python C{str} | C{unicode} (VCF file paths)
        for recalibration
    @type known_sites_recalibration: list
    @ivar annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
        Python C{tuple} of
        Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
    @type annotation_resources_dict: dict
    @ivar truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
    @type truth_sensitivity_filter_level_indel: str
    @ivar truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
    @type truth_sensitivity_filter_level_snp: str
    @ivar vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
    @type vqsr_skip_indel: bool, None
    @ivar vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
    @type vqsr_skip_snp: bool, None
    @ivar vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_indel_list: list
    @ivar vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
    @type vqsr_annotations_snp_list: list
    @ivar exclude_intervals_list: Python C{list} of Python C{str} (intervals) to exclude from the analysis
    @type exclude_intervals_list: list
    @ivar include_intervals_list: Python C{list} of Python C{str} (intervals) to include in the analysis
    @type include_intervals_list: list
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
    drms_name_process_cohort = 'variant_calling_process_cohort'
    drms_name_split_cohort = 'variant_calling_split_cohort'

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{VariantCallingGATK} object from a UNIX-style configuration file via the
        C{Configuration} class.

        @param config_path: UNIX-style configuration file
        @type config_path: str | unicode
        @return: C{VariantCallingGATK}
        @rtype: VariantCallingGATK
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{VariantCallingGATK} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{VariantCallingGATK}
        @rtype: VariantCallingGATK
        """

        assert isinstance(configuration, Configuration)

        variant_calling = cls(configuration=configuration)

        # A "bsf.analyses.variant_calling.VariantCallingGATK" section specifies defaults for this Analysis.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        variant_calling.set_configuration(variant_calling.configuration, section=section)

        return variant_calling

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 replicate_grouping=False, bwa_genome_db=None, comparison_path=None,
                 cohort_name=None, accessory_cohort_gvcfs=None, skip_mark_duplicates=False,
                 known_sites_discovery=None, known_sites_realignment=None, known_sites_recalibration=None,
                 annotation_resources_dict=None,
                 truth_sensitivity_filter_level_indel=None,
                 truth_sensitivity_filter_level_snp=None,
                 vqsr_skip_indel=None, vqsr_skip_snp=None,
                 vqsr_resources_indel_dict=None, vqsr_resources_snp_dict=None,
                 vqsr_annotations_indel_list=None, vqsr_annotations_snp_list=None,
                 exclude_intervals_list=None,
                 include_intervals_list=None,
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
        @type drms_list: list
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{list} objects of C{Sample} objects
        @type comparisons: dict
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        @param replicate_grouping: Group individual C{PairedReads} objects for processing or run them separately
        @type replicate_grouping: bool
        @param bwa_genome_db: Genome sequence file path with BWA index
        @type bwa_genome_db: str | unicode
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode
        @param cohort_name: Cohort name
        @type cohort_name: str
        @param accessory_cohort_gvcfs: Python C{list} of Python C{str} | C{unicode} (GVCF file path) objects
        @type accessory_cohort_gvcfs: list
        @param skip_mark_duplicates: Skip the Picard MarkDuplicates step
        @type skip_mark_duplicates: bool
        @param known_sites_discovery: VCF file path for variant discovery via The Haplotype Caller or Unified Genotyper
        @type known_sites_discovery: str | unicode
        @param known_sites_realignment: Python C{list} of Python C{str} | C{unicode} (VCF file paths)
            for realignment
        @type known_sites_realignment:list
        @param known_sites_recalibration: Python C{list} of Python C{str} | C{unicode} (VCF file paths)
            for recalibration
        @type known_sites_recalibration: list
        @param annotation_resources_dict: Python C{dict} of Python C{str} (annotation resource name) key and
            Python C{tuple} of
            Python C{str} (file path) and Python C{list} of Python C{str} (annotation) value data
        @type annotation_resources_dict: dict
        @param truth_sensitivity_filter_level_indel: Truth sensitivity filter level for INDELs
        @type truth_sensitivity_filter_level_indel: str
        @param truth_sensitivity_filter_level_snp: Truth sensitivity filter level for SNPs
        @type truth_sensitivity_filter_level_snp: str
        @param vqsr_skip_indel: Skip the Variant Quality Score Recalibration on INDELs
        @type vqsr_skip_indel: bool, None
        @param vqsr_skip_snp: Skip the Variant Quality Score Recalibration on SNPs
        @type vqsr_skip_snp: bool, None
        @param vqsr_resources_indel_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_indel_dict: dict
        @param vqsr_resources_snp_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_snp_dict: dict
        @param vqsr_annotations_indel_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_indel_list: list
        @param vqsr_annotations_snp_list: Python C{list} of Python C{str} (variant annotation) objects
        @type vqsr_annotations_snp_list: list
        @param exclude_intervals_list: Python C{list} of Python C{str} (intervals) to exclude from the analysis
        @type exclude_intervals_list: list
        @param include_intervals_list: Python C{list} of Python C{str} (intervals) to include in the analysis
        @type include_intervals_list: list
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
        """

        super(VariantCallingGATK, self).__init__(
            configuration=configuration,
            project_name=project_name, genome_version=genome_version,
            input_directory=input_directory, output_directory=output_directory,
            project_directory=project_directory, genome_directory=genome_directory,
            e_mail=e_mail, debug=debug, drms_list=drms_list,
            collection=collection, comparisons=comparisons, samples=samples)

        # Sub-class specific ...

        self.replicate_grouping = replicate_grouping

        if bwa_genome_db:
            self.bwa_genome_db = bwa_genome_db
        else:
            self.bwa_genome_db = str()

        if comparison_path:
            self.comparison_path = comparison_path
        else:
            self.comparison_path = str()

        if cohort_name:
            self.cohort_name = cohort_name
        else:
            self.cohort_name = str()

        if accessory_cohort_gvcfs:
            self.accessory_cohort_gvcfs = accessory_cohort_gvcfs
        else:
            self.accessory_cohort_gvcfs = list()

        self.skip_mark_duplicates = skip_mark_duplicates

        if known_sites_discovery:
            self.known_sites_discovery = known_sites_discovery
        else:
            self.known_sites_discovery = str()

        if known_sites_realignment:
            self.known_sites_realignment = known_sites_realignment
        else:
            self.known_sites_realignment = list()

        if known_sites_recalibration:
            self.known_sites_recalibration = known_sites_recalibration
        else:
            self.known_sites_recalibration = list()

        if annotation_resources_dict:
            self.annotation_resources_dict = annotation_resources_dict
        else:
            self.annotation_resources_dict = dict()

        if truth_sensitivity_filter_level_indel:
            self.truth_sensitivity_filter_level_indel = truth_sensitivity_filter_level_indel
        else:
            self.truth_sensitivity_filter_level_indel = str()

        if truth_sensitivity_filter_level_snp:
            self.truth_sensitivity_filter_level_snp = truth_sensitivity_filter_level_snp
        else:
            self.truth_sensitivity_filter_level_snp = str()

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

        if vqsr_resources_indel_dict:
            self.vqsr_resources_indel_dict = vqsr_resources_indel_dict
        else:
            self.vqsr_resources_indel_dict = dict()

        if vqsr_resources_snp_dict:
            self.vqsr_resources_snp_dict = vqsr_resources_snp_dict
        else:
            self.vqsr_resources_snp_dict = dict()

        if vqsr_annotations_indel_list:
            self.vqsr_annotations_indel_list = vqsr_annotations_indel_list
        else:
            self.vqsr_annotations_indel_list = list()

        if vqsr_annotations_snp_list:
            self.vqsr_annotations_snp_list = vqsr_annotations_snp_list
        else:
            self.vqsr_annotations_snp_list = list()

        if exclude_intervals_list:
            self.exclude_intervals_list = exclude_intervals_list
        else:
            self.exclude_intervals_list = list()

        if include_intervals_list:
            self.include_intervals_list = include_intervals_list
        else:
            self.include_intervals_list = list()

        if downsample_to_fraction:
            self.downsample_to_fraction = downsample_to_fraction
        else:
            self.downsample_to_fraction = str()

        if gatk_bundle_version:
            self.gatk_bundle_version = gatk_bundle_version
        else:
            self.gatk_bundle_version = str()

        if snpeff_genome_version:
            self.snpeff_genome_version = snpeff_genome_version
        else:
            self.snpeff_genome_version = str()

        if classpath_gatk:
            self.classpath_gatk = classpath_gatk
        else:
            self.classpath_gatk = str()

        if classpath_picard:
            self.classpath_picard = classpath_picard
        else:
            self.classpath_picard = str()

        if classpath_snpeff:
            self.classpath_snpeff = classpath_snpeff
        else:
            self.classpath_snpeff = str()

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
        """

        super(VariantCallingGATK, self).set_configuration(configuration=configuration, section=section)

        if configuration.config_parser.has_option(section=section, option='replicate_grouping'):
            self.replicate_grouping = configuration.config_parser.getboolean(
                section=section,
                option='replicate_grouping')

        # Get the genome database.

        if configuration.config_parser.has_option(section=section, option='bwa_genome_db'):
            self.bwa_genome_db = configuration.config_parser.get(
                section=section,
                option='bwa_genome_db')

        # Read a comparison file.

        # if configuration.config_parser.has_option(section=section, option='comparison_path'):
        #     self.comparison_path = configuration.config_parser.get(section=section, option='comparison_path')
        # Use the sample annotation sheet instead of a separate comparison file.
        if configuration.config_parser.has_option(section=section, option='sas_file'):
            self.comparison_path = configuration.config_parser.get(
                section=section,
                option='sas_file')

        # Get the cohort name.

        if configuration.config_parser.has_option(section=section, option='cohort_name'):
            self.cohort_name = configuration.config_parser.get(
                section=section,
                option='cohort_name')

        # Comma-separated list of GVCF files from accessory cohorts
        # that should be used in the recalibration procedure.

        if configuration.config_parser.has_option(section=section, option='accessory_cohort_gvcfs'):
            for file_path in configuration.config_parser.get(
                    section=section,
                    option='accessory_cohort_gvcfs').split(','):
                self.accessory_cohort_gvcfs.append(file_path.strip())  # Strip white space around commas.

        # Get the skip mark duplicates option.

        if configuration.config_parser.has_option(section=section, option='skip_mark_duplicates'):
            self.skip_mark_duplicates = configuration.config_parser.getboolean(
                section=section,
                option='skip_mark_duplicates')

        # Get the truth sensitivity filter level for INDELs.

        if configuration.config_parser.has_option(section=section, option='truth_sensitivity_filter_level_indel'):
            self.truth_sensitivity_filter_level_indel = configuration.config_parser.get(
                section=section,
                option='truth_sensitivity_filter_level_indel')

        # Get the truth sensitivity filter level for SNPs.

        if configuration.config_parser.has_option(section=section, option='truth_sensitivity_filter_level_snp'):
            self.truth_sensitivity_filter_level_snp = configuration.config_parser.get(
                section=section,
                option='truth_sensitivity_filter_level_snp')

        # Get the flag for skipping the Variant Quality Score Recalibration (VQSR) for INDELs.

        if configuration.config_parser.has_option(section=section, option='vqsr_skip_indel'):
            self.vqsr_skip_indel = configuration.config_parser.getboolean(
                section=section,
                option='vqsr_skip_indel')

        # Get the flag for skipping the Variant Quality Score Recalibration (VQSR) for SNPs.

        if configuration.config_parser.has_option(section=section, option='vqsr_skip_snp'):
            self.vqsr_skip_snp = configuration.config_parser.getboolean(
                section=section,
                option='vqsr_skip_snp')

        # Get the list of annotations for the Variant Quality Score Recalibration (VQSR) for INDELs.

        if configuration.config_parser.has_option(section=section, option='vqsr_annotations_indel'):
            for annotation in configuration.config_parser.get(
                    section=section,
                    option='vqsr_annotations_indel').split(','):
                self.vqsr_annotations_indel_list.append(annotation.strip())  # Strip white space around commas.

        # Get the list of annotations for the Variant Quality Score Recalibration (VQSR) for SNPs.

        if configuration.config_parser.has_option(section=section, option='vqsr_annotations_snp'):
            for annotation in configuration.config_parser.get(
                    section=section,
                    option='vqsr_annotations_snp').split(','):
                self.vqsr_annotations_snp_list.append(annotation.strip())  # Strip white space around commas.

        self._read_vqsr_configuration(vqsr_resources_dict=self.vqsr_resources_indel_dict, variation_type='indel')
        self._read_vqsr_configuration(vqsr_resources_dict=self.vqsr_resources_snp_dict, variation_type='snp')

        # Read additionally requested annotation resources for the GATK AnnotateVariants step.

        # Python dict of Python str (annotation resource name) key and
        # Python tuple of
        # Python str (file path) and Python list of Python str (annotation) value data.

        if configuration.config_parser.has_option(section=section, option='annotation_resources'):
            for annotation_resource in configuration.config_parser.get(
                    section=section,
                    option='annotation_resources').split(','):
                annotation_resource = annotation_resource.strip()  # Strip white space around commas.
                resource_section = string.join(words=(annotation_resource, 'resource'), sep='_')
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

        if configuration.config_parser.has_option(section=section, option='known_sites_discovery'):
            self.known_sites_discovery = configuration.config_parser.get(
                section=section,
                option='known_sites_discovery')

        # Comma-separated list of VCF files with known variant sites for the
        # GATK RealignerTargetCreator and IndelRealigner steps.

        if configuration.config_parser.has_option(section=section, option='known_sites_realignment'):
            for file_path in configuration.config_parser.get(
                    section=section,
                    option='known_sites_realignment').split(','):
                self.known_sites_realignment.append(file_path.strip())  # Strip white space around commas.

        # Comma-separated list of VCF files with known variant sites for the
        # GATK BaseRecalibrator and PrintReads steps.

        if configuration.config_parser.has_option(section=section, option='known_sites_recalibration'):
            for file_path in configuration.config_parser.get(
                    section=section,
                    option='known_sites_recalibration').split(','):
                self.known_sites_recalibration.append(file_path.strip())  # Strip white space around commas.

        # Get the list of intervals to exclude.

        if configuration.config_parser.has_option(section=section, option='exclude_intervals'):
            exclude_intervals = configuration.config_parser.get(
                section=section,
                option='exclude_intervals')
            # For comma-separated interval lists split into components on commas, strip white space
            # and push them onto the list individually.
            self.exclude_intervals_list.extend(map(lambda x: x.strip(), exclude_intervals.split(',')))

        # Get the list of intervals to include.

        if configuration.config_parser.has_option(section=section, option='include_intervals'):
            include_intervals = configuration.config_parser.get(
                section=section,
                option='include_intervals')
            # For comma-separated interval lists split into components on commas, strip white space and
            # push them onto the list individually.
            self.include_intervals_list.extend(map(lambda x: x.strip(), include_intervals.split(',')))

        # Get the down-sample to fraction information.

        if configuration.config_parser.has_option(section=section, option='downsample_to_fraction'):
            self.downsample_to_fraction = configuration.config_parser.get(
                section=section,
                option='downsample_to_fraction')

        # Get the GATK bundle version.

        if configuration.config_parser.has_option(section=section, option='gatk_bundle_version'):
            self.gatk_bundle_version = configuration.config_parser.get(
                section=section,
                option='gatk_bundle_version')

        # Get the snpEff genome version.

        if configuration.config_parser.has_option(section=section, option='snpeff_genome_version'):
            self.snpeff_genome_version = configuration.config_parser.get(
                section=section,
                option='snpeff_genome_version')

        # Get the Genome Analysis Tool Kit (GATK) Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_gatk'):
            self.classpath_gatk = configuration.config_parser.get(
                section=section,
                option='classpath_gatk')

        # Get the Picard tools Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_picard'):
            self.classpath_picard = configuration.config_parser.get(
                section=section,
                option='classpath_picard')

        # Get the snpEff tool Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_snpeff'):
            self.classpath_snpeff = configuration.config_parser.get(
                section=section,
                option='classpath_snpeff')

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
        """

        sas = SampleAnnotationSheet.from_file_path(file_path=comparison_path)

        for row_dict in sas.row_dicts:
            self.add_sample(sample=self.collection.get_sample_from_row_dict(row_dict=row_dict))

    def _read_vqsr_configuration(self, vqsr_resources_dict, variation_type=None):
        """Private method to read variant quality score recalibration (VQSR) configuration information.

        @param vqsr_resources_dict: Python C{dict} of Python C{str} (resource name) and Python C{dict} values
        @type vqsr_resources_dict: dict
        @param variation_type: Variation type I{indel} or I{snp}
        @type variation_type: str
        """

        if variation_type not in ('indel', 'snp'):
            raise Exception("Variation type has to be 'indel' or 'snp', not {!r}.".format(variation_type))

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        resource_option = string.join(words=('vqsr_resources', variation_type), sep='_')
        if config_parser.has_option(section=config_section, option=resource_option):
            for resource in config_parser.get(section=config_section, option=resource_option).split(','):
                resource = resource.strip()
                resource_section = string.join(words=('vqsr', variation_type, resource), sep='_')
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

    def run(self):
        """Run this C{VariantCallingGATK} analysis.
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
            if interval[-10:] == '.intervals' or interval[-14:] == '.interval_list':
                # For Picard-style interval lists prepend the current directory if necessary.
                if not os.path.isabs(interval):
                    interval = os.path.join(os.path.realpath(os.path.curdir), interval)
                if os.path.exists(interval):
                    temporary_list.append(interval)
                else:
                    raise Exception('Exclude intervals file {!r} does not exist.'.format(interval))
        self.exclude_intervals_list = temporary_list

        # Include intervals

        temporary_list = list()
        for interval in self.include_intervals_list:
            if interval[-10:] == '.intervals' or interval[-14:] == '.interval_list':
                # For Picard-style interval lists prepend the current directory if necessary.
                if not os.path.isabs(interval):
                    interval = os.path.join(os.path.realpath(os.path.curdir), interval)
                if os.path.exists(interval):
                    temporary_list.append(interval)
                else:
                    raise Exception('Include intervals file {!r} does not exist.'.format(interval))
        self.include_intervals_list = temporary_list

        # The comparison file can be absolute, in the same directory or in the project directory.
        if not os.path.isabs(self.comparison_path) and not os.path.exists(self.comparison_path):
            self.comparison_path = Default.get_absolute_path(
                file_path=self.comparison_path,
                default_path=self.project_directory)

        # Real comparisons would be required for somatic mutation calling.
        self._read_comparisons(comparison_path=self.comparison_path)

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # bsf_run_bwa.py script.

        drms_align_lane = DRMS.from_analysis(
            name=self.drms_name_align_lane,
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(drms_align_lane)

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_process_lane Runnable.

        drms_process_lane = DRMS.from_analysis(
            name=self.drms_name_process_lane,
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(drms_process_lane)

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_process_sample Runnable.

        drms_process_sample = DRMS.from_analysis(
            name=self.drms_name_process_sample,
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(drms_process_sample)

        # Initialise a Distributed Resource Management System (DRMS) object for the
        # variant_calling_process_cohort Runnable.

        drms_process_cohort = DRMS.from_analysis(
            name=self.drms_name_process_cohort,
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(drms_process_cohort)

        drms_split_cohort = DRMS.from_analysis(
            name=self.drms_name_split_cohort,
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(drms_split_cohort)

        vc_process_cohort_dependencies = list()
        vc_process_cohort_replicates = list()

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
                    bwa_mem.arguments.append(string.join(words=reads1, sep=','))
                elif len(reads1) and len(reads2):
                    bwa_mem.arguments.append(string.join(words=reads1, sep=','))
                    bwa_mem.arguments.append(string.join(words=reads2, sep=','))
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

                run_bwa = Executable.from_analysis(
                    name=string.join(words=(drms_align_lane.name, replicate_key), sep='_'),
                    program='bsf_run_bwa.py',
                    analysis=self)
                drms_align_lane.add_executable(executable=run_bwa)

                # Only submit this Executable if the final result file does not exist.
                if (os.path.exists(
                        os.path.join(self.genome_directory, file_path_align_lane['aligned_md5']))
                    and os.path.getsize(
                        os.path.join(self.genome_directory, file_path_align_lane['aligned_md5']))):
                    run_bwa.submit = False
                # Check also for existence of a new-style Runnable status file.
                if os.path.exists(os.path.join(
                        drms_align_lane.work_directory,
                        string.join(words=(drms_align_lane.name, replicate_key, 'completed.txt'), sep='_'))):
                    run_bwa.submit = False

                    # Set run_bwa options.

                run_bwa.add_option_long(key='pickler_path', value=pickler_path)
                run_bwa.add_option_long(key='debug', value=str(self.debug))

                prefix_lane = string.join(words=(drms_process_lane.name, replicate_key), sep='_')

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
                    realigner_targets=prefix_lane + '_realigner.intervals',
                    realigned_bam=prefix_lane + '_realigned.bam',
                    realigned_bai=prefix_lane + '_realigned.bai',
                    recalibration_table_pre=prefix_lane + '_recalibration_pre.table',
                    recalibration_table_post=prefix_lane + '_recalibration_post.table',
                    recalibration_plot=prefix_lane + '_recalibration_report.pdf',
                    recalibrated_bam=prefix_lane + '_recalibrated.bam',
                    recalibrated_bai=prefix_lane + '_recalibrated.bai',
                    alignment_summary_metrics=prefix_lane + '_alignment_summary_metrics.tsv')

                # Lane-specific Runnable

                runnable_process_lane = Runnable(
                    name=prefix_lane,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_dict=file_path_dict_lane,
                    debug=self.debug)
                self.add_runnable(runnable=runnable_process_lane)

                # Run the Picard MarkDuplicates step, unless configured to skip it.

                if self.skip_mark_duplicates:
                    file_path_dict_lane['duplicates_marked_bam'] = file_path_dict_lane['aligned_bam']
                    file_path_dict_lane['duplicates_marked_bai'] = file_path_dict_lane['aligned_bai']
                    file_path_dict_lane['duplicates_marked_md5'] = file_path_dict_lane['aligned_md5']
                else:
                    java_process = RunnableStep(
                        name='picard_mark_duplicates',
                        program='java',
                        sub_command=Command(command=str()),
                        obsolete_file_path_list=[
                            # It may not be the best idea to remove the aligned BAM file from the previous
                            # lane-specific alignment step here. For the moment, keep pipeline steps independent
                            # from each other.
                            file_path_dict_lane['aligned_bam'],
                            file_path_dict_lane['aligned_bai'],
                            file_path_dict_lane['aligned_md5']
                        ])
                    runnable_process_lane.add_runnable_step(runnable_step=java_process)

                    java_process.add_switch_short(
                        key='d64')
                    java_process.add_option_short(
                        key='jar',
                        value=os.path.join(self.classpath_picard, 'MarkDuplicates.jar'))
                    java_process.add_switch_short(
                        key='Xmx6G')
                    java_process.add_option_pair(
                        key='-Djava.io.tmpdir',
                        value=file_path_dict_lane['temporary_directory'])

                    sub_command = java_process.sub_command
                    sub_command.add_option_pair(key='INPUT', value=file_path_dict_lane['aligned_bam'])
                    sub_command.add_option_pair(key='OUTPUT', value=file_path_dict_lane['duplicates_marked_bam'])
                    sub_command.add_option_pair(key='METRICS_FILE', value=file_path_dict_lane['duplicate_metrics'])
                    # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
                    # as otherwise, optical duplicates could not be detected. This is a consequence of using
                    # Illumina2bam rather than Picard ExtractIlluminaBarcodes, IlluminaBasecallsToFastq and
                    # IlluminaBasecallsToSam.
                    # See BioStar post: http://www.biostars.org/p/12538/
                    # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    sub_command.add_option_pair(key='READ_NAME_REGEX',
                                                value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*')
                    sub_command.add_option_pair(key='TMP_DIR', value=file_path_dict_lane['temporary_directory'])
                    sub_command.add_option_pair(key='VERBOSITY', value='WARNING')
                    sub_command.add_option_pair(key='QUIET', value='false')
                    sub_command.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')
                    sub_command.add_option_pair(key='COMPRESSION_LEVEL', value='5')
                    sub_command.add_option_pair(key='MAX_RECORDS_IN_RAM', value='4000000')
                    sub_command.add_option_pair(key='CREATE_INDEX', value='true')
                    sub_command.add_option_pair(key='CREATE_MD5_FILE', value='true')

                # Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

                java_process = RunnableStep(
                    name='gatk_realigner_target_creator',
                    program='java',
                    sub_command=Command(command=str()))
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_long(key='analysis_type', value='RealignerTargetCreator')
                sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    sub_command.add_option_long(key='excludeIntervals', value=interval)
                for interval in self.include_intervals_list:
                    sub_command.add_option_long(key='intervals', value=interval)
                for file_path in self.known_sites_realignment:
                    sub_command.add_option_long(key='known', value=file_path)
                sub_command.add_option_long(key='input_file', value=file_path_dict_lane['duplicates_marked_bam'])
                sub_command.add_option_long(key='out', value=file_path_dict_lane['realigner_targets'])

                # Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

                java_process = RunnableStep(
                    name='gatk_indel_realigner',
                    program='java',
                    sub_command=Command(command=str()),
                    obsolete_file_path_list=[
                        file_path_dict_lane['duplicates_marked_bam'],
                        file_path_dict_lane['duplicates_marked_bai'],
                        file_path_dict_lane['duplicates_marked_md5']
                    ])
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_long(key='analysis_type', value='IndelRealigner')
                sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    sub_command.add_option_long(key='excludeIntervals', value=interval)
                for interval in self.include_intervals_list:
                    sub_command.add_option_long(key='intervals', value=interval)
                for file_path in self.known_sites_realignment:
                    sub_command.add_option_long(key='knownAlleles', value=file_path)
                sub_command.add_option_long(key='input_file', value=file_path_dict_lane['duplicates_marked_bam'])
                sub_command.add_option_long(key='targetIntervals', value=file_path_dict_lane['realigner_targets'])
                sub_command.add_option_long(key='out', value=file_path_dict_lane['realigned_bam'])

                # Run the GATK BaseRecalibrator step as a first-pass walker for the GATK PrintReads step.

                java_process = RunnableStep(
                    name='gatk_base_recalibrator_pre',
                    program='java',
                    sub_command=Command(command=str()))
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_long(key='analysis_type', value='BaseRecalibrator')
                sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    sub_command.add_option_long(key='excludeIntervals', value=interval)
                for interval in self.include_intervals_list:
                    sub_command.add_option_long(key='intervals', value=interval)
                for file_path in self.known_sites_recalibration:
                    sub_command.add_option_long(key='knownSites', value=file_path)
                sub_command.add_option_long(key='input_file', value=file_path_dict_lane['realigned_bam'])
                sub_command.add_option_long(key='out', value=file_path_dict_lane['recalibration_table_pre'])

                # Run the GATK BaseRecalibrator on-the-fly recalibration step to generate plots.

                java_process = RunnableStep(
                    name='gatk_base_recalibrator_post',
                    program='java',
                    sub_command=Command(command=str()))
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_long(key='analysis_type', value='BaseRecalibrator')
                sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    sub_command.add_option_long(key='excludeIntervals', value=interval)
                for interval in self.include_intervals_list:
                    sub_command.add_option_long(key='intervals', value=interval)
                for file_path in self.known_sites_recalibration:
                    sub_command.add_option_long(key='knownSites', value=file_path)
                sub_command.add_option_long(key='BQSR', value=file_path_dict_lane['recalibration_table_pre'])
                sub_command.add_option_long(key='input_file', value=file_path_dict_lane['realigned_bam'])
                sub_command.add_option_long(key='out', value=file_path_dict_lane['recalibration_table_post'])

                # Run the GATK AnalyzeCovariates step to create a recalibration plot.

                java_process = RunnableStep(
                    name='gatk_analyze_covariates',
                    program='java',
                    sub_command=Command(command=str()))
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_long(key='analysis_type', value='AnalyzeCovariates')
                sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    sub_command.add_option_long(key='excludeIntervals', value=interval)
                for interval in self.include_intervals_list:
                    sub_command.add_option_long(key='intervals', value=interval)
                sub_command.add_option_long(key='afterReportFile',
                                            value=file_path_dict_lane['recalibration_table_post'])
                sub_command.add_option_long(key='beforeReportFile',
                                            value=file_path_dict_lane['recalibration_table_pre'])
                sub_command.add_option_long(key='plotsReportFile',
                                            value=file_path_dict_lane['recalibration_plot'])
                # sub_command.add_option_long(key='logging_level', value='DEBUG')

                # Run the GATK PrintReads step as second-pass walker after the BaseRecalibrator step.

                java_process = RunnableStep(
                    name='gatk_print_reads',
                    program='java',
                    sub_command=Command(command=str()),
                    obsolete_file_path_list=[
                        file_path_dict_lane['realigned_bam'],
                        file_path_dict_lane['realigned_bai']
                    ])
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_long(key='analysis_type', value='PrintReads')
                sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
                if self.downsample_to_fraction:
                    sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
                for interval in self.exclude_intervals_list:
                    sub_command.add_option_long(key='excludeIntervals', value=interval)
                for interval in self.include_intervals_list:
                    sub_command.add_option_long(key='intervals', value=interval)
                sub_command.add_option_long(key='input_file', value=file_path_dict_lane['realigned_bam'])
                sub_command.add_option_long(key='BQSR', value=file_path_dict_lane['recalibration_table_pre'])
                sub_command.add_option_long(key='out', value=file_path_dict_lane['recalibrated_bam'])

                # Run the Picard CollectAlignmentSummaryMetrics step.

                java_process = RunnableStep(
                    name='picard_collect_alignment_summary_metrics',
                    program='java',
                    sub_command=Command(command=str()))
                runnable_process_lane.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_picard, 'CollectAlignmentSummaryMetrics.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_lane['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_pair(key='INPUT', value=file_path_dict_lane['recalibrated_bam'])
                sub_command.add_option_pair(key='OUTPUT', value=file_path_dict_lane['alignment_summary_metrics'])
                sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS')
                sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE')
                sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY')
                sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP')
                sub_command.add_option_pair(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
                sub_command.add_option_pair(key='TMP_DIR', value=file_path_dict_lane['temporary_directory'])
                sub_command.add_option_pair(key='VERBOSITY', value='WARNING')
                sub_command.add_option_pair(key='QUIET', value='false')
                sub_command.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')
                sub_command.add_option_pair(key='COMPRESSION_LEVEL', value='5')
                sub_command.add_option_pair(key='MAX_RECORDS_IN_RAM', value='4000000')
                sub_command.add_option_pair(key='CREATE_INDEX', value='true')
                sub_command.add_option_pair(key='CREATE_MD5_FILE', value='true')

                # Create an Executable for the Runnable processing the lane.

                executable_process_lane = Executable.from_analysis_runnable(
                    analysis=self,
                    runnable_name=runnable_process_lane.name)
                drms_process_lane.add_executable(executable=executable_process_lane)

                executable_process_lane.dependencies.append(run_bwa.name)

                # Set dependencies for the next stage.
                vc_process_sample_dependencies.append(executable_process_lane.name)
                # Add the result of the variant_calling_process_lane Runnable.
                vc_process_sample_replicates.append(file_path_dict_lane['recalibrated_bam'])

            # Step 2: Process per sample.
            #
            #   Picard MergeSamFiles
            #   Picard MarkDuplicates
            #   GATK RealignerTargetCreator
            #   GATK IndelRealigner
            #   Picard CollectAlignmentSummaryMetrics
            #   GATK HaplotypeCaller

            prefix_sample = string.join(words=(drms_process_sample.name, sample.name), sep='_')

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
                raw_variants_gvcf_vcf=prefix_sample + '_raw_variants_gvcf.vcf',
                raw_variants_gvcf_idx=prefix_sample + '_raw_variants_gvcf.vcf.idx')

            # Sample-specific Runnable

            runnable_process_sample = Runnable(
                name=prefix_sample,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                file_path_dict=file_path_dict_sample,
                debug=self.debug)
            self.add_runnable(runnable=runnable_process_sample)

            # Run the Picard MergeSamFiles step.

            java_process = RunnableStep(
                name='picard_merge_sam_files',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_sample.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_picard, 'MergeSamFiles.jar'))
            java_process.add_switch_short(
                key='Xmx6G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_sample['temporary_directory'])

            sub_command = java_process.sub_command
            for file_path in vc_process_sample_replicates:
                sub_command.add_option_pair(key='INPUT', value=file_path)
            sub_command.add_option_pair(key='OUTPUT', value=file_path_dict_sample['merged_bam'])
            sub_command.add_option_pair(key='USE_THREADING', value='true')
            sub_command.add_option_pair(key='COMMENT', value='Merged from the following files:')
            for file_path in vc_process_sample_replicates:
                sub_command.add_option_pair(key='COMMENT', value=file_path)
            sub_command.add_option_pair(key='TMP_DIR', value=file_path_dict_sample['temporary_directory'])
            sub_command.add_option_pair(key='VERBOSITY', value='WARNING')
            sub_command.add_option_pair(key='QUIET', value='false')
            sub_command.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')
            sub_command.add_option_pair(key='COMPRESSION_LEVEL', value='5')
            sub_command.add_option_pair(key='MAX_RECORDS_IN_RAM', value='4000000')
            sub_command.add_option_pair(key='CREATE_INDEX', value='true')
            sub_command.add_option_pair(key='CREATE_MD5_FILE', value='true')

            # Run the Picard MarkDuplicates step, unless configured to skip it.
            # Optical duplicates should already have been flagged in the lane-specific processing step.

            if self.skip_mark_duplicates:
                file_path_dict_sample['duplicates_marked_bam'] = file_path_dict_sample['merged_bam']
                file_path_dict_sample['duplicates_marked_bai'] = file_path_dict_sample['merged_bai']
                file_path_dict_sample['duplicates_marked_md5'] = file_path_dict_sample['merged_md5']
            else:
                java_process = RunnableStep(
                    name='picard_mark_duplicates',
                    program='java',
                    sub_command=Command(command=str()),
                    obsolete_file_path_list=[
                        file_path_dict_sample['merged_bam'],
                        file_path_dict_sample['merged_bai'],
                        file_path_dict_sample['merged_md5']
                    ])
                runnable_process_sample.add_runnable_step(runnable_step=java_process)

                java_process.add_switch_short(
                    key='d64')
                java_process.add_option_short(
                    key='jar',
                    value=os.path.join(self.classpath_picard, 'MarkDuplicates.jar'))
                java_process.add_switch_short(
                    key='Xmx6G')
                java_process.add_option_pair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict_sample['temporary_directory'])

                sub_command = java_process.sub_command
                sub_command.add_option_pair(key='INPUT', value=file_path_dict_sample['merged_bam'])
                sub_command.add_option_pair(key='OUTPUT', value=file_path_dict_sample['duplicates_marked_bam'])
                sub_command.add_option_pair(key='METRICS_FILE', value=file_path_dict_sample['duplicate_metrics'])
                # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
                # as otherwise optical duplicates could not be detected. This is a consequence of using Illumina2bam
                # rather than Picard ExtractIlluminaBarcodes, IlluminaBasecallsToFastq and IlluminaBasecallsToSam.
                # See BioStar post: http://www.biostars.org/p/12538/
                # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                sub_command.add_option_pair(key='READ_NAME_REGEX',
                                            value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*')
                sub_command.add_option_pair(key='TMP_DIR', value=file_path_dict_sample['temporary_directory'])
                sub_command.add_option_pair(key='VERBOSITY', value='WARNING')
                sub_command.add_option_pair(key='QUIET', value='false')
                sub_command.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')
                sub_command.add_option_pair(key='COMPRESSION_LEVEL', value='5')
                sub_command.add_option_pair(key='MAX_RECORDS_IN_RAM', value='4000000')
                sub_command.add_option_pair(key='CREATE_INDEX', value='true')
                sub_command.add_option_pair(key='CREATE_MD5_FILE', value='true')

            # Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

            java_process = RunnableStep(
                name='gatk_realigner_target_creator',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_sample.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx6G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_sample['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='RealignerTargetCreator')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            if self.downsample_to_fraction:
                sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            for file_path in self.known_sites_realignment:
                sub_command.add_option_long(key='known', value=file_path)
            sub_command.add_option_long(key='input_file', value=file_path_dict_sample['duplicates_marked_bam'])
            sub_command.add_option_long(key='out', value=file_path_dict_sample['realigner_targets'])

            # Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

            java_process = RunnableStep(
                name='gatk_indel_realigner',
                program='java',
                sub_command=Command(command=str()),
                obsolete_file_path_list=[
                    file_path_dict_sample['duplicates_marked_bam'],
                    file_path_dict_sample['duplicates_marked_bai'],
                    file_path_dict_sample['duplicates_marked_md5']
                ])
            runnable_process_sample.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx6G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_sample['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='IndelRealigner')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            if self.downsample_to_fraction:
                sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            for file_path in self.known_sites_realignment:
                sub_command.add_option_long(key='knownAlleles', value=file_path)
            sub_command.add_option_long(key='input_file', value=file_path_dict_sample['duplicates_marked_bam'])
            sub_command.add_option_long(key='targetIntervals', value=file_path_dict_sample['realigner_targets'])
            sub_command.add_option_long(key='out', value=file_path_dict_sample['realigned_bam'])
            # For debugging only.
            # sub_command.add_option_long(key='logging_level', value='DEBUG')

            # Run the Picard CollectAlignmentSummaryMetrics step.

            java_process = RunnableStep(
                name='picard_collect_alignment_summary_metrics',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_sample.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_picard, 'CollectAlignmentSummaryMetrics.jar'))
            java_process.add_switch_short(
                key='Xmx6G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_sample['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_pair(key='INPUT', value=file_path_dict_sample['realigned_bam'])
            sub_command.add_option_pair(key='OUTPUT', value=file_path_dict_sample['alignment_summary_metrics'])
            sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS')
            sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE')
            sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY')
            sub_command.add_option_pair(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP')
            sub_command.add_option_pair(key='REFERENCE_SEQUENCE', value=self.bwa_genome_db)
            sub_command.add_option_pair(key='TMP_DIR', value=file_path_dict_sample['temporary_directory'])
            sub_command.add_option_pair(key='VERBOSITY', value='WARNING')
            sub_command.add_option_pair(key='QUIET', value='false')
            sub_command.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')
            sub_command.add_option_pair(key='COMPRESSION_LEVEL', value='5')
            sub_command.add_option_pair(key='MAX_RECORDS_IN_RAM', value='4000000')
            sub_command.add_option_pair(key='CREATE_INDEX', value='true')
            sub_command.add_option_pair(key='CREATE_MD5_FILE', value='true')

            # Run the GATK HaplotypeCaller per sample.

            java_process = RunnableStep(
                name='gatk_haplotype_caller',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_sample.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx8G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_sample['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='HaplotypeCaller')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            if self.downsample_to_fraction:
                sub_command.add_option_long(key='downsample_to_fraction', value=self.downsample_to_fraction)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            # TODO: The number of threads should be configurable.
            # sub_command.add_option_long(key='num_cpu_threads_per_data_thread', value='1')
            sub_command.add_option_long(key='pair_hmm_implementation', value='VECTOR_LOGLESS_CACHING')
            sub_command.add_option_long(key='genotyping_mode', value='DISCOVERY')
            sub_command.add_option_long(key='standard_min_confidence_threshold_for_emitting', value='10')
            sub_command.add_option_long(key='standard_min_confidence_threshold_for_calling', value='30')
            sub_command.add_option_long(key='emitRefConfidence', value='GVCF')
            if self.known_sites_discovery:
                sub_command.add_option_long(key='dbsnp', value=self.known_sites_discovery)
            sub_command.add_option_long(key='input_file', value=file_path_dict_sample['realigned_bam'])
            sub_command.add_option_long(key='out', value=file_path_dict_sample['raw_variants_gvcf_vcf'])
            # Parameter to pass to the VCF/BCF IndexCreator
            sub_command.add_option_long(key='variant_index_type', value='LINEAR')
            sub_command.add_option_long(key='variant_index_parameter', value='128000')

            # Create an Executable for processing the sample.

            executable_process_sample = Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_process_sample.name)
            drms_process_sample.add_executable(executable=executable_process_sample)

            executable_process_sample.dependencies.extend(vc_process_sample_dependencies)

            # Record dependencies for the next stage.
            vc_process_cohort_dependencies.append(executable_process_sample.name)
            # Add the result of the variant_calling_process_sample Runnable.
            vc_process_cohort_replicates.append(file_path_dict_sample['raw_variants_gvcf_vcf'])

        # Step 3: Process per cohort.
        #
        #   GATK CombineGVCFs
        #   GATK GenotypeGVCFs
        #   GATK VariantRecalibrator for SNPs
        #   GATK VariantRecalibrator for INDELs
        #   GATK ApplyRecalibration for SNPs
        #   GATK ApplyRecalibration for INDELs

        prefix_cohort = string.join(words=(drms_process_cohort.name, self.cohort_name), sep='_')

        file_path_dict_cohort = dict(
            temporary_directory=prefix_cohort + '_temporary',
            # Combined GVCF file for the cohort defined in this project.
            combined_gvcf_vcf=prefix_cohort + '_combined_gvcf.vcf',
            combined_gvcf_idx=prefix_cohort + '_combined_gvcf.vcf.idx',
            # Temporary GVCF file with other cohorts merged in to facilitate recalibration.
            temporary_gvcf_vcf=prefix_cohort + '_temporary_gvcf.vcf',
            temporary_gvcf_idx=prefix_cohort + '_temporary_gvcf.vcf.idx',
            genotyped_raw_vcf=prefix_cohort + '_genotyped_raw_snp_raw_indel.vcf',
            genotyped_raw_idx=prefix_cohort + '_genotyped_raw_snp_raw_indel.vcf.idx',
            recalibrated_snp_raw_indel_vcf=prefix_cohort + '_recalibrated_snp_raw_indel.vcf',
            recalibrated_snp_raw_indel_idx=prefix_cohort + '_recalibrated_snp_raw_indel.vcf.idx',
            recalibrated_snp_recalibrated_indel_vcf=prefix_cohort + '_recalibrated_snp_recalibrated_indel.vcf',
            recalibrated_snp_recalibrated_indel_idx=prefix_cohort + '_recalibrated_snp_recalibrated_indel.vcf.idx',
            multi_sample_vcf=prefix_cohort + '_multi_sample.vcf',
            multi_sample_idx=prefix_cohort + '_multi_sample.vcf.idx',
            snpeff_vcf=prefix_cohort + '_snpeff.vcf',
            snpeff_idx=prefix_cohort + '_snpeff.vcf.idx',
            snpeff_stats=prefix_cohort + '_snpeff_summary.html',
            annotated_vcf=prefix_cohort + '_annotated.vcf',
            annotated_idx=prefix_cohort + '_annotated.vcf.idx',
            recalibration_indel=prefix_cohort + '_recalibration_indel.recal',
            recalibration_snp=prefix_cohort + '_recalibration_snp.recal',
            tranches_indel=prefix_cohort + '_recalibration_indel.tranches',
            tranches_snp=prefix_cohort + '_recalibration_snp.tranches',
            plots_indel=prefix_cohort + '_recalibration_indel.R',
            plots_snp=prefix_cohort + '_recalibration_snp.R')

        # Cohort-specific Runnable

        runnable_process_cohort = Runnable(
            name=prefix_cohort,
            code_module='bsf.runnables.generic',
            working_directory=self.genome_directory,
            file_path_dict=file_path_dict_cohort,
            debug=self.debug)
        self.add_runnable(runnable=runnable_process_cohort)

        # Hierarchical merging of samples before GenotypeGVCFs processing.
        #
        # First, merge all samples from this analysis project into a project-specific cohort.
        # Ideally, an VariantCallingGATK-specific sub-class of the SampleAnnotationSheet extends it with
        # another column 'cohort', so that samples of the same project could be merged into one or multiple cohorts.
        #
        # Second, allow merging with other cohorts to facilitate recalibration of project with too few samples.
        # This would generate a temporary cohort file on which all recalibration steps are run on.
        #
        # Third, it would be good to run another GATK SelectVariants step to select all samples from this cohort,
        # so that a clean multi-sample VCF file can be delivered.
        #
        # The cohorts to merge in need to be configurable and it would be essential,
        # to check for sample name clashes beforehand. How should this be done?
        # Should sample annotation sheets be read or can the combined GVCF file be read in
        # to extract the actual sample names?

        # Run the GATK CombineGVCFs step for the cohort defined in this project.

        java_process = RunnableStep(
            name='gatk_combine_gvcfs',
            program='java',
            sub_command=Command(command=str()))
        runnable_process_cohort.add_runnable_step(runnable_step=java_process)

        java_process.add_switch_short(
            key='d64')
        java_process.add_option_short(
            key='jar',
            value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_switch_short(
            key='Xmx4G')
        java_process.add_option_pair(
            key='-Djava.io.tmpdir',
            value=file_path_dict_cohort['temporary_directory'])

        sub_command = java_process.sub_command
        sub_command.add_option_long(key='analysis_type', value='CombineGVCFs')
        sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
        for interval in self.exclude_intervals_list:
            sub_command.add_option_long(key='excludeIntervals', value=interval)
        for interval in self.include_intervals_list:
            sub_command.add_option_long(key='intervals', value=interval)
        for file_path in vc_process_cohort_replicates:
            sub_command.add_option_long(key='variant', value=file_path)
        sub_command.add_option_long(key='out', value=file_path_dict_cohort['combined_gvcf_vcf'])

        # Run an additional GATK CombineGVCFs step to merge into a super-cohort.

        if len(self.accessory_cohort_gvcfs):
            java_process = RunnableStep(
                name='gatk_combine_gvcfs_accessory',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx4G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_cohort['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='CombineGVCFs')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            for file_path in self.accessory_cohort_gvcfs:
                sub_command.add_option_long(key='variant', value=file_path)
            sub_command.add_option_long(key='variant', value=file_path_dict_cohort['combined_gvcf_vcf'])
            sub_command.add_option_long(key='out', value=file_path_dict_cohort['temporary_gvcf_vcf'])

        # Run the GATK GenotypeGVCFs step.

        java_process = RunnableStep(
            name='gatk_genotype_gvcfs',
            program='java',
            sub_command=Command(command=str()))
        runnable_process_cohort.add_runnable_step(runnable_step=java_process)

        java_process.add_switch_short(
            key='d64')
        java_process.add_option_short(
            key='jar',
            value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_switch_short(
            key='Xmx6G')
        java_process.add_option_pair(
            key='-Djava.io.tmpdir',
            value=file_path_dict_cohort['temporary_directory'])

        sub_command = java_process.sub_command
        sub_command.add_option_long(key='analysis_type', value='GenotypeGVCFs')
        sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
        for interval in self.exclude_intervals_list:
            sub_command.add_option_long(key='excludeIntervals', value=interval)
        for interval in self.include_intervals_list:
            sub_command.add_option_long(key='intervals', value=interval)
        if self.known_sites_discovery:
            sub_command.add_option_long(key='dbsnp', value=self.known_sites_discovery)
        if len(self.accessory_cohort_gvcfs):
            sub_command.add_option_long(key='variant', value=file_path_dict_cohort['temporary_gvcf_vcf'])
        else:
            sub_command.add_option_long(key='variant', value=file_path_dict_cohort['combined_gvcf_vcf'])
        sub_command.add_option_long(key='out', value=file_path_dict_cohort['genotyped_raw_vcf'])

        # Run the VQSR procedure on SNPs.

        if self.vqsr_skip_snp:
            file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'] = file_path_dict_cohort['genotyped_raw_vcf']
            file_path_dict_cohort['recalibrated_snp_raw_indel_idx'] = file_path_dict_cohort['genotyped_raw_idx']
        else:

            # Run the GATK VariantRecalibrator for SNPs.

            java_process = RunnableStep(
                name='gatk_variant_recalibrator_snp',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx8G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_cohort['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='VariantRecalibrator')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            sub_command.add_option_long(key='mode', value='SNP')
            for resource in self.vqsr_resources_snp_dict.keys():
                resource_option = 'resource:{},known={},training={},truth={},prior={}'. \
                    format(resource,
                           self.vqsr_resources_snp_dict[resource]['known'],
                           self.vqsr_resources_snp_dict[resource]['training'],
                           self.vqsr_resources_snp_dict[resource]['truth'],
                           self.vqsr_resources_snp_dict[resource]['prior'])
                sub_command.add_option_long(
                    key=resource_option,
                    value=self.vqsr_resources_snp_dict[resource]['file_path'])
            for annotation in self.vqsr_annotations_snp_list:
                sub_command.add_option_long(key='use_annotation', value=annotation)
            sub_command.add_option_long(key='input', value=file_path_dict_cohort['genotyped_raw_vcf'])
            sub_command.add_option_long(key='recal_file', value=file_path_dict_cohort['recalibration_snp'])
            sub_command.add_option_long(key='tranches_file', value=file_path_dict_cohort['tranches_snp'])
            sub_command.add_option_long(key='rscript_file', value=file_path_dict_cohort['plots_snp'])

            # Run the GATK ApplyRecalibration step for SNPs.

            java_process = RunnableStep(
                name='gatk_apply_recalibration_snp',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx4G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_cohort['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='ApplyRecalibration')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            sub_command.add_option_long(key='mode', value='SNP')
            sub_command.add_option_long(key='input', value=file_path_dict_cohort['genotyped_raw_vcf'])
            sub_command.add_option_long(key='recal_file', value=file_path_dict_cohort['recalibration_snp'])
            sub_command.add_option_long(key='tranches_file', value=file_path_dict_cohort['tranches_snp'])
            sub_command.add_option_long(key='out', value=file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'])
            # The lodCutoff (VQSLOD score) filter is not applied for the moment.
            if self.truth_sensitivity_filter_level_snp:
                sub_command.add_option_long(key='ts_filter_level', value=self.truth_sensitivity_filter_level_snp)

        # Run the VQSR procedure on INDELs.

        if self.vqsr_skip_indel:
            file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf'] = \
                file_path_dict_cohort['recalibrated_snp_raw_indel_vcf']
            file_path_dict_cohort['recalibrated_snp_recalibrated_indel_idx'] = \
                file_path_dict_cohort['recalibrated_snp_raw_indel_idx']
        else:

            # Run the GATK VariantRecalibrator for INDELs.

            java_process = RunnableStep(
                name='gatk_variant_recalibrator_indel',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx8G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_cohort['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='VariantRecalibrator')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            sub_command.add_option_long(key='mode', value='INDEL')
            for resource in self.vqsr_resources_indel_dict.keys():
                resource_option = 'resource:{},known={},training={},truth={},prior={}'. \
                    format(resource,
                           self.vqsr_resources_indel_dict[resource]['known'],
                           self.vqsr_resources_indel_dict[resource]['training'],
                           self.vqsr_resources_indel_dict[resource]['truth'],
                           self.vqsr_resources_indel_dict[resource]['prior'])
                sub_command.add_option_long(
                    key=resource_option,
                    value=self.vqsr_resources_indel_dict[resource]['file_path'])
            for annotation in self.vqsr_annotations_indel_list:
                sub_command.add_option_long(key='use_annotation', value=annotation)
            sub_command.add_option_long(key='maxGaussians', value='4')  # TODO: Would be good to have this configurable.
            sub_command.add_option_long(key='input', value=file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'])
            sub_command.add_option_long(key='recal_file', value=file_path_dict_cohort['recalibration_indel'])
            sub_command.add_option_long(key='tranches_file', value=file_path_dict_cohort['tranches_indel'])
            sub_command.add_option_long(key='rscript_file', value=file_path_dict_cohort['plots_indel'])

            # Run the GATK ApplyRecalibration step for INDELs.

            java_process = RunnableStep(
                name='gatk_apply_recalibration_indel',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx4G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_cohort['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='ApplyRecalibration')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)
            sub_command.add_option_long(key='mode', value='INDEL')
            sub_command.add_option_long(key='input', value=file_path_dict_cohort['recalibrated_snp_raw_indel_vcf'])
            sub_command.add_option_long(key='recal_file', value=file_path_dict_cohort['recalibration_indel'])
            sub_command.add_option_long(key='tranches_file', value=file_path_dict_cohort['tranches_indel'])
            sub_command.add_option_long(
                key='out',
                value=file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf'])
            # The lodCutoff (VQSLOD score) filter is not applied for the moment.
            if self.truth_sensitivity_filter_level_indel:
                sub_command.add_option_long(key='ts_filter_level', value=self.truth_sensitivity_filter_level_indel)

        # In case accessory GVCF files have been used, re-create a multi-sample VCF file with just the samples
        # in this cohort.

        if len(self.accessory_cohort_gvcfs):
            java_process = RunnableStep(
                name='gatk_select_variants_cohort',
                program='java',
                sub_command=Command(command=str()))
            runnable_process_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx4G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_cohort['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='SelectVariants')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)

            sub_command.add_option_long(
                key='variant',
                value=file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf'])
            sub_command.add_option_long(key='out', value=file_path_dict_cohort['multi_sample_vcf'])
            for sample in self.samples:
                sub_command.add_option_long(key='sample_name', value=sample.name)
            sub_command.add_switch_long(key='excludeNonVariants')
        else:
            file_path_dict_cohort['multi_sample_vcf'] = \
                file_path_dict_cohort['recalibrated_snp_recalibrated_indel_vcf']
            file_path_dict_cohort['multi_sample_idx'] = \
                file_path_dict_cohort['recalibrated_snp_recalibrated_indel_idx']

        # Run the snpEff tool for functional variant annotation.

        java_process = RunnableStep(
            name='snpeff',
            program='java',
            sub_command=Command(command='eff'))
        runnable_process_cohort.add_runnable_step(runnable_step=java_process)

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

        # Run the GATK VariantAnnotator

        java_process = RunnableStep(
            name='gatk_variant_annotator',
            program='java',
            sub_command=Command(command=str()))
        runnable_process_cohort.add_runnable_step(runnable_step=java_process)

        java_process.add_switch_short(
            key='d64')
        java_process.add_option_short(
            key='jar', value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_switch_short(
            key='Xmx4G')
        java_process.add_option_pair(
            key='-Djava.io.tmpdir',
            value=file_path_dict_cohort['temporary_directory'])

        sub_command = java_process.sub_command
        sub_command.add_option_long(key='analysis_type', value='VariantAnnotator')
        sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
        for interval in self.exclude_intervals_list:
            sub_command.add_option_long(key='excludeIntervals', value=interval)
        for interval in self.include_intervals_list:
            sub_command.add_option_long(key='intervals', value=interval)
        if self.known_sites_discovery:
            sub_command.add_option_long(key='dbsnp', value=self.known_sites_discovery)

        # Add annotation resources and their corresponding expression options.
        for annotation_resource in self.annotation_resources_dict.keys():
            assert isinstance(annotation_resource, str)
            if len(self.annotation_resources_dict[annotation_resource][0]) \
                    and len(self.annotation_resources_dict[annotation_resource][1]):
                sub_command.add_option_long(
                    key=string.join(words=('resource', annotation_resource), sep=':'),
                    value=self.annotation_resources_dict[annotation_resource][0])
                for annotation in self.annotation_resources_dict[annotation_resource][1]:
                    assert isinstance(annotation, str)
                    sub_command.add_option_long(
                        key='expression',
                        value=string.join(words=(annotation_resource, annotation), sep='.'))

        sub_command.add_option_long(key='variant', value=file_path_dict_cohort['multi_sample_vcf'])
        # The AlleleBalanceBySample annotation does not seem to work in either GATK 3.1-1 or GATK 3.2-0.
        # sub_command.add_option_long(key='annotation', value='AlleleBalanceBySample')
        sub_command.add_option_long(key='annotation', value='SnpEff')
        sub_command.add_option_long(key='snpEffFile', value=file_path_dict_cohort['snpeff_vcf'])
        sub_command.add_option_long(key='out', value=file_path_dict_cohort['annotated_vcf'])

        # Create an Executable for processing the cohort.

        executable_process_cohort = Executable.from_analysis_runnable(
            analysis=self,
            runnable_name=runnable_process_cohort.name)
        drms_process_cohort.add_executable(executable=executable_process_cohort)

        executable_process_cohort.dependencies.extend(vc_process_cohort_dependencies)

        # Re-process and split the cohort by sample.

        for sample in self.samples:

            prefix_split = string.join(words=(drms_split_cohort.name, sample.name), sep='_')

            file_path_dict_split = dict(
                temporary_directory=prefix_split + '_temporary',
                sample_vcf=prefix_split + '.vcf',
                sample_idx=prefix_split + '.vcf.idx',
                sample_tsv=prefix_split + '.tsv')

            runnable_split_cohort = Runnable(
                name=prefix_split,
                code_module='bsf.runnables.generic',
                working_directory=self.genome_directory,
                file_path_dict=file_path_dict_split,
                debug=self.debug)
            self.add_runnable(runnable=runnable_split_cohort)

            # Run the GATK SelectVariants step to split multi-sample VCF files into one per sample.

            java_process = RunnableStep(
                name='gatk_select_variants',
                program='java',
                sub_command=Command(command=str()))
            runnable_split_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx2G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_split['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='SelectVariants')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)

            sub_command.add_option_long(key='variant', value=file_path_dict_cohort['annotated_vcf'])
            sub_command.add_option_long(key='out', value=file_path_dict_split['sample_vcf'])
            sub_command.add_option_long(key='sample_name', value=sample.name)
            sub_command.add_switch_long(key='excludeNonVariants')

            # Run the GATK VariantsToTable step.

            java_process = RunnableStep(
                name='gatk_variants_to_table',
                program='java',
                sub_command=Command(command=str()))
            runnable_split_cohort.add_runnable_step(runnable_step=java_process)

            java_process.add_switch_short(
                key='d64')
            java_process.add_option_short(
                key='jar',
                value=os.path.join(self.classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_switch_short(
                key='Xmx2G')
            java_process.add_option_pair(
                key='-Djava.io.tmpdir',
                value=file_path_dict_split['temporary_directory'])

            sub_command = java_process.sub_command
            sub_command.add_option_long(key='analysis_type', value='VariantsToTable')
            sub_command.add_option_long(key='reference_sequence', value=self.bwa_genome_db)
            for interval in self.exclude_intervals_list:
                sub_command.add_option_long(key='excludeIntervals', value=interval)
            for interval in self.include_intervals_list:
                sub_command.add_option_long(key='intervals', value=interval)

            sub_command.add_option_long(key='variant', value=file_path_dict_split['sample_vcf'])
            sub_command.add_option_long(key='out', value=file_path_dict_split['sample_tsv'])
            sub_command.add_switch_long(key='allowMissingData')
            sub_command.add_switch_long(key='showFiltered')
            # Set of standard VCF fields.
            sub_command.add_option_long(key='fields', value='CHROM')
            sub_command.add_option_long(key='fields', value='POS')
            sub_command.add_option_long(key='fields', value='ID')
            sub_command.add_option_long(key='fields', value='REF')
            sub_command.add_option_long(key='fields', value='ALT')
            sub_command.add_option_long(key='fields', value='QUAL')
            sub_command.add_option_long(key='fields', value='FILTER')
            #
            sub_command.add_option_long(key='fields', value='VQSLOD')
            sub_command.add_option_long(key='fields', value='AF')
            # GATK Haplotype Caller genotype fields: GT:AD:DP:GQ:PL
            sub_command.add_option_long(key='genotypeFields', value='GT')
            sub_command.add_option_long(key='genotypeFields', value='AD')
            sub_command.add_option_long(key='genotypeFields', value='DP')
            sub_command.add_option_long(key='genotypeFields', value='GQ')
            sub_command.add_option_long(key='genotypeFields', value='PL')
            # Set of snpEff fields.
            sub_command.add_option_long(key='fields', value='SNPEFF_EFFECT')
            sub_command.add_option_long(key='fields', value='SNPEFF_IMPACT')
            sub_command.add_option_long(key='fields', value='SNPEFF_FUNCTIONAL_CLASS')
            sub_command.add_option_long(key='fields', value='SNPEFF_CODON_CHANGE')
            sub_command.add_option_long(key='fields', value='SNPEFF_AMINO_ACID_CHANGE')
            sub_command.add_option_long(key='fields', value='SNPEFF_GENE_NAME')
            sub_command.add_option_long(key='fields', value='SNPEFF_GENE_BIOTYPE')
            sub_command.add_option_long(key='fields', value='SNPEFF_TRANSCRIPT_ID')
            sub_command.add_option_long(key='fields', value='SNPEFF_EXON_ID')

            # Automatically add all fields defined for the Variant Annotator resources, above.
            for annotation_resource in self.annotation_resources_dict.keys():
                assert isinstance(annotation_resource, str)
                if len(self.annotation_resources_dict[annotation_resource][0]) \
                        and len(self.annotation_resources_dict[annotation_resource][1]):
                    for annotation in self.annotation_resources_dict[annotation_resource][1]:
                        assert isinstance(annotation, str)
                        sub_command.add_option_long(
                            key='fields',
                            value=string.join(words=(annotation_resource, annotation), sep='.'))

            # Create an Executable for splitting the cohort.

            executable_split_cohort = Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_split_cohort.name)
            drms_split_cohort.add_executable(executable=executable_split_cohort)

            executable_split_cohort.dependencies.append(executable_process_cohort.name)

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

        output += '<h2>Sample and Replicate Level</h2>\n'

        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr>\n'
        output += '<th>Sample</th>\n'
        output += '<th>Variants VCF file</th>\n'
        output += '<th>Variants TSV file</th>\n'
        output += '<th>Aligned BAM file</th>\n'
        output += '<th>Aligned BAI file</th>\n'
        output += '<th>Replicate</th>\n'
        output += '<th>Duplicate Metrics</th>\n'
        output += '<th>Alignment Summary Metrics</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        # Group via UCSC super tracks.

        track_output += 'track Alignments\n'
        track_output += 'shortLabel Alignments\n'
        track_output += 'longLabel BWA NGS read alignments\n'
        track_output += 'visibility show\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            runnable_process_sample = self.runnable_dict[
                string.join(words=(self.drms_name_process_sample, sample.name), sep='_')]
            assert isinstance(runnable_process_sample, Runnable)
            runnable_split_cohort = self.runnable_dict[
                string.join(words=(self.drms_name_split_cohort, sample.name), sep='_')]
            assert isinstance(runnable_split_cohort, Runnable)

            track_output += 'track {}_alignments\n'. \
                format(sample.name)
            track_output += 'type bam\n'
            track_output += 'shortLabel {}_alignments\n'. \
                format(sample.name)
            track_output += 'longLabel {} BWA NGS read alignments\n'. \
                format(sample.name)
            track_output += 'bigDataUrl ./{}\n'. \
                format(runnable_process_sample.file_path_dict['realigned_bam'])
            track_output += 'visibility squish\n'
            # track_output += 'html {}\n'.format()

            # Common optional settings.

            track_output += 'color {}\n'. \
                format('0,0,0')

            # Compressed Sequence Alignment track settings.

            # None so far.

            # Composite track settings.

            track_output += 'parent Alignments\n'
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
            output += '<td></td>\n'  # Replicate
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
                    string.join(words=(self.drms_name_process_lane, replicate_key), sep='_')]
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
            string.join(words=(self.drms_name_process_cohort, self.cohort_name), sep='_')]
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
