#! /usr/bin/env python
#
# BSF Python script to run MuTect after the GATK pipeline.
#
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


from argparse import ArgumentParser
import os
from pickle import Pickler, HIGHEST_PROTOCOL
import string

from Bio.BSF import Analysis, Command, Default, DRMS, Executable
from Bio.BSF.annotation import SampleAnnotationSheet


def _read_comparisons(analysis, cmp_file):
    sas = SampleAnnotationSheet.read_from_file(file_path=cmp_file, name='MuTect Comparisons')

    for row_dict in sas.row_dicts:

        key = str()
        comparison_groups = list()

        for prefix in 'Normal', 'Tumor':
            group_name, group_samples = analysis.collection.get_Samples_from_row_dict(
                row_dict=row_dict, prefix=prefix)
            if group_name and len(group_samples):
                key += group_name
                key += '__'
                # key = string.join(words=(key, group_name), sep='__')
                comparison_groups.append((group_name, group_samples))
                # Also expand each Python list of Sample objects to get all those Sample objects
                # that this Analysis needs considering.
                for sample in group_samples:
                    if analysis.debug > 1:
                        print '  {} Sample name: {!r} file_path:{!r}'.format(prefix, sample.name, sample.file_path)
                        # print sample.trace(1)
                    analysis.add_Sample(sample=sample)

        analysis.comparisons[key[:-2]] = comparison_groups


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = ArgumentParser(
    description='BSF script to run MuTect.')

argument_parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--stage',
    help='limit job submission to a particular Analysis stage',
    dest='stage',
    required=False)

argument_parser.add_argument(
    'configuration',
    help='configuration (*.ini) file path')

name_space = argument_parser.parse_args()

# TODO: Load the Sample Annotation sheet and comparison sheets.

analysis = Analysis.from_config_file(config_file=name_space.configuration)

if name_space.debug:
    analysis.debug = name_space.debug

# Call the Analysis.run method here to properly initialise the Analysis object.
analysis.run()

# Get global defaults.

default = Default.get_global_default()

config_parser = analysis.configuration.config_parser
config_section = analysis.configuration.section_from_instance(analysis)

if config_parser.has_option(section=config_section, option='gatk_bundle_version'):
    gatk_bundle_version = config_parser.get(section=config_section, option='gatk_bundle_version')
else:
    raise Exception("A 'VariantCallingGATK' analysis requires a gatk_bundle_version configuration option.")

if config_parser.has_option(section=config_section, option='classpath_gatk'):
    classpath_gatk = config_parser.get(section=config_section, option='classpath_gatk')
else:
    classpath_gatk = default.classpath_gatk

# TODO: Should be configurable via Default.
classpath_mutect = '/cm/shared/apps/muTect/1.1.5'
classpath_indel_genotyper = '/home/mschuster/src'

if config_parser.has_option(section=config_section, option='classpath_snpeff'):
    classpath_snpeff = config_parser.get(section=config_section, option='classpath_snpeff')
else:
    classpath_snpeff = default.classpath_snpeff

bwa_genome_db = str(config_parser.get(section=config_section, option='bwa_genome_db'))
if not os.path.isabs(bwa_genome_db):
    bwa_genome_db = os.path.join(
        Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                     genome_version=analysis.genome_version),
        bwa_genome_db)

snpeff_genome_version = config_parser.get(section=config_section, option='snpeff_genome_version')

include_intervals_list = list()
if config_parser.has_option(section=config_section, option='include_intervals'):
    annotation_resources = config_parser.get(section=config_section, option='include_intervals')
    include_intervals_list.extend(annotation_resources.split(','))

exclude_intervals_list = list()
if config_parser.has_option(section=config_section, option='exclude_intervals'):
    annotation_resources = config_parser.get(section=config_section, option='exclude_intervals')
    exclude_intervals_list.extend(annotation_resources.split(','))

# Single VCF file of known sites for the MuTect step.
known_sites_discovery = str()
if config_parser.has_option(section=config_section, option='known_sites_discovery'):
    known_sites_discovery = str(config_parser.get(section=config_section, option='known_sites_discovery'))
    if not os.path.isabs(known_sites_discovery):
        known_sites_discovery = os.path.join(
            Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                         genome_version=analysis.genome_version),
            known_sites_discovery)

known_sites_somatic = str()
if config_parser.has_option(section=config_section, option='known_sites_somatic'):
    known_sites_somatic = str(config_parser.get(section=config_section, option='known_sites_somatic'))
    if not os.path.isabs(known_sites_somatic):
        raise Exception("Need absolute path to 'known_sites_somatic' not {!r}.".format(known_sites_somatic))

# Read additionally requested annotation resources for the GATK AnnotateVariants step.

# Python dict of Python str (annotation resource name) key and
# Python tuple of
# Python str (file path) and Python list of Python str (annotation) value data.
annotation_resources_dict = dict()

if config_parser.has_option(section=config_section, option='annotation_resources'):
    for annotation_resource in \
            config_parser.get(section=config_section, option='annotation_resources').split(','):
        resource_section = string.join(words=(annotation_resource, 'resource'), sep='_')
        if config_parser.has_section(section=resource_section):
            annotation_list = list()
            if config_parser.has_option(section=resource_section, option='file_path'):
                file_path = str(config_parser.get(section=resource_section, option='file_path'))
                if not os.path.isabs(file_path):
                    file_path = os.path.join(
                        Default.absolute_gatk_bundle(
                            gatk_bundle_version=gatk_bundle_version,
                            genome_version=analysis.genome_version),
                        file_path)
            else:
                raise Exception(
                    "Missing configuration option 'file_path' in configuration section {!r}.".
                    format(resource_section))
            if config_parser.has_option(section=resource_section, option='annotations'):
                annotation_list.extend(
                    config_parser.get(section=resource_section, option='annotations').split(','))
            else:
                raise Exception(
                    "Missing configuration option 'annotations' in configuration section {!r}.".
                    format(resource_section))
            # Create a dict key and a tuple of a Python str and Python list.
            if annotation_resource not in annotation_resources_dict:
                annotation_resources_dict[annotation_resource] = file_path, annotation_list
        else:
            raise Exception(
                'Missing configuration section {!r} declared in option annotation_resources {!r}.'.
                format(resource_section,
                       config_parser.get(section=config_section, option='annotation_resources')))

# Initialise a Distributed Resource Management System (DRMS) object for the
# bsf_run_variant_calling_somatic.py script.

vc_run_somatic_drms = DRMS.from_Analysis(
    name='variant_calling_somatic',
    work_directory=analysis.genome_directory,
    analysis=analysis)
analysis.drms_list.append(vc_run_somatic_drms)

_read_comparisons(analysis=analysis, cmp_file=config_parser.get(section=config_section, option='cmp_file'))

key_list = analysis.comparisons.keys()
key_list.sort(cmp=lambda x, y: cmp(x, y))

for key in key_list:

    # The list of samples must contain exactly one normal and one tumor samples.
    if len(analysis.comparisons[key]) != 2:
        continue

    prefix_somatic = string.join((vc_run_somatic_drms.name, key), sep='_')

    # Somatic variant calling-specific file paths

    file_path_somatic = dict(
        # TODO: the prefix_somatic is everything that is needed here.
        temporary_directory=prefix_somatic + '_temporary',
        mutect_vcf=prefix_somatic + '_mutations.vcf',
        mutect_idx=prefix_somatic + '_mutations.vcf.idx',
        mutect_out=prefix_somatic + '_mutations.txt',
        mutect_wig=prefix_somatic + '_mutations.wig',
        indel_vcf=prefix_somatic + '_indels.vcf',
        indel_idx=prefix_somatic + '_indels.vcf.idx',
        indel_bed=prefix_somatic + '_indels.bed',
        indel_vrb=prefix_somatic + '_indels.txt',
        combined_vcf=prefix_somatic + '_combined.vcf',
        combined_idx=prefix_somatic + '_combined.vcf.idx',
        snpeff_vcf=prefix_somatic + '_snpeff.vcf',
        snpeff_idx=prefix_somatic + '_snpeff.vcf.idx',
        snpeff_stats=prefix_somatic + '_snpeff_summary.html',
        annotated_vcf=prefix_somatic + '_annotated.vcf',
        annotated_idx=prefix_somatic + '_annotated.vcf.idx',
        annotated_csv=prefix_somatic + '_annotated.csv')

    # Somatic variant calling-specific pickler_dict

    pickler_dict_somatic = dict(
        file_path_dict=file_path_somatic,
        prefix=vc_run_somatic_drms.name,
        comparison_key=key)

    # Run the MuTect analysis

    java_process = Executable(name='mutect',
                              program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_mutect, 'muTect.jar'))
    java_process.add_SwitchShort(key='Xmx4G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_somatic['temporary_directory'])

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='MuTect')
    sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
    for interval in exclude_intervals_list:
        sub_command.add_OptionLong(key='excludeIntervals', value=interval)
    for interval in include_intervals_list:
        sub_command.add_OptionLong(key='intervals', value=interval)
    if known_sites_discovery:
        sub_command.add_OptionLong(key='dbsnp', value=known_sites_discovery)
    if known_sites_somatic:
        sub_command.add_OptionLong(key='cosmic', value=known_sites_somatic)

    sub_command.add_OptionLong(
        key='input_file:normal',
        value='variant_calling_process_sample_{}_realigned.bam'.format(analysis.comparisons[key][0][1][0].name))
    sub_command.add_OptionLong(
        key='input_file:tumor',
        value='variant_calling_process_sample_{}_realigned.bam'.format(analysis.comparisons[key][-1][1][0].name))

    sub_command.add_OptionLong(key='out', value=file_path_somatic['mutect_out'])
    sub_command.add_OptionLong(key='vcf', value=file_path_somatic['mutect_vcf'])
    sub_command.add_OptionLong(key='coverage_file', value=file_path_somatic['mutect_wig'])

    pickler_dict_somatic[java_process.name] = java_process

    # Run the Indel Genotyper analysis
    # Note that the Indel Genotyper is based on a much older GATK version.

    java_process = Executable(name='indel_genotyper',
                              program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_indel_genotyper,
                                                               'IndelGenotyper.36.3336-GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx4G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_somatic['temporary_directory'])

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='IndelGenotyperV2')
    sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
    for interval in exclude_intervals_list:
        sub_command.add_OptionLong(key='excludeIntervals', value=interval)
    for interval in include_intervals_list:
        sub_command.add_OptionLong(key='intervals', value=interval)
    # Not supported by the old GATK version behind the Somatic Indel Genotyper
    # MESSAGE: --DBSNP (-D) argument currently does not support VCF.
    # To use dbSNP in VCF format, please use -B:dbsnp,vcf <filename>.
    # if known_sites_discovery:
    # sub_command.add_OptionLong(key='DBSNP', value=known_sites_discovery)

    sub_command.add_SwitchLong(key='somatic')
    sub_command.add_OptionLong(
        key='input_file:normal',
        value='variant_calling_process_sample_{}_realigned.bam'.format(analysis.comparisons[key][0][1][0].name))
    sub_command.add_OptionLong(
        key='input_file:tumor',
        value='variant_calling_process_sample_{}_realigned.bam'.format(analysis.comparisons[key][-1][1][0].name))

    sub_command.add_OptionLong(key='out', value=file_path_somatic['indel_vcf'])
    sub_command.add_OptionLong(key='bedOutput', value=file_path_somatic['indel_bed'])
    sub_command.add_OptionLong(key='verboseOutput', value=file_path_somatic['indel_vrb'])
    # Extend the window size to get around a bug in the IndelGenotyperV2? Sigh.
    # ##### ERROR MESSAGE: Invalid command line: Argument window_size has a bad value:
    #  Read HWI-ST181_0391:7:2214:8748:86539#737C: out of coverage window bounds.
    #  Probably window is too small, so increase the value of the window_size argument.
    # ##### ERROR Read length=100; cigar=100M; start=3833506; end=3833605;
    #  window start (after trying to accomodate the read)=3833405; window end=3833604

    sub_command.add_OptionLong(key='window_size', value='1000')

    pickler_dict_somatic[java_process.name] = java_process

    # Run the GATK Combine Variants analysis

    java_process = Executable(name='gatk_combine_variants',
                              program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_somatic['temporary_directory'])

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='CombineVariants')
    sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
    for interval in exclude_intervals_list:
        sub_command.add_OptionLong(key='excludeIntervals', value=interval)
    for interval in include_intervals_list:
        sub_command.add_OptionLong(key='intervals', value=interval)

    # TODO: Should this use the option --assumeIdenticalSamples to just concatenate the VCFs?
    sub_command.add_OptionLong(key='variant', value=file_path_somatic['mutect_vcf'])
    sub_command.add_OptionLong(key='variant', value=file_path_somatic['indel_vcf'])
    sub_command.add_OptionLong(key='out', value=file_path_somatic['combined_vcf'])

    pickler_dict_somatic[java_process.name] = java_process

    # Run the snpEff tool for functional variant annotation.

    java_process = Executable(name='snpeff',
                              program='java',
                              sub_command=Command(command='eff'))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_snpeff, 'snpEff.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_somatic['temporary_directory'])
    java_process.stdout_path = file_path_somatic['snpeff_vcf']

    sub_command = java_process.sub_command
    sub_command.add_SwitchShort(key='download')
    sub_command.add_OptionShort(key='o', value='gatk')
    sub_command.add_OptionShort(key='stats', value=file_path_somatic['snpeff_stats'])
    sub_command.add_OptionShort(key='config', value=os.path.join(classpath_snpeff, 'snpEff.config'))

    sub_command.arguments.append(snpeff_genome_version)
    sub_command.arguments.append(file_path_somatic['combined_vcf'])

    pickler_dict_somatic[java_process.name] = java_process

    # Run the GATK Variant Annotator

    java_process = Executable(name='gatk_variant_annotator',
                              program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_somatic['temporary_directory'])

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='VariantAnnotator')
    sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
    for interval in exclude_intervals_list:
        sub_command.add_OptionLong(key='excludeIntervals', value=interval)
    for interval in include_intervals_list:
        sub_command.add_OptionLong(key='intervals', value=interval)
    if known_sites_discovery:
        sub_command.add_OptionLong(key='dbsnp', value=known_sites_discovery)

    # Add annotation resources and their corresponding expression options.
    for annotation_resource in annotation_resources_dict.keys():
        if len(annotation_resources_dict[annotation_resource][0]) \
                and len(annotation_resources_dict[annotation_resource][1]):
            sub_command.add_OptionLong(
                key=string.join(words=('resource', annotation_resource), sep=':'),
                value=annotation_resources_dict[annotation_resource][0])
            for annotation in annotation_resources_dict[annotation_resource][1]:
                sub_command.add_OptionLong(
                    key='expression',
                    value=string.join(words=(annotation_resource, annotation), sep='.'))

    sub_command.add_OptionLong(key='variant', value=file_path_somatic['combined_vcf'])
    # TODO: Test whether annotation AlleleBalanceBySample works in GATK 3.2.
    sub_command.add_OptionLong(key='annotation', value='AlleleBalanceBySample')
    sub_command.add_OptionLong(key='annotation', value='SnpEff')
    sub_command.add_OptionLong(key='snpEffFile', value=file_path_somatic['snpeff_vcf'])
    sub_command.add_OptionLong(key='out', value=file_path_somatic['annotated_vcf'])

    pickler_dict_somatic[java_process.name] = java_process

    # Run the GATK VariantsToTable step.

    java_process = Executable(name='gatk_variants_to_table',
                              program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx4G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_somatic['temporary_directory'])

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='VariantsToTable')
    sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
    for interval in exclude_intervals_list:
        sub_command.add_OptionLong(key='excludeIntervals', value=interval)
    for interval in include_intervals_list:
        sub_command.add_OptionLong(key='intervals', value=interval)

    sub_command.add_OptionLong(key='variant', value=file_path_somatic['annotated_vcf'])
    sub_command.add_OptionLong(key='out', value=file_path_somatic['annotated_csv'])
    sub_command.add_SwitchLong(key='allowMissingData')
    sub_command.add_SwitchLong(key='showFiltered')
    # Set of standard VCF fields.
    sub_command.add_OptionLong(key='fields', value='CHROM')
    sub_command.add_OptionLong(key='fields', value='POS')
    sub_command.add_OptionLong(key='fields', value='ID')
    sub_command.add_OptionLong(key='fields', value='REF')
    sub_command.add_OptionLong(key='fields', value='ALT')
    sub_command.add_OptionLong(key='fields', value='QUAL')
    sub_command.add_OptionLong(key='fields', value='FILTER')
    #
    sub_command.add_OptionLong(key='fields', value='AF')
    # MuTect genotype fields  GT:AD:BQ:DP:FA
    sub_command.add_OptionLong(key='genotypeFields', value='GT')
    sub_command.add_OptionLong(key='genotypeFields', value='AD')
    sub_command.add_OptionLong(key='genotypeFields', value='BQ')
    sub_command.add_OptionLong(key='genotypeFields', value='DP')
    sub_command.add_OptionLong(key='genotypeFields', value='FA')
    # Indel Genotyper genotype fields GT:GQ
    sub_command.add_OptionLong(key='genotypeFields', value='GQ')
    # Set of snpEff fields.
    sub_command.add_OptionLong(key='fields', value='SNPEFF_EFFECT')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_IMPACT')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_FUNCTIONAL_CLASS')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_CODON_CHANGE')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_AMINO_ACID_CHANGE')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_GENE_NAME')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_GENE_BIOTYPE')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_TRANSCRIPT_ID')
    sub_command.add_OptionLong(key='fields', value='SNPEFF_EXON_ID')

    # Automatically add all fields defined for the Variant Annotator resources, above.
    for annotation_resource in annotation_resources_dict.keys():
        if len(annotation_resources_dict[annotation_resource][0]) \
                and len(annotation_resources_dict[annotation_resource][1]):
            for annotation in annotation_resources_dict[annotation_resource][1]:
                sub_command.add_OptionLong(
                    key='fields',
                    value=string.join(words=(annotation_resource, annotation), sep='.'))

    pickler_dict_somatic[java_process.name] = java_process

    # Write the Pickler dict file for somatic variant calling.

    pickler_path = os.path.join(analysis.genome_directory, prefix_somatic + '.pkl')
    pickler_file = open(pickler_path, 'wb')
    pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
    pickler.dump(obj=pickler_dict_somatic)
    pickler_file.close()

    # Create a BSF Executable for somatic variant calling.

    vc_run_somatic = Executable.from_Analysis(
        name=prefix_somatic,
        program='bsf_run_variant_calling_somatic.py',
        analysis=analysis)
    vc_run_somatic_drms.add_Executable(vc_run_somatic)

    # Only submit this Executable if the final result file does not exist.
    if (os.path.exists(
            os.path.join(analysis.genome_directory, file_path_somatic['annotated_csv']))
        and os.path.getsize(
                os.path.join(analysis.genome_directory, file_path_somatic['annotated_csv']))):
        vc_run_somatic.submit = False

    # TODO: Set dependencies on sample-level processing in case this gets moved into the VariantCalling module.

    # Set variant_calling_run_process_lane options.

    vc_run_somatic.add_OptionLong(key='pickler_path', value=pickler_path)
    vc_run_somatic.add_OptionLong(key='debug', value=str(analysis.debug))

# Submit all Executable objects of all Distributed Resource Management System objects.

submit = 0

for drms in analysis.drms_list:

    if name_space.stage:
        if name_space.stage == drms.name:
            submit += 1
        else:
            continue

    drms.submit(debug=analysis.debug)

    if analysis.debug:
        print repr(drms)
        print drms.trace(1)

if name_space.stage:
    if name_space.stage == 'report':
        analysis.report()
        pass
    elif not submit:
        name_list = [drms.name for drms in analysis.drms_list]
        name_list.append('report')
        print 'Valid Analysis stages are: {!r}'.format(name_list)

print 'Somatic Variant Calling Analysis'
print 'Project name:      ', analysis.project_name
print 'Genome version:    ', analysis.genome_version
print 'Input directory:   ', analysis.input_directory
print 'Output directory:  ', analysis.output_directory
print 'Project directory: ', analysis.project_directory
print 'Genome directory:  ', analysis.genome_directory
