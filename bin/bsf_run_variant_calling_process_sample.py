#! /usr/bin/env python
#
# BSF Python script to process samples for the GATK pipeline.
#
#   Picard MergeSamFiles
#   Picard MarkDuplicates
#   GATK RealignerTargetCreator
#   GATK IndelRealigner
#   Picard CollectAlignmentSummaryMetrics
#   GATK HaplotypeCaller
#
#
# Copyright 2013 Michael K. Schuster
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


import argparse
import errno
import os.path
from pickle import Unpickler
import shutil

from Bio.BSF import Command, Default, Executable, Runnable

# Set the environment consistently.

os.environ['LANG'] = 'C'

# Get global defaults.

default = Default.get_global_default()

# Parse the arguments.

parser = argparse.ArgumentParser(description='BSF Runner for the Burrows Wheeler Aligner (BWA).')

parser.add_argument('--debug', required=False, type=int,
                    help='debug level')

parser.add_argument('--pickler_path', required=True,
                    help='file path to a Python Pickler file.')

args = parser.parse_args()

# Unpickle the file into a Python dict object.

pickler_file = open(args.pickler_path, 'rb')

unpickler = Unpickler(file=pickler_file)

pickler_dict = unpickler.load()

pickler_file.close()

# Picker dictionary entries that must exist.

key = 'prefix'
if key in pickler_dict and pickler_dict[key]:
    prefix = pickler_dict[key]
    if prefix[-1:] != '_':
        prefix += '_'
else:
    prefix = str()

key = 'sample_key'
if key in pickler_dict and pickler_dict[key]:
    sample_key = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

key = 'classpath_gatk'
if key in pickler_dict and pickler_dict[key]:
    classpath_gatk = pickler_dict[key]
else:
    classpath_gatk = default.classpath_gatk

key = 'classpath_picard'
if key in pickler_dict and pickler_dict[key]:
    classpath_picard = pickler_dict[key]
else:
    classpath_picard = default.classpath_picard

key = 'path_reference_sequence'
if key in pickler_dict and pickler_dict[key]:
    path_reference_sequence = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

key = 'replicate_file_list'
if key in pickler_dict and pickler_dict[key]:
    replicate_file_list = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

key = 'known_sites_realignment'
if key in pickler_dict and pickler_dict[key]:
    known_sites_realignment = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

# The known_sites_recalibration parameter is not required here.

key = 'known_sites_discovery'
if key in pickler_dict and pickler_dict[key]:
    known_sites_discovery = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, sample_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

path_merged_bam = '{}{}_merged.bam'.format(prefix, sample_key)
path_merged_bai = '{}{}_merged.bai'.format(prefix, sample_key)
path_merged_md5 = '{}{}_merged.bam.md5'.format(prefix, sample_key)
path_duplicates_marked_bam = '{}{}_duplicates_marked.bam'.format(prefix, sample_key)
path_duplicates_marked_bai = '{}{}_duplicates_marked.bai'.format(prefix, sample_key)
path_duplicates_marked_md5 = '{}{}_duplicates_marked.bam.md5'.format(prefix, sample_key)
path_realigner_targets = '{}{}_realigner.intervals'.format(prefix, sample_key)
path_realigned_bam = '{}{}_realigned.bam'.format(prefix, sample_key)
path_realigned_bai = '{}{}_realigned.bai'.format(prefix, sample_key)
path_alignment_summary_metrics = '{}{}_alignment_summary_metrics.csv'.format(prefix, sample_key)
path_raw_variants = '{}{}_raw_variants.vcf'.format(prefix, sample_key)

# Run the Picard MergeSamFiles step.

if not os.path.exists(path_duplicates_marked_bam):
    java_process = Executable(name='picard_merge_sam_files', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'MergeSamFiles.jar'))
    java_process.add_SwitchShort(key='Xmx6G')

    sub_command = java_process.sub_command
    for replicate_file in replicate_file_list:
        sub_command.add_OptionPair(key='INPUT', value=replicate_file)
    sub_command.add_OptionPair(key='OUTPUT', value=path_merged_bam)
    sub_command.add_OptionPair(key='COMMENT', value='Merged from the following files:')
    for replicate_file in replicate_file_list:
        sub_command.add_OptionPair(key='COMMENT', value=replicate_file[:-4])
    sub_command.add_OptionPair(key='TMP_DIR', value=path_temporary)
    sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
    sub_command.add_OptionPair(key='QUIET', value='false')
    sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
    sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
    sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
    sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
    sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the Picard MergeSamFiles step.'
        raise Exception(message)

# Run the Picard MarkDuplicates step.
# Optical duplicates should already have been flagged in the lane-specific processing step.

if not os.path.exists(path_duplicates_marked_bam):
    java_process = Executable(name='picard_mark_duplicates', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'MarkDuplicates.jar'))
    java_process.add_SwitchShort(key='Xmx6G')

    sub_command = java_process.sub_command
    sub_command.add_OptionPair(key='INPUT', value=path_merged_bam)
    sub_command.add_OptionPair(key='OUTPUT', value=path_duplicates_marked_bam)
    sub_command.add_OptionPair(key='METRICS_FILE', value='{}{}_duplicate_metrics.csv'.format(prefix, sample_key))
    # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
    # as otherwise optical duplicates could not be detected.
    # See BioStar post: http://www.biostars.org/p/12538/
    # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.
    # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.
    sub_command.add_OptionPair(key='READ_NAME_REGEX', value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.')
    sub_command.add_OptionPair(key='TMP_DIR', value=path_temporary)
    sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
    sub_command.add_OptionPair(key='QUIET', value='false')
    sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
    sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
    sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
    sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
    sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the Picard MarkDuplicates step.'
        raise Exception(message)

# if args.debug < 1:
#     os.remove(path_merged_bam)
#     os.remove(path_merged_bai)
#     os.remove(path_merged_md5)

# Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

if not os.path.exists(path_realigner_targets):
    java_process = Executable(name='gatk_realigner_target_creator', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='RealignerTargetCreator')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    for file_path in known_sites_realignment:
        sub_command.add_OptionLong(key='known', value=file_path)
    sub_command.add_OptionLong(key='input_file', value=path_duplicates_marked_bam)
    sub_command.add_OptionLong(key='out', value=path_realigner_targets)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK RealignerTargetCreator step.'

# Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

if not os.path.exists(path_realigned_bam):
    java_process = Executable(name='gatk_indel_realigner', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='IndelRealigner')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    for file_path in known_sites_realignment:
        sub_command.add_OptionLong(key='knownAlleles', value=file_path)
    sub_command.add_OptionLong(key='input_file', value=path_duplicates_marked_bam)
    sub_command.add_OptionLong(key='targetIntervals', value=path_realigner_targets)
    sub_command.add_OptionLong(key='out', value=path_realigned_bam)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK IndelRealigner step.'
        raise Exception(message)

# if args.debug < 1:
#     os.remove(path_duplicates_marked_bam)
#     os.remove(path_duplicates_marked_bai)
#     os.remove(path_duplicates_marked_md5)

# The GATK BaseRecalibrator step has already been run on a per lane basis.

# Run the Picard CollectAlignmentSummaryMetrics step.

if not os.path.exists(path_alignment_summary_metrics):
    java_process = Executable(name='picard_collect_alignment_summary_metrics', program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'CollectAlignmentSummaryMetrics.jar'))
    java_process.add_SwitchShort(key='Xmx6G')

    sub_command = java_process.sub_command
    sub_command.add_OptionPair(key='INPUT', value=path_realigned_bam)
    sub_command.add_OptionPair(key='OUTPUT', value=path_alignment_summary_metrics)
    sub_command.add_OptionPair(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS')
    sub_command.add_OptionPair(key='REFERENCE_SEQUENCE', value=path_reference_sequence)
    sub_command.add_OptionPair(key='TMP_DIR', value=path_temporary)
    sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
    sub_command.add_OptionPair(key='QUIET', value='false')
    sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
    sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
    sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
    sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
    sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the Picard CollectAlignmentSummaryMetrics step.'
        raise Exception(message)

# Run the GATK HaplotypeCaller per sample.

if not os.path.exists(path_raw_variants):
    java_process = Executable(name='gatk_haplotype_caller', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='HaplotypeCaller')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    sub_command.add_OptionLong(key='genotyping_mode', value='DISCOVERY')
    sub_command.add_OptionLong(key='standard_min_confidence_threshold_for_emitting', value='10')
    sub_command.add_OptionLong(key='standard_min_confidence_threshold_for_calling', value='30')
    sub_command.add_OptionLong(key='emitRefConfidence', value='GVCF')
    for file_path in known_sites_discovery:
        sub_command.add_OptionLong(key='dbsnp', value=file_path)
    sub_command.add_OptionLong(key='input_file', value=path_realigned_bam)
    sub_command.add_OptionLong(key='out', value=path_raw_variants)
    # Parameter to pass to the VCF/BCF IndexCreator
    sub_command.add_OptionLong(key='variant_index_type', value='LINEAR')
    sub_command.add_OptionLong(key='variant_index_parameter', value='128000')

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK HaplotypeCaller step.'
        raise Exception(message)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
