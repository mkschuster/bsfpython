#! /usr/bin/env python
#
# BSF Python script to process lanes for the GATK pipeline.
#
#   Picard MarkDuplicates
#   GATK RealignerTargetCreator
#   GATK IndelRealigner
#   GATK BaseRecalibrator first-pass
#   GATK BaseRecalibrator second-pass
#   GATK AnalyzeCovariates
#   GATK PrintReads
#   Picard CollectAlignmentSummaryMetrics
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

parser = argparse.ArgumentParser(
    description='BSF Runner for post-processing alignments prior to GATK variant calling.')

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

key = 'replicate_key'
if key in pickler_dict and pickler_dict[key]:
    replicate_key = pickler_dict[key]
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

key = 'path_replicate'
if key in pickler_dict and pickler_dict[key]:
    path_replicate = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

key = 'path_reference_sequence'
if key in pickler_dict and pickler_dict[key]:
    path_reference_sequence = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

key = 'known_sites_realignment'
# Mandatory key with optional value; an empty Python list is fine.
# if key in pickler_dict and pickler_dict[key]:
if key in pickler_dict:
    known_sites_realignment = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

key = 'known_sites_recalibration'
# Mandatory key with optional value; an empty Python list is fine.
# if key in pickler_dict and pickler_dict[key]:
if key in pickler_dict:
    known_sites_recalibration = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, replicate_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

path_aligned_bam = str(path_replicate)
# path_aligned_bam = '{}{}.bam'.format(prefix, replicate_key)
# path_aligned_bai = '{}{}.bai'.format(prefix, replicate_key)
# path_aligned_md5 = '{}{}.bam.md5'.format(prefix, replicate_key)
path_duplicates_marked_bam = '{}{}_duplicates_marked.bam'.format(prefix, replicate_key)
path_duplicates_marked_bai = '{}{}_duplicates_marked.bai'.format(prefix, replicate_key)
path_duplicates_marked_md5 = '{}{}_duplicates_marked.bam.md5'.format(prefix, replicate_key)
path_realigner_targets = '{}{}_realigner.intervals'.format(prefix, replicate_key)
path_realigned_bam = '{}{}_realigned.bam'.format(prefix, replicate_key)
path_realigned_bai = '{}{}_realigned.bai'.format(prefix, replicate_key)
path_recalibration_table_pre = '{}{}_recalibration_pre.table'.format(prefix, replicate_key)
path_recalibration_table_post = '{}{}_recalibration_post.table'.format(prefix, replicate_key)
path_recalibration_plot = '{}{}_recalibration_report.pdf'.format(prefix, replicate_key)
path_recalibrated_bam = '{}{}_recalibrated.bam'.format(prefix, replicate_key)
path_recalibrated_bai = '{}{}_recalibrated.bai'.format(prefix, replicate_key)
path_alignment_summary_metrics = '{}{}_alignment_summary_metrics.csv'.format(prefix, replicate_key)

# Run the Picard MarkDuplicates step.

if not os.path.exists(path_duplicates_marked_bam):
    java_process = Executable(name='picard_mark_duplicates', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'MarkDuplicates.jar'))
    java_process.add_SwitchShort(key='Xmx6G')

    sub_command = java_process.sub_command
    sub_command.add_OptionPair(key='INPUT', value=path_aligned_bam)
    sub_command.add_OptionPair(key='OUTPUT', value=path_duplicates_marked_bam)
    sub_command.add_OptionPair(key='METRICS_FILE', value='{}{}_duplicate_metrics.csv'.format(prefix, replicate_key))
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

# TODO: Not sure the aligned BAM file from the previous step should be deleted here.
# if args.debug < 1:
#     os.remove(path_aligned_bam)
#     os.remove(path_aligned_bai)
#     os.remove(path_aligned_md5)

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
    java_process.add_SwitchShort(key='Xmx8G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='IndelRealigner')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    for file_path in known_sites_realignment:
        sub_command.add_OptionLong(key='knownAlleles', value=file_path)
    sub_command.add_OptionLong(key='input_file', value=path_duplicates_marked_bam)
    sub_command.add_OptionLong(key='targetIntervals', value=path_realigner_targets)
    sub_command.add_OptionLong(key='out', value=path_realigned_bam)
    # TODO: For debugging only.
    sub_command.add_OptionLong(key='logging_level', value='DEBUG')

    child_return_code = Runnable.run(executable=java_process, max_loop_counter=1)

    if child_return_code:
        message = 'Could not complete the GATK IndelRealigner step.'
        raise Exception(message)

# TODO: Delete once this works.
# if args.debug < 1:
#     os.remove(path_duplicates_marked_bam)
#     os.remove(path_duplicates_marked_bai)
#     os.remove(path_duplicates_marked_md5)

# Run the GATK BaseRecalibrator step as a first-pass walker for the GATK PrintReads step.

if not os.path.exists(path_recalibration_table_pre):
    java_process = Executable(name='gatk_base_recalibrator', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='BaseRecalibrator')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    for file_path in known_sites_recalibration:
        sub_command.add_OptionLong(key='knownSites', value=file_path)
    sub_command.add_OptionLong(key='input_file', value=path_realigned_bam)
    sub_command.add_OptionLong(key='out', value=path_recalibration_table_pre)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK BaseRecalibrator step.'
        raise Exception(message)

# Run the GATK BaseRecalibrator on-the-fly recalibration step to generate plots.

if not os.path.exists(path_recalibration_table_post):
    java_process = Executable(name='gatk_base_recalibrator_post', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='BaseRecalibrator')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    for file_path in known_sites_recalibration:
        sub_command.add_OptionLong(key='knownSites', value=file_path)
    sub_command.add_OptionLong(key='BQSR', value=path_recalibration_table_pre)
    sub_command.add_OptionLong(key='input_file', value=path_realigned_bam)
    sub_command.add_OptionLong(key='out', value=path_recalibration_table_post)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK BaseRecalibrator post step.'
        raise Exception(message)

# Run the GATK AnalyzeCovariates step to create a recalibration plot.

if not os.path.exists(path_recalibration_plot):
    java_process = Executable(name='gatk_analyze_covariates', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='AnalyzeCovariates')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    sub_command.add_OptionLong(key='afterReportFile', value=path_recalibration_table_post)
    sub_command.add_OptionLong(key='beforeReportFile', value=path_recalibration_table_pre)
    sub_command.add_OptionLong(key='plotsReportFile', value=path_recalibration_plot)
    # sub_command.add_OptionLong(key='logging_level', value='DEBUG')

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK AnalyzeCovariates step.'
        raise Exception(message)

# Run the GATK PrintReads step as second-pass walker after the BaseRecalibrator step.

if not os.path.exists(path_recalibrated_bam):
    java_process = Executable(name='gatk_print_reads', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='PrintReads')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    sub_command.add_OptionLong(key='input_file', value=path_realigned_bam)
    sub_command.add_OptionLong(key='BQSR', value=path_recalibration_table_pre)
    sub_command.add_OptionLong(key='out', value=path_recalibrated_bam)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK PrintReads step.'
        raise Exception(message)

# TODO: Delete once this works.
# if args.debug < 1:
#     os.remove(path_realigned_bam)
#     os.remove(path_realigned_bai)

# Run the Picard CollectAlignmentSummaryMetrics step.

if not os.path.exists(path_alignment_summary_metrics):
    java_process = Executable(name='picard_collect_alignment_summary_metrics', program='java',
                              sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'CollectAlignmentSummaryMetrics.jar'))
    java_process.add_SwitchShort(key='Xmx6G')

    sub_command = java_process.sub_command
    sub_command.add_OptionPair(key='INPUT', value=path_recalibrated_bam)
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

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
