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

from Bio.BSF import Runnable


def run_executable(key):
    """Run an Executable defined in the Pickler dict.

    :param key: Key for the Executable
    :type key: str
    :return: Nothing
    :rtype: None
    """

    executable = pickler_dict[key]
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_picard_merge_sam_files():
    """Run the Picard MergeSamFiles step.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['merged_bam']):
        return

    run_executable(key='picard_merge_sam_files')


def run_picard_mark_duplicates():
    """Run the Picard MarkDuplicates step.

    Optical duplicates should already have been flagged in the lane-specific processing step.
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['duplicates_marked_bam']):
        return

    run_picard_merge_sam_files()

    # The Picard MarkDuplicates step may be skipped.
    if not 'picard_mark_duplicates' in pickler_dict:
        return

    run_executable(key='picard_mark_duplicates')

    if args.debug < 1:
        for file_key in ('merged_bam', 'merged_bai', 'merged_md5'):
            if os.path.exists(pickler_dict['file_path_dict'][file_key]):
                os.remove(pickler_dict['file_path_dict'][file_key])


def run_gatk_realigner_target_creator():
    """Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['realigner_targets']):
        return

    run_picard_mark_duplicates()
    run_executable(key='gatk_realigner_target_creator')


def run_gatk_indel_realigner():
    """Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['realigned_bam']):
        return

    run_gatk_realigner_target_creator()
    run_executable(key='gatk_indel_realigner')

    if args.debug < 1:
        # Remove file from the previous Picard MarkDuplicates step and additionally,
        # if the Picard MarkDuplicates step was skipped, remove the merged BAM files here.
        for file_key in ('duplicates_marked_bam', 'duplicates_marked_bai', 'duplicates_marked_md5',
                         'merged_bam', 'merged_bai', 'merged_md5'):
            if os.path.exists(pickler_dict['file_path_dict'][file_key]):
                os.remove(pickler_dict['file_path_dict'][file_key])


def run_picard_collect_alignment_summary_metrics():
    """Run the Picard CollectAlignmentSummaryMetrics step.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['alignment_summary_metrics']):
        return

    run_gatk_indel_realigner()
    run_executable(key='picard_collect_alignment_summary_metrics')


def run_gatk_haplotype_caller():
    """Run the GATK HaplotypeCaller per sample.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['raw_variants_gvcf_idx']):
        return

    run_picard_collect_alignment_summary_metrics()
    run_executable(key='gatk_haplotype_caller')


# Set the environment consistently.

os.environ['LANG'] = 'C'

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

# Create a temporary directory.

path_temporary = pickler_dict['file_path_dict']['temporary_sample']

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Run the chain of executables back up the function hierarchy so that
# dependencies on temporarily created files become simple to manage.

run_gatk_haplotype_caller()

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
