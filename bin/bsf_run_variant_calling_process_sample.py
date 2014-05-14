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


def run_picard_merge_sam_files(pickler_dict):
    """Run the Picard MergeSamFiles step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['merged_bam']):
        return

    executable = pickler_dict['picard_merge_sam_files']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        message = "Could not complete the '{}' step.".format(executable.name)
        raise Exception(message)


def run_picard_mark_duplicates(pickler_dict):
    """Run the Picard MarkDuplicates step.

    Optical duplicates should already have been flagged in the lane-specific processing step.
    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    run_picard_merge_sam_files(pickler_dict=pickler_dict)

    if os.path.exists(pickler_dict['file_path_dict']['duplicates_marked_bam']):
        return

    executable = pickler_dict['picard_mark_duplicates']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        message = "Could not complete the '{}' step.".format(executable.name)
        raise Exception(message)

    # if args.debug < 1:
    #     os.remove(pickler_dict['file_path_dict']['merged_bam'])
    #     os.remove(pickler_dict['file_path_dict']['merged_bai'])
    #     os.remove(pickler_dict['file_path_dict']['merged_md5'])


def run_gatk_realigner_target_creator(pickler_dict):
    """Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    run_picard_mark_duplicates(pickler_dict=pickler_dict)

    if os.path.exists(pickler_dict['file_path_dict']['realigner_targets']):
        return

    executable = pickler_dict['gatk_realigner_target_creator']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        message = "Could not complete the '{}' step.".format(executable.name)
        raise Exception(message)


def run_gatk_indel_realigner(pickler_dict):
    """Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    run_gatk_realigner_target_creator(pickler_dict=pickler_dict)

    if os.path.exists(pickler_dict['file_path_dict']['realigned_bam']):
        return

    executable = pickler_dict['gatk_indel_realigner']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        message = "Could not complete the '{}' step.".format(executable.name)
        raise Exception(message)

    # if args.debug < 1:
    #     os.remove(pickler_dict['file_path_dict']['duplicates_marked_bam'])
    #     os.remove(pickler_dict['file_path_dict']['duplicates_marked_bai'])
    #     os.remove(pickler_dict['file_path_dict']['duplicates_marked_md5'])


def run_picard_collect_alignment_summary_metrics(pickler_dict):
    """Run the Picard CollectAlignmentSummaryMetrics step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    run_gatk_indel_realigner(pickler_dict=pickler_dict)

    if os.path.exists(pickler_dict['file_path_dict']['alignment_summary_metrics']):
        return

    executable = pickler_dict['picard_collect_alignment_summary_metrics']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        message = "Could not complete the '{}' step.".format(executable.name)
        raise Exception(message)


def run_gatk_haplotype_caller(pickler_dict):
    """Run the GATK HaplotypeCaller per sample.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    run_picard_collect_alignment_summary_metrics(pickler_dict=pickler_dict)

    if os.path.exists(pickler_dict['file_path_dict']['raw_variants']):
        return

    executable = pickler_dict['gatk_haplotype_caller']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        message = "Could not complete the '{}' step.".format(executable.name)
        raise Exception(message)


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

run_gatk_haplotype_caller(pickler_dict=pickler_dict)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
