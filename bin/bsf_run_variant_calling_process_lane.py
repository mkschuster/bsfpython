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

from Bio.BSF import Runnable


def run_picard_mark_duplicates(pickler_dict):
    """Run the Picard MarkDuplicates step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['duplicates_marked_bam']):
        return

    # The Picard MarkDuplicates step may be skipped.
    if not 'picard_mark_duplicates' in pickler_dict:
        return

    executable = pickler_dict['picard_mark_duplicates']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))

    # It may not be the best idea to remove the aligned BAM file from the previous lane-specific alignment step here.
    # For the moment, keep pipeline steps independent from each other.
    # if args.debug < 1:
    #     if os.path.exists(path=pickler_dict['file_path_dict']['aligned_bam']):
    #         os.remove(pickler_dict['file_path_dict']['aligned_bam'])
    #     if os.path.exists(path=pickler_dict['file_path_dict']['aligned_bai']):
    #         os.remove(pickler_dict['file_path_dict']['aligned_bai'])
    #     if os.path.exists(path=pickler_dict['file_path_dict']['aligned_md5']):
    #         os.remove(pickler_dict['file_path_dict']['aligned_md5'])


def run_gatk_realigner_target_creator(pickler_dict):
    """Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['realigner_targets']):
        return

    run_picard_mark_duplicates(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_realigner_target_creator']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_indel_realigner(pickler_dict):
    """Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['realigned_bam']):
        return

    run_gatk_realigner_target_creator(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_indel_realigner']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))

    if args.debug < 1:
        if os.path.exists(path=pickler_dict['file_path_dict']['duplicates_marked_bam']):
            os.remove(pickler_dict['file_path_dict']['duplicates_marked_bam'])
        if os.path.exists(path=pickler_dict['file_path_dict']['duplicates_marked_bai']):
            os.remove(pickler_dict['file_path_dict']['duplicates_marked_bai'])
        if os.path.exists(path=pickler_dict['file_path_dict']['duplicates_marked_md5']):
            os.remove(pickler_dict['file_path_dict']['duplicates_marked_md5'])


def run_gatk_base_recalibrator_pre(pickler_dict):
    """Run the GATK BaseRecalibrator step as a first-pass walker for the GATK PrintReads step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['recalibration_table_pre']):
        return

    run_gatk_indel_realigner(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_base_recalibrator_pre']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_base_recalibrator_post(pickler_dict):
    """Run the GATK BaseRecalibrator on-the-fly recalibration step to generate plots.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['recalibration_table_post']):
        return

    run_gatk_base_recalibrator_pre(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_base_recalibrator_post']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_analyze_covariates(pickler_dict):
    """Run the GATK AnalyzeCovariates step to create a recalibration plot.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['recalibration_plot']):
        return

    run_gatk_base_recalibrator_post(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_analyze_covariates']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_print_reads(pickler_dict):
    """Run the GATK PrintReads step as second-pass walker after the BaseRecalibrator step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['recalibrated_bam']):
        return

    run_gatk_analyze_covariates(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_print_reads']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))

    if args.debug < 1:
        if os.path.exists(path=pickler_dict['file_path_dict']['realigned_bam']):
            os.remove(pickler_dict['file_path_dict']['realigned_bam'])
        if os.path.exists(path=pickler_dict['file_path_dict']['realigned_bai']):
            os.remove(pickler_dict['file_path_dict']['realigned_bai'])


def run_picard_collect_alignment_summary_metrics(pickler_dict):
    """Run the Picard CollectAlignmentSummaryMetrics step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['alignment_summary_metrics']):
        return

    run_gatk_print_reads(pickler_dict=pickler_dict)

    executable = pickler_dict['picard_collect_alignment_summary_metrics']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


# Set the environment consistently.

os.environ['LANG'] = 'C'

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

# Create a temporary directory.

path_temporary = pickler_dict['file_path_dict']['temporary_lane']

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Run the chain of executables back up the function hierarchy so that
# dependencies on temporarily created files become simple to manage.

run_picard_collect_alignment_summary_metrics(pickler_dict=pickler_dict)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
