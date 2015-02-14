"""bsf.runnables.variant_calling_process_lane

A package of classes and methods to run variant calling processes.
"""

#
# Copyright 2013 - 2014 Michael K. Schuster
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
import os
import shutil

from bsf import Runnable


def run_picard_mark_duplicates(runnable):
    """Run the I{Picard MarkDuplicates} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['duplicates_marked_bam']):
        return

    # The Picard MarkDuplicates step may be skipped.
    if not 'picard_mark_duplicates' in runnable.executable_dict:
        return

    runnable.run_executable(name='picard_mark_duplicates')

    # It may not be the best idea to remove the aligned BAM file from the previous lane-specific alignment step here.
    # For the moment, keep pipeline steps independent from each other.
    # if runnable.debug < 1:
    #     for file_key in ('aligned_bam', 'aligned_bai', 'aligned_md5'):
    #         if os.path.exists(runnable.file_path_dict[file_key]):
    #             os.remove(runnable.file_path_dict[file_key])


def run_gatk_realigner_target_creator(runnable):
    """Run the I{GATK RealignerTargetCreator} step as the first-pass walker for the I{GATK IndelRealigner} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['realigner_targets']):
        return

    run_picard_mark_duplicates(runnable=runnable)
    runnable.run_executable(name='gatk_realigner_target_creator')


def run_gatk_indel_realigner(runnable):
    """Run the I{GATK IndelRealigner} step as a second-pass walker after the I{GATK RealignerTargetCreator} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['realigned_bam']):
        return

    run_gatk_realigner_target_creator(runnable=runnable)
    runnable.run_executable(name='gatk_indel_realigner')

    if runnable.debug < 1:
        for file_key in ('duplicates_marked_bam', 'duplicates_marked_bai', 'duplicates_marked_md5'):
            if os.path.exists(runnable.file_path_dict[file_key]):
                os.remove(runnable.file_path_dict[file_key])


def run_gatk_base_recalibrator_pre(runnable):
    """Run the I{GATK BaseRecalibrator} step as a first-pass walker for the I{GATK PrintReads} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibration_table_pre']):
        return

    run_gatk_indel_realigner(runnable=runnable)
    runnable.run_executable(name='gatk_base_recalibrator_pre')


def run_gatk_base_recalibrator_post(runnable):
    """Run the I{GATK BaseRecalibrator} on-the-fly recalibration step to generate plots.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibration_table_post']):
        return

    run_gatk_base_recalibrator_pre(runnable=runnable)
    runnable.run_executable(name='gatk_base_recalibrator_post')


def run_gatk_analyze_covariates(runnable):
    """Run the I{GATK AnalyzeCovariates} step to create a recalibration plot.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibration_plot']):
        return

    run_gatk_base_recalibrator_post(runnable=runnable)
    runnable.run_executable(name='gatk_analyze_covariates')


def run_gatk_print_reads(runnable):
    """Run the I{GATK PrintReads} step as second-pass walker after the I{GATK BaseRecalibrator} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibrated_bam']):
        return

    run_gatk_analyze_covariates(runnable=runnable)
    runnable.run_executable(name='gatk_print_reads')

    if runnable.debug < 1:
        for file_key in ('realigned_bam', 'realigned_bai'):
            if os.path.exists(runnable.file_path_dict[file_key]):
                os.remove(runnable.file_path_dict[file_key])


def run_picard_collect_alignment_summary_metrics(runnable):
    """Run the I{Picard CollectAlignmentSummaryMetrics} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['alignment_summary_metrics']):
        return

    run_gatk_print_reads(runnable=runnable)
    runnable.run_executable(name='picard_collect_alignment_summary_metrics')


def run(runnable):
    """Run the the C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    # Create a temporary directory.

    path_temporary = runnable.file_path_dict['temporary_directory']

    if not os.path.isdir(path_temporary):
        try:
            os.makedirs(path_temporary)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Run the chain of executables back up the function hierarchy so that
    # dependencies on temporarily created files become simple to manage.

    run_picard_collect_alignment_summary_metrics(runnable=runnable)

    # Remove the temporary directory and everything within it.

    shutil.rmtree(path=path_temporary, ignore_errors=False)

    # Job done.
