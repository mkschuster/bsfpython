"""Bio.BSF.Runnables.variant_calling_process_sample

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

from Bio.BSF import Runnable


def run_picard_merge_sam_files(runnable):
    """Run the Picard MergeSamFiles step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['merged_bam']):
        return

    runnable.run_executable(name='picard_merge_sam_files')


def run_picard_mark_duplicates(runnable):
    """Run the Picard MarkDuplicates step.

    Optical duplicates should already have been flagged in the lane-specific processing step.
    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['duplicates_marked_bam']):
        return

    run_picard_merge_sam_files(runnable=runnable)

    # The Picard MarkDuplicates step may be skipped.
    if not 'picard_mark_duplicates' in runnable.executable_dict:
        return

    runnable.run_executable(name='picard_mark_duplicates')

    if runnable.debug < 1:
        for file_key in ('merged_bam', 'merged_bai', 'merged_md5'):
            if os.path.exists(runnable.file_path_dict[file_key]):
                os.remove(runnable.file_path_dict[file_key])


def run_gatk_realigner_target_creator(runnable):
    """Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['realigner_targets']):
        return

    run_picard_mark_duplicates(runnable=runnable)
    runnable.run_executable(name='gatk_realigner_target_creator')


def run_gatk_indel_realigner(runnable):
    """Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['realigned_bam']):
        return

    run_gatk_realigner_target_creator(runnable=runnable)
    runnable.run_executable(name='gatk_indel_realigner')

    if runnable.debug < 1:
        # Remove file from the previous Picard MarkDuplicates step and additionally,
        # if the Picard MarkDuplicates step was skipped, remove the merged BAM files here.
        for file_key in ('duplicates_marked_bam', 'duplicates_marked_bai', 'duplicates_marked_md5',
                         'merged_bam', 'merged_bai', 'merged_md5'):
            if os.path.exists(runnable.file_path_dict[file_key]):
                os.remove(runnable.file_path_dict[file_key])


def run_picard_collect_alignment_summary_metrics(runnable):
    """Run the Picard CollectAlignmentSummaryMetrics step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['alignment_summary_metrics']):
        return

    run_gatk_indel_realigner(runnable=runnable)
    runnable.run_executable(name='picard_collect_alignment_summary_metrics')


def run_gatk_haplotype_caller(runnable):
    """Run the GATK HaplotypeCaller per sample.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['raw_variants_gvcf_idx']):
        return

    run_picard_collect_alignment_summary_metrics(runnable=runnable)
    runnable.run_executable(name='gatk_haplotype_caller')


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

    run_gatk_haplotype_caller(runnable=runnable)

    # Remove the temporary directory and everything within it.

    shutil.rmtree(path=path_temporary, ignore_errors=False)

    # Job done.
