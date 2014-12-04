"""bsf.runnables.variant_calling_process_cohort

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


sample_names = list()


def run_gatk_combine_gvcfs(runnable):
    """Run the I{GATK CombineGVCFs} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['combined_gvcf_idx']):
        return

    runnable.run_executable(name='gatk_combine_gvcfs')


def run_gatk_combine_gvcfs_accessory(runnable):
    """Run the I{GATK CombineGVCFs} step on accessory GVCF files.

    It is only required for hierarchically merging samples before GenotypeGVCFs,
    if accessory samples for recalibration need processing.
    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['temporary_gvcf_idx']):
        return

    run_gatk_combine_gvcfs(runnable=runnable)
    if 'gatk_combine_gvcfs_accessory' in runnable.executable_dict:
        runnable.run_executable(name='gatk_combine_gvcfs_accessory')


def run_gatk_genotype_gvcfs(runnable):
    """Run the I{GATK GenotypeGVCFs} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['genotyped_raw_idx']):
        return

    run_gatk_combine_gvcfs_accessory(runnable=runnable)
    runnable.run_executable(name='gatk_genotype_gvcfs')


def run_gatk_variant_recalibrator_snp(runnable):
    """Run the I{GATK VariantRecalibrator} for I{SNPs}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibration_snp']):
        return

    run_gatk_genotype_gvcfs(runnable=runnable)
    runnable.run_executable(name='gatk_variant_recalibrator_snp')


def run_gatk_variant_recalibrator_indel(runnable):
    """Run the I{GATK VariantRecalibrator} for I{INDELs}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibration_indel']):
        return

    run_gatk_variant_recalibrator_snp(runnable=runnable)
    runnable.run_executable(name='gatk_variant_recalibrator_indel')


def run_gatk_apply_recalibration_snp(runnable):
    """Run the I{GATK ApplyRecalibration} step for I{SNPs}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibrated_snp_raw_indel_idx']):
        return

    run_gatk_variant_recalibrator_indel(runnable=runnable)
    runnable.run_executable(name='gatk_apply_recalibration_snp')


def run_gatk_apply_recalibration_indel(runnable):
    """Run the I{GATK ApplyRecalibration} step for I{INDELs}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['recalibrated_snp_recalibrated_indel_idx']):
        return

    run_gatk_apply_recalibration_snp(runnable=runnable)
    runnable.run_executable(name='gatk_apply_recalibration_indel')


def run_gatk_select_variants_cohort(runnable):
    """Run the I{GATK SelectVariants} step only, if accessory GVCF files have been used.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['multi_sample_idx']):
        return

    run_gatk_apply_recalibration_indel(runnable=runnable)
    if 'gatk_select_variants_cohort' in runnable.executable_dict:
        runnable.run_executable(name='gatk_select_variants_cohort')


def run_snpeff(runnable):
    """Run the I{snpEff} tool.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['snpeff_vcf']):
        return

    run_gatk_select_variants_cohort(runnable=runnable)
    runnable.run_executable(name='snpeff')


def run_gatk_variant_annotator(runnable):
    """Run the I{GATK VariantAnnotator} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    if os.path.exists(runnable.file_path_dict['annotated_idx']):
        return

    run_snpeff(runnable=runnable)
    runnable.run_executable(name='gatk_variant_annotator')


def run_gatk_select_variants(runnable):
    """Run the I{GATK SelectVariants} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    # Check that all files of this function have already been created.

    complete = True
    for sample_name in sample_names:
        if not os.path.exists(runnable.file_path_dict['sample_vcf_' + sample_name]):
            complete = False
            break
    if complete:
        return

    run_gatk_variant_annotator(runnable=runnable)

    # The GATK SelectVariants step has to be run for each sample name separately.

    for sample_name in sample_names:
        if os.path.exists(runnable.file_path_dict['sample_vcf_' + sample_name]):
            continue
        runnable.run_executable(name='gatk_select_variants_sample_' + sample_name)


def run_gatk_variants_to_table(runnable):
    """Run the I{GATK VariantsToTable} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    # Check that all files of this function have already been created.

    complete = True
    for sample_name in sample_names:
        if not os.path.exists(runnable.file_path_dict['sample_csv_' + sample_name]):
            complete = False
            break
    if complete:
        return

    run_gatk_select_variants(runnable=runnable)

    # The GATK SelectVariants step has to be run for each sample name separately.

    for sample_name in sample_names:
        if os.path.exists(runnable.file_path_dict['sample_csv_' + sample_name]):
            continue
        runnable.run_executable(name='gatk_variants_to_table_sample_' + sample_name)


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

    # Get all sample names, from file_path_dict keys that start with 'sample_vcf_'.

    keys = runnable.file_path_dict.keys()
    keys.sort(cmp=lambda x, y: cmp(x.name, y.name))

    for key in keys:
        if key[:11] == 'sample_vcf_':
            sample_names.append(key[11:])

    # Run the chain of executables back up the function hierarchy so that
    # dependencies on temporarily created files become simple to manage.

    run_gatk_variants_to_table(runnable=runnable)

    # Remove the temporary directory and everything within it.

    shutil.rmtree(path=path_temporary, ignore_errors=False)

    # Job done.
