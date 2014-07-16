#! /usr/bin/env python
#
# BSF Python script to process cohorts for the GATK pipeline.
#
#   GATK CombineGVCFs
#   GATK GenotypeGVCFs
#   GATK VariantRecalibrator for SNPs
#   GATK VariantRecalibrator for INDELs
#   GATK ApplyRecalibration for SNPs
#   GATK ApplyRecalibration for INDELs
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


def run_gatk_combine_gvcfs(pickler_dict):
    """Run the GATK CombineGVCFs step.

    It is only required for hierarchically merging samples before GenotypeGVCFs,
    if too many samples need processing.
    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['gvcf_combined']):
        return

    executable = pickler_dict['gatk_combine_gvcfs']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_genotype_gvcfs(pickler_dict):
    """Run the GATK GenotypeGVCFs step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['gvcf_genotyped_raw']):
        return

    run_gatk_combine_gvcfs(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_genotype_gvcfs']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_variant_recalibrator_snp(pickler_dict):
    """Run the GATK VariantRecalibrator for SNPs.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['recalibration_snp']):
        return

    run_gatk_genotype_gvcfs(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_variant_recalibrator_snp']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_variant_recalibrator_indel(pickler_dict):
    """Run the GATK VariantRecalibrator for INDELs.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['recalibration_indel']):
        return

    run_gatk_variant_recalibrator_snp(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_variant_recalibrator_indel']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_apply_recalibration_snp(pickler_dict):
    """Run the GATK ApplyRecalibration step for SNPs.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['gvcf_recalibrated_snp_raw_indel']):
        return

    run_gatk_variant_recalibrator_indel(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_apply_recalibration_snp']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_apply_recalibration_indel(pickler_dict):
    """Run the GATK ApplyRecalibration step for INDELs.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['gvcf_recalibrated_snp_recalibrated_indel']):
        return

    run_gatk_apply_recalibration_snp(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_apply_recalibration_indel']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_snpeff(pickler_dict):
    """Run the snpEff tool.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['snpeff_annotated']):
        return

    run_gatk_apply_recalibration_indel(pickler_dict=pickler_dict)

    executable = pickler_dict['snpeff']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_variant_annotator(pickler_dict):
    """Run the GATK VariantAnnotator step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['gvcf_annotated']):
        return

    run_snpeff(pickler_dict=pickler_dict)

    executable = pickler_dict['gatk_variant_annotator']
    child_return_code = Runnable.run(executable=executable)

    if child_return_code:
        raise Exception('Could not complete the {}!r step.'.format(executable.name))


def run_gatk_select_variants(pickler_dict):
    """Run the GATK SelectVariants step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    # Check that all files of this function have already been created.

    complete = True
    for sample_name in pickler_dict['sample_names']:
        if not os.path.exists(pickler_dict['file_path_dict']['sample_vcf_' + sample_name]):
            complete = False
            break
    if complete:
        return

    run_gatk_variant_annotator(pickler_dict=pickler_dict)

    # The GATK SelectVariants step has to be run for each sample name separately.

    for sample_name in pickler_dict['sample_names']:

        if os.path.exists(pickler_dict['file_path_dict']['sample_vcf_' + sample_name]):
            continue

        executable = pickler_dict['gatk_select_variants_sample_' + sample_name]
        child_return_code = Runnable.run(executable=executable)

        if child_return_code:
            raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_gatk_variants_to_table(pickler_dict):
    """Run the GATK VariantsToTable step.

    :param pickler_dict: Pickler dict
    :type pickler_dict: dict
    :return: Nothing
    :rtype: None
    """

    # Check that all files of this function have already been created.

    complete = True
    for sample_name in pickler_dict['sample_names']:
        if not os.path.exists(pickler_dict['file_path_dict']['sample_csv_' + sample_name]):
            complete = False
            break
    if complete:
        return

    run_gatk_select_variants(pickler_dict=pickler_dict)

    # The GATK SelectVariants step has to be run for each sample name separately.

    for sample_name in pickler_dict['sample_names']:

        if os.path.exists(pickler_dict['file_path_dict']['sample_csv_' + sample_name]):
            continue

        executable = pickler_dict['gatk_variants_to_table_sample_' + sample_name]
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

path_temporary = pickler_dict['file_path_dict']['temporary_cohort']

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Run the chain of executables back up the function hierarchy so that
# dependencies on temporarily created files become simple to manage.

run_gatk_variants_to_table(pickler_dict=pickler_dict)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
