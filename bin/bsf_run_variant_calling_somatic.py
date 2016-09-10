#! /usr/bin/env python
#
# BSF Python script to run MuTect.
#
#
# Copyright 2013 - 2016 Michael K. Schuster
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
import os
from pickle import Unpickler
import shutil


def run_executable(key):
    """Run an Executable defined in the Pickler dict.

    :param key: Key for the Executable
    :type key: str
    :return: Nothing
    :rtype: None
    """

    executable = pickler_dict[key]
    child_return_code = executable.run()

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(executable.name))


def run_mutect():
    """Run MuTect.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['mutect_vcf']) \
            and os.path.getsize(pickler_dict['file_path_dict']['mutect_vcf']):
        return

    run_executable(key='mutect')


def run_indel_genotyper():
    """Run the Indel Genotyper

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['indel_idx']) \
            and os.path.getsize(pickler_dict['file_path_dict']['indel_idx']):
        return

    run_mutect()
    run_executable(key='indel_genotyper')


def run_gatk_combine_variants():
    """Run the GATK Combine Variants


    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['combined_idx']) \
            and os.path.getsize(pickler_dict['file_path_dict']['combined_idx']):
        return

    run_indel_genotyper()
    run_executable(key='gatk_combine_variants')

    # if args.debug < 1:
    #     for file_key in ('mutect_vcf', 'mutect_idx', 'indel_vcf', 'indel_idx'):
    #         if os.path.exists(pickler_dict['file_path_dict'][file_key]):
    #             os.remove(pickler_dict['file_path_dict'][file_key])


def run_snpeff():
    """Run the snpEff tool.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['snpeff_vcf']) \
            and os.path.getsize(pickler_dict['file_path_dict']['snpeff_vcf']):
        return

    run_gatk_combine_variants()
    run_executable(key='snpeff')


def run_gatk_variant_annotator():
    """Run the GATK VariantAnnotator step.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['annotated_idx']) \
            and os.path.getsize(pickler_dict['file_path_dict']['annotated_idx']):
        return

    run_snpeff()
    run_executable(key='gatk_variant_annotator')

    # if args.debug < 1:
    #     for file_key in ('combined_vcf', 'combined_idx'):
    #         if os.path.exists(pickler_dict['file_path_dict'][file_key]):
    #             os.remove(pickler_dict['file_path_dict'][file_key])


def run_gatk_variant_to_table():
    """Run the GATK VariantsToTable step.

    :return: Nothing
    :rtype: None
    """

    if os.path.exists(pickler_dict['file_path_dict']['annotated_tsv']) \
            and os.path.getsize(pickler_dict['file_path_dict']['annotated_tsv']):
        return

    run_gatk_variant_annotator()
    run_executable(key='gatk_variants_to_table')


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

parser = argparse.ArgumentParser(
    description='BSF Runner for running the MuTect somatic variant caller.')

parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

parser.add_argument(
    '--pickler_path',
    help='file path to a Python Pickler file',
    required=True)

args = parser.parse_args()

# Unpickle the file into a Python dict object.

pickler_file = open(args.pickler_path, 'rb')

unpickler = Unpickler(file=pickler_file)

pickler_dict = unpickler.load()

pickler_file.close()

# Create a temporary directory.

path_temporary = pickler_dict['file_path_dict']['temporary_directory']

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# Run the chain of executables back up the function hierarchy so that
# dependencies on temporarily created files become simple to manage.

run_gatk_variant_to_table()

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
