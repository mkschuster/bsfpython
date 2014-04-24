#! /usr/bin/env python
#
# BSF Python script to process cohorts for the GATK pipeline.
#
#   Picard MarkDuplicates
#   GATK RealignerTargetCreator
#   GATK IndelRealigner
#   GATK BaseRecalibrator
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

key = 'cohort_key'
if key in pickler_dict and pickler_dict[key]:
    cohort_key = pickler_dict[key]
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

# The known_sites_recalibration parameter is not required here.

key = 'known_sites_discovery'
if key in pickler_dict and pickler_dict[key]:
    known_sites_discovery = pickler_dict[key]
else:
    message = "The pickler dict in the pickler file needs to contain a '{}' entry.".format(key)
    raise Exception(message)

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, cohort_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

path_combined_gvcf = '{}{}_combined_gvcf.vcf'.format(prefix, cohort_key)
path_genotyped_gvcf = '{}{}_genotyped_gvcf.vcf'.format(prefix, cohort_key)

# Run the GATK CombineGVCFs step.

if not os.path.exists(path_combined_gvcf):
    java_process = Executable(name='gatk_combine_gvcfs', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='CombineGVCFs')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    sub_command.add_OptionLong(key='variant', value=path_combined_gvcf)
    sub_command.add_OptionLong(key='out', value=path_combined_gvcf)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK CombineGVCFs step.'
        raise Exception(message)

# Run the GATK GenotypeGVCFs step.

if not os.path.exists(path_combined_gvcf):
    java_process = Executable(name='gatk_combine_gvcfs', program='java', sub_command=Command(command=str()))
    java_process.add_SwitchShort(key='d64')
    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
    java_process.add_SwitchShort(key='Xmx6G')
    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=path_temporary)

    sub_command = java_process.sub_command
    sub_command.add_OptionLong(key='analysis_type', value='GenotypeGVCFs')
    sub_command.add_OptionLong(key='reference_sequence', value=path_reference_sequence)
    sub_command.add_OptionLong(key='dbsnp', value=known_sites_discovery)
    for file_path in replicate_file_list:
        sub_command.add_OptionLong(key='variant', value=file_path)
    sub_command.add_OptionLong(key='out', value=path_combined_gvcf)

    child_return_code = Runnable.run(executable=java_process)

    if child_return_code:
        message = 'Could not complete the GATK GenotypeGVCFs step.'
        raise Exception(message)

# Run the GATK

# TODO: Complete the VQSLOD step.

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
