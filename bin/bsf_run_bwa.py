#! /usr/bin/env python
#
# BSF Python script to run the Burrows Wheeler Aligner BWA.
#
#
# Copyright 2013 - 2018 Michael K. Schuster
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
import pickle
import shutil

from bsf.process import Command, Executable
from bsf.standards import JavaClassPath

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

unpickler = pickle.Unpickler(file=pickler_file)

pickler_dict = unpickler.load()

pickler_file.close()

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
    raise Exception('The pickler dict in the pickler file needs to contain key {!r}.'.format(key))

key = 'classpath_picard'
if key in pickler_dict and pickler_dict[key]:
    classpath_picard = pickler_dict[key]
else:
    classpath_picard = JavaClassPath.get_picard()

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, replicate_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

path_fastq_1 = "{}{}_R1.fastq".format(prefix, replicate_key)
path_fastq_2 = "{}{}_R2.fastq".format(prefix, replicate_key)
path_aligned_sam = "{}{}_aligned.sam".format(prefix, replicate_key)
path_cleaned_sam = "{}{}_cleaned.sam".format(prefix, replicate_key)
path_header_sam = "{}{}_header.sam".format(prefix, replicate_key)
path_temporary_sam = "{}{}_temporary.sam".format(prefix, replicate_key)
path_sorted_bam = "{}{}.bam".format(prefix, replicate_key)

# SAM header lines that need propagating around FASTQ files. Sigh!
sam_header_pg = list()
""" @type sam_header_pg: list[str | unicode] """
sam_header_rg = list()
""" @type sam_header_rg: list[str | unicode] """

# Run BWA to produce an aligned SAM file.

run_bwa = pickler_dict['bwa_executable']
assert isinstance(run_bwa, Executable)
run_bwa.stdout_path = path_aligned_sam

# Check if the run_bwa sub-command got a BAM file to align in mem mode.
# If so, convert it first into FASTQ file(s).
# TODO: Comma-separated lists of FASTQ or BAM files are not yet supported,
# but they may not be supported by mem, either.
# TODO: This code only works if there is only one read group @RG per BAM file.
# It would be good to run SamToFastq to get one FASTQ file (pair) per read group,
# align them separately and set the correct read group for each.
# Picard writes out a file per platform unit PU, with special characters replaced by underscores and
# '_1.fastq', as well as '_2.fastq' appended.
# See function makeFileNameSafe
# https://github.com/samtools/htsjdk/blob/1.114/src/java/htsjdk/samtools/util/IOUtil.java#L719
# TODO: Once this works, it should be turned into a generic method useful in the ChIPSeq and RNASeq analyses, too.

if run_bwa.sub_command.program == 'mem' and run_bwa.sub_command.arguments[1][-4:] == '.bam':

    # Propagate SAM header lines @PG and @RG into the final BAM file.

    samtools = Executable(
        name='samtools_view',
        program='samtools',
        sub_command=Command(program='view'),
        stdout_path=path_temporary_sam)

    samtools_view = samtools.sub_command
    samtools_view.add_switch_short(key='H')
    samtools_view.arguments.append(run_bwa.sub_command.arguments[1])

    child_return_code = samtools.run()

    if child_return_code:
        raise Exception('Could not complete the {!r} step on the BAM file for the replicate.'.format(samtools.name))

    sam_temporary_handle = open(path_temporary_sam, 'r')
    for line in sam_temporary_handle:
        if line[:3] == '@PG':
            sam_header_pg.append(line.rstrip())
        if line[:3] == '@RG':
            sam_header_rg.append(line.rstrip())
    sam_temporary_handle.close()
    os.remove(path_temporary_sam)

    # At this stage, the SAM @PG and @RG lines are stored internally.
    # Now run Picard SamToFastq to convert.

    java_process = Executable(name='sam_to_fastq', program='java', sub_command=Command())
    java_process.add_switch_short(key='d64')
    java_process.add_switch_short(key='server')
    java_process.add_switch_short(key='Xmx4G')
    java_process.add_option_pair(key='-Djava.io.tmpdir', value=path_temporary)

    picard_process = java_process.sub_command
    picard_process.add_option_short(key='jar', value=os.path.join(classpath_picard, 'picard.jar'))
    picard_process.sub_command = Command(program='SamToFastq')

    sam_to_fastq = picard_process.sub_command
    sam_to_fastq.add_option_pair(key='INPUT', value=run_bwa.sub_command.arguments[1])
    sam_to_fastq.add_option_pair(key='FASTQ', value=path_fastq_1)
    sam_to_fastq.add_option_pair(key='SECOND_END_FASTQ', value=path_fastq_2)
    sam_to_fastq.add_option_pair(key='INCLUDE_NON_PF_READS', value='false')  # TODO: Make this configurable.
    sam_to_fastq.add_option_pair(key='TMP_DIR', value=path_temporary)
    sam_to_fastq.add_option_pair(key='VERBOSITY', value='WARNING')
    sam_to_fastq.add_option_pair(key='QUIET', value='false')
    sam_to_fastq.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')

    child_return_code = java_process.run()

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(java_process.name))

    # The Read Group (@RG) option needs inserting at the BWA stage,
    # otherwise reads are not associated with a read group in the aligned SAM file.

    # TODO: So far, this module only works with one read group per BAM file.
    if len(sam_header_rg) > 1:
        raise Exception(
            "BAM file {!r} contains more than one read group line, which is not supported, yet.\n"
            "RGs: {!r}".format(run_bwa.sub_command.arguments[1], sam_header_rg))

    run_bwa.sub_command.add_option_short(key='R', value=sam_header_rg[0].rstrip().replace("\t", "\\t"))

    # After the SamToFastq conversion, the second and third arguments in the BWA sub-command need replacing.

    run_bwa.sub_command.arguments = [run_bwa.sub_command.arguments[0]]
    run_bwa.sub_command.arguments.append(path_fastq_1)
    if os.path.getsize(filename=path_fastq_2):
        run_bwa.sub_command.arguments.append(path_fastq_2)

child_return_code = run_bwa.run()

if child_return_code:
    raise Exception('Could not complete the {!r} step.'.format(run_bwa.name))

# Remove the temporary, uncompressed FASTQ files if they exist.
if os.path.exists(path=path_fastq_1):
    os.remove(path_fastq_1)
if os.path.exists(path=path_fastq_2):
    os.remove(path_fastq_2)

# Run Picard CleanSam to convert the aligned SAM file into a cleaned SAM file.

java_process = Executable(name='clean_sam', program='java', sub_command=Command())
java_process.add_switch_short(key='d64')
java_process.add_switch_short(key='server')
java_process.add_switch_short(key='Xmx4G')
java_process.add_option_pair(key='-Djava.io.tmpdir', value=path_temporary)

picard_process = java_process.sub_command
picard_process.add_option_short(key='jar', value=os.path.join(classpath_picard, 'picard.jar'))
picard_process.sub_command = Command(program='CleanSam')

clean_sam = picard_process.sub_command
clean_sam.add_option_pair(key='INPUT', value=path_aligned_sam)
clean_sam.add_option_pair(key='OUTPUT', value=path_cleaned_sam)
clean_sam.add_option_pair(key='TMP_DIR', value=path_temporary)
clean_sam.add_option_pair(key='VERBOSITY', value='WARNING')
clean_sam.add_option_pair(key='QUIET', value='false')
clean_sam.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')

child_return_code = java_process.run()

if child_return_code:
    raise Exception('Could not complete the {!r} step.'.format(java_process.name))

# Retrofit the stored SAM header lines if required.

if len(sam_header_pg) or len(sam_header_rg):

    samtools = Executable(
        name='samtools_view',
        program='samtools',
        sub_command=Command(program='view'),
        stdout_path=path_temporary_sam)

    samtools_view = samtools.sub_command
    samtools_view.add_switch_short(key='H')
    samtools_view.add_switch_short(key='S')
    samtools_view.arguments.append(path_cleaned_sam)

    child_return_code = samtools.run()

    if child_return_code:
        raise Exception('Could not complete the {!r} step on the SAM file after CleanSam.'.format(samtools.name))

    sam_temporary_handle = open(path_temporary_sam, 'r')
    sam_header_handle = open(path_header_sam, 'w')
    for line in sam_temporary_handle:
        if line[:3] == '@PG':
            # Insert all @PG lines before this one, then clear the list so that no further insertion is possible.
            for line_pg in sam_header_pg:
                sam_header_handle.write(line_pg + "\n")
            sam_header_pg = list()
            """ @type sam_header_pg: list[str | unicode] """
        sam_header_handle.write(line)

    # Add remaining @PG lines it they were not present in the cleaned SAM file.

    for line_pg in sam_header_pg:
        sam_header_handle.write(line_pg + "\n")
    sam_header_pg = list()
    """ @type sam_header_pg: list[str | unicode] """

    sam_header_handle.close()
    sam_temporary_handle.close()

    os.remove(path_temporary_sam)

    # Run Picard ReplaceSamHeader.

    java_process = Executable(name='replace_sam_header', program='java', sub_command=Command())
    java_process.add_switch_short(key='d64')
    java_process.add_switch_short(key='server')
    java_process.add_switch_short(key='Xmx4G')
    java_process.add_option_pair(key='-Djava.io.tmpdir', value=path_temporary)

    picard_process = java_process.sub_command
    picard_process.add_option_short(key='jar', value=os.path.join(classpath_picard, 'picard.jar'))
    picard_process.sub_command = Command(program='ReplaceSamHeader')

    replace_sam_header = picard_process.sub_command
    replace_sam_header.add_option_pair(key='INPUT', value=path_cleaned_sam)
    replace_sam_header.add_option_pair(key='HEADER', value=path_header_sam)
    replace_sam_header.add_option_pair(key='OUTPUT', value=path_temporary_sam)
    replace_sam_header.add_option_pair(key='TMP_DIR', value=path_temporary)
    replace_sam_header.add_option_pair(key='VERBOSITY', value='WARNING')
    replace_sam_header.add_option_pair(key='QUIET', value='false')
    replace_sam_header.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')

    child_return_code = java_process.run()

    if child_return_code:
        raise Exception('Could not complete the {!r} step.'.format(java_process.name))

    # Remove the now redundant cleaned SAM file.
    os.remove(path_cleaned_sam)
    # Rename the temporary SAM file that has the new header into the cleaned SAM file with the original header.
    os.rename(path_temporary_sam, path_cleaned_sam)
    # Remove the now redundant header file.
    os.remove(path_header_sam)

# Remove the aligned SAM file.

if os.path.exists(path=path_aligned_sam):
    os.remove(path_aligned_sam)

# Run Picard SortSam to convert the cleaned SAM file into a coordinate sorted BAM file.

java_process = Executable(name='sort_sam', program='java', sub_command=Command())
java_process.add_switch_short(key='d64')
java_process.add_switch_short(key='server')
java_process.add_switch_short(key='Xmx4G')
java_process.add_option_pair(key='-Djava.io.tmpdir', value=path_temporary)

picard_process = java_process.sub_command
picard_process.add_option_short(key='jar', value=os.path.join(classpath_picard, 'picard.jar'))
picard_process.sub_command = Command(program='SortSam')

sort_sam = picard_process.sub_command
sort_sam.add_option_pair(key='INPUT', value=path_cleaned_sam)
sort_sam.add_option_pair(key='OUTPUT', value=path_sorted_bam)
sort_sam.add_option_pair(key='SORT_ORDER', value='coordinate')
sort_sam.add_option_pair(key='TMP_DIR', value=path_temporary)
sort_sam.add_option_pair(key='VERBOSITY', value='WARNING')
sort_sam.add_option_pair(key='QUIET', value='false')
sort_sam.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')
sort_sam.add_option_pair(key='COMPRESSION_LEVEL', value='9')
sort_sam.add_option_pair(key='MAX_RECORDS_IN_RAM', value='2000000')
sort_sam.add_option_pair(key='CREATE_INDEX', value='true')
sort_sam.add_option_pair(key='CREATE_MD5_FILE', value='true')

child_return_code = java_process.run()

if child_return_code:
    raise Exception('Could not complete the {!r} step.'.format(java_process.name))

# Remove the cleaned SAM file.

if os.path.exists(path=path_cleaned_sam):
    os.remove(path_cleaned_sam)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Write a status file.

status_path = "{}{}_completed.txt".format(prefix, replicate_key)
open(status_path, 'w').close()

# Job done.
