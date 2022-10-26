#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
#
#  Biomedical Sequencing Facility (BSF), part of the genomics core facility
#  of the Research Center for Molecular Medicine (CeMM) of the
#  Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
#  This file is part of BSF Python.
#
#  BSF Python is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BSF Python is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.
#
#
#  BSF Python script to run Tophat.
#
import errno
import os
import pickle
import shutil
from argparse import ArgumentParser
from typing import List

from bsf.connector import ConnectorFile
from bsf.process import Command, Executable, RunnableStep
from bsf.standards import SafeFileName, JavaArchive


def run_picard_sam_to_fastq(input_path, temporary_path):
    """Convert a BAM file into a pair of FASTQ files for each read group (@RG).

    :param input_path: A BAM file path.
    :type input_path: str
    :param temporary_path: A temporary directory path.
    :type temporary_path: str
    :return: A Python :py:class:`list` object of :py:class:`tuple` objects of
        Python :py:class:`str` (file path read 1) and Python :py:class:`str` (file path read2) objects.
    :rtype: list[(str, str)]
    """

    file_paths = list()

    # Propagate SAM header lines @PG and @RG into the final BAM file.
    # Create a temporary SAM file of the same name to store the SAM header.

    path_temporary_sam = os.path.basename(input_path)
    path_temporary_sam = path_temporary_sam.replace('.bam', '.sam')

    executable_samtools = Executable(
        name='samtools_view',
        program='samtools',
        sub_command=Command(program='view'),
        stdout=ConnectorFile(file_path=path_temporary_sam, file_mode='wt'))

    samtools_view = executable_samtools.sub_command
    samtools_view.add_switch_short(key='H')
    samtools_view.arguments.append(input_path)

    local_exception_str_list = executable_samtools.run()

    if local_exception_str_list:
        raise Exception('\n'.join(local_exception_str_list))

    # SAM header lines that need propagating around FASTQ files. Sigh!
    sam_header_pg: List[str] = list()
    sam_header_rg: List[str] = list()

    with open(file=path_temporary_sam, mode='rt') as _input_file:
        for line_str in _input_file:
            if line_str.startswith('@PG'):
                sam_header_pg.append(line_str.rstrip())
            if line_str.startswith('@RG'):
                sam_header_rg.append(line_str.rstrip())

    os.remove(path_temporary_sam)

    # At this stage, the SAM @PG and @RG lines are stored internally.
    # Now run Picard SamToFastq to convert.

    executable_java = Executable(name='sam_to_fastq', program='java', sub_command=Command())
    executable_java.add_switch_short(key='server')
    executable_java.add_switch_short(key='Xmx4G')
    executable_java.add_switch_short(key='XX:+UseG1GC')
    executable_java.add_option_pair_short(key='XX:ActiveProcessorCount', value='8')
    executable_java.add_option_pair_short(key='XX:CICompilerCount', value='2')
    executable_java.add_option_pair_short(key='XX:ParallelGCThreads', value='8')
    executable_java.add_option_pair_short(key='Djava.io.tmpdir', value=path_temporary)

    picard_process = executable_java.sub_command
    picard_process.add_option_short(key='jar', value=java_archive_picard)
    picard_process.sub_command = Command(program='SamToFastq')

    sam_to_fastq = picard_process.sub_command
    sam_to_fastq.add_option_pair(key='INPUT', value=input_path)
    sam_to_fastq.add_option_pair(key='OUTPUT_PER_RG', value='true')
    sam_to_fastq.add_option_pair(key='OUTPUT_DIR', value=temporary_path)
    sam_to_fastq.add_option_pair(key='INCLUDE_NON_PF_READS', value='false')  # TODO: Make this configurable.
    sam_to_fastq.add_option_pair(key='TMP_DIR', value=path_temporary)
    sam_to_fastq.add_option_pair(key='VERBOSITY', value='WARNING')
    sam_to_fastq.add_option_pair(key='QUIET', value='false')
    sam_to_fastq.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')

    local_exception_str_list = executable_java.run()

    if local_exception_str_list:
        raise Exception('\n'.join(local_exception_str_list))

    platform_unit = str()
    for line_str in sam_header_rg:
        for field in line_str.rstrip().split():
            if field.startswith('PU:'):
                platform_unit = field[3:]

        file_name_prefix = SafeFileName.get_safe_file_name(file_name=platform_unit)
        file_name_1 = os.path.join(path_temporary, file_name_prefix + '_1.fastq')
        file_name_2 = os.path.join(path_temporary, file_name_prefix + '_2.fastq')

        if os.path.exists(file_name_2):
            file_paths.append((file_name_1, file_name_2))
        else:
            file_paths.append((file_name_1, None))

    return file_paths


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = ArgumentParser(
    description='BSF Runner for running the Tuxedo Tophat application.')

argument_parser.add_argument(
    '--pickler_path',
    help='file path to a Python Pickler file',
    required=True)

name_space = argument_parser.parse_args()

# Unpickle the file into a Python dict object.

with open(file=name_space.pickler_path, mode='rb') as input_file:
    pickler_dict = pickle.Unpickler(input_file).load()

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

key = 'java_archive_picard'
if key in pickler_dict and pickler_dict[key]:
    java_archive_picard = pickler_dict[key]
else:
    java_archive_picard = JavaArchive.get_picard()

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, replicate_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

runnable_step_tophat: RunnableStep = pickler_dict['runnable_step']
assert isinstance(runnable_step_tophat, RunnableStep)

# print('Executable before conversion')
# sys.stdout.writelines(runnable_step_tophat.trace(level=1))

# Check the list of file paths in the second and third arguments for FASTQ versus BAM files.

new_file_paths_1 = list()
new_file_paths_2 = list()
old_file_paths_1 = runnable_step_tophat.arguments[1].split(',')
old_file_paths_2 = runnable_step_tophat.arguments[2].split(',')
temporary_files = list()

# Tophat does not like empty compressed FASTQ files.
# Since empty GNU Zip files seem to have 20 bytes,
# only append a file if its size is equal to or greater than 1024 bytes.

minimum_size = 1024

for i in range(0, len(old_file_paths_1)):
    if old_file_paths_1[i].endswith('.bam'):
        # This file needs converting.
        for file_path_1, file_path_2 in run_picard_sam_to_fastq(input_path=old_file_paths_1[i],
                                                                temporary_path=path_temporary):
            if file_path_2 and os.path.exists(file_path_2):
                if (os.path.getsize(file_path_2) >= minimum_size) and (os.path.getsize(file_path_1) >= minimum_size):
                    new_file_paths_1.append(file_path_1)
                    new_file_paths_2.append(file_path_2)
                temporary_files.append(file_path_1)
                temporary_files.append(file_path_2)
            else:
                if os.path.getsize(file_path_1) > minimum_size:
                    new_file_paths_1.append(file_path_1)
                temporary_files.append(file_path_1)
    else:
        file_path_1 = old_file_paths_1[i]
        try:
            # The second list may be shorter (i.e., empty) in case read 2 files are not defined.
            file_path_2 = old_file_paths_2[i]
        except IndexError:
            file_path_2 = None

        if file_path_2 and os.path.exists(file_path_2):
            if (os.path.getsize(file_path_2) >= minimum_size) and (os.path.getsize(file_path_1) >= minimum_size):
                # Only append the pair if both files are larger or equal to 1024 bytes.
                # Empty GNU Zip files teem to have 20 bytes.
                new_file_paths_1.append(file_path_1)
                new_file_paths_2.append(file_path_2)
        else:
            if os.path.getsize(file_path_1) > minimum_size:
                new_file_paths_1.append(file_path_1)

runnable_step_tophat.arguments[1] = ','.join(new_file_paths_1)

if len(new_file_paths_2):
    runnable_step_tophat.arguments[2] = ','.join(new_file_paths_2)
else:
    # If the list of arguments is now empty truncate it to just two (0: genome index, 1: R1 FASTQ files).
    runnable_step_tophat.arguments = runnable_step_tophat.arguments[:2]

# print('Executable after conversion')
# sys.stdout.writelines(runnable_step_tophat.trace(level=1))

exception_str_list = runnable_step_tophat.run()

if exception_str_list:
    raise Exception('\n'.join(exception_str_list))

# Clean up temporary files.

for file_path in temporary_files:
    if os.path.exists(file_path):
        os.remove(file_path)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
