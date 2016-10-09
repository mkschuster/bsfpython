#! /usr/bin/env python
#
# BSF Python script to run Tophat.
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
import re
import shutil

from bsf.process import Command, Executable
from bsf.standards import Default


def run_picard_sam_to_fastq(input_path, temporary_path):
    """Convert a BAM file into a pair of FASTQ files for each read group (@RG).

    @param input_path: File path to the BAM file
    @type input_path: str | unicode
    @param temporary_path: File path to a temporary directory
    @type temporary_path: str | unicode
    @return: Python list of tuple objects of file paths of read1 and read2
    @rtype: list
    """

    file_paths = list()

    # Propagate SAM header lines @PG and @RG into the final BAM file.
    # Create a temporary SAM file of the same name to store the SAM header.

    path_temporary_sam = os.path.basename(input_path)
    path_temporary_sam = path_temporary_sam.replace('.bam', '.sam')

    samtools = Executable(
        name='samtools_view',
        program='samtools',
        sub_command=Command(program='view'),
        stdout_path=path_temporary_sam)

    samtools_view = samtools.sub_command
    samtools_view.add_switch_short(key='H')
    samtools_view.arguments.append(input_path)

    child_return_code_local = samtools.run()

    if child_return_code_local:
        raise Exception('Could not complete the samtools view step on the BAM file for the replicate.')

    # SAM header lines that need propagating around FASTQ files. Sigh!
    sam_header_pg = list()
    sam_header_rg = list()

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

    picard_process = java_process.sub_command
    picard_process.add_option_short(key='jar', value=os.path.join(classpath_picard, 'picard.jar'))
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

    child_return_code_local = java_process.run()

    if child_return_code_local:
        raise Exception('Could not complete the Picard SamToFastq step.')

    platform_unit = str()
    for line in sam_header_rg:
        assert isinstance(line, basestring)
        for field in line.rstrip().split():
            if field.startswith('PU:'):
                platform_unit = field[3:]

        file_name_prefix = re.sub(pattern='[^0-9A-Za-z_-]', repl='_', string=platform_unit)
        file_name_1 = os.path.join(path_temporary, file_name_prefix + '_1.fastq')
        file_name_2 = os.path.join(path_temporary, file_name_prefix + '_2.fastq')

        if os.path.exists(path=file_name_2):
            file_paths.append((file_name_1, file_name_2))
        else:
            file_paths.append((file_name_1, None))

    return file_paths


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Get global defaults.

default = Default.get_global_default()

# Parse the arguments.

parser = argparse.ArgumentParser(
    description='BSF Runner for running the Tuxedo Tophat application.')

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
    classpath_picard = default.classpath_picard

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, replicate_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

run_tophat = pickler_dict['tophat_executable']
assert isinstance(run_tophat, Executable)

# Check the list of file paths in the second and third arguments for FASTQ versus BAM files.

new_file_paths_1 = list()
new_file_paths_2 = list()
old_file_paths_1 = run_tophat.arguments[1].split(',')
old_file_paths_2 = run_tophat.arguments[2].split(',')
temporary_files = list()

for i in range(0, len(old_file_paths_1)):
    if old_file_paths_1[i][-4:] == '.bam':
        # This file needs converting.
        for file_path_1, file_path_2 in run_picard_sam_to_fastq(input_path=old_file_paths_1[i],
                                                                temporary_path=path_temporary):
            if file_path_2 and os.path.exists(path=file_path_2):
                new_file_paths_1.append(file_path_1)
                new_file_paths_2.append(file_path_2)
                temporary_files.append(file_path_1)
                temporary_files.append(file_path_2)
            else:
                new_file_paths_1.append(file_path_1)
                temporary_files.append(file_path_1)
    else:
        new_file_paths_1.append(old_file_paths_1[i])
        if old_file_paths_2[i] and os.path.exists(path=old_file_paths_2[i]):
            new_file_paths_2.append(old_file_paths_2[i])

run_tophat.arguments[1] = ','.join(new_file_paths_1)

if len(new_file_paths_2):
    run_tophat.arguments[2] = ','.join(new_file_paths_2)
else:
    # If the list of arguments is now empty truncate it to just two.
    run_tophat.arguments = run_tophat.arguments[:2]

child_return_code = run_tophat.run()

if child_return_code:
    raise Exception('Could not complete the Tophat step.')

# Clean up temporary files.

for file_path in temporary_files:
    if os.path.exists(path=file_path):
        os.remove(file_path)

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
