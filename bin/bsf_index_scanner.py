#! /usr/bin/env python
#
# BSF Python script to list multiplexing indices (barcodes).
#
#
# Copyright 2012 - 2013 Andreas Schoenegger
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

from __future__ import print_function

import argparse
import collections
import datetime
import os
import platform
import subprocess
import sys
import threading

import pysam

import bsf.process

# Set the environment consistently.

os.environ['LANG'] = 'C'

max_thread_joins = 5
thread_join_timeout = 5

# Parse the arguments.

parser = argparse.ArgumentParser(
    description='BSF barcode index scanner.')

parser.add_argument('--debug', required=False, type=int,
                    help='debug level')

parser.add_argument('--input_file', help='file path to a BAM file')
parser.add_argument('--output_file', help='file path to a TSV file')

name_space = parser.parse_args()

# TODO: dual indices not supported
# TODO: GNU Zip FASTQ files could be read directly: http://docs.python.org/2/library/gzip.html

file_type = os.path.splitext(name_space.input_file)[-1]  # supported formats: .fastq/.sam/.bam

log_after_x_processed_reads = 100000

# Global variable

barcode_dict = dict()

alignment_file = pysam.AlignmentFile(name_space.input_file, 'rb', check_sq=False)
alignment_iterator = alignment_file.fetch(until_eof=True)
""" @type alignment_iterator: pysam.IteratorRowAll """
alignment_counter = 0
for aligned_segment in alignment_iterator:
    """ @type aligned_segment: pysam.AlignedSegment """
    barcode = aligned_segment.get_tag(tag='BC')
    if barcode not in barcode_dict:
        barcode_dict[barcode] = 1
    else:
        barcode_dict[barcode] += 1

    alignment_counter += 1
    if not alignment_counter % log_after_x_processed_reads:
        print('Processed {} reads'.format(alignment_counter))

alignment_file.close()

alignment_file = open(name_space.output_file, 'w')
alignment_file.write('Barcode\tCount\n')

for barcode in sorted(barcode_dict):
    alignment_file.write(barcode + '\t' + str(barcode_dict[barcode]) + '\n')

alignment_file.close()

exit()


def parse_sam_format(sam_file_handle):
    """Parses SAM format columns.
    This function identifies the column with the barcode based on its suffix BC:Z:
    and increments the respective count in a Python C{dict}.

    @param sam_file_handle: File handle
    @type sam_file_handle: FileIO
    @return: Python C{dict} of barcode key and count value data
    @rtype: dict
    """

    for line in sam_file_handle:
        assert isinstance(line, str)
        if line.startswith('@'):
            continue
        columns = line.rstrip().split('\t')
        # Find the column with the BC:Z: tag.
        # The first 11 columns are fixed.
        for i in range(10, len(columns)):
            if columns[i].startswith('BC:Z:'):
                barcode_str = columns[i][5:]
                if barcode_str in barcode_dict:
                    barcode_dict[barcode_str] += 1
                else:
                    barcode_dict[barcode_str] = 1
                break


def parse_sam_file(input_filename):
    """
    Parse a SAM file, processing barcodes in the 12th column.

    @param input_filename: File path
    @type input_filename: str | unicode
    @return: Python C{dict} of barcode key and count value data
    @rtype: dict
    """

    # TODO: breaks if barcode not in 12th column or tag has different name (expecting BC:Z:....)

    print('Processing SAM file (expecting barcode in the 12th column of the SAM file as a tag, '
          'like BC:Z:GGAACC (Picard - IlluminaToSam creates these unmapped BAM files)')

    input_file = open(input_filename)
    barcodes = dict()

    i = 0
    for line in input_file:
        # ignore SAM header (starts with @)
        if not line.startswith('@'):
            # split each line at tabs, extract 12 column, split barcode tag
            bc = line.split('\t')[11].split(':')[2]

            # store barcode in a dict
            if bc not in barcodes:
                barcodes[bc] = 0

            barcodes[bc] += 1

        # log message
        i += 1
        if i % log_after_x_processed_reads == 0:
            print('Read ' + str(i) + ' lines.')

    input_file.close()

    return barcodes


# converts a BAM file to SAM, runs parse_sam_file on the temporary SAM file
def parse_bam_file(input_filename):
    """
    Parse a BAM file.
    Internally converts a BAM file to a SAM file and subsequently runs parse_sam_file over it.

    @param input_filename: File path
    @type input_filename: str | unicode
    @return: Python C{dict} of barcode key and count value data
    @rtype: dict
    """

    if platform.system() == 'Windows':
        raise Exception('BAM files not supported on Windows')

    print('Processing BAM file, running samtools view <bamfile> > <tmp_samfile> (samtools required in PATH!)')

    tmp_samfile = os.path.splitext(input_filename)[0] + '.tmp.sam'

    exit_code = os.system('samtools view -h ' + input_filename + ' > ' + tmp_samfile)

    if not exit_code == 0:
        raise Exception('Error executing samtools. Check if samtools is in your PATH (e.g. which samtools)')

    barcodes = parse_sam_file(tmp_samfile)

    os.remove(tmp_samfile)

    return barcodes


# parses a FASTQ file, processes barcode in header line of each FASTQ entry
# TODO: breaks if barcode has /1 or /2 at the end (I think FASTQ files created by Picard from PE BAM files have these)
def parse_fastq_file(input_filename):
    print('Processing FASTQ file (expecting barcode in the header)')

    input_file = open(input_filename)
    barcodes = dict()

    i = 0
    for line in input_file:
        if i % 4 == 0:
            bc = line.rstrip().split(':').pop()
            if bc not in barcodes:
                barcodes[bc] = 0

            barcodes[bc] += 1

        i += 1
        if i % log_after_x_processed_reads == 0:
            print('Read ' + str(i) + ' lines.')

    input_file.close()
    return barcodes


if file_type == '.fastq':
    barcodes_dict = parse_fastq_file(name_space.input_file)
elif file_type == '.sam':
    file_handle = open(name_space.input_file, 'r')
    parse_sam_format(sam_file_handle=file_handle)
    file_handle.close()
elif file_type == '.bam':
    executable = bsf.process.Executable(
        name='samtools_view',
        program='samtools',
        sub_command=bsf.process.Command(program='view'))
    sub_command = executable.sub_command
    sub_command.add_option_short(key='f', value='64')  # Select only the first read of a pair.
    sub_command.add_option_short(key='F', value='512')  # Suppress all reads that fail vendor quality filtering.
    sub_command.add_switch_short(key='h')
    sub_command.arguments.append(name_space.input_file)

    on_posix = 'posix' in sys.builtin_module_names

    child_process = subprocess.Popen(args=executable.command_list(),
                                     bufsize=-1,
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=False,
                                     close_fds=on_posix)

    # Two threads, thread_out and thread_err reading STDOUT and STDERR, respectively,
    # should make sure that buffers are not filling up.

    thread_lock = threading.Lock()

    thread_out = threading.Thread(
        target=parse_sam_format,
        kwargs={'sam_file_handle': child_process.stdout, })
    thread_out.daemon = True  # Thread dies with the program.
    thread_out.start()

    thread_err = threading.Thread(
        target=bsf.process.Executable.process_stderr,
        kwargs={
            'stderr_handle': child_process.stderr,
            'thread_lock': thread_lock,
            'stderr_path': executable.stderr_path,
            'debug': 0,
        })
    thread_err.daemon = True  # Thread dies with the program.
    thread_err.start()

    # Wait for the child process to finish.

    child_return_code = child_process.wait()

    thread_join_counter = 0

    while thread_out.is_alive() and thread_join_counter < max_thread_joins:
        thread_lock.acquire(True)
        # print('[{}] Waiting for STDOUT processor to finish.'.format(
        #     datetime.datetime.now().isoformat()))
        thread_lock.release()

        thread_out.join(timeout=thread_join_timeout)
        thread_join_counter += 1

    thread_join_counter = 0

    while thread_err.is_alive() and thread_join_counter < max_thread_joins:
        thread_lock.acquire(True)
        # print('[{}] Waiting for STDERR processor to finish.'.format(
        #     datetime.datetime.now().isoformat()))
        thread_lock.release()

        thread_err.join(timeout=thread_join_timeout)
        thread_join_counter += 1

    if child_return_code > 0:
        if name_space.debug > 0:
            print('[{}] Child process {!r} failed with exit code {}'.format(
                datetime.datetime.now().isoformat(), executable.name, +child_return_code))
    elif child_return_code < 0:
        if name_space.debug > 0:
            print('[{}] Child process {!r} received signal {}.'.format(
                datetime.datetime.now().isoformat(), executable.name, -child_return_code))
    else:
        if name_space.debug > 0:
            print('[{}] Child process {!r} completed successfully {}.'.format(
                datetime.datetime.now().isoformat(), executable.name, +child_return_code))
else:
    raise Exception('Unsupported file format.')

# write report
# outfile = open('unmatched_barcode_report.csv', 'w')
# Write dict sorted by highest number of occurrence
for barcode in collections.OrderedDict(sorted(barcode_dict.items(), reverse=True, key=lambda t: t[1])):
    message = barcode + ';' + str(barcode_dict[barcode])

    # outfile.write(message + '\n')
    print(message)

# Check again for Illumina Barcodes:
# outfile.close()

print('--------')

illumina_adapters = [
    'ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT',
    'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT', 'GGCTAC', 'CTTGTA',
    'AGTCAA', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTCCGC', 'GTGAAA',
    'GTGGCC', 'GTTTCG', 'CGTACG', 'GAGTGG', 'ACTGAT', 'ATTCCT'
]

# outfile = open('unmatched_barcode_report.illumina.csv', 'w')

for adapter in illumina_adapters:
    assert isinstance(adapter, str)
    count = 0
    if adapter in barcode_dict:
        count = barcode_dict[adapter]

    message = adapter + ';' + str(count)

    # outfile.write(message + '\n')
    print(message)

# outfile.close()
