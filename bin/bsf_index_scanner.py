#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
#
# BSF Python script to list multiplexing indices (barcodes).
#
#
# Copyright 2012 - 2013 Andreas Schoenegger
# Copyright 2013 - 2019 Michael K. Schuster
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
#

from __future__ import print_function

import argparse
import collections
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

argument_parser = argparse.ArgumentParser(
    description='BSF barcode index scanner.')

argument_parser.add_argument(
    '--debug',
    default=0,
    required=False,
    type=int,
    help='debug level')

argument_parser.add_argument('--input_file', help='file path to a BAM file')
argument_parser.add_argument('--output_file', help='file path to a TSV file')

name_space = argument_parser.parse_args()

# TODO: dual indices not supported
# TODO: GNU Zip FASTQ files could be read directly: http://docs.python.org/2/library/gzip.html

file_type = os.path.splitext(name_space.input_file)[-1]  # supported formats: .fastq/.sam/.bam

log_after_x_processed_reads = 100000

# Global variable

barcode_dict = dict()
""" @type barcode_dict: dict[str , int] """

with pysam.AlignmentFile(name_space.input_file, 'rb', check_sq=False) as alignment_file:
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

with open(file=name_space.output_file, mode='wt') as output_file:
    output_file.write('Barcode\tCount\n')

    for barcode in sorted(barcode_dict):
        output_file.write(barcode + '\t' + str(barcode_dict[barcode]) + '\n')

exit()


def parse_sam_format(file_path):
    """Parses SAM format columns.
    This function identifies the column with the barcode based on its suffix BC:Z:
    and increments the respective count in a Python C{dict}.

    @param file_path: File path
    @type file_path: str | unicode
    @return:
    @rtype:
    """
    with open(file=file_path, mode='rt') as input_file:
        for line_str in input_file:
            if line_str.startswith('@'):
                continue
            columns = line_str.rstrip().split('\t')
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

    barcodes = dict()
    with open(file=input_filename, mode='rt') as input_file:
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

    return barcodes


# converts a BAM file to SAM, runs parse_sam_file on the temporary SAM file
def parse_bam_file(input_filename):
    """Parse a BAM file.

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
def parse_fastq_file(file_path):
    print('Processing FASTQ file (expecting barcode in the header)')

    with open(file=file_path, mode='rt') as input_file:
        i = 0
        for line_str in input_file:
            if i % 4 == 0:
                bc = line_str.rstrip().split(':').pop()
                if bc not in barcode_dict:
                    barcode_dict[bc] = 0

                barcode_dict[bc] += 1

            i += 1
            if i % log_after_x_processed_reads == 0:
                print('Read ' + str(i) + ' lines.')


if file_type == '.fastq':
    parse_fastq_file(file_path=name_space.input_file)
elif file_type == '.sam':
    parse_sam_format(file_path=name_space.input_file)
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

    child_process = subprocess.Popen(
        args=executable.command_list(),
        bufsize=-1,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=False,
        close_fds=on_posix,
        text=True)

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
            # 'stderr_path': executable.stderr_path,
            'debug': 0,
        })
    thread_err.daemon = True  # Thread dies with the program.
    thread_err.start()

    # Wait for the child process to finish.

    child_return_code = child_process.wait()

    thread_join_counter = 0

    while thread_out.is_alive() and thread_join_counter < max_thread_joins:
        # thread_lock.acquire(True)
        # print('[{}] Waiting for STDOUT processor to finish.'.format(
        #     datetime.datetime.now().isoformat()))
        # thread_lock.release()

        thread_out.join(timeout=thread_join_timeout)
        thread_join_counter += 1

    thread_join_counter = 0

    while thread_err.is_alive() and thread_join_counter < max_thread_joins:
        # thread_lock.acquire(True)
        # print('[{}] Waiting for STDERR processor to finish.'.format(
        #     datetime.datetime.now().isoformat()))
        # thread_lock.release()

        thread_err.join(timeout=thread_join_timeout)
        thread_join_counter += 1

    if name_space.debug > 0:
        executable.evaluate_return_code(return_code=child_return_code)
else:
    raise Exception('Unsupported file format.')

# write report
# with open(file='unmatched_barcode_report.csv', mode='wt') as output_file:
# Write dict sorted by highest number of occurrence
for barcode in collections.OrderedDict(sorted(barcode_dict.items(), reverse=True, key=lambda item: item[1])):
    message = barcode + ';' + str(barcode_dict[barcode])

    # output_file.write(message + '\n')
    print(message)

# Check again for Illumina Barcodes:

print('--------')

illumina_adapters = [
    'ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT',
    'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT', 'GGCTAC', 'CTTGTA',
    'AGTCAA', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTCCGC', 'GTGAAA',
    'GTGGCC', 'GTTTCG', 'CGTACG', 'GAGTGG', 'ACTGAT', 'ATTCCT',
]
""" @type illumina_adapters: list[str] """

# with open(file='unmatched_barcode_report.illumina.csv', mode='wt') as output_file:

for adapter in illumina_adapters:
    count = 0
    if adapter in barcode_dict:
        count = barcode_dict[adapter]

    message = adapter + ';' + str(count)

    # output_file.write(message + '\n')
    print(message)
