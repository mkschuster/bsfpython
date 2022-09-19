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
#  BSF Python script to extract the index read sequence (BC) and quality scores (QT) of an unaligned BAM file
#  into a separate GNU Zip-compressed FASTQ file.
#
import os
import queue
import re
import threading
import warnings
from argparse import ArgumentParser
from gzip import GzipFile
from queue import Queue
from typing import Dict, List

import pysam
from pysam import AlignmentFile


def write_gzip_file(task_gzip_file, task_fifo_queue):
    """Write items from a or to a :py:class:`gzip.GzipFile` object.

    :param task_gzip_file: A :py:class:`gzip.GzipFile` object.
    :type task_gzip_file: gzip.GzipFile
    :param task_fifo_queue: A :py:class:`queue.Queue` object.
    :type task_fifo_queue: queue.Queue
    """
    while True:
        line_str: str = task_fifo_queue.get()
        task_gzip_file.write(line_str)
        task_fifo_queue.task_done()


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = ArgumentParser(
    description='BSF utility to extract index sequences of an unaligned BAM file into a FASTQ file.')

argument_parser.add_argument(
    '--input-path',
    dest='input_path',
    help='BAM file path',
    required=True)

argument_parser.add_argument(
    '--output-path',
    default='.',
    dest='output_path',
    help='Output directory path',
    required=False)

argument_parser.add_argument(
    '--npf-reads',
    action='store_true',
    dest='npf_reads',
    help='include reads not passing vendor quality filters')

name_space = argument_parser.parse_args()

vendor_filter = not name_space.npf_reads

alignment_file = AlignmentFile(filename=name_space.input_path, mode='rb', check_sq=False)

alignment_header_dict = alignment_file.header.to_dict()

# Open a GzipFile object for each ReadGroup on the basis of @RG PU entries.

# This procedure only works for query name sorted SAM files.

hd_dict = alignment_header_dict['HD']

if 'SO' in hd_dict:
    if hd_dict['SO'] not in ('queryname', 'unsorted'):
        raise Exception("Can only work on 'queryname' or 'unsorted' BAM files")
else:
    warnings.warn("Could not find a 'SO' tag in the '@HD' line.")

gzip_file_dict: Dict[str, List[GzipFile]] = dict()

fifo_queue_dict: Dict[str, List[Queue]] = dict()

rg_dict: Dict[str, str]
for rg_dict in alignment_header_dict['RG']:
    # The makeFileNameSafe() method of htsjdk.samtools.util.IOUtil uses the following pattern:
    # [\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]
    platform_unit = re.sub(
        pattern='[\\s!"#$%&\'()*/:;<=>?@\\[\\]\\\\^`{|}~]',
        repl='_',
        string=rg_dict['PU'])

    if rg_dict['ID'] not in gzip_file_dict:
        gzip_file_dict[rg_dict['ID']] = list()
    else:
        warnings.warn('ReadGroup ID already present in gzip_file_dict: ' + rg_dict['ID'])

    if rg_dict['ID'] not in fifo_queue_dict:
        fifo_queue_dict[rg_dict['ID']] = list()
    else:
        warnings.warn('ReadGroup ID already present in fifo_queue_dict: ' + rg_dict['ID'])

    gzip_file_list = gzip_file_dict[rg_dict['ID']]

    fifo_queue_list = fifo_queue_dict[rg_dict['ID']]

    for suffix in ('1', '2', 'i'):
        # Create a GzipFile.
        gzip_file = GzipFile(
            filename=os.path.join(name_space.output_path, platform_unit + '_' + suffix + '.fastq.gz'),
            mode='wb',
            compresslevel=9)
        gzip_file_list.append(gzip_file)

        # Create a FIFO Queue with 100 items maximum.
        fifo_queue = queue.Queue(maxsize=100)
        fifo_queue_list.append(fifo_queue)

        # Create a daemon Thread and start it.
        read_thread = threading.Thread(
            target=write_gzip_file,
            kwargs={
                'task_gzip_file': gzip_file,
                'task_fifo_queue': fifo_queue,
            })
        read_thread.daemon = True
        read_thread.start()

for aligned_segment in alignment_file.fetch(until_eof=True):
    if vendor_filter and aligned_segment.is_qcfail:
        continue

    # Assign the AlignedSegment to its ReadGroup-specific GzipFile.
    if aligned_segment.is_read1:
        fifo_queue = fifo_queue_dict[aligned_segment.get_tag('RG')][0]
    else:
        fifo_queue = fifo_queue_dict[aligned_segment.get_tag('RG')][1]

    fifo_queue.put(
        '@' + aligned_segment.query_name + '\n' +
        aligned_segment.query_sequence + '\n' +
        '+\n' +
        pysam.libcutils.array_to_qualitystring(aligned_segment.query_qualities) + '\n'
    )

    if aligned_segment.has_tag('BC'):
        # Assign the AlignedSegment to its ReadGroup-specific GzipFile.
        fifo_queue = fifo_queue_dict[aligned_segment.get_tag('RG')][2]

        fifo_queue.put(
            '@' + aligned_segment.query_name + '\n' +
            aligned_segment.get_tag('BC') + '\n' +
            '+\n' +
            aligned_segment.get_tag('QT') + '\n'
        )

for read_group_id in fifo_queue_dict:
    for fifo_queue in fifo_queue_dict[read_group_id]:
        fifo_queue.join()

for read_group_id in gzip_file_dict:
    for gzip_file in gzip_file_dict[read_group_id]:
        gzip_file.close()

alignment_file.close()
