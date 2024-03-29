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
"""The :py:mod:`bin.bsf_utility_sam_bc_qt_to_fastq` module is a script to
extract the index read sequence (:literal:`BC`) and quality scores (:literal:`QT`) of an
unaligned :emphasis:`BAM` file into a separate :emphasis:`GNU Zip`-compressed :emphasis:`FASTQ` file.
"""

import os
import queue
import sys
import threading
import warnings
from argparse import ArgumentParser
from gzip import GzipFile
from queue import Queue
from typing import Optional

import pysam
from pysam import AlignmentFile

from bsf.standards import SafeFileName


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


def run(
        input_path: str,
        output_path: str,
        npf_reads: Optional[bool] = None) -> int:
    """Run function.

    :param input_path: An input BAM file path.
    :type input_path: str
    :param output_path: An output directory path.
    :type output_path: str
    :param npf_reads: Request including reads not passing vendor quality filters.
    :type npf_reads: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    # Set the environment consistently.

    os.environ['LANG'] = 'C'

    vendor_filter = not npf_reads

    alignment_file = AlignmentFile(filename=input_path, mode='rb', check_sq=False)

    alignment_header_dict = alignment_file.header.to_dict()

    # Open a GzipFile object for each ReadGroup on the basis of @RG PU entries.

    # This procedure only works for query name sorted SAM files.

    hd_dict = alignment_header_dict['HD']

    if 'SO' in hd_dict:
        if hd_dict['SO'] not in ('queryname', 'unsorted'):
            raise Exception("Can only work on 'queryname' or 'unsorted' BAM files.")
    else:
        warnings.warn("Could not find an 'SO' tag in the '@HD' line.", UserWarning)

    gzip_file_dict: dict[str, list[GzipFile]] = dict()

    fifo_queue_dict: dict[str, list[Queue]] = dict()

    rg_dict: dict[str, str]
    for rg_dict in alignment_header_dict['RG']:
        platform_unit = SafeFileName.get_safe_file_name(file_name=rg_dict['PU'])

        if rg_dict['ID'] not in gzip_file_dict:
            gzip_file_dict[rg_dict['ID']] = list()
        else:
            warnings.warn(f'The ReadGroup ID {rg_dict["ID"]!r} is already present in the gzip_file_dict.', UserWarning)

        if rg_dict['ID'] not in fifo_queue_dict:
            fifo_queue_dict[rg_dict['ID']] = list()
        else:
            warnings.warn(f'The ReadGroup ID {rg_dict["ID"]!r} is already present in the fifo_queue_dict.', UserWarning)

        gzip_file_list = gzip_file_dict[rg_dict['ID']]

        fifo_queue_list = fifo_queue_dict[rg_dict['ID']]

        for suffix in ('1', '2', 'i'):
            # Create a GzipFile.
            gzip_file = GzipFile(
                filename=os.path.join(output_path, platform_unit + '_' + suffix + '.fastq.gz'),
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

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='BSF utility to extract index sequences of an unaligned BAM file into a FASTQ file.')

    argument_parser.add_argument(
        '--input-path',
        required=True,
        help='input BAM file path')

    argument_parser.add_argument(
        '--output-path',
        default='.',
        help='output directory path')

    argument_parser.add_argument(
        '--npf-reads',
        action='store_true',
        help='include reads not passing vendor quality filters')

    name_space = argument_parser.parse_args()

    return run(input_path=name_space.input_path, output_path=name_space.output_path, npf_reads=name_space.npf_reads)


if __name__ == '__main__':
    sys.exit(main())
