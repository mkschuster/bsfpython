#! /usr/bin/env python
#
# BSF Python script to extract the index read sequence (BC) and quality scores (QT) of an unaligned BAM file
# into a separate GNU Zip-compressed FASTQ file.
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

from __future__ import print_function

import argparse
import gzip
import os
import re

import pysam

# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = argparse.ArgumentParser(
    description='BSF utility to extract index sequences of an unaligned BAM file into a FASTQ file.')

argument_parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--input-path',
    dest='input_path',
    help='BED file path',
    required=True)

argument_parser.add_argument(
    '--output-path',
    default='.',
    dest='output_path',
    help='Picard-style intervals file path',
    required=False)

name_space = argument_parser.parse_args()

alignment_file = pysam.AlignmentFile(name_space.input_path, 'rb', check_sq=False)

# Open a GzipFile object for each ReadGroup on the basis of @RG PU entries.

gzip_file_dict = dict()
""" @type file_dict: dict[str, gzip.GzipFile] """

for read_group_dict in alignment_file.header['RG']:
    """ @type read_group_dict: dict[str, str] """
    # The makeFileNameSafe() method of htsjdk.samtools.util.IOUtil uses the following pattern:
    # [\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]
    platform_unit = re.sub(
        pattern='[\\s!"#$%&\'()*/:;<=>?@\\[\\]\\\\^`{|}~]',
        repl='_',
        string=read_group_dict['PU'])

    gzip_file_dict[read_group_dict['ID']] = gzip.GzipFile(
        filename=os.path.join(name_space.output_path, platform_unit + '_i.fastq.gz'),
        mode='wb',
        compresslevel=9)

for aligned_segment in alignment_file:
    """ @type aligned_segment: pysam.AlignedSegment """
    if aligned_segment.has_tag('BC'):
        # Assign the AlignedSegment to its ReadGroup-specific GzipFile.
        gzip_file = gzip_file_dict[aligned_segment.get_tag('RG')]

        gzip_file.write('@' + aligned_segment.query_name + '\n')
        gzip_file.write(aligned_segment.get_tag('BC') + '\n')
        gzip_file.write('+\n')
        gzip_file.write(aligned_segment.get_tag('QT') + '\n')

for gzip_file in gzip_file_dict.itervalues():
    gzip_file.close()

alignment_file.close()
