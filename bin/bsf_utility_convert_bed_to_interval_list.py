#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# BSF Python script to convert a BED file into a Picard interval list file.
#
#
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
import os

# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = argparse.ArgumentParser(
    description='BSF utility to convert a BED file into a Picard interval list file.')

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
    dest='output_path',
    help='Picard-style intervals file path',
    required=False)

argument_parser.add_argument(
    '--dictionary',
    help='sequence dictionary i.e. SAM header file path',
    required=True)

argument_parser.add_argument(
    '--ucsc-to-grc',
    action='store_true',
    dest='ucsc_to_grc',
    help="convert UCSC-style 'chr' prefixed sequence names to GRC reference names")

name_space = argument_parser.parse_args()

if name_space.output_path is None or not name_space.output_path:
    # If the output-path option is missing, construct the output path from the BED input path.
    output_path = name_space.input_path
    assert isinstance(output_path, (str, unicode))
    if output_path.endswith('.bed'):
        output_path = output_path[:-3] + 'interval_list'
    else:
        raise Exception('The input-path option does not specify a *.bed file and the output-path option is missing.')
else:
    output_path = name_space.output_path

output_file = open(output_path, 'wb')

# Read the SAM header dictionary and copy it to the output file.

input_file = open(name_space.dictionary, 'rb')

# To retain the order of sequence regions as defined in the sequence dictionary (i.e. SAM header),
# build a Python list of sequence region names and a Python dict of Python str (sequence region name) key
# and Python list of Python list (interval) objects.

sequence_name_dict = dict()
sequence_name_list = list()

for line in input_file:
    output_file.write(line)
    if line.startswith('@SQ'):
        for sam_field in line.split('\t'):
            if sam_field[:3] == 'SN:':
                sequence_name_dict[sam_field[3:]] = list()
                sequence_name_list.append(sam_field[3:])

input_file.close()

# Read the BED file.

input_file = open(name_space.input_path, 'rb')

for line in input_file:

    bed_fields = line.strip().split()

    if bed_fields[0] == 'browser':
        continue
    if bed_fields[0] == 'track':
        continue

    interval_fields = list()

    if name_space.ucsc_to_grc:
        if bed_fields[0] == 'chrM':  # The 'chrM' needs converting into 'MT'.
            interval_fields.append('MT')
        elif bed_fields[0][:3] == 'chr':  # Sequence with or without 'chr'.
            interval_fields.append(bed_fields[0][3:])
        else:
            interval_fields.append(bed_fields[0])
    else:
        interval_fields.append(bed_fields[0])

    interval_fields.append(str(int(bed_fields[1]) + 1))  # Start: the BED format is half-open, zero-based.
    interval_fields.append(bed_fields[2])  # End
    interval_fields.append('+')  # Strand
    interval_fields.append(bed_fields[3])  # Name

    sequence_name_dict[interval_fields[0]].append(interval_fields)

input_file.close()

# Write interval lines in the order of sequence region names in the sequence dictionary (i.e. SAM header) file.

for sequence_name in sequence_name_list:
    assert isinstance(sequence_name, str)
    interval_list = sequence_name_dict[sequence_name]
    assert isinstance(interval_list, list)
    print('Sequence name: {} lines: {}'.format(sequence_name, len(interval_list)))
    # Sort numerically on the sequence region start field.
    interval_list.sort(key=lambda item: int(item[1]))

    for interval_fields in interval_list:
        output_file.write('\t'.join(interval_fields) + "\n")

output_file.close()
