#! /usr/bin/env python
#
# BSF Python script to convert a BED file into a Picard intervals file.
#
#
# Copyright 2014 Michael K. Schuster
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
import os
import string


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = argparse.ArgumentParser(
    description='BSF utility to convert a BED file into a Picard intervals file.')

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
    required=True)

argument_parser.add_argument(
    '--dictionary',
    help='sequence dictionary i.e. SAM header file path',
    required=True)

argument_parser.add_argument(
    '--remove-chr-prefix',
    action='store_true',
    dest='no_prefix',
    help="remove a UCSC-style 'chr' prefix from the sequence name")

name_space = argument_parser.parse_args()

output_file = open(name_space.output_path, 'wb')

# Read the SAM header dictionary and copy to the output file.
input_file = open(name_space.dictionary, 'rb')

for line in input_file:
    output_file.write(line)

input_file.close()

input_file = open(name_space.input_path, 'rb')

for line in input_file:

    bed_fields = line.strip().split()

    if bed_fields[0] == 'browser':
        continue
    if bed_fields[0] == 'track':
        continue

    interval_fields = list()

    if name_space.no_prefix and bed_fields[0][:3] == 'chr':  # Sequence with or without 'chr'.
        interval_fields.append(bed_fields[0][3:])
    else:
        interval_fields.append(bed_fields[0])

    interval_fields.append(str(int(bed_fields[1]) + 1))  # Start: the BED format is half-open, zero-based.
    interval_fields.append(bed_fields[2])  # End
    interval_fields.append('+')  # Strand
    interval_fields.append(bed_fields[3])  # Name

    output_file.write(string.join(words=interval_fields, sep="\t") + "\n")

input_file.close()

output_file.close()
