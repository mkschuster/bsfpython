#! /usr/bin/env python
#
# BSF Python script to convert a Picard sequence dictionary based on a
# SAM header file into a UCSC chromosome sizes file.
#
# The SAM file header specifies in sequence (@SQ) lines
# sequence names (SN:) and length (LN:) fields,
# while the UCSC chromosome sizes file just contains sequence name
# and length separated by tabs.
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
import os

# Set the environment consistently.

os.environ['LANG'] = 'C'

# Parse the arguments.

argument_parser = argparse.ArgumentParser(
    description='BSF utility to convert a Picard sequence dictionary (SAM header) '
                'into a UCSC chromosome sizes file.')

argument_parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--input_path',
    help='file path to a Picard sequence dictionary file.',
    required=True)

argument_parser.add_argument(
    '--output_path',
    help='file path to a UCSC chromosome sizes file.',
    required=True)

name_space = argument_parser.parse_args()

input_file = open(name_space.input_path, 'r')
output_file = open(name_space.output_path, 'w')

for line in input_file:
    assert isinstance(line, str)
    if not line.startswith('@SQ'):
        continue
    columns = line.rstrip().split("\t")
    sequence_name = str()
    sequence_length = str()

    # Find the column with the SN: tag.
    for column in columns:
        if column.startswith('SN:'):
            sequence_name = column[3:]
            break

    # Find the column with the LN: tag.
    for column in columns:
        if column.startswith('LN:'):
            sequence_length = column[3:]
            break

    output_file.write("\t".join((sequence_name, sequence_length)) + "\n")

input_file.close()
output_file.close()
