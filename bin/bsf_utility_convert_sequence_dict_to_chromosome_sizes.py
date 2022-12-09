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
"""The :py:mod:`bin.bsf_utility_convert_sequence_dict_to_chromosome_sizes` module is a script to
convert a Picard sequence dictionary based on a SAM header file into a UCSC chromosome sizes file.

The SAM file header specifies in sequence (:literal:`@SQ`) lines
sequence names (:literal:`SN`) and length (:literal:`LN`) fields,
while the UCSC chromosome sizes file just contains sequence name
and length separated by tabs.
"""

import os
import sys
from argparse import ArgumentParser


def run(
        input_path: str,
        output_path: str) -> int:
    """Run function.

    :param input_path: A Picard sequence dictionary file path.
    :type input_path: str
    :param output_path: A UCSC chromosome sizes file path.
    :type output_path: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    # Set the environment consistently.

    os.environ['LANG'] = 'C'

    # Parse the arguments.

    with open(file=output_path, mode='wt') as output_text_io:
        with open(file=input_path, mode='rt') as input_text_io:
            for line_str in input_text_io:
                if not line_str.startswith('@SQ'):
                    continue

                column_list = line_str.rstrip().split("\t")
                sequence_name = str()
                sequence_length = str()

                # Find the column with the SN: tag.
                for column_str in column_list:
                    if column_str.startswith('SN:'):
                        sequence_name = column_str[3:]
                        break

                # Find the column with the LN: tag.
                for column_str in column_list:
                    if column_str.startswith('LN:'):
                        sequence_length = column_str[3:]
                        break

                print(sequence_name, sequence_length, sep='/t', file=output_text_io)

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='BSF utility to convert a Picard sequence dictionary (SAM header) '
                    'into a UCSC chromosome sizes file.')

    argument_parser.add_argument(
        '--input-path',
        required=True,
        help='Picard sequence dictionary file path')

    argument_parser.add_argument(
        '--output-path',
        required=True,
        help='UCSC chromosome sizes file path')

    name_space = argument_parser.parse_args()

    return run(input_path=name_space.input_path, output_path=name_space.output_path)


if __name__ == '__main__':
    sys.exit(main())
