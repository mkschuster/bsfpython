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
"""The :py:mod:`bin.bsf_utility_convert_bed_to_interval_list` module is a script to
convert a BED file into a Picard interval list file.
"""

import os
import sys
from argparse import ArgumentParser
from typing import Optional


def run(
        dictionary_path: str,
        input_path: str,
        output_path: Optional[str] = None,
        drop_names: Optional[bool] = None,
        ucsc_to_grc: Optional[bool] = None) -> int:
    """Run function.

    :param dictionary_path: A sequence dictionary (i.e., SAM header) file path.
    :type dictionary_path: str
    :param input_path: An input BED file path.
    :type input_path: str
    :param output_path: An output Picard-style intervals file path.
    :type output_path: str | None
    :param drop_names: Request dropping probe names.
    :type drop_names: bool | None
    :param ucsc_to_grc: Request converting UCSC-style :emphasis:`chr`-prefixed sequence names to GRC reference names.
    :type ucsc_to_grc: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    # Set the environment consistently.

    os.environ['LANG'] = 'C'

    if output_path:
        output_path = output_path
    else:
        # If the output-path option is missing, construct the output path from the BED input path.
        output_path = input_path
        if output_path.endswith('.bed'):
            output_path = output_path[:-3] + 'interval_list'
        else:
            raise Exception('The --input-path option does not specify a BED file '
                            'and the --output-path option is missing.')

    # To retain the order of sequence regions as defined in the sequence dictionary (i.e. SAM header),
    # build a Python list of sequence region names and a Python dict of Python str (sequence region name) key
    # and Python list of Python list (interval) objects.

    sequence_name_dict: dict[str, list[list[str]]] = dict()
    sequence_name_list: list[str] = list()

    with open(file=output_path, mode='wt') as output_text_io:
        # Read the SAM header dictionary and copy it to the output file.
        with open(file=dictionary_path, mode='rt') as input_text_io:
            for line_str in input_text_io:
                output_text_io.write(line_str)
                if line_str.startswith('@SQ'):
                    for sam_field in line_str.split('\t'):
                        if sam_field.startswith('SN:'):
                            sequence_name_dict[sam_field[3:]] = list()
                            sequence_name_list.append(sam_field[3:])

        # Read the BED file.
        with open(file=input_path, mode='rt') as input_text_io:
            for line_str in input_text_io:
                bed_fields = line_str.strip().split()

                if bed_fields[0] == 'browser':
                    continue
                if bed_fields[0] == 'track':
                    continue

                interval_fields = list()

                if ucsc_to_grc:
                    if bed_fields[0] == 'chrM':  # The 'chrM' needs converting into 'MT'.
                        interval_fields.append('MT')
                    elif bed_fields[0] == 'chrUn_gl000228':
                        interval_fields.append('GL000228.1')
                    elif bed_fields[0].startswith('chr'):  # Sequence with or without 'chr'.
                        interval_fields.append(bed_fields[0][3:])
                    else:
                        interval_fields.append(bed_fields[0])
                else:
                    interval_fields.append(bed_fields[0])

                interval_fields.append(str(int(bed_fields[1]) + 1))  # Start: the BED format is half-open, zero-based.
                interval_fields.append(bed_fields[2])  # End
                interval_fields.append('+')  # Strand
                if len(bed_fields) >= 4 and not drop_names:
                    interval_fields.append(bed_fields[3])  # Name

                sequence_name_dict[interval_fields[0]].append(interval_fields)

        # Write interval lines in the order of sequence region names in the sequence dictionary (i.e. SAM header) file.
        for sequence_name in sequence_name_list:
            interval_list = sequence_name_dict[sequence_name]
            print('Sequence name: {} lines: {}'.format(sequence_name, len(interval_list)))
            # Sort numerically on the sequence region start field.
            interval_list.sort(key=lambda item: int(item[1]))

            for interval_fields in interval_list:
                output_text_io.write('\t'.join(interval_fields) + "\n")

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='BSF utility to convert a BED file into a Picard interval list file.')

    argument_parser.add_argument(
        '--dictionary-path',
        required=True,
        help='sequence dictionary (i.e., SAM header) file path')

    argument_parser.add_argument(
        '--input-path',
        required=True,
        help='BED file path')

    argument_parser.add_argument(
        '--output-path',
        help='Picard-style intervals file path')

    argument_parser.add_argument(
        '--drop-names',
        action='store_true',
        help='drop probe names')

    argument_parser.add_argument(
        '--ucsc-to-grc',
        action='store_true',
        help="convert UCSC-style 'chr'-prefixed sequence names to GRC reference names")

    name_space = argument_parser.parse_args()

    return run(
        dictionary_path=name_space.dictionary_path,
        input_path=name_space.input_path,
        output_path=name_space.output_path,
        drop_names=name_space.drop_names,
        ucsc_to_grc=name_space.ucsc_to_grc)


if __name__ == '__main__':
    sys.exit(main())
