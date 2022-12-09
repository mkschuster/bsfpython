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
"""The :py:mod:`bin.bsf_utility_convert_vcf_to_23andMe` module is a script to
convert a VCF file into the simple :emphasis:`23andMe` format.

The :emphasis:`23andMe` format consists of an :emphasis:`NCBI dbSNP` reference SNP identifier (rsid),
the sequence name, the position and the genotype for each allele (diploid).
"""

import sys

from argparse import ArgumentParser


def run(
        input_path: str,
        output_path: str) -> int:
    """Run function.

    :param input_path: An input VCF file path.
    :type input_path: str
    :param output_path: An output :emphasis:`23andMe` file path.
    :type output_path: str
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    with open(file=output_path, mode='wt') as output_text_io:
        with open(file=input_path, mode='rt') as input_text_io:
            for line_str in input_text_io:
                # Ignore comment lines.
                if line_str.startswith('#'):
                    continue

                # Split VCF lines on tabs into VCF fields.
                vcf_field_list = line_str.rstrip().split('\t')

                # The genotype (GT) is in field 10, the REF in field 4 and ALT in field 5.

                genotype_index = 0

                info_field_list = vcf_field_list[8].split(':')

                for i in range(0, len(info_field_list)):
                    if info_field_list[i] == 'GT':
                        genotype_index = i
                        break

                sample_field_list = vcf_field_list[9].split(':')

                if sample_field_list[genotype_index] == '0/0':
                    alleles = vcf_field_list[3] + vcf_field_list[3]
                elif sample_field_list[genotype_index] == '0/1':
                    alleles = vcf_field_list[3] + vcf_field_list[4]
                elif sample_field_list[genotype_index] == '1/1':
                    alleles = vcf_field_list[4] + vcf_field_list[4]
                elif sample_field_list[genotype_index] == './.':
                    continue
                else:
                    alleles = '..'
                    print(f'Unexpected genotype {vcf_field_list[9]!r} in line {line_str!r}')

                # The identifier (ID) is in field 2,
                # the chromosome (CHROM) in field 0 and the position (POS) in field 1.

                print(vcf_field_list[2], vcf_field_list[0], vcf_field_list[1], alleles, sep='\t', file=output_text_io)

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='BSF utility to convert from VCF to 23andMe format.')

    argument_parser.add_argument(
        '--input-path',
        required=True,
        help='VCF file path')

    argument_parser.add_argument(
        '--output-path',
        required=True,
        help='23andMe file path')

    name_space = argument_parser.parse_args()

    return run(input_path=name_space.input_path, output_path=name_space.output_path)


if __name__ == '__main__':
    sys.exit(main())
