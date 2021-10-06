#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  Copyright 2013 - 2021 Michael K. Schuster
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
#  BSF Python script to convert a VCF file into the simple 23andMe format.
#
#  The 23andMe format consists of a dbSNP reference SNP identifier (rsid),
#  the sequence name, the position and the genotype for each allele (diploid).
#
from argparse import ArgumentParser

argument_parser = ArgumentParser(
    description='BSF utility to convert from VCF to 23andMe format.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--input',
    help='VCF file path',
    required=True)

argument_parser.add_argument(
    '--output',
    help='23andMe file path',
    required=True)

name_space = argument_parser.parse_args()

with open(file=name_space.output, mode='wt') as output_file:
    with open(file=name_space.input, mode='rt') as input_file:
        for line_str in input_file:
            # Ignore comment lines.
            if line_str.startswith('#'):
                continue

            # Split VCF lines on tabs into VCF fields.
            vcf_field_list = line_str.rstrip().split('\t')
            alleles = '..'

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
                print('Unexpected genotype {!r} in line: {}'.format(vcf_field_list[9], line_str))

            # The identifier (ID) is in field 2,  the chromosome (CHROM) in field 0 and the position (POS) in field 1.
            output_file.write('\t'.join((vcf_field_list[2], vcf_field_list[0], vcf_field_list[1], alleles)) + '\n')
