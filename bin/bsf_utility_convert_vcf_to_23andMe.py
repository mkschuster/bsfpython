#! /usr/bin/env python
#
# BSF Python script to convert a VCF file into the simple 23andMe format.
#
# The 23andMe format consists of a dbSNP reference SNP identifier (rsid),
# the sequence name, the position and the genotype for each allele (diploid).
#
#
# Copyright 2013 - 2016 Michael K. Schuster
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

argument_parser = argparse.ArgumentParser(
    description='BSF utility to convert from VCF to 23andMe format.')

argument_parser.add_argument(
    '--debug',
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

input_fh = open(name_space.input, 'r')
output_fh = open(name_space.output, 'w')

for line in input_fh:

    # Ignore comment lines.
    if line.startswith('#'):
        continue

    # Split VCF lines on tabs into VCF fields.
    vcf_fields = line.rstrip().split('\t')
    alleles = '..'

    # The genotype (GT) is in field 10, the REF in field 4 and ALT in field 5.

    genotype_index = 0

    info_fields = vcf_fields[8].split(':')

    for i in range(0, len(info_fields)):
        if info_fields[i] == 'GT':
            genotype_index = i
            break

    sample_fields = vcf_fields[9].split(':')

    if sample_fields[genotype_index] == '0/0':
        alleles = vcf_fields[3] + vcf_fields[3]
    elif sample_fields[genotype_index] == '0/1':
        alleles = vcf_fields[3] + vcf_fields[4]
    elif sample_fields[genotype_index] == '1/1':
        alleles = vcf_fields[4] + vcf_fields[4]
    elif sample_fields[genotype_index] == './.':
        continue
    else:
        print('Unexpected genotype {!r} in line: {}'.format(vcf_fields[9], line))

    # The identifier (ID) is in field 2,  the chromosome (CHROM) in field 0 and the position (POS) in field 1.
    output_fh.write('\t'.join((vcf_fields[2], vcf_fields[0], vcf_fields[1], alleles)) + '\n')

input_fh.close()
output_fh.close()
