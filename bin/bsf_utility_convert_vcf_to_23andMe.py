#! /usr/bin/env python
#
# BSF Python script to convert a VCF file into the simple 23andMe format.
#
# The 23andMe format consists of a dbSNP reference SNP identifier (rsid),
# the sequence name, the position and the genotype for each allele (diploid).
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
import string


parser = argparse.ArgumentParser(description='BSF Converter from VCF to 23andMe format.')

parser.add_argument('--debug', required=False, type=int,
                    help='debug level')

parser.add_argument('--input', required=True,
                    help='file path to a VCF file.')

parser.add_argument('--output', required=True,
                    help='file path to a 23andMe file.')

args = parser.parse_args()

input_fh = open(args.input, 'r')
output_fh = open(args.output, 'w')

for line in input_fh:

    # Ignore comment lines.
    if line.startswith('#'):
        continue

    # Split VCF lines on tabs into VCF fields.
    vcf_fields = string.split(s=line.rstrip(), sep='\t')
    alleles = '..'

    # The genotype (GT) is in field 10, the REF in field 4 and ALT in field 5.

    genotype_index = 0

    info_fields = string.split(s=vcf_fields[8], sep=':')

    for i in range(0, len(info_fields) - 1):
        if info_fields[i] == 'GT':
            genotype_index = i
            break

    sample_fields = string.split(vcf_fields[9], sep=':')

    if sample_fields[genotype_index] == '0/0':
        alleles = vcf_fields[3] + vcf_fields[3]
    elif sample_fields[genotype_index] == '0/1':
        alleles = vcf_fields[3] + vcf_fields[4]
    elif sample_fields[genotype_index] == '1/1':
        alleles = vcf_fields[4] + vcf_fields[4]
    elif sample_fields[genotype_index] == './.':
        continue
    else:
        print 'Unexpected genotype {!r} in line: {}'.format(vcf_fields[9], line)

    # The identifier (ID) is in field 2,  the chromosome (CHROM) in field 0 and the position (POS) in field 1.
    output_fh.write(string.join(words=(vcf_fields[2], vcf_fields[0], vcf_fields[1], alleles), sep='\t') + '\n')

input_fh.close()
output_fh.close()
