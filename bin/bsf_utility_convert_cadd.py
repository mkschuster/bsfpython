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
"""The :py:mod:`bin.bsf_utility_convert_cadd` module is a script to
convert a CADD file into a VCF file that can be used by the GATK :literal:`VariantAnnotator`.
"""

import os
import sys
from argparse import ArgumentParser
from subprocess import Popen, PIPE
from threading import Lock
from typing import IO

from bsf.connector import StandardOutputStream
from bsf.process import Executable


def process_stdout(input_file_handle: IO, thread_lock: Lock, output_file_path: str, contig_list: list[str]) -> None:
    """Process :emphasis:`standard output` via a separate :py:class:`threading.Thread` object.

    :param input_file_handle: Standard input stream
    :type input_file_handle: IO
    :param thread_lock: A Python :py:class:`threading.Lock` object
    :type thread_lock: Lock
    :param output_file_path: Output file path
    :type output_file_path: str
    :param contig_list: A Python :py:class:`list` object of Python :py:class:`str` (contig name) objects
    :type contig_list: list[str]
    """
    thread_lock.acquire(True)
    output_binary_io = open(file=output_file_path, mode='wb')
    output_process = Popen(
        args=['bgzip'],
        stdin=PIPE,
        stdout=output_binary_io,
        stderr=None,
        text=True)
    thread_lock.release()

    # Write the new VCF header declaring the file format and the two CADD INFO fields.

    thread_lock.acquire(True)
    output_process.stdin.write('##fileformat=VCFv4.1\n')
    output_process.stdin.write('##INFO=<ID=CADDr,Number=A,Type=Float,'
                               'Description="Combined Annotation Dependent Depletion (CADD) '
                               'raw combined SVM score (C-score) on an arbitrary scale '
                               'depending on the annotations used.">\n')
    output_process.stdin.write('##INFO=<ID=CADDs,Number=A,Type=Float,'
                               'Description="Combined Annotation Dependent Depletion (CADD) '
                               'PHRED-scaled combined SVM score (scaled C-score) ranging from 1 to 99, '
                               'based on the rank of each variant relative to all possible '
                               '8.6 billion substitutions in the human reference genome.">\n')
    thread_lock.release()

    for vcf_line in input_file_handle:
        if vcf_line[:6] == '#Chrom':
            # Replace the CADD column header line with the standard VCF line.
            thread_lock.acquire(True)
            output_process.stdin.writelines(contig_list)
            output_process.stdin.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            thread_lock.release()
            continue
        if vcf_line[:1] == '#':
            # Copy all comment lines from the CADD file.
            thread_lock.acquire(True)
            output_process.stdin.write(vcf_line)
            thread_lock.release()
            continue
        fields = vcf_line.split()
        # GRCh37 has one ambiguity base (M) at position 3:60830534, which GATK does not like.
        if fields[2] not in ('A', 'C', 'G', 'T'):
            continue
        new_line = '\t'.join((
            fields[0],  # CHROM
            fields[1],  # POS
            '.',  # ID (missing value)
            fields[2],  # REF
            fields[3],  # ALT
            '.',  # QUAL (missing value)
            '.',  # FILTER (missing value)
            'CADDr={};CADDs={}'.format(fields[4], fields[5])  # INFO
        ))
        # Since this is the only thread processing STDOUT and writing to this file,
        # the locking may not be required after all.
        thread_lock.acquire(True)
        output_process.stdin.write(new_line + '\n')
        thread_lock.release()

    thread_lock.acquire(True)
    # Close STDIN of the output process to signal end of stream, then wait for it to exit.
    output_process.stdin.close()
    output_return_code = output_process.wait()
    print('Output process return code:', repr(output_return_code))
    output_binary_io.close()
    thread_lock.release()

    return


def run(
        input_path: str,
        output_path: str,
        reference_vcf: str) -> int:
    """Run function.

    :param input_path: A CADD format bgzip compressed input file path.
    :type input_path: str
    :param output_path: A :py:class:`bsf.analysis.Stage` name.
    :type output_path: str
    :param reference_vcf: A reference VCF file path.
    :type reference_vcf: str
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    # Set the environment consistently.

    os.environ['LANG'] = 'C'

    # Read the contig annotation from the reference file.

    contig_list: list[str] = list()

    with open(file=reference_vcf, mode='rt') as input_text_io:
        for line_str in input_text_io:
            if line_str.startswith('##contig'):
                contig_list.append(line_str)
            if line_str.startswith('##reference'):
                contig_list.append(line_str)
            if line_str.startswith('#CHROM'):
                break

    # Now the complicated part, the sub process for reading via bgzip.

    stdin_executable = Executable(
        program='bgzip',
        stdout=StandardOutputStream(
            thread_callable=process_stdout,
            thread_kwargs={'output_file_path': output_path, 'contig_list': contig_list}))
    stdin_executable.add_switch_long(key='decompress')
    stdin_executable.add_switch_long(key='stdout')
    stdin_executable.arguments.append(input_path)

    exception_str_list = stdin_executable.run()

    if exception_str_list:
        exception_str_list.append('Command list representation: ' + repr(stdin_executable.command_list()))
        raise Exception('\n'.join(exception_str_list))

    print('All done.')

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
        help='CADD format bgzip-compressed input file path')

    argument_parser.add_argument(
        '--output-path',
        required=True,
        help='bgzip-compressed VCF output file path')

    argument_parser.add_argument(
        '--reference-vcf',
        required=True,
        help='reference VCF file path')

    name_space = argument_parser.parse_args()

    return run(
        input_path=name_space.input_path,
        output_path=name_space.output_path,
        reference_vcf=name_space.reference_vcf)


if __name__ == '__main__':
    sys.exit(main())
