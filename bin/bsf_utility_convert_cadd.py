#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# BSF Python script to convert a CADD file into a VCF file that can be used by the
# GATK VariantAnnotator.
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
import datetime
import os
import subprocess
import sys
import threading

import bsf.process

on_posix = 'posix' in sys.builtin_module_names
max_thread_joins = 10
thread_join_timeout = 10


def process_stdout(input_file_handle, thread_lock, output_file_path, debug=0):
    thread_lock.acquire(True)
    output_file = open(output_file_path, 'wb')
    output_process = subprocess.Popen(
        args=['bgzip'],
        bufsize=-1,
        stdin=subprocess.PIPE,
        stdout=output_file,
        stderr=None,
        shell=False,
        close_fds=on_posix)
    thread_lock.release()

    if debug > 0:
        pass  # Just make debug, which has to be part of the function interface, used.

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
    output_file.close()
    thread_lock.release()

    return


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
    '--input-path',
    dest='input_path',
    help='file path to a CADD format bgzip cpmpressed input file',
    required=True)

argument_parser.add_argument(
    '--output-path',
    dest='output_path',
    help='file path to a bgzip compressed VCF output file',
    required=True)

argument_parser.add_argument(
    '--reference-vcf',
    dest='reference_vcf',
    help="file path to a reference VCF to read 'contig' annotation",
    required=True)

name_space = argument_parser.parse_args()

# Read the contig annotation from the reference file.

contig_list = list()

with open(name_space.reference_vcf, 'rt') as input_file:
    for line_str in input_file:
        if line_str.startswith('##contig'):
            contig_list.append(line_str)
        if line_str.startswith('##reference'):
            contig_list.append(line_str)
        if line_str.startswith('#CHROM'):
            break

# Now the complicated part, the sub process for reading via bgzip.

stdin_command = bsf.process.Command(program='bgzip')
stdin_command.add_switch_long(key='decompress')
stdin_command.add_switch_long(key='stdout')
stdin_command.arguments.append(name_space.input_path)

input_process = subprocess.Popen(
    args=stdin_command.command_list(),
    bufsize=-1,
    stdin=None,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    shell=False,
    close_fds=on_posix)

# Two threads, thread_out and thread_err reading STDOUT and STDERR, respectively,
# should make sure that buffers are not filling up.

input_thread_lock = threading.Lock()

input_thread_out = threading.Thread(
    target=process_stdout,
    kwargs={
        'input_file_handle': input_process.stdout,
        'thread_lock': input_thread_lock,
        'output_file_path': name_space.output_path,
        'debug': name_space.debug,
    })
input_thread_out.daemon = True  # Thread dies with the program.
input_thread_out.start()

input_thread_err = threading.Thread(
    target=bsf.process.Executable.process_stderr,
    kwargs={
        'stderr_handle': input_process.stderr,
        'thread_lock': input_thread_lock,
        'stderr_path': None,
        'debug': name_space.debug,
    })
input_thread_err.daemon = True  # Thread dies with the program.
input_thread_err.start()

# Wait for the child process to finish.

input_return_code = input_process.wait()

thread_join_counter = 0

while input_thread_out.is_alive() and thread_join_counter < max_thread_joins:
    input_thread_lock.acquire(True)
    if name_space.debug > 0:
        print('[{}] Waiting for STDOUT processor to finish.'.format(datetime.datetime.now().isoformat()))
    input_thread_lock.release()

    input_thread_out.join(timeout=thread_join_timeout)
    thread_join_counter += 1

thread_join_counter = 0

while input_thread_err.is_alive() and thread_join_counter < max_thread_joins:
    input_thread_lock.acquire(True)
    if name_space.debug > 0:
        print('[{}] Waiting for STDERR processor to finish.'.format(datetime.datetime.now().isoformat()))
    input_thread_lock.release()

    input_thread_err.join(timeout=thread_join_timeout)
    thread_join_counter += 1

print('All done.')
