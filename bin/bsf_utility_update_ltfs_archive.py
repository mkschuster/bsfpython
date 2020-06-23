#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
#
# BSF Python script to collect and update LTFS archive contents.
#
# A directory tree (directory-path) is scanned for files matching a wildcard pattern
# (--pattern), by default *_log.txt and updated into an existing (--file-path) content file.
# The LTFS content file lists run folder names and a comma-separated list of LTO barcode labels
# in tab-separated value (TSV) format.
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
import fnmatch
import os
from argparse import ArgumentParser

argument_parser = ArgumentParser(
    description='LTFS archive collection and update script.')

argument_parser.add_argument(
    'directory_path',
    help='directory path',
    type=str)

argument_parser.add_argument(
    '--pattern',
    default='*_log.txt',
    help='wildcard pattern (glob) [*_log.txt]',
    type=str)

argument_parser.add_argument(
    '--file-path',
    dest='file_path',
    help='Existing LTO content file path to update',
    required=True,
    type=str)

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

name_space = argument_parser.parse_args()

# The LTFS dict uses IRF keys and list value data of LTFS tapes.
ltfs_dict = dict()
""" @type ltfs_dict: dict[str, list[str]] """

# Read the initial LTFS content file that needs updating.

if os.path.exists(name_space.file_path):
    with open(file=name_space.file_path, mode='rt') as input_file:
        for line_str in input_file:
            line_list = line_str.split()
            if line_list[0] not in ltfs_dict:
                ltfs_dict[line_list[0]] = list()
            tape_list = ltfs_dict[line_list[0]]
            for tape_id in line_list[1].split(','):
                if tape_id not in tape_list:
                    tape_list.append(tape_id)

# Iterate through the directory path to find new LTFS log files.

for file_path, directory_name_list, file_name_list in os.walk(top=name_space.directory_path, topdown=True):
    if name_space.debug > 0:
        print('file_path:', file_path)
        print('directory_name_list:', directory_name_list)
        print('file_name_list:', file_name_list)
        print()

    for file_name in file_name_list:
        if not fnmatch.fnmatch(file_name, name_space.pattern):
            continue

        # BS0048L6_log.txt
        tape_name = file_name[:-8]
        with open(file=os.path.join(file_path, file_name), mode='rt') as input_file:
            for line_str in input_file:
                # ILT30505I Copy /scratch/lab_bsf/archive_irf/BS0048L6/BSF_0648_H7WC2BBXY_3.bam to
                #   /mnt/ltfs/BSF_0648_H7WC2BBXY_3.bam
                if not line_str.startswith('ILT30505I'):
                    continue

                line_list = line_str.split()
                base_name = line_list[2].split('/')[-1]

                if name_space.debug > 0:
                    print('File name:', base_name)

                if base_name not in ltfs_dict:
                    ltfs_dict[base_name] = list()
                tape_list = ltfs_dict[base_name]
                if tape_name not in tape_list:
                    tape_list.append(tape_name)

# Write the updated LTFS content file.

with open(file=name_space.file_path, mode='wt') as output_file:
    for irf_file_name in sorted(ltfs_dict):
        print(irf_file_name, ','.join(sorted(ltfs_dict[irf_file_name])), sep='\t', file=output_file)
