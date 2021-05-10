#!/usr/bin/env python3.8
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
#  BSF Python script to collect and update MD5 sums.
#
#  A directory tree (directory-path) is scanned for files matching a regular expression pattern
#  (--pattern), by default \\.md5 and updated into an existing (--file-path) md5sum file.
#  Picard-style MD5 file that only contain the MD5 digest are reformatted to obey the md5sum format.
#
import os
import re
import warnings
from argparse import ArgumentParser
from typing import Dict, Tuple


def split_md5sum_line(md5sum_str):
    """Split a line of a md5sum file into check_sum, check_mode (i.e. ' ' for text or '*' for binary) and file name.

    @param md5sum_str: md5sum string
    @type md5sum_str: str
    @return: Tuple of check_sum, check_mode and file_name
    @rtype: (str, str, str)
    """
    md5sum_str = md5sum_str.strip()

    if ' ' in md5sum_str:
        _index_int = md5sum_str.index(' ')

        # The MD5 check sum lies up until the index location.
        _check_sum = md5sum_str[:_index_int]

        # The check mode marker for text (' ') or binary ('*') lies after the index location.
        _check_mode = md5sum_str[_index_int + 1:_index_int + 2]

        # The file name excluding the marker follows the index location + 1.
        _file_name = md5sum_str[_index_int + 2:]
    else:
        _check_sum = md5sum_str
        _check_mode = '*'
        _file_name = ''

    return _check_sum, _check_mode, _file_name


md5_dict: Dict[str, Tuple[str, str]] = dict()


def add_md5_entry(entry_file_name, entry_check_sum, entry_check_mode):
    """Add an dm5sum entry to the global md5sum Python C{dict} object.

    @param entry_file_name: File name
    @type entry_file_name: str
    @param entry_check_sum: Check sum
    @type entry_check_sum: str
    @param entry_check_mode: Check mode
    @type entry_check_mode: str
    """
    if entry_file_name in md5_dict and md5_dict[entry_file_name][0] != entry_check_sum:
        warnings.warn('Non-matching check sum ' + repr(entry_check_sum) + ' for file name: ' +
                      repr(entry_file_name))
    else:
        md5_dict[entry_file_name] = (entry_check_sum, entry_check_mode)

    return


argument_parser = ArgumentParser(
    description='MD5 sum collection and update script.')

argument_parser.add_argument(
    'directory_path',
    help='directory path',
    type=str)

argument_parser.add_argument(
    '--pattern',
    default='\\.md5',
    help='Regular expression pattern [\\.md5]',
    type=str)

argument_parser.add_argument(
    '--file-path',
    dest='file_path',
    help='Existing MD5 sum file path to update',
    required=True,
    type=str)

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

name_space = argument_parser.parse_args()

re_pattern = re.compile(pattern=name_space.pattern)

# Read the initial MD5 sum file that needs updating.

if os.path.exists(name_space.file_path):
    with open(file=name_space.file_path, mode='rt') as input_file:
        for line_str in input_file:
            md5_check_sum, md5_check_mode, md5_file_name = split_md5sum_line(md5sum_str=line_str)

            if name_space.debug > 0:
                print('md5_check_sum:', md5_check_sum,
                      'md5_check_mode:', md5_check_mode,
                      'md5_file_name:', md5_file_name)

            if md5_file_name:
                add_md5_entry(
                    entry_check_sum=md5_check_sum,
                    entry_file_name=md5_file_name,
                    entry_check_mode=md5_check_mode)
            else:
                raise Exception('The md5sum file ' + repr(name_space.file_path) +
                                " does not obey the standard 'MD5SUM *file_path' format.")

# Iterate through the directory path to find new MD5 sum files.

for file_path, directory_name_list, file_name_list in os.walk(top=name_space.directory_path, topdown=True):
    if name_space.debug > 0:
        print('file_path:', file_path)
        print('directory_name_list:', directory_name_list)
        print('file_name_list:', file_name_list)
        print()

    for file_name in file_name_list:
        if re.search(pattern=re_pattern, string=file_name) is None:
            if name_space.debug > 0:
                print('Excluding:', file_name)
            continue

        with open(file=os.path.join(file_path, file_name), mode='rt') as input_file:
            for line_str in input_file:
                md5_check_sum, md5_check_mode, md5_file_name = split_md5sum_line(md5sum_str=line_str)

                if not md5_file_name:
                    # In case the md5sum file does not specify a file name (e.g. Picard MD5 sum files),
                    # use the file name without its '.md5' suffix.
                    if file_name.endswith('.md5'):
                        md5_file_name = file_name[:-4]
                    else:
                        raise Exception('Unexpected suffix of (Picard) md5sum file: ' + repr(file_name))

                md5_file_name = os.path.basename(md5_file_name)

                add_md5_entry(
                    entry_check_sum=md5_check_sum,
                    entry_file_name=md5_file_name,
                    entry_check_mode=md5_check_mode)

# Write the updated MD5 sum file.

with open(file=name_space.file_path, mode='wt') as output_file:
    for md5_file_name in sorted(md5_dict):
        md5_check_sum = md5_dict[md5_file_name][0]
        md5_check_mode = md5_dict[md5_file_name][1]
        # Adjust the mode to binary for certain files.
        for suffix in ('.bam', '.gz'):
            if md5_file_name.endswith(suffix):
                md5_check_mode = '*'

        line_str = ''.join((md5_check_sum, ' ', md5_check_mode, md5_file_name))
        if name_space.debug > 0:
            print(line_str)
        output_file.write(line_str + '\n')
