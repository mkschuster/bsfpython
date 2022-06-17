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
#
#  BSF Python script to collect and update MD5 sums.
#
#  A directory tree (directory-path) is scanned for file names matching a regular expression pattern
#  (--file-pattern), by default \\.md5 and updated into an existing (--file-path) md5sum file.
#  Picard-style MD5 file that only contain the MD5 digest are reformatted to obey the md5sum format.
#
import logging
import os
import re
from argparse import ArgumentParser

from bsf.md5sum import MD5Sum, MD5SumArchive


def process_md5_files(md5sum_archive, directory_path, file_pattern):
    """Process individual MD5 sum files.

    :param md5sum_archive: A :py:class:`bsf.md5sum.MD5SumArchive` object.
    :type md5sum_archive: MD5SumArchive
    :param directory_path: A directory file path.
    :type directory_path: str
    :param file_pattern: A file pattern regular expression.
    :type file_pattern: str
    """
    re_pattern = re.compile(pattern=file_pattern)

    for directory_path, directory_name_list, file_name_list in os.walk(top=directory_path):
        logging.debug("directory_path: '%s'", directory_path)
        logging.debug("directory_name_list: '%s'", directory_name_list)
        logging.debug("file_name_list: '%s'", file_name_list)

        for file_name in file_name_list:
            if re_pattern.search(string=file_name) is None:
                logging.debug("Excluding: '%s'", file_name)
                continue

            with open(file=os.path.join(directory_path, file_name), mode='rt') as text_io:
                for line_str in text_io:
                    md5sum = MD5Sum.from_line(md5sum_str=line_str)

                    if not md5sum.file_path:
                        # In case the md5sum file does not specify a file name (e.g., Picard MD5 sum files),
                        # use the file name without its '.md5' suffix.
                        if file_name.endswith('.md5'):
                            md5sum.file_path = file_name[:-4]
                        else:
                            raise Exception('Unexpected suffix of (Picard) md5sum file: ' + repr(file_name))

                    # Archive just the base name of the file path.
                    md5sum.file_path = os.path.basename(md5sum.file_path)

                    md5sum_archive.add_md5sum(md5sum=md5sum)

    return


argument_parser = ArgumentParser(
    description='Update an MD5 sum archive file')

argument_parser.add_argument(
    'directory_path',
    help='directory path',
    type=str)

argument_parser.add_argument(
    '--file-pattern',
    default='\\.md5$',
    dest='file_pattern',
    help='File name regular expression pattern [\\.md5$]',
    type=str)

argument_parser.add_argument(
    '--file-path',
    dest='file_path',
    help='MD5 sum archive file path',
    required=True,
    type=str)

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

name_space = argument_parser.parse_args()

if name_space.debug:
    if name_space.debug > 1:
        logging.basicConfig(level=logging.DEBUG)
    elif name_space.debug > 0:
        logging.basicConfig(level=logging.INFO)

archive_object = MD5SumArchive.from_file_path(file_path=name_space.file_path)

process_md5_files(
    md5sum_archive=archive_object,
    directory_path=name_space.directory_path,
    file_pattern=name_space.file_pattern)

archive_object.to_file_path()
