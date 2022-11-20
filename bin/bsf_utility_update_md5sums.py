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
"""The :py:mod:`bin.bsf_utility_update_md5sums` module is a script to collect and update MD5 sums.

A directory tree (:literal:`directory_path`) is scanned for file names matching a regular expression pattern
(:literal:`--file-pattern`), by default :literal:`\\.md5` and updated into an existing
(:literal:`--file-path`) md5sum file.
Picard-style MD5 file that only contain the MD5 digest are reformatted to obey the GNU md5sum format.
"""

import logging
import os
import re
from argparse import ArgumentParser

from bsf.md5sum import MD5Sum, MD5SumArchive


def process_md5_files(md5sum_archive: MD5SumArchive, directory_path: str, file_pattern: str) -> None:
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
        logging.debug('directory_path: %r', directory_path)
        logging.debug('directory_name_list: %r', directory_name_list)
        logging.debug('file_name_list: %r', file_name_list)

        for file_name in file_name_list:
            if re_pattern.search(string=file_name) is None:
                logging.debug('Excluding file_name: %r', file_name)
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
                            raise Exception(f'The (Picard) md5sum file prefix is unsupported: {file_name!r}')

                    # Archive just the base name of the file path.
                    md5sum.file_path = os.path.basename(md5sum.file_path)

                    md5sum_archive.add_md5sum(md5sum=md5sum)

    return


argument_parser = ArgumentParser(
    description='Update an MD5 sum archive file')

argument_parser.add_argument(
    'directory_path',
    help='directory path')

argument_parser.add_argument(
    '--file-pattern',
    default='\\.md5$',
    dest='file_pattern',
    help='File name regular expression pattern [\\.md5$]')

argument_parser.add_argument(
    '--file-path',
    dest='file_path',
    help='MD5 sum archive file path',
    required=True)

argument_parser.add_argument(
    '--logging-level',
    choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
    default='INFO',
    dest='logging_level',
    help='Logging level [INFO]',
    required=False)

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

archive_object = MD5SumArchive.from_file_path(file_path=name_space.file_path)

process_md5_files(
    md5sum_archive=archive_object,
    directory_path=name_space.directory_path,
    file_pattern=name_space.file_pattern)

archive_object.to_file_path()
