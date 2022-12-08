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
"""The :py:mod:`bin.bsf_utility_check_md5sums` module is a script to check files against the repository of MD5 sums.
"""

import logging
import os
from argparse import ArgumentParser
from tempfile import NamedTemporaryFile

from bsf.md5sum import MD5SumArchive
from bsf.process import Executable

argument_parser = ArgumentParser(
    description='Update an MD5 sum archive file')

argument_parser.add_argument(
    '--logging-level',
    choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
    default='INFO',
    dest='logging_level',
    help='Logging level [INFO]',
    required=False)

argument_parser.add_argument(
    '--file-path',
    dest='file_path',
    help='MD5 sum archive file path',
    required=True)

argument_parser.add_argument(
    'file_paths',
    help='one or more file paths',
    nargs='+')

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

md5sum_archive = MD5SumArchive.from_file_path(file_path=name_space.file_path)

text_io = NamedTemporaryFile(mode='wt', suffix='.md5', delete=False)

for file_path in name_space.file_paths:
    file_name = os.path.basename(file_path)
    if file_name in md5sum_archive.md5sum_dict:
        md5sum = md5sum_archive.md5sum_dict[file_name]
        md5_file_path = os.path.normpath(os.path.join(os.getcwd(), file_path))
        print(md5sum.check_sum + ' ' + md5sum.check_mode + md5_file_path, file=text_io)

text_io.close()

executable = Executable(name='md5sum', program='md5sum')
executable.add_switch_long(key='check')
executable.arguments.append(text_io.name)
message_list = executable.run()
if message_list:
    print(repr(message_list))

os.remove(text_io.name)
print('All done.')
