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
#  BSF Python utility script to collect processed run folder statistics.
#
import os
import re
import stat
from argparse import ArgumentParser

from bsf.standards import Configuration, StandardFilePath

argument_parser = ArgumentParser(
    description='Collect statistics from a processed run folder after de-multiplexing.')

argument_parser.add_argument(
    'input',
    help='Processed run folder directory')

name_space = argument_parser.parse_args()

prf_path = name_space.input

prf_path = Configuration.get_absolute_path(
    file_path=prf_path,
    default_path=StandardFilePath.get_sequences(absolute=True))

if not os.path.exists(prf_path):
    raise Exception('Could not find processed run folder directory {!r}'.format(prf_path))

for file_name in os.listdir(prf_path):
    file_path = os.path.join(prf_path, file_name)
    mode = os.stat(file_path).st_mode
    match = re.search(pattern=r'^(\d+)$', string=file_name)
    if stat.S_ISDIR(mode) and match:
        # This is the lane directory. Should be changed from 1 to L001 for CASAVA compatibility ...
        for file_name_2 in os.listdir(file_path):
            file_path_2 = os.path.join(file_path, file_name_2)
            mode = os.stat(file_path_2).st_mode
            match = re.search(pattern=r'([^.]+).([^.]+).output.metrics.txt', string=file_name_2)
            if stat.S_ISREG(mode) and match:
                with open(file=file_path_2, mode='rt') as metrics_file:
                    for line_str in metrics_file:
                        if not line_str:
                            continue

                        if line_str.startswith('#'):
                            continue

                        values = line_str.split()

                        print('Line:', ' '.join(values))
