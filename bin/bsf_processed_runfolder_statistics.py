#! /usr/bin/env python
#
# BSF Python script to collect processed run folder statistics.
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

import argparse
import os
import re
import stat

from bsf.standards import Default


parser = argparse.ArgumentParser(description='Collect statistics from a processed run folder after de-multiplexing.')

parser.add_argument('input',
                    help='Processed run folder directory')

args = parser.parse_args()

prf_path = args.input

if not os.path.isabs(prf_path):
    prf_path = os.path.join(Default.absolute_sequences(), prf_path)

if not os.path.exists(path=prf_path):
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
                metrics_file = open(file_path_2, 'r')
                for line in metrics_file:

                    if not line:
                        continue

                    match = re.search(pattern='^#', string=line)

                    if match:
                        continue

                    values = line.split()

                    print 'Line: {}'.format(' '.join(values))
