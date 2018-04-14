#! /usr/bin/env python
#
# BSF Python utility script to report the number of samples per lane.
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

from __future__ import print_function

import argparse
import os

from bsf.analyses.illumina_to_bam_tools import LibraryAnnotationSheet

argument_parser = argparse.ArgumentParser(
    description='Count samples per lane based on BamIndexDecoder library annotation files.')

argument_parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--directory',
    help='directory of IlluminaToBam tools BamIndexDecoder library annotation files',
    required=False,
    type=str)

argument_parser.add_argument(
    '--ascending',
    action='store_true',
    help='sort flow cells in ascending order rather than in descending by default')

name_space = argument_parser.parse_args()

print('Lane,Sample Number,Comment')

file_name_list = os.listdir(name_space.directory)
file_name_list.sort(cmp=lambda x, y: cmp(x, y))

if not name_space.ascending:
    file_name_list.reverse()

for file_name in file_name_list:
    if file_name[-14:] == '_libraries.csv':
        sas = LibraryAnnotationSheet.from_file_path(file_path=os.path.join(name_space.directory, file_name))

        lane_dict = dict()
        for row_dict in sas.row_dicts:
            if row_dict['lane'] in lane_dict:
                lane_dict[row_dict['lane']] += 1
            else:
                lane_dict[row_dict['lane']] = 1

        for i in range(1, 9):
            key = str(i)
            lane_name = '_'.join((file_name[:-15], key))
            if key in lane_dict:
                print(lane_name + ',' + str(lane_dict[key]))
            else:
                print(lane_name + ',1,not annotated')
