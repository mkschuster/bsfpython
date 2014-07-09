#! /usr/bin/env python
#
# BSF Python script to summarise an Illumina Run Folder.
#
#
# Copyright 2013 Michael K. Schuster
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
import datetime
import os.path
from stat import *
import string

from Bio.BSF import Default
from Bio.BSF.Illumina import RunFolder


parser = argparse.ArgumentParser(description='Summarise an Illumina Run Folder.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('file_path',
                    help='File path to an Illumina Run Folder')

args = parser.parse_args()

file_path = str(args.file_path)
file_path = os.path.expanduser(path=file_path)
file_path = os.path.expandvars(path=file_path)

if not os.path.isabs(file_path):
    file_path = os.path.join(Default.absolute_runs_illumina(), file_path)

irf = RunFolder.from_file_path(file_path=file_path)

run_format = list()

data_read1 = irf.run_parameters.get_read1()
if data_read1 and data_read1 != '0':
    run_format.append('{}B'.format(data_read1))

index_read1 = irf.run_parameters.get_index_read1()
if index_read1 and index_read1 != '0':
    run_format.append('{}I'.format(index_read1))

index_read2 = irf.run_parameters.get_index_read2()
if index_read2 and index_read2 != '0':
    run_format.append('{}I'.format(index_read2))

data_read2 = irf.run_parameters.get_read2()
if data_read2 and data_read2 != '0':
    run_format.append('{}B'.format(data_read2))

print 'Flow Cell Identifier: {}_{}'.format(irf.run_parameters.get_experiment_name(),
                                           irf.run_parameters.get_flow_cell_barcode())
print 'Read Structure: {}'.format(string.join(run_format, sep=' + '))

file_path_start = os.path.join(file_path, 'runParameters.xml')
if os.path.exists(file_path_start):
    print 'Start date:     {}'.format(datetime.date.fromtimestamp(os.stat(file_path_start).st_mtime))

file_path_end = os.path.join(file_path, 'RTAComplete.txt')
if os.path.exists(file_path_end):
    print 'End date:       {}'.format(datetime.date.fromtimestamp(os.stat(file_path_end).st_mtime))

print 'Experiment:     {}'.format(irf.run_parameters.get_experiment_name())
print 'Flow Cell:      {}'.format(irf.run_parameters.get_flow_cell_barcode())

position = irf.run_parameters.get_position()
if position:
    print 'Position:       {}'.format(irf.run_parameters.get_position())

print 'Run Identifier: {}'.format(irf.run_information.run_identifier)
