#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
#
# BSF Python script to summarise an Illumina Run Folder.
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

import bsf.illumina

parser = argparse.ArgumentParser(description='Summarise an Illumina Run Folder.')

parser.add_argument(
    '--debug',
    help='Debug level',
    required=False,
    type=int)

parser.add_argument(
    '--check',
    action='store_true',
    help='check for completeness',
    required=False)

parser.add_argument(
    'file_path',
    help='File path to an Illumina Run Folder')

name_space = parser.parse_args()

file_path = bsf.illumina.RunFolder.absolute_file_path(name=name_space.file_path)
if file_path is None:
    raise Exception("Could not resolve the --file-name value to a valid Illumina Run Folder location.")

irf = bsf.illumina.RunFolder.from_file_path(file_path=file_path)

if not irf.run_parameters.get_experiment_name:
    raise Exception("No experiment name set in the Illumina Run Folder configuration.")

print('Flow Cell Identifier: ', '_'.join((irf.run_parameters.get_experiment_name,
                                          irf.run_parameters.get_flow_cell_barcode)))
print('Read Structure:       ', ' + '.join(irf.run_information.get_read_structure_list))

file_path_start = os.path.join(file_path, 'Config')
if os.path.exists(file_path_start):
    print('Start date:           ', datetime.date.fromtimestamp(os.stat(file_path_start).st_mtime))

file_path_end = os.path.join(file_path, 'RTAComplete.txt')
if os.path.exists(file_path_end):
    print('End date:             ', datetime.date.fromtimestamp(os.stat(file_path_end).st_mtime))

print('Experiment:           ', irf.run_parameters.get_experiment_name)
print('Flow Cell:            ', irf.run_parameters.get_flow_cell_barcode)

position = irf.run_parameters.get_position
if position:
    print('Position:             ', irf.run_parameters.get_position)

print('Run Identifier:       ', irf.run_information.run_identifier)
print('Application Name:     ', irf.run_parameters.get_application_name)
print('Application Version:  ', irf.run_parameters.get_application_version)
print('RTA Version:          ', irf.run_parameters.get_real_time_analysis_version)

flow_cell_type = irf.run_parameters.get_flow_cell_type
if flow_cell_type:
    print('Flow Cell Type:       ', flow_cell_type)

index_type = irf.run_parameters.get_index_type
if index_type:
    print('Index Type:           ', index_type)

pe_type = irf.run_parameters.get_pe_type
if pe_type:
    print('Paired-end Type:      ', pe_type)

sbs_type = irf.run_parameters.get_sbs_type
if sbs_type:
    print('SBS Type:             ', sbs_type)

if name_space.check:
    irf.check(debug=name_space.debug)
