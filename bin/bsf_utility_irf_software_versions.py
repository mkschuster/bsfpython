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
#  BSF Python utility script to report software versions for Illumina Run Folders.
#
import os
import re
from argparse import ArgumentParser
from xml.etree.ElementTree import ParseError

from bsf.annotation import AnnotationSheet
from bsf.illumina import RunFolder

argument_parser = ArgumentParser(
    description='List software version of Illumina Run Folders.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--directory',
    help='directory of illumina run folders',
    required=True,
    type=str)

argument_parser.add_argument(
    '--output-file',
    dest='output_file',
    help='output (*.csv) file path',
    required=True,
    type=str)

argument_parser.add_argument(
    '--ascending',
    action='store_true',
    help='sort flow cells in ascending order rather than in descending by default')

name_space = argument_parser.parse_args()

file_name_list = os.listdir(name_space.directory)
file_name_list.sort()

if not name_space.ascending:
    file_name_list.reverse()

field_names = [
    'experiment',
    'experiment_name',
    'run_identifier',
    'application_name',
    'application_version',
    'rta_version',
    'picard_read_structure',
    'lane_count',
    'keep_intensities',
]

annotation_sheet = AnnotationSheet(
    file_path=name_space.output_file,
    header=True,
    field_names=field_names)

# Illumina Run Folders obey a pattern and additionally directories with just Sequence Analysis Viewer (SAV)
# information should also be allowed.

irf_pattern = re.compile(pattern=r'^[0-9]{6,6}_.*(?:_sav)?$')

for file_name in file_name_list:
    if name_space.debug:
        print('File name:', repr(file_name))

    # Process just entries that obey the Illumina Run Folder pattern.
    match = re.search(pattern=irf_pattern, string=file_name)
    if not match:
        print('No match:', repr(file_name))
        continue

    file_path = os.path.join(name_space.directory, file_name)
    if not (os.path.exists(os.path.join(file_path, 'runParameters.xml')) or
            os.path.exists(os.path.join(file_path, 'RunParameters.xml'))):
        print('Directory ' + repr(file_name) + ' not an Illumina Run Folder')
        continue

    # Temporarily catch IOError and xml.etree.ElementTree.ParseError exceptions
    # that result from a broken FhGFS file system.
    try:
        irf = RunFolder.from_file_path(file_path=file_path)
    except IOError as exception:
        if name_space.debug:
            print('\t'.join((file_name, '?', '?', '?')))
        continue
    except ParseError:
        if name_space.debug:
            print('\t'.join((file_name, '?', '?', '?')))
        continue

    annotation_sheet.row_dicts.append({
        'experiment': irf.run_parameters.get_experiment_name,
        'experiment_name': '_'.join((irf.run_parameters.get_experiment_name,
                                     irf.run_parameters.get_flow_cell_barcode)),
        'run_identifier': irf.run_parameters.get_run_identifier,
        'application_name': irf.run_parameters.get_application_name,
        'application_version': irf.run_parameters.get_application_version,
        'rta_version': irf.run_parameters.get_real_time_analysis_version,
        'picard_read_structure': irf.run_information.get_picard_read_structure,
        'lane_count': irf.run_information.flow_cell_layout.lane_count,
        'keep_intensities': irf.run_parameters.get_keep_intensities,
    })

annotation_sheet.to_file_path()
