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
#  BSF Python utility script to reverse barcodes in a library annotation sheet.
#
from argparse import ArgumentParser

from Bio.Seq import reverse_complement

from bsf.analyses.illumina_to_bam_tools import LibraryAnnotationSheet

argument_parser = ArgumentParser(
    description='BSF Python utility script to reverse complement a Library Annotation Sheet.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='Debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--input-path',
    dest='input_path',
    help='Input Library Annotation Sheet',
    required=True,
    type=str)

argument_parser.add_argument(
    '--output-path',
    dest='output_path',
    help='Output Library Annotation Sheet',
    required=True,
    type=str)

argument_parser.add_argument(
    '--lane-list',
    dest='lane_list',
    help='Comma-separated list of lanes to reverse complement',
    required=False,
    type=str)

argument_parser.add_argument(
    '--barcode_1',
    action='store_true',
    dest='barcode_1',
    help='Reverse complement barcode 1')

argument_parser.add_argument(
    '--barcode_2',
    action='store_true',
    dest='barcode_2',
    help='Reverse complement barcode 2')

name_space = argument_parser.parse_args()

if name_space.lane_list:
    lane_list = name_space.lane_list.split(',')
else:
    lane_list = list()

library_annotation_sheet = LibraryAnnotationSheet.from_file_path(
    file_path=name_space.input_path)

for row_dict in library_annotation_sheet.row_dicts:
    if not lane_list or row_dict['lane'] in lane_list:
        if row_dict['barcode_sequence_1'] and name_space.barcode_1:
            row_dict['barcode_sequence_1'] = reverse_complement(sequence=row_dict['barcode_sequence_1'])

        if row_dict['barcode_sequence_2'] and name_space.barcode_2:
            row_dict['barcode_sequence_2'] = reverse_complement(sequence=row_dict['barcode_sequence_2'])

library_annotation_sheet.file_path = name_space.output_path
library_annotation_sheet.to_file_path(adjust_field_names=True)
