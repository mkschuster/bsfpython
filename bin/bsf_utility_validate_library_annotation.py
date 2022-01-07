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
#  BSF Python utility script to validate a Library Annotation Sheet file
#  used by the BamIndexDecoder Analysis.
#
import warnings
from argparse import ArgumentParser

from bsf.analyses.illumina_to_bam_tools import LibraryAnnotationSheet

argument_parser = ArgumentParser(
    description='BSF Python utility script to validate Library Annotation Sheet files.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='Debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--mode',
    default='high',
    help='HiSeq run mode i.e. "high" (high-output) or "rapid" (rapid run) or "miseq" for a MiSeq run',
    required=False,
    type=str)

argument_parser.add_argument(
    'library_path',
    help='library annotation sheet (*.csv) file path')

name_space = argument_parser.parse_args()

library_annotation_sheet = LibraryAnnotationSheet.from_file_path(
    file_path=name_space.library_path)

if name_space.mode == 'high':
    lanes = 8
elif name_space.mode == 'rapid':
    lanes = 2
elif name_space.mode == 'miseq':
    lanes = 1
elif name_space.mode == 'nextseq':
    lanes = 4
else:
    raise Exception("Unknown output mode " + name_space.mode)

messages = library_annotation_sheet.validate(lanes=lanes)

if messages:
    warnings.warn('\n' + messages)
