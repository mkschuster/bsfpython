#! /usr/bin/env python
#
# BSF Python script to check a library annotation sheet used by
# the BamIndexDecoder module.
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
import warnings


from Bio.BSF.Data import BamIndexDecoderSheet

parser = argparse.ArgumentParser(description='BSF library annotation sheet checking script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('file_path',
                    help='Library annotation file (*.csv)')

args = parser.parse_args()

annotation_sheet = BamIndexDecoderSheet.read_from_file(file_path=args.file_path)

messages = annotation_sheet.validate()

if messages:
    warnings.warn('\n' + messages)
