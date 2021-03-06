#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
#
# BSF Python driver script for the conversion of BAM files into FASTQ format.
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
from argparse import ArgumentParser

from bsf.analyses import RunBamToFastq

argument_parser = ArgumentParser(
    description=RunBamToFastq.name + ' driver script.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='Debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--stage',
    help='Limit job submission to a particular Analysis stage',
    required=False)

argument_parser.add_argument(
    'configuration',
    help='Configuration file (*.ini)')

name_space = argument_parser.parse_args()

# Create a RunBamToFastq BSF Analysis and run it.

analysis = RunBamToFastq.from_config_file_path(config_path=name_space.configuration)

if name_space.debug:
    analysis.debug = name_space.debug

# Do the work.

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage)

print(analysis.name)
print('Project name:      ', analysis.project_name)
print('Input directory:   ', analysis.input_directory)
print('Output directory:  ', analysis.output_directory)
print('Project directory: ', analysis.project_directory)
print('Genome directory:  ', analysis.genome_directory)
