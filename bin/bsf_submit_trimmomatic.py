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
#  BSF Python script to drive the Trimmomatic analysis pipeline.
#
import sys
from argparse import ArgumentParser

from bsf.analyses.trimmomatic import Trimmomatic
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=Trimmomatic.name + ' driver script.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--stage',
    help='limit job submission to a particular Analysis stage',
    required=False,
    type=str)

argument_parser.add_argument(
    '--configuration',
    default=Configuration.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--project-name',
    dest='project_name',
    help='project name',
    required=False,
    type=str)

argument_parser.add_argument(
    '--sas-file',
    dest='sas_file',
    help='sample annotation sheet (*.csv) file path',
    required=False,
    type=str)

name_space = argument_parser.parse_args()

# This analysis requires either a non-default --configuration argument or a
# --project-name and --sas-file argument.

if name_space.configuration == Configuration.global_file_path:
    if name_space.project_name is None:
        raise Exception("argument --project-name is required if --configuration is not set")
    if name_space.sas_file is None:
        raise Exception("argument --sas-file is required if --configuration is not set")

# Create a Trimmomatic Analysis and run it.

analysis = Trimmomatic.from_config_file_path(config_path=name_space.configuration)

if name_space.debug:
    analysis.debug = name_space.debug

if name_space.project_name:
    analysis.project_name = name_space.project_name

if name_space.sas_file:
    analysis.sas_file = name_space.sas_file

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage)

print(analysis.name)
print('Project name:      ', analysis.project_name)
print('Input directory:   ', analysis.input_directory)
print('Output directory:  ', analysis.output_directory)
print('Project directory: ', analysis.project_directory)
print('Genome directory:  ', analysis.genome_directory)

if analysis.debug >= 2:
    print(repr(analysis), 'final trace:')
    sys.stdout.writelines(analysis.trace(level=1))
