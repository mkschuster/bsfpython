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
#  BSF Python script to drive the Trimmomatic analysis pipeline.
#
import logging
from argparse import ArgumentParser

from bsf.analyses.trimmomatic import Trimmomatic
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=Trimmomatic.name + ' driver script.')

argument_parser.add_argument(
    '--dry-run',
    action='store_false',
    default=True,
    dest='drms_submit',
    help='dry run',
    required=False)

argument_parser.add_argument(
    '--logging-level',
    choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
    default='INFO',
    dest='logging_level',
    help='Logging level [INFO]',
    required=False)

argument_parser.add_argument(
    '--stage',
    help='limit job submission to a particular Analysis stage',
    required=False)

argument_parser.add_argument(
    '--configuration',
    default=Configuration.global_file_path,
    help='configuration (*.ini) file path',
    required=False)

argument_parser.add_argument(
    '--project-name',
    dest='project_name',
    help='project name',
    required=False)

argument_parser.add_argument(
    '--sas-file',
    dest='sas_file',
    help='sample annotation sheet (*.csv) file path',
    required=False)

name_space = argument_parser.parse_args()

# This analysis requires either a non-default --configuration argument or a
# --project-name and --sas-file argument.

if name_space.configuration == Configuration.global_file_path:
    if not name_space.project_name:
        raise Exception("argument --project-name is required if --configuration is not set")
    if not name_space.sas_file:
        raise Exception("argument --sas-file is required if --configuration is not set")

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

analysis = Trimmomatic.from_config_file_path(config_path=name_space.configuration)

if name_space.project_name:
    analysis.project_name = name_space.project_name

if name_space.sas_file:
    analysis.sas_file = name_space.sas_file

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

print(analysis.name)
print('Project name:      ', analysis.project_name)
print('Input directory:   ', analysis.input_directory)
print('Project directory: ', analysis.project_directory)
print('Genome directory:  ', analysis.genome_directory)
if name_space.stage == 'report':
    print('Report URL:        ', analysis.get_html_report_url())
