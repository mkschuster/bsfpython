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
#  BSF Python script to drive the Picard ExtractIlluminaRunFolder analysis.
#
import logging
from argparse import ArgumentParser

from bsf.analyses.picard import ExtractIlluminaRunFolder
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=ExtractIlluminaRunFolder.name + ' driver script.')

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
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False)

argument_parser.add_argument(
    '--configuration',
    default=Configuration.global_file_path,
    help='configuration (*.ini) file path',
    required=False)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode i.e., high (high-output) or rapid (rapid run)',
    required=False)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing of an incomplete Illumina Run Folder')

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

analysis = ExtractIlluminaRunFolder.from_config_file_path(config_path=name_space.configuration)

if name_space.irf:
    analysis.run_directory = name_space.irf

if name_space.mode:
    if name_space.mode == 'high':
        analysis.lanes = 8
    elif name_space.mode == 'rapid':
        analysis.lanes = 2
    elif name_space.mode == 'miseq':
        analysis.lanes = 1
    elif name_space.mode == 'nextseq':
        analysis.lanes = 4
    elif name_space.mode == 'novaseq':
        analysis.lanes = 4
    else:
        raise Exception(f'The mode option {name_space.mode!r} is not supported.')

if name_space.force:
    analysis.force = name_space.force

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

print(analysis.name)
print('Project name:           ', analysis.project_name)
print('Project directory:      ', analysis.project_directory)
print('Illumina run directory: ', analysis.run_directory)
