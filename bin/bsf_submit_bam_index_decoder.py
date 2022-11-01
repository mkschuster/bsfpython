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
#  BSF Python script to drive the IlluminaToBamTools BamIndexDecoder analysis.
#
import logging
import os
from argparse import ArgumentParser

from bsf.analyses.illumina_to_bam_tools import BamIndexDecoder
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=BamIndexDecoder.name + ' driver script.')

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
    help='project name i.e. flow cell identifier',
    required=False)

argument_parser.add_argument(
    '--library-path',
    dest='library_path',
    help='Library annotation sheet file path',
    required=False)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode i.e. high (high-output) or rapid (rapid run)',
    required=False)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing even if library annotation sheet validation fails',
    required=False)

name_space = argument_parser.parse_args()

# This analysis requires either a non-default --configuration argument or a
# --project-name. The --library-path argument can be worked out on the basis
# of the --project-name.

if name_space.configuration == Configuration.global_file_path:
    if not name_space.project_name:
        raise Exception('The argument --project-name is required if --configuration is not set.')

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

analysis = BamIndexDecoder.from_config_file_path(config_path=name_space.configuration)

if name_space.project_name:
    project_name: str = name_space.project_name
    if project_name.endswith('.ini'):
        raise Exception('The --project-name option should not be a configuration (*.ini) file.')

    analysis.project_name = project_name

if name_space.mode:
    if name_space.mode == 'high':
        analysis.lanes = 8
    elif name_space.mode == 'rapid':
        analysis.lanes = 2
    elif name_space.mode == 'miseq':
        analysis.lanes = 1
    elif name_space.mode == 'nextseq':
        analysis.lanes = 4
    else:
        raise Exception(f'The --mode option {name_space.mode!r} is not supported.')

if name_space.force:
    analysis.force = name_space.force

if name_space.library_path:
    analysis.library_path = name_space.library_path

# If a library file has not been defined so far, check,
# if a standard library file i.e. PROJECT_NAME_libraries.csv exists in the current directory.

if not analysis.library_path:
    library_path = '_'.join((analysis.project_name, 'libraries.csv'))
    if os.path.exists(library_path):
        analysis.library_path = library_path

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

print(analysis.name)
print('Project name:         ', analysis.project_name)
print('Project directory:    ', analysis.project_directory)
print('Sequences directory:  ', analysis.sequences_directory)
print('Experiment directory: ', analysis.get_experiment_directory)
