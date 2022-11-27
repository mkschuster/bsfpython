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
#  BSF Python script to restore an Illumina Run Folder (IRF) from a magnetic tape library.
#
import logging
from argparse import ArgumentParser

from bsf.analyses.illumina_run_folder import IlluminaRunFolderRestore
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=IlluminaRunFolderRestore.name + ' driver script.')

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
    default=Configuration.get_global_file_path(),
    help=f'configuration (*.ini) file path [{Configuration.get_global_file_path()!s}]',
    required=False)

argument_parser.add_argument(
    '--project-name',
    dest='project_name',
    help='project name i.e. instrument run identifier',
    required=True)

argument_parser.add_argument(
    '--archive-directory',
    dest='archive_directory',
    help='archive directory',
    required=True)

argument_parser.add_argument(
    '--extract-intensities',
    action='store_true',
    dest='extract_intensities',
    help='extract cluster intensity (*.cif) files',
    required=False)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing even if a run folder exists already',
    required=False)

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

analysis = IlluminaRunFolderRestore.from_config_file_path(config_path=name_space.configuration)

if name_space.project_name:
    project_name: str = name_space.project_name
    if project_name.endswith('.ini'):
        raise Exception('The --project-name option should not be a configuration (*.ini) file.')

    analysis.project_name = project_name

if name_space.archive_directory:
    analysis.archive_directory = name_space.archive_directory

if name_space.extract_intensities:
    analysis.extract_intensities = name_space.extract_intensities

if name_space.force:
    analysis.force = name_space.force

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

print(analysis.name)
print('Project name:           ', analysis.project_name)
print('Project directory:      ', analysis.project_directory)
print('Illumina run directory: ', analysis.illumina_directory)
print('Archive directory:      ', analysis.archive_directory)
