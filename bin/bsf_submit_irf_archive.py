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
#  BSF Python script to archive an Illumina Run Folder (IRF) from a magnetic tape library.
#
import sys
from argparse import ArgumentParser

from bsf.analyses.illumina_run_folder import IlluminaRunFolderArchive
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=IlluminaRunFolderArchive.name + ' driver script.')

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
    help='project name i.e. instrument run identifier',
    required=False,
    type=str)

argument_parser.add_argument(
    '--archive-directory',
    dest='archive_directory',
    help='archive directory',
    required=False,
    type=str)

argument_parser.add_argument(
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing even if a run folder exists already',
    required=False)

name_space = argument_parser.parse_args()

# Create a BSF IlluminaRunFolderRestore analysis, run and submit it.

analysis = IlluminaRunFolderArchive.from_config_file_path(
    config_path=name_space.configuration)

# Set arguments that override the configuration file.

if name_space.debug:
    analysis.debug = name_space.debug

if name_space.project_name:
    analysis.project_name = name_space.project_name

if name_space.archive_directory:
    analysis.archive_directory = name_space.archive_directory

if name_space.irf:
    analysis.run_directory = name_space.irf

if name_space.force:
    analysis.force = name_space.force

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage)

print(analysis.name)
print('Project name:           ', analysis.project_name)
print('Project directory:      ', analysis.project_directory)
print('Illumina run directory: ', analysis.run_directory)
print('Archive directory:      ', analysis.archive_directory)

if analysis.debug >= 2:
    print('{!r} final trace:'.format(analysis))
    sys.stdout.writelines(analysis.trace(level=1))
