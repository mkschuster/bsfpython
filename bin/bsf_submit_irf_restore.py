#! /usr/bin/env python
#
# BSF Python script to restore an Illumina Run Folder (IRF) from a magnetic tape library.
#
#
# Copyright 2013 - 2016 Michael K. Schuster
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

from argparse import ArgumentParser

from bsf.analyses.illumina_run_folder import IlluminaRunFolderRestore
from bsf.standards import Default


argument_parser = ArgumentParser(
    description='IlluminaRunFolderRestore Analysis driver script.')

argument_parser.add_argument(
    '--debug',
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
    default=Default.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--project-name',
    dest='project_name',
    help='project name i.e. instrument run identifier',
    required=True,
    type=str)

argument_parser.add_argument(
    '--archive-directory',
    dest='archive_directory',
    help='archive directory',
    required=True,
    type=str)

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

# Create a BSF IlluminaRunFolderRestore analysis, run and submit it.

irf_restore = IlluminaRunFolderRestore.from_config_file_path(config_path=name_space.configuration)

# Set arguments that override the configuration file.

if name_space.debug:
    irf_restore.debug = name_space.debug

if name_space.project_name:
    irf_restore.project_name = name_space.project_name

if name_space.archive_directory:
    irf_restore.archive_directory = name_space.archive_directory

if name_space.extract_intensities:
    irf_restore.extract_intensities = name_space.extract_intensities

if name_space.force:
    irf_restore.force = name_space.force

irf_restore.run()
irf_restore.check_state()
irf_restore.submit(name=name_space.stage)

print 'IlluminaRunFolderRestore Analysis'
print 'Project name:           ', irf_restore.project_name
print 'Project directory:      ', irf_restore.project_directory
print 'Illumina run directory: ', irf_restore.illumina_directory
print 'Archive directory:      ', irf_restore.archive_directory
