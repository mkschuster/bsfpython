#! /usr/bin/env python
#
# BSF Python script to archive an Illumina Run Folder (IRF) from a magnetic tape library.
#
#
# Copyright 2013 - 2015 Michael K. Schuster
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

from bsf import Default
from bsf.analyses.illumina_run_folder import IlluminaRunFolderArchive


argument_parser = ArgumentParser(
    description='IlluminaRunFolderArchive Analysis driver script.')

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

irf_archive = IlluminaRunFolderArchive.from_config_file_path(config_path=name_space.configuration)

# Set arguments that override the configuration file.

if name_space.debug:
    irf_archive.debug = name_space.debug

if name_space.project_name:
    irf_archive.project_name = name_space.project_name

if name_space.archive_directory:
    irf_archive.archive_directory = name_space.archive_directory

if name_space.irf:
    irf_archive.run_directory = name_space.irf

if name_space.force:
    irf_archive.force = name_space.force

irf_archive.run()
irf_archive.submit(drms_name=name_space.stage)

print 'IlluminaRunFolderArchive Analysis'
print 'Project name:           ', irf_archive.project_name
print 'Project directory:      ', irf_archive.project_directory
print 'Illumina run directory: ', irf_archive.run_directory
print 'Archive directory:      ', irf_archive.archive_directory