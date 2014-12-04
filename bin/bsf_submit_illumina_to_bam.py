#! /usr/bin/env python
#
# BSF Python script to drive the IlluminaToBamTools IlluminaToBam analysis.
#
#
# Copyright 2014 Michael K. Schuster
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
from bsf.analyses.illumina_to_bam_tools import IlluminaToBam


argument_parser = ArgumentParser(
    description='IlluminaToBamTools Illumina2bam analysis driver script.')

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
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--configuration',
    default=Default.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing of an incomplete Illumina Run Folder')

name_space = argument_parser.parse_args()

# Create a BSF IlluminaToBam analysis, run and submit it.

itb = IlluminaToBam.from_config_file_path(config_path=name_space.configuration)

# Set arguments that override the configuration file.

if name_space.debug:
    itb.debug = name_space.debug

if name_space.irf:
    itb.illumina_run_folder = name_space.irf

if name_space.force:
    itb.force = name_space.force

# Do the work.

itb.run()
itb.submit(drms_name=name_space.stage)

print 'IlluminaToBamTools IlluminaToBam Analysis'
print 'Project name:         ', itb.project_name
print 'Project directory:    ', itb.project_directory
print 'Illumina Run Folder:  ', itb.illumina_run_folder
print 'Experiment directory: ', itb.experiment_directory
