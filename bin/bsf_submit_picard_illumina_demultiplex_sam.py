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
#  BSF Python script to drive the Picard IlluminaDemultiplexSam analysis.
#
import sys
from argparse import ArgumentParser

from bsf.analyses.picard import IlluminaDemultiplexSam
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description=IlluminaDemultiplexSam.name + ' driver script.')

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

# FIXME: For the moment, the IRF needs passing in.
argument_parser.add_argument(
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--project-name',
    dest='project_name',
    help='project name i.e. flow cell identifier',
    required=False,
    type=str)

argument_parser.add_argument(
    '--library-path',
    dest='library_path',
    help='Library annotation sheet file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--lane-list',
    default='',
    dest='lane_list',
    help='Comma-separated list of lanes to process',
    required=False,
    type=str)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode i.e. high (high-output) or rapid (rapid run)',
    required=False,
    type=str)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing even if library annotation sheet validation fails')

name_space = argument_parser.parse_args()

# This analysis requires either a non-default --configuration argument or a
# --project-name. The --library-path argument can be worked out on the basis
# of the --project-name.

if name_space.configuration == Configuration.global_file_path:
    if name_space.project_name is None and name_space.irf is None:
        raise Exception("argument --project-name or --irf are required if --configuration is not set")

# Create a IlluminaDemultiplexSam analysis, run and submit it.

analysis = IlluminaDemultiplexSam.from_config_file_path(config_path=name_space.configuration)

# Set arguments that override the configuration file.

if name_space.debug:
    analysis.debug = name_space.debug

if name_space.irf:
    analysis.run_directory = name_space.irf

if name_space.project_name:
    analysis.project_name = name_space.project_name

if name_space.lane_list:
    analysis.lane_list = name_space.lane_list.split(',')

# FIXME: Remove the mode option.
# Since the script requires an IRF option, the lanes can be read from the bsf.illumina.RunFolder.
# TODO: Synchronise to other de-multiplexing modules?
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
        raise Exception("Unknown output mode " + name_space.mode)

if name_space.force:
    analysis.force = name_space.force

if name_space.library_path:
    analysis.library_path = name_space.library_path

# Do the work.

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage)

print(analysis.name)
print('Project name:         ', analysis.project_name)
print('Project directory:    ', analysis.project_directory)
print('Sequences directory:  ', analysis.sequences_directory)
print('Experiment directory: ', analysis.get_experiment_directory)

if analysis.debug >= 2:
    print(repr(analysis), 'final trace:')
    sys.stdout.writelines(analysis.trace(level=1))
