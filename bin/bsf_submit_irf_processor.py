#! /usr/bin/env python
#
# BSF Python script to process an Illumina Run Folder (IRF) after sequencing and
# drive the IlluminaToBamTools IlluminaToBam and BamIndexDecoder analyses.
#
#
# Copyright 2013 - 2018 Michael K. Schuster
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

from __future__ import print_function

import argparse
import datetime
import os
import time

from bsf.analyses.illumina_to_bam_tools import BamIndexDecoder, IlluminaToBam
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.standards import Configuration

argument_parser = argparse.ArgumentParser(
    description='IlluminaToBamTools Illumina2bam and BamIndexDecoder Analysis driver script.')

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

# argument_parser.add_argument(
#     '--archive-directory',
#     dest='archive_directory',
#     help='archive directory',
#     required=False,
#     type=str)
#
# argument_parser.add_argument(
#     '--project-name',
#     dest='project_name',
#     help='project name i.e. flow cell identifier',
#     required=False,
#     type=str)
#
# argument_parser.add_argument(
#     '--extract-intensities',
#     action='store_true',
#     dest='extract_intensities',
#     help='extract cluster intensity (*.cif) files',
#     required=False)
#
# argument_parser.add_argument(
#     '--force',
#     action='store_true',
#     help='force expanding a run folder even if it exists already')

argument_parser.add_argument(
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--library-path',
    dest='library_path',
    help='library annotation sheet file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--configuration',
    default=Configuration.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode i.e. "high" (high-output) or "rapid" (rapid run) or "miseq" for a MiSeq run',
    required=False,
    type=str)

argument_parser.add_argument(
    '--no-validation',
    action='store_true',
    dest='no_validation',
    help='force processing even if library annotation sheet validation fails')

argument_parser.add_argument(
    '--loop',
    action='store_true',
    help='loop until a RTAComplete.txt file has been copied by the Illumina Real-Time Analysis (RTA) software',
    required=False)

argument_parser.add_argument(
    '--interval',
    help='loop interval [s]',
    default=120,
    type=int)

name_space = argument_parser.parse_args()

# if name_space.archive_directory:
#
#     # Extract the Illumina Run Folder first.
#
#     irf_restore = IlluminaRunFolderRestore.from_config_file_path(config_path=name_space.configuration)
#
#     if name_space.debug:
#         irf_restore.debug = name_space.debug
#
#     if name_space.project_name:
#         irf_restore.project_name = name_space.project_name
#
#     if name_space.archive_directory:
#         irf_restore.archive_directory = name_space.archive_directory
#
#     if name_space.extract_intensities:
#         irf_restore.extract_intensities = name_space.extract_intensities
#
#     if name_space.force:
#         irf_restore.force = name_space.force
#
#     irf_restore.run()
#     irf_restore.submit(name=name_space.stage)
#
#     print(irf_restore.name)
#     print('Project name:           ', irf_restore.project_name)
#     print('Project directory:      ', irf_restore.project_directory)
#     print('Illumina run directory: ', irf_restore.illumina_directory)
#     print('Archive directory:      ', irf_restore.archive_directory)
# else:
#     irf_restore = None

# Create a BSF IlluminaToBam analysis, run and submit it.

analysis_itb = IlluminaToBam.from_config_file_path(config_path=name_space.configuration)
""" @type analysis_itb: bsf.analyses.illumina_to_bam_tools.IlluminaToBam """

# Set arguments that override the configuration file.

if name_space.debug:
    assert isinstance(name_space.debug, int)
    analysis_itb.debug = name_space.debug

if name_space.irf:
    assert isinstance(name_space.irf, (str, unicode))
    analysis_itb.run_directory = name_space.irf

# if irf_restore is not None:
#     # If the IlluminaRunFolderRestore has not run, the run folder is not complete.
#     itb.run_directory = irf_restore.get_run_directory_path
#     itb.force = True

# Do the work.

if name_space.loop:
    assert isinstance(name_space.loop, bool)
    # If the --loop option has been set, wait until RTAComplete.txt has been copied.
    loop_counter = int(1)
    while True:
        print('[{}] Loop {}:'.format(datetime.datetime.now().isoformat(), loop_counter))
        loop_counter += int(1)
        try:
            analysis_itb.run()
        except RunFolderNotComplete as exception:
            print(exception)
        else:
            print('Illumina Run Folder seems complete.')
            break
        time.sleep(name_space.interval)
else:
    analysis_itb.run()

analysis_itb.check_state()
analysis_itb.submit(name=name_space.stage)

print(analysis_itb.name)
print('Project name:           ', analysis_itb.project_name)
print('Project directory:      ', analysis_itb.project_directory)
print('Illumina run directory: ', analysis_itb.run_directory)
print('Experiment directory:   ', analysis_itb.get_experiment_directory)

if analysis_itb.debug >= 2:
    print(repr(analysis_itb), 'final trace:')
    print(analysis_itb.trace(level=1))

# Create a BSF BamIndexDecoder analysis, run and submit it.

analysis_bid = BamIndexDecoder.from_config_file_path(config_path=name_space.configuration)
""" @type analysis_bid: bsf.analyses.illumina_to_bam_tools.BamIndexDecoder """

# Transfer the project name from the IlluminaToBam to the BamIndexDecoder analysis.

analysis_bid.project_name = analysis_itb.project_name

# Set arguments that override the configuration file.

if name_space.debug:
    assert isinstance(name_space.debug, int)
    analysis_bid.debug = name_space.debug

if name_space.mode:
    assert isinstance(name_space.mode, str)
    if name_space.mode == 'high':
        analysis_bid.lanes = 8
    elif name_space.mode == 'rapid':
        analysis_bid.lanes = 2
    elif name_space.mode == 'miseq':
        analysis_bid.lanes = 1
    elif name_space.mode == 'nextseq':
        analysis_bid.lanes = 4
    else:
        raise Exception("Unknown output mode " + name_space.mode)
else:
    analysis_bid.lanes = \
        RunFolder.from_file_path(file_path=analysis_itb.run_directory).run_information.flow_cell_layout.lane_count

if name_space.no_validation:
    assert isinstance(name_space.no_validation, bool)
    analysis_bid.force = name_space.no_validation

if name_space.library_path:
    assert isinstance(name_space.library_path, (str, unicode))
    analysis_bid.library_path = name_space.library_path

# If a library file has not been defined so far, check,
# if a standard library file i.e. PROJECT_NAME_libraries.csv exists in the current directory.

if not analysis_bid.library_path:
    library_path = '_'.join((analysis_bid.project_name, 'libraries.csv'))
    if os.path.exists(path=library_path):
        analysis_bid.library_path = library_path

if analysis_bid.library_path:
    # Do the work if, at this stage, a library file has been set.

    analysis_bid.run()
    analysis_bid.check_state()
    analysis_bid.submit(name=name_space.stage)

    print('')
    print(analysis_bid.name)
    print('Project name:         ', analysis_bid.project_name)
    print('Project directory:    ', analysis_bid.project_directory)
    print('Sequences directory:  ', analysis_bid.sequences_directory)
    print('Experiment directory: ', analysis_bid.get_experiment_directory)

if analysis_bid.debug >= 2:
    print(repr(analysis_bid), 'final trace:')
    print(analysis_bid.trace(level=1))
