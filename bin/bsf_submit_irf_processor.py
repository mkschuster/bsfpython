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
#  BSF Python script to process an Illumina Run Folder (IRF) after sequencing and
#  drive the IlluminaToBamTools IlluminaToBam and BamIndexDecoder analyses.
#
import datetime
import logging
import os
import time
from argparse import ArgumentParser

from bsf.analyses.illumina_to_bam_tools import IlluminaToBam, BamIndexDecoder
from bsf.analyses.picard import IlluminaMultiplexSam, IlluminaDemultiplexSam
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.standards import Configuration

argument_parser = ArgumentParser(
    description='Illumina Run Folder processor driver script.')

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
    '--illumina2bam',
    action='store_true',
    help='Use Illumina2bam rather than Picard tools',
    required=False)

# argument_parser.add_argument(
#     '--archive-directory',
#     dest='archive_directory',
#     help='archive directory',
#     required=False)
#
# argument_parser.add_argument(
#     '--project-name',
#     dest='project_name',
#     help='project name (i.e., flow cell identifier)',
#     required=False)
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
#     help='force expanding a run folder even if it exists already',
#     required=False)

argument_parser.add_argument(
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False)

argument_parser.add_argument(
    '--library-path',
    dest='library_path',
    help='library annotation sheet file path',
    required=False)

argument_parser.add_argument(
    '--configuration',
    default=Configuration.get_global_file_path(),
    help=f'configuration (*.ini) file path [{Configuration.get_global_file_path()!s}]',
    required=False)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode (i.e., "high" (high-output) or "rapid" (rapid run) or "miseq" for a MiSeq run)',
    required=False)

argument_parser.add_argument(
    '--compression',
    help='Picard (Zlib) compression level [9]',
    required=False,
    type=int)

argument_parser.add_argument(
    '--no-validation',
    action='store_true',
    dest='no_validation',
    help='force processing even if library annotation sheet validation fails',
    required=False)

argument_parser.add_argument(
    '--loop',
    action='store_true',
    help='loop until a RTAComplete.txt file has been copied by the Illumina Real-Time Analysis (RTA) software',
    required=False)

argument_parser.add_argument(
    '--interval',
    help='loop interval [s]',
    default=120,
    required=False,
    type=int)

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

# if name_space.archive_directory:
#
#     # Extract the Illumina Run Folder first.
#
#     irf_restore = IlluminaRunFolderRestore.from_config_file_path(config_path=name_space.configuration)
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
#     irf_restore.submit(name=name_space.stage, drms_submit=name_space.drms_submit)
#
#     print(irf_restore.name)
#     print('Project name:           ', irf_restore.project_name)
#     print('Project directory:      ', irf_restore.project_directory)
#     print('Illumina run directory: ', irf_restore.illumina_directory)
#     print('Archive directory:      ', irf_restore.archive_directory)
# else:
#     irf_restore = None

if name_space.illumina2bam:
    analysis_itb = IlluminaToBam.from_config_file_path(config_path=name_space.configuration)

    if name_space.irf:
        analysis_itb.run_directory = name_space.irf

    # if irf_restore is not None:
    #     # If the IlluminaRunFolderRestore has not run, the run folder is not complete.
    #     itb.run_directory = irf_restore.get_run_directory_path
    #     itb.force = True

    if name_space.loop:
        # If the --loop option has been set, wait until RTAComplete.txt has been copied.
        loop_counter = 1
        while True:
            print('[{}] Loop {}:'.format(datetime.datetime.now().isoformat(), loop_counter))
            loop_counter += 1
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
    analysis_itb.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

    print(analysis_itb.name)
    print('Project name:           ', analysis_itb.project_name)
    print('Project directory:      ', analysis_itb.project_directory)
    print('Illumina run directory: ', analysis_itb.run_directory)
    print('Experiment directory:   ', analysis_itb.get_experiment_directory)

    analysis_bid = BamIndexDecoder.from_config_file_path(config_path=name_space.configuration)

    # Transfer the project name from the IlluminaToBam to the BamIndexDecoder analysis.

    analysis_bid.project_name = analysis_itb.project_name

    if name_space.mode:
        if name_space.mode == 'high':
            analysis_bid.lanes = 8
        elif name_space.mode == 'rapid':
            analysis_bid.lanes = 2
        elif name_space.mode == 'miseq':
            analysis_bid.lanes = 1
        elif name_space.mode == 'nextseq':
            analysis_bid.lanes = 4
        elif name_space.mode == 'novaseq':
            analysis_bid.lanes = 4
        else:
            raise Exception(f'The --mode option {name_space.mode!r} is not supported.')
    else:
        analysis_bid.lanes = RunFolder.from_file_path(
            file_path=analysis_itb.run_directory).run_information.flow_cell_layout.lane_count

    if name_space.no_validation:
        analysis_bid.force = name_space.no_validation

    if name_space.library_path:
        analysis_bid.library_path = name_space.library_path

    # If a library file has not been defined so far, check,
    # if a standard library file (i.e., PROJECT_NAME_libraries.csv) exists in the current directory.

    if not analysis_bid.library_path:
        library_path = '_'.join((analysis_bid.project_name, 'libraries.csv'))
        if os.path.exists(library_path):
            analysis_bid.library_path = library_path

    if analysis_bid.library_path:
        # Do the work if, at this stage, a library file has been set.

        analysis_bid.run()
        analysis_bid.check_state()
        analysis_bid.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

        print('')
        print(analysis_bid.name)
        print('Project name:         ', analysis_bid.project_name)
        print('Project directory:    ', analysis_bid.project_directory)
        print('Sequences directory:  ', analysis_bid.sequences_directory)
        print('Experiment directory: ', analysis_bid.get_experiment_directory)
else:
    analysis_ims = IlluminaMultiplexSam.from_config_file_path(config_path=name_space.configuration)

    if name_space.irf:
        analysis_ims.run_directory = name_space.irf

    if name_space.compression is not None:
        analysis_ims.compression_level = name_space.compression

    if name_space.loop:
        # If the --loop option has been set, wait until RTAComplete.txt has been copied.
        loop_counter = 1
        while True:
            print('[{}] Loop {}:'.format(datetime.datetime.now().isoformat(), loop_counter))
            loop_counter += 1
            try:
                analysis_ims.run()
            except RunFolderNotComplete as exception:
                print(exception)
            else:
                print('Illumina Run Folder seems complete.')
                break
            time.sleep(name_space.interval)
    else:
        analysis_ims.run()

    analysis_ims.check_state()
    analysis_ims.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

    print(analysis_ims.name)
    print('Project name:           ', analysis_ims.project_name)
    print('Project directory:      ', analysis_ims.project_directory)
    print('Illumina run directory: ', analysis_ims.run_directory)
    print('Experiment directory:   ', analysis_ims.get_experiment_directory)

    analysis_ids = IlluminaDemultiplexSam.from_config_file_path(config_path=name_space.configuration)

    # Transfer the project_name and run_directory from the IlluminaMultiplexSam to the
    # IlluminaDemultiplexSam analysis.

    analysis_ids.project_name = analysis_ims.project_name
    analysis_ids.run_directory = analysis_ims.run_directory

    if name_space.compression is not None:
        analysis_ids.compression_level = name_space.compression

    if name_space.mode:
        if name_space.mode == 'high':
            analysis_ids.lanes = 8
        elif name_space.mode == 'rapid':
            analysis_ids.lanes = 2
        elif name_space.mode == 'miseq':
            analysis_ids.lanes = 1
        elif name_space.mode == 'nextseq':
            analysis_ids.lanes = 4
        elif name_space.mode == 'novaseq':
            analysis_ids.lanes = 4
        else:
            raise Exception(f'The --mode option {name_space.mode!r} is not supported.')
    else:
        analysis_ids.lanes = RunFolder.from_file_path(
            file_path=analysis_ims.run_directory).run_information.flow_cell_layout.lane_count

    if name_space.no_validation:
        analysis_ids.force = name_space.no_validation

    if name_space.library_path:
        analysis_ids.library_path = name_space.library_path

    # If a library file has not been defined so far, check,
    # if a standard library file (i.e., PROJECT_NAME_libraries.csv) exists in the current directory.

    if not analysis_ids.library_path:
        library_path = '_'.join((analysis_ids.project_name, 'libraries.csv'))
        if os.path.exists(library_path):
            analysis_ids.library_path = library_path

    if analysis_ids.library_path:
        # Do the work if, at this stage, a library file has been set.

        analysis_ids.run()
        analysis_ids.check_state()
        analysis_ids.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

        print('')
        print(analysis_ids.name)
        print('Project name:         ', analysis_ids.project_name)
        print('Project directory:    ', analysis_ids.project_directory)
        print('Sequences directory:  ', analysis_ids.sequences_directory)
        print('Experiment directory: ', analysis_ids.get_experiment_directory)
