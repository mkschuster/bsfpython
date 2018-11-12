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
import sys
import time

import bsf.analyses.illumina_to_bam_tools
import bsf.analyses.picard
import bsf.illumina
import bsf.standards

argument_parser = argparse.ArgumentParser(
    description='Illumina Run Folder processor driver script.')

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
    '--illumina2bam',
    action='store_true',
    help='Use Illumina2bam rather than Picard tools')

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
    default=bsf.standards.Configuration.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode i.e. "high" (high-output) or "rapid" (rapid run) or "miseq" for a MiSeq run',
    required=False,
    type=str)

argument_parser.add_argument(
    '--compression',
    help='Picard (Zlib) compression level [9]',
    required=False,
    type=int)

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

if name_space.illumina2bam:
    # Create an IlluminaToBam analysis, run and submit it.

    analysis_itb = bsf.analyses.illumina_to_bam_tools.IlluminaToBam.from_config_file_path(
        config_path=name_space.configuration)
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
            except bsf.illumina.RunFolderNotComplete as exception:
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
        sys.stdout.writelines(analysis_itb.trace(level=1))

    # Create a BamIndexDecoder analysis, run and submit it.

    analysis_bid = bsf.analyses.illumina_to_bam_tools.BamIndexDecoder.from_config_file_path(
        config_path=name_space.configuration)
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
        elif name_space.mode == 'novaseq':
            analysis_bid.lanes = 4
        else:
            raise Exception("Unknown output mode " + name_space.mode)
    else:
        analysis_bid.lanes = \
            bsf.illumina.RunFolder.from_file_path(
                file_path=analysis_itb.run_directory).run_information.flow_cell_layout.lane_count

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
        sys.stdout.writelines(analysis_bid.trace(level=1))
else:
    # Create an IlluminaMultiplexSam analysis, run and submit it.

    analysis_ims = bsf.analyses.picard.IlluminaMultiplexSam.from_config_file_path(config_path=name_space.configuration)
    """ @type analysis_ims: bsf.analyses.picard.IlluminaMultiplexSam """

    # Set arguments that override the configuration file.

    if name_space.debug:
        assert isinstance(name_space.debug, int)
        analysis_ims.debug = name_space.debug

    if name_space.irf:
        assert isinstance(name_space.irf, (str, unicode))
        analysis_ims.run_directory = name_space.irf

    if name_space.compression is not None:
        analysis_ims.compression_level = name_space.compression

    # Do the work.

    if name_space.loop:
        assert isinstance(name_space.loop, bool)
        # If the --loop option has been set, wait until RTAComplete.txt has been copied.
        loop_counter = int(1)
        while True:
            print('[{}] Loop {}:'.format(datetime.datetime.now().isoformat(), loop_counter))
            loop_counter += int(1)
            try:
                analysis_ims.run()
            except bsf.illumina.RunFolderNotComplete as exception:
                print(exception)
            else:
                print('Illumina Run Folder seems complete.')
                break
            time.sleep(name_space.interval)
    else:
        analysis_ims.run()

    analysis_ims.check_state()
    analysis_ims.submit(name=name_space.stage)

    print(analysis_ims.name)
    print('Project name:           ', analysis_ims.project_name)
    print('Project directory:      ', analysis_ims.project_directory)
    print('Illumina run directory: ', analysis_ims.run_directory)
    print('Experiment directory:   ', analysis_ims.get_experiment_directory)

    if analysis_ims.debug >= 2:
        print(repr(analysis_ims), 'final trace:')
        sys.stdout.writelines(analysis_ims.trace(level=1))

    # Create a IlluminaDemultiplexSam analysis, run and submit it.

    analysis_ids = bsf.analyses.picard.IlluminaDemultiplexSam.from_config_file_path(
        config_path=name_space.configuration)
    """ @type analysis_ids: bsf.analyses.picard.IlluminaDemultiplexSam """

    # Transfer the project_name and run_directory from the IlluminaMultiplexSam to the
    # IlluminaDemultiplexSam analysis.

    analysis_ids.project_name = analysis_ims.project_name
    analysis_ids.run_directory = analysis_ims.run_directory

    # Set arguments that override the configuration file.

    if name_space.debug:
        assert isinstance(name_space.debug, int)
        analysis_ids.debug = name_space.debug

    if name_space.compression is not None:
        analysis_ids.compression_level = name_space.compression

    if name_space.mode:
        assert isinstance(name_space.mode, str)
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
            raise Exception("Unknown output mode " + name_space.mode)
    else:
        analysis_ids.lanes = \
            bsf.illumina.RunFolder.from_file_path(
                file_path=analysis_ims.run_directory).run_information.flow_cell_layout.lane_count

    if name_space.no_validation:
        assert isinstance(name_space.no_validation, bool)
        analysis_ids.force = name_space.no_validation

    if name_space.library_path:
        assert isinstance(name_space.library_path, (str, unicode))
        analysis_ids.library_path = name_space.library_path

    # If a library file has not been defined so far, check,
    # if a standard library file i.e. PROJECT_NAME_libraries.csv exists in the current directory.

    if not analysis_ids.library_path:
        library_path = '_'.join((analysis_ids.project_name, 'libraries.csv'))
        if os.path.exists(path=library_path):
            analysis_ids.library_path = library_path

    if analysis_ids.library_path:
        # Do the work if, at this stage, a library file has been set.

        analysis_ids.run()
        analysis_ids.check_state()
        analysis_ids.submit(name=name_space.stage)

        print('')
        print(analysis_ids.name)
        print('Project name:         ', analysis_ids.project_name)
        print('Project directory:    ', analysis_ids.project_directory)
        print('Sequences directory:  ', analysis_ids.sequences_directory)
        print('Experiment directory: ', analysis_ids.get_experiment_directory)

    if analysis_ids.debug >= 2:
        print(repr(analysis_ids), 'final trace:')
        sys.stdout.writelines(analysis_ids.trace(level=1))
