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
"""The :py:mod:`bin.bsf_submit_irf_processor` module is a script to drive the
`Picard <http://broadinstitute.github.io/picard/>`_
:py:class:`bsf.analyses.picard.IlluminaMultiplexSam` and
:py:class:`bsf.analyses.picard.IlluminaDemultiplexSam` analyses or if the :literal:`--illumina2bam` option is set, the
`Illumina2bam <https://github.com/wtsi-npg/illumina2bam>`_
:py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBam` and
:py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` analyses.

This analysis requires either a positional :literal:`configuration` argument or
a :literal:`--irf` argument.
"""

import datetime
import logging
import os
import sys
import time
from argparse import ArgumentParser
from typing import Optional

from bsf.analyses.illumina_to_bam_tools import IlluminaToBam, BamIndexDecoder
from bsf.analyses.picard import IlluminaMultiplexSam, IlluminaDemultiplexSam
from bsf.illumina import RunFolder, RunFolderNotComplete
from bsf.standards import Configuration


def run(
        configuration_path: str,
        stage_name: Optional[str] = None,
        drms_submit: Optional[bool] = None,
        irf_path: Optional[str] = None,
        library_path: Optional[str] = None,
        mode: Optional[str] = None,
        illumina2bam: Optional[bool] = None,
        no_validation: Optional[bool] = None,
        compression_level: Optional[int] = None,
        loop: Optional[bool] = None,
        interval: Optional[float] = None) -> int:
    """Run function.

    :param configuration_path: A UNIX configuration (*.ini) file path.
    :type configuration_path: str
    :param stage_name: A :py:class:`bsf.analysis.Stage` name.
    :type stage_name: str | None
    :param drms_submit: Request submitting into the DRMS.
    :type drms_submit: bool | None
    :param irf_path: An :emphasis:`Illumina Run Folder` path.
    :type irf_path: str | None
    :param library_path: A library path.
    :type library_path: str | None
    :param mode: A sequencer mode.
    :type mode: str | None
    :param illumina2bam: Request processing with Illumina2bam rather than Picard tools.
    :type illumina2bam: bool | None
    :param no_validation: Request skipping library annotation sheet validation.
    :type no_validation: bool | None
    :param compression_level: A Zlib compression level.
    :type compression_level: int | None
    :param loop: Request looping until :literal:`RTAComplete.txt` is available.
    :type loop: bool | None
    :param interval: A looping interval in seconds.
    :type interval: float | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if configuration_path == Configuration.get_global_file_path():
        if not irf_path:
            raise Exception('The --irf argument is required if configuration is not set.')

    # if name_space.archive_directory:
    #
    #     # Extract the Illumina Run Folder first.
    #
    #     irf_restore: IlluminaRunFolderRestore = IlluminaRunFolderRestore.from_config_file_path(
    #         config_path=name_space.configuration)
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

    if illumina2bam:
        analysis_itb: IlluminaToBam = IlluminaToBam.from_config_file_path(config_path=configuration_path)

        if irf_path:
            analysis_itb.run_directory = irf_path

        # if irf_restore is not None:
        #     # If the IlluminaRunFolderRestore has not run, the run folder is not complete.
        #     itb.run_directory = irf_restore.get_run_directory_path
        #     itb.force = True

        if loop:
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
                time.sleep(interval)
        else:
            analysis_itb.run()

        analysis_itb.check_state()
        analysis_itb.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis_itb.name)
        print('Project name:           ', analysis_itb.project_name)
        print('Project directory:      ', analysis_itb.project_directory)
        print('Illumina run directory: ', analysis_itb.run_directory)
        print('Experiment directory:   ', analysis_itb.get_experiment_directory)

        analysis_bid: BamIndexDecoder = BamIndexDecoder.from_config_file_path(config_path=configuration_path)

        # Transfer the project name from the IlluminaToBam to the BamIndexDecoder analysis.

        analysis_bid.project_name = analysis_itb.project_name

        if mode:
            if mode == 'high':
                analysis_bid.lanes = 8
            elif mode == 'rapid':
                analysis_bid.lanes = 2
            elif mode == 'miseq':
                analysis_bid.lanes = 1
            elif mode == 'nextseq':
                analysis_bid.lanes = 4
            elif mode == 'novaseq':
                analysis_bid.lanes = 4
            else:
                raise Exception(f'The --mode option {mode!r} is not supported.')
        else:
            analysis_bid.lanes = RunFolder.from_file_path(
                file_path=analysis_itb.run_directory).run_information.flow_cell_layout.lane_count

        if no_validation:
            analysis_bid.force = no_validation

        if library_path:
            analysis_bid.library_path = library_path

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
            analysis_bid.submit(name=stage_name, drms_submit=drms_submit)

            print('')
            print(analysis_bid.name)
            print('Project name:         ', analysis_bid.project_name)
            print('Project directory:    ', analysis_bid.project_directory)
            print('Sequences directory:  ', analysis_bid.sequences_directory)
            print('Experiment directory: ', analysis_bid.get_experiment_directory)
    else:
        analysis_ims: IlluminaMultiplexSam = IlluminaMultiplexSam.from_config_file_path(config_path=configuration_path)

        if irf_path:
            analysis_ims.run_directory = irf_path

        if compression_level is not None:
            analysis_ims.compression_level = compression_level

        if loop:
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
                time.sleep(interval)
        else:
            analysis_ims.run()

        analysis_ims.check_state()
        analysis_ims.submit(name=stage_name, drms_submit=drms_submit)

        print(analysis_ims.name)
        print('Project name:           ', analysis_ims.project_name)
        print('Project directory:      ', analysis_ims.project_directory)
        print('Illumina run directory: ', analysis_ims.run_directory)
        print('Experiment directory:   ', analysis_ims.get_experiment_directory)

        analysis_ids: IlluminaDemultiplexSam = IlluminaDemultiplexSam.from_config_file_path(
            config_path=configuration_path)

        # Transfer the project_name and run_directory from the IlluminaMultiplexSam to the
        # IlluminaDemultiplexSam analysis.

        analysis_ids.project_name = analysis_ims.project_name
        analysis_ids.run_directory = analysis_ims.run_directory

        if compression_level is not None:
            analysis_ids.compression_level = compression_level

        if mode:
            if mode == 'high':
                analysis_ids.lanes = 8
            elif mode == 'rapid':
                analysis_ids.lanes = 2
            elif mode == 'miseq':
                analysis_ids.lanes = 1
            elif mode == 'nextseq':
                analysis_ids.lanes = 4
            elif mode == 'novaseq':
                analysis_ids.lanes = 4
            else:
                raise Exception(f'The --mode option {mode!r} is not supported.')
        else:
            analysis_ids.lanes = RunFolder.from_file_path(
                file_path=analysis_ims.run_directory).run_information.flow_cell_layout.lane_count

        if no_validation:
            analysis_ids.force = no_validation

        if library_path:
            analysis_ids.library_path = library_path

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
            analysis_ids.submit(name=stage_name, drms_submit=drms_submit)

            print('')
            print(analysis_ids.name)
            print('Project name:         ', analysis_ids.project_name)
            print('Project directory:    ', analysis_ids.project_directory)
            print('Sequences directory:  ', analysis_ids.sequences_directory)
            print('Experiment directory: ', analysis_ids.get_experiment_directory)

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Illumina Run Folder processor driver script.')

    argument_parser.add_argument(
        '--dry-run',
        action='store_false',
        help='dry run',
        dest='drms_submit')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--stage-name',
        help='limit job submission to a particular analysis stage')

    argument_parser.add_argument(
        '--illumina2bam',
        action='store_true',
        help='Use Illumina2bam rather than Picard tools')

    # argument_parser.add_argument(
    #     '--archive-directory',
    #     help='archive directory')
    #
    # argument_parser.add_argument(
    #     '--project-name',
    #     help='project name (i.e., flow cell identifier)')
    #
    # argument_parser.add_argument(
    #     '--extract-intensities',
    #     action='store_true',
    #     help='extract cluster intensity (*.cif) files')
    #
    # argument_parser.add_argument(
    #     '--force',
    #     action='store_true',
    #     help='force expanding a run folder even if it exists already')

    argument_parser.add_argument(
        '--irf',
        help='Illumina Run Folder name or file path',
        dest='irf_path')

    argument_parser.add_argument(
        '--library-path',
        help='library annotation sheet (*.csv) file path')

    argument_parser.add_argument(
        '--mode',
        help='HiSeq run mode (i.e., "high" (high-output) or "rapid" (rapid run) or "miseq" for a MiSeq run)')

    argument_parser.add_argument(
        '--compression',
        type=int,
        help='Picard (Zlib) compression level [9]')

    argument_parser.add_argument(
        '--no-validation',
        action='store_true',
        help='force processing even if library annotation sheet validation fails')

    argument_parser.add_argument(
        '--loop',
        action='store_true',
        help='loop until a RTAComplete.txt file has been copied by the Illumina Real-Time Analysis (RTA) software')

    argument_parser.add_argument(
        '--interval',
        type=float,
        default=120.0,
        help='loop interval in seconds [120.0]')

    argument_parser.add_argument(
        'configuration',
        nargs='?',
        default=Configuration.get_global_file_path(),
        help=f'configuration (*.ini) file path [{Configuration.get_global_file_path()!s}]')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return run(
        configuration_path=name_space.configuration,
        stage_name=name_space.stage_name,
        drms_submit=name_space.drms_submit,
        irf_path=name_space.irf_path,
        library_path=name_space.library_path,
        mode=name_space.mode,
        illumina2bam=name_space.illumina2bam,
        no_validation=name_space.no_validation,
        compression_level=name_space.compression,
        loop=name_space.loop,
        interval=name_space.interval)


if __name__ == '__main__':
    sys.exit(main())
