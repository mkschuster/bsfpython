#! /usr/bin/env python
#
# BSF Python script to process an Illumina Run Folder (IRF) after sequencing and
# drive the IlluminaToBamTools IlluminaToBam and BamIndexDecoder analyses.
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
import datetime
import os
import string
import time

from bsf import Default
from bsf.analyses.illumina_to_bam_tools import BamIndexDecoder, IlluminaToBam, IlluminaRunFolderNotComplete


argument_parser = ArgumentParser(
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
    default=Default.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--mode',
    help='HiSeq run mode i.e. high (high-output) or rapid (rapid run)',
    required=False,
    type=str)

argument_parser.add_argument(
    '--loop',
    action='store_true',
    help='loop until a RTAComplete.txt file has been copied by the Illumina Real-Time Analysis (RTA) software')

argument_parser.add_argument(
    '--interval',
    help='loop interval [s]',
    default=120,
    type=int)

name_space = argument_parser.parse_args()

# Create a BSF IlluminaToBam analysis, run and submit it.

itb = IlluminaToBam.from_config_file_path(config_path=name_space.configuration)

# Set arguments that override the configuration file.

if name_space.debug:
    itb.debug = name_space.debug

if name_space.irf:
    itb.illumina_run_folder = name_space.irf

# Do the work.

if name_space.loop:
    # If the --loop option has been set, wait until RTAComplete.txt has been copied.
    loop_counter = int(1)
    while 1:
        print '[{}] Loop {}:'.format(datetime.datetime.now().isoformat(), loop_counter)
        loop_counter += int(1)
        try:
            itb.run()
        except IlluminaRunFolderNotComplete as exception:
            print exception
        else:
            print 'Illumina Run Folder seems complete.'
            break
        time.sleep(name_space.interval)
else:
    itb.run()

itb.submit(drms_name=name_space.stage)

print 'IlluminaToBamTools IlluminaToBam Analysis'
print 'Project name:         ', itb.project_name
print 'Project directory:    ', itb.project_directory
print 'Illumina Run Folder:  ', itb.illumina_run_folder
print 'Experiment directory: ', itb.experiment_directory

# Create a BSF BamIndexDecoder analysis, run and submit it.

bid = BamIndexDecoder.from_config_file_path(config_path=name_space.configuration)

# Transfer the project name from the IlluminaToBam to the BamIndexDecoder analysis.

bid.project_name = itb.project_name

# Set arguments that override the configuration file.

if name_space.debug:
    bid.debug = name_space.debug

if name_space.mode:
    if name_space.mode == 'high':
        bid.lanes = int(8)
    elif name_space.mode == 'rapid':
        bid.lanes = int(2)
    else:
        raise Exception("Unknown output mode " + name_space.mode)

if name_space.library_path:
    bid.library_path = name_space.library_path

# If a library file has not been defined so far, check,
# if a standard library file i.e. PROJECT_NAME_libraries.csv exists in the current directory.

if not bid.library_path:
    library_path = string.join(words=(bid.project_name, 'libraries.csv'), sep='_')
    if os.path.exists(path=library_path):
        bid.library_path = library_path

if bid.library_path:

    # Do the work if, at this stage, a library file has been set.

    bid.run()
    bid.submit(drms_name=name_space.stage)

    print ''
    print 'IlluminaToBamTools BamIndexDecoder Analysis'
    print 'Project name:         ', bid.project_name
    print 'Project directory:    ', bid.project_directory
    print 'Sequences directory:  ', bid.sequences_directory
    print 'Experiment directory: ', bid.experiment_directory
