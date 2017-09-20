#! /usr/bin/env python
#
# BSF Python script to drive the Picard CollectHiSeqXPfFailMetrics analysis.
#
#
# Copyright 2013 - 2017 Michael K. Schuster
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

import argparse

from bsf.analyses.picard import CollectHiSeqXPfFailMetrics
from bsf.standards import Default

argument_parser = argparse.ArgumentParser(
    description='Picard CollectHiSeqXPfFailMetrics Analysis driver script.')

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
    '--irf',
    help='Illumina Run Folder name or file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--force',
    action='store_true',
    help='force processing of an incomplete Illumina Run Folder')

name_space = argument_parser.parse_args()

# This analysis requires either a non-default --configuration argument or a
# --irf argument.

if name_space.configuration == Default.global_file_path:
    if name_space.irf is None:
        raise Exception("argument --irf is required if --configuration is not set")

# Create a CollectHiSeqXPfFailMetrics analysis, run and submit it.

analysis = CollectHiSeqXPfFailMetrics.from_config_file_path(config_path=name_space.configuration)
""" @type analysis: bsf.analyses.picard.CollectHiSeqXPfFailMetrics """

# Set arguments that override the configuration file.

if name_space.debug:
    assert isinstance(name_space.debug, int)
    analysis.debug = name_space.debug

if name_space.irf:
    assert isinstance(name_space.irf, (str, unicode))
    analysis.run_directory = name_space.irf

if name_space.force:
    assert isinstance(name_space.force, bool)
    analysis.force = name_space.force

# Do the work.

analysis.run()
analysis.check_state()
analysis.submit(name=name_space.stage)

print 'Picard CollectHiSeqXPfFailMetrics Analysis'
print 'Project name:           ', analysis.project_name
print 'Project directory:      ', analysis.project_directory
print 'Illumina run directory: ', analysis.run_directory

if analysis.debug >= 2:
    print '{!r} final trace:'.format(analysis)
    print analysis.trace(level=1)
