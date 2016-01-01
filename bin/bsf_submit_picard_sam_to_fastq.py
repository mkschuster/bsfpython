#! /usr/bin/env python
#
# BSF Python script to drive the Picard SamToFastq analysis.
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

import argparse

from bsf import Default
from bsf.analyses.picard import SamToFastq


argument_parser = argparse.ArgumentParser(
    description='Picard SamToFastq analysis driver script.')

argument_parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--stage',
    dest='stage',
    help='limit job submission to a particular Analysis stage',
    required=False)

argument_parser.add_argument(
    '--configuration',
    default=Default.global_file_path,
    help='configuration (*.ini) file path',
    required=False,
    type=str)

argument_parser.add_argument(
    '--project-name',
    dest='project_name',
    help='project name',
    required=True,
    type=str)

argument_parser.add_argument(
    '--sas-file',
    dest='sas_file',
    help='sample annotation sheet (*.csv) file path',
    required=True,
    type=unicode)

name_space = argument_parser.parse_args()

# Create a Picard SamToFastq analysis and run it.

stf = SamToFastq.from_config_file_path(config_path=name_space.configuration)

if name_space.debug:
    assert isinstance(name_space.debug, int)
    stf.debug = name_space.debug

if name_space.project_name:
    assert isinstance(name_space.project_name, str)
    stf.project_name = name_space.project_name

if name_space.sas_file:
    assert isinstance(name_space.sas_file, (str, unicode))
    stf.sas_file = name_space.sas_file

annotation_sheet = stf.run()
stf.submit(drms_name=name_space.stage)

print 'Picard SamToFastq Analysis'
print 'Project name:      ', stf.project_name
print 'Genome version:    ', stf.genome_version
print 'Input directory:   ', stf.input_directory
print 'Output directory:  ', stf.output_directory
print 'Project directory: ', stf.project_directory
print 'Genome directory:  ', stf.genome_directory

if stf.debug >= 2:
    print '{!r} final trace:'.format(stf)
    print stf.trace(level=1)
