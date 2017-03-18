#! /usr/bin/env python
#
# BSF Python script to drive the STAR Aligner analysis.
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

from bsf.analyses.star_aligner import StarAligner
from bsf.standards import Default


argument_parser = argparse.ArgumentParser(
    description='STAR aligner analysis driver script.')

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
    required=False,
    type=str)

argument_parser.add_argument(
    '--sas-file',
    dest='sas_file',
    help='sample annotation sheet (*.csv) file path',
    required=False,
    type=unicode)

name_space = argument_parser.parse_args()

# This analysis requires either a non-default --configuration argument or a
# --project-name and --sas-file argument.

if name_space.configuration == Default.global_file_path:
    if name_space.project_name is None:
        raise Exception("argument --project-name is required if --configuration is not set")
    if name_space.sas_file is None:
        raise Exception("argument --sas-file is required if --configuration is not set")

# Create a STAR Aligner analysis and run it.

star_aligner = StarAligner.from_config_file_path(config_path=name_space.configuration)

if name_space.debug:
    assert isinstance(name_space.debug, int)
    star_aligner.debug = name_space.debug

if name_space.project_name:
    assert isinstance(name_space.project_name, str)
    star_aligner.project_name = name_space.project_name

if name_space.sas_file:
    assert isinstance(name_space.sas_file, basestring)
    star_aligner.sas_file = name_space.sas_file

annotation_sheet = star_aligner.run()
star_aligner.check_state()
star_aligner.submit(name=name_space.stage)

print 'STAR Aligner Analysis'
print 'Project name:      ', star_aligner.project_name
print 'Genome version:    ', star_aligner.genome_version
print 'Input directory:   ', star_aligner.input_directory
print 'Output directory:  ', star_aligner.output_directory
print 'Project directory: ', star_aligner.project_directory
print 'Genome directory:  ', star_aligner.genome_directory

if star_aligner.debug >= 2:
    print '{!r} final trace:'.format(star_aligner)
    print star_aligner.trace(level=1)
