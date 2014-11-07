#! /usr/bin/env python
#
# BSF Python script to drive the ChIPSeq analysis.
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

from Bio.BSF.Analyses.ChIPSeq import ChIPSeq


argument_parser = ArgumentParser(
    description='ChIPSeq analysis driver script.')

argument_parser.add_argument(
    '--debug',
    help='Debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--stage',
    help='Limit job submission to a particular Analysis stage',
    required=False,
    type=str)

argument_parser.add_argument(
    'configuration',
    help='Configuration file (*.ini)',
    type=str)

arguments = argument_parser.parse_args()

# Create a BSF BamIndexDecoder analysis, run and submit it.

chipseq = ChIPSeq.from_config_file(config_file=arguments.configuration)

# Set arguments that override the configuration file.

if arguments.debug:
    chipseq.debug = arguments.debug

# Do the work.

chipseq.run()
chipseq.submit(drms_name=arguments.stage)

print 'ChIPSeq Analysis'
print 'Project name:      ', chipseq.project_name
print 'Input directory:   ', chipseq.input_directory
print 'Project directory: ', chipseq.project_directory
