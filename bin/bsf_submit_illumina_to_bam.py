#! /usr/bin/env python
#
# BSF Python script to drive the IlluminaToBamTools IlluminaToBam analysis.
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

from Bio.BSF.Analyses.IlluminaToBamTools import IlluminaToBam


argument_parser = ArgumentParser(
    description='IlluminaToBamTools Illumina2bam analysis driver script.')

argument_parser.add_argument('--debug', required=False, type=int,
                             help='Debug level')

argument_parser.add_argument('--stage', required=False, type=str,
                             help='Limit job submission to a particular Analysis stage')

argument_parser.add_argument('configuration', type=str,
                             help='Configuration file (*.ini)')

arguments = argument_parser.parse_args()

# Create a BSF IlluminaToBam analysis, run and submit it.

itb = IlluminaToBam.from_config_file(config_file=arguments.configuration)

# Set arguments that override the configuration file.

if arguments.debug:
    itb.debug = arguments.debug

# Do the work.

itb.run()
itb.submit(drms_name=arguments.stage)

print 'IlluminaToBamTools IlluminaToBam Analysis'
print 'Project name:      ', itb.project_name
print 'Input directory:   ', itb.input_directory
print 'Project directory: ', itb.project_directory
