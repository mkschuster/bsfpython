#! /usr/bin/env python
#
# BSF Python script to drive the Variant Calling analysis pipeline.
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

from bsf.analyses.variant_calling import VariantCallingGATK


argument_parser = argparse.ArgumentParser(
    description='Variant Calling analysis driver script.')

argument_parser.add_argument(
    '--debug',
    required=False,
    type=int,
    help='debug level')

argument_parser.add_argument(
    '--stage',
    required=False,
    help='limit job submission to a particular Analysis stage')

argument_parser.add_argument(
    'configuration',
    help='configuration (*.ini) file path')

name_space = argument_parser.parse_args()

# Create a BSF Variant Calling analysis and run it.

variant_calling = VariantCallingGATK.from_config_file_path(config_path=name_space.configuration)

if name_space.debug:
    variant_calling.debug = name_space.debug

variant_calling.run()
variant_calling.submit(drms_name=name_space.stage)

print 'Variant Calling Analysis'
print 'Project name:      ', variant_calling.project_name
print 'Genome version:    ', variant_calling.genome_version
print 'Input directory:   ', variant_calling.input_directory
print 'Output directory:  ', variant_calling.output_directory
print 'Project directory: ', variant_calling.project_directory
print 'Genome directory:  ', variant_calling.genome_directory
