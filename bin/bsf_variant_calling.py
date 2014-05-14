#! /usr/bin/env python
#
# BSF Python script to drive the Variant Calling analysis pipeline.
#
#
# Copyright 2013 Michael K. Schuster
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

from Bio.BSF.Analyses.VariantCalling import VariantCallingGATK


parser = argparse.ArgumentParser(description='Variant Calling analysis driver script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('--stage', required=False,
                    help='Limit job submission to a particular Analysis stage')

parser.add_argument('configuration',
                    help='Configuration file (*.ini)')

args = parser.parse_args()

# Create a BSF Variant Calling analysis and run it.

variant_calling = VariantCallingGATK.from_config_file(config_file=args.configuration)

if args.debug:
    variant_calling.debug = args.debug

variant_calling.run()

# Submit all Executable objects of all Distributed Resource Management System objects.

submit = 0

for drms in variant_calling.drms_list:

    if args.stage:
        if args.stage == drms.name:
            submit += 1
        else:
            continue

    drms.submit(debug=variant_calling.debug)

    if variant_calling.debug:
        print repr(drms)
        print drms.trace(1)

if args.stage:
    if args.stage == 'report':
        variant_calling.report()
    elif not submit:
        name_list = [drms.name for drms in variant_calling.drms_list]
        name_list.append('report')
        print 'Valid Analysis stages are: {!r}'.format(name_list)

print 'Variant Calling Analysis'
print 'Project name:      ', variant_calling.project_name
print 'Genome version:    ', variant_calling.genome_version
print 'Input directory:   ', variant_calling.input_directory
print 'Output directory:  ', variant_calling.output_directory
print 'Project directory: ', variant_calling.project_directory
print 'Genome directory:  ', variant_calling.genome_directory
