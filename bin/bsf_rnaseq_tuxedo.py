#! /usr/bin/env python
#
# BSF Python script to drive the Tuxedo suite-based RNA-Seq analysis pipeline.
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

from Bio.BSF.Analyses.RNASeq import Tuxedo


parser = argparse.ArgumentParser(description='RNA-Seq analysis driver script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('--stage', dest='stage', required=False,
                    help='Limit job submission to a particular Analysis stage')

parser.add_argument('configuration',
                    help='Configuration file (*.ini)')

args = parser.parse_args()

# Create a BSF ChIPSeq analysis and run it.

tuxedo = Tuxedo.from_config_file(config_file=args.configuration)

if args.debug:
    tuxedo.debug = args.debug

tuxedo.run()

# Submit all Executable objects of all Distributed Resource Management System objects.

submit = 0

for drms in tuxedo.drms_list:

    if args.stage:
        if args.stage == drms.name:
            submit += 1
        else:
            continue

    drms.submit(debug=tuxedo.debug)

    if tuxedo.debug:
        print repr(drms)
        print drms.trace(1)

if args.stage:
    if args.stage == 'report':
        tuxedo.report()
        pass
    elif not submit:
        name_list = [drms.name for drms in tuxedo.drms_list]
        name_list.append('report')
        print 'Valid Analysis stages are: {!r}'.format(name_list)

print 'RNA-Seq Analysis'
print 'Project name:      ', tuxedo.project_name
print 'Genome version:    ', tuxedo.genome_version
print 'Input directory:   ', tuxedo.input_directory
print 'Output directory:  ', tuxedo.output_directory
print 'Project directory: ', tuxedo.project_directory
print 'Genome directory:  ', tuxedo.genome_directory
