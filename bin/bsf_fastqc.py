#! /usr/bin/env python
#
# BSF Python script to run FastQC over an Illumina run folder.
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

from Bio.BSF.Analyses import RunFastQC


parser = argparse.ArgumentParser(description='FastQC driver script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('--stage', required=False,
                    help='Limit job submission to a particular Analysis stage')

parser.add_argument('configuration',
                    help='Configuration file (*.ini)')

args = parser.parse_args()

fastqc = RunFastQC.from_config_file(config_file=args.configuration)

if args.debug:
    fastqc.debug = args.debug

fastqc.run()

# Submit all Executable objects of all Distributed Resource Management System objects.

submit = 0

for drms in fastqc.drms_list:

    if args.stage:
        if args.stage == drms.name:
            submit += 1
        else:
            continue

    drms.submit(debug=fastqc.debug)

    if fastqc.debug:
        print repr(drms)
        print drms.trace(1)

if args.stage:
    if args.stage == 'report':
        fastqc.report()
    elif not submit:
        name_list = [drms.name for drms in fastqc.drms_list]
        name_list.append('report')
        print 'Valid Analysis stages are: {!r}'.format(name_list)

print 'FastQC Analysis'
print 'Project name:      ', fastqc.project_name
print 'Input directory:   ', fastqc.input_directory
print 'Output directory:  ', fastqc.output_directory
print 'Project directory: ', fastqc.project_directory
