#! /usr/bin/env python
#
# BSF Python script to run FastQC over an Illumina run folder.
#
#
# Copyright 2013 - 2018 Michael K. Schuster
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

from __future__ import print_function

import argparse

from bsf.analyses import RunFastQC

parser = argparse.ArgumentParser(description='FastQC driver script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('--stage', required=False,
                    help='Limit job submission to a particular Analysis stage')

parser.add_argument('configuration',
                    help='Configuration file (*.ini)')

args = parser.parse_args()

fastqc = RunFastQC.from_config_file_path(config_path=args.configuration)

if args.debug:
    fastqc.debug = args.debug

fastqc.run()

# Submit all Executable objects of all Stage objects.

submit = 0

for stage in fastqc.stage_list:

    if args.stage:
        if args.stage == stage.name:
            submit += 1
        else:
            continue

    stage.submit(debug=fastqc.debug)

    if fastqc.debug:
        print(repr(stage))
        print(stage.trace(1))

if args.stage:
    if args.stage == 'report':
        fastqc.report()
    elif not submit:
        name_list = [stage.name for stage in fastqc.stage_list]
        name_list.append('report')
        print('Valid Analysis stages are:', repr(name_list))

print(fastqc.name)
print('Project name:      ', fastqc.project_name)
print('Input directory:   ', fastqc.input_directory)
print('Output directory:  ', fastqc.output_directory)
print('Project directory: ', fastqc.project_directory)
