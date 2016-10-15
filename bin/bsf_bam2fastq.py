#! /usr/bin/env python
#
# BSF Python driver script for the conversion of BAM files into FASTQ format.
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

from bsf.analyses import RunBamToFastq


parser = argparse.ArgumentParser(description='BAM To FASTQ analysis driver script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('--stage', required=False,
                    help='Limit job submission to a particular Analysis stage')

parser.add_argument('configuration',
                    help='Configuration file (*.ini)')

args = parser.parse_args()

# Create a RunBamToFastq BSF Analysis and run it.

bam2fastq = RunBamToFastq.from_config_file_path(config_path=args.configuration)

if args.debug:
    bam2fastq.debug = args.debug

bam2fastq.run()

# Submit all Executable objects of all Stage objects.

submit = 0

for stage in bam2fastq.stage_list:

    if args.stage:
        if args.stage == stage.name:
            submit += 1
        else:
            continue

    stage.submit(debug=bam2fastq.debug)

    if bam2fastq.debug:
        print repr(stage)
        print stage.trace(1)

if args.stage:
    if args.stage == 'report':
        pass
    elif not submit:
        name_list = [stage.name for stage in bam2fastq.stage_list]
        name_list.append('report')
        print 'Valid Analysis stages are: {!r}'.format(name_list)

print 'RunBamToFastq Analysis'
print 'Project name:      ', bam2fastq.project_name
print 'Genome version:    ', bam2fastq.genome_version
print 'Input directory:   ', bam2fastq.input_directory
print 'Output directory:  ', bam2fastq.output_directory
print 'Project directory: ', bam2fastq.project_directory
print 'Genome directory:  ', bam2fastq.genome_directory
