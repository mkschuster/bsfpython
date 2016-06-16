#! /usr/bin/env python
#
# BSF Python script to drive the ChIP-Seq analysis pipeline.
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

from bsf.analyses import ChIPSeq


parser = argparse.ArgumentParser(description='ChIP-Seq analysis driver script.')

parser.add_argument('--debug', required=False, type=int,
                    help='Debug level')

parser.add_argument('--stage', required=False,
                    help='Limit job submission to a particular Analysis stage')

parser.add_argument('configuration',
                    help='Configuration file (*.ini)')

args = parser.parse_args()

# Create a BSF ChIPSeq analysis and run it.

chipseq = ChIPSeq.from_config_file_path(config_path=args.configuration)

if args.debug:
    chipseq.debug = args.debug

chipseq.run()
chipseq.check_state()

# Submit all Executable objects of all Distributed Resource Management System objects.

submit = 0

for drms in chipseq.drms_list:

    if args.stage:
        if args.stage == drms.name:
            submit += 1
        else:
            continue

    drms.submit(debug=chipseq.debug)

    if chipseq.debug:
        print repr(drms)
        print drms.trace(1)

if args.stage:
    if args.stage == 'report':
        chipseq.report()
    elif not submit:
        name_list = [drms.name for drms in chipseq.drms_list]
        name_list.append('report')
        print 'Valid Analysis stages are: {!r}'.format(name_list)

print 'ChIP-Seq Analysis'
print 'Project name:      ', chipseq.project_name
print 'Genome version:    ', chipseq.genome_version
print 'Input directory:   ', chipseq.input_directory
print 'Output directory:  ', chipseq.output_directory
print 'Project directory: ', chipseq.project_directory
print 'Genome directory:  ', chipseq.genome_directory
