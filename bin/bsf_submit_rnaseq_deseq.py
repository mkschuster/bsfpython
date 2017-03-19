#! /usr/bin/env python
#
# BSF Python script to drive the DESeq-based RNA-Seq analysis pipeline.
#
#
# Copyright 2013 - 2017 Michael K. Schuster
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

from bsf.analyses.rna_seq import DESeq


argument_parser = argparse.ArgumentParser(
    description='RNA-Seq analysis driver script.')

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
    'configuration',
    help='configuration (*.ini) file path')

name_space = argument_parser.parse_args()

# Create a DESeq analysis and run it.

deseq = DESeq.from_config_file_path(config_path=name_space.configuration)

if name_space.debug:
    deseq.debug = name_space.debug

deseq.run()
deseq.check_state()
deseq.submit(name=name_space.stage)

print 'RNA-Seq Analysis'
print 'Project name:      ', deseq.project_name
print 'Genome version:    ', deseq.genome_version
print 'Input directory:   ', deseq.input_directory
print 'Output directory:  ', deseq.output_directory
print 'Project directory: ', deseq.project_directory
print 'Genome directory:  ', deseq.genome_directory

if deseq.debug >= 2:
    print '{!r} final trace:'.format(deseq)
    print deseq.trace(level=1)
