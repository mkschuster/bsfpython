#! /usr/bin/env python
#
# BSF Python script to submit Illumina Run Folder archive processes.
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


import argparse
import os

from bsf import Default, DRMS, Executable


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Get global defaults.

default = Default.get_global_default()

# Parse the arguments.

parser = argparse.ArgumentParser(
    description='BSF Runner for post-processing alignments prior to GATK variant calling.')

parser.add_argument('--debug', required=False, type=int,
                    help='debug level')

parser.add_argument('--force', required=False,
                    help='force archiving')

parser.add_argument('irf_path',
                    help='file path to an Illumina Run Folder')

args = parser.parse_args()

irf_path = os.path.abspath(path=args.irf_path)

if not os.path.exists(path=irf_path):
    raise Exception('Could not find Illumina Run Folder: {!r}'.format(irf_path))

drms = DRMS(name='irf_archive', work_directory=os.path.dirname(irf_path))
drms.set_default(default=default)

executable = Executable(name='irf_archive_' + os.path.basename(irf_path),
                        program='bsf_run_irf_archive.bash')
drms.add_executable(executable=executable)

executable.arguments.append(os.path.basename(irf_path))

print "Resulting DRMS:"
print drms.trace(level=1)

# TODO: Submit this Executable.
# TODO: This needs some basic configuration for all Illumina Run Folders that are archived.
# Maybe we should have a concept of a configuration directory from which those could be read?
# The .bsfpython.ini file may not be specific enough?
# What about including default configurations for bsf.analyses modules in the .bsfpython.ini file? Implemented.
