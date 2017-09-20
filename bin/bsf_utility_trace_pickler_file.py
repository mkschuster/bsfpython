#! /usr/bin/env python
#
# BSF Python wrapper script that traces a Runnable.
# The script un-pickles the Runnable object from a file path
# and calls Runnable.trace().
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

from bsf import Runnable

argument_parser = argparse.ArgumentParser(
    description='Trace a pickled Runnable.')

argument_parser.add_argument(
    '--format',
    choices=['list', 'str'],
    default='str',
    help='Output format Python str or Python list [str]',
    required=False)

argument_parser.add_argument(
    'pickler_path',
    help='file path to a pickled Runnable object')

name_space = argument_parser.parse_args()

runnable = Runnable.from_pickler_path(file_path=name_space.pickler_path)

print runnable.trace(level=1)

for runnable_step in runnable.runnable_step_list:
    if name_space.format == 'list':
        print 'RunnableStep command list:', runnable_step.command_list()
        print
    elif name_space.format == 'str':
        print 'RunnableStep command str:', runnable_step.command_str()
        print
