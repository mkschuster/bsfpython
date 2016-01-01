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

from argparse import ArgumentParser

from bsf import Executable, Runnable


argument_parser = ArgumentParser(
    description='Trace a pickled BSF Runnable.')

argument_parser.add_argument(
    '--pickler-path',
    dest='pickler_path',
    help='file path to a pickled Runnable object',
    required=True)

arguments = argument_parser.parse_args()

runnable = Runnable.from_pickler_path(file_path=arguments.pickler_path)

print runnable.trace(level=1)

keys = runnable.executable_dict.keys()
keys.sort()

for key in keys:
    executable = runnable.executable_dict[key]
    assert isinstance(executable, Executable)
    print 'Executable key:     {}'.format(key)
    print 'Executable command: {}'.format(executable.command_str())
