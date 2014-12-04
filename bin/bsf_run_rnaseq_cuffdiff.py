#! /usr/bin/env python
#
# BSF Python script to run Cuffdiff.
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
import errno
import os
from pickle import Unpickler
import shutil

from bsf import Default, Executable, Runnable


# Set the environment consistently.

os.environ['LANG'] = 'C'

# Get global defaults.

default = Default.get_global_default()

# Parse the arguments.

parser = argparse.ArgumentParser(
    description='BSF Runner for the Tuxedo Cuffdiff application.')

parser.add_argument('--debug', required=False, type=int,
                    help='debug level')

parser.add_argument('--pickler_path', required=True,
                    help='file path to a Python Pickler file.')

args = parser.parse_args()

# Unpickle the file into a Python dict object.

pickler_file = open(args.pickler_path, 'rb')

unpickler = Unpickler(file=pickler_file)

pickler_dict = unpickler.load()

pickler_file.close()

key = 'prefix'
if key in pickler_dict and pickler_dict[key]:
    prefix = pickler_dict[key]
    if prefix[-1:] != '_':
        prefix += '_'
else:
    prefix = str()

key = 'replicate_key'
if key in pickler_dict and pickler_dict[key]:
    replicate_key = pickler_dict[key]
else:
    raise Exception("The pickler dict in the pickler file needs to contain key {!r}.".format(key))

# Create a temporary directory.

path_temporary = "{}{}_temporary".format(prefix, replicate_key)

if not os.path.isdir(path_temporary):
    try:
        os.makedirs(path_temporary)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

run_cuffdiff = pickler_dict['cuffdiff_executable']
assert isinstance(run_cuffdiff, Executable)

child_return_code = Runnable.run(executable=run_cuffdiff)

if child_return_code:
    raise Exception('Could not complete the Cuffdiff step.')

# Remove the temporary directory and everything within it.

shutil.rmtree(path=path_temporary, ignore_errors=False)

# Job done.
