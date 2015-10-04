#! /usr/bin/env python
#
# BSF Python utility script to append a UUID to a string.
#
#
# Copyright 2013 - 2015 Michael K. Schuster
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

"""
Generate UUIDs for uses such as keeping public sequencing project URLs
just a bit more private.
"""

from argparse import ArgumentParser
import uuid

argument_parser = ArgumentParser(
    description='Append a UUID to a string.')

argument_parser.add_argument(
    'input',
    help="input string, can be empty by supplying ''",
    type=str)

name_space = argument_parser.parse_args()

# By passing in '' as an empty string, a leading '_' needs stripping.

uuid_string = '_'.join((name_space.input, uuid.uuid4().hex)).lstrip('_')

print "UUID: {}".format(uuid_string)
