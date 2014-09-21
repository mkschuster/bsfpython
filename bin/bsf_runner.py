#! /usr/bin/env python
#
# BSF Python wrapper script that runs BSF Runnable modules.
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

from argparse import ArgumentParser
import importlib

from Bio.BSF import Analysis


argument_parser = ArgumentParser(description='Generic BSF runner script.')

argument_parser.add_argument('--debug', required=False, type=int,
                             help='Debug level')

argument_parser.add_argument('--runnable_name', required=True,
                             help='BSF Runnable name')

argument_parser.add_argument('--pickler_path', required=True,
                             help='File path to the pickled Analysis object.')

arguments = argument_parser.parse_args()

analysis = Analysis.from_picker_file(file_path=arguments.pickler_path)

if not arguments.runnable_name in analysis.runnable_dict:
    raise Exception(
        "A Runnable with name {!r} has not been defined in Analysis project name {!r}.".
        format(arguments.runnable_name), analysis.project_name)

runnable = analysis.runnable_dict[arguments.runnable_name]

module = importlib.import_module(name=runnable.code_module)

module.run(runnable=runnable)
