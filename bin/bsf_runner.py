#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
#
# BSF Python wrapper script that runs a Runnable module.
# The script un-pickles the Runnable object from a file path,
# loads the specific processing logic from a runnables code module
# and finally, calls Runnable.run(). Exit codes are defined in the
# accessory code module.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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
#
import importlib
from argparse import ArgumentParser

from bsf.procedure import Runnable

argument_parser = ArgumentParser(
    description='Generic BSF runner script.')

argument_parser.add_argument(
    '--pickler-path',
    dest='pickler_path',
    help='file path to a pickled Runnable object',
    required=True)

name_space = argument_parser.parse_args()

runnable = Runnable.from_pickler_path(file_path=name_space.pickler_path)

module_type = importlib.import_module(name=runnable.code_module)

run_function = getattr(module_type, 'run')
run_function(runnable=runnable)
