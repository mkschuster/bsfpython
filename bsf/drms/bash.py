# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
#
#  Biomedical Sequencing Facility (BSF), part of the genomics core facility
#  of the Research Center for Molecular Medicine (CeMM) of the
#  Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
#  This file is part of BSF Python.
#
#  BSF Python is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BSF Python is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.
#
"""GNU Bourne-Again Shell (Bash) DRMS module.

A package of methods supporting the GNU Bourne-Again Shell (Bash) as
Distributed Resource Management System (DRMS) module.
"""
import os
from typing import List

from bsf.connector import ConnectorFile


def submit(stage, debug=0):
    """Submit each C{bsf.process.Executable} of a C{bsf.analysis.Stage}.

    Submits each C{bsf.process.Executable} by writing a GNU Bourne-Again Shell (BASH) script
    into the C{bsf.analysis.Stage.work_directory}.

    @param stage: C{bsf.analysis.Stage}
    @type stage: bsf.analysis.Stage
    @param debug: Debug level
    @type debug: int
    """

    output_list: List[str] = list()

    output_list.append('#!/usr/bin/env bash\n')
    output_list.append('\n')

    if debug > 0:
        output_list.append('# BSF-Python debug mode: ' + repr(debug) + '\n')
        output_list.append('\n')

    for executable in stage.executable_list:
        if not executable.submit:
            output_list.append('# ')
        output_list.append(executable.command_str())
        if isinstance(executable.stdout, ConnectorFile):
            output_list.append(' 1>' + executable.stdout.file_path)
        if isinstance(executable.stderr, ConnectorFile):
            output_list.append(' 2>' + executable.stderr.file_path)
        output_list.append('\n')
        output_list.append('\n')

    script_path = os.path.join(stage.working_directory, 'bsfpython_bash_' + repr(stage.name) + '.bash')
    with open(file=script_path, mode='wt') as script_file:
        script_file.writelines(output_list)

    return


def check_state(stage, debug=0):
    """Check the state of each C{bsf.process.Executable} of a C{bsf.analysis.Stage}.

    @param stage: C{bsf.analysis.Stage}
    @type stage: bsf.analysis.Stage
    @param debug: Debug level
    @type debug: int
    """

    if stage:
        pass

    if debug:
        pass

    return
