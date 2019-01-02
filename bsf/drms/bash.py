# -*- coding: utf-8 -*-
"""GNU Bourne-Again Shell (Bash) DRMS module

A package of methods supporting the GNU Bourne-Again Shell (Bash) as
Distributed Resource Management System (DRMS) module.
"""
#  Copyright 2013 - 2019 Michael K. Schuster
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

import os


def submit(stage, debug=0):
    """Submit each C{bsf.process.Executable} of a C{bsf.Stage}.

    Submits each C{bsf.process.Executable} by writing a GNU Bourne-Again Shell (BASH) script
    into the C{bsf.Stage.work_directory}.

    @param stage: C{bsf.Stage}
    @type stage: bsf.Stage
    @param debug: Debug level
    @type debug: int
    @return:
    @rtype:
    """

    output = str()
    output += '#! /usr/bin/env bash\n'
    output += '\n'

    if debug > 0:
        output += '# BSF-Python debug mode: ' + repr(debug) + '\n'
        output += '\n'

    for executable in stage.executable_list:
        if not executable.submit:
            output += '# '
        output += executable.command_str()
        if executable.stdout_path:
            output += ' 1>' + executable.stdout_path
        if executable.stderr_path:
            output += ' 2>' + executable.stderr_path
        output += '\n'
        output += '\n'

    script_path = os.path.join(stage.working_directory, 'bsfpython_bash_' + repr(stage.name) + '.bash')
    script_file = open(script_path, 'w')
    script_file.write(output)
    script_file.close()

    return


def check_state(stage, debug=0):
    """Check the state of each C{bsf.process.Executable} of a C{bsf.Stage}.

    @param stage: C{bsf.Stage}
    @type stage: bsf.Stage
    @param debug: Debug level
    @type debug: int
    @return:
    @rtype:
    """

    if stage:
        pass

    if debug:
        pass

    return
