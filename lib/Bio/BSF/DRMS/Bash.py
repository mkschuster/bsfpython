"""Bio.BSF.DRMS.Bash

A package of methods supporting the GNU Bourne-Again Shell (Bash).
"""

#
# Copyright 2013 Michael K. Schuster
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


import os


def submit(drms, debug=0):
    """Submit C{Executable} objects by writing a GNU Bourne-Again Shell (BASH) script into the C{DRMS.work_directory}.

    @param drms: Distributed Resource Management System (C{DRMS})
    @type drms: DRMS
    @param debug: Debug level
    @type debug: int
    """

    output = str()
    output += "#! /bin/bash\n"
    output += "\n"

    if debug > 0:
        output += "# BSF-Python debug mode: {}\n".format(debug)
        output += "\n"

    for executable in drms.executables:
        if not executable.submit:
            output += '# '
        output += executable.command_str()
        if executable.stdout_path:
            output += " 1>{}".format(executable.stdout_path)
        if executable.stderr_path:
            output += " 2>{}".format(executable.stderr_path)
        output += "\n"
        output += "\n"

    script_path = os.path.join(drms.work_directory, 'bsfpython_bash_{}.bash'.format(drms.name))
    script_file = open(name=script_path, mode='w')
    script_file.write(output)
    script_file.close()
