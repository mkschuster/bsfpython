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
"""The :py:mod:`bsf.drms.bash` module supports the GNU Bourne-Again Shell (Bash) as
Distributed Resource Management System (DRMS) module.
"""
import os
from typing import Optional

from bsf.analysis import Stage
from bsf.connector import ConnectorFile


def submit(stage, drms_submit: Optional[bool] = None, write_script: Optional[bool] = None) -> None:
    """Submit each :py:class:`bsf.process.Executable` object of a :py:class:`bsf.analysis.Stage` object.

    Submits each :py:class:`bsf.process.Executable` object by writing a
    :emphasis:`GNU Bourne-Again Shell` (BASH) script
    into the :py:attr:`bsf.analysis.Stage.work_directory`.

    :param stage: A :py:class:`bsf.analysis.Stage` object.
    :type stage: bsf.analysis.Stage
    :param drms_submit: Submit to the DRMS.
    :type drms_submit: bool | None
    :param write_script: Write a :py:class:`bsf.analysis.Stage`-specific GNU Bash script.
    :type write_script: bool | None
    """
    # Since write_script is required for compatibility of the function interface with the
    # SLURM and SGE modules, it always needs overriding here.
    if not write_script:
        write_script = True

    output_list: list[str] = list()

    output_list.append('#!/usr/bin/env bash\n')
    output_list.append('\n')

    for executable in stage.executable_list:
        if not (executable.submit and drms_submit):
            output_list.append('# ')
        output_list.append(executable.command_str())
        if isinstance(executable.stdout, ConnectorFile):
            output_list.append(' 1>' + executable.stdout.file_path)
        if isinstance(executable.stderr, ConnectorFile):
            output_list.append(' 2>' + executable.stderr.file_path)
        output_list.append('\n')
        output_list.append('\n')

    if write_script:
        with open(
                file=os.path.join(stage.working_directory, 'bsfpython_bash_' + repr(stage.name) + '.bash'),
                mode='wt') as output_text_io:
            output_text_io.writelines(output_list)

    return


def check_state(stage):
    """Check the state of each :py:class:`bsf.process.Executable` object of a :py:class:`bsf.analysis.Stage` object.

    :param stage: A :py:class:`bsf.analysis.Stage` object.
    :type stage: Stage
    """
    if stage:
        pass

    return
