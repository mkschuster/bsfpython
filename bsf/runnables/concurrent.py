# -*- coding: utf-8 -*-
"""Generic ConcurrentRunnable module.

A package of classes and methods to run C{bsf.process.RunnableStep} objects of a C{bsf.procedure.ConcurrentRunnable}.
Empty status files keep track of completed C{bsf.process.RunnableStep} objects and allow for
restarting of the C{bsf.procedure.ConcurrentRunnable} object processing.
"""
#  Copyright 2013 - 2021 Michael K. Schuster
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

from bsf.procedure import ConcurrentRunnable
from bsf.process import RunnableStep, run_executables


def run(runnable):
    """Run the the C{bsf.procedure.ConcurrentRunnable}.

    @param runnable: C{bsf.procedure.ConcurrentRunnable}
    @type runnable: ConcurrentRunnable
    """

    def _run_consecutively(runnable_step_list):
        """Run a Python C{list} of C{bsf.process.RunnableStep} objects consecutively.

        @param runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
        @type runnable_step_list: list[RunnableStep]
        """
        _exception = runnable.run_consecutively(runnable_step_list=runnable_step_list)

        if _exception is not None:
            # If an Exception occurred, remove the ConcurrentRunnable-specific cache and temporary directories
            # and everything within them and create a status file that indicates failure of the whole
            # ConcurrentRunnable.
            runnable.cache_directory_remove()
            runnable.temporary_directory_remove()
            runnable.runnable_status_file_create(success=False)
            raise _exception

        return

    # If the ConcurrentRunnable status file exists, there is nothing to do and
    # this ConcurrentRunnable should not have been submitted in the first place.

    if os.path.exists(runnable.runnable_status_file_path(success=True, absolute=False)):
        return

    # Create and populate a ConcurrentRunnable-specific cache directory if it does not exist already.

    runnable.cache_directory_create()

    # Create a ConcurrentRunnable-specific temporary directory if it does not exist already.

    runnable.temporary_directory_create()

    # Set conventional environment variables.

    os.environ['TMPDIR'] = runnable.temporary_directory_path(absolute=True)
    os.environ['TEMP'] = runnable.temporary_directory_path(absolute=True)
    os.environ['TMP'] = runnable.temporary_directory_path(absolute=True)

    # Run the RunnableStep objects on the pre-run list.

    _run_consecutively(runnable_step_list=runnable.runnable_step_list_prologue)

    # Run the RunnableStep objects on the concurrent list.

    exception_str_list = run_executables(runnable.runnable_step_list_concurrent, debug=runnable.debug)

    if exception_str_list:
        # If an exception occurred, remove the ConcurrentRunnable-specific cache and temporary directories
        # and everything within them and create a status file that indicates failure of the whole ConcurrentRunnable.
        runnable.cache_directory_remove()
        runnable.temporary_directory_remove()
        runnable.runnable_status_file_create(success=False)
        raise Exception('\n'.join(exception_str_list))
    else:
        # If no Exception occurred, delete the list of file paths that the bsf.process.RunnableStep
        # declared to be obsolete now.
        for runnable_step in runnable.runnable_step_list_concurrent:
            runnable_step.remove_obsolete_file_paths()

    # Run the RunnableStep object on the post-run list.

    _run_consecutively(runnable_step_list=runnable.runnable_step_list_epilogue)

    # Upon success, remove the ConcurrentRunnable-specific cache and temporary directories and everything
    # within them and create a status file that indicates completion of the whole ConcurrentRunnable.

    runnable.cache_directory_remove()
    runnable.temporary_directory_remove()
    runnable.runnable_status_file_remove()
    runnable.runnable_status_file_create(success=True)

    return
