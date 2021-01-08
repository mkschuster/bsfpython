# -*- coding: utf-8 -*-
"""Generic ConsecutiveRunnable module.

A package of classes and methods to run C{bsf.process.RunnableStep} objects of a C{bsf.procedure.ConsecutiveRunnable}.
Empty status files keep track of completed C{bsf.process.RunnableStep} objects and allow for
restarting of the C{bsf.procedure.ConsecutiveRunnable} object processing.
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

from bsf.procedure import ConsecutiveRunnable


def run(runnable):
    """Run the the C{bsf.procedure.ConsecutiveRunnable}.

    @param runnable: C{bsf.procedure.ConsecutiveRunnable}
    @type runnable: ConsecutiveRunnable
    """

    # If the ConsecutiveRunnable status file exists, there is nothing to do and
    # this ConsecutiveRunnable should not have been submitted in the first place.

    if os.path.exists(runnable.runnable_status_file_path(success=True, absolute=False)):
        return

    # Create and populate a ConsecutiveRunnable-specific cache directory if it does not exist already.

    runnable.cache_directory_create()

    # Create a ConsecutiveRunnable-specific temporary directory if it does not exist already.

    runnable.temporary_directory_create()

    # Set conventional environment variables.

    os.environ['TMPDIR'] = runnable.temporary_directory_path(absolute=True)
    os.environ['TEMP'] = runnable.temporary_directory_path(absolute=True)
    os.environ['TMP'] = runnable.temporary_directory_path(absolute=True)

    exception = runnable.run()

    # Irrespective of failure ...

    # Remove the ConsecutiveRunnable-specific cache directory and everything within it.

    runnable.cache_directory_remove()

    # Remove the ConsecutiveRunnable-specific temporary directory and everything within it.

    runnable.temporary_directory_remove()

    runnable.runnable_status_file_remove()

    if exception is not None:
        runnable.runnable_status_file_create(success=False)
        raise exception
    else:
        runnable.runnable_status_file_create(success=True)

    return
