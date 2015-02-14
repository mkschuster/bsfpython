"""bsf.runnables.rnaseq_cuffquant

A package of classes and methods to run Cuffquant of the Cufflinks package.
"""

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


import errno
import os
import shutil

from bsf import Runnable


def run_cuffquant(runnable):
    """Run the I{Cuffquant} step.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    # if os.path.exists(runnable.file_path_dict['annotated_idx']):
    #     return

    runnable.run_executable(name='cuffquant')


def run(runnable):
    """Run the the C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    """

    # Create a temporary directory.

    path_temporary = runnable.file_path_dict['temporary_directory']

    if not os.path.isdir(path_temporary):
        try:
            os.makedirs(path_temporary)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Run the chain of executables back up the function hierarchy so that
    # dependencies on temporarily created files become simple to manage.

    run_cuffquant(runnable=runnable)

    # Remove the temporary directory and everything within it.

    shutil.rmtree(path=path_temporary, ignore_errors=False)

    # Job done.
