"""bsf.runnables.generic

A package of classes and methods to run RunnableStep objects of a Runnable.
Empty status files keep track of completed RunnableStep objects and allow for
restarting of the Runnable object processing.
"""

#
# Copyright 2013 - 2016 Michael K. Schuster
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


import datetime
import errno
import os
import shutil

from bsf import Runnable
from bsf.process import RunnableStep


def _runnable_step_remove_obsolete_file_paths(runnable_step):
    """Remove the list of file path objects that the RunnableStep declared to be obsolete.

    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @return: Nothing
    @rtype: None
    """
    if runnable_step is None:
        return

    for file_path in runnable_step.obsolete_file_path_list:
        assert isinstance(file_path, (str, unicode))
        if os.path.exists(file_path):
            os.remove(file_path)

    return


def _runnable_step_status_file_path(runnable, runnable_step, success=True):
    """Get the status file path for a C{RunnableStep} of a C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @param success: Successful completion
    @type success: bool
    @return: Status file path
    @rtype: str
    """
    if success:
        return '_'.join((runnable.name, runnable_step.name, 'completed.txt'))
    else:
        return '_'.join((runnable.name, runnable_step.name, 'failed.txt'))


def _runnable_step_status_file_create(runnable, runnable_step, success=True):
    """Create an empty status file for a C{RunnableStep} of a C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @param success: Successful completion
    @type success: bool
    @return: Nothing
    @rtype: None
    """
    if runnable_step is None:
        return

    status_path = _runnable_step_status_file_path(runnable=runnable, runnable_step=runnable_step, success=success)
    open(status_path, 'w').close()

    return


def _runnable_step_status_file_remove(runnable, runnable_step):
    """Remove the status file for a C{RunnableStep} of a C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @return: Nothing
    @rtype: None
    """

    if runnable_step is None:
        return

    # Automatically remove both status files, successful or not.

    status_path = _runnable_step_status_file_path(runnable=runnable, runnable_step=runnable_step, success=True)
    if os.path.exists(status_path):
        os.remove(status_path)

    status_path = _runnable_step_status_file_path(runnable=runnable, runnable_step=runnable_step, success=False)
    if os.path.exists(status_path):
        os.remove(status_path)

    return


def run(runnable):
    """Run the the C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @return: Nothing
    @rtype: None
    """

    # If the Runnable status file exists, there is nothing to do and
    # this Runnable should not have been submitted in the first place.

    if os.path.exists(runnable.get_relative_status_path):
        return

    # Create a Runnable-specific cache directory if it does not exist already.

    cache_directory_path = runnable.get_absolute_cache_directory_path
    if not os.path.isdir(cache_directory_path):
        try:
            os.makedirs(cache_directory_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Copy files from the cache dictionary into the cache directory.

    for key in runnable.cache_path_dict:
        source_path = runnable.cache_path_dict[key]
        target_path = os.path.join(runnable.get_absolute_cache_directory_path, os.path.basename(source_path))

        runnable_step = RunnableStep(name='_'.join((runnable.name, 'cache', key)), program='cp')
        runnable_step.add_switch_short(key='p')
        runnable_step.arguments.append(source_path)
        runnable_step.arguments.append(target_path)

        child_return_code = runnable_step.run()

        if child_return_code != 0:
            # Remove the Runnable-specific cache directory and everything within it.
            if os.path.exists(cache_directory_path):
                shutil.rmtree(path=cache_directory_path, ignore_errors=False)
            # Raise an Exception.
            if child_return_code > 0:
                raise Exception('[{}] Child process {}_{} failed with return code {}'.
                                format(datetime.datetime.now().isoformat(),
                                       runnable.name,
                                       runnable_step.name,
                                       +child_return_code))
            elif child_return_code < 0:
                raise Exception('[{}] Child process {}_{} received signal {}.'.
                                format(datetime.datetime.now().isoformat(),
                                       runnable.name,
                                       runnable_step.name,
                                       -child_return_code))

    # Create a Runnable-specific temporary directory if it does not exist already.

    temporary_directory_path = runnable.get_relative_temporary_directory_path
    if not os.path.isdir(temporary_directory_path):
        try:
            os.makedirs(temporary_directory_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Check the Python list of RunnableStep objects in reverse order to see what has completed already.
    # If a RunnableStep is complete, it will be the first one on the list to become the previous RunnableStep.

    runnable_step_list = list()
    for runnable_step in reversed(runnable.runnable_step_list):
        assert isinstance(runnable_step, RunnableStep)
        runnable_step_list.append(runnable_step)
        if os.path.exists(path=_runnable_step_status_file_path(runnable=runnable, runnable_step=runnable_step)):
            break
    runnable_step_list.reverse()

    # Work through the list of RunnableStep objects in logical order. Keep track of the previous RunnableStep.

    child_return_code = 0
    runnable_step_current = None
    runnable_step_previous = None

    for runnable_step_current in runnable_step_list:
        assert isinstance(runnable_step_current, RunnableStep)

        # Check for a RunnableStep-specific status file.
        status_path = _runnable_step_status_file_path(
            runnable=runnable, runnable_step=runnable_step_current, success=True)
        if os.path.exists(path=status_path):
            # If a status file exists, this RunnableStep is complete.
            # Set it as the previous RunnableStep and continue with the next RunnableStep.
            runnable_step_previous = runnable_step_current
            continue

        # Do the work.

        child_return_code = runnable_step_current.run()

        # Upon failure, break out of this loop without deleting obsolete files or altering status files at this point.

        if child_return_code != 0:
            break

        # Delete the list of file paths that the current RunnableStep declared to be obsolete now.

        _runnable_step_remove_obsolete_file_paths(runnable_step=runnable_step_current)

        # Create an empty status file upon success.

        _runnable_step_status_file_create(runnable=runnable, runnable_step=runnable_step_current, success=True)

        # Remove the status file of the previous RunnableStep.

        _runnable_step_status_file_remove(runnable=runnable, runnable_step=runnable_step_previous)

        # Finally, make the current RunnableStep the previous RunnableStep.

        runnable_step_previous = runnable_step_current

    # Irrespective of failure ...

    # Remove the Runnable-specific cache directory and everything within it.

    if os.path.exists(cache_directory_path):
        shutil.rmtree(path=cache_directory_path, ignore_errors=False)

    # Remove the Runnable-specific temporary directory and everything within it.

    if os.path.exists(temporary_directory_path):
        shutil.rmtree(path=temporary_directory_path, ignore_errors=False)

    if child_return_code == 0:
        # Create a status file that indicates completion for the whole Runnable.

        open(runnable.get_relative_status_path, 'w').close()

        # Remove the status file of the previous RunnableStep.

        _runnable_step_status_file_remove(runnable=runnable, runnable_step=runnable_step_previous)

        # Job done.

        return
    else:
        _runnable_step_status_file_create(runnable=runnable, runnable_step=runnable_step_current, success=False)

        # Upon failure, create a RunnableStep-specific status file showing failure and raise an Exception.

        if child_return_code > 0:
            raise Exception('[{}] Child process {}_{} failed with return code {}'.
                            format(datetime.datetime.now().isoformat(),
                                   runnable.name,
                                   runnable_step_current.name,
                                   +child_return_code))
        elif child_return_code < 0:
            raise Exception('[{}] Child process {}_{} received signal {}.'.
                            format(datetime.datetime.now().isoformat(),
                                   runnable.name,
                                   runnable_step_current.name,
                                   -child_return_code))

        return
