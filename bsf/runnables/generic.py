# -*- coding: utf-8 -*-
"""Generic Runnables module

A package of classes and methods to run C{bsf.process.RunnableStep} objects of a C{bsf.Runnable}.
Empty status files keep track of completed C{bsf.process.RunnableStep} objects and allow for
restarting of the C{bsf.Runnable} object processing.
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

import datetime
import errno
import os
import shutil

import bsf.process


def run(runnable):
    """Run the the C{bsf.Runnable}.

    @param runnable: C{bsf.Runnable}
    @type runnable: bsf.Runnable
    @return:
    @rtype:
    """

    # If the Runnable status file exists, there is nothing to do and
    # this Runnable should not have been submitted in the first place.

    if os.path.exists(runnable.runnable_status_file_path(success=True, absolute=False)):
        return

    # Create a Runnable-specific cache directory if it does not exist already.

    if runnable.cache_path_dict:
        cache_directory_path = runnable.cache_directory_path(absolute=True)
        if not os.path.isdir(cache_directory_path):
            try:
                os.makedirs(cache_directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        # Copy files from the cache dictionary into the cache directory.

        for key in runnable.cache_path_dict:
            source_path = runnable.cache_path_dict[key]
            target_path = os.path.join(cache_directory_path, os.path.basename(source_path))

            runnable_step = bsf.process.RunnableStep(name='_'.join((runnable.name, 'cache', key)), program='cp')
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
    else:
        cache_directory_path = None

    # Create a Runnable-specific temporary directory if it does not exist already.

    temporary_directory_path = runnable.temporary_directory_path(absolute=False)
    if not os.path.isdir(temporary_directory_path):
        try:
            os.makedirs(temporary_directory_path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Check the Python list of bsf.process.RunnableStep objects in reverse order to see what has completed already.
    # If a bsf.process.RunnableStep is complete, it will be the first one on the list to become the
    # previous bsf.process.RunnableStep.

    runnable_step_list = list()
    """ @type runnable_step_list: list[bsf.process.RunnableStep] """
    for runnable_step in reversed(runnable.runnable_step_list):
        runnable_step_list.append(runnable_step)
        if os.path.exists(path=runnable.runnable_step_status_file_path(
                runnable_step=runnable_step,
                success=True)):
            break
    runnable_step_list.reverse()

    # Work through the list of bsf.process.RunnableStep objects in logical order.
    # Keep track of the previous bsf.process.RunnableStep.

    child_return_code = 0
    runnable_step_current = None
    runnable_step_previous = None

    for runnable_step_current in runnable_step_list:
        # Check for a bsf.process.RunnableStep-specific status file.
        if os.path.exists(path=runnable.runnable_step_status_file_path(
                runnable_step=runnable_step_current,
                success=True)):
            # If a status file exists, this RunnableStep is complete.
            # Set it as the previous RunnableStep and continue with the next RunnableStep.
            runnable_step_previous = runnable_step_current
            continue

        # Do the work.

        child_return_code = runnable_step_current.run()

        # Upon failure, break out of this loop without deleting obsolete files or altering status files at this point.

        if child_return_code != 0:
            break

        # Delete the list of file paths that the current bsf.process.RunnableStep declared to be obsolete now.

        runnable_step_current.remove_obsolete_file_paths()

        # Create an empty status file upon success.

        runnable.runnable_step_status_file_create(runnable_step=runnable_step_current, success=True)

        # Remove the status file of the previous bsf.process.RunnableStep, if it has been defined at this stage.

        if runnable_step_previous is not None:
            runnable.runnable_step_status_file_remove(runnable_step=runnable_step_previous)

        # Finally, make the current RunnableStep the previous RunnableStep.

        runnable_step_previous = runnable_step_current

    # Irrespective of failure ...

    # Remove the Runnable-specific cache directory and everything within it.

    if cache_directory_path and os.path.exists(cache_directory_path):
        shutil.rmtree(path=cache_directory_path, ignore_errors=False)

    # Remove the Runnable-specific temporary directory and everything within it.

    if os.path.exists(temporary_directory_path):
        shutil.rmtree(path=temporary_directory_path, ignore_errors=False)

    if child_return_code == 0:
        # Upon success, create a Runnable-specific status file that indicates completion for the whole Runnable.

        runnable.runnable_status_file_remove()
        runnable.runnable_status_file_create(success=True)

        # Remove the status file of the previous RunnableStep, if it has been defined at this stage.

        if runnable_step_previous is not None:
            runnable.runnable_step_status_file_remove(runnable_step=runnable_step_previous)
    else:
        # Upon failure, create a RunnableStep-specific status file showing failure and raise an Exception.

        if runnable_step_current is not None:
            runnable.runnable_step_status_file_create(runnable_step=runnable_step_current, success=False)

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
