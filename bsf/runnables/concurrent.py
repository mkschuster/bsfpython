# -*- coding: utf-8 -*-
"""Generic ConcurrentRunnable module

A package of classes and methods to run C{bsf.process.RunnableStep} objects of a C{bsf.procedure.ConcurrentRunnable}.
Empty status files keep track of completed C{bsf.process.RunnableStep} objects and allow for
restarting of the C{bsf.procedure.ConcurrentRunnable} object processing.
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
import errno
import os
import subprocess

import bsf.procedure


def run(runnable):
    """Run the the C{bsf.procedure.ConcurrentRunnable}.

    @param runnable: C{bsf.procedure.ConcurrentRunnable}
    @type runnable: bsf.procedure.ConcurrentRunnable
    @return:
    @rtype:
    """

    def _run_consecutively(runnable_step_list):
        """Run a Python C{list} of C{bsf.process.RunnableStep} objects consecutively.

        @param runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
        @type runnable_step_list: list[bsf.process.RunnableStep]
        @return:
        @rtype:
        """
        _exception = runnable.run_consecutively(runnable_step_list=runnable_step_list)

        if _exception is not None:
            # Remove the ConcurrentRunnable-specific cache directory and everything within it.
            runnable.cache_directory_remove()
            # Remove the ConcurrentRunnable-specific temporary directory and everything within it.
            runnable.temporary_directory_remove()
            runnable.runnable_status_file_create(success=False)
            raise _exception

        return

    def _map_connector(_connector):
        """Map a connector to a file handle.

        @param _connector: C{Connector} or sub-class thereof.
        @type _connector: procedure.Connector
        @return: File handle
        @rtype: file | subprocess.PIPE
        """
        if isinstance(_connector, bsf.procedure.ConnectorSink):
            return open('/dev/null', 'wb')

        if isinstance(_connector, bsf.procedure.ConnectorFile):
            if isinstance(_connector, bsf.procedure.ConnectorPipeNamed):
                # A named pipe needs creating before it can be opened.
                if not os.path.exists(_connector.file_path):
                    os.mkfifo(_connector.file_path)
            return open(_connector.file_path, _connector.file_mode)

        if isinstance(_connector, bsf.procedure.ConnectorPipe):
            return subprocess.PIPE

        if isinstance(_connector, bsf.procedure.ConnectorProcess):
            for _sub_process in runnable.sub_process_list:
                if _sub_process.runnable_step.name == _connector.name:
                    if _sub_process.sub_process is None:
                        raise Exception(
                            'Sub-process ' + repr(_sub_process.runnable_step.name) + ' not initialised, yet.')
                    return getattr(_sub_process.sub_process, _connector.connection)

        return None

    # If the ConcurrentRunnable status file exists, there is nothing to do and
    # this ConcurrentRunnable should not have been submitted in the first place.

    if os.path.exists(runnable.runnable_status_file_path(success=True, absolute=False)):
        return

    # Create and populate a ConcurrentRunnable-specific cache directory if it does not exist already.

    runnable.cache_directory_create()

    # Create a ConcurrentRunnable-specific temporary directory if it does not exist already.

    runnable.temporary_directory_create()

    # Now, process the RunnableStep objects on the pre-run list.

    _run_consecutively(runnable_step_list=runnable.runnable_step_list_pre)

    for sub_process in runnable.sub_process_list:
        # Create and assign sub-processes.
        try:
            sub_process.sub_process = subprocess.Popen(
                args=sub_process.runnable_step.command_list(),
                stdin=_map_connector(_connector=sub_process.stdin),
                stdout=_map_connector(_connector=sub_process.stdout),
                stderr=_map_connector(_connector=sub_process.stderr),
                close_fds=True,
                shell=False)
        except OSError as exception:
            if exception.errno == errno.ENOENT:
                raise Exception(
                    repr(sub_process.runnable_step) + ' ' +
                    repr(sub_process.runnable_step.program) +
                    ' could not be found.')
            else:
                # Re-raise the Exception object.
                raise exception

    # At this stage all subprocess.Popen objects should have been created.
    # Now, wait for all child processes to complete.

    exception_str_list = list()
    for sub_process in runnable.sub_process_list:
        child_return_code = sub_process.sub_process.wait()

        if child_return_code > 0:
            exception_str_list.append(
                bsf.process.get_timestamp() +
                ' Child process ' + repr(runnable.name) + ' ' + repr(sub_process.runnable_step.name) +
                ' failed with return code ' +
                repr(+child_return_code) + '.')
        elif child_return_code < 0:
            exception_str_list.append(
                bsf.process.get_timestamp() +
                ' Child process ' + repr(runnable.name) + ' ' + repr(sub_process.runnable_step.name) +
                ' received signal ' +
                repr(-child_return_code) + '.')
        else:
            # Delete the list of file paths that the bsf.process.RunnableStep declared to be obsolete now.

            sub_process.runnable_step.remove_obsolete_file_paths()

    if exception_str_list:
        runnable.cache_directory_remove()
        runnable.temporary_directory_remove()
        runnable.runnable_status_file_create(success=False)
        raise Exception('\n'.join(exception_str_list))

    # Finally, process the RunnableStep object on the post-run list.

    _run_consecutively(runnable_step_list=runnable.runnable_step_list_post)

    # Remove the ConcurrentRunnable-specific cache directory and everything within it.

    runnable.cache_directory_remove()

    # Remove the ConcurrentRunnable-specific temporary directory and everything within it.

    runnable.temporary_directory_remove()

    # Upon success, create a ConcurrentRunnable-specific status file that indicates completion
    # for the whole ConcurrentRunnable.
    runnable.runnable_status_file_remove()
    runnable.runnable_status_file_create(success=True)

    return
