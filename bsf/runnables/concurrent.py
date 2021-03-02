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
import errno
import os
from io import TextIOWrapper
from subprocess import Popen, DEVNULL, PIPE
from threading import Lock, Thread

from bsf.connector import *
from bsf.procedure import ConcurrentRunnable
from bsf.process import Executable, RunnableStep, get_timestamp


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
            # Remove the ConcurrentRunnable-specific cache directory and everything within it.
            runnable.cache_directory_remove()
            # Remove the ConcurrentRunnable-specific temporary directory and everything within it.
            runnable.temporary_directory_remove()
            runnable.runnable_status_file_create(success=False)
            raise _exception

        return

    def _map_connector(_connector):
        """Map a connector to a file handle.

        @param _connector: C{bsf.connector.Connector} or sub-class thereof.
        @type _connector: Connector
        @return: File handle
        @rtype: TextIOWrapper | PIPE
        """
        if isinstance(_connector, ElectronicSink):
            return DEVNULL

        if isinstance(_connector, ConnectorFile):
            if isinstance(_connector, ConnectorPipeNamed):
                # A named pipe needs creating before it can be opened.
                if not os.path.exists(_connector.file_path):
                    os.mkfifo(_connector.file_path)
            return open(file=_connector.file_path, mode=_connector.file_mode)

        if isinstance(_connector, ConnectorPipe):
            return PIPE

        if isinstance(_connector, StandardStream):
            return PIPE

        if isinstance(_connector, ConcurrentProcess):
            for _runnable_step in runnable.runnable_step_list_concurrent:
                if _runnable_step.name == _connector.name:
                    if _runnable_step.sub_process is None:
                        raise Exception(
                            'Sub-process ' + repr(_runnable_step.name) + ' not initialised, yet.')
                    return getattr(_runnable_step.sub_process, _connector.connection)
            else:
                raise Exception('Could not find a suitable RunnableStep for Connector.name ' + repr(_connector.name))

        return None

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

    # Now, process the RunnableStep objects on the pre-run list.

    _run_consecutively(runnable_step_list=runnable.runnable_step_list_pre)

    thread_lock = Lock()

    for runnable_step in runnable.runnable_step_list_concurrent:
        # Per RunnableStep, up to three threads, writing to STDIN, as well as reading from STDOUT and STDERR,
        # should make sure that buffers are not filling up. If STDOUT or STDERR connectors are not defined,
        # defaults need be set to avoid the sub-process from blocking.

        if runnable_step.stdout is None:
            runnable_step.stdout = StandardOutputStream()

        if runnable_step.stderr is None:
            runnable_step.stderr = StandardErrorStream()

        # Create and assign sub-processes.
        try:
            runnable_step.sub_process = Popen(
                args=runnable_step.command_list(),
                bufsize=1,
                stdin=_map_connector(_connector=runnable_step.stdin),
                stdout=_map_connector(_connector=runnable_step.stdout),
                stderr=_map_connector(_connector=runnable_step.stderr),
                text=True)
        except OSError as exception:
            if exception.errno == errno.ENOENT:
                raise Exception(
                    repr(runnable_step) + ' ' +
                    repr(runnable_step.program) +
                    ' could not be found.')
            else:
                # Re-raise the Exception object.
                raise exception

        for attribute in ('stdin', 'stdout', 'stderr'):
            connector = getattr(runnable_step, attribute)
            """ @type connector: Connector """

            if isinstance(connector, StandardInputStream):
                if connector.thread_callable is None:
                    # If a specific STDIN callable is not defined, run bsf.process.Executable.process_stdin().
                    pass
                else:
                    connector.thread = Thread(
                        target=connector.thread_callable,
                        args=[runnable_step.sub_process.stdin, thread_lock, runnable.debug],
                        kwargs=connector.thread_kwargs)

            if isinstance(connector, StandardOutputStream):
                if connector.thread_callable is None:
                    # If a specific STDOUT callable is not defined, run bsf.process.Executable.process_stdout().
                    connector.thread = Thread(
                        target=Executable.process_stdout,
                        args=[runnable_step.sub_process.stdout, thread_lock, runnable.debug],
                        kwargs={'stdout_path': connector.file_path})
                else:
                    connector.thread = Thread(
                        target=connector.thread_callable,
                        args=[runnable_step.sub_process.stdout, thread_lock, runnable.debug],
                        kwargs=connector.thread_kwargs)

            if isinstance(connector, StandardErrorStream):
                if connector.thread_callable is None:
                    # If a specific STDERR callable is not defined, run bsf.process.Executable.process_stderr().
                    connector.thread = Thread(
                        target=Executable.process_stderr,
                        args=[runnable_step.sub_process.stderr, thread_lock, runnable.debug],
                        kwargs={'stderr_path': connector.file_path})
                else:
                    connector.thread = Thread(
                        target=connector.thread_callable,
                        args=[runnable_step.sub_process.stderr, thread_lock, runnable.debug],
                        kwargs=connector.thread_kwargs)

            if isinstance(connector, StandardStream) and connector.thread:
                connector.thread.daemon = True
                connector.thread.start()

    # At this stage all subprocess.Popen objects should have been created.
    # Now, wait for all child processes to complete.

    exception_str_list = list()
    for runnable_step in runnable.runnable_step_list_concurrent:
        child_return_code = runnable_step.sub_process.wait()

        # First, join all standard stream processing threads.

        for attribute in ('stdin', 'stdout', 'stderr'):
            connector = getattr(runnable_step, attribute)
            """ @type connector: Connector """
            if isinstance(connector, StandardStream) and connector.thread:
                thread_join_counter = 0
                while connector.thread.is_alive() and thread_join_counter < connector.thread_joins:
                    if runnable.debug > 0:
                        thread_lock.acquire(True)
                        print(get_timestamp(), 'Waiting for ' + repr(attribute) + ' processor to finish.')
                        thread_lock.release()

                    connector.thread.join(timeout=connector.thread_timeout)
                    thread_join_counter += 1

        # Second, inspect the child process' return code.

        if child_return_code > 0:
            # Child return code.
            exception_str_list.append(
                get_timestamp() +
                ' Child process ' + repr(runnable.name) + ' ' + repr(runnable_step.name) +
                ' failed with return code ' +
                repr(+child_return_code) + '.')
        elif child_return_code < 0:
            # Child signal.
            exception_str_list.append(
                get_timestamp() +
                ' Child process ' + repr(runnable.name) + ' ' + repr(runnable_step.name) +
                ' received signal ' +
                repr(-child_return_code) + '.')
        else:
            # Upon success, delete the list of file paths that the bsf.process.RunnableStep
            # declared to be obsolete now.

            runnable_step.remove_obsolete_file_paths()

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
