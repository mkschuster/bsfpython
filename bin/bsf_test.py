#! /usr/bin/env python
#
# BSF Python script to test library functions.
#
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
import re
import subprocess
import sys
import threading

import bsf.argument
from bsf.process import Executable


def bsf_test_argument():
    print 'Parsed via bsf.argument.Argument.from_key_value()'
    for key, value in (
            ('-switch_short', None),
            ('--switch_long', None),
            ('-option', 'short'),
            ('--option', 'long'),
            ('option', 'pair'),
            ('option', 'pair space'),
            ('-option_pair=short', None),
            ('--option_pair=long', None),
            ('--option=pair space', None)):
        argument = bsf.argument.Argument.from_key_value(key=key, value=value)
        print argument, argument.get_str(), argument.get_list()

    print

    print 'Instantiated directly via class method'
    argument_list = [
        bsf.argument.Switch(key='key'),
        bsf.argument.SwitchLong(key='key'),
        bsf.argument.SwitchShort(key='key'),
        # Single values
        bsf.argument.Option(key='key', value='value'),
        bsf.argument.OptionLong(key='key', value='value'),
        bsf.argument.OptionShort(key='key', value='value'),
        bsf.argument.OptionPair(key='key', value='value'),
        bsf.argument.OptionPairShort(key='key', value='value'),
        bsf.argument.OptionPairLong(key='key', value='value'),
        bsf.argument.OptionMulti(key='key', value='value'),
        bsf.argument.OptionMultiLong(key='key', value='value'),
        bsf.argument.OptionMultiShort(key='key', value='value'),
        bsf.argument.OptionMultiPair(key='key', value='value'),
        bsf.argument.OptionMultiPairLong(key='key', value='value'),
        bsf.argument.OptionMultiPairShort(key='key', value='value'),
        # Multiple values
        bsf.argument.Option(key='key', value='value1 value2'),
        bsf.argument.OptionLong(key='key', value='value1 value2'),
        bsf.argument.OptionShort(key='key', value='value1 value2'),
        bsf.argument.OptionPair(key='key', value='value1 value2'),
        bsf.argument.OptionPairShort(key='key', value='value1 value2'),
        bsf.argument.OptionPairLong(key='key', value='value1 value2'),
        bsf.argument.OptionMulti(key='key', value='value1 value2'),
        bsf.argument.OptionMultiLong(key='key', value='value1 value2'),
        bsf.argument.OptionMultiShort(key='key', value='value1 value2'),
        bsf.argument.OptionMultiPair(key='key', value='value1 value2'),
        bsf.argument.OptionMultiPairLong(key='key', value='value1 value2'),
        bsf.argument.OptionMultiPairShort(key='key', value='value1 value2'),
    ]

    for argument in argument_list:
        print argument, argument.get_str(), argument.get_list()

    print

    return


def bsf_test_thread(file_handle, thread_lock, debug, executable):
    """Thread callable to parse the SLURM process identifier and set it in the C{Executable}.

    @param file_handle: File handle (i.e. pipe)
    @type file_handle: file
    @param thread_lock: Thread lock
    @type thread_lock: threading.lock
    @param debug:
    @param executable: C{Executable}
    @type executable: bsf.process.Executable
    @return:
    """
    for line in file_handle:
        match = re.search(pattern=r'Submitted batch job (\d+)', string=line)
        thread_lock.acquire(True)
        if debug > 0:
            print 'Line:', line
        if match:
            executable.process_identifier = match.group(1)
        else:
            print 'Could not parse the process identifier from the SLURM sbatch response line', line
            executable.process_identifier = '9999999'
            executable.process_name = 'testing successful'
        thread_lock.release()

    return


def bsf_test_subprocess(
        maximum_attempts=3,
        max_thread_joins=10,
        thread_join_timeout=10,
        debug=1):
    """Test the Python subprocess module.

    @param maximum_attempts: Maximum number of attempts to run
    @type maximum_attempts: int
    @param max_thread_joins: Maximum number of attempts to join threads
    @type max_thread_joins: int
    @param thread_join_timeout: Thread join timeout
    @type thread_join_timeout: int
    @param debug: Debug level
    @type debug: int
    @return: Child return code
    @rtype: int
    """
    executable = Executable(name='bsf_test', program='ls')
    executable.add_switch_short(key='a')
    executable.add_switch_short(key='l')
    # Conclusion:
    # The following re-directions only work, if the Popen() call is allowed to use the shell i.e. Popen(shell=True)
    # '1>bsf_test_stdout_redirect.txt'
    # '2>bsf_test_stderr_redirect.txt'
    # However, if the shell is used, the command has to be specified as a single string rather than a list.
    # stdout_path = "bsf_test_stdout.txt"
    stderr_path = "bsf_test_stderr.txt"

    on_posix = 'posix' in sys.builtin_module_names

    child_return_code = 0
    attempt_counter = 0

    while attempt_counter < maximum_attempts:
        child_process = subprocess.Popen(
            args=executable.command_list(),
            bufsize=4096,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False,
            close_fds=on_posix)

        thread_lock = threading.Lock()

        # thread_out = Thread(
        #     target=Executable.process_stdout,
        #     kwargs={
        #         'stdout_handle': child_process.stdout,
        #         'thread_lock': thread_lock,
        #         'stdout_path': stdout_path,
        #         'debug': debug,
        #     })
        thread_out = threading.Thread(
            target=bsf_test_thread,
            kwargs={
                'file_handle': child_process.stdout,
                'thread_lock': thread_lock,
                'debug': debug,
                'executable': executable,
            })
        thread_out.daemon = True  # Thread dies with the program.
        thread_out.start()

        thread_err = threading.Thread(
            target=Executable.process_stderr,
            kwargs={
                'stderr_handle': child_process.stderr,
                'thread_lock': thread_lock,
                'stderr_path': stderr_path,
                'debug': debug,
            })
        thread_err.daemon = True  # Thread dies with the program.
        thread_err.start()

        # Wait for the child process to finish.

        child_return_code = child_process.wait()

        thread_join_counter = 0

        while thread_out.is_alive() and thread_join_counter < max_thread_joins:
            thread_lock.acquire(True)
            if debug > 0:
                print '[{}] Waiting for STDOUT processor to finish.'. \
                    format(datetime.datetime.now().isoformat())
            thread_lock.release()

            thread_out.join(timeout=thread_join_timeout)
            thread_join_counter += 1

        thread_join_counter = 0

        while thread_err.is_alive() and thread_join_counter < max_thread_joins:
            thread_lock.acquire(True)
            if debug > 0:
                print '[{}] Waiting for STDERR processor to finish.'. \
                    format(datetime.datetime.now().isoformat())
            thread_lock.release()

            thread_err.join(timeout=thread_join_timeout)
            thread_join_counter += 1

        if child_return_code > 0:
            if debug > 0:
                print '[{}] Child process {!r} failed with exit code {}'. \
                    format(datetime.datetime.now().isoformat(), executable.name, +child_return_code)
            attempt_counter += 1
        elif child_return_code < 0:
            if debug > 0:
                print '[{}] Child process {!r} received signal {}.'. \
                    format(datetime.datetime.now().isoformat(), executable.name, -child_return_code)
        else:
            if debug > 0:
                print '[{}] Child process {!r} completed successfully {}.'. \
                    format(datetime.datetime.now().isoformat(), executable.name, +child_return_code)
            break
    else:
        if debug > 0:
            print '[{}] Runnable {!r} exceeded the maximum retry counter {}.'. \
                format(datetime.datetime.now().isoformat(), executable.name, maximum_attempts)

    print executable.trace(level=0)

    return child_return_code


bsf_test_argument()
