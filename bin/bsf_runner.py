#! /usr/bin/env python
#
# BSF Python wrapper script that runs binary programs submitted into a
# Distributed Resource Management System (DRMS) as a Python subprocess.
# The wrapper allows for systematic capturing and parsing of
# standard output and standard error streams, as well as the return value.
#
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

import argparse
import datetime
import importlib
from subprocess import PIPE, Popen
import sys
from threading import Lock, Thread


max_loop_counter = 10
max_thread_counter = 10
threading_timeout = 10.0

on_posix = 'posix' in sys.builtin_module_names

parser = argparse.ArgumentParser(description='BSF runner script.')

parser.add_argument('--runnable', required=True,
                    help='BSF Runnable name')

(args, arguments) = parser.parse_known_args()

module = importlib.import_module(name='Bio.BSF.Runnable.' + args.runnable)

# TODO: How to set the binary executable as first option needs re-thinking.
# A BSF Executable, Command and subordinate command infrastructure exists.
# This module would need to remove the bsf_runner.py, as well as all options specific to the
# BSF Runner script from the arguments list.

arguments.insert(0, module.executable)

loop_counter = 0
child_return_code = 0

while loop_counter < max_loop_counter:

    child_process = Popen(args=arguments,
                          bufsize=1,
                          stdin=PIPE,
                          stdout=PIPE,
                          stderr=PIPE,
                          shell=False,
                          close_fds=on_posix)

    # Two threads, thread_out and thread_err reading STDOUT and STDERR, respectively,
    # should make sure that buffers are not filling up.

    thread_lock = Lock()

    thread_out = Thread(target=module.process_stdout,
                        kwargs={'stdout_handle': child_process.stdout, 'lock': thread_lock})
    thread_out.daemon = True  # Thread dies with the program.
    thread_out.start()

    thread_err = Thread(target=module.process_stderr,
                        kwargs={'stderr_handle': child_process.stderr, 'lock': thread_lock})
    thread_err.daemon = True  # Thread dies with the program.
    thread_err.start()

    # Wait for the child process to finish.

    child_return_code = child_process.wait()

    thread_counter = 0

    while thread_out.is_alive and thread_counter < max_thread_counter:

        thread_lock.acquire(True)
        print '[{}] Waiting for STDOUT processor to finish.'.format(datetime.datetime.now().isoformat())
        thread_lock.release()

        thread_out.join(threading_timeout)
        thread_counter += 1

    thread_counter = 0

    while thread_err.is_alive and thread_counter < max_thread_counter:

        thread_lock.acquire(True)
        print '[{}] Waiting for STDERR processor to finish.'.format(datetime.datetime.now().isoformat())
        thread_lock.release()

        thread_err.join(threading_timeout)
        thread_counter += 1

    # Close all pipes - not sure if it's required.
    child_process.stdin.close()
    child_process.stdout.close()
    child_process.stderr.close()

    if child_return_code > 0:
        print 'Child process {} returned exit code {}'. \
            format(module.executable, +child_return_code)
        loop_counter += 1
    elif child_return_code < 0:
        print 'Child process {} received signal {}.'. \
            format(module.executable, -child_return_code)
    else:
        print 'Child process {} completed successfully {}.'. \
            format(module.executable, +child_return_code)
        break

else:
    print 'This BSF runner exceeded the maximum re-run counter {}.' \
        .format(max_loop_counter)

if child_return_code < 0:
    # This process has been sent a signal.
    # Reset the return code to 1.
    child_return_code = 1

sys.exit(child_return_code)
