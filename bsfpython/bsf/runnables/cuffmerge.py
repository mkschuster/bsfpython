"""bsf.runnables.cuffmerge

A package of functions supporting the BSF Runner script.
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


executable = 'cuffmerge'  # TODO: This is only needed until a better solution has been found.

# TODO: Both, STDOUT and STDERR of the child are not in the STDOUT of the BSF Runner.


def process_stdout(stdout_handle, lock):
    """BSF Runner function to process STDOUT from the child process.

    :param stdout_handle: The STDOUT file handle
    :type stdout_handle: file
    :param lock: A Python threading.Lock object
    :type lock: Lock
    :return: Nothing
    :rtype: None
    """

    lock.acquire(True)
    print '[{}] Started BSF Runner STDOUT processor in module {}.'. \
        format(datetime.datetime.now().isoformat(), __name__)
    lock.release()

    for line in stdout_handle:
        lock.acquire(True)
        print '[{}] STDOUT: {}'.format(datetime.datetime.now().isoformat(), line.rstrip())
        lock.release()

    lock.acquire(True)
    print '[{}] STDOUT: EOF on pipe.'.format(datetime.datetime.now().isoformat())
    lock.release()


def process_stderr(stderr_handle, lock):
    """BSF Runner function to process STDERR from the child process.

    :param stderr_handle: The STDERR file handle
    :type stderr_handle: file
    :param lock: A Python threading.Lock object
    :type lock: Lock
    :return: Nothing
    :rtype: None
    """

    lock.acquire(True)
    print '[{}] Started BSF Runner STDERR processor in module {}.'. \
        format(datetime.datetime.now().isoformat(), __name__)
    lock.release()

    for line in stderr_handle:
        lock.acquire(True)
        print '[{}] STDERR: {}'.format(datetime.datetime.now().isoformat(), line.rstrip())
        lock.release()

    lock.acquire(True)
    print '[{}] STDERR: EOF on pipe.'.format(datetime.datetime.now().isoformat())
    lock.release()
