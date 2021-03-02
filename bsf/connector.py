# -*- coding: utf-8 -*-
"""Connector module.

A package of classes and methods modelling inter-process connectors.
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
from threading import Thread

__all__ = \
    'Connector', 'ConnectorFile', 'ConnectorPipe', 'ConnectorPipeNamed', \
    'ConcurrentProcess', 'ElectronicSink', \
    'StandardStream', 'StandardInputStream', 'StandardOutputStream', 'StandardErrorStream'


class Connector(object):
    """The C{Connector} class represents an abstract super-class of inter-process connections.

    Attributes:
    """

    pass


class ConnectorFile(Connector):
    """The C{ConnectorFile} class represents a C{file} connection.

    Attributes:
    @ivar file_path: File path
    @type file_path: str
    @ivar file_mode: File mode
    @type file_mode: str
    """

    def __init__(self, file_path, file_mode):
        """Initialise a C{ConnectorFile} object.

        @param file_path: File path
        @type file_path: str
        @param file_mode: File mode
        @type file_mode: str
        """
        super(ConnectorFile, self).__init__()

        self.file_path = file_path
        self.file_mode = file_mode

        return


class ConnectorPipe(Connector):
    """The C{ConnectorPipe} class represents an abstract pipe.

    Attributes:
    """

    pass


class ConnectorPipeNamed(ConnectorFile):
    """The C{ConnectorPipeNamed} class represents a named pipe.

    Attributes:
    """

    pass


class ConcurrentProcess(Connector):
    """The C{ConcurrentProcess} class represents a pipe to or from a concurrent process.

    Attributes:
    @ivar name: C{bsf.process.Executable.name}
    @type name: str
    @ivar connection: Connection type 'stdin', 'stdout' or 'stderr'
    @type connection: str
    """

    def __init__(self, name, connection):
        """Initialise a C{ConcurrentProcess} object.

        @param name: C{bsf.process.Executable.name}
        @type name: str
        @param connection: Connection type 'stdin', 'stdout' or 'stderr'
        @type connection: str
        """
        super(ConcurrentProcess, self).__init__()

        self.name = name
        self.connection = connection

        return


class ElectronicSink(Connector):
    """The C{ElectronicSink} class represents a C{file} connection to /dev/null.

    Attributes:
    """

    pass


class StandardStream(Connector):
    """The C{StandardStream} class represents a standard stream processed via C{threading.Thread}.

    Standard streams (i.e. STDIN, STDOUT, STDERR) require processing via a C{threading.Thread} to prevent buffers from
    filling up and subsequently sub-processes (C{subprocess.Popen}) from blocking.
    Attributes:
    @ivar file_path: File path
    @type file_path: str | None
    @ivar thread_callable: C{Callable} object for C{threading.Thread.target}
    @type thread_callable: object | None
    @ivar thread_kwargs: Python C{dict} of keyword arguments for C{threading.Thread.kwargs}
    @type thread_kwargs: dict[str, object] | None
    @ivar thread_joins: Maximum number of attempts calling C{threading.Thread.join}
    @type thread_joins: int
    @ivar thread_timeout: Timeout in seconds for calling C{threading.Thread.join}
    @type thread_timeout: int
    @ivar thread: C{threading.Thread} object
    @type thread: Thread | None
    """

    def __init__(
            self,
            file_path=None,
            thread_callable=None,
            thread_kwargs=None,
            thread_joins=10,
            thread_timeout=10,
            thread=None):
        """Initialise a C{StandardStream} object.

        @param file_path: File path
        @type file_path: str | None
        @param thread_callable: C{Callable} object for C{threading.Thread.target}
        @type thread_callable: object | None
        @param thread_kwargs: Python C{dict} of keyword arguments for C{threading.Thread.kwargs}
        @type thread_kwargs: dict[str, object] | None
        @param thread_joins: Maximum number of attempts calling C{threading.Thread.join}
        @type thread_joins: int
        @param thread_timeout: Timeout in seconds for calling C{threading.Thread.join}
        @type thread_timeout: int
        @param thread: C{threading.Thread} object
        @type thread: Thread | None
        """
        self.file_path = file_path
        self.thread_callable = thread_callable
        self.thread_kwargs = thread_kwargs
        self.thread_joins = thread_joins
        self.thread_timeout = thread_timeout
        self.thread = thread

        return


class StandardInputStream(StandardStream):
    """The C{StandardInputStream} class represents a STDIN stream processed via C{threading.Thread}.

    Attributes:
    """

    pass


class StandardOutputStream(StandardStream):
    """The C{StandardOutputStream} class represents a STDOUT stream processed via C{threading.Thread}.

    Attributes:
    """

    pass


class StandardErrorStream(StandardStream):
    """The C{StandardErrorStream} class represents a STDERR stream processed via C{threading.Thread}.

    Attributes:
    """

    pass
