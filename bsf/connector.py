# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
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
"""The :py:mod:`bsf.connector` module provides classes modelling inter-process connectors.
"""
from threading import Thread
from typing import Any, Callable, Optional

__all__ = \
    'Connector', 'ConnectorFile', 'ConnectorPipe', 'ConnectorPipeNamed', \
    'ConcurrentProcess', 'ElectronicSink', \
    'StandardStream', 'StandardInputStream', 'StandardOutputStream', 'StandardErrorStream'


class Connector(object):
    """The :py:class:`bsf.connector.Connector` class represents an abstract super-class of inter-process connections.
    """

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}()'


class ConnectorFile(Connector):
    """The :py:class:`bsf.connector.ConnectorFile` class represents a file.

    :ivar file_path: A file path.
    :type file_path: str
    :ivar file_mode: A file mode.
    :type file_mode: str
    """

    def __init__(self, file_path: str, file_mode: str) -> None:
        """Initialise a :py:class:`bsf.connector.ConnectorFile` object.

        :param file_path: A file path.
        :type file_path: str
        :param file_mode: A file mode.
        :type file_mode: str
        """
        super(ConnectorFile, self).__init__()

        self.file_path = file_path
        self.file_mode = file_mode

        return

    def __repr__(self) -> str:
        return \
            f'{self.__class__.__name__}(' \
            f'file_path={self.file_path!r}, ' \
            f'file_mode={self.file_mode!r})'


class ConnectorPipe(Connector):
    """The :py:class:`bsf.connector.ConnectorPipe` class represents an abstract pipe.
    """

    pass


class ConnectorPipeNamed(ConnectorFile):
    """The :py:class:`bsf.connector.ConnectorPipeNamed` class represents a named pipe.
    """

    pass


class ConcurrentProcess(Connector):
    """The :py:class:`bsf.connector.ConcurrentProcess` class represents a pipe to or from a concurrent process.

    :ivar name: A name matching a :py:attr:`bsf.process.Executable.name` attribute of the
        concurrent process to connect to.
    :type name: str
    :ivar connection: A connection type :literal:`stdin`, :literal:`stdout` or :literal:`stderr`,
        representing the :py:data:`sys.stdin`, :py:data:`sys.stdout` and :py:data:`sys.stderr` streams.
    :type connection: str
    """

    def __init__(self, name: str, connection: str) -> None:
        """Initialise a :py:class:`bsf.connector.ConcurrentProcess` object.

        :param name: A name matching a :py:attr:`bsf.process.Executable.name` attribute of the
            concurrent process to connect to.
        :type name: str
        :param connection: A connection type :literal:`stdin`, :literal:`stdout` or :literal:`stderr`,
            representing the :py:data:`sys.stdin`, :py:data:`sys.stdout` and :py:data:`sys.stderr` streams.
        :type connection: str
        """
        super(ConcurrentProcess, self).__init__()

        self.name = name
        self.connection = connection

        return

    def __repr__(self) -> str:
        return \
            f'{self.__class__.__name__}(' \
            f'name={self.name!r}, ' \
            f'connection={self.connection!r})'


class ElectronicSink(Connector):
    """The :py:class:`bsf.connector.ElectronicSink` class represents a :py:class:`io.FileIO` connection to
    :py:data:`os.devnull` via a Python :py:class:`subprocess.Popen` object and :py:data:`subprocess.DEVNULL`.
    """

    pass


class StandardStream(Connector):
    """The :py:class:`bsf.connector.StandardStream` class represents a standard stream processed via a
    Python :py:class:`threading.Thread` object.

    Standard streams (i.e., :py:data:`sys.stdin`, :py:data:`sys.stdout` and :py:data:`sys.stderr`)
    require processing via a
    Python :py:class:`threading.Thread` object to prevent buffers from filling up and subsequently,
    Python :py:class:`subprocess.Popen` sub-process from blocking.

    :ivar file_path: A file path.
    :type file_path: str | None
    :ivar thread_callable: A :py:class:`typing.Callable` object for the :py:attr:`threading.Thread.target` attribute.
    :type thread_callable: (...) -> Any | None
    :ivar thread_kwargs: A Python :py:class:`dict` object of keyword arguments for the
        :py:attr:`threading.Thread.kwargs` attribute.
    :type thread_kwargs: dict[str, Any] | None
    :ivar thread_joins: A maximum number of attempts calling the :py:meth:`threading.Thread.join` method.
    :type thread_joins: int
    :ivar thread_timeout: A timeout in seconds for calling the :py:meth:`threading.Thread.join` method.
    :type thread_timeout: int
    :ivar thread: A :py:class:`threading.Thread` object.
    :type thread: Thread | None
    """

    def __init__(
            self,
            file_path: Optional[str] = None,
            thread_callable: Optional[Callable] = None,
            thread_kwargs: Optional[dict[str, Any]] = None,
            thread_joins: int = 10,
            thread_timeout: int = 10,
            thread: Optional[Thread] = None) -> None:
        """Initialise a :py:class:`bsf.connector.StandardStream` object.

        :param file_path: A file path.
        :type file_path: str | None
        :param thread_callable: A :py:class:`Callable` object for :py:attr:`threading.Thread.target` attribute.
        :type thread_callable: (...) -> Any | None
        :param thread_kwargs: A Python :py:class:`dict` object of keyword arguments for the
            :py:attr:`threading.Thread.kwargs` attribute.
        :type thread_kwargs: dict[str, Any] | None
        :param thread_joins: A maximum number of attempts calling the :py:meth:`threading.Thread.join` method.
        :type thread_joins: int
        :param thread_timeout: A timeout in seconds for calling the :py:meth:`threading.Thread.join` method.
        :type thread_timeout: int
        :param thread: A :py:class:`threading.Thread` object.
        :type thread: Thread | None
        """
        self.file_path = file_path
        self.thread_callable = thread_callable
        self.thread_kwargs = thread_kwargs
        self.thread_joins = thread_joins
        self.thread_timeout = thread_timeout
        self.thread = thread

        return

    def __repr__(self) -> str:
        return \
            f'{self.__class__.__name__}(' \
            f'file_path={self.file_path!r}, ' \
            f'thread_callable={self.thread_callable!r}, ' \
            f'thread_kwargs={self.thread_kwargs!r}, ' \
            f'thread_joins={self.thread_joins!r}, ' \
            f'thread_timeout={self.thread_timeout!r}, ' \
            f'thread={self.thread!r})'


class StandardInputStream(StandardStream):
    """The :py:class:`bsf.connector.StandardInputStream` class represents a :py:data:`sys.stdin` stream processed
    via a Python :py:class:`threading.Thread` object.
    """

    pass


class StandardOutputStream(StandardStream):
    """The :py:class:`bsf.connector.StandardOutputStream` class represents a :py:data:`sys.stdout` stream processed
    via a Python :py:class:`threading.Thread` object.
    """

    pass


class StandardErrorStream(StandardStream):
    """The :py:class:`bsf.connector.StandardErrorStream` class represents a :py:data:`sys.stderr` stream processed
    via a Python :py:class:`threading.Thread` object.
    """

    pass
