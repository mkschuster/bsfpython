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
"""The :py:mod:`bsf.procedure` module provides classes modelling procedures.
"""
import errno
import os
import pickle
import shutil
from typing import Optional, TypeVar

from bsf.process import RunnableStep
from bsf.standards import Configuration

RunnableType = TypeVar(name='RunnableType', bound='Runnable')


class FilePath(object):
    """The :py:class:`bsf.procedure.FilePath` class represents formalised file path information for the
    :py:class:`bsf.procedure.Runnable` class.

    Each :py:class:`bsf.procedure.Runnable` class is expected to define its corresponding
    :py:class:`bsf.procedure.FilePath` subclass.

    :ivar prefix: A file path prefix.
    :type prefix: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.procedure.FilePath` object.

        :param prefix: A file path prefix.
        :type prefix: str
        """
        self.prefix = prefix

        return


class Runnable(object):
    """The :py:class:`bsf.procedure.Runnable` class represents one or more :py:class:`bsf.process.Executable` objects
    for the :literal:`Runner` script.

    A :py:class:`bsf.procedure.Runnable` holds all information to run one or more :py:class:`bsf.process.Executable`
    objects through the runner script in the :py:attr:`bsf.procedure.Runnable.runner_script` attribute.
    It can be thought of a GNU Bash script that executes as set of :py:class:`bsf.process.RunnableStep` objects
    reflecting commands of a GNU Bash script.

    :cvar runner_script: A global :literal:`Runner` script name.
    :type runner_script: str
    :ivar name: A name.
    :type name: str
    :ivar working_directory: A working directory for writing a :py:class:`pickle.Pickler` file.
    :type working_directory: str
    :ivar code_module: A name of a module, usually under the :py:mod:`bsf.runnables` module
        that implements the logic required to run :py:class:`bsf.process.Executable` objects via the
        :py:attr:`bsf.procedure.Runnable.runner_script`.
    :type code_module: str
    :ivar cache_directory: A cache directory path.
    :type cache_directory: str | None
    :ivar cache_path_dict: A Python :py:class:`dict` object of
        Python :py:class:`str` (name) key and
        Python :py:class:`str` (file_path) value data of files that will be copied into the
        directory path in the :py:attr:`bsf.procedure.Runnable.cache_directory` attribute.
    :type cache_path_dict: dict[str, str]
    """

    runner_script = 'bsf_runner.py'

    def __init__(
            self,
            name: str,
            working_directory: str,
            code_module: str,
            cache_directory: Optional[str] = None,
            cache_path_dict: Optional[dict[str, str]] = None) -> None:
        """Initialise a :py:class:`bsf.procedure.Runnable` object.

        :param name: A name.
        :type name: str
        :param working_directory: A working directory for writing a Python :py:class:`pickle.Pickler` file.
        :type working_directory: str
        :param code_module: A name of a module, usually under the :py:mod:`bsf.runnables` module
            that implements the logic required to run :py:class:`bsf.process.Executable` objects via the
            :py:attr:`bsf.procedure.Runnable.runner_script`.
        :type code_module: str
        :param cache_directory: A cache directory path.
        :type cache_directory: str | None
        :param cache_path_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (name) key and
            Python :py:class:`str` (file_path) value data of files that will be copied into the
            directory path in the :py:attr:`bsf.procedure.Runnable.cache_directory` attribute.
        :type cache_path_dict: dict[str, str] | None
        """

        super(Runnable, self).__init__()

        self.name = name
        self.working_directory = working_directory
        self.code_module = code_module
        self.cache_directory = cache_directory

        if cache_path_dict is None:
            self.cache_path_dict = dict()
        else:
            self.cache_path_dict = cache_path_dict

        return

    def trace(self, level: int = 1) -> list[str]:
        """Trace a :py:class:`bsf.procedure.Runnable` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: list[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        return str_list

    @property
    def pickler_path(self) -> str:
        """Get the Python :py:class:`pickle.Pickler` file path.

        :return: A Python :py:class:`pickle.Pickler` file path.
        :rtype: str
        """
        return os.path.join(self.working_directory, '.'.join((self.name, 'pkl')))

    def to_pickler_path(self) -> None:
        """Write this :py:class:`bsf.procedure.Runnable` as a Python :py:class:`pickle.Pickler` file into the
        working directory.
        """
        with open(file=self.pickler_path, mode='wb') as output_binary_io:
            pickler = pickle.Pickler(file=output_binary_io, protocol=pickle.HIGHEST_PROTOCOL)
            pickler.dump(self)

        return

    @classmethod
    def from_pickler_path(cls, file_path: str) -> RunnableType:
        """Create a :py:class:`bsf.procedure.Runnable` from a Python :py:class:`pickle.Pickler` file via
        Python :py:class:`pickle.Unpickler`.

        :param file_path: A Python :py:class:`pickle.Pickler` file path.
        :type file_path: str
        :return: A :py:class:`bsf.procedure.Runnable` object.
        :rtype: Runnable
        """
        with open(file=file_path, mode='rb') as input_binary_io:
            runnable: Runnable = pickle.Unpickler(input_binary_io).load()

        # Did the Unpickler really return a Runnable object?
        assert isinstance(runnable, Runnable)

        return runnable

    def cache_directory_path(self, absolute: bool = False) -> str:
        """Get the absolute or relative cache directory path of a :py:class:`bsf.procedure.Runnable` object.

        If the :py:attr:`bsf.procedure.Runnable.cache_directory` attribute is not defined,
        :py:attr:`bsf.procedure.Runnable.working_directory` will be prepended.
        Since the relative cache directory path includes the :py:attr:`bsf.procedure.Runnable.name`,
        the directory is :py:class:`bsf.procedure.Runnable`-specific.
        (i.e., :py:attr:`bsf.procedure.Runnable.cache_directory`\\/:py:attr:`bsf.procedure.Runnable.name`\\_cache or
        :py:attr:`bsf.procedure.Runnable.working_directory`\\/:py:attr:`bsf.procedure.Runnable.name`\\_cache)

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: An absolute or relative cache directory path.
        :rtype: str
        """
        directory_name = '_'.join((self.name, 'cache'))

        if absolute:
            if self.cache_directory:
                default_path = self.cache_directory
            else:
                default_path = self.working_directory
            return Configuration.get_absolute_path(
                file_path=directory_name,
                default_path=default_path)
        else:
            return directory_name

    def cache_directory_create(self) -> None:
        """Create a cache directory, if it does not exist and populate it from the cache dictionary.

        In case of an error during the copy, the entire cache directory is removed before an Exception is re-raised.
        """
        if not self.cache_path_dict:
            return

        cache_directory_path = self.cache_directory_path(absolute=True)

        # Create a Runnable-specific cache directory if it does not exist already.

        if not os.path.isdir(cache_directory_path):
            try:
                os.makedirs(cache_directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise exception

        # Copy files from the cache dictionary into the cache directory.

        for key in self.cache_path_dict:
            try:
                shutil.copy2(src=self.cache_path_dict[key], dst=cache_directory_path)
            except OSError as exception:
                shutil.rmtree(path=cache_directory_path, ignore_errors=False)
                raise exception

        return

    def cache_directory_remove(self) -> None:
        """Remove the cache directory and its contents, if it exists.
        """
        cache_directory_path = self.cache_directory_path(absolute=True)

        # Remove the Runnable-specific cache directory and everything within it.

        if os.path.exists(cache_directory_path):
            shutil.rmtree(path=cache_directory_path, ignore_errors=False)

        return

    def get_cache_file_path(self, file_path: str, absolute: bool = False) -> str:
        """Get the absolute or relative cache file path for a file path.

        :param file_path: A default file path.
        :type file_path: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: An absolute or relative cache file path.
        :rtype: str
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        return os.path.join(self.cache_directory_path(absolute=absolute), file_name)

    def temporary_directory_path(self, absolute: bool = False) -> str:
        """Get the absolute or relative temporary directory path of a :py:class:`bsf.procedure.Runnable` object.

        The absolute path prepends the :py:attr:`bsf.procedure.Runnable.working_directory`,
        the relative uses just the :py:attr:`bsf.procedure.Runnable.name` attribute and a :literal:`_temporary` suffix.

        :param absolute: Get an absolute or relative file path.
        :type absolute: bool
        :return: An absolute or relative temporary directory path.
        :rtype: str
        """
        directory_name = '_'.join((self.name, 'temporary'))

        if absolute:
            return os.path.join(self.working_directory, directory_name)
        else:
            return directory_name

    def temporary_directory_create(self) -> None:
        """Create a temporary directory, if it does not exist.
        """
        temporary_directory_path = self.temporary_directory_path(absolute=False)

        if not os.path.isdir(temporary_directory_path):
            try:
                os.makedirs(temporary_directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise exception

        return

    def temporary_directory_remove(self) -> None:
        """Remove the temporary directory and its contents, if it exists.
        """
        temporary_directory_path = self.temporary_directory_path(absolute=False)

        if os.path.exists(temporary_directory_path):
            shutil.rmtree(path=temporary_directory_path, ignore_errors=False)

        return

    def runnable_status_file_path(self, success: bool = True, absolute: bool = False) -> str:
        """Get the status file path for a :py:class:`bsf.procedure.Runnable` object.

        :param success: Get a status file for :emphasis:`successful` or :emphasis:`failed` completion.
        :type success: bool
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A status file path.
        :rtype: str
        """
        if success:
            file_name = '_'.join((self.name, 'completed.txt'))
        else:
            file_name = '_'.join((self.name, 'failed.txt'))

        if absolute:
            return os.path.join(self.working_directory, file_name)
        else:
            return file_name

    def runnable_status_file_create(self, success: bool = True) -> None:
        """Create an empty status file for a :py:class:`bsf.procedure.Runnable` object.

        This method is mainly used by the :py:mod:`bsf.runnable.consecutive` module and related ones.

        :param success: Get a status file for :emphasis:`successful` or :emphasis:`failed` completion.
        :type success: bool
        """
        status_path = self.runnable_status_file_path(success=success)
        open(file=status_path, mode='wt').close()

        return

    def runnable_status_file_remove(self) -> None:
        """Remove the status file for a :py:class:`bsf.procedure.Runnable` object.

        This method is mainly used by the :py:mod:`bsf.runnable.consecutive` module and related ones.
        """
        # Automatically remove both status files, successful or not.

        status_path = self.runnable_status_file_path(success=True)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise exception

        status_path = self.runnable_status_file_path(success=False)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise exception

        return

    def runnable_step_status_file_path(
            self,
            runnable_step: Optional[RunnableStep] = None,
            success: bool = True) -> Optional[str]:
        """Get the status file path for a :py:class:`bsf.process.RunnableStep` object of a
        :py:class:`bsf.procedure.Runnable` object.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        :param success: Get a status file for :emphasis:`successful` or :emphasis:`failed` completion.
        :type success: bool
        :return: A status file path.
        :rtype: str | None
        """
        if runnable_step is None:
            return

        if success:
            return '_'.join((self.name, runnable_step.name, 'completed.txt'))
        else:
            return '_'.join((self.name, runnable_step.name, 'failed.txt'))

    def runnable_step_status_file_create(
            self,
            runnable_step: Optional[RunnableStep] = None,
            success: bool = True) -> None:
        """Create an empty status file for a :py:class:`bsf.process.RunnableStep` object of a
        :py:class:`bsf.procedure.Runnable` object.

        This method is mainly used by the :py:mod:`bsf.runnable.consecutive` module and related ones.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        :param success: Get a status file for :emphasis:`successful` or :emphasis:`failed` completion.
        :type success: bool
        """
        if runnable_step is None:
            return

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=success)
        open(file=status_path, mode='wt').close()

        return

    def runnable_step_status_file_remove(self, runnable_step: Optional[RunnableStep]) -> None:
        """Remove the status file for a :py:class:`bsf.process.RunnableStep` object of a
        :py:class:`bsf.procedure.Runnable` object.

        This method is mainly used by the :py:mod:`bsf.runnable.consecutive` module and related ones.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        # Automatically remove both status files, successful or not.

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=True)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise exception

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=False)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise exception

        return

    def run_consecutively(self, runnable_step_list: list[RunnableStep]) -> Optional[Exception]:
        """Run a Python :py:class:`list` object of :py:class:`bsf.process.RunnableStep` objects through a
         :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable_step_list: A Python :py:class:`list` object of :py:class:`bsf.process.RunnableStep` objects.
        :type runnable_step_list: list[RunnableStep]
        :return: A Python :py:class:`Exception` in case of a :py:class:`bsf.process.RunnableStep` failure.
        :rtype: Exception | None
        """
        # Check the Python list of bsf.process.RunnableStep objects in reverse order to see what has completed already.
        # If a bsf.process.RunnableStep is complete, it will be the first one on the list to become the
        # previous bsf.process.RunnableStep.

        new_runnable_step_list: list[RunnableStep] = list()
        for runnable_step in reversed(runnable_step_list):
            new_runnable_step_list.append(runnable_step)
            if os.path.exists(self.runnable_step_status_file_path(
                    runnable_step=runnable_step,
                    success=True)):
                break
        new_runnable_step_list.reverse()

        # Work through the list of bsf.process.RunnableStep objects in logical order.
        # Keep track of the previous bsf.process.RunnableStep.

        exception_str_list = None
        runnable_step_current = None
        runnable_step_previous = None

        for runnable_step_current in new_runnable_step_list:
            # Check for a bsf.process.RunnableStep-specific status file.
            if os.path.exists(self.runnable_step_status_file_path(
                    runnable_step=runnable_step_current,
                    success=True)):
                # If a status file exists, this RunnableStep is complete.
                # Set it as the previous RunnableStep and continue with the next RunnableStep.
                runnable_step_previous = runnable_step_current
                continue

            # Do the work.

            exception_str_list = runnable_step_current.run()

            # Upon failure, break out of this loop without deleting obsolete files or altering status files
            # at this point.

            if exception_str_list:
                break

            # Delete the list of file paths that the current bsf.process.RunnableStep declared to be obsolete now.

            runnable_step_current.remove_obsolete_file_paths()

            # Create an empty status file upon success.

            self.runnable_step_status_file_create(runnable_step=runnable_step_current, success=True)

            # Remove the status file of the previous bsf.process.RunnableStep, if it has been defined at this stage.

            if runnable_step_previous is not None:
                self.runnable_step_status_file_remove(runnable_step=runnable_step_previous)

            # Finally, make the current RunnableStep the previous RunnableStep.

            runnable_step_previous = runnable_step_current

        if not exception_str_list:
            # Upon success, remove the status file of the previous RunnableStep, if it has been defined at this stage.

            if runnable_step_previous is not None:
                self.runnable_step_status_file_remove(runnable_step=runnable_step_previous)
        else:
            # Upon failure, create a RunnableStep-specific status file showing failure and return an Exception object.

            if runnable_step_current is not None:
                self.runnable_step_status_file_create(runnable_step=runnable_step_current, success=False)

            return Exception('\n'.join(exception_str_list))

        return None


class ConsecutiveRunnable(Runnable):
    """The :py:class:`bsf.procedure.ConsecutiveRunnable` class represents a procedure of consecutively running
    :py:class:`bsf.process.RunnableStep` or subclasses thereof.

    :ivar runnable_step_list: A Python :py:class:`list` object of :py:class:´bsf.process.RunnableStep` objects.
    :type runnable_step_list: list[RunnableStep]
    """

    def __init__(
            self,
            name: str,
            working_directory: str,
            code_module: str = 'bsf.runnables.consecutive',
            cache_directory: Optional[str] = None,
            cache_path_dict: Optional[dict[str, str]] = None,
            runnable_step_list: Optional[list[RunnableStep]] = None) -> None:
        """Initialise a :py:class:`bsf.ConsecutiveRunnable` object.

        :param name: A name.
        :type name: str
        :param working_directory: A working directory for writing a Python :py:class:`pickle.Pickler` file.
        :type working_directory: str
        :param code_module: A name of a module, usually under the :py:mod:`bsf.runnables` module
            that implements the logic required to run :py:class:`bsf.process.Executable` objects via the
            :py:attr:`bsf.procedure.Runnable.runner_script` consecutively
        :type code_module: str
        :param cache_directory: A cache directory path.
        :type cache_directory: str | None
        :param cache_path_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (name) key and
          Python :py:class:`str` (file_path) value data of files that will be copied into the
          directory path in the :py:attr:`bsf.procedure.Runnable.cache_directory` attribute.
        :type cache_path_dict: dict[str, str] | None
        :param runnable_step_list: A Python :py:class:`list` object of :py:class:`bsf.process.RunnableStep` objects.
        :type runnable_step_list: list[RunnableStep] | None
        """
        super(ConsecutiveRunnable, self).__init__(
            name=name,
            working_directory=working_directory,
            code_module=code_module,
            cache_directory=cache_directory,
            cache_path_dict=cache_path_dict)

        if runnable_step_list is None:
            self.runnable_step_list = list()
        else:
            self.runnable_step_list = runnable_step_list

        return

    def trace(self, level: int = 1) -> list[str]:
        """Trace a :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: list[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  runnable_step_list: {!r}\n'.format(indent, self.runnable_step_list))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        str_list.append('{}  Python list of RunnableStep objects:\n'.format(indent))
        for runnable_step in self.runnable_step_list:
            str_list.extend(runnable_step.trace(level=level + 1))

        return str_list

    def add_runnable_step(self, runnable_step: Optional[RunnableStep]) -> None:
        """Convenience method to facilitate initialising, adding and returning a
        :py:class:`bsf.process.RunnableStep` object.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list.append(runnable_step)

        return

    def run(self) -> Optional[Exception]:
        """Convenience method to run a :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :return: A Python :py:class:`Exception` in case of a :py:class:`bsf.process.RunnableStep` failure.
        :rtype: Exception | None
        """
        return self.run_consecutively(runnable_step_list=self.runnable_step_list)


class ConcurrentRunnable(Runnable):
    """The :py:class:`bsf.procedure.ConcurrentRunnable` class represents a procedure of concurrently running
    :py:class:`bsf.process.RunnableStep` objects or subclasses thereof.

    :ivar runnable_step_list_prologue: A Python :py:class:`list` object of
        :py:class:`bsf.process.RunnableStep` objects run as prologue.
    :type runnable_step_list_prologue: list[RunnableStep]
    :ivar runnable_step_list_concurrent: A Python :py:class:`list` object of
        :py:class:`bsf.process.RunnableStep` objects run concurrently.
    :type runnable_step_list_concurrent: list[RunnableStep]
    :ivar runnable_step_list_epilogue: A Python :py:class:`list` object of
        :py:class:`bsf.process.RunnableStep` objects run as epilogue.
    :type runnable_step_list_epilogue: list[RunnableStep]
    """

    def __init__(
            self,
            name: str,
            working_directory: str,
            code_module: str = 'bsf.runnables.concurrent',
            cache_directory: Optional[str] = None,
            cache_path_dict: Optional[dict[str, str]] = None,
            runnable_step_list_prologue: Optional[list[RunnableStep]] = None,
            runnable_step_list_concurrent: Optional[list[RunnableStep]] = None,
            runnable_step_list_epilogue: Optional[list[RunnableStep]] = None) -> None:
        """Initialise a :py:class:`bsf.ConcurrentRunnable` object.

        :param name: A name.
        :type name: str
        :param working_directory: A working directory for writing a Python :py:class:`pickle.Pickler` file.
        :type working_directory: str
        :param code_module: A name of a module, usually under the :py:mod:`bsf.runnables` module
            that implements the logic required to run :py:class:`bsf.process.Executable` objects via the
            :py:attr:`bsf.procedure.Runnable.runner_script` concurrently.
        :type code_module: str
        :param cache_directory: A cache directory path.
        :type cache_directory: str | None
        :param cache_path_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (name) key and
            Python :py:class:`str` (file_path) value data of files that will be copied into the
            directory path in the :py:attr:`bsf.procedure.Runnable.cache_directory` attribute.
        :type cache_path_dict: dict[str, str] | None
        :param runnable_step_list_prologue: A Python :py:class:`list` object of
            :py:class:´bsf.process.RunnableStep` objects run as prologue.
        :type runnable_step_list_prologue: list[RunnableStep] | None
        :param runnable_step_list_concurrent: A Python :py:class:`list` object of
            :py:class:`bsf.process.RunnableStep` objects run concurrently.
        :type runnable_step_list_concurrent: list[RunnableStep] | None
        :param runnable_step_list_epilogue: A Python :py:class:`list` object of
            :py:class:`bsf.process.RunnableStep` objects run as epilogue.
        :type runnable_step_list_epilogue: list[RunnableStep] | None
        """
        super(ConcurrentRunnable, self).__init__(
            name=name,
            working_directory=working_directory,
            code_module=code_module,
            cache_directory=cache_directory,
            cache_path_dict=cache_path_dict)

        if runnable_step_list_prologue is None:
            self.runnable_step_list_prologue = list()
        else:
            self.runnable_step_list_prologue = runnable_step_list_prologue

        if runnable_step_list_concurrent is None:
            self.runnable_step_list_concurrent = list()
        else:
            self.runnable_step_list_concurrent = runnable_step_list_concurrent

        if runnable_step_list_epilogue is None:
            self.runnable_step_list_epilogue = list()
        else:
            self.runnable_step_list_epilogue = runnable_step_list_epilogue

        return

    def trace(self, level: int = 1) -> list[str]:
        """Trace a :py:class:`bsf.procedure.ConcurrentRunnable` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: list[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  runnable_step_list_prologue: {!r}\n'.format(indent, self.runnable_step_list_prologue))
        str_list.append('{}  runnable_step_list_concurrent: {!r}\n'.format(indent, self.runnable_step_list_concurrent))
        str_list.append('{}  runnable_step_list_epilogue: {!r}\n'.format(indent, self.runnable_step_list_epilogue))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        str_list.append('{}  Python list of pre RunnableStep objects:\n'.format(indent))
        for runnable_step in self.runnable_step_list_prologue:
            str_list.extend(runnable_step.trace(level=level + 1))

        str_list.append('{}  Python list of concurrent RunnableStep objects:\n'.format(indent))
        for runnable_step in self.runnable_step_list_concurrent:
            str_list.extend(runnable_step.trace(level=level + 1))

        str_list.append('{}  Python list of postRunnableStep objects:\n'.format(indent))
        for runnable_step in self.runnable_step_list_epilogue:
            str_list.extend(runnable_step.trace(level=level + 1))

        return str_list

    def add_runnable_step_prologue(self, runnable_step: Optional[RunnableStep]) -> None:
        """Convenience method to add a :py:class:`bsf.process.RunnableStep` object to the prologue list.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list_prologue.append(runnable_step)

        return

    def add_runnable_step_epilogue(self, runnable_step: Optional[RunnableStep]) -> None:
        """Convenience method to add a :py:class:`bsf.process.RunnableStep` object to the epilogue list.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list_epilogue.append(runnable_step)

        return

    def add_runnable_step(self, runnable_step: Optional[RunnableStep]) -> None:
        """Convenience method to add a :py:class:`bsf.process.RunnableStep` object to the concurrent list.

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list_concurrent.append(runnable_step)

        return
