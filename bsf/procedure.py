# -*- coding: utf-8 -*-
"""Procedure module.

A package of classes and methods modelling procedures.
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
import pickle
import shutil

from bsf.process import RunnableStep
from bsf.standards import Configuration


class FilePath(object):
    """The C{bsf.procedure.FilePath} class represents formalised file path information for the
    C{bsf.procedure.Runnable} class.

    Each C{bsf.procedure.Runnable} class is expected to define its corresponding C{bsf.procedure.FilePath} sub-class.
    Attributes:
    @ivar prefix: File path prefix
    @type prefix: str
    #ivar temporary_directory: Temporary directory path
    #type temporary_directory: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.procedure.FilePath}.

        @param prefix: File path prefix
        @type prefix: str
        """
        self.prefix = prefix

        return


class Runnable(object):
    """The C{bsf.procedure.Runnable} class represents one or more C{bsf.process.Executable} objects
    for the I{Runner} script.

    A C{bsf.procedure.Runnable} holds all information to run one or more C{bsf.process.Executable} objects through the
    C{bsf.procedure.Runnable.runner_script}. It can be thought of a GNU Bash script that executes as set of
    C{bsf.process.RunnableStep} objects reflecting commands of a GNU Bash script.

    Attributes:
    @cvar runner_script: Name of the I{Runner} script
    @type runner_script: str
    @ivar name: Name
    @type name: str
    @ivar working_directory: Working directory to write C{pickle.Pickler} files
    @type working_directory: str
    @ivar code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
        C{bsf.process.Executable} objects via the C{bsf.procedure.Runnable.runner_script}.
    @type code_module: str
    @ivar cache_directory: Cache directory
    @type cache_directory: str | None
    @ivar cache_path_dict: Python C{dict} of Python C{str} (name) key and
        Python C{str} (file_path) value data of files that will be copied into the
        C{bsf.procedure.Runnable.cache_directory}
    @type cache_path_dict: dict[str, str]
    @ivar debug: Debug level
    @type debug: int
    """

    runner_script = 'bsf_runner.py'

    def __init__(
            self,
            name,
            working_directory,
            code_module,
            cache_directory=None,
            cache_path_dict=None,
            debug=0):
        """Initialise a C{bsf.procedure.Runnable}.

        @param name: Name
        @type name: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.Executable} objects via the C{bsf.procedure.Runnable.runner_script}
        @type code_module: str
        @param cache_directory: Cache directory
        @type cache_directory: str | None
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the
            C{bsf.procedure.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str] | None
        @param debug: Integer debugging level
        @type debug: int
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

        if debug is None:
            self.debug = 0
        else:
            assert isinstance(debug, int)
            self.debug = debug

        return

    def trace(self, level=1):
        """Trace a C{bsf.procedure.Runnable}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  debug: {!r}\n'.format(indent, self.debug))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        return str_list

    @property
    def pickler_path(self):
        """Get the Python C{pickle.Pickler} file path.

        @return: Python C{pickle.Pickler} file path
        @rtype: str
        """
        return os.path.join(self.working_directory, '.'.join((self.name, 'pkl')))

    def to_pickler_path(self):
        """Write this C{bsf.procedure.Runnable} as a Python C{pickle.Pickler} file into the working directory.
        """
        with open(file=self.pickler_path, mode='wb') as output_file:
            pickler = pickle.Pickler(file=output_file, protocol=pickle.HIGHEST_PROTOCOL)
            pickler.dump(self)

        return

    @classmethod
    def from_pickler_path(cls, file_path):
        """Create a C{bsf.procedure.Runnable} from a Python C{pickle.Pickler} file via Python C{pickle.Unpickler}.

        @param file_path: File path to a Python C{pickle.Pickler} file
        @type file_path: str
        @return: C{bsf.procedure.Runnable}
        @rtype: Runnable
        """
        with open(file=file_path, mode='rb') as input_file:
            runnable = pickle.Unpickler(input_file).load()
            """ @type runnable: Runnable """

        # Did the Unpickler really return a Runnable object?
        assert isinstance(runnable, Runnable)

        return runnable

    def cache_directory_path(self, absolute=False):
        """Get the absolute or relative cache directory path of a C{bsf.procedure.Runnable}.

        If C{bsf.procedure.Runnable.cache_directory} is not defined, C{bsf.procedure.Runnable.working_directory}
        will be prepended.
        Since the relative cache directory path includes the C{bsf.procedure.Runnable.name},
        the directory is C{bsf.procedure.Runnable}-specific.
        (i.e. C{bsf.procedure.Runnable.cache_directory}/C{bsf.procedure.Runnable.name}_cache or
        C{bsf.procedure.Runnable.working_directory}/C{bsf.procedure.Runnable.name}_cache)

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Absolute or relative cache directory path
        @rtype: str
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

    def cache_directory_create(self):
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
                    raise

        # Copy files from the cache dictionary into the cache directory.

        for key in self.cache_path_dict:
            try:
                shutil.copy2(src=self.cache_path_dict[key], dst=cache_directory_path)
            except OSError:
                shutil.rmtree(path=cache_directory_path, ignore_errors=False)
                raise

        return

    def cache_directory_remove(self):
        """Remove the cache directory and its contents, if it exists.
        """
        cache_directory_path = self.cache_directory_path(absolute=True)

        # Remove the Runnable-specific cache directory and everything within it.

        if os.path.exists(cache_directory_path):
            shutil.rmtree(path=cache_directory_path, ignore_errors=False)

        return

    def get_cache_file_path(self, file_path, absolute=False):
        """Get the absolute or relative cache file path for a file path.

        @param file_path: Default file path
        @type file_path: str
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Absolute or relative cache file path
        @rtype: str
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        return os.path.join(self.cache_directory_path(absolute=absolute), file_name)

    def temporary_directory_path(self, absolute=False):
        """Get the absolute or relative temporary directory path of a C{bsf.procedure.Runnable}.

        The absolute path prepends the C{bsf.procedure.Runnable.working_directory},
        the relative just C{bsf.procedure.Runnable.name}_temporary.
        @param absolute: Absolute or relative file path
        @type absolute: bool
        @return: Absolute or relative temporary directory path
        @rtype: str
        """
        directory_name = '_'.join((self.name, 'temporary'))

        if absolute:
            return os.path.join(self.working_directory, directory_name)
        else:
            return directory_name

    def temporary_directory_create(self):
        """Create a temporary directory, if it does not exist.
        """
        temporary_directory_path = self.temporary_directory_path(absolute=False)

        if not os.path.isdir(temporary_directory_path):
            try:
                os.makedirs(temporary_directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return

    def temporary_directory_remove(self):
        """Remove the temporary directory and its contents, if it exists.
        """
        temporary_directory_path = self.temporary_directory_path(absolute=False)

        if os.path.exists(temporary_directory_path):
            shutil.rmtree(path=temporary_directory_path, ignore_errors=False)

        return

    def runnable_status_file_path(self, success=True, absolute=False):
        """Get the status file path for a C{bsf.procedure.Runnable}.

        @param success: Successful completion
        @type success: bool
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Status file path
        @rtype: str
        """
        if success:
            file_name = '_'.join((self.name, 'completed.txt'))
        else:
            file_name = '_'.join((self.name, 'failed.txt'))

        if absolute:
            return os.path.join(self.working_directory, file_name)
        else:
            return file_name

    def runnable_status_file_create(self, success=True):
        """Create an empty status file for a C{bsf.procedure.Runnable}.

        This method is mainly used by C{bsf.runnable.consecutive} and related modules.

        @param success: Successful completion
        @type success: bool
        """
        status_path = self.runnable_status_file_path(success=success)
        open(file=status_path, mode='wt').close()

        return

    def runnable_status_file_remove(self):
        """Remove the status file for a C{bsf.procedure.Runnable}.

        This method is mainly used by C{bsf.runnable.consecutive} and related modules.
        """
        # Automatically remove both status files, successful or not.

        status_path = self.runnable_status_file_path(success=True)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        status_path = self.runnable_status_file_path(success=False)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        return

    def runnable_step_status_file_path(self, runnable_step=None, success=True):
        """Get the status file path for a C{bsf.process.RunnableStep} of a C{bsf.procedure.Runnable}.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep | None
        @param success: Successful completion
        @type success: bool
        @return: Status file path
        @rtype: str | None
        """
        if runnable_step is None:
            return

        if success:
            return '_'.join((self.name, runnable_step.name, 'completed.txt'))
        else:
            return '_'.join((self.name, runnable_step.name, 'failed.txt'))

    def runnable_step_status_file_create(self, runnable_step=None, success=True):
        """Create an empty status file for a C{bsf.process.RunnableStep} of a C{bsf.procedure.Runnable}.

        This method is mainly used by C{bsf.runnable.consecutive} and related modules.

        @param runnable_step: C{bsf.process.RunnableStep} | None
        @type runnable_step: RunnableStep
        @param success: Successful completion
        @type success: bool
        """
        if runnable_step is None:
            return

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=success)
        open(file=status_path, mode='wt').close()

        return

    def runnable_step_status_file_remove(self, runnable_step):
        """Remove the status file for a C{bsf.process.RunnableStep} of a C{bsf.procedure.Runnable}.

        This method is mainly used by C{bsf.runnable.consecutive} and related modules.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        # Automatically remove both status files, successful or not.

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=True)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=False)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        return

    def run_consecutively(self, runnable_step_list):
        """Run a Python C{list} of C{bsf.process.RunnableStep} objects through a C{bsf.procedure.ConsecutiveRunnable}.

        @param runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
        @type runnable_step_list: list[RunnableStep]
        @return: Exception in case of a RunnableStep failure
        @rtype: Exception | None
        """
        # Check the Python list of bsf.process.RunnableStep objects in reverse order to see what has completed already.
        # If a bsf.process.RunnableStep is complete, it will be the first one on the list to become the
        # previous bsf.process.RunnableStep.

        new_runnable_step_list = list()
        """ @type new_runnable_step_list: list[RunnableStep] """
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

            Exception('\n'.join(exception_str_list))

        return None


class ConsecutiveRunnable(Runnable):
    """The C{bsf.procedure.ConsecutiveRunnable} represents a procedure of consecutively running
    C{bsf.process.RunnableStep} or sub-classes thereof.

    Attributes:
    @ivar runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
    @type runnable_step_list: list[RunnableStep]
    """

    def __init__(
            self,
            name,
            working_directory,
            code_module='bsf.runnables.consecutive',
            cache_directory=None,
            cache_path_dict=None,
            debug=0,
            runnable_step_list=None):
        """Initialise a C{bsf.ConsecutiveRunnable}.

        @param name: Name
        @type name: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.RunnableStep} objects via the C{bsf.procedure.Runnable.runner_script} consecutively
        @type code_module: str
        @param cache_directory: Cache directory
        @type cache_directory: str | None
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the
            C{bsf.procedure.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str] | None
        @param debug: Integer debugging level
        @type debug: int
        @param runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
        @type runnable_step_list: list[RunnableStep]
        """
        super(ConsecutiveRunnable, self).__init__(
            name=name,
            working_directory=working_directory,
            code_module=code_module,
            cache_directory=cache_directory,
            cache_path_dict=cache_path_dict,
            debug=debug)

        if runnable_step_list is None:
            self.runnable_step_list = list()
        else:
            self.runnable_step_list = runnable_step_list

        return

    def trace(self, level=1):
        """Trace a C{bsf.procedure.ConsecutiveRunnable}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  runnable_step_list: {!r}\n'.format(indent, self.runnable_step_list))
        str_list.append('{}  debug: {!r}\n'.format(indent, self.debug))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        str_list.append('{}  Python list of RunnableStep objects:\n'.format(indent))
        for runnable_step in self.runnable_step_list:
            str_list.extend(runnable_step.trace(level=level + 1))

        return str_list

    def add_runnable_step(self, runnable_step):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.process.RunnableStep}.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list.append(runnable_step)

        return

    def run(self):
        """Convenience function to run a C{bsf.procedure.ConsecutiveRunnable}.

        @return: Exception in case of a RunnableStep failure
        @rtype: Exception | None
        """
        return self.run_consecutively(runnable_step_list=self.runnable_step_list)


class ConcurrentRunnable(Runnable):
    """The C{bsf.procedure.ConcurrentRunnable} represents a procedure of concurrently running
    C{bsf.process.RunnableStep} or sub-classes thereof.

    Attributes:
    @ivar runnable_step_list_prologue: Python C{list} of C{bsf.process.RunnableStep} objects run as prologue
    @type runnable_step_list_prologue: list[RunnableStep] | None
    @ivar runnable_step_list_concurrent: Python C{list} of C{bsf.process.RunnableStep} objects run concurrently
    @type runnable_step_list_concurrent: list[RunnableStep] | None
    @ivar runnable_step_list_epilogue: Python C{list} of C{bsf.process.RunnableStep} objects run as epilogue
    @type runnable_step_list_epilogue: list[RunnableStep] | None
    """

    def __init__(
            self,
            name,
            working_directory,
            code_module='bsf.runnables.concurrent',
            cache_directory=None,
            cache_path_dict=None,
            debug=0,
            runnable_step_list_prologue=None,
            runnable_step_list_concurrent=None,
            runnable_step_list_epilogue=None):
        """Initialise a C{bsf.ConcurrentRunnable}.

        @param name: Name
        @type name: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.RunnableStep} objects via the C{bsf.procedure.Runnable.runner_script} concurrently
        @type code_module: str
        @param cache_directory: Cache directory
        @type cache_directory: str | None
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the
            C{bsf.procedure.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str] | None
        @param debug: Integer debugging level
        @type debug: int
        @param runnable_step_list_prologue: Python C{list} of C{bsf.process.RunnableStep} objects run as prologue
        @type runnable_step_list_prologue: list[RunnableStep] | None
        @param runnable_step_list_concurrent: Python C{list} of C{bsf.process.RunnableStep} objects run concurrently
        @type runnable_step_list_concurrent: list[RunnableStep] | None
        @param runnable_step_list_epilogue: Python C{list} of C{bsf.process.RunnableStep} objects run as epilogue
        @type runnable_step_list_epilogue: list[RunnableStep] | None
        """
        super(ConcurrentRunnable, self).__init__(
            name=name,
            working_directory=working_directory,
            code_module=code_module,
            cache_directory=cache_directory,
            cache_path_dict=cache_path_dict,
            debug=debug)

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

    def trace(self, level=1):
        """Trace a C{bsf.procedure.ConcurrentRunnable}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  runnable_step_list_prologue: {!r}\n'.format(indent, self.runnable_step_list_prologue))
        str_list.append('{}  runnable_step_list_concurrent: {!r}\n'.format(indent, self.runnable_step_list_concurrent))
        str_list.append('{}  runnable_step_list_epilogue: {!r}\n'.format(indent, self.runnable_step_list_epilogue))
        str_list.append('{}  debug: {!r}\n'.format(indent, self.debug))

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

    def add_runnable_step_prologue(self, runnable_step):
        """Convenience method to add a C{bsf.process.RunnableStep} to the prologue list.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list_prologue.append(runnable_step)

        return

    def add_runnable_step_epilogue(self, runnable_step):
        """Convenience method to add a C{bsf.process.RunnableStep} to the epilogue list.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list_epilogue.append(runnable_step)

        return

    def add_runnable_step(self, runnable_step):
        """Convenience method to add a C{bsf.process.RunnableStep} to the concurrent list.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep | None
        """
        if runnable_step is None:
            return

        self.runnable_step_list_concurrent.append(runnable_step)

        return
