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
import pickle
import shutil

import bsf.process
import bsf.standards


class FilePath(object):
    """The C{bsf.procedure.FilePath} class represents formalised file path information for the
    C{bsf.procedure.Runnable} class.

    Each C{bsf.procedure.Runnable} class is expected to define its corresponding C{bsf.procedure.FilePath} sub-class.
    Attributes:
    @ivar prefix: File path prefix
    @type prefix: str | unicode
    #ivar temporary_directory: Temporary directory path
    #type temporary_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.procedure.FilePath}.

        @param prefix: File path prefix
        @type prefix: str | unicode
        @return:
        @rtype:
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
    @type runner_script: str | unicode
    @ivar name: Name
    @type name: str
    @ivar code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
        C{bsf.process.Executable} objects via the C{bsf.procedure.Runnable.runner_script}.
    @type code_module: str
    @ivar cache_directory: Cache directory
    @type cache_directory: str | unicode | None
    @ivar cache_path_dict: Python C{dict} of Python C{str} (name) key and
        Python C{str} (file_path) value data of files that will be copied into the
        C{bsf.procedure.Runnable.cache_directory}
    @type cache_path_dict: dict[str, str | unicode]
    @ivar file_path_object: C{bsf.procedure.FilePath}
    @type file_path_object: bsf.procedure.FilePath
    @ivar working_directory: Working directory to write C{pickle.Pickler} files
    @type working_directory: str | unicode | None
    @ivar debug: Debug level
    @type debug: int
    """

    runner_script = 'bsf_runner.py'

    def __init__(
            self,
            name,
            code_module,
            working_directory,
            cache_directory=None,
            cache_path_dict=None,
            file_path_object=None,
            debug=0):
        """Initialise a C{bsf.procedure.Runnable}.

        @param name: Name
        @type name: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.Executable} objects via the C{bsf.procedure.Runnable.runner_script}
        @type code_module: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str | unicode
        @param cache_directory: Cache directory
        @type cache_directory: str | unicode | None
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the
            C{bsf.procedure.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str | unicode]
        @param file_path_object: C{bsf.procedure.FilePath}
        @type file_path_object: bsf.procedure.FilePath
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        super(Runnable, self).__init__()

        self.name = name  # Can be None.
        self.code_module = code_module  # Can be None.
        self.working_directory = working_directory
        self.cache_directory = cache_directory

        if cache_path_dict is None:
            self.cache_path_dict = dict()
        else:
            self.cache_path_dict = cache_path_dict

        if file_path_object is None:
            self.file_path_object = FilePath(prefix='default_file_path')
        else:
            self.file_path_object = file_path_object

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
        @rtype: list[str | unicode]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  file_path_object: {!r}\n'.format(indent, self.file_path_object))
        str_list.append('{}  debug: {!r}\n'.format(indent, self.debug))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        return str_list

    @property
    def pickler_path(self):
        """Get the Python C{pickle.Pickler} file path.

        @return: Python C{pickle.Pickler} file path
        @rtype: str | unicode
        """
        return os.path.join(self.working_directory, '.'.join((self.name, 'pkl')))

    def to_pickler_path(self):
        """Write this C{bsf.procedure.Runnable} as a Python C{pickle.Pickler} file into the working directory.

        @return:
        @rtype:
        """
        with open(self.pickler_path, 'wb') as output_file:
            pickler = pickle.Pickler(file=output_file, protocol=pickle.HIGHEST_PROTOCOL)
            pickler.dump(obj=self)

        return

    @classmethod
    def from_pickler_path(cls, file_path):
        """Create a C{bsf.procedure.Runnable} from a Python C{pickle.Pickler} file via Python C{pickle.Unpickler}.

        @param file_path: File path to a Python C{pickle.Pickler} file
        @type file_path: str | unicode
        @return: C{bsf.procedure.Runnable}
        @rtype: bsf.procedure.Runnable
        """
        with open(file_path, 'rb') as input_file:
            runnable = pickle.Unpickler(file=input_file).load()
            """ @type runnable: bsf.procedure.Runnable """

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
        @rtype: str | unicode
        """
        directory_name = '_'.join((self.name, 'cache'))

        if absolute:
            if self.cache_directory:
                default_path = self.cache_directory
            else:
                default_path = self.working_directory
            return bsf.standards.Configuration.get_absolute_path(
                file_path=directory_name,
                default_path=default_path)
        else:
            return directory_name

    def cache_directory_create(self):
        """Create a cache directory, if it does not exist and populate it from the cache dictionary.

        In case of an error during the copy, the entire cache directory is removed before an Exception is re-raised.
        @return:
        @rtype:
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

        @return:
        @rtype:
        """
        cache_directory_path = self.cache_directory_path(absolute=True)

        # Remove the Runnable-specific cache directory and everything within it.

        if os.path.exists(cache_directory_path):
            shutil.rmtree(path=cache_directory_path, ignore_errors=False)

        return

    def get_cache_file_path(self, file_path, absolute=False):
        """Get the absolute or relative cache file path for a file path.

        @param file_path: Default file path
        @type file_path: str | unicode
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Absolute or relative cache file path
        @rtype: str | unicode
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
        @rtype: str | unicode
        """
        directory_name = '_'.join((self.name, 'temporary'))

        if absolute:
            return os.path.join(self.working_directory, directory_name)
        else:
            return directory_name

    def temporary_directory_create(self):
        """Create a temporary directory, if it does not exist.

        @return:
        @rtype:
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

        @return:
        @rtype:
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
        @rtype: str | unicode
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

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @param success: Successful completion
        @type success: bool
        @return:
        @rtype:
        """
        status_path = self.runnable_status_file_path(success=success)
        open(status_path, 'wt').close()

        return

    def runnable_status_file_remove(self):
        """Remove the status file for a C{bsf.procedure.Runnable}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @return:
        @rtype:
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
        @type runnable_step: bsf.process.RunnableStep | None
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

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @param runnable_step: C{bsf.process.RunnableStep} | None
        @type runnable_step: bsf.process.RunnableStep
        @param success: Successful completion
        @type success: bool
        @return:
        @rtype:
        """
        if runnable_step is None:
            return

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=success)
        open(status_path, 'wt').close()

        return

    def runnable_step_status_file_remove(self, runnable_step):
        """Remove the status file for a C{bsf.process.RunnableStep} of a C{bsf.procedure.Runnable}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep | None
        @return:
        @rtype:
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
        @type runnable_step_list: list[bsf.process.RunnableStep]
        @return: Exception in case of a RunnableStep failure
        @rtype: Exception | None
        """
        # Check the Python list of bsf.process.RunnableStep objects in reverse order to see what has completed already.
        # If a bsf.process.RunnableStep is complete, it will be the first one on the list to become the
        # previous bsf.process.RunnableStep.

        new_runnable_step_list = list()
        """ @type new_runnable_step_list: list[bsf.process.RunnableStep] """
        for runnable_step in reversed(runnable_step_list):
            new_runnable_step_list.append(runnable_step)
            if os.path.exists(path=self.runnable_step_status_file_path(
                    runnable_step=runnable_step,
                    success=True)):
                break
        new_runnable_step_list.reverse()

        # Work through the list of bsf.process.RunnableStep objects in logical order.
        # Keep track of the previous bsf.process.RunnableStep.

        child_return_code = 0
        runnable_step_current = None
        runnable_step_previous = None

        for runnable_step_current in new_runnable_step_list:
            # Check for a bsf.process.RunnableStep-specific status file.
            if os.path.exists(path=self.runnable_step_status_file_path(
                    runnable_step=runnable_step_current,
                    success=True)):
                # If a status file exists, this RunnableStep is complete.
                # Set it as the previous RunnableStep and continue with the next RunnableStep.
                runnable_step_previous = runnable_step_current
                continue

            # Do the work.

            child_return_code = runnable_step_current.run()

            # Upon failure, break out of this loop without deleting obsolete files or altering status files
            # at this point.

            if child_return_code != 0:
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

        if child_return_code == 0:
            # Upon success, remove the status file of the previous RunnableStep, if it has been defined at this stage.

            if runnable_step_previous is not None:
                self.runnable_step_status_file_remove(runnable_step=runnable_step_previous)
        else:
            # Upon failure, create a RunnableStep-specific status file showing failure and return an Exception object.

            if runnable_step_current is not None:
                self.runnable_step_status_file_create(runnable_step=runnable_step_current, success=False)

            if child_return_code > 0:
                return Exception(
                    bsf.process.get_timestamp() +
                    ' Child process ' + repr(self.name) + ' ' + repr(runnable_step_current.name) +
                    ' failed with return code ' +
                    repr(+child_return_code) + '.')
            elif child_return_code < 0:
                return Exception(
                    bsf.process.get_timestamp() +
                    ' Child process ' + repr(self.name) + ' ' + repr(runnable_step_current.name) +
                    ' received signal ' +
                    repr(-child_return_code) + '.')

        return None


class ConsecutiveRunnable(Runnable):
    """The C{bsf.procedure.ConsecutiveRunnable} represents a procedure of consecutively running
    C{bsf.process.RunnableStep} or sub-classes thereof.

    Attributes:
    @ivar runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
    @type runnable_step_list: list[bsf.process.RunnableStep]
    """

    def __init__(
            self,
            name,
            code_module,
            working_directory,
            cache_directory=None,
            cache_path_dict=None,
            file_path_object=None,
            debug=0,
            runnable_step_list=None):
        """Initialise a C{bsf.ConsecutiveRunnable}.

        @param name: Name
        @type name: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.RunnableStep} objects via the C{bsf.procedure.Runnable.runner_script} consecutively
        @type code_module: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str | unicode
        @param cache_directory: Cache directory
        @type cache_directory: str | unicode | None
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the
            C{bsf.procedure.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str | unicode]
        @param file_path_object: C{bsf.procedure.FilePath}
        @type file_path_object: bsf.procedure.FilePath
        @param debug: Integer debugging level
        @type debug: int
        @param runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
        @type runnable_step_list: list[bsf.process.RunnableStep]
        @return:
        @rtype:
        """
        super(ConsecutiveRunnable, self).__init__(
            name=name,
            code_module=code_module,
            working_directory=working_directory,
            cache_directory=cache_directory,
            cache_path_dict=cache_path_dict,
            file_path_object=file_path_object,
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
        @rtype: list[str | unicode]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name: {!r}\n'.format(indent, self.name))
        str_list.append('{}  code_module: {!r}\n'.format(indent, self.code_module))
        str_list.append('{}  working_directory: {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict))
        str_list.append('{}  file_path_object: {!r}\n'.format(indent, self.file_path_object))
        str_list.append('{}  runnable_step_list: {!r}\n'.format(indent, self.runnable_step_list))
        str_list.append('{}  debug: {!r}\n'.format(indent, self.debug))

        str_list.append('{}  Python dict of Python str (cache path) objects:\n'.format(indent))
        for key in sorted(self.cache_path_dict):
            str_list.append('{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key]))

        str_list.append('{}  Python list of RunnableStep objects:\n'.format(indent))
        for runnable_step in self.runnable_step_list:
            str_list.extend(runnable_step.trace(level=level + 1))

        return str_list

    def add_runnable_step(self, runnable_step=None):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.process.RunnableStep}.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep | None
        @return: C{bsf.process.RunnableStep}
        @rtype: bsf.process.RunnableStep
        """
        if runnable_step is None:
            return

        assert isinstance(runnable_step, bsf.process.RunnableStep)

        self.runnable_step_list.append(runnable_step)

        return runnable_step

    def run(self):
        """Convenience function to run a C{bsf.procedure.ConsecutiveRunnable}.

        @return: Exception in case of a RunnableStep failure
        @rtype: Exception | None
        """
        return self.run_consecutively(runnable_step_list=self.runnable_step_list)


class Connector(object):
    """The C{Connector} class represents an abstract super-class of inter-process connections.

    Attributes:
    """

    pass


class ConnectorFile(Connector):
    """The C{ConnectorFile} class represents a C{file} connection.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_mode: File mode
    @type file_mode: str
    """

    def __init__(self, file_path, file_mode):
        """Initialise a C{ConnectorFile} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_mode: File mode
        @type file_mode: str
        @return:
        @rtype:
        """
        super(ConnectorFile, self).__init__()

        self.file_path = file_path
        self.file_mode = file_mode

        return


class ConnectorPipe(Connector):
    """The C{ConnectorPipe} class represents an abstract pipe to a downstream sub-process.

    Attributes:
    """

    def __init__(self):
        """Initialise a C{ConnectorPipe} object.

        @return:
        @rtype:
        """
        super(ConnectorPipe, self).__init__()

        return


class ConnectorPipeNamed(ConnectorFile):
    """The C{ConnectorPipeNamed} class represents a named pipe.

    Attributes:
    """

    pass


class ConnectorProcess(Connector):
    """The C{ConnectorProcess} class represents a concrete pipe to or from a sub-process.

    Attributes:
    """

    def __init__(self, name, connection):
        """Initialise a C{ConnectorProcess} object.

        @param name: C{bsf.process.Executable.name}
        @type name: str
        @param connection: Connection type stdin, stdout or stderr.
        @type connection: str
        """
        super(ConnectorProcess, self).__init__()

        self.name = name
        self.connection = connection

        return


class ConnectorSink(Connector):
    """The C{ConnectorSink} class represents a C{file} connection to /dev/null.

    Attributes:
    """

    pass


class ConnectorStandardStream(Connector):
    """The C{ConnectorStandardStream} class represents a standard stream processed via C{threading.Thread}.

    Standard streams (i.e. STDIN, STDOUT, STDERR) require processing via a C{threading.Thread} to prevent buffers from
    filling up and subsequently sub-processes (C{subprocess.Popen}) from blocking.
    Attributes:
    @ivar thread_callable: Callable for C{threading.Thread.target}
    @type thread_callable: object | None
    @ivar thread_kwargs: Python C{dict} of keyword arguments for C{threading.Thread.kwargs}
    @type thread_kwargs: dict[str, object] | None
    @ivar thread_joins: Maximum number of attempts calling C{threading.Thread.join}
    @type thread_joins: int
    @ivar thread_timeout: Timeout in seconds for calling C{threading.Thread.join}
    @type thread_timeout: int
    @ivar thread: C{threading.Thread} object | None
    @type thread: threading.Thread
    """

    def __init__(self, thread_callable=None, thread_kwargs=None, thread_joins=10, thread_timeout=10, thread=None):
        """Initialise a C{bsf.procedure.ConnectorStandardStream} object.

        @param thread_callable: Callable for C{threading.Thread.target}
        @type thread_callable: object | None
        @param thread_kwargs: Python C{dict} of keyword arguments for C{threading.Thread.kwargs}
        @type thread_kwargs: dict[str, object] | None
        @param thread_joins: Maximum number of attempts calling C{threading.Thread.join}
        @type thread_joins: int
        @param thread_timeout: Timeout in seconds for calling C{threading.Thread.join}
        @type thread_timeout: int
        @param thread: C{threading.Thread} object
        @type thread: threading.Thread | None
        @return:
        @rtype:
        """
        self.thread_callable = thread_callable
        self.thread_kwargs = thread_kwargs
        self.thread_joins = thread_joins
        self.thread_timeout = thread_timeout
        self.thread = thread

        return


class ConnectorStandardInput(ConnectorStandardStream):

    pass


class ConnectorStandardOutput(ConnectorStandardStream):

    pass


class ConnectorStandardError(ConnectorStandardStream):

    pass


class SubProcess(object):
    """The C{SubProcess} class represents one C{bsf.process.RunnableStep} to be run via C{subprocess.Popen}.

    Attributes:
    @ivar runnable_step: C{bsf.process.RunnableStep} or sub-class thereof
    @type runnable_step: bsf.process.RunnableStep
    @ivar stdin: C{Connector} for STDIN
    @type stdin: Connector
    @ivar stdout: C{Connector} for STDOUT
    @type stdout: Connector
    @ivar stderr: C{Connector} for STDERR
    @type stderr Connector
    """

    def __init__(self, runnable_step, stdin=None, stdout=None, stderr=None, sub_process=None):
        """Initialise a C{SubProcess} object.

        @param runnable_step: C{bsf.process.RunnableStep} or sub-class thereof
        @type runnable_step: bsf.process.RunnableStep
        @param stdin: C{Connector} for STDIN
        @type stdin: Connector | None
        @param stdout: C{Connector} for STDOUT
        @type stdout: Connector | None
        @param stderr: C{Connector} for STDERR
        @type stderr Connector | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: subprocess.Popen | None
        @return:
        @rtype:
        """
        super(SubProcess, self).__init__()

        self.runnable_step = runnable_step
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        self.sub_process = sub_process

        return


class ConcurrentRunnable(Runnable):
    """The C{bsf.procedure.ConcurrentRunnable} represents a procedure of concurrently running
    C{bsf.process.RunnableStep} or sub-classes thereof.

    Attributes:
    @ivar runnable_step_list_pre: Python C{list} of C{bsf.process.RunnableStep} object to pre-run
    @type runnable_step_list_pre: list[bsf.process.RunnableStep] | None
    @ivar sub_process_list: Python C{list} of C{SubProcess} objects
    @type sub_process_list: list[SubProcess]
    @ivar runnable_step_list_post: Python C{list} of C{bsf.process.RunnableStep} object to post-run
    @type runnable_step_list_post: list[bsf.process.RunnableStep] | None
    """

    def __init__(
            self,
            name,
            code_module,
            working_directory,
            cache_directory=None,
            cache_path_dict=None,
            file_path_object=None,
            debug=0,
            runnable_step_list_pre=None,
            sub_process_list=None,
            runnable_step_list_post=None):
        """Initialise a C{bsf.ConcurrentRunnable}.

        @param name: Name
        @type name: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.RunnableStep} objects via the C{bsf.procedure.Runnable.runner_script} concurrently
        @type code_module: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str | unicode
        @param cache_directory: Cache directory
        @type cache_directory: str | unicode | None
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the
            C{bsf.procedure.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str | unicode] | None
        @param file_path_object: C{bsf.procedure.FilePath}
        @type file_path_object: bsf.procedure.FilePath | None
        @param debug: Integer debugging level
        @type debug: int
        @param runnable_step_list_pre: Python C{list} of C{bsf.process.RunnableStep} object to pre-run
        @type runnable_step_list_pre: list[bsf.process.RunnableStep] | None
        @param sub_process_list: Python C{list} of C{SubProcess} objects
        @type sub_process_list: list[SubProcess] | None
        @param runnable_step_list_post: Python C{list} of C{bsf.process.RunnableStep} object to post-run
        @type runnable_step_list_post: list[bsf.process.RunnableStep] | None
        @return:
        @rtype:
        """
        super(ConcurrentRunnable, self).__init__(
            name=name,
            code_module=code_module,
            working_directory=working_directory,
            cache_directory=cache_directory,
            cache_path_dict=cache_path_dict,
            file_path_object=file_path_object,
            debug=debug)

        if runnable_step_list_pre is None:
            self.runnable_step_list_pre = list()
        else:
            self.runnable_step_list_pre = runnable_step_list_pre

        if sub_process_list is None:
            self.sub_process_list = list()
        else:
            self.sub_process_list = sub_process_list

        if runnable_step_list_post is None:
            self.runnable_step_list_post = list()
        else:
            self.runnable_step_list_post = runnable_step_list_post

        return

    def add_runnable_step_pre(self, runnable_step):
        """Convenience method to add a C{bsf.process.RunnableStep} to the pre-run list.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep
        @return: C{bsf.process.RunnableStep}
        @rtype: bsf.process.RunnableStep
        """
        if runnable_step is None:
            return

        assert isinstance(runnable_step, bsf.process.RunnableStep)

        self.runnable_step_list_pre.append(runnable_step)

        return runnable_step

    def add_runnable_step_post(self, runnable_step):
        """Convenience method to add a C{bsf.process.RunnableStep} to the post-run list.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep
        @return: C{bsf.process.RunnableStep}
        @rtype: bsf.process.RunnableStep
        """
        if runnable_step is None:
            return

        assert isinstance(runnable_step, bsf.process.RunnableStep)

        self.runnable_step_list_post.append(runnable_step)

        return runnable_step

    def add_sub_process(self, sub_process=None):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.procedure.SubProcess}.

        @param sub_process: C{bsf.procedure.SubProcess}
        @type sub_process: bsf.procedure.SubProcess | None
        @return: C{bsf.procedure.SubProcess}
        @rtype: bsf.procedure.SubProcess
        """
        if sub_process is None:
            return

        assert isinstance(sub_process, bsf.procedure.SubProcess)

        self.sub_process_list.append(sub_process)

        return sub_process
