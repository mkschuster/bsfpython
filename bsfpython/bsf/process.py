"""bsf.process

A package of classes and methods modelling processes.
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
import errno
import os
import shutil
import sys
import time
import warnings
from subprocess import Popen, PIPE
from threading import Lock, Thread

from bsf.argument import Argument, SwitchLong, SwitchShort, OptionLong, OptionShort, OptionPair
from bsf.standards import Configuration


class Command(object):
    """C{Command} class representing one program, its options and arguments and possibly
    another subordinate C{Command}.

    Attributes:
    @ivar program: Program
    @type program: str
    @ivar options: Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
        C{Argument} objects
    @type options: dict[Argument.key, list[Argument]]
    @ivar arguments: Python C{list} of Python C{str} (program argument) objects
    @type arguments: list[str]
    @ivar sub_command: Subordinate C{Command}
    @type sub_command: Command
    """

    def __init__(
            self,
            program=None,
            options=None,
            arguments=None,
            sub_command=None):
        """Initialise a C{Command} object.

        @param program: Program
        @type program: str
        @param options: Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str]
        @param sub_command: Subordinate C{Command} object
        @type sub_command: Command
        @return:
        @rtype:
        """

        super(Command, self).__init__()

        self.program = program  # Can be None.

        if options is None:
            self.options = dict()
        else:
            self.options = options

        if arguments is None:
            self.arguments = list()
        else:
            self.arguments = arguments

        self.sub_command = sub_command  # Can be None.

        return

    def trace(self, level):
        """Trace a C{Command} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  program:            {!r}\n'. \
            format(indent, self.program)

        # List all options

        output += '{}  options:\n'.format(indent)

        for key in self.options.keys():
            assert isinstance(key, str)
            output += '{}    key: {!r} Argument objects:\n'.format(indent, key)
            for argument in self.options[key]:
                assert isinstance(argument, Argument)
                output += argument.trace(level=level + 2)

        # List all arguments

        output += '{}  arguments:\n'.format(indent)

        i = 0
        for argument in self.arguments:
            assert isinstance(argument, str)
            output += '{}    {:2d}: {!r}\n'.format(indent, i, argument)
            i += 1

        if self.sub_command:
            output += self.sub_command.trace(level=level + 1)

        return output

    def add_argument(self, argument, override):
        """Add an C{Argument} or one of its sub-classes.

        @param argument: C{Argument}
        @type argument: Argument
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        assert isinstance(argument, Argument)
        assert isinstance(override, bool)

        if not override and argument.key in self.options:
            warnings.warn(
                'Adding an Argument with key {!r} that exits already in Command.program {!r}.'.format(
                    argument.key,
                    self.program),
                UserWarning)

        if argument.key in self.options:
            arguments_list = self.options[argument.key]
            assert isinstance(arguments_list, list)
        else:
            arguments_list = list()
            self.options[argument.key] = arguments_list

        arguments_list.append(argument)

        return

    def add_switch_long(self, key, override=False):
        """Initialise and add a C{SwitchLong} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=SwitchLong(key=key), override=override)

    def add_switch_short(self, key, override=False):
        """Initialise and add a C{SwitchShort} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=SwitchShort(key=key), override=override)

    def add_option_long(self, key, value, override=False):
        """Initialise and add an C{OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=OptionLong(key=key, value=value), override=override)

    def add_option_short(self, key, value, override=False):
        """Initialise and add an C{OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=OptionShort(key=key, value=value), override=override)

    def add_option_pair(self, key, value, override=False):
        """Initialise and add an C{OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_argument(self, argument, override):
        """Set an C{Argument} or one of its sub-classes.

        @param argument: C{Argument}
        @type argument: Argument
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        assert isinstance(argument, Argument)
        assert isinstance(override, bool)

        if not override and argument.key in self.options:
            warnings.warn(
                'Setting an Argument with key {!r} that exits already in Command.program {!r}.'.format(
                    argument.key,
                    self.program),
                UserWarning)

        self.options[argument.key] = [argument]

        return

    def set_switch_long(self, key, override=False):
        """Initialise and set a C{SwitchLong} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=SwitchLong(key=key), override=override)

    def set_switch_short(self, key, override=False):
        """Initialise and set a C{SwitchShort} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=SwitchShort(key=key), override=override)

    def set_option_long(self, key, value, override=False):
        """Initialise and set an C{OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @type value: str | unicode
        @return:
        @rtype:
        """
        return self.set_argument(argument=OptionLong(key=key, value=value), override=override)

    def set_option_short(self, key, value, override=False):
        """Initialise and set an C{OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=OptionShort(key=key, value=value), override=override)

    def set_option_pair(self, key, value, override=False):
        """Initialise and set an C{OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{Command} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section, defaults to instance class
        @type section: str
        @return:
        @rtype:
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            warnings.warn(
                'Section {!r} not defined in Configuration files: {!r}'.format(section, configuration.file_path_list),
                UserWarning)

            return

        # The configuration section is available.

        for option in configuration.config_parser.options(section=section):
            self.add_argument(
                argument=Argument.from_key_value(
                    key=option,
                    value=configuration.config_parser.get(
                        section=section,
                        option=option)),
                override=False)

        return

    def command_list(self):
        """Assemble the command line from program, options and arguments.

        @return: Python C{list} of program, options, switches and arguments
        @rtype: list[str | unicode]
        """

        command_line = list()

        if self.program:
            command_line.append(self.program)

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            assert isinstance(key, str)
            options_list = self.options[key]
            assert isinstance(options_list, list)
            for argument in options_list:
                assert isinstance(argument, Argument)
                if isinstance(argument, SwitchLong):
                    command_line.append('--{}'.format(argument.key))
                elif isinstance(argument, SwitchShort):
                    command_line.append('-{}'.format(argument.key))
                elif isinstance(argument, OptionLong):
                    command_line.append('--{}'.format(argument.key))
                    if argument.value:
                        command_line.append(str(argument.value))
                elif isinstance(argument, OptionShort):
                    command_line.append('-{}'.format(argument.key))
                    if argument.value:
                        command_line.append(str(argument.value))
                elif isinstance(argument, OptionPair):
                    command_line.append('{}={}'.format(argument.key, argument.value))
                else:
                    warnings.warn(
                        'Unexpected object {!r} in Command.options dict.'.format(argument),
                        UserWarning)

        # Add all arguments.

        for argument in self.arguments:
            assert isinstance(argument, str)
            command_line.append(str(argument))

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line.extend(self.sub_command.command_list())

        return command_line

    def command_str(self):
        """Assemble the command line from program, options, switches and arguments.

        @return: A Python C{str} of program, options, switches and arguments
        @rtype: str
        """

        command_line = str()

        if self.program:
            command_line += self.program

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            assert isinstance(key, str)
            options_list = self.options[key]
            assert isinstance(options_list, list)
            for argument in options_list:
                assert isinstance(argument, Argument)
                if isinstance(argument, SwitchLong):
                    command_line += ' --{}'.format(argument.key)
                elif isinstance(argument, SwitchShort):
                    command_line += ' -{}'.format(argument.key)
                elif isinstance(argument, OptionLong):
                    command_line += ' --{} {}'.format(argument.key, argument.value)
                elif isinstance(argument, OptionShort):
                    command_line += ' -{} {}'.format(argument.key, argument.value)
                elif isinstance(argument, OptionPair):
                    command_line += ' {}={}'.format(argument.key, argument.value)
                else:
                    warnings.warn(
                        'Unexpected object {!r} in Command.options dict.'.format(argument),
                        UserWarning)

        # Add all arguments.

        for argument in self.arguments:
            assert isinstance(argument, str)
            command_line += ' '
            command_line += argument

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line += ' '
            command_line += self.sub_command.command_str()

        return command_line


class Executable(Command):
    """The Executable class represents an executable program,
    its options and arguments.

    Attributes:
    @ivar name: Name in the context of a DRMS dependency
    @type name: str
    @ivar program: Program (executable or full file path)
    @type program: str
    @ivar options: Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
        C{Argument} objects
    @type options: dict[Argument.key, list[Argument]]
    @ivar arguments: Python C{list} of Python C{str} or C{unicode} (argument) objects
    @type arguments: list[str | unicode]
    @ivar sub_command: Subordinate Command
    @type sub_command: bsf.process.Command
    @ivar stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
    @type stdout_path: str | unicode
    @ivar stderr_path: Standard error (STDERR) redirection in Bash (2>word)
    @type stderr_path: str | unicode
    @ivar dependencies: Python C{list} of C{Executable.name} properties in the
        context of C{DRMS} dependencies
    @type dependencies: list[Executable.name]
    @ivar hold: Hold on job scheduling
    @type hold: str
    @ivar submit: Submit the Executable into the DRMS
    @type submit: bool
    @ivar maximum_attempts: Maximum number of attempts to run this C{Executable}
    @type maximum_attempts: int
    @ivar process_identifier: Process identifier
    @type process_identifier: str
    @ivar process_name: Process name
    @type process_name: str
    """

    @staticmethod
    def process_stream(file_type, file_handle, thread_lock, file_path=None, debug=0):
        """C{Executable} function to process I{STDOUT} or I{STDERR} from the child process as a thread.

        @param file_type: File handle type I{STDOUT} or I{STDERR}
        @type file_type: str
        @param file_handle: The I{STDOUT} or I{STDERR} file handle
        @type file_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: thread.lock
        @param file_path: I{STDOUT} file path
        @type file_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        @raise Exception: The file_type has to be either I{STDOUT} or I{STDERR}
        """

        if file_type not in ('STDOUT', 'STDERR'):
            raise Exception('The file_type has to be either STDOUT or STDERR.')

        thread_lock.acquire(True)
        if debug > 0:
            print '[{}] Started Runner {} processor in module {}.'. \
                format(datetime.datetime.now().isoformat(), file_type, __name__)
        output_file = None
        if file_path:
            output_file = open(file_path, 'w')
            if debug > 0:
                print '[{}] Opened {} file {!r}.'. \
                    format(datetime.datetime.now().isoformat(), file_type, file_path)
        thread_lock.release()

        for line in file_handle:
            thread_lock.acquire(True)
            if output_file:
                output_file.write(line)
            else:
                print '[{}] {}: {}'.format(datetime.datetime.now().isoformat(), file_type, line.rstrip())
            thread_lock.release()

        thread_lock.acquire(True)
        if debug > 0:
            print '[{}] Received EOF on {} pipe.'.format(datetime.datetime.now().isoformat(), file_type)
        if output_file:
            output_file.close()
            if debug > 0:
                print '[{}] Closed {} file {!r}.'. \
                    format(datetime.datetime.now().isoformat(), file_type, file_path)
        thread_lock.release()

        return

    @staticmethod
    def process_stdout(stdout_handle, thread_lock, stdout_path=None, debug=0):
        """C{Executable} function to process I{STDOUT} from the child process as a thread.

        @param stdout_handle: The I{STDOUT} file handle
        @type stdout_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: thread.lock
        @param stdout_path: I{STDOUT} file path
        @type stdout_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        return Executable.process_stream(
            file_type='STDOUT',
            file_handle=stdout_handle,
            thread_lock=thread_lock,
            file_path=stdout_path,
            debug=debug)

    @staticmethod
    def process_stderr(stderr_handle, thread_lock, stderr_path=None, debug=0):
        """C{Executable} function to process I{STDERR} from the child process as a thread.

        @param stderr_handle: The I{STDERR} file handle
        @type stderr_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: thread.lock
        @param stderr_path: I{STDERR} file path
        @type stderr_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        return Executable.process_stream(
            file_type='STDERR',
            file_handle=stderr_handle,
            thread_lock=thread_lock,
            file_path=stderr_path,
            debug=debug)

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            maximum_attempts=1,
            process_identifier=None,
            process_name=None):
        """Initialise an C{Executable} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param maximum_attempts: Maximum number of attempts to run this C{Executable}
        @type maximum_attempts: int
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @return:
        @rtype:
        """

        super(Executable, self).__init__(
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command)

        self.name = name  # Can be None.

        if stderr_path is None:
            self.stderr_path = str()
        else:
            self.stderr_path = stderr_path

        if stdout_path is None:
            self.stdout_path = str()
        else:
            self.stdout_path = stdout_path

        if dependencies is None:
            self.dependencies = list()
        else:
            self.dependencies = dependencies

        if hold is None:
            self.hold = str()
        else:
            self.hold = hold

        if submit is None:
            self.submit = True
        else:
            assert isinstance(submit, bool)
            self.submit = submit

        if maximum_attempts is None:
            self.maximum_attempts = int(x=1)
        else:
            self.maximum_attempts = maximum_attempts

        if process_identifier is None:
            self.process_identifier = str()
        else:
            self.process_identifier = process_identifier

        if process_name is None:
            self.process_name = str()
        else:
            self.process_name = process_name

        return

    def trace(self, level):
        """Trace an C{Executable} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  name:               {!r}\n'. \
            format(indent, self.name)
        output += '{}  stdout:             {!r}\n'. \
            format(indent, self.stdout_path)
        output += '{}  stderr:             {!r}\n'. \
            format(indent, self.stderr_path)
        output += '{}  hold:               {!r}\n'. \
            format(indent, self.hold)
        output += '{}  submit:             {!r}\n'. \
            format(indent, self.submit)
        output += '{}  maximum_attempts: {!r}\n'. \
            format(indent, self.maximum_attempts)
        output += '{}  process_identifier: {!r}\n'. \
            format(indent, self.process_identifier)
        output += '{}  process_name:       {!r}\n'. \
            format(indent, self.process_name)

        # List all dependencies.

        output += '{}  dependencies:\n'.format(indent)

        i = 0
        for dependency in self.dependencies:
            assert isinstance(dependency, str)
            output += '{}    {:2d} {!r}\n'.format(indent, i, dependency)
            i += 1

        # Trace the Command super-class.

        output += super(Executable, self).trace(level=level + 1)

        return output

    def command_list(self):
        """Assemble the command line from program, options and arguments.

        @return: Python C{list} of program, options and arguments
        @rtype: list[str | unicode]
        """

        command = list()

        command.extend(super(Executable, self).command_list())

        # The stdout_path and stderr_path gets appended in specific modules.

        return command

    def command_str(self):
        """Assemble the command line from program, options, switches and arguments.

        @return: A Python C{str} of program, options, switches and arguments
        @rtype: str
        """

        command = str()

        command += super(Executable, self).command_str()

        # The stdout_path and stderr_path gets appended in specific modules.

        return command

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run an C{Executable} object via the Python C{subprocess.Popen} class.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """
        on_posix = 'posix' in sys.builtin_module_names

        child_return_code = 0
        attempt_counter = 0

        while attempt_counter < self.maximum_attempts:

            child_process = Popen(
                args=self.command_list(),
                bufsize=0,
                stdin=PIPE,
                stdout=PIPE,
                stderr=PIPE,
                shell=False,
                close_fds=on_posix)

            # Two threads, thread_out and thread_err reading STDOUT and STDERR, respectively,
            # should make sure that buffers are not filling up.

            thread_lock = Lock()

            thread_out = Thread(
                target=Executable.process_stdout,
                kwargs=dict(
                    stdout_handle=child_process.stdout,
                    thread_lock=thread_lock,
                    stdout_path=self.stdout_path,
                    debug=debug))
            thread_out.daemon = True  # Thread dies with the program.
            thread_out.start()

            thread_err = Thread(
                target=Executable.process_stderr,
                kwargs=dict(
                    stderr_handle=child_process.stderr,
                    thread_lock=thread_lock,
                    stderr_path=self.stderr_path,
                    debug=debug))
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
                        format(datetime.datetime.now().isoformat(), self.name, +child_return_code)
                attempt_counter += 1
            elif child_return_code < 0:
                if debug > 0:
                    print '[{}] Child process {!r} received signal {}.'. \
                        format(datetime.datetime.now().isoformat(), self.name, -child_return_code)
            else:
                if debug > 0:
                    print '[{}] Child process {!r} completed successfully {}.'. \
                        format(datetime.datetime.now().isoformat(), self.name, +child_return_code)
                break

        else:
            if debug > 0:
                print '[{}] Runnable {!r} exceeded the maximum retry counter {}.' \
                    .format(datetime.datetime.now().isoformat(), self.name, self.maximum_attempts)

        return child_return_code

    def evaluate_return_code(self, return_code):
        """Evaluate a return code from the run method.

        @param return_code: Return code
        @type return_code: int
        @return:
        @rtype:
        """

        if return_code > 0:
            print '[{}] Child process {!r} failed with return code {}'. \
                format(datetime.datetime.now().isoformat(), self.name, +return_code)
        elif return_code < 0:
            print '[{}] Child process {!r} received signal {}.'. \
                format(datetime.datetime.now().isoformat(), self.name, -return_code)
        else:
            print '[{}] Child process {!r} completed with return code {}.'. \
                format(datetime.datetime.now().isoformat(), self.name, +return_code)

        return


class RunnableStep(Executable):
    """The C{RunnableStep} represents a step in a C{Runnable} class.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{RunnableStep}
    @type obsolete_file_path_list: list[str | unicode]
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None):
        """Initialise a C{RunnableStep} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @return:
        @rtype:
        """

        super(RunnableStep, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name)

        if obsolete_file_path_list is None:
            self.obsolete_file_path_list = list()
        else:
            self.obsolete_file_path_list = obsolete_file_path_list

        return

    def trace(self, level=1):
        """Trace a C{RunnableStep} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  obsolete_file_path_list: {!r}\n'.format(indent, self.obsolete_file_path_list)
        output += super(RunnableStep, self).trace(level=level + 1)

        return output


class RunnableStepJava(RunnableStep):
    """The C{RunnableStepJava} class represents a C{RunnableStep} with all peculiarities of Java programs.

    Attributes:
    None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_maximum=None,
            java_jar_path=None):
        """Create a C{RunnableStep} for a Java program.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @return: C{RunnableStep}
        @rtype: RunnableStep
        """

        super(RunnableStepJava, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        # JavaVM command
        if self.program is None:
            self.program = 'java'

        # JavaVM options
        if 'd64' not in self.options:
            self.add_switch_short(key='d64')

        if 'server' not in self.options:
            self.add_switch_short(key='server')

        if java_heap_maximum and java_heap_maximum not in self.options:
            self.add_switch_short(key=java_heap_maximum)

        if '-Djava.io.tmpdir' not in self.options:
            self.add_option_pair(key='-Djava.io.tmpdir', value=java_temporary_path)

        if self.sub_command is None:
            # The Picard command line interface is a bit broken, as the -jar option needs to come last,
            # just before the Picard command. GATK does this slightly better with the --analysis_type option.
            # To be on the safe side, an empty sub command is required to separate the -jar option from the
            # other JavaVM options.
            self.sub_command = Command()
            if java_jar_path is not None:
                self.sub_command.add_option_short(key='jar', value=java_jar_path)

        return


class RunnableStepPicard(RunnableStepJava):
    """The C{RunnableStepPicard} class represents a C{RunnableStepJava} specific to Picard tools.

    Attributes:
    None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_maximum=None,
            java_jar_path=None,
            picard_classpath=None,
            picard_command=None):
        """Create a C{RunnableStep} for a Picard algorithm.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @param picard_classpath: Picard class path
        @type picard_classpath: str | unicode
        @param picard_command: Picard command
        @type picard_command: str
        @return: C{RunnableStep}
        @rtype: RunnableStep
        """

        super(RunnableStepPicard, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path,
            java_heap_maximum=java_heap_maximum,
            java_jar_path=java_jar_path)

        # Set the Picard classpath and the Picard Java archive.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(picard_classpath, 'picard.jar'))

        # The Picard algorithm is then another sub-command.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command(program=picard_command)

        return

    def add_picard_option(self, key, value, override=False):
        """Add an option to the Picard command.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class RunnableStepLink(RunnableStep):
    """The C{RunnableStepLink} represents a step in a C{Runnable} class.

    Attributes:
    @ivar source_path: Source path
    @type source_path: str | unicode
    @ivar target_path: Target path
    @type target_path: str | unicode
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            source_path=None,
            target_path=None):
        """Initialise a C{RunnableStepLink} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param source_path: Source path
        @type source_path: str | unicode
        @param target_path: Target path
        @type target_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepLink, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if source_path is None:
            self.source_path = str()
        else:
            self.source_path = source_path

        if target_path is None:
            self.target_path = str()
        else:
            self.target_path = target_path

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """

        if self.source_path and self.target_path and not os.path.exists(self.target_path):
            try:
                os.symlink(self.source_path, self.target_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return 0


class RunnableStepMakeDirectory(RunnableStep):
    """The C{RunnableStepMakeDirectory} represents a step in a C{Runnable} class.

    Attributes:
    @ivar directory_path: Directory path
    @type directory_path: str | unicode
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            directory_path=None):
        """Initialise a C{RunnableStepMakeDirectory} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param directory_path: Directory path
        @type directory_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepMakeDirectory, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if directory_path is None:
            self.directory_path = str()
        else:
            self.directory_path = directory_path

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """

        if self.directory_path and not os.path.isdir(self.directory_path):
            try:
                os.makedirs(self.directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return 0


class RunnableStepMove(RunnableStep):
    """The C{RunnableStepMove} represents a step in a C{Runnable} class.

    Attributes:
    @ivar source_path: Source path
    @type source_path: str | unicode
    @ivar target_path: Target path
    @type target_path: str | unicode
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            source_path=None,
            target_path=None):
        """Initialise a C{RunnableStepMove} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param source_path: Source path
        @type source_path: str | unicode
        @param target_path: Target path
        @type target_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepMove, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if source_path is None:
            self.source_path = str()
        else:
            self.source_path = source_path

        if target_path is None:
            self.target_path = str()
        else:
            self.target_path = target_path

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """

        if self.source_path and self.target_path:
            # os.rename(self.source_path, self.target_path)
            shutil.move(src=self.source_path, dst=self.target_path)

        return 0


class RunnableStepSleep(RunnableStep):
    """The C{RunnableStepSleep} represents a step in a C{Runnable} class.

    Attributes:
    @ivar sleep_time: Sleep time in seconds
    @type sleep_time: float
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            sleep_time=None):
        """Initialise a C{RunnableStepSleep} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: bsf.process.Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param sleep_time: Sleep time in seconds
        @type sleep_time: float
        @return:
        @rtype:
        """

        super(RunnableStepSleep, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if sleep_time is None:
            self.sleep_time = float()
        else:
            assert isinstance(sleep_time, float)
            self.sleep_time = sleep_time

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """

        time.sleep(self.sleep_time)

        return 0
