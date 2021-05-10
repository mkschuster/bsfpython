# -*- coding: utf-8 -*-
"""Process module.

A package of classes and methods modelling processes.
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
import datetime
import errno
import os
import shutil
import stat
import sys
import time
import warnings
from io import IOBase, TextIOWrapper
from subprocess import Popen, PIPE, DEVNULL
from threading import Lock, Thread
from typing import List

from bsf.argument import *
from bsf.connector import *
from bsf.standards import Configuration


def get_timestamp():
    """Get the current time stamp in ISO format.

    @return: ISO format time stamp
    @rtype: str
    """
    return '[' + datetime.datetime.now().isoformat() + ']'


def map_connector(connector=None, executable_list=None):
    """Map a C{bsf.connector.Connector} object to a file handle.

    @param connector: C{bsf.connector.Connector} object or sub-class thereof.
    @type connector: Connector | None
    @param executable_list: Python C{list} object of C{Executable} objects.
    @type executable_list: list[Executable] | None
    @return: File handle object
    @rtype: IOBase | DEVNULL | PIPE | None
    """
    if isinstance(connector, ElectronicSink):
        return DEVNULL

    if isinstance(connector, ConnectorFile):
        if isinstance(connector, ConnectorPipeNamed):
            # A named pipe needs creating before it can be opened.
            if not os.path.exists(connector.file_path):
                os.mkfifo(connector.file_path)
        return open(file=connector.file_path, mode=connector.file_mode)

    if isinstance(connector, ConnectorPipe):
        return PIPE

    if isinstance(connector, StandardStream):
        return PIPE

    if isinstance(connector, ConcurrentProcess):
        if executable_list is None:
            raise Exception('A ConcurrentProcess Connector requires a list of Executable objects to connect to.')

        for executable in executable_list:
            if executable.name == connector.name:
                if executable.sub_process is None:
                    raise Exception(
                        'Sub-process ' + repr(executable.name) + ' not initialised, yet.')
                return getattr(executable.sub_process, connector.connection)
        else:
            raise Exception('Could not find a suitable Executable.name for Connector.name ' + repr(connector.name))

    return None


def run_executables(executable_list, debug=0):
    """Run a Python C{list} of C{Executable} objects concurrently.

    @param executable_list: Python C{list} of C{Executable} objects
    @type executable_list: list[Executable]
    @param debug: Integer debugging level
    @type debug: int
    @return: Python C{list} of Python C{str} (exception) objects
    @rtype: list[str] | None
    """
    thread_lock = Lock()

    for executable in executable_list:
        # Per RunnableStep, up to three threads, writing to STDIN, as well as reading from STDOUT and STDERR,
        # should make sure that buffers are not filling up. If STDOUT or STDERR connectors are not defined,
        # defaults need be set to avoid the sub-process from blocking.

        if executable.stdout is None:
            executable.stdout = StandardOutputStream()

        if executable.stderr is None:
            executable.stderr = StandardErrorStream()

        # Create and assign sub-processes.
        try:
            executable.sub_process = Popen(
                args=executable.command_list(),
                bufsize=1,
                stdin=map_connector(connector=executable.stdin, executable_list=executable_list),
                stdout=map_connector(connector=executable.stdout, executable_list=executable_list),
                stderr=map_connector(connector=executable.stderr, executable_list=executable_list),
                text=True)
        except OSError as exception:
            if exception.errno == errno.ENOENT:
                raise Exception(
                    "For Executable.name {!r} Executable.program {!r} could not be found.".format(
                        executable.name, executable.program))
            else:
                # Re-raise the Exception object.
                raise exception

        for attribute in ('stdin', 'stdout', 'stderr'):
            connector: Connector = getattr(executable, attribute)

            if isinstance(connector, StandardInputStream):
                if connector.thread_callable is None:
                    # If a specific STDIN callable is not defined, run bsf.process.Executable.process_stdin().
                    pass
                else:
                    connector.thread = Thread(
                        target=connector.thread_callable,
                        args=[executable.sub_process.stdin, thread_lock, debug],
                        kwargs=connector.thread_kwargs)

            if isinstance(connector, StandardOutputStream):
                if connector.thread_callable is None:
                    # If a specific STDOUT callable is not defined, run bsf.process.Executable.process_stdout().
                    connector.thread = Thread(
                        target=Executable.process_stdout,
                        args=[executable.sub_process.stdout, thread_lock, debug],
                        kwargs={'stdout_path': connector.file_path})
                else:
                    connector.thread = Thread(
                        target=connector.thread_callable,
                        args=[executable.sub_process.stdout, thread_lock, debug],
                        kwargs=connector.thread_kwargs)

            if isinstance(connector, StandardErrorStream):
                if connector.thread_callable is None:
                    # If a specific STDERR callable is not defined, run bsf.process.Executable.process_stderr().
                    connector.thread = Thread(
                        target=Executable.process_stderr,
                        args=[executable.sub_process.stderr, thread_lock, debug],
                        kwargs={'stderr_path': connector.file_path})
                else:
                    connector.thread = Thread(
                        target=connector.thread_callable,
                        args=[executable.sub_process.stderr, thread_lock, debug],
                        kwargs=connector.thread_kwargs)

            if isinstance(connector, StandardStream) and connector.thread:
                connector.thread.daemon = True
                connector.thread.start()

    # At this stage all subprocess.Popen objects should have been created.
    # Now, wait for all child processes to complete.

    exception_str_list = list()
    for executable in executable_list:
        child_return_code = executable.sub_process.wait()

        # First, join all standard stream processing threads.

        for attribute in ('stdin', 'stdout', 'stderr'):
            connector: Connector = getattr(executable, attribute)
            if isinstance(connector, StandardStream) and connector.thread:
                thread_join_counter = 0
                while connector.thread.is_alive() and thread_join_counter < connector.thread_joins:
                    if debug > 0:
                        thread_lock.acquire(True)
                        print("{} Waiting for Executable.name {!r} {!r} processor to finish.".format(
                            get_timestamp(),
                            executable.name,
                            attribute))
                        thread_lock.release()

                    connector.thread.join(timeout=connector.thread_timeout)
                    thread_join_counter += 1

        # Second, inspect the child process' return code.

        if child_return_code > 0:
            # Child return code.
            exception_str_list.append("{} Child process Executable.name {!r} failed with return code {!r}.".format(
                get_timestamp(), executable.name, +child_return_code))
        elif child_return_code < 0:
            # Child signal.
            exception_str_list.append("{} Child process Executable.name {!r} received signal {!r}.".format(
                get_timestamp(), executable.name, -child_return_code))

    return exception_str_list


class Command(object):
    """The C{bsf.process.Command} class represents one program, its options and arguments.

    A C{bsf.process.Command} object can possibly contain another subordinate C{bsf.process.Command} object.

    @ivar name: Name
    @type name: str
    @ivar program: Program
    @type program: str
    @ivar options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
        Python C{list} value objects of C{bsf.argument.Argument} objects
    @type options: dict[Argument.key, list[Argument]]
    @ivar arguments: Python C{list} of Python C{str} (program argument) objects
    @type arguments: list[str]
    @ivar sub_command: Subordinate C{bsf.process.Command} object
    @type sub_command: Command | None
    """

    def __init__(
            self,
            name=None,
            program=None,
            options=None,
            arguments=None,
            sub_command=None):
        """Initialise a C{bsf.process.Command} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        """

        super(Command, self).__init__()

        if name is None:
            self.name = str()
        else:
            self.name = name

        if program is None:
            self.program = str()
        else:
            self.program = program

        if options is None:
            self.options = dict()
        else:
            self.options = options

        if arguments is None:
            self.arguments = list()
        else:
            self.arguments = arguments

        self.sub_command = sub_command

        return

    def trace(self, level):
        """Trace a C{bsf.process.Command} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name:               {!r}\n'.format(indent, self.name))
        str_list.append('{}  program:            {!r}\n'.format(indent, self.program))

        # List all options

        str_list.append('{}  options:\n'.format(indent))

        for key in sorted(self.options):
            str_list.append('{}    key: {!r} Argument objects:\n'.format(indent, key))
            for argument in self.options[key]:
                str_list.extend(argument.trace(level=level + 2))

        # List all arguments

        str_list.append('{}  arguments:\n'.format(indent))

        i = 0
        for argument in self.arguments:
            str_list.append('{}    {:2d}: {!r}\n'.format(indent, i, argument))
            i += 1

        if self.sub_command is not None:
            str_list.extend(self.sub_command.trace(level=level + 1))

        return str_list

    def add_argument(self, argument, override):
        """Add a C{bsf.argument.Argument} object or one of its sub-classes.

        @param argument: C{bsf.argument.Argument} object
        @type argument: Argument
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        if not override and argument.key in self.options:
            warnings.warn(
                'Adding an Argument with key ' + repr(argument.key) +
                ' that exits already in Command.program ' + repr(self.program) + '.',
                UserWarning)

        if argument.key not in self.options:
            self.options[argument.key] = list()

        arguments_list = self.options[argument.key]
        arguments_list.append(argument)

        return

    def add_switch_long(self, key, override=False):
        """Initialise and add a C{bsf.argument.SwitchLong} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=SwitchLong(key=key), override=override)

    def add_switch_short(self, key, override=False):
        """Initialise and add a C{bsf.argument.SwitchShort} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=SwitchShort(key=key), override=override)

    def add_option_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionLong(key=key, value=value), override=override)

    def add_option_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionShort(key=key, value=value), override=override)

    def add_option_pair(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionPair(key=key, value=value), override=override)

    def add_option_pair_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPairShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionPairShort(key=key, value=value), override=override)

    def add_option_pair_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPairLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionPairLong(key=key, value=value), override=override)

    def add_option_multi(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMulti} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMulti(key=key, value=value), override=override)

    def add_option_multi_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiLong(key=key, value=value), override=override)

    def add_option_multi_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiShort(key=key, value=value), override=override)

    def add_option_multi_pair(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiPair(key=key, value=value), override=override)

    def add_option_multi_pair_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPairLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiPairLong(key=key, value=value), override=override)

    def add_option_multi_pair_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPairShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiPairShort(key=key, value=value), override=override)

    def set_argument(self, argument, override):
        """Set a C{bsf.argument.Argument} objects or one of its sub-classes.

        @param argument: C{bsf.argument.Argument} object
        @type argument: Argument
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        if not override and argument.key in self.options:
            warnings.warn(
                'Setting an Argument with key ' + repr(argument.key) +
                ' that exits already in Command.program ' + repr(self.program) + '.',
                UserWarning)

        self.options[argument.key] = [argument]

        return

    def set_switch_long(self, key, override=False):
        """Initialise and set a C{bsf.argument.SwitchLong} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=SwitchLong(key=key), override=override)

    def set_switch_short(self, key, override=False):
        """Initialise and set a C{bsf.argument.SwitchShort} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=SwitchShort(key=key), override=override)

    def set_option_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionLong(key=key, value=value), override=override)

    def set_option_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionShort(key=key, value=value), override=override)

    def set_option_pair(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_option_pair_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPairShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionPairShort(key=key, value=value), override=override)

    def set_option_pair_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPairLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionPairLong(key=key, value=value), override=override)

    def set_option_multi(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMulti} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMulti(key=key, value=value), override=override)

    def set_option_multi_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiLong(key=key, value=value), override=override)

    def set_option_multi_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiShort(key=key, value=value), override=override)

    def set_option_multi_pair(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiPair(key=key, value=value), override=override)

    def set_option_multi_pair_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPairLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiPairLong(key=key, value=value), override=override)

    def set_option_multi_pair_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPairShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiPairShort(key=key, value=value), override=override)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.process.Command} object via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration} object
        @type configuration: Configuration
        @param section: Configuration file section, defaults to instance class
        @type section: str
        """
        if not configuration.config_parser.has_section(section=section):
            warnings.warn(
                'Section ' + repr(section) + ' not defined in configuration files:\n' +
                repr(configuration.from_file_path_list),
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
        @rtype: list[str]
        """
        command_line: List[str] = list()

        if self.program:
            command_line.append(self.program)

        # Add all options and switches in alphabetical order.

        for argument_key in sorted(self.options):
            for argument in self.options[argument_key]:
                command_line.extend(argument.get_list())

        # Add all arguments.

        for argument in self.arguments:
            command_line.append(argument)

        # Expand a subordinate command, if defined.

        if self.sub_command is not None:
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

        for argument_key in sorted(self.options):
            for argument in self.options[argument_key]:
                command_line += ' ' + argument.get_str()

        # Add all arguments.

        for argument in self.arguments:
            command_line += ' '
            command_line += argument

        # Expand a subordinate command, if defined.

        if self.sub_command is not None:
            command_line += ' '
            command_line += self.sub_command.command_str()

        return command_line


class Executable(Command):
    """The C{bsf.process.Executable} class represents one C{bsf.process.Command} as UNIX process.

    @ivar stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
    @type stdin: Connector | None
    @ivar stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
    @type stdout: Connector | None
    @ivar stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
    @type stderr: Connector | None
    @ivar dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
        in the context of C{bsf.analysis.Stage} dependencies
    @type dependencies: list[Executable.name]
    @ivar hold: Hold on job scheduling
    @type hold: str | None
    @ivar submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
    @type submit: bool
    @ivar maximum_attempts: Maximum number of attempts to run this C{bsf.process.Executable} object
    @type maximum_attempts: int
    @ivar process_identifier: Process identifier
    @type process_identifier: str | None
    @ivar process_name: Process name
    @type process_name: str | None
    @ivar sub_process: C{subprocess.Popen} object
    @type sub_process: Popen | None
    """

    @staticmethod
    def process_stream(file_handle, thread_lock, debug, file_type, file_path=None):
        """Process a I{STDOUT} or I{STDERR} text stream from the child process as a thread.

        If a file_path was provided, a corresponding Python C{file} will be opened in text mode,
        if not, C{sys.stdout} or C{sys.stderr} will be used according to the I{file_type}.
        If a debug level was set, diagnostic output will be printed to C{sys.stdout}, as well as the output stream.
        @param file_handle: The I{STDOUT} or I{STDERR} C{io.TextIOWrapper} object
        @type file_handle: TextIOWrapper
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: Lock
        @param debug: Debug level
        @type debug: int
        @param file_type: File handle type I{STDOUT} or I{STDERR}
        @type file_type: str
        @param file_path: I{STDOUT} file path
        @type file_path: str | None
        @raise Exception: If file_type is neither I{STDOUT} nor I{STDERR}
        """
        if file_type not in ('STDOUT', 'STDERR'):
            raise Exception('The file_type has to be either STDOUT or STDERR.')

        thread_lock.acquire(True)
        if debug > 0:
            print(get_timestamp(),
                  'Started Runner ' + repr(file_type) + ' processor in module ' + repr(__name__) + '.',
                  file=sys.stdout, flush=True)
        output_file = None
        if file_path:
            output_file = open(file=file_path, mode='wt')
            if debug > 0:
                print(get_timestamp(), 'Opened ' + repr(file_type) + ' file ' + repr(file_path) + '.',
                      file=sys.stdout, flush=True)
        elif file_type == 'STDOUT':
            output_file = sys.stdout
        elif file_type == 'STDERR':
            output_file = sys.stderr
        thread_lock.release()

        for line_str in file_handle:
            thread_lock.acquire(True)
            if debug > 0:
                print(get_timestamp(), file_type + ': ' + line_str.rstrip(),
                      file=output_file, flush=True)
            else:
                print(line_str.rstrip(), file=output_file, flush=True)
            thread_lock.release()

        thread_lock.acquire(True)
        if debug > 0:
            print(get_timestamp(), 'Received EOF on ' + repr(file_type) + ' pipe.',
                  file=sys.stdout, flush=True)
        if file_path:
            output_file.close()
            if debug > 0:
                print(get_timestamp(), 'Closed ' + repr(file_type) + ' file ' + repr(file_path) + '.',
                      file=sys.stdout, flush=True)
        thread_lock.release()

        return

    @staticmethod
    def process_stdout(stdout_handle, thread_lock, debug, stdout_path=None):
        """Process I{STDOUT} from the child process as a thread.

        @param stdout_handle: The I{STDOUT} C{io.TextIOWrapper} object
        @type stdout_handle: TextIOWrapper
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: Lock
        @param debug: Debug level
        @type debug: int
        @param stdout_path: I{STDOUT} file path
        @type stdout_path: str | None
        """
        return Executable.process_stream(
            file_handle=stdout_handle,
            thread_lock=thread_lock,
            debug=debug,
            file_type='STDOUT',
            file_path=stdout_path)

    @staticmethod
    def process_stderr(stderr_handle, thread_lock, debug, stderr_path=None):
        """Process I{STDERR} from the child process as a thread.

        @param stderr_handle: The I{STDERR} C{io.TextIOWrapper} object
        @type stderr_handle: TextIOWrapper
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: Lock
        @param debug: Debug level
        @type debug: int
        @param stderr_path: I{STDERR} file path
        @type stderr_path: str | None
        """
        return Executable.process_stream(
            file_handle=stderr_handle,
            thread_lock=thread_lock,
            debug=debug,
            file_type='STDERR',
            file_path=stderr_path)

    def __init__(
            self,
            name=None,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            maximum_attempts=1,
            process_identifier=None,
            process_name=None,
            sub_process=None):
        """Initialise a C{bsf.process.Executable} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param maximum_attempts: Maximum number of attempts to run this C{bsf.process.Executable} object
        @type maximum_attempts: int
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        """
        # Further constrain the name instance variable in the Executable class.

        if not name:
            raise Exception(
                'The Executable class requires a non-empty name option ' + repr(name) +
                ' for program ' + repr(program) + '.')

        super(Executable, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command)

        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr

        if dependencies is None:
            self.dependencies = list()
        else:
            self.dependencies = dependencies

        self.hold = hold

        if submit is None:
            self.submit = True
        else:
            assert isinstance(submit, bool)
            self.submit = submit

        if maximum_attempts is None:
            self.maximum_attempts = 1
        else:
            self.maximum_attempts = maximum_attempts

        self.process_identifier = process_identifier
        self.process_name = process_name

        self.sub_process = sub_process

        return

    def trace(self, level):
        """Trace a C{bsf.process.Executable} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  stdin:              {!r}\n'.format(indent, self.stdin))
        str_list.append('{}  stdout:             {!r}\n'.format(indent, self.stdout))
        str_list.append('{}  stderr:             {!r}\n'.format(indent, self.stderr))
        str_list.append('{}  hold:               {!r}\n'.format(indent, self.hold))
        str_list.append('{}  submit:             {!r}\n'.format(indent, self.submit))
        str_list.append('{}  maximum_attempts:   {!r}\n'.format(indent, self.maximum_attempts))
        str_list.append('{}  process_identifier: {!r}\n'.format(indent, self.process_identifier))
        str_list.append('{}  process_name:       {!r}\n'.format(indent, self.process_name))

        # List all dependencies.

        str_list.append('{}  dependencies:\n'.format(indent))

        i = 0
        for dependency in self.dependencies:
            str_list.append('{}    {:2d} {!r}\n'.format(indent, i, dependency))
            i += 1

        # Trace the Command super-class.

        str_list.extend(super(Executable, self).trace(level=level + 1))

        return str_list

    def run(self, debug=0):
        """Run a C{bsf.process.Executable} object via the Python C{subprocess.Popen} class.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        return run_executables(executable_list=[self], debug=debug)


class RunnableStep(Executable):
    """The C{bsf.process.RunnableStep} class represents one C{bsf.process.Executable} object
    in a C{bsf.procedure.Runnable} object.

    @ivar obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
        after successfully completing C{bsf.process.RunnableStep.run}
    @type obsolete_file_path_list: list[str]
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None):
        """Initialise a C{bsf.process.RunnableStep} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        """
        super(RunnableStep, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process)

        if obsolete_file_path_list is None:
            self.obsolete_file_path_list = list()
        else:
            self.obsolete_file_path_list = obsolete_file_path_list

        return

    def trace(self, level=1):
        """Trace a C{bsf.process.RunnableStep} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  obsolete_file_path_list: {!r}\n'.format(indent, self.obsolete_file_path_list))
        str_list.extend(super(RunnableStep, self).trace(level=level + 1))

        return str_list

    def remove_obsolete_file_paths(self):
        """Remove file paths of the C{bsf.process.RunnableStep.obsolete_file_path_list} Python C{list} object.

        This method is mainly used by C{bsf.runnable.consecutive} and related modules.
        """
        if self is None:
            return

        for file_path in self.obsolete_file_path_list:
            if os.path.exists(file_path):
                os.remove(file_path)

        return


class RunnableStepChangeMode(RunnableStep):
    """The C{bsf.process.RunnableStepChangeMode} class represents a step changing file access mode.

    @ivar file_path: File path
    @type file_path: str | None
    @ivar mode_directory: Comma-separated list of C{stat} constant names defining directory access modes
    @type mode_directory: str | None
    @ivar mode_file: Comma-separated list of C{stat} constant names defining file access modes
    @type mode_file: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            file_path=None,
            mode_directory=None,
            mode_file=None):
        """Initialise a C{bsf.process.RunnableStepChangeMode} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param file_path: File path
        @type file_path: str | None
        @param mode_directory: Comma-separated list of C{stat} constant names defining directory access modes
        @type mode_directory: str | None
        @param mode_file: Comma-separated list of C{stat} constant names defining file access modes
        @type mode_file: str | None
        """
        super(RunnableStepChangeMode, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.file_path = file_path
        self.mode_directory = mode_directory
        self.mode_file = mode_file

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepChangeMode} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """

        def convert_mode(mode_str):
            """Private function to convert the mode string into a mode integer.

            @param mode_str: Mode string
            @type mode_str: str
            @return: Mode integer
            @rtype: int
            """
            mode_int = 0
            for permission_str in filter(lambda x: x != '', map(lambda x: x.strip(), mode_str.split(','))):
                # The default of 0 does not modify the mode in case the attribute does not exist.
                mode_int |= getattr(stat, permission_str.upper(), 0)

            return mode_int

        # If the file path is not defined, all further efforts are futile.
        if not self.file_path:
            return None

        # Convert comma-separated directory permission constants to an integer.
        if self.mode_directory is None:
            int_mode_directory = None
        else:
            int_mode_directory = convert_mode(mode_str=self.mode_directory)

        # Convert comma-separated file permission constants to an integer.
        if self.mode_file is None:
            int_mode_file = None
        else:
            int_mode_file = convert_mode(mode_str=self.mode_file)

        # Change the mode of directories and files simultaneously, but only if the mode is not None.

        for file_path, directory_name_list, file_name_list in os.walk(top=self.file_path, topdown=True):
            # Change the mode of the file_path.
            if int_mode_directory is not None:
                os.chmod(file_path, int_mode_directory)
            # Change the mode of each directory name.
            # This is redundant, because each sub-directory will also appear as a file_path once.
            # if int_mode_directory is not None:
            #     for directory_name in directory_name_list:
            #         os.chmod(os.path.join(file_path, directory_name), int_mode_directory)
            # Change the mode of each file name.
            if int_mode_file is not None:
                for file_name in file_name_list:
                    os.chmod(os.path.join(file_path, file_name), int_mode_file)

        return None


class RunnableStepCopy(RunnableStep):
    """The C{bsf.process.RunnableStepCopy} class represents a step copying files.

    @ivar source_path: Source path
    @type source_path: str | None
    @ivar target_path: Target path
    @type target_path: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            source_path=None,
            target_path=None):
        """Initialise a C{bsf.process.RunnableStepCopy} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param source_path: Source path
        @type source_path: str | None
        @param target_path: Target path
        @type target_path: str | None
        """

        super(RunnableStepCopy, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.source_path = source_path
        self.target_path = target_path

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepLink} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """

        if self.source_path and self.target_path and not os.path.exists(self.target_path):
            # Copy data and all stat info ("cp -p src dst").
            shutil.copy2(src=self.source_path, dst=self.target_path)

        return None


class RunnableStepJava(RunnableStep):
    """The C{bsf.process.RunnableStepJava} class represents peculiarities of a Java program.
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_minimum=None,
            java_heap_maximum=None,
            java_jar_path=None):
        """Initialise a C{bsf.process.RunnableStepJava} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | None
        @param java_heap_minimum: Java heap minimum size (-Xms option)
        @type java_heap_minimum: str | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java Archive (JAR) file path
        @type java_jar_path: str | None
        """

        super(RunnableStepJava, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        # JavaVM name
        if not self.name:
            self.name = 'java'

        # JavaVM command
        if not self.program:
            self.program = 'java'

        # JavaVM options

        if 'server' not in self.options:
            self.add_switch_short(key='server')

        if java_heap_minimum and java_heap_minimum not in self.options:
            self.add_switch_short(key=java_heap_minimum)

        if java_heap_maximum and java_heap_maximum not in self.options:
            self.add_switch_short(key=java_heap_maximum)

        if 'XX:+UseG1GC' not in self.options:
            self.add_switch_short(key='XX:+UseG1GC')

        if 'XX:ActiveProcessorCount' not in self.options:
            self.add_option_pair_short(key='XX:ActiveProcessorCount', value='8')

        if 'XX:CICompilerCount' not in self.options:
            self.add_option_pair_short(key='XX:CICompilerCount', value='2')

        if 'XX:ParallelGCThreads' not in self.options:
            self.add_option_pair_short(key='XX:ParallelGCThreads', value='8')

        if java_temporary_path and 'Djava.io.tmpdir' not in self.options:
            self.add_option_pair_short(key='Djava.io.tmpdir', value=java_temporary_path)

        if self.sub_command is None:
            # The Picard command line interface is a bit broken, as the -jar option needs to come last,
            # just before the Picard command. GATK does this slightly better with the --analysis_type option.
            # To be on the safe side, an empty sub command is required to separate the -jar option from the
            # other JavaVM options.
            self.sub_command = Command()
            if java_jar_path:
                self.sub_command.add_option_short(key='jar', value=java_jar_path)

        return


class RunnableStepPicard(RunnableStepJava):
    """The C{bsf.process.RunnableStepPicard} class represents a Picard tool program.
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            java_temporary_path=None,
            java_heap_minimum=None,
            java_heap_maximum=None,
            java_jar_path=None,
            picard_command=None):
        """Initialise a C{bsf.process.RunnableStepPicard} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | None
        @param java_heap_minimum: Java heap minimum size (-Xms option)
        @type java_heap_minimum: str | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java Archive (JAR) file path
        @type java_jar_path: str | None
        @param picard_command: Picard command
        @type picard_command: str | None
        """
        super(RunnableStepPicard, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path,
            java_heap_minimum=java_heap_minimum,
            java_heap_maximum=java_heap_maximum,
            java_jar_path=java_jar_path)

        # The Picard algorithm is then another sub-command.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command(name=picard_command, program=picard_command)

        return

    def add_picard_option(self, key, value, override=False):
        """Add a C{bsf.argument.OptionPair} object to a C{bsf.process.RunnableStepPicard} object.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} objects without warning
        @type override: bool
        """
        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class RunnableStepLink(RunnableStep):
    """The C{bsf.process.RunnableStepLink} class represents a step creating a symbolic link.

    @ivar source_path: Source path
    @type source_path: str | None
    @ivar target_path: Target path
    @type target_path: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            source_path=None,
            target_path=None):
        """Initialise a C{bsf.process.RunnableStepLink} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param source_path: Source path
        @type source_path: str | None
        @param target_path: Target path
        @type target_path: str | None
        """
        super(RunnableStepLink, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.source_path = source_path
        self.target_path = target_path

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepLink} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.source_path and self.target_path and not os.path.exists(self.target_path):
            try:
                os.symlink(self.source_path, self.target_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return None


class RunnableStepMakeDirectory(RunnableStep):
    """The C{bsf.process.RunnableStepMakeDirectory} class represents a step creating a directory.

    @ivar directory_path: Directory path
    @type directory_path: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            directory_path=None):
        """Initialise a C{bsf.process.RunnableStepMakeDirectory} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param directory_path: Directory path
        @type directory_path: str | None
        """
        super(RunnableStepMakeDirectory, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.directory_path = directory_path

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepMakeDirectory} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.directory_path and not os.path.isdir(self.directory_path):
            try:
                os.makedirs(self.directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return None


class RunnableStepMakeNamedPipe(RunnableStep):
    """The C{bsf.process.RunnableStepMakeNamedPipe} class represents a step creating a named pipe.

    @ivar file_path: Named pipe file path
    @type file_path: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            file_path=None):
        """Initialise a C{bsf.process.RunnableStepMakeNamedPipe} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param file_path: Named pipe file path
        @type file_path: str | None
        """
        super(RunnableStepMakeNamedPipe, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.file_path = file_path

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepMakeNamedPipe} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.file_path and not os.path.exists(self.file_path):
            try:
                os.mkfifo(self.file_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return None


class RunnableStepMove(RunnableStep):
    """The C{bsf.process.RunnableStepMove} class represents a step moving a directory or file.

    @ivar source_path: Source path
    @type source_path: str | None
    @ivar target_path: Target path
    @type target_path: str | None
    @ivar use_shutil: Use the Python shutil module rather than a (POSIX) utility
    @type use_shutil: bool | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            source_path=None,
            target_path=None,
            use_shutil=None):
        """Initialise a C{bsf.process.RunnableStepMove} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param source_path: Source path
        @type source_path: str | None
        @param target_path: Target path
        @type target_path: str | None
        @param use_shutil: Use the Python shutil module rather than a (POSIX) utility
        @type use_shutil: bool | None
        """
        super(RunnableStepMove, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.source_path = source_path
        self.target_path = target_path
        self.use_shutil = use_shutil

        if not self.use_shutil:
            self.program = 'mv'
            self.arguments = [self.source_path, self.target_path]

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepMove} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.source_path and self.target_path:
            if self.use_shutil:
                shutil.move(src=self.source_path, dst=self.target_path)
                return None
            else:
                return super(RunnableStepMove, self).run(debug=debug)


class RunnableStepRemoveDirectory(RunnableStep):
    """The C{bsf.process.RunnableStepRemoveDirectory} class represents a step removing a directory.

    @ivar directory_path: Directory path
    @type directory_path: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            directory_path=None):
        """Initialise a C{bsf.process.RunnableStepRemoveDirectory} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param directory_path: Directory path
        @type directory_path: str | None
        """
        super(RunnableStepRemoveDirectory, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.directory_path = directory_path

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableRemoveDirectory} object.

        FileNotFoundError and OSError ENOTEMPTY Exceptions are caught.
        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.directory_path and not os.path.isdir(self.directory_path):
            try:
                os.rmdir(self.directory_path)
            except FileNotFoundError:
                pass
            except OSError as exception:
                if exception.errno != errno.ENOTEMPTY:
                    raise

        return None


class RunnableStepRemoveFile(RunnableStep):
    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            file_path=None):
        """Initialise a C{bsf.process.RunnableStepRemoveFile} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param file_path: File path
        @type file_path: str | None
        """
        super(RunnableStepRemoveFile, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.file_path = file_path

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepRemoveFile} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        os.remove(path=self.file_path)

        return None


class RunnableStepRemoveTree(RunnableStep):
    """The C{bsf.process.RunnableStepRemoveTree} class represents a step removing a file system tree.

    @ivar file_path: File path
    @type file_path: str | None
    @ivar ignore_errors: Ignore errors
    @type ignore_errors: bool
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            file_path=None,
            ignore_errors=None):
        """Initialise a C{bsf.process.RunnableStepRemoveTree} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param file_path: File path
        @type file_path: str | None
        """
        super(RunnableStepRemoveTree, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.file_path = file_path
        self.ignore_errors = ignore_errors

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableRemoveTree} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.ignore_errors:
            ignore_errors = True
        else:
            ignore_errors = False

        shutil.rmtree(path=self.file_path, ignore_errors=ignore_errors)

        return None


class RunnableStepSleep(RunnableStep):
    """The C{bsf.process.RunnableStepSleep} class represents a step sleeping the process.

    @ivar sleep_time: Sleep time in seconds
    @type sleep_time: float | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            sleep_time=None):
        """Initialise a C{bsf.process.RunnableStepSleep} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str]
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param sleep_time: Sleep time in seconds
        @type sleep_time: float | None
        """
        super(RunnableStepSleep, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.sleep_time = sleep_time

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepSleep} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.sleep_time is not None:
            time.sleep(self.sleep_time)

        return None


class RunnableStepSetEnvironment(RunnableStep):
    """The C{bsf.process.RunnableStepSetEnvironment} class represents a step setting the process environment.

    @ivar key: Environment key
    @type key: str
    @ivar value: Environment value
    @type value: str
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            key=None,
            value=None):
        """Initialise a C{bsf.process.RunnableStepSetEnvironment} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str]
        @param sub_command: Subordinate C{bsf.process.Command} object
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector} object
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector} object
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector} object
        @type stderr: Connector | None
        @param dependencies: Python C{list} of Python C{str} (C{bsf.process.Executable.name}) objects
            in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit this C{bsf.process.Executable} object during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen} object
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of Python C{str} file path objects that can be removed
            after successfully completing C{bsf.process.RunnableStep.run}
        @type obsolete_file_path_list: list[str] | None
        @param key: Environment key
        @type key: str
        """
        super(RunnableStepSetEnvironment, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.key = key
        self.value = value

        return

    def run(self, debug=0):
        """Run a C{bsf.process.RunnableStepSetEnvironment} object.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        if self.key is not None:
            os.environ[self.key] = self.value

        return None
