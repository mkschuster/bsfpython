# -*- coding: utf-8 -*-
"""Process module.

A package of classes and methods modelling processes.
"""
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

from bsf.argument import *
from bsf.connector import *
from bsf.standards import Configuration


def get_timestamp():
    """Get the current time stamp in ISO format.

    @return: ISO format time stamp
    @rtype: str
    """
    return '[' + datetime.datetime.now().isoformat() + ']'


class Command(object):
    """The C{bsf.process.Command} class represents one program, its options and arguments.

    A C{bsf.process.Command} can possibly contain another subordinate C{bsf.process.Command}.

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar program: Program
    @type program: str
    @ivar options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
        Python C{list} value objects of C{bsf.argument.Argument} objects
    @type options: dict[Argument.key, list[Argument]]
    @ivar arguments: Python C{list} of Python C{str} (program argument) objects
    @type arguments: list[str]
    @ivar sub_command: Subordinate C{bsf.process.Command}
    @type sub_command: Command | None
    """

    def __init__(
            self,
            name=None,
            program=None,
            options=None,
            arguments=None,
            sub_command=None):
        """Initialise a C{bsf.process.Command}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
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
        """Trace a C{bsf.process.Command}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

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
        """Add a C{bsf.argument.Argument} or one of its sub-classes.

        @param argument: C{bsf.argument.Argument}
        @type argument: Argument
        @param override: Override existing C{bsf.argument.Argument} without warning
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
        """Initialise and add a C{bsf.argument.SwitchLong}.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=SwitchLong(key=key), override=override)

    def add_switch_short(self, key, override=False):
        """Initialise and add a C{bsf.argument.SwitchShort}.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=SwitchShort(key=key), override=override)

    def add_option_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionLong(key=key, value=value), override=override)

    def add_option_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionShort(key=key, value=value), override=override)

    def add_option_pair(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionPair(key=key, value=value), override=override)

    def add_option_pair_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionPairShort(key=key, value=value), override=override)

    def add_option_pair_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionPairLong(key=key, value=value), override=override)

    def add_option_multi(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMulti}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMulti(key=key, value=value), override=override)

    def add_option_multi_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiLong(key=key, value=value), override=override)

    def add_option_multi_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiShort(key=key, value=value), override=override)

    def add_option_multi_pair(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiPair(key=key, value=value), override=override)

    def add_option_multi_pair_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiPairLong(key=key, value=value), override=override)

    def add_option_multi_pair_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.add_argument(argument=OptionMultiPairShort(key=key, value=value), override=override)

    def set_argument(self, argument, override):
        """Set a C{bsf.argument.Argument} or one of its sub-classes.

        @param argument: C{bsf.argument.Argument}
        @type argument: Argument
        @param override: Override existing C{bsf.argument.Argument} without warning
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
        """Initialise and set a C{bsf.argument.SwitchLong}.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=SwitchLong(key=key), override=override)

    def set_switch_short(self, key, override=False):
        """Initialise and set a C{bsf.argument.SwitchShort}.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=SwitchShort(key=key), override=override)

    def set_option_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionLong(key=key, value=value), override=override)

    def set_option_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionShort(key=key, value=value), override=override)

    def set_option_pair(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_option_pair_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionPairShort(key=key, value=value), override=override)

    def set_option_pair_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionPairLong(key=key, value=value), override=override)

    def set_option_multi(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMulti}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMulti(key=key, value=value), override=override)

    def set_option_multi_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiLong(key=key, value=value), override=override)

    def set_option_multi_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiShort(key=key, value=value), override=override)

    def set_option_multi_pair(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiPair(key=key, value=value), override=override)

    def set_option_multi_pair_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiPairLong(key=key, value=value), override=override)

    def set_option_multi_pair_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.set_argument(argument=OptionMultiPairShort(key=key, value=value), override=override)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.process.Command} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
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
        command_line = list()
        """ @type command_line: list[str] """

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

    Attributes:
    @ivar stdin: Standard input I{STDIN} C{bsf.connector.Connector}
    @type stdin: Connector | None
    @ivar stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
    @type stdout: Connector | None
    @ivar stderr: Standard error I{STDERR} C{bsf.connector.Connector}
    @type stderr: Connector | None
    @ivar dependencies: Python C{list} of C{bsf.process.Executable.name} properties in the
        context of C{bsf.analysis.Stage} dependencies
    @type dependencies: list[Executable.name]
    @ivar hold: Hold on job scheduling
    @type hold: str | None
    @ivar submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
    @type submit: bool
    @ivar maximum_attempts: Maximum number of attempts to run this C{bsf.process.Executable}
    @type maximum_attempts: int
    @ivar process_identifier: Process identifier
    @type process_identifier: str | None
    @ivar process_name: Process name
    @type process_name: str | None
    @ivar sub_process: C{subprocess.Popen}
    @type sub_process: Popen | None
    """

    @staticmethod
    def process_stream(file_handle, thread_lock, debug, file_type, file_path=None):
        """Process a I{STDOUT} or I{STDERR} text stream from the child process as a thread.

        If a file_path was provided, a corresponding Python C{file} will be opened in text mode,
        if not, C{sys.stdout} or C{sys.stderr} will be used according to the I{file_type}.
        If a debug level was set, diagnostic output will be printed to C{sys.stdout}, as well as the output stream.
        @param file_handle: The I{STDOUT} or I{STDERR} C{io.TextIOWrapper} file handle
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

        @param stdout_handle: The I{STDOUT} C{io.TextIOWrapper} file handle
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

        @param stderr_handle: The I{STDERR} C{io.TextIOWrapper} file handle
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
        """Initialise a C{bsf.process.Executable}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param maximum_attempts: Maximum number of attempts to run this C{bsf.process.Executable}
        @type maximum_attempts: int
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
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
        """Trace a C{bsf.process.Executable}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.Executable} via the Python C{subprocess.Popen} class.

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

        def _map_connector(_connector):
            """Map a connector to a file handle.

            @param _connector: C{bsf.connector.Connector} or sub-class thereof.
            @type _connector: Connector
            @return: File handle
            @rtype: IOBase | PIPE | DEVNULL
            """
            if isinstance(_connector, ElectronicSink):
                return DEVNULL

            if isinstance(_connector, ConnectorFile):
                if isinstance(_connector, ConnectorPipeNamed):
                    # A named pipe needs creating before it can be opened.
                    if not os.path.exists(_connector.file_path):
                        os.mkfifo(_connector.file_path)
                return open(file=_connector.file_path, mode=_connector.file_mode)

            if isinstance(_connector, ConnectorPipe):
                return PIPE

            if isinstance(_connector, StandardStream):
                return PIPE

            if isinstance(_connector, ConcurrentProcess):
                # A bsf.connector.ConcurrentProcess object
                raise Exception('A ' + type(_connector).__name__ +
                                ' object cannot be set for a ' + type(self).__name__ +
                                ' object.')

            return None

        # Up to three threads, writing to STDIN, as well as reading from STDOUT and STDERR,
        # should make sure that buffers are not filling up. If STDOUT or STDERR connectors are not defined,
        # the defaults need be set to avoid the process from blocking.

        if self.stdout is None:
            self.stdout = StandardOutputStream()

        if self.stderr is None:
            self.stderr = StandardErrorStream()

        child_return_code = 0
        attempt_counter = 0

        thread_lock = Lock()

        while attempt_counter < self.maximum_attempts:
            attempt_counter += 1

            try:
                self.sub_process = Popen(
                    args=self.command_list(),
                    bufsize=0,
                    stdin=_map_connector(_connector=self.stdin),
                    stdout=_map_connector(_connector=self.stdout),
                    stderr=_map_connector(_connector=self.stderr),
                    close_fds='posix' in sys.builtin_module_names,
                    shell=False,
                    text=True)
            except OSError as exception:
                if exception.errno == errno.ENOENT:
                    raise Exception(repr(self) + ' ' + repr(self.program) + ' could not be found.')
                else:
                    # Re-raise the Exception object.
                    raise exception

            # Start threads for standard stream processing to prevent buffers from filling up and
            # sub-processes from blocking.

            for attribute in ('stdin', 'stdout', 'stderr'):
                connector = getattr(self, attribute)
                """ @type connector: Connector """

                if isinstance(connector, StandardInputStream):
                    if connector.thread_callable is None:
                        # If a specific STDIN callable is not defined, run bsf.process.Executable.process_stdin().
                        pass
                    else:
                        connector.thread = Thread(
                            target=connector.thread_callable,
                            args=[self.sub_process.stdin, thread_lock, debug],
                            kwargs=connector.thread_kwargs)

                if isinstance(connector, StandardOutputStream):
                    if connector.thread_callable is None:
                        # If a specific STDOUT callable is not defined, run bsf.process.Executable.process_stdout().
                        connector.thread = Thread(
                            target=Executable.process_stdout,
                            args=[self.sub_process.stdout, thread_lock, debug],
                            kwargs={'stdout_path': connector.file_path})
                    else:
                        connector.thread = Thread(
                            target=connector.thread_callable,
                            args=[self.sub_process.stdout, thread_lock, debug],
                            kwargs=connector.thread_kwargs)

                if isinstance(connector, StandardErrorStream):
                    if connector.thread_callable is None:
                        # If a specific STDERR callable is not defined, run bsf.process.Executable.process_stderr().
                        connector.thread = Thread(
                            target=Executable.process_stderr,
                            args=[self.sub_process.stderr, thread_lock, debug],
                            kwargs={'stderr_path': connector.file_path})
                    else:
                        connector.thread = Thread(
                            target=connector.thread_callable,
                            args=[self.sub_process.stderr, thread_lock, debug],
                            kwargs=connector.thread_kwargs)

                if isinstance(connector, StandardStream) and connector.thread:
                    connector.thread.daemon = True
                    connector.thread.start()

            child_return_code = self.sub_process.wait()

            # First, join all standard stream processing threads.

            for attribute in ('stdin', 'stdout', 'stderr'):
                connector = getattr(self, attribute)
                """ @type connector: Connector """
                if isinstance(connector, StandardStream) and connector.thread:
                    thread_join_counter = 0
                    while connector.thread.is_alive() and thread_join_counter < connector.thread_joins:
                        if debug > 0:
                            thread_lock.acquire(True)
                            print(get_timestamp(),
                                  'Waiting for ' + repr(attribute) + ' processor to finish.')
                            thread_lock.release()

                        connector.thread.join(timeout=connector.thread_timeout)
                        thread_join_counter += 1

            # Second, inspect the child process' return code.

            if child_return_code > 0:
                if debug > 0:
                    print(get_timestamp(),
                          'Child process ' + repr(self.name) +
                          ' failed with exit code ' + repr(+child_return_code) + '.')
            elif child_return_code < 0:
                if debug > 0:
                    print(get_timestamp(),
                          'Child process ' + repr(self.name) +
                          ' received signal ' + repr(-child_return_code) + '.')
            else:
                if debug > 0:
                    print(get_timestamp(),
                          'Child process ' + repr(self.name) +
                          ' completed successfully ' + repr(+child_return_code) + '.')
                break
        else:
            if debug > 0:
                print(get_timestamp(),
                      'Runnable ' + repr(self.name) +
                      ' exceeded the maximum retry counter ' + repr(self.maximum_attempts) + '.')

        return child_return_code

    def evaluate_return_code(self, return_code):
        """Evaluate a return code from the run method.

        @param return_code: Return code
        @type return_code: int
        """
        if return_code > 0:
            print(get_timestamp(),
                  'Child process ' + repr(self.name) +
                  ' failed with return code ' + repr(+return_code) + '.')
        elif return_code < 0:
            print(get_timestamp(),
                  'Child process ' + repr(self.name) +
                  ' received signal ' + repr(-return_code) + '.')
        else:
            print(get_timestamp(),
                  'Child process ' + repr(self.name) +
                  ' completed with return code ' + repr(+return_code) + '.')

        return


class RunnableStep(Executable):
    """The C{bsf.process.RunnableStep} represents one C{bsf.process.Executable} in a C{bsf.procedure.Runnable}.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{bsf.process.RunnableStep}
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
        """Initialise a C{bsf.process.RunnableStep}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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
        """Trace a C{bsf.process.RunnableStep}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  obsolete_file_path_list: {!r}\n'.format(indent, self.obsolete_file_path_list))
        str_list.extend(super(RunnableStep, self).trace(level=level + 1))

        return str_list

    def remove_obsolete_file_paths(self):
        """Remove file paths on the C{bsf.process.RunnableStep.obsolete_file_path_list} Python C{list}.

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

    Attributes:
    @ivar file_path: File path
    @type file_path: str | None
    @ivar mode_directory: Directory access mode according to C{stat}
    @type mode_directory: str | None
    @ivar mode_file: File access mode for files according to C{stat}
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
        """Initialise a C{bsf.process.RunnableStepChangeMode}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param file_path: File path
        @type file_path: str | None
        @param mode_directory: Directory access mode according to C{stat}
        @type mode_directory: str | None
        @param mode_file: File access mode for files according to C{stat}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepChangeMode}.

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

        def convert_mode(mode_str):
            """Private function to convert the mode string into a mode integer.

            @param mode_str: Mode string
            @type mode_str: str
            @return: Mode integer
            @rtype: int
            """
            mode_int = 0
            for permission_bit in filter(
                    lambda x: x != '',
                    map(
                        lambda x: x.strip(),
                        mode_str.split(','))):
                if permission_bit.upper() in permission_dict:
                    mode_int |= permission_dict[permission_bit.upper()]

            return mode_int

        # If the file path is not defined, all further efforts are futile.
        if not self.file_path:
            return 0

        # Use a dictionary to map string literals to integers defined in the stat module rather than
        # evaluating code directly, which can be rather dangerous.

        permission_dict = {
            'S_ISUID': stat.S_ISUID,
            'S_ISGID': stat.S_ISGID,
            # System V file locking enforcement.
            'S_ENFMT': stat.S_ENFMT,
            # Sticky bit.
            'S_ISVTX': stat.S_ISVTX,
            'S_IREAD': stat.S_IREAD,
            'S_IWRITE': stat.S_IWRITE,
            # UNIX V7 synonyms.
            'S_IEXEC': stat.S_IEXEC,
            'S_IRWXU': stat.S_IRWXU,
            # User permissions.
            'S_IRUSR': stat.S_IRUSR,
            'S_IWUSR': stat.S_IWUSR,
            'S_IXUSR': stat.S_IXUSR,
            'S_IRWXG': stat.S_IRWXG,
            # Group permissions.
            'S_IRGRP': stat.S_IRGRP,
            'S_IWGRP': stat.S_IWGRP,
            'S_IXGRP': stat.S_IXGRP,
            'S_IRWXO': stat.S_IRWXO,
            # Other permissions.
            'S_IROTH': stat.S_IROTH,
            'S_IWOTH': stat.S_IWOTH,
            'S_IXOTH': stat.S_IXOTH,
        }

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

        return 0


class RunnableStepCopy(RunnableStep):
    """The C{bsf.process.RunnableStepCopy} represents a step copying files.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepCopy}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepLink}.

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
            # Copy data and all stat info ("cp -p src dst").
            shutil.copy2(src=self.source_path, dst=self.target_path)

        return 0


class RunnableStepJava(RunnableStep):
    """The C{bsf.process.RunnableStepJava} class represents peculiarities of a Java program.

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
        """Initialise a C{bsf.process.RunnableStepJava}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | None
        @param java_heap_minimum: Java heap minimum size (-Xms option)
        @type java_heap_minimum: str | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java archive file path
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
        if 'd64' not in self.options:
            self.add_switch_short(key='d64')

        if 'server' not in self.options:
            self.add_switch_short(key='server')

        if java_heap_minimum and java_heap_minimum not in self.options:
            self.add_switch_short(key=java_heap_minimum)

        if java_heap_maximum and java_heap_maximum not in self.options:
            self.add_switch_short(key=java_heap_maximum)

        if java_temporary_path and '-Djava.io.tmpdir' not in self.options:
            self.add_option_pair(key='-Djava.io.tmpdir', value=java_temporary_path)

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
            picard_classpath=None,
            picard_command=None):
        """Initialise a C{bsf.process.RunnableStepPicard}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | None
        @param java_heap_minimum: Java heap minimum size (-Xms option)
        @type java_heap_minimum: str | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | None
        @param picard_classpath: Picard class path
        @type picard_classpath: str | None
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

        # Set the Picard classpath and the Picard Java archive.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(picard_classpath, 'picard.jar'))

        # The Picard algorithm is then another sub-command.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command(name=picard_command, program=picard_command)

        return

    def add_picard_option(self, key, value, override=False):
        """Add a C{bsf.argument.OptionPair} to a C{bsf.process.RunnableStepPicard}.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        """
        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class RunnableStepLink(RunnableStep):
    """The C{bsf.process.RunnableStepLink} represents a step creating a symbolic link.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepLink}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepLink}.

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
    """The C{bsf.process.RunnableStepMakeDirectory} represents a step creating a directory.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepMakeDirectory}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepMakeDirectory}.

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


class RunnableStepMakeNamedPipe(RunnableStep):
    """The C{bsf.process.RunnableStepMakeNamedPipe} represents a step creating a named pipe.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepMakeNamedPipe}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepMakeNamedPipe}.

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
        if self.file_path and not os.path.exists(self.file_path):
            try:
                os.mkfifo(self.file_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return 0


class RunnableStepMove(RunnableStep):
    """The C{bsf.process.RunnableStepMove} class represents a step moving a directory or file.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepMove}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param source_path: Source path
        @type source_path: str | None
        @param target_path: Target path
        @type target_path: str | None
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

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepMove}.

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
            shutil.move(src=self.source_path, dst=self.target_path)

        return 0


class RunnableStepSleep(RunnableStep):
    """The C{bsf.process.RunnableStepSleep} class represents a step sleeping the process.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepSleep}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str]
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepSleep}.

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
        if self.sleep_time is not None:
            time.sleep(self.sleep_time)

        return 0


class RunnableStepSetEnvironment(RunnableStep):
    """The C{bsf.process.RunnableStepSetEnvironment} class represents a step setting the process environment.

    Attributes:
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
        """Initialise a C{bsf.process.RunnableStepSetEnvironment}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[Argument.key, list[Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str]
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.process.RunnableStepSetEnvironment}.

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
        if self.key is not None:
            os.environ[self.key] = self.value

        return 0
