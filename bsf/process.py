# -*- coding: utf-8 -*-
"""Process module

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
from __future__ import print_function

import datetime
import errno
import os
import shutil
import stat
import subprocess
import sys
import threading
import time
import warnings

import bsf.argument
import bsf.standards


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
    @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]]
    @ivar arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
    @type arguments: list[str | unicode]
    @ivar sub_command: Subordinate C{bsf.process.Command}
    @type sub_command: bsf.process.Command | None
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
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @return:
        @rtype:
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
        @rtype: list[str | unicode]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str | unicode] """

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
        @type argument: bsf.argument.Argument
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        assert isinstance(argument, bsf.argument.Argument)
        assert isinstance(override, bool)

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
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.SwitchLong(key=key), override=override)

    def add_switch_short(self, key, override=False):
        """Initialise and add a C{bsf.argument.SwitchShort}.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.SwitchShort(key=key), override=override)

    def add_option_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionLong(key=key, value=value), override=override)

    def add_option_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionShort(key=key, value=value), override=override)

    def add_option_pair(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionPair(key=key, value=value), override=override)

    def add_option_pair_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionPairShort(key=key, value=value), override=override)

    def add_option_pair_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionPairLong(key=key, value=value), override=override)

    def add_option_multi(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMulti}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionMulti(key=key, value=value), override=override)

    def add_option_multi_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionMultiLong(key=key, value=value), override=override)

    def add_option_multi_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionMultiShort(key=key, value=value), override=override)

    def add_option_multi_pair(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionMultiPair(key=key, value=value), override=override)

    def add_option_multi_pair_long(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionMultiPairLong(key=key, value=value), override=override)

    def add_option_multi_pair_short(self, key, value, override=False):
        """Initialise and add a C{bsf.argument.OptionMultiPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.add_argument(argument=bsf.argument.OptionMultiPairShort(key=key, value=value), override=override)

    def set_argument(self, argument, override):
        """Set a C{bsf.argument.Argument} or one of its sub-classes.

        @param argument: C{bsf.argument.Argument}
        @type argument: bsf.argument.Argument
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        assert isinstance(argument, bsf.argument.Argument)
        assert isinstance(override, bool)

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
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.SwitchLong(key=key), override=override)

    def set_switch_short(self, key, override=False):
        """Initialise and set a C{bsf.argument.SwitchShort}.

        @param key: Key
        @type key: str
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.SwitchShort(key=key), override=override)

    def set_option_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionLong(key=key, value=value), override=override)

    def set_option_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionShort(key=key, value=value), override=override)

    def set_option_pair(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionPair(key=key, value=value), override=override)

    def set_option_pair_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionPairShort(key=key, value=value), override=override)

    def set_option_pair_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionPairLong(key=key, value=value), override=override)

    def set_option_multi(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMulti}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionMulti(key=key, value=value), override=override)

    def set_option_multi_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionMultiLong(key=key, value=value), override=override)

    def set_option_multi_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionMultiShort(key=key, value=value), override=override)

    def set_option_multi_pair(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPair}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionMultiPair(key=key, value=value), override=override)

    def set_option_multi_pair_long(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPairLong}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionMultiPairLong(key=key, value=value), override=override)

    def set_option_multi_pair_short(self, key, value, override=False):
        """Initialise and set a C{bsf.argument.OptionMultiPairShort}.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{bsf.argument.Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=bsf.argument.OptionMultiPairShort(key=key, value=value), override=override)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.process.Command} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section, defaults to instance class
        @type section: str
        @return:
        @rtype:
        """
        assert isinstance(configuration, bsf.standards.Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            warnings.warn(
                'Section ' + repr(section) + ' not defined in configuration files:\n' +
                repr(configuration.from_file_path_list),
                UserWarning)

            return

        # The configuration section is available.

        for option in configuration.config_parser.options(section=section):
            self.add_argument(
                argument=bsf.argument.Argument.from_key_value(
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
        @rtype: str | unicode
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
    @ivar stdin_callable: Standard input I{STDIN} callable
    @type stdin_callable: object | None
    @ivar stdin_kwargs: Standard input I{STDIN} keyword arguments
    @type stdin_kwargs: dict[str, object]
    @ivar stdout_callable: Standard output I{STDOUT} callable
    @type stdout_callable: object | None
    @ivar stdout_kwargs: Standard output I{STDOUT} keyword arguments
    @type stdout_kwargs: dict[str, object]
    @ivar stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
    @type stdout_path: str | unicode | None
    @ivar stderr_callable: Standard error I{STDERR} callable
    @type stderr_callable: object | None
    @ivar stderr_kwargs: Standard error I{STDERR} keyword arguments
    @type stderr_kwargs: dict[str, object]
    @ivar stderr_path: Standard error (STDERR) redirection in Bash (2>word)
    @type stderr_path: str | unicode | None
    @ivar dependencies: Python C{list} of C{bsf.process.Executable.name} properties in the
        context of C{bsf.Stage} dependencies
    @type dependencies: list[bsf.process.Executable.name]
    @ivar hold: Hold on job scheduling
    @type hold: str | None
    @ivar submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
    @type submit: bool
    @ivar maximum_attempts: Maximum number of attempts to run this C{bsf.process.Executable}
    @type maximum_attempts: int
    @ivar process_identifier: Process identifier
    @type process_identifier: str | None
    @ivar process_name: Process name
    @type process_name: str | None
    """

    @staticmethod
    def process_stream(file_handle, thread_lock, debug, file_type, file_path=None):
        """Process a I{STDOUT} or I{STDERR} text stream from the child process as a thread.

        If a file_path was provided, a corresponding Python C{file} will be opened in text mode,
        if not C{sys.stdout} or C{sys.stderr} will be used according to the I{file_type}.
        If a debug level was set, diagnostic output will be printed to C{sys.stdout}, as well as the output stream.
        @param file_handle: The I{STDOUT} or I{STDERR} file handle
        @type file_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: threading.Lock
        @param debug: Debug level
        @type debug: int
        @param file_type: File handle type I{STDOUT} or I{STDERR}
        @type file_type: str
        @param file_path: I{STDOUT} file path
        @type file_path: str | unicode | None
        @return:
        @rtype:
        @raise Exception: If file_type is neither I{STDOUT} nor I{STDERR}
        """
        def get_timestamp():
            return '[' + datetime.datetime.now().isoformat() + ']'

        if file_type not in ('STDOUT', 'STDERR'):
            raise Exception('The file_type has to be either STDOUT or STDERR.')

        thread_lock.acquire(True)
        if debug > 0:
            print(get_timestamp(),
                  'Started Runner ' + repr(file_type) + ' processor in module ' + repr(__name__) + '.',
                  file=sys.stdout)
        output_file = None
        if file_path:
            output_file = open(file_path, 'wt')
            if debug > 0:
                print(get_timestamp(), 'Opened ' + repr(file_type) + ' file ' + repr(file_path) + '.', file=sys.stdout)
        elif file_type == 'STDOUT':
            output_file = sys.stdout
        elif file_type == 'STDERR':
            output_file = sys.stderr
        thread_lock.release()

        for line_str in file_handle:
            thread_lock.acquire(True)
            if debug > 0:
                print(get_timestamp(), file_type + ': ' + line_str.rstrip(), file=output_file)
            else:
                output_file.write(line_str)
            thread_lock.release()

        thread_lock.acquire(True)
        if debug > 0:
            print(get_timestamp(), 'Received EOF on ' + repr(file_type) + ' pipe.', file=sys.stdout)
        if file_path:
            output_file.close()
            if debug > 0:
                print(get_timestamp(), 'Closed ' + repr(file_type) + ' file ' + repr(file_path) + '.', file=sys.stdout)
        thread_lock.release()

        return

    @staticmethod
    def process_stdout(stdout_handle, thread_lock, debug, stdout_path=None):
        """Process I{STDOUT} from the child process as a thread.

        @param stdout_handle: The I{STDOUT} file handle
        @type stdout_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: threading.Lock
        @param debug: Debug level
        @type debug: int
        @param stdout_path: I{STDOUT} file path
        @type stdout_path: str | unicode | None
        @return:
        @rtype:
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

        @param stderr_handle: The I{STDERR} file handle
        @type stderr_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: threading.Lock
        @param debug: Debug level
        @type debug: int
        @param stderr_path: I{STDERR} file path
        @type stderr_path: str | unicode | None
        @return:
        @rtype:
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
            stdin_callable=None,
            stdin_kwargs=None,
            stdout_callable=None,
            stdout_kwargs=None,
            stdout_path=None,
            stderr_callable=None,
            stderr_kwargs=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            maximum_attempts=1,
            process_identifier=None,
            process_name=None):
        """Initialise a C{bsf.process.Executable}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdin_callable: Standard input I{STDIN} callable
        @type stdin_callable: object | None
        @param stdin_kwargs: Standard input I{STDIN} keyword arguments
        @type stdin_kwargs: dict[str, object] | None
        @param stdout_callable: Standard output I{STDOUT} callable
        @type stdout_callable: object | None
        @param stdout_kwargs: Standard output I{STDOUT} keyword arguments
        @type stdout_kwargs: dict[str, object] | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_callable: Standard error I{STDERR} callable
        @type stderr_callable: object | None
        @param stderr_kwargs: Standard error I{STDERR} keyword arguments
        @type stderr_kwargs: dict[str, object] | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param maximum_attempts: Maximum number of attempts to run this C{bsf.process.Executable}
        @type maximum_attempts: int
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @return:
        @rtype:
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

        self.stdin_callable = stdin_callable
        self.stdin_kwargs = stdin_kwargs

        self.stdout_callable = stdout_callable
        self.stdout_kwargs = stdout_kwargs
        self.stdout_path = stdout_path

        self.stderr_callable = stderr_callable
        self.stderr_kwargs = stderr_kwargs
        self.stderr_path = stderr_path

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

        return

    def trace(self, level):
        """Trace a C{bsf.process.Executable}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str | unicode]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  stdin_callable:     {!r}\n'.format(indent, self.stdin_callable))
        str_list.append('{}  stdout_callable:    {!r}\n'.format(indent, self.stdout_callable))
        str_list.append('{}  stdout_path:        {!r}\n'.format(indent, self.stdout_path))
        str_list.append('{}  stderr_callable:    {!r}\n'.format(indent, self.stderr_callable))
        str_list.append('{}  stderr_path:        {!r}\n'.format(indent, self.stderr_path))
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
        on_posix = 'posix' in sys.builtin_module_names

        child_return_code = 0
        attempt_counter = 0

        while attempt_counter < self.maximum_attempts:
            attempt_counter += 1

            try:
                child_process = subprocess.Popen(
                    args=self.command_list(),
                    bufsize=0,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=False,
                    close_fds=on_posix)
            except OSError as exception:
                if exception.errno == errno.ENOENT:
                    raise Exception('Executable ' + repr(self.program) + ' could not be found.')
                else:
                    # Re-raise the Exception object.
                    raise exception

            # Up to three threads, writing to STDIN, as well as reading from STDOUT and STDERR,
            # should make sure that buffers are not filling up.

            thread_lock = threading.Lock()

            # Start a STDIN thread.

            if self.stdin_callable is None:
                thread_in = None
            else:
                thread_in = threading.Thread(
                    target=self.stdin_callable,
                    args=[child_process.stdin, thread_lock, debug],
                    kwargs=self.stdin_kwargs)
                thread_in.daemon = True  # Thread dies with the program.
                thread_in.start()

            # Start a STDOUT thread.

            if self.stdout_callable is None:
                # If a specific STDOUT callable is not defined, run bsf.process.Executable.process_stdout().
                thread_out = threading.Thread(
                    target=bsf.process.Executable.process_stdout,
                    args=[child_process.stdout, thread_lock, debug],
                    kwargs={'stdout_path': self.stdout_path})
            else:
                thread_out = threading.Thread(
                    target=self.stdout_callable,
                    args=[child_process.stdout, thread_lock, debug],
                    kwargs=self.stdout_kwargs)

            thread_out.daemon = True  # Thread dies with the program.
            thread_out.start()

            # Start a STDERR thread.

            if self.stderr_callable is None:
                # If a specific STDERR callable is not defined, run bsf.process.Executable.process_stderr().
                thread_err = threading.Thread(
                    target=bsf.process.Executable.process_stderr,
                    args=[child_process.stderr, thread_lock, debug],
                    kwargs={'stderr_path': self.stderr_path})
            else:
                thread_err = threading.Thread(
                    target=self.stderr_callable,
                    args=[child_process.stderr, thread_lock, debug],
                    kwargs=self.stderr_kwargs)

            thread_err.daemon = True  # Thread dies with the program.
            thread_err.start()

            # Wait for the child process to finish.

            child_return_code = child_process.wait()

            thread_join_counter = 0

            if thread_in is not None:
                while thread_in.is_alive() and thread_join_counter < max_thread_joins:
                    if debug > 0:
                        thread_lock.acquire(True)
                        print('[' + datetime.datetime.now().isoformat() + ']',
                              'Waiting for STDIN processor to finish.')
                        thread_lock.release()

                    thread_in.join(timeout=thread_join_timeout)
                    thread_join_counter += 1

            thread_join_counter = 0

            while thread_out.is_alive() and thread_join_counter < max_thread_joins:
                if debug > 0:
                    thread_lock.acquire(True)
                    print('[' + datetime.datetime.now().isoformat() + ']',
                          'Waiting for STDOUT processor to finish.')
                    thread_lock.release()

                thread_out.join(timeout=thread_join_timeout)
                thread_join_counter += 1

            thread_join_counter = 0

            while thread_err.is_alive() and thread_join_counter < max_thread_joins:
                if debug > 0:
                    thread_lock.acquire(True)
                    print('[' + datetime.datetime.now().isoformat() + ']',
                          'Waiting for STDERR processor to finish.')
                    thread_lock.release()

                thread_err.join(timeout=thread_join_timeout)
                thread_join_counter += 1

            if child_return_code > 0:
                if debug > 0:
                    print('[' + datetime.datetime.now().isoformat() + ']',
                          'Child process ' + repr(self.name) +
                          ' failed with exit code ' + repr(+child_return_code) + '.')
            elif child_return_code < 0:
                if debug > 0:
                    print('[' + datetime.datetime.now().isoformat() + ']',
                          'Child process ' + repr(self.name) +
                          ' received signal ' + repr(-child_return_code) + '.')
            else:
                if debug > 0:
                    print('[' + datetime.datetime.now().isoformat() + ']',
                          'Child process ' + repr(self.name) +
                          ' completed successfully ' + repr(+child_return_code) + '.')
                break
        else:
            if debug > 0:
                print('[' + datetime.datetime.now().isoformat() + ']',
                      'Runnable ' + repr(self.name) +
                      ' exceeded the maximum retry counter ' + repr(self.maximum_attempts) + '.')

        return child_return_code

    def evaluate_return_code(self, return_code):
        """Evaluate a return code from the run method.

        @param return_code: Return code
        @type return_code: int
        @return:
        @rtype:
        """
        if return_code > 0:
            print('[' + datetime.datetime.now().isoformat() + ']',
                  'Child process ' + repr(self.name) +
                  ' failed with return code ' + repr(+return_code) + '.')
        elif return_code < 0:
            print('[' + datetime.datetime.now().isoformat() + ']',
                  'Child process ' + repr(self.name) +
                  ' received signal ' + repr(-return_code) + '.')
        else:
            print('[' + datetime.datetime.now().isoformat() + ']',
                  'Child process ' + repr(self.name) +
                  ' completed with return code ' + repr(+return_code) + '.')

        return


class RunnableStep(Executable):
    """The C{bsf.process.RunnableStep} represents one C{bsf.process.Executable} in a C{bsf.procedure.Runnable}.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{bsf.process.RunnableStep}
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
        """Initialise a C{bsf.process.RunnableStep}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
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
        """Trace a C{bsf.process.RunnableStep}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str | unicode]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  obsolete_file_path_list: {!r}\n'.format(indent, self.obsolete_file_path_list))
        str_list.extend(super(RunnableStep, self).trace(level=level + 1))

        return str_list

    def remove_obsolete_file_paths(self):
        """Remove file paths on the C{bsf.process.RunnableStep.obsolete_file_path_list} Python C{list}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @return:
        @rtype:
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
    @type file_path: str | unicode | None
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
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
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
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param file_path: File path
        @type file_path: str | unicode | None
        @param mode_directory: Directory access mode according to C{stat}
        @type mode_directory: str | None
        @param mode_file: File access mode for files according to C{stat}
        @type mode_file: str | None
        @return:
        @rtype:
        """
        super(RunnableStepChangeMode, self).__init__(
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
    @type source_path: str | unicode | None
    @ivar target_path: Target path
    @type target_path: str | unicode | None
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
        """Initialise a C{bsf.process.RunnableStepCopy}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param source_path: Source path
        @type source_path: str | unicode | None
        @param target_path: Target path
        @type target_path: str | unicode | None
        @return:
        @rtype:
        """

        super(RunnableStepCopy, self).__init__(
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
        """Initialise a C{bsf.process.RunnableStepJava}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode | None
        @return:
        @rtype:
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
        """Initialise a C{bsf.process.RunnableStepPicard}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode | None
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str | None
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode | None
        @param picard_classpath: Picard class path
        @type picard_classpath: str | unicode | None
        @param picard_command: Picard command
        @type picard_command: str | None
        @return:
        @rtype:
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
        @return:
        @rtype:
        """
        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class RunnableStepLink(RunnableStep):
    """The C{bsf.process.RunnableStepLink} represents a step creating a symbolic link.

    Attributes:
    @ivar source_path: Source path
    @type source_path: str | unicode | None
    @ivar target_path: Target path
    @type target_path: str | unicode | None
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
        """Initialise a C{bsf.process.RunnableStepLink}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param source_path: Source path
        @type source_path: str | unicode | None
        @param target_path: Target path
        @type target_path: str | unicode | None
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
    @type directory_path: str | unicode | None
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
        """Initialise a C{bsf.process.RunnableStepMakeDirectory}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param directory_path: Directory path
        @type directory_path: str | unicode | None
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
    @type file_path: str | unicode | None
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
            file_path=None):
        """Initialise a C{bsf.process.RunnableStepMakeNamedPipe}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param file_path: Named pipe file path
        @type file_path: str | unicode | None
        @return:
        @rtype:
        """
        super(RunnableStepMakeNamedPipe, self).__init__(
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
    @type source_path: str | unicode | None
    @ivar target_path: Target path
    @type target_path: str | unicode | None
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
        """Initialise a C{bsf.process.RunnableStepMove}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param source_path: Source path
        @type source_path: str | unicode | None
        @param target_path: Target path
        @type target_path: str | unicode | None
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
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            obsolete_file_path_list=None,
            sleep_time=None):
        """Initialise a C{bsf.process.RunnableStepSleep}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param sleep_time: Sleep time in seconds
        @type sleep_time: float | None
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
            stdout_path=None,
            stderr_path=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
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
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode | None
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param key: Environment key
        @type key: str
        @return:
        @rtype:
        """
        super(RunnableStepSetEnvironment, self).__init__(
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
