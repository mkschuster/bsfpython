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
"""The :py:mod:`bsf.argument` module provides classes to model :py:class:`bsf.process.Command` arguments.

    - :py:class:`bsf.argument.SwitchLong` (:literal:`--key`),
    - :py:class:`bsf.argument.SwitchShort` (:literal:`-k`),
    - :py:class:`bsf.argument.OptionLong` (:literal:`--key value`),
    - :py:class:`bsf.argument.OptionShort` (:literal:`-k value`),
    - :py:class:`bsf.argument.OptionPair` (:literal:`key=value`),
    - :py:class:`bsf.argument.OptionPairLong` (:literal:`--key=value`)
    - :py:class:`bsf.argument.OptionPairShort` (:literal:`-key=value`)
"""
import sys

from typing import List, Type

__all__ = \
    'Argument', \
    'Switch', 'SwitchLong', 'SwitchShort', \
    'Option', 'OptionLong', 'OptionShort', \
    'OptionPair', 'OptionPairLong', 'OptionPairShort', \
    'OptionMulti', 'OptionMultiLong', 'OptionMultiShort', \
    'OptionMultiPair', 'OptionMultiPairLong', 'OptionMultiPairShort'


class Argument(object):
    """The :py:class:`bsf.argument.Argument` class represents an argument obeying a
    :literal:`key` schema.

    :ivar key: A key.
    :type key: str
    :ivar value: A value.
    :type value: str | None
    """

    @classmethod
    def from_key_value(cls, key, value):
        """Create a :py:class:`bsf.argument.Argument` object from a key and value argument pair.

        If the key starts with a valid :literal:`bsf.argument` class name, the corresponding object is initialised from
        the key, which is separated by a dot from the class name and, if applicable for the class, the value.

        :literal:`class.key = value`

        This is particularly useful for STAR option :literal:`--outSAMattributes NH HI NM MD XS AS` that contains
        a list of SAM tags that need splitting via class :py:class:`bsf.argument.OptionMultiPairLong`
        upon command line passing as list.

        :literal:`OptionMultiPairLong.outSAMattributes = NH HI NM MD XS AS`

        Legacy functionality:
        If the key starts with one or two hyphen and has no value associated, create a
        :py:class:`bsf.argument.SwitchShort` or :py:class:`bsf.argument.SwitchLong`, respectively.
        If the key is additionally associated with a value, create a :py:class:`bsf.argument.OptionShort` or
        a :py:class:`bsf.argument.OptionLong`, respectively.

        :param key: A key.
        :type key: str
        :param value: A value.
        :type value: str
        :return: A :py:class:`bsf.argument.Argument` object.
        :rtype: Argument
        """
        # Test for a >>>class.key = value<<< INI file option construct.
        if '.' in key:
            key_list = key.split(sep='.')
            if key_list[0] in __all__:
                argument_class: Type = getattr(sys.modules[__name__], key_list[0])

                if argument_class.__name__ in ('Argument', 'Switch', 'SwitchLong', 'SwitchShort'):
                    return argument_class(key=key_list[1])
                else:
                    return argument_class(key=key_list[1], value=value)

        # Legacy functionality below.

        if key.startswith('--'):
            # Long switch, option or option pair ...
            str_position = key.find('=')
            if str_position != -1:
                argument = OptionPairLong(key=key[:str_position].lstrip('-'), value=key[str_position + 1:])
            elif value:
                argument = OptionLong(key=key.lstrip('-'), value=value)
            else:
                argument = SwitchLong(key=key.lstrip('-'))
        elif key.startswith('-'):
            # Short switch, option or option pair ...
            str_position = key.find('=')
            if str_position != -1:
                argument = OptionPairShort(key=key[:str_position].lstrip('-'), value=key[str_position + 1:])
            elif value:
                argument = OptionShort(key=key.lstrip('-'), value=value)
            else:
                argument = SwitchShort(key=key.lstrip('-'))
        else:
            # Must be an option pair ...
            argument = OptionPair(key=key, value=value)

        return argument

    def __init__(self, key):
        """Initialise a :py:class:`bsf.argument.Argument` object.

        :param key: A key.
        :type key: str
        """
        super(Argument, self).__init__()

        if not key:
            raise Exception('The key argument has to be defined.')

        self.key = key
        self.value = None

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.argument.Argument` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  key: {!r}\n'.format(indent, self.key))

        return str_list

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return self.key

    def get_list(self):
        """Get the list representation as Python :py:class:`list`.

        :return: A Python :py:class:`list` representation.
        :rtype: list[str]
        """
        return [self.get_str()]


class Switch(Argument):
    """The :py:class:`bsf.argument.Switch` class represents an argument obeying a
    :literal:`--key` or :literal:`-k` schema.
    """
    pass


class SwitchLong(Switch):
    """The :py:class:`bsf.argument.SwitchLong` class represents an argument obeying a
    :literal:`--key` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`bsf.argument.Switch.get_str` to prepend two dashes
        (i.e., :literal:`'--key'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '--' + super(SwitchLong, self).get_str()


class SwitchShort(Switch):
    """The :py:class:`bsf.argument.SwitchShort` class represents an argument obeying a
    :literal:`-k` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`bsf.argument.Switch.get_str` to prepend one dash
        (i.e., :literal:`'-key'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '-' + super(SwitchShort, self).get_str()


class Option(Switch):
    """The :py:class:`bsf.argument.Option` class represents arguments obeying a
    :literal:`--key value` or :literal:`-k value` schema.

    :ivar value: A value.
    :type value: str
    """

    def __init__(self, key, value):
        """Initialise a :py:class:`bsf.argument.Option` object.

        :param key: A key.
        :type key: str
        :param value: A value.
        :type value: str
        """
        super(Option, self).__init__(key=key)

        if not value:
            raise Exception('The value argument has to be defined.')

        self.value = value

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.argument.Option` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  key: {!r}\n'.format(indent, self.key))
        str_list.append('{}  value:   {!r}\n'.format(indent, self.value))

        # Do not call (i.e., trace) the super class, as the separation of key and value becomes quite confusing.
        # str_list.extend(super(Option, self).trace(level=level + 1))

        return str_list

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`bsf.argument.Switch.get_str` to join key and value with space
        (i.e., :literal:`'key value'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return ' '.join((self.key, self.value))

    def get_list(self):
        """Get the list representation as Python :py:class:`list`.

        Overrides method :py:meth:`bsf.argument.Switch.get_list` to split the Python :py:class:`str` representation on
        white space only once (i.e., :literal:`['key', 'value']`).

        :return: A Python :py:class:`list` representation.
        :rtype: list[str]
        """
        return self.get_str().split(None, 1)


class OptionLong(Option):
    """The :py:class:`bsf.argument.OptionLong` class represents an argument obeying a
    :literal:`--key value` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`bsf.argument.Option.get_str` to prepend two dashes
        (i.e., :literal:`'--key value'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '--' + super(OptionLong, self).get_str()


class OptionShort(Option):
    """The :py:class:`bsf.argument.OptionShort` class represents an argument obeying a
    :literal:`-k value` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`baf.argument.Option.get_str` to prepend one dash
        (i.e., :literal:`'-key value'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '-' + super(OptionShort, self).get_str()


class OptionPair(Option):
    """The :py:class:`bsf.argument.OptionPair` class represents an argument obeying a
    :literal:`key=value` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`baf.argument.Option.get_str` to combine key and value with an equal sign.
        (i.e., :literal:`'key=value'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '='.join((self.key, self.value))

    def get_list(self):
        """Get the list representation as Python :py:class:`list`.

        Overrides method :py:meth:`bsf.argument.Option.get_list` to avoid splitting by white space
        (i.e, :literal:`['key=value']`).

        :return: A Python :py:class:`list` representation.
        :rtype: list[str]
        """
        return [self.get_str()]


class OptionPairLong(OptionPair):
    """The :py:class:`bsf.argument.OptionPairLong` class represents an argument obeying a
    :literal:`--key=value` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`bsf.argument.OptionPair.get_str` to prepend two dashes
        (i.e., :literal:`'--key=value'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '--' + super(OptionPairLong, self).get_str()


class OptionPairShort(OptionPair):
    """The :py:class:`bsf.argument.OptionPairShort` class represents an argument obeying a
    :literal:`-key=value` schema.
    """

    def get_str(self):
        """Get the string representation as Python :py:class:`str`.

        Overrides method :py:meth:`bsf.argument.OptionPair.get_str` to prepend one dash
        (i.e., :literal:`'-key=value'`).

        :return: A Python :py:class:`str` representation.
        :rtype: str
        """
        return '-' + super(OptionPairShort, self).get_str()


class OptionMulti(Option):
    """The :py:class:`bsf.argument.Option` class represents arguments obeying a
    :literal:`--key value1 value2` or :literal:`-k value1 value2` schema.
    """

    def get_list(self):
        """Get the list representation as Python :py:class:`list` object.

        Overrides method :py:meth:`bsf.argument.Option.get_list` to split on white space
        (i.e., :literal:`['key', 'value1', 'value2', '...']`).

        This method supports programs like the STAR aligner.

        :return: A Python :py:class:`list` representation.
        :rtype: list[str]
        """
        return self.get_str().split()


class OptionMultiLong(OptionMulti, OptionLong):
    """The :py:class:`bsf.argument.OptionMultiLong` class represents arguments obeying a
    :literal:`--key value1 value2` schema.

    The order of inheritance (:py:class:`bsf.argument.OptionMulti`, :py:class:`bsf.argument.OptionLong`)
    ascertains that
    :py:meth:`bsf.argument.OptionMulti.get_list` gets called before
    :py:meth:`bsf.argument.OptionLong.get_list`.
    """

    pass


class OptionMultiShort(OptionMulti, OptionShort):
    """The :py:class:`bsf.argument.OptionMultiLong` class represents arguments obeying a
    :literal:`-k value1 value2` schema.

    The order of inheritance (:py:class:`bsf.argument.OptionMulti`, :py:class:`bsf.argument.OptionShort`)
    ascertains that
    :py:meth:`bsf.argument.OptionMulti.get_list` gets called before
    :py:meth:`bsf.argument.OptionShort.get_list`.
    """

    pass


class OptionMultiPair(OptionMulti, OptionPair):
    """The :py:class:`bsf.argument.OptionMultiLong` class represents arguments obeying a
    :literal:`key=value1 value2` schema.

    The order of inheritance (:py:class:`bsf.argument.OptionMulti`, :py:class:`bsf.argument.OptionPair`)
    ascertains that
    :py:meth:`bsf.argument.OptionMulti.get_list` gets called before
    :py:meth:`bsf.argument.OptionPair.get_list`.
    """

    pass


class OptionMultiPairLong(OptionMulti, OptionPairLong):
    """The :py:class:`bsf.argument.OptionMultiPairLong` class represents arguments obeying a
    :literal:`--key=value1 value2` schema.

    The order of inheritance (:py:class:`bsf.argument.OptionMulti`, :py:class:`bsf.argument.OptionPairLong`)
    ascertains that
    :py:meth:`bsf.argument.OptionMulti.get_list` gets called before
    :py:meth:`bsf.argument.OptionPairLong.get_list`.
    """

    pass


class OptionMultiPairShort(OptionMulti, OptionPairShort):
    """The :py:class:`bsf.argument.OptionMultiPairShort` class represents arguments obeying a
    :literal:`-key=value1 value2` schema.

    The order of inheritance (:py:class:`bsf.argument.OptionMulti`, :py:class:`bsf.argument.OptionPairShort`)
    ascertains that
    :py:meth:`bsf.argument.OptionMulti.get_list` gets called before
    :py:meth:`bsf.argument.OptionPairShort.get_list`.
    """

    pass
