"""bsf.argument

A package of classes and methods to model C{bsf.process.Command} arguments, i.e.
C{bsf.argument.SwitchLong} (--key),
C{bsf.argument.SwitchShort} (-k),
C{bsf.argument.OptionLong} (--key value),
C{bsf.argument.OptionShort} (-k value),
C{bsf.argument.OptionPair} (key=value),
C{bsf.argument.OptionPairLong} (--key=value) and
C{bsf.argument.OptionPairShort} (-key=value).
"""


#
# Copyright 2013 - 2018 Michael K. Schuster
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


class Argument(object):
    """The C{bsf.argument.Argument} class represents an argument obeying a I{key} schema.

    Attributes:
    @ivar key: Key
    @type key: str | unicode
    """

    @classmethod
    def from_key_value(cls, key, value):
        """Create a C{bsf.argument.Argument} from a key and value argument pair.

        If the key starts with one or two hyphen and has no value associated, create a
        C{bsf.argument.SwitchShort} or C{bsf.argument.SwitchLong}, respectively.
        If the key is additionally associated with a value, create a C{bsf.argument.OptionShort} or
        a C{bsf.argument.OptionLong}, respectively.
        @param key: Key
        @type key: str | unicode
        @param value: Value
        @type value: str | unicode
        @return: C{bsf.argument.Argument}
        @rtype: bsf.argument.Argument
        """
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
        """Initialise a C{bsf.argument.Argument}.

        @param key: Key
        @type key: str | unicode
        @return:
        @rtype:
        """
        super(Argument, self).__init__()

        assert isinstance(key, (str, unicode))

        if not key:
            raise Exception('The key argument has to be defined.')

        self.key = key

        return

    def trace(self, level):
        """Trace a C{bsf.argument.Argument}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """
        indent = '  ' * level

        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  key: {!r}\n'.format(indent, self.key)

        return output

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}

        @return: String representation
        @rtype: str | unicode
        """
        return self.key

    def get_list(self):
        """Get the list representation as Python C{list}.

        @return: List representation
        @rtype: list[str | unicode]
        """
        return [self.get_str()]


class Switch(Argument):
    """The C{bsf.argument.Switch} class represents an argument obeying a I{--key} or I{-k} schema.

    Attributes:
    """
    pass


class SwitchLong(Switch):
    """The C{bsf.argument.SwitchLong} class represents an argument obeying a I{--key} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}

        Overrides method Switch.get_str() to prepend two dashes ('--key').
        @return: String representation
        @rtype: str | unicode
        """
        return '--' + super(SwitchLong, self).get_str()


class SwitchShort(Switch):
    """The C{bsf.argument.SwitchShort} class represents an argument obeying a I{-k} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}

        Overrides method Switch.get_str() to prepend one dash ('-key').
        @return: String representation
        @rtype: str | unicode
        """
        return '-' + super(SwitchShort, self).get_str()


class Option(Switch):
    """The C{bsf.argument.Option} class represents arguments obeying a I{--key value} or I{-k value} schema.

    Attributes:
    @ivar value: Value
    @type value: str | unicode
    """

    def __init__(self, key, value):
        """Initialise a C{bsf.argument.Option}.

        @param key: Key
        @type key: str | unicode
        @param value: Value
        @type value: str | unicode
        @return:
        @rtype:
        """
        super(Option, self).__init__(key=key)

        assert isinstance(value, (str, unicode))

        if not value:
            raise Exception('The value argument has to be defined.')

        self.value = value

        return

    def trace(self, level):
        """Trace a C{bsf.argument.Option}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """
        indent = '  ' * level

        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  key: {!r}\n'.format(indent, self.key)
        output += '{}  value:   {!r}\n'.format(indent, self.value)

        # Do not call (i.e. trace) the super class, as the separation of key and value becomes quite confusing.
        # output += super(Option, self).trace(level=level + 1)

        return output

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}

        Overrides method Switch.get_str() to join key and value with space ('key value1 value2 ...').
        @return: String representation
        @rtype: str | unicode
        """
        return ' '.join((self.key, self.value))

    def get_list(self):
        """Get the list representation as Python C{list}.

        Overrides method Switch.get_list() to split the Python str representation on
        white space only once (['key', 'value1 value2 ...']).
        @return: List representation
        @rtype: list[str | unicode]
        """
        return self.get_str().split(None, 1)


class OptionLong(Option):
    """The C{bsf.argument.OptionLong} class represents an argument obeying a I{--key value} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}

        Overrides method Option.get_str() to prepend two dashes ('--key value1 value2 ...').
        @return: String representation
        @rtype: str | unicode
        """
        return '--' + super(OptionLong, self).get_str()


class OptionShort(Option):
    """The C{bsf.argument.OptionShort} class represents an argument obeying a I{-k value} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}

        Overrides method Option.get_str() to prepend one dash ('-key value1 value2 ...').
        @return: String representation
        @rtype: str | unicode
        """
        return '-' + super(OptionShort, self).get_str()


class OptionPair(Option):
    """The C{bsf.argument.OptionPair} class represents an argument obeying a I{KEY=VALUE} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}.

        @return: String representation
        @rtype: str | unicode
        """
        return '='.join((self.key, self.value))

    def get_list(self):
        """Get the list representation as Python C{list}.

        Overrides method Option.get_list() to avoid splitting by white space (['key=value1 value2 ...']).
        @return: List representation
        @rtype: list[str | unicode]
        """
        return [self.get_str()]


class OptionPairLong(OptionPair):
    """The C{bsf.argument.OptionPairLong} class represents an argument obeying a I{--KEY=VALUE} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}.

        Overrides method OptionPair.get_str() to prepend two dashes ('--key=value1 value2 ...').
        @return: String representation
        @rtype: str | unicode
        """
        return '--' + super(OptionPairLong, self).get_str()


class OptionPairShort(OptionPair):
    """The C{bsf.argument.OptionPairShort} class represents an argument obeying a I{-KEY=VALUE} schema.

    Attributes:
    """

    def get_str(self):
        """Get the string representation as Python C{str} or C{unicode}.

        Overrides method OptionPair.get_str() to prepend one dash ('-key=value1 value2 ...').
        @return: String representation
        @rtype: str | unicode
        """
        return '-' + super(OptionPairShort, self).get_str()


class OptionMulti(Option):
    """The C{bsf.argument.Option} class represents arguments obeying a I{--key value1 value2}
    or I{-k value1 value2} schema.

    Attributes:
    """

    def get_list(self):
        """Get the list representation as Python C{list}.

        Override method Option.get_list() to split on white space (['key', 'value1', 'value2', '...']).
        This method supports programs like the STAR aligner.
        @return: List representation
        @rtype: list[str | unicode]
        """
        return self.get_str().split()


class OptionMultiLong(OptionMulti, OptionLong):
    """The C{bsf.argument.OptionMultiLong} class represents arguments obeying a I{--key value1 value2} or schema.

    The order of inheritance (OptionMulti, OptionLong) ascertains that OptionMulti.get_list()
    gets called before OptionLong.get_list().
    Attributes:
    """
    pass


class OptionMultiShort(OptionMulti, OptionShort):
    """The C{bsf.argument.OptionMultiLong} class represents arguments obeying a I{-k value1 value2} schema.

    The order of inheritance (OptionMulti, OptionShort) ascertains that OptionMulti.get_list()
    gets called before OptionShort.get_list().
    Attributes:
    """
    pass


class OptionMultiPair(OptionMulti, OptionPair):
    """The C{bsf.argument.OptionMultiLong} class represents arguments obeying a I{key=value1 value2} schema.

    The order of inheritance (OptionMulti, OptionPair) ascertains that OptionMulti.get_list()
    gets called before OptionPair.get_list().
    Attributes:
    """
    pass


class OptionMultiPairLong(OptionMulti, OptionPairLong):
    """The C{bsf.argument.OptionMultiPairLong} class represents arguments obeying a I{--key=value1 value2} schema.

    The order of inheritance (OptionMulti, OptionPairLong) ascertains that OptionMulti.get_list()
    gets called before OptionPairLong.get_list().
    Attributes:
    """
    pass


class OptionMultiPairShort(OptionMulti, OptionPairShort):
    """The C{bsf.argument.OptionMultiPairShort} class represents arguments obeying a I{-key=value1 value2} schema.

    The order of inheritance (OptionMulti, OptionPairShort) ascertains that OptionMulti.get_list()
    gets called before OptionPairShort.get_list().
    Attributes:
    """
    pass
