"""Bio.BSF.Argument

A package of classes and methods to model Command arguments,
i.e. SwitchLong (--key), SwitchShort (-k),
OptionLong (--key value), OptionShort (-k value), OptionPair (key=value),
SubCommand (key) and Argument (key).
"""

#
# Copyright 2013 Michael K. Schuster
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
    """The C{Argument} class represents an argument obeying a C{key} schema.

    Attributes:
    @ivar key: Key
    @type key: str
    """

    @classmethod
    def from_key_value(cls, key, value):
        """Create an C{Argument} or sub-class object from a key and value argument pair.

        If the key starts with one or two hyphen and has no value associated, create a
        C{SwitchShort} or C{SwitchLong} object, respectively. If the key is additionally associated
        with a value, create an C{OptionShort} or C{OptionLong} object, respectively.
        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @return: C{Argument} object or a sub-class thereof
        @rtype: Argument
        """

        if key.startswith('--'):

            # Long option or switch ...

            if value:
                argument = OptionLong(key=key.lstrip('-'), value=str(value))
            else:
                argument = SwitchLong(key=key.lstrip('-'))

        elif key.startswith('-'):

            # Short option or switch ...

            if value:
                argument = OptionShort(key=key.lstrip('-'), value=str(value))
            else:
                argument = SwitchShort(key=key.lstrip('-'))

        else:

            # Must be an option pair ...

            argument = OptionPair(key=str(key), value=str(value))

        return argument

    def __init__(self, key):
        """Initialise an C{Argument} object.

        @param key: Key
        @type key: str
        """

        assert isinstance(key, basestring)

        if not key:
            raise Exception('The key argument has to be defined.')

        self.key = key

    def trace(self, level):
        """Trace an C{Argument} object.

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


class Switch(Argument):
    """The C{Switch} class represents an argument obeying a C{--key} or C{-k} schema.

    Attributes:
    @ivar is_long: GNU-style long switch (C{--key})
    @type is_long: bool
    """

    def __init__(self, key, is_long=False):
        """Initialise a C{Switch} object.

        @param key: Key
        @type key: str
        @param is_long: GNU-style long switch (C{--key})
        @type is_long: bool
        """

        super(Switch, self).__init__(key=key)

        assert isinstance(is_long, bool)

        self.is_long = is_long

    def trace(self, level):
        """Trace a C{Switch} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  is_long: {!r}\n'.format(indent, self.is_long)

        output += super(Switch, self).trace(level=level + 1)

        return output


class SwitchLong(Switch):
    """The C{SwitchLong} class represents an argument obeying a C{--key} schema.

    Attributes:
    """

    def __init__(self, key):
        """Initialise a C{SwitchLong} object.

        @param key: Key
        @type key: str
        """

        super(SwitchLong, self).__init__(key=key, is_long=True)

    def trace(self, level):
        """Trace a C{SwitchLong} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(SwitchLong, self).trace(level=level + 1)

        return output


class SwitchShort(Switch):
    """The C{SwitchShort} class represents an argument obeying a C{-k} schema.

    Attributes:
    """

    def __init__(self, key):
        """Initialise a C{SwitchShort} object.

        @param key: Key
        @type key: str
        """

        super(SwitchShort, self).__init__(key=key, is_long=False)

    def trace(self, level):
        """Trace a C{SwitchShort} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(SwitchShort, self).trace(level=level + 1)

        return output


class Option(Switch):
    """The C{Option} class represents arguments obeying a C{--key value} or C{-k value} schema.

    Attributes:
    @ivar is_pair: A C{KEY=VALUE} pair
    @type is_pair: bool
    @ivar value: Value
    @type value: str | unicode
    """

    def __init__(self, key, value, is_long=False, is_pair=False):
        """Initialise an C{Option} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode | None
        @param is_long: GNU-style long option (C{--key value})
        @type is_long: bool
        @param is_pair: A C{KEY=VALUE} pair
        @type is_pair: bool
        """

        assert isinstance(value, basestring)
        assert isinstance(is_pair, bool)

        super(Option, self).__init__(key=key, is_long=is_long)

        if not value:
            raise Exception('The value argument has to be defined.')

        self.value = value
        self.is_pair = is_pair

    def trace(self, level):
        """Trace an C{Option} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  is_pair: {!r}\n'.format(indent, self.is_pair)
        output += '{}  value:   {!r}\n'.format(indent, self.value)

        output += super(Option, self).trace(level=level + 1)

        return output


class OptionLong(Option):
    """The C{OptionLong} class represents an argument obeying a C{--key value} schema.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise an C{OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        """

        super(OptionLong, self).__init__(key=key, value=value, is_long=True)

    def trace(self, level):
        """Trace an C{OptionLong} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(OptionLong, self).trace(level=level + 1)

        return output


class OptionShort(Option):
    """The C{OptionShort} class represents an argument obeying a C{-k value} schema.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise an C{OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        """

        super(OptionShort, self).__init__(key=key, value=value, is_long=False)

    def trace(self, level):
        """Trace an C{OptionShort} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(OptionShort, self).trace(level=level + 1)

        return output


class OptionPair(Option):
    """The C{OptionPair} class represents an argument obeying a C{KEY=VALUE} schema.

    Although the C{KEY=VALUE} expressions could be added as simple arguments, the benefit
    of having a specific sub-class is the tracking of the key in the C{Command.options}
    dictionary.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise an C{OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        """

        super(OptionPair, self).__init__(key=key, value=value, is_pair=True)

    def trace(self, level):
        """Trace an C{OptionPair} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(OptionPair, self).trace(level=level + 1)

        return output
