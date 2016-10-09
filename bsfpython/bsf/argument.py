"""bsf.argument

A package of classes and methods to model C{bsf.process.Command} arguments,
i.e. C{bsf.argument.SwitchLong} (--key), C{bsf.argument.SwitchShort} (-k),
C{bsf.argument.OptionLong} (--key value), C{bsf.argument.OptionShort} (-k value),
C{bsf.argument.OptionPair} (key=value), C{bsf.argument.SubCommand} (key) and
C{bsf.argument.Argument} (key).
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


class Argument(object):
    """The C{bsf.argument.Argument} class represents an argument obeying a C{key} schema.

    Attributes:
    @ivar key: Key
    @type key: str
    """

    @classmethod
    def from_key_value(cls, key, value):
        """Create a C{bsf.argument.Argument} or sub-class object from a key and value argument pair.

        If the key starts with one or two hyphen and has no value associated, create a
        C{bsf.argument.SwitchShort} or C{bsf.argument.SwitchLong} object, respectively.
        If the key is additionally associated with a value, create a C{bsf.argument.OptionShort} or
        a C{bsf.argument.OptionLong} object, respectively.
        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @return: C{bsf.argument.Argument} object or a sub-class thereof
        @rtype: bsf.argument.Argument
        """

        if key.startswith('--'):
            # Long option or switch ...
            if value:
                argument = OptionLong(key=key.lstrip('-'), value=value)
            else:
                argument = SwitchLong(key=key.lstrip('-'))
        elif key.startswith('-'):
            # Short option or switch ...
            if value:
                argument = OptionShort(key=key.lstrip('-'), value=value)
            else:
                argument = SwitchShort(key=key.lstrip('-'))
        else:
            # Must be an option pair ...
            argument = OptionPair(key=key, value=value)

        return argument

    def __init__(self, key):
        """Initialise a C{bsf.argument.Argument} object.

        @param key: Key
        @type key: str
        @return:
        @rtype:
        """

        assert isinstance(key, basestring)

        if not key:
            raise Exception('The key argument has to be defined.')

        super(Argument, self).__init__()

        self.key = key

        return

    def trace(self, level):
        """Trace a C{bsf.argument.Argument} object.

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
    """The C{bsf.argument.Switch} class represents an argument obeying a C{--key} or C{-k} schema.

    Attributes:
    @ivar is_long: GNU-style long switch (C{--key})
    @type is_long: bool
    """

    def __init__(self, key, is_long=False):
        """Initialise a C{bsf.argument.Switch} object.

        @param key: Key
        @type key: str
        @param is_long: GNU-style long switch (C{--key})
        @type is_long: bool
        @return:
        @rtype:
        """

        assert isinstance(is_long, bool)

        super(Switch, self).__init__(key=key)

        self.is_long = is_long

        return

    def trace(self, level):
        """Trace a C{bsf.argument.Switch} object.

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
    """The C{bsf.argument.SwitchLong} class represents an argument obeying a C{--key} schema.

    Attributes:
    """

    def __init__(self, key):
        """Initialise a C{bsf.argument.SwitchLong} object.

        @param key: Key
        @type key: str
        @return:
        @rtype:
        """

        super(SwitchLong, self).__init__(key=key, is_long=True)

        return

    def trace(self, level):
        """Trace a C{bsf.argument.SwitchLong} object.

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
    """The C{bsf.argument.SwitchShort} class represents an argument obeying a C{-k} schema.

    Attributes:
    """

    def __init__(self, key):
        """Initialise a C{bsf.argument.SwitchShort} object.

        @param key: Key
        @type key: str
        @return:
        @rtype:
        """

        super(SwitchShort, self).__init__(key=key, is_long=False)

        return

    def trace(self, level):
        """Trace a C{bsf.argument.SwitchShort} object.

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
    """The C{bsf.argument.Option} class represents arguments obeying a C{--key value} or C{-k value} schema.

    Attributes:
    @ivar is_pair: A I{KEY=VALUE} pair
    @type is_pair: bool
    @ivar value: Value
    @type value: str | unicode
    """

    def __init__(self, key, value, is_long=False, is_pair=False):
        """Initialise a C{bsf.argument.Option} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode | None
        @param is_long: GNU-style long option (C{--key value})
        @type is_long: bool
        @param is_pair: A I{KEY=VALUE} pair
        @type is_pair: bool
        @return:
        @rtype:
        """

        assert isinstance(value, basestring)
        assert isinstance(is_pair, bool)

        super(Option, self).__init__(key=key, is_long=is_long)

        if not value:
            raise Exception('The value argument has to be defined.')

        self.value = value
        self.is_pair = is_pair

        return

    def trace(self, level):
        """Trace a C{bsf.argument.Option} object.

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
    """The C{bsf.argument.OptionLong} class represents an argument obeying a C{--key value} schema.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise a C{bsf.argument.OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @return:
        @rtype:
        """

        super(OptionLong, self).__init__(key=key, value=value, is_long=True)

        return

    def trace(self, level):
        """Trace a C{bsf.argument.OptionLong} object.

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
    """The C{bsf.argument.OptionShort} class represents an argument obeying a C{-k value} schema.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise a C{bsf.argument.OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @return:
        @rtype:
        """

        super(OptionShort, self).__init__(key=key, value=value, is_long=False)

        return

    def trace(self, level):
        """Trace a C{bsf.argument.OptionShort} object.

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
    """The C{bsf.argument.OptionPair} class represents an argument obeying a I{KEY=VALUE} schema.

    Although the I{KEY=VALUE} expressions could be added as simple arguments, the benefit
    of having a specific sub-class is the tracking of the key in the C{bsf.process.Command.options}
    dictionary.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise a C{bsf.argument.OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @return:
        @rtype:
        """

        super(OptionPair, self).__init__(key=key, value=value, is_pair=True)

        return

    def trace(self, level):
        """Trace a C{bsf.argument.OptionPair} object.

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
