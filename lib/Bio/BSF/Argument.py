"""Bio.BSF.Argument

A package of classes and methods to model Bio.BSF.Command arguments,
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
    """Bio.BSF.Argument.Argument class for an argument obeying a "key" schema.

    Attributes:
    :ivar key: Key
    :type key: str
    """

    @classmethod
    def from_key_value(cls, key, value):

        """Create a Bio.BSF.Argument.Argument or sub-class object from a key and value argument pair.

        If the key starts with one or two hyphen and has no value associated, create a
        SwitchShort or SwitchLong object, respectively. If the key is additionally associated
        with a value, create an OptionShort or OptionLong object, respectively.
        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Bio.BSF.Argument object or a sub-class thereof
        :rtype: Option, Switch
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

        """Initialise a Bio.BSF.Argument.Argument object.

        :param key: Key
        :type key: str
        :return: Nothing
        :rtype: None
        """

        assert isinstance(key, basestring)

        if not key:
            message = 'The key argument has to be defined.'
            raise Exception(message)

        self.key = key

    def trace(self, level):

        """Trace a Bio.BSF.Argument.Argument object.

        :param self: Bio.BSF.Argument.Argument
        :type self: Argument
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  key: {!r}\n'.format(indent, self.key)

        return output


class Switch(Argument):
    """Bio.BSF.Argument.Switch class, for an argument obeying a "--key" or "-k" schema.

    Attributes:
    :ivar is_long: GNU-style long switch (--key)
    :type is_long: bool
    """

    def __init__(self, key, is_long=False):
        """Initialise a Bio.BSF.Argument.Switch object.

        :param key: Key
        :type key: str
        :param is_long: GNU-style long switch (--key)
        :type is_long: bool
        :return: Nothing
        :rtype: None
        """

        super(Switch, self).__init__(key=key)

        assert isinstance(is_long, bool)

        self.is_long = is_long

    def trace(self, level):
        """Trace a Bio.BSF.Argument.Switch object.

        :param self: Bio.BSF.Argument.Switch
        :type self: Switch
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  is_long: {!r}\n'.format(indent, self.is_long)

        output += super(Switch, self).trace(level=level + 1)

        return output


class SwitchLong(Switch):
    """Bio.BSF.Argument.SwitchLong class, for an argument obeying a "--key" schema.

    Attributes:
    """

    def __init__(self, key):
        """Initialise a Bio.BSF.Argument.SwitchLong object.

        :param key: Key
        :type key: str
        :return: Nothing
        :rtype: None
        """

        super(SwitchLong, self).__init__(key=key, is_long=True)

    def trace(self, level):
        """Trace a Bio.BSF.Argument.SwitchLong object.

        :param self: Bio.BSF.Argument.SwitchLong
        :type self: SwitchLong
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(SwitchLong, self).trace(level=level + 1)

        return output


class SwitchShort(Switch):
    """Bio.BSF.Argument.SwitchShort class, for an argument obeying a "-k" schema.

    Attributes:
    """

    def __init__(self, key):
        """Initialise a Bio.BSF.Argument.SwitchShort object.

        :param key: Key
        :type key: str
        :return: Nothing
        :rtype: None
        """

        super(SwitchShort, self).__init__(key=key, is_long=False)

    def trace(self, level):
        """Trace a Bio.BSF.Argument.SwitchShort object.

        :param self: Bio.BSF.Argument.SwitchShort
        :type self: SwitchShort
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(SwitchShort, self).trace(level=level + 1)

        return output


class Option(Switch):
    """Bio.BSF.Argument.Option class for arguments obeying a "--key value" or "-k value" schema.

    Attributes:
    :ivar is_pair: A KEY=VALUE pair
    :type is_pair: bool
    :ivar value: Value
    :type value: str
    """

    def __init__(self, key, value, is_long=False, is_pair=False):
        """Initialise a Bio.BSF.Argument.Option object.

        :param key: Key
        :type key: str
        :param value: Value
        :type value: str, None
        :param is_long: GNU-style long option (--key value)
        :type is_long: bool
        :param is_pair: A KEY=VALUE pair
        :type is_pair: bool
        :return: Nothing
        :rtype: None
        """

        assert isinstance(value, basestring)
        assert isinstance(is_pair, bool)

        super(Option, self).__init__(key=key, is_long=is_long)

        if not value:
            message = 'The value argument has to be defined.'
            raise Exception(message)

        self.value = value
        self.is_pair = is_pair

    def trace(self, level):
        """Trace a Bio.BSF.Argument.Option object.

        :param self: Bio.BSF.Argument.Option
        :type self: Option
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  is_pair: {!r}\n'.format(indent, self.is_pair)
        output += '{}  value:   {!r}\n'.format(indent, self.value)

        output += super(Option, self).trace(level=level + 1)

        return output


class OptionLong(Option):
    """Bio.BSF.Argument.OptionLong class, representing an argument obeying a "--key value" schema.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise a Bio.BSF.Argument.OptionLong object.

        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Nothing
        :rtype: None
        """

        super(OptionLong, self).__init__(key=key, value=value, is_long=True)

    def trace(self, level):
        """Trace a Bio.BSF.Argument.OptionLong object.

        :param self: Bio.BSF.Argument.OptionLong
        :type self: OptionLong
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(OptionLong, self).trace(level=level + 1)

        return output


class OptionShort(Option):
    """Bio.BSF.Argument.OptionShort class, representing an argument obeying a "-k value" schema.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise a Bio.BSF.Argument.OptionShort object.

        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Nothing
        :rtype: None
        """

        super(OptionShort, self).__init__(key=key, value=value, is_long=False)

    def trace(self, level):
        """Trace a Bio.BSF.Argument.OptionShort object.

        :param self: Bio.BSF.Argument.OptionShort
        :type self: OptionShort
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(OptionShort, self).trace(level=level + 1)

        return output


class OptionPair(Option):
    """Bio.BSF.Argument.OptionPair class, for an argument obeying a KEY=VALUE schema.

    Although the KEY=VALUE expressions could be added as simple arguments, the benefit
    of having a specific sub-class is the tracking of the key in the Bio.BSF:Command.options
    dictionary.

    Attributes:
    """

    def __init__(self, key, value):
        """Initialise a Bio.BSF.Argument.OptionPair object.

        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Nothing
        :rtype: None
        """

        super(OptionPair, self).__init__(key=key, value=value, is_pair=True)

    def trace(self, level):
        """Trace a Bio.BSF.Argument.OptionPair object.

        :param self: Bio.BSF.Argument.OptionPair
        :type self: OptionPair
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)

        output += super(OptionPair, self).trace(level=level + 1)

        return output
