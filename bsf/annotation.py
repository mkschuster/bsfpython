# -*- coding: utf-8 -*-
"""Annotation module.

A package of classes and methods modelling Comma-Separated Value (CSV) and Tab-Separated Value (TSV) annotation files.
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
import re
import warnings
from csv import DictReader, DictWriter
from io import TextIOWrapper
from typing import Callable, Dict, List, Optional


class AnnotationSheet(object):
    """The C{bsf.annotation.AnnotationSheet} class represents comma-separated value (CSV) files.

    This class is a bit unusual in that values of instance variables around the
    C{csv.DictReader} and C{csv.DictWriter} classes can be initialised from a corresponding set of class variables.
    Therefore, subclasses with fixed defaults can be defined, while generic objects can be initialised directly
    via the C{bsf.annotation.AnnotationSheet.__init__} method without the need to create a
    C{bsf.annotation.AnnotationSheet} subclass in code.

    @cvar _field_names: Python C{list} if field (column) names
    @type _field_names: list[str]
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{Callable} value data
    @type _test_methods: Dict[str, List[Callable[[int, Dict[str, str], str], str]]]
    @ivar file_path: File path
    @type file_path: str | None
    @ivar file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type file_type: str
    @ivar name: Name
    @type name: str | None
    @ivar field_names: Python C{list} of Python C{str} (field name) objects
    @type field_names: list[str]
    @ivar test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{Callable} value data
    @type test_methods: Dict[str, List[Callable[[int, Dict[str, str], str], str]]]
    @ivar row_dicts: Python C{list} of Python C{dict} objects
    @type row_dicts: list[dict[str, str]]
    @ivar _csv_reader_file: C{io.TextIO}
    @type _csv_reader_file: io.TextIO | None
    @ivar _csv_reader_object: C{csv.DictReader}
    @type _csv_reader_object: DictReader | None
    @ivar _csv_writer_file: C{io.TextIO}
    @type _csv_writer_file: io.TextIO | None
    @ivar _csv_writer_object: C{csv.DictWriter}
    @type _csv_writer_object: DictWriter | None
    """

    # Regular expression for non-alphanumeric characters
    _regular_expression_non_alpha = re.compile(pattern='\\W')

    # Regular expression for non-numeric characters
    _regular_expression_non_numeric = re.compile(pattern='\\D')

    # Regular expression for non-sequence characters
    _regular_expression_non_sequence = re.compile(pattern='[^ACGTacgt]')

    # Regular expression for non-ambiguous sequence characters
    _regular_expression_non_ambiguous_sequence = re.compile(pattern='[^ACGTacgtWSMKRYwsmkryBDHVbdhvNn]')

    # Regular expression for multiple underscore characters
    _regular_expression_multiple_underscore = re.compile(pattern='_{2,}')

    # File type (i.e. 'excel' or 'excel-tab' defined in the csv.Dialect class)
    _file_type = 'excel'

    # Header line exists
    _header_line = True

    # Python list of Python str (field name) objects
    _field_names = list()

    # Python dict of Python str (field name) key data and
    # Python list of Python typing.Callable value data
    _test_methods: Dict[str, List[Callable[[int, Dict[str, str], str], str]]] = dict()

    # Python dict of (boolean state) Python str objects and Python bool value objects.
    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False
    }

    @classmethod
    def check_column_value(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Check for a column name and return its associated value, if any.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value, which also implies requiring the column_name
        @type require_value: bool
        @return: Python C{tuple} of Python C{str} (warning message) and Python C{str} (column value)
        @rtype: (str, str)
        """
        message = str()

        if column_name not in row_dict:
            if require_column or require_value:
                message += 'Column name {!r} not in dict for row {}.\n'.format(column_name, row_number)
            return message, None

        if row_dict[column_name]:
            return message, row_dict[column_name]
        else:
            if require_value:
                message += 'Column name {!r} without value in row {}.\n'.format(column_name, row_number)
            return message, None

    @classmethod
    def _check_alphanumeric(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as I{alphanumeric}.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only alphanumeric characters.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value, which also implies requiring the column_name.
        @type require_value: bool
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=require_column,
            require_value=require_value)

        if column_value is not None:
            match = re.search(pattern=cls._regular_expression_non_alpha, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with non-alphanumeric characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_alphanumeric(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{alphanumeric}.

        Neither the column nor the value needs existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_alphanumeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name)

    @classmethod
    def check_alphanumeric_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{alphanumeric}.

        The column, but not the value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_alphanumeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_alphanumeric_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{alphanumeric}.

        Both, the column and value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_alphanumeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def _check_numeric(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as I{numeric}.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only numeric characters.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value, which also implies requiring the column_name.
        @type require_value: bool
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=require_column,
            require_value=require_value)

        if column_value is not None:
            match = re.search(pattern=cls._regular_expression_non_numeric, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with non-numeric characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_numeric(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{numeric}.

        Neither the column nor the value needs existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_numeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=False,
            require_value=False)

    @classmethod
    def check_numeric_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{numeric}.

        The column, but not the value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_numeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_numeric_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{numeric}.

        Both, the column and value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_numeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def _check_sequence(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as I{IUPAC sequence}.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only valid IUPAC sequence characters.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value
        @type require_value: bool
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=require_column,
            require_value=require_value)

        if column_value is not None:
            match = re.search(pattern=cls._regular_expression_non_sequence, string=row_dict[column_name])
            if match:
                messages += 'Field {!r} in row {} contains a IUPAC sequence {!r} ' \
                            'with illegal characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_sequence(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{IUPAC sequence}.

        Neither the column nor the value needs existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=False,
            require_value=False)

    @classmethod
    def check_sequence_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{IUPAC sequence}.

        The column, but not the value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_sequence_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{IUPAC sequence}.

        Both, the column and value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def _check_ambiguous_sequence(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as I{IUPAC ambiguous sequence}.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only valid IUPAC ambiguity sequence characters.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value
        @type require_value: bool
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=require_column,
            require_value=require_value)

        if column_value is not None:
            match = re.search(pattern=cls._regular_expression_non_ambiguous_sequence, string=row_dict[column_name])
            if match:
                messages += 'Field {!r} in row {} contains a IUPAC ambiguous sequence {!r} ' \
                            'with illegal characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_ambiguous_sequence(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{IUPAC ambiguous sequence}.

        Neither the column nor the value needs existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_ambiguous_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=False,
            require_value=False)

    @classmethod
    def check_ambiguous_sequence_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{IUPAC ambiguous sequence}.

        The column, but not the value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_ambiguous_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_ambiguous_sequence_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{IUPAC ambiguous sequence}.

        Both, the column and value need existing.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls._check_ambiguous_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def check_underscore_leading(cls, row_number, row_dict, column_name):
        """Validate a particular column value for I{leading underscore characters}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has no leading underscore.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name)

        if column_value is not None:
            if column_value[:1] == '_':
                messages += 'Column {!r} in row {} contains a value {!r} with a leading underscore character.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_underscore_trailing(cls, row_number, row_dict, column_name):
        """Validate a particular column value for I{trailing underscore characters}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has no trailing underscore.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name)

        if column_value is not None:
            if column_value[-1:] == '_':
                messages += 'Column {!r} in row {} contains a value {!r} with a trailing underscore character.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_underscore_multiple(cls, row_number, row_dict, column_name):
        """Validate a particular column value for I{multiple underscore characters}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has not multiple underscore adjacent to each other.

        @param row_number: The row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column_value(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name)

        if column_value is not None:
            match = re.search(pattern=cls._regular_expression_multiple_underscore, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with multiple underscore characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def from_file_path(cls, file_path=None, file_type=None, name=None):
        """Construct a C{bsf.annotation.AnnotationSheet} from a comma-separated value (CSV) file.

        This method reads the whole CSV file at once and stores a Python C{list} of Python C{dict} objects
        representing each row. For large files, the C{bsf.annotation.AnnotationSheet.csv_reader_open},
        C{bsf.annotation.AnnotationSheet.csv_reader_next} and
        C{bsf.annotation.AnnotationSheet.csv_reader_close} methods should be called explicitly.
        @param file_path: File path
        @type file_path: str
        @param file_type: File type (i.e. I{excel} or I{excel_tab} defined in the C{csv.Dialect} class)
        @type file_type: str
        @param name: Name
        @type name: str
        @return: C{bsf.annotation.AnnotationSheet}
        @rtype: AnnotationSheet
        """
        annotation_sheet = cls(file_path=file_path, file_type=file_type, name=name)

        annotation_sheet.csv_reader_open()

        for row_dict in annotation_sheet._csv_reader_object:
            annotation_sheet.row_dicts.append(row_dict)

        annotation_sheet.csv_reader_close()

        return annotation_sheet

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            header=None,
            field_names=None,
            test_methods: Optional[Dict[str, List[Callable[[int, Dict[str, str], str], str]]]] = None,
            row_dicts=None):
        """Initialise a C{bsf.annotation.AnnotationSheet} object.

        @param file_path: File path
        @type file_path: str | None
        @param file_type: File type (i.e. I{excel} or I{excel_tab} defined in the C{csv.Dialect} class)
        @type file_type: str | None
        @param name: Name | None
        @type name: str | None
        @param header: Header line
        @type header: bool | None
        @param field_names: Python C{list} of Python C{str} (field name) objects
        @type field_names: list[str] | None
        @param test_methods: Python C{dict} of Python C{str} (field name) key data and
            Python C{list} of Python C{Callable} value data
        @type test_methods: Optional[Dict[str, List[Callable[[int, Dict[str, str], str], str]]]]
        @param row_dicts: Python C{list} of Python C{dict} objects
        @type row_dicts: list[dict[str, str]] | None
        """
        super(AnnotationSheet, self).__init__()

        self.file_path = file_path

        if file_type is None:
            # Copy the class variable.
            self.file_type = str(self._file_type)
        else:
            self.file_type = file_type

        self.name = name

        if header is None:
            # Copy the class variable.
            self.header = bool(self._header_line)
        else:
            self.header = header

        if field_names is None:
            # Copy the class variable.
            self.field_names = list(self._field_names)
        else:
            self.field_names = field_names

        if test_methods is None:
            # Copy the class variable.
            self.test_methods = self._test_methods.copy()
        else:
            self.test_methods = test_methods

        if row_dicts is None:
            self.row_dicts = list()
        else:
            self.row_dicts = row_dicts

        self._csv_reader_file: Optional[TextIOWrapper] = None
        self._csv_reader_object: Optional[DictReader] = None
        self._csv_writer_file: Optional[TextIOWrapper] = None
        self._csv_writer_object: Optional[DictWriter] = None

        return

    def trace(self, level=1):
        """Trace a C{bsf.annotation.SampleAnnotationSheet} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  file_path:    {!r}\n'.format(indent, self.file_path))
        str_list.append('{}  file_type:    {!r}\n'.format(indent, self.file_type))
        str_list.append('{}  name:         {!r}\n'.format(indent, self.name))
        str_list.append('{}  header:       {!r}\n'.format(indent, self.header))
        str_list.append('{}  field_names:  {!r}\n'.format(indent, self.field_names))
        str_list.append('{}  test_methods: {!r}\n'.format(indent, self.test_methods))
        str_list.append('{}  row_dicts:    {}\n'.format(indent, self.row_dicts))

        str_list.append('{}  Row dictionaries {}:\n'.format(indent, self.row_dicts))
        for row_dict in self.row_dicts:
            str_list.append('{}    {!r}\n'.format(indent, row_dict))

        return str_list

    def csv_reader_open(self):
        """Open a Comma-Separated Value (CSV) file linked to a C{bsf.annotation.AnnotationSheet} object for reading
        and initialise a Python C{csv.DictReader} object.
        """
        if not self.file_path:
            raise Exception('Cannot read an AnnotationSheet without a valid file_name.')

        # Although the AnnotationSheet is initialised with an empty Python list object,
        # the DictReader really needs None to automatically populate the fieldnames instance variable.
        # However, the DictReader should always do so, if a header line is expected since
        # otherwise, the header line would get repeated as data line.

        if self.field_names and not self.header:
            csv_field_names = self.field_names
        else:
            csv_field_names = None

        if self.file_type:
            csv_file_type = self.file_type
        else:
            csv_file_type = None

        # For Python2.7, the open() function has to use the binary 'b' flag.
        # For Python3, the open() function has to use newline=''.
        self._csv_reader_file = open(file=self.file_path, mode='r', newline='')
        self._csv_reader_object = DictReader(
            f=self._csv_reader_file,
            fieldnames=csv_field_names,
            dialect=csv_file_type)

        # Automatically set the field names from the DictReader,
        # if the field_names list is empty and if possible.

        if self._csv_reader_object.fieldnames and not self.field_names:
            self.field_names.extend(self._csv_reader_object.fieldnames)

        return

    def csv_reader_next(self):
        """Read the next line of a CSV file linked to a C{bsf.annotation.AnnotationSheet} object.

        @return: Python C{dict} of column key and row value data
        @rtype: dict[str, str]
        """
        return self._csv_reader_object.next()

    def csv_reader_close(self):
        """Close a Comma-Separated Value (CSV) file linked to a C{bsf.annotation.AnnotationSheet} object for reading.
        """
        self._csv_reader_object = None
        self._csv_reader_file.close()
        self._csv_reader_file = None

        return

    def csv_writer_open(self):
        """Open a Comma-Separated Value (CSV) file linked to a C{bsf.annotation.AnnotationSheet} object for writing,
        initialise a Python C{csv.DictWriter} object and write the header line if one has been defined.
        """
        if not self.file_path:
            raise Exception('Cannot write an AnnotationSheet without a valid file_name.')

        if not self.field_names:
            raise Exception('A csv.DictWriter object requires a Python list of field_names.')

        if self.file_type:
            csv_file_type = self.file_type
        else:
            csv_file_type = None

        # For Python2.7, the open() function has to use the binary 'b' flag.
        # For Python3, the open() function has to use newline=''.
        self._csv_writer_file = open(file=self.file_path, mode='w', newline='')
        self._csv_writer_object = DictWriter(
            f=self._csv_writer_file,
            fieldnames=self.field_names,
            dialect=csv_file_type)

        if self.header:
            self._csv_writer_object.writeheader()

        return

    def csv_writer_next(self, row_dict):
        """Write the next line of a CSV file linked to a C{bsf.annotation.AnnotationSheet} object.

        @param row_dict: Row Python C{dict}
        @type row_dict: dict[str, str]
        """
        self._csv_writer_object.writerow(rowdict=row_dict)

        return

    def csv_writer_close(self):
        """Close a Comma-Separated Value (CSV) file linked to a C{bsf.annotation.AnnotationSheet} object for writing.
        """
        self._csv_writer_object = None
        self._csv_writer_file.close()
        self._csv_writer_file = None

        return

    def get_boolean(self, row_dict, key):
        """Get the Boolean state of a cell of an C{AnnotationSheet}.

        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str]
        @param key: Key
        @param key: str
        @return: Boolean
        @rtype: bool | None
        """
        if key in row_dict:
            value = row_dict[key].lower()
            if value in self._boolean_states:
                return self._boolean_states[value]
            else:
                raise ValueError('Value ' + repr(row_dict[key]) + ' in field ' + repr(key) +
                                 ' of AnnotationSheet ' + repr(self.name) + ' is not a boolean.')

    def sort(self):
        """Sort a C{bsf.annotation.AnnotationSheet}.

        This method has to implemented in the subclass,
        as it requires information about field-specific sorting.
        """
        warnings.warn(
            'Sorting of AnnotationSheet objects has to implemented in the sub-class.',
            UserWarning)

        return

    def validate(self):
        """Validate a C{bsf.annotation.AnnotationSheet}.

        @return: Warning messages
        @rtype: str
        """
        messages = str()
        row_number = 0

        for row_dict in self.row_dicts:
            row_number += 1
            for field_name in self.field_names:
                if field_name in self.test_methods:
                    for class_method_pointer in self.test_methods[field_name]:
                        # Only validate fields, for which instructions exist in the test_methods dict.
                        messages += class_method_pointer(row_number, row_dict, field_name)

        return messages

    def adjust_field_names(self):
        """Adjust the Python C{list} of Python C{str} field names to the keys used in Python C{dict} (row) objects.
        """
        field_names = list()

        for row_dict in self.row_dicts:
            for key in row_dict:
                if key not in field_names:
                    field_names.append(key)

        del self.field_names[:]

        self.field_names.extend(sorted(field_names))

        return

    def to_file_path(self, adjust_field_names=None):
        """Write a C{bsf.annotation.AnnotationSheet} to a file path.

        @param adjust_field_names: Clear and adjust the Python C{list} of Python C{str} field name objects
        @type adjust_field_names: bool
        """
        if adjust_field_names:
            self.adjust_field_names()

        self.csv_writer_open()

        for row_dict in self.row_dicts:
            self.csv_writer_next(row_dict=row_dict)

        self.csv_writer_close()

        return
