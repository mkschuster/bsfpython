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
"""The :py:mod:`bsf.annotation` module provides classes modelling
Comma-Separated Value (CSV) and
Tab-Separated Value (TSV) annotation files.
"""
import logging
import re
from csv import DictReader, DictWriter
from io import TextIOWrapper
from typing import Callable, Dict, List, Optional

module_logger = logging.getLogger(name=__name__)


class AnnotationSheet(object):
    """The :py:class:`bsf.annotation.AnnotationSheet` class represents comma-separated value (CSV) files.

    This class is a bit unusual in that values of instance variables around the
    :py:class:`csv.DictReader` and :py:class:`csv.DictWriter` classes can be initialised from a corresponding set
    of class variables.
    Therefore, subclasses with fixed defaults can be defined, while generic objects can be initialised directly
    via the :py:meth:`bsf.annotation.AnnotationSheet.__init__` method without the need to create a
    :py:class:`bsf.annotation.AnnotationSheet` subclass in code.

    :ivar file_path: A file path.
    :type file_path: str | None
    :ivar file_type: A file type (i.e., :literal:`excel` or :literal:`excel-tab`
        defined in the :py:class:`csv.Dialect` class).
    :type file_type: str
    :ivar name: A name.
    :type name: str | None
    :ivar field_names: A Python :py:class:`list` object of Python :py:class:`str` (field name) objects.
    :type field_names: list[str]
    :ivar test_methods: A Python :py:class:`dict` of Python :py:class:`str` (field name) key and
        Python :py:class:`list` of Python :py:class:`Callable` value objects.
    :type test_methods: Dict[str, List[Callable[[int, Dict[str, str], str], str]]]
    :ivar row_dicts: A Python :py:class:`list` object of Python :py:class:`dict` (row) objects.
    :type row_dicts: list[dict[str, str]]
    """

    # Regular expression for non-alphanumeric characters
    _regular_expression_non_alpha = re.compile(pattern=r'\W')

    # Regular expression for non-numeric characters
    _regular_expression_non_numeric = re.compile(pattern=r'\D')

    # Regular expression for non-sequence characters
    _regular_expression_non_sequence = re.compile(pattern=r'[^ACGTacgt]')

    # Regular expression for non-ambiguous sequence characters
    _regular_expression_non_ambiguous_sequence = re.compile(pattern=r'[^ACGTacgtWSMKRYwsmkryBDHVbdhvNn]')

    # Regular expression for multiple underscore characters
    _regular_expression_multiple_underscore = re.compile(pattern=r'_{2,}')

    # File type (i.e., 'excel' or 'excel-tab' defined in the csv.Dialect class)
    _file_type = 'excel'

    # Header line exists
    _header_line = True

    # Python list of Python str (field name) objects
    _field_names: List[str] = list()

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

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :param require_column: Require the column_name to be defined in the row_dict.
        :type require_column: bool
        :param require_value: Require a value, which also implies requiring the column_name.
        :type require_value: bool
        :return: A Python :py:class:`tuple` of Python :py:class:`str` (warning message) and
            Python :py:class:`str` (column value) objects.
        :rtype: (str, str) | (str, None)
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
        """Validate a particular column value as :literal:`alphanumeric`.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only alphanumeric characters.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :param require_column: Require the column_name to be defined in the row_dict.
        :type require_column: bool
        :param require_value: Require a value, which also implies requiring the column_name.
        :type require_value: bool
        :return: Warning messages.
        :rtype: str
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
        """Validate a particular column value as :literal:`alphanumeric`.

        Neither the column nor the value needs existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_alphanumeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name)

    @classmethod
    def check_alphanumeric_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`alphanumeric`.

        The column, but not the value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_alphanumeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_alphanumeric_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`alphanumeric`.

        Both, the column and value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_alphanumeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def _check_numeric(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as :literal:`numeric`.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only numeric characters.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :param require_column: Require the column_name to be defined in the row_dict.
        :type require_column: bool
        :param require_value: Require a value, which also implies requiring the column_name.
        :type require_value: bool
        :return: Warning messages.
        :rtype: str
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
        """Validate a particular column value as :literal:`numeric`.

        Neither the column nor the value needs existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_numeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=False,
            require_value=False)

    @classmethod
    def check_numeric_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`numeric`.

        The column, but not the value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_numeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_numeric_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`numeric`.

        Both, the column and value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_numeric(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def _check_sequence(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as :literal:`IUPAC sequence`.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only valid IUPAC sequence characters.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :param require_column: Require the column_name to be defined in the row_dict.
        :type require_column: bool
        :param require_value: Require a value.
        :type require_value: bool
        :return: Warning messages.
        :rtype: str
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
        """Validate a particular column value as :literal:`IUPAC sequence`.

        Neither the column nor the value needs existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=False,
            require_value=False)

    @classmethod
    def check_sequence_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`IUPAC sequence`.

        The column, but not the value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_sequence_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`IUPAC sequence`.

        Both, the column and value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def _check_ambiguous_sequence(cls, row_number, row_dict, column_name, require_column=False, require_value=False):
        """Validate a particular column value as :literal:`IUPAC ambiguous sequence`.

        If the particular column name key exists in the row dictionary and if it has
        an associated value, it must contain only valid IUPAC ambiguity sequence characters.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :param require_column: Require the column_name to be defined in the row_dict.
        :type require_column: bool
        :param require_value: Require a value.
        :type require_value: bool
        :return: Warning messages.
        :rtype: str
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
        """Validate a particular column value as :literal:`IUPAC ambiguous sequence`.

        Neither the column nor the value needs existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_ambiguous_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=False,
            require_value=False)

    @classmethod
    def check_ambiguous_sequence_column(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`IUPAC ambiguous sequence`.

        The column, but not the value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_ambiguous_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_ambiguous_sequence_value(cls, row_number, row_dict, column_name):
        """Validate a particular column value as :literal:`IUPAC ambiguous sequence`.

        Both, the column and value need existing.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
        """
        return cls._check_ambiguous_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def check_underscore_leading(cls, row_number, row_dict, column_name):
        """Validate a particular column value for :literal:`leading underscore characters`.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has no leading underscore.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
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
        """Validate a particular column value for :literal:`trailing underscore characters`.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has no trailing underscore.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
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
        """Validate a particular column value for :literal:`multiple underscore characters`.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has not multiple underscore adjacent to each other.

        :param row_number: A current row number for warning messages.
        :type row_number: int
        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param column_name: A column name.
        :type column_name: str
        :return: Warning messages.
        :rtype: str
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
        """Construct a :py:class:`bsf.annotation.AnnotationSheet` from a comma-separated value (CSV) file.

        This method reads the whole CSV file at once and stores a Python :py:class:`list` of
        Python :py:class:`dict` objects representing each row.
        For large files, the :py:meth:`bsf.annotation.AnnotationSheet.csv_reader_open`,
        :py:meth:`bsf.annotation.AnnotationSheet.csv_reader_next` and
        :py:meth:`bsf.annotation.AnnotationSheet.csv_reader_close` methods should be called explicitly.

        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`excel` or :literal:`excel_tab` defined in the
            :py:class:`csv.Dialect` class).
        :type file_type: str | None
        :param name: A name.
        :type name: str | None
        :return: A :py:class:`bsf.annotation.AnnotationSheet` object.
        :rtype: AnnotationSheet
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
        """Initialise a :py:class:`bsf.annotation.AnnotationSheet` object.

        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`excel` or :literal:`excel_tab` defined in the
            :py:class:`csv.Dialect` class).
        :type file_type: str | None
        :param name: A name.
        :type name: str | None
        :param header: A header line exists.
        :type header: bool | None
        :param field_names: A Python :py:class:`list` object of Python :py:class:`str` (field name) objects.
        :type field_names: list[str] | None
        :param test_methods: A Python :py:class:`dict` of Python :py:class:`str` (field name) key and
        Python :py:class:`list` of Python :py:class:`Callable` value objects.
        :type test_methods: Optional[Dict[str, List[Callable[[int, Dict[str, str], str], str]]]]
        :param row_dicts: A Python :py:class:`list` object of Python :py:class:`dict` (row) objects.
        :type row_dicts: list[dict[str, str]] | None
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

        self._csv_reader_text_io: Optional[TextIOWrapper] = None
        self._csv_reader_object: Optional[DictReader] = None
        self._csv_writer_text_io: Optional[TextIOWrapper] = None
        self._csv_writer_object: Optional[DictWriter] = None

        return

    def trace(self, level=1):
        """Trace a :py:class:`bsf.annotation.SampleAnnotationSheet` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
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
        """Open a Comma-Separated Value (CSV) file linked to a :py:class:`bsf.annotation.AnnotationSheet` object
        for reading and initialise a Python :py:class:`csv.DictReader` object.
        """
        if not self.file_path:
            raise Exception("Cannot read an AnnotationSheet without a valid 'file_name'.")

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
        self._csv_reader_text_io = open(file=self.file_path, mode='r', newline='')
        self._csv_reader_object = DictReader(
            f=self._csv_reader_text_io,
            fieldnames=csv_field_names,
            dialect=csv_file_type)

        # Automatically set the field names from the DictReader,
        # if the field_names list is empty and if possible.

        if self._csv_reader_object.fieldnames and not self.field_names:
            self.field_names.extend(self._csv_reader_object.fieldnames)

        return

    def csv_reader_close(self):
        """Close a Comma-Separated Value (CSV) file linked to a :py:class:`bsf.annotation.AnnotationSheet` object
        for reading.
        """
        self._csv_reader_object = None
        self._csv_reader_text_io.close()
        self._csv_reader_text_io = None

        return

    def csv_writer_open(self):
        """Open a Comma-Separated Value (CSV) file linked to a :py:class:`bsf.annotation.AnnotationSheet` object
        for writing, initialise a Python :py:class:`csv.DictWriter` object and write the header line
        if one has been defined.
        """
        if not self.file_path:
            raise Exception("Cannot write an AnnotationSheet without a valid 'file_name'.")

        if not self.field_names:
            raise Exception("A csv.DictWriter object requires a Python list of 'field_names'.")

        if self.file_type:
            csv_file_type = self.file_type
        else:
            csv_file_type = None

        # For Python2.7, the open() function has to use the binary 'b' flag.
        # For Python3, the open() function has to use newline=''.
        self._csv_writer_text_io = open(file=self.file_path, mode='w', newline='')
        self._csv_writer_object = DictWriter(
            f=self._csv_writer_text_io,
            fieldnames=self.field_names,
            dialect=csv_file_type)

        if self.header:
            self._csv_writer_object.writeheader()

        return

    def csv_writer_write_row(self, row_dict):
        """Write the next row of a CSV file linked to a :py:class:`bsf.annotation.AnnotationSheet` object.

        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictWriter` object.
        :type row_dict: dict[str, str]
        """
        self._csv_writer_object.writerow(rowdict=row_dict)

        return

    def csv_writer_close(self):
        """Close a Comma-Separated Value (CSV) file linked to a :py:class:`bsf.annotation.AnnotationSheet` object
        for writing.
        """
        self._csv_writer_object = None
        self._csv_writer_text_io.close()
        self._csv_writer_text_io = None

        return

    def get_boolean(self, row_dict, key):
        """Get the Boolean state of a cell of an :py:class:`AnnotationSheet` object.

        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:class:`csv.DictReader` object.
        :type row_dict: dict[str, str]
        :param key: A key.
        :param key: str
        :return: A Python :py:class:`bool` object.
        :rtype: bool
        """
        if key in row_dict:
            value = row_dict[key].lower()
            if value in self._boolean_states:
                return self._boolean_states[value]
            else:
                raise ValueError(
                    f'Value {row_dict[key]!r} in field {key!r} of AnnotationSheet {self.name!r} is not a boolean.')

    def sort(self):
        """Sort a :py:class:`bsf.annotation.AnnotationSheet` object.

        This method has to implemented in the subclass,
        as it requires information about field-specific sorting.
        """
        module_logger.warning(
            'Sorting of AnnotationSheet objects has to implemented in the %r sub-class.',
            self.__class__.__name__)

        return

    def validate(self):
        """Validate a :py:class:`bsf.annotation.AnnotationSheet` object.

        :return: Warning messages.
        :rtype: str
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
        """Adjust the Python :py:class:`list` object of Python :py:class:`str` (field name) objects to the keys used in
        Python :py:class:`dict` (row) objects.
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
        """Write a :py:class:`bsf.annotation.AnnotationSheet` object to a file path.

        :param adjust_field_names: Clear and adjust the Python :py:class:`list` of
            Python :py:class:`str` (field name) objects.
        :type adjust_field_names: bool | None
        """
        if adjust_field_names:
            self.adjust_field_names()

        self.csv_writer_open()

        for row_dict in self.row_dicts:
            self.csv_writer_write_row(row_dict=row_dict)

        self.csv_writer_close()

        return
