"""bsf.annotation

A package of classes and methods modelling Comma-Separated Value (CSV) and Tab-Separated Value (TSV) annotation files.
"""

#
# Copyright 2013 - 2014 Michael K. Schuster
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


import csv
import re
import warnings


class AnnotationSheet(object):
    """The C{AnnotationSheet} class represents comma-separated value (CSV) files.

    This class is a bit unusual in that values of instance variables around the
    C{csv.DictReader} and C{csv.DictWriter} classes can be initialised from a corresponding set of class variables.
    Therefore, sub-classes with fixed defaults can be defined, while generic objects can be initialised directly
    via the C{AnnotationSheet.__init__} method without the need to create an C{AnnotationSheet} sub-class in code.

    Attributes:
    @cvar _regular_expression_non_alpha: Regular expression for non-alphanumeric characters
    @type _regular_expression_non_alpha: __Regex
    @cvar _regular_expression_non_numeric: Regular expression for non-numeric characters
    @type _regular_expression_non_numeric: __Regex
    @cvar _regular_expression_non_sequence: Regular expression for non-sequence characters
    @type _regular_expression_non_sequence: __Regex
    @cvar _regular_expression_multiple_underscore: Regular expression for multiple underscore characters
    @type _regular_expression_multiple_underscore: __Regex
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar field_names: Python C{list} of Python C{str} (field name) objects
    @type field_names: list
    @ivar test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type test_methods: dict
    @ivar row_dicts: Python C{list} of Python C{dict} objects
    @type row_dicts: list
    """

    _regular_expression_non_alpha = re.compile(pattern='\W')
    _regular_expression_non_numeric = re.compile(pattern='\D')
    _regular_expression_non_sequence = re.compile(pattern='[^ACGTacgt]')
    _regular_expression_multiple_underscore = re.compile(pattern='_{2,}')

    _file_type = 'excel'
    _header_line = True
    _field_names = list()
    _test_methods = dict()

    @classmethod
    def check_column(cls, row_number, row_dict, column_name, require_column=True, require_value=True):
        """Check for a column name and return its associated value, if any.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value
        @type require_value: bool
        @return: Python C{tuple} of warning message and column value
        @rtype: tuple
        """

        message = str()

        if column_name not in row_dict:
            if require_column:
                message += 'Column name {!r} not in dict for row {}.\n'.format(column_name, row_number)
            return message, None

        if row_dict[column_name]:
            return message, row_dict[column_name]
        else:
            if require_value:
                message += 'Column name {!r} without value in row {}.\n'.format(column_name, row_number)
            return message, None

    @classmethod
    def check_alphanumeric(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{alphanumeric}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value contains only alphanumeric characters.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            match = re.search(pattern=cls._regular_expression_non_alpha, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with non-alphanumeric characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_numeric(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{numeric}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value contains only numeric characters.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            match = re.search(pattern=cls._regular_expression_non_numeric, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with non-numeric characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_sequence(cls, row_number, row_dict, column_name, require_column=True, require_value=True):
        """Validate a particular column value as I{sequence}.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @param require_column: Require the column_name to be defined in the row_dict
        @type require_column: bool
        @param require_value: Require a value
        @type require_value: bool
        @return: Warning messages
        @rtype: str
        """
        messages, column_value = cls.check_column(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=require_column,
            require_value=require_value)

        if column_value:
            match = re.search(pattern=cls._regular_expression_non_sequence, string=row_dict[column_name])
            if match:
                messages += 'Field {!r} in row {} contains a sequence {!r} with illegal characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_sequence_mandatory(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{mandatory sequence}.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls.check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=True)

    @classmethod
    def check_sequence_optional(cls, row_number, row_dict, column_name):
        """Validate a particular column value as I{optional sequence}.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """
        return cls.check_sequence(
            row_number=row_number,
            row_dict=row_dict,
            column_name=column_name,
            require_column=True,
            require_value=False)

    @classmethod
    def check_underscore_leading(cls, row_number, row_dict, column_name):
        """Validate a particular column value for I{leading underscore characters}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has no leading underscore.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            if column_value[:1] == '_':
                messages += 'Column {!r} in row {} contains a value {!r} with a leading underscore character.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_underscore_trailing(cls, row_number, row_dict, column_name):
        """Validate a particular column value for I{trailing underscore characters}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has no trailing underscore.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            if column_value[-1:] == '_':
                messages += 'Column {!r} in row {} contains a value {!r} with a trailing underscore character.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_underscore_multiple(cls, row_number, row_dict, column_name):
        """Validate a particular column value for I{multiple underscore characters}.

        Check that the particular column name key exists in the row dictionary and that
        its associated value has not multiple underscore adjacent to each other.

        @param row_number: Row number for warning messages
        @type row_number: int
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param column_name: Column name
        @type column_name: str
        @return: Warning messages
        @rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            match = re.search(pattern=cls._regular_expression_multiple_underscore, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with multiple underscore characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def from_file_path(cls, file_path=None, file_type=None, name=None):
        """Construct an C{AnnotationSheet} from a comma-separated value (CSV) file.

        This method reads the whole CSV file at once and stores a Python C{list} of Python C{dict} objects
        representing each row. For large files, the C{AnnotationSheet.csv_reader_open},
        C{AnnotationSheet.csv_reader_next} and C{AnnotationSheet.csv_reader_close} methods should be called explicitly.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (i.e. I{excel} or I{excel_tab} defined in the C{csv.Dialect} class)
        @type file_type: str
        @return: C{AnnotationSheet}
        @rtype: AnnotationSheet
        """

        annotation_sheet = cls(file_path=file_path, file_type=file_type, name=name)

        annotation_sheet.csv_reader_open()

        for row_dict in annotation_sheet._csv_reader_object:
            annotation_sheet.row_dicts.append(row_dict)

        annotation_sheet.csv_reader_close()

        return annotation_sheet

    def __init__(self, file_path=None, file_type=None, name=None, header=None, field_names=None, test_methods=None,
                 row_dicts=None):
        """Initialise an C{AnnotationSheet} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (i.e. I{excel} or I{excel_tab} defined in the C{csv.Dialect} class)
        @type file_type: str
        @param name: Name
        @type name: str
        @param header: Header line
        @type header: bool
        @param field_names: Python C{list} of Python C{str} (field name) objects
        @type field_names: list
        @param test_methods: Python C{dict} of Python C{str} (field name) key data and
            Python C{list} of Python C{function} value data
        @type test_methods: dict
        @param row_dicts: Python C{list} of Python C{dict} objects
        @type row_dicts: list
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = self._file_type

        if name:
            self.name = name
        else:
            self.name = str()

        if header is None:
            self.header = self._header_line
        else:
            self.header = header

        if field_names:
            self.field_names = field_names
        else:
            self.field_names = self._field_names

        if test_methods:
            self.test_methods = test_methods
        else:
            self.test_methods = self._test_methods

        if row_dicts:
            self.row_dicts = row_dicts
        else:
            self.row_dicts = list()

        self._csv_reader_file = None
        self._csv_reader_object = None

        self._csv_writer_file = None
        self._csv_writer_object = None

    def csv_reader_open(self):
        """Open a Comma-Separated Value (CSV) file linked to an C{AnnotationSheet} object for reading
        and initialise a Python C{csv.DictReader} object.
        """

        # Although the AnnotationSheet is initialised with an empty Python list object,
        # the DictReader really needs None to automatically populate the fieldnames instance variable.
        # However, the DictReader should always do so, if a header line is expected since
        # otherwise, the header line would get repeated as data line.

        if len(self.field_names) and not self.header:
            csv_field_names = self.field_names
        else:
            csv_field_names = None

        if len(self.file_type):
            csv_file_type = self.file_type
        else:
            csv_file_type = None

        self._csv_reader_file = open(name=self.file_path, mode='rb')
        self._csv_reader_object = csv.DictReader(
            f=self._csv_reader_file,
            fieldnames=csv_field_names,
            dialect=csv_file_type)

        # Automatically set the field names from the DictReader,
        # if the field_names list is empty and if possible.

        if len(self._csv_reader_object.fieldnames) and not len(self.field_names):
            self.field_names.extend(self._csv_reader_object.fieldnames)

    def csv_reader_next(self):
        """Read the next line of a CSV file linked to an C{AnnotationSheet} object.

        @return: Python C{dict} of column key and row value data
        @rtype: dict
        """

        return self._csv_reader_object.next()

    def csv_reader_close(self):
        """Close a Comma-Separated Value (CSV) file linked to an C{AnnotationSheet} object for reading.
        """

        self._csv_reader_object = None
        self._csv_reader_file.close()
        self._csv_reader_file = None

    def csv_writer_open(self):
        """Open a Comma-Separated Value (CSV) file linked to an C{AnnotationSheet} object for writing,
        initialise a Python C{csv.DictWriter} object and write the header line if one has been defined.
        """

        if not len(self.field_names):
            raise Exception("A csv.DictWriter object requires a Python list of field_names.")

        if len(self.file_type):
            csv_file_type = self.file_type
        else:
            csv_file_type = None

        self._csv_writer_file = open(name=self.file_path, mode='wb')
        self._csv_writer_object = csv.DictWriter(
            f=self._csv_writer_file,
            fieldnames=self.field_names,
            dialect=csv_file_type)

        if self.header:
            self._csv_writer_object.writeheader()

    def csv_writer_next(self, row_dict):
        """Write the next line of a CSV file linked to an C{AnnotationSheet} object.

        @param row_dict: Row Python C{dict}
        @type row_dict: dict
        """

        self._csv_writer_object.writerow(rowdict=row_dict)

    def csv_writer_close(self):
        """Close a Comma-Separated Value (CSV) file linked to an C{AnnotationSheet} object for writing.
        """

        self._csv_writer_object = None
        self._csv_writer_file.close()
        self._csv_writer_file = None

    def sort(self):
        """Sort an C{AnnotationSheet}.

        This method has to implemented in the sub-class,
        as it requires information about field-specific sorting.
        """

        warnings.warn(
            'Sorting of AnnotationSheet objects has to implemented in the sub-class.',
            UserWarning)

    def validate(self):
        """Validate an C{AnnotationSheet}.

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
                        messages += class_method_pointer(
                            row_number=row_number,
                            row_dict=row_dict,
                            column_name=field_name
                        )

        return messages

    def write_to_file(self):
        """Write an C{AnnotationSheet} to a file path.
        """
        # TODO: This method should be renamed to to_file_path.
        self.csv_writer_open()

        for row_dict in self.row_dicts:
            self.csv_writer_next(row_dict=row_dict)

        self.csv_writer_close()


class BamIndexDecoderSheet(AnnotationSheet):
    """The C{BamIndexDecoderSheet} class represents a Tab-Separated Value (TSV) table of
    library information for the C{IlluminaToBamTools.BamIndexDecoder} C{Analysis}.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    """

    _file_type = 'excel-tab'

    # The field names are defined in the IndexDecoder.java source file.
    # https://github.com/wtsi-npg/illumina2bam/blob/devel/src/uk/ac/sanger/npg/picard/IndexDecoder.java

    _header_line = True

    _field_names = [
        'barcode_sequence',
        'barcode_name',
        'library_name',
        'sample_name',
        'description'
    ]

    _test_methods = dict()


class LibraryAnnotationSheet(AnnotationSheet):
    """The C{LibraryAnnotationSheet} class represents a Comma-Separated Value (CSV) table of
    library information for the C{IlluminaToBamTools.BamIndexDecoder} C{Analysis}.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'lane',  # Lane number
        'barcode_sequence_1',  # Index read sequence 1
        'barcode_sequence_2',  # Index read sequence 2
        'sample_name',  # Sample name (alphanumeric including '_' characters)
        'library_name',  # Library name (alphanumeric including '_' characters)
        'library_size'  # Library size (numeric)
    ]

    _test_methods = dict(
        lane=[
            AnnotationSheet.check_alphanumeric
        ],
        barcode_sequence_1=[
            AnnotationSheet.check_sequence_optional
        ],
        barcode_sequence_2=[
            AnnotationSheet.check_sequence_optional
        ],
        sample_name=[
            AnnotationSheet.check_alphanumeric,
            AnnotationSheet.check_underscore_leading,
            AnnotationSheet.check_underscore_trailing,
            AnnotationSheet.check_underscore_multiple
        ],
        library_name=[
            AnnotationSheet.check_alphanumeric,
            AnnotationSheet.check_underscore_leading,
            AnnotationSheet.check_underscore_trailing,
            AnnotationSheet.check_underscore_multiple
        ]
    )

    def validate(self, lanes=8):
        """
        Validate a C{LibraryAnnotationSheet}.

        @param lanes: Number of lanes to validate
        @type lanes: int
        @return: Warning messages
        @rtype: str
        """

        messages = str()

        # Check the header line via the pre-defined field names.

        for index in range(0, len(self._field_names)):
            if not self.field_names[index]:
                messages += 'Column with name {!r} is missing from the header line.\n'. \
                    format(self._field_names[index])

            if self.field_names[index] != self._field_names[index]:
                messages += 'Column name {!r} in the header line does not match template {!r}.\n'. \
                    format(self.field_names[index], self._field_names[index])

        # Validate the field values for alphanumeric or sequence grade in the context of the
        # AnnotationSheet super-class.

        messages += super(LibraryAnnotationSheet, self).validate()

        lane_index = dict()

        row_number = 0

        for row_dict in self.row_dicts:

            row_number += 1

            # Check that all required fields are defined.

            if row_dict['lane'] in lane_index:
                barcode_dict, sample_dict, library_name = lane_index[row_dict['lane']]
            else:
                barcode_dict = dict()
                sample_dict = dict()
                library_name = str(row_dict['library_name'])
                lane_index[row_dict['lane']] = (barcode_dict, sample_dict, library_name)

            barcode_sequence = str()

            if 'barcode_sequence_1' in row_dict and row_dict['barcode_sequence_1']:
                barcode_sequence += row_dict['barcode_sequence_1']
            else:
                barcode_sequence += '-NoIndex-'

            if 'barcode_sequence_2' in row_dict and row_dict['barcode_sequence_2']:
                barcode_sequence += row_dict['barcode_sequence_2']
            else:
                barcode_sequence += '-NoIndex-'

            if barcode_sequence in barcode_dict:
                messages += 'Barcode sequence {!r} from row {} duplicated in row {}.\n'. \
                    format(barcode_sequence, barcode_dict[barcode_sequence], row_number)
            else:
                barcode_dict[barcode_sequence] = row_number

            sample_name = str()
            sample_name += row_dict['sample_name']

            if sample_name in sample_dict:
                messages += 'Sample name {!r} from row {} duplicated in row {}.\n'. \
                    format(sample_name, sample_dict[sample_name], row_number)
            else:
                sample_dict[sample_name] = row_number

            if library_name != row_dict['library_name']:
                messages += 'Library name {!r} in row {} does not match previous name {!r}.\n'. \
                    format(row_dict['library_name'], row_number, library_name)

        for lane_number in range(0 + 1, lanes + 1):
            lane_string = str(lane_number)

            # Check that all lanes have annotation.
            if lane_string not in lane_index:
                messages += 'No annotation for lane number {!r}.\n'.format(lane_number)
                continue

            barcode_dict, sample_dict, library_name = lane_index[lane_string]

            # Check that all or none of the rows has barcode sequence 1 or 2 populated.
            no_index_1 = 0
            no_index_2 = 0
            for key in barcode_dict.keys():
                if key[:9] == '-NoIndex-':
                    no_index_1 += 1
                if key[-9:] == '-NoIndex-':
                    no_index_2 += 1

            if not (no_index_1 == 0 or no_index_1 == len(barcode_dict)):
                messages += 'Some empty barcode_sequence_1 fields in lane {}.\n'.format(lane_number)
            if not (no_index_2 == 0 or no_index_2 == len(barcode_dict)):
                messages += 'Some empty barcode_sequence_2 fields in lane {}.\n'.format(lane_number)

            # Check that all barcode sequences have the same length.
            # This test also finds cases of missing sequences tested for above.
            key_list = barcode_dict.keys()
            key_length = len(key_list[0])
            for key in key_list[1:]:
                if len(key) != key_length:
                    messages += 'Mismatching barcode sequence lengths in lane {}.\n'.format(lane_number)

        return messages


class SampleAnnotationSheet(AnnotationSheet):
    """The C{SampleAnnotationSheet} class represents a Comma-Separated Value (CSV) table of sample information
    after running the C{IlluminaToBamTools.BamIndexDecoder} C{Analysis}.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'ProcessedRunFolder', 'Project', 'Sample',
        'Reads1', 'File1', 'Reads2', 'File2',
        'LibrarySize', 'Barcode1', 'Barcode2'
    ]

    _test_methods = dict(
        ProcessedRunFolder=[
            AnnotationSheet.check_alphanumeric
        ],
        Project=[
            AnnotationSheet.check_alphanumeric
        ],
        Sample=[
            AnnotationSheet.check_alphanumeric
        ],
        Reads1=[
            AnnotationSheet.check_alphanumeric
        ],
        Reads2=[
            AnnotationSheet.check_alphanumeric
        ],
        Barcode1=[
            AnnotationSheet.check_sequence_optional
        ],
        Barcode2=[
            AnnotationSheet.check_sequence_optional
        ]
    )


class ChIPSeqDiffBindSheet(AnnotationSheet):
    """ChIP-Seq Bioconductor DiffBind annotation sheet class.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'SampleID', 'Tissue', 'Factor', 'Condition', 'Treatment', 'Replicate',
        'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat'
    ]

    _test_methods = dict(
        SampleID=[
            AnnotationSheet.check_alphanumeric
        ],
        Tissue=[
            AnnotationSheet.check_alphanumeric
        ],
        Factor=[
            AnnotationSheet.check_alphanumeric
        ],
        Condition=[
            AnnotationSheet.check_alphanumeric
        ],
        Treatment=[
            AnnotationSheet.check_alphanumeric
        ],
        Replicate=[
            AnnotationSheet.check_numeric
        ],
        ControlID=[
            AnnotationSheet.check_alphanumeric
        ],
        PeakCaller=[
            AnnotationSheet.check_alphanumeric
        ],
        PeakFormat=[
            AnnotationSheet.check_alphanumeric
        ]
    )

    def sort(self):
        """Sort by I{Tissue}, I{Factor}, I{Condition}, I{Treatment} and I{Replicate} columns.
        """

        self.row_dicts.sort(
            cmp=lambda x, y:
            cmp(x['Tissue'], y['Tissue']) or
            cmp(x['Factor'], y['Factor']) or
            cmp(x['Condition'], y['Condition']) or
            cmp(x['Treatment'], y['Treatment']) or
            cmp(int(x['Replicate']), int(y['Replicate'])))

    def write_to_file(self):
        """Write a C{ChIPSeqDiffBindSheet} to a file.
        """

        # Override the method from the super-class to automatically sort before writing to a file.

        self.sort()
        super(ChIPSeqDiffBindSheet, self).write_to_file()


class TuxedoSamplePairSheet(AnnotationSheet):
    """The C{TuxedoSamplePairSheet} class represents C{Sample} pairs defined by
    the C{bsf_rnaseq_process_cuffdiff.R} script.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    """

    _file_type = "excel-tab"

    _field_names = ["V1", "V2"]

    _test_methods = dict()
