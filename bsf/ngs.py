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
"""The :py:mod:`bsf.ngs` module provides classes modelling next-generation sequencing (NGS) data directories and files.
"""
import logging
import os
import re
from typing import Callable, Dict, List, Optional
from weakref import ReferenceType

from bsf.annotation import AnnotationSheet
from bsf.standards import Configuration

module_logger = logging.getLogger(name=__name__)


class NextGenerationBase(object):
    """The :py:class:`bsf.ngs.NextGenerationBase` class represents a super-class for
    Next Generation Sequencing (NGS) objects.

    The :py:attr:`bsf.ngs.NextGenerationBase.file_type` attribute accepts the following values:
        - :literal:`CASAVA`: FASTQ file after post-processing with CASAVA
        - :literal:`External`: other data files

    :ivar name: A name.
    :type name: str | None
    :ivar file_path: A file path.
    :type file_path: str | None
    :ivar file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
    :type file_type: str | None
    :ivar annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
        Python :py:class:`list` value objects of Python :py:class:`str` objects.
    :type annotation_dict: dict[str, list[str]] | None
    """

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None):
        """Initialise a :py:class:`bsf.ngs.Reads` object.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        """
        super(NextGenerationBase, self).__init__()

        self.name = name
        self.file_path = file_path
        self.file_type = file_type

        if annotation_dict is None:
            self.annotation_dict = dict()
        else:
            self.annotation_dict = annotation_dict

        return

    def __eq__(self, other):
        """Test :py:class:`bsf.ngs.NextGenerationBase` objects for equality.

        :param other: A :py:class:`bsf.ngs.NextGenerationBase` object.
        :type other: NextGenerationBase
        :return: :py:const:`True` if equal, :py:const:`False` otherwise.
        :rtype: bool
        """
        assert isinstance(other, NextGenerationBase)

        if self is other:
            return True

        return \
            self.name == other.name \
            and self.file_path == other.file_path \
            and self.file_type == other.file_type \
            and self.annotation_dict == other.annotation_dict

    def __repr__(self):
        """Return a printable representation of a :py:class:`bsf.ngs.NextGenerationBase` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'NextGenerationBase(' \
            f'name={self.name!r} ' \
            f'file_path={self.file_path!r} ' \
            f'file_type={self.file_type!r} ' \
            f'annotation_dict={self.annotation_dict!r}' \
            f')'

    def add_annotation(self, key, value):
        """Add an annotation value under a key.

        :param key: An annotation key.
        :type key: str
        :param value: An annotation value.
        :type value: str
        """
        if self.annotation_dict is None:
            self.annotation_dict = dict()

        if key not in self.annotation_dict:
            self.annotation_dict[key] = list()

        value_list = self.annotation_dict[key]

        if value not in value_list:
            value_list.append(value)

        return

    def process_annotation(self, row_dict, key_list, prefix=None):
        """Process annotation from a Python :py:class:`dict` object of row entries of a
        Python :py:mod:`csv` module object.

        :param row_dict: A Python :py:class:`dict` object of row entries of a Python :py:mod:`csv` module object.
        :type row_dict: dict[str, str]
        :param key_list: A Python :py:class:`list` object of Python :py:class:`str` (key) objects in the row.
        :type key_list: list[str]
        :param prefix: An optional configuration prefix
            (e.g., :literal:`[Control] ReadsN`, :literal:`[Treatment] ReadsN`, :literal:`[Point N] ReadsN`, ...).
        :type prefix: str | None
        """
        pattern = self.__class__.__name__

        if prefix:
            pattern = ' '.join((prefix, pattern))

        for key in list(key_list):
            if key.startswith(pattern):
                key_list.remove(key)
                # Strip an eventual prefix and the class name from the annotation key.
                annotation_key = key[len(pattern):].strip()
                annotation_value = row_dict[key]

                # Exclude empty annotation keys and empty annotation values.
                if annotation_key and annotation_value:
                    self.add_annotation(key=annotation_key, value=annotation_value)

        return


class Reads(NextGenerationBase):
    """The :py:class:`bsf.ngs.Reads` class represents a file of Next-Generation Sequencing (NGS) reads.

    Typically, a :py:class:`bsf.ngs.Reads` object represents a FASTQ or unmapped BAM file.

    :ivar barcode: A barcode used for sample multiplexing.
    :type barcode: str | None
    :ivar lane: A lane number.
    :type lane: str | None
    :ivar read: A read number (e.g., :literal:`R1`, :literal:`R2`, ...).
    :type read: str | None
    :ivar chunk: A chunk number (e.g., :literal:`001`, :literal:`002`, ...).
    :type chunk: str | None
    :ivar weak_reference_paired_reads: A :py:class:`weakref.ReferenceType` object pointing at a
        :py:class:`bsf.ngs.PairedReads` object.
    :type weak_reference_paired_reads: ReferenceType | None
    """

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a :py:class:`bsf.ngs.Reads` object from a file path.

        For a :literal:`file_type` :literal:`CASAVA`, :py:attr:`bsf.ngs.Reads.file_path` obeys a
        :literal:`SampleName_Index_Lane_Read_Chunk` schema, so that :py:attr:`bsf.ngs.Reads.name`,
        :py:attr:`bsf.ngs.Reads.barcode`, :py:attr:`bsf.ngs.Reads.lane`, :py:class:`bsf.ngs.Reads.read` and
        :py:attr:`bsf.ngs.Reads.chunk` can be populated automatically.

        For :literal:`file_type` :literal:`External`, the instance variables need populating manually.

        :param file_path: A file path.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str
        :return: A :py:class:`bsf.ngs.Reads` object.
        :rtype: Reads
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_type == 'CASAVA':
            # CASAVA Reads obey a SampleName_Index_Lane_Read_Chunk.fastq.gz schema.
            # Recursively split off the extension, before splitting by underline character.

            extension = '.'
            while extension:
                file_name, extension = os.path.splitext(file_path)

            if isinstance(file_name, bytes):
                file_name = file_name.decode()

            component_list = file_name.split('_')
            # Since SampleName can contain underscores, list components need assigning from the end
            # i.e., via negative indices.

            return cls(
                name='_'.join(component_list[:-4]),  # Exclude the last four components
                file_path=file_path,
                file_type=file_type,
                annotation_dict=None,
                barcode=component_list[-4],
                lane=component_list[-3],
                read=component_list[-2],
                chunk=component_list[-1],
                weak_reference_paired_reads=None)
        else:
            raise Exception(f'Unsupported file_type: {file_type!r}')

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            barcode=None,
            lane=None,
            read=None,
            chunk=None,
            weak_reference_paired_reads=None):
        """Initialise a :py:class:`bsf.ngs.Reads` object.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        :param barcode: A barcode used for sample multiplexing.
        :type barcode: str | None
        :param lane: A lane number.
        :type lane: str | None
        :param read: A read number (e.g., :literal:`R1`, :literal:`R2`, ...).
        :type read: str | None
        :param chunk: A chunk number (e.g., :literal:`001`, :literal:`002`, ...).
        :type chunk: str | None
        :param weak_reference_paired_reads: A :py:class:`weakref.ReferenceType` object pointing at a
            :py:class:`bsf.ngs.PairedReads` object.
        :type weak_reference_paired_reads: ReferenceType | None
        """
        super(Reads, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        self.barcode = barcode
        self.lane = lane
        self.read = read
        self.chunk = chunk
        self.weak_reference_paired_reads = weak_reference_paired_reads

        return

    def __eq__(self, other):
        """Test :py:class:`bsf.ngs.Reads` objects for equality.

        :param other: A :py:class:`bsf.ngs.Reads` object.
        :type other: Reads
        :return: :py:const:`True` if equal, :py:const:`False` otherwise.
        :rtype: bool
        """
        assert isinstance(other, Reads)

        if self is other:
            return True

        return \
            super(Reads, self).__eq__(other=other) \
            and self.barcode == other.barcode \
            and self.lane == other.lane \
            and self.read == other.read \
            and self.chunk == other.chunk

    def __bool__(self):
        """Test a :py:class:`bsf.ngs.Reads` object for non-zero.

        :return: :py:const:`True` if non-zero (i.e., file_path or name are meaningfully defined),
            otherwise :py:const:`False`.
        :rtype: bool
        """
        return bool(self.file_path and self.name)

    def __repr__(self):
        """Return a printable representation of a :py:class:`bsf.ngs.Reads` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'Reads(' \
            f'name={self.name!r} ' \
            f'file_path={self.file_path!r} ' \
            f'file_type={self.file_type!r} ' \
            f'annotation_dict={self.annotation_dict!r} ' \
            f'barcode={self.barcode!r} ' \
            f'lane={self.lane!r} ' \
            f'read={self.read!r} ' \
            f'chunk={self.chunk!r}' \
            f')'

    def trace(self, level):
        """Trace a :py:class:`bsf.ngs.Reads` object.

        :param level: Indentation level.
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append(f'{indent}{self.__class__!r}\n')
        str_list.append(f'{indent}  weak_reference_paired_reads: {self.weak_reference_paired_reads!r}\n')
        str_list.append(f'{indent}  name:      {self.name!r}\n')
        str_list.append(f'{indent}  file_path: {self.file_path!r}\n')
        str_list.append(f'{indent}  file_type: {self.file_type!r}\n')
        str_list.append(f'{indent}  annotation_dict:\n')
        for annotation_key in sorted(self.annotation_dict):
            str_list.append(f'{indent}    {annotation_key!r} {self.annotation_dict[annotation_key]!r}\n')
        str_list.append(f'{indent}  barcode:   {self.barcode!r}\n')
        str_list.append(f'{indent}  lane:      {self.lane!r}\n')
        str_list.append(f'{indent}  read:      {self.read!r}\n')
        str_list.append(f'{indent}  chunk:     {self.chunk!r}\n')

        return str_list

    def match(self, reads):
        """Match :py:class:`bsf.ngs.Reads` objects.

        Two :py:class:`bsf.ngs.Reads` objects are identical, if all their instance variables match.

        :param reads: A second :py:class:`bsf.ngs.Reads` object.
        :type reads: Reads
        :return: :py:const:`True` if both objects match, :py:const:`False` otherwise.
        :rtype: bool
        """
        assert isinstance(reads, Reads)

        # Quick test first - if the objects are identical, the rest has to match.

        if self is reads:
            return True

        if self.name != reads.name:
            return False

        if self.file_path != reads.file_path:
            return False

        if self.file_type != reads.file_type:
            return False

        if self.annotation_dict != reads.annotation_dict:
            return False

        if self.barcode != reads.barcode:
            return False

        if self.lane != reads.lane:
            return False

        if self.read != reads.read:
            return False

        if self.chunk != reads.chunk:
            return False

        return True

    def match_paired(self, reads):
        """Test that two :py:class:`bsf.ngs.Reads` objects are compatible to form a
        :py:class:`bsf.ngs.PairedReads` object.

        For a :literal:`file_type` of :literal:`CASAVA`,
        all :py:class:`bsf.ngs.Reads` instance variables need to be equal, except
        :py:attr:`bsf.ngs.Reads.read` (i.e., :literal:`R1` or :literal:`R2`) and
        :py:attr:`bsf.ngs.Reads.file_path`, which also contains the :py:attr:`bsf.ngs.Reads.read` information.
        Hence, matching criteria need relaxing with respect to method :py:meth:`bsf.ngs.Reads.__eq__`.

        :param reads: A second :py:class:`bsf.ngs.Reads` object.
        :type reads: Reads
        :return: :py:const:`True` if both objects match, :py:const:`False` otherwise.
        :rtype: bool
        """
        assert isinstance(reads, Reads)

        if self.file_type == 'CASAVA':
            if self.name != reads.name:
                return False

            # Do not compare the 'file_path' instance variable.
            if self.file_path != reads.file_path:
                pass

            if self.file_type != reads.file_type:
                return False

            if self.barcode != reads.barcode:
                return False

            if self.lane != reads.lane:
                return False

            # Do not compare the 'read' instance variable.
            if self.read != reads.read:
                pass

            if self.chunk != reads.chunk:
                return False

            return True
        else:
            # It is difficult to match PairedReads outside of CASAVA conventions.

            return True


class PairedReads(NextGenerationBase):
    """The :py:class:`bsf.ngs.PairedReads` class represents a pair of :py:class:`bsf.ngs.Reads` objects.

    For the :py:attr:`bsf.ngs.Reads.file_type` :literal:`CASAVA`, this represents a read pair
    (i.e., :literal:`R1` and :literal:`R2`).

    :ivar reads_1: A first :py:class:`bsf.ngs.Reads` object.
    :type reads_1: Reads | None
    :ivar reads_2: A second :py:class:`bsf.ngs.Reads` object.
    :type reads_2: Reads | None
    :ivar exclude: Exclude from processing.
    :type exclude: bool | None
    :ivar index_1: An index sequence 1.
    :type index_1: str | None
    :ivar index_2: An index sequence 2.
    :type index_2: str | None
    :ivar read_group: A SAM read group (:literal:`@RG`) information.
    :type read_group: str | None
    :ivar weak_reference_sample: A :py:class:`weakref.ReferenceType` object pointing at a
        :py:class:`bsf.ngs.Sample` object.
    :type weak_reference_sample: ReferenceType | None
    """

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            reads_1=None,
            reads_2=None,
            exclude=None,
            index_1=None,
            index_2=None,
            read_group=None,
            weak_reference_sample=None):
        """Initialise a :py:class:`bsf.ngs.PairedReads` object.

        For the :py:attr:`bsf.ngs.Reads.file_type` :literal:`CASAVA`, the reads object will be
        automatically assigned on the basis of the :py:attr:`bsf.ngs.Reads.read`
        instance variable (i.e., :literal:`R1` or :literal:`R2`).
        Upon initialisation of this :py:class:`bsf.ngs.PairedReads` object,
        weak references (:py:class:`weakref.ReferenceType`) are set in the :py:class:`bsf.ngs.Reads` objects.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        :param reads_1: A first :py:class:`bsf.ngs.Reads` object.
        :type reads_1: Reads | None
        :param reads_2: A second :py:class:`bsf.ngs.Reads` object.
        :type reads_2: Reads | None
        :param index_1: An index sequence 1.
        :type index_1: str | None
        :param index_2: An index sequence 2.
        :type index_2: str | None
        :param exclude: Exclude from processing.
        :type exclude: bool | None
        :param read_group: A SAM read group (:literal:`@RG`) information.
        :type read_group: str | None
        :param weak_reference_sample: A :py:class:`weakref.ReferenceType` object pointing at a
            :py:class:`bsf.ngs.Sample` object.
        :type weak_reference_sample: ReferenceType | None
        :raise Exception: For :py:attr:`bsf.ngs.Reads.file_type} :literal:`CASAVA`, :literal:`R1` or :literal:`R2`
            must be set in the :py:class:`bsf.ngs.Reads` object.
        """
        if (reads_1 and reads_2) and (not reads_1.match_paired(reads=reads_2)):
            raise Exception('The Reads objects are not compatible to form a PairedReads object.')

        super(PairedReads, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        self.reads_1 = None
        self.reads_2 = None

        if reads_1:
            if reads_1.file_type == 'CASAVA':
                if reads_1.read == 'R1':
                    self.reads_1 = reads_1
                    reads_1.weak_reference_paired_reads = ReferenceType(self)
                elif reads_1.read == 'R2':
                    self.reads_2 = reads_1
                    reads_1.weak_reference_paired_reads = ReferenceType(self)
                else:
                    raise Exception(f'Unexpected reads_1.read value {reads_1.read!r} for CASAVA.')
            else:
                # Other file types go here ...
                self.reads_1 = reads_1
                reads_1.weak_reference_paired_reads = ReferenceType(self)

        if reads_2:
            if reads_2.file_type == 'CASAVA':
                if reads_2.read == 'R1':
                    self.reads_1 = reads_2
                    reads_2.weak_reference_paired_reads = ReferenceType(self)
                elif reads_2.read == 'R2':
                    self.reads_2 = reads_2
                    reads_2.weak_reference_paired_reads = ReferenceType(self)
                else:
                    raise Exception(f'Unexpected reads_2.read value {reads_2.read!r} for CASAVA.')
            else:
                # Other file types go here ...
                self.reads_2 = reads_2
                reads_2.weak_reference_paired_reads = ReferenceType(self)

        self.exclude = exclude
        self.index_1 = index_1
        self.index_2 = index_2
        self.read_group = read_group

        self.weak_reference_sample = weak_reference_sample

        return

    def __eq__(self, other):
        """Test :py:class:`bsf.ngs.PairedReads` objects for equality.

        :param other: A :py:class:`bsf.ngs.PairedReads` object.
        :type other: PairedReads
        :return: :py:const:`True` if equal, :py:const:`False` otherwise.
        :rtype: bool
        """
        assert isinstance(other, PairedReads)

        if self is other:
            return True

        return \
            super(PairedReads, self).__eq__(other=other) \
            and self.reads_1 == other.reads_1 \
            and self.reads_2 == other.reads_2 \
            and self.exclude == other.exclude \
            and self.index_1 == other.index_1 \
            and self.index_2 == other.index_2 \
            and self.read_group == other.read_group

    def __repr__(self):
        """Return a printable representation of a :py:class:`bsf.ngs.PairedReads` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'PairedReads(' \
            f'name={self.name!r} ' \
            f'file_path={self.file_path!r} ' \
            f'file_type={self.file_type!r} ' \
            f'annotation_dict={self.annotation_dict!r} ' \
            f'reads_1={self.reads_1!r} ' \
            f'reads_2={self.reads_2!r} ' \
            f'exclude={self.exclude!r} ' \
            f'index_1={self.index_1!r} ' \
            f'index_2={self.index_2!r} ' \
            f'read_group={self.read_group!r}' \
            f')'

    def trace(self, level):
        """Trace a :py:class:`bsf.ngs.PairedReads` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append(f'{indent}{self.__class__!r}\n')
        str_list.append(f'{indent}  weak_reference_sample: {self.weak_reference_sample!r}\n')
        str_list.append(f'{indent}  name:       {self.name!r}\n')
        str_list.append(f'{indent}  file_path:  {self.file_path!r}\n')
        str_list.append(f'{indent}  file_type:  {self.file_type!r}\n')
        str_list.append(f'{indent}  annotation_dict:\n')
        for annotation_key in sorted(self.annotation_dict):
            str_list.append(f'{indent}    {annotation_key!r} {self.annotation_dict[annotation_key]!r}\n')
        str_list.append(f'{indent}  reads_1:    {self.reads_1!r}\n')
        str_list.append(f'{indent}  reads_2:    {self.reads_2!r}\n')
        str_list.append(f'{indent}  index_1:    {self.index_1!r}\n')
        str_list.append(f'{indent}  index_2:    {self.index_2!r}\n')
        str_list.append(f'{indent}  exclude:    {self.exclude!r}\n')
        str_list.append(f'{indent}  read_group: {self.read_group!r}\n')

        if self.reads_1 is not None:
            str_list.extend(self.reads_1.trace(level=level + 1))
        if self.reads_2 is not None:
            str_list.extend(self.reads_2.trace(level=level + 1))

        return str_list

    def add_reads(self, reads):
        """Add a :py:class:`bsf.ngs.Reads` object.

        For a :py:attr:`bsf.ngs.Reads.file_type` :literal:`CASAVA` the :py:class:`bsf.ngs.Reads` object can be
        automatically assigned on the basis of the :py:attr:`bsf.ngs.Reads.read` instance variable
        (i.e., :literal:`R1` or :literal:`R2`).

        See also the :py:meth:`bsf.ngs.Sample.add_reads` method.

        :param reads: A :py:class:`bsf.ngs.Reads` object.
        :type reads: Reads
        :return: :py:const:`True` upon success, :py:const:`False` otherwise.
        :rtype: bool
        """
        assert isinstance(reads, Reads)

        if self.reads_1 is not None:
            if not self.reads_1.match_paired(reads=reads):
                # logger.warning(
                #     f'The Reads object is not compatible with the PairedReads.reads_1 object. '
                #     f'reads: {reads!r} '
                #     f'paired_reads.reads_1: {self.reads_1!r}'
                # )
                return False

            if self.reads_1.file_type == 'CASAVA':
                if reads.read == 'R1':
                    raise Exception(
                        'PairedReads.reads_1 has already been defined.\n'
                        f'  reads_1: {self.reads_1.file_path!r}\n'
                        f'  reads:   {reads.file_path!r}'
                    )
                elif reads.read == 'R2':
                    self.reads_2 = reads
                    reads.weak_reference_paired_reads = ReferenceType(self)
                    return True
                else:
                    raise Exception(f'Unexpected Reads.read value {reads.read!r} for CASAVA.')

        if self.reads_2 is not None:
            if not self.reads_2.match_paired(reads=reads):
                # logger.warning(
                #     f'The Reads object is not compatible with the PairedReads.reads_2 object. '
                #     f'reads: {reads!r} '
                #     f'paired_reads.reads_2: {self.reads_2!r}'
                # )
                return False

            if self.reads_2.file_type == 'CASAVA':
                if reads.read == 'R1':
                    self.reads_1 = reads
                    reads.weak_reference_paired_reads = ReferenceType(self)
                    return True
                elif reads.read == 'R2':
                    raise Exception(
                        'PairedReads.reads_2 has already been defined.\n'
                        f'  reads_2: {self.reads_2.file_path!r}\n'
                        f'  reads:   {reads.file_path!r}'
                    )
                else:
                    raise Exception(f'Unexpected Reads.read value {reads.read!r} for CASAVA.')

        return False

    def match(self, paired_reads):
        """Match :py:class:`bsf.ngs.PairedReads` objects.

        :param paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
        :type paired_reads: PairedReads
        :return: :py:const:`True` if both objects match, :py:const:`False` otherwise.
        :rtype: bool
        """
        if self is paired_reads:
            return True

        if self.index_1 != paired_reads.index_1:
            return False

        if self.index_2 != paired_reads.index_2:
            return False

        if self.read_group != paired_reads.read_group:
            return False

        if not self.reads_1.match(reads=paired_reads.reads_1):
            return False

        if (self.reads_2 is not None) and (paired_reads.reads_2 is not None):
            if not self.reads_2.match(reads=paired_reads.reads_2):
                return False

        return True

    def get_name(self, full=False):
        """Get the name of a :py:class:`bsf.ngs.PairedReads` object.

        For the :py:attr:`bsf.ngs.Reads.file_type` :literal:`CASAVA` the name is a concatenation of the
        :py:attr:`bsf.ngs.Reads.name`, :py:attr:`bsf.ngs.Reads.barcode` and
        :py:attr:`bsf.ngs.Reads.lane` instance variables, preferentially derived from the first
        :py:class:`bsf.ngs.Reads` object in the :py:class:`bsf.ngs.PairedReads` object.

        If the :literal:`full` parameter is set, :py:attr:`bsf.ngs.Reads.read` and
        :py:attr:`bsf.ngs.Reads.chunk` are also added.

        :param full: Return the full name including read and chunk information.
        :type full: bool
        :return: A name.
        :rtype: str
        """
        if self.reads_1 is not None:
            reads = self.reads_1
        elif self.reads_2 is not None:
            reads = self.reads_2
        else:
            return

        if reads.file_type == 'CASAVA':
            if full:
                name = '_'.join((reads.name, reads.barcode, reads.lane, reads.read, reads.chunk))
            else:
                name = '_'.join((reads.name, reads.barcode, reads.lane))
        else:
            name = reads.name

        return name


class Sample(NextGenerationBase):
    """The :py:class:`bsf.ngs.Sample` class represents a Next-Generation Sequencing sample.

    It consists of one or more :py:class:`bsf.ngs.PairedReads` objects as (biological or technical) replicates
    that result from the same flow cell.

    :cvar default_name: A default name.
    :type default_name: str
    :ivar paired_reads_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.PairedReads` objects.
    :type paired_reads_list: list[PairedReads]
    :ivar weak_reference_project: A :py:class:`weakref.ReferenceType` object pointing at a
        :py:class:`bsf.ngs.Project` object.
    :type weak_reference_project: ReferenceType | None
    """

    default_name = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a :py:class:`bsf.ngs.Sample` object from a file path.

        For a :literal:`file_type` :literal:`CASAVA` the name is automatically populated,
        while :py:class:`bsf.ngs.PairedReads` objects are automatically discovered.

        See also the :py:meth:`bsf.ngs.Project.from_file_path` method.

        :param file_path: A file path.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str
        :return: A :py:class:`bsf.ngs.Sample` object.
        :rtype: Sample
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_type == 'CASAVA':
            # CASAVA Samples obey a "Sample_name" schema.

            component_list = file_name.split('_')

            sample = cls(
                file_path=file_path,
                file_type=file_type,
                name='_'.join(component_list[1:]))

            # Automatically discover CASAVA Reads files ...

            for file_name in os.listdir(sample.file_path):
                if not file_name.endswith('fastq.gz'):
                    continue
                file_path = os.path.join(sample.file_path, file_name)
                if not os.path.isfile(file_path):
                    continue
                sample.add_reads(reads=Reads.from_file_path(file_path=file_path, file_type=file_type))
        else:
            raise Exception(f'Unsupported file_type: {file_type!r}')

        return sample

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            paired_reads_list=None,
            weak_reference_project=None):
        """Initialise a :py:class:`bsf.ngs.Sample` object.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        :param paired_reads_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.PairedReads` objects.
        :type paired_reads_list: list[PairedReads] | None
        :param weak_reference_project: A :py:class:`weakref.ReferenceType` object pointing at a
            :py:class:`bsf.ngs.Project` object.
        :type weak_reference_project: ReferenceType | None
        """
        super(Sample, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        if paired_reads_list is None:
            self.paired_reads_list = list()
        else:
            self.paired_reads_list = paired_reads_list
            # Set this bsf.ngs.Sample as a weak reference for each bsf.ngs.PairedReads object on the list.
            for paired_reads_object in paired_reads_list:
                paired_reads_object.weak_reference_sample = ReferenceType(self)

        self.weak_reference_project = weak_reference_project

        return

    def __repr__(self):
        """Return a printable representation of a :py:class:`bsf.ngs.Sample` object.

        :return: A printable representation.
        :rtype: str
        """
        return \
            f'Sample(' \
            f'name={self.name!r} ' \
            f'file_path={self.file_path!r} ' \
            f'file_type={self.file_type!r} ' \
            f'annotation_dict={self.annotation_dict!r} ' \
            f'paired_reads_list={self.paired_reads_list!r}' \
            f')'

    def trace(self, level):
        """Trace a :py:class:`bsf.ngs.Sample` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append(f'{indent}{self.__class__!r}\n')
        str_list.append(f'{indent}  weak_reference_project: {self.weak_reference_project!r}\n')
        str_list.append(f'{indent}  name:      {self.name!r}\n')
        str_list.append(f'{indent}  file_path: {self.file_path!r}\n')
        str_list.append(f'{indent}  file_type: {self.file_type!r}\n')
        str_list.append(f'{indent}  annotation_dict:\n')
        for annotation_key in sorted(self.annotation_dict):
            str_list.append(f'{indent}    {annotation_key!r} {self.annotation_dict[annotation_key]!r}\n')
        str_list.append(f'{indent}  paired_reads:\n')
        for paired_reads in self.paired_reads_list:
            str_list.extend(paired_reads.trace(level=level + 1))

        return str_list

    def match(self, sample):
        """Match :py:class:`bsf.ngs.Sample` objects.

        :param sample: A :py:class:`bsf.ngs.Sample` object.
        :type sample: Sample
        :return: :py:const:`True` if both objects match, :py:const:`False` otherwise.
        :rtype: bool
        """
        # When Sample objects get merged, the file_path does no longer make sense.
        # Hence, cope with empty file paths here ...

        if self is sample:
            return True

        if (self.file_path and sample.file_path) and (self.file_path != sample.file_path):
            return False

        if self.file_type != sample.file_type:
            return False

        if self.name != sample.name:
            return False

        # Match Python list objects of PairedReads objects ...

        for paired_reads_1 in self.paired_reads_list:
            match = False
            for paired_reads_2 in sample.paired_reads_list:
                if paired_reads_1.match(paired_reads=paired_reads_2):
                    match = True
                    break

            if not match:
                # No match for this particular PairedReads object.
                return False

        return True

    def add_paired_reads(self, paired_reads):
        """Add a :py:class:`bsf.ngs.PairedReads` object.

        This method checks, whether a matching :py:class:`bsf.ngs.PairedReads` object is already present in this
        :py:class:`bsf.ngs.Sample` object.

        :param paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
        :type paired_reads: PairedReads
        """
        assert isinstance(paired_reads, PairedReads)

        # Iterate through the Python list of PairedReads objects.
        # The PairedReads object must not match.

        for old_paired_reads in self.paired_reads_list:
            if old_paired_reads.match(paired_reads=paired_reads):
                break
        else:
            # None of the existing PairedReads objects has matched,
            # so add this one to the Sample object.
            self.paired_reads_list.append(paired_reads)
            paired_reads.weak_reference_sample = ReferenceType(self)

        return

    def add_reads(self, reads):
        """Add a :py:class:`bsf.ngs.Reads` object.

        See also the :py:meth:`bsf.ngs.Sample.from_file_path` method for CASAVA.

        :param reads: A :py:class:`bsf.ngs.Reads` object.
        :type reads: Reads
        """
        assert isinstance(reads, Reads)

        # Iterate through the Python list of PairedReads objects
        # The sample name must be identical.
        # The chunk must be identical
        # The read must fit

        for paired_reads in self.paired_reads_list:
            if paired_reads.add_reads(reads=reads):
                break
        else:
            # The Reads object could not be added.
            # Create a new PairedReads object initialised
            # with the Reads object and add it to the Sample.
            self.add_paired_reads(paired_reads=PairedReads(reads_1=reads))

        return

    def get_all_paired_reads(self, replicate_grouping, exclude=False, full=False):
        """Get all :py:class:`bsf.ngs.PairedReads` objects of a :py:class:`bsf.ngs.Sample` grouped or un-grouped.

        A :py:class:`bsf.ngs.Sample` object can hold several :py:class:`bsf.ngs.PairedReads` objects
        (i.e., SAM Read Group (:literal:`@RG`) entries) that have been sequenced on different lanes of the same or
        even another flow cell.

        The :py:class:`bsf.ngs.PairedReads` objects will therefore differ in
        :py:attr:`bsf.ngs.Reads.name`,
        :py:attr:`bsf.ngs.Reads.barcode` or
        :py:attr:`bsf.ngs.Reads.lane` information.
        Depending on the :literal:`replicate_grouping` parameter,
        :py:class:`bsf.ngs.PairedReads` objects can be returned as a group or separately.
        However, :py:class:`bsf.ngs.PairedReads` objects that share the
        :py:attr:`bsf.ngs.Reads.name`,
        :py:attr:`bsf.ngs.Reads.barcode` and
        :py:attr:`bsf.ngs.Reads.lane`,
        but only differ in their :py:attr:`bsf.ngs.Reads.chunk` number are always grouped together.

        :param replicate_grouping: Group all :py:class:`bsf.ngs.PairedReads` objects of a
            :py:class:`bsf.ngs.Sample` object or list them individually.
        :type replicate_grouping: bool
        :param exclude: Exclude on the basis of the :py:attr:`bsf.ngs.PairedReads.exclude` attribute.
        :type exclude: bool
        :param full: Return the full name including read and chunk information.
        :type full: bool
        :return: A Python :py:class:`dict` object of Python :py:class:`str` (sensible replicate name) key and
            Python :py:class:`list` object of :py:class:`bsf.ngs.PairedReads` objects value data.
        :rtype: dict[str, list[PairedReads]]
        """
        paired_reads_dict: Dict[str, List[PairedReads]] = dict()

        for paired_reads in self.paired_reads_list:
            if exclude and paired_reads.exclude:
                # Skip excluded PairedReads.
                continue

            if replicate_grouping:
                # If grouped, use the Sample.name instance variable as key so that all PairedReads
                # of this Sample end up on the same Python list.
                key = self.name
            else:
                # If un-grouped use the PairedReads.get_name() method as key so that each
                # PairedReads of this Sample ends up on a separate Python list.
                key = paired_reads.get_name(full=full)

            if not key:
                continue

            if key not in paired_reads_dict:
                paired_reads_dict[key] = list()

            paired_reads_dict[key].append(paired_reads)

        return paired_reads_dict

    def is_excluded(self):
        """Is this :py:class:`bsf.ngs.Sample` object excluded,
        because all its :py:class:`bsf.ngs.PairedReads` objects are?

        :return: :py:const:`True` if the :py:class:`bsf.ngs.Sample` object is excluded since all
            :py:class:`bsf.ngs.PairedReads` objects are excluded.
        :rtype: bool
        """
        for paired_read in self.paired_reads_list:
            if not paired_read.exclude:
                return False
        else:
            return True


class Project(NextGenerationBase):
    """The :py:class:`bsf.ngs.Project` class represents a Next-Generation Sequencing Project
    consisting of one or more :py:class:`bsf.ngs.Sample` objects.

    :cvar default_name: A default name.
    :type default_name: str
    :ivar sample_dict: A Python :py:class:`dict` object of
        :py:attr:`bsf.ngs.Sample.name` key objects and
        :py:class:`bsf.ngs.Sample` value objects.
    :type sample_dict: dict[str, Sample]
    :ivar weak_reference_prf: A :py:class:`weakref.ReferenceType` object pointing at a
        :py:class:`bsf.ngs.ProcessedRunFolder` object.
    :type weak_reference_prf: ReferenceType | None
    """

    default_name = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a :py:class:`bsf.ngs.Project` object from a file path.

        See also the :py:meth:`bsf.ngs.ProcessedRunFolder.from_file_path` method for CASAVA.

        :param file_path: A file path.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str
        :return: A :py:class:`bsf.ngs.Project` object.
        :rtype: Project
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_type == 'CASAVA':
            # CASAVA Projects obey a "Project_name" schema.

            component_list = file_name.split('_')

            project = cls(
                file_path=file_path,
                file_type=file_type,
                name='_'.join(component_list[1:]))

            # Automatically discover CASAVA Sample directories ...

            re_pattern = re.compile(pattern=r'^Sample_(.*)$')
            for file_name in os.listdir(project.file_path):
                file_path = os.path.join(project.file_path, file_name)
                re_match = re_pattern.search(string=file_name)
                if os.path.isdir(file_path) and re_match is not None:
                    if re_match.group(1) in project.sample_dict:
                        raise Exception(f'Sample with name {re_match.group(1)!r} already exists.')
                    else:
                        project.add_sample(sample=Sample.from_file_path(file_path=file_path, file_type=file_type))
        else:
            raise Exception(f'Unsupported file_type: {file_type!r}')

        return project

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            sample_dict=None,
            weak_reference_prf=None):
        """Initialise a :py:class:`bsf.ngs.Project` object.

        For a :literal:`file_type` :literal:`CASAVA` the name is automatically populated,
        while :py:class:`bsf.ngs.Sample` objects are automatically discovered.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        :param sample_dict: A Python :py:class:`dict` object of
            :py:attr:`bsf.ngs.Sample.name` key objects and
            :py:class:`bsf.ngs.Sample´ value objects.
        :type sample_dict: dict[str, Sample] | None
        :param weak_reference_prf: A :py:class:`weakref.ReferenceType` object pointing at a
            :py:class:`bsf.ngs.ProcessedRunFolder` object.
        :type weak_reference_prf: ReferenceType | None
        :raise Exception: If :py:attr:`bsf.ngs.Sample.name` values are not unique for
            :literal:`file_type` :literal:`CASAVA`
        """
        super(Project, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        if sample_dict is None:
            self.sample_dict = dict()
        else:
            self.sample_dict = sample_dict
            # Set this bsf.ngs.Project as weak reference for each bsf.ngs.Sample in the dict.
            for sample in self.sample_dict.values():
                sample.weak_reference_project = ReferenceType(self)

        self.weak_reference_prf = weak_reference_prf

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.ngs.Project` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append(f'{indent}{self.__class__!r}\n')
        str_list.append(f'{indent}  weak_reference_prf: {self.weak_reference_prf!r}\n')
        str_list.append(f'{indent}  file_path: {self.file_path!r}\n')
        str_list.append(f'{indent}  file_type: {self.file_type!r}\n')
        str_list.append(f'{indent}  name:      {self.name!r}\n')
        str_list.append(f'{indent}  annotation_dict:\n')
        for annotation_key in sorted(self.annotation_dict):
            str_list.append(f'{indent}    {annotation_key!r} {self.annotation_dict[annotation_key]!r}\n')
        str_list.append(f'{indent}  sample_dict:\n')
        for sample_name in sorted(self.sample_dict):
            str_list.extend(self.sample_dict[sample_name].trace(level=level + 1))

        return str_list

    def add_sample(self, sample):
        """Add a :py:class:`bsf.ngs.Sample` object to a :py:class:`bsf.ngs.Project` object and set
        a weak reference in to :py:class:`bsf.ngs.Sample` object back to the :py:class:`bsf.ngs.Project` object.

        :param sample: A :py:class:`bsf.ngs.Sample` object.
        :type sample: Sample
        :return: A :py:class:`bsf.ngs.Sample` object.
        :rtype: Sample
        """
        assert isinstance(sample, Sample)
        # Delete an eventual bsf.ngs.Sample stored under the same bsf.ngs.Sample.name to remove the weak reference.
        self.del_sample(name=sample.name)
        self.sample_dict[sample.name] = sample
        sample.weak_reference_project = ReferenceType(self)

        return sample

    def del_sample(self, name):
        """Delete a :py:class:`bsf.ngs.Sample` object from a :py:class:`bsf.ngs.Project` object and
        clear the weak reference, if it points back at the :py:class:`bsf.ngs.Project` object.

        :param name: A :py:class:`bsf.ngs.Sample.name` object.
        :type name: str
        :return: A :py:class:`bsf.ngs.Sample` object.
        :rtype: Sample
        """
        if name in self.sample_dict:
            sample = self.sample_dict[name]
            del self.sample_dict[name]
            if (sample.weak_reference_project is not None) and (sample.weak_reference_project() is self):
                sample.weak_reference_project = None
            return sample
        else:
            return

    def get_all_samples(self, exclude=False):
        """Get an ordered Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.

        :param exclude: Exclude :py:class:`bsf.ngs.Sample` objects on the basis of the
            :py:attr:`bsf.ngs.PairedReads.exclude` attribute.
        :type exclude: bool
        :return: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :rtype: list[Sample]
        """
        sample_list: List[Sample] = list()

        for sample_name in sorted(self.sample_dict):
            sample = self.sample_dict[sample_name]
            if not (exclude and sample.is_excluded()):
                sample_list.append(sample)

        return sample_list


class ProcessedRunFolder(NextGenerationBase):
    """The :py:class:`bsf.ngs.ProcessedRunFolder` class represents an Illumina Run Folder after processing with CASAVA.

    :cvar default_name: A default name.
    :type default_name: str
    :ivar prefix: A prefix.
    :type prefix: str | None
    :ivar flow_cell: A flow cell identifier.
    :type flow_cell: str | None
    :ivar version: A version number.
    :type version: str | None
    :ivar project_dict: A Python :py:class:`dict` object of
        :py:attr:`bsf.ngs.Project.name` key objects and
        :py:class:`bsf.ngs.Project` value objects.
    :type project_dict: dict[str, Project]
    :ivar weak_reference_collection: A :py:class:`weakref.ReferenceType` object pointing at a
        :py:class:`bsf.ngs.Collection` object.
    :type weak_reference_collection: ReferenceType | None
    """

    default_name = 'Default'

    @staticmethod
    def guess_file_type(file_path):
        """Guess the :literal:`file_type` of a :py:class:`bsf.ngs.ProcessedRunFolder` object on the basis of the
        :literal:`file_path`.

        CASAVA :py:class:`bsf.ngs.ProcessedRunFolder` objects obey a :literal:`Prefix_FCID_CASAVA182`
        schema.

        The following prefixes are currently in use:
            - :literal:`BSF_` Biomedical Sequencing Facility
            - :literal:`NGS_` Kaan Boztug group
            - :literal:`MUW_` Medical University Vienna
            - :literal:`SET_` Robert Kralovics group

        :param file_path: A file path.
        :type file_path: str
        :return: A file type (i.e., :literal:`CASAVA` or :literal:`External`).
        :rtype: str
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if isinstance(file_name, bytes):
            file_name = file_name.decode()

        component_list = file_name.split('_')

        if component_list[-1].startswith('CASAVA'):
            return 'CASAVA'
        else:
            return 'External'

    @classmethod
    def from_file_path(cls, file_path, file_type=None):
        """Construct a :py:class:`bsf.ngs.ProcessedRunFolder` object from a file path.

        For the :literal:`file_type` :literal:`CASAVA`, the :py:attr:`bsf.ngs.ProcessedRunFolder.name`,
        :py:attr:`bsf.ngs.ProcessedRunFolder.prefix`,
        :py:attr:`bsf.ngs.ProcessedRunFolder.flow_cell` and :py:attr:`bsf.ngs.ProcessedRunFolder.version`
        instance variables can be automatically parsed from the :literal:`file_path`, while
        :py:class:`bsf.ngs.Project` objects can be automatically discovered.

        :param file_path: A file path.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :return: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
        :rtype: ProcessedRunFolder
        """
        # Try to determine the file_type if not explicitly specified.

        if not file_type or file_type == 'Automatic':
            file_type = ProcessedRunFolder.guess_file_type(file_path=file_path)

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_type == 'CASAVA':
            # CASAVA Processed Run Folders obey a "Prefix_FCID_CASAVA182"
            # schema. The following prefixes are currently in use:
            # -- BSF_ Biomedical Sequencing Facility
            # -- NGS_ Kaan Boztug group
            # -- MUW_ Medical University Vienna
            # -- SET_ Robert Kralovics group

            if isinstance(file_name, bytes):
                file_name = file_name.decode()

            component_list = file_name.split('_')

            prf = cls(
                file_path=file_path,
                file_type=file_type,
                name=file_name,
                prefix=component_list[0],
                flow_cell=component_list[1],
                version=component_list[2])

            # Automatically discover CASAVA Project directories ...

            re_pattern = re.compile(pattern=r'^Project_(.*)$')
            for file_name in os.listdir(prf.file_path):
                file_path = os.path.join(prf.file_path, file_name)
                re_match = re_pattern.search(string=file_name)
                if os.path.isdir(file_path) and re_match is not None:
                    if re_match.group(1) in prf.project_dict:
                        raise Exception(f'Project with name {re_match.group(1)!r} already exists.')
                    else:
                        prf.add_project(project=Project.from_file_path(file_path=file_path, file_type=file_type))
        elif file_type == 'External':
            # Create a new, minimal ProcessedRunFolder.
            prf = cls(file_path=file_path, file_type=file_type, name=file_name)
        else:
            raise Exception(f'Unsupported file_type: {file_type!r}')

        return prf

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            prefix=None,
            flow_cell=None,
            version=None,
            project_dict=None,
            weak_reference_collection=None):
        """Initialise a :py:class:`bsf.ngs.ProcessedRunFolder` object.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        :param prefix: A prefix.
        :type prefix: str | None
        :param flow_cell: A flow cell identifier.
        :type flow_cell: str | None
        :param version: A version.
        :type version: str | None
        :param project_dict: A Python :py:class:`dict` object of
            :py:attr:`bsf.ngs.Project.name` key objects and
            :py:class:`bsf.ngs.Project` value objects.
        :type project_dict: dict[str, Project] | None
        :param weak_reference_collection: A :py:class:`weakref.ReferenceType` object pointing at a
            :py:class:`bsf.ngs.Collection` object.
        :type weak_reference_collection: ReferenceType | None
        :raise Exception: If :py:attr:`bsf.ngs.Project.name` attribute values are not unique for
            file_type :literal:`CASAVA`
        """
        super(ProcessedRunFolder, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        self.prefix = prefix
        self.flow_cell = flow_cell
        self.version = version

        if project_dict is None:
            self.project_dict = dict()
        else:
            self.project_dict = project_dict
            # Set this bsf.ngs.ProcessedRunFolder as weak reference for each bsf.ngs.Project in the dict.
            for project in self.project_dict.values():
                project.weak_reference_prf = ReferenceType(self)

        self.weak_reference_collection = weak_reference_collection

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.ngs.ProcessedRunFolder` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append(f'{indent}{self.__class__!r}\n')
        str_list.append(f'{indent}  weak_reference_collection: {self.weak_reference_collection!r}\n')
        str_list.append(f'{indent}  file_path: {self.file_path!r}\n')
        str_list.append(f'{indent}  file_type: {self.file_type!r}\n')
        str_list.append(f'{indent}  name:      {self.name!r}\n')
        str_list.append(f'{indent}  prefix:    {self.prefix!r}\n')
        str_list.append(f'{indent}  flow_cell: {self.flow_cell!r}\n')
        str_list.append(f'{indent}  version:   {self.version!r}\n')
        str_list.append(f'{indent}  annotation_dict:\n')
        for annotation_key in sorted(self.annotation_dict):
            str_list.append(f'{indent}    {annotation_key!r} {self.annotation_dict[annotation_key]!r}\n')
        str_list.append(f'{indent}  project_dict:\n')
        for project_name in sorted(self.project_dict):
            str_list.extend(self.project_dict[project_name].trace(level=level + 1))

        return str_list

    def add_project(self, project):
        """Add a :py:class:`bsf.ngs.Project` object to the :py:class:`bsf.ngs.ProcessedRunFolder` object and set
        a weak reference in the :py:class:`bsf.ngs.Project` object to the
        :py:class:`bsf.ngs.ProcessedRunFolder` object.

        :param project: A :py:class:`bsf.ngs.Project` object.
        :type project: Project
        :return: A :py:class:`bsf.ngs.Project` object.
        :rtype: Project
        """
        assert isinstance(project, Project)
        # Delete an eventual bsf.ngs.Project that is stored under the same bsf.ngs.Project.name.
        self.del_project(name=project.name)
        self.project_dict[project.name] = project
        project.weak_reference_prf = ReferenceType(self)

        return project

    def del_project(self, name):
        """Delete a :py:class:`bsf.ngs.Project` object from am :py:class:`bsf.ngs.ProcessedRunFolder` object and
        clear the weak reference, if it points back at the :py:class:`bsf.ngs.ProcessedRunFolder` object.

        :param name: A :py:attr:`bsf.ngs.Project.name` attribute.
        :type name: str
        :return: A :py:class:`bsf.ngs.Project` object.
        :rtype: Project
        """
        if name in self.project_dict:
            project = self.project_dict[name]
            del self.project_dict[name]
            if (project.weak_reference_prf is not None) and (project.weak_reference_prf() is self):
                project.weak_reference_prf = None
            return project
        else:
            return

    def get_all_projects(self):
        """Get an ordered Python :py:class:`list` object of :py:class:`bsf.ngs.Project` objects.

        :return: A Python :py:class:`list` object of :py:class:`bsf.ngs.Project` objects.
        :rtype: list[Project]
        """
        project_list: List[Project] = list()

        for project_name in sorted(self.project_dict):
            project_list.append(self.project_dict[project_name])

        return project_list


class Collection(NextGenerationBase):
    """The :py:class:`bsf.ngs.Collection` class represents a collection of
    one or more :py:class:`bsf.ngs.ProcessedRunFolder` objects.

    :cvar default_name: A default name.
    :type default_name: str
    :ivar processed_run_folder_dict: A Python :py:class:`dict` object of
        :py:attr:`bsf.ngs.ProcessedRunFolder.name` key objects and
        :py:class:`bsf.ngs.ProcessedRunFolder` value objects.
    :type processed_run_folder_dict: dict[str, ProcessedRunFolder]
    :ivar sample_group_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (group name) key objects and
            second-level Python :py:class:`dict` value objects of
            :py:attr:`bsf.ngs.Sample.name` key objects and
            :py:class:`bsf.ngs.Sample` value objects.
    :type sample_group_dict: dict[str, list[Sample]]
    """

    default_name = 'Default'

    # Python dict to map from an int (status value) to a (private) function.
    _reads_function_dict = None
    _paired_reads_function_dict = None

    @classmethod
    def from_sas_path(cls, file_path, file_type, name, sas_path, sas_prefix=None):
        """Construct a :py:class:`bsf.ngs.Collection` from a
        :py:class:`bsf.ngs.SampleAnnotationSheet` object file path.

        :param file_path: A file path as prefix for relative file paths in the
            :py:class:`bsf.ngs.SampleAnnotationSheet` object.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str
        :param name: A name.
        :type name: str
        :param sas_path: A :py:class:`bsf.ngs.SampleAnnotationSheet` object file path.
        :type sas_path: str
        :param sas_prefix: A prefix to columns in a :py:class:`bsf.ngs.SampleAnnotationSheet` object
            (e.g., :literal:`[Control] Sample`, :literal:`[Treatment] Sample`, ...).
        :type sas_prefix: str | None
        :return: A :py:class:`bsf.ngs.Collection` object.
        :rtype: Collection
        """
        return cls.from_sas(
            file_path=file_path,
            file_type=file_type,
            name=name,
            sas=SampleAnnotationSheet.from_file_path(file_path=sas_path),
            sas_prefix=sas_prefix)

    @classmethod
    def from_sas(cls, file_path, file_type, name, sas, sas_prefix=None):
        """Construct a :py:class:`bsf.ngs.Collection` object from a :py:class:`bsf.ngs.SampleAnnotationSheet` object.

        This method creates :py:class:`bsf.ngs.Reads`, :py:class:`bsf.ngs.PairedReads`, :py:class:`bsf.ngs.Sample`,
        :py:class:`bsf.ngs.Project` and :py:class:`bsf.ngs.ProcessedRunFolder` objects from a
        :py:class:`bsf.ngs.SampleAnnotationSheet` row Python :py:class:`dict` object.

        If object-specific keys or their corresponding values are not
        available from the Python row dict, new objects will be created
        with a default key.

        Sample Annotation Sheet format:
            - FileType: Data object :literal:`file_type` (i.e., :literal:`CASAVA` or :literal:`External`),
                defaults to :literal:`External`.
                The FileType :literal:`CASAVA` allows for auto-discovery of
                :py:class:`bsf.ngs.ProcessedRunFolder` objects.
            - ProcessedRunFolder Name: A :py:attr:`bsf.ngs.ProcessedRunFolder.file_path`,
                can be automatically registered.
            - Project Name: A :py:attr:`bsf.ngs.Project.name` attribute.
            - Sample Name: A :py:attr:`bsf.ngs.Sample.name` attribute.
            - Reads1 File: A :py:attr:`bsf.ngs.Reads.file_name` attribute.
                Subjected to :py:func:`os.path.expanduser` (i.e., on UNIX ~ or ~user) and
                :py:func:`os.path.expandvars` (i.e., on UNIX ${NAME} or $NAME).
                If still relative, the :py:attr:`bsf.ngs.Collection.file_path` is prepended.
            - Reads1 Name: A :py:attr:`bsf.ngs.Reads.name` attribute
            - Reads2 File: Same as Reads1 File
            - Reads2 Name: Same as Reads1 Name
            - Group: :py:class:`bsf.ngs.Sample` objects can be grouped for further analysis
                (e.g., RNA-Seq or ChIP-Seq experiments).

        :param file_path: A file path as prefix for relative file paths in the
            :py:class:`bsf.ngs.SampleAnnotationSheet` object.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str
        :param name: A name.
        :type name: str
        :param sas: A :py:class:`bsf.ngs.SampleAnnotationSheet` object.
        :type sas: SampleAnnotationSheet
        :param sas_prefix: A prefix to columns in a :py:class:`bsf.ngs.SampleAnnotationSheet` object
            (e.g., :literal:`[Control] Sample`, :literal:`[Treatment] Sample`, ...).
        :type sas_prefix: str | None
        :return: A :py:class:`bsf.ngs.Collection` object.
        :rtype: Collection
        """
        current_prf: Optional[ProcessedRunFolder] = None
        current_project: Optional[Project] = None
        current_sample: Optional[Sample] = None
        current_paired_reads: Optional[PairedReads] = None

        def process_file_type():
            """Private function to get file type information.

            A :literal:`[Prefix] FileType` key is optional, its value defaults to :literal:`Automatic`.

            :return: A file type.
            :rtype: str
            """
            _key = 'File Type'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                # The key exists ...
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    # ... and has a meaningful value ...
                    file_type_new = _value
                else:
                    # ... and has no meaningful value ...
                    file_type_new = 'Automatic'
            else:
                # The key does not exists ...
                file_type_new = 'Automatic'

            return file_type_new

        def process_processed_run_folder(collection, prf):
            """Private function to get or create a :py:class:`bsf.ngs.ProcessedRunFolder` object.

            A :literal:`[Prefix] ProcessedRunFolder Name` key is optional, its value defaults to :literal:`Default`.

            :param collection: A :py:class:`bsf.ngs.Collection` object.
            :type collection: Collection
            :param prf: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
            :type prf: ProcessedRunFolder | None
            :return: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
            :rtype: ProcessedRunFolder
            """
            _key = 'ProcessedRunFolder Name'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                # The key exists ...
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    # ... and has a meaningful value ...
                    if _value in collection.processed_run_folder_dict:
                        # ..., which exists in the dict of ProcessedRunFolder objects.
                        prf = collection.processed_run_folder_dict[_value]
                    else:
                        # ..., which does not exist in the dict of ProcessedRunFolder objects.
                        # Try to automatically discover a ProcessedRunFolder.
                        prf = collection.get_processed_run_folder(file_path=_value, file_type=file_type)

            if prf is None:
                if ProcessedRunFolder.default_name in collection.processed_run_folder_dict:
                    prf = collection.processed_run_folder_dict[ProcessedRunFolder.default_name]
                else:
                    prf = collection.add_processed_run_folder(
                        prf=ProcessedRunFolder(name=ProcessedRunFolder.default_name, file_type=file_type))

            prf.process_annotation(row_dict=row_dict, key_list=key_list, prefix=sas_prefix)

            return prf

        def process_project(prf, project):
            """Private function to get or create a :py:class:`bsf.ngs.Project` object.

            A :literal:`[Prefix] Project Name` key is optional, its value defaults to :literal:`Default`.

            :param prf: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
            :type prf: ProcessedRunFolder
            :param project: A :py:class:`bsf.ngs.Project` object.
            :type project: Project
            :return: A :py:class:`bsf.ngs.Project` object.
            :rtype: Project
            """
            _key = 'Project Name'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                # The key exists ...
                key_list.remove(_key)  # Remove the 'Name' key.
                _value = row_dict[_key]
                if _value:
                    # ... and has a meaningful value ...
                    if _value in prf.project_dict:
                        # ..., which exists in the dict of Project objects.
                        project = prf.project_dict[_value]
                    else:
                        # ..., which does not exist in the dict of Project objects.
                        # Create a new Project.
                        project = prf.add_project(project=Project(name=_value, file_type=file_type))

            if project is None:
                if Project.default_name in prf.project_dict:
                    project = prf.project_dict[Project.default_name]
                else:
                    project = prf.add_project(project=Project(name=Project.default_name, file_type=file_type))

            project.process_annotation(row_dict=row_dict, key_list=key_list, prefix=sas_prefix)

            return project

        def process_sample(project, sample):
            """Private function to get or create a :py:class:`bsf.ngs.Sample` object.

            A :literal:`[Prefix] Sample Name` key is optional, its value defaults to :literal:`Default`.

            :param project: A :py:class:`bsf.ngs.Project` object.
            :type project: Project
            :param sample: A :py:class:`bsf.ngs.Sample` object.
            :type sample: Sample | None
            :return: A :py:class:`bsf.ngs.Sample` object.
            :rtype: Sample
            """
            module_logger.debug(
                f'process_sample('
                f'project={project!r}, '
                f'sample={sample!r}'
                f')'
            )

            _key = 'Sample Name'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                # The key exists ...
                key_list.remove(_key)  # Remove the 'Name' key.
                _value = row_dict[_key]
                if _value:
                    # ... and has a meaningful value ...
                    if _value in project.sample_dict:
                        # ..., which exists in the dict of Sample objects.
                        sample = project.sample_dict[_value]
                    else:
                        # ..., which does not exist in the dict of Sample objects.
                        # Create a new Sample.
                        sample = project.add_sample(sample=Sample(name=_value, file_type=file_type))

            # If no Sample is defined create a default one.
            if sample is None:
                if Sample.default_name in project.sample_dict:
                    sample = project.sample_dict[Sample.default_name]
                else:
                    sample = project.add_sample(sample=Sample(name=Sample.default_name, file_type=file_type))

            sample.process_annotation(row_dict=row_dict, key_list=key_list, prefix=sas_prefix)

            module_logger.debug(
                f'process_sample final: '
                f'sample={sample!r}'
            )

            return sample

        def process_reads(reads, suffix, default_path):
            """Get or create a first or second :py:class:`bsf.ngs.Reads` object.

            A :literal:`[Prefix] Reads{suffix} Name` key and :literal:`[Prefix] Reads{suffix} File` key are optional,
            in which case the default is a :literal:`None` object.

            :param reads: A current :py:class:`bsf.ngs.Reads` object that may get replaced upon encountering a new
                :literal:`[Prefix] ReadsN Name` key.
            :type reads: Reads | None
            :param suffix: A read suffix (i.e., :literal:`1` or :literal:`2`).
            :type suffix: str
            :param default_path: A default file path.
            :type default_path: str
            :return: A :py:class:`bsf.ngs.Reads` object.
            :rtype: Reads | None
            """

            def process_reads_status(new, old):
                """Calculate a comparison status.
                    - 1 << 0 (1): new defined
                    - 1 << 1 (2): old defined
                    - 1 << 2 (4): new == old

                :param new: A new :py:class:`str` object.
                :type new: str | None
                :param old: An old :py:class:`str` object.
                :type old: str | None
                :return: A comparison status.
                :rtype: int
                """
                status = 0

                if new:
                    status |= 1 << 0  # 1

                if old:
                    status |= 1 << 1  # 2

                if (status == 3) and (new == old):
                    # Both arguments are defined (3) and have identical values.
                    status |= 1 << 2  # 4

                return status

            def process_reads_new(_reads, _reads_file, _reads_name, _file_type):
                """Private function to process (i.e., initialise) a new :py:class:`bsf.ngs.Reads` object.

                :param _reads: A current :py:class:`bsf.ngs.Reads` object.
                :type _reads: Reads | None
                :param _reads_file: A file path.
                :type _reads_file: str
                :param _reads_name: A name.
                :type _reads_name: str
                :param _file_type: A file type.
                :type _file_type: str
                :return: A new :py:class:`bsf.ngs.Reads` object.
                :rtype: Reads | None
                """
                module_logger.debug(
                    f'process_reads_new('
                    f'_reads={_reads!r}, '
                    f'_reads_file={_reads_file!r}, '
                    f'_reads_name={_reads_name!r}, '
                    f'_file_type={_file_type!r}'
                    f')'
                )

                return Reads(name=_reads_name, file_path=_reads_file, file_type=_file_type)

            def process_reads_new_file(_reads, _reads_file, _reads_name, _file_type):
                """Private function to process (i.e., initialise) a new
                :py:class:`bsf.ngs.Reads` object with file_path only.

                :param _reads: A current :py:class:`bsf.ngs.Reads` object.
                :type _reads: Reads | None
                :param _reads_file: A file path.
                :type _reads_file: str
                :param _reads_name: A name.
                :type _reads_name: str
                :param _file_type: A file type.
                :type _file_type: str
                :return: A new :py:class:`bsf.ngs.Reads` object.
                :rtype: Reads | None
                """
                module_logger.debug(
                    f'process_reads_new_file('
                    f'_reads={_reads!r}, '
                    f'_reads_file={_reads_file!r}, '
                    f'_reads_name={_reads_name!r}, '
                    f'_file_type={_file_type!r}'
                    f')'
                )

                return Reads(name=None, file_path=_reads_file, file_type=_file_type)

            def process_reads_new_name(_reads, _reads_file, _reads_name, _file_type):
                """Private function to process (i.e., initialise) a new
                :py:class:`bsf.ngs.Reads` object with :py:attr:`bsf.ngs.Reads.name` only.

                :param _reads: A current :py:class:`bsf.ngs.Reads` object.
                :type _reads: Reads | None
                :param _reads_file: A file path.
                :type _reads_file: str
                :param _reads_name: A name.
                :type _reads_name: str
                :param _file_type: A file type.
                :type _file_type: str
                :return: A new :py:class:`bsf.ngs.Reads` object.
                :rtype: Reads | None
                """
                module_logger.debug(
                    f'process_reads_new_name('
                    f'_reads={_reads!r}, '
                    f'_reads_file={_reads_file!r}, '
                    f'_reads_name={_reads_name!r}, '
                    f'_file_type={_file_type!r}'
                    f')'
                )

                return Reads(name=_reads_name, file_path=None, file_type=_file_type)

            def process_reads_old(_reads, _reads_file, _reads_name, _file_type):
                """Private function to process (i.e. complete) an old :py:class:`bsf.ngs.Reads` object.

                :param _reads: A current :py:class:`bsf.ngs.Reads` object.
                :type _reads: Reads | None
                :param _reads_file: A file path.
                :type _reads_file: str
                :param _reads_name: A name.
                :type _reads_name: str
                :param _file_type: A file type.
                :type _file_type: str
                :return: A new :py:class:`bsf.ngs.Reads` object.
                :rtype: Reads | None
                """
                module_logger.debug(
                    f'process_reads_old('
                    f'_reads={_reads!r}, '
                    f'_reads_file={_reads_file!r}, '
                    f'_reads_name={_reads_name!r}, '
                    f'_file_type={_file_type!r}'
                    f')'
                )

                if _reads is None:
                    return

                if _reads_file:
                    _reads.file_path = _reads_file

                if _reads_name:
                    _reads.name = _reads_name

                if _file_type:
                    _reads.file_type = _file_type

                return _reads

            module_logger.debug(
                f'process_reads('
                f'reads={reads!r}, '
                f'suffix={suffix!r}, '
                f'default_path={default_path!r}'
                f')'
            )

            # Cache the Python dict of function pointers as populating the dict is most likely expensive.

            if cls._reads_function_dict is None:
                cls._reads_function_dict = {
                    # File | Name
                    # ION ION  (Identical|Old|New)
                    # 000 000  (0) -> old Reads (nothing new)
                    0: process_reads_old,
                    # 000 001  (1) -> new Reads (new name)
                    1: process_reads_new,
                    # 000 010  (2) -> old Reads (old name)
                    2: process_reads_old,
                    # 000 011  (3) -> new Reads (new name)
                    3: process_reads_new,
                    # 000 100  (4) -> impossible (old and new names not defined, but identical)
                    # 000 101  (5) -> impossible (old name not defined but identical)
                    # 000 110  (6) -> impossible (new name not defined but identical)
                    # 000 111  (7) -> old Reads (identical old and new name)
                    7: process_reads_old,
                    # 001 000  (8) -> new Reads (new file)
                    8: process_reads_new,
                    # 001 001  (9) -> new Reads (new file and new name)
                    9: process_reads_new,
                    # 001 010 (10) -> old Reads (new file - old name)
                    10: process_reads_old,
                    # 001 011 (11) -> new Reads (new file and new name - old name)
                    11: process_reads_new,
                    # 001 100 (12) -> impossible (old and new names not defined, but identical)
                    # 001 101 (13) -> impossible (old name not defined, but identical)
                    # 001 110 (14) -> impossible (new name not defined, but identical)
                    # 001 111 (15) -> new Reads file (new file - identical old and new name)
                    15: process_reads_new_file,
                    # 010 000 (16) -> old Reads (old file)
                    16: process_reads_old,
                    # 010 001 (17) -> old Reads (old file, new name)
                    17: process_reads_old,
                    # 010 010 (18) -> old Reads (old file, old name)
                    18: process_reads_old,
                    # 010 011 (19) -> new Reads (new name)
                    19: process_reads_new,
                    # 010 100 (20) -> impossible (names not defined, but identical)
                    # 010 101 (21) -> impossible (old name not defined but identical)
                    # 010 110 (22) -> impossible (new name not defined but identical)
                    # 010 111 (23) -> old Reads (old file, old name, new identical name)
                    23: process_reads_old,
                    # 011 000 (24) -> new Reads (new file)
                    24: process_reads_new,
                    # 011 001 (25) -> new Reads (new file and new name)
                    25: process_reads_new,
                    # 011 010 (26) -> new Reads (new file)
                    26: process_reads_new,
                    # 011 011 (27) -> new Reads (new file and new name)
                    27: process_reads_new,
                    # 011 100 (28) -> impossible (names not defined, but supposedly identical)
                    # 011 101 (29) -> impossible (old name not defined, but supposedly identical)
                    # 011 110 (30) -> impossible (new name not defined, but supposedly identical)
                    # 011 111 (31) -> new Reads file (new file, identical names)
                    31: process_reads_new_file,
                    # 100 XXX (32-39) -> impossible (no files, but supposedly identical)
                    # 101 XXX (40-47) -> impossible (no old file, but supposedly identical)
                    # 110 XXX (48-55) -> impossible (no new file, but supposedly identical)
                    # 111 000 (56) -> old Reads (same old and new files)
                    56: process_reads_old,
                    # 111 001 (57) -> new Reads name (same file, but new name)
                    57: process_reads_new_name,
                    # 111 010 (58) -> old Reads (identical file, old name)
                    58: process_reads_old,
                    # 111 011 (59) -> new Reads name (identical file, but new name)
                    59: process_reads_new_name,
                    # 111 100 (60) -> impossible (names not defined, but supposedly identical)
                    # 111 101 (61) -> impossible (old name not defined, but supposedly identical)
                    # 111 110 (62) -> impossible (new name not defined, but supposedly identical)
                    # 111 111 (63) -> old Reads (all values defined and identical)
                    63: process_reads_old,
                }

            reads_file: Optional[str] = None
            reads_name: Optional[str] = None

            # Reads{suffix} File
            _key = f'Reads{suffix} File'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                # The key exists ...
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    reads_file = Configuration.get_absolute_path(
                        file_path=_value,
                        default_path=default_path)

            # Reads{suffix} Name
            _key = f'Reads{suffix} Name'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                # The key exists ...
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    reads_name = _value

            reads_status = 0
            if reads is None:
                reads_status |= process_reads_status(new=reads_file, old=None) << 3
                reads_status |= process_reads_status(new=reads_name, old=None)
            else:
                reads_status |= process_reads_status(new=reads_file, old=reads.file_path) << 3
                reads_status |= process_reads_status(new=reads_name, old=reads.name)

            reads_new = cls._reads_function_dict[reads_status](
                _reads=reads,
                _reads_file=reads_file,
                _reads_name=reads_name,
                _file_type=file_type)

            module_logger.debug(
                f'process_reads final: '
                f'suffix: {suffix!r} '
                f'Reads: {reads!r} '
                f'reads_file: {reads_file!r} '
                f'reads_name: {reads_name!r} '
                f'reads_status: {reads_status!r} '
                f'reads_new: {reads_new!r}'
            )

            return reads_new

        def process_paired_reads(sample, paired_reads, default_path):
            """Get or create a :py:class:`bsf.ngs.PairedReads` object.

            The :literal:`[Prefix] PairedReads Exclude` key is optional,
            in which case the default is Python :py:class:`bool` :py:const:`False`.
            The :literal:`[Prefix] PairedReads Index 1`, :literal:`[Prefix] PairedReads Index 2` and
            :literal:`[Prefix] PairedReads ReadGroup` keys are optional,
            in which case the default is an empty Python :py:class:`str` object.

            :param sample: A :py:class:`bsf.ngs.Sample` object.
            :type sample: Sample
            :param paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
            :type paired_reads: PairedReads | None
            :param default_path: A default file path.
            :type default_path: str
            :return: A :py:class:`bsf.ngs.PairedReads` object.
            :rtype: PairedReads | None
            """

            def process_paired_reads_status(new, old):
                """Calculate a comparison status.
                    - 2**0 1: new defined
                    - 2**1 2: old defined
                    - 2**2 4: new is args_1

                :param new: A new :py:class:`bsf.ngs.Reads` object.
                :type new: Reads | None
                :param old: An old :py:class:`bsf.ngs.Reads` object.
                :type old: Reads | None
                :return: A comparison status.
                :rtype: int
                """
                status = 0

                if new is not None:
                    status |= 1 << 0  # 1

                if old is not None:
                    status |= 1 << 1  # 2

                if (status == 3) and (new is old):
                    status |= 1 << 2  # 4

                return status

            def process_paired_reads_new(_sample, _paired_reads, _file_type, _reads_1, _reads_2):
                """Private function to process (i.e. initialise) a new :py:class:`bsf.ngs.PairedReads` object.

                :param _sample: A :py:class:`bsf.ngs.Sample` object.
                :type _sample: Sample
                :param _paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
                :type _paired_reads: PairedReads | None
                :param _file_type: A file type.
                :type _file_type: str
                :param _reads_1: A first :py:class:`bsf.ngs.Reads` object.
                :type _reads_1: Reads | None
                :param _reads_2: A second :py:class:`bsf.ngs.Reads` object.
                :type _reads_2: Reads | None
                :return: A :py:class:`bsf.ngs.PairedReads` object.
                :rtype: PairedReads | None
                """
                module_logger.debug(
                    f'process_paired_reads_new('
                    f'_sample.name={_sample.name!r}, '
                    f'_pairedReads={_paired_reads!r}, '
                    f'_file_type={_file_type!r}, '
                    f'_reads_1={_reads_1!r}, '
                    f'_reads_2={_reads_2!r}'
                    f')'
                )

                _paired_reads = PairedReads(file_type=_file_type, reads_1=_reads_1, reads_2=_reads_2)

                _sample.add_paired_reads(paired_reads=_paired_reads)

                return _paired_reads

            def process_paired_reads_new_1(_sample, _paired_reads, _file_type, _reads_1, _reads_2):
                """Private function to process (i.e. initialise) a new :py:class:`bsf.ngs.PairedReads` object
                with only :py:attr:`bsf.ngs.PairedReads.reads_1`.

                :param _sample: A :py:class:`bsf.ngs.Sample` object.
                :type _sample: Sample
                :param _paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
                :type _paired_reads: PairedReads | None
                :param _file_type: A file type.
                :type _file_type: str
                :param _reads_1: A first :py:class:`bsf.ngs.Reads` object.
                :type _reads_1: Reads | None
                :param _reads_2: A second :py:class:`bsf.ngs.Reads` object.
                :type _reads_2: Reads | None
                :return: A :py:class:`bsf.ngs.PairedReads` object.
                :rtype: PairedReads | None
                """
                module_logger.debug(
                    f'process_paired_reads_new_1('
                    f'_sample.name={_sample.name!r}, '
                    f'_pairedReads={_paired_reads!r}, '
                    f'_reads_1={_reads_1!r}, '
                    f'_reads_2={_reads_2!r}'
                    f')'
                )

                _paired_reads = PairedReads(file_type=_file_type, reads_1=_reads_1, reads_2=None)

                _sample.add_paired_reads(paired_reads=_paired_reads)

                return _paired_reads

            def process_paired_reads_new_2(_sample, _paired_reads, _file_type, _reads_1, _reads_2):
                """Private function to process (i.e. initialise) a new :py:class:`bsf.ngs.PairedReads` object
                with only :py:attr:`bsf.ngs.PairedReads.reads_2`.

                :param _sample: A :py:class:`bsf.ngs.Sample` object.
                :type _sample: Sample
                :param _paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
                :type _paired_reads: PairedReads | None
                :param _file_type: A file type.
                :type _file_type: str
                :param _reads_1: A first :py:class:`bsf.ngs.Reads` object.
                :type _reads_1: Reads | None
                :param _reads_2: A second :py:class:`bsf.ngs.Reads` object.
                :type _reads_2: Reads | None
                :return: A :py:class:`bsf.ngs.PairedReads` object.
                :rtype: PairedReads | None
                """
                module_logger.debug(
                    f'process_paired_reads_new_2('
                    f'_sample.name={_sample.name!r}, '
                    f'_paired_reads={_paired_reads!r}, '
                    f'_reads_1={_reads_1!r}, '
                    f'_reads_2={_reads_2!r}'
                    f')'
                )

                _paired_reads = PairedReads(file_type=_file_type, reads_1=None, reads_2=_reads_2)

                _sample.add_paired_reads(paired_reads=_paired_reads)

                return _paired_reads

            def process_paired_reads_old(_sample, _paired_reads, _file_type, _reads_1, _reads_2):
                """Private function to process (i.e. complete) an old :py:class:`bsf.ngs.PairedReads` object.

                :param _sample: A :py:class:`bsf.ngs.Sample` object.
                :type _sample: Sample
                :param _paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
                :type _paired_reads: PairedReads | None
                :param _file_type: A file type.
                :type _file_type: str
                :param _reads_1: A first :py:class:`bsf.ngs.Reads` object.
                :type _reads_1: Reads | None
                :param _reads_2: A second :py:class:`bsf.ngs.Reads` object.
                :type _reads_2: Reads | None
                :return: A :py:class:`bsf.ngs.PairedReads` object.
                :rtype: PairedReads | None
                """
                module_logger.debug(
                    f'process_paired_reads_old('
                    f'_sample.name={_sample.name!r}, '
                    f'_paired_reads={_paired_reads!r}, '
                    f'_reads_1={_reads_1!r}, '
                    f'_reads_2={_reads_2!r}'
                    f')'
                )

                if _paired_reads is None:
                    return

                if _file_type:
                    _paired_reads.file_type = _file_type

                if _reads_1 is not None:
                    _paired_reads.reads_1 = _reads_1

                if _reads_2 is not None:
                    _paired_reads.reads_2 = _reads_2

                return _paired_reads

            module_logger.debug(
                f'process_paired_reads('
                f'sample.name={sample.name!r}, '
                f'paired_reads={paired_reads!r}, '
                f'default_path={default_path!r}'
                f')'
            )

            # Cache the Python dict of function pointers as populating the dict is most likely expensive.

            if cls._paired_reads_function_dict is None:
                cls._paired_reads_function_dict = {
                    # Reads1 | Reads2
                    # ION ION  (Identical|Old|New)
                    # 000 000  (0) -> old PairedReads (nothing new)
                    0: process_paired_reads_old,
                    # 000 001  (1) -> new PairedReads (new Reads2 object)
                    1: process_paired_reads_new,
                    # 000 010  (2) -> old PairedReads (old Reads2 object)
                    2: process_paired_reads_old,
                    # 000 011  (3) -> new PairedReads (old and new Reads2 objects, not identical)
                    3: process_paired_reads_new,
                    # 000 100  (4) -> impossible (old and new Reads2 objects not defined, but supposedly identical)
                    # 000 101  (5) -> impossible (old Reads2 object not defined but supposedly identical)
                    # 000 110  (6) -> impossible (new Reads2 object not defined but supposedly identical)
                    # 000 111  (7) -> old PairedReads (identical old and new Reads2 object)
                    7: process_paired_reads_old,
                    # 001 000  (8) -> new PairedReads (new Reads1 object)
                    8: process_paired_reads_new,
                    # 001 001  (9) -> new PairedReads (new Reads1 and new Reads2 objects)
                    9: process_paired_reads_new,
                    # 001 010 (10) -> old PairedReads (new Reads1 and old Reads2 objects)
                    10: process_paired_reads_old,
                    # 001 011 (11) -> new PairedReads (new Reads1 and new Reads2 - old Reads2 objects)
                    11: process_paired_reads_new,
                    # 001 100 (12) -> impossible (old and new Reads2 not defined, but supposedly identical)
                    # 001 101 (13) -> impossible (old Reads2 not defined, but supposedly identical)
                    # 001 110 (14) -> impossible (new Reads2 not defined, but supposedly identical)
                    # 001 111 (15) -> new PairedReads1 (new Reads1 - identical old and new Reads2 objects)
                    15: process_paired_reads_new_1,
                    # 010 000 (16) -> old PairedReads (old Reads1 object)
                    16: process_paired_reads_old,
                    # 010 001 (17) -> old PairedReads (old Reads1 object - new Reads2 object)
                    17: process_paired_reads_old,
                    # 010 010 (18) -> old PairedReads (old Reads1 object - old Reads2 object)
                    18: process_paired_reads_old,
                    # 010 011 (19) -> new PairedReads (new Reads2 object)
                    19: process_paired_reads_new,
                    # 010 100 (20) -> impossible (Reads2 object not defined, but supposedly identical)
                    # 010 101 (21) -> impossible (old Reads2 object not defined but supposedly identical)
                    # 010 110 (22) -> impossible (new Reads2 object not defined but supposedly identical)
                    # 010 111 (23) -> old PairedReads (old Reads1, old Reads2, new identical Reads2 object)
                    23: process_paired_reads_old,
                    # 011 000 (24) -> new PairedReads (new Reads1 object)
                    24: process_paired_reads_new,
                    # 011 001 (25) -> new PairedReads (new Reads1 and new Reads2 objects)
                    25: process_paired_reads_new,
                    # 011 010 (26) -> new PairedReads (new Reads1 object)
                    26: process_paired_reads_new,
                    # 011 011 (27) -> new PairedReads (new Reads1 and new Reads2 objects)
                    27: process_paired_reads_new,
                    # 011 100 (28) -> impossible (Reads2 not defined, but supposedly identical)
                    # 011 101 (29) -> impossible (old Reads2 not defined, but supposedly identical)
                    # 011 110 (30) -> impossible (new Reads2 not defined, but supposedly identical)
                    # 011 111 (31) -> new PairedReads1 (different Reads1, identical Reads2)
                    # NOTE: In process_reads() an Exception gets raised.
                    31: process_paired_reads_new_1,
                    # 100 XXX (32-39) -> impossible (no Reads1, but supposedly identical)
                    # 101 XXX (40-47) -> impossible (no old Reads1, but supposedly identical)
                    # 110 XXX (48-55) -> impossible (no new Reads1, but supposedly identical)
                    # 111 000 (56) -> old PairedReads (same old and new Reads1)
                    56: process_paired_reads_old,
                    # 111 001 (57) -> New PairedReads2 (same Reads1, but new Reads2)
                    # NOTE: In process_reads() an Exception gets raised.
                    57: process_paired_reads_new_2,
                    # 111 010 (58) -> old PairedReads (identical Reads1, old Reads2)
                    58: process_paired_reads_old,
                    # 111 011 (59) -> new PairedReads2 (identical Reads1, but different Reads2 objects)
                    59: process_paired_reads_new_2,
                    # 111 100 (60) -> impossible (Reads2 objects not defined, but supposedly identical)
                    # 111 101 (61) -> impossible (old Reads2 not defined, but supposedly identical)
                    # 111 110 (62) -> impossible (new Reads2 not defined, but supposedly identical)
                    # 111 111 (63) -> old PairedReads (all objects defined and identical)
                    63: process_paired_reads_old,
                }

            paired_reads_status = 0

            if paired_reads is None:
                reads_1 = process_reads(reads=None, suffix='1', default_path=default_path)
                reads_2 = process_reads(reads=None, suffix='2', default_path=default_path)
                paired_reads_status |= process_paired_reads_status(new=reads_1, old=None) << 3
                paired_reads_status |= process_paired_reads_status(new=reads_2, old=None)
            else:
                reads_1 = process_reads(reads=paired_reads.reads_1, suffix='1', default_path=default_path)
                reads_2 = process_reads(reads=paired_reads.reads_2, suffix='2', default_path=default_path)
                paired_reads_status |= process_paired_reads_status(new=reads_1, old=paired_reads.reads_1) << 3
                paired_reads_status |= process_paired_reads_status(new=reads_2, old=paired_reads.reads_2)

            paired_reads_new = cls._paired_reads_function_dict[paired_reads_status](
                _sample=sample,
                _paired_reads=paired_reads,
                _file_type=file_type,
                _reads_1=reads_1,
                _reads_2=reads_2)

            if paired_reads_new is None:
                # If a PairedReads object is not defined at this stage it cannot be annotated.
                # Remove all PairedReads strings from the key list.
                _key = 'PairedReads'
                if sas_prefix:
                    _key = sas_prefix + ' ' + _key

                for _component in list(key_list):
                    if _component.startswith(_key):
                        key_list.remove(_component)
                        if row_dict[_component]:
                            module_logger.warning(
                                f'Dismissed PairedReads annotation '
                                f'key: {_component} '
                                f'value: {row_dict[_component]!r}')

                return

            _key = 'PairedReads Exclude'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    paired_reads_new.exclude = sas.get_boolean(row_dict=row_dict, key=_key)

            _key = 'PairedReads Index 1'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    paired_reads_new.index_1 = _value

            _key = 'PairedReads Index 2'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                key_list.remove(_key)
                _value = row_dict[_key]
                if _value:
                    paired_reads_new.index_2 = _value

            _key = 'PairedReads ReadGroup'
            if sas_prefix:
                _key = sas_prefix + ' ' + _key

            if _key in row_dict:
                _value = row_dict[_key]
                key_list.remove(_key)
                if _value:
                    paired_reads_new.read_group = _value

            paired_reads_new.process_annotation(row_dict=row_dict, key_list=key_list, prefix=sas_prefix)

            module_logger.debug(
                f'process_paired_reads final: '
                f'Sample.name: {sample.name!r} '
                f'Reads_1: {reads_1!r} '
                f'Reads_2: {reads_2!r} '
                f'status: {paired_reads_status!r} '
                f'paired_reads_new: {paired_reads_new!r}'
            )

            return paired_reads_new

        assert isinstance(sas, SampleAnnotationSheet)

        current_collection = cls(file_path=file_path, file_type=file_type, name=name)

        for row_dict in sas.row_dicts:
            module_logger.debug(f'from_sas: row_dict={row_dict!r}')

            # Generate a Python list of key objects since the list is subsequently modified.
            key_list = [key for key in row_dict]
            file_type = process_file_type()
            current_prf = process_processed_run_folder(collection=current_collection, prf=current_prf)
            current_project = process_project(prf=current_prf, project=current_project)
            current_sample = process_sample(project=current_project, sample=current_sample)
            current_paired_reads = process_paired_reads(
                sample=current_sample,
                paired_reads=current_paired_reads,
                default_path=current_collection.file_path)

            if len(key_list):
                module_logger.warning(f'Unexpected keys in SampleAnnotationSheet: {key_list!r} Row: {row_dict!r}')

        # Quench empty default objects that are a consequence of empty lines in the sample annotation sheet.
        # Use a less efficient list() constructor here to allow for modification of the dict object while iterating.

        for _prf in list(current_collection.processed_run_folder_dict.values()):
            for _project in list(_prf.project_dict.values()):
                for _sample in list(_project.sample_dict.values()):
                    if _sample.name == Sample.default_name and not len(_sample.paired_reads_list):
                        _project.del_sample(name=_sample.name)
                if _project.name == Project.default_name and not len(_project.sample_dict):
                    _prf.del_project(name=_project.name)
            if _prf.name == ProcessedRunFolder.default_name and not _prf.project_dict:
                current_collection.del_processed_run_folder(name=_prf.name)

        # Group Sample objects on the basis of 'Sample Group' annotation.

        for _prf in current_collection.processed_run_folder_dict.values():
            for _project in _prf.project_dict.values():
                for _sample in _project.sample_dict.values():
                    if 'Group' in _sample.annotation_dict:
                        for _sample_group_name in _sample.annotation_dict['Group']:
                            if _sample_group_name not in current_collection.sample_group_dict:
                                current_collection.sample_group_dict[_sample_group_name] = list()
                            sample_list = current_collection.sample_group_dict[_sample_group_name]
                            if _sample not in sample_list:
                                sample_list.append(_sample)

        return current_collection

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            processed_run_folder_dict=None,
            sample_group_dict=None):
        """Initialise a :py:class:`bsf.ngs.Collection` object.

        :param name: A name.
        :type name: str | None
        :param file_path: A file path.
        :type file_path: str | None
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :param annotation_dict: A Python :py:class:`dict` object for annotation of Python :py:class:`str` key and
            Python :py:class:`list` value objects of Python :py:class:`str` objects.
        :type annotation_dict: dict[str, list[str]] | None
        :param processed_run_folder_dict: A Python :py:class:`dict` object of
            :py:attr:`bsf.ngs.ProcessedRunFolder.name` key objects and
            :py:class:`bsf.ngs.ProcessedRunFolder` value objects.
        :type processed_run_folder_dict: dict[str, ProcessedRunFolder] | None
        :param sample_group_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (group name) key objects and
            second-level Python :py:class:`dict` value objects of
            :py:attr:`bsf.ngs.Sample.name` key objects and
            :py:class:`bsf.ngs.Sample` value objects.
        :type sample_group_dict: dict[str, list[Sample]] | None
        """
        super(Collection, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        if processed_run_folder_dict is None:
            self.processed_run_folder_dict = dict()
        else:
            self.processed_run_folder_dict = processed_run_folder_dict
            # Set this bsf.ngs.Collection as weak reference for each bsf.ngs.ProcessedRunFolder in the dict.
            for prf in self.processed_run_folder_dict.values():
                prf.weak_reference_collection = ReferenceType(self)

        if sample_group_dict is None:
            self.sample_group_dict = dict()
        else:
            self.sample_group_dict = sample_group_dict

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.ngs.Collection` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append(f'{indent}{self.__class__!r}\n')
        str_list.append(f'{indent}  name:      {self.name!r}\n')
        str_list.append(f'{indent}  file_path: {self.file_path!r}\n')
        str_list.append(f'{indent}  file_type: {self.file_type!r}\n')
        str_list.append(f'{indent}  annotation_dict:\n')
        for annotation_key in sorted(self.annotation_dict):
            str_list.append(f'{indent}    {annotation_key!r} {self.annotation_dict[annotation_key]!r}\n')
        str_list.append(f'{indent}  processed_run_folder_dict:\n')
        for prf_name in sorted(self.processed_run_folder_dict):
            str_list.extend(self.processed_run_folder_dict[prf_name].trace(level=level + 1))
        str_list.append(f'{indent}  sample_group_dict:\n')
        for sample_group_name in sorted(self.sample_group_dict):
            str_list.append(f'{indent}    group: {sample_group_name!r}\n')
            # List all Sample objects of this Python list object.
            for sample in self.sample_group_dict[sample_group_name]:
                str_list.append(f'{indent}      Sample name: {sample.name!r} file_path: {sample.file_path!r}\n')

        return str_list

    def add_processed_run_folder(self, prf):
        """Add a :py:class:`bsf.ngs.ProcessedRunFolder` object to the :py:class:`bsf.ngs.Collection` object and
        set a weak reference in the :py:class:`bsf.ngs.ProcessedRunFolder` object back to the
        :py:class:`bsf.ngs.Collection` object.

        :param prf: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
        :type prf: ProcessedRunFolder
        :return: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
        :rtype: ProcessedRunFolder
        """
        assert isinstance(prf, ProcessedRunFolder)
        # Delete an eventual bsf.ngs.ProcessedRunFolder that may exist under the same bsf.ngs.ProcessedRunFolder.name.
        self.del_processed_run_folder(name=prf.name)
        self.processed_run_folder_dict[prf.name] = prf
        prf.weak_reference_collection = ReferenceType(self)

        return prf

    def del_processed_run_folder(self, name):
        """Delete a :py:class:`bsf.ngs.ProcessedRunFolder` object from a :py:class:`bsf.ngs.Collection` object and
        clear the weak reference, if it points back at the :py:class:`bsf.ngs.Collection` object.

        :param name: A :py:attr:`bsf.ngs.ProcessedRunFolder.name` attribute.
        :type name: str
        :return: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
        :rtype: ProcessedRunFolder
        """
        if name in self.processed_run_folder_dict:
            prf = self.processed_run_folder_dict[name]
            del self.processed_run_folder_dict[name]
            if (prf.weak_reference_collection is not None) and (prf.weak_reference_collection() is self):
                prf.weak_reference_collection = None
            return prf
        else:
            return

    def get_processed_run_folder(self, file_path, file_type=None):
        """Get a :py:class:`bsf.ngs.ProcessedRunFolder` object by file path.

        If the :py:class:`bsf.ngs.ProcessedRunFolder` object does not exist in the
        :py:class:`bsf.ngs.Collection` object, it will be automatically discovered and added.

        :param file_path: A file path.
        :type file_path: str
        :param file_type: A file type (i.e., :literal:`CASAVA`, :literal:`External` or :literal:`Automatic`).
        :type file_type: str | None
        :return: A :py:class:`bsf.ngs.ProcessedRunFolder` object.
        :rtype: ProcessedRunFolder
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_name in self.processed_run_folder_dict:
            return self.processed_run_folder_dict[file_name]
        else:
            return self.add_processed_run_folder(
                prf=ProcessedRunFolder.from_file_path(
                    file_path=Configuration.get_absolute_path(
                        file_path=file_path,
                        default_path=self.file_path),
                    file_type=file_type))

    def get_all_processed_run_folders(self):
        """Get an ordered Python :py:class:`list` object of :py:class:`bsf.ngs.ProcessedRunFolder` objects.

        :return: A Python :py:class:`list` object of :py:class:`bsf.ngs.ProcessedRunFolder` objects.
        :rtype: list[ProcessedRunFolder]
        """
        processed_run_folder_list: List[ProcessedRunFolder] = list()

        for processed_run_folder_name in sorted(self.processed_run_folder_dict):
            processed_run_folder_list.append(self.processed_run_folder_dict[processed_run_folder_name])

        return processed_run_folder_list

    def get_all_projects(self):
        """Get an ordered Python :py:class:`list` object of :py:class:`bsf.ngs.Project` objects.

        :return: A Python :py:class:`list` object of :py:class:`bsf.ngs.Project` objects.
        :rtype: list[Project]
        """
        project_list: List[Project] = list()

        for processed_run_folder in self.get_all_processed_run_folders():
            project_list.extend(processed_run_folder.get_all_projects())

        return project_list

    def get_all_samples(self, exclude=False):
        """Get an ordered Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.

        :param exclude: Exclude :py:class:`bsf.ngs.Sample` objects on the basis of
            :py:attr:`bsf.ngs.PairedReads.exclude` attributes.
        :type exclude: bool
        :return: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :rtype: list[Sample]
        """
        sample_list: List[Sample] = list()

        for project in self.get_all_projects():
            sample_list.extend(project.get_all_samples(exclude=exclude))

        return sample_list

    def get_sample_from_row_dict(self, row_dict, prefix=None):
        """Get a :py:class:`bsf.ngs.Sample` object from a :py:class:`bsf.ngs.SampleAnnotationSheet` row
        Python :py:class:`dict` object.

        Look-up a hierarchy of :py:class:`bsf.ngs.ProcessedRunFolder`, :py:class:`bsf.ngs.Project` and
        :py:class:`bsf.ngs.Sample` objects based on a :py:class:`bsf.ngs.SampleAnnotationSheet` row dictionary.
        :py:class:`bsf.ngs.ProcessedRunFolder` objects of file type :literal:`CASAVA` can be
        automatically discovered and registered.
        Return the corresponding :py:class:`bsf.ngs.Sample` object.

        :param row_dict: A :py:class:`bsf.ngs.SampleAnnotationSheet` row Python :py:class:`dict` object.
        :type row_dict: dict[str, str]
        :param prefix: An optional configuration prefix.
            (e.g., :literal:`[Control] Sample`, :literal:`[Treatment] Sample`, :literal:`[Point N] Sample`, ...)
        :type prefix: str | None
        :return: A :py:class:`bsf.ngs.Sample` object.
        :rtype: Sample
        """
        # NOTE: For the moment, the row_dict has to include keys for 'ProcessedRunFolder',
        # 'Project' and 'Sample'. Maybe, this method could search for a Sample name in
        # the Collection object. However, the problem is that a Collection contains
        # complete, auto-registered RunFolder objects. Thus, Sample names (typically 1, 2, 3, ...)
        # may clash. Therefore, it is probably best to stick to the three fields for
        # unambiguous sample resolution.

        key = 'ProcessedRunFolder Name'
        if prefix:
            key = prefix + ' ' + key

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = ProcessedRunFolder.default_name

        # The Collection.get_processed_run_folder method can automatically register
        # ProcessedRunFolder objects of file type 'CASAVA'.

        prf = self.get_processed_run_folder(file_path=value)

        key = 'Project Name'
        if prefix:
            key = prefix + ' ' + key

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Project.default_name

        project = prf.project_dict[value]

        key = 'Sample Name'
        if prefix:
            key = prefix + ' ' + key

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Sample.default_name

        sample = project.sample_dict[value]

        return sample

    def get_samples_from_row_dict(self, row_dict, prefix=None):
        """Get a Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects and a
        Python :py:class:`str` object of the Group column value
        as a Python :py:class:`tuple` object from a :py:class:`bsf.ngs.SampleAnnotationSheet` row
        Python :py:class:`dict` object.

        :param row_dict: A comparison CSV file row Python :py:class:`dict` object.
        :type row_dict: dict[str, str]
        :param prefix: An optional configuration prefix.
            (e.g., :literal:`[Control] Sample`, :literal:`[Treatment] Sample`, :literal:`[Point N] Sample`, ...)
        :type prefix: str | None
        :return: A Python :py:class:`tuple` object of
            a Python :py:class:`str` object of :literal:`[Prefix] Group` column value and
            a Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :rtype: (str, list[Sample])
        """
        sample_list: List[Sample] = list()

        value = str()

        key = 'Group'
        if prefix:
            key = prefix + ' ' + key

        # Does the key exist and is its value defined?
        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Extend the Python list with all Sample objects of this group.
            if value in self.sample_group_dict:
                sample_list.extend(self.sample_group_dict[value])

        key = 'Sample Name'
        if prefix:
            key = prefix + ' ' + key

        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Append the Sample object, if defined.

            sample = self.get_sample_from_row_dict(row_dict=row_dict, prefix=prefix)

            if sample:
                sample_list.append(sample)

        return value, sample_list

    def to_sas(self, file_path=None, name=None):
        """Convert a :py:class:`bsf.ngs.Collection` object into a :py:class:`bsf.ngs.SampleAnnotationSheet` object.

        :param file_path: A file path.
        :type file_path: str | None
        :param name: A name.
        :type name: str | None
        :return: A :py:class:`bsf.ngs.SampleAnnotationSheet` object.
        :rtype: SampleAnnotationSheet
        """
        sas = SampleAnnotationSheet(file_path=file_path, name=name)

        # Scan the Collection and its contained objects for additional (annotation) field names.

        for prf in self.processed_run_folder_dict.values():
            if prf.annotation_dict:  # not None and not empty
                for prf_annotation_key in prf.annotation_dict.keys():
                    prf_annotation_field = ' '.join(('ProcessedRunFolder', prf_annotation_key))
                    if prf_annotation_field not in sas.field_names:
                        sas.field_names.append(prf_annotation_field)
            for project in prf.project_dict.values():
                if project.annotation_dict:  # not None and not empty
                    for project_annotation_key in project.annotation_dict.keys():
                        project_annotation_field = ' '.join(('Project', project_annotation_key))
                        if project_annotation_field not in sas.field_names:
                            sas.field_names.append(project_annotation_field)
                for sample in project.sample_dict.values():
                    if sample.annotation_dict:  # not None and not empty
                        for sample_annotation_key in sample.annotation_dict.keys():
                            sample_annotation_field = ' '.join(('Sample', sample_annotation_key))
                            if sample_annotation_field not in sas.field_names:
                                sas.field_names.append(sample_annotation_field)
                    for paired_reads in sample.paired_reads_list:
                        if paired_reads.annotation_dict:  # not None and not empty
                            for paired_reads_annotation_key in paired_reads.annotation_dict.keys():
                                paired_reads_annotation_field = ' '.join(('PairedReads', paired_reads_annotation_key))
                                if paired_reads_annotation_field not in sas.field_names:
                                    sas.field_names.append(paired_reads_annotation_field)

        # At this stage all annotation keys should be added. Partition and sort the list of field names.

        field_names = sorted(sas.field_names)  # Create a new sorted list.

        del sas.field_names[:]  # Clear the existing list.

        sas.field_names.extend([item for item in field_names if item.startswith('File Type')])
        field_names = [item for item in field_names if not item.startswith('File Type')]
        sas.field_names.extend([item for item in field_names if item.startswith('ProcessedRunFolder')])
        field_names = [item for item in field_names if not item.startswith('ProcessedRunFolder')]
        sas.field_names.extend([item for item in field_names if item.startswith('Project')])
        field_names = [item for item in field_names if not item.startswith('Project')]
        sas.field_names.extend([item for item in field_names if item.startswith('Sample')])
        field_names = [item for item in field_names if not item.startswith('Sample')]
        sas.field_names.extend([item for item in field_names if item.startswith('PairedReads')])
        field_names = [item for item in field_names if not item.startswith('PairedReads')]
        sas.field_names.extend([item for item in field_names if item.startswith('Reads1')])
        field_names = [item for item in field_names if not item.startswith('Reads1')]
        sas.field_names.extend([item for item in field_names if item.startswith('Reads2')])
        field_names = [item for item in field_names if not item.startswith('Reads2')]
        sas.field_names.extend(field_names)

        # Finally, construct the SampleAnnotationSheet.

        def row_dict_add(_row_dict, key, value):
            """Private function to add a :literal:`key` and :literal:`value` pair to a
            Python :py:class:`dict` object representing a row.

            If the :literal:`key` exists in the current Python :py:class:`dict` object,
            the Python :py:class:`dict` object is pushed onto the Python :py:class:`list` object of (row)
            Python :py:class:`dict` objects and a new one is created before the key and value pair is added.

            :param _row_dict: A row Python :py:class:`dict` object.
            :type _row_dict: dict[str, str]
            :param key: A key.
            :type key: str
            :param value: A value.
            :type value: str
            :return: A row Python :py:class:`dict` object.
            :rtype: dict[str, str]
            """
            if _row_dict is None:
                # If the row dict does not exist, create it.
                _row_dict = dict()

            if key in _row_dict:
                # If the key already exists in the row dict, append it and create a new one.
                sas.row_dicts.append(_row_dict)
                _row_dict = dict()

            _row_dict[key] = value

            return _row_dict

        def row_dict_complete(_row_dict):
            """Private function to complete a Python :py:class:`dict` object representing a row.

            :param _row_dict: A row Python :py:class:`dict` object.
            :type _row_dict: dict[str, str]
            """
            if row_dict:
                sas.row_dicts.append(_row_dict)

            return dict()

        row_dict: Optional[Dict[str, str]] = None

        for prf_name in sorted(self.processed_run_folder_dict):
            prf = self.processed_run_folder_dict[prf_name]
            row_dict = row_dict_add(_row_dict=row_dict, key='ProcessedRunFolder Name', value=prf.name)
            if prf.annotation_dict:  # not None and not empty
                for prf_annotation_key in sorted(prf.annotation_dict):
                    prf_annotation_list = prf.annotation_dict[prf_annotation_key]
                    prf_annotation_field = ' '.join(('ProcessedRunFolder', prf_annotation_key))
                    for annotation in prf_annotation_list:
                        row_dict = row_dict_add(_row_dict=row_dict, key=prf_annotation_field, value=annotation)
            for project_name in sorted(prf.project_dict):
                project = prf.project_dict[project_name]
                row_dict = row_dict_add(_row_dict=row_dict, key='Project Name', value=project.name)
                if project.annotation_dict:  # not None and not empty
                    for project_annotation_key in sorted(project.annotation_dict):
                        project_annotation_list = project.annotation_dict[project_annotation_key]
                        project_annotation_field = ' '.join(('Project', project_annotation_key))
                        for annotation in project_annotation_list:
                            row_dict = row_dict_add(_row_dict=row_dict, key=project_annotation_field, value=annotation)
                for sample_name in sorted(project.sample_dict):
                    sample = project.sample_dict[sample_name]
                    row_dict = row_dict_add(_row_dict=row_dict, key='Sample Name', value=sample.name)
                    if sample.annotation_dict:  # not None and not empty
                        for sample_annotation_key in sorted(sample.annotation_dict):
                            sample_annotation_list = sample.annotation_dict[sample_annotation_key]
                            sample_annotation_field = ' '.join(('Sample', sample_annotation_key))
                            for annotation in sample_annotation_list:
                                row_dict = row_dict_add(
                                    _row_dict=row_dict,
                                    key=sample_annotation_field,
                                    value=annotation)
                    for paired_reads in sample.paired_reads_list:
                        if paired_reads.exclude is None:
                            paired_reads.exclude = False
                        row_dict = row_dict_add(
                            _row_dict=row_dict,
                            key='PairedReads Exclude',
                            value=repr(paired_reads.exclude).upper())
                        row_dict = row_dict_add(
                            _row_dict=row_dict,
                            key='PairedReads Index 1',
                            value=paired_reads.index_1)
                        row_dict = row_dict_add(
                            _row_dict=row_dict,
                            key='PairedReads Index 2',
                            value=paired_reads.index_2)
                        row_dict = row_dict_add(
                            _row_dict=row_dict,
                            key='PairedReads ReadGroup',
                            value=paired_reads.read_group)
                        if paired_reads.reads_1 is not None:
                            row_dict = row_dict_add(
                                _row_dict=row_dict,
                                key='Reads1 Name',
                                value=paired_reads.reads_1.name)
                            row_dict = row_dict_add(
                                _row_dict=row_dict,
                                key='Reads1 File',
                                value=paired_reads.reads_1.file_path)
                        if paired_reads.reads_2 is not None:
                            row_dict = row_dict_add(
                                _row_dict=row_dict,
                                key='Reads2 Name',
                                value=paired_reads.reads_2.name)
                            row_dict = row_dict_add(
                                _row_dict=row_dict,
                                key='Reads2 File',
                                value=paired_reads.reads_2.file_path)
                        if paired_reads.annotation_dict:  # not None and not empty
                            for paired_reads_annotation_key in sorted(paired_reads.annotation_dict):
                                paired_reads_annotation_list = paired_reads.annotation_dict[paired_reads_annotation_key]
                                paired_reads_annotation_field = ' '.join(('PairedReads', paired_reads_annotation_key))
                                for annotation in paired_reads_annotation_list:
                                    row_dict = row_dict_add(
                                        _row_dict=row_dict,
                                        key=paired_reads_annotation_field,
                                        value=annotation)
                        # Complete the row dict at the end of an object so that fields higher up the hierarchy
                        # are not filled in.
                        row_dict = row_dict_complete(_row_dict=row_dict)
                    row_dict = row_dict_complete(_row_dict=row_dict)
                row_dict = row_dict_complete(_row_dict=row_dict)
            row_dict = row_dict_complete(_row_dict=row_dict)
        row_dict = row_dict_complete(_row_dict=row_dict)

        return sas

    def to_sas_path(self, file_path=None, name=None):
        """Write as a :py:class:`bsf.ngs.SampleAnnotationSheet` to a file path.

        :param file_path: A file path.
        :type file_path: str | None
        :param name: A name.
        :type name: str | None
        """
        sas = self.to_sas(file_path=file_path, name=name)
        sas.to_file_path()

        return


class SampleGroup(object):
    """The :py:class:`bsf.ngs.SampleGroup` class represents a named Python :py:class:`list` of
    :py:class:`bsf.ngs.Sample` objects.

    The grouping is usually defined in a sample annotation sheet.

    :ivar name: A name.
    :type name: str | None
    :ivar sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
    :type sample_list: list[Sample]
    """

    # Sample and PairedReads objects from different ProcessRunFolder objects
    # (i.e. flow cells) could bear the same name, leading to problems with SGE job names.
    # This would need further re-thinking.

    def __init__(
            self,
            name=None,
            sample_list=None):
        """Initialise a :py:class:`bsf.ngs.SampleGroup` object.

        :param name: A name.
        :type name: str | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects
        :type sample_list: list[Sample] | None
        """
        super(SampleGroup, self).__init__()

        self.name = name

        if sample_list is None:
            self.sample_list = list()
        else:
            self.sample_list = sample_list

        return

    def add_sample(self, sample):
        """Add a :py:class:`bsf.ngs.Sample` object.

        :param sample: A :py:class:`bsf.ngs.Sample` object.
        :type sample: Sample
        """
        assert isinstance(sample, Sample)
        if sample not in self.sample_list:
            self.sample_list.append(sample)

        return

    def is_excluded(self):
        """Is this :py:class:`bsf.ngs.SampleGroup` object excluded,
        because all its :py:class:`ngs.bsf.Sample` objects are?

        :return: :py:const:`True` if this :py:class:`bsf.ngs.SampleGroup` is excluded
            as all :py:class:`bsf.ngs.Sample` objects are excluded, :py:const:`False` otherwise.
        :rtype: bool
        """
        for sample in self.sample_list:
            if not sample.is_excluded():
                return False
        else:
            return True

    def get_all_paired_reads(self, replicate_grouping):
        """Get all :py:class:`bsf.ngs.PairedReads` objects of a :`py:class:`bsf.ngs.SampleGroup`.

        For the moment, replicates are :py:class:`bsf.ngs.PairedReads` objects that do not share
        anything but the chunk number.

        :param replicate_grouping: Group all :py:class:`bsf.ngs.PairedReads` objects of a
            :py:class:`bsf.ngs.Sample` object or list them individually.
        :type replicate_grouping: bool
        :return: A Python :py:class:`dict` object  of
            Python :py:class:`str` (replicate name) key and
            Python :py:class:`list` value objects of :py:class:`bsf.ngs.PairedReads` objects value data.
        :rtype: dict[str, list[PairedReads]]
        """
        group_dict: Dict[str, List[PairedReads]] = dict()

        for sample in self.sample_list:
            for paired_reads_name, paired_reads_list in \
                    sample.get_all_paired_reads(replicate_grouping=replicate_grouping).items():
                if paired_reads_name not in group_dict:
                    group_dict[paired_reads_name] = list()

                # Add PairedReads objects one-by-one and check if they are not already there.

                for paired_reads in paired_reads_list:
                    if paired_reads not in group_dict[paired_reads_name]:
                        group_dict[paired_reads_name].append(paired_reads)

        return group_dict


class SampleAnnotationSheet(AnnotationSheet):
    """The :py:class:`bsf.ngs.SampleAnnotationSheet` class represents a Comma-Separated Value (CSV) table of
    sample information after running the
    :py:class:`bsf.analyses.illumina_to_bam_tools.IlluminaToBamTools.BamIndexDecoder` object.
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'File Type',
        'ProcessedRunFolder Name',
        'Project Name',
        'Project Size',
        'Sample Name',
        'PairedReads Exclude',
        'PairedReads Index 1',
        'PairedReads Index 2',
        'PairedReads ReadGroup',
        'Reads1 Name',
        'Reads1 File',
        'Reads2 Name',
        'Reads2 File',
    ]

    _test_methods: Dict[str, List[Callable[[int, Dict[str, str], str], str]]] = {
        # File Type
        'ProcessedRunFolder Name': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Project Name': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Project Size': [
            AnnotationSheet.check_numeric,
        ],
        'Sample Name': [
            AnnotationSheet.check_alphanumeric,
        ],
        'PairedReads Index 1': [
            AnnotationSheet.check_sequence,
        ],
        'PairedReads Index 2': [
            AnnotationSheet.check_sequence,
        ],
        # PairedReads ReadGroup
        'Reads1 Name': [
            AnnotationSheet.check_alphanumeric,
        ],
        # Reads1 File
        'Reads2 Name': [
            AnnotationSheet.check_alphanumeric,
        ],
        # Reads2 File
    }
