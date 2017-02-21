"""bsf.ngs

A package of classes and methods modelling next-generation sequencing data directories and files.
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


import os
import re
import stat
import warnings
import weakref

from bsf.annotation import AnnotationSheet


class NextGenerationBase(object):
    """The C{bsf.ngs.NextGenerationBase} class represents a super-class for Next Generation Sequencing (NGS) objects.

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @type file_type: str
    @ivar annotation_dict: Python C{dict} for annotation of Python C{str} key and
        Python C{list} of Python C{str} value data
    @type annotation_dict: dict[str, list[str]]
    """

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None):
        """Initialise a C{bsf.ngs.Reads} object.

        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @return:
        @rtype:
        """
        super(NextGenerationBase, self).__init__()

        if name is None:
            self.name = str()
        else:
            self.name = name

        if file_path is None:
            self.file_path = str()
        else:
            self.file_path = file_path

        if file_type is None:
            self.file_type = str()
        else:
            self.file_type = file_type

        if annotation_dict is None:
            self.annotation_dict = dict()
        else:
            self.annotation_dict = annotation_dict

        return

    def __eq__(self, other):
        """Test C{bsf.ngs.NextGenerationBase} objects for equality.

        @param other: C{bsf.ngs.NextGenerationBase}
        @type other: bsf.ngs.NextGenerationBase
        @return: C{True} if equal, C{False} otherwise
        @rtype: bool
        """
        assert isinstance(other, NextGenerationBase)

        if self is other:
            return True

        return self.name == other.name \
            and self.file_path == other.file_path \
            and self.file_type == other.file_type \
            and self.annotation_dict == other.annotation_dict

    def __ne__(self, other):
        """Test C{bsf.ngs.NextGenerationBase} objects for inequality.

        @param other: C{bsf.ngs.NextGenerationBase}
        @type other: bsf.ngs.NextGenerationBase
        @return: C{True} if unequal, C{False} otherwise
        @rtype: bool
        """
        assert isinstance(other, NextGenerationBase)

        if self is not other:
            return True

        return self.name != other.name \
            or self.file_path != other.file_path \
            or self.file_type != other.file_type \
            or self.annotation_dict != other.annotation_dict

    def add_annotation(self, key, value):
        """Add an annotation value under a key.

        @param key: Annotation key
        @type key: str
        @param value: Annotation value
        @type value: str
        @return:
        @rtype:
        """

        if self.annotation_dict is None:
            self.annotation_dict = dict()

        if key not in self.annotation_dict:
            self.annotation_dict[key] = list()

        value_list = self.annotation_dict[key]

        if value not in value_list:
            value_list.append(value)

        return

    def process_annotation(self, row_dict, key_list, prefix):
        """Process annotation from a Python C{dict} of row entries of a Python C{csv} object.

        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] ReadsN', '[Treatment] ReadsN', '[Point N] ReadsN', ...)
        @type prefix: str
        @return:
        @rtype:
        """
        re_pattern = re.compile(pattern=' '.join((prefix, self.__class__.__name__)).lstrip())
        for key in list(key_list):
            # Only match the pattern at the start of the string.
            re_match = re.match(pattern=re_pattern, string=key)
            if re_match is not None:
                key_list.remove(re_match.string)
                # Capture the string from the end of the match to the end of the string and strip white space.
                annotation_key = re_match.string[re_match.end(0):].strip()
                if annotation_key and row_dict[re_match.string]:
                    # Exclude empty keys and values.
                    self.add_annotation(key=annotation_key, value=row_dict[re_match.string])

        return


class Reads(NextGenerationBase):
    """The C{bsf.ngs.Reads} class represents a file of Next-Generation Sequencing (NGS) reads.

    Typically, a C{Reads} object represents a FASTQ or unmapped BAM file.

    Attributes:
    @ivar barcode: Barcode used for sample multiplexing
    @type barcode: str
    @ivar lane: Lane number
    @type lane: str
    @ivar read: Read number (e.g. I{R1}, I{R2}, ...)
    @type read: str
    @ivar chunk: Chunk number (e.g. I{001}, I{002}, ...)
    @type chunk: str
    @ivar weak_reference_paired_reads: C{weakref.ReferenceType} pointing at a C{bsf.ngs.PairedReads} object
    @type weak_reference_paired_reads: weakref.ReferenceType
    """

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.ngs.Reads} object from a file path.

        For a I{file_type} I{CASAVA}, C{bsf.ngs.Reads.file_path} obeys a I{SampleName_Index_Lane_Read_Chunk} schema,
        so that C{bsf.ngs.Reads.name}, C{bsf.ngs.Reads.barcode}, C{bsf.ngs.Reads.lane}, C{bsf.ngs.Reads.read} and
        C{bsf.ngs.Reads.chunk} can be populated automatically.
        For I{file_type} I{External}, the attributes need populating manually.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.ngs.Reads} object
        @rtype: bsf.ngs.Reads
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_type == 'CASAVA':
            # CASAVA Reads obey a SampleName_Index_Lane_Read_Chunk schema.

            component_list = file_name.split('.')  # Split file extensions

            component_list[0] = component_list[0].split('_')  # Split components

            reads = cls(
                name='_'.join(component_list[0][:-4]),  # Exclude the SampleName
                file_path=file_path,
                file_type=file_type,
                annotation_dict=None,
                barcode=component_list[0][-4],
                lane=component_list[0][-3],
                read=component_list[0][-2],
                chunk=component_list[0][-1],
                weak_reference_paired_reads=None)
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return reads

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
        """Initialise a C{bsf.ngs.Reads} object.

        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @param barcode: Barcode used for sample multiplexing
        @type barcode: str
        @param lane: Lane number
        @type lane: str
        @param read: Read number (e.g. I{R1}, I{R2}, ...)
        @type read: str
        @param chunk: Chunk number (e.g. I{001}, I{002}, ...)
        @type chunk:str
        @param weak_reference_paired_reads: C{weakref.ReferenceType} pointing at a C{bsf.ngs.PairedReads} object
        @type weak_reference_paired_reads: weakref.ReferenceType
        @return:
        @rtype:
        """

        super(Reads, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        if barcode is None:
            self.barcode = str()
        else:
            self.barcode = barcode

        if lane is None:
            self.lane = str()
        else:
            self.lane = lane

        if read is None:
            self.read = str()
        else:
            self.read = read

        if chunk is None:
            self.chunk = str()
        else:
            self.chunk = chunk

        self.weak_reference_paired_reads = weak_reference_paired_reads  # Can be None.

        return

    def __eq__(self, other):
        """Test C{bsf.ngs.Reads} objects for equality.

        @param other: C{bsf.ngs.Reads}
        @type other: bsf.ngs.Reads
        @return: C{True} if equal, C{False} otherwise
        @rtype: bool
        """
        assert isinstance(other, Reads)

        return self is other \
            or super(Reads, self).__eq__(other=other) \
            and self.barcode == other.barcode \
            and self.lane == other.lane \
            and self.read == other.read \
            and self.chunk == other.chunk

    def __ne__(self, other):
        """Test C{bsf.ngs.Reads} objects for inequality.

        @param other: C{bsf.ngs.Reads}
        @type other: bsf.ngs.Reads
        @return: C{True} if unequal, C{False} otherwise
        @rtype: bool
        """
        assert isinstance(other, Reads)

        return self is not other \
            or super(Reads, self).__ne__(other=other) \
            or self.barcode != other.barcode \
            or self.lane != other.lane \
            or self.read != other.read \
            or self.chunk != other.chunk

    def __nonzero__(self):
        """Test C{bsf.ngs.Reads} objects for non-zero.

        @return: C{True} if non-zero, i.e file_path or name are meaningfully defined.
        @rtype: bool
        """
        if self.file_path and self.name:
            return True
        else:
            return False

    def trace(self, level):
        """Trace a C{bsf.ngs.Reads} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_paired_reads: {!r}\n'.format(indent, self.weak_reference_paired_reads)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  annotation_dict:\n'.format(indent, self.annotation_dict)
        for annotation_key in self.annotation_dict.keys():
            output += '{}    {!r} {!r}\n'.format(indent, annotation_key, self.annotation_dict[annotation_key])
        output += '{}  barcode:   {!r}\n'.format(indent, self.barcode)
        output += '{}  lane:      {!r}\n'.format(indent, self.lane)
        output += '{}  read:      {!r}\n'.format(indent, self.read)
        output += '{}  chunk:     {!r}\n'.format(indent, self.chunk)

        return output

    def match(self, reads):
        """Match C{bsf.ngs.Reads} objects.

        Two C{bsf.ngs.Reads} objects are identical, if all their instance variables match.
        @param reads: Second C{bsf.ngs.Reads} object
        @type reads: bsf.ngs.Reads
        @return: True if both objects match, False otherwise
        @rtype: bool
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
        """Match paired C{bsf.ngs.Reads} objects, by relaxing matching criteria.

        @param reads: Second C{bsf.ngs.Reads} object
        @type reads: bsf.ngs.Reads
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
        """

        assert isinstance(reads, Reads)

        if self.file_type == 'CASAVA':

            # All CASAVA Reads attributes need to be equal,
            # with the exception of the read (R1 or R2) and
            # the file_path, which also contains the reads
            # information.

            if self.name != reads.name:
                return False

            if self.file_type != reads.file_type:
                return False

            if self.barcode != reads.barcode:
                return False

            if self.lane != reads.lane:
                return False

            # Do not compare the read instance variable.

            if self.chunk != reads.chunk:
                return False

            return True
        else:
            # It is difficult to match PairedReads outside of CASAVA conventions.
            warnings.warn(
                'Matching of paired Reads objects for file_type other than CASAVA not implemented yet.',
                UserWarning)

            return True


class PairedReads(NextGenerationBase):
    """The C{bsf.ngs.PairedReads} class represents a pair of C{bsf.ngs.Reads}.

    For the C{bsf.ngs.Reads.file_type} I{CASAVA} this represents a read pair (i.e. I{R1} and I{R2}).
    Attributes:
    @ivar reads_1: First C{bsf.ngs.Reads} object
    @type reads_1: bsf.ngs.Reads
    @ivar reads_2: Second C{bsf.ngs.Reads} object
    @type reads_2: bsf.ngs.Reads
    @ivar exclude: Exclude from processing
    @type exclude: bool
    @ivar index_1: Index sequence 1
    @type index_1: str
    @ivar index_2: Index sequence 2
    @type index_2: str
    @ivar read_group: SAM read group (@RG) information
    @type read_group: str
    @ivar weak_reference_sample: C{weakref.ReferenceType} to a C{bsf.ngs.Sample}
    @type weak_reference_sample: weakref.ReferenceType
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
        """Initialise a C{bsf.ngs.PairedReads} object.

        For the C{bsf.ngs.Reads.file_type} I{CASAVA} the reads object will be
        automatically assigned on the basis of the C{bsf.ngs.Reads.read}
        attribute (i.e. I{R1} or I{R2}).
        Upon initialisation of this C{bsf.ngs.PairedReads} object, weak references (C{weakref.ReferenceType})
        are set in the C{bsf.ngs.Reads} objects.
        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @param reads_1: First C{bsf.ngs.Reads} object
        @type reads_1: bsf.ngs.Reads
        @param reads_2: Second C{bsf.ngs.Reads} object
        @type reads_2: bsf.ngs.Reads
        @param index_1: Index sequence 1
        @type index_1: str
        @param index_2: Index sequence 2
        @type index_2: str
        @param exclude: Exclude from processing
        @type exclude: bool
        @param read_group: SAM read group (@RG) information
        @type read_group: str
        @param weak_reference_sample: C{weakref.ReferenceType} pointing at a C{bsf.ngs.Sample}
        @type weak_reference_sample: weakref.ReferenceType
        @return:
        @rtype:
        @raise Exception: For C{bsf.ngs.Reads.file_type} I{CASAVA}, I{R1} or I{R2} must be set in the
            C{bsf.ngs.Reads} object.
        """

        if (reads_1 and reads_2) and (not reads_1.match_paired(reads=reads_2)):
            raise Exception('The Reads objects do not match.')

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
                    reads_1.weak_reference_paired_reads = weakref.ReferenceType(self)
                elif reads_1.read == 'R2':
                    self.reads_2 = reads_1
                    reads_1.weak_reference_paired_reads = weakref.ReferenceType(self)
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads_1.read))
            else:
                # Other file types go here ...
                self.reads_1 = reads_1
                reads_1.weak_reference_paired_reads = weakref.ReferenceType(self)

        if reads_2:
            if reads_2.file_type == 'CASAVA':
                if reads_2.read == 'R1':
                    self.reads_1 = reads_2
                    reads_2.weak_reference_paired_reads = weakref.ReferenceType(self)
                elif reads_2.read == 'R2':
                    self.reads_2 = reads_2
                    reads_2.weak_reference_paired_reads = weakref.ReferenceType(self)
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads_2.read))
            else:
                # Other file types go here ...
                self.reads_2 = reads_2
                reads_2.weak_reference_paired_reads = weakref.ReferenceType(self)

        if exclude is None:
            self.exclude = False
        else:
            self.exclude = exclude

        if index_1 is None:
            self.index_1 = str()
        else:
            self.index_1 = index_1

        if index_2 is None:
            self.index_2 = str()
        else:
            self.index_2 = index_2

        if read_group is None:
            self.read_group = str()
        else:
            self.read_group = read_group

        self.weak_reference_sample = weak_reference_sample  # Can be None.

        return

    def __eq__(self, other):
        """Test C{bsf.ngs.PairedReads} objects for equality.

        @param other: C{bsf.ngs.PairedReads}
        @type other: bsf.ngs.PairedReads
        @return: C{True} if equal, C{False} otherwise
        @rtype: bool
        """
        assert isinstance(other, PairedReads)

        return self is other \
            or super(PairedReads, self).__eq__(other=other) \
            and self.reads_1 == other.reads_1 \
            and self.reads_2 == other.reads_2 \
            and self.exclude == other.exclude \
            and self.index_1 == other.index_1 \
            and self.index_2 == other.index_2 \
            and self.read_group == other.read_group

    def __ne__(self, other):
        """Test C{bsf.ngs.PairedReads} objects for inequality.

        @param other: C{bsf.ngs.PairedReads}
        @type other: bsf.ngs.PairedReads
        @return: C{True} if unequal, C{False} otherwise
        @rtype: bool
        """
        assert isinstance(other, PairedReads)

        return self is not other \
            or super(PairedReads, self).__ne__(other=other) \
            or self.reads_1 != other.reads_1 \
            or self.reads_2 != other.reads_2 \
            or self.exclude != other.exclude \
            or self.index_1 != other.index_1 \
            or self.index_2 != other.index_2 \
            or self.read_group != other.read_group

    def trace(self, level):
        """Trace a C{bsf.ngs.PairedReads} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_sample: {!r}\n'.format(indent, self.weak_reference_sample)
        output += '{}  name:       {!r}\n'.format(indent, self.name)
        output += '{}  file_path:  {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type:  {!r}\n'.format(indent, self.file_type)
        output += '{}  annotation_dict:\n'.format(indent, self.annotation_dict)
        for annotation_key in self.annotation_dict.keys():
            output += '{}    {!r} {!r}\n'.format(indent, annotation_key, self.annotation_dict[annotation_key])
        output += '{}  reads_1:    {!r}\n'.format(indent, self.reads_1)
        output += '{}  reads_2:    {!r}\n'.format(indent, self.reads_2)
        output += '{}  index_1:    {!r}\n'.format(indent, self.index_1)
        output += '{}  index_2:    {!r}\n'.format(indent, self.index_2)
        output += '{}  exclude:    {!r}\n'.format(indent, self.exclude)
        output += '{}  read_group: {!r}\n'.format(indent, self.read_group)

        if isinstance(self.reads_1, Reads):
            output += self.reads_1.trace(level + 1)
        if isinstance(self.reads_2, Reads):
            output += self.reads_2.trace(level + 1)

        return output

    def add_reads(self, reads):
        """Add a C{bsf.ngs.Reads} object.

        For a C{bsf.ngs.Reads.file_type} I{CASAVA} the C{bsf.ngs.Reads} object can be automatically
        assigned on the basis of the C{bsf.ngs.Reads.read} attribute (i.e. I{R1} or I{R2}).
        @param reads: C{bsf.ngs.Reads}
        @type reads: bsf.ngs.Reads
        @return: C{True} upon success, C{False} otherwise
        @rtype: bool
        """

        # For CASAVA projects, reads are automatically added to either
        # reads_1 or reads_2 according to the file name.
        # Returns True upon success, False otherwise.

        assert isinstance(reads, Reads)

        if self.reads_1:
            if not self.reads_1.match_paired(reads=reads):
                return False

            if self.reads_1.file_type == 'CASAVA':
                if reads.read == 'R1':
                    raise Exception(
                        'PairedReads reads_1 has already been defined.\n'
                        '  reads_1: {!r}\n'
                        '  reads:   {!r}'.format(self.reads_1.file_path, reads.file_path))
                elif reads.read == 'R2':
                    self.reads_2 = reads
                    reads.weak_reference_paired_reads = weakref.ReferenceType(self)
                    return True
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads.read))
            else:
                # Other file types go here ...
                warnings.warn(
                    'Method not implemented for file types other than CASAVA.',
                    UserWarning)

        if self.reads_2:
            if not self.reads_2.match_paired(reads=reads):
                return False

            if self.reads_2.file_type == 'CASAVA':
                if reads.read == 'R1':
                    self.reads_1 = reads
                    reads.weak_reference_paired_reads = weakref.ReferenceType(self)
                    return True
                elif reads.read == 'R2':
                    raise Exception(
                        'PairedReads reads_2 has already been defined.\n'
                        '  reads_2: {!r}\n'
                        '  reads:   {!r}'.format(self.reads_2.file_path, reads.file_path))
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads.read))
            else:
                # Other file types go here ...
                warnings.warn(
                    'Method not implemented for file types other than CASAVA.',
                    UserWarning)

        return False

    def match(self, paired_reads):
        """Match C{bsf.ngs.PairedReads} objects.

        @param paired_reads: C{bsf.ngs.PairedReads}
        @type paired_reads: bsf.ngs.PairedReads
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
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

        if (self.reads_2 and paired_reads.reads_2) and not self.reads_2.match(reads=paired_reads.reads_2):
            return False

        return True

    def get_name(self, full=False):
        """Get the name of a C{bsf.ngs.PairedReads} object.

        For the C{bsf.ngs.Reads.file_type} I{CASAVA} the name is a concatenation of the
        C{bsf.ngs.Reads.name}, C{bsf.ngs.Reads.barcode} and C{bsf.ngs.Reads.lane} attributes,
        preferentially derived from the first C{bsf.ngs.Reads} object in the
        C{bsf.ngs.PairedReads} object.
        If the I{full} parameter is set, C{bsf.ngs.Reads.read} and C{bsf.ngs.Reads.chunk} are also added.
        @param full: Return the full name including read and chunk information
        @type full: bool
        @return: Name
        @rtype: str
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
    """The C{bsf.ngs.Sample} class represents a Next-Generation Sequencing sample.

    It consists of one or more C{bsf.ngs.PairedReads} objects as (biological or technical) replicates
    that result from the same flow cell.
    Attributes:
    @cvar default_name: Default key
    @type default_name: str
    @ivar paired_reads_list: Python C{list} of C{bsf.ngs.PairedReads} objects
    @type paired_reads_list: list[bsf.ngs.PairedReads]
    @ivar weak_reference_project: C{weakref.ReferenceType} pointing at a C{bsf.ngs.Project} object
    @type weak_reference_project: weakref.ReferenceType
    """

    default_name = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.ngs.Sample} object from a file path.

        For a I{file_type} I{CASAVA} the name is automatically populated,
        while C{bsf.ngs.PairedReads} objects are automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.ngs.Sample}
        @rtype: bsf.ngs.Sample
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
                file_mode = os.stat(file_path).st_mode
                if not stat.S_ISREG(file_mode):
                    continue
                sample.add_reads(reads=Reads.from_file_path(file_path=file_path, file_type=file_type))
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return sample

    @classmethod
    def from_samples(cls, sample1, sample2):
        """Create a merged C{bsf.ngs.Sample} from two C{bsf.ngs.Sample} objects.

        @param sample1: C{bsf.ngs.Sample}
        @type sample1: bsf.ngs.Sample
        @param sample2: C{bsf.ngs.Sample}
        @type sample2: bsf.ngs.Sample
        @return: C{bsf.ngs.Sample}
        @rtype: bsf.ngs.Sample
        """

        assert isinstance(sample1, Sample)
        assert isinstance(sample2, Sample)

        if sample1.name != sample2.name:
            warnings.warn(
                'Merged Sample objects {!r} and {!r} should have the same name.'.format(sample1.name, sample2.name),
                UserWarning)

        # A file_path does not make sense for merged Sample objects.

        sample = cls(file_type=sample1.file_type, name=sample1.name)

        # Merge the PairedReads objects from both Sample objects,
        # but check, if the PairedReads objects are not already there.
        # TODO: This method assumes that the PairedReads object are at the same address.
        # The test with 'in' will not work for distinct PairedReads objects that point to the same files.

        for paired_reads in sample1.paired_reads_list:
            if paired_reads not in sample.paired_reads_list:
                sample.add_paired_reads(paired_reads=paired_reads)

        for paired_reads in sample2.paired_reads_list:
            if paired_reads not in sample.paired_reads_list:
                sample.add_paired_reads(paired_reads=paired_reads)

        return sample

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            paired_reads_list=None,
            weak_reference_project=None):
        """Initialise a C{bsf.ngs.Sample} object.

        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @param paired_reads_list: Python C{list} of C{bsf.ngs.PairedReads} objects
        @type paired_reads_list: list[bsf.ngs.PairedReads]
        @param weak_reference_project: C{weakref.ReferenceType} pointing at a C{bsf.ngs.Project} object
        @type weak_reference_project: weakref.ReferenceType
        @return:
        @rtype:
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
                paired_reads_object.weak_reference_sample = weakref.ReferenceType(self)

        self.weak_reference_project = weak_reference_project  # Can be None.

        return

    def trace(self, level):
        """Trace a C{bsf.ngs.Sample} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_project: {!r}\n'.format(indent, self.weak_reference_project)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  annotation_dict:\n'.format(indent, self.annotation_dict)
        for annotation_key in self.annotation_dict.keys():
            output += '{}    {!r} {!r}\n'.format(indent, annotation_key, self.annotation_dict[annotation_key])
        output += '{}  paired_reads:\n'.format(indent)
        for paired_reads in self.paired_reads_list:
            output += paired_reads.trace(level + 1)

        return output

    def match(self, sample):
        """Match C{bsf.ngs.Sample} objects.

        @param sample: C{bsf.ngs.Sample}
        @type sample: bsf.ngs.Sample
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
        """

        # When Sample objects get merged, the file_path does no longer make sense.
        # Hence cope with empty file paths here ...

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
        """Add a C{bsf.ngs.PairedReads} object.

        This method checks, whether a matching C{bsf.ngs.PairedReads} object is already present in this
        C{bsf.ngs.Sample} object.
        @param paired_reads: C{bsf.ngs.PairedReads}
        @type paired_reads: bsf.ngs.PairedReads
        @return:
        @rtype:
        """

        assert isinstance(paired_reads, PairedReads)

        # Iterate through the Python list of PairedReads objects.
        # The PairedReads object must not match.

        for old_paired_reads in self.paired_reads_list:
            assert isinstance(old_paired_reads, PairedReads)
            if old_paired_reads.match(paired_reads=paired_reads):
                break
        else:
            # None of the existing PairedReads objects has matched,
            # so add this one to the Sample object.
            self.paired_reads_list.append(paired_reads)
            paired_reads.weak_reference_sample = weakref.ReferenceType(self)

        return

    def add_reads(self, reads):
        """Add a C{bsf.ngs.Reads} object.

        @param reads: C{bsf.ngs.Reads}
        @type reads: bsf.ngs.Reads
        @return:
        @rtype:
        """

        assert isinstance(reads, Reads)

        # Iterate through the Python list of PairedReads objects
        # The sample name must be identical.
        # The chunk must be identical
        # The read must fit

        for paired_reads in self.paired_reads_list:
            assert isinstance(paired_reads, PairedReads)
            if paired_reads.add_reads(reads=reads):
                break
        else:
            # The Reads object could not be added.
            # Create a new PairedReads object initialised
            # with the Reads object and add it to the Sample.
            self.add_paired_reads(paired_reads=PairedReads(reads_1=reads))

        return

    def get_all_paired_reads(self, replicate_grouping, exclude=False, full=False):
        """Get all C{bsf.ngs.PairedReads} of a C{bsf.ngs.Sample} grouped or un-grouped.

        A C{bsf.ngs.Sample} can hold several C{bsf.ngs.PairedReads} (i.e. biological or technical
        replicates) that have been sequenced on different lanes of the same flow cell.
        The C{bsf.ngs.PairedReads} will therefore differ in I{name}, I{barcode} or I{lane} information.
        Depending on the I{replicate_grouping} parameter they can be returned as a group or separately.
        C{bsf.ngs.PairedReads} that share the I{name}, I{barcode} and I{lane}, but differ in I{chunk} number
        are always grouped together.
        @param replicate_grouping: Group all C{bsf.ngs.PairedReads} of a C{bsf.ngs.Sample} or
            list them individually
        @type replicate_grouping: bool
        @param exclude: Exclude on the basis of C{bsf.ngs.PairedReads.exclude}
        @type exclude: bool
        @param full: Return the full name including read and chunk information
        @type full: bool
        @return: Python C{dict} of Python C{str} (sensible replicate name) key and
            Python C{list} of C{bsf.ngs.PairedReads} value data
        @rtype: dict[str, list[bsf.ngs.PairedReads]]
        """

        paired_reads_dict = dict()

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


class Project(NextGenerationBase):
    """The C{bsf.ngs.Project} class represents a Next-Generation Sequencing Project
    consisting of one or more C{bsf.ngs.Sample} objects.

    Attributes:

    @cvar default_name: Default key
    @type default_name: str
    @ivar sample_dict: Python C{dict} of C{bsf.ngs.Sample.name} key objects and C{bsf.ngs.Sample} value objects
    @type sample_dict: dict[bsf.ngs.Sample.name, bsf.ngs.Sample]
    @ivar weak_reference_prf: C{weakref.ReferenceType} pointing at a C{bsf.ngs.ProcessedRunFolder} object
    @type weak_reference_prf: weakref.ReferenceType
    """

    default_name = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.ngs.Project} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.ngs.Project}
        @rtype: bsf.ngs.Project
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
                file_mode = os.stat(file_path).st_mode
                re_match = re_pattern.search(string=file_name)
                if stat.S_ISDIR(file_mode) and re_match is not None:
                    if re_match.group(1) in project.sample_dict:
                        raise Exception(
                            'Sample with name {!r} already exists.'.format(re_match.group(1)))
                    else:
                        project.add_sample(sample=Sample.from_file_path(file_path=file_path, file_type=file_type))
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return project

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            sample_dict=None,
            weak_reference_prf=None):
        """Initialise a C{bsf.ngs.Project} object.

        For a I{file_type} I{CASAVA} the name is automatically populated,
        while C{bsf.ngs.Sample} objects are automatically discovered.
        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @param sample_dict: Python C{dict} of C{bsf.ngs.Sample.name} key objects and C{bsf.ngs.Sample} value objects
        @type sample_dict: dict[bsf.ngs.Sample.name, bsf.ngs.Sample]
        @param weak_reference_prf: C{weakref.ReferenceType} pointing at a C{bsf.ngs.ProcessedRunFolder} object
        @type weak_reference_prf: weakref.ReferenceType
        @raise Exception: If C{bsf.ngs.Sample.name} values are not unique for I{file_type} I{CASAVA}
        @return:
        @rtype:
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
            for sample in self.sample_dict.itervalues():
                assert isinstance(sample, Sample)
                sample.weak_reference_project = weakref.ReferenceType(self)

        self.weak_reference_prf = weak_reference_prf  # Can be None.

        return

    def trace(self, level):
        """Trace a C{bsf.ngs.Project} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_prf: {!r}\n'.format(indent, self.weak_reference_prf)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  annotation_dict:\n'.format(indent, self.annotation_dict)
        for annotation_key in self.annotation_dict.keys():
            output += '{}    {!r} {!r}\n'.format(indent, annotation_key, self.annotation_dict[annotation_key])
        output += '{}  sample_dict:\n'.format(indent)
        for sample_name in self.sample_dict.keys():
            output += self.sample_dict[sample_name].trace(level + 1)

        return output

    def add_sample(self, sample):
        """Add a C{bsf.ngs.Sample} object to a C{bsf.ngs.Project} object and set
        a weak reference in to C{bsf.ngs.Sample} object back to the C{bsf.ngs.Project} object.

        @param sample: C{bsf.ngs.Sample}
        @type sample: bsf.ngs.Sample
        @return: C{bsf.ngs.Sample}
        @rtype: bsf.ngs.Sample
        """

        assert isinstance(sample, Sample)
        # Delete an eventual bsf.ngs.Sample stored under the same bsf.ngs.Sample.name to remove the weak reference.
        self.del_sample(name=sample.name)
        self.sample_dict[sample.name] = sample
        sample.weak_reference_project = weakref.ReferenceType(self)

        return sample

    def del_sample(self, name):
        """Delete a C{bsf.ngs.Sample} object from a C{bsf.ngs.Project} object and
        clear the weak reference, if it points back at the C{bsf.ngs.Project} object.

        @param name: C{bsf.ngs.Sample.name}
        @type name: str
        @return: C{bsf.ngs.Sample}
        @rtype: bsf.ngs.Sample
        """

        if name in self.sample_dict:
            sample = self.sample_dict[name]
            del self.sample_dict[name]
            if sample.weak_reference_project() is self:
                sample.weak_reference_project = None
            return sample
        else:
            return

    def get_all_samples(self):
        """Get an ordered Python C{list} of C{bsf.ngs.Sample} objects.

        @return: Python C{list} of C{bsf.ngs.Sample} objects
        @rtype: list[bsf.ngs.Sample]
        """

        sample_list = list()

        sample_name_list = self.sample_dict.keys()
        sample_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for sample_name in sample_name_list:
            sample_list.append(self.sample_dict[sample_name])

        return sample_list


class ProcessedRunFolder(NextGenerationBase):
    """The C{bsf.ngs.ProcessedRunFolder} class represents an Illumina Run Folder after processing with CASAVA.

    Attributes:
    @cvar default_name: Default key
    @type default_name: str
    @ivar prefix: Prefix
    @type prefix: str
    @ivar flow_cell: Flow cell identifier
    @type flow_cell: str
    @ivar version: Version number
    @type version: str
    @ivar project_dict: Python C{dict} of C{bsf.ngs.Project.name} key objects and C{bsf.ngs.Project} value objects
    @type project_dict: dict[bsf.ngs.Project.name, bsf.ngs.Project]
    @ivar weak_reference_collection: C{weakref.ReferenceType} pointing at a C{bsf.ngs.Collection} object
    @type weak_reference_collection: weakref.ReferenceType
    """

    default_name = 'Default'

    @staticmethod
    def guess_file_type(file_path):
        """Guess the I{file_type} of a C{bsf.ngs.ProcessedRunFolder} on the basis of the I{file_path}.

        CASAVA C{bsf.ngs.ProcessedRunFolder} objects obey a I{Prefix_FCID_CASAVA182}
        schema. The following prefixes are currently in use:
            - BSF_ Biomedical Sequencing Facility
            - NGS_ Kaan Boztug group
            - MUW_ Medical University Vienna
            - SET_ Robert Kralovics group
        @param file_path: File path
        @type file_path: str | unicode
        @return: File type (i.e. I{CASAVA} or I{External})
        @rtype: str
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        component_list = file_name.split('_')

        if component_list[-1].startswith('CASAVA'):
            return 'CASAVA'
        else:
            return 'External'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.ngs.ProcessedRunFolder} object from a file path.

        For the I{file_type} I{CASAVA}, the C{bsf.ngs.ProcessedRunFolder.name}, C{bsf.ngs.ProcessedRunFolder.prefix},
        C{bsf.ngs.ProcessedRunFolder.flow_cell} and C{bsf.ngs.ProcessedRunFolder.version}
        attributes can be automatically parsed from the I{file_path}, while
        C{bsf.ngs.Project} objects can be automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.ngs.ProcessedRunFolder}
        @rtype: bsf.ngs.ProcessedRunFolder
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
                file_mode = os.stat(file_path).st_mode
                re_match = re_pattern.search(string=file_name)
                if stat.S_ISDIR(file_mode) and re_match is not None:
                    if re_match.group(1) in prf.project_dict:
                        raise Exception(
                            'Project with name {!r} already exists.'.format(re_match.group(1)))
                    else:
                        prf.add_project(project=Project.from_file_path(file_path=file_path, file_type=file_type))
        elif file_type == 'External':
            # Create a new, minimal ProcessedRunFolder.
            prf = cls(file_path=file_path, file_type=file_type, name=file_name)
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

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
        """Initialise a C{bsf.ngs.ProcessedRunFolder} object.

        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @param prefix: Prefix
        @type prefix: str
        @param flow_cell: Flow cell identifier
        @type flow_cell: str
        @param version: Version
        @type version: str
        @param project_dict: Python C{dict} of C{bsf.ngs.Project.name} key objects and C{bsf.ngs.Project} value objects
        @type project_dict: dict[bsf.ngs.Project.name, bsf.ngs.Project]
        @param weak_reference_collection: C{weakref.ReferenceType} pointing at a C{bsf.ngs.Collection} object
        @type weak_reference_collection: weakref.ReferenceType
        @return:
        @rtype:
        @raise Exception: If C{bsf.ngs.Project.name} values are not unique for file_type I{CASAVA}
        """

        super(ProcessedRunFolder, self).__init__(
            name=name,
            file_path=file_path,
            file_type=file_type,
            annotation_dict=annotation_dict
        )

        if prefix is None:
            self.prefix = str()
        else:
            self.prefix = prefix

        if flow_cell is None:
            self.flow_cell = str()
        else:
            self.flow_cell = flow_cell

        if version is None:
            self.version = str()
        else:
            self.version = version

        if project_dict is None:
            self.project_dict = dict()
        else:
            self.project_dict = project_dict
            # Set this bsf.ngs.ProcessedRunFolder as weak reference for each bsf.ngs.Project in the dict.
            for project in self.project_dict.itervalues():
                assert isinstance(project, Project)
                project.weak_reference_prf = weakref.ReferenceType(self)

        self.weak_reference_collection = weak_reference_collection  # Can be None.

        return

    def trace(self, level):
        """Trace a C{bsf.ngs.ProcessedRunFolder} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_collection: {!r}\n'.format(indent, self.weak_reference_collection)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  prefix:    {!r}\n'.format(indent, self.prefix)
        output += '{}  flow_cell: {!r}\n'.format(indent, self.flow_cell)
        output += '{}  version:   {!r}\n'.format(indent, self.version)
        output += '{}  annotation_dict:\n'.format(indent, self.annotation_dict)
        for annotation_key in self.annotation_dict.keys():
            output += '{}    {!r} {!r}\n'.format(indent, annotation_key, self.annotation_dict[annotation_key])
        output += '{}  project_dict:\n'.format(indent)
        for project_name in self.project_dict.keys():
            output += self.project_dict[project_name].trace(level + 1)

        return output

    def add_project(self, project):
        """Add a C{bsf.ngs.Project} object to the C{bsf.ngs.ProcessedRunFolder} object and set
        a weak reference in the C{bsf.ngs.Project} object to the C{bsf.ngs.ProcessedRunFolder} object.

        @param project: C{bsf.ngs.Project}
        @type project: bsf.ngs.Project
        @return: C{bsf.ngs.Project}
        @rtype: bsf.ngs.Project
        """

        assert isinstance(project, Project)
        # Delete an eventual bsf.ngs.Project that is stored under the same bsf.ngs.Project.name.
        self.del_project(name=project.name)
        self.project_dict[project.name] = project
        project.weak_reference_prf = weakref.ReferenceType(self)

        return project

    def del_project(self, name):
        """Delete a C{bsf.ngs.Project} object from am C{bsf.ngs.ProcessedRunFolder} object and
        clear the weak reference, if it points back at the C{bsf.ngs.ProcessedRunFolder} object.

        @param name: C{bsf.ngs.Project.name}
        @type name: str
        @return: C{bsf.ngs.Project}
        @rtype: bsf.ngs.Project
        """

        if name in self.project_dict:
            project = self.project_dict[name]
            del self.project_dict[name]
            if project.weak_reference_prf() is self:
                project.weak_reference_prf = None
            return project
        else:
            return

    def get_all_projects(self):
        """Get an ordered Python C{list} of C{bsf.ngs.Project} objects.

        @return: A Python C{list} of C{bsf.ngs.Project} objects
        @rtype: list[bsf.ngs.Project]
        """

        project_list = list()

        project_name_list = self.project_dict.keys()
        project_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for project_name in project_name_list:
            project_list.append(self.project_dict[project_name])

        return project_list


class Collection(NextGenerationBase):
    """The C{bsf.ngs.Collection} class represents a collection of
    one or more C{bsf.ngs.ProcessedRunFolder} objects.

    Attributes:
    @cvar default_name: Default key
    @type default_name: str
    @ivar processed_run_folder_dict: Python C{dict} of C{bsf.ngs.ProcessedRunFolder.name} key objects and
        C{bsf.ngs.ProcessedRunFolder} value objects
    @type processed_run_folder_dict: dict[bsf.ngs.ProcessedRunFolder.name, bsf.ngs.ProcessedRunFolder]
    @ivar sample_group_dict: Python C{dict} of Python C{str} (group) and C{bsf.ngs.Sample} objects
    @type sample_group_dict: dict[str, list[bsf.ngs.Sample]]
    """

    default_name = 'Default'

    @classmethod
    def from_sas_path(cls, file_path, file_type, name, sas_path, sas_prefix=None):
        """Construct a C{bsf.ngs.Collection} from a C{bsf.ngs.SampleAnnotationSheet} file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param sas_path: C{bsf.ngs.SampleAnnotationSheet} file path
        @type sas_path: str | unicode
        @param sas_prefix: Optional column header prefix
            (e.g. '[Control ]Sample', '[Treatment ]Sample', ...)
        @type sas_prefix: str
        @return: C{bsf.ngs.Collection}
        @rtype: bsf.ngs.Collection
        """

        return cls.from_sas(
            file_path=file_path,
            file_type=file_type,
            name=name,
            sas=SampleAnnotationSheet.from_file_path(file_path=sas_path),
            sas_prefix=sas_prefix)

    @classmethod
    def from_sas(cls, file_path, file_type, name, sas, sas_prefix=None):
        """Construct a C{bsf.ngs.Collection} from a C{bsf.ngs.SampleAnnotationSheet}.

        This method creates C{bsf.ngs.Reads}, C{bsf.ngs.PairedReads}, C{bsf.ngs.Sample}, C{bsf.ngs.Project} and
        C{bsf.ngs.ProcessedRunFolder} objects from a C{bsf.ngs.SampleAnnotationSheet} row Python C{dict} object.
        If object-specific keys or their corresponding values are not
        available from the Python row dict, new objects will be created
        with a default key.
        Sample Annotation Sheet format:
            - FileType: Data object I{file_type} (i.e. I{CASAVA} or I{External}), defaults to I{External}.
                The FileType I{CASAVA} allows for auto-discovery of C{bsf.ngs.ProcessedRunFolder} objects.
            - ProcessedRunFolder Name: C{bsf.ngs.ProcessedRunFolder} file_path, can be automatically registered.
            - Project Name: C{bsf.ngs.Project} name
            - Sample Name: C{bsf.ngs.Sample} name
            - Reads1 File: C{bsf.ngs.Reads.file_name} instance variable.
                Subjected to C{os.path.expanduser} and C{os.path.expandvars}.
                If still relative, the C{bsf.ngs.Collection.file_path} is prepended.
            - Reads1 Name: C{bsf.ngs.Reads.name} instance variable
            - Reads2 File: Same as Reads1 File
            - Reads2 Name: Same as Reads1 Name
            - Group: C{bsf.ngs.Sample} objects can be grouped for further analysis in
                e.g. RNA-Seq or ChIP-Seq experiments.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param sas: C{bsf.ngs.SampleAnnotationSheet}
        @type sas: bsf.ngs.SampleAnnotationSheet
        @param sas_prefix: Optional column header prefix
            (e.g. '[Control ]Sample', '[Treatment ]Sample', '[Point N ]Sample', ...)
        @type sas_prefix: str
        @return: C{bsf.ngs.Collection}
        @rtype: bsf.ngs.Collection
        """
        assert isinstance(sas, SampleAnnotationSheet)

        collection = cls(file_path=file_path, file_type=file_type, name=name)

        if sas_prefix is None:
            sas_prefix = str()

        prf = None
        project = None
        sample = None
        paired_reads = None

        for row_dict in sas.row_dicts:
            key_list = row_dict.keys()
            file_type = collection._process_file_type(row_dict=row_dict, key_list=key_list, prefix=sas_prefix)

            prf = collection._process_processed_run_folder(
                prf=prf,
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix,
                file_type=file_type)

            project = collection._process_project(
                project=project,
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix,
                file_type=file_type,
                prf=prf)

            sample = collection._process_sample(
                sample=sample,
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix,
                file_type=file_type,
                project=project)

            paired_reads = collection._process_paired_reads(
                paired_reads=paired_reads,
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix,
                file_type=file_type,
                sample=sample)

            # Optionally group the Sample objects.

            key = '{} Group'.format(sas_prefix).lstrip()

            if key in key_list:
                key_list.remove(key)

            if key in row_dict and row_dict[key]:
                value = row_dict[key]
                if value not in collection.sample_group_dict:
                    collection.sample_group_dict[value] = list()
                sample_list = collection.sample_group_dict[value]
                if sample not in sample_list:
                    sample_list.append(sample)

            if len(key_list):
                warnings.warn("Unexpected keys in sample annotation sheet: {!r}".format(key_list), UserWarning)

        # Quench empty default objects that are a consequence of empty lines in the sample annotation sheet.

        for prf_name in collection.processed_run_folder_dict.keys():
            prf = collection.processed_run_folder_dict[prf_name]
            assert isinstance(prf, ProcessedRunFolder)
            for project_name in prf.project_dict.keys():
                project = prf.project_dict[project_name]
                assert isinstance(project, Project)
                for sample_name in project.sample_dict.keys():
                    sample = project.sample_dict[sample_name]
                    assert isinstance(sample, Sample)
                    if sample.name == Sample.default_name and not len(sample.paired_reads_list):
                        project.del_sample(name=sample.name)
                if project.name == Project.default_name and not len(project.sample_dict):
                    prf.del_project(name=project.name)
            if prf.name == ProcessedRunFolder.default_name and not prf.project_dict:
                collection.del_processed_run_folder(name=prf.name)

        return collection

    def __init__(
            self,
            name=None,
            file_path=None,
            file_type=None,
            annotation_dict=None,
            processed_run_folder_dict=None,
            sample_group_dict=None):
        """Initialise a C{bsf.ngs.Collection} object.

        @param name: Name
        @type name: str
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param annotation_dict: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation_dict: dict[str, list[str]]
        @param processed_run_folder_dict: Python C{dict} of C{bsf.ngs.ProcessedRunFolder.name} key objects and
            C{bsf.ngs.ProcessedRunFolder} value objects
        @type processed_run_folder_dict: dict[bsf.ngs.ProcessedRunFolder.name, bsf.ngs.ProcessedRunFolder]
        @param sample_group_dict: Python C{dict} of Python C{str} (group name) key objects and
            second-level Python C{dict} value objects of C{bsf.ngs.Sample.name} key objects and
            C{bsf.ngs.Sample} value objects
        @type sample_group_dict: dict[str, list[bsf.ngs.Sample]]
        @return:
        @rtype:
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
            for prf in self.processed_run_folder_dict.itervalues():
                prf.weak_reference_collection = weakref.ReferenceType(self)

        if sample_group_dict is None:
            self.sample_group_dict = dict()
        else:
            self.sample_group_dict = sample_group_dict

        return

    def trace(self, level):
        """Trace a C{bsf.ngs.Collection} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  annotation_dict:\n'.format(indent, self.annotation_dict)
        for annotation_key in self.annotation_dict.keys():
            output += '{}    {!r} {!r}\n'.format(indent, annotation_key, self.annotation_dict[annotation_key])
        output += '{}  processed_run_folder_dict:\n'.format(indent)
        for prf_name in self.processed_run_folder_dict.keys():
            output += self.processed_run_folder_dict[prf_name].trace(level + 1)
        output += '{}  sample_group_dict:\n'.format(indent)
        sample_group_name_list = self.sample_group_dict.keys()
        sample_group_name_list.sort(cmp=lambda x, y: cmp(x, y))
        for sample_group_name in sample_group_name_list:
            output += '{}    group: {!r}\n'.format(indent, sample_group_name)
            # List all Sample objects of this Python list object.
            for sample in self.sample_group_dict[sample_group_name]:
                output += '{}      Sample name: {!r} file_path: {!r}\n'.format(indent, sample.name, sample.file_path)

        return output

    def add_processed_run_folder(self, prf):
        """Add a C{bsf.ngs.ProcessedRunFolder} object to the C{bsf.ngs.Collection} object and
        set a weak reference in the C{bsf.ngs.ProcessedRunFolder} back to the C{bsf.ngs.Collection}.

        @param prf: C{bsf.ngs.ProcessedRunFolder}
        @type prf: bsf.ngs.ProcessedRunFolder
        @return: C{bsf.ngs.ProcessedRunFolder}
        @rtype: bsf.ngs.ProcessedRunFolder
        """

        assert isinstance(prf, ProcessedRunFolder)
        # Delete an eventual bsf.ngs.ProcessedRunFolder that may exist under the same bsf.ngs.ProcessedRunFolder.name.
        self.del_processed_run_folder(name=prf.name)
        self.processed_run_folder_dict[prf.name] = prf
        prf.weak_reference_collection = weakref.ReferenceType(self)

        return prf

    def del_processed_run_folder(self, name):
        """Delete a C{bsf.ngs.ProcessedRunFolder} object from a C{bsf.ngs.Collection} object and
        clear the weak reference, if it points back at the C{bsf.ngs.Collection} object.

        @param name: C{bsf.ngs.ProcessedRunFolder.name}
        @type name: str
        @return: C{bsf.ngs.ProcessedRunFolder}
        @rtype: bsf.ngs.ProcessedRunFolder
        """

        if name in self.processed_run_folder_dict:
            prf = self.processed_run_folder_dict[name]
            del self.processed_run_folder_dict[name]
            if prf.weak_reference_collection() is self:
                prf.weak_reference_collection = None
            return prf
        else:
            return

    def get_processed_run_folder(self, file_path, file_type=None):
        """Get a C{bsf.ngs.ProcessedRunFolder} by file path.

        If the C{bsf.ngs.ProcessedRunFolder} does not exist in the C{bsf.ngs.Collection} object,
        it will be automatically discovered and added.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.ngs.ProcessedRunFolder}
        @rtype: bsf.ngs.ProcessedRunFolder
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_name in self.processed_run_folder_dict:
            return self.processed_run_folder_dict[file_name]
        else:
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if not os.path.isabs(file_path):
                file_path = os.path.join(self.file_path, file_path)

            return self.add_processed_run_folder(
                prf=ProcessedRunFolder.from_file_path(file_path=file_path, file_type=file_type))

    def get_all_processed_run_folders(self):
        """Get an ordered Python C{list} of C{bsf.ngs.ProcessedRunFolder} objects.

        @return: Python C{list} of C{bsf.ngs.ProcessedRunFolder} objects
        @rtype: list[bsf.ngs.ProcessedRunFolder]
        """

        processed_run_folder_list = list()

        processed_run_folder_name_list = self.processed_run_folder_dict.keys()
        processed_run_folder_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for processed_run_folder_name in processed_run_folder_name_list:
            processed_run_folder_list.append(self.processed_run_folder_dict[processed_run_folder_name])

        return processed_run_folder_list

    def get_all_projects(self):
        """Get an ordered Python C{list} of C{bsf.ngs.Project} objects.

        @return: Python C{list} of C{bsf.ngs.Project} objects
        @rtype: list[bsf.ngs.Project]
        """

        project_list = list()

        for processed_run_folder in self.get_all_processed_run_folders():
            project_list.extend(processed_run_folder.get_all_projects())

        return project_list

    def get_all_samples(self):
        """Get an ordered Python C{list} of C{bsf.ngs.Sample} objects.

        @return: Python C{list} of C{bsf.ngs.Sample} objects
        @rtype: list[bsf.ngs.Sample]
        """

        sample_list = list()

        for project in self.get_all_projects():
            sample_list.extend(project.get_all_samples())

        return sample_list

    @staticmethod
    def _process_file_type(row_dict, key_list, prefix):
        """Get file type information.

        A 'I{[Prefix] FileType}' key is optional, its value defaults to I{Automatic}.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] FileType', '[Treatment] FileType', '[Point N]  FileType', ...)
        @type prefix: str
        @return: File type
        @rtype: str
        """

        key = '{} File Type'.format(prefix).lstrip()

        # If present, delete the key from the key list.
        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            file_type = row_dict[key]
        else:
            file_type = 'Automatic'

        return file_type

    def _process_processed_run_folder(self, prf, row_dict, key_list, prefix, file_type):
        """Get or create a C{bsf.ngs.ProcessedRunFolder}.

        A 'I{[Prefix] ProcessedRunFolder Name}' key is optional, its value defaults to I{Default}.
        @param prf: Current C{bsf.ngs.ProcessedRunFolder} that may get replaced upon encountering a new
            'I{[Prefix] ProcessedRunFolder Name}' key
        @type prf: bsf.ngs.ProcessedRunFolder
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] ProcessedRunFolder', '[Treatment] ProcessedRunFolder',
            '[Point N] ProcessedRunFolder', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.ngs.ProcessedRunFolder}
        @rtype: bsf.ngs.ProcessedRunFolder
        """

        new_prf = None

        key = '{} ProcessedRunFolder Name'.format(prefix).lstrip()

        if key in row_dict:
            # The key exists ...
            key_list.remove(key)
            if row_dict[key]:
                # ... and has a meaningful value ...
                value = row_dict[key]
                if value in self.processed_run_folder_dict:
                    # ..., which exists in the dict of ProcessedRunFolder objects.
                    new_prf = self.processed_run_folder_dict[value]
                else:
                    # ..., which does not exist in the dict of ProcessedRunFolder objects.
                    # Try to automatically discover a ProcessedRunFolder.
                    new_prf = self.get_processed_run_folder(file_path=value, file_type=file_type)
            else:
                # ... and has no meaningful value ...
                if prf is not None:
                    # ..., but an old PRF exists.
                    new_prf = prf

        if new_prf is None:
            if ProcessedRunFolder.default_name in self.processed_run_folder_dict:
                new_prf = self.processed_run_folder_dict[ProcessedRunFolder.default_name]
            else:
                new_prf = self.add_processed_run_folder(
                    prf=ProcessedRunFolder(name=ProcessedRunFolder.default_name, file_type=file_type))

        new_prf.process_annotation(row_dict=row_dict, key_list=key_list, prefix=prefix)

        return new_prf

    @staticmethod
    def _process_project(project, row_dict, key_list, prefix, file_type, prf):
        """Get or create a C{bsf.ngs.Project}.

        A 'I{[Prefix] Project Name}' key is optional, its value defaults to I{Default}.
        @param project: Current C{bsf.ngs.Project} that may get replaced upon encountering a new
            'I{[Prefix] Project Name}' key
        @type project: bsf.ngs.Project
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Project', '[Treatment] Project', '[Point N] Project', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param prf: C{bsf.ngs.ProcessedRunFolder}
        @type prf: bsf.ngs.ProcessedRunFolder
        @return: C{bsf.ngs.Project}
        @rtype: bsf.ngs.Project
        """

        new_project = None

        key = '{} Project Name'.format(prefix).lstrip()

        if key in row_dict:
            # The key exists ...
            key_list.remove(key)  # Remove the 'Name' key.
            if row_dict[key]:
                # ... and has a meaningful value ...
                value = row_dict[key]
                if value in prf.project_dict:
                    # ..., which exists in the dict of Project objects.
                    new_project = prf.project_dict[value]
                else:
                    # ..., which does not exist in the dict of Project objects.
                    # Create a new Project.
                    new_project = prf.add_project(project=Project(name=value, file_type=file_type))
            else:
                # ... and has no meaningful value ...
                if project is not None:
                    # ..., but a current Project exists.
                    new_project = project

        if new_project is None:
            if Project.default_name in prf.project_dict:
                new_project = prf.project_dict[Project.default_name]
            else:
                new_project = prf.add_project(project=Project(name=Project.default_name, file_type=file_type))

        new_project.process_annotation(row_dict=row_dict, key_list=key_list, prefix=prefix)

        return new_project

    @staticmethod
    def _process_sample(sample, row_dict, key_list, prefix, file_type, project):
        """Get or create a C{bsf.ngs.Sample}.

        A 'I{[Prefix] Sample Name}' key is optional, its value defaults to I{Default}.
        @param sample: Current C{bsf.ngs.Sample} that may get replaced upon encountering a new
            'I{[Prefix] Sample Name}' key
        @type sample: bsf.ngs.Sample
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param project: C{bsf.ngs.Project}
        @type project: bsf.ngs.Project
        @return: C{bsf.ngs.Sample}
        @rtype: bsf.ngs.Sample
        """

        new_sample = None

        key = '{} Sample Name'.format(prefix).lstrip()

        if key in row_dict:
            key_list.remove(key)  # Remove the 'Name' key.
            # The key exists ...
            if row_dict[key]:
                # ... and has a meaningful value ...
                value = row_dict[key]
                if value in project.sample_dict:
                    # ..., which exists in the dict of Sample objects.
                    new_sample = project.sample_dict[value]
                else:
                    # ..., which does not exist in the dict of Sample objects.
                    # Create a new Sample.
                    new_sample = project.add_sample(sample=Sample(name=value, file_type=file_type))
            else:
                # ... and has no meaningful value ...
                if sample is not None:
                    # ..., but a current Sample exists.
                    new_sample = sample

        if new_sample is None:
            if Sample.default_name in project.sample_dict:
                new_sample = project.sample_dict[Sample.default_name]
            else:
                new_sample = project.add_sample(sample=Sample(name=Sample.default_name, file_type=file_type))

        new_sample.process_annotation(row_dict=row_dict, key_list=key_list, prefix=prefix)

        return new_sample

    def _process_reads(self, reads, row_dict, key_list, prefix, file_type, suffix):
        """Get or create a C{bsf.ngs.Reads} object.

        A 'I{[Prefix] Reads{suffix} Name}' key and 'I{[Prefix] Reads{suffix} File}' key are optional,
        in which case the default is a C{None} object.
        @param reads: Current C{bsf.ngs.Reads} that may get replaced upon encountering a new
            'I{[Prefix] ReadsN Name}' key
        @type reads: bsf.ngs.Reads
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] ReadsN', '[Treatment] ReadsN', '[Point N] ReadsN', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param suffix: The read suffix (i.e. I{1} or I{2})
        @type suffix: str
        @return: C{bsf.ngs.Reads}
        @rtype: bsf.ngs.Reads
        """

        if reads is None:
            reads = Reads()

        new_reads = Reads(file_type=file_type)

        # Pre-process the Reads.file_path instance variable.

        key = '{} Reads{} File'.format(prefix, suffix).lstrip()
        if key in key_list:
            key_list.remove(key)

        is_new_file_path = False
        if key in row_dict and row_dict[key]:
            new_reads.file_path = row_dict[key]
            new_reads.file_path = os.path.expanduser(new_reads.file_path)
            new_reads.file_path = os.path.expandvars(new_reads.file_path)
            if not os.path.isabs(new_reads.file_path):
                new_reads.file_path = os.path.join(self.file_path, new_reads.file_path)
            new_reads.file_path = os.path.normpath(new_reads.file_path)
            # Check for a non-matching, i.e. new "file_path" instance variable.
            if new_reads.file_path != reads.file_path:
                is_new_file_path = True

        # Pre-process the Reads.name instance variable.

        key = '{} Reads{} Name'.format(prefix, suffix).lstrip()
        if key in key_list:
            key_list.remove(key)

        is_new_name = False
        if key in row_dict and row_dict[key]:
            new_reads.name = row_dict[key]
            # Check for a non-matching i.e. new "name" instance variable.
            if new_reads.name != reads.name:
                is_new_name = True

        if is_new_file_path and is_new_name:
            # All is well. Just return the new Reads object.
            pass
        elif is_new_file_path:
            # A new file_path, but check the name.
            if new_reads.name:
                if reads.name and new_reads.name == reads.name:
                    # The old Reads object has the same Reads.name instance variable.
                    raise Exception("Encountered new Reads.file_path {} -> {}, but the same Reads.name {}.".format(
                        reads.file_path, new_reads.file_path, new_reads.name))
                else:
                    # Set the Reads.name instance variable.
                    reads.name = new_reads.name
                    new_reads = reads
        elif is_new_name:
            # A new name, but check the file_path.
            if new_reads.file_path:
                if reads.file_path and new_reads.file_path == reads.file_path:
                    # The old Reads object has the same Reads.file_path instance variable.
                    raise Exception("Encountered new Reads.name {} -> {}, but same reads.file_path {}.".format(
                        reads.name, new_reads.name, new_reads.file_path))
                else:
                    # Set the Reads.file_path instance variable.
                    reads.file_path = new_reads.file_path
                    new_reads = reads
        else:
            # Nothing new, just return the old Reads object.
            new_reads = reads

        return new_reads

    # Taken from ConfigParser.RawConfigParser.getboolean()

    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False}

    def _process_paired_reads(self, paired_reads, row_dict, key_list, prefix, file_type, sample):
        """Get or create a C{bsf.ngs.PairedReads} object.

        The 'I{[Prefix] PairedReads Exclude}' key is optional, in which case the default is Python C{bool} I{False}.
        The 'I{[Prefix] PairedReads Index 1}', 'I{[Prefix] PairedReads Index 2}' and 'I{[Prefix] PairedReads ReadGroup}'
        keys are optional, in which case the default is an empty Python C{str} object.
        @param paired_reads: C{PairedReads}
        @type paired_reads: PairedReads
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] PairedReads ReadGroup',
            '[Treatment] PairedReads ReadGroup',
            '[Point N] PairedReads ReadGroup', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @return: C{PairedReads}
        @rtype: PairedReads
        """

        def _is_new_reads(reads_1, reads_2):
            """Test whether a C{bsf.ngs.Reads} object is new.

            To test for object equality is not good enough here. A new C{bsf.ngs.Reads} object is encountered
            upon mismatches in the C{bsf.ngs.Reads.file_path} or C{bsf.ngs.Reads.name} instance variables.
            @param reads_1: C{bsf.ngs.Reads}
            @type reads_1: Reads
            @param reads_2: C{bsf.ngs.Reads}
            @type reads_2: Reads
            @return: True if C{Reads} is new
            @rtype: bool
            """
            if reads_1 is None and reads_2 is None:
                # Both Reads objects are not defined, so, nothing new.
                is_new_reads = False
            elif not (reads_1 or reads_2):
                is_new_reads = False
            elif reads_1 is None:
                # At least reads_2 must be defined.
                is_new_reads = True
            elif reads_2 is None:
                # At least reads_1 must be defined.
                is_new_reads = True
            elif reads_1.name != reads_2.name or reads_1.file_path != reads_2.file_path:
                is_new_reads = True
            else:
                is_new_reads = False

            return is_new_reads

        # PairedReads objects have no name instance variable. Thus, a new PairedReads object has to begin,
        # when a new Reads() object is encountered.

        if paired_reads is None:
            paired_reads = PairedReads()

        new_paired_reads = PairedReads(
            file_type=file_type,
            reads_1=self._process_reads(
                reads=paired_reads.reads_1,
                row_dict=row_dict,
                key_list=key_list,
                prefix=prefix,
                file_type=file_type,
                suffix='1'),
            reads_2=self._process_reads(
                reads=paired_reads.reads_2,
                row_dict=row_dict,
                key_list=key_list,
                prefix=prefix,
                file_type=file_type,
                suffix='2'))

        if _is_new_reads(new_paired_reads.reads_1, paired_reads.reads_1) or \
                _is_new_reads(new_paired_reads.reads_2, paired_reads.reads_2):
            sample.add_paired_reads(paired_reads=new_paired_reads)
        else:
            new_paired_reads = paired_reads

        key = '{} PairedReads Exclude'.format(prefix).lstrip()
        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            if row_dict[key].lower() not in Collection._boolean_states:
                raise ValueError('Value in field {!r} is not a boolean: {!r}'.format(key, row_dict[key]))
            new_paired_reads.exclude = Collection._boolean_states[row_dict[key].lower()]

        key = '{} PairedReads Index 1'.format(prefix).lstrip()
        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            new_paired_reads.index_1 = row_dict[key]

        key = '{} PairedReads Index 2'.format(prefix).lstrip()
        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            new_paired_reads.index_2 = row_dict[key]

        key = '{} PairedReads ReadGroup'.format(prefix).lstrip()
        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            new_paired_reads.read_group = row_dict[key]

        new_paired_reads.process_annotation(row_dict=row_dict, key_list=key_list, prefix=prefix)

        return new_paired_reads

    def get_sample_from_row_dict(self, row_dict, prefix=None):
        """Get a Sample from a C{bsf.ngs.SampleAnnotationSheet} row Python C{dict}.

        Look-up a hierarchy of C{bsf.ngs.ProcessedRunFolder}, C{bsf.ngs.Project} and C{bsf.ngs.Sample} objects
        based on a C{bsf.ngs.SampleAnnotationSheet} row dictionary.
        C{bsf.ngs.ProcessedRunFolder} objects of file type I{CASAVA} can be
        automatically discovered and registered.
        Return the corresponding C{bsf.ngs.Sample}.
        @param row_dict: C{bsf.ngs.SampleAnnotationSheet} row Python C{dict}
        @type row_dict: dict[str, str | unicode]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample')
        @type prefix: str
        @return: C{bsf.ngs.Sample}
        @rtype: bsf.ngs.Sample
        """

        # NOTE: For the moment, the row_dict has to include keys for 'ProcessedRunFolder',
        # 'Project' and 'Sample'. Maybe, this method could search for a Sample name in
        # the Collection object. However, the problem is that a Collection contains
        # complete, auto-registered RunFolder objects. Thus, Sample names (typically 1, 2, 3, ...)
        # may clash. Therefore, it is probably best to stick to the three fields for
        # unambiguous sample resolution.

        if not prefix:
            prefix = str()

        key = '{} ProcessedRunFolder Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = ProcessedRunFolder.default_name

        # The Collection.get_processed_run_folder method can automatically register
        # ProcessedRunFolder objects of file type 'CASAVA'.

        prf = self.get_processed_run_folder(file_path=value)

        key = '{} Project Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Project.default_name

        project = prf.project_dict[value]

        key = '{} Sample Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Sample.default_name

        sample = project.sample_dict[value]

        return sample

    def get_samples_from_row_dict(self, row_dict, prefix=None):
        """Get a Python C{list} of C{bsf.ngs.Sample} objects and a Python C{str} of the Group column value
        as a Python C{tuple} from a C{bsf.ngs.SampleAnnotationSheet} row Python C{dict}.

        @param row_dict: Comparison CSV file row Python C{dict}
        @type row_dict: dict[str, str | unicode]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        @return: Python C{tuple} of Python C{str} of '[Prefix] Group' column value and
            Python C{list} of C{bsf.ngs.Sample} objects
        @rtype: (str, list[bsf.ngs.Sample])
        """

        sample_list = list()
        value = str()

        if not prefix:
            prefix = str()

        key = '{} Group'.format(prefix).lstrip()

        # Does the key exist and is its value defined?
        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Extend the Python list with all Sample objects of this group.
            if value in self.sample_group_dict:
                sample_list.extend(self.sample_group_dict[value])

        key = '{} Sample Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Append the Sample object, if defined.

            sample = self.get_sample_from_row_dict(row_dict=row_dict, prefix=prefix)

            if sample:
                sample_list.append(sample)

        return value, sample_list

    def to_sas(self, file_path=None, name=None):
        """Convert a C{bsf.ngs.Collection} into a SampleAnnotationSheet object.

        @param file_path: File path
        @type file_path: str | unicode
        @param name: Name
        @type name: str
        @return: SampleAnnotationSheet object
        @rtype: SampleAnnotationSheet
        """

        sas = SampleAnnotationSheet(file_path=file_path, name=name)

        # Scan the Collection and its contained objects for additional (annotation) field names.

        for prf in self.processed_run_folder_dict.itervalues():
            assert isinstance(prf, ProcessedRunFolder)
            if prf.annotation_dict is not None:
                for prf_annotation_key in prf.annotation_dict.iterkeys():
                    prf_annotation_field = ' '.join(('ProcessedRunFolder', prf_annotation_key))
                    if prf_annotation_field not in sas.field_names:
                        sas.field_names.append(prf_annotation_field)
            for project in prf.project_dict.itervalues():
                assert isinstance(project, Project)
                if project.annotation_dict is not None:
                    for project_annotation_key in project.annotation_dict.iterkeys():
                        project_annotation_field = ' '.join(('Project', project_annotation_key))
                        if project_annotation_field not in sas.field_names:
                            sas.field_names.append(project_annotation_field)
                for sample in project.sample_dict.itervalues():
                    assert isinstance(sample, Sample)
                    if sample.annotation_dict is not None:
                        for sample_annotation_key in sample.annotation_dict.iterkeys():
                            sample_annotation_field = ' '.join(('Sample', sample_annotation_key))
                            if sample_annotation_field not in sas.field_names:
                                sas.field_names.append(sample_annotation_field)
                    for paired_reads in sample.paired_reads_list:
                        assert isinstance(paired_reads, PairedReads)
                        if paired_reads.annotation_dict is not None:
                            for paired_reads_annotation_key in paired_reads.annotation_dict.iterkeys():
                                paired_reads_annotation_field = ' '.join(('PairedReads', paired_reads_annotation_key))
                                if paired_reads_annotation_field not in sas.field_names:
                                    sas.field_names.append(paired_reads_annotation_field)

        # At this stage all annotation keys should be added. Partition and sort the list of field names.

        field_names = list(sas.field_names)
        field_names.sort(cmp=lambda x, y: cmp(x, y))

        sas.field_names = list()
        sas.field_names += filter(lambda x: x.startswith('File Type'), field_names)
        field_names = filter(lambda x: not x.startswith('File Type'), field_names)
        sas.field_names += filter(lambda x: x.startswith('ProcessedRunFolder'), field_names)
        field_names = filter(lambda x: not x.startswith('ProcessedRunFolder'), field_names)
        sas.field_names += filter(lambda x: x.startswith('Project'), field_names)
        field_names = filter(lambda x: not x.startswith('Project'), field_names)
        sas.field_names += filter(lambda x: x.startswith('Sample'), field_names)
        field_names = filter(lambda x: not x.startswith('Sample'), field_names)
        sas.field_names += filter(lambda x: x.startswith('PairedReads'), field_names)
        field_names = filter(lambda x: not x.startswith('PairedReads'), field_names)
        sas.field_names += filter(lambda x: x.startswith('Reads1'), field_names)
        field_names = filter(lambda x: not x.startswith('Reads1'), field_names)
        sas.field_names += filter(lambda x: x.startswith('Reads2'), field_names)
        field_names = filter(lambda x: not x.startswith('Reads2'), field_names)
        sas.field_names += field_names

        # Finally, construct the SampleAnnotationSheet.

        prf_name_list = self.processed_run_folder_dict.keys()
        prf_name_list.sort(cmp=lambda x, y: cmp(x, y))
        for prf_name in prf_name_list:
            prf = self.processed_run_folder_dict[prf_name]
            assert isinstance(prf, ProcessedRunFolder)
            sas.row_dicts.append({
                'ProcessedRunFolder Name': prf.name,
            })
            if prf.annotation_dict is not None:
                prf_annotation_key_list = prf.annotation_dict.keys()
                prf_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                for prf_annotation_key in prf_annotation_key_list:
                    prf_annotation_list = prf.annotation_dict[prf_annotation_key]
                    prf_annotation_field = ' '.join(('ProcessedRunFolder', prf_annotation_key))
                    for annotation in prf_annotation_list:
                        sas.row_dicts.append({prf_annotation_field: annotation})
            project_name_list = prf.project_dict.keys()
            project_name_list.sort(cmp=lambda x, y: cmp(x, y))
            for project_name in project_name_list:
                project = prf.project_dict[project_name]
                assert isinstance(project, Project)
                sas.row_dicts.append({
                    'Project Name': project.name,
                })
                if project.annotation_dict is not None:
                    project_annotation_key_list = project.annotation_dict.keys()
                    project_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                    for project_annotation_key in project_annotation_key_list:
                        project_annotation_list = project.annotation_dict[project_annotation_key]
                        project_annotation_field = ' '.join(('Project', project_annotation_key))
                        for annotation in project_annotation_list:
                            sas.row_dicts.append({project_annotation_field: annotation})
                sample_name_list = project.sample_dict.keys()
                sample_name_list.sort(cmp=lambda x, y: cmp(x, y))
                for sample_name in sample_name_list:
                    sample = project.sample_dict[sample_name]
                    assert isinstance(sample, Sample)
                    sas.row_dicts.append({
                        'Sample Name': sample.name
                    })
                    if sample.annotation_dict is not None:
                        sample_annotation_key_list = sample.annotation_dict.keys()
                        sample_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                        for sample_annotation_key in sample_annotation_key_list:
                            sample_annotation_list = sample.annotation_dict[sample_annotation_key]
                            sample_annotation_field = ' '.join(('Sample', sample_annotation_key))
                            for annotation in sample_annotation_list:
                                sas.row_dicts.append({sample_annotation_field: annotation})
                    for paired_reads in sample.paired_reads_list:
                        assert isinstance(paired_reads, PairedReads)
                        row_dict = {
                            'PairedReads Exclude': '{}'.format(paired_reads.exclude).lower(),
                            'PairedReads Index 1': paired_reads.index_1,
                            'PairedReads Index 2': paired_reads.index_2,
                            'PairedReads ReadGroup': paired_reads.read_group,
                        }
                        if paired_reads.reads_1 is not None:
                            row_dict['Reads1 Name'] = paired_reads.reads_1.name
                            row_dict['Reads1 File'] = paired_reads.reads_1.file_path
                        if paired_reads.reads_2 is not None:
                            row_dict['Reads2 Name'] = paired_reads.reads_2.name
                            row_dict['Reads2 File'] = paired_reads.reads_2.file_path
                        sas.row_dicts.append(row_dict)
                        if paired_reads.annotation_dict is not None:
                            paired_reads_annotation_key_list = paired_reads.annotation_dict.keys()
                            paired_reads_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                            for paired_reads_annotation_key in paired_reads_annotation_key_list:
                                paired_reads_annotation_list = paired_reads.annotation_dict[paired_reads_annotation_key]
                                paired_reads_annotation_field = ' '.join(('PairedReads', paired_reads_annotation_key))
                                for annotation in paired_reads_annotation_list:
                                    sas.row_dicts.append({paired_reads_annotation_field: annotation})

        return sas

    def to_sas_path(self, file_path=None, name=None):
        """Write as a SampleAnnotationSheet to a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param name: Name
        @type name: str
        @return:
        @rtype:
        """

        sas = self.to_sas(file_path=file_path, name=name)
        sas.to_file_path()

        return


class SampleGroup(object):
    """The C{bsf.ngs.SampleGroup} class represents a group of C{bsf.ngs.Sample} objects.

    The grouping is usually defined in a sample annotation sheet.
    Attributes:
    @ivar name: Name
    @type name: str
    @ivar sample_list: Python C{list} of C{bsf.ngs.Sample} objects
    @type sample_list: list[bsf.ngs.Sample]
    """

    # TODO: The SampleGroup class is currently not in use.
    # Sample and PairedReads objects from different ProcessRunFolder objects
    # (i.e. flow cells) could bear the same name, leading to problems with SGE job names.
    # This would need further re-thinking.

    def __init__(
            self,
            name=None,
            sample_list=None):
        """Initialise a C{bsf.ngs.SampleGroup} object.

        @param name: Name
        @type name: str
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @return:
        @rtype:
        """

        super(SampleGroup, self).__init__()

        if name is None:
            self.name = str()
        else:
            self.name = name

        if sample_list is None:
            self.sample_list = list()
        else:
            self.sample_list = sample_list

        return

    def add_sample(self, sample):
        """Add a C{bsf.ngs.Sample} object.

        @param sample: C{bsf.ngs.Sample}
        @type sample: bsf.ngs.Sample
        @return:
        @rtype:
        """

        assert isinstance(sample, Sample)
        if sample not in self.sample_list:
            self.sample_list.append(sample)

        return

    def get_all_paired_reads(self, replicate_grouping):
        """Get all C{bsf.ngs.PairedReads} objects of a C{bsf.ngs.SampleGroup}.

        For the moment, replicates are C{bsf.ngs.PairedReads} objects that do not share
        anything but the chunk number.
        @param replicate_grouping: Group all C{bsf.ngs.PairedReads} objects of a C{bsf.ngs.Sample} or
            list them individually
        @type replicate_grouping: bool
        @return: Python C{dict} of Python C{str} (replicate name) key and
            Python C{list} value objects of C{bsf.ngs.PairedReads} objects value data
        @rtype: dict[str, list[bsf.ngs.PairedReads]]
        """

        group_dict = dict()

        for sample in self.sample_list:
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping)
            for paired_reads_name in paired_reads_dict.keys():
                if paired_reads_name not in group_dict:
                    group_dict[paired_reads_name] = list()

                # Add PairedReads objects one-by-one and check if they are not already there.

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if paired_reads not in group_dict[paired_reads_name]:
                        group_dict[paired_reads_name].append(paired_reads)

        return group_dict


class SampleAnnotationSheet(AnnotationSheet):
    """The C{bsf.ngs.SampleAnnotationSheet} class represents a Comma-Separated Value (CSV) table of sample information
    after running the C{bsf.analyses.illumina_to_bam_tools.IlluminaToBamTools.BamIndexDecoder} C{bsf.Analysis}.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the C{csv.Dialect} class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list[str]
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict[str, list[function]]
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

    _test_methods = {
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
