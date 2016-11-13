"""bsf.data

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


class Reads(object):
    """The C{bsf.data.Reads} class represents a file of Next-Generation Sequencing (NGS) reads,
    such as a FASTQ or unmapped BAM file.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar barcode: Barcode used for sample multiplexing
    @type barcode: str
    @ivar lane: Lane number
    @type lane: str
    @ivar read: Read number (e.g. I{R1}, I{R2}, ...)
    @type read: str
    @ivar chunk: Chunk number (e.g. I{001}, I{002}, ...)
    @type chunk: str
    @ivar weak_reference_paired_reads: C{weakref.ReferenceType} pointing at a C{bsf.data.PairedReads} object
    @type weak_reference_paired_reads: weakref.ReferenceType
    """

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.data.Reads} object from a file path.

        For a I{file_type} I{CASAVA}, C{bsf.data.Reads.file_path} obeys a I{SampleName_Index_Lane_Read_Chunk} schema,
        so that C{bsf.data.Reads.name}, C{bsf.data.Reads.barcode}, C{bsf.data.Reads.lane}, C{bsf.data.Reads.read} and
        C{bsf.data.Reads.chunk} can be populated automatically.
        For I{file_type} I{External}, the attributes need populating manually.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.data.Reads} object
        @rtype: bsf.data.Reads
        """

        if file_type == 'CASAVA':
            file_path = os.path.normpath(file_path)
            file_name = os.path.basename(file_path)

            # CASAVA Reads obey a SampleName_Index_Lane_Read_Chunk schema.
            components = file_name.split('.')
            components[0] = components[0].split('_')

            name = '_'.join(components[0][:-4])
            barcode = components[0][-4]
            lane = components[0][-3]
            read = components[0][-2]
            chunk = components[0][-1]

            reads = cls(file_path=file_path, file_type=file_type, name=name,
                        barcode=barcode, lane=lane, read=read, chunk=chunk)
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return reads

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            barcode=None,
            lane=None,
            read=None,
            chunk=None,
            weak_reference_paired_reads=None):
        """Initialise a C{bsf.data.Reads} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param barcode: Barcode used for sample multiplexing
        @type barcode: str
        @param lane: Lane number
        @type lane: str
        @param read: Read number (e.g. I{R1}, I{R2}, ...)
        @type read: str
        @param chunk: Chunk number (e.g. I{001}, I{002}, ...)
        @type chunk:str
        @param weak_reference_paired_reads: C{weakref.ReferenceType} pointing at a C{bsf.data.PairedReads} object
        @type weak_reference_paired_reads: weakref.ReferenceType
        @return:
        @rtype:
        """

        super(Reads, self).__init__()

        if file_path is None:
            self.file_path = str()
        else:
            self.file_path = file_path

        if file_type is None:
            self.file_type = str()
        else:
            self.file_type = file_type

        if name is None:
            self.name = str()
        else:
            self.name = name

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

    def trace(self, level):
        """Trace a C{bsf.data.Reads} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_paired_reads: {!r}\n'.format(indent, self.weak_reference_paired_reads)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  barcode:   {!r}\n'.format(indent, self.barcode)
        output += '{}  lane:      {!r}\n'.format(indent, self.lane)
        output += '{}  read:      {!r}\n'.format(indent, self.read)
        output += '{}  chunk:     {!r}\n'.format(indent, self.chunk)

        return output

    def match(self, reads):
        """Match C{bsf.data.Reads} objects.

        Two C{bsf.data.Reads} objects are identical, if all their instance variables match.
        @param reads: Second C{bsf.data.Reads} object
        @type reads: bsf.data.Reads
        @return: True if both objects match, False otherwise
        @rtype: bool
        """

        assert isinstance(reads, Reads)

        # Quick test first - if the objects are identical, the rest has to match.

        if self == reads:
            return True

        if self.file_path != reads.file_path:
            return False

        if self.file_type != reads.file_type:
            return False

        if self.name != reads.name:
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
        """Match paired C{bsf.data.Reads} objects, by relaxing matching criteria.

        @param reads: Second C{bsf.data.Reads} object
        @type reads: bsf.data.Reads
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
        """

        assert isinstance(reads, Reads)

        if self.file_type == 'CASAVA':

            # All CASAVA Reads attributes need to be equal,
            # with the exception of the read (R1 or R2) and
            # the file_path, which also contains the reads
            # information.

            if self.file_type != reads.file_type:
                return False

            if self.name != reads.name:
                return False

            if self.barcode != reads.barcode:
                return False

            if self.lane != reads.lane:
                return False

            if self.chunk != reads.chunk:
                return False

            return True
        else:
            # It is difficult to match PairedReads outside of CASAVA conventions.
            warnings.warn(
                'Matching of paired Reads objects for file_type other than CASAVA not implemented yet.',
                UserWarning)

            return True


class PairedReads(object):
    """The C{bsf.data.PairedReads} class represents a pair of C{bsf.data.Reads} objects
    representing a read pair (i.e. I{R1} and I{R2}).

    Attributes:
    @ivar reads1: First C{bsf.data.Reads} object
    @type reads1: bsf.data.Reads
    @ivar reads2: Second C{bsf.data.Reads} object
    @type reads2: bsf.data.Reads
    @ivar index_1: Index sequence 1
    @type index_1: str
    @ivar index_2: Index sequence 2
    @type index_2: str
    @ivar annotation: Python C{dict} for annotation of Python C{str} key and
        Python C{list} of Python C{str} value data
    @type annotation: dict[str, list[str]]
    @ivar exclude: Exclude from processing
    @type exclude: bool
    @ivar read_group: SAM read group (@RG) information
    @type read_group: str
    @ivar weak_reference_sample: C{weakref.ReferenceType} to a C{bsf.data.Sample}
    @type weak_reference_sample: weakref.ReferenceType
    """

    def __init__(
            self,
            reads1=None,
            reads2=None,
            annotation=None,
            exclude=None,
            index_1=None,
            index_2=None,
            read_group=None,
            weak_reference_sample=None):
        """Initialise a C{bsf.data.PairedReads} object.

        For the C{bsf.data.Reads.file_type} I{CASAVA} the reads object will be
        automatically assigned on the basis of the C{bsf.data.Reads.read}
        attribute (i.e. I{R1} or I{R2}).
        @param reads1: First C{bsf.data.Reads} object
        @type reads1: bsf.data.Reads
        @param reads2: Second C{bsf.data.Reads} object
        @type reads2: bsf.data.Reads
        @param index_1: Index sequence 1
        @type index_1: str
        @param index_2: Index sequence 2
        @type index_2: str
        @param annotation: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation: dict[str, list[str]]
        @param exclude: Exclude from processing
        @type exclude: bool
        @param read_group: SAM read group (@RG) information
        @type read_group: str
        @param weak_reference_sample: C{weakref.ReferenceType} pointing at a C{bsf.data.Sample}
        @type weak_reference_sample: weakref.ReferenceType
        @return:
        @rtype:
        @raise Exception: For C{bsf.data.Reads.file_type} I{CASAVA}, I{R1} or I{R2} must be set in the
            C{bsf.dta.Reads} object.
        """

        if (reads1 and reads2) and (not reads1.match_paired(reads=reads2)):
            raise Exception('The Reads objects do not match.')

        super(PairedReads, self).__init__()

        self.reads1 = None
        self.reads2 = None

        if reads1:

            if reads1.file_type == 'CASAVA':
                if reads1.read == 'R1':
                    self.reads1 = reads1
                    reads1.weak_reference_paired_reads = weakref.ref(self)
                elif reads1.read == 'R2':
                    self.reads2 = reads1
                    reads1.weak_reference_paired_reads = weakref.ref(self)
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads1.read))
            else:
                # Other file types go here ...
                self.reads1 = reads1
                reads1.weak_reference_paired_reads = weakref.ref(self)

        if reads2:

            if reads2.file_type == 'CASAVA':
                if reads2.read == 'R1':
                    self.reads1 = reads2
                    reads2.weak_reference_paired_reads = weakref.ref(self)
                elif reads2.read == 'R2':
                    self.reads2 = reads2
                    reads2.weak_reference_paired_reads = weakref.ref(self)
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads2.read))
            else:
                # Other file types go here ...
                self.reads2 = reads2
                reads2.weak_reference_paired_reads = weakref.ref(self)

        if annotation is None:
            self.annotation = dict()
        else:
            self.annotation = annotation

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

    def trace(self, level):
        """Trace a C{bsf.data.PairedReads} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_sample: {!r}\n'.format(indent, self.weak_reference_sample)
        output += '{}  reads1: {!r}\n'.format(indent, self.reads1)
        output += '{}  reads2: {!r}\n'.format(indent, self.reads2)
        output += '{}  annotation:\n'.format(indent, self.annotation)
        for key in self.annotation.keys():
            output += '{}    {!r} {!r}\n'.format(indent, key, self.annotation[key])

        output += '{}  index_1: {!r}\n'.format(indent, self.index_1)
        output += '{}  index_2: {!r}\n'.format(indent, self.index_2)
        output += '{}  exclude: {!r}\n'.format(indent, self.exclude)
        output += '{}  read_group: {!r}\n'.format(indent, self.read_group)

        if isinstance(self.reads1, Reads):
            output += self.reads1.trace(level + 1)
        if isinstance(self.reads2, Reads):
            output += self.reads2.trace(level + 1)

        return output

    def add_reads(self, reads):
        """Add a C{bsf.data.Reads} object.

        For a C{bsf.data.Reads.file_type} I{CASAVA} the C{bsf.data.Reads} object can be automatically
        assigned on the basis of the C{bsf.data.Reads.read} attribute (i.e. I{R1} or I{R2}).
        @param reads: C{bsf.data.Reads}
        @type reads: bsf.data.Reads
        @return: C{True} upon success, C{False} otherwise
        @rtype: bool
        """

        # For CASAVA projects, reads are automatically added to either
        # reads1 or reads2 according to the file name.
        # Returns True upon success, False otherwise.

        assert isinstance(reads, Reads)

        if self.reads1:

            if not self.reads1.match_paired(reads=reads):
                return False

            if self.reads1.file_type == 'CASAVA':
                if reads.read == 'R1':
                    raise Exception(
                        'PairedReads reads1 has already been defined.\n'
                        '  reads1: {!r}\n'
                        '  reads:  {!r}'.format(self.reads1.file_path, reads.file_path))
                elif reads.read == 'R2':
                    self.reads2 = reads
                    reads.weak_reference_paired_reads = weakref.ref(self)
                    return True
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads.read))
            else:
                # Other file types go here ...
                warnings.warn(
                    'Method not implemented for file types other than CASAVA.',
                    UserWarning)

        if self.reads2:

            if not self.reads2.match_paired(reads=reads):
                return False

            if self.reads2.file_type == 'CASAVA':
                if reads.read == 'R1':
                    self.reads1 = reads
                    reads.weak_reference_paired_reads = weakref.ref(self)
                    return True
                elif reads.read == 'R2':
                    raise Exception(
                        'PairedReads reads2 has already been defined.\n'
                        '  reads2: {!r}\n'
                        '  reads:  {!r}'.format(self.reads2.file_path, reads.file_path))
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads.read))
            else:
                # Other file types go here ...
                warnings.warn(
                    'Method not implemented for file types other than CASAVA.',
                    UserWarning)

        return False

    def match(self, paired_reads):
        """Match C{bsf.data.PairedReads} objects.

        @param paired_reads: C{bsf.data.PairedReads}
        @type paired_reads: bsf.data.PairedReads
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

        if not self.reads1.match(reads=paired_reads.reads1):
            return False

        if (self.reads2 and paired_reads.reads2) and not self.reads2.match(reads=paired_reads.reads2):
            return False

        return True

    def get_name(self, full=False):
        """Get the name of a C{bsf.data.PairedReads} object.

        For the C{bsf.data.Reads.file_type} I{CASAVA} the name is a concatenation of the
        C{bsf.data.Reads.name}, C{bsf.data.Reads.barcode} and C{bsf.data.Reads.lane} attributes,
        preferentially derived from the first C{bsf.data.Reads} object in the
        C{bsf.data.PairedReads} object.
        If the I{full} parameter is set, C{bsf.data.Reads.read} and C{bsf.data.Reads.chunk} are also added.
        @param full: Return the full name including read and chunk information
        @type full: bool
        @return: Name
        @rtype: str
        """

        if self.reads1 is not None:
            reads = self.reads1
        elif self.reads2 is not None:
            reads = self.reads2
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


class Sample(object):
    """The C{bsf.data.Sample} class represents a Next-Generation Sequencing sample that
    consists of one or more C{bsf.data.PairedReads} objects as (biological or technical) replicates
    that result from the same flow cell.

    Attributes:
    @cvar default_key: Default key
    @type default_key: str
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @ivar name: Name
    @type name: str
    @ivar annotation: Python C{dict} for annotation of Python C{str} key and
        Python C{list} of Python C{str} value data
    @type annotation: dict[str, list[str]]
    @ivar paired_reads_list: Python C{list} of C{bsf.data.PairedReads} objects
    @type paired_reads_list: list[bsf.data.PairedReads]
    @ivar weak_reference_project: C{weakref.ReferenceType} pointing at a C{bsf.data.Project} object
    @type weak_reference_project: weakref.ReferenceType
    """

    default_key = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.data.Sample} object from a file path.

        For a I{file_type} I{CASAVA} the name is automatically populated,
        while C{bsf.data.PairedReads} objects are automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.data.Sample}
        @rtype: bsf.data.Sample
        """

        if file_type == 'CASAVA':
            file_path = os.path.normpath(file_path)
            file_name = os.path.basename(file_path)

            # CASAVA Samples obey a "Sample_name" schema.

            components = file_name.split('_')

            name = '_'.join(components[1:])

            sample = cls(file_path=file_path, file_type=file_type, name=name)

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
        """Create a merged C{bsf.data.Sample} from two C{bsf.data.Sample} objects.

        @param sample1: C{bsf.data.Sample}
        @type sample1: bsf.data.Sample
        @param sample2: C{bsf.data.Sample}
        @type sample2: bsf.data.Sample
        @return: C{bsf.data.Sample}
        @rtype: bsf.data.Sample
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
            file_path=None,
            file_type=None,
            name=None,
            annotation=None,
            paired_reads_list=None,
            weak_reference_project=None):
        """Initialise a C{bsf.data.Sample} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @param name: Name
        @type name: str
        @param annotation: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation: dict[str, list[str]]
        @param paired_reads_list: Python C{list} of C{bsf.data.PairedReads} objects
        @type paired_reads_list: list[bsf.data.PairedReads]
        @param weak_reference_project: C{weakref.ReferenceType} pointing at a C{bsf.data.Project} object
        @type weak_reference_project: weakref.ReferenceType
        @return:
        @rtype:
        """

        super(Sample, self).__init__()

        if file_path is None:
            self.file_path = str()
        else:
            self.file_path = file_path

        if file_type is None:
            self.file_type = str()
        else:
            self.file_type = file_type

        if name is None:
            self.name = str()
        else:
            self.name = name

        if annotation is None:
            self.annotation = dict()
        else:
            self.annotation = annotation

        if paired_reads_list is None:
            self.paired_reads_list = list()
        else:
            self.paired_reads_list = paired_reads_list
            # Setting the Sample object as a weak reference is problematic as it modifies the
            # PairedReads objects on the list.
            for paired_reads_object in paired_reads_list:
                paired_reads_object.weak_reference_sample = weakref.ref(self)

        self.weak_reference_project = weak_reference_project  # Can be None.

        return

    def trace(self, level):
        """Trace a C{bsf.data.Sample} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_project: {!r}\n'.format(indent, self.weak_reference_project)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  annotation:\n'.format(indent, self.annotation)
        for key in self.annotation.keys():
            output += '{}    {!r} {!r}\n'.format(indent, key, self.annotation[key])

        output += '{}  paired_reads:\n'.format(indent)

        for paired_reads in self.paired_reads_list:
            output += paired_reads.trace(level + 1)

        return output

    def match(self, sample):
        """Match C{bsf.data.Sample} objects.

        @param sample: C{bsf.data.Sample}
        @type sample: bsf.data.Sample
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

        for paired_reads1 in self.paired_reads_list:

            match = False

            for paired_reads2 in sample.paired_reads_list:

                if paired_reads1.match(paired_reads=paired_reads2):
                    match = True
                    break

            if not match:
                # No match for this particular PairedReads object.
                return False

        return True

    def add_annotation(self, key, value):
        """Add an annotation value under a key.

        @param key: Annotation key
        @type key: str
        @param value: Annotation value
        @type value: str
        @return:
        @rtype:
        """

        if key in self.annotation:
            value_list = self.annotation[key]
        else:
            value_list = list()
            self.annotation[key] = value_list

        if value not in value_list:
            value_list.append(value)

        return

    def add_paired_reads(self, paired_reads):
        """Add a C{bsf.data.PairedReads} object.

        This method checks, whether a matching C{bsf.data.PairedReads} object is already present in this
        C{bsf.data.Sample} object.
        @param paired_reads: C{bsf.data.PairedReads}
        @type paired_reads: bsf.data.PairedReads
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
            paired_reads.weak_reference_sample = weakref.ref(self)

        return

    def add_reads(self, reads):
        """Add a C{bsf.data.Reads} object.

        @param reads: C{bsf.data.Reads}
        @type reads: bsf.data.Reads
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
            self.add_paired_reads(paired_reads=PairedReads(reads1=reads))

        return

    def get_all_paired_reads(self, replicate_grouping, exclude=False, full=False):
        """Get all C{bsf.data.PairedReads} objects of a C{bsf.data.Sample} grouped or un-grouped.

        A C{bsf.data.Sample} object can hold several C{bsf.data.PairedReads} objects (i.e. biological or technical
        replicates) that have been sequenced on different lanes of the same flow cell.
        The C{bsf.data.PairedReads} objects will therefore differ in I{name}, I{barcode} or I{lane} information.
        Depending on the I{replicate_grouping} parameter they can be returned as a group or separately.
        C{bsf.data.PairedReads} objects that share the I{name}, I{barcode} and I{lane}, but differ in I{chunk} number
        are always grouped together.
        @param replicate_grouping: Group all C{bsf.data.PairedReads} objects of a C{bsf.data.Sample} or
            list them individually
        @type replicate_grouping: bool
        @param exclude: Exclude on the basis of C{bsf.data.PairedReads.exclude}
        @type exclude: bool
        @param full: Return the full name including read and chunk information
        @type full: bool
        @return: Python C{dict} of Python C{str} (sensible replicate name) key and
            Python C{list} object of C{bsf.data.PairedReads} objects value data
        @rtype: dict[str, list[bsf.data.PairedReads]]
        """

        groups = dict()

        for paired_reads in self.paired_reads_list:

            if exclude and paired_reads.exclude:
                # Skip excluded PairedReads objects.
                continue

            if replicate_grouping:
                # If grouped, use the Sample.name instance variable as key so that all PairedReds objects
                # of this Sample object end up on the same Python list object.
                key = self.name
            else:
                # If un-grouped use the PairedReads.get_name() method as key so that each
                # PairedReads object of this Sample object ends up on a separate Python list
                # object.
                key = paired_reads.get_name(full=full)

            if not key:
                continue

            if key not in groups:
                groups[key] = list()

            groups[key].append(paired_reads)

        # Return all values of the reads dictionary, which should be
        # a Python list of Python list objects of PairedReads objects.

        return groups


class Project(object):
    """The C{bsf.data.Project} class represents a Next-Generation Sequencing Project
    consisting of one or more C{bsf.data.Sample} objects.

    Attributes:

    @cvar default_key: Default key
    @type default_key: str
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar annotation: Python C{dict} for annotation of Python C{str} key and
        Python C{list} of Python C{str} value data
    @type annotation: dict[str, list[str]]
    @ivar samples: Python C{dict} of C{bsf.data.Sample.name} key objects and C{bsf.data.Sample} value objects
    @type samples: dict[bsf.data.Sample.name, bsf.data.Sample]
    @ivar weak_reference_prf: C{weakref.ReferenceType} pointing at a C{bsf.data.ProcessedRunFolder} object
    @type weak_reference_prf: weakref.ReferenceType
    """

    default_key = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.data.Project} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.data.Project}
        @rtype: bsf.data.Project
        """

        if file_type == 'CASAVA':
            file_path = os.path.normpath(file_path)
            file_name = os.path.basename(file_path)

            # CASAVA Projects obey a "Project_name" schema.

            components = file_name.split('_')

            name = '_'.join(components[1:])

            project = cls(file_path=file_path, file_type=file_type, name=name)

            # Automatically discover CASAVA Sample directories ...

            re_pattern = re.compile(pattern=r'^Sample_(.*)$')
            for file_name in os.listdir(project.file_path):
                file_path = os.path.join(project.file_path, file_name)
                file_mode = os.stat(file_path).st_mode
                re_match = re_pattern.search(string=file_name)
                if stat.S_ISDIR(file_mode) and re_match is not None:
                    if re_match.group(1) in project.samples:
                        raise Exception(
                            'Sample with name {!r} already exists.'.format(re_match.group(1)))
                    else:
                        project.add_sample(sample=Sample.from_file_path(file_path=file_path, file_type=file_type))
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return project

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            annotation=None,
            samples=None,
            weak_reference_prf=None):
        """Initialise a C{bsf.data.Project} object.

        For a I{file_type} I{CASAVA} the name is automatically populated,
        while C{bsf.data.Sample} objects are automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param annotation: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation: dict[str, list[str]]
        @param samples: Python C{dict} of C{bsf.data.Sample.name} key objects and C{bsf.data.Sample} value objects
        @type samples: dict[bsf.data.Sample.name, bsf.data.Sample]
        @param weak_reference_prf: C{weakref.ReferenceType} pointing at a C{bsf.data.ProcessedRunFolder} object
        @type weak_reference_prf: weakref.ReferenceType
        @raise Exception: If C{bsf.data.Sample.name} values are not unique for I{file_type} I{CASAVA}
        @return:
        @rtype:
        """

        super(Project, self).__init__()

        if file_path is None:
            self.file_path = str()
        else:
            self.file_path = file_path

        if file_type is None:
            self.file_type = str()
        else:
            self.file_type = file_type

        if name is None:
            self.name = str()
        else:
            self.name = name

        if annotation is None:
            self.annotation = dict()
        else:
            self.annotation = annotation

        if samples is None:
            self.samples = dict()
        else:
            self.samples = samples

        self.weak_reference_prf = weak_reference_prf  # Can be None.

        return

    def trace(self, level):
        """Trace a C{bsf.data.Project} object.

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
        output += '{}  annotation:\n'.format(indent, self.annotation)
        for key in self.annotation.keys():
            output += '{}    {!r} {!r}\n'.format(indent, key, self.annotation[key])

        output += '{}  samples:\n'.format(indent)

        for key in self.samples.keys():
            output += self.samples[key].trace(level + 1)

        return output

    def add_annotation(self, key, value):
        """Add an annotation value under a key.

        @param key: Annotation key
        @type key: str
        @param value: Annotation value
        @type value: str
        @return:
        @rtype:
        """

        if key in self.annotation:
            value_list = self.annotation[key]
        else:
            value_list = list()
            self.annotation[key] = value_list

        if value not in value_list:
            value_list.append(value)

        return

    def add_sample(self, sample):
        """Add a C{bsf.data.Sample} object to a C{bsf.data.Project} object and set
        a weak reference in to C{bsf.data.Sample} object back to the C{bsf.data.Project} object.

        @param sample: C{bsf.data.Sample}
        @type sample: bsf.data.Sample
        @return: C{bsf.data.Sample}
        @rtype: bsf.data.Sample
        """

        assert isinstance(sample, Sample)
        self.samples[sample.name] = sample
        sample.weak_reference_project = weakref.ref(self)

        return sample

    def del_sample(self, name):
        """Delete a C{bsf.data.Sample} object from am C{bsf.data.Project} object and
        clear the weak reference, if it points back at the C{bsf.data.Project} object.

        @param name: C{bsf.data.Sample.name}
        @type name: str
        @return: C{bsf.data.Sample}
        @rtype: bsf.data.Sample
        """

        if name in self.samples:
            sample = self.samples[name]
            del self.samples[name]
            if sample.weak_reference_project() == self:
                sample.weak_reference_project = None
            return sample
        else:
            return

    def get_all_samples(self):
        """Get an ordered Python C{list} of C{bsf.data.Sample} objects.

        @return: Python C{list} of C{bsf.data.Sample} objects
        @rtype: list[bsf.data.Sample]
        """

        samples = list()

        keys = self.samples.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            samples.append(self.samples[key])

        return samples


class ProcessedRunFolder(object):
    """The C{bsf.data.ProcessedRunFolder} class represents an (Illumina)
    run folder after processing with CASAVA.

    Attributes:
    @cvar default_key: Default key
    @type default_key: str
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA.
        I{External}: other data files.
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar prefix: Prefix
    @type prefix: str
    @ivar flow_cell: Flow cell identifier
    @type flow_cell: str
    @ivar version: Version number
    @type version: str
    @ivar annotation: Python C{dict} for annotation of Python C{str} key and
        Python C{list} of Python C{str} value data
    @type annotation: dict[str, list[str]]
    @ivar projects: Python C{dict} of C{bsf.data.Project.name} key objects and C{bsf.data.Project} value objects
    @type projects: dict[bsf.data.Project.name, bsf.data.Project]
    @ivar weak_reference_collection: C{weakref.ReferenceType} pointing at a C{bsf.data.Collection} object
    @type weak_reference_collection: weakref.ReferenceType
    """

    default_key = 'Default'

    @staticmethod
    def guess_file_type(file_path):
        """Guess the I{file_type} of a C{bsf.data.ProcessedRunFolder} on the basis of the I{file_path}.

        CASAVA C{bsf.data.ProcessedRunFolder} objects obey a I{Prefix_FCID_CASAVA182}
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

        components = file_name.split('_')

        if components[-1].startswith('CASAVA'):
            return 'CASAVA'
        else:
            return 'External'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{bsf.data.ProcessedRunFolder} object from a file path.

        For the I{file_type} I{CASAVA}, the C{bsf.data.ProcessedRunFolder.name}, C{bsf.data.ProcessedRunFolder.prefix},
        C{bsf.data.ProcessedRunFolder.flow_cell} and C{bsf.data.ProcessedRunFolder.version}
        attributes can be automatically parsed from the I{file_path}, while
        C{bsf.data.Project} objects can be automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.data.ProcessedRunFolder}
        @rtype: bsf.data.ProcessedRunFolder
        """

        # Try to determine the file_type if not explicitly specified.

        if not file_type or file_type == 'Automatic':
            file_type = ProcessedRunFolder.guess_file_type(file_path=file_path)

        if file_type == 'CASAVA':
            file_path = os.path.normpath(file_path)
            file_name = os.path.basename(file_path)

            # CASAVA Processed Run Folders obey a "Prefix_FCID_CASAVA182"
            # schema. The following prefixes are currently in use:
            # -- BSF_ Biomedical Sequencing Facility
            # -- NGS_ Kaan Boztug group
            # -- MUW_ Medical University Vienna
            # -- SET_ Robert Kralovics group

            components = file_name.split('_')

            name = '_'.join(components[:])
            prefix = components[0]
            flow_cell = components[1]
            version = components[2]

            prf = cls(file_path=file_path, file_type=file_type, name=name,
                      prefix=prefix, flow_cell=flow_cell, version=version)

            # Automatically discover CASAVA Project directories ...

            re_pattern = re.compile(pattern=r'^Project_(.*)$')
            for file_name in os.listdir(prf.file_path):
                file_path = os.path.join(prf.file_path, file_name)
                file_mode = os.stat(file_path).st_mode
                re_match = re_pattern.search(string=file_name)
                if stat.S_ISDIR(file_mode) and re_match is not None:
                    if re_match.group(1) in prf.projects:
                        raise Exception(
                            'Project with name {!r} already exists.'.format(re_match.group(1)))
                    else:
                        prf.add_project(project=Project.from_file_path(file_path=file_path, file_type=file_type))
        elif file_type == 'External':
            # Create a new, minimal ProcessedRunFolder.
            file_path = os.path.normpath(file_path)
            file_name = os.path.basename(file_path)

            prf = cls(file_path=file_path, file_type=file_type, name=file_name)
        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return prf

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            prefix=None,
            flow_cell=None,
            version=None,
            annotation=None,
            projects=None,
            weak_reference_collection=None):
        """Initialise a C{bsf.data.ProcessedRunFolder} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str
        @param name: Name
        @type name: str
        @param prefix: Prefix
        @type prefix: str
        @param flow_cell: Flow cell identifier
        @type flow_cell: str
        @param version: Version
        @type version: str
        @param annotation: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation: dict[str, list[str]]
        @param projects: Python C{dict} of C{bsf.data.Project.name} key objects and C{bsf.data.Project} value objects
        @type projects: dict[bsf.data.Project.name, bsf.data.Project]
        @param weak_reference_collection: C{weakref.ReferenceType} pointing at a C{bsf.data.Collection} object
        @type weak_reference_collection: weakref.ReferenceType
        @return:
        @rtype:
        @raise Exception: If C{bsf.data.Project.name} values are not unique for file_type I{CASAVA}
        """

        super(ProcessedRunFolder, self).__init__()

        if file_path is None:
            self.file_path = str()
        else:
            self.file_path = file_path

        if file_type is None:
            self.file_type = str()
        else:
            self.file_type = file_type

        if name is None:
            self.name = str()
        else:
            self.name = name

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

        if annotation is None:
            self.annotation = dict()
        else:
            self.annotation = annotation

        if projects is None:
            self.projects = dict()
        else:
            self.projects = projects

        self.weak_reference_collection = weak_reference_collection  # Can be None.

        return

    def trace(self, level):
        """Trace a C{bsf.data.ProcessedRunFolder} object.

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
        output += '{}  annotation:\n'.format(indent, self.annotation)
        for key in self.annotation.keys():
            output += '{}    {!r} {!r}\n'.format(indent, key, self.annotation[key])

        output += '{}  projects:\n'.format(indent)

        for key in self.projects.keys():
            output += self.projects[key].trace(level + 1)

        return output

    def add_annotation(self, key, value):
        """Add an annotation value under a key.

        @param key: Annotation key
        @type key: str
        @param value: Annotation value
        @type value: str
        @return:
        @rtype:
        """

        if key in self.annotation:
            value_list = self.annotation[key]
        else:
            value_list = list()
            self.annotation[key] = value_list

        if value not in value_list:
            value_list.append(value)

        return

    def add_project(self, project):
        """Add a C{bsf.data.Project} object to the C{bsf.data.ProcessedRunFolder} object and set
        a weak reference in the C{bsf.data.Project} object to the C{bsf.data.ProcessedRunFolder} object.

        @param project: C{bsf.data.Project}
        @type project: bsf.data.Project
        @return: C{bsf.data.Project}
        @rtype: bsf.data.Project
        """

        assert isinstance(project, Project)
        self.projects[project.name] = project
        project.weak_reference_prf = weakref.ref(self)

        return project

    def del_project(self, name):
        """Delete a C{bsf.data.Project} object from am C{bsf.data.ProcessedRunFolder} object and
        clear the weak reference, if it points back at the C{bsf.data.ProcessedRunFolder} object.

        @param name: C{bsf.data.Project.name}
        @type name: str
        @return: C{bsf.data.Project}
        @rtype: bsf.data.Project
        """

        if name in self.projects:
            project = self.projects[name]
            del self.projects[name]
            if project.weak_reference_prf() == self:
                project.weak_reference_prf = None
            return project
        else:
            return

    def get_all_projects(self):
        """Get an ordered Python C{list} of C{bsf.data.Project} objects.

        @return: A Python C{list} of C{bsf.data.Project} objects
        @rtype: list[bsf.data.Project]
        """

        projects = list()

        keys = self.projects.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            projects.append(self.projects[key])

        return projects


class Collection(object):
    """The C{bsf.data.Collection} class represents a collection of
    one or more C{bsf.data.ProcessedRunFolder} objects.

    Attributes:
    @cvar default_key: Default key
    @type default_key: str
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar annotation: Python C{dict} for annotation of Python C{str} key and
        Python C{list} of Python C{str} value data
    @type annotation: dict[str, list[str]]
    @ivar processed_run_folders: Python C{dict} of C{bsf.data.ProcessedRunFolder.name} key objects and
        C{bsf.data.ProcessedRunFolder} value objects
    @type processed_run_folders: dict[bsf.data.ProcessedRunFolder.name, bsf.data.ProcessedRunFolder]
    @ivar sample_groups: Python C{dict} of Python C{str} (group) and C{bsf.data.Sample} objects
    @type sample_groups: dict[str, list[bsf.data.Sample]]
    """

    default_key = 'Default'

    @classmethod
    def from_sas_path(cls, file_path, file_type, name, sas_path, sas_prefix=None):
        """Construct a C{bsf.data.Collection} from a C{bsf.data.SampleAnnotationSheet} file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param sas_path: C{bsf.data.SampleAnnotationSheet} file path
        @type sas_path: str | unicode
        @param sas_prefix: Optional column header prefix
            (e.g. '[Control ]Sample', '[Treatment ]Sample', ...)
        @type sas_prefix: str
        @return: C{bsf.data.Collection}
        @rtype: bsf.data.Collection
        """

        return cls.from_sas(
            file_path=file_path,
            file_type=file_type,
            name=name,
            sas=SampleAnnotationSheet.from_file_path(file_path=sas_path),
            sas_prefix=sas_prefix)

    @classmethod
    def from_sas(cls, file_path, file_type, name, sas, sas_prefix=None):
        """Construct a C{bsf.data.Collection} from a C{bsf.data.SampleAnnotationSheet}.

        This method creates C{bsf.data.Reads}, C{bsf.data.PairedReads}, C{bsf.data.Sample}, C{bsf.data.Project} and
        C{bsf.data.ProcessedRunFolder} objects from a C{bsf.data.SampleAnnotationSheet} row Python C{dict} object.
        If object-specific keys or their corresponding values are not
        available from the Python row dict, new objects will be created
        with a default key.
        Sample Annotation Sheet format:
            - FileType: Data object I{file_type} (i.e. I{CASAVA} or I{External}), defaults to I{External}.
                The FileType I{CASAVA} allows for auto-discovery of C{bsf.data.ProcessedRunFolder} objects.
            - ProcessedRunFolder Name: C{bsf.data.ProcessedRunFolder} file_path, can be automatically registered.
            - Project Name: C{bsf.data.Project} name
            - Sample Name: C{bsf.data.Sample} name
            - Reads1 File: C{bsf.data.Reads.file_name} instance variable.
                Subjected to C{os.path.expanduser} and C{os.path.expandvars}.
                If still relative, the C{bsf.data.Collection.file_path} is prepended.
            - Reads1 Name: C{bsf.data.Reads.name} instance variable
            - Reads2 File: Same as Reads1 File
            - Reads2 Name: Same as Reads1 Name
            - Group: C{bsf.data.Sample} objects can be grouped for further analysis in
                e.g. RNA-Seq or ChIP-Seq experiments.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param sas: C{bsf.data.SampleAnnotationSheet}
        @type sas: bsf.data.SampleAnnotationSheet
        @param sas_prefix: Optional column header prefix
            (e.g. '[Control ]Sample', '[Treatment ]Sample', '[Point N ]Sample', ...)
        @type sas_prefix: str
        @return: C{bsf.data.Collection}
        @rtype: bsf.data.Collection
        """
        assert isinstance(sas, SampleAnnotationSheet)

        collection = cls(file_path=file_path, file_type=file_type, name=name)

        if sas_prefix is None:
            sas_prefix = str()

        prf = None
        project = None
        sample = None

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

            reads1 = collection._process_reads(
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix,
                file_type=file_type,
                suffix='1')

            reads2 = collection._process_reads(
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix,
                file_type=file_type,
                suffix='2')

            pr_exclude, pr_index_1, pr_index_2, pr_read_group = collection._process_paired_reads(
                row_dict=row_dict,
                key_list=key_list,
                prefix=sas_prefix)

            # If none of the Reads objects has been defined, the Sample
            # may have been automatically loaded from a CASAVA ProcessedRunFolder.

            if reads1 or reads2:
                sample.add_paired_reads(paired_reads=PairedReads(
                    reads1=reads1,
                    reads2=reads2,
                    exclude=pr_exclude,
                    index_1=pr_index_1,
                    index_2=pr_index_2,
                    read_group=pr_read_group))

            # Optionally group the Sample objects.

            key = '{} Group'.format(sas_prefix).lstrip()

            if key in key_list:
                key_list.remove(key)

            if key in row_dict and row_dict[key]:

                value = row_dict[key]

                if value in collection.sample_groups:
                    sample_group = collection.sample_groups[value]
                else:
                    sample_group = list()
                    collection.sample_groups[value] = sample_group

                if sample not in sample_group:
                    sample_group.append(sample)

            if len(key_list):
                warnings.warn("Unexpected keys in sample annotation sheet: {!r}".format(key_list), UserWarning)

        # Quench empty default objects that are a consequence of empty lines in the sample annotation sheet.

        for prf_key in collection.processed_run_folders.keys():
            prf = collection.processed_run_folders[prf_key]
            assert isinstance(prf, ProcessedRunFolder)
            for project_key in prf.projects.keys():
                project = prf.projects[project_key]
                assert isinstance(project, Project)
                for sample_key in project.samples.keys():
                    sample = project.samples[sample_key]
                    assert isinstance(sample, Sample)
                    if sample.name == Sample.default_key and not len(sample.paired_reads_list):
                        project.del_sample(name=sample.name)
                if project.name == Project.default_key and not len(project.samples):
                    prf.del_project(name=project.name)
            if prf.name == ProcessedRunFolder.default_key and not prf.projects:
                collection.del_processed_run_folder(name=prf.name)

        return collection

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            annotation=None,
            processed_run_folders=None,
            sample_groups=None):
        """Initialise a C{bsf.data.Collection} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param annotation: Python C{dict} for annotation of Python C{str} key and
            Python C{list} of Python C{str} value data
        @type annotation: dict[str, list[str]]
        @param processed_run_folders: Python C{dict} of C{bsf.data.ProcessedRunFolder.name} key objects and
            C{bsf.data.ProcessedRunFolder} value objects
        @type processed_run_folders: dict[bsf.data.ProcessedRunFolder.name, bsf.data.ProcessedRunFolder]
        @param sample_groups: Python C{dict} of Python C{str} (group name) key objects and
            second-level Python C{dict} value objects of C{bsf.data.Sample.name} key objects and
            C{bsf.data.Sample} value objects
        @type sample_groups: dict[str, dict[bsf.data.Sample.name, bsf.data.Sample]]
        @return:
        @rtype:
        """

        super(Collection, self).__init__()

        if file_path is None:
            self.file_path = str()
        else:
            self.file_path = file_path

        if file_type is None:
            self.file_type = str()
        else:
            self.file_type = file_type

        if name is None:
            self.name = str()
        else:
            self.name = name

        if annotation is None:
            self.annotation = dict()
        else:
            self.annotation = annotation

        if processed_run_folders is None:
            self.processed_run_folders = dict()
        else:
            self.processed_run_folders = processed_run_folders

        if sample_groups is None:
            self.sample_groups = dict()
        else:
            self.sample_groups = sample_groups

        return

    def trace(self, level):
        """Trace a C{bsf.data.Collection} object.

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

        output += '{}  annotation:\n'.format(indent, self.annotation)
        for key in self.annotation.keys():
            output += '{}    {!r} {!r}\n'.format(indent, key, self.annotation[key])

        output += '{}  processed_run_folders:\n'.format(indent)
        for key in self.processed_run_folders.keys():
            output += self.processed_run_folders[key].trace(level + 1)

        output += '{}  sample_groups:\n'.format(indent)

        keys = self.sample_groups.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            output += '{}    group: {!r}\n'.format(indent, key)
            # List all Sample objects of this Python list object.
            for sample in self.sample_groups[key]:
                output += '{}      Sample name: {!r} file_path: {!r}\n'.format(indent, sample.name, sample.file_path)

        return output

    def add_processed_run_folder(self, prf):
        """Add a C{bsf.data.ProcessedRunFolder} object to the C{bsf.data.Collection} object and
        set a weak reference in the C{bsf.data.ProcessedRunFolder} back to the C{bsf.data.Collection}.

        @param prf: C{bsf.data.ProcessedRunFolder}
        @type prf: bsf.data.ProcessedRunFolder
        @return: C{bsf.data.ProcessedRunFolder}
        @rtype: bsf.data.ProcessedRunFolder
        """

        assert isinstance(prf, ProcessedRunFolder)
        self.processed_run_folders[prf.name] = prf
        prf.weak_reference_collection = weakref.ref(self)

        return prf

    def del_processed_run_folder(self, name):
        """Delete a C{bsf.data.ProcessedRunFolder} object from am C{bsf.data.Collection} object and
        clear the weak reference, if it points back at the C{bsf.data.Collection} object.

        @param name: C{bsf.data.ProcessedRunFolder.name}
        @type name: str
        @return: C{bsf.data.ProcessedRunFolder}
        @rtype: bsf.data.ProcessedRunFolder
        """

        if name in self.processed_run_folders:
            prf = self.processed_run_folders[name]
            del self.processed_run_folders[name]
            if prf.weak_reference_collection() == self:
                prf.weak_reference_collection = None
            return prf
        else:
            return

    def get_processed_run_folder(self, file_path, file_type=None):
        """Get a C{bsf.data.ProcessedRunFolder} by file path.

        If the C{bsf.data.ProcessedRunFolder} does not exist in the C{bsf.data.Collection} object,
        it will be automatically discovered and added.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{bsf.data.ProcessedRunFolder}
        @rtype: bsf.data.ProcessedRunFolder
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        if file_name in self.processed_run_folders:
            return self.processed_run_folders[file_name]
        else:
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if not os.path.isabs(file_path):
                file_path = os.path.join(self.file_path, file_path)

            return self.add_processed_run_folder(
                prf=ProcessedRunFolder.from_file_path(file_path=file_path, file_type=file_type))

    def get_all_processed_run_folders(self):
        """Get an ordered Python C{list} of C{bsf.data.ProcessedRunFolder} objects.

        @return: Python C{list} of C{bsf.data.ProcessedRunFolder} objects
        @rtype: list[bsf.data.ProcessedRunFolder]
        """

        processed_run_folders = list()

        keys = self.processed_run_folders.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            processed_run_folders.append(self.processed_run_folders[key])

        return processed_run_folders

    def get_all_projects(self):
        """Get an ordered Python C{list} of C{bsf.data.Project} objects.

        @return: Python C{list} of C{bsf.data.Project} objects
        @rtype: list[bsf.data.Project]
        """

        projects = list()

        for processed_run_folder in self.get_all_processed_run_folders():
            projects.extend(processed_run_folder.get_all_projects())

        return projects

    def get_all_samples(self):
        """Get an ordered Python C{list} of C{bsf.data.Sample} objects.

        @return: Python C{list} of C{bsf.data.Sample} objects
        @rtype: list[bsf.data.Sample]
        """

        samples = list()

        for project in self.get_all_projects():
            samples.extend(project.get_all_samples())

        return samples

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
        """Get or create a C{bsf.data.ProcessedRunFolder}.

        A 'I{[Prefix] ProcessedRunFolder Name}' key is optional, its value defaults to I{Default}.
        @param prf: Current C{bsf.data.ProcessedRunFolder} that may get replaced upon encountering a new
            'I{[Prefix] ProcessedRunFolder Name}' key
        @type prf: bsf.data.ProcessedRunFolder
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
        @return: C{bsf.data.ProcessedRunFolder}
        @rtype: bsf.data.ProcessedRunFolder
        """

        new_prf = None

        key = '{} ProcessedRunFolder Name'.format(prefix).lstrip()

        if key in row_dict:
            # The key exists ...
            if row_dict[key]:
                # ... and has a meaningful value ...
                value = row_dict[key]
                if value in self.processed_run_folders:
                    # ..., which exists in the dict of ProcessedRunFolder objects.
                    new_prf = self.processed_run_folders[value]
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
            if ProcessedRunFolder.default_key in self.processed_run_folders:
                new_prf = self.processed_run_folders[ProcessedRunFolder.default_key]
            else:
                new_prf = self.add_processed_run_folder(
                    prf=ProcessedRunFolder(name=ProcessedRunFolder.default_key, file_type=file_type))

        # Look for additional ProcessedRunFolder keys that provide further annotation.

        re_pattern = re.compile(pattern='{} ProcessedRunFolder'.format(prefix).lstrip())
        for key1 in row_dict.keys():
            # Only search for the pattern at the start of the string.
            re_match = re.match(pattern=re_pattern, string=key1)
            if re_match is not None:
                # If present, delete the key from the key list, which includes the 'Name' key.
                if re_match.string in key_list:
                    key_list.remove(re_match.string)
                # Capture the string from the end of the match to the end of the string and strip white space.
                key2 = re_match.string[re_match.end(0):].strip()
                if key2 and key2 != 'Name':
                    # Exclude empty key strings and the 'Name' key that is not an annotation as such.
                    if row_dict[re_match.string]:
                        # Exclude empty fields.
                        new_prf.add_annotation(key=key2, value=row_dict[re_match.string])

        return new_prf

    @staticmethod
    def _process_project(project, row_dict, key_list, prefix, file_type, prf):
        """Get or create a C{bsf.data.Project}.

        A 'I{[Prefix] Project Name}' key is optional, its value defaults to I{Default}.
        @param project: Current C{bsf.data.Project} that may get replaced upon encountering a new
            'I{[Prefix] Project Name}' key
        @type project: bsf.data.Project
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Project', '[Treatment] Project', '[Point N] Project', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param prf: C{bsf.data.ProcessedRunFolder}
        @type prf: bsf.data.ProcessedRunFolder
        @return: C{bsf.data.Project}
        @rtype: bsf.data.Project
        """

        new_project = None

        key = '{} Project Name'.format(prefix).lstrip()

        if key in row_dict:
            # The key exists ...
            if row_dict[key]:
                # ... and has a meaningful value ...
                value = row_dict[key]
                if value in prf.projects:
                    # ..., which exists in the dict of Project objects.
                    new_project = prf.projects[value]
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
            if Project.default_key in prf.projects:
                new_project = prf.projects[Project.default_key]
            else:
                new_project = prf.add_project(project=Project(name=Project.default_key, file_type=file_type))

        # Look for additional Project keys that provide additional annotation.

        re_pattern = re.compile(pattern='{} Project'.format(prefix).lstrip())
        for key1 in row_dict.keys():
            # Only search for the pattern at the start of the string.
            re_match = re.match(pattern=re_pattern, string=key1)
            if re_match is not None:
                # If present, delete the key from the key list, which includes the 'Name' key.
                if re_match.string in key_list:
                    key_list.remove(re_match.string)
                # Capture the string from the end of the match to the end of the string and strip white space.
                key2 = re_match.string[re_match.end(0):].strip()
                if key2 and key2 != 'Name':
                    # Exclude empty key strings and the 'Name' key that is not an annotation as such.
                    if row_dict[re_match.string]:
                        # Exclude empty fields.
                        new_project.add_annotation(key=key2, value=row_dict[re_match.string])

        return new_project

    @staticmethod
    def _process_sample(sample, row_dict, key_list, prefix, file_type, project):
        """Get or create a C{bsf.data.Sample}.

        A 'I{[Prefix] Sample Name}' key is optional, its value defaults to I{Default}.
        @param sample: Current C{bsf.data.Sample} that may get replaced upon encountering a new
            'I{[Prefix] Sample Name}' key
        @type sample: bsf.data.Sample
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param project: C{bsf.data.Project}
        @type project: bsf.data.Project
        @return: C{bsf.data.Sample}
        @rtype: bsf.data.Sample
        """

        new_sample = None

        key = '{} Sample Name'.format(prefix).lstrip()

        if key in row_dict:
            # The key exists ...
            if row_dict[key]:
                # ... and has a meaningful value ...
                value = row_dict[key]
                if value in project.samples:
                    # ..., which exists in the dict of Sample objects.
                    new_sample = project.samples[value]
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
            if Sample.default_key in project.samples:
                new_sample = project.samples[Sample.default_key]
            else:
                new_sample = project.add_sample(sample=Sample(name=Sample.default_key, file_type=file_type))

        # Look for additional Sample keys that provide additional annotation.

        re_pattern = re.compile(pattern='{} Sample'.format(prefix).lstrip())
        for key1 in row_dict.keys():
            # Only search for the pattern at the start of the string.
            re_match = re.match(pattern=re_pattern, string=key1)
            if re_match is not None:
                # If present, delete the key from the key list, which includes the 'Name' key.
                if re_match.string in key_list:
                    key_list.remove(re_match.string)
                # Capture the string from the end of the match to the end of the string and strip white space.
                key2 = re_match.string[re_match.end(0):].strip()
                if key2 and key2 != 'Name':
                    # Exclude empty key strings and the 'Name' key that is not an annotation as such.
                    if row_dict[re_match.string]:
                        # Exclude empty fields.
                        new_sample.add_annotation(key=key2, value=row_dict[re_match.string])

        return new_sample

    def _process_reads(self, row_dict, key_list, prefix, file_type, suffix):
        """Get or create a C{bsf.data.Reads} object.

        A 'I{[Prefix] Reads{suffix} Name}' key and 'I{[Prefix] Reads{suffix} File}' key are optional,
        in which case the default is a C{None} object.
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
        @return: C{bsf.data.Reads}
        @rtype: bsf.data.Reads
        """

        key_file = '{} Reads{} File'.format(prefix, suffix).lstrip()
        key_name = '{} Reads{} Name'.format(prefix, suffix).lstrip()

        if key_file in key_list:
            key_list.remove(key_file)
        if key_name in key_list:
            key_list.remove(key_name)

        if key_file in row_dict and row_dict[key_file]:
            file_path = row_dict[key_file]
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if not os.path.isabs(file_path):
                file_path = os.path.join(self.file_path, file_path)

            file_path = os.path.normpath(file_path)

            if key_name in row_dict and row_dict[key_name]:
                reads = Reads(name=row_dict[key_name], file_path=file_path, file_type=file_type)
            else:
                reads = None
        else:
            reads = None

        return reads

    # Taken from ConfigParser.RawConfigParser.getboolean()

    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False}

    @staticmethod
    def _process_paired_reads(row_dict, key_list, prefix):
        """Get or create a C{bsf.data.PairedReads} object.

        The 'I{[Prefix] PairedReads Exclude}' key is optional, in which case the default is Python C{bool} I{False}.
        The 'I{[Prefix] PairedReads Index 1}', 'I{[Prefix] PairedReads Index 2}' and 'I{[Prefix] PairedReads ReadGroup}'
        keys are optional, in which case the default is an empty Python C{str} object.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict[str, str | unicode]
        @param key_list: A Python C{list} of Python C{str} (key) objects in the row
        @type key_list: list[str]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] PairedReads ReadGroup',
            '[Treatment] PairedReads ReadGroup',
            '[Point N] PairedReads ReadGroup', ...)
        @type prefix: str
        @return: Python C{tuple} of Python C{bool} (exclude flag), Python C{str} (read group), Python C{str} (index 1)
            and Python C{str} (index 2)
        @rtype: tuple[bool, str, str, str]
        """

        key = '{} PairedReads Exclude'.format(prefix).lstrip()

        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            if row_dict[key].lower() not in Collection._boolean_states:
                raise ValueError('Value in field {!r} is not a boolean: {!r}'.format(key, row_dict[key]))
            exclude = Collection._boolean_states[row_dict[key].lower()]
        else:
            exclude = False

        key = '{} PairedReads Index 1'.format(prefix).lstrip()

        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            index_1 = row_dict[key]
        else:
            index_1 = str()

        key = '{} PairedReads Index 2'.format(prefix).lstrip()

        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            index_2 = row_dict[key]
        else:
            index_2 = str()

        key = '{} PairedReads ReadGroup'.format(prefix).lstrip()

        if key in key_list:
            key_list.remove(key)

        if key in row_dict and row_dict[key]:
            read_group = row_dict[key]
        else:
            read_group = str()

        return exclude, index_1, index_2, read_group

    def get_sample_from_row_dict(self, row_dict, prefix=None):
        """Get a Sample from a C{bsf.data.SampleAnnotationSheet} row Python C{dict}.

        Look-up a hierarchy of C{bsf.data.ProcessedRunFolder}, C{bsf.data.Project} and C{bsf.data.Sample} objects
        based on a C{bsf.data.SampleAnnotationSheet} row dictionary.
        C{bsf.data.ProcessedRunFolder} objects of file type I{CASAVA} can be
        automatically discovered and registered.
        Return the corresponding C{bsf.data.Sample}.
        @param row_dict: C{bsf.data.SampleAnnotationSheet} row Python C{dict}
        @type row_dict: dict[str, str | unicode]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample')
        @type prefix: str
        @return: C{bsf.data.Sample}
        @rtype: bsf.data.Sample
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
            value = ProcessedRunFolder.default_key

        # The Collection.get_processed_run_folder method can automatically register
        # ProcessedRunFolder objects of file type 'CASAVA'.

        prf = self.get_processed_run_folder(file_path=value)

        key = '{} Project Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Project.default_key

        project = prf.projects[value]

        key = '{} Sample Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Sample.default_key

        sample = project.samples[value]

        return sample

    def get_samples_from_row_dict(self, row_dict, prefix=None):
        """Get a Python C{list} of C{bsf.data.Sample} objects and a Python C{str} of the Group column value
        as a Python C{tuple} from a C{bsf.data.SampleAnnotationSheet} row Python C{dict}.

        @param row_dict: Comparison CSV file row Python C{dict}
        @type row_dict: dict[str, str | unicode]
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        @return: Python C{tuple} of Python C{str} of '[Prefix] Group' column value and
            Python C{list} of C{bsf.data.Sample} objects
        @rtype: (str, list[bsf.data.Sample])
        """

        samples = list()
        value = str()

        if not prefix:
            prefix = str()

        key = '{} Group'.format(prefix).lstrip()

        # Does the key exist and is its value defined?
        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Extend the Python list with all Sample objects of this group.
            if value in self.sample_groups:
                samples.extend(self.sample_groups[value])

        key = '{} Sample Name'.format(prefix).lstrip()

        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Append the Sample object, if defined.

            sample = self.get_sample_from_row_dict(row_dict=row_dict, prefix=prefix)

            if sample:
                samples.append(sample)

        return value, samples

    def to_sas(self, file_path=None, name=None):
        """Convert a C{bsf.data.Collection} into a SampleAnnotationSheet object.

        @param file_path: File path
        @type file_path: str | unicode
        @param name: Name
        @type name: str
        @return: SampleAnnotationSheet object
        @rtype: SampleAnnotationSheet
        """

        sas = SampleAnnotationSheet(file_path=file_path, name=name)

        # Scan the Collection and its contained objects for additional (annotation) field names.

        for prf in self.processed_run_folders.itervalues():
            assert isinstance(prf, ProcessedRunFolder)
            for prf_annotation_key in prf.annotation.iterkeys():
                prf_annotation_field = ' '.join(('ProcessedRunFolder', prf_annotation_key))
                if prf_annotation_field not in sas.field_names:
                    sas.field_names.append(prf_annotation_field)
            for project in prf.projects.itervalues():
                assert isinstance(project, Project)
                for project_annotation_key in project.annotation.iterkeys():
                    project_annotation_field = ' '.join(('Project', project_annotation_key))
                    if project_annotation_field not in sas.field_names:
                        sas.field_names.append(project_annotation_field)
                for sample in project.samples.itervalues():
                    assert isinstance(sample, Sample)
                    for sample_annotation_key in sample.annotation.iterkeys():
                        sample_annotation_field = ' '.join(('Sample', sample_annotation_key))
                        if sample_annotation_field not in sas.field_names:
                            sas.field_names.append(sample_annotation_field)
                    for paired_reads in sample.paired_reads_list:
                        assert isinstance(paired_reads, PairedReads)
                        for paired_reads_annotation_key in paired_reads.annotation.iterkeys():
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

        prf_name_list = self.processed_run_folders.keys()
        prf_name_list.sort(cmp=lambda x, y: cmp(x, y))
        for prf_name in prf_name_list:
            prf = self.processed_run_folders[prf_name]
            assert isinstance(prf, ProcessedRunFolder)
            sas.row_dicts.append({
                'ProcessedRunFolder Name': prf.name,
            })
            prf_annotation_key_list = prf.annotation.keys()
            prf_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
            for prf_annotation_key in prf_annotation_key_list:
                prf_annotation_list = prf.annotation[prf_annotation_key]
                prf_annotation_field = ' '.join(('ProcessedRunFolder', prf_annotation_key))
                for annotation in prf_annotation_list:
                    sas.row_dicts.append({prf_annotation_field: annotation})
            project_name_list = prf.projects.keys()
            project_name_list.sort(cmp=lambda x, y: cmp(x, y))
            for project_name in project_name_list:
                project = prf.projects[project_name]
                assert isinstance(project, Project)
                sas.row_dicts.append({
                    'Project Name': project.name,
                })
                project_annotation_key_list = project.annotation.keys()
                project_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                for project_annotation_key in project_annotation_key_list:
                    project_annotation_list = project.annotation[project_annotation_key]
                    project_annotation_field = ' '.join(('Project', project_annotation_key))
                    for annotation in project_annotation_list:
                        sas.row_dicts.append({project_annotation_field: annotation})
                sample_name_list = project.samples.keys()
                sample_name_list.sort(cmp=lambda x, y: cmp(x, y))
                for sample_name in sample_name_list:
                    sample = project.samples[sample_name]
                    assert isinstance(sample, Sample)
                    sas.row_dicts.append({
                        'Sample Name': sample.name
                    })
                    sample_annotation_key_list = sample.annotation.keys()
                    sample_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                    for sample_annotation_key in sample_annotation_key_list:
                        sample_annotation_list = sample.annotation[sample_annotation_key]
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
                        if paired_reads.reads1 is not None:
                            row_dict['Reads1 Name'] = paired_reads.reads1.name
                            row_dict['Reads1 File'] = paired_reads.reads1.file_path
                        if paired_reads.reads2 is not None:
                            row_dict['Reads2 Name'] = paired_reads.reads2.name
                            row_dict['Reads2 File'] = paired_reads.reads2.file_path
                        sas.row_dicts.append(row_dict)
                        paired_reads_annotation_key_list = paired_reads.annotation.keys()
                        paired_reads_annotation_key_list.sort(cmp=lambda x, y: cmp(x, y))
                        for paired_reads_annotation_key in paired_reads_annotation_key_list:
                            paired_reads_annotation_list = paired_reads.annotation[paired_reads_annotation_key]
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
    """The C{bsf.data.SampleGroup} class represents a group of C{bsf.data.Sample} objects.

    The grouping is usually defined in a sample annotation sheet.
    Attributes:
    @ivar name: Name
    @type name: str
    @ivar samples: Python C{list} of C{bsf.data.Sample} objects
    @type samples: list[bsf.data.Sample]
    """

    # TODO: The SampleGroup class is currently not in use.
    # Sample and PairedReads objects from different ProcessRunFolder objects
    # (i.e. flow cells) could bear the same name, leading to problems with SGE job names.
    # This would need further re-thinking.

    def __init__(
            self,
            name=None,
            samples=None):
        """Initialise a C{bsf.data.SampleGroup} object.

        @param name: Name
        @type name: str
        @param samples: Python C{list} of C{bsf.data.Sample} objects
        @type samples: list[bsf.data.Sample]
        @return:
        @rtype:
        """

        super(SampleGroup, self).__init__()

        if name is None:
            self.name = str()
        else:
            self.name = name

        if samples is None:
            self.samples = list()
        else:
            self.samples = samples

        return

    def add_sample(self, sample):
        """Add a C{bsf.data.Sample} object.

        @param sample: C{bsf.data.Sample}
        @type sample: bsf.data.Sample
        @return:
        @rtype:
        """

        assert isinstance(sample, Sample)
        if sample not in self.samples:
            self.samples.append(sample)

        return

    def get_all_paired_reads(self, replicate_grouping):
        """Get all C{bsf.data.PairedReads} objects of a C{bsf.data.SampleGroup}.

        For the moment, replicates are C{bsf.data.PairedReads} objects that do not share
        anything but the chunk number.
        @param replicate_grouping: Group all C{bsf.data.PairedReads} objects of a C{bsf.data.Sample} or
            list them individually
        @type replicate_grouping: bool
        @return: Python C{dict} of Python C{str} (replicate name) key objects and
            Python C{list} value objects of C{bsf.data.PairedReads} objects value data
        @rtype: dict[str, list[bsf.data.PairedReads]]
        """

        group_dict = dict()

        for sample in self.samples:
            replicate_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping)
            for replicate_key in replicate_dict.keys():
                if replicate_key not in group_dict:
                    group_dict[replicate_key] = list()

                # group_dict[replicate_key].extend(replicate_dict[replicate_key])

                # Add PairedReads objects one-by-one and check if they are not already there.

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads not in group_dict[replicate_key]:
                        group_dict[replicate_key].append(paired_reads)

        return group_dict


class SampleAnnotationSheet(AnnotationSheet):
    """The C{bsf.data.SampleAnnotationSheet} class represents a Comma-Separated Value (CSV) table of sample information
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
