"""bsf.data

A package of classes and methods modelling data directories and files.
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
import os
import re
from stat import *
import string
import warnings
import weakref

from bsf.annotation import SampleAnnotationSheet


class Reads(object):
    """The C{Reads} class represents a file of Next-Generation Sequencing (NGS) reads,
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
    @ivar weak_reference_paired_reads: Weak Reference to a C{PairedReads} object
    @type weak_reference_paired_reads: PairedReads
    """

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{Reads} object from a file path.

        For a I{file_type} I{CASAVA}, C{Reads.file_path} obeys a I{Sample_name_Index_Lane_Read_Chunk} schema,
        so that C{Reads.name}, C{Reads.barcode}, C{Reads.lane}, C{Reads.read} and C{Reads.chunk}
        can be populated automatically.
        For I{file_type} I{External}, the attributes need to be populated manually.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{Reads} object
        @rtype: Reads
        """

        if file_type == 'CASAVA':

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Reads obey a Sample_name_Index_Lane_Read_Chunk schema.

            components = file_name.split('.')
            components[0] = components[0].split('_')

            name = string.join(words=components[0][:-4], sep='_')
            barcode = components[0][-4]
            lane = components[0][-3]
            read = components[0][-2]
            chunk = components[0][-1]

            reads = cls(file_path=file_path, file_type=file_type, name=name,
                        barcode=barcode, lane=lane, read=read, chunk=chunk)

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return reads

    def __init__(self, file_path=None, file_type=None, name=None,
                 barcode=None, lane=None, read=None, chunk=None,
                 weak_reference_paired_reads=None):
        """Initialise a C{Reads} object.

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
        @param weak_reference_paired_reads: Weak Reference to a C{PairedReads} object
        @type weak_reference_paired_reads: PairedReads
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if barcode:
            self.barcode = barcode
        else:
            self.barcode = str()

        if lane:
            self.lane = lane
        else:
            self.lane = str()

        if read:
            self.read = read
        else:
            self.read = str()

        if chunk:
            self.chunk = chunk
        else:
            self.chunk = str()

        if weak_reference_paired_reads:
            self.weak_reference_paired_reads = weak_reference_paired_reads
        else:
            self.weak_reference_paired_reads = None

    def trace(self, level):
        """Trace a C{Reads} object.

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
        """Match C{Reads} objects.

        Two C{Reads} objects are identical, if all their instance variables match.
        @param reads: Second C{Reads} object
        @type reads: Reads
        @return: True if both objects match, False otherwise
        @rtype: bool
        """

        assert isinstance(reads, Reads)

        # Quick test first - if the objects are identical, the rest has to match.

        if self == reads:
            return True

        if not self.file_path == reads.file_path:
            return False

        if not self.file_type == reads.file_type:
            return False

        if not self.name == reads.name:
            return False

        if not self.barcode == reads.barcode:
            return False

        if not self.lane == reads.lane:
            return False

        if not self.read == reads.read:
            return False

        if not self.chunk == reads.chunk:
            return False

        return True

    def match_paired(self, reads):
        """Match paired C{Reads} objects, by relaxing matching criteria.

        @param reads: Second C{Reads} object
        @type reads: Reads
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
        """

        assert isinstance(reads, Reads)

        if self.file_type == 'CASAVA':

            # All CASAVA Reads attributes need to be equal,
            # with the exception of the read (R1 or R2) and
            # the file_path, which also contains the reads
            # information.

            if not self.file_type == reads.file_type:
                return False

            if not self.name == reads.name:
                return False

            if not self.barcode == reads.barcode:
                return False

            if not self.lane == reads.lane:
                return False

            if not self.chunk == reads.chunk:
                return False

            return True
        else:
            # It is difficult to match PairedReads outside of CASAVA conventions.
            warnings.warn(
                'Matching of paired Reads objects for file_type other than CASAVA not implemented yet.',
                UserWarning)

            return True


class PairedReads(object):
    """The C{PairedReads} class represents a pair of C{Reads} objects
    representing a read pair (i.e. I{R1} and I{R2}).

    Attributes:
    @ivar reads1: First C{Reads} object
    @type reads1: Reads
    @ivar reads2: Second C{Reads} object
    @type reads2: Reads
    @ivar read_group: SAM read group (@RG) information
    @type read_group: str
    @ivar weak_reference_sample: Weak Reference to a C{Sample}
    @type weak_reference_sample: Sample
    """

    def __init__(self, reads1=None, reads2=None, read_group=None, weak_reference_sample=None):
        """Initialise a C{PairedReads} object.

        For the C{Reads.file_type} I{CASAVA} the reads object will be
        automatically assigned on the basis of the C{Reads.read}
        attribute (i.e. I{R1} or I{R2}).
        @param reads1: First C{Reads} object
        @type reads1: Reads
        @param reads2: Second C{Reads} object
        @type reads2: Reads
        @param read_group: SAM read group (@RG) information
        @type read_group: str
        @param weak_reference_sample: Weak Reference to a C{Sample}
        @type weak_reference_sample: Sample
        @raise Exception: For C{Reads.file_type} I{CASAVA}, I{R1} or I{R2} must be set in the
        C{Reads} object.
        """

        if (reads1 and reads2) and (not reads1.match_paired(reads=reads2)):
            raise Exception('The Reads objects do not match.')

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

        if read_group:
            self.read_group = read_group
        else:
            self.read_group = str()

        if weak_reference_sample:
            self.weak_reference_sample = weak_reference_sample
        else:
            self.weak_reference_sample = None

    def trace(self, level):
        """Trace a C{PairedReads} object.

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
        output += '{}  read_group: {!r}\n'.format(indent, self.read_group)

        if isinstance(self.reads1, Reads):
            output += self.reads1.trace(level + 1)
        if isinstance(self.reads2, Reads):
            output += self.reads2.trace(level + 1)

        return output

    def add_reads(self, reads):
        """Add a C{Reads} object.

        For a C{Reads.file_type} I{CASAVA} the C{Reads} object can be automatically
        assigned on the basis of the C{Reads.read} attribute (i.e. I{R1} or I{R2}).
        @param reads: C{Reads}
        @type reads: Reads
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
                        '  reads:  {!r}'.
                        format(self.reads1.file_path, reads.file_path))
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
                        '  reads:  {!r}'.
                        format(self.reads2.file_path, reads.file_path))
                else:
                    raise Exception('Unknown Reads read attribute {!r}.'.format(reads.read))
            else:
                # Other file types go here ...
                warnings.warn(
                    'Method not implemented for file types other than CASAVA.',
                    UserWarning)

        return False

    def match(self, paired_reads):
        """Match C{PairedReads} objects.

        @param paired_reads: C{PairedReads}
        @type paired_reads: PairedReads
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
        """

        if not self.read_group == paired_reads.read_group:
            return False

        if not self.reads1.match(reads=paired_reads.reads1):
            return False

        if (self.reads2 and paired_reads.reads2) and not self.reads2.match(reads=paired_reads.reads2):
            return False

        return True

    def get_name(self, full=False):
        """Get the name of a C{PairedReads} object.

        For the C{Reads.file_type} I{CASAVA} the name is a concatenation of the
        C{Reads.name}, C{Reads.barcode} and C{Reads.lane} attributes, preferentially derived from
        the first C{Reads} object in the C{PairedReads} object.
        If the I{full} parameter is set, C{Reads.read} and C{Reads.chunk} are also added.
        @param full: Return the full name including read and chunk information
        @type full: bool
        @return: Name
        @rtype: str
        """

        if self.reads1:
            reads = self.reads1
        elif self.reads2:
            reads = self.reads2
        else:
            return

        if reads.file_type == 'CASAVA':
            if full:
                name = string.join(words=(reads.name, reads.barcode, reads.lane, reads.read, reads.chunk), sep='_')
            else:
                name = string.join(words=(reads.name, reads.barcode, reads.lane), sep='_')
        else:
            name = reads.name

        return name


class Sample(object):
    """The C{Sample} class represents a Next-Generation Sequencing sample that
    consists of one or more C{PairedReads} objects as (biological or technical) replicates
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
    @ivar paired_reads_list: Python C{list} of C{PairedReads} objects
    @type paired_reads_list: list
    @ivar weak_reference_project: Weak reference to a C{Project} object
    @type weak_reference_project: Project
    """

    default_key = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{Sample} object from a file path.

        For a I{file_type} I{CASAVA} the name is automatically populated,
        while C{PairedReads} objects are automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: Sample object
        @rtype: Sample
        """

        if file_type == 'CASAVA':

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Samples obey a "Sample_name" schema.

            components = file_name.split('_')

            name = string.join(words=components[1:], sep='_')

            sample = cls(file_path=file_path, file_type=file_type, name=name)

            # Automatically discover CASAVA Reads files ...

            for file_name in os.listdir(sample.file_path):
                file_path = os.path.join(sample.file_path, file_name)
                mode = os.stat(file_path).st_mode
                match = re.search(pattern=r'fastq.gz$', string=file_name)
                if S_ISREG(mode) and match:
                    sample.add_reads(reads=Reads.from_file_path(file_path=file_path, file_type=file_type))

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return sample

    @classmethod
    def from_samples(cls, sample1, sample2):
        """Create a merged C{Sample} from two C{Sample} objects.

        @param sample1: C{Sample}
        @type sample1: Sample
        @param sample2: C{Sample}
        @type sample2: Sample
        @return: C{Sample}
        @rtype: Sample
        """

        assert isinstance(sample1, Sample)
        assert isinstance(sample2, Sample)

        if sample1.name != sample2.name:
            warnings.warn(
                'Merged Sample objects {!r} and {!r} should have the same name.'.
                format(sample1.name, sample2.name),
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

    def __init__(self, file_path=None, file_type=None, name=None, paired_reads_list=None, weak_reference_project=None):
        """Initialise a C{Sample} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @param name: Name
        @type name: str
        @param paired_reads_list: Python C{list} of C{PairedReads} objects
        @type paired_reads_list: list
        @param weak_reference_project: Weak Reference to a C{Project} object
        @type weak_reference_project: Project
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if paired_reads_list:
            self.paired_reads_list = paired_reads_list
            # Setting the Sample object as a weak reference is problematic as it modifies the
            # PairedReads objects on the list.
            for paired_reads_object in paired_reads_list:
                paired_reads_object.weak_reference_sample = weakref.ref(self)
        else:
            self.paired_reads_list = list()

        if weak_reference_project:
            self.weak_reference_project = weak_reference_project
        else:
            self.weak_reference_project = None

    def trace(self, level):
        """Trace a C{Sample} object.

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
        output += '{}  paired_reads:\n'.format(indent)

        for paired_reads in self.paired_reads_list:
            output += paired_reads.trace(level + 1)

        return output

    def match(self, sample):
        """Match C{Sample} objects.

        @param sample: C{Sample}
        @type sample: Sample
        @return: C{True} if both objects match, C{False} otherwise
        @rtype: bool
        """

        # When Sample objects get merged, the file_path does no longer make sense.
        # Hence cope with empty file paths here ...

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

    def add_paired_reads(self, paired_reads):
        """Add a C{PairedReads} object.

        This method checks, whether a matching C{PairedReads} object is already present in this C{Sample) object.
        @param paired_reads: C{PairedReads}
        @type paired_reads: PairedReads
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

    def add_reads(self, reads):
        """Add a C{Reads} object.

        @param reads: C{Reads}
        @type reads: Reads
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

    def get_all_paired_reads(self, replicate_grouping, full=False):
        """Get all C{PairedReads} objects of a C{Sample} grouped or un-grouped.

        A C{Sample} object can hold several C{PairedReads} objects (i.e. biological or technical
        replicates) that have been sequenced on different lanes of the same flow-cell.
        The C{PairedReads} objects will therefore differ in I{name}, I{barcode} or I{lane} information.
        Depending on the I{replicate_grouping} parameter they can be returned as a group or separately.
        C{PairedReads} objects that share the I{name}, I{barcode} and I{lane}, but differ in I{chunk} number
        are always grouped together.
        @param replicate_grouping: Group all C{PairedReads} objects of a C{Sample} or
            list them individually
        @type replicate_grouping: bool
        @param full: Return the full name including read and chunk information
        @type full: bool
        @return: Python C{dict} of Python C{str} (sensible replicate name) key and
            Python C{list} object of C{PairedReads} objects value data
        @rtype: dict
        """

        groups = dict()

        for paired_reads in self.paired_reads_list:

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
    """The C{Project} class represents a Next-Generation Sequencing Project
    consisting of one or more C{Sample} objects.

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
    @ivar samples: Python C{dict} of C{Sample} objects
    @type samples: dict
    @ivar weak_reference_prf: Weak Reference to a C{ProcessedRunFolder} object
    @type weak_reference_prf: ProcessedRunFolder
    """

    default_key = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{Project} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{Project}
        @rtype: Project
        """

        if file_type == 'CASAVA':

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Projects obey a "Project_name" schema.

            components = file_name.split('_')

            name = string.join(words=components[1:], sep='_')

            project = cls(file_path=file_path, file_type=file_type, name=name)

            # Automatically discover CASAVA Sample directories ...

            for file_name in os.listdir(project.file_path):
                file_path = os.path.join(project.file_path, file_name)
                mode = os.stat(file_path).st_mode
                match = re.search(pattern=r'^Sample_(.*)$', string=file_name)
                if S_ISDIR(mode) and match:
                    if match.group(1) in project.samples:
                        raise Exception(
                            'Sample with name {!r} already exists.'.format(match.group(1)))
                    else:
                        project.add_sample(sample=Sample.from_file_path(file_path=file_path, file_type=file_type))

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return project

    def __init__(self, file_path=None, file_type=None, name=None,
                 samples=None, weak_reference_prf=None):
        """Initialise a C{Project} object.

        For a I{file_type} I{CASAVA} the name is automatically populated,
        while C{Sample} objects are automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param samples: A Python C{dict} of C{Sample} objects
        @type samples: dict
        @param weak_reference_prf: Weak Reference to a C{ProcessedRunFolder} object
        @type weak_reference_prf: ProcessedRunFolder
        @raise Exception: If C{Sample.name} values are not unique for I{file_type} I{CASAVA}
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if samples:
            self.samples = samples
        else:
            self.samples = dict()

        if weak_reference_prf:
            self.weak_reference_prf = weak_reference_prf
        else:
            self.weak_reference_prf = None

    def trace(self, level):
        """Trace a C{Project} object.

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
        output += '{}  samples:\n'.format(indent)

        for key in self.samples.keys():
            output += self.samples[key].trace(level + 1)

        return output

    def add_sample(self, sample):
        """Add a C{Sample}.

        @param sample: C{Sample}
        @type sample: Sample
        """

        if sample:
            assert isinstance(sample, Sample)
            self.samples[sample.name] = sample
            sample.weak_reference_project = weakref.ref(self)

    def get_all_samples(self):
        """Get an ordered Python C{list} of C{Sample} objects.

        @return: Python C{list} of C{Sample} objects
        @rtype: list
        """

        samples = list()

        keys = self.samples.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            samples.append(self.samples[key])

        return samples


class ProcessedRunFolder(object):
    """The C{ProcessedRunFolder} class represents an (Illumina)
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
    @ivar flow_cell: Flow-cell identifier
    @type flow_cell: str
    @ivar version: Version number
    @type version: str
    @ivar projects: Python C{dict} of C{Project} objects
    @type projects: dict
    @ivar weak_reference_collection: Weak Reference to a C{Collection} object
    @type weak_reference_collection: Collection
    """

    default_key = 'Default'

    @staticmethod
    def guess_file_type(file_path):
        """Guess the I{file_type} of a C{ProcessedRunFolder} on the basis of the I{file_path}.

        CASAVA C{ProcessedRunFolder} objects obey a I{Prefix_FCID_CASAVA182}
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

        file_name = os.path.basename(file_path.rstrip('/ '))

        components = file_name.split('_')

        match = re.match(r'^CASAVA', components[-1])

        if match:
            return 'CASAVA'
        else:
            return 'External'

    @classmethod
    def from_file_path(cls, file_path, file_type):
        """Construct a C{ProcessedRunFolder} object from a file path.

        For the I{file_type} I{CASAVA}, the C{ProcessedRunFolder.name}, C{ProcessedRunFolder.prefix},
        C{ProcessedRunFolder.flow_cell} and C{ProcessedRunFolder.version}
        attributes can be automatically parsed from the I{file_path}, while
        C{Project} objects can be automatically discovered.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{ProcessedRunFolder}
        @rtype: ProcessedRunFolder
        """

        # Try to determine the file_type if not explicitly specified.

        if not file_type or file_type == 'Automatic':
            file_type = ProcessedRunFolder.guess_file_type(file_path)

        if file_type == 'CASAVA':

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Processed Run Folders obey a "Prefix_FCID_CASAVA182"
            # schema. The following prefixes are currently in use:
            # -- BSF_ Biomedical Sequencing Facility
            # -- NGS_ Kaan Boztug group
            # -- MUW_ Medical University Vienna
            # -- SET_ Robert Kralovics group

            components = file_name.split('_')

            name = string.join(words=components[:], sep='_')
            prefix = components[0]
            flow_cell = components[1]
            version = components[2]

            prf = cls(file_path=file_path, file_type=file_type, name=name,
                      prefix=prefix, flow_cell=flow_cell, version=version)

            # Automatically discover CASAVA Project directories ...

            for file_name in os.listdir(prf.file_path):
                file_path = os.path.join(prf.file_path, file_name)
                mode = os.stat(file_path).st_mode
                match = re.search(pattern=r'^Project_(.*)$', string=file_name)
                if S_ISDIR(mode) and match:
                    if match.group(1) in prf.projects:
                        raise Exception(
                            'Project with name {!r} already exists.'.format(match.group(1)))
                    else:
                        prf.add_project(project=Project.from_file_path(file_path=file_path, file_type=file_type))

        elif file_type == 'External':

            # Create a new, minimal ProcessedRunFolder.

            name = os.path.basename(file_path.rstrip('/ '))

            prf = cls(file_path=file_path, file_type=file_type, name=name)

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return prf

    def __init__(self, file_path=None, file_type=None, name=None,
                 prefix=None, flow_cell=None, version=None,
                 projects=None, weak_reference_collection=None):
        """Initialise a C{ProcessedRunFolder} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str
        @param name: Name
        @type name: str
        @param prefix: Prefix
        @type prefix: str
        @param flow_cell: Flow-cell identifier
        @type flow_cell: str
        @param version: Version
        @type version: str
        @param projects: Python C{dict} of C{Project} objects
        @type projects: dict
        @param weak_reference_collection: Weak Reference to a C{Collection} object
        @type weak_reference_collection: Collection
        @raise Exception: If C{Project.name} values are not unique for file_type I{CASAVA}
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if prefix:
            self.prefix = prefix
        else:
            self.prefix = str()

        if flow_cell:
            self.flow_cell = flow_cell
        else:
            self.flow_cell = str()

        if version:
            self.version = version
        else:
            self.version = str()

        if projects:
            self.projects = projects
        else:
            self.projects = dict()

        if weak_reference_collection:
            self.weak_reference_collection = weak_reference_collection
        else:
            self.weak_reference_collection = None

    def trace(self, level):
        """Trace a C{ProcessedRunFolder} object.

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
        output += '{}  projects:\n'.format(indent)

        for key in self.projects.keys():
            output += self.projects[key].trace(level + 1)

        return output

    def add_project(self, project):
        """Add a C{Project} object.

        @param project: C{Project}
        @type project: Project
        """

        if project:
            assert isinstance(project, Project)
            self.projects[project.name] = project
            project.weak_reference_prf = weakref.ref(self)

    def get_all_projects(self):
        """Get an ordered Python C{list} of C{Project} objects.

        @return: A Python C{list} of C{Project} objects
        @rtype: list
        """

        projects = list()

        keys = self.projects.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            projects.append(self.projects[key])

        return projects


class Collection(object):
    """The C{Collection} class represents a collection of
    one or more C{ProcessedRunFolder} objects.

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
    @ivar processed_run_folders: Python C{dict} of C{ProcessedRunFolder} objects
    @type processed_run_folders: dict
    @ivar sample_groups: Python C{dict} of C{Sample} objects
    @type sample_groups: dict
    """

    default_key = 'Default'

    @classmethod
    def from_sas_path(cls, file_path, file_type, name, sas_path, sas_prefix=None):
        """Construct a C{Collection} from a sample annotation sheet.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param sas_path: C{SampleAnnotationSheet} file path
        @type sas_path: str | unicode
        @param sas_prefix: Optional column header prefix
            (e.g. '[Control ]Sample', '[Treatment ]Sample', ...)
        @type sas_prefix: str
        @return: C{Collection}
        @rtype: Collection
        """

        collection = cls(file_path=file_path, file_type=file_type, name=name)

        sas = SampleAnnotationSheet.from_file_path(file_path=sas_path)

        for row_dict in sas.row_dicts:
            collection._process_row_dict(row_dict=row_dict, prefix=sas_prefix)

        return collection

    def __init__(self, file_path=None, file_type=None, name=None,
                 processed_run_folders=None, sample_groups=None):
        """Initialise a C{Collection} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External}, ...)
        @type file_type: str
        @param name: Name
        @type name: str
        @param processed_run_folders: Python C{dict} of C{ProcessedRunFolder} objects
        @type processed_run_folders: dict
        @param sample_groups: Python C{dict} of group name key data and second-level Python C{dict} value data
            Second-level Python C{dict} of C{Sample.name} key data and C{Sample} object value data
        @type sample_groups: dict
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if processed_run_folders:
            self.processed_run_folders = processed_run_folders
        else:
            self.processed_run_folders = dict()

        if sample_groups:
            self.sample_groups = sample_groups
        else:
            self.sample_groups = dict()

    def trace(self, level):
        """Trace a C{Collection} object.

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
        """Add a C{ProcessedRunFolder} object.

        @param prf: C{ProcessedRunFolder}
        @type prf: ProcessedRunFolder
        """

        if prf:
            assert isinstance(prf, ProcessedRunFolder)
            self.processed_run_folders[prf.name] = prf
            prf.weak_reference_collection = weakref.ref(self)

    def get_processed_run_folder(self, file_path, file_type=None):
        """Get a C{ProcessedRunFolder} by file path.

        If the C{ProcessedRunFolder} does not exist in the C{Collection} object,
        it will be automatically discovered and added.
        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type
        @type file_type: str
        @return: C{ProcessedRunFolder}
        @rtype: ProcessedRunFolder
        """

        file_name = os.path.basename(file_path.rstrip('/ '))

        if file_name in self.processed_run_folders:
            return self.processed_run_folders[file_name]
        else:
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if not os.path.isabs(file_path):
                file_path = os.path.join(self.file_path, file_path)

            prf = ProcessedRunFolder.from_file_path(file_path=file_path, file_type=file_type)

            self.add_processed_run_folder(prf=prf)

            return prf

    def get_all_processed_run_folders(self):
        """Get an ordered Python C{list} of C{ProcessedRunFolder} objects.

        @return: Python C{list} of C{ProcessedRunFolder} objects
        @rtype: list
        """

        processed_run_folders = list()

        keys = self.processed_run_folders.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            processed_run_folders.append(self.processed_run_folders[key])

        return processed_run_folders

    def get_all_projects(self):
        """Get an ordered Python C{list} of C{Project} objects.

        @return: Python C{list} of C{Project} objects
        @rtype: list
        """

        projects = list()

        for processed_run_folder in self.get_all_processed_run_folders():
            projects.extend(processed_run_folder.get_all_projects())

        return projects

    def get_all_samples(self):
        """Get an ordered Python C{list} of C{Sample} objects.

        @return: Python C{list} of C{Sample} objects
        @rtype: list
        """

        samples = list()

        for project in self.get_all_projects():
            samples.extend(project.get_all_samples())

        return samples

    def _process_row_dict(self, row_dict, prefix=None):
        """Private method to read a hierarchy of data objects.

        This method creates C{Reads}, C{PairedReads}, C{Sample}, C{Project} and C{ProcessedRunFolder}
        objects from a C{SampleAnnotationSheet} row Python C{dict} object.
        If object-specific keys or their corresponding values are not
        available from the Python row dict, new objects will be created
        with a default key.
        Sample Annotation Sheet format:
            - FileType: Data object I{file_type} (i.e. I{CASAVA} or I{External}), defaults to I{External}.
                The FileType I{CASAVA} allows for auto-discovery of C{ProcessedRunFolder} objects.
            - ProcessedRunFolder: C{ProcessedRunFolder} file_path, can be automatically registered.
            - Project: C{Project} name
            - Sample: C{Sample} name
            - File1: C{Reads.file_name} instance variable.
                Subjected to C{os.path.expanduser} and C{os.path.expandvars}.
                If still relative, the C{Collection.file_path} is prepended.
            - Reads1: C{Reads.name} instance variable
            - File2: Same as File1
            - Reads2: Same as Reads1
            - Group: C{Sample} objects can be grouped for further analysis in
                e.g. RNA-Seq or ChIP-Seq experiments.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        """

        if not prefix:
            prefix = str()

        file_type = self._process_file_type(row_dict=row_dict, prefix=prefix)
        prf = self._process_processed_run_folder(row_dict=row_dict, prefix=prefix, file_type=file_type)
        project = self._process_project(row_dict=row_dict, prefix=prefix, file_type=file_type, prf=prf)
        sample = self._process_sample(row_dict=row_dict, prefix=prefix, file_type=file_type, project=project)
        reads1 = self._process_reads(row_dict=row_dict, prefix=prefix, file_type=file_type, suffix='1')
        reads2 = self._process_reads(row_dict=row_dict, prefix=prefix, file_type=file_type, suffix='2')
        read_group = self._process_read_group(row_dict=row_dict, prefix=prefix)

        # If none of the Reads objects has been defined, the Sample
        # may have been automatically loaded from a CASAVA ProcessedRunFolder.

        if reads1 or reads2:
            sample.add_paired_reads(paired_reads=PairedReads(reads1=reads1, reads2=reads2, read_group=read_group))

        # Optionally group the Sample objects.

        key = '{} Group'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:

            value = row_dict[key]

            if value in self.sample_groups:
                sample_group = self.sample_groups[value]
            else:
                sample_group = list()
                self.sample_groups[value] = sample_group

            if sample not in sample_group:
                sample_group.append(sample)

    def _process_file_type(self, row_dict, prefix):
        """Get file type information.

        A 'I{[Prefix] FileType}' key is optional, its value defaults to I{Automatic}.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] FileType', '[Treatment] FileType', '[Point N]  FileType', ...)
        @type prefix: str
        @return: File type
        @rtype: str
        """

        key = '{} FileType'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            file_type = row_dict[key]
        else:
            file_type = 'Automatic'

        return file_type

    def _process_processed_run_folder(self, row_dict, prefix, file_type):
        """Get or create a C{ProcessedRunFolder}.

        A 'I{[Prefix] ProcessedRunFolder}' key is optional, its value defaults to I{Default}.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] ProcessedRunFolder', '[Treatment] ProcessedRunFolder',
            '[Point N] ProcessedRunFolder', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @return: C{ProcessedRunFolder}
        @rtype: ProcessedRunFolder
        """

        key = '{} ProcessedRunFolder'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
            if value in self.processed_run_folders:
                prf = self.processed_run_folders[value]
            else:
            #    prf = ProcessedRunFolder(name=value, file_type=file_type)
            #    self.add_processed_run_folder(prf=prf)
                # Try to automatically discover a ProcessedRunFolder.
                prf = self.get_processed_run_folder(file_path=row_dict[key], file_type=file_type)
        elif ProcessedRunFolder.default_key in self.processed_run_folders:
            prf = self.processed_run_folders[ProcessedRunFolder.default_key]
        else:
            # No key or value therefore create a default object.
            prf = ProcessedRunFolder(name=ProcessedRunFolder.default_key, file_type=file_type)
            self.add_processed_run_folder(prf=prf)

        return prf

    def _process_project(self, row_dict, prefix, file_type, prf):
        """Get or create a C{Project}.

        A 'I{[Prefix] Project}' key is optional, its value defaults to I{Default}.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Project', '[Treatment] Project', '[Point N] Project', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param prf: C{ProcessedRunFolder}
        @type prf: ProcessedRunFolder
        @return: C{Project}
        @rtype: Project
        """

        key = '{} Project'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
            if value in prf.projects:
                project = prf.projects[value]
            else:
                project = Project(name=value, file_type=file_type)
                prf.add_project(project=project)
        elif Project.default_key in prf.projects:
            project = prf.projects[Project.default_key]
        else:
            project = Project(name=Project.default_key, file_type=file_type)
            prf.add_project(project=project)

        return project

    def _process_sample(self, row_dict, prefix, file_type, project):
        """Get or create a C{Sample}.

        A 'I{[Prefix] Sample}' key is optional, its value defaults to I{Default}.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param project: C{Project}
        @type project: Project
        @return: C{Sample}
        @rtype: Sample
        """

        key = '{} Sample'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
            if value in project.samples:
                sample = project.samples[value]
            else:
                sample = Sample(name=value, file_type=file_type)
                project.add_sample(sample=sample)
        elif Sample.default_key in project.samples:
            sample = project.samples[Sample.default_key]
        else:
            sample = Sample(name=Sample.default_key, file_type=file_type)
            project.add_sample(sample=sample)

        return sample

    def _process_reads(self, row_dict, prefix, file_type, suffix):
        """Get or create a C{Reads} object.

        A 'I{[Prefix] Reads{suffix}}' key is optional, in which case the default is a C{None} object.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] ReadsN', '[Treatment] ReadsN', '[Point N] ReadsN', ...)
        @type prefix: str
        @param file_type: File type
        @type file_type: str
        @param suffix: The read suffix (i.e. I{1} or I{2})
        @type suffix: str
        @return: C{Reads}
        @rtype: Reads
        """

        key1 = '{} File{}'.format(prefix, suffix).lstrip(' ')

        if key1 in row_dict and row_dict[key1]:
            file_path = row_dict[key1]
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if not os.path.isabs(file_path):
                file_path = os.path.join(self.file_path, file_path)

            key2 = '{} Reads{}'.format(prefix, suffix).lstrip(' ')

            if key2 in row_dict and row_dict[key2]:
                reads = Reads(name=row_dict[key2], file_path=file_path, file_type=file_type)
            else:
                reads = None
        else:
            reads = None

        return reads

    def _process_read_group(self, row_dict, prefix):
        """Get or create a read group.

        A 'I{[Prefix] ReadGroup}' key is optional, in which case the default is an empty Python C{str) object.
        @param row_dict: A Python C{dict} of row entries of a Python C{csv} object
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] ReadGroup', '[Treatment] ReadGroup', '[Point N] ReadGroup', ...)
        @type prefix: str
        @return: Read group string
        @rtype: str
        """

        key = '{} ReadGroup'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            read_group = row_dict[key]
        else:
            read_group = str()

        return read_group

    def get_sample_from_row_dict(self, row_dict, prefix=None):
        """Get a Sample from a C{SampleAnnotationSheet} row Python C{dict}.

        Look-up a hierarchy of C{ProcessedRunFolder}, C{Project} and C{Sample} objects
        based on a C{SampleAnnotationSheet} row dictionary.
        C{ProcessedRunFolder} objects of file type I{CASAVA} can be
        automatically discovered and registered.
        Return the corresponding C{Sample}.
        @param row_dict: C{SampleAnnotationSheet} row Python C{dict}
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample')
        @type prefix: str
        @return: C{Sample}
        @rtype: Sample
        """

        # NOTE: For the moment, the row_dict has to include keys for 'ProcessedRunFolder',
        # 'Project' and 'Sample'. Maybe, this method could search for a Sample name in
        # the Collection object. However, the problem is that a Collection contains
        # complete, auto-registered RunFolder objects. Thus, Sample names (typically 1, 2, 3, ...)
        # may clash. Therefore, it is probably best to stick to the three fields for
        # unambiguous sample resolution.

        if not prefix:
            prefix = str()

        key = '{} ProcessedRunFolder'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = ProcessedRunFolder.default_key

        # The Collection.get_processed_run_folder method can automatically register
        # ProcessedRunFolder objects of file type 'CASAVA'.

        prf = self.get_processed_run_folder(file_path=value)

        key = '{} Project'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Project.default_key

        project = prf.projects[value]

        key = '{} Sample'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = Sample.default_key

        sample = project.samples[value]

        return sample

    def get_samples_from_row_dict(self, row_dict, prefix=None):
        """Get a Python C{list} of C{Sample} objects and a Python C{str} of the Group column value as a Python C{tuple}
        from a C{SampleAnnotationSheet} row Python C{dict}.

        @param row_dict: Comparison CSV file row Python C{dict}
        @type row_dict: dict
        @param prefix: Optional configuration prefix
            (e.g. '[Control] Sample', '[Treatment] Sample', '[Point N] Sample', ...)
        @type prefix: str
        @return: Python C{tuple} of Python C{str} of '[Prefix] Group' column value and
            Python C{list} of C{Sample} objects
        @rtype: tuple
        """

        samples = list()
        value = str()

        if not prefix:
            prefix = str()

        key = '{} Group'.format(prefix).lstrip(' ')

        # Does the key exist and is its value defined?
        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Extend the Python list with all Sample objects of this group.
            if value in self.sample_groups:
                samples.extend(self.sample_groups[value])

        key = '{} Sample'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Append the Sample object, if defined.

            sample = self.get_sample_from_row_dict(row_dict=row_dict, prefix=prefix)

            if sample:
                samples.append(sample)

        return value, samples


class SampleGroup(object):
    """The C{SampleGroup} class represents a group of C{Sample} objects.

    The grouping is usually defined in a sample annotation sheet.
    Attributes:
    @ivar name: Name
    @type name: str
    @ivar samples: Python C{list} of C{Sample} objects
    @type samples: list
    """

    # TODO: The SampleGroup class is currently not in use.
    # Sample and PairedReads objects from different ProcessRunFolder objects
    # (i.e. flow-cells) could bear the same name, leading to problems with SGE job names.
    # This would need further re-thinking.

    def __init__(self, name=None, samples=None):
        """Initialise a C{SampleGroup} object.

        @param name: Name
        @type name: str
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        """

        if name:
            self.name = name
        else:
            self.name = str()

        if samples:
            self.samples = samples
        else:
            self.samples = list()

    def add_sample(self, sample):
        """Add a C{Sample} object.

        @param sample: C{Sample}
        @type sample: Sample
        """

        if sample not in self.samples:
            self.samples.append(sample)

    def get_all_paired_reads(self, replicate_grouping):
        """Get all C{PairedReads} objects of a C{SampleGroup}.

        For the moment, replicates are C{PairedReads} objects that do not share
        anything but the chunk number.
        @param replicate_grouping: Group all C{PairedReads} objects of a C{Sample} or list them individually
        @type replicate_grouping: bool
        @return: Python C{dict} of Python C{str} (replicate name) key and
            Python C{list} object of C{PairedReads} objects value data
        @rtype: dict
        """

        groups = dict()

        for sample in self.samples:

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping)

            for replicate_key in replicate_dict.keys():

                if replicate_key not in groups:
                    groups[replicate_key] = list()

                # groups[replicate_key].extend(replicate_dict[replicate_key])

                # Add PairedReads objects one-by-one and check if they are not already there.

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads not in groups[replicate_key]:
                        groups[replicate_key].append(paired_reads)

        return groups
