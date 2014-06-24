"""Bio.BSF.Data

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


import csv
import os
import re
from stat import *
import string
import warnings
import weakref


class Reads(object):
    """BSF Next-Generation Sequencing Reads class.

    The BSF Reads class represents a file of Next-Generation Sequencing reads,
    like a FASTQ or unmapped BAM file.

    Attributes:
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type
                     CASAVA: FASTQ file after post-processing with CASAVA
                     External: other data files
    :type file_type: str
    :ivar name: Name
    :type name: str
    :ivar barcode: Barcode used for sample multiplexing
    :type barcode: str
    :ivar lane: Lane number
    :type lane: str
    :ivar read: Read number (e.g. R1, R2, ...)
    :type read: str
    :ivar chunk: Chunk number (e.g. 001, 002, ...)
    :type chunk: str
    :ivar weak_reference_paired_reads: Weak Reference to a BSF PairedReads object
    :type weak_reference_paired_reads: PairedReads
    :return: Nothing
    :rtype: None
    """

    @classmethod
    def from_file_path(cls, file_path, file_type):

        """Construct a BSF Reads object from a file path.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type
        :type file_type: str
        :return: BSF Reads object
        :rtype: Reads
        """

        if file_type == 'CASAVA':

            # Since os.path.basename returns an empty string at trailing slashes,
            # use string.rstrip to remove them.

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Reads obey a Sample_name_Index_Lane_Read_Chunk schema.

            components = file_name.split('.')
            components[0] = components[0].split('_')

            name = string.join(components[0][:-4], '_')
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

        """Initialise a BSF Reads object.

        For the file_type 'CASAVA', the name, barcode, lane, read and chunk
        attributes can be automatically parsed from the file_path.
        For file_type 'External' the attributes need to be populated manually.
        :param self: BSF Reads
        :type self: Reads
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External, ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param barcode: Barcode used for sample multiplexing
        :type barcode: str
        :param lane: Lane number
        :type lane: str
        :param read: Read number (e.g. R1, R2, ...)
        :type read: str
        :param chunk: Chunk number (e.g. 001, 002, ...)
        :type chunk:str
        :param weak_reference_paired_reads: Weak Reference to a BSF PairedReads object
        :type weak_reference_paired_reads: PairedReads
        :return: Nothing
        :rtype: None
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

        """Trace a BSF Reads object.

        :param self: BSF Reads
        :type self: Reads
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
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

        """Match BSF Reads objects.

        :param self: First BSF Reads object
        :type self: Reads
        :param reads: Second BSF Reads object
        :type reads: Reads
        :return: True if both objects match, False otherwise
        :rtype: bool
        """

        assert isinstance(reads, Reads)

        if self.file_type == 'CASAVA':

            # All CASAVA Reads attributes need to be equal,
            # with the exception of the read (R1 or R2) and
            # the file_path, which also contains the reads
            # information.

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
        else:
            # It is difficult to match PairedReads outside of CASAVA conventions.
            warnings.warn(
                'Matching of paired BSF Reads objects for file_type other than CASAVA not implemented yet.',
                UserWarning)

            return True

    def match_paired(self, reads):

        """Match paired BSF Reads objects, by relaxing matching criteria.

        :param self: First BSF Reads object
        :type self: Reads
        :param reads: Second BSF Reads object
        :type reads: Reads
        :return: True if both objects match, False otherwise
        :rtype: bool
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
                'Matching of paired BSF Reads objects for file_type other than CASAVA not implemented yet.',
                UserWarning)

            return True


class PairedReads(object):
    """BSF Next-Generation Sequencing Paired Reads class.

    The BSF PairedReads class represents a pair of BSF Reads objects
    representing a read pair (i.e. R1 and R2).

    Attributes:
    :ivar reads1: First BSF Reads
    :type reads1: Reads
    :ivar reads2: Second BSF Reads
    :type reads2: Reads
    :ivar weak_reference_sample: Weak Reference to a BSF Sample
    :type weak_reference_sample: Sample
    """

    def __init__(self, reads1=None, reads2=None, weak_reference_sample=None):

        """Initialise a BSF PairedReads object.

        For the file_type 'CASAVA' the reads object will be
        automatically assigned on the basis of the Reads.read
        attribute (i.e. R1 or R2).
        :param self: BSF PairedReads
        :type self: PairedReads
        :param reads1: First BSF Reads
        :type reads1: Reads
        :param reads2: Second BSF Reads
        :type reads2: Reads
        :param weak_reference_sample: Weak Reference to a BSF Sample
        :type weak_reference_sample: Sample
        :return: Nothing
        :rtype: None
        :raise Exception: For file_type 'CASAVA', 'R1' or 'R2' must be set in the
        Reads object.
        """

        if (reads1 and reads2) and (not reads1.match_paired(reads=reads2)):
            raise Exception('The BSF Reads objects do not match.')

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
                    raise Exception('Unknown BSF Reads read attribute {!r}.'.format(reads1.read))
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
                    raise Exception('Unknown BSF Reads read attribute {!r}.'.format(reads2.read))
            else:
                # Other file types go here ...
                self.reads2 = reads2
                reads2.weak_reference_paired_reads = weakref.ref(self)

        if weak_reference_sample:
            self.weak_reference_sample = weak_reference_sample
        else:
            self.weak_reference_sample = None

    def trace(self, level):

        """Trace a BSF PairedReads object.

        :param self: BSF PairedReads
        :type self: PairedReads
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_sample: {!r}\n'.format(indent, self.weak_reference_sample)
        output += '{}  reads1: {!r}\n'.format(indent, self.reads1)
        output += '{}  reads2: {!r}\n'.format(indent, self.reads2)

        if isinstance(self.reads1, Reads):
            output += self.reads1.trace(level + 1)
        if isinstance(self.reads2, Reads):
            output += self.reads2.trace(level + 1)

        return output

    def add_Reads(self, reads):

        """Add a BSF Reads object.

        For a file_type 'CASAVA' the BSF Reads object can be automatically
        assigned on the basis of the Reads.read attribute (i.e. R1 or R2).
        :param self: BSF PairedReads
        :type self: PairedReads
        :param reads: BSF Reads
        :type reads: Reads
        :return: True upon success, False otherwise
        :rtype: bool
        """

        # For CASAVA projects, reads are automatically added to either
        # reads1 or reads2 according to the file name.
        # Returns True upon success, False otherwise.

        if not isinstance(reads, Reads):
            return False

        if self.reads1:

            if not self.reads1.match_paired(reads=reads):
                return False

            if self.reads1.file_type == 'CASAVA':
                if reads.read == 'R1':
                    raise Exception(
                        'BSF PairedReads reads1 has already been defined.\n'
                        '  reads1: {!r}\n'
                        '  reads:  {!r}'.
                        format(self.reads1.file_path, reads.file_path))
                elif reads.read == 'R2':
                    self.reads2 = reads
                    reads.weak_reference_paired_reads = weakref.ref(self)
                    return True
                else:
                    raise Exception('Unknown BSF Reads read attribute {!r}.'.format(reads.read))
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
                        'BSF PairedReads reads2 has already been defined.\n'
                        '  reads2: {!r}\n'
                        '  reads:  {!r}'.
                        format(self.reads2.file_path, reads.file_path))
                else:
                    raise Exception('Unknown BSF Reads read attribute {!r}.'.format(reads.read))
            else:
                # Other file types go here ...
                warnings.warn(
                    'Method not implemented for file types other than CASAVA.',
                    UserWarning)

        return False

    def match(self, paired_reads):

        """Match BSF PairedReads objects.

        :param self: BSF PairedReads
        :type self: PairedReads
        :param paired_reads: BSF PairedReads
        :type paired_reads: PairedReads
        :return: True if both objects match, False otherwise
        :rtype: bool
        """

        if not self.reads1.match(reads=paired_reads.reads1):
            return False

        if (self.reads2 and paired_reads.reads2) and not self.reads2.match(reads=paired_reads.reads2):
            return False

        return True

    def get_name(self, full=False):

        """Get the name of a BSF PairedReads object.

        For the file_type 'CASAVA' the name is a concatenation of the
        name, barcode and lane attributes, preferentially derived from
        the first BSF Reads object in the BSF PairedReads object.
        If the full parameter is set, read and chunk are also added.
        :param self: BSF PairedReads
        :type self: PairedReads
        :param full: Return the full name including read and chunk information.
        :type full: bool
        :return: Name string
        :rtype: str
        """

        if self.reads1:
            reads = self.reads1
        elif self.reads2:
            reads = self.reads2
        else:
            return

        if reads.file_type == 'CASAVA':
            if full:
                name = string.join([reads.name, reads.barcode, reads.lane, reads.read, reads.chunk], '_')
            else:
                name = string.join([reads.name, reads.barcode, reads.lane], '_')
        else:
            name = reads.name

        return name


class Sample(object):
    """BSF Sample class.

    The BSF Sample class represents a Next-Generation Sequencing sample that
    consists of one or more BSF PairedReads objects as (biological or technical) replicates
    that result from the same flow cell.

    Attributes:
    :cvar default_key: Default key
    :type default_key: str
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type
                     CASAVA: FASTQ file after post-processing with CASAVA
                     External: other data files
    :ivar name: Name
    :type name: str
    :ivar paired_reads_list: Python list of BSF PairedReads objects
    :type paired_reads_list: list
    :ivar weak_reference_project: Weak reference to a BSF Project object
    :type weak_reference_project: Project
    """

    default_key = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):

        """Construct a BSF Sample objects from a file path.

        For a file_type 'CASAVA' the name is automatically populated,
        while BSF PairedReads objects are automatically discovered.
        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type
        :type file_type: str
        :return: BSF Sample object
        :rtype: Sample
        """

        if file_type == 'CASAVA':

            # Since os.path.basename returns an empty string at trailing slashes,
            # use string.rstrip to remove them.

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Samples obey a "Sample_name" schema.

            components = file_name.split('_')

            name = string.join(components[1:], '_')

            sample = cls(file_path=file_path, file_type=file_type, name=name)

            # Automatically discover CASAVA Reads files ...

            for file_name in os.listdir(sample.file_path):
                file_path = os.path.join(sample.file_path, file_name)
                mode = os.stat(file_path).st_mode
                match = re.search(pattern=r'fastq.gz$', string=file_name)
                if S_ISREG(mode) and match:
                    sample.add_Reads(reads=Reads.from_file_path(file_path=file_path, file_type=file_type))

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return sample

    @classmethod
    def from_Samples(cls, sample1, sample2):

        """Create a merged BSF Sample from two BSF Sample objects.

        :param cls: Class
        :type cls: Class
        :param sample1: BSF Sample
        :type sample1: Sample
        :param sample2: BSF Sample
        :type sample2: Sample
        :return: BSF Sample object
        :rtype: Sample
        """

        assert isinstance(sample1, Sample)
        assert isinstance(sample2, Sample)

        if sample1.name != sample2.name:
            warnings.warn(
                'Merged BSF Sample objects {!r} and {!r} should have the same name.'.
                format(sample1.name, sample2.name),
                UserWarning)

        # A file_path does not make sense for merged BSF Sample objects.

        sample = cls(file_type=sample1.file_type, name=sample1.name)

        # Merge the BSF PairedReads objects from both BSF Sample objects,
        # but check, if the BSF PairedReads objects are not already there.

        for paired_reads in sample1.paired_reads_list:
            if paired_reads not in sample.paired_reads_list:
                sample.add_PairedReads(paired_reads=paired_reads)

        for paired_reads in sample2.paired_reads_list:
            if paired_reads not in sample.paired_reads_list:
                sample.add_PairedReads(paired_reads=paired_reads)

        return sample

    def __init__(self, file_path=None, file_type=None, name=None, paired_reads_list=None, weak_reference_project=None):

        """Initialise a BSF Sample object.

        :param self: BSF Sample
        :type self: Sample
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type
        :type file_type: str
        :param name: Name
        :type name: str
        :param paired_reads_list: Python list of BSF PairedReads objects
        :type paired_reads_list: list
        :param weak_reference_project: Weak Reference to a BSF Project object
        :type weak_reference_project: Project
        :return: Nothing
        :rtype: None
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
            # BSF PairedReads objects on the list.
            for paired_reads_object in paired_reads_list:
                paired_reads_object.weak_reference_sample = weakref.ref(self)
        else:
            self.paired_reads_list = list()

        if weak_reference_project:
            self.weak_reference_project = weak_reference_project
        else:
            self.weak_reference_project = None

    def trace(self, level):

        """Trace a BSF Sample object.

        :param self: BSF Sample
        :type self: Sample
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
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

        """Match BSF Sample objects.

        :param self: BSF Sample
        :type self: Sample
        :param sample: BSF Sample
        :type sample: Sample
        :return: True if both objects match, False otherwise
        :rtype: bool
        """

        # When BSF Sample objects get merged, the file_path does no longer make sense.
        # Hence cope with empty file paths here ...

        if (self.file_path and sample.file_path) and (self.file_path != sample.file_path):
            return False

        if self.file_type != sample.file_type:
            return False

        if self.name != sample.name:
            return False

        # Match Python list objects of BSF PairedReads objects ...

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

    def add_PairedReads(self, paired_reads):

        """Add a BSF Paired Reads object.

        :param self: BSF Sample
        :type self: Sample
        :param paired_reads: BSF PairedReads
        :type paired_reads: PairedReads
        :return: Nothing
        :rtype: None
        """

        if paired_reads:
            assert isinstance(paired_reads, PairedReads)
            self.paired_reads_list.append(paired_reads)
            paired_reads.weak_reference_sample = weakref.ref(self)

    def add_Reads(self, reads):

        """Add a BSF Reads object.

        :param self: BSF Sample
        :type self: Sample
        :param reads: BSF Reads
        :type reads: Reads
        :return: Nothing
        :rtype: None
        """

        assert isinstance(reads, Reads)

        # Iterate through the Python list of BSF PairedReads objects
        # The sample name must be identical.
        # The chunk must be identical
        # The read must fit

        for paired_reads in self.paired_reads_list:
            if paired_reads.add_Reads(reads=reads):
                break
        else:
            # The BSF Reads object could not be added.
            # Create a new BSF PairedReads object initialised
            # with the BSF Reads object and add it to the BSF Sample.
            self.add_PairedReads(paired_reads=PairedReads(reads1=reads))

    def get_all_PairedReads(self, replicate_grouping, full=False):

        """Get all BSF PairedReads objects of a BSF Sample grouped or un-grouped.

        A BSF Sample object can hold several BSF PairedReads objects (i.e. biological or technical
        replicates) that have been sequenced on different lanes of the same flow-cell.
        The BSF PairedReads objects will therefore differ in name, barcode or lane information.
        Depending on the replicate_grouping parameter they can be returned as a group or separately.
        BSF PairedReads objects that share the name, barcode and lane, but differ in chunk number
        are always grouped together.
        :param self: BSF Sample
        :type self: Sample
        :param replicate_grouping: Group all BSF PairedReads objects of a BSF Sample or
         list them individually
        :type replicate_grouping: bool
        :param full: Return the full name including read and chunk information.
        :type full: bool
        :return: Python dict of Python str (sensible replicate name) key and
         Python list object of BSF PairedReads objects value data
        :rtype: dict
        """

        groups = dict()

        for paired_reads in self.paired_reads_list:

            if replicate_grouping:
                # If grouped, use the Bio.BSF.Data.Sample.name as key so that all PairedReds objects
                # of this BSF Sample object end up on the same Python list object.
                key = self.name
            else:
                # If un-grouped use the Bio.BSF.Data.PairedReads.get_name() as key so that each
                # BSF PairedReads object of this BSF Sample object ends up on a separate Python list
                # object.
                key = paired_reads.get_name(full=full)

            if not key:
                continue

            if key not in groups:
                groups[key] = list()

            groups[key].append(paired_reads)

        # Return all values of the reads dictionary, which should be
        # a Python list of Python list objects of BSF PairedReads objects.

        return groups


class Project(object):
    """BSF Project class.

    The BSF Project class represents a BSF Next-Generation Sequencing Project
    consisting of one or more BSF Sample objects.

    Attributes:

    :cvar default_key: Default key
    :type default_key: str
    :ivar file_path: File path.
    :type file_path: str, unicode
    :ivar file_type: File type.
                     CASAVA: FASTQ file after post-processing with CASAVA
                     External: other data files
    :type file_type: str
    :ivar name: Name.
    :type name: str
    :ivar samples: Python dict of BSF Sample objects
    :type samples: dict
    :ivar weak_reference_prf: Weak Reference to a BSF Processed Run Folder object
    :type weak_reference_prf: ProcessedRunFolder
    """

    default_key = 'Default'

    @classmethod
    def from_file_path(cls, file_path, file_type):

        """Construct a BSF Project object from a file path.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type
        :type file_type: str
        :return: BSF Project object
        :rtype: Project
        """

        if file_type == 'CASAVA':

            # Since os.path.basename returns an empty string at trailing slashes,
            # use string.rstrip to remove them.

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Projects obey a "Project_name" schema.

            components = file_name.split('_')

            name = string.join(components[1:], '_')

            project = cls(file_path=file_path, file_type=file_type, name=name)

            # Automatically discover CASAVA Sample directories ...

            for file_name in os.listdir(project.file_path):
                file_path = os.path.join(project.file_path, file_name)
                mode = os.stat(file_path).st_mode
                match = re.search(pattern=r'^Sample_(.*)$', string=file_name)
                if S_ISDIR(mode) and match:
                    if match.group(1) in project.samples:
                        raise Exception(
                            'BSF Sample with name {!r} already exists.'.format(match.group(1)))
                    else:
                        project.add_Sample(sample=Sample.from_file_path(file_path=file_path, file_type=file_type))

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return project

    def __init__(self, file_path=None, file_type=None, name=None,
                 samples=None, weak_reference_prf=None):

        """Initialise a BSF Project object.

        For a file_type 'CASAVA' the name is automatically populated,
        while BSF Sample objects are automatically discovered.
        :param self: BSF Project
        :type self: Project
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External, ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param samples: A Python dict of BSF Sample objects
        :type samples: dict
        :param weak_reference_prf: Weak Reference to a BSF Processed Run Folder object
        :type weak_reference_prf: ProcessedRunFolder
        :return: Nothing
        :rtype: None
        :raise Exception: For file_type 'CASAVA' BSF Sample name members
        must be unique.
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

        """Trace a BSF Project object.

        :param self: BSF Project
        :type self: Project
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  weak_reference_prf: {!r}'.format(indent, self.weak_reference_prf)
        output += '{}  file_path: {!r}\n'.format(indent, self.file_path)
        output += '{}  file_type: {!r}\n'.format(indent, self.file_type)
        output += '{}  name:      {!r}\n'.format(indent, self.name)
        output += '{}  samples:\n'.format(indent)

        for key in self.samples.keys():
            output += self.samples[key].trace(level + 1)

        return output

    def add_Sample(self, sample):

        """Add a BSF Sample.

        :param self: BSF Project
        :type self: Project
        :param sample: BSF Sample
        :type sample: Sample
        :return: Nothing
        :rtype: None
        """

        if sample:
            assert isinstance(sample, Sample)
            self.samples[sample.name] = sample
            sample.weak_reference_project = weakref.ref(self)

    def get_all_Samples(self):

        """Get an ordered Python list of BSF Sample objects.

        :param self: BSF Project
        :type self: Project
        :return: Python list of BSF Sample objects
        :rtype: list
        """

        samples = list()

        keys = self.samples.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            samples.append(self.samples[key])

        return samples


class ProcessedRunFolder(object):
    """BSF (CASAVA-) Processed Run Folder class.

    The BSF ProcessedRunFolder class represents an (Illumina)
    run folder after processing with CASAVA.

    Attributes:
    :cvar default_key: Default key
    :type default_key: str
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type
                     CASAVA: FASTQ file after post-processing with CASAVA.
                     External: other data files.
    :type file_type: str
    :ivar name: Name
    :type name: str
    :ivar prefix: Prefix.
    :type prefix: str
    :ivar flow_cell: Flow-cell identifier.
    :type flow_cell: str
    :ivar version: Version number
    :type version: str
    :ivar projects: Python dict of BSF Project objects
    :type projects: dict
    :ivar weak_reference_collection: Weak Reference to a BSF Collection object
    :type weak_reference_collection: Collection
    """

    default_key = 'Default'

    @staticmethod
    def guess_file_type(file_path):

        """Guess the file_type of a BSF ProcessedRunFolder on the basis of the file_path.

        CASAVA Processed Run Folder objects obey a "Prefix_FCID_CASAVA182"
        schema. The following prefixes are currently in use:
        BSF_ Biomedical Sequencing Facility
        NGS_ Kaan Boztug group
        MUW_ Medical University Vienna
        SET_ Robert Kralovics group
        :param file_path: File path
        :type file_path: str, unicode
        :return: File type (i.e. 'CASAVA' or 'External')
        :rtype: str
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

        """Construct a BSF ProcessedRunFolder object form a file path.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type
        :type file_type: str
        :return: BSF ProcessedRunFolder object
        :rtype: ProcessedRunFolder
        """

        # Try to determine the file_type if not explicitly specified.

        if not file_type or file_type == 'Automatic':
            file_type = ProcessedRunFolder.guess_file_type(file_path)

        if file_type == 'CASAVA':

            # Since os.path.basename returns an empty string at trailing slashes,
            # use string.rstrip to remove them.

            file_name = os.path.basename(file_path.rstrip('/ '))

            # CASAVA Processed Run Folders obey a "Prefix_FCID_CASAVA182"
            # schema. The following prefixes are currently in use:
            # -- BSF_ Biomedical Sequencing Facility
            # -- NGS_ Kaan Boztug group
            # -- MUW_ Medical University Vienna
            # -- SET_ Robert Kralovics group

            components = file_name.split('_')

            name = string.join(components[:], '_')
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
                            'BSF Project with name {!r} already exists.'.format(match.group(1)))
                    else:
                        prf.add_Project(project=Project.from_file_path(file_path=file_path, file_type=file_type))

        elif file_type == 'External':

            # Create a new, minimal BSF ProcessedRunFolder.

            name = os.path.basename(file_path.rstrip('/ '))

            prf = cls(file_path=file_path, file_type=file_type, name=name)

        else:
            raise Exception('Unsupported file_type {!r}.'.format(file_type))

        return prf

    def __init__(self, file_path=None, file_type=None, name=None,
                 prefix=None, flow_cell=None, version=None,
                 projects=None, weak_reference_collection=None):

        """Initialise a BSF ProcessedRunFolder object.

        For the file_type 'CASAVA', the name, prefix, flow_cell and version
        attributes can be automatically parsed from the file_path, while
        BSF Project objects can be automatically discovered.
        :param self: BSF ProcessedRunFolder
        :type self: ProcessedRunFolder
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External or Automatic)
        :type file_type: str
        :param name: Name
        :type name: str
        :param prefix: Prefix
        :type prefix: str
        :param flow_cell: Flow-cell identifier
        :type flow_cell: str
        :param version: Version
        :type version: str
        :param projects: Python dict of BSF Project objects
        :type projects: dict
        :param weak_reference_collection: Weak Reference to a BSF Collection object
        :type weak_reference_collection: Collection
        :return: Nothing
        :rtype: None
        :raise Exception: For file_type 'CASAVA' BSF Project
        members name have to be unique.
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

        """Trace a BSF ProcessedRunFolder object.

        :param self: BSF ProcessedRunFolder
        :type self: ProcessedRunFolder
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
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

    def add_Project(self, project):

        """Add a BSF Project object.

        :param self: BSF ProcessedRunFolder
        :type self: ProcessedRunFolder
        :param project: BSF Project
        :type project: Project
        :return: Nothing
        :rtype: None
        """

        if project:
            assert isinstance(project, Project)
            self.projects[project.name] = project
            project.weak_reference_prf = weakref.ref(self)

    def get_all_Projects(self):

        """Get an ordered Python list of BSF Project objects.

        :param self: BSF ProcessedRunFolder
        :type self: ProcessedRunFolder
        :return: A Python list of BSF Project objects
        :rtype: list
        """

        projects = list()

        keys = self.projects.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            projects.append(self.projects[key])

        return projects


class Collection(object):
    """BSF Collection class.

    The BSF Collection class represents a collection of
    one or more BSF ProcessedRunFolder objects.

    Attributes:
    :cvar default_key: Default key
    :type default_key: str
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type
                     CASAVA: FASTQ file after post-processing with CASAVA
                     External: other data files
    :type file_type: str
    :ivar name: Name
    :type name: str
    :ivar processed_run_folders: Python dict of BSF ProcessedRunFolder objects
    :type processed_run_folders: dict
    :ivar sample_groups: Python dict of BSF Sample objects
    :type sample_groups: dict
    """

    default_key = 'Default'

    @classmethod
    def from_sas(cls, file_path, file_type, name, sas_file, sas_prefix=None):

        """Construct a BSF Collection from a sample annotation sheet.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External, ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param sas_file: SampleAnnotationSheet file path
        :type sas_file: str, unicode
        :param sas_prefix: Optional column header prefix
        (e.g. '[Control ]Sample', '[Treatment ]Sample', ...)
        :type sas_prefix: str
        :return: BSF Collection object
        :rtype: Collection
        """

        collection = cls(file_path=file_path, file_type=file_type, name=name)

        sas = SampleAnnotationSheet(file_path=sas_file)
        collection.read_SampleAnnotationSheet(sas=sas, prefix=sas_prefix)

        return collection

    def __init__(self, file_path=None, file_type=None, name=None,
                 processed_run_folders=None, sample_groups=None):

        """Initialise a BSF Collection object.

        :param self: BSF Collection
        :type self: Collection
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External, ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param processed_run_folders: Python dict of BSF ProcessedRunFolder objects
        :type processed_run_folders: dict
        :param sample_groups: Python dict of group name key data and second-level Python dict value data.
        Second-level Python dict of BSF Sample name key data and BSF Sample object value data.
        :type sample_groups: dict
        :return: Nothing
        :rtype: None
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

        """Trace a BSF Collection object.

        :param self: BSF Collection
        :type self: Collection
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
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

            output += '{}    {!r}\n'.format(indent, key)

            # List all BSF Sample objects of this Python list object.

            for sample in self.sample_groups[key]:
                output += '{}      {!r}: {!r}\n'.format(indent, sample.name, sample.file_path)

        return output

    def add_ProcessedRunFolder(self, prf):

        """Add a BSF ProcessedRunFolder object.

        :param self: BSF Collection
        :type self: Collection
        :param prf: BSF ProcessedRunFolder
        :type prf: ProcessedRunFolder
        :return: Nothing
        :rtype: None
        """

        if prf:
            assert isinstance(prf, ProcessedRunFolder)
            self.processed_run_folders[prf.name] = prf
            prf.weak_reference_collection = weakref.ref(self)

    def get_ProcessedRunFolder(self, file_path, file_type=None):

        """Get a BSF ProcessedRunFolder by file path.

        If the BSF ProcessedRunFolder does not exist in the BSF Collection object,
        it will be automatically discovered and added.
        :param self: BSF Collection
        :type self: Collection
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type
        :type file_type: str
        :return: BSF ProcessedRunFolder
        :rtype: ProcessedRunFolder
        """

        # Since os.path.basename returns an empty string at trailing slashes,
        # use string.rstrip to remove them.

        file_name = os.path.basename(file_path.rstrip('/ '))

        if file_name in self.processed_run_folders:
            return self.processed_run_folders[file_name]
        else:
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if not os.path.isabs(file_path):
                file_path = os.path.join(self.file_path, file_path)

            prf = ProcessedRunFolder.from_file_path(file_path=file_path, file_type=file_type)

            self.add_ProcessedRunFolder(prf=prf)

            return prf

    def get_all_ProcessedRunFolders(self):

        """Get an ordered Python list of BSF ProcessedRunFolder objects.

        :param self: BSF Collection
        :type self: Collection
        :return: Python list of BSF ProcessedRunFolder objects
        :rtype: list
        """

        processed_run_folders = list()

        keys = self.processed_run_folders.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            processed_run_folders.append(self.processed_run_folders[key])

        return processed_run_folders

    def get_all_Projects(self):

        """Get an ordered Python list of BSF Project objects.

        :param self: BSF Collection
        :type self: Collection
        :return: Python list of BSF Project objects
        :rtype: list
        """

        projects = list()

        for processed_run_folder in self.get_all_ProcessedRunFolders():
            projects.extend(processed_run_folder.get_all_Projects())

        return projects

    def get_all_Samples(self):

        """Get an ordered Python list of BSF Sample objects.

        :param self: BSF Collection
        :type self: Collection
        :return: Python list of BSF Sample objects
        :rtype: list
        """

        samples = list()

        for project in self.get_all_Projects():
            samples.extend(project.get_all_Samples())

        return samples

    def read_SampleAnnotationSheet(self, sas, prefix=None):

        """Read a BSF SampleAnnotationSheet from a file.

        This method imports a hierarchy of BSF ProcessedRunFolder, BSF Project, BSF Sample,
        BSF PairedReads and BSF Reads objects into a BSF Collection.

        :param self: BSF Collection
        :type self: Collection
        :param sas: BSF SampleAnnotationSheet
        :type sas: SampleAnnotationSheet
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] Sample', '[Treatment] Sample')
        :type prefix: str
        :return: Nothing
        :rtype: None
        """

        sas.csv_reader_open()

        for row_dict in sas._csv_reader:
            self._process_row_dict(row_dict=row_dict, prefix=prefix)

        sas.csv_reader_close()

    def _process_row_dict(self, row_dict, prefix=None):

        """Private method to read a hierarchy of BSF data objects.

        This method BSF Reads, BSF PairedReads, BSF Sample, BSF Project and BSF ProcessedRunFolder
        objects from a BSF SampleAnnotationSheet row Python dict object.
        If object-specific keys or their corresponding values are not
        available from the Python row dict, new objects will be created
        with a default key.
        Sample Annotation Sheet format:
         - FileType: BSF Data object file_type (i.e. 'CASAVA' or 'External'), defaults to 'External'.
                     The FileType 'CASAVA' allows for auto-discovery of BSF ProcessedRunFolder objects.
         - ProcessedRunFolder: BSF ProcessedRunFolder file_path, can be automatically registered.
         - Project: BSF Project name
         - Sample: BSF Sample name
         - File1: BSF Reads file_name instance variable.
                  Subjected to os.path.expanduser and os.path.expandvars.
                  If still relative, the BSF Collection file_path is prepended.
         - Reads1: BSF Reads name instance variable
         - File2: Same as File1
         - Reads2: Same as Reads1
         - Group: BSF Sample objects can be grouped for further analysis in
                  e.g. RNA-Seq or ChIP-Seq experiments.
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] Sample', '[Treatment] Sample')
        :type prefix: str
        :return: Nothing
        :rtype: None
        """

        # Get or create a BSF ProcessedRunFolder.
        # A '[Prefix ]ProcessedRunFolder' key is optional,
        # its value defaults to 'Default'.

        if not prefix:
            prefix = str()

        file_type = self._process_file_type(row_dict=row_dict, prefix=prefix)
        prf = self._process_processed_run_folder(row_dict=row_dict, prefix=prefix, file_type=file_type)
        project = self._process_project(row_dict=row_dict, prefix=prefix, file_type=file_type, prf=prf)
        sample = self._process_sample(row_dict=row_dict, prefix=prefix, file_type=file_type, project=project)
        reads1 = self._process_reads(row_dict=row_dict, prefix=prefix, file_type=file_type, suffix='1')
        reads2 = self._process_reads(row_dict=row_dict, prefix=prefix, file_type=file_type, suffix='2')

        # If none of the BSF Reads objects has been defined, the BSF Sample
        # may have been automatically loaded from a CASAVA ProcessedRunFolder.

        if reads1 or reads2:
            sample.add_PairedReads(paired_reads=PairedReads(reads1=reads1, reads2=reads2))

        # Optionally group the BSF Sample objects.

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

        A '[Prefix] FileType' key is optional, its value defaults to 'Automatic'.
        :param self: BSF Collection
        :type self: Collection
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] ProcessedRunFolder', '[Treatment] ProcessedRunFolder')
        :type prefix: str
        :return: File type
        :rtype: str
        """

        key = '{} FileType'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            file_type = row_dict[key]
        else:
            file_type = 'Automatic'

        return file_type

    def _process_processed_run_folder(self, row_dict, prefix, file_type):

        """Get or create a BSF ProcessedRunFolder.

        A '[Prefix] ProcessedRunFolder' key is optional, its value defaults to 'Default'.
        :param self: BSF Collection
        :type self: Collection
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] ProcessedRunFolder', '[Treatment] ProcessedRunFolder')
        :type prefix: str
        :param file_type: File type
        :type file_type: str
        :return: BSF ProcessedRunFolder
        :rtype: ProcessedRunFolder
        """

        key = '{} ProcessedRunFolder'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
            if value in self.processed_run_folders:
                prf = self.processed_run_folders[value]
            else:
            #    prf = ProcessedRunFolder(name=value,
            #                             file_type=file_type)
            #    self.add_ProcessedRunFolder(prf=prf)
                # Try to automatically discover a BSF ProcessedRunFolder.
                prf = self.get_ProcessedRunFolder(file_path=row_dict[key], file_type=file_type)
        elif ProcessedRunFolder.default_key in self.processed_run_folders:
            prf = self.processed_run_folders[ProcessedRunFolder.default_key]
        else:
            # No key or value therefore create a default object.
            prf = ProcessedRunFolder(name=ProcessedRunFolder.default_key, file_type=file_type)
            self.add_ProcessedRunFolder(prf=prf)

        return prf

    def _process_project(self, row_dict, prefix, file_type, prf):

        """Get or create a BSF Project.

        A '[Prefix] Project' key is optional, its value defaults to 'Default'.
        :param self: BSF Collection
        :type self: Collection
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] Project', '[Treatment] Project')
        :type prefix: str
        :param file_type: File type
        :type file_type: str
        :param prf: BSF ProcessedRunFolder
        :type prf: ProcessedRunFolder
        :return: BSF Project
        :rtype: Project
        """

        key = '{} Project'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
            if value in prf.projects:
                project = prf.projects[value]
            else:
                project = Project(name=value, file_type=file_type)
                prf.add_Project(project=project)
        elif Project.default_key in prf.projects:
            project = prf.projects[Project.default_key]
        else:
            project = Project(name=Project.default_key, file_type=file_type)
            prf.add_Project(project=project)

        return project

    def _process_sample(self, row_dict, prefix, file_type, project):

        """Get or create a BSF Sample.

        A '[Prefix] Sample' key is optional, its value defaults to 'Default'.
        :param self: BSF Collection
        :type self: Collection
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] Sample', '[Treatment] Sample')
        :type prefix: str
        :param file_type: File type
        :type file_type: str
        :param project: BSF Project
        :type project: Project
        :return: BSF Sample
        :rtype: Sample
        """

        key = '{} Sample'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
            if value in project.samples:
                sample = project.samples[value]
            else:
                sample = Sample(name=value, file_type=file_type)
                project.add_Sample(sample=sample)
        elif Sample.default_key in project.samples:
            sample = project.samples[Sample.default_key]
        else:
            sample = Sample(name=Sample.default_key, file_type=file_type)
            project.add_Sample(sample=sample)

        return sample

    def _process_reads(self, row_dict, prefix, file_type, suffix):

        # Get or create a BSF Reads object.
        # The '[Prefix ]Reads{suffix}' and '[Prefix ]File{suffix}' keys are optional,
        # in which case the default is a None object.

        """Get or create a BSF Reads object.

        A '[Prefix] Reads{suffix}' key is optional, in which case the default is a None object.
        :param self: BSF Collection
        :type self: Collection
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] ReadsN', '[Treatment] ReadsN')
        :type prefix: str
        :param file_type: File type
        :type file_type: str
        :param suffix: The read suffix (i.e. 1 or 2)
        :type suffix: str
        :return: BSF Reads
        :rtype: Reads
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

    def get_Sample_from_row_dict(self, row_dict, prefix=None):

        """Get a BSF Sample from a BSF SampleAnnotationSheet row Python dict.

        Look-up a hierarchy of BSF ProcessedRunFolder, BSF Project and BSF Sample objects
        based on a BSF SampleAnnotationSheet row dictionary.
        BSF ProcessedRunFolder objects of file type 'CASAVA' can be
        automatically discovered and registered.
        Return the corresponding BSF Sample.
        :param row_dict: BSF SampleAnnotationSheet row Python dict
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] Sample', '[Treatment] Sample')
        :type prefix: str
        :return: BSF Sample
        :rtype: Sample
        """

        # NOTE: For the moment, the row_dict has to include keys for 'ProcessedRunFolder',
        # 'Project' and 'Sample'. Maybe, this method could search for a BSF Sample name in
        # the BSF Collection object. However, the problem is that a BSF Collection contains
        # complete, auto-registered BSF RunFolder objects. Thus, Sample names (typically 1, 2, 3, ...)
        # may clash. Therefore, it is probably best to stick to the three fields for
        # unambiguous sample resolution.

        if not prefix:
            prefix = str()

        key = '{} ProcessedRunFolder'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]
        else:
            value = ProcessedRunFolder.default_key

        # The Collection.get_ProcessedRunFolder method can automatically register
        # BSF ProcessedRunFolder objects of file type 'CASAVA'.

        prf = self.get_ProcessedRunFolder(file_path=value)

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

    def get_Samples_from_row_dict(self, row_dict, prefix=None):

        """Get a Python list of BSF Sample objects and a Python str of the Group column value as a Python tuple
        from a BSF SampleAnnotationSheet row Python dict.

        :param row_dict: Comparison CSV file row Python dict
        :type row_dict: dict
        :param prefix: Optional configuration prefix
                       (e.g. '[Control] Sample', '[Treatment] Sample')
        :type prefix: str
        :return: Tuple of Python str of '[Prefix] Group' column value and
         Python list of BSF Sample objects
        :rtype: tuple
        """

        samples = list()
        value = str()

        if not prefix:
            prefix = str()

        key = '{} Group'.format(prefix).lstrip(' ')

        # Does the key exist and is its value defined?
        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Extend the Python list with all BSF Sample objects of this group.
            if value in self.sample_groups:
                samples.extend(self.sample_groups[value])

        key = '{} Sample'.format(prefix).lstrip(' ')

        if key in row_dict and row_dict[key]:
            value = row_dict[key]

            # Append the BSF Sample object, if defined.

            sample = self.get_Sample_from_row_dict(row_dict=row_dict, prefix=prefix)

            if sample:
                samples.append(sample)

        return value, samples


class AnnotationSheet(object):
    """BSF Annotation Sheet class.

    The BSF AnnotationSheet class represents comma-separated value (CSV) files.

    Attributes:
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type (e.g. ...)
    :type file_type: str
    :ivar name: Name
    :type name: str
    :ivar field_names: Python list of field names
    :type field_names: list
    :ivar row_dicts: Python list of Python dict objects
    :type row_dicts: list
    """

    non_alpha_expression = re.compile(pattern='\W')
    non_numeric_expression = re.compile(pattern='\D')
    non_sequence_expression = re.compile(pattern='[^ACGTacgt]')

    @classmethod
    def check_column(cls, row_number, row_dict, column_name, require_column=True, require_value=True):

        """Check for a column name and return its associated value, if any.

        :param cls: Class
        :type cls: Class
        :param row_number: Row number for warning messages
        :type row_number: int
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param column_name: Column name
        :type column_name: str
        :param require_column: Require the column_name to be defined in teh row_dict
        :type require_column: bool
        :param require_value: Require a value
        :type require_value: bool
        :return: Tuple of warning message and column value
        :rtype: tuple
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

        """Validate a particular column value as alphanumeric.

         Check that the particular column name key exists in the row dictionary and that
         its associated value contains only alphanumeric characters.

        :param cls: Class
        :type cls: Class
        :param row_number: Row number for warning messages
        :type row_number: int
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param column_name: Column name
        :type column_name: str
        :return: Warning messages
        :rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            match = re.search(pattern=cls.non_alpha_expression, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with non-alphanumeric characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_numeric(cls, row_number, row_dict, column_name):

        """Validate a particular column value as numeric.

         Check that the particular column name key exists in the row dictionary and that
         its associated value contains only numeric characters.

        :param cls: Class
        :type cls: Class
        :param row_number: Row number for warning messages
        :type row_number: int
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param column_name: Column name
        :type column_name: str
        :return: Warning messages
        :rtype: str
        """

        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name)

        if column_value:
            match = re.search(pattern=cls.non_numeric_expression, string=row_dict[column_name])
            if match:
                messages += 'Column {!r} in row {} contains a value {!r} with non-numeric characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_sequence(cls, row_number, row_dict, column_name, require_column=True, require_value=True):

        """Validate a particular column value as sequence.

        :param cls: Class
        :type cls: Class
        :param row_number: Row number for warning messages
        :type row_number: int
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param column_name: Column name
        :type column_name: str
        :param require_column: Require the column_name to be defined in teh row_dict
        :type require_column: bool
        :param require_value: Require a value
        :type require_value: bool
        :return: Warning messages
        :rtype: str
        """
        messages, column_value = cls.check_column(row_number=row_number, row_dict=row_dict, column_name=column_name,
                                                  require_column=require_column, require_value=require_value)

        if column_value:
            match = re.search(pattern=cls.non_sequence_expression, string=row_dict[column_name])
            if match:
                messages += 'Field {!r} in row {} contains a sequence {!r} with illegal characters.\n'. \
                    format(column_name, row_number, row_dict[column_name])

        return messages

    @classmethod
    def check_sequence_mandatory(cls, row_number, row_dict, column_name):

        """Validate a particular column value as mandatory sequence.

        :param cls: Class
        :type cls: Class
        :param row_number: Row number for warning messages
        :type row_number: int
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param column_name: Column name
        :type column_name: str
        :return: Warning messages
        :rtype: str
        """
        return cls.check_sequence(row_number=row_number, row_dict=row_dict, column_name=column_name,
                                  require_column=True, require_value=True)

    @classmethod
    def check_sequence_optional(cls, row_number, row_dict, column_name):

        """Validate a particular column value as optional sequence.

        :param cls: Class
        :type cls: Class
        :param row_number: Row number for warning messages
        :type row_number: int
        :param row_dict: A Python dict of row entries of a Python csv object
        :type row_dict: dict
        :param column_name: Column name
        :type column_name: str
        :return: Warning messages
        :rtype: str
        """
        return cls.check_sequence(row_number=row_number, row_dict=row_dict, column_name=column_name,
                                  require_column=True, require_value=False)

    @classmethod
    def read_from_file(cls, file_path=None, file_type=None, name=None):

        """Construct an Annotation Sheet from a comma-separated value (CSV) file.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. ...)
        :type file_type: str
        :return: BSF Annotation Sheet
        :rtype: AnnotationSheet
        """

        annotation_sheet = cls(file_path=file_path, file_type=file_type, name=name)

        annotation_sheet.csv_reader_open()

        for row_dict in annotation_sheet._csv_reader:
            annotation_sheet.row_dicts.append(row_dict)

        annotation_sheet.csv_reader_close()

        return annotation_sheet

    def __init__(self, file_path=None, file_type=None, name=None, field_names=None, row_dicts=None):

        """Initialise a BSF AnnotationSheet object.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param field_names: Python list of field names
        :type field_names: list
        :param row_dicts: Python list of Python dict objects
        :type row_dicts: list
        :return: Nothing
        :rtype: None
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

        if field_names:
            self.field_names = field_names
        else:
            self.field_names = list()

        if row_dicts:
            self.row_dicts = row_dicts
        else:
            self.row_dicts = list()

    def csv_reader_open(self, header=True):

        """Open a Comma-Separated Value (CSV) file for reading and initialise a Python csv.DictWriter object.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :param header: A header line is present and populates the filed_names Python list.
        :type header: bool
        :return: Nothing
        :rtype: None
        """

        # Although the AnnotationSheet is initialised with an empty Python list object,
        # the DictReader really needs None to automatically populate the fieldnames instance variable.
        # However, the DictReader should always do so, if a header line is expected since
        # otherwise, the header line would get repeated as data line.

        if len(self.field_names) and not header:
            field_names = self.field_names
        else:
            field_names = None

        self._csv_reader_file = open(name=self.file_path, mode='rb')
        self._csv_reader = csv.DictReader(f=self._csv_reader_file, fieldnames=field_names)

        # Automatically set the field names from the DictReader,
        # if the field_names list is empty and if possible.

        if len(self._csv_reader.fieldnames) and not len(self.field_names):
            self.field_names.extend(self._csv_reader.fieldnames)

    def csv_reader_next(self):

        """Read the next line of a CSV file.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Python dict of column key and row value data
        :rtype: dict
        """

        return self._csv_reader.next()

    def csv_reader_close(self):

        """Close a Comma-Separated Value (CSV) file for reading.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Nothing
        :rtype: None
        """

        self._csv_reader = None
        self._csv_reader_file.close()
        self._csv_reader_file = None

    def read_tsv(self):

        """Read a Tab-Separated Value (TSV) file.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Nothing
        :rtype: None
        """

        self._tsv_file = open(name=self.file_path, mode='r')
        self._tsv_rows = list()

        for line in self._tsv_file:

            # Exclude comment lines.

            if re.search(pattern=r'\s*#', string=line):
                continue

            self._tsv_rows.append(line.split())

        self._tsv_file.close()

    def csv_writer_open(self):

        """Open a Comma-Separated Value (CSV) file for writing and initialise a Python csv.DictWriter object.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Nothing
        :rtype: None
        """

        self._csv_writer_file = open(name=self.file_path, mode='wb')
        self._csv_writer = csv.DictWriter(f=self._csv_writer_file, fieldnames=self.field_names)
        self._csv_writer.writeheader()

    def csv_writer_next(self, row_dict):

        """Write the next line of a CSV file.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :param row_dict: Row dictionary
        :type row_dict: dict
        :return: Nothing
        :rtype: None
        """

        self._csv_writer.writerow(rowdict=row_dict)

    def csv_writer_close(self):

        """Close a Comma-Separated Value (CSV) file for writing.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Nothing
        :rtype: None
        """

        self._csv_writer = None
        self._csv_writer_file.close()

    def sort(self):

        """Sort a BSF Annotation Sheet.
        This method has to implemented in the sub-class,
        as it requires information about field-specific sorting.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Nothing
        :rtype: None
        """

        warnings.warn(
            'Sorting of BSF Annotation Sheets has to implemented in the sub-class.',
            UserWarning)

    def validate(self, test_methods):

        """Validate a BSF Annotation Sheet.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :param test_methods: Dict of column names and associated test methods,
        i.e. AnnotationSheet.check_alphanumeric, AnnotationSheet.check_sequence, ...
        :return: Warning messages
        :rtype: str
        """

        messages = str()
        row_number = 0

        for row_dict in self.row_dicts:
            row_number += 1
            for field_name in self.field_names:
                if field_name in test_methods:
                    # Only validate fields, for which instructions exist in the tests dict.
                    messages += test_methods[field_name](row_number=row_number,
                                                         row_dict=row_dict,
                                                         column_name=field_name)

        return messages

    def write_to_file(self):

        """Write a BSF Annotation Sheet to a file.

        :param self: BSF AnnotationSheet
        :type self: AnnotationSheet
        :return: Nothing
        :rtype: None
        """

        self.csv_writer_open()

        for row_dict in self.row_dicts:
            self.csv_writer_next(row_dict=row_dict)

        self.csv_writer_close()


class BamIndexDecoderSheet(AnnotationSheet):
    """
    BSF BamIndexDecoder Sheet class.

    The BSF BamIndexDecoder Sheet class represents a comma-separated value (CSV) table of
    library information for the Illumina2bam BamIndexDecoder.
    """

    # Column header names.

    field_names = ['lane', 'barcode_sequence_1', 'barcode_sequence_2', 'sample_name', 'library_name']

    # Python dict to correlate columns with validation methods.

    test_methods = dict(lane=AnnotationSheet.check_alphanumeric,
                        barcode_sequence_1=AnnotationSheet.check_sequence_mandatory,
                        barcode_sequence_2=AnnotationSheet.check_sequence_optional,
                        sample_name=AnnotationSheet.check_alphanumeric,
                        library_name=AnnotationSheet.check_alphanumeric)

    @classmethod
    def read_from_file(cls, file_path=None, file_type=None, name=None):

        """Construct a BamIndexDecoder Sheet from a comma-separated value (CSV) file.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. ...)
        :type file_type: str
        :return: BSF BamIndexDecoder Sheet
        :rtype: BamIndexDecoderSheet
        """

        return super(BamIndexDecoderSheet, cls).read_from_file(file_path=file_path)

    def __init__(self, file_path=None, file_type=None, name=None, row_dicts=None):

        """Initialise a BSF BamIndexDecoder Sheet object.

        :param self: BSF BamIndexDecoder Sheet
        :type self: BamIndexDecoderSheet
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param row_dicts: Python list of Python dict objects
        :type row_dicts: list
        :return: Nothing
        :rtype: None
        """

        super(BamIndexDecoderSheet, self).__init__(file_path=file_path,
                                                   file_type=file_type,
                                                   name=name,
                                                   field_names=BamIndexDecoderSheet.field_names,
                                                   row_dicts=row_dicts)

    def validate(self, test_methods=None):

        """
        Validate a BSF BamIndexDecoder Sheet.

        :param test_methods:
        :type test_methods: dict
        :return: Warning messages
        :rtype: str
        """

        messages = str()

        # Check the header line via the pre-defined field names.

        for index in range(0, len(BamIndexDecoderSheet.field_names)):
            if not self.field_names[index]:
                messages += 'Column with name {!r} is missing from the header line.\n'. \
                    format(BamIndexDecoderSheet.field_names[index])

            if self.field_names[index] != BamIndexDecoderSheet.field_names[index]:
                messages += 'Column name {!r} in the header line does not match template {!r}.\n'. \
                    format(self.field_names[index], BamIndexDecoderSheet.field_names[index])

        # Validate the field values for alphanumeric or sequence grade in the context of the
        # BSF Annotation Sheet super-class.

        if not test_methods:
            test_methods = BamIndexDecoderSheet.test_methods

        messages += super(BamIndexDecoderSheet, self).validate(test_methods=test_methods)

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

            if 'barcode_sequence_2' in row_dict and row_dict['barcode_sequence_2']:
                barcode_sequence += row_dict['barcode_sequence_2']

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

        # TODO: The lane number needs to be configurable. Maybe this needs checking at the point where the
        # Annotation Sheet gets used?
        for lane_number in range(1, 9):
            lane_string = str(lane_number)
            if lane_string not in lane_index:
                messages += 'No annotation for lane number {!r}.\n'.format(lane_number)

        return messages


class SampleAnnotationSheet(AnnotationSheet):
    pass


class SampleGroup(object):
    """BSF Sample Group class.

    The BSF Sample Group class represents a group of BSF Sample objects.

    The grouping is usually defined in a sample annotation sheet.
    Attributes:
    :ivar name: Name
    :type name: str
    :ivar samples: Python list of BSF Sample objects
    :type samples: list
    """

    # TODO: The BSF SampleGroup class is currently not in use.
    # BSF Sample and BSF PairedReads objects from different BSF ProcessRunFolder objects
    # (i.e. flow-cells) could bear the same name, leading to problems with SGE job names.
    # This would need further re-thinking.

    def __init__(self, name=None, samples=None):

        """Initialise a BSF SampleGroup object.

        :param self: BSF SampleGroup
        :type self: SampleGroup
        :param name: Name
        :type name: str
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :return: Nothing
        :rtype: None
        """

        if name:
            self.name = name
        else:
            self.name = str()

        if samples:
            self.samples = samples
        else:
            self.samples = list()

    def add_Sample(self, sample):

        """Add a BSF Sample object.

        :param self: BSF SampleGroup
        :type self: SampleGroup
        :param sample: BSF Sample
        :type sample: Sample
        :return: Nothing
        :rtype: None
        """

        if sample not in self.samples:
            self.samples.append(sample)

    def get_all_PairedReads(self, replicate_grouping):

        """Get all BSF PairedReads objects of a BSF SampleGroup.

        For the moment, replicates are BSF PairedReads objects that do not share
        anything but the chunk number.
        :param self: BSF SampleGroup
        :type self: SampleGroup
        :param replicate_grouping: Group all BSF PairedReads objects of a BSF Sample or list them individually
        :type replicate_grouping: bool
        :return: Python dict of Python str (replicate name) key and
         Python list object of BSF PairedReads objects value data
        :rtype: dict
        """

        groups = dict()

        for sample in self.samples:

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            for replicate_key in replicate_dict.keys():

                if replicate_key not in groups:
                    groups[replicate_key] = list()

                # groups[replicate_key].extend(replicate_dict[replicate_key])

                # Add BSF PairedReads objects one-by-one and check if they are not already there.

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads not in groups[replicate_key]:
                        groups[replicate_key].append(paired_reads)

        return groups
