#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
#
# BSF Python script to support IBM Linear Tape File System (LTFS) operations.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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
#

from __future__ import print_function

import argparse
import os
import xml.etree.ElementTree

cartridge_dict = {
    # IBM LTO Ultrium Cartridge Label Specification (Revision 6)
    # Part Number 19P0034
    # EC - M10321
    # http://www-01.ibm.com/support/docview.wss?uid=ssg1S7000429
    # 'L1': # Generation 1 Type A (100 GB)
    # 'LA': # Generation 1 Type B (50 GB)
    # 'LB': # Generation 1 Type C (30 GB)
    # 'LC': # Generation 1 Type D (10 GB)
    # 'L2': # Generation 2 Type A (200 GB)
    # 'L3': # Generation 3 Type A (400 GB)
    # 'L4': # Generation 4 Type A (600 GB)
    'L5': 1391601152 * 1024,  # Generation 5 Type A (1500 GB)
    'L6': 2351648768 * 1024,  # Generation 6 Type A (2500 GB)
}

si_prefix_dict = {
    0: 'B',
    1: 'kiB',
    2: 'MiB',
    3: 'GiB',
    4: 'TiB',
    5: 'PiB',
    6: 'EiB',
}


def convert_into_readable(integer_bytes):
    """Convert an integer number of bytes into a readable number with SI prefixes.

    @param integer_bytes: Integer number of bytes
    @type integer_bytes: int
    @return: Python C{tuple} of Python C{float} readable bytes and Python C{str} (SI) unit.
    @rtype: (float, str)
    """
    si_index = 0
    readable_bytes = float(integer_bytes)
    while readable_bytes > 1024.0:
        readable_bytes /= 1024.0
        si_index += 1

    return readable_bytes, si_prefix_dict[si_index]


class LinearTapeFileSystemDirectory(object):
    """The C{LinearTapeFileSystemDirectory} class represents a source (and target) directory in the copy process.

    Attributes:
    @ivar source_path: Source directory path
    @type source_path: str | unicode | None
    @ivar target_path: Target directory path
    @type target_path: str | unicode | None
    @ivar source_specification: Source specification pattern including wildcard characters
    @type source_specification: str | unicode | None
    @ivar source_file_path_list: Python C{list} of source file paths in the directory
    @type source_file_path_list: list[str | unicode]
    """

    def __init__(
            self,
            source_path=None,
            target_path=None,
            source_specification=None,
            source_file_path_list=None):
        """Initialise a C{LinearTapeFileSystemDirectory} object.

        @param source_path: Source directory path
        @type source_path: str | unicode | None
        @param target_path: Target directory path
        @type target_path: str | unicode | None
        @param source_specification: Source specification pattern including wildcard characters
        @type source_specification: str | unicode | None
        @param source_file_path_list: Python C{list} of source file paths in the directory
        @type source_file_path_list: list[str | unicode] | None
        @return:
        @rtype:
        """
        super(LinearTapeFileSystemDirectory, self).__init__()

        self.source_path = source_path
        self.target_path = target_path
        self.source_specification = source_specification

        if source_file_path_list is None:
            self.source_file_path_list = list()
        else:
            self.source_file_path_list = source_file_path_list

        return

    def add_source_file_path(self, source_file_path):
        """Add a source file path to a C{LinearTapeFileSystemDirectory} object.

        @param source_file_path: Source file path
        @type source_file_path: str | unicode
        @return:
        @rtype:
        """
        if source_file_path is None:
            return

        if source_file_path not in self.source_file_path_list:
            self.source_file_path_list.append(source_file_path)

        return


class LinearTapeFileSystemCopy(object):
    """The C{LinearTapeFileSystemCopy} class represents one Linear Tape File System Copy (ltfscp) process.

    Attributes:
    @ivar total_buffer_size: Total buffer size
    @type total_buffer_size: str | None
    @ivar buffer_size: Buffer size per thread
    @type buffer_size: str | None
    @ivar log_level: Log level
    @type log_level: str | None
    @ivar recursive: Recursive processing
    @type recursive: bool | None
    @ivar sparse: Support sparse files
    @type sparse: bool | None
    @ivar default_target_path: Default target path for all C{LinearTapeFileSystemDirectory} objects
    @type default_target_path: str | unicode | None
    @ivar ltfs_directory_dict: Python C{dict} of Python C{str} directory path and
        C{LinearTapeFileSystemDirectory} objects
    @type ltfs_directory_dict: dict[str | unicode, LinearTapeFileSystemDirectory]

    XML batch file format
    http://www-01.ibm.com/support/knowledgecenter/STZMZN/com.ibm.storage.hollywood.doc/ltfs_reference_lcp_cli.html

    <?xml version="1.0" encoding="UTF-8"?>
    <ltfscpspec version="1.0">
     <params>
      [parameters]
     </params>
     <data>
      <file>
       <srcpath>path</srcpath>
       <dstpath>path</dstpath>
       <srcspec>expression</srcspec>
       <sf>filename</sf>
       <sf rename="newname">filename</sf>
      </file>
      <file>
       ...
      </file>
     </data>
    </ltfscpspec>
    """

    def __init__(
            self,
            total_buffer_size=None,
            buffer_size=None,
            log_level=None,
            recursive=None,
            sparse=None,
            default_target_path=None,
            ltfs_directory_dict=None):
        """Initialise a C{LinearTapeFileSystemCopy} object.

        @param total_buffer_size: Total buffer size
        @type total_buffer_size: str | None
        @param buffer_size: Buffer size per thread
        @type buffer_size: str | None
        @param log_level: Log level
        @type log_level: str | None
        @param recursive: Recursive processing
        @type recursive: bool | None
        @param sparse: Support sparse files
        @type sparse: bool | None
        @param default_target_path: Default target path for all C{LinearTapeFileSystemDirectory} objects
        @type default_target_path: str | unicode | None
        @param ltfs_directory_dict: Python C{dict} of Python C{str} or C{unicode} directory path and
            C{LinearTapeFileSystemDirectory} objects
        @type ltfs_directory_dict: dict[str | unicode, LinearTapeFileSystemDirectory] | None
        @return:
        @rtype:
        """
        super(LinearTapeFileSystemCopy, self).__init__()

        self.total_buffer_size = total_buffer_size
        self.buffer_size = buffer_size
        self.log_level = log_level
        self.recursive = recursive
        self.sparse = sparse
        self.default_target_path = default_target_path

        if ltfs_directory_dict is None:
            self.ltfs_directory_dict = dict()
        else:
            self.ltfs_directory_dict = ltfs_directory_dict

        return

    def get_or_add_ltfs_directory(self, ltfs_directory):
        """Add a C{LinearTapeFileSystemDirectory} object to the C{LinearTapeFileSystemCopy} object,
        if it does not exist already.

        @param ltfs_directory: C{LinearTapeFileSystemDirectory}
        @type ltfs_directory: LinearTapeFileSystemDirectory
        @return: C{LinearTapeFileSystemDirectory}
        @rtype: LinearTapeFileSystemDirectory
        """
        if ltfs_directory is None:
            return

        if ltfs_directory.source_path in self.ltfs_directory_dict:
            return self.ltfs_directory_dict[ltfs_directory.source_path]
        else:
            self.ltfs_directory_dict[ltfs_directory.source_path] = ltfs_directory
            return ltfs_directory

    def get_or_add_source_directory(self, source_path, source_specification=None):
        """Add a source directory to the C{LinearTapeFileSystemCopy} object and return the
        corresponding C{LinearTapeFileSystemDirectory} object.

        @param source_path: Source directory path
        @type source_path: str | unicode
        @param source_specification: Source specification
        @type source_specification: str | unicode | None
        @return: C{LinearTapeFileSystemDirectory}
        @rtype: LinearTapeFileSystemDirectory
        """
        if source_path in self.ltfs_directory_dict:
            return self.ltfs_directory_dict[source_path]
        else:
            ltfs_directory = LinearTapeFileSystemDirectory(
                source_path=source_path,
                source_specification=source_specification)
            self.ltfs_directory_dict[ltfs_directory.source_path] = ltfs_directory
            return ltfs_directory

    def add_source_file_path(self, source_path):
        """Add a source file path to a C{LinearTapeFileSystemCopy} object.

        @param source_path: Source file path
        @type source_path: str | unicode
        @return:
        @rtype:
        """
        source_path = os.path.normpath(source_path)
        source_directory = os.path.dirname(source_path)
        source_file_path = os.path.basename(source_path)

        ltfs_directory = self.get_or_add_source_directory(source_path=source_directory)
        ltfs_directory.add_source_file_path(source_file_path=source_file_path)

        return

    def get_element_tree(self):
        """Get the C{xml.etree.ElementTree.ElementTree} representation of the C{LinearTapeFileSystemCopy} object.

        @return: XML Element Tree
        @rtype: xml.etree.ElementTree.ElementTree
        """
        ltfs_specification = xml.etree.ElementTree.Element(tag='ltfscpspec', attrib={'version': '1.0'})

        ltfs_element_tree = xml.etree.ElementTree.ElementTree(element=ltfs_specification)

        ltfs_parameters = xml.etree.ElementTree.Element(tag='params')
        ltfs_specification.append(element=ltfs_parameters)

        if self.total_buffer_size:
            param_total_buffer_size = xml.etree.ElementTree.Element(tag='total-buffer-size')
            param_total_buffer_size.text = self.total_buffer_size
            ltfs_parameters.append(element=param_total_buffer_size)

        if self.buffer_size:
            param_buffer_size = xml.etree.ElementTree.Element(tag='buffer-size')
            param_buffer_size.text = self.buffer_size
            ltfs_parameters.append(element=param_buffer_size)

        if self.log_level:
            param_log_level = xml.etree.ElementTree.Element(tag='loglevel')
            param_log_level.text = self.log_level
            ltfs_parameters.append(element=param_log_level)

        if self.recursive:
            param_recursive = xml.etree.ElementTree.Element(tag='recursive')
            param_recursive.text = 'enable'
            ltfs_parameters.append(element=param_recursive)

        if self.sparse:
            param_sparse = xml.etree.ElementTree.Element(tag='sparse')
            param_sparse.text = 'enable'
            ltfs_parameters.append(element=param_sparse)

        ltfs_data = xml.etree.ElementTree.Element(tag='data')
        ltfs_specification.append(element=ltfs_data)

        for source_path in sorted(self.ltfs_directory_dict):
            ltfs_directory = self.ltfs_directory_dict[source_path]

            ltfs_file = xml.etree.ElementTree.Element(tag='file')
            ltfs_data.append(element=ltfs_file)

            ltfs_file_source_path = xml.etree.ElementTree.Element(tag='srcpath')
            ltfs_file_source_path.text = ltfs_directory.source_path
            ltfs_file.append(element=ltfs_file_source_path)

            ltfs_file_destination_path = xml.etree.ElementTree.Element(tag='dstpath')
            if ltfs_directory.target_path:
                ltfs_file_destination_path.text = ltfs_directory.target_path
            else:
                ltfs_file_destination_path.text = self.default_target_path
            ltfs_file.append(element=ltfs_file_destination_path)

            if ltfs_directory.source_specification:
                ltfs_file_source_specification = xml.etree.ElementTree.Element(tag='srcspec')
                ltfs_file_source_specification.text = ltfs_directory.source_specification
                ltfs_file.append(element=ltfs_file_source_specification)

            ltfs_directory.source_file_path_list.sort()

            for source_file_path in ltfs_directory.source_file_path_list:
                ltfs_source_file = xml.etree.ElementTree.Element(tag='sf')
                ltfs_source_file.text = source_file_path
                ltfs_file.append(element=ltfs_source_file)

        return ltfs_element_tree

    def write_batch_file(self, file_path):
        """Write a Linear Tape File System Copy tool batch file.

        @param file_path: File path
        @type file_path: str | unicode
        @return:
        @rtype:
        """
        ltfs_element_tree = self.get_element_tree()
        ltfs_element_tree.write(file_or_filename=file_path)


argument_parser = argparse.ArgumentParser(
    description='Linear Tape File System Copy tool driver script.')

argument_parser.add_argument(
    'cartridge',
    help='cartridge barcode',
    type=str)

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--total-buffer-size',
    default='4G',
    dest='total_buffer_size',
    help='The default total buffer size is 400M. The maximum total buffer size is 4G. [4G]',
    required=False,
    type=str)

argument_parser.add_argument(
    '--buffer-size',
    default='1G',
    dest='buffer_size',
    help='The default buffer size is 128K. The maximum buffer size is 1G. [1G]',
    required=False,
    type=str)

argument_parser.add_argument(
    '--no-sparse',
    action='store_false',
    dest='sparse',
    help='do not support sparse files')

argument_parser.add_argument(
    '--recursive',
    action='store_true',
    dest='recursive',
    help='process recursively')

argument_parser.add_argument(
    '--log-level',
    choices=['URGENT', 'WARNING', 'INFO', 'DEBUG'],
    default='INFO',
    dest='log_level',
    help='log level [INFO]',
    required=False,
    type=str)

argument_parser.add_argument(
    '--target-path',
    default='/mnt/ltfs',
    dest='target_path',
    help='default target path [/mnt/ltfs]',
    required=False,
    type=str)

argument_parser.add_argument(
    '--source-specification',
    default='',
    dest='source_specification',
    help='source specification pattern []',
    required=False,
    type=str)

argument_parser.add_argument(
    '--mounted-cartridge',
    action='store_true',
    dest='mounted_cartridge',
    help='calculate remaining space for a mounted cartridge')

name_space = argument_parser.parse_args()

linear_tape_file_system_copy = LinearTapeFileSystemCopy(
    total_buffer_size=name_space.total_buffer_size,
    buffer_size=name_space.buffer_size,
    log_level=name_space.log_level,
    recursive=name_space.recursive,
    sparse=name_space.sparse,
    default_target_path=name_space.target_path)

cartridge_code = name_space.cartridge
""" @type cartridge_code: str """

if cartridge_code.endswith('.txt'):
    # In case shell completion is used, remove the trailing *.txt suffix.
    cartridge_code = name_space.cartridge[:-4]

if cartridge_code[-2:] not in cartridge_dict:
    raise Exception('The cartridge barcode media characters ' + repr(cartridge_code[-2:]) +
                    ' are currently not supported.')

if name_space.mounted_cartridge:
    # If a cartridge is mounted, stat the virtual file system.
    ltfs_statvfs_result = os.statvfs(linear_tape_file_system_copy.default_target_path)
    ltfs_free_bytes = ltfs_statvfs_result.f_bavail * ltfs_statvfs_result.f_bsize
    ltfs_total_bytes = ltfs_statvfs_result.f_blocks * ltfs_statvfs_result.f_bsize
else:
    ltfs_free_bytes = cartridge_dict[cartridge_code[-2:]]
    ltfs_total_bytes = cartridge_dict[cartridge_code[-2:]]

# Summarize the sizes of all source files.

total_size = 0

with open(file=cartridge_code + '.txt', mode='rt') as input_file:
    for main_file_path in input_file:
        main_file_path = main_file_path.strip()
        main_file_path = os.path.normpath(main_file_path)
        ltfs_stat_result = os.stat(main_file_path)
        total_size += ltfs_stat_result.st_size
        linear_tape_file_system_copy.add_source_file_path(source_path=main_file_path)

if ltfs_free_bytes <= total_size:
    difference_bytes = total_size - ltfs_free_bytes
    difference_readable, difference_unit = convert_into_readable(integer_bytes=difference_bytes)
    print('LTFS total size {:,d} LTFS free size {:,d} Total file size {:,d} Exceeding {:,d} ({:0.1f} {})'.format(
        ltfs_total_bytes, ltfs_free_bytes, total_size, difference_bytes, difference_readable, difference_unit))
else:
    difference_bytes = ltfs_free_bytes - total_size
    difference_readable, difference_unit = convert_into_readable(integer_bytes=difference_bytes)
    print('LTFS total size {:,d} LTFS free size {:,d} Total file size {:,d} Remaining {:,d} ({:0.1f} {})'.format(
        ltfs_total_bytes, ltfs_free_bytes, total_size, difference_bytes, difference_readable, difference_unit))

linear_tape_file_system_copy.write_batch_file(file_path=cartridge_code + '.xml')
