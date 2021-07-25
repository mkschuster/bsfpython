# -*- coding: utf-8 -*-
"""MD5 sum module.

A package of classes and methods modelling a repository of GNU md5sum information.
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
import logging
import os


class MD5Sum(object):
    """The C{MD5Sum} object represents one GNU md5sum entry.

    @ivar file_path: File path
    @type file_path: str
    @ivar check_sum: Check sum
    @type check_sum: str
    @ivar check_mode: Check mode
    @type check_mode: str
    """

    def __init__(self, file_path, check_sum, check_mode):
        """Initialise a C{MD5Sum} object.

        @param file_path: File path
        @type file_path: str
        @param check_sum: Check sum
        @type check_sum: str
        @param check_mode: Check mode
        @type check_mode: str
        """
        self.file_path = file_path
        self.check_sum = check_sum
        self.check_mode = check_mode

        return

    @classmethod
    def from_line(cls, md5sum_str):
        """Create a C{MD5Sum} object from a GNU md5sum line.

        @param md5sum_str: GNU md5sum line
        @type md5sum_str: str
        @return: C{MD5Sum} object
        @rtype: MD5Sum
        """
        md5sum_str = md5sum_str.strip()

        if ' ' in md5sum_str:
            index_int = md5sum_str.index(' ')

            # The MD5 check sum lies up until the index location.
            check_sum = md5sum_str[:index_int]

            # The check mode marker for text (' ') or binary ('*') lies after the index location.
            check_mode = md5sum_str[index_int + 1:index_int + 2]

            # The file name excluding the marker follows the index location + 1.
            file_path = md5sum_str[index_int + 2:]
        else:
            check_sum = md5sum_str
            check_mode = '*'
            file_path = ''

        return cls(file_path=file_path, check_sum=check_sum, check_mode=check_mode)

    def to_line(self):
        """Convert a C{MD5Sum} object to a GNU md5sum line.

        @return: GNU md5sum line
        @rtype: str
        """
        return ' '.join((self.check_sum, self.check_mode + self.file_path))


class MD5SumArchive(object):
    """The C{MD5SumArchive} models a file of GNU md5sum entries.

    @ivar file_path: File path
    @type file_path: str
    @ivar md5sum_dict: Python C{dict} of file name Python C{str} key and C{MD5Sum} value objects
    @type md5sum_dict: dict[str, MD5Sum]
    """

    def __init__(self, file_path=None, md5sum_dict=None):
        """Initialise a C{MD5SumArchive} object.

        @param file_path: File path of the MD5 sum repository
        @type file_path: str | None
        @param md5sum_dict: Python C{dict} of file name Python C{str} key and C{MD5Sum} value objects
        @type md5sum_dict: dict[str, MD5Sum] | None
        """
        self.file_path = file_path

        if md5sum_dict is None:
            self.md5sum_dict = dict()
        else:
            self.md5sum_dict = md5sum_dict

        return

    def add_md5sum(self, md5sum):
        """Add a C{MD5Sum} object to the C{MD5SumArchive} object.

        @param md5sum: C{MD5Sum} object
        @type md5sum: MD5Sum
        @return: True upon success, False otherwise
        @rtype: bool
        """
        # Check if a file_path instance variable is set.
        if not md5sum.file_path:
            logging.warning("Missing file_path instance variable for check sum '%s'", md5sum.check_sum)
            return False

        # Check if no other check_sum value is set for this file_path.
        if md5sum.file_path in self.md5sum_dict and self.md5sum_dict[md5sum.file_path].check_sum != md5sum.check_sum:
            logging.warning("Non-matching check sum '%s' for file path '%s'", md5sum.check_sum, md5sum.file_path)
            return False

        self.md5sum_dict[md5sum.file_path] = md5sum

        return True

    def add_md5sum_line(self, md5sum_str):
        """Split a GNU md5sum line into its components and add it to the C{MD5SumArchive} object.

        @param md5sum_str: GNU md5sum line
        @type md5sum_str: str
        @return: True upon success, False otherwise
        @rtype: bool
        """
        return self.add_md5sum(md5sum=MD5Sum.from_line(md5sum_str=md5sum_str))

    def read_md5sum_archive(self, file_path=None):
        """Read a MD5 sum archive file.

        @param file_path: File path
        @type file_path: str | None
        """
        if not file_path:
            file_path = self.file_path

        if os.path.exists(file_path):
            with open(file=file_path, mode='rt') as input_file:
                for line_str in input_file:
                    md5sum = MD5Sum.from_line(md5sum_str=line_str)

                    logging.debug("md5_file_path:  '%s'", md5sum.file_path)
                    logging.debug("md5_check_sum:  '%s'", md5sum.check_sum)
                    logging.debug("md5_check_mode: '%s'", md5sum.check_mode)

                    if not self.add_md5sum(md5sum=md5sum):
                        raise Exception('The md5sum archive file ' + repr(file_path) +
                                        " does not obey the standard 'MD5SUM *file_path' format.")
        else:
            logging.warning("The md5sum archive does not exists: '%s'", file_path)

        return

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{MD5SumArchive} from a file path.

        @param file_path: File path
        @type file_path: str
        @return: C{MD5SumArchive}
        @rtype: MD5SumArchive
        """
        logging.debug("Reading MD5SumArchive from file path: '%s'", file_path)

        if not file_path:
            raise Exception('Require a valid file_path option.')

        md5sum_archive = cls(file_path=file_path)
        md5sum_archive.read_md5sum_archive()

        return md5sum_archive

    def to_file_path(self, file_path=None):
        """Write a C{MD5SumArchive} object to a file path.

        @param file_path: File path
        @type file_path: str | None
        """
        if not file_path:
            file_path = self.file_path

        logging.debug("Writing MD5SumArchive to file path: '%s'", file_path)

        with open(file=file_path, mode='wt') as text_io:
            for md5_file_name in sorted(self.md5sum_dict):
                md5sum = self.md5sum_dict[md5_file_name]

                # Adjust the mode to binary for certain files.
                for suffix in ('.bam', '.gz'):
                    if md5sum.file_path.endswith(suffix):
                        md5sum.check_mode = '*'

                print(md5sum.to_line(), file=text_io)

        return
