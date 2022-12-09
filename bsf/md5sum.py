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
"""The :py:mod:`bsf.md5sum` module supports the :emphasis:`GNU md5sum` utility.
"""
import logging
import os
import re
from argparse import ArgumentParser
from tempfile import NamedTemporaryFile
from typing import Optional, TypeVar

from bsf.process import Executable

MD5SumType = TypeVar(name='MD5SumType', bound='MD5Sum')
MD5SumArchiveType = TypeVar(name='MD5SumArchiveType', bound='MD5SumArchive')

module_logger = logging.getLogger(name=__name__)


class MD5Sum(object):
    """The :py:class:`bsf.md5sum.MD5Sum` object represents one :emphasis:`GNU md5sum` entry.

    :ivar file_path: A file path.
    :type file_path: str
    :ivar check_sum: A check sum.
    :type check_sum: str
    :ivar check_mode: A check mode.
    :type check_mode: str
    """

    def __init__(self, file_path: str, check_sum: str, check_mode: str) -> None:
        """Initialise a :py:class:`bsf.md5sum.MD5Sum` object.

        :param file_path: A file path.
        :type file_path: str
        :param check_sum: A check sum.
        :type check_sum: str
        :param check_mode: A check mode.
        :type check_mode: str
        """
        self.file_path = file_path
        self.check_sum = check_sum
        self.check_mode = check_mode

        return

    @classmethod
    def from_line(cls, md5sum_str: str) -> MD5SumType:
        """Create a :py:class:`bsf.md5sum.MD5Sum` object from a :emphasis:`GNU md5sum` line.

        :param md5sum_str: A :emphasis:`GNU md5sum` line.
        :type md5sum_str: str
        :return: A :py:class:`bsf.md5sum.MD5Sum` object.
        :rtype: MD5Sum
        """
        md5sum_str = md5sum_str.strip()

        if ' ' in md5sum_str:
            index_int = md5sum_str.index(' ')

            # The MD5 checksum lies up until the index location.
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

    def to_line(self) -> str:
        """Convert a :py:class:`bsf.md5sum.MD5Sum` object to a :emphasis:`GNU md5sum` line.

        :return: A :emphasis:`GNU md5sum` line.
        :rtype: str
        """
        return ' '.join((self.check_sum, self.check_mode + self.file_path))


class MD5SumArchive(object):
    """The :py:class:`bsf.md5sum.MD5SumArchive` models a file of :emphasis:`GNU md5sum` entries.

    :ivar file_path: A MD5 checksum archive file path.
    :type file_path: str
    :ivar md5sum_dict: A Python :py:class:`dict` object of
        Python :py:class:`str` (file name) key and
        :py:class:`bsf.md5sum.MD5Sum` value objects.
    :type md5sum_dict: dict[str, MD5Sum]
    """

    def __init__(self, file_path: Optional[str] = None, md5sum_dict: Optional[dict[str, MD5Sum]] = None) -> None:
        """Initialise a :py:class:`bsf.md5sum.MD5SumArchive` object.

        :param file_path: A MD5 checksum archive file path.
        :type file_path: str | None
        :param md5sum_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (file name) key and
            :py:class:`bsf.md5sum.MD5Sum` value objects.
        :type md5sum_dict: dict[str, MD5Sum] | None
        """
        self.file_path = file_path

        if md5sum_dict is None:
            self.md5sum_dict = dict()
        else:
            self.md5sum_dict = md5sum_dict

        return

    def add_md5sum(self, md5sum: MD5Sum) -> bool:
        """Add a :py:class:`bsf.md5sum.MD5Sum` object to the :py:class:`bsf.md5sum.MD5SumArchive` object.

        :param md5sum: A :py:class:`bsf.md5sum.MD5Sum` object.
        :type md5sum: MD5Sum
        :return: :py:const:`True` upon success, :py:const:`False` otherwise.
        :rtype: bool
        """
        # Check if a file_path instance variable is set.
        if not md5sum.file_path:
            module_logger.warning("Missing 'file_path' instance variable for 'check_sum' %r.", md5sum.check_sum)
            return False

        # Check if a check_sum instance variable is set.
        if not md5sum.check_sum:
            module_logger.warning("Missing 'check_sum' instance variable for 'file_path' %r.", md5sum.file_path)
            return False

        # Check if no other check_sum value is set for this file_path.
        if md5sum.file_path in self.md5sum_dict and self.md5sum_dict[md5sum.file_path].check_sum != md5sum.check_sum:
            module_logger.warning('Non-matching check sum %r for file path %r.', md5sum.check_sum, md5sum.file_path)
            return False

        self.md5sum_dict[md5sum.file_path] = md5sum

        return True

    def add_md5sum_line(self, md5sum_str: str) -> bool:
        """Split a GNU :literal:`md5sum` line into its components and add it to the
        :py:class:`bsf.md5sum.MD5SumArchive` object.

        :param md5sum_str: A :emphasis:`GNU md5sum` line.
        :type md5sum_str: str
        :return: :py:const:`True` upon success, :py:const:`False` otherwise.
        :rtype: bool
        """
        return self.add_md5sum(md5sum=MD5Sum.from_line(md5sum_str=md5sum_str))

    def read_md5sum_archive(self, file_path: Optional[str] = None) -> None:
        """Read an MD5 checksum archive file.

        :param file_path: File path
        :type file_path: str | None
        """
        if not file_path:
            file_path = self.file_path

        if os.path.exists(file_path):
            with open(file=file_path, mode='rt') as input_text_io:
                for line_str in input_text_io:
                    md5sum = MD5Sum.from_line(md5sum_str=line_str)

                    module_logger.debug(
                        'read_md5sum_archive: MD5Sum.file_path: %r MD5Sum.check_sum: %r MD5Sum.check_mode: %r',
                        md5sum.file_path, md5sum.check_sum, md5sum.check_mode)

                    if not self.add_md5sum(md5sum=md5sum):
                        raise Exception(
                            f"The MD5SumArchive file {file_path!r} "
                            "does not obey the standard '<MD5_SUM> *<file_path>' format.")
        else:
            module_logger.warning('The MD5SumArchive.file_path %r does not exist.', file_path)

        return

    @classmethod
    def from_file_path(cls, file_path: str) -> MD5SumArchiveType:
        """Create a :py:class:`bsf.md5sum.MD5SumArchive` from a file path.

        :param file_path: A file path.
        :type file_path: str
        :return: A :py:class:`bsf.md5sum.MD5SumArchive` object.
        :rtype: MD5SumArchive
        """
        module_logger.debug('Reading MD5SumArchive from file path: %r', file_path)

        if not file_path:
            raise Exception("Require a valid 'file_path' option.")

        md5sum_archive = cls(file_path=file_path)
        md5sum_archive.read_md5sum_archive()

        return md5sum_archive

    def to_file_path(self, file_path: Optional[str] = None) -> None:
        """Write a :py:class:`bsf.md5sum.MD5SumArchive` object to a file path.

        :param file_path: A file path.
        :type file_path: str | None
        """
        if not file_path:
            file_path = self.file_path

        module_logger.debug('Writing MD5SumArchive to file path: %r', file_path)

        with open(file=file_path, mode='wt') as output_text_io:
            for md5_file_name in sorted(self.md5sum_dict):
                md5sum = self.md5sum_dict[md5_file_name]

                # Adjust the mode to binary for certain files.
                for suffix in ('.bam', '.gz'):
                    if md5sum.file_path.endswith(suffix):
                        md5sum.check_mode = '*'

                print(md5sum.to_line(), file=output_text_io)

        return

    @classmethod
    def console_update(
            cls,
            directory_path: Optional[str] = None,
            file_path: Optional[str] = None,
            file_pattern: Optional[str] = None) -> int:
        """Console function to collect and update a :py:class:`bsf.md5sum.MD5SumArchive` object
        with a collection of :emphasis:`GNU md5sum` files.

        A directory tree (:literal:`directory_path`) is scanned for file names matching a regular expression pattern
        (:literal:`file_pattern`), by default :literal:`\\.md5` and updated into an existing
        (:literal:`file_path`) :py:class:`bsf.md5sum.MD5SumArchive` file.

        Picard-style MD5 files that only contain the MD5 digest are reformatted to obey the
        :emphasis:`GNU md5sum` format.

        :param directory_path: A :emphasis:`GNU md5sum` directory path.
        :type directory_path: str | None
        :param file_path: A :emphasis:`GNU md5sum` archive file path.
        :type file_path: str | None
        :param file_pattern: A file name regular expression pattern.
        :type file_pattern: str | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        md5sum_archive: MD5SumArchive = MD5SumArchive.from_file_path(file_path=file_path)

        re_pattern = re.compile(pattern=file_pattern)

        for directory_path, directory_name_list, file_name_list in os.walk(top=directory_path):
            logging.debug('directory_path: %r', directory_path)
            logging.debug('directory_name_list: %r', directory_name_list)
            logging.debug('file_name_list: %r', file_name_list)

            for file_name in file_name_list:
                if re_pattern.search(string=file_name) is None:
                    logging.debug('Excluding file_name: %r', file_name)
                    continue

                with open(file=os.path.join(directory_path, file_name), mode='rt') as text_io:
                    for line_str in text_io:
                        md5sum = MD5Sum.from_line(md5sum_str=line_str)

                        if not md5sum.file_path:
                            # In case the md5sum file does not specify a file name (e.g., Picard MD5 checksum files),
                            # use the file name without its '.md5' suffix.
                            if file_name.endswith('.md5'):
                                md5sum.file_path = file_name[:-4]
                            else:
                                raise Exception(f'The (Picard) MD5 file suffix is unsupported: {file_name!r}')

                        # Archive just the base name of the file path.
                        md5sum.file_path = os.path.basename(md5sum.file_path)

                        md5sum_archive.add_md5sum(md5sum=md5sum)

        md5sum_archive.to_file_path()

        return 0

    @classmethod
    def entry_point_update(cls) -> int:
        """Console entry point to collect and update a :py:class:`bsf.md5sum.MD5SumArchive` object
        with a collection of :emphasis:`GNU md5sum` files.

        A directory tree (:literal:`directory_path`) is scanned for file names matching a regular expression pattern
        (:literal:`--file-pattern`), by default :literal:`\\.md5` and updated into an existing
        (:literal:`--file-path`) :py:class:`bsf.md5sum.MD5SumArchive` file.

        Picard-style MD5 files that only contain the MD5 digest are reformatted to obey the
        :emphasis:`GNU md5sum` format.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description='Update a GNU md5sum archive file')

        argument_parser.add_argument(
            '--logging-level',
            default='WARNING',
            choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
            help='logging level [WARNING]')

        argument_parser.add_argument(
            '--file-path',
            required=True,
            help='GNU md5sum archive file path')

        argument_parser.add_argument(
            '--file-pattern',
            default=r'\.md5$',
            help=r'file name regular expression pattern [\.md5$]')

        argument_parser.add_argument(
            'directory_path',
            help='GNU md5sum directory path')

        name_space = argument_parser.parse_args()

        if name_space.logging_level:
            logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
            logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

            logging.basicConfig(level=name_space.logging_level)

        return cls.console_update(
            directory_path=name_space.directory_path,
            file_path=name_space.file_path,
            file_pattern=name_space.file_pattern)

    @classmethod
    def console_check(
            cls,
            md5sum_archive_path: Optional[str] = None,
            file_path_list: Optional[list[str]] = None) -> int:
        """Console function to check files against a :py:class:`bsf.md5sum.MD5SumArchive`
        :emphasis:`GNU md5sum` file.

        :param md5sum_archive_path: A :emphasis:`GNU md5sum` file path.
        :type md5sum_archive_path: str | None
        :param file_path_list: A Python :py:class:`list`´object of Python :py:class:`str` (file path) objects.
        :type file_path_list: list[str] | None
        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        md5sum_archive = MD5SumArchive.from_file_path(file_path=md5sum_archive_path)

        text_io = NamedTemporaryFile(mode='wt', suffix='.md5', delete=False)

        for file_path in file_path_list:
            file_name = os.path.basename(file_path)
            if file_name in md5sum_archive.md5sum_dict:
                md5sum = md5sum_archive.md5sum_dict[file_name]
                md5_file_path = os.path.normpath(os.path.join(os.getcwd(), file_path))
                print(md5sum.check_sum + ' ' + md5sum.check_mode + md5_file_path, file=text_io)

        text_io.close()

        executable = Executable(name='md5sum', program='md5sum')
        executable.add_switch_long(key='check')
        executable.arguments.append(text_io.name)

        message_list = executable.run()

        if message_list:
            print(repr(message_list))

        os.remove(text_io.name)

        print('All done.')

        return 0

    @classmethod
    def entry_point_check(cls) -> int:
        """Console entry point to check files against a :py:class:`bsf.md5sum.MD5SumArchive`
        :emphasis:`GNU md5sum` file.

        :return: A :py:class:`SystemExit` status value.
        :rtype: int
        """
        argument_parser = ArgumentParser(
            description='Check one or more files against a GNU md5sum archive file')

        argument_parser.add_argument(
            '--logging-level',
            default='WARNING',
            choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
            help='logging level [WARNING]')

        argument_parser.add_argument(
            '--md5sum-path',
            required=True,
            help='GNU md5sum archive file path')

        argument_parser.add_argument(
            'file_path',
            nargs='+',
            help='one or more file paths',
            dest='file_path_list')

        name_space = argument_parser.parse_args()

        if name_space.logging_level:
            logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
            logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

            logging.basicConfig(level=name_space.logging_level)

        return cls.console_check(
            md5sum_archive_path=name_space.md5sum_path,
            file_path_list=name_space.file_path_list)
