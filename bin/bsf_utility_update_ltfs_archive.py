#!/usr/bin/env python
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
"""The :py:mod:`bin.bsf_utility_update_ltfs_archive` module is a script to
collect and update :emphasis:`Linear Tape File System` (LTFS) archive contents.

A directory tree (:literal:`directory_path`) is scanned for files matching a wildcard pattern
(:literal:`--pattern`), by default :literal:`*_log.txt` and updated into an existing
(:literal:`--file-path`) content file.
The LTFS content file lists run folder names and a comma-separated list of LTO barcode labels
in :emphasis:`tab-separated value` (TSV) format.
"""

import json
import logging
import os
import re
import sys
from argparse import ArgumentParser
from typing import Optional


def read_ltfs_archive(archive_file_path: str) -> dict[str, list[str]]:
    """Read a :emphasis:`Linear Tape File System` archive file in :emphasis:`tab-separated value` (TSV) format.

    :param archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
    :type archive_file_path: str
    :return: A Python :py:class:`dict` object of
        Python :py:class:`str` (file name) key and
        Python :py:class:`list` object of
        Python :py:class:`str` (cartridge barcode) objects.
    :rtype: dict[str, list[str]]
    """
    ltfs_dict: dict[str, list[str]] = dict()

    if os.path.exists(archive_file_path):
        with open(file=archive_file_path, mode='rt') as text_io:
            for line_str in text_io:
                line_list = line_str.split()
                if line_list[0] not in ltfs_dict:
                    ltfs_dict[line_list[0]] = list()
                tape_list = ltfs_dict[line_list[0]]
                for tape_id in line_list[1].split(','):
                    if tape_id not in tape_list:
                        tape_list.append(tape_id)

    return ltfs_dict


def write_ltfs_archive(archive_file_path: str, ltfs_dict: dict[str, list[str]]) -> None:
    """Write an :emphasis:`Linear Tape File System` archive file in :emphasis:`tab-separated value` (TSV) format.

    :param archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
    :type archive_file_path: str
    :param ltfs_dict: A Python :py:class:`dict` object of
        Python :py:class:`str` (file name) key and
        Python :py:class:`list` object of
        Python :py:class:`str` (cartridge barcode) objects.
    :type ltfs_dict: dict[str, list[str]]
    """
    with open(file=archive_file_path, mode='wt') as text_io:
        for irf_file_name in sorted(ltfs_dict):
            print(irf_file_name, ','.join(sorted(ltfs_dict[irf_file_name])), sep='\t', file=text_io)

    return


def process_json_files(ltfs_dict: dict[str, list[str]], directory_path: str, json_pattern: str) -> None:
    """Process JSON files with LTFS :emphasis:`virtual extended attributes` (VEA).

    :param ltfs_dict: A Python :py:class:`dict` object of
        Python :py:class:`str` (file name) key and
        Python :py:class:`list` object of
        Python :py:class:`str` (cartridge barcode) objects.
    :type ltfs_dict: dict[str, list[str]]
    :param directory_path: A directory path.
    :type directory_path: str
    :param json_pattern: A regular expression pattern matching JSON file names.
    :type json_pattern: str | None
    """
    re_pattern = re.compile(pattern=json_pattern)

    for base_path, directory_name_list, file_name_list in os.walk(top=directory_path):
        logging.debug('directory_path: %r', base_path)
        logging.debug('directory_name_list: %r', directory_name_list)
        logging.debug('file_name_list: %r', file_name_list)

        for file_name in file_name_list:
            re_match = re_pattern.search(string=file_name)
            if re_match is None:
                logging.debug('Excluding file_name: %r', file_name)
                continue

            # The tape name could be parsed from the file name via a regular expression, or better, retrieved from
            # the corresponding LTFS virtual extended attribute (VEA).
            # volume_name = re_match.group(1)
            with open(file=os.path.join(base_path, file_name), mode='rt') as text_io:
                vea_dict = json.load(fp=text_io)

                volume_name = vea_dict['volume']['ltfs.volumeName']

                for vea_object_dict in vea_dict['objects']:
                    file_name = vea_object_dict['file.path']

                    if file_name not in ltfs_dict:
                        ltfs_dict[file_name] = list()

                    tape_list = ltfs_dict[file_name]

                    if volume_name not in tape_list:
                        tape_list.append(volume_name)

    return


def process_ltfscp_files(ltfs_dict: dict[str, list[str]], directory_path: str, ltfscp_pattern: str) -> None:
    """Process :emphasis:`LTFSCP` output files, particularly the :literal:`ILT30505I` informational entry.

    :param ltfs_dict: A Python :py:classy:`dict` object of
        Python :py:class:`str` (file name) key and
        Python :py:class:`list` object of
        Python :py:class:`str` (cartridge barcode) objects.
    :type ltfs_dict: dict[str, list[str]]
    :param directory_path: A directory path.
    :type directory_path: str
    :param ltfscp_pattern: A regular expression pattern matching :emphasis:`LTFSCP` output file names.
    :type ltfscp_pattern: str
    """
    re_pattern = re.compile(pattern=ltfscp_pattern)

    for base_path, directory_name_list, file_name_list in os.walk(top=directory_path):
        logging.debug('directory_path: %r', base_path)
        logging.debug('directory_name_list: %r', directory_name_list)
        logging.debug('file_name_list: %r', file_name_list)

        for file_name in file_name_list:
            re_match = re_pattern.search(string=file_name)
            if re_match is None:
                logging.debug('Excluding file_name: %r', file_name)
                continue

            # Parse the volume name from the file name (e.g., BS0048L6_log.txt)
            volume_name = re_match.group(1)

            with open(file=os.path.join(base_path, file_name), mode='rt') as input_file:
                for line_str in input_file:
                    # ILT30505I Copy /scratch/lab_bsf/archive_irf/BS0048L6/BSF_0648_H7WC2BBXY_3.bam to
                    #   /mnt/ltfs/BSF_0648_H7WC2BBXY_3.bam
                    if not line_str.startswith('ILT30505I'):
                        continue

                    line_list = line_str.split()
                    base_name = line_list[2].split('/')[-1]

                    logging.debug('File name: %r', base_name)

                    if base_name not in ltfs_dict:
                        ltfs_dict[base_name] = list()

                    tape_list = ltfs_dict[base_name]

                    if volume_name not in tape_list:
                        tape_list.append(volume_name)

    return


def run(
        directory_path: Optional[str] = None,
        archive_file_path: Optional[str] = None,
        json_pattern: Optional[str] = None,
        ltfscp_pattern: Optional[str] = None) -> int:
    """Run function.

    :param directory_path: A :emphasis:`Linear Tape File System` output directory path.
    :type directory_path: str | None
    :param archive_file_path: A :emphasis:`Linear Tape File System` archive file path.
    :type archive_file_path: str | None
    :param json_pattern: A regular expression pattern matching JSON file names.
    :type json_pattern: str | None
    :param ltfscp_pattern: A regular expression pattern matching :emphasis:`LTFSCP` output file names.
    :type ltfscp_pattern: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    # Read the initial LTFS content file that needs updating.
    ltfs_dict = read_ltfs_archive(archive_file_path=archive_file_path)

    # Process JSON files if a pattern was provided.
    if json_pattern:
        process_json_files(ltfs_dict=ltfs_dict, directory_path=directory_path, json_pattern=json_pattern)

    # Process LTFSCP files if a pattern was provided.
    if ltfscp_pattern:
        process_ltfscp_files(ltfs_dict=ltfs_dict, directory_path=directory_path, ltfscp_pattern=ltfscp_pattern)

    # Write the updated LTFS content file back to its original location.
    write_ltfs_archive(archive_file_path=archive_file_path, ltfs_dict=ltfs_dict)

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Update an LTFS archive file')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--json-pattern',
        default=r'^([0-9A-Z]{6,6}L[5-6])\.json$',
        help=r'regular expression matching JSON file names [^([0-9A-Z]{6,6}L[5-6])\.json$]')

    argument_parser.add_argument(
        '--ltfscp-pattern',
        default=r'^([0-9A-Z]{6,6}L[5-6])_log\.txt$',
        help=r'regular expression matching LTFSCP file names [^([0-9A-Z]{6,6}L[5-6])_log\.txt$]')

    argument_parser.add_argument(
        '--file-path',
        help='LTFS archive file path')

    argument_parser.add_argument(
        'directory_path',
        help='LTFS output directory path')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return run(
        directory_path=name_space.directory_path,
        archive_file_path=name_space.file_path,
        json_pattern=name_space.json_pattern,
        ltfscp_pattern=name_space.ltfscp_pattern)


if __name__ == '__main__':
    sys.exit(main())
