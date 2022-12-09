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
"""The :py:mod:`bin.bsf_utility_fix_sample_links` module is a script to
find and fix broken symbolic links pointing from the BSF samples to the
BSF sequences directory. These symbolic links are set by the
:py:class:`bsf.analyses.illumina_to_bam_tools.BamIndexDecoder` analysis
in case a flow-cell lane does not need de-multiplexing.
"""

import errno
import os
import stat
import sys
from argparse import ArgumentParser
from typing import Optional


def scan_directory(samples_path: str, sequences_path: str, commit: Optional[bool] = None) -> int:
    """Scan a directory for symbolic links or directories recursively.

    :param samples_path: A samples archive directory path
    :type samples_path: str
    :param sequences_path: A sequence archive directory path
    :type sequences_path: str
    :param commit: Commit changes to the file system
    :type commit: bool | None
    :return: A :py:class:`SystemExit` status value
    :rtype: int
    """
    for file_name in os.listdir(samples_path):
        file_path = os.path.join(samples_path, file_name)
        file_mode = os.lstat(file_path).st_mode
        if stat.S_ISDIR(file_mode):
            # For a directory, recurse further down the tree.
            scan_directory(samples_path=file_path, sequences_path=sequences_path, commit=commit)
        elif stat.S_ISLNK(file_mode):
            # For a link, evaluate the link.
            source_path = os.readlink(file_path)
            if not os.path.exists(source_path):
                print(f'source {source_path!r} target {file_name!r} of old symbolic link')
                # Split the entire source path to get the last two components of the source path including
                # the lanes-specific BAM file and the flow cell directory.
                source_path_list = list()
                while 1:
                    source_path, component = os.path.split(source_path)
                    if component:
                        source_path_list.append(component)
                    else:
                        if source_path:
                            source_path_list.append(source_path)
                        break
                source_path_list.reverse()
                source_path = os.path.relpath(
                    os.path.join(sequences_path, *source_path_list[-2:]),
                    samples_path)
                print(f'source {source_path!r} target {file_name!r} of new symbolic link')

                # Check that the directory is writable.
                directory_mode = os.lstat(samples_path)

                if directory_mode[stat.ST_MODE] & stat.S_IWUSR:
                    # The old symbolic link need deleting before a new one can be set.
                    # Secure against race conditions in case the link has already been deleted
                    # before calling os.remove() or added before calling os.symlink().
                    if commit:
                        try:
                            os.remove(file_path)
                        except OSError as exception:
                            if exception.errno != errno.ENOENT:
                                raise exception

                        try:
                            os.symlink(source_path, file_path)
                        except OSError as exception:
                            if exception.errno != errno.EEXIST:
                                raise exception
                else:
                    print('No write permission for link in directory:', repr(samples_path))

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Fix symbolic links between the BSF samples and BSF sequences directory.')

    argument_parser.add_argument(
        '--dry-run',
        action='store_false',
        help='dry run',
        dest='commit')

    argument_parser.add_argument(
        '--samples-path',
        help='samples archive directory path')

    argument_parser.add_argument(
        '--sequences-path',
        help='sequence archive directory path')

    name_space = argument_parser.parse_args()

    return scan_directory(
        samples_path=name_space.samples_path,
        sequences_path=name_space.sequences_path,
        commit=name_space.commit)


if __name__ == '__main__':
    sys.exit(main())
