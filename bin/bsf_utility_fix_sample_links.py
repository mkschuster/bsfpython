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
#
#  BSF Python utility script to find and fix broken symbolic links pointing from the BSF samples to the
#  BSF sequences directory. These symbolic links are set by the bsf.analyses.illumina_to_bam_tools.BamIndexDecoder
#  analysis in case a flow-cell lane does not need de-multiplexing.
#
import errno
import os
import stat
from argparse import ArgumentParser


def scan_directory(directory_path, commit=None):
    """Scan a directory for symbolic links or directories recursively.

    :param directory_path: A directory path.
    :type directory_path: str
    :param commit: Commit changes to the file system.
    :type commit: bool | None
    """
    for file_name in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file_name)
        file_mode = os.lstat(file_path).st_mode
        if stat.S_ISDIR(file_mode):
            # For a directory, recurse further down the tree.
            scan_directory(directory_path=file_path)
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
                    os.path.join(name_space.sequences_path, *source_path_list[-2:]),
                    directory_path)
                print(f'source {source_path!r} target {file_name!r} of new symbolic link')

                # Check that the directory is writable.
                directory_mode = os.lstat(directory_path)

                if directory_mode[stat.ST_MODE] & stat.S_IWUSR:
                    # The old symbolic link need deleting before a new one can be set.
                    # Secure against race conditions in case the link has already been deleted
                    # before calling os.remove() or added before calling os.symlink().
                    if commit:
                        try:
                            os.remove(file_path)
                        except OSError as exception:
                            if exception.errno != errno.ENOENT:
                                raise

                        try:
                            os.symlink(source_path, file_path)
                        except OSError as exception:
                            if exception.errno != errno.EEXIST:
                                raise
                else:
                    print('No write permission for link in directory:', repr(directory_path))

    return


argument_parser = ArgumentParser(
    description='Fix symbolic links between the BSF samples and BSF sequences directory.')

argument_parser.add_argument(
    '--dry-run',
    action='store_false',
    default=True,
    dest='commit',
    help='dry run',
    required=False)

argument_parser.add_argument(
    '--samples-path',
    dest='samples_path',
    help='BSF samples directory path',
    required=False)

argument_parser.add_argument(
    '--sequences-path',
    dest='sequences_path',
    help='BSF sequences directory path',
    required=False)

name_space = argument_parser.parse_args()

scan_directory(directory_path=name_space.samples_path, commit=name_space.commit)
