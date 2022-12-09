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
"""The :py:mod:`bin.bsf_utility_count_samples_per_lane` module is a script to report the number of samples per lane.
"""

import os
import sys
from argparse import ArgumentParser
from typing import Optional

from bsf.analyses.illumina_to_bam_tools import LibraryAnnotationSheet


def run(
        directory_path: str,
        ascending_order: Optional[bool] = None) -> int:
    """Run function.

    :param directory_path: A UNIX configuration (*.ini) file path.
    :type directory_path: str
    :param ascending_order: Request submitting into the DRMS.
    :type ascending_order: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    print('Lane,Sample Number,Comment')

    file_name_list = os.listdir(directory_path)
    file_name_list.sort()

    if not ascending_order:
        file_name_list.reverse()

    for file_name in file_name_list:
        if file_name[-14:] == '_libraries.csv':
            sas: LibraryAnnotationSheet = LibraryAnnotationSheet.from_file_path(
                file_path=os.path.join(directory_path, file_name))

            lane_dict = dict()
            for row_dict in sas.row_dict_list:
                if row_dict['lane'] in lane_dict:
                    lane_dict[row_dict['lane']] += 1
                else:
                    lane_dict[row_dict['lane']] = 1

            for i in range(1, 9):
                key = str(i)
                lane_name = '_'.join((file_name[:-15], key))
                if key in lane_dict:
                    print(lane_name + ',' + str(lane_dict[key]))
                else:
                    print(lane_name + ',1,not annotated')

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Count samples per lane based on BamIndexDecoder library annotation files.')

    argument_parser.add_argument(
        '--directory-path',
        help='directory of BamIndexDecoder library annotation files')

    argument_parser.add_argument(
        '--ascending-order',
        action='store_true',
        help='sort flow cells in ascending order rather than in descending by default')

    name_space = argument_parser.parse_args()

    return run(directory_path=name_space.directory_path, ascending_order=name_space.ascending_order)


if __name__ == '__main__':
    sys.exit(main())
