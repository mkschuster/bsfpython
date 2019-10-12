#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""Sample Annotation Sheet Runnables module.

A package of classes and methods to prune a sample annotation sheet.

Reads objects not associated with (FASTQ or BAM) files get deleted from their corresponding PairedReads objects.
A Reads object is retained, if a file with a size equal to or larger than the configured minimum file size or
an (empty) status file (*.status) exists in the file system. Files smaller then the configured minimum file size
are automatically deleted.
"""
#  Copyright 2013 - 2019 Michael K. Schuster
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
import argparse
import os
import sys

import bsf.argument
import bsf.ngs
import bsf.procedure


def _convert_collection(file_path_old, file_path_new, minimum_size, drop_read_1=False, drop_read_2=False, debug=0):
    """Private function to convert a collection into a Sample Annotation Sheet.

    @param file_path_old: Old sample annotation sheet file path
    @type file_path_old: str
    @param file_path_new: New sample annotation sheet file path
    @type file_path_new: str
    @param minimum_size: Minimum file size of the FASTQ file to keep it
    @type minimum_size: int
    @param drop_read_1: Drop read 1
    @type drop_read_1: bool
    @param drop_read_2: Drop read 2
    @type drop_read_2: bool
    @param debug: Debug level
    @type debug: int
    @return:
    @rtype:
    """
    collection = bsf.ngs.Collection.from_sas_path(
        file_path='',
        file_type='',
        name='picard_sam_to_fastq',
        sas_path=file_path_old,
        debug=debug)

    if debug > 0:
        print('Initial Collection:')
        sys.stdout.writelines(collection.trace(level=1))

    for prf in collection.processed_run_folder_dict.values():
        for project in prf.project_dict.values():
            for sample in project.sample_dict.values():
                new_paired_reads_list = list()
                """ @type new_paired_reads_list: list[bsf.ngs.PairedReads] """
                for paired_reads in sample.paired_reads_list:
                    paired_reads_keep = False
                    if paired_reads.reads_1 is not None:
                        if os.path.exists(paired_reads.reads_1.file_path):
                            if os.path.getsize(paired_reads.reads_1.file_path) >= minimum_size:
                                paired_reads_keep = True
                            else:
                                os.remove(paired_reads.reads_1.file_path)
                                paired_reads.reads_1 = None
                        elif os.path.exists(paired_reads.reads_1.file_path + '.truncated'):
                            paired_reads_keep = True
                        else:
                            # The PairedReads object does not have a meaningful Reads object in reads1.
                            paired_reads.reads_1 = None
                    if drop_read_1:
                        paired_reads.reads_1 = None
                    if paired_reads.reads_2 is not None:
                        if os.path.exists(paired_reads.reads_2.file_path):
                            if os.path.getsize(paired_reads.reads_2.file_path) >= minimum_size:
                                paired_reads_keep = True
                            else:
                                os.remove(paired_reads.reads_2.file_path)
                                paired_reads.reads_2 = None
                        elif os.path.exists(paired_reads.reads_2.file_path + '.truncated'):
                            paired_reads_keep = True
                        else:
                            # The PairedReads object does not have a meaningful Reads object in reads2.
                            paired_reads.reads_2 = None
                    if drop_read_2:
                        paired_reads.reads_2 = None
                    if paired_reads.reads_1 is None and paired_reads.reads_2 is None:
                        paired_reads_keep = False
                    if paired_reads_keep:
                        new_paired_reads_list.append(paired_reads)
                sample.paired_reads_list = new_paired_reads_list
                # The Sample object could have lost all its PairedReads objects.
                # The PairedReads objects may no longer have the correct weak reference to their Sample.

    collection.to_sas_path(name='picard_sam_to_fastq', file_path=file_path_new)

    if debug > 0:
        print('Final Collection:')
        sys.stdout.writelines(collection.trace(level=1))

    return


def run(runnable):
    """Run the the C{bsf.procedure.ConsecutiveRunnable}.

    @param runnable: C{bsf.procedure.ConsecutiveRunnable}
    @type runnable: bsf.procedure.ConsecutiveRunnable
    @return:
    @rtype:
    """

    def run_get_value(key):
        """Get the value of the first C{bsf.argument.OptionLong} object registered under a key
        in the first C{bsf.process.RunnableStep} object of this C{bsf.procedure.ConsecutiveRunnable} object.

        @param key: C{bsf.argument.OptionLong} key
        @type key: str
        @return: C{bsf.argument.OptionLong} value
        @rtype: str | None
        """
        if key in runnable_step.options and runnable_step.options[key]:
            return runnable_step.options[key][0].value
        else:
            return None

    # If the ConsecutiveRunnable status file exists, there is nothing to do and
    # this ConsecutiveRunnable should not have been submitted in the first place.

    if os.path.exists(runnable.runnable_status_file_path(success=True, absolute=False)):
        return

    # Do the work.

    runnable_step = runnable.runnable_step_list[0]

    _convert_collection(
        file_path_old=run_get_value(key='sas_path_old'),
        file_path_new=run_get_value(key='sas_path_new'),
        minimum_size=int(run_get_value(key='minimum_size')),
        drop_read_1=bool(run_get_value(key='drop_read_1')),
        drop_read_2=bool(run_get_value(key='drop_read_2')),
        debug=runnable.debug)

    runnable_step.remove_obsolete_file_paths()

    # Upon success, create a ConsecutiveRunnable-specific status file that indicates completion
    # for the whole ConsecutiveRunnable.

    runnable.runnable_status_file_remove()
    runnable.runnable_status_file_create(success=True)

    return


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--debug',
        help='Debug level',
        required=False,
        type=int)

    argument_parser.add_argument(
        '--old-sas-path',
        dest='old_sas_path',
        help='Old sample annotation sheet file path',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--new-sas-path',
        dest='new_sas_path',
        help='New sample annotation sheet file path',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--minimum-size',
        default=1024,
        dest='minimum_size',
        help='Minimum file size',
        required=False,
        type=int)

    argument_parser.add_argument(
        '--drop-read-1',
        action='store_true',
        dest='drop_read_1',
        help='Drop read 1',
        required=False)

    argument_parser.add_argument(
        '--drop-read-2',
        action='store_true',
        dest='drop_read_2',
        help='Drop read 2',
        required=False)

    name_space = argument_parser.parse_args()

    _convert_collection(
        file_path_old=name_space.old_sas_path,
        file_path_new=name_space.new_sas_path,
        minimum_size=name_space.minimum_size,
        drop_read_1=name_space.drop_read_1,
        drop_read_2=name_space.drop_read_2,
        debug=name_space.debug)
