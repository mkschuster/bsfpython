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
"""NGS Collection module.

A package of classes and methods to prune a sample annotation sheet.

All Reads objects not associated with (FASTQ or BAM) files get deleted from their corresponding PairedReads objects.
A Reads object is retained, if a file with a size equal to or larger than the configured minimum file size or
an (empty) status file (*.truncated) exists in the file system. Files smaller than the configured minimum file size
are automatically deleted.
"""
import os
import sys
from argparse import ArgumentParser
from subprocess import Popen

from bsf.connector import Connector
from bsf.ngs import Collection
from bsf.process import Command, Executable, RunnableStep


class RunnableStepCollectionPruneFastq(RunnableStep):
    """The C{bsf.executables.collection.RunnableStepCollectionPruneFastq} class prunes a sample annotation sheet.

    @ivar file_path_old: Old Sample Annotation Sheet file path
    @type file_path_old: str | None
    @ivar file_path_new: New Sample Annotation Sheet file path
    @type file_path_new: str | None
    @ivar minimum_size: Minimum FASTQ file size
    @type minimum_size: int
    @ivar drop_read_1: Drop read 1
    @type drop_read_1: bool
    @ivar drop_read_2: Drop read 2
    @type drop_read_2: bool
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            file_path_old=None,
            file_path_new=None,
            minimum_size=1024,
            drop_read_1=False,
            drop_read_2=False):
        """Initialise a C{bsf.executables.collection.RunnableStepCollectionPruneFastq}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param file_path_old: Old Sample Annotation Sheet file path
        @type file_path_old: str | None
        @param file_path_new: New Sample Annotation Sheet file path
        @type file_path_new: str | None
        @param minimum_size: Minimum FASTQ file size
        @type minimum_size: int
        @param drop_read_1: Drop read 1
        @type drop_read_1: bool
        @param drop_read_2: Drop read 2
        @type drop_read_2: bool
        """
        super(RunnableStepCollectionPruneFastq, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.file_path_old = file_path_old
        self.file_path_new = file_path_new
        self.minimum_size = minimum_size
        self.drop_read_1 = drop_read_1
        self.drop_read_2 = drop_read_2

        return

    def run(self, debug=0):
        """Run a C{bsf.executables.bed.RunnableStepRescaleScore}.

        @param debug: Debug level
        @type debug: int
        @return: Python C{list} of Python C{str} (exception) objects
        @rtype: list[str] | None
        """
        collection = Collection.from_sas_path(
            file_path='',
            file_type='',
            name='sample_annotation_sheet',
            sas_path=self.file_path_old,
            sas_prefix=None)

        if debug > 0:
            print('Initial Collection:')
            sys.stdout.writelines(collection.trace(level=1))

        for prf in collection.processed_run_folder_dict.values():
            for project in prf.project_dict.values():
                for sample in project.sample_dict.values():
                    new_paired_reads_list = list()
                    for paired_reads in sample.paired_reads_list:
                        paired_reads_keep = False
                        if paired_reads.reads_1 is not None:
                            if os.path.exists(paired_reads.reads_1.file_path):
                                if os.path.getsize(paired_reads.reads_1.file_path) >= self.minimum_size:
                                    paired_reads_keep = True
                                else:
                                    os.remove(paired_reads.reads_1.file_path)
                                    paired_reads.reads_1 = None
                            elif os.path.exists(paired_reads.reads_1.file_path + '.truncated'):
                                paired_reads_keep = True
                            else:
                                # The PairedReads object does not have a meaningful Reads object in reads1.
                                paired_reads.reads_1 = None
                        if self.drop_read_1:
                            paired_reads.reads_1 = None
                        if paired_reads.reads_2 is not None:
                            if os.path.exists(paired_reads.reads_2.file_path):
                                if os.path.getsize(paired_reads.reads_2.file_path) >= self.minimum_size:
                                    paired_reads_keep = True
                                else:
                                    os.remove(paired_reads.reads_2.file_path)
                                    paired_reads.reads_2 = None
                            elif os.path.exists(paired_reads.reads_2.file_path + '.truncated'):
                                paired_reads_keep = True
                            else:
                                # The PairedReads object does not have a meaningful Reads object in reads2.
                                paired_reads.reads_2 = None
                        if self.drop_read_2:
                            paired_reads.reads_2 = None
                        if paired_reads.reads_1 is None and paired_reads.reads_2 is None:
                            paired_reads_keep = False
                        if paired_reads_keep:
                            new_paired_reads_list.append(paired_reads)
                    sample.paired_reads_list = new_paired_reads_list
                    # The Sample object could have lost all its PairedReads objects.
                    # The PairedReads objects may no longer have the correct weak reference to their Sample.

        collection.to_sas_path(name='picard_sam_to_fastq', file_path=self.file_path_new)

        if debug > 0:
            print('Final Collection:')
            sys.stdout.writelines(collection.trace(level=1))

        return None


if __name__ == '__main__':
    argument_parser = ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--debug',
        default=0,
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

    runnable_step = RunnableStepCollectionPruneFastq(
        name='prune_sample_annotation_sheet',
        file_path_old=name_space.old_sas_path,
        file_path_new=name_space.new_sas_path,
        minimum_size=name_space.minimum_size,
        drop_read_1=name_space.drop_read_1,
        drop_read_2=name_space.drop_read_2)

    runnable_step.run(debug=name_space.debug)
