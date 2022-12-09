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
"""The :py:mod:`bsf.executables.collection` module provides classes to prune a sample annotation sheet represented by
a :py:class:`bsf.ngs.Collection` object.

All :py:class:`bsf.ngs.Reads` objects not associated with (FASTQ or BAM) files get deleted from their corresponding
:py:class:`bsf.ngs.PairedReads` objects.
A :py:class:`bsf.ngs.Reads` object is retained, if a file with a size equal to or larger than the
configured minimum file size or an (empty) status file (:literal:`*.truncated`) exists in the file system.
Files smaller than the configured minimum file size are automatically deleted.
"""
import logging
import os
from argparse import ArgumentParser
from subprocess import Popen

from bsf.connector import Connector
from bsf.ngs import Collection
from bsf.process import Command, RunnableStep

module_logger = logging.getLogger(name=__name__)


class RunnableStepCollectionPruneFastq(RunnableStep):
    """The :py:class:`bsf.executables.collection.RunnableStepCollectionPruneFastq` class prunes a
    sample annotation sheet.

    :ivar file_path_old: An old Sample Annotation Sheet file path.
    :type file_path_old: str | None
    :ivar file_path_new: A new Sample Annotation Sheet file path.
    :type file_path_new: str | None
    :ivar minimum_size: A minimum FASTQ file size.
    :type minimum_size: int
    :ivar drop_read_1: Request dropping of read 1.
    :type drop_read_1: bool
    :ivar drop_read_2: Request dropping of read 2.
    :type drop_read_2: bool
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
        """Initialise a :py:class:`bsf.executables.collection.RunnableStepCollectionPruneFastq` object.

        :param name: A name.
        :type name: str | None
        :param program: A program.
        :type program: str | None
        :param options: A Python :py:class:`dict` object of
            Python :py:class:`str` (:py:attr:`bsf.argument.Argument.key`) key and
            Python :py:class:`list` value objects of :py:class:`bsf.argument.Argument` objects.
        :type options: dict[str, list[Argument]] | None
        :param arguments: A Python :py:class:`list` object of Python :py:class:`str` (program argument) objects.
        :type arguments: list[str] | None
        :param sub_command: A subordinate :py:class:`bsf.process.Command` object.
        :type sub_command: Command | None
        :param stdin: A standard input :literal:`STDIN` :py:class:`bsf.connector.Connector` object.
        :type stdin: Connector | None
        :param stdout: A standard output :literal:`STDOUT` :py:class:`bsf.connector.Connector` object.
        :type stdout: Connector | None
        :param stderr: A standard error :literal:`STDERR` :py:class:`bsf.connector.Connector` object.
        :type stderr: Connector | None
        :param dependencies: A Python :py:class:`list` object of
            Python :py:class:`str` (:py:attr:`bsf.process.Executable.name`) objects
            in the context of :py:class:`bsf.analysis.Stage` dependencies.
        :type dependencies: list[str] | None
        :param hold: Request a hold on job scheduling.
        :type hold: bool | None
        :param submit: Request the submission via the :py:meth:`bsf.analysis.Stage.submit` method.
        :type submit: bool
        :param process_identifier: A process identifier.
        :type process_identifier: str | None
        :param process_name: A process name.
        :type process_name: str | None
        :param sub_process: A :py:class:`subprocess.Popen` object.
        :type sub_process: Popen | None
        :param obsolete_file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (file path) objects
            that can be removed after successfully completing the :py:meth:`bsf.process.RunnableStep.run` method.
        :type obsolete_file_path_list: list[str] | None
        :param file_path_old: An old Sample Annotation Sheet file path.
        :type file_path_old: str | None
        :param file_path_new: A new Sample Annotation Sheet file path.
        :type file_path_new: str | None
        :param minimum_size: A minimum FASTQ file size.
        :type minimum_size: int
        :param drop_read_1: Request dropping of read 1.
        :type drop_read_1: bool
        :param drop_read_2: Request dropping of read 2.
        :type drop_read_2: bool
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

    def run(self):
        """Run a :py:class:`bsf.executables.collection.RunnableStepCollectionPruneFastq` object.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (exception) objects.
        :rtype: list[str] | None
        """
        collection = Collection.from_sas_path(
            file_path='',
            file_type='',
            name='sample_annotation_sheet',
            sas_path=self.file_path_old,
            sas_prefix=None)

        module_logger.debug('Initial Collection: %r', collection)

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

        module_logger.debug('Final Collection: %r', collection)

        return None


if __name__ == '__main__':
    argument_parser = ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--logging-level',
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
        default='INFO',
        dest='logging_level',
        help='Logging level [INFO]',
        required=False)

    argument_parser.add_argument(
        '--old-sas-path',
        dest='old_sas_path',
        help='Old sample annotation sheet file path',
        required=True)

    argument_parser.add_argument(
        '--new-sas-path',
        dest='new_sas_path',
        help='New sample annotation sheet file path',
        required=True)

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

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    runnable_step = RunnableStepCollectionPruneFastq(
        name='prune_sample_annotation_sheet',
        file_path_old=name_space.old_sas_path,
        file_path_new=name_space.new_sas_path,
        minimum_size=name_space.minimum_size,
        drop_read_1=name_space.drop_read_1,
        drop_read_2=name_space.drop_read_2)

    runnable_step.run()
