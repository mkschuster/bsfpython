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
"""The :py:mod:`bsf.executables.bed` module provides classes and functions to rescale
Browser Extensible Data (BED) scores from 0 to 1000.
"""
import logging
from argparse import ArgumentParser
from subprocess import Popen
from typing import Optional

from bsf.connector import Connector
from bsf.process import Command, RunnableStep


class RunnableStepRescaleScore(RunnableStep):
    """The :py:class:`bsf.executables.bed.RunnableStepRescaleScore` class rescales the BED score.

    :ivar file_path_old: An old BED file path.
    :type file_path_old: str | None
    :ivar file_path_new: A new BED file path.
    :type file_path_new: str | None
    :ivar keep_header_lines: Request keeping :literal:`browser` and :literal:`track` lines.
    :type keep_header_lines: bool | None
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
            keep_header_lines=None):
        """Initialise a :py:class:`bsf.executables.bed.RunnableStepRescaleScore` object.

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
        :param file_path_old: An old BED file path.
        :type file_path_old: str | None
        :param file_path_new: A new BED file path.
        :type file_path_new: str | None
        :param keep_header_lines: Request keeping :literal:`browser` and :literal:`track` lines.
        :type keep_header_lines: bool | None
        """
        super(RunnableStepRescaleScore, self).__init__(
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
        self.keep_header_lines = keep_header_lines

        return

    def run(self):
        """Run a :py:class:`bsf.executables.bed.RunnableStepRescaleScore` object.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (exception) objects.
        :rtype: list[str] | None
        """
        # BED variables:
        #
        # chrom
        # chromStart
        # chromEnd
        # name
        # score
        # strand
        # thickStart
        # thickEnd
        # itemRgb
        # blockCount
        # blockSizes
        # blockStarts

        # First pass, get the maximum and minimum scores.
        maximum_score: Optional[float] = None
        minimum_score: Optional[float] = None
        integer_score: bool = True

        with open(file=self.file_path_old, mode='rt') as input_text_io:
            for line in input_text_io:
                if line.startswith('#') or line.startswith('browser') or line.startswith('track'):
                    continue

                line_list = line.split(sep='\t')

                if len(line_list) >= 5:
                    # Skip scores that are undefined (i.e., '.').
                    if line_list[4] == '.':
                        continue

                    # Test for integer scores.
                    score = float(line_list[4])

                    if integer_score and score != round(score):
                        integer_score = False

                    if maximum_score is None:
                        maximum_score = score
                    else:
                        maximum_score = max(maximum_score, score)

                    if minimum_score is None:
                        minimum_score = score
                    else:
                        minimum_score = min(minimum_score, score)

        upper_score = min(1000.0, maximum_score)
        lower_score = max(0.0, minimum_score)

        scale_factor = (maximum_score - minimum_score) / (upper_score - lower_score) * 1000.0

        # Second pass, rewrite the scaled BED file.
        with open(file=self.file_path_new, mode='wt') as output_text_io:
            with open(file=self.file_path_old, mode='rt') as input_text_io:
                for line in input_text_io:
                    if line.startswith('#') or line.startswith('browser:') or line.startswith('track:'):
                        if self.keep_header_lines:
                            print(line, file=output_text_io)
                        continue
                    line_list = line.split(sep='\t')
                    if len(line_list) >= 5 and line_list[4] != '.':
                        score = float(line_list[4])
                        line_list[4] = str(int(round(score * scale_factor + score - 0.0)))
                    print(line_list, sep='\t', file=output_text_io)

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
        '--keep-header-lines',
        action='store_true',
        default=False,
        dest='keep_header_lines',
        help='Keep header (browser and track) lines',
        required=False)

    argument_parser.add_argument(
        '--old-bed-path',
        dest='bed_vcf_path',
        help='Old (input) BED file path',
        required=True)

    argument_parser.add_argument(
        '--new-bed-path',
        dest='new_bed_path',
        help='New (output) BED file path',
        required=True)

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    runnable_step = RunnableStepRescaleScore(
        name='rescale_bed_score',
        file_path_old=name_space.old_bed_path,
        file_path_new=name_space.new_bed_path,
        keep_header_lines=name_space.keep_header_lines)

    runnable_step.run()
