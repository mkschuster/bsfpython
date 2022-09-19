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
"""The :py:mod:`bsf.analyses.bowtie` module provides classes and methods supporting the Bowtie aligners.
"""
import os

from bsf.analyses.aligner import Aligner, FilePathAlign
from bsf.analysis import Stage
from bsf.connector import ConnectorFile
from bsf.procedure import ConcurrentRunnable
from bsf.process import RunnableStep
from bsf.standards import StandardFilePath


class Bowtie1(Aligner):
    """The :py:class:`bsf.analyses.bowtie.Bowtie1` class represents the logic to run the
    `Bowtie1 <http://bowtie-bio.sourceforge.net/index.shtml>`_ short read aligner.
    """
    name = 'Bowtie1 Analysis'
    prefix = 'bowtie1'

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a Bowtie1-specific :py:class:`bsf.process.RunnableStep` object to the
        :py:class:`bsf.procedure.ConcurrentRunnable` object.

        :param runnable_align: A :py:class:`bsf.procedure.ConcurrentRunnable` object.
        :type runnable_align: ConcurrentRunnable
        :param stage_align: A :py:class:`bsf.analysis.Stage` object.
        :type stage_align: Stage
        :param file_path_1: A :literal:`FASTQ` file path 1.
        :type file_path_1: str | None
        :param file_path_2: A :literal:`FASTQ` file path 2.
        :type file_path_2: str | None
        """
        file_path_align = FilePathAlign(prefix=runnable_align.name)

        runnable_step = RunnableStep(
            name='bowtie1',
            program='bowtie',
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'))
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        runnable_step.arguments.append(self.genome_index)

        # For Bowtie1 unpaired reads are an argument, paired reads come with options -1 <m1> and -2 <m2>.
        if file_path_1 and not file_path_2:
            runnable_step.arguments.append(file_path_1)
        else:
            runnable_step.add_option_short(key='1', value=file_path_1)
            runnable_step.add_option_short(key='2', value=file_path_2)

        runnable_step.arguments.append(file_path_align.aligned_sam)

        return

    def run(self):
        """Run a :py:class:`bsf.analyses.bowtie.Bowtie1` object.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Get the sample annotation sheet.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception(f'Sample annotation sheet {self.sas_file!r} does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[Bowtie1.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation sheet in the current working directory.')

        # The Bowtie1 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                StandardFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        super(Bowtie1, self).run()

        return


class Bowtie2(Aligner):
    """The :py:class:`bsf.analyses.bowtie.Bowtie2` class represents the logic to run the
    `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ short read aligner.
    """
    name = 'Bowtie2 Analysis'
    prefix = 'bowtie2'

    sam_attributes_to_retain_list = [
        # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output

        # XS:i:<N> Alignment score for the best-scoring alignment found other than the alignment reported.
        'XS',
        # YS:i:<N> Alignment score for opposite mate in the paired-end alignment.
        'YS',
        # XN:i:<N> The number of ambiguous bases in the reference covering this alignment.
        'XN',
        # XM:i:<N> The number of mismatches in the alignment.
        'XM',
        # XO:i:<N> The number of gap opens, for both read and reference gaps, in the alignment.
        'XO',
        # XG:i:<N> The number of gap extensions, for both read and reference gaps, in the alignment.
        'XG',
        # YF:Z:<S> String indicating reason why the read was filtered out.
        'YF',
        # YT:Z:<S> Value of UU indicates the read was not part of a pair.
        'YT',
    ]

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a Bowtie2-specific :py:class:`bsf.process.RunnableStep` object to the
        :py:class:`bsf.procedure.ConcurrentRunnable` object.

        :param runnable_align: A :py:class:`bsf.procedure.ConcurrentRunnable` object.
        :type runnable_align: ConcurrentRunnable
        :param stage_align: A :py:class:`bsf.analysis.Stage` object.
        :type stage_align: Stage
        :param file_path_1: A :literal:`FASTQ` file path 1.
        :type file_path_1: str | None
        :param file_path_2: A :literal:`FASTQ` file path 2.
        :type file_path_2: str | None
        """
        file_path_align = FilePathAlign(prefix=runnable_align.name)

        runnable_step = RunnableStep(
            name='bowtie2',
            program='bowtie2',
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'))
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_option_short(key='S', value=file_path_align.aligned_sam)
        runnable_step.add_option_short(key='x', value=self.genome_index)
        runnable_step.add_option_long(key='threads', value=str(stage_align.threads))

        # NOTE: The following options are properties of the Sample,
        # PairedReads and Reads objects.
        # runnable_step.add_switch_short(key='q')
        # runnable_step.add_switch_short(key='phred33')

        if file_path_1 and not file_path_2:
            runnable_step.add_option_short(key='U', value=file_path_1)
        else:
            runnable_step.add_option_short(key='1', value=file_path_1)
            runnable_step.add_option_short(key='2', value=file_path_2)

        return

    def run(self):
        """Run a :py:class:`bsf.analyses.bowtie.Bowtie2` object.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Get the sample annotation sheet.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception(f'Sample annotation sheet {self.sas_file!r} does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[Bowtie2.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation sheet in the current working directory.')

        # The Bowtie2 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                StandardFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        super(Bowtie2, self).run()

        return
