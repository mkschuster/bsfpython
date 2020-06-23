# -*- coding: utf-8 -*-
"""Bowtie Analysis module.

A package of classes and methods supporting Bowtie alignment analyses.
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
import os

from bsf.analyses.aligner import Aligner, FilePathAlign
from bsf.analysis import Stage
from bsf.connector import ConnectorFile
from bsf.procedure import ConcurrentRunnable
from bsf.process import RunnableStep
from bsf.standards import FilePath as StandardsFilePath


class Bowtie1(Aligner):
    """The C{bsf.analyses.bowtie.Bowtie1} class represents the logic to run the Bowtie1 aligner.

    Attributes:
    """
    name = 'Bowtie1 Analysis'
    prefix = 'bowtie1'

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a Bowtie1-specific C{bsf.process.RunnableStep} to the C{bsf.procedure.ConcurrentRunnable}.

        @param runnable_align: C{bsf.procedure.ConcurrentRunnable}
        @type runnable_align: ConcurrentRunnable
        @param stage_align: C{bsf.analysis.Stage}
        @type stage_align: Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | None
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
        """Run a C{bsf.analyses.bowtie.Bowtie1} analysis.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Get the sample annotation sheet.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[Bowtie1.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # The Bowtie1 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                StandardsFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        super(Bowtie1, self).run()

        return


class Bowtie2(Aligner):
    """The C{bsf.analyses.bowtie.Bowtie2} class represents the logic to run the Bowtie2 aligner.

    Attributes:
    """
    name = 'Bowtie2 Analysis'
    prefix = 'bowtie2'

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a Bowtie2-specific C{bsf.process.RunnableStep} to the C{bsf.procedure.ConcurrentRunnable}.

        @param runnable_align: C{bsf.procedure.ConcurrentRunnable}
        @type runnable_align: ConcurrentRunnable
        @param stage_align: C{bsf.analysis.Stage}
        @type stage_align: Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | None
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
        """Run a C{bsf.analyses.bowtie.Bowtie2} analysis.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Get the sample annotation sheet.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception('Sample annotation file ' + repr(self.sas_file) + ' does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[Bowtie2.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # The Bowtie2 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                StandardsFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        super(Bowtie2, self).run()

        return
