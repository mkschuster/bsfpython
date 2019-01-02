# -*- coding: utf-8 -*-
"""Bowtie Analysis module

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

from __future__ import print_function

import os

import bsf
import bsf.analyses.aligner
import bsf.process


class Bowtie1(bsf.analyses.aligner.Aligner):
    """The C{bsf.analyses.bowtie.Bowtie1} class represents the logic to run the Bowtie1 aligner.

    Attributes:
    """
    name = 'Bowtie1 Analysis'
    prefix = 'bowtie1'

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a Bowtie1-specific C{bsf.process.RunnableStep} to the C{bsf.Runnable}.

        @param runnable_align: C{bsf.Runnable}
        @type runnable_align: bsf.Runnable
        @param stage_align: C{bsf.Stage}
        @type stage_align: bsf.Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | unicode | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | unicode | None
        @return:
        @rtype:
        """
        file_path_align = runnable_align.file_path_object
        """ @type file_path_align bsf.analyses.aligner.FilePathAlign """

        runnable_step = runnable_align.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='bowtie1',
                program='bowtie'))
        """ @type runnable_step: bsf.process.RunnableStep """
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

        @return:
        @rtype:
        """
        # The Bowtie1 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                bsf.standards.FilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        super(Bowtie1, self).run()

        return


class Bowtie2(bsf.analyses.aligner.Aligner):
    """The C{bsf.analyses.bowtie.Bowtie2} class represents the logic to run the Bowtie2 aligner.

    Attributes:
    """
    name = 'Bowtie2 Analysis'
    prefix = 'bowtie2'

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a Bowtie2-specific C{bsf.process.RunnableStep} to the C{bsf.Runnable}.

        @param runnable_align: C{bsf.Runnable}
        @type runnable_align: bsf.Runnable
        @param stage_align: C{bsf.Stage}
        @type stage_align: bsf.Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | unicode | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | unicode | None
        @return:
        @rtype:
        """
        file_path_align = runnable_align.file_path_object
        """ @type file_path_align bsf.analyses.aligner.FilePathAlign """

        runnable_step = runnable_align.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='bowtie2',
                program='bowtie2'))
        """ @type runnable_step: bsf.process.RunnableStep """
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

        @return:
        @rtype:
        """
        # The Bowtie2 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                bsf.standards.FilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        super(Bowtie2, self).run()

        return
