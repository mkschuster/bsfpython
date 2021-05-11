# -*- coding: utf-8 -*-
"""BWA Analysis module.

A package of classes and methods supporting Burrows-Wheeler Aligner (BWA) analyses.
"""
#  Copyright 2013 - 2021 Michael K. Schuster
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
from bsf.process import Command, RunnableStep
from bsf.standards import StandardFilePath


class MaximalExactMatches(Aligner):
    """The C{bsf.analyses.bwa.MaximalExactMatches} class represents the BWA Maximal Exact Matches (MEM) algorithm.

    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar sam_attributes_to_retain_list: A Python C{list} of aligner-specific, private SAM tags (i.e. X*, Y*, z*)
        that should be retained by Picard MergeBamAlignment
    @type sam_attributes_to_retain_list: list[str]
    """
    name = 'BWA Maximal Exact Matches Analysis'
    prefix = 'bwa_mem'

    sam_attributes_to_retain_list = [
        # Number of best hits
        'X0',
        # Number of suboptimal hits found by BWA
        'X1',
        # Number of ambiguous bases in the reference
        'XN',
        # Number of mismatches in the alignment
        'XM',
        # Number of gap opens
        'XO',
        # Number of gap extensions
        'XG',
        # Type: Unique/Repeat/N/Mate-sw
        'XT',
        # Alternative hits; format: (chr,pos,CIGAR,NM;)*
        'XA',
        # Suboptimal alignment score
        'XS',
        # Support from forward / reverse alignment
        'XF',
        # Number of supporting seeds
        'XE',
    ]

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add a BWA MEM-specific C{bsf.process.RunnableStep} to the C{bsf.procedure.ConcurrentRunnable}.

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
            name='bwa_mem',
            program='bwa',
            sub_command=Command(name='mem', program='mem'),
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'))
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        sub_command = runnable_step.sub_command
        # -t [1] Number of threads
        sub_command.add_option_short(key='t', value=str(stage_align.threads))
        # Output errors only.
        sub_command.add_option_short(key='v', value='1')
        # Mark shorter split hits as secondary (for Picard compatibility).
        sub_command.add_switch_short(key='M')
        #  Use soft clipping CIGAR operation for supplementary alignments.
        sub_command.add_switch_short(key='Y')
        # Output file [standard output]
        sub_command.add_option_short(key='o', value=file_path_align.aligned_sam)
        # -H header lines
        # -R @RG line
        sub_command.arguments.append(self.genome_index)
        sub_command.arguments.append(file_path_1)
        if file_path_2:
            sub_command.arguments.append(file_path_2)

        return

    def run(self):
        """Run a C{bsf.analyses.bwa.MaximalExactMatches} analysis.
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
            self.sas_file = self.get_annotation_file(prefix_list=[MaximalExactMatches.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # BWA requires the genome.fasta file.
        if not self.genome_index:
            self.genome_index = StandardFilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bwa')

        super(MaximalExactMatches, self).run()

        return
