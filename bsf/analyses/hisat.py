# -*- coding: utf-8 -*-
"""HISAT Analysis module.

A package of classes and methods supporting the HISAT graph-based alignment of next generation sequencing reads
to a population of genomes.

Project:  https://ccb.jhu.edu/software/hisat2/index.shtml
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

from bsf.analyses.aligner import Aligner, FilePathAlign as AlignerFilePathAlign
from bsf.connector import ConnectorFile
from bsf.ngs import Collection, Sample
from bsf.process import RunnableStep
from bsf.standards import Configuration, StandardFilePath


class FilePathAlign(AlignerFilePathAlign):
    """The C{bsf.analyses.hisat.FilePathAlign} class models file paths at the alignment stage.

    Attributes:
    @ivar summary_txt: Alignment summary file
    @type summary_txt: str
    """

    def __init__(self, prefix):
        """Initialise a C{FilePathAlign} object.

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        self.summary_txt = os.path.join(prefix, '_'.join((prefix, 'summary.txt')))

        return


class Hisat2(Aligner):
    """The C{bsf.analyses.hisat.Hisat2} class represents the logic to run a (short read) aligner.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar rna_strand: mRNA strand (i.e. F, R, FR or RF)
    @type rna_strand: str
    """

    name = 'HISAT2 Analysis'
    prefix = 'hisat2'

    @classmethod
    def get_file_path_align(cls, paired_reads_name):
        """Get a C{FilePathAlign} object from this or a sub-class.

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: C{FilePathAlign} or sub-class object
        @rtype: FilePathAlign
        """
        return FilePathAlign(prefix=cls.get_prefix_align(paired_reads_name=paired_reads_name))

    def __init__(
            self,
            configuration=None,
            project_name=None,
            genome_version=None,
            input_directory=None,
            output_directory=None,
            project_directory=None,
            genome_directory=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            genome_fasta=None,
            genome_index=None,
            skip_mark_duplicates=None,
            java_archive_picard=None,
            rna_strand=None):
        """Initialise a C{bsf.analyses.hisat.Hisat2}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str | None
        @param genome_version: Genome version
        @type genome_version: str | None
        @param input_directory: C{bsf.analysis.Analysis}-wide input directory
        @type input_directory: str | None
        @param output_directory: C{bsf.analysis.Analysis}-wide output directory
        @type output_directory: str | None
        @param project_directory: C{bsf.analysis.Analysis}-wide project directory,
            normally under the C{bsf.analysis.Analysis}-wide output directory
        @type project_directory: str | None
        @param genome_directory: C{bsf.analysis.Analysis}-wide genome directory,
            normally under the C{bsf.analysis.Analysis}-wide project directory
        @type genome_directory: str | None
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str | None
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.analysis.Stage} objects
        @type stage_list: list[bsf.analysis.Stage] | None
        @param collection: C{bsf.ngs.Collection}
        @type collection: Collection | None
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[Sample] | None
        @param genome_fasta: Genome FASTA file
        @type genome_fasta: str | None
        @param genome_index: Genome index
        @type genome_index: str | None
        @param skip_mark_duplicates: Mark duplicates
        @type skip_mark_duplicates: bool | None
        @param java_archive_picard: Picard tools Java Archive (JAR) file path
        @type java_archive_picard: str | None
        @param rna_strand: mRNA strand (i.e. F, R, FR or RF)
        @type rna_strand: str
        """
        super(Hisat2, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            e_mail=e_mail,
            debug=debug,
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list,
            genome_fasta=genome_fasta,
            genome_index=genome_index,
            skip_mark_duplicates=skip_mark_duplicates,
            java_archive_picard=java_archive_picard)

        # Sub-class specific ...

        self.rna_strand = rna_strand

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.hisat.Hisat2} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """
        super(Hisat2, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'rna_strand'
        if configuration.config_parser.has_option(section=section, option=option):
            self.rna_strand = configuration.config_parser.get(section=section, option=option)

        return

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add one or more HISAT2-specific C{bsf.process.RunnableStep} objects to the
        C{bsf.procedure.ConcurrentRunnable}.

        @param runnable_align: C{bsf.procedure.ConcurrentRunnable}
        @type runnable_align: bsf.procedure.ConcurrentRunnable
        @param stage_align: C{bsf.analysis.Stage}
        @type stage_align: bsf.analysis.Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | None
        """
        file_path_align = FilePathAlign(prefix=runnable_align.name)

        runnable_step = RunnableStep(
            name='HISAT2',
            program='hisat2',
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'))
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        self.set_runnable_step_configuration(runnable_step=runnable_step)

        runnable_step.add_option_short(key='S', value=file_path_align.aligned_sam)
        runnable_step.add_option_short(key='x', value=self.genome_index)
        runnable_step.add_option_long(key='threads', value=str(stage_align.threads))

        if file_path_1 and not file_path_2:
            runnable_step.add_option_short(key='U', value=file_path_1)
        else:
            runnable_step.add_option_short(key='1', value=file_path_1)
            runnable_step.add_option_short(key='2', value=file_path_2)

        # TODO: The --rna-strandness option would require attributes from the sample.
        # This function interface thus needs to change to accommodate for both Sample and PairedReads objects.
        if self.rna_strand:
            runnable_step.add_option_long(key='rna-strandness', value=self.rna_strand)

        # For cufflinks compatible alignment reporting (XS:A:[+|-]).
        runnable_step.add_switch_long(key='dta-cufflinks')

        # Paired-end options, if --no-spliced-alignment is set.
        # --minins <int>
        # --maxins <int>

        # --fr/--rf/--ff

        # --summary-file
        runnable_step.add_option_long(key='summary-file', value=file_path_align.summary_txt)
        # --new-summary
        runnable_step.add_switch_long(key='new-summary')

        # --rg-id <text>
        # --rg <text>

        # Enable memory-mapped (i.e. shared) indices.
        runnable_step.add_switch_long(key='mm')

        return

    def run(self):
        """Run a C{bsf.analyses.hisat.Hisat2} analysis.
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
            self.sas_file = self.get_annotation_file(prefix_list=[Hisat2.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation file in the current working directory.')

        # The Hisat2 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                StandardFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='hisat2'),
                self.genome_version)

        if self.rna_strand and self.rna_strand not in ('F', 'R', 'FR', 'RF'):
            raise Exception("The 'rna_strand' configuration option has to be 'F', 'R', 'FR', 'RF' or ''.")

        super(Hisat2, self).run()

        return
