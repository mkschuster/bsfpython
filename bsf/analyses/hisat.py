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
"""The :py:mod:`bsf.analyses.hisat` module provides classes and methods supporting the
`HISAT <https://daehwankimlab.github.io/hisat2/>`_ graph-based alignment of
next generation sequencing reads to a population of genomes.
"""
import os

from bsf.analyses.aligner import Aligner, FilePathAlign as AlignerFilePathAlign
from bsf.analysis import Stage
from bsf.connector import ConnectorFile
from bsf.ngs import Collection, Sample
from bsf.procedure import ConcurrentRunnable
from bsf.process import RunnableStep
from bsf.standards import Configuration, StandardFilePath, Transcriptome


class FilePathAlign(AlignerFilePathAlign):
    """The :py:class:`bsf.analyses.hisat.FilePathAlign` class models file paths at the alignment stage.

    :ivar summary_txt: An alignment summary file path.
    :type summary_txt: str
    """

    def __init__(self, prefix):
        """Initialise a :py:class:`bsf.analyses.hisat.FilePathAlign` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        self.summary_txt = os.path.join(prefix, '_'.join((prefix, 'summary.txt')))

        return


class Hisat2(Aligner):
    """The :py:class:`bsf.analyses.hisat.Hisat2` class represents the logic to run a (short read) aligner.

    :ivar transcriptome_version: A transcriptome version.
    :type transcriptome_version: str | None
    :ivar transcriptome_gtf: A transcriptome annotation GTF file path.
    :type transcriptome_gtf: str | None
    :ivar transcriptome_index: A transcriptome index directory path.
    :type transcriptome_index: str
    :ivar rna_strand: An mRNA strand (i.e., F, R, FR or RF).
    :type rna_strand: str
    :ivar threads_number: A number of threads.
    :type threads_number: int | None
    """

    name = 'HISAT2 Analysis'
    prefix = 'hisat2'

    sam_attributes_to_retain_list = [
        # http://daehwankimlab.github.io/hisat2/manual/

        # ZS:i:<N> : Alignment score for the best-scoring alignment found other than the alignment reported.
        'ZS',

        # YS:i:<N> : Alignment score for opposite mate in the paired-end alignment.
        'YS',

        # XN:i:<N> : The number of ambiguous bases in the reference covering this alignment.
        'XN',

        # XM:i:<N> : The number of mismatches in the alignment.
        'XM',

        # XO:i:<N> : The number of gap opens, for both read and reference gaps, in the alignment.
        'XO',

        # XG:i:<N> : The number of gap extensions, for both read and reference gaps, in the alignment.
        'XG',

        # YF:Z:<S> : String indicating reason why the read was filtered out.
        'YF',

        # YT:Z:<S> : Value of UU indicates the read was not part of a pair.
        'YT',

        # XS:A:<A> : Values of + and - indicate the read is mapped to transcripts on sense and anti-sense strands,
        # respectively.
        'XS',

        # Zs:Z:<S> : When the alignment of a read involves SNPs that are in the index,
        # this option is used to indicate where exactly the read involves the SNPs.
        'Zs',
    ]

    @classmethod
    def get_file_path_align(cls, paired_reads_name):
        """Get a :py:class:`bsf.analyses.hisat.FilePathAlign` object from this or a subclass.

        :param paired_reads_name: A :py:attr:`bsf.ngs.PairedReads.name` attribute.
        :type paired_reads_name: str
        :return: A :py:class:`bsf.analyses.hisat.FilePathAlign` object or subclass thereof.
        :rtype: FilePathAlign
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
            report_style_path=None,
            report_header_path=None,
            report_footer_path=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            genome_fasta=None,
            genome_index=None,
            transcriptome_version=None,
            transcriptome_gtf=None,
            transcriptome_index=None,
            skip_mark_duplicates=None,
            java_archive_picard=None,
            rna_strand=None,
            threads_number=None):
        """Initialise a :py:class:`bsf.analyses.hisat.Hisat2` object.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration | None
        :param project_name: A project name.
        :type project_name: str | None
        :param genome_version: A genome assembly version.
        :type genome_version: str | None
        :param input_directory: An input directory path.
        :type input_directory: str | None
        :param output_directory: An output directory path.
        :type output_directory: str | None
        :param project_directory: A project directory path, normally under the output directory path.
        :type project_directory: str | None
        :param genome_directory: A genome directory path, normally under the project directory path.
        :type genome_directory: str | None
        :param report_style_path: Report :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: Report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: Report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a UCSC Genome Browser Track Hub.
        :type e_mail: str | None
        :param debug: An integer debugging level.
        :type debug: int | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param genome_fasta: A genome FASTA file path.
        :type genome_fasta: str | None
        :param genome_index: A genome index file path.
        :type genome_index: str | None
        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str | None
        :param transcriptome_gtf: A transcriptome annotation GTF file path.
        :type transcriptome_gtf: str | None
        :param transcriptome_index: A transcriptome index directory path.
        :type transcriptome_index: str
        :param skip_mark_duplicates: Request skipping the Picard :literal:`MarkDuplicates` step.
        :type skip_mark_duplicates: bool | None
        :param java_archive_picard: A Picard tools Java Archive (JAR) file path.
        :type java_archive_picard: str | None
        :param rna_strand: An mRNA strand (i.e., F, R, FR or RF).
        :type rna_strand: str
        :param threads_number: A number of threads.
        :type threads_number: int | None
        """
        super(Hisat2, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            report_style_path=report_style_path,
            report_header_path=report_header_path,
            report_footer_path=report_footer_path,
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

        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf = transcriptome_gtf
        self.transcriptome_index = transcriptome_index

        self.rna_strand = rna_strand
        self.threads_number = threads_number

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a :py:class:`bsf.analyses.hisat.Hisat2` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(Hisat2, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index = configuration.config_parser.get(section=section, option=option)

        option = 'rna_strand'
        if configuration.config_parser.has_option(section=section, option=option):
            self.rna_strand = configuration.config_parser.get(section=section, option=option)

        option = 'threads_number'
        if configuration.config_parser.has_option(section=section, option=option):
            self.threads_number = configuration.config_parser.getint(section=section, option=option)

        return

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add one or more HISAT2-specific :py:class:`bsf.process.RunnableStep` objects to the
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
            name='HISAT2',
            program='hisat2',
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'))
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        self.set_runnable_step_configuration(runnable_step=runnable_step)

        runnable_step.add_option_short(key='S', value=file_path_align.aligned_sam)

        if self.transcriptome_index:
            runnable_step.add_option_short(key='x', value=self.transcriptome_index)
        else:
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

        # Enable memory-mapped (i.e., shared) indices.
        runnable_step.add_switch_long(key='mm')

        return

    def run(self):
        """Run a :py:class:`bsf.analyses.hisat.Hisat2` object.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a valid 'transcriptome_version' configuration option.")

        # The Hisat2 genome index is quite peculiar as it has to be the index file path without file extensions.

        if self.transcriptome_version:
            if not self.transcriptome_index:
                self.transcriptome_index = os.path.join(
                    StandardFilePath.get_resource_transcriptome_index(
                        transcriptome_version=self.transcriptome_version,
                        transcriptome_index='hisat2'),
                    self.transcriptome_version)
        else:
            if not self.genome_index:
                self.genome_index = os.path.join(
                    StandardFilePath.get_resource_genome_index(
                        genome_version=self.genome_version,
                        genome_index='hisat2'),
                    self.genome_version)

        if self.rna_strand and self.rna_strand not in ('F', 'R', 'FR', 'RF'):
            raise Exception("The 'rna_strand' configuration option has to be 'F', 'R', 'FR', 'RF' or ''.")

        if not self.threads_number:
            self.threads_number = 1

        super(Hisat2, self).run()

        return
