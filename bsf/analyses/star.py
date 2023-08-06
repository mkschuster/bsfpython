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
"""The :py:mod:`bsf.analyses.star` module provides classes and methods supporting the
`Spliced Transcripts Alignment to a Reference (STAR) <https://github.com/alexdobin/STAR>`_ by Alexander Dobin.
"""
import os
from typing import Optional

from bsf.analyses.aligner import Aligner, \
    FilePathAlign as AlignerFilePathAlign, \
    FilePathSample, \
    FilePathSummary as AlignerFilePathSummary
from bsf.analysis import Stage
from bsf.connector import ConnectorFile
from bsf.ngs import Collection, Sample
from bsf.procedure import ConcurrentRunnable, ConsecutiveRunnable
from bsf.process import RunnableStep
from bsf.standards import Configuration, StandardFilePath, Transcriptome


class FilePathAlign(AlignerFilePathAlign):
    """The :py:class:`bsf.analyses.star.FilePathAlign` class models file paths at the alignment stage.

    :ivar aligned_sam: Aligned sequence alignment map (SAM) file path.
    :type aligned_sam: str
    :ivar splice_junctions_tsv: Splice junctions :emphasis:`Tab-Separated Value` (TSV) file path.
    :type splice_junctions_tsv: str
    :ivar star_prefix: STAR :literal:`outFileNamePrefix` file path.
    :type star_prefix: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.star.FilePathAlign` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        # Override the aligned_sam instance variable of the super-class.
        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'Aligned.out.sam')))

        self.splice_junctions_tsv = os.path.join(prefix, '_'.join((prefix, 'SJ.out.tab')))

        self.star_prefix = os.path.join(prefix, '_'.join((prefix, '')))

        return


class FilePathSummary(AlignerFilePathSummary):
    """The :py:class:`bsf.analyses.star.FilePathSummary` class models file paths at the summary stage.
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.star.FilePathSummary` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathSummary, self).__init__(prefix=prefix)

        self.alignment_read_group_pdf = prefix + '_star_alignment_read_group.pdf'
        self.alignment_read_group_png = prefix + '_star_alignment_read_group.png'
        self.alignment_sample_png = prefix + '_star_alignment_sample.png'
        self.alignment_sample_pdf = prefix + '_star_alignment_sample.pdf'
        self.junction_fraction_read_group_pdf = prefix + '_star_junction_fraction_read_group.pdf'
        self.junction_fraction_read_group_png = prefix + '_star_junction_fraction_read_group.png'
        self.junction_fraction_sample_pdf = prefix + '_star_junction_fraction_sample.pdf'
        self.junction_fraction_sample_png = prefix + '_star_junction_fraction_sample.png'
        self.junction_number_read_group_pdf = prefix + '_star_junction_number_read_group.pdf'
        self.junction_number_read_group_png = prefix + '_star_junction_number_read_group.png'
        self.junction_number_sample_pdf = prefix + '_star_junction_number_sample.pdf'
        self.junction_number_sample_png = prefix + '_star_junction_number_sample.png'
        self.mapped_fraction_read_group_pdf = prefix + '_star_mapped_fraction_read_group.pdf'
        self.mapped_fraction_read_group_png = prefix + '_star_mapped_fraction_read_group.png'
        self.mapped_fraction_sample_png = prefix + '_star_mapped_fraction_sample.png'
        self.mapped_fraction_sample_pdf = prefix + '_star_mapped_fraction_sample.pdf'
        self.mapped_number_read_group_pdf = prefix + '_star_mapped_number_read_group.pdf'
        self.mapped_number_read_group_png = prefix + '_star_mapped_number_read_group.png'
        self.mapped_number_sample_png = prefix + '_star_mapped_number_sample.png'
        self.mapped_number_sample_pdf = prefix + '_star_mapped_number_sample.pdf'
        self.table_read_group_tsv = prefix + '_star_table_read_group.tsv'
        self.table_sample_tsv = prefix + '_star_table_sample.tsv'

        return


class Star(Aligner):
    """The :py:class:`bsf.analyses.star.Star` class models the
    `Spliced Transcripts Alignment to a Reference (STAR) <https://github.com/alexdobin/STAR>`_ by Alexander Dobin.

    :ivar transcriptome_version: A transcriptome version.
    :type transcriptome_version: str | None
    :ivar transcriptome_gtf: A transcriptome annotation GTF file path.
    :type transcriptome_gtf: str | None
    :ivar transcriptome_index: A transcriptome index directory path.
    :type transcriptome_index: str | None
    :ivar two_pass_mapping: Basic two-pass mapping
    :type two_pass_mapping: str | None
    """

    name = 'STAR Analysis'
    prefix = 'star'

    sam_attributes_to_retain_list = [
        # The nM SAM tag indicates the number of mismatches (excluding Ns) per read pair.
        'nM',
        # The uT SAM tag indicates the reason for not mapping a read.
        'uT',
    ]

    @classmethod
    def get_file_path_align(cls, paired_reads_name: str) -> FilePathAlign:
        """Get a :py:class:`bsf.analyses.star.FilePathAlign` object from this or a subclass.

        :param paired_reads_name: A :py:attr:`bsf.ngs.PairedReads.name` attribute.
        :type paired_reads_name: str
        :return: A bsf.analyses.star.FilePathAlign` object or subclass thereof.
        :rtype: FilePathAlign
        """
        return FilePathAlign(prefix=cls.get_prefix_align(paired_reads_name=paired_reads_name))

    @classmethod
    def get_file_path_summary(cls) -> FilePathSummary:
        """Get a bsf.analyses.star.FilePathSummary` object from this or a subclass.

        :return: A bsf.analyses.star.FilePathSummary` object or subclass thereof.
        :rtype: FilePathSummary
        """
        return FilePathSummary(prefix=cls.get_prefix_summary())

    def __init__(
            self,
            configuration: Optional[Configuration] = None,
            project_name: Optional[str] = None,
            genome_version: Optional[str] = None,
            input_directory: Optional[str] = None,
            output_directory: Optional[str] = None,
            project_directory: Optional[str] = None,
            genome_directory: Optional[str] = None,
            report_style_path: Optional[str] = None,
            report_header_path: Optional[str] = None,
            report_footer_path: Optional[str] = None,
            e_mail: Optional[str] = None,
            stage_list: Optional[list[Stage]] = None,
            collection: Optional[Collection] = None,
            sample_list: Optional[list[Sample]] = None,
            transcriptome_version: Optional[str] = None,
            transcriptome_gtf: Optional[str] = None,
            transcriptome_index: Optional[str] = None,
            two_pass_mapping: Optional[str] = None,
            skip_mark_duplicates: Optional[bool] = None,
            java_archive_picard: Optional[str] = None):
        """Initialise a :py:class:`bsf.analyses.star.Star` object.

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
        :param report_style_path: A report style :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: A report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: A report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param e_mail: An e-mail address for a :emphasis:`UCSC Genome Browser Track Hub`.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str | None
        :param transcriptome_gtf: A transcriptome annotation GTF file path.
        :type transcriptome_gtf: str | None
        :param transcriptome_index: A transcriptome index directory path.
        :type transcriptome_index: str | None
        :param two_pass_mapping: Basic two-pass mapping
        :type two_pass_mapping: str | None
        :param skip_mark_duplicates: Request skipping the Picard :literal:`MarkDuplicates` step.
        :type skip_mark_duplicates: bool | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        """
        super(Star, self).__init__(
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
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list,
            genome_fasta=None,
            genome_index=None,
            skip_mark_duplicates=skip_mark_duplicates,
            java_archive_picard=java_archive_picard)

        # Sub-class specific ...

        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf = transcriptome_gtf
        self.transcriptome_index = transcriptome_index
        self.two_pass_mapping = two_pass_mapping

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.star.Star` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(Star, self).set_configuration(configuration=configuration, section=section)

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

        option = 'two_pass_mapping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.two_pass_mapping = configuration.config_parser.get(section=section, option=option)

        return

    def add_runnable_step_aligner(
            self,
            runnable_align: ConcurrentRunnable,
            stage_align: Stage,
            file_path_1: Optional[str],
            file_path_2: Optional[str]) -> None:
        """Add one or more STAR-specific :py:class:`bsf.process.RunnableStep` objects to the
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

        # Run STAR

        runnable_step = RunnableStep(
            name='STAR',
            program='STAR',
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'))
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        self.set_runnable_step_configuration(runnable_step=runnable_step)

        runnable_step.add_option_long(key='runThreadN', value=str(stage_align.threads))
        runnable_step.add_option_long(key='genomeDir', value=self.transcriptome_index)
        runnable_step.add_option_long(key='outFileNamePrefix', value=file_path_align.star_prefix)
        runnable_step.add_option_multi_long(key='outSAMunmapped', value='Within KeepPairs')
        if self.two_pass_mapping == 'basic':
            runnable_step.add_option_long(key='twopassMode', value='Basic')
        # NOTE: The STAR command line interface is seriously broken,
        # as the readFilesIn option requires two values.
        # Hence, use class bsf.argument.OptionMultiLong via wrapper Command.add_option_multi_long().
        if file_path_2 is None:
            runnable_step.add_option_long(
                key='readFilesIn',
                value=file_path_1)
        else:
            runnable_step.add_option_multi_long(
                key='readFilesIn',
                value=' '.join((file_path_1, file_path_2)))
        if file_path_1.endswith('fastq.gz'):
            runnable_step.add_option_long(key='readFilesCommand', value='zcat')

        # Run GNU Zip over the rather large splice junction table.

        runnable_step = RunnableStep(
            name='gzip',
            program='gzip')
        runnable_align.add_runnable_step_epilogue(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='best')
        runnable_step.arguments.append(file_path_align.splice_junctions_tsv)

        return

    def add_runnable_step_sample(self, runnable_sample: ConsecutiveRunnable, stage_sample: Stage) -> None:
        """Add one or more STAR-specific :py:class:`bsf.process.RunnableStep` objects
        to the :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable_sample: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :type runnable_sample: ConsecutiveRunnable
        :param stage_sample: A :py:class:`bsf.analysis.Stage` object.
        :type stage_sample: Stage
        """
        # NOTE: This method a copy of Tophat2.add_runnable_step_sample().
        file_path_sample = FilePathSample(prefix=runnable_sample.name)

        # This requires kentUtils to automatically convert via wigToBigWig.
        runnable_step = RunnableStep(
            name='bam2wig',
            program='bam2wig.py',
            obsolete_file_path_list=[file_path_sample.sample_wig])
        runnable_step.add_option_long(key='input-file', value=file_path_sample.sample_bam)
        runnable_step.add_option_long(
            key='chromSize',
            value=StandardFilePath.get_resource_genome_fasta_index(genome_version=self.genome_version))
        runnable_step.add_option_long(key='out-prefix', value=file_path_sample.prefix_prefix)

        runnable_sample.add_runnable_step(runnable_step=runnable_step)

        runnable_step = RunnableStep(
            name='bigwig_info',
            program='bigWigInfo',
            stdout=ConnectorFile(file_path=file_path_sample.sample_bwi, file_mode='wt'))
        runnable_step.arguments.append(file_path_sample.sample_bw)

        runnable_sample.add_runnable_step(runnable_step=runnable_step)

        return

    def add_runnable_step_summary(self, runnable_summary: ConsecutiveRunnable, stage_summary: Stage) -> None:
        """Add one or more STAR-specific :py:class:`bsf.process.RunnableStep` objects
        to the :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable_summary: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :type runnable_summary: ConsecutiveRunnable
        :param stage_summary: A :py:class:`bsf.analysis.Stage` object.
        :type stage_summary: Stage
        """
        file_path_summary = self.get_file_path_summary()

        runnable_step = RunnableStep(
            name='star_summary',
            program='bsf_star_summary.R')
        runnable_summary.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_option_long(key='prefix', value=f'{self.get_stage_name_summary()}_star')
        runnable_step.add_option_long(key='sample-table-path', value=file_path_summary.read_group_to_sample_tsv)

        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.star.Star` object.

        Although STAR can directly count reads according to its splice junction database,
        more than one read group may need aligning so that the count tables had to be combined.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # STAR requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception(f"A {self.name!s} requires a valid 'transcriptome_version' configuration option.")

        if not self.transcriptome_index:
            self.transcriptome_index = StandardFilePath.get_resource_transcriptome_index(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='star')

        if not self.transcriptome_gtf:
            # FIXME: The transcriptome_gtf is currently not used.
            self.transcriptome_gtf = StandardFilePath.get_resource_transcriptome_gtf(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='none',
                basic=True,
                absolute=True)

        if not self.two_pass_mapping:
            self.two_pass_mapping = 'none'

        two_pass_mapping_tuple = ('none', 'basic', 'full')
        if self.two_pass_mapping not in two_pass_mapping_tuple:
            raise Exception(f"The 'two_pass_mapping' option {self.two_pass_mapping!r} "
                            f"is not a member of {two_pass_mapping_tuple!r}.")

        super(Star, self).run()

        return

    def report(self) -> None:
        """Create a :literal:`XHTML 1.0` report and a :emphasis:`UCSC Genome Browser Track Hub`.
        """
        # Create a symbolic link containing the project name and a UUID.
        # This code only needs the public URL.
        link_path = self.create_public_project_link()

        str_list: list[str] = list()

        str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
        str_list.append('\n')

        str_list.extend(self.get_html_genome(genome_version=self.genome_version))
        str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
        str_list.append('\n')

        str_list.append('<h2 id="alignment_visualisation">Alignment Visualisation</h2>\n')
        str_list.append('\n')

        str_list.append('<p id="ucsc_track_hub">')
        str_list.append('Alignments can be visualised by attaching the ')
        str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
        str_list.append('.\n')
        str_list.append('Upon following the link, a project-specific track configuration section ')
        str_list.append('<strong>' + self.project_name + '</strong> ')
        str_list.append('gets added to the UCSC Genome Browser display. By default, all tracks are turned off. ')
        str_list.append('While all tracks can be switched on directly from the configuration section, ')
        str_list.append('especially for larger projects, it may be better to activate individual tracks by ')
        str_list.append('following the track category label, first.\n')
        str_list.append('</p>\n')
        str_list.append('\n')

        str_list.append('<h2 id="qc_plots">QC Plots</h2>\n')
        str_list.append('\n')
        str_list.append('<table id="qc_table">\n')
        str_list.append('<thead>\n')
        str_list.append('<tr>\n')
        str_list.append('<th>Sample</th>\n')
        str_list.append('<th>Read Group</th>\n')
        str_list.append('<th>Metrics</th>\n')
        str_list.append('</tr>\n')
        str_list.append('</thead>\n')
        str_list.append('<tbody>\n')

        file_path_summary = self.get_file_path_summary()

        # Alignment Summary Plots
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.alignment_sample_pdf + '">')
        str_list.append('<img alt="Mapped - Sample"')
        str_list.append(' src="' + file_path_summary.alignment_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.alignment_read_group_pdf + '">')
        str_list.append('<img alt="Mapped - Read Group"')
        str_list.append(' src="' + file_path_summary.alignment_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Mapped</td>\n')
        str_list.append('</tr>\n')

        # Mapped Fraction Plots
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.mapped_fraction_sample_pdf + '">')
        str_list.append('<img alt="Mapped Fraction - Sample"')
        str_list.append(' src="' + file_path_summary.mapped_fraction_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.mapped_fraction_read_group_pdf + '">')
        str_list.append('<img alt="Mapped Fraction - Read Group"')
        str_list.append(' src="' + file_path_summary.mapped_fraction_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Mapped Fraction</td>\n')
        str_list.append('</tr>\n')

        # Mapped Number Plots
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.mapped_number_sample_pdf + '">')
        str_list.append('<img alt="Mapped Number - Sample"')
        str_list.append(' src="' + file_path_summary.mapped_number_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.mapped_number_read_group_pdf + '">')
        str_list.append('<img alt="Mapped Number - Read Group"')
        str_list.append(' src="' + file_path_summary.mapped_number_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Mapped Number</td>\n')
        str_list.append('</tr>\n')

        # Junction Fraction Plots
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.junction_fraction_sample_pdf + '">')
        str_list.append('<img alt="Junction Fraction - Sample"')
        str_list.append(' src="' + file_path_summary.junction_fraction_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.junction_fraction_read_group_pdf + '">')
        str_list.append('<img alt="Junction Fraction - Read Group"')
        str_list.append(' src="' + file_path_summary.junction_fraction_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Junction Fraction</td>\n')
        str_list.append('</tr>\n')

        # Junction Number Plots
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.junction_number_sample_pdf + '">')
        str_list.append('<img alt="Junction Number - Sample"')
        str_list.append(' src="' + file_path_summary.junction_number_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.junction_number_read_group_pdf + '">')
        str_list.append('<img alt="Junction Number - Read Group"')
        str_list.append(' src="' + file_path_summary.junction_number_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Junction Number</td>\n')
        str_list.append('</tr>\n')

        # Summary Tables
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.table_sample_tsv + '">')
        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.table_read_group_tsv + '">')
        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Summary</td>\n')
        str_list.append('</tr>\n')

        # The table rows below are a copy of the Aligner.report() method.

        # Alignment Summary Plot
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_summary_sample_pdf + '">')
        str_list.append('<img alt="Alignment Summary - Sample"')
        str_list.append(' src="' + file_path_summary.pasm_summary_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_summary_read_group_pdf + '">')
        str_list.append('<img alt="Alignment Summary - Read Group"')
        str_list.append(' src="' + file_path_summary.pasm_summary_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Alignment Summary</td>\n')
        str_list.append('</tr>\n')

        # Alignment Summary Plot Mapped Numbers
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_numbers_sample_pdf + '">')
        str_list.append('<img alt="Mapped Numbers - Sample"')
        str_list.append(' src="' + file_path_summary.pasm_numbers_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_numbers_read_group_pdf + '">')
        str_list.append('<img alt="Mapped Numbers - Read Group"')
        str_list.append(' src="' + file_path_summary.pasm_numbers_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Mapped Numbers</td>\n')
        str_list.append('</tr>\n')

        # Alignment Summary Plot Mapped Fractions
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_fractions_sample_pdf + '">')
        str_list.append('<img alt="Mapped Fractions - Sample"')
        str_list.append(' src="' + file_path_summary.pasm_fractions_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_fractions_read_group_pdf + '">')
        str_list.append('<img alt="Mapped Fractions - Read Group"')
        str_list.append(' src="' + file_path_summary.pasm_fractions_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Mapped Fractions</td>\n')
        str_list.append('</tr>\n')

        # Alignment Summary Plot Strand Balance
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_strand_balance_sample_pdf + '">')
        str_list.append('<img alt="Strand Balance - Sample"')
        str_list.append(' src="' + file_path_summary.pasm_strand_balance_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_strand_balance_read_group_pdf + '">')
        str_list.append('<img alt="Strand Balance - Read Group"')
        str_list.append(' src="' + file_path_summary.pasm_strand_balance_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Strand Balance</td>\n')
        str_list.append('</tr>\n')

        # Alignment Summary Metrics Tables
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_table_sample_tsv + '">')
        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_table_read_group_tsv + '">')
        str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Alignment Summary</td>\n')
        str_list.append('</tr>\n')

        if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_levels_pdf)):
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pdm_levels_pdf + '">')
            if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_levels_png)):
                str_list.append('<img alt="Duplication Levels - Sample"')
                str_list.append(' src="' + file_path_summary.pdm_levels_png + '"')
                str_list.append(' height="100" width="100" />')
            else:
                str_list.append('PDF')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Duplication Levels</td>\n')
            str_list.append('</tr>\n')

        if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_fractions_pdf)):
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pdm_fractions_pdf + '">')
            if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_fractions_png)):
                str_list.append('<img alt="Duplication Percentage - Sample"')
                str_list.append(' src="' + file_path_summary.pdm_fractions_png + '"')
                str_list.append(' height="100" width="100" />')
            else:
                str_list.append('PDF')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Duplication Percentage</td>\n')
            str_list.append('</tr>\n')

        if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_table_tsv)):
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pdm_table_tsv + '">')
            str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Duplication Summary</td>\n')
            str_list.append('</tr>\n')

        str_list.append('</tbody>\n')
        str_list.append('</table>\n')
        str_list.append('\n')

        # Add the sample table.
        str_list.extend(self.report_html_sample())

        # Save the HTML report.
        self.report_to_file(content=str_list)

        # Save a UCSC Track Hub with alignment and coverage composite tracks.
        self.ucsc_hub_to_file(content=self.report_hub_alignment() + self.report_hub_coverage())

        return
