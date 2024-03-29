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
"""The :py:mod:`bsf.analyses.aligner` module provides classes supporting Aligner analyses.
"""
import logging
import os
from typing import Optional

from bsf.analysis import Analysis, Stage
from bsf.annotation import AnnotationSheet
from bsf.ngs import Collection, PairedReads, Sample
from bsf.procedure import FilePath, ConcurrentRunnable, ConsecutiveRunnable
from bsf.process import RunnableStepMakeDirectory, RunnableStepMakeNamedPipe, RunnableStepPicard, \
    RunnableStepMove, RunnableStep, RunnableStepLink
from bsf.standards import SafeFileName, Configuration, StandardFilePath, JavaArchive

module_logger = logging.getLogger(name=__name__)


class FilePathAlign(FilePath):
    """The :py:class:`bsf.analyses.aligner.FilePathAlign` class models file paths at the alignment stage.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar stderr_txt: Text file to capture :literal:`STDERR` of the aligner
    :type stderr_txt: str
    :ivar stdout_txt: Text file to capture :literal:`STDOUT` of the aligner
    :type stdout_txt: str
    :ivar aligned_sam: An aligned sequence alignment map (SAM) file path.
    :type aligned_sam: str
    :ivar cleaned_bam: A cleaned sequence alignment map (SAM) file path.
    :type cleaned_bam: str
    :ivar aligned_bam: An aligned binary alignment map (BAM) file path.
    :type aligned_bam: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.aligner.FilePathAlign` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.stderr_txt = os.path.join(prefix, prefix + '_stderr.txt')
        self.stdout_txt = os.path.join(prefix, prefix + '_stdout.txt')

        self.aligned_sam = os.path.join(prefix, prefix + '_aligned.sam')
        self.cleaned_bam = os.path.join(prefix, prefix + '_cleaned.bam')
        self.aligned_bam = os.path.join(prefix, prefix + '_aligned.bam')

        return


class FilePathReadGroup(FilePath):
    """The :py:class:`bsf.analyses.aligner.FilePathReadGroup` class models file paths at the
    read group processing stage.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar merged_bam: A merged binary alignment map (BAM) file path.
    :type merged_bam: str
    :ivar sorted_bai: A sorted binary alignment map index (BAI) file path.
    :type sorted_bai: str
    :ivar sorted_bam: A sorted binary alignment map (BAM) file path.
    :type sorted_bam: str
    :ivar sorted_md5: A sorted binary alignment map (BAM) MD5 checksum file path.
    :type sorted_md5: str
    :ivar read_group_bai: A read group binary alignment map index (BAI) file path alias.
    :type read_group_bai: str
    :ivar read_group_bam: A read group binary alignment map (BAM) file path alias.
    :type read_group_bam: str
    :ivar read_group_md5: A read group binary alignment map (BAM) MD5 checksum file path alias.
    :type read_group_md5: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.aligner.FilePathReadGroup` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.merged_bam = os.path.join(prefix, prefix + '_merged.bam')

        self.sorted_bai = os.path.join(prefix, prefix + '.bai')
        self.sorted_bam = os.path.join(prefix, prefix + '.bam')
        self.sorted_md5 = os.path.join(prefix, prefix + '.bam.md5')

        self.read_group_bai = self.sorted_bai
        self.read_group_bam = self.sorted_bam
        self.read_group_md5 = self.sorted_md5

        return


class FilePathSample(FilePath):
    """The :py:class:`bsf.analyses.aligner.FilePathSample` class models file paths at the sample processing stage.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar merged_bam: A merged binary alignment map (BAM) file path.
    :type merged_bam: str
    :ivar duplicate_bam: A duplicate-marked binary alignment map (BAM) file path.
    :type duplicate_bam: str
    :ivar duplicate_metrics_tsv: A Picard Duplicate Metrics TSV file path.
    :type duplicate_metrics_tsv: str
    :ivar sample_bai_link_source: A symbolic link source of the binary alignment map index (BAI) file path.
    :type sample_bai_link_source: str
    :ivar sample_bai_link_target: A symbolic link target of the binary alignment map index (BAI) file path.
    :type sample_bai_link_target: str
    :ivar alignment_summary_metrics_tsv: A Picard Alignment Summary Metrics TSV file path.
    :type alignment_summary_metrics_tsv: str
    :ivar sample_bai: A sample binary alignment map index (BAI) file path.
    :type sample_bai: str
    :ivar sample_bam: A sample binary alignment map (BAM) file path.
    :type sample_bam: str
    :ivar sample_md5: A sample binary alignment map (BAM) MD5 checksum file path.
    :type sample_md5: str
    :ivar sample_wig: A sample wiggle file path.
    :type sample_wig: str
    :ivar sample_bw: A sample bigWig file path.
    :type sample_bw: str
    :ivar sample_bwi: A sample bigWig info file path.
    :type sample_bwi: str
    :ivar prefix_prefix: A double prefix for RSeQC :literal:`bam2wig.py`.
    :type prefix_prefix: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.aligner.FilePathSample` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathSample, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.merged_bam = os.path.join(prefix, prefix + '_merged.bam')

        self.duplicate_bam = os.path.join(prefix, prefix + '_duplicate.bam')

        self.duplicate_metrics_tsv = os.path.join(prefix, prefix + '_duplicate_metrics.tsv')

        self.sample_bai = os.path.join(prefix, prefix + '.bai')
        self.sample_bam = os.path.join(prefix, prefix + '.bam')
        self.sample_md5 = os.path.join(prefix, prefix + '.bam.md5')
        self.sample_wig = os.path.join(prefix, prefix + '.wig')
        self.sample_bw = os.path.join(prefix, prefix + '.bw')
        self.sample_bwi = os.path.join(prefix, prefix + '_bwi.txt')

        self.sample_bai_link_source = prefix + '.bai'
        self.sample_bai_link_target = os.path.join(prefix, prefix + '.bam.bai')

        self.alignment_summary_metrics_tsv = os.path.join(prefix, '_'.join((prefix, 'alignment_summary_metrics.tsv')))

        self.prefix_prefix = os.path.join(prefix, prefix)

        return


class FilePathSummary(FilePath):
    """The :py:class:`bsf.analyses.aligner.FilePathSummary` class models file paths at the summary stage.
    They are mainly based on `Picard Alignment Summary Metrics
    <http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics>`_ and
    `Picard Duplication Metrics
    <http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics>`_ reports.

    :ivar output_directory: An output directory path.
    :type output_directory: str
    :ivar pasm_summary_read_group_pdf: A Picard Alignment Summary Metrics per read group PDF plot file path.
    :type pasm_summary_read_group_pdf: str
    :ivar pasm_summary_read_group_png: A Picard Alignment Summary Metrics per read group PNG plot file path.
    :type pasm_summary_read_group_png: str
    :ivar pasm_summary_sample_pdf: A Picard Alignment Summary Metrics per sample PDF plot file path.
    :type pasm_summary_sample_pdf: str
    :ivar pasm_summary_sample_png: A Picard Alignment Summary Metrics per sample PNG plot file path.
    :type pasm_summary_sample_png: str
    :ivar pasm_numbers_read_group_pdf: A Picard Alignment Summary Metrics numbers per read group PDF plot file path.
    :type pasm_numbers_read_group_pdf: str
    :ivar pasm_numbers_read_group_png: A Picard Alignment Summary Metrics numbers per read group PNG plot file path.
    :type pasm_numbers_read_group_png: str
    :ivar pasm_numbers_sample_pdf: A Picard Alignment Summary Metrics numbers per sample PDF plot file path.
    :type pasm_numbers_sample_pdf: str
    :ivar pasm_numbers_sample_png: A Picard Alignment Summary Metrics numbers per sample PNG plot file path.
    :type pasm_numbers_sample_png: str
    :ivar pasm_fractions_read_group_pdf: A Picard Alignment Summary Metrics fractions per read group
        PDF plot file path.
    :type pasm_fractions_read_group_pdf: str
    :ivar pasm_fractions_read_group_png: A Picard Alignment Summary Metrics fractions per read group
        PNG plot file path.
    :type pasm_fractions_read_group_png: str
    :ivar pasm_fractions_sample_pdf: A Picard Alignment Summary Metrics fractions per sample PDF plot file path.
    :type pasm_fractions_sample_pdf: str
    :ivar pasm_fractions_sample_png: A Picard Alignment Summary Metrics fractions per sample PNG plot file path.
    :type pasm_fractions_sample_png: str
    :ivar pasm_strand_balance_read_group_pdf: A Picard Alignment Summary Metrics strand balance per read group
        PDF plot file path.
    :type pasm_strand_balance_read_group_pdf: str
    :ivar pasm_strand_balance_read_group_png: A Picard Alignment Summary Metrics strand balance per read group
        PNG plot file path.
    :type pasm_strand_balance_read_group_pdf: str
    :ivar pasm_strand_balance_sample_pdf: A Picard Alignment Summary Metrics strand balance per sample
        PDF plot file path.
    :type pasm_strand_balance_sample_pdf: str
    :ivar pasm_strand_balance_sample_png: A Picard Alignment Summary Metrics strand balance per sample
        PDF plot file path.
    :type pasm_strand_balance_sample_png: str
    :ivar pasm_table_read_group_tsv: A Picard Alignment Summary Metrics table per read group file path.
    :type pasm_table_read_group_tsv: str
    :ivar pasm_table_sample_tsv: A Picard Alignment Summary Metrics TSV table per sample file path.
    :type pasm_table_sample_tsv: str
    :ivar pdm_levels_pdf: A Picard Duplication Metrics levels per sample PDF plot file path.
    :type pdm_levels_pdf: str
    :ivar pdm_levels_png: A Picard Duplication Metrics levels per sample PNG plot file path.
    :type pdm_levels_png: str
    :ivar pdm_fractions_pdf: A Picard Duplication Metrics percentage per sample PDF plot file path.
    :type pdm_fractions_pdf: str
    :ivar pdm_fractions_png: A Picard Duplication Metrics percentage per sample PNG plot file path.
    :type pdm_fractions_png: str
    :ivar pdm_table_tsv: A Picard Duplication Metrics TSV table per sample file path.
    :type pdm_table_tsv: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.aligner.FilePathSummary` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathSummary, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.read_group_to_sample_tsv = prefix + '_read_group_to_sample.tsv'

        # Picard Alignment Summary Metrics (PASM)

        self.pasm_summary_read_group_pdf = prefix + '_pasm_summary_read_group.pdf'
        self.pasm_summary_read_group_png = prefix + '_pasm_summary_read_group.png'

        self.pasm_summary_sample_pdf = prefix + '_pasm_summary_sample.pdf'
        self.pasm_summary_sample_png = prefix + '_pasm_summary_sample.png'

        self.pasm_numbers_read_group_pdf = prefix + '_pasm_numbers_read_group.pdf'
        self.pasm_numbers_read_group_png = prefix + '_pasm_numbers_read_group.png'

        self.pasm_numbers_sample_pdf = prefix + '_pasm_numbers_sample.pdf'
        self.pasm_numbers_sample_png = prefix + '_pasm_numbers_sample.png'

        self.pasm_fractions_read_group_pdf = prefix + '_pasm_fractions_read_group.pdf'
        self.pasm_fractions_read_group_png = prefix + '_pasm_fractions_read_group.png'

        self.pasm_fractions_sample_pdf = prefix + '_pasm_fractions_sample.pdf'
        self.pasm_fractions_sample_png = prefix + '_pasm_fractions_sample.png'

        self.pasm_strand_balance_read_group_pdf = prefix + '_pasm_strand_balance_read_group.pdf'
        self.pasm_strand_balance_read_group_png = prefix + '_pasm_strand_balance_read_group.png'

        self.pasm_strand_balance_sample_pdf = prefix + '_pasm_strand_balance_sample.pdf'
        self.pasm_strand_balance_sample_png = prefix + '_pasm_strand_balance_sample.png'

        self.pasm_table_read_group_tsv = prefix + '_pasm_table_read_group.tsv'

        self.pasm_table_sample_tsv = prefix + '_pasm_table_sample.tsv'

        # Picard Duplication Metrics (PDM)

        self.pdm_levels_pdf = prefix + '_pdm_levels.pdf'
        self.pdm_levels_png = prefix + '_pdm_levels.png'

        self.pdm_fractions_pdf = prefix + '_pdm_fractions.pdf'
        self.pdm_fractions_png = prefix + '_pdm_fractions.png'

        self.pdm_table_tsv = prefix + '_pdm_table.tsv'

        return


class Aligner(Analysis):
    """The :py:class:`bsf.analyses.aligner.Aligner` class represents the logic to run a (short read) aligner.

    :cvar name: :py:attr:`bsf.analysis.Analysis.name` that should be overridden by subclasses.
    :type name: str
    :cvar prefix: :py:attr:`bsf.analysis.Analysis.prefix` that should be overridden by subclasses.
    :type prefix: str
    :ivar genome_fasta: A genome :emphasis:`FASTA` file path.
    :type genome_fasta: str | None
    :ivar genome_index: A genome index file path.
    :type genome_index: str | None
    :ivar skip_mark_duplicates: Request skipping the Picard :literal:`MarkDuplicates` step.
    :type skip_mark_duplicates: bool | None
    :ivar java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
    :type java_archive_picard: str | None
    """

    name = 'Aligner Analysis'
    prefix = 'aligner'

    sam_attributes_to_retain_list: list[str] = []

    @classmethod
    def get_stage_name_align(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'align'))

    @classmethod
    def get_stage_name_read_group(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_stage_name_sample(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'sample'))

    @classmethod
    def get_stage_name_summary(cls) -> str:
        """Get a particular :py:attr:`bsf.analysis.Stage.name` attribute.

        :return: A :py:attr:`bsf.analysis.Stage.name` attribute.
        :rtype: str
        """
        return '_'.join((cls.prefix, 'summary'))

    @classmethod
    def get_prefix_align(cls, paired_reads_name) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param paired_reads_name: A :py:attr:`bsf.ngs.PairedReads.name` attribute.
        :type paired_reads_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_align(), paired_reads_name))

    @classmethod
    def get_prefix_read_group(cls, read_group_name) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_prefix_sample(cls, sample_name) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return '_'.join((cls.get_stage_name_sample(), sample_name))

    @classmethod
    def get_prefix_summary(cls) -> str:
        """Get a Python :py:class:`str` prefix representing a :py:class:`bsf.procedure.Runnable` object.

        :return: A Python :py:class:`str` (prefix)  object representing a :py:class:`bsf.procedure.Runnable` object.
        :rtype: str
        """
        return cls.get_stage_name_summary()

    @classmethod
    def get_file_path_align(cls, paired_reads_name) -> FilePathAlign:
        """Get a :py:class:`bsf.analyses.aligner.FilePathAlign` object from this or a subclass.

        :param paired_reads_name: A :py:attr:`bsf.ngs.PairedReads.name` attribute.
        :type paired_reads_name: str
        :return: A :py:class:`bsf.analyses.aligner.FilePathAlign` object or subclass thereof.
        :rtype: FilePathAlign
        """
        return FilePathAlign(prefix=cls.get_prefix_align(paired_reads_name=paired_reads_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name) -> FilePathReadGroup:
        """Get a :py:class:`bsf.analyses.aligner.FilePathReadGroup` object from this or a subclass.

        :param read_group_name: A read group name.
        :type read_group_name: str
        :return: A :py:class:`bsf.analyses.aligner.FilePathReadGroup` object or subclass thereof.
        :rtype: FilePathReadGroup
        """
        return FilePathReadGroup(prefix=cls.get_prefix_read_group(read_group_name=read_group_name))

    @classmethod
    def get_file_path_sample(cls, sample_name) -> FilePathSample:
        """Get a :py:class:`bsf.analyses.aligner.FilePathSample` object from this or a subclass.

        :param sample_name: A :py:attr:`bsf.ngs.Sample.name` attribute.
        :type sample_name: str
        :return: A :py:class:`bsf.analyses.aligner.FilePathSample` object or subclass thereof.
        :rtype: FilePathSample
        """
        return FilePathSample(prefix=cls.get_prefix_sample(sample_name=sample_name))

    @classmethod
    def get_file_path_summary(cls) -> FilePathSummary:
        """Get a :py:class:`bsf.analyses.aligner.FilePathSummary` object from this or a subclass.

        :return: A :py:class:`bsf.analyses.aligner.FilePathSummary` object or subclass thereof.
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
            genome_fasta: Optional[str] = None,
            genome_index: Optional[str] = None,
            skip_mark_duplicates: Optional[bool] = None,
            java_archive_picard: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.analyses.aligner.Aligner` object.

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
        :param genome_fasta: A genome :emphasis:`FASTA` file path.
        :type genome_fasta: str | None
        :param genome_index: A genome index file path.
        :type genome_index: str | None
        :param skip_mark_duplicates: Request skipping the Picard :literal:`MarkDuplicates` step.
        :type skip_mark_duplicates: bool | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        """
        super(Aligner, self).__init__(
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
            sample_list=sample_list)

        # Sub-class specific ...

        self.genome_fasta = genome_fasta
        self.genome_index = genome_index
        self.skip_mark_duplicates = skip_mark_duplicates
        self.java_archive_picard = java_archive_picard

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.aligner.Aligner` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(Aligner, self).set_configuration(configuration=configuration, section=section)

        config_parser = configuration.config_parser

        # Get the genome database.

        option = 'genome_fasta'
        if config_parser.has_option(section=section, option=option):
            self.genome_fasta = config_parser.get(section=section, option=option)

        option = 'genome_index'
        if config_parser.has_option(section=section, option=option):
            self.genome_index = config_parser.get(section=section, option=option)

        option = 'skip_mark_duplicates'
        if config_parser.has_option(section=section, option=option):
            self.skip_mark_duplicates = config_parser.getboolean(section=section, option=option)

        option = 'java_archive_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.java_archive_picard = configuration.config_parser.get(section=section, option=option)

        return

    def add_runnable_step_aligner(
            self,
            runnable_align: ConcurrentRunnable,
            stage_align: Stage,
            file_path_1: Optional[str],
            file_path_2: Optional[str]) -> None:
        """Add one or more Aligner-specific :py:class:`bsf.process.RunnableStep` objects to the
        :py:class:`bsf.procedure.ConcurrentRunnable` object.

        :param runnable_align: A :py:class:`bsf.procedure.ConcurrentRunnable` object.
        :type runnable_align: ConcurrentRunnable
        :param stage_align: A :py:class:`bsf.analysis.Stage` object.
        :type stage_align: Stage
        :param file_path_1: A :emphasis:`FASTQ` file path 1.
        :type file_path_1: str | None
        :param file_path_2: A :emphasis:`FASTQ` file path 2.
        :type file_path_2: str | None
        """
        return

    def add_runnable_step_sample(self, runnable_sample: ConsecutiveRunnable, stage_sample: Stage) -> None:
        """Add one or more Aligner-specific :py:class:`bsf.process.RunnableStep` objects
        to the :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable_sample: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :type runnable_sample: ConsecutiveRunnable
        :param stage_sample: A :py:class:`bsf.analysis.Stage` object.
        :type stage_sample: Stage
        """
        return

    def add_runnable_step_summary(self, runnable_summary: ConsecutiveRunnable, stage_summary: Stage) -> None:
        """Add one or more Aligner-specific :py:class:`bsf.process.RunnableStep` objects
        to the :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable_summary: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :type runnable_summary: ConsecutiveRunnable
        :param stage_summary: A :py:class:`bsf.analysis.Stage` object.
        :type stage_summary: Stage
        """
        return

    def run(self) -> None:
        """Run a :py:class:`bsf.analyses.aligner.Aligner` object.

        - For each :py:class:`bsf.ngs.Sample` object, which :py:class:`bsf.ngs.PairedReads` objects
          have not been excluded:

          - For each :py:class:`bsf.ngs.PairedReads` object, which has not been excluded:

            - Run an aligner and :emphasis:`Picard CleanSam` concurrently.

          - For each :py:class:`bsf.ngs.Sample`-specific unaligned BAM file:

            - For a single aligned BAM file per unaligned BAM file, the aligned BAM file is just renamed.
            - For more than one aligned BAM file per unaligned BAM file,
              all aligned BAM files are merged via :emphasis:`Picard MergeSamFiles` in query name order.
            - If an unaligned BAM file is available, the merged aligned BAM files are annotated
              with information from the unaligned BAM file via :emphasis:`Picard MergeBamAlignment`.
            - If an unaligned BAM file is not available, the merged aligned BAM file is just renamed.

          - If duplicates should *not* be marked:

            - For a single merged aligned BAM file, just rename the single read group-specific aligned BAM file.
            - For more than one read group-specific merged aligned BAM file,
              run :emphasis:`Picard MergeSamFile` to simply merge them into one sample-specific aligned BAM file.

          - If duplicates should be marked,
            run :emphasis:`Picard MarkDuplicates` on all read group-specific merged aligned BAM files at once.
          - Run :emphasis:`Picard SortSam` by coordinate.
          - Create a symbolic link from the Picard-style <file_name>.bai file to a
            samtools-style <file_name>.bam.bai file.
          - Run the :emphasis:`Picard CollectAlignmentSummaryMetrics` analysis.

        - Run the :emphasis:`bsfR` alignment summary script.
        - Run aligner-specific :emphasis:`bsfR` summary scripts.
        """

        def run_read_comparisons() -> None:
            """Private function to read a :py:class:`bsf.annotation.AnnotationSheet` specifying comparisons
            from a CSV file path.

            This implementation just adds all :py:class:`bsf.ngs.Sample` objects from the
            :py:attr:`bsf.analysis.Analysis.collection` instance variable
            (i.e., :py:class:`bsf.ngs.Collection` object) to the
            :py:attr:`bsf.analysis.Analysis.sample_list` instance variable.
            """
            self.sample_list.extend(self.collection.get_all_samples())

            return

        def run_get_unmapped_bam(_paired_reads: PairedReads) -> tuple[str, Optional[str]]:
            """Get the unmapped BAM file annotation of a :py:class:`bsf.ngs.PairedReads` object.

            :param _paired_reads: A :py:class:`bsf.ngs.PairedReads` object.
            :type _paired_reads: PairedReads
            :return: A Python :py:class:`tuple` object of
                Python :py:class:`str` (base name) and Python :py:class:`str` (unmapped BAM file path) objects.
            :rtype: (str, str | None)
            """
            if 'BAM File' in _paired_reads.annotation_dict and _paired_reads.annotation_dict['BAM File']:
                # Consider only the first list component.
                bam_name, bam_extension = os.path.splitext(
                    os.path.basename(_paired_reads.annotation_dict['BAM File'][0]))

                bam_name = SafeFileName.get_safe_file_name(file_name=bam_name)

                return bam_name, _paired_reads.annotation_dict['BAM File'][0]
            else:
                module_logger.warning("PairedReads object %r without 'BAM File' annotation.", _paired_reads.get_name())
                return _paired_reads.get_name(), None

        # Start of the run() method body.

        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Get the sample annotation sheet before calling the run() method of the bsf.analysis.Analysis super-class.

        if self.sas_file:
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)
            if not os.path.exists(self.sas_file):
                raise Exception(f'The sample annotation sheet {self.sas_file!r} does not exist.')
        else:
            self.sas_file = self.get_annotation_file(prefix_list=[self.prefix], suffix='samples.csv')
            if not self.sas_file:
                raise Exception('No suitable sample annotation sheet in the current working directory.')

        # The Aligner analysis requires a genome version.

        if not self.genome_version:
            raise Exception(f"A {self.name!s} requires a 'genome_version' configuration option.")

        super(Aligner, self).run()

        if not self.genome_fasta:
            # Set the genome_index to None, to point to the genome directory that contains the
            # Picard sequence dictionary.
            self.genome_fasta = StandardFilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index=None)

        if not self.java_archive_picard:
            self.java_archive_picard = JavaArchive.get_picard()
            if not self.java_archive_picard:
                raise Exception(f"A {self.name!s} requires a 'java_archive_picard' configuration option.")

        run_read_comparisons()

        stage_align = self.get_stage(name=self.get_stage_name_align())
        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())
        stage_sample = self.get_stage(name=self.get_stage_name_sample())
        stage_summary = self.get_stage(name=self.get_stage_name_summary())

        file_path_summary = self.get_file_path_summary()

        # Create an annotation sheet linking sample name and read group name, which is required for the
        # summary script.

        annotation_sheet = AnnotationSheet(
            file_path=os.path.join(self.genome_directory, file_path_summary.read_group_to_sample_tsv),
            file_type='excel-tab',
            name='star_read_group',
            field_name_list=['sample', 'read_group'])

        runnable_sample_list: list[ConsecutiveRunnable] = list()

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            module_logger.debug('Sample.name: %r', sample.name)
            module_logger.log(logging.DEBUG - 2, 'Sample: %r', sample)

            # To run Picard MergeBamAlignment, all alignments from a BAM file need merging into one.

            unmapped_bam_file_dict: dict[str, tuple[Optional[str], list[ConcurrentRunnable]]] = dict()

            runnable_read_group_list: list[ConsecutiveRunnable] = list()

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            # Record whether the sample has UMIs on the basis of the PairedReads Structure annotation.
            sample_umis = False

            for paired_reads_name in sorted(paired_reads_dict):
                if not paired_reads_dict[paired_reads_name]:
                    # Skip PairedReads.name keys, which PairedReads objects have all been excluded.
                    continue

                if len(paired_reads_dict[paired_reads_name]) > 1:
                    raise Exception('Cannot align and process more than one PairedReads object at a time.')

                # Record the Runnable for the Picard MergeBamAlignment step.
                paired_reads = paired_reads_dict[paired_reads_name][0]

                # Get the file paths for Reads1 and Reads2 and check for FASTQ files.
                if paired_reads.reads_1 is None:
                    raise Exception(f'A {self.name!s} requires a Reads1 object.')
                else:
                    if paired_reads.reads_1.file_path is not None:
                        if not (paired_reads.reads_1.file_path.endswith('.fastq') or
                                paired_reads.reads_1.file_path.endswith('.fastq.gz')):
                            raise Exception(f'A {self.name!s} requires a (GNU Zip compressed) FASTQ file.')
                    file_path_1 = paired_reads.reads_1.file_path

                if paired_reads.reads_2 is None:
                    file_path_2 = None
                else:
                    if paired_reads.reads_2.file_path is not None:
                        if not (paired_reads.reads_2.file_path.endswith('.fastq') or
                                paired_reads.reads_2.file_path.endswith('.fastq.gz')):
                            raise Exception(f'A {self.name!s} requires a (GNU Zip compressed) FASTQ file.')
                    file_path_2 = paired_reads.reads_2.file_path

                # Check the PairedReads Structure annotation for UMIs (i.e., M).

                if 'Structure' in paired_reads.annotation_dict:
                    for structure_str in paired_reads.annotation_dict['Structure']:
                        if 'M' in structure_str:
                            sample_umis = True

                ###################
                # Alignment Stage #
                ###################

                annotation_sheet.row_dict_list.append({'sample': sample.name, 'read_group': paired_reads.get_name()})

                # Create a Runnable and Executable for alignment and processing.

                file_path_align = self.get_file_path_align(paired_reads_name=paired_reads_name)

                # Create a new ConcurrentRunnable to run the aligner and Picard CleanSam, concurrently.

                runnable_align = self.add_runnable(
                    runnable=ConcurrentRunnable(
                        name=self.get_prefix_align(paired_reads_name=paired_reads_name),
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory))
                self.set_stage_runnable(
                    stage=stage_align,
                    runnable=runnable_align)

                # Make a directory via RunnableStepMakeDirectory.

                runnable_step = RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=file_path_align.output_directory)
                runnable_align.add_runnable_step_prologue(runnable_step=runnable_step)

                # Make named pipes via RunnableStepMakeNamedPipe.

                runnable_step = RunnableStepMakeNamedPipe(
                    name='make_fifo_aligned_sam',
                    file_path=file_path_align.aligned_sam)
                runnable_align.add_runnable_step_prologue(runnable_step=runnable_step)

                runnable_step = RunnableStepMakeNamedPipe(
                    name='make_fifo_cleaned_bam',
                    file_path=file_path_align.cleaned_bam)
                runnable_align.add_runnable_step_prologue(runnable_step=runnable_step)

                # Start the programs in reverse order so that they do not block while opening their named pipe(s).

                # Run Picard SortSam to convert the cleaned BAM file into a query name-sorted BAM file.
                # Not every aligner may return reads in the order of the unaligned BAM (FASTQ) files.

                runnable_step = RunnableStepPicard(
                    name='picard_sort_sam',
                    obsolete_file_path_list=[
                        file_path_align.cleaned_bam,
                    ],
                    java_temporary_path=runnable_align.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='SortSam')
                runnable_align.add_runnable_step(runnable_step=runnable_step)

                # Required Arguments
                # INPUT []
                runnable_step.add_picard_option(key='INPUT', value=file_path_align.cleaned_bam)
                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_align.aligned_bam)
                # SORT_ORDER []
                runnable_step.add_picard_option(key='SORT_ORDER', value='queryname')

                # Optional Common Arguments
                # COMPRESSION_LEVEL [5]
                # CREATE_INDEX [false]
                # CREATE_MD5_FILE [false]
                # MAX_RECORDS_IN_RAM [500000]
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                # QUIET [false]
                # REFERENCE_SEQUENCE [null]
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_align.temporary_directory_path(absolute=False))
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]
                # VALIDATION_STRINGENCY [STRICT]
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

                # Run Picard CleanSam to split alignments at sequence boundaries and
                # set mapping quality of unmapped reads to 0.

                runnable_step = RunnableStepPicard(
                    name='picard_clean_sam',
                    obsolete_file_path_list=[
                        file_path_align.aligned_sam,
                    ],
                    java_temporary_path=runnable_align.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx2G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='CleanSam')
                runnable_align.add_runnable_step(runnable_step=runnable_step)

                # Required Arguments
                # INPUT []
                runnable_step.add_picard_option(key='INPUT', value=file_path_align.aligned_sam)
                # OUTPUT []
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_align.cleaned_bam)

                # Optional Common Arguments
                # COMPRESSION_LEVEL [5]
                # CREATE_INDEX [false]
                # CREATE_MD5_FILE [false]
                # MAX_RECORDS_IN_RAM [500000]
                # QUIET [false]
                # REFERENCE_SEQUENCE [null]
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_align.temporary_directory_path(absolute=False))
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]
                # VALIDATION_STRINGENCY [STRICT]
                runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='0')
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

                # Start the Aligner.

                self.add_runnable_step_aligner(
                    runnable_align=runnable_align,
                    stage_align=stage_align,
                    file_path_1=file_path_1,
                    file_path_2=file_path_2)

                # Record the alignment Runnable under the unmapped BAM file name
                # to merge the sub-alignments of all (trimmed) FASTQ files that
                # resulted from an unmapped BAM file.
                bam_file_name, bam_file_path = run_get_unmapped_bam(_paired_reads=paired_reads)

                if bam_file_name in unmapped_bam_file_dict:
                    unmapped_bam_file_dict[bam_file_name][1].append(runnable_align)
                else:
                    unmapped_bam_file_dict[bam_file_name] = (bam_file_path, [runnable_align])

            ###############################
            # Read Group Processing Stage #
            ###############################

            # Merge all aligned BAM files of each read group of an unaligned BAM file.

            for bam_file_name, (bam_file_path, runnable_align_list) in unmapped_bam_file_dict.items():
                # Create a Runnable and Executable object for merging each read group.

                file_path_read_group = self.get_file_path_read_group(read_group_name=bam_file_name)

                runnable_read_group = self.add_runnable(
                    runnable=ConsecutiveRunnable(
                        name=self.get_prefix_read_group(read_group_name=bam_file_name),
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory))
                executable_read_group = self.set_stage_runnable(
                    stage=stage_read_group,
                    runnable=runnable_read_group)

                # Record the dependency for all alignment Runnable objects.
                for runnable_align in runnable_align_list:
                    executable_read_group.dependencies.append(runnable_align.name)

                runnable_read_group_list.append(runnable_read_group)

                runnable_step = RunnableStepMakeDirectory(
                    name='make_directory',
                    directory_path=file_path_read_group.output_directory)
                runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                if len(runnable_align_list) == 1:
                    # For a single aligned BAM file per unaligned BAM file,
                    # the aligned BAM file is just renamed.
                    runnable_align = runnable_align_list[0]
                    # NOTE: Since the Runnable.name already contains the prefix, FilePathAlign() has to be used.
                    file_path_align = FilePathAlign(prefix=runnable_align.name)
                    runnable_step = RunnableStepMove(
                        name='move_aligned_bam',
                        source_path=file_path_align.aligned_bam,
                        target_path=file_path_read_group.merged_bam)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)
                else:
                    # For more than one aligned BAM file per unaligned BAM file,
                    # all aligned BAM files are merged via Picard MergeSamFiles in query name order.
                    runnable_step = RunnableStepPicard(
                        name='picard_merge_sam_files',
                        java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_picard,
                        picard_command='MergeSamFiles')
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    # Required Arguments
                    # INPUT Required
                    for runnable_align in runnable_align_list:
                        # NOTE: Since the Runnable.name already contains the prefix, FilePathAlign() has to be used.
                        file_path_align = FilePathAlign(prefix=runnable_align.name)
                        runnable_step.add_picard_option(
                            key='INPUT',
                            value=file_path_align.aligned_bam,
                            override=True)
                        runnable_step.obsolete_file_path_list.append(file_path_align.aligned_bam)
                    # OUTPUT Required
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.merged_bam)

                    # Optional Tool Arguments
                    # ASSUME_SORTED [false]
                    # COMMENT [null]
                    # INTERVALS [null]
                    # MERGE_SEQUENCE_DICTIONARIES [false]
                    # SORT_ORDER [coordinate]
                    runnable_step.add_picard_option(key='SORT_ORDER', value='queryname')
                    # USE_THREADING [false]
                    runnable_step.add_picard_option(key='USE_THREADING', value='true')

                    # Optional Common Arguments
                    # COMPRESSION_LEVEL [5]
                    # CREATE_INDEX [false]
                    # CREATE_MD5_FILE [false]
                    # MAX_RECORDS_IN_RAM [500000]
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                    # QUIET [false]
                    # REFERENCE_SEQUENCE [null]
                    # TMP_DIR [null]
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_read_group.temporary_directory_path(absolute=False))
                    # USE_JDK_DEFLATER [false]
                    # USE_JDK_INFLATER [false]
                    # VALIDATION_STRINGENCY [STRICT]
                    # VERBOSITY [INFO]
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

                if bam_file_path:
                    # If an unaligned BAM file is available, the merged aligned BAM files are annotated
                    # with information from the unaligned BAM file via Picard MergeBamAlignment.
                    runnable_step = RunnableStepPicard(
                        name='picard_merge_bam_alignment',
                        obsolete_file_path_list=[
                            file_path_read_group.merged_bam,
                        ],
                        java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx12G',
                        java_jar_path=self.java_archive_picard,
                        picard_command='MergeBamAlignment')
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    # Required Arguments
                    # OUTPUT []
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.sorted_bam)
                    # REFERENCE_SEQUENCE [null]
                    runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.genome_fasta)
                    # UNMAPPED_BAM []
                    runnable_step.add_picard_option(key='UNMAPPED_BAM', value=bam_file_path)

                    # Optional Tool Arguments
                    # ADD_MATE_CIGAR [true]
                    # ALIGNED_BAM [null]
                    runnable_step.add_picard_option(key='ALIGNED_BAM', value=file_path_read_group.merged_bam)
                    # ALIGNED_READS_ONLY [false]
                    # ALIGNER_PROPER_PAIR_FLAGS [false]
                    # ATTRIBUTES_TO_REMOVE [null]
                    # ATTRIBUTES_TO_RETAIN [null]
                    for sam_attribute in self.sam_attributes_to_retain_list:
                        runnable_step.add_picard_option(key='ATTRIBUTES_TO_RETAIN', value=sam_attribute, override=True)
                    # ATTRIBUTES_TO_REVERSE [OQ, U2]
                    # ATTRIBUTES_TO_REVERSE_COMPLEMENT [E2, SQ]
                    # CLIP_ADAPTERS [true]
                    # CLIP_OVERLAPPING_READS [true]
                    # EXPECTED_ORIENTATIONS [null]
                    # HARD_CLIP_OVERLAPPING_READS [false]
                    # INCLUDE_SECONDARY_ALIGNMENTS [true]
                    # IS_BISULFITE_SEQUENCE [false]
                    # JUMP_SIZE [null]
                    # MATCHING_DICTIONARY_TAGS [M5, LN]
                    # MAX_INSERTIONS_OR_DELETIONS [1]
                    runnable_step.add_picard_option(key='MAX_INSERTIONS_OR_DELETIONS', value='-1')
                    # MIN_UNCLIPPED_BASES [32]
                    # PAIRED_RUN [true] obsolete!
                    # PRIMARY_ALIGNMENT_STRATEGY [BestMapq]
                    # PROGRAM_GROUP_COMMAND_LINE [null]
                    # PROGRAM_GROUP_NAME [null]
                    # PROGRAM_GROUP_VERSION [null]
                    # PROGRAM_RECORD_ID [null]
                    # READ1_ALIGNED_BAM [null]
                    # READ1_TRIM [0]
                    # READ2_ALIGNED_BAM [null]
                    # READ2_TRIM [0]
                    # SORT_ORDER [coordinate]
                    runnable_step.add_picard_option(key='SORT_ORDER', value='queryname')
                    # UNMAP_CONTAMINANT_READS [false]
                    # UNMAPPED_READ_STRATEGY [DO_NOT_CHANGE]

                    # Optional Common Arguments
                    # ADD_PG_TAG_TO_READS [true]
                    # COMPRESSION_LEVEL [5]
                    # CREATE_INDEX [false]
                    # CREATE_MD5_FILE [false]
                    # MAX_RECORDS_IN_RAM [500000]
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                    # QUIET [false]
                    # TMP_DIR [null]
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_read_group.temporary_directory_path(absolute=False))
                    # USE_JDK_DEFLATER [false]
                    # USE_JDK_INFLATER [false]
                    # VALIDATION_STRINGENCY [STRICT]
                    # VERBOSITY [INFO]
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                else:
                    # If an unaligned BAM file is not available, the merged aligned BAM file is just renamed.
                    runnable_step = RunnableStepMove(
                        name='move_merged_bam',
                        source_path=file_path_read_group.merged_bam,
                        target_path=file_path_read_group.sorted_bam)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

            ###########################
            # Sample Processing Stage #
            ###########################

            # For more than one unaligned BAM file per Sample,
            # the aligned BAM files need merging into Sample-specific ones.

            file_path_sample = self.get_file_path_sample(sample_name=sample.name)

            runnable_sample = self.add_runnable(
                runnable=ConsecutiveRunnable(
                    name=self.get_prefix_sample(sample_name=sample.name),
                    working_directory=self.genome_directory,
                    cache_directory=self.cache_directory))
            executable_sample = self.set_stage_runnable(
                stage=stage_sample,
                runnable=runnable_sample)

            runnable_sample_list.append(runnable_sample)

            # Add dependencies on Runnable objects of the read group processing stage.

            for runnable_read_group in runnable_read_group_list:
                executable_sample.dependencies.append(runnable_read_group.name)

            runnable_step = RunnableStepMakeDirectory(
                name='make_directory',
                directory_path=file_path_sample.output_directory)
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            if self.skip_mark_duplicates:
                # If duplicates should *not* be marked, run Picard MergeSamFile to simply merge all
                # read group-specific aligned BAM files into one sample-specific aligned BAM file or
                # just rename a single read group-specific aligned BAM file.
                if len(runnable_read_group_list) == 1:
                    # For a single read group-specific aligned BAM file, just move (i.e., rename) the file.
                    runnable_read_group = runnable_read_group_list[0]
                    # NOTE: Since the Runnable.name already contains the prefix, FilePathReadGroup() has to be used.
                    file_path_read_group = FilePathReadGroup(prefix=runnable_read_group.name)
                    runnable_step = RunnableStepMove(
                        name='move_bam',
                        source_path=file_path_read_group.sorted_bam,
                        target_path=file_path_sample.duplicate_bam)
                    runnable_sample.add_runnable_step(runnable_step=runnable_step)
                else:
                    # For more than one read group-specific merged aligned BAM file,
                    # run Picard MergeSamFile to simply merge them into one sample-specific aligned BAM file.
                    runnable_step = RunnableStepPicard(
                        name='picard_merge_sam_files',
                        java_temporary_path=runnable_sample.temporary_directory_path(absolute=False),
                        java_heap_maximum='Xmx6G',
                        java_jar_path=self.java_archive_picard,
                        picard_command='MergeSamFiles')
                    runnable_sample.add_runnable_step(runnable_step=runnable_step)

                    # Required Arguments
                    # INPUT Required
                    for runnable_read_group in runnable_read_group_list:
                        # NOTE: Since the Runnable.name already contains the prefix, FilePathReadGroup() has to be used.
                        file_path_read_group = FilePathReadGroup(prefix=runnable_read_group.name)
                        runnable_step.add_picard_option(
                            key='INPUT',
                            value=file_path_read_group.sorted_bam,
                            override=True)
                        runnable_step.obsolete_file_path_list.append(file_path_read_group.sorted_bam)
                    # OUTPUT Required
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_sample.duplicate_bam)

                    # Optional Tool Arguments
                    # ASSUME_SORTED [false]
                    # COMMENT [null]
                    # INTERVALS [null]
                    # MERGE_SEQUENCE_DICTIONARIES [false]
                    # SORT_ORDER [coordinate]
                    # USE_THREADING [false]
                    runnable_step.add_picard_option(key='USE_THREADING', value='true')

                    # Optional Common Arguments
                    # COMPRESSION_LEVEL [5]
                    # CREATE_INDEX [false]
                    # CREATE_MD5_FILE [false]
                    # MAX_RECORDS_IN_RAM [500000]
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                    # QUIET [false]
                    # REFERENCE_SEQUENCE [null]
                    # TMP_DIR [null]
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_sample.temporary_directory_path(absolute=False))
                    # USE_JDK_DEFLATER [false]
                    # USE_JDK_INFLATER [false]
                    # VALIDATION_STRINGENCY [STRICT]
                    # VERBOSITY [INFO]
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
            else:
                # If duplicates should be marked,
                # run Picard MarkDuplicates on all read group-specific merged aligned BAM files at once.
                runnable_step = RunnableStepPicard(
                    name='picard_mark_duplicates',
                    java_temporary_path=runnable_sample.temporary_directory_path(absolute=False),
                    java_heap_maximum='Xmx6G',
                    java_jar_path=self.java_archive_picard,
                    picard_command='MarkDuplicates')
                runnable_sample.add_runnable_step(runnable_step=runnable_step)

                # Required Arguments
                # INPUT Required
                for runnable_read_group in runnable_read_group_list:
                    # NOTE: Since the Runnable.name already contains the prefix, FilePathReadGroup() has to be used.
                    file_path_read_group = FilePathReadGroup(prefix=runnable_read_group.name)
                    runnable_step.add_picard_option(key='INPUT', value=file_path_read_group.sorted_bam, override=True)
                    runnable_step.obsolete_file_path_list.append(file_path_read_group.sorted_bam)
                # METRICS_FILE Required
                runnable_step.add_picard_option(key='METRICS_FILE', value=file_path_sample.duplicate_metrics_tsv)
                # OUTPUT Required
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_sample.duplicate_bam)

                # Optional Tool Arguments
                # ASSUME_SORT_ORDER [null]
                # ASSUME_SORTED Deprecated
                # BARCODE_TAG [null]
                if sample_umis:
                    runnable_step.add_picard_option(key='BARCODE_TAG', value='RX')
                # CLEAR_DT [true]
                # COMMENT [null]
                # DUPLEX_UMI [false]
                # DUPLICATE_SCORING_STRATEGY [SUM_OF_BASE_QUALITIES]
                # FLOW_EFFECTIVE_QUALITY_THRESHOLD [15]
                # FLOW_MODE [false]
                # FLOW_Q_IS_KNOWN_END [false]
                # FLOW_QUALITY_SUM_STRATEGY [false]
                # FLOW_SKIP_FIRST_N_FLOWS [0]
                # MAX_FILE_HANDLES_FOR_READ_ENDS_MAP [8000]
                # MAX_OPTICAL_DUPLICATE_SET_SIZE [300000]
                # MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP Obsolete [50000]
                # MOLECULAR_IDENTIFIER_TAG [null]
                # OPTICAL_DUPLICATE_PIXEL_DISTANCE [100]
                # NOTE: Should be 2500 for patterned flow cells.
                runnable_step.add_picard_option(key='OPTICAL_DUPLICATE_PIXEL_DISTANCE', value='3000')
                # PROGRAM_GROUP_COMMAND_LINE [null]
                # PROGRAM_GROUP_NAME [MarkDuplicates]
                # PROGRAM_GROUP_VERSION [null]
                # PROGRAM_RECORD_ID [MarkDuplicates]
                # READ_NAME_REGEX [.]
                # READ_ONE_BARCODE_TAG [null]
                # READ_TWO_BARCODE_TAG [null]
                # REMOVE_DUPLICATES [false]
                # REMOVE_SEQUENCING_DUPLICATES [false]
                # SORTING_COLLECTION_SIZE_RATIO [0.25]
                # TAG_DUPLICATE_SET_MEMBERS [false]
                # TAGGING_POLICY [DontTag]
                # UNPAIRED_END_UNCERTAINTY [0]
                # USE_END_IN_UNPAIRED_READS [false]
                # USE_UNPAIRED_CLIPPED_END [false]

                # Optional Common Arguments
                # ADD_PG_TAG_TO_READS [true]
                # COMPRESSION_LEVEL [5]
                # CREATE_INDEX [false]
                # CREATE_MD5_FILE [false]
                # MAX_RECORDS_IN_RAM [500000]
                runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
                # QUIET [false]
                # REFERENCE_SEQUENCE [null]
                # TMP_DIR [null]
                runnable_step.add_picard_option(
                    key='TMP_DIR',
                    value=runnable_sample.temporary_directory_path(absolute=False))
                # USE_JDK_DEFLATER [false]
                # USE_JDK_INFLATER [false]
                # VALIDATION_STRINGENCY [STRICT]
                # VERBOSITY [INFO]
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

            # Finally, run Picard SortSam by coordinate.

            runnable_step = RunnableStepPicard(
                name='picard_sort_sam',
                obsolete_file_path_list=[
                    file_path_sample.duplicate_bam,
                ],
                java_temporary_path=runnable_sample.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx6G',
                java_jar_path=self.java_archive_picard,
                picard_command='SortSam')
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            # Required Arguments
            # INPUT []
            runnable_step.add_picard_option(key='INPUT', value=file_path_sample.duplicate_bam)
            # OUTPUT []
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_sample.sample_bam)
            # SORT_ORDER []
            runnable_step.add_picard_option(key='SORT_ORDER', value='coordinate')

            # Optional Common Arguments
            # COMPRESSION_LEVEL [5]
            runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='9')
            # CREATE_INDEX [false]
            runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
            # CREATE_MD5_FILE [false]
            runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')
            # MAX_RECORDS_IN_RAM [500000]
            runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='2000000')
            # QUIET [false]
            # REFERENCE_SEQUENCE [null]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_sample.temporary_directory_path(absolute=False))
            # USE_JDK_INFLATER [false]
            # USE_JDK_DEFLATER [false]
            # VALIDATION_STRINGENCY [STRICT]
            # VERBOSITY [INFO]
            runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

            # Create a symbolic link from the Picard-style <file_name>.bai file to a
            # samtools-style <file_name>.bam.bai file.

            runnable_step = RunnableStepLink(
                name='link',
                source_path=file_path_sample.sample_bai_link_source,
                target_path=file_path_sample.sample_bai_link_target)
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            # Run the Picard CollectAlignmentSummaryMetrics analysis.

            runnable_step = RunnableStepPicard(
                name='picard_collect_alignment_summary_metrics',
                java_temporary_path=runnable_sample.temporary_directory_path(absolute=False),
                java_heap_maximum='Xmx6G',
                java_jar_path=self.java_archive_picard,
                picard_command='CollectAlignmentSummaryMetrics')
            runnable_sample.add_runnable_step(runnable_step=runnable_step)

            # Required Arguments
            # INPUT []
            runnable_step.add_picard_option(key='INPUT', value=file_path_sample.sample_bam)
            # OUTPUT []
            runnable_step.add_picard_option(key='OUTPUT', value=file_path_sample.alignment_summary_metrics_tsv)

            # Optional Tool Arguments
            # ADAPTER_SEQUENCE [
            #   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT,
            #   AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG,
            #   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT,
            #   AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG,
            #   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT,
            #   AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
            # ]
            # ASSUME_SORTED [true]
            # COLLECT_ALIGNMENT_INFORMATION [true]
            # EXPECTED_PAIR_ORIENTATIONS [FR]
            # HISTOGRAM_FILE []
            # IS_BISULFITE_SEQUENCED [false]
            # MAX_INSERT_SIZE [100000]
            # METRIC_ACCUMULATION_LEVEL [ALL_READS]
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='SAMPLE', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='LIBRARY', override=True)
            runnable_step.add_picard_option(key='METRIC_ACCUMULATION_LEVEL', value='READ_GROUP', override=True)
            # REFERENCE_SEQUENCE [null]
            runnable_step.add_picard_option(key='REFERENCE_SEQUENCE', value=self.genome_fasta)
            # STOP_AFTER [0]

            # Optional Common Arguments
            # COMPRESSION_LEVEL [5]
            # CREATE_INDEX [false]
            # CREATE_MD5_FILE [false]
            # MAX_RECORDS_IN_RAM [500000]
            # QUIET [false]
            # TMP_DIR [null]
            runnable_step.add_picard_option(
                key='TMP_DIR',
                value=runnable_sample.temporary_directory_path(absolute=False))
            # USE_JDK_DEFLATER [false]
            # USE_JDK_INFLATER [false]
            # VALIDATION_STRINGENCY [STRICT]
            # VERBOSITY [INFO]
            runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')

            # Add aligner-specific RunnableStep objects.

            self.add_runnable_step_sample(runnable_sample=runnable_sample, stage_sample=stage_sample)

        #################
        # Summary Stage #
        #################

        # Write the AnnotationSheet to disk.

        annotation_sheet.to_file_path()

        # Create a Runnable and Executable for the summary.

        runnable_summary = self.add_runnable(
            runnable=ConsecutiveRunnable(
                name=self.get_prefix_summary(),
                working_directory=self.genome_directory,
                cache_directory=self.cache_directory))
        executable_summary = self.set_stage_runnable(
            stage=stage_summary,
            runnable=runnable_summary)

        # Add dependencies on Runnable objects of the sample stage.
        for runnable_sample in runnable_sample_list:
            executable_summary.dependencies.append(runnable_sample.name)

        # Add a RunnableStep to summarise and plot the Picard Alignment Summary Metrics reports.

        runnable_step = RunnableStep(
            name='picard_alignment_summary_metrics',
            program='picard_reports_alignment_summary.R')
        runnable_summary.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_option_long(key='prefix', value=f'{self.get_stage_name_summary()}_pasm')
        runnable_step.add_option_long(
            key='file-pattern',
            value=f'^{self.get_stage_name_sample()}_(.*)_alignment_summary_metrics.tsv$')

        # Add a RunnableStep to summarise and plot the Picard Duplication Metrics reports.

        runnable_step = RunnableStep(
            name='picard_duplication_metrics',
            program='picard_reports_duplication.R')
        runnable_summary.add_runnable_step(runnable_step=runnable_step)

        runnable_step.add_option_long(key='prefix', value=f'{self.get_stage_name_summary()}_pdm')
        runnable_step.add_option_long(
            key='file-pattern',
            value=f'^{self.get_stage_name_sample()}_(.*)_duplicate_metrics.tsv$')

        # Add aligner-specific RunnableStep objects.

        self.add_runnable_step_summary(runnable_summary=runnable_summary, stage_summary=stage_summary)

        return

    def report_html_sample(self) -> list[str]:
        """Create a :literal:`XHTML 1.0` sample table.

        :return: A Python :py:class:`list` object of Python :py:class:`str` objects.
        :rtype: list[str]
        """
        str_list: list[str] = list()

        str_list.append('<h2 id="sample_section">Sample Table</h2>\n')
        str_list.append('\n')
        str_list.append('<table id="sample_table">\n')
        str_list.append('<thead>\n')
        str_list.append('<tr>\n')
        str_list.append('<th>Sample</th>\n')
        str_list.append('<th>BAM</th>\n')
        str_list.append('<th>BAI</th>\n')
        str_list.append('<th>MD5</th>\n')
        str_list.append('</tr>\n')
        str_list.append('</thead>\n')
        str_list.append('<tbody>\n')

        for sample in self.sample_list:
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            file_path_sample = self.get_file_path_sample(sample_name=sample.name)

            str_list.append('<tr>\n')
            # Sample
            str_list.append('<td class="left">' + sample.name + '</td>\n')
            # BAM
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_sample.sample_bam + '">')
            str_list.append('<abbr title="Binary Alignment/Map">BAM</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            # BAI
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_sample.sample_bai + '">')
            str_list.append('<abbr title="Binary Alignment/Map Index">BAI</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            # MD5
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_sample.sample_md5 + '">')
            str_list.append('<abbr title="Message Digest 5 Checksum">MD5</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('</tr>\n')

        str_list.append('</tbody>\n')
        str_list.append('</table>\n')
        str_list.append('\n')

        return str_list

    def report_hub_alignment(self) -> list[str]:
        """Create an :literal:`Alignment` UCSC Composite Track.

        :return: A Python :py:class:`list` object of Python :py:class:`str` objects.
        :rtype: list[str]
        """
        str_list: list[str] = list()

        # Common track settings
        str_list.append('track alignment\n')
        str_list.append('type bam\n')
        str_list.append('shortLabel Alignment\n')
        str_list.append('longLabel ' + self.name + ' Alignment\n')
        # str_list.append('html ...\n')
        str_list.append('visibility hide\n')
        # Composite track settings
        str_list.append('compositeTrack on\n')
        str_list.append('allButtonPair on\n')  # Has to be "off" to allow for configuration via a matrix.
        str_list.append('centerLabelsDense on\n')
        # str_list.append('dragAndDrop subTracks\n')
        str_list.append('\n')

        # Sample-specific tracks

        for sample in self.sample_list:
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            file_path_sample = self.get_file_path_sample(sample_name=sample.name)

            # Common track settings
            str_list.append('  track ' + sample.name + '_alignment\n')
            str_list.append('  type bam\n')
            str_list.append('  shortLabel ' + '_'.join((sample.name, self.prefix, 'alignment')) + '\n')
            str_list.append('  longLabel ' + ' '.join((sample.name, self.name, 'Alignment')) + '\n')
            str_list.append('  bigDataUrl ' + file_path_sample.sample_bam + '\n')
            # str_list.append('  html ...\n')
            str_list.append('  visibility dense\n')
            # Common optional track settings
            # str_list.append('  color 0,0,0\n')
            # bam/cram - Compressed Sequence Alignment track settings
            # Composite track settings
            str_list.append('  parent alignment on\n')
            str_list.append('  centerLabelsDense on\n')
            # str_list.append('  dragAndDrop subTracks\n')
            str_list.append('  \n')

        return str_list

    def report_hub_coverage(self) -> list[str]:
        """Create a :literal:`Coverage` UCSC Composite Track.

        :return: A Python :py:class:`list` of Python :py:class:`str` objects.
        :rtype: list[str]
        """
        str_list: list[str] = list()

        # Common track settings
        str_list.append('track coverage\n')
        str_list.append('type bigWig\n')
        str_list.append('shortLabel Coverage\n')
        str_list.append('longLabel ' + self.name + ' Coverage\n')
        # str_list.append('html ...\n')
        str_list.append('visibility hide\n')
        # Composite track settings
        str_list.append('compositeTrack on\n')
        str_list.append('allButtonPair on\n')  # Has to be "off" to allow for configuration via a matrix.
        str_list.append('centerLabelsDense on\n')
        # str_list.append('dragAndDrop subTracks\n')
        str_list.append('\n')

        # Sample-specific tracks

        for sample in self.sample_list:
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            file_path_sample = self.get_file_path_sample(sample_name=sample.name)

            # Common track settings
            str_list.append('  track ' + sample.name + '_coverage\n')
            str_list.append('  ' + self.ucsc_hub_bigwig_info_signal_range(file_path=file_path_sample.sample_bwi))
            str_list.append('  shortLabel ' + sample.name + '_coverage\n')
            str_list.append('  longLabel ' + sample.name + ' STAR RNA-seq alignment coverage\n')
            str_list.append('  bigDataUrl ' + file_path_sample.sample_bw + '\n')
            # str_list.append('  html ...\n')
            str_list.append('  visibility full\n')
            # Common optional track settings
            str_list.append('  color 0,0,0\n')
            # bigWig - Signal graphing track settings
            str_list.append('  alwaysZero on\n')
            str_list.append('  autoScale on\n')
            str_list.append('  graphTypeDefault bar\n')
            str_list.append('  maxHeightPixels 100:60:20\n')
            # str_list.append('  maxWindowToQuery 10000000\n')
            # str_list.append('  smoothingWindow 5\n')
            # str_list.append('  transformFunc NONE\n')
            # str_list.append('  viewLimits 0:45\n')
            # str_list.append('  viewLimitsMax 0:50\n')
            # str_list.append('  windowingFunction maximum\n')
            # str_list.append('  yLineMark <#>\n')
            # str_list.append('  yLineOnOff on \n')
            # str_list.append('  gridDefault on\n')
            # Composite track settings
            str_list.append('  parent coverage on\n')
            str_list.append('  centerLabelsDense on\n')
            # str_list.append('  dragAndDrop subTracks\n')
            str_list.append('  \n')

        return str_list

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

        # Alignment Summary Plot Absolute Mapped
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_numbers_sample_pdf + '">')
        str_list.append('<img alt="Absolute Mapped - Sample"')
        str_list.append(' src="' + file_path_summary.pasm_numbers_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_numbers_read_group_pdf + '">')
        str_list.append('<img alt="Absolute Mapped - Read Group"')
        str_list.append(' src="' + file_path_summary.pasm_numbers_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Absolute Mapped</td>\n')
        str_list.append('</tr>\n')

        # Alignment Summary Plot Percentage Mapped
        str_list.append('<tr>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_fractions_sample_pdf + '">')
        str_list.append('<img alt="Percentage Mapped - Sample"')
        str_list.append(' src="' + file_path_summary.pasm_fractions_sample_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="center">')
        str_list.append('<a href="' + file_path_summary.pasm_fractions_read_group_pdf + '">')
        str_list.append('<img alt="Percentage Mapped - Read Group"')
        str_list.append(' src="' + file_path_summary.pasm_fractions_read_group_png + '"')
        str_list.append(' height="100" width="100" />')
        str_list.append('</a>')
        str_list.append('</td>\n')
        str_list.append('<td class="left">Percentage Mapped</td>\n')
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

        if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_levels_png)):
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pdm_levels_pdf + '">')
            str_list.append('<img alt="Duplication Levels - Sample"')
            str_list.append(' src="' + file_path_summary.pdm_levels_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Duplication Levels</td>\n')
            str_list.append('</tr>\n')

        if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdm_fractions_png)):
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pdm_fractions_pdf + '">')
            str_list.append('<img alt="Duplication Percentage - Sample"')
            str_list.append(' src="' + file_path_summary.pdm_fractions_png + '"')
            str_list.append(' height="100" width="100" />')
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

        # Save a UCSC Track Hub with an alignment composite track.
        self.ucsc_hub_to_file(content=self.report_hub_alignment())

        return
