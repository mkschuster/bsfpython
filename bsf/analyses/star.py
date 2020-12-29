# -*- coding: utf-8 -*-
"""STAR Analysis module.

A package of classes and methods supporting the Spliced Transcripts Alignment
to a Reference (STAR) by Alexander Dobin.

Project:  https://github.com/alexdobin/STAR
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

from bsf.analyses.aligner import Aligner, \
    FilePathAlign as AlignerFilePathAlign, \
    FilePathSummary as AlignerFilePathSummary
from bsf.analysis import Stage
from bsf.connector import ConnectorFile
from bsf.ngs import Collection, Sample
from bsf.process import RunnableStep
from bsf.standards import Configuration, StandardFilePath, Transcriptome


class FilePathAlign(AlignerFilePathAlign):
    """The C{bsf.analyses.star.FilePathAlign} class models file paths at the alignment stage.

    Attributes:
    @ivar aligned_sam: Aligned sequence alignment map (SAM) file path
    @type aligned_sam: str
    @ivar splice_junctions_tsv: Splice junctions tab-separated value (TSV) file path
    @type splice_junctions_tsv: str
    @ivar star_prefix: STAR outFileNamePrefix file path
    @type star_prefix: str
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.star.FilePathAlign} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        # Override the aligned_sam instance variable of the super-class.
        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'Aligned.out.sam')))

        self.splice_junctions_tsv = os.path.join(prefix, '_'.join((prefix, 'SJ.out.tab')))

        self.star_prefix = os.path.join(prefix, '_'.join((prefix, '')))

        return


class FilePathSummary(AlignerFilePathSummary):
    """The C{bsf.analyses.star.FilePathSummary} class models file paths at the summary stage.

    Attributes:
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.star.FilePathSummary} object

        @param prefix: Prefix
        @type prefix: str
        """
        super(FilePathSummary, self).__init__(prefix=prefix)

        self.alignment_read_group_pdf = prefix + '_alignment_read_group.pdf'
        self.alignment_read_group_png = prefix + '_alignment_read_group.png'
        self.alignment_sample_png = prefix + '_alignment_sample.png'
        self.alignment_sample_pdf = prefix + '_alignment_sample.pdf'
        self.junction_fraction_read_group_pdf = prefix + '_junction_fraction_read_group.pdf'
        self.junction_fraction_read_group_png = prefix + '_junction_fraction_read_group.png'
        self.junction_fraction_sample_pdf = prefix + '_junction_fraction_sample.pdf'
        self.junction_fraction_sample_png = prefix + '_junction_fraction_sample.png'
        self.junction_number_read_group_pdf = prefix + '_junction_number_read_group.pdf'
        self.junction_number_read_group_png = prefix + '_junction_number_read_group.png'
        self.junction_number_sample_pdf = prefix + '_junction_number_sample.pdf'
        self.junction_number_sample_png = prefix + '_junction_number_sample.png'
        self.mapped_fraction_read_group_pdf = prefix + '_mapped_fraction_read_group.pdf'
        self.mapped_fraction_read_group_png = prefix + '_mapped_fraction_read_group.png'
        self.mapped_fraction_sample_png = prefix + '_mapped_fraction_sample.png'
        self.mapped_fraction_sample_pdf = prefix + '_mapped_fraction_sample.pdf'
        self.mapped_number_read_group_pdf = prefix + '_mapped_number_read_group.pdf'
        self.mapped_number_read_group_png = prefix + '_mapped_number_read_group.png'
        self.mapped_number_sample_png = prefix + '_mapped_number_sample.png'
        self.mapped_number_sample_pdf = prefix + '_mapped_number_sample.pdf'
        self.table_read_group_tsv = prefix + '_table_read_group.tsv'
        self.table_sample_tsv = prefix + '_table_sample.tsv'

        return


class Star(Aligner):
    """STAR C{bsf.analyses.aligner.Aligner} sub-class.

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar sam_attributes_to_retain_list: A Python C{list} of aligner-specific, private SAM tags (i.e. X*, Y*, z*)
        that should be retained by Picard MergeBamAlignment
    @type sam_attributes_to_retain_list: list[str]
    @ivar transcriptome_index: Transcriptome index directory path
    @type transcriptome_index: str | None
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str | None
    @ivar transcriptome_gtf: Transcriptome annotation GTF file path
    @type transcriptome_gtf: str | None
    @ivar two_pass_mapping: Basic two-pass mapping
    @type two_pass_mapping: str | None
    @ivar java_archive_picard: Picard tools Java Archive (JAR) file path
    @type java_archive_picard: str | None
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
    def get_file_path_align(cls, paired_reads_name):
        """Get a C{FilePathAlign} object from this or a sub-class.

        @param paired_reads_name: C{bsf.ngs.PairedReads.name}
        @type paired_reads_name: str
        @return: C{FilePathAlign} or sub-class object
        @rtype: FilePathAlign
        """
        return FilePathAlign(prefix=cls.get_prefix_align(paired_reads_name=paired_reads_name))

    @classmethod
    def get_file_path_summary(cls):
        """Get a C{FilePathSummary} object from this or a sub-class.

        @return: C{FilePathSummary} or sub-class object
        @rtype: FilePathSummary
        """
        return FilePathSummary(prefix=cls.get_prefix_summary())

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
            transcriptome_index=None,
            transcriptome_version=None,
            transcriptome_gtf=None,
            two_pass_mapping=None,
            skip_mark_duplicates=None,
            java_archive_picard=None):
        """Initialise a C{bsf.analyses.star.Star} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.analysis.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.analysis.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.analysis.Analysis}-wide project directory,
            normally under the C{bsf.analysis.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.analysis.Analysis}-wide genome directory,
            normally under the C{bsf.analysis.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.analysis.Stage} objects
        @type stage_list: list[Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[Sample]
        @param transcriptome_index: Transcriptome index directory path
        @type transcriptome_index: str
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @param transcriptome_gtf: Transcriptome annotation GTF file path
        @type transcriptome_gtf: str | None
        @param two_pass_mapping: Basic two-pass mapping
        @type two_pass_mapping: str | None
        @param skip_mark_duplicates: Mark duplicates
        @type skip_mark_duplicates: bool | None
        @param java_archive_picard: Picard tools Java Archive (JAR) file path
        @type java_archive_picard: str | None
        """
        super(Star, self).__init__(
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
            genome_fasta=None,
            genome_index=None,
            skip_mark_duplicates=skip_mark_duplicates,
            java_archive_picard=java_archive_picard)

        # Sub-class specific ...

        self.transcriptome_index = transcriptome_index
        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf = transcriptome_gtf
        self.two_pass_mapping = two_pass_mapping

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.star.Star} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """
        super(Star, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf = configuration.config_parser.get(section=section, option=option)

        option = 'two_pass_mapping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.two_pass_mapping = configuration.config_parser.get(section=section, option=option)

        return

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add one or more STAR-specific C{bsf.process.RunnableStep} objects to the C{bsf.procedure.ConcurrentRunnable}.

        @param runnable_align: C{bsf.procedure.ConcurrentRunnable}
        @type runnable_align: bsf.procedure.ConcurrentRunnable
        @param stage_align: C{bsf.analysis.Stage}
        @type stage_align: Stage
        @param file_path_1: FASTQ file path 1
        @type file_path_1: str | None
        @param file_path_2: FASTQ file path 2
        @type file_path_2: str | None
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
        runnable_align.add_runnable_step_post(runnable_step=runnable_step)

        runnable_step.add_switch_long(key='best')
        runnable_step.arguments.append(file_path_align.splice_junctions_tsv)

        return

    def add_runnable_step_summary(self, runnable_summary, stage_summary):
        """Add one or more STAR-specific C{bsf.process.RunnableStep} objects to the C{bsf.procedure.Runnable}.

        @param runnable_summary: C{bsf.procedure.ConsecutiveRunnable}
        @type runnable_summary: bsf.procedure.ConsecutiveRunnable
        @param stage_summary: C{bsf.analysis.Stage}
        @type stage_summary: Stage
        """
        runnable_step = RunnableStep(
            name='star_summary',
            program='bsf_star_summary.R')
        runnable_summary.add_runnable_step(runnable_step=runnable_step)

        return

    def run(self):
        """Run this C{bsf.analyses.star.Star} analysis.

        Although STAR can directly count reads according to its splice junction database,
        more than one read group may need aligning so that the count tables had to be combined.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # STAR requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A ' + self.name + " requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a valid 'transcriptome_version' configuration option.")

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
            raise ('The ' + self.name + " 'two_pass_mapping' option can only take values " +
                   repr(two_pass_mapping_tuple) + ' not ' + repr(self.two_pass_mapping) + '.')

        super(Star, self).run()

        return

    def report(self):
        """Create a report.
        """

        def report_html():
            """Private function to create a HTML report.
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str] """

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
            str_list.append('<a href="' + file_path_summary.pasm_alignment_sample_pdf + '">')
            str_list.append('<img alt="Alignment Summary - Sample"')
            str_list.append(' src="' + file_path_summary.pasm_alignment_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pasm_alignment_read_group_pdf + '">')
            str_list.append('<img alt="Alignment Summary - Read Group"')
            str_list.append(' src="' + file_path_summary.pasm_alignment_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Alignment Summary</td>\n')
            str_list.append('</tr>\n')

            # Alignment Summary Plot Absolute Mapped
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pasm_absolute_sample_pdf + '">')
            str_list.append('<img alt="Absolute Mapped - Sample"')
            str_list.append(' src="' + file_path_summary.pasm_absolute_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pasm_absolute_read_group_pdf + '">')
            str_list.append('<img alt="Absolute Mapped - Read Group"')
            str_list.append(' src="' + file_path_summary.pasm_absolute_read_group_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Absolute Mapped</td>\n')
            str_list.append('</tr>\n')

            # Alignment Summary Plot Percentage Mapped
            str_list.append('<tr>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pasm_percentage_sample_pdf + '">')
            str_list.append('<img alt="Percentage Mapped - Sample"')
            str_list.append(' src="' + file_path_summary.pasm_percentage_sample_png + '"')
            str_list.append(' height="100" width="100" />')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.pasm_percentage_read_group_pdf + '">')
            str_list.append('<img alt="Percentage Mapped - Read Group"')
            str_list.append(' src="' + file_path_summary.pasm_percentage_read_group_png + '"')
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

            if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdsm_levels_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.pdsm_levels_sample_pdf + '">')
                str_list.append('<img alt="Duplication Levels - Sample"')
                str_list.append(' src="' + file_path_summary.pdsm_levels_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Duplication Levels</td>\n')
                str_list.append('</tr>\n')

            if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdsm_percentage_sample_png)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.pdsm_percentage_sample_pdf + '">')
                str_list.append('<img alt="Duplication Percentage - Sample"')
                str_list.append(' src="' + file_path_summary.pdsm_percentage_sample_png + '"')
                str_list.append(' height="100" width="100" />')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('<td class="center">')
                str_list.append('</td>\n')
                str_list.append('<td class="left">Duplication Percentage</td>\n')
                str_list.append('</tr>\n')

            if os.path.exists(os.path.join(self.genome_directory, file_path_summary.pdsm_table_sample_tsv)):
                str_list.append('<tr>\n')
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_summary.pdsm_table_sample_tsv + '">')
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

            self.report_to_file(content=str_list)

            return

        def report_hub():
            """Private function to create a UCSC Track Hub.
            """

            str_list = list()
            """ @type str_list: list[str] """

            # Group via UCSC super tracks.

            str_list.append('track alignment\n')
            str_list.append('type bam\n')
            str_list.append('shortLabel Alignment\n')
            str_list.append('longLabel ' + self.name + ' Alignment\n')
            str_list.append('visibility hide\n')
            str_list.append('compositeTrack on\n')
            str_list.append('allButtonPair on\n')  # Has to be off to allow for configuration via a matrix.
            str_list.append('centerLabelsDense on\n')
            str_list.append('\n')

            # Sample-specific tracks

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                file_path_sample = self.get_file_path_sample(sample_name=sample.name)

                #
                # Add a trackDB entry for each Tophat accepted_hits.bam file.
                #
                # Common settings
                str_list.append('  track ' + sample.name + '_alignment\n')
                str_list.append('  type bam\n')
                str_list.append('  shortLabel ' + '_'.join((sample.name, self.prefix, 'alignment')) + '\n')
                str_list.append('  longLabel ' + ' '.join((sample.name, self.name, 'Alignment')) + '\n')
                str_list.append('  bigDataUrl ' + file_path_sample.sample_bam + '\n')
                # str_list.append('  html ...\n')
                str_list.append('  visibility dense\n')

                # Common optional settings
                # str_list.append('  color 0,0,0\n')

                # bam/cram - Compressed Sequence Alignment track settings
                # None

                # Composite track settings
                str_list.append('  parent alignment on\n')
                str_list.append('  centerLabelsDense on\n')
                str_list.append('  \n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return
