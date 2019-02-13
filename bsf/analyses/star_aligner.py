# -*- coding: utf-8 -*-
"""STAR Analysis module

A package of classes and methods supporting the spliced Transcripts Alignment
to a Reference (STAR) aligner by Alexander Dobin.

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

from __future__ import print_function

import os

import bsf
import bsf.analyses.aligner
import bsf.annotation
import bsf.process
import bsf.standards


class FilePathAlign(bsf.analyses.aligner.FilePathAlign):
    """The C{bsf.analyses.star_aligner.FilePathAlign} class models file paths at the alignment stage.

    Attributes:
    @ivar aligned_sam: Aligned sequence alignment map (SAM) file path
    @type aligned_sam: str | unicode
    @ivar splice_junctions_tsv: Splice junctions tab-separated value (TSV) file path
    @type splice_junctions_tsv: str | unicode
    @ivar star_prefix: STAR Aligner outFileNamePrefix file path
    @type star_prefix: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.star_aligner.FilePathAlign} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        # Override the aligned_sam instance variable of the super-class.
        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'Aligned.out.sam')))

        self.splice_junctions_tsv = os.path.join(prefix, '_'.join((prefix, 'SJ.out.tab')))

        self.star_prefix = os.path.join(prefix, '_'.join((prefix, '')))

        return


class FilePathSummary(bsf.analyses.aligner.FilePathSummary):
    """The C{bsf.analyses.star_aligner.FilePathSummary} class models file paths at the summary stage.

    Attributes:
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.star_aligner.FilePathSummary} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype:
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
        self.table_read_group = prefix + '_table_read_group.tsv'
        self.table_sample = prefix + '_table_sample.tsv'

        return


class StarAligner(bsf.analyses.aligner.Aligner):
    """STAR Aligner C{bsf.analyses.aligner.Aligner} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar index_directory: Genome directory with STAR indices
    @type index_directory: str | unicode | None
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str | None
    @ivar transcriptome_gtf: GTF file path of transcriptome annotation
    @type transcriptome_gtf: str | unicode | None
    @ivar two_pass_mapping: Basic two-pass mapping
    @type two_pass_mapping: bool | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode | None
    """

    name = 'STAR Aligner Analysis'
    prefix = 'star_aligner'

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
            index_directory=None,
            transcriptome_version=None,
            transcriptome_gtf=None,
            two_pass_mapping=None,
            classpath_picard=None):
        """Initialise a C{bsf.analyses.star_aligner.StarAligner} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @param index_directory: Genome directory with STAR indices
        @type index_directory: str | unicode
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @param transcriptome_gtf: GTF file path of transcriptome annotation
        @type transcriptome_gtf: str | unicode | None
        @param two_pass_mapping: Basic two-pass mapping
        @type two_pass_mapping: str | unicode | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode | None
        """
        super(StarAligner, self).__init__(
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
            classpath_picard=classpath_picard)

        # Sub-class specific ...

        self.index_directory = index_directory
        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf = transcriptome_gtf
        self.two_pass_mapping = two_pass_mapping

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.star_aligner.StarAligner} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(StarAligner, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'index_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.index_directory = configuration.config_parser.get(section=section, option=option)

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
        """Add one or more STAR-specific C{bsf.process.RunnableStep} objects to the C{bsf.Runnable}.

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
        """ @type file_path_align FilePathAlign """

        # Run the STAR Aligner

        runnable_step = runnable_align.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='STAR',
                program='STAR'))
        """ @type runnable_step: bsf.process.RunnableStep """
        self.set_runnable_step_configuration(runnable_step=runnable_step)
        runnable_step.add_option_long(key='runThreadN', value=str(stage_align.threads))
        runnable_step.add_option_long(key='genomeDir', value=self.index_directory)
        runnable_step.add_option_long(key='outFileNamePrefix', value=file_path_align.star_prefix)
        if self.two_pass_mapping == 'basic':
            runnable_step.add_option_long(key='twopassMode', value='Basic')
        # NOTE: The STAR aligner command line interface is seriously broken,
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

        runnable_step = runnable_align.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='gzip',
                program='gzip'))
        """ @type runnable_step: bsf.process.RunnableStep """
        runnable_step.add_switch_long(key='best')
        runnable_step.arguments.append(file_path_align.splice_junctions_tsv)

        return

    def add_runnable_step_summary(self, runnable_summary, stage_summary):
        """Add one or more STAR-specific C{bsf.process.RunnableStep} objects to the C{bsf.Runnable}.

        @param runnable_summary: C{bsf.Runnable}
        @type runnable_summary: bsf.Runnable
        @param stage_summary: C{bsf.Stage}
        @type stage_summary: bsf.Stage
        @return:
        @rtype:
        """
        runnable_summary.add_runnable_step(
            runnable_step=bsf.process.RunnableStep(
                name='summary',
                program='bsf_star_aligner_summary.R'))
        """ @type runnable_step: bsf.process.RunnableStep """

        return

    def run(self):
        """Run this C{bsf.analyses.star_aligner.StarAligner} analysis.

        Although the STAR aligner can directly count reads according to its splice junction database,
        more than one read group may need aligning so that the count tables had to be combined.
        @return:
        @rtype:
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # The STAR Aligner requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A ' + self.name + " requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.Analysis super-class.

        if not self.genome_version:
            self.genome_version = bsf.standards.Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.index_directory:
            self.index_directory = bsf.standards.FilePath.get_resource_transcriptome_index(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='star')

        if not self.transcriptome_gtf:
            # FIXME: The transcriptome_gtf is currently not used.
            self.transcriptome_gtf = bsf.standards.FilePath.get_resource_transcriptome_gtf(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='star')

        if not self.two_pass_mapping:
            self.two_pass_mapping = 'none'

        two_pass_mapping_tuple = ('none', 'basic', 'full')
        if self.two_pass_mapping not in two_pass_mapping_tuple:
            raise ('The ' + self.name + " 'two_pass_mapping' option can only take values " +
                   repr(two_pass_mapping_tuple) + ' not ' + repr(self.two_pass_mapping) + '.')

        super(StarAligner, self).run()

        return

    def report(self):
        """Create a report.

        @return:
        @rtype:
        """

        def report_html():
            """Private function to create a HTML report.

            @return:
            @rtype:
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # This code only needs the public URL.

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str | unicode] """

            str_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.extend(self.get_html_genome(genome_version=self.genome_version))
            str_list.extend(self.get_html_transcriptome(transcriptome_version=self.transcriptome_version))
            str_list.append('\n')

            str_list.append('<p id="ucsc_track_hub">')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('</p>\n')
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

                runnable_sample = self.runnable_dict[self.get_prefix_sample(sample_name=sample.name)]
                file_path_sample = runnable_sample.file_path_object
                """ @type file_path_sample: bsf.analyses.aligner.FilePathSample """

                str_list.append('<tr>\n')
                # Sample
                str_list.append('<td class="left">' + sample.name + '</td>\n')
                # BAM
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_sample.merged_bam + '">')
                str_list.append('<abbr title="Binary Alignment/Map">BAM</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # BAI
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_sample.merged_bai + '">')
                str_list.append('<abbr title="Binary Alignment/Map Index">BAI</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                # MD5
                str_list.append('<td class="center">')
                str_list.append('<a href="' + file_path_sample.merged_md5 + '">')
                str_list.append('<abbr title="Message Digest 5 Checksum">MD5</abbr>')
                str_list.append('</a>')
                str_list.append('</td>\n')
                str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
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

            runnable_summary = self.runnable_dict[self.get_prefix_summary()]
            file_path_summary = runnable_summary.file_path_object
            """ @type file_path_summary: FilePathSummary """

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
            str_list.append('<a href="' + file_path_summary.table_sample + '">')
            str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="center">')
            str_list.append('<a href="' + file_path_summary.table_read_group + '">')
            str_list.append('<abbr title="Tab-Separated Value">TSV</abbr>')
            str_list.append('</a>')
            str_list.append('</td>\n')
            str_list.append('<td class="left">Summary</td>\n')
            str_list.append('</tr>\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            self.report_to_file(content=str_list)

            return

        def report_hub():
            """Private function to create a UCSC Track Hub.

            @return:
            @rtype:
            """

            str_list = list()
            """ @type str_list: list[str | unicode] """

            # Group via UCSC super tracks.

            str_list.append('track Alignments\n')
            str_list.append('shortLabel Alignments\n')
            str_list.append('longLabel STAR alignments\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignments\n')
            str_list.append('\n')

            # Sample-specific tracks

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=True)

                if not paired_reads_dict:
                    # Skip Sample objects, which PairedReads objects have all been excluded.
                    continue

                runnable_sample = self.runnable_dict[self.get_prefix_sample(sample_name=sample.name)]
                file_path_sample = runnable_sample.file_path_object
                """ @type file_path_sample: bsf.analyses.aligner.FilePathSample """

                #
                # Add a trackDB entry for each Tophat accepted_hits.bam file.
                #
                # Common settings
                str_list.append('track ' + sample.name + '_alignments\n')
                str_list.append('type bam\n')
                str_list.append('shortLabel ' + sample.name + '_alignments\n')
                str_list.append('longLabel ' + sample.name + ' STAR alignments\n')
                str_list.append('bigDataUrl ' + file_path_sample.merged_bam + '\n')
                # str_list.append('html ...\n')
                str_list.append('visibility dense\n')

                # Common optional settings
                str_list.append('color 0,0,0\n')

                # bam/cram - Compressed Sequence Alignment track settings
                # None

                # Composite track settings
                str_list.append('parent Alignments\n')
                str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return
