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
"""The :py:mod:`bsf.analyses.tophat` module provides classes and methods supporting the Tophat aligner.
"""
import os

from bsf.analyses.aligner import Aligner, \
    FilePathAlign as AlignerFilePathAlign, \
    FilePathSample
from bsf.connector import ConnectorFile
from bsf.executables.bed import RunnableStepRescaleScore
from bsf.procedure import ConcurrentRunnable
from bsf.process import RunnableStep, RunnableStepLink
from bsf.standards import Configuration, Index, StandardFilePath, Transcriptome


class FilePathAlign(AlignerFilePathAlign):
    """The :py:class:`bsf.analyses.tophat.FilePathAlign` class models file paths at the alignment stage.

    :ivar aligned_sam: An aligned sequence alignment map (SAM) file path.
    :type aligned_sam: str
    :ivar unaligned_bam: An unaligned binary alignment map (BAM) file path.
    :type unaligned_bam: str
    :ivar align_summary_txt_link_source: An alignment summary link source file pth.
    :type align_summary_txt_link_source: str
    :ivar align_summary_txt_link_target: An alignment summary link target file path.
    :type align_summary_txt_link_target: str
    :ivar deletions_raw: A raw Tophat2 deletions (BED) file path.
    :type deletions_raw: str
    :ivar deletions_bed: A deletions (BED) file path.
    :type deletions_bed: str
    :ivar deletions_bb: A deletions (bigBED) fil epath.
    :type deletions_bb: str
    :ivar insertions_raw: A raw Tophat2 insertions (BED) file path.
    :type insertions_raw: str
    :ivar insertions_bed: An insertions (BED) file path.
    :type insertions_bed: str
    :ivar insertions_bb: An insertions (bigBED) file path.
    :type insertions_bb: str
    :ivar junctions_raw: A raw Tophat2 junctions (BED) file path.
    :type junctions_raw: str
    :ivar junctions_bed: A junctions (BED) file path.
    :type junctions_bed: str
    :ivar junctions_bb: A junctions (bigBED) file path.
    :type junctions_bb: str
    """

    def __init__(self, prefix: str) -> None:
        """Initialise a :py:class:`bsf.analyses.tophat.FilePathAlign` object.

        :param prefix: A Python :py:class:`str` prefix representing a :py:attr:`bsf.procedure.Runnable.name` attribute.
        :type prefix: str
        """
        super(FilePathAlign, self).__init__(prefix=prefix)

        # Override the aligned_sam instance variable of the super-class.
        self.aligned_sam = os.path.join(prefix, 'accepted_hits.sam')

        self.unaligned_bam = os.path.join(prefix, 'unmapped.bam')

        self.align_summary_txt_link_source = 'align_summary.txt'
        self.align_summary_txt_link_target = os.path.join(prefix, '_'.join((prefix, 'align_summary.txt')))

        self.deletions_raw = os.path.join(prefix, 'deletions.bed')
        self.deletions_bed = os.path.join(prefix, '_'.join((prefix, 'deletions.bed')))
        self.deletions_bb = os.path.join(prefix, '_'.join((prefix, 'deletions.bb')))

        self.insertions_raw = os.path.join(prefix, 'insertions.bed')
        self.insertions_bed = os.path.join(prefix, '_'.join((prefix, 'insertions.bed')))
        self.insertions_bb = os.path.join(prefix, '_'.join((prefix, 'insertions.bb')))

        self.junctions_raw = os.path.join(prefix, 'junctions.bed')
        self.junctions_bed = os.path.join(prefix, '_'.join((prefix, 'junctions.bed')))
        self.junctions_bb = os.path.join(prefix, '_'.join((prefix, 'junctions.bb')))

        return


class Tophat2(Aligner):
    """The :py:class:`bsf.analyses.tophat.Tophat2` models the Tophat2 spliced aligner.

    :ivar insert_size: An insert size.
    :type insert_size: int | None
    :ivar insert_size_sd: An insert size standard deviation.
    :type insert_size_sd: int | None
    :ivar read_length: A read length.
    :type read_length: int | None
    :ivar library_type: A Tophat2 library type
        (i.e., :literal:`fr-unstranded`, :literal:`fr-firststrand` or :literal:`fr-secondstrand`).
    :type library_type: str | None
    :ivar transcriptome_version: A transcriptome version.
    :type transcriptome_version: str | None
    :ivar transcriptome_gtf: A transcriptome annotation GTF file path.
    :type transcriptome_gtf: str | None
    :ivar transcriptome_index: A transcriptome index directory path.
    :type transcriptome_index: str
    :ivar threads_number: A number of threads.
    :type threads_number: int | None
    """

    name = 'Tophat2 Analysis'
    prefix = 'tophat2'

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

    @classmethod
    def get_file_path_align(cls, paired_reads_name: str) -> FilePathAlign:
        """Get a :py:class:`bsf.analyses.tophat.FilePathAlign` object from this or a subclass.

        :param paired_reads_name: A :py:attr:`bsf.ngs.PairedReads.name` attribute.
        :type paired_reads_name: str
        :return: A :py:class:`bsf.analyses.tophat.FilePathAlign` object or a subclass thereof.
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
            stage_list=None,
            collection=None,
            sample_list=None,
            skip_mark_duplicates=None,
            java_archive_picard=None,
            insert_size=None,
            insert_size_sd=None,
            read_length=None,
            library_type=None,
            transcriptome_version=None,
            transcriptome_gtf=None,
            transcriptome_index=None,
            threads_number=None):
        """Initialise a :py:class:`bsf.analyses.tophat.Tophat2` object.

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
        :param skip_mark_duplicates: Request skipping the Picard :literal:`MarkDuplicates` step.
        :type skip_mark_duplicates: bool | None
        :param java_archive_picard: A Picard tools :emphasis:`Java Archive` (JAR) file path.
        :type java_archive_picard: str | None
        :param insert_size: An insert size.
        :type insert_size: int | None
        :param insert_size_sd: An insert size standard deviation.
        :type insert_size_sd: int | None
        :param read_length: A read length.
        :type read_length: int | None
        :param library_type: A Tophat2 library type
            (i.e., :literal:`fr-unstranded`, :literal:`fr-firststrand` or :literal:`fr-secondstrand`).
        :type library_type: str | None
        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str | None
        :param transcriptome_gtf: A transcriptome annotation GTF file path.
        :type transcriptome_gtf: str | None
        :param transcriptome_index: A transcriptome index directory path.
        :type transcriptome_index: str
        :param threads_number: A number of threads.
        :type threads_number: int | None
        """
        super(Tophat2, self).__init__(
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

        self.insert_size = insert_size
        self.insert_size_sd = insert_size_sd
        self.read_length = read_length
        self.library_type = library_type

        self.transcriptome_version = transcriptome_version
        self.transcriptome_gtf = transcriptome_gtf
        self.transcriptome_index = transcriptome_index

        self.threads_number = threads_number

        return

    def set_configuration(self, configuration: Configuration, section: str) -> None:
        """Set instance variables of a :py:class:`bsf.analyses.tophat.Tophat2` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        super(Tophat2, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'insert_size'
        if configuration.config_parser.has_option(section=section, option=option):
            self.insert_size = configuration.config_parser.getint(section=section, option=option)

        option = 'insert_std_dev'
        if configuration.config_parser.has_option(section=section, option=option):
            self.insert_size_sd = configuration.config_parser.getint(section=section, option=option)

        option = 'read_length'
        if configuration.config_parser.has_option(section=section, option=option):
            self.read_length = configuration.config_parser.getint(section=section, option=option)

        option = 'library_type'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_type = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index = configuration.config_parser.get(section=section, option=option)

        option = 'threads_number'
        if configuration.config_parser.has_option(section=section, option=option):
            self.threads_number = configuration.config_parser.getint(section=section, option=option)

        return

    def add_runnable_step_aligner(self, runnable_align, stage_align, file_path_1, file_path_2):
        """Add one or more Tophat2-specific :py:class:`bsf.process.RunnableStep` objects
        to the :py:class:`bsf.procedure.ConcurrentRunnable` object.

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

        # Run Tophat2

        runnable_step = RunnableStep(
            name='tophat2',
            program='tophat',
            stdout=ConnectorFile(file_path=file_path_align.stdout_txt, file_mode='wt'),
            stderr=ConnectorFile(file_path=file_path_align.stderr_txt, file_mode='wt'),
            obsolete_file_path_list=[
                # The Tophat2 unmapped.bam is not needed after Picard MergeBamAlignment has run.
                file_path_align.unaligned_bam,
            ])
        runnable_align.add_runnable_step(runnable_step=runnable_step)

        self.set_runnable_step_configuration(runnable_step=runnable_step)

        # Set Tophat2 options.

        runnable_step.add_option_long(key='GTF', value=self.transcriptome_gtf)

        if self.transcriptome_index:
            runnable_step.add_option_long(key='transcriptome-index', value=self.transcriptome_index)

        runnable_step.add_option_long(key='output-dir', value=file_path_align.output_directory)

        runnable_step.add_option_long(key='num-threads', value=str(self.threads_number))

        # NOTE: In case a single Tophat2 thread is requested, the final BAM file is moved to "accepted_hits.bam",
        # which is not compatible with the named pipe of the same name. The tophat process needs to open
        # the named pipe fpr writing. The --no-convert-bam option implies that the Tophat2-internal samtools view
        # command writes to the named pipe.
        runnable_step.add_switch_long(key='no-convert-bam')

        # Picard MergeBamAlignments requires the SQ dictionary in FASTA order.
        # Without the --keep-fasta-order option, Tophat2 sorts the @SQ entries of the SAM header lexicographically.
        runnable_step.add_switch_long(key='keep-fasta-order')

        # Picard CleanSam expects a coordinate-sorted file, because the SAM header declares it.
        # Irrespective of the --no-sort-bam option, Tophat2 always declares the SAM header to be coordinate-sorted.
        # runnable_step.add_switch_long(key='no-sort-bam')

        # TODO: These really are properties of the Reads, PairedReads or Sample objects rather than an Analysis.
        # To query PairedReads or Sample annotation these objects would need passing into the
        # Aligner.add_runnable_step_aligner() method. The stage_align option could possibly be dropped, because
        if self.insert_size and self.read_length:
            runnable_step.add_option_long(key='mate-inner-dist', value=str(self.insert_size - 2 * self.read_length))

        if self.insert_size_sd:
            runnable_step.add_option_long(key='mate-std-dev', value=str(self.insert_size_sd))

        if self.library_type:
            runnable_step.add_option_long(key='library-type', value=self.library_type)

        # The TopHat coverage search finds additional 'GT-AG' introns, but is only recommended for
        # short reads (< 45 bp) and small read numbers (<= 10 M).
        # TODO: This option should possibly become configurable per sample.
        runnable_step.add_switch_long(key='no-coverage-search')

        # Set Tophat2 arguments.

        runnable_step.arguments.append(self.genome_index)
        runnable_step.arguments.append(file_path_1)
        if file_path_2 is not None:
            runnable_step.arguments.append(file_path_2)

        # Link the align_summary_txt file.

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStepLink(
                name='link_align_summary_txt',
                source_path=file_path_align.align_summary_txt_link_source,
                target_path=file_path_align.align_summary_txt_link_target))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStepRescaleScore(
                name='rescale_deletions',
                file_path_old=file_path_align.deletions_raw,
                file_path_new=file_path_align.deletions_bed,
                keep_header_lines=False))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStep(
                name='sort_deletions',
                program='sortBed',
                arguments=[file_path_align.deletions_bed, file_path_align.deletions_bed]))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStep(
                name='bed_to_bb_deletions',
                program='bedToBigBed',
                arguments=[
                    file_path_align.deletions_bed,
                    StandardFilePath.get_resource_genome_fasta_index(genome_version=self.genome_version),
                    file_path_align.deletions_bb],
                obsolete_file_path_list=[file_path_align.deletions_bed]))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStepRescaleScore(
                name='rescale_insertions',
                file_path_old=file_path_align.insertions_raw,
                file_path_new=file_path_align.insertions_bed,
                keep_header_lines=False))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStep(
                name='sort_insertions',
                program='sortBed',
                arguments=[file_path_align.insertions_bed, file_path_align.insertions_bed]))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStep(
                name='bed_to_bb_insertions',
                program='bedToBigBed',
                arguments=[
                    file_path_align.insertions_bed,
                    StandardFilePath.get_resource_genome_fasta_index(genome_version=self.genome_version),
                    file_path_align.insertions_bb],
                obsolete_file_path_list=[file_path_align.insertions_bed]))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStepRescaleScore(
                name='rescale_junctions',
                file_path_old=file_path_align.junctions_raw,
                file_path_new=file_path_align.junctions_bed,
                keep_header_lines=False))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStep(
                name='sort_junctions',
                program='sortBed',
                arguments=[file_path_align.junctions_bed, file_path_align.junctions_bed]))

        runnable_align.add_runnable_step_epilogue(
            runnable_step=RunnableStep(
                name='bed_to_bb_junctions',
                program='bedToBigBed',
                arguments=[
                    file_path_align.junctions_bed,
                    StandardFilePath.get_resource_genome_fasta_index(genome_version=self.genome_version),
                    file_path_align.junctions_bb],
                obsolete_file_path_list=[file_path_align.junctions_bed]))

        return

    def add_runnable_step_sample(self, runnable_sample, stage_sample):
        """Add one or more Tophat2-specific :py:class:`bsf.process.RunnableStep` objects
        to a :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable_sample: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :type runnable_sample: ConsecutiveRunnable
        :param stage_sample: A :py:class:`bsf.analysis.Stage` object.
        :type stage_sample: Stage
        """
        # NOTE: This method a copy of STAR.add_runnable_step_sample().
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

    def run(self):
        """Run a :py:class:`bsf.analyses.tophat.Tophat2` object.
        """
        # Check for the project name already here,
        # since the super class method has to be called later.
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Tophat2 requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_version' configuration option.")

        # Get the genome version before calling the run() method of the bsf.analysis.Analysis super-class.

        if not self.genome_version:
            self.genome_version = Transcriptome.get_genome(
                transcriptome_version=self.transcriptome_version)

        if not self.genome_version:
            raise Exception(f"A {self.name!s} requires a valid 'transcriptome_version' configuration option.")

        # The Bowtie2 genome index is quite peculiar as it has to be the index file path without file extensions.

        if not self.genome_index:
            self.genome_index = os.path.join(
                StandardFilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bowtie2'),
                self.genome_version)

        if not self.genome_fasta:
            self.genome_fasta = StandardFilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        # Define a reference transcriptome index directory or a GTF file path.

        if self.transcriptome_index:
            # Check if the transcriptome_index is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_index = self.configuration.get_absolute_path(
                file_path=self.transcriptome_index,
                default_path=StandardFilePath.get_resource_transcriptome(
                    transcriptome_version=None,
                    absolute=True))

            if not os.path.isdir(self.transcriptome_index):
                raise Exception(f'The reference transcriptome index directory {self.transcriptome_index!r} '
                                f'does not exist.')

            transcriptome_prefix = os.path.basename(self.transcriptome_index)

            # Does an indices_for_TopHat directory exist?
            transcriptome_index = os.path.join(
                self.transcriptome_index,
                Index.get(option='tophat2'))
            if os.path.isdir(transcriptome_index):
                self.transcriptome_index = transcriptome_index

            # Finally, set the transcriptome GTF file path.
            # The tophat --transcript-index process puts a GFF file into the index directory
            # that really is a GTF file. A symbolic link to a GTF file is needed to make the
            # process cuffdiff script work.
            # For the moment, use the symbolic link in the indices_for_TopHat directory.

            self.transcriptome_gtf = os.path.join(
                self.transcriptome_index,
                '.'.join((transcriptome_prefix, 'gtf')))

            if not os.path.exists(self.transcriptome_gtf):
                raise Exception(f'The reference transcriptome GTF file {self.transcriptome_gtf!r} does not exist.')
        elif self.transcriptome_gtf:
            # Check, if transcriptome_gtf_path is absolute and if not,
            # prepend the default transcriptome directory.
            self.transcriptome_gtf = self.configuration.get_absolute_path(
                file_path=self.transcriptome_gtf,
                default_path=StandardFilePath.get_resource_transcriptome(
                    transcriptome_version=self.transcriptome_version,
                    absolute=True))

            if not os.path.exists(self.transcriptome_gtf):
                raise Exception(f'The reference transcriptome GTF file {self.transcriptome_gtf!r} does not exist.')
        else:
            # Neither was provided, automatically discover on the basis of the transcriptome version.
            self.transcriptome_index = os.path.join(
                StandardFilePath.get_resource_transcriptome_index(
                    transcriptome_version=self.transcriptome_version,
                    transcriptome_index='tophat2'),
                self.transcriptome_version,  # TopHat puts the transcriptome index into a subdirectory.
                self.transcriptome_version)

            self.transcriptome_gtf = StandardFilePath.get_resource_transcriptome_gtf(
                transcriptome_version=self.transcriptome_version,
                transcriptome_index='tophat2',
                basic=True)

            if not os.path.exists(self.transcriptome_gtf):
                raise Exception(f'The reference transcriptome GTF file path {self.transcriptome_gtf!r} does not exist.')

        if not self.transcriptome_gtf:
            raise Exception(f"A {self.name!s} requires a 'transcriptome_index' or 'transcriptome_gtf' "
                            "configuration option.")

        if not self.threads_number:
            self.threads_number = 1

        super(Tophat2, self).run()

        return
