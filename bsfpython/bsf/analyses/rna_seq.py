"""bsf.analyses.rna_seq

A package of classes and methods supporting RNA-Seq analyses.
"""

#
# Copyright 2013 - 2017 Michael K. Schuster
#
# Biomedical Sequencing Facility (BSF), part of the genomics core facility
# of the Research Center for Molecular Medicine (CeMM) of the
# Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
# This file is part of BSF Python.
#
# BSF Python is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BSF Python is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.


import errno
import os
from pickle import Pickler, HIGHEST_PROTOCOL
import re

from bsf import Analysis, FilePath, Runnable
from bsf.analyses.star_aligner import StarAligner
from bsf.annotation import AnnotationSheet, TuxedoSamplePairSheet
from bsf.executables import TopHat
from bsf.process import Executable, RunnableStep, RunnableStepLink
from bsf.standards import Default


class FilePathTophat(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathTophat} models files in a sample-specific TopHat directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathTophat} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathTophat, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.accepted_hits_bam = os.path.join(prefix, 'accepted_hits.bam')
        self.accepted_hits_bam_link_source = 'accepted_hits.bam'
        self.accepted_hits_bam_link_target = os.path.join(prefix, prefix + '_accepted_hits.bam')
        self.accepted_hits_bai = os.path.join(prefix, 'accepted_hits.bam.bai')
        self.accepted_hits_bai_link_source = 'accepted_hits.bam.bai'
        self.accepted_hits_bai_link_target = os.path.join(prefix, prefix + '_accepted_hits.bam.bai')
        self.accepted_hits_bw = os.path.join(prefix, 'accepted_hits.bw')
        self.align_summary = os.path.join(prefix, 'align_summary.txt')
        self.deletions_bed = os.path.join(prefix, 'deletions.bed')
        self.deletions_bb = os.path.join(prefix, 'deletions.bb')
        self.insertions_bed = os.path.join(prefix, 'insertions.bed')
        self.insertions_bb = os.path.join(prefix, 'insertions.bb')
        self.junctions_bed = os.path.join(prefix, 'junctions.bed')
        self.junctions_bb = os.path.join(prefix, 'junctions.bb')
        self.prep_reads_info = os.path.join(prefix, 'prep_reads.info')
        self.unmapped_bam = os.path.join(prefix, 'unmapped.bam')
        self.unmapped_bam_link_source = 'unmapped.bam'
        self.unmapped_bam_link_target = os.path.join(prefix, prefix + '_unmapped.bam')

        return


class FilePathCufflinks(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCufflinks} models files in a sample-specific Cufflinks directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCufflinks} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCufflinks, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.transcripts_gtf = os.path.join(prefix, 'transcripts.gtf')
        self.temporary_gene_prediction = os.path.join(prefix, 'transcripts_gene_prediction.tsv')
        self.temporary_big_gene_prediction = os.path.join(prefix, 'transcripts_big_gene_prediction.tsv')
        self.temporary_sorted_tsv = os.path.join(prefix, 'transcripts_sorted.tsv')
        self.transcripts_bb = os.path.join(prefix, 'transcripts.bb')
        self.transcripts_gtf_link_source = 'transcripts.gtf'
        self.transcripts_gtf_link_target = os.path.join(prefix, prefix + '_transcripts.gtf')
        self.transcripts_bb_link_source = 'transcripts.bb'
        self.transcripts_bb_link_target = os.path.join(prefix, prefix + '_transcripts.bb')
        # TODO: Some files are missing.

        return


class FilePathCuffmerge(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffmerge} models files in a comparison-specific Cuffmerge directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCufflinks} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffmerge, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.assembly_txt = '_'.join((prefix, 'assembly.txt'))
        self.merged_gtf = os.path.join(prefix, 'merged.gtf')
        self.temporary_gene_prediction = os.path.join(prefix, 'merged_gene_prediction.tsv')
        self.temporary_big_gene_prediction = os.path.join(prefix, 'merged_big_gene_prediction.tsv')
        self.temporary_sorted_tsv = os.path.join(prefix, 'merged_sorted.tsv')
        self.merged_bb = os.path.join(prefix, 'merged.bb')
        self.merged_gtf_link_source = 'merged.gtf'
        self.merged_gtf_link_target = os.path.join(prefix, prefix + '_merged.gtf')
        self.merged_bb_link_source = 'merged.bb'
        self.merged_bb_link_target = os.path.join(prefix, prefix + '_merged.bb')

        return


class FilePathCuffquant(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffquant} models files in a sample-specific Cuffquant directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCuffquant} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffquant, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.abundances = os.path.join(prefix, 'abundances.cxb')
        # TODO: Remove after rearranging _create_tophat_...() and _create_cuffmerge_...() into run().
        # self.tophat_accepted_hits = os.path.join('_'.join(('rnaseq_tophat', paired_reads_name)), 'accepted_hits.bam')
        self.tophat_accepted_hits = ''

        return


class FilePathCuffnorm(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffnorm} models files in a comparison-specific Cuffnorm directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCuffnorm} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffnorm, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class FilePathCuffdiff(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathCuffdiff} models files in a comparison-specific Cuffdiff directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathCuffdiff} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathCuffdiff, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class Tuxedo(Analysis):
    """Tuxedo RNASeq C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_run_tophat: C{bsf.Stage.name} for the run Tophat stage
    @type stage_name_run_tophat: str
    @cvar stage_name_process_tophat: C{bsf.Stage.name} for the process Tophat stage
    @type stage_name_process_tophat: str
    @cvar stage_name_run_cufflinks: C{bsf.Stage.name} for the run Cufflinks stage
    @type stage_name_run_cufflinks: str
    @cvar stage_name_process_cufflinks: C{bsf.Stage.name} for the process Cufflinks stage
    @type stage_name_process_cufflinks: str
    @cvar stage_name_run_cuffmerge: C{bsf.Stage.name} for the Cuffmerge stage
    @type stage_name_run_cuffmerge: str
    @cvar stage_name_run_cuffquant: C{bsf.Stage.name} for the Cuffquant stage
    @type stage_name_run_cuffquant: str
    @cvar stage_name_run_cuffnorm: C{bsf.Stage.name} for the Cuffnorm stage
    @type stage_name_run_cuffnorm: str
    @cvar stage_name_run_cuffdiff: C{bsf.Stage.name} for the Cuffdiff stage
    @type stage_name_run_cuffdiff: str
    @cvar stage_name_process_cuffdiff: C{bsf.Stage.name} for the process Cuffdiff stage
    @type stage_name_process_cuffdiff: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar cmp_file: Comparison file
    @type cmp_file: str | unicode
    @ivar genome_index_path: Bowtie genome index path
    @type genome_index_path: str | unicode
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode
    @ivar transcriptome_version: Transcriptome version
    @type transcriptome_version: str
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | unicode
    @ivar transcriptome_gtf_path: Reference transcriptome GTF file path
    @type transcriptome_gtf_path: str | unicode
    @ivar mask_gtf_path: GTF file path to mask transcripts
    @type mask_gtf_path: str | unicode
    @ivar multi_read_correction: Apply multi-read correction
    @type multi_read_correction: bool
    @ivar library_type: Library type
        Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
    @type library_type: str
    @ivar novel_transcripts: Assemble novel transcripts
    @type novel_transcripts: bool
    @ivar no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
    @type no_length_correction: bool
    """

    name = 'RNA-seq Analysis'
    prefix = 'rnaseq'

    # Replicate stage
    stage_name_run_tophat = '_'.join((prefix, 'run_tophat'))
    stage_name_process_tophat = '_'.join((prefix, 'process_tophat'))
    stage_name_run_cufflinks = '_'.join((prefix, 'run_cufflinks'))
    stage_name_process_cufflinks = '_'.join((prefix, 'process_cufflinks'))

    # Comparison stage
    stage_name_run_cuffmerge = '_'.join((prefix, 'cuffmerge'))
    stage_name_run_cuffquant = '_'.join((prefix, 'cuffquant'))
    stage_name_run_cuffnorm = '_'.join((prefix, 'cuffnorm'))
    stage_name_run_cuffdiff = '_'.join((prefix, 'cuffdiff'))
    stage_name_process_cuffdiff = '_'.join((prefix, 'process_cuffdiff'))

    @classmethod
    def get_prefix_rnaseq_run_tophat(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: Replicate key
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_tophat, paired_reads_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffmerge(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffmerge, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffquant(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffquant, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffnorm(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffnorm, comparison_name))

    @classmethod
    def get_prefix_rnaseq_run_cuffdiff(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_run_cuffdiff, comparison_name))

    @classmethod
    def get_prefix_rnaseq_process_cuffdiff(cls, comparison_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param comparison_name: Comparison key
        @type comparison_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_process_cuffdiff, comparison_name))

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
            comparisons=None,
            sample_list=None,
            replicate_grouping=False,
            cmp_file=None,
            genome_index_path=None,
            genome_fasta_path=None,
            transcriptome_version=None,
            transcriptome_index_path=None,
            transcriptome_gtf_path=None,
            mask_gtf_path=None,
            multi_read_correction=None,
            library_type=None,
            novel_transcripts=None,
            no_length_correction=False):
        """Initialise a C{bsf.analyses.rna_seq.Tuxedo} object.

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
        @param comparisons: Python C{dict} of Python C{str} (comparison name) key objects and
            Python C{tuple} value objects of C{bsf.ngs.Sample.name} and Python C{list} of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, (bsf.ngs.Sample.name, list[bsf.ngs.Sample])]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @param genome_index_path: Bowtie genome index path
        @type genome_index_path: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode
        @param mask_gtf_path: GTF file path to mask transcripts
        @type mask_gtf_path: str | unicode
        @param multi_read_correction: Apply multi-read correction
        @type multi_read_correction: bool
        @param library_type: Library type
            Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
        @type library_type: str
        @param novel_transcripts: Assemble novel transcripts
        @type novel_transcripts: bool
        @param no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
        @type no_length_correction: bool
        @return:
        @rtype:
        """

        super(Tuxedo, self).__init__(
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
            comparisons=comparisons,
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if cmp_file is None:
            self.cmp_file = str()
        else:
            self.cmp_file = cmp_file

        if genome_index_path is None:
            self.genome_index_path = str()
        else:
            self.genome_index_path = genome_index_path

        if genome_fasta_path is None:
            self.genome_fasta_path = str()
        else:
            self.genome_fasta_path = genome_fasta_path

        if transcriptome_version is None:
            self.transcriptome_version = str()
        else:
            self.transcriptome_version = transcriptome_version

        if transcriptome_index_path is None:
            self.transcriptome_index_path = str()
        else:
            self.transcriptome_index_path = transcriptome_index_path

        if transcriptome_gtf_path is None:
            self.transcriptome_gtf_path = str()
        else:
            self.transcriptome_gtf_path = transcriptome_gtf_path

        if mask_gtf_path is None:
            self.mask_gtf_path = str()
        else:
            self.mask_gtf_path = mask_gtf_path

        if multi_read_correction is None:
            self.multi_read_correction = False
        else:
            assert isinstance(multi_read_correction, bool)
            self.multi_read_correction = multi_read_correction

        if library_type is None:
            self.library_type = str()
        else:
            self.library_type = library_type

        if novel_transcripts is None:
            self.novel_transcripts = True
        else:
            self.novel_transcripts = novel_transcripts

        if no_length_correction is None:
            self.no_length_correction = False
        else:
            assert isinstance(no_length_correction, bool)
            self.no_length_correction = no_length_correction

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.rna_seq.Tuxedo} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(Tuxedo, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cmp_file = configuration.config_parser.get(section=section, option=option)

        option = 'genome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_version = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'mask_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mask_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'multi_read_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.multi_read_correction = configuration.config_parser.getboolean(section=section, option=option)

        option = 'library_type'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_type = configuration.config_parser.get(section=section, option=option)

        option = 'novel_transcripts'
        if configuration.config_parser.has_option(section=section, option=option):
            self.novel_transcripts = configuration.config_parser.getboolean(section=section, option=option)

        option = 'no_length_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.no_length_correction = configuration.config_parser.getboolean(section=section, option=option)

        return

    def _read_comparisons(self, cmp_file):
        """Read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

        All C{bsf.ngs.Sample} objects referenced in a comparison are added from the C{bsf.ngs.Collection} to the
        C{bsf.Analysis} object.

            - Column headers for CASAVA folders:
                - Treatment/Control/Point N ProcessedRunFolder:
                    - CASAVA processed run folder name or
                    - C{bsf.Analysis.input_directory} by default
                - Treatment/Control/Point N Project:
                    - CASAVA Project name or
                    - C{bsf.Analysis.project_name} by default
                - Treatment/Control/Point N Sample:
                    - CASAVA Sample name, no default
            - Column headers for independent samples:
                - Treatment/Control/Point N Sample:
                - Treatment/Control/Point N Reads:
                - Treatment/Control/Point N File:
        @param cmp_file: Comparisons file path
        @type cmp_file: str | unicode
        @return:
        @rtype:
        """

        if self.debug > 1:
            print '{!r} method _read_comparisons:'.format(self)

        regular_expression = re.compile(pattern='\W')

        annotation_sheet = AnnotationSheet.from_file_path(file_path=cmp_file)

        # TODO: Adjust by introducing a new class RNASeqComparisonSheet(AnnotationSheet) in this module?
        for row_dict in annotation_sheet.row_dicts:
            # In addition to defining samples, allow also the definition of groups in comparison files.
            # If the row dictionary has a 'Group' key, then the Sample in the same row gets added to the group.
            # So, 'ProcessedRunFolder', 'Project', 'Sample', 'Group' defines the groups, while ...
            # 'Control Group','Treatment Group' defines a comparison, as does ...
            # 'Control Group','Treatment ProcessedRunFolder','Treatment Project','Treatment Sample'

            # Get Sample objects for classical 'Control' and 'Treatment' keys,
            # before looking up a series of 'Point N' keys.
            i = -2
            key = str()
            comparison_groups = list()
            """ @type comparison_groups: list[(str, list[bsf.ngs.Sample])] """
            while True:
                i += 1
                if i == -1:
                    prefix = 'Control'
                elif i == 0:
                    prefix = 'Treatment'
                else:
                    prefix = 'Point {}'.format(i)
                # Get Sample objects for 'Point N' keys for as long as they are defined.
                # The Collection.get_sample_from_row_dict method can return one or more Sample objects,
                # depending on 'Group' or 'Sample' column entries.
                # In RNA-Seq experiments, entire pools of Sample objects (replicates) are compared with each other.
                group_name, group_samples = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                if group_name and len(group_samples):
                    key += group_name
                    key += '__'
                    comparison_groups.append((group_name, group_samples))
                    # Also expand each Python list of bsf.ngs.Sample objects to get all those bsf.ngs.Sample objects
                    # that this bsf.Analysis needs considering.
                    for sample in group_samples:
                        if self.debug > 1:
                            print '  {} Sample name: {!r} file_path: {!r}'.format(prefix, sample.name, sample.file_path)
                            # print sample.trace(1)
                        self.add_sample(sample=sample)
                elif i < 1:
                    # A Control and Treatment prefix is not required.
                    continue
                else:
                    # Only break if there is no further 'Point N' prefix.
                    break

            # For a successful comparison, more than one Sample (pool) has to be defined.

            if len(comparison_groups) < 2:
                if self.debug > 1:
                    print 'Comparison line with less than two Sample or Group keys. {!r}'. \
                        format(row_dict)
                continue

            if 'Comparison Name' in row_dict and row_dict['Comparison Name']:
                # For ridiculously large comparisons involving loads of groups or samples a comparison name can be
                # explicitly specified. Any non-word characters get replaced by underscore characters.
                key = re.sub(pattern=regular_expression, repl='_', string=row_dict['Comparison Name'])
            else:
                # Truncate the last '__' separator off the key string.
                key = key.rstrip('_')

            self.comparisons[key] = comparison_groups
            """ @type comparisons: dic[str, (str, list[bsf.ngs.Sample])] """

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.Tuxedo} analysis.
        @return:
        @rtype:
        """

        super(Tuxedo, self).run()

        # Get global defaults.

        default = Default.get_global_default()

        # Method configuration with regards to Cuffquant and Cuffdiff.
        run_cuffquant_before_cuffdiff = False

        # Tuxedo requires a genome version.

        if not self.genome_version:
            raise Exception('A Tuxedo analysis requires a genome_version configuration option.')

        # Tuxedo requires a transcriptome version.

        if not self.transcriptome_version:
            raise Exception('A Tuxedo analysis requires a transcriptome_version configuration option.')

        if self.cmp_file:
            # A comparison file path has been provided.
            if self.cmp_file == '*groups*':
                # The special file name *groups* creates the comparisons on the basis of an
                # all-against-all group comparison.
                # Without a comparison file path, simply add all Sample objects from the Collection.
                self.sample_list.extend(self.collection.get_all_samples())
                # Create a global comparison by adding all sample groups.
                comparison_groups = list()
                for key in self.collection.sample_group_dict.keys():
                    comparison_groups.append((key, self.collection.sample_group_dict[key]))
                # Sort the list of comparison groups by group name.
                comparison_groups.sort(cmp=lambda x, y: cmp(x[0], y[0]))
                self.comparisons['global'] = comparison_groups
            elif self.cmp_file == '*samples*':
                # The special file name *samples* creates the comparisons on the basis of an
                # all-against-all sample comparison.
                # Without a comparison file path, simply add all Sample objects from the Collection.
                self.sample_list.extend(self.collection.get_all_samples())
                # Create a global comparison by adding all samples under their sample name as group name.
                comparison_groups = list()
                for sample in self.sample_list:
                    # Add a tuple of group (i.e. sample) name and a Python list of the Sample object.
                    comparison_groups.append((sample.name, [sample]))
                # Sort the list of comparison groups by Sample.name.
                comparison_groups.sort(cmp=lambda x, y: cmp(x[0], y[0]))
                self.comparisons['global'] = comparison_groups
            else:
                # A comparison file path has been provided.
                # Expand an eventual user part i.e. on UNIX ~ or ~user and
                # expand any environment variables i.e. on UNIX ${NAME} or $NAME
                # Check if an absolute path has been provided, if not,
                # automatically prepend standard directory paths.
                self.cmp_file = os.path.expanduser(path=self.cmp_file)
                self.cmp_file = os.path.expandvars(path=self.cmp_file)
                if not os.path.isabs(self.cmp_file) and not os.path.exists(path=self.cmp_file):
                    self.cmp_file = os.path.join(self.project_directory, self.cmp_file)
                # Read and process the comparison file, which includes adding only those sample objects,
                # which are referenced in a comparison.
                self._read_comparisons(cmp_file=self.cmp_file)
        else:
            # Without a comparison file path, simply add all Sample objects from the Collection.
            # This means that only the initial pipeline stages, but not the comparison stage gets run.
            self.sample_list.extend(self.collection.get_all_samples())

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Sample objects are only put into the
        # bsf.Analysis.samples list by the _read_comparisons method.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        # Get the genome index

        if not self.genome_index_path:
            self.genome_index_path = os.path.join(
                Default.absolute_genome_resource(genome_version=self.genome_version),
                default.indices['bowtie2'],
                self.genome_version)

        if not self.genome_fasta_path:
            self.genome_fasta_path = self.genome_index_path + '.fa'

        # FIXME: For the moment use the FAI file instead of the UCSC chromosome sizes file.
        genome_sizes_path = self.genome_fasta_path + '.fai'

        if not os.path.exists(self.genome_fasta_path):
            raise Exception("Genome FASTA file path {!r} does not exists.".format(self.genome_fasta_path))

        # Define a reference transcriptome index directory or a GTF file path.

        if self.transcriptome_index_path:
            # Check if the transcriptome_index_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_index_path = Default.get_absolute_path(
                file_path=self.transcriptome_index_path,
                default_path=Default.absolute_transcriptome_resource(transcriptome_version=''))

            if not os.path.isdir(self.transcriptome_index_path):
                raise Exception("Reference transcriptome index directory {!r} does not exist.".
                                format(self.transcriptome_index_path))

            transcriptome_index = os.path.basename(self.transcriptome_index_path)

            # Does an indices_for_TopHat directory exist?
            transcriptome_index_path = os.path.join(self.transcriptome_index_path, default.indices['tophat'])
            if os.path.isdir(transcriptome_index_path):
                self.transcriptome_index_path = transcriptome_index_path

            # Finally, set the transcriptome GTF file path.
            # The tophat --transcript-index process puts a GFF file into the index directory
            # that really is a GTF file. A symbolic link to a GTF file is needed to make the
            # process cuffdiff script work.
            # For the moment, use the symbolic link in the indices_for_TopHat directory.

            self.transcriptome_gtf_path = os.path.join(
                self.transcriptome_index_path,
                '.'.join((transcriptome_index, 'gtf')))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        elif self.transcriptome_gtf_path:
            # Check if transcriptome_gtf_path is absolute and if not,
            # prepend the default transcriptomes directory.
            self.transcriptome_gtf_path = Default.get_absolute_path(
                file_path=self.transcriptome_gtf_path,
                default_path=Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        else:
            # Neither was provided, automatically discover on the basis of the transcriptome version.
            self.transcriptome_index_path = os.path.join(
                Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version),
                default.indices['tophat'],
                self.transcriptome_version,  # TopHat puts the transcriptome index into a sub directory.
                self.transcriptome_version)

            self.transcriptome_gtf_path = os.path.join(
                Default.absolute_transcriptome_resource(transcriptome_version=self.transcriptome_version),
                default.indices['tophat'],
                self.transcriptome_version) + '.gtf'
            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file path {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        if not self.transcriptome_gtf_path:
            raise Exception("Reference transcriptome GTF file not defined.")

        # Read configuration options.

        # TODO: Move the ConfigParser code.
        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        stage_run_tophat = self.get_stage(name=self.stage_name_run_tophat)
        stage_process_tophat = self.get_stage(name=self.stage_name_process_tophat)
        stage_run_cufflinks = self.get_stage(name=self.stage_name_run_cufflinks)
        stage_process_cufflinks = self.get_stage(name=self.stage_name_process_cufflinks)
        stage_run_cuffmerge = self.get_stage(name=self.stage_name_run_cuffmerge)
        stage_run_cuffquant = self.get_stage(name=self.stage_name_run_cuffquant)
        stage_run_cuffnorm = self.get_stage(name=self.stage_name_run_cuffnorm)
        stage_run_cuffdiff = self.get_stage(name=self.stage_name_run_cuffdiff)
        stage_process_cuffdiff = self.get_stage(name=self.stage_name_process_cuffdiff)

        runnable_run_cufflinks_list = list()
        """ @type runnable_run_cufflinks_list: list[Runnable] """

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                # Create a Tophat Runnable per paired_reads_name.

                # TODO: Activate the new code once the bsf_run_rnaseq_tophat.py script has been retired.

                # prefix_run_tophat = '_'.join((stage_run_tophat.name, paired_reads_name))

                # The output directory deviates from the prefix_run_tophat that itself is based on
                # stage_run_tophat.name. Both rnaseq_run_tophat and rnaseq_process_tophat processes should
                # use the same rnaseq_tophat directory.
                prefix_tophat = '_'.join(('rnaseq_tophat', paired_reads_name))

                file_path_tophat = FilePathTophat(prefix=prefix_tophat)

                # runnable_run_tophat = self.add_runnable(
                #         runnable=Runnable(
                #                 name=self.get_prefix_rnaseq_run_tophat(paired_reads_name=paired_reads_name),
                #                 code_module='bsf.runnables.generic',
                #                 working_directory=self.genome_directory,
                #                 file_path_object=file_path_tophat,
                #                 debug=self.debug))
                # executable_run_tophat = self.set_stage_runnable(
                #         stage=stage_run_tophat,
                #         runnable=runnable_run_tophat)
                #
                # Create a new Tophat RunnableStep.

                # tophat = runnable_run_tophat.add_runnable_step(
                #     runnable_step=RunnableStep(
                #     name='tophat2',
                #     program='tophat2'))

                tophat = TopHat(
                    name='_'.join(('rnaseq_tophat', paired_reads_name)),
                    analysis=self)

                # Set tophat options.

                tophat.add_option_long(
                    key='GTF',
                    value=self.transcriptome_gtf_path)
                if self.transcriptome_index_path:
                    tophat.add_option_long(
                        key='transcriptome-index',
                        value=self.transcriptome_index_path)
                tophat.add_option_long(
                    key='output-dir',
                    value=file_path_tophat.output_directory)
                tophat.add_option_long(
                    key='num-threads',
                    value=str(stage_run_tophat.threads))
                # TODO: These really are properties of the Reads, PairedReads or Sample objects rather than an Analysis.
                # TODO: Move the ConfigParser code.
                if config_parser.has_option(section=config_section, option='insert_size'):
                    insert_size = config_parser.getint(section=config_section, option='insert_size')
                    read_length = config_parser.getint(section=config_section, option='read_length')
                    mate_inner_dist = insert_size - 2 * read_length
                    tophat.add_option_long(
                        key='mate-inner-dist',
                        value=str(mate_inner_dist))
                # TODO: Move the ConfigParser code.
                if config_parser.has_option(section=config_section, option='mate-std-dev'):
                    tophat.add_option_long(
                        key='mate-std-dev',
                        value=config_parser.getint(section=config_section, option='mate-std-dev'))
                if self.library_type:
                    tophat.add_option_long(
                        key='library-type',
                        value=self.library_type)
                # The TopHat coverage search finds additional "GT-AG" introns, but is only recommended for
                # short reads (< 45 bp) and small read numbers (<= 10 M).
                # TODO: This option should possibly become configurable per sample.
                tophat.add_switch_long(key='no-coverage-search')
                # TODO: Set -rg-* options to back fill the read group from Illumina2bam.

                # Set rnaseq_tophat arguments.

                tophat.arguments.append(self.genome_index_path)

                # Set rnaseq_tophat arguments for reads1 and reads2.

                reads1 = list()
                reads2 = list()

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads_1:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads2.append(paired_reads.reads_2.file_path)

                tophat.arguments.append(','.join(reads1))
                tophat.arguments.append(','.join(reads2))

                # Create a new rnaseq_run_tophat Executable.
                # TODO: The following code block is required as long as the bsf_run_rnaseq_tophat.py script
                # has not been retired.

                pickler_dict_run_tophat = {
                    'prefix': stage_run_tophat.name,
                    'replicate_key': paired_reads_name,
                    'tophat_executable': tophat,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(stage_run_tophat.name, paired_reads_name))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_run_tophat)
                pickler_file.close()

                executable_run_tophat = stage_run_tophat.add_executable(
                    executable=Executable(
                        name='_'.join((stage_run_tophat.name, paired_reads_name)),
                        program='bsf_run_rnaseq_tophat.py'))
                # Set dependencies on previous Runnable or Executable objects.
                # None.

                # Set rnaseq_run_tophat options.
                self.set_command_configuration(command=executable_run_tophat)
                executable_run_tophat.add_option_long(key='pickler_path', value=pickler_path)
                executable_run_tophat.add_option_long(key='debug', value=str(self.debug))

                # Only submit this Executable if the "align_summary.txt" file does not exist.
                file_path_temporary = os.path.join(self.genome_directory, file_path_tophat.align_summary)
                if os.path.exists(file_path_temporary) and os.path.getsize(file_path_temporary) > 0:
                    executable_run_tophat.submit = False

                # TODO: End of code block.

                # Create a new rnaseq_process_tophat Executable.

                executable_process_tophat = stage_process_tophat.add_executable(
                    executable=Executable(
                        name='_'.join((stage_process_tophat.name, paired_reads_name)),
                        program='bsf_rnaseq_process_tophat2.sh'))
                # Set dependencies on previous Runnable or Executable objects.
                executable_process_tophat.dependencies.append(executable_run_tophat.name)

                # Set rnaseq_process_tophat options.

                # Set rnaseq_process_tophat arguments.
                self.set_command_configuration(command=executable_process_tophat)
                executable_process_tophat.arguments.append(file_path_tophat.output_directory)
                executable_process_tophat.arguments.append(genome_sizes_path)

                # Only submit this Executable if the "accepted_hits.bam.bai" file does not exist.
                file_path_temporary = os.path.join(
                    self.genome_directory,
                    file_path_tophat.accepted_hits_bai_link_source)
                if os.path.exists(file_path_temporary) and os.path.getsize(file_path_temporary):
                    executable_process_tophat.submit = False

                # TODO: Switch from an external Bash script to a set of Runnable and RunnableStep objects.
                # Since the Bash script includes Perl code to reset the BED score field to 0, rather than
                # re-scale it properly, it would be good to write a new bsf.runnables.process_tophat module
                # to implement this in Python code.

                # Create a process_tophat Runnable.

                # file_path_process_tophat = FilePathProcessTophat(prefix=prefix_tophat)
                #
                # runnable_process_tophat = self.add_runnable(runnable=Runnable(
                #     name='rnaseq_process_tophat',
                #     code_module='bsf.runnables.generic',
                #     working_directory=self.genome_directory,
                #     file_path_object=file_path_process_tophat,
                #     debug=self.debug))

                # Create a Cufflinks Runnable per paired_reads_name.

                prefix_run_cufflinks = '_'.join((stage_run_cufflinks.name, paired_reads_name))

                prefix_cufflinks = '_'.join(('rnaseq_cufflinks', paired_reads_name))
                # The output directory deviates from the prefix_run_cufflinks that itself is based on
                # stage_run_cufflinks.name. Both rnaseq_run_cufflinks and rnaseq_process_cufflinks processes
                # should use the same rnaseq_cufflinks directory.

                file_path_cufflinks = FilePathCufflinks(prefix=prefix_cufflinks)

                runnable_run_cufflinks = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_run_cufflinks,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        file_path_object=file_path_cufflinks,
                        debug=self.debug))
                executable_run_cufflinks = self.set_stage_runnable(
                    stage=stage_run_cufflinks,
                    runnable=runnable_run_cufflinks)
                # Set dependencies for subsequent Runnable or Executable objects.
                runnable_run_cufflinks_list.append(runnable_run_cufflinks)
                # Set dependencies on previous Runnable or Executable objects.
                executable_run_cufflinks.dependencies.append(executable_run_tophat.name)

                # Create a new Cufflinks RunnableStep.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='cufflinks',
                        program='cufflinks'))
                """ @type runnable_step: RunnableStep """

                # Cufflinks has a GTF option, in which case it will not assemble
                # novel transcripts and a GTF-guide option in which case it will
                # assemble novel transcripts.

                if self.novel_transcripts:
                    runnable_step.add_option_long(
                        key='GTF-guide',
                        value=self.transcriptome_gtf_path)
                else:
                    runnable_step.add_option_long(
                        key='GTF',
                        value=self.transcriptome_gtf_path)

                if self.mask_gtf_path:
                    runnable_step.add_option_long(
                        key='mask-file',
                        value=self.mask_gtf_path)

                runnable_step.add_option_long(
                    key='frag-bias-correct',
                    value=self.genome_fasta_path)

                if self.multi_read_correction:
                    runnable_step.add_switch_long(
                        key='multi-read-correct')

                if self.library_type:
                    runnable_step.add_option_long(
                        key='library-type',
                        value=self.library_type)

                # Cufflinks has a --library-norm-method option, but only one option (classic-fpkm) seems supported.

                if self.no_length_correction:
                    runnable_step.add_switch_long(
                        key='no-length-correction')

                runnable_step.add_option_long(
                    key='output-dir',
                    value=file_path_cufflinks.output_directory)

                runnable_step.add_option_long(
                    key='num-threads',
                    value=str(stage_run_cufflinks.threads))

                runnable_step.add_switch_long(
                    key='quiet')

                runnable_step.add_switch_long(
                    key='no-update-check')

                # Set Cufflinks arguments.

                runnable_step.arguments.append(file_path_tophat.accepted_hits_bam)

                # Convert the resulting transcripts GTF file into a UCSC genePred file.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='gtf_to_gp',
                        program='gtfToGenePred'))
                """ @type runnable_step: RunnableStep """
                runnable_step.add_switch_short(key='genePredExt')
                runnable_step.arguments.append(file_path_cufflinks.transcripts_gtf)
                runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)

                # Convert the UCSC genePred into a UCSC bigGenePred file.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='gp_to_bgp',
                        program='genePredToBigGenePred',
                        obsolete_file_path_list=[
                            file_path_cufflinks.temporary_gene_prediction,
                        ]))
                """ @type runnable_step: RunnableStep """
                runnable_step.arguments.append(file_path_cufflinks.temporary_gene_prediction)
                runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)

                # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bed_sort',
                        program='bedSort',
                        obsolete_file_path_list=[
                            file_path_cufflinks.temporary_big_gene_prediction,
                        ]))
                """ @type runnable_step: RunnableStep """
                runnable_step.arguments.append(file_path_cufflinks.temporary_big_gene_prediction)
                runnable_step.arguments.append(file_path_cufflinks.temporary_sorted_tsv)

                # Convert the sorted UCSC bigGenePred into a bigBed file.

                runnable_step = runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bgp_to_bb',
                        program='bedToBigBed',
                        obsolete_file_path_list=[
                            file_path_cufflinks.temporary_sorted_tsv,
                        ]))
                """ @type runnable_step: RunnableStep """
                # FIXME: It would be good to allow options with and without an equal sign.
                # i.e. --type=bed12+8 versus --type bed12+8, which does not work.
                # How could this work with configuration.ini files?
                # FIXME: This needs to be configurable.
                # Either another instance variable is required or via a generic configuration option
                # for the RunnableStep.
                # TODO: The location of the autoSQL file needs to be configurable.
                runnable_step.add_switch_short(key='as=/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
                runnable_step.add_switch_short(key='tab')
                runnable_step.add_switch_short(key='type=bed12+8')
                runnable_step.arguments.append(file_path_cufflinks.temporary_sorted_tsv)
                runnable_step.arguments.append(genome_sizes_path)
                runnable_step.arguments.append(file_path_cufflinks.transcripts_bb)

                # Add a symbolic link for the transcripts bigBed file, that includes a sample name prefix.
                runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='link_transcripts_bb',
                        source_path=file_path_cufflinks.transcripts_bb_link_source,
                        target_path=file_path_cufflinks.transcripts_bb_link_target))

                # Add a symbolic link for the transcripts GTF file, that includes a sample name prefix.
                runnable_run_cufflinks.add_runnable_step(
                    runnable_step=RunnableStepLink(
                        name='link_transcripts_gtf',
                        source_path=file_path_cufflinks.transcripts_gtf_link_source,
                        target_path=file_path_cufflinks.transcripts_gtf_link_target))

        # Create one process_cufflinks Executable to process all sub-directories.

        if len(runnable_run_cufflinks_list):
            executable_process_cufflinks = stage_process_cufflinks.add_executable(
                executable=Executable(
                    name=stage_process_cufflinks.name,
                    program='bsf_rnaseq_process_cufflinks.R'))
            # Set dependencies on previous Runnable or Executable objects.
            for runnable_run_cufflinks in runnable_run_cufflinks_list:
                executable_process_cufflinks.dependencies.append(runnable_run_cufflinks.name)

            # Set process_cufflinks options.
            self.set_command_configuration(command=executable_process_cufflinks)
            executable_process_cufflinks.add_option_long(
                key='gtf-reference',
                value=self.transcriptome_gtf_path)
            executable_process_cufflinks.add_option_long(
                key='genome-version',
                value=self.genome_version)

        # The Cuffmerge process generates temporary files in the working directory that
        # have the same name for each comparison. If more than one Cuffmerge process runs at the same time,
        # contention occurs where more than one process writes to the same file and one process deletes files
        # expected to exist by another process.
        # Circumvent such a situation by introducing dependencies on previous Cuffmerge processes. Sigh.
        # TODO: Report this to the Cufflinks author.
        previous_cuffmerge_name = str()

        comparison_keys = self.comparisons.keys()
        comparison_keys.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_key in comparison_keys:
            # TODO: Should the comparison prefix also include the project name or number?
            prefix_cuffmerge = self.get_prefix_rnaseq_run_cuffmerge(comparison_name=comparison_key)

            file_path_cuffmerge = FilePathCuffmerge(prefix=prefix_cuffmerge)

            runnable_run_cuffmerge = self.add_runnable(
                runnable=Runnable(
                    name=prefix_cuffmerge,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_object=file_path_cuffmerge,
                    debug=self.debug))
            executable_run_cuffmerge = self.set_stage_runnable(
                stage=stage_run_cuffmerge,
                runnable=runnable_run_cuffmerge)
            # Set a dependency on the previous Cuffmerge process to avoid file contention.
            if previous_cuffmerge_name:
                executable_run_cuffmerge.dependencies.append(previous_cuffmerge_name)
            previous_cuffmerge_name = executable_run_cuffmerge.name

            # Create a new Cuffmerge RunnableStep.

            runnable_step_cuffmerge = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='cuffmerge',
                    program='cuffmerge'))
            """ @type runnable_step_cuffmerge: RunnableStep """

            # Set rnaseq_cuffmerge options.

            runnable_step_cuffmerge.add_option_long(
                key='output-dir',
                value=file_path_cuffmerge.output_directory)
            runnable_step_cuffmerge.add_option_long(
                key='num-threads',
                value=str(stage_run_cuffmerge.threads))
            runnable_step_cuffmerge.add_option_long(
                key='ref-gtf',
                value=self.transcriptome_gtf_path)
            runnable_step_cuffmerge.add_option_long(
                key='ref-sequence',
                value=self.genome_fasta_path)

            # Set rnaseq_cuffmerge arguments.

            # Create an assembly manifest file to merge all replicates of each Sample object.
            # This file requires an absolute path, because the working directory is not set at the stage of
            # job submission.

            assembly_path = os.path.join(self.genome_directory, file_path_cuffmerge.assembly_txt)
            assembly_file = open(assembly_path, 'w')

            # Process rnaseq_cuffmerge and rnaseq_cuffdiff arguments in parallel.

            # Create a Python list of Python list objects of Cuffquant abundances per comparison group.
            cuffdiff_cuffnorm_abundances = list()
            """ @type cuffdiff_cuffnorm_abundances: list[list[str | unicode]] """
            cuffdiff_cuffnorm_alignments = list()
            """ @type cuffdiff_cuffnorm_alignments: list[list[str | unicode]] """
            cuffdiff_cuffnorm_dependencies = list()
            """ @type cuffdiff_cuffnorm_dependencies: list[str] """
            cuffdiff_cuffnorm_labels = list()
            """ @type cuffdiff_cuffnorm_labels: list[str] """

            for group_name, group_samples in self.comparisons[comparison_key]:
                """ @type group_name: str """
                """ @type group_samples: list[Sample] """
                per_group_abundances_list = list()
                """ @type per_group_abundances_list: list[str | unicode] """
                per_group_alignments_list = list()
                """ @type per_group_alignments_list: list[str | unicode] """
                # Count samples that remain after removing excluded PairedReads objects.
                sample_count = 0

                for sample in group_samples:
                    """ @type sample: bsf.ngs.Sample """
                    paired_reads_dict = sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping,
                        exclude=True)

                    paired_reads_name_list = paired_reads_dict.keys()
                    if len(paired_reads_name_list):
                        sample_count += 1
                    else:
                        # Skip Sample objects, which PairedReads objects have all been excluded.
                        continue
                    paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                    for paired_reads_name in paired_reads_name_list:
                        # Add the Cufflinks assembled transcripts to the Cuffmerge manifest.

                        transcripts_path = os.path.join(
                            self.genome_directory,
                            '_'.join(('rnaseq_cufflinks', paired_reads_name)),
                            'transcripts.gtf')
                        assembly_file.write(transcripts_path + '\n')

                        # Wait for each Cufflinks replicate to finish, before Cuffmerge can run.

                        executable_run_cuffmerge.dependencies.append(
                            '_'.join((self.stage_name_run_cufflinks, paired_reads_name)))

                        # Create a Cuffquant Runnable per comparison (comparison_key) and replicate (paired_reads_name)
                        # on the basis of the Cuffmerge GTF file.

                        prefix_cuffquant = '_'.join((stage_run_cuffquant.name, comparison_key, paired_reads_name))

                        file_path_cuffquant = FilePathCuffquant(prefix=prefix_cuffquant)
                        # TODO: Remove after rearranging this method into the run method.
                        file_path_cuffquant.tophat_accepted_hits = \
                            os.path.join('_'.join(('rnaseq_tophat', paired_reads_name)), 'accepted_hits.bam')

                        runnable_run_cuffquant = self.add_runnable(
                            runnable=Runnable(
                                name=prefix_cuffquant,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                file_path_object=file_path_cuffquant,
                                debug=self.debug))
                        executable_run_cuffquant = self.set_stage_runnable(
                            stage=stage_run_cuffquant,
                            runnable=runnable_run_cuffquant)

                        # Each Cuffquant process depends on Cuffmerge.

                        executable_run_cuffquant.dependencies.append(executable_run_cuffmerge.name)

                        # Create a new cuffquant RunnableStep.

                        runnable_step_cuffquant = runnable_run_cuffquant.add_runnable_step(
                            runnable_step=RunnableStep(
                                name='cuffquant',
                                program='cuffquant'))
                        """ @type runnable_step_cuffquant: RunnableStep """

                        # Set Cuffquant options.

                        runnable_step_cuffquant.add_option_long(
                            key='output-dir',
                            value=file_path_cuffquant.output_directory)
                        runnable_step_cuffquant.add_option_long(
                            key='num-threads',
                            value=str(stage_run_cuffquant.threads))
                        if self.mask_gtf_path:
                            runnable_step_cuffquant.add_option_long(
                                key='mask-file',
                                value=self.mask_gtf_path)
                        runnable_step_cuffquant.add_option_long(
                            key='frag-bias-correct',
                            value=self.genome_fasta_path)
                        if self.multi_read_correction:
                            runnable_step_cuffquant.add_switch_long(
                                key='multi-read-correct')
                        if self.library_type:
                            runnable_step_cuffquant.add_option_long(
                                key='library-type',
                                value=self.library_type)
                        if self.no_length_correction:
                            runnable_step_cuffquant.add_switch_long(
                                key='no-length-correction')
                        runnable_step_cuffquant.add_switch_long(
                            key='quiet')
                        runnable_step_cuffquant.add_switch_long(
                            key='no-update-check')

                        # Set Cuffquant arguments.
                        # Add the Cuffmerge GTF file and the TopHat BAM file as Cuffquant arguments.

                        runnable_step_cuffquant.arguments.append(file_path_cuffmerge.merged_gtf)
                        runnable_step_cuffquant.arguments.append(file_path_cuffquant.tophat_accepted_hits)

                        # Add the Cuffquant abundances file to the Cuffdiff list.

                        per_group_abundances_list.append(file_path_cuffquant.abundances)

                        # Add the TopHat BAM file to the Cuffdiff alignments list.

                        per_group_alignments_list.append(file_path_cuffquant.tophat_accepted_hits)

                        # Add the Cuffquant Runnable process name to the Cuffdiff and Cuffnorm dependencies list.

                        cuffdiff_cuffnorm_dependencies.append(executable_run_cuffquant.name)

                if sample_count:
                    cuffdiff_cuffnorm_labels.append(group_name)
                    cuffdiff_cuffnorm_abundances.append(per_group_abundances_list)
                    cuffdiff_cuffnorm_alignments.append(per_group_alignments_list)

            assembly_file.close()

            # Add the assembly manifest file as Cuffmerge argument.
            runnable_step_cuffmerge.arguments.append(assembly_path)

            # Convert the resulting merged GTF file into a UCSC genePred file.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='gtf_to_gp',
                    program='gtfToGenePred'))
            """ @type runnable_step: RunnableStep """
            runnable_step.add_switch_short(key='genePredExt')
            runnable_step.arguments.append(file_path_cuffmerge.merged_gtf)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)

            # Convert the UCSC genePred into a UCSC bigGenePred file.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='gp_to_bgp',
                    program='genePredToBigGenePred',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.temporary_gene_prediction,
                    ]))
            """ @type runnable_step: RunnableStep """
            runnable_step.arguments.append(file_path_cuffmerge.temporary_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)

            # Run bedSort on the UCSC bigGenePred file to sort field 1 in lexicographic mode and 2 in numeric mode.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='bed_sort',
                    program='bedSort',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.temporary_big_gene_prediction,
                    ]))
            """ @type runnable_step: RunnableStep """
            runnable_step.arguments.append(file_path_cuffmerge.temporary_big_gene_prediction)
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)

            # Convert the sorted UCSC bigGenePred into a bigBed file.

            runnable_step = runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStep(
                    name='bgp_to_bb',
                    program='bedToBigBed',
                    obsolete_file_path_list=[
                        file_path_cuffmerge.temporary_sorted_tsv,
                    ]))
            """ @type runnable_step: RunnableStep """
            # FIXME: It would be good to allow options with and without an equal sign.
            # i.e. --type=bed12+8 versus --type bed12+8, which does not work.
            # How could this work with configuration.ini files?
            # FIXME: This needs to be configurable.
            # Either another instance variable is required or via a generic configuration option
            # for the RunnableStep.
            # TODO: The location of the autoSQL file needs to be configurable.
            runnable_step.add_switch_short(key='as=/scratch/lab_bsf/resources/UCSC/bigGenePred.as')
            runnable_step.add_switch_short(key='tab')
            runnable_step.add_switch_short(key='type=bed12+8')
            runnable_step.arguments.append(file_path_cuffmerge.temporary_sorted_tsv)
            runnable_step.arguments.append(genome_sizes_path)
            runnable_step.arguments.append(file_path_cuffmerge.merged_bb)

            # Add a symbolic link for the merged bigBed file, that includes a comparison name prefix.
            runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStepLink(
                    name='link_merged_bb',
                    source_path=file_path_cuffmerge.merged_bb_link_source,
                    target_path=file_path_cuffmerge.merged_bb_link_target))

            # Add a symbolic link for the merged GTF file, that includes a comparison name prefix.
            runnable_run_cuffmerge.add_runnable_step(
                runnable_step=RunnableStepLink(
                    name='link_merged_gtf',
                    source_path=file_path_cuffmerge.merged_gtf_link_source,
                    target_path=file_path_cuffmerge.merged_gtf_link_target))

            # Create a Cuffnorm Runnable per comparison.

            prefix_run_cuffnorm = '_'.join((stage_run_cuffnorm.name, comparison_key))

            file_path_run_cuffnorm = FilePathCuffnorm(prefix=prefix_run_cuffnorm)

            runnable_run_cuffnorm = self.add_runnable(
                runnable=Runnable(
                    name=prefix_run_cuffnorm,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_object=file_path_run_cuffnorm,
                    debug=self.debug))
            executable_run_cuffnorm = self.set_stage_runnable(
                stage=stage_run_cuffnorm,
                runnable=runnable_run_cuffnorm)

            executable_run_cuffnorm.dependencies.extend(cuffdiff_cuffnorm_dependencies)

            # Create a new Cuffnorm RunnableStep.

            runnable_step_cuffnorm = runnable_run_cuffnorm.add_runnable_step(
                runnable_step=RunnableStep(
                    name='cuffnorm',
                    program='cuffnorm'))
            """ @type runnable_step: RunnableStep """

            # Set Cuffnorm options.

            runnable_step_cuffnorm.add_option_long(
                key='output-dir',
                value=file_path_run_cuffnorm.output_directory)
            runnable_step_cuffnorm.add_option_long(
                key='labels',
                value=','.join(cuffdiff_cuffnorm_labels))
            runnable_step_cuffnorm.add_option_long(
                key='num-threads',
                value=str(stage_run_cuffnorm.threads))
            if self.library_type:
                runnable_step_cuffnorm.add_option_long(
                    key='library-type',
                    value=self.library_type)
            runnable_step_cuffnorm.add_switch_long(
                key='quiet')
            runnable_step_cuffnorm.add_switch_long(
                key='no-update-check')

            # Add the Cuffmerge GTF file as first Cuffnorm argument.
            runnable_step_cuffnorm.arguments.append(file_path_cuffmerge.merged_gtf)

            # Add the Cuffquant abundances files per point as Cuffnorm arguments.
            for per_group_abundances_list in cuffdiff_cuffnorm_abundances:
                runnable_step_cuffnorm.arguments.append(','.join(per_group_abundances_list))

            # Create a Cuffdiff Runnable per comparison.

            prefix_run_cuffdiff = '_'.join((stage_run_cuffdiff.name, comparison_key))

            file_path_run_cuffdiff = FilePathCuffdiff(prefix=prefix_run_cuffdiff)

            runnable_run_cuffdiff = self.add_runnable(
                runnable=Runnable(
                    name=prefix_run_cuffdiff,
                    code_module='bsf.runnables.generic',
                    working_directory=self.genome_directory,
                    file_path_object=file_path_run_cuffdiff,
                    debug=self.debug))
            executable_run_cuffdiff = self.set_stage_runnable(
                stage=stage_run_cuffdiff,
                runnable=runnable_run_cuffdiff)

            if run_cuffquant_before_cuffdiff:
                # Add all executable_run_cuffquant dependencies to the executable_run_cuffdiff process.
                executable_run_cuffdiff.dependencies.extend(cuffdiff_cuffnorm_dependencies)
            else:
                # Add the executable_run_cuffmerge dependency to the executable_run_cuffdiff process.
                executable_run_cuffdiff.dependencies.append(executable_run_cuffmerge.name)

            # Create a new Cuffdiff RunnableStep.

            runnable_step_cuffdiff = runnable_run_cuffdiff.add_runnable_step(
                runnable_step=RunnableStep(
                    name='cuffdiff',
                    program='cuffdiff'))
            """ @type runnable_step: RunnableStep """

            # Set Cuffdiff options.

            runnable_step_cuffdiff.add_option_long(
                key='output-dir',
                value=file_path_run_cuffdiff.output_directory)
            runnable_step_cuffdiff.add_option_long(
                key='labels',
                value=','.join(cuffdiff_cuffnorm_labels))
            runnable_step_cuffdiff.add_option_long(
                key='num-threads',
                value=str(stage_run_cuffdiff.threads))
            if self.mask_gtf_path:
                runnable_step_cuffdiff.add_option_long(
                    key='mask-file',
                    value=self.mask_gtf_path)
            runnable_step_cuffdiff.add_option_long(
                key='frag-bias-correct',
                value=self.genome_fasta_path)
            if self.multi_read_correction:
                runnable_step_cuffdiff.add_switch_long(
                    key='multi-read-correct')
            if self.library_type:
                runnable_step_cuffdiff.add_option_long(
                    key='library-type',
                    value=self.library_type)
            if self.no_length_correction:
                runnable_step_cuffdiff.add_switch_long(
                    key='no-length-correction')
            runnable_step_cuffdiff.add_switch_long(
                key='quiet')
            runnable_step_cuffdiff.add_switch_long(
                key='no-update-check')

            # Add the Cuffmerge GTF file as first Cuffdiff argument.
            runnable_step_cuffdiff.arguments.append(file_path_cuffmerge.merged_gtf)

            # Cuffdiff seems to have a problem with Cuffquant abundances files in that the isoforms.count_tracking
            # files show ridiculously low numbers such as 1e-319 for some splice variants. Usually, other splice
            # variants in the same cluster seem fine.
            if run_cuffquant_before_cuffdiff:
                # Add the Cuffquant abundances files per comparison group as Cuffdiff arguments.
                for per_group_abundances_list in cuffdiff_cuffnorm_abundances:
                    runnable_step_cuffdiff.arguments.append(','.join(per_group_abundances_list))
            else:
                # Add the TopHat BAM files per comparison group as Cuffdiff arguments.
                for per_group_alignments_list in cuffdiff_cuffnorm_alignments:
                    runnable_step_cuffdiff.arguments.append(','.join(per_group_alignments_list))

            # Create a new rnaseq_process_cuffdiff Executable.

            executable_process_cuffdiff = stage_process_cuffdiff.add_executable(
                executable=Executable(
                    name='_'.join((stage_process_cuffdiff.name, comparison_key)),
                    program='bsf_rnaseq_process_cuffdiff.R'))
            executable_process_cuffdiff.dependencies.append(executable_run_cuffdiff.name)

            # Set rnaseq_process_cuffdiff options.
            self.set_command_configuration(command=executable_process_cuffdiff)
            executable_process_cuffdiff.add_option_long(
                key='comparison-name',
                value=comparison_key)
            executable_process_cuffdiff.add_option_long(
                key='gtf-assembly',
                value=file_path_cuffmerge.merged_gtf)
            executable_process_cuffdiff.add_option_long(
                key='gtf-reference',
                value=self.transcriptome_gtf_path)
            executable_process_cuffdiff.add_option_long(
                key='genome-version',
                value=self.genome_version)

            # Set rnaseq_process_cuffdiff arguments.

            # None so far.

        return

    def report(self):
        """Create a C{bsf.analyses.rna_seq.Tuxedo} report in HTML format and a UCSC Genome Browser Track Hub.
        @return:
        @rtype:
        """

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        # This code only needs the public URL.

        output_hub = str()

        # Write a HTML document.

        output_html = str()

        output_html += '<h1 id="{}_analysis">{} {}</h1>\n'.format(self.prefix, self.project_name, self.name)
        output_html += '\n'

        # TopHat and Cufflinks table.

        output_html += '<h2 id="transcriptome_browsing">Transcriptome Browsing</h2>\n'
        output_html += '\n'

        output_html += '<h3 id="read_alignments">Read Alignments</h3>\n'
        output_html += '\n'

        output_html += '<p id ="tophat">\n'
        # http://tophat.cbcb.umd.edu/manual.html
        output_html += '<strong><a href="http://ccb.jhu.edu/software/tophat/index.shtml">TopHat</a></strong> '
        output_html += 'aligns RNA-Seq reads to a genome in order to identify '
        output_html += 'exon-exon splice junctions. It is built on the ultra fast\n'
        output_html += 'short read mapping program\n'
        output_html += '<strong>'
        output_html += '<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie 2</a>'
        output_html += '</strong>.\n'
        output_html += '</p>\n'
        output_html += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = {
            'db': self.genome_version,
            'hubUrl': '{}/{}/rnaseq_hub.txt'.format(Default.url_absolute_projects(), link_name),
        }

        output_html += '<p id="track_hub">\n'
        output_html += 'View TopHat <strong>read alignments</strong> tracks for each sample\n'
        output_html += 'in their genomic context via the project-specific\n'
        output_html += 'UCSC Genome Browser Track Hub <a href="{}">{}</a>.\n'.format(
            self.ucsc_track_url(options_dict=options_dict),
            self.project_name)
        output_html += '</p>\n'
        output_html += '\n'

        output_html += '<p>\n'
        output_html += '<a href="rnaseq_tophat_alignment_summary.pdf"><img ' \
                       'alt="TopHat Alignment Summary" id="tophat_alignment_summary_img"' \
                       'src="rnaseq_tophat_alignment_summary.png" ' \
                       'height="80" ' \
                       'width="80" ' \
                       '/>' \
                       '</a>\n'
        output_html += 'Alignment summary statistics <a href="rnaseq_tophat_alignment_summary.tsv">TSV</a>\n'
        output_html += '</p>\n'

        output_html += '<h3 id="alignment_events">Splice Junctions, Insertions and Deletions</h3>\n'
        output_html += '\n'

        output_html += '<p>\n'
        output_html += 'TopHat reports <strong>splice junctions</strong> on the basis of RNA-Seq\n'
        output_html += 'read alignments in UCSC BED track format.\n'
        output_html += 'Each junction consists of two connected BED blocks,\n'
        output_html += 'where each block is as long as the maximal overhang\n'
        output_html += 'of any read spanning the junction. The score is\n'
        output_html += 'the number of alignments spanning the junction.\n'
        output_html += 'UCSC BED tracks of <strong>insertions</strong> and\n'
        output_html += '<strong>deletions</strong> are also reported by TopHat.\n'
        output_html += '</p>\n'

        output_html += '<p>\n'
        output_html += 'View the corresponding TopHat tracks for junctions, deletions and insertions\n'
        output_html += 'for each sample in their genomic context via the project-specific\n'
        output_html += 'UCSC Genome Browser Track Hub <a href="{}">{}</a>.\n'.format(
            self.ucsc_track_url(options_dict=options_dict),
            self.project_name)
        output_html += '</p>\n'
        output_html += '\n'

        # output += '<p>\n'
        # output += 'Follow the links below to attach\n'
        # output += 'Tophat junction, deletion and insertion annotation to the\n'
        # output += 'UCSC Genome Browser. Since each file needs transferring to\n'
        # output += 'the UCSC site, subsequent pages will take some time to load.\n'
        # output += '</p>\n'

        output_html += '<h2 id="gene_expression_profiles">Gene Expression Profiles</h2>\n'
        output_html += '\n'

        output_html += '<p id="cufflinks">\n'
        # http://cufflinks.cbcb.umd.edu/howitworks.html
        output_html += '<strong><a href="http://cole-trapnell-lab.github.io/cufflinks/">Cufflinks</a></strong>\n'
        output_html += 'assembles aligned RNA-Seq reads into transcripts,\n'
        output_html += 'estimates their abundances, and tests for differential\n'
        output_html += 'expression and regulation transcriptome-wide.\n'
        output_html += 'It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of\n'
        output_html += 'transcripts. Cufflinks then estimates the relative abundances of these transcripts based\n'
        output_html += 'on how many reads support each one, taking into account biases in library preparation '
        output_html += 'protocols.\n'
        output_html += '</p>\n'

        output_html += '<p>\n'
        output_html += 'The Cufflinks <strong>assembled transcripts</strong> can be attached to the \n'
        output_html += 'UCSC Genome Browser, by following the "Transcript Assembly" links\n'
        output_html += 'below.\n'
        output_html += 'The isoforms.fpkm_tracking and genes.fpkm_tracking files\n'
        output_html += 'contain the estimated isoform or gene expression values in the generic\n'
        # http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
        output_html += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                       'fpkm-tracking-format">FPKM Tracking format</a>.\n'
        output_html += 'The isoforms.count_tracking and genes.count_tracking files\n'
        output_html += 'contain the scaled isoform or gene count values in the generic\n'
        output_html += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                       'count-tracking-format">Count Tracking format</a>.\n'
        output_html += '</p>\n'

        output_html += '<p>\n'
        output_html += 'Please see a more detailed description of\n'
        # http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output
        output_html += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                       'output-formats-used-in-the-cufflinks-suite">Cufflinks output</a>.\n'
        output_html += '</p>\n'

        output_html += '<table id="gene_expression_table">\n'
        output_html += '<thead>\n'
        output_html += '<tr>\n'
        output_html += '<th>Sample</th>\n'
        output_html += '<th>Assembled Transcripts</th>\n'
        output_html += '<th>Gene FPKM</th>\n'
        output_html += '<th>Transcript FPKM</th>\n'
        output_html += '<th>Genes (Symbols)</th>\n'
        output_html += '<th>Isoforms (Symbols)</th>\n'
        output_html += '<th>Aligned BAM file</th>\n'
        output_html += '<th>Aligned BAI file</th>\n'
        output_html += '<th>Unaligned BAM file</th>\n'
        output_html += '</tr>\n'
        output_html += '</thead>\n'
        output_html += '<tbody>\n'

        # Group via UCSC super tracks.

        output_hub += 'track Alignments\n'
        output_hub += 'shortLabel Alignments\n'
        output_hub += 'longLabel TopHat RNA-Seq read alignments\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group alignments\n'
        output_hub += '\n'

        output_hub += 'track Assemblies\n'
        output_hub += 'shortLabel Assemblies\n'
        output_hub += 'longLabel Cuffmerge transcript structures\n'
        output_hub += 'visibility full\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group assemblies\n'
        output_hub += '\n'

        output_hub += 'track Coverage\n'
        output_hub += 'shortLabel Coverage\n'
        output_hub += 'longLabel TopHat RNA-Seq alignment coverage\n'
        output_hub += 'visibility full\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group coverage\n'
        output_hub += '\n'

        output_hub += 'track Deletions\n'
        output_hub += 'shortLabel Deletions\n'
        output_hub += 'longLabel TopHat RNA-Seq deletions\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group alignments\n'
        output_hub += '\n'

        output_hub += 'track Insertions\n'
        output_hub += 'shortLabel Insertions\n'
        output_hub += 'longLabel TopHat RNA-Seq insertions\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group alignments\n'
        output_hub += '\n'

        output_hub += 'track Junctions\n'
        output_hub += 'shortLabel Junctions\n'
        output_hub += 'longLabel TopHat RNA-Seq splice junctions\n'
        output_hub += 'visibility show\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group alignments\n'
        output_hub += '\n'

        output_hub += 'track Transcripts\n'
        output_hub += 'shortLabel Transcripts\n'
        output_hub += 'longLabel Cufflinks transcript structures\n'
        output_hub += 'visibility show\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group transcripts\n'
        output_hub += '\n'

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                # TopHat produces accepted_hits.bam, deletions.bb,
                # insertions.bb and junctions.bb files.

                #
                # Add a trackDB entry for each accepted_hits.bam file.
                #

                # Common trackDb settings.

                output_hub += 'track {}_alignments\n'. \
                    format(paired_reads_name)
                output_hub += 'type bam\n'
                output_hub += 'shortLabel {}_alignments\n'. \
                    format(paired_reads_name)
                output_hub += 'longLabel {} TopHat RNA-Seq read alignments\n'. \
                    format(paired_reads_name)
                output_hub += 'bigDataUrl rnaseq_tophat_{}/accepted_hits.bam\n'. \
                    format(paired_reads_name)
                output_hub += 'visibility dense\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                output_hub += 'color {}\n'. \
                    format('0,0,0')

                # Compressed Sequence Alignment track settings.

                # None so far.

                # Composite track settings.

                output_hub += 'parent Alignments\n'
                output_hub += '\n'

                #
                # Add a trackDB entry for each accepted_hits.bw file.
                #

                # Common trackDB settings.

                output_hub += 'track {}_coverage\n'. \
                    format(paired_reads_name)
                # TODO: The bigWig type must declare the expected signal range.
                # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                output_hub += 'type bigWig\n'
                output_hub += 'shortLabel {}_coverage\n'. \
                    format(paired_reads_name)
                output_hub += 'longLabel {} TopHat RNA-Seq alignment coverage\n'. \
                    format(paired_reads_name)
                output_hub += 'bigDataUrl rnaseq_tophat_{}/accepted_hits.bw\n'. \
                    format(paired_reads_name)
                output_hub += 'visibility full\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                output_hub += 'color {}\n'. \
                    format('0,0,0')

                # bigWig - Signal graphing track settings.

                output_hub += 'alwaysZero on\n'
                output_hub += 'autoScale on\n'
                output_hub += 'graphTypeDefault bar\n'
                output_hub += 'maxHeightPixels 100:60:20\n'
                # track_output += 'maxWindowToQuery 10000000\n'
                # track_output += 'smoothingWindow 5\n'
                # track_output += 'transformFunc NONE\n'
                # track_output += 'viewLimits 0:45\n'
                # track_output += 'viewLimitsMax 0:50\n'
                # track_output += 'windowingFunction maximum\n'
                # track_output += 'yLineMark <#>\n'
                # track_output += 'yLineOnOff on \n'
                # track_output += 'gridDefault on\n'

                # Composite track settings.

                output_hub += 'parent Coverage\n'
                output_hub += 'centerLabelsDense off\n'
                output_hub += '\n'

                #
                # Add a trackDB entry for each deletions.bb file.
                #

                output_hub += 'track {}_deletions\n'. \
                    format(paired_reads_name)
                output_hub += 'type bigBed\n'
                output_hub += 'shortLabel {}_deletions\n'. \
                    format(paired_reads_name)
                output_hub += 'longLabel {} TopHat RNA-Seq deletions\n'. \
                    format(paired_reads_name)
                output_hub += 'bigDataUrl rnaseq_tophat_{}/deletions.bb\n'. \
                    format(paired_reads_name)
                output_hub += 'visibility hide\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                output_hub += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                output_hub += 'parent Deletions\n'
                output_hub += '\n'

                # Insertions

                output_hub += 'track insertions_{}\n'. \
                    format(paired_reads_name)
                output_hub += 'type bigBed\n'
                output_hub += 'shortLabel {}_insertions\n'. \
                    format(paired_reads_name)
                output_hub += 'longLabel {} TopHat RNA-Seq insertions\n'. \
                    format(paired_reads_name)
                output_hub += 'bigDataUrl rnaseq_tophat_{}/insertions.bb\n'. \
                    format(paired_reads_name)
                output_hub += 'visibility hide\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                output_hub += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                output_hub += 'parent Insertions\n'
                output_hub += '\n'

                # Junctions

                output_hub += 'track {}_junctions\n'. \
                    format(paired_reads_name)
                output_hub += 'type bigBed\n'
                output_hub += 'shortLabel {}_junctions\n'. \
                    format(paired_reads_name)
                output_hub += 'longLabel {} TopHat RNA-Seq splice junctions\n'. \
                    format(paired_reads_name)
                output_hub += 'bigDataUrl rnaseq_tophat_{}/junctions.bb\n'. \
                    format(paired_reads_name)
                output_hub += 'visibility pack\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                output_hub += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                output_hub += 'parent Junctions\n'
                output_hub += '\n'

                # Transcripts

                output_hub += 'track {}_transcripts\n'. \
                    format(paired_reads_name)
                output_hub += 'type bigGenePred\n'
                output_hub += 'shortLabel {}_transcripts\n'. \
                    format(paired_reads_name)
                output_hub += 'longLabel {} Cufflinks transcript assembly\n'. \
                    format(paired_reads_name)
                output_hub += 'bigDataUrl rnaseq_cufflinks_{}/transcripts.bb\n'. \
                    format(paired_reads_name)
                output_hub += 'visibility hide\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                output_hub += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                output_hub += 'parent Transcripts\n'
                output_hub += '\n'

                # Cufflinks produces genes.fpkm_tracking, isoforms.fpkm_tracking,
                # skipped.gtf and transcripts.gtf.

                prefix = 'rnaseq_cufflinks_{}'.format(paired_reads_name)

                output_html += '<tr>\n'
                output_html += '<td class="left">' \
                               '{}' \
                               '</td>\n'.format(paired_reads_name)
                output_html += '<td class="center">' \
                               '<a href="{}/{}_transcripts.gtf">Transcript Assembly</a>' \
                               '</td>\n'.format(prefix, prefix)
                output_html += '<td class="center">' \
                               '<a href="{}/genes.fpkm_tracking">Genes FPKM</a>' \
                               '</td>\n'.format(prefix)
                output_html += '<td class="center">' \
                               '<a href="{}/isoforms.fpkm_tracking">Isoforms FPKM</a>' \
                               '</td>\n'.format(prefix)
                output_html += '<td class="center">' \
                               '<a href="{}/{}_genes_fpkm_tracking.tsv">Genes (Symbols)</a>' \
                               '</td>\n'.format(prefix, prefix)
                output_html += '<td class="center">' \
                               '<a href="{}/{}_isoforms_fpkm_tracking.tsv">Isoforms (Symbols)</a>' \
                               '</td>\n'.format(prefix, prefix)
                output_html += '<td class="center">' \
                               '<a href="{}/rnaseq_tophat_{}_accepted_hits.bam">Aligned BAM</a>' \
                               '</td>\n'.format(prefix, paired_reads_name)
                output_html += '<td class="center">' \
                               '<a href="{}/rnaseq_tophat_{}_accepted_hits.bam.bai">Aligned BAI</a>' \
                               '</td>\n'.format(prefix, paired_reads_name)
                output_html += '<td class="center">' \
                               '<a href="{}/rnaseq_tophat_{}_unaligned.bam">Unaligned BAM</a>' \
                               '</td>\n'.format(prefix, paired_reads_name)
                output_html += '</tr>\n'

        output_html += '</tbody>\n'
        output_html += '</table>\n'
        output_html += '\n'

        # Cuffdiff produces cds_exp.diff, gene_exp.diff, isoform_exp.diff
        # promoters.diff, splicing.diff and tss_group_exp.diff amongst many others.

        output_html += '<h2 id="differential_expression">Differential Expression</h2>\n'
        output_html += '\n'

        output_html += '<p id="cuffdiff">\n'
        output_html += '<strong><a href="http://cufflinks.cbcb.umd.edu/howitworks.html#diff">Cuffdiff</a></strong>\n'
        output_html += 'finds significant changes in transcript\n'
        output_html += 'expression, splicing, and promoter use.'
        output_html += '</p>\n'
        output_html += '\n'

        output_html += '<h3 id="all_genes">All Genes</h3>\n'

        output_html += '<table id="differential_expression_table">\n'
        output_html += '<thead>\n'
        output_html += '<tr>\n'
        output_html += '<th>Comparison</th>\n'
        output_html += '<th>Samples</th>\n'
        output_html += '<th>Replicates</th>\n'
        output_html += '<th>Coding Sequences</th>\n'
        output_html += '<th>Genes</th>\n'
        output_html += '<th>Isoforms</th>\n'
        output_html += '<th>Promoters</th>\n'
        output_html += '<th>Splicing</th>\n'
        output_html += '<th>Transcription Start Sites</th>\n'
        output_html += '<th>Gene FPKM Replicates</th>\n'
        output_html += '<th>Gene Count Replicates</th>\n'
        output_html += '<th>Isoform FPKM Replicates</th>\n'
        output_html += '<th>Isoform Count Replicates</th>\n'
        output_html += '</tr>\n'
        output_html += '</thead>\n'
        output_html += '<tbody>\n'

        comparison_keys = self.comparisons.keys()
        comparison_keys.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_key in comparison_keys:
            # Assemblies Super Track

            output_hub += 'track {}_assembly\n'. \
                format(comparison_key)
            output_hub += 'type bigGenePred\n'
            output_hub += 'shortLabel {}_assembly\n'. \
                format(comparison_key)
            output_hub += 'longLabel {} Cufflinks transcript assembly\n'. \
                format(comparison_key)
            output_hub += 'bigDataUrl rnaseq_cuffmerge_{}/merged.bb\n'. \
                format(comparison_key)
            output_hub += 'visibility pack\n'
            # 'html' is missing from the common settings.

            # Common optional settings.

            output_hub += 'color {}\n'. \
                format('0,0,0')

            # Composite track settings.

            output_hub += 'parent Assemblies\n'
            output_hub += '\n'

            prefix = 'rnaseq_process_cuffdiff_{}'.format(comparison_key)

            # Link to comparison-specific symbolic links in the directory after cummeRbund processing.

            output_html += '<tr>\n'
            output_html += '<td class="left">' \
                           '{}' \
                           '</td>\n'.format(comparison_key)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_samples.tsv">Samples</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_replicates.tsv">Replicates</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_cds_exp_diff.tsv">Coding Sequences</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_exp_diff.tsv"><strong>Genes</strong></a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_exp_diff.tsv">Isoforms</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_promoters_diff.tsv">Promoters</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_splicing_diff.tsv">Splicing</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_tss_group_exp_diff.tsv">Transcription Start Sites</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_fpkm_replicates.tsv">Gene FPKM Replicates</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_counts_replicates.tsv">Gene Count Replicates</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_fpkm_replicates.tsv">Isoform FPKM Replicates</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_counts_replicates.tsv">Isoform Count Replicates</a>' \
                           '</td>\n'.format(prefix, prefix)
            output_html += '</tr>\n'

            # Read sample pair information if available.

            sample_pair_path = os.path.join(
                self.genome_directory,
                prefix,
                '_'.join((prefix, 'sample_pairs.tsv')))

            if os.path.exists(sample_pair_path):

                sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                for row_dict in sample_pair_sheet.row_dicts:
                    output_html += '<tr>\n'
                    output_html += '<td>' \
                                   '' \
                                   '</td>\n'  # Comparison
                    output_html += '<td class="left" colspan="3">' \
                                   '<strong>{}</strong> versus <strong>{}</strong>' \
                                   '</td>\n'.format(row_dict['V1'], row_dict['V2'])  # Sample
                    output_html += '<td class="center">' \
                                   '<a href="{}/{}_{}_{}_genes_diff.tsv"><strong>Genes</strong></a>' \
                                   '</td>\n'.format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '<td class="center">' \
                                   '<a href="{}/{}_{}_{}_isoforms_diff.tsv">Isoforms</a>' \
                                   '</td>\n'.format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '<td class="left" colspan="5">' \
                                   '' \
                                   '</td>\n'
                    output_html += '</tr>\n'

        output_html += '</tbody>\n'
        output_html += '</table>\n'
        output_html += '\n'

        output_html += '<h3 id="significant_genes">Significant Genes</h3>\n'

        output_html += '<table id="significant_genes_table">\n'
        output_html += '<thead>\n'
        output_html += '<tr>\n'
        output_html += '<th>Comparison</th>\n'
        output_html += '<th>Genes</th>\n'
        output_html += '<th>Isoforms</th>\n'
        output_html += '</tr>\n'
        output_html += '</thead>\n'
        output_html += '<tbody>\n'

        for comparison_key in comparison_keys:
            prefix = 'rnaseq_process_cuffdiff_{}'.format(comparison_key)

            output_html += '<tr>\n'
            output_html += '<td class="left">' \
                           '{}' \
                           '</td>\n'.format(comparison_key)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_significance_matrix.pdf">' \
                           '<img ' \
                           'alt="Significance Matrix Plot - Genes - {}" ' \
                           'src="{}/{}_genes_significance_matrix.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_significance_matrix.pdf">' \
                           '<img ' \
                           'alt="Significance Matrix Plot - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_significance_matrix.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)
            output_html += '</tr>\n'

        output_html += '</tbody>\n'
        output_html += '</table>\n'

        # Show cummeRbund quality plots.

        output_html += '<h2 id="quality_plots">Quality Plots</h2>\n'
        output_html += '\n'

        output_html += '<p>\n'
        output_html += '</p>\n'
        output_html += '\n'

        output_html += '<table id="quality_plots_table">\n'
        output_html += '<thead>\n'
        output_html += '<tr>\n'
        output_html += '<th>Comparison</th>\n'
        output_html += '<th>Dispersion Plot - Genes</th>\n'
        output_html += '<th>Dispersion Plot - Isoforms</th>\n'
        output_html += '<th>Squared Coefficient of Variation - Genes</th>\n'
        output_html += '<th>Squared Coefficient of Variation - Isoforms</th>\n'
        output_html += '<th>Density Plot without Replicates - Genes</th>\n'
        output_html += '<th>Density Plot with Replicates - Genes</th>\n'
        output_html += '<th>Density Plot without Replicates - Isoforms</th>\n'
        output_html += '<th>Density Plot with Replicates - Isoforms</th>\n'
        output_html += '<th>Box Plot without Replicates - Genes</th>\n'
        output_html += '<th>Box Plot with Replicates - Genes</th>\n'
        output_html += '<th>Box Plot without Replicates - Isoforms</th>\n'
        output_html += '<th>Box Plot with Replicates - Isoforms</th>\n'
        output_html += '<th>Scatter Matrix Plot - Genes</th>\n'
        output_html += '<th>Scatter Matrix Plot - Isoforms</th>\n'
        output_html += '<th>Dendrogram Plot</th>\n'
        output_html += '<th>Volcano Matrix Plot - Genes</th>\n'
        output_html += '<th>Multidimensional Scaling Plot - Genes</th>\n'
        output_html += '<th>Principal Component Analysis Plot - Genes</th>\n'
        output_html += '</tr>\n'
        output_html += '</thead>\n'
        output_html += '<tbody>\n'

        for comparison_key in comparison_keys:
            prefix = 'rnaseq_process_cuffdiff_{}'.format(comparison_key)

            output_html += '<tr>\n'
            output_html += '<td class="left">{}</td>\n'.format(comparison_key)

            # Dispersion Plots for Genes and Isoforms

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_dispersion.pdf">' \
                           '<img ' \
                           'alt="Dispersion Plot - Genes - {}" ' \
                           'src="{}/{}_genes_dispersion.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)
            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_dispersion.pdf">' \
                           '<img ' \
                           'alt="Dispersion Plot - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_dispersion.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms

            if os.path.exists(
                    path=os.path.join(self.genome_directory, '{}/{}_genes_scv.png'.format(prefix, prefix))):
                output_html += '<td class="center">' \
                               '<a href="{}/{}_genes_scv.pdf">' \
                               '<img ' \
                               'alt="Squared Coefficient of Variation (SCV) - Genes - {}" ' \
                               'src="{}/{}_genes_scv.png" ' \
                               'height="80" ' \
                               'width="80" ' \
                               '/>' \
                               '</a>' \
                               '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)
            else:
                output_html += '<td class="center"></td>\n'

            if os.path.exists(
                    path=os.path.join(self.genome_directory, '{}/{}_isoforms_scv.png'.format(prefix, prefix))):
                output_html += '<td class="center">' \
                               '<a href="{}/{}_isoforms_scv.pdf">' \
                               '<img ' \
                               'alt="Squared Coefficient of Variation (SCV) - Isoforms - {}" ' \
                               'src="{}/{}_isoforms_scv.png" ' \
                               'height="80" ' \
                               'width="80" ' \
                               '/>' \
                               '</a>' \
                               '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)
            else:
                output_html += '<td class="center"></td>\n'

            # Density Plots for Genes without and with Replicates

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_density_wo_replicates.pdf">' \
                           '<img ' \
                           'alt="Density Plot without Replicates - Genes- {}" ' \
                           'src="{}/{}_genes_density_wo_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_density_w_replicates.pdf">' \
                           '<img ' \
                           'alt="Density Plot with Replicates - Genes - {}" ' \
                           'src="{}/{}_genes_density_w_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Density Plots for Isoforms without and with Replicates

            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_density_wo_replicates.pdf">' \
                           '<img ' \
                           'alt="Density Plot without Replicates - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_density_wo_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_density_w_replicates.pdf">' \
                           '<img ' \
                           'alt="Density Plot with Replicates - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_density_w_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Box Plots for Genes without and with Replicates

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_box_wo_replicates.pdf">' \
                           '<img ' \
                           'alt="Box Plot without Replicates - Genes - {}" ' \
                           'src="{}/{}_genes_box_wo_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_box_w_replicates.pdf">' \
                           '<img ' \
                           'alt="Box Plot with Replicates - Genes - {}" ' \
                           'src="{}/{}_genes_box_w_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Box Plots for Isoforms with and without Replicates

            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_box_wo_replicates.pdf">' \
                           '<img ' \
                           'alt="Box Plot without Replicates - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_box_wo_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_box_w_replicates.pdf">' \
                           '<img ' \
                           'alt="Box Plot with Replicates - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_box_w_replicates.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Scatter Matrix Plot for Genes and Isoforms

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_scatter_matrix.pdf">' \
                           '<img ' \
                           'alt="Scatter Matrix Plot - Genes - {}" ' \
                           'src="{}/{}_genes_scatter_matrix.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            output_html += '<td class="center">' \
                           '<a href="{}/{}_isoforms_scatter_matrix.pdf">' \
                           '<img ' \
                           'alt="Scatter Matrix Plot - Isoforms - {}" ' \
                           'src="{}/{}_isoforms_scatter_matrix.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Dendrogram Plot for Genes

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_dendrogram.pdf">' \
                           '<img ' \
                           'alt="Dendrogram Plot - Genes - {}" ' \
                           'src="{}/{}_genes_dendrogram.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Volcano Matrix Plot for Genes

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_volcano_matrix.pdf">' \
                           '<img ' \
                           'alt="Volcano Matrix Plot - Genes - {}" ' \
                           'src="{}/{}_genes_volcano_matrix.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            # Multidimensional Scaling Plot for Genes

            if os.path.exists(
                    path=os.path.join(self.genome_directory, '{}/{}_genes_mds.png'.format(prefix, prefix))):
                output_html += '<td class="center">' \
                               '<a href="{}/{}_genes_mds.pdf">' \
                               '<img ' \
                               'alt="Multidimensional Scaling Plot - Genes - {}" ' \
                               'src="{}/{}_genes_mds.png" ' \
                               'height="80" ' \
                               'width="80" ' \
                               '/>' \
                               '</a>' \
                               '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)
            else:
                output_html += '<td></td>\n'

            # Principal Component Analysis Plot for Genes

            output_html += '<td class="center">' \
                           '<a href="{}/{}_genes_pca.pdf">' \
                           '<img ' \
                           'alt="Principal Component Analysis Plot - Genes - {}" ' \
                           'src="{}/{}_genes_pca.png" ' \
                           'height="80" ' \
                           'width="80" ' \
                           '/>' \
                           '</a>' \
                           '</td>\n'.format(prefix, prefix, comparison_key, prefix, prefix)

            output_html += '</tr>\n'

            # Read sample pair information if available.

            sample_pair_path = os.path.join(self.genome_directory, prefix, '_'.join((prefix, 'sample_pairs.tsv')))

            if os.path.exists(sample_pair_path):

                sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                for row_dict in sample_pair_sheet.row_dicts:
                    output_html += '<tr>\n'

                    output_html += '<td class="left"></td>\n'
                    output_html += '<td  class="left" colspan="10">' \
                                   '<strong>{}</strong> versus <strong>{}</strong>' \
                                   '</td>\n'.format(row_dict['V1'], row_dict['V2'])

                    output_html += '<td class="center">'
                    output_html += '<a href="{}/{}_{}_{}_genes_scatter.pdf">'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '<img ' \
                                   'alt="Scatter Plot on genes {} versus {}" ' \
                                   'src="{}/{}_{}_{}_genes_scatter.png" ' \
                                   'height="80" ' \
                                   'width="80" ' \
                                   '/>'. \
                        format(row_dict['V1'], row_dict['V2'], prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '</a>'
                    output_html += '</td>\n'

                    output_html += '<td class="center"></td>\n'

                    output_html += '<td class="center">'
                    output_html += '<a href="{}/{}_{}_{}_maplot.pdf">'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '<img ' \
                                   'alt="M vs A Plot on genes {} versus {}" ' \
                                   'src="{}/{}_{}_{}_maplot.png" ' \
                                   'height="80" ' \
                                   'width="80" ' \
                                   '/>'. \
                        format(row_dict['V1'], row_dict['V2'], prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '</a>'
                    output_html += '</td>\n'

                    output_html += '<td class="center">'
                    output_html += '<a href="{}/{}_{}_{}_genes_volcano.pdf">'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '<img ' \
                                   'alt="Volcano Plot on genes {} versus {}" ' \
                                   'src="{}/{}_{}_{}_genes_volcano.png" ' \
                                   'height="80" ' \
                                   'width="80" ' \
                                   '/>'. \
                        format(row_dict['V1'], row_dict['V2'], prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output_html += '</a>'
                    output_html += '</td>\n'

                    output_html += '<td  class="center"colspan="4"></td>\n'

                    output_html += '</tr>\n'

        output_html += '</tbody>\n'
        output_html += '</table>\n'
        output_html += '\n'

        self.report_to_file(content=output_html)
        self.ucsc_hub_to_file(content=output_hub)

        return


class FilePathDESeq(FilePath):
    """The C{bsf.analyses.rna_seq.FilePathDESeq} models files in a comparison-specific DESeq directory.
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.rna_seq.FilePathDESeq} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathDESeq, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class DESeq(Analysis):
    """DESeq RNASeq C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_deseq: C{bsf.Stage.name} for the run DESeq stage
    @type stage_name_deseq: str
    """

    name = 'RNA-seq Analysis'
    prefix = 'rnaseq'

    # Replicate stage
    stage_name_deseq = '_'.join((prefix, 'deseq'))

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
            comparisons=None,
            sample_list=None,
            replicate_grouping=False,
            cmp_file=None,
            genome_fasta_path=None,
            transcriptome_gtf_path=None,
            transcriptome_index_path=None,
            design_formula=None,
            mask_gtf_path=None,
            multi_read_correction=None,
            library_type=None,
            no_length_correction=False):
        """Initialise a C{bsf.analyses.rna_seq.DESeq} object.

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
        @param comparisons: Python C{dict} of Python C{str} (comparison name) key objects and
            Python C{tuple} value objects of C{bsf.ngs.Sample.name} and Python C{list} of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, (bsf.ngs.Sample.name, list[bsf.ngs.Sample])]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param transcriptome_gtf_path: Reference transcriptome GTF file path
        @type transcriptome_gtf_path: str | unicode
        @param transcriptome_index_path: Tophat transcriptome index path
        @type transcriptome_index_path: str | unicode
        @param design_formula: Design formula
        @type design_formula: str | unicode
        @param mask_gtf_path: GTF file path to mask transcripts
        @type mask_gtf_path: str | unicode
        @param multi_read_correction: Apply multi-read correction
        @type multi_read_correction: bool
        @param library_type: Library type
            Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
        @type library_type: str
        @param no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
        @type no_length_correction: bool
        @return:
        @rtype:
        """

        super(DESeq, self).__init__(
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
            comparisons=comparisons,
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if cmp_file is None:
            self.cmp_file = str()
        else:
            self.cmp_file = cmp_file

        if genome_fasta_path is None:
            self.genome_fasta_path = str()
        else:
            self.genome_fasta_path = genome_fasta_path

        if transcriptome_gtf_path is None:
            self.transcriptome_gtf_path = str()
        else:
            self.transcriptome_gtf_path = transcriptome_gtf_path

        if transcriptome_index_path is None:
            self.transcriptome_index_path = str()
        else:
            self.transcriptome_index_path = transcriptome_index_path

        if design_formula is None:
            self.design_formula = str()
        else:
            self.design_formula = design_formula

        if mask_gtf_path is None:
            self.mask_gtf_path = str()
        else:
            self.mask_gtf_path = mask_gtf_path

        if multi_read_correction is None:
            self.multi_read_correction = False
        else:
            assert isinstance(multi_read_correction, bool)
            self.multi_read_correction = multi_read_correction

        if library_type is None:
            self.library_type = str()
        else:
            self.library_type = library_type

        if no_length_correction is None:
            self.no_length_correction = False
        else:
            assert isinstance(no_length_correction, bool)
            self.no_length_correction = no_length_correction

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.rna_seq.Tuxedo} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(DESeq, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cmp_file = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index_path = configuration.config_parser.get(section=section, option=option)

        option = 'design_formula'
        if configuration.config_parser.has_option(section=section, option=option):
            self.design_formula = configuration.config_parser.get(section=section, option=option)

        option = 'mask_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.mask_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'multi_read_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.multi_read_correction = configuration.config_parser.getboolean(section=section, option=option)

        option = 'library_type'
        if configuration.config_parser.has_option(section=section, option=option):
            self.library_type = configuration.config_parser.get(section=section, option=option)

        option = 'no_length_correction'
        if configuration.config_parser.has_option(section=section, option=option):
            self.no_length_correction = configuration.config_parser.getboolean(section=section, option=option)

        return

    def run(self):
        """Run this C{bsf.analyses.rna_seq.Tuxedo} analysis.
        @return:
        @rtype:
        """

        super(DESeq, self).run()

        # Tuxedo requires a genome version.

        if not self.genome_version:
            raise Exception('A Tuxedo analysis requires a genome_version configuration option.')

        # For DESeq, all samples need adding to the Analysis regardless.
        # TODO: Sort out the comparison issue.
        for sample in self.collection.get_all_samples():
            self.add_sample(sample=sample)

        prefix = 'rnaseq_deseq_global'

        comparison_directory = os.path.join(self.genome_directory, prefix)

        if not os.path.isdir(comparison_directory):
            try:
                os.makedirs(comparison_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        annotation_sheet = AnnotationSheet(
            file_path=os.path.join(comparison_directory, prefix + '_comparisons.tsv'),
            file_type='excel-tab',
            name='DESeq Annotation',
            header=True)

        # Re-index the sample group dict by sample name.
        sample_dict = dict()
        for (group_name, sample_list) in self.collection.sample_group_dict.iteritems():
            for sample in sample_list:
                if sample.name not in sample_dict:
                    sample_dict[sample.name] = list()
                sample_dict[sample.name].append(group_name)

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            row_dict = dict()
            # Set the group from the original Sample Annotation Sheet.
            row_dict['sample'] = sample.name
            row_dict['group'] = sample_dict[sample.name][0]
            row_dict['bam_path'] = StarAligner.get_prefix_star_aligner_merge(sample_name=sample.name) + '.bam'
            row_dict['bai_path'] = StarAligner.get_prefix_star_aligner_merge(sample_name=sample.name) + '.bai'
            # Set additional columns from the Sample Annotation Sheet prefixed with "Sample DESeq ".
            for annotation_key in filter(lambda x: x.startswith('DESeq '), sample.annotation_dict.keys()):
                row_dict[annotation_key[6:]] = sample.annotation_dict[annotation_key][0]

            # for paired_reads_name in paired_reads_name_list:
                # Create a ??? Runnable per paired_reads_name.

                # prefix_run_tophat = '_'.join((stage_run_tophat.name, paired_reads_name))

                # file_path_deseq = FilePathDESeq(prefix=prefix)

            annotation_sheet.row_dicts.append(row_dict)
        # TODO: Create a directory for the comparison.
        # TODO: Write a data frame defining the experiment for DESeq2.

        # Process all row_dict objects to get the superset of field names.
        for row_dict in annotation_sheet.row_dicts:
            for key in row_dict.iterkeys():
                if key not in annotation_sheet.field_names:
                    annotation_sheet.field_names.append(key)

        # Write the Annotation Sheet to disk.
        annotation_sheet.to_file_path()

        return
