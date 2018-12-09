"""bsf.analyses.chip_seq

A package of classes and methods supporting ChIP-Seq analyses.
"""

#
# Copyright 2013 - 2018 Michael K. Schuster
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


from __future__ import print_function

import errno
import os
import pickle
import sys
import warnings

import bsf
import bsf.annotation
import bsf.executables
import bsf.ngs
import bsf.process
import bsf.standards


class ChIPSeqComparison(object):
    """ChIP-Seq comparison annotation sheet.

    Attributes:
    @ivar c_name: Control name
    @type c_name: str | None
    @ivar t_name: Treatment name
    @type t_name: str | None
    @ivar c_samples: Python C{list} of control C{bsf.ngs.Sample} objects
    @type c_samples: list[bsf.ngs.Sample]
    @ivar t_samples: Python C{list} of treatment C{bsf.ngs.Sample} objects
    @type t_samples: list[bsf.ngs.Sample]
    @ivar factor: ChIP factor
    @type factor: str
    @ivar tissue: Tissue
    @type tissue: str
    @ivar condition: Condition
    @type condition: str
    @ivar treatment: Treatment
    @type treatment: str
    @ivar replicate: replicate number
    @type replicate: int
    @ivar diff_bind: Run the DiffBind analysis
    @type diff_bind: bool
    """

    def __init__(
            self,
            c_name,
            t_name,
            c_samples,
            t_samples,
            factor,
            tissue=None,
            condition=None,
            treatment=None,
            replicate=None,
            diff_bind=None):
        """Initialise a C{bsf.analyses.chip_seq.ChIPSeqComparison}.

        @param c_name: Control name
        @type c_name: str | None
        @param t_name: Treatment name
        @type t_name: str | None
        @param c_samples: Python C{list} of control C{bsf.ngs.Sample} objects
        @type c_samples: list[bsf.ngs.Sample] | None
        @param t_samples: Python C{list} of treatment C{bsf.ngs.Sample} objects
        @type t_samples: list[bsf.ngs.Sample] | None
        @param factor: ChIP factor
        @type factor: str | None
        @param tissue: Tissue
        @type tissue: str | None
        @param condition: Condition
        @type condition: str | None
        @param treatment: Treatment
        @type treatment: str | None
        @param replicate: replicate number
        @type replicate: int | None
        @param diff_bind: Run the DiffBind analysis
        @type diff_bind: bool | None
        @return:
        @rtype:
        """

        super(ChIPSeqComparison, self).__init__()

        # Condition', 'Treatment', 'Replicate',
        # 'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat'

        self.c_name = c_name
        self.t_name = t_name

        if c_samples is None:
            self.c_samples = list()
        else:
            self.c_samples = c_samples

        if t_samples is None:
            self.t_samples = list()
        else:
            self.t_samples = t_samples

        if factor is None:
            self.factor = str()
        else:
            self.factor = factor

        if tissue is None:
            self.tissue = str()
        else:
            self.tissue = tissue

        if condition is None:
            self.condition = str()
        else:
            self.condition = condition

        if treatment is None:
            self.treatment = str()
        else:
            self.treatment = treatment

        if replicate is None:
            self.replicate = 0
        else:
            self.replicate = replicate

        if diff_bind is None:
            self.diff_bind = True
        else:
            self.diff_bind = diff_bind

        return


class FilePathAlignment(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chip_seq.FilePathAlignment} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathAlignment, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'aligned.sam')))
        self.cleaned_sam = os.path.join(prefix, '_'.join((prefix, 'cleaned.sam')))
        self.aligned_bam = os.path.join(prefix, prefix + '.bam')
        self.aligned_bai = os.path.join(prefix, prefix + '.bai')
        self.aligned_md5 = os.path.join(prefix, prefix + '.bam.md5')
        self.aligned_bai_link_source = prefix + '.bai'
        self.aligned_bai_link_target = os.path.join(prefix, prefix + '.bam.bai')

        return


class FilePathPeakCalling(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chip_seq.FilePathPeakCalling} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """

        super(FilePathPeakCalling, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.name_prefix = os.path.join(prefix, prefix)  # MACS2 --name option

        self.control_bdg = os.path.join(prefix, '_'.join((prefix, 'control_lambda.bdg')))
        self.control_bw = os.path.join(prefix, '_'.join((prefix, 'control_lambda.bw')))
        self.treatment_bdg = os.path.join(prefix, '_'.join((prefix, 'treat_pileup.bdg')))
        self.treatment_bw = os.path.join(prefix, '_'.join((prefix, 'treat_pileup.bw')))
        self.comparison_bdg = os.path.join(prefix, '_'.join((prefix, 'bdgcmp.bdg')))
        self.comparison_bw = os.path.join(prefix, '_'.join((prefix, 'bdgcmp.bw')))

        # self.peaks_bed = os.path.join(prefix, '_'.join((prefix, 'peaks.bed')))
        # self.peaks_bb = os.path.join(prefix, '_'.join((prefix, 'peaks.bb')))
        self.peaks_narrow = os.path.join(prefix, '_'.join((prefix, 'peaks.narrowPeak')))
        self.narrow_peaks_bed = os.path.join(prefix, '_'.join((prefix, 'narrow_peaks.bed')))
        self.narrow_peaks_bb = os.path.join(prefix, '_'.join((prefix, 'narrow_peaks.bb')))
        self.summits_bed = os.path.join(prefix, '_'.join((prefix, 'summits.bed')))
        self.summits_bb = os.path.join(prefix, '_'.join((prefix, 'summits.bb')))
        self.peaks_xls = os.path.join(prefix, '_'.join((prefix, 'peaks.xls')))
        self.model_r = os.path.join(prefix, '_'.join((prefix, 'model.r')))

        return


class FilePathDiffBind(bsf.FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chip_seq.FilePathDiffBind} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathDiffBind, self).__init__(prefix=prefix)

        self.output_directory = prefix
        self.sample_annotation_sheet = prefix + '_samples.csv'
        self.correlation_read_counts_png = os.path.join(prefix, prefix + '_correlation_read_counts.png')
        self.correlation_peak_caller_score_png = os.path.join(prefix, prefix + '_correlation_peak_caller_score.png')
        self.correlation_analysis_png = os.path.join(prefix, prefix + '_correlation_analysis.png')
        self.pca_png = os.path.join(prefix, prefix + '_pca_plot.png')
        self.contrasts_csv = os.path.join(prefix, prefix + '_contrasts.csv')

        return


class FilePathDiffBindContrast(bsf.FilePath):
    def __init__(self, prefix, group_1, group_2):
        """Initialise a C{bsf.analyses.chip_seq.FilePathDiffBind} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @param group_1: Group 1
        @type group_1: str
        @param group_2: Group 2
        @type group_2: str
        @return:
        @rtype
        """
        super(FilePathDiffBindContrast, self).__init__(prefix=prefix)

        suffix = group_1 + '__' + group_2

        self.ma_plot_png = os.path.join(prefix, prefix + '_ma_plot_' + suffix + '.png')
        self.scatter_plot_png = os.path.join(prefix, prefix + '_scatter_plot_' + suffix + '.png')
        self.pca_plot_png = os.path.join(prefix, prefix + '_pca_plot_' + suffix + '.png')
        self.box_plot_png = os.path.join(prefix, prefix + '_box_plot_' + suffix + '.png')
        self.dba_report_csv = os.path.join(prefix, prefix + '_report_' + suffix + '.csv')

        return


class ChIPSeqDiffBindSheet(bsf.annotation.AnnotationSheet):
    """ChIP-Seq Bioconductor DiffBind annotation sheet class.

    Attributes:
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'SampleID',
        'Tissue',
        'Factor',
        'Condition',
        'Treatment',
        'Replicate',
        'bamReads',
        'bamControl',
        'ControlID',
        'Peaks',
        'PeakCaller',
        'PeakFormat',
    ]

    _test_methods = {
        'SampleID': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'Tissue': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'Factor': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'Condition': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'Treatment': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'Replicate': [
            bsf.annotation.AnnotationSheet.check_numeric,
        ],
        'ControlID': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'PeakCaller': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
        'PeakFormat': [
            bsf.annotation.AnnotationSheet.check_alphanumeric,
        ],
    }

    def sort(self):
        """Sort by columns I{Tissue}, I{Factor}, I{Condition}, I{Treatment} and I{Replicate}.

        @return:
        @rtype:
        """
        self.row_dicts.sort(
            cmp=lambda x, y:
            cmp(x['Tissue'], y['Tissue']) or
            cmp(x['Factor'], y['Factor']) or
            cmp(x['Condition'], y['Condition']) or
            cmp(x['Treatment'], y['Treatment']) or
            cmp(int(x['Replicate']), int(y['Replicate'])))

        return

    def to_file_path(self, adjust_field_names=None):
        """Write a C{bsf.analyses.ChIPSeqDiffBindSheet} to a file.

        @param adjust_field_names: Clear and adjust the Python C{list} of Python C{str} field name objects
        @type adjust_field_names: bool
        @return:
        @rtype:
        """

        # Override the method from the super-class to automatically sort before writing to a file.

        self.sort()
        super(ChIPSeqDiffBindSheet, self).to_file_path(adjust_field_names=adjust_field_names)

        return


class ChIPSeq(bsf.Analysis):
    """The C{bsf.analyses.chip_seq.ChIPSeq} class represents the logic to run a ChIP-Seq-specific C{bsf.Analysis}.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_alignment: C{bsf.Stage.name} for the alignment stage
    @type stage_name_alignment: str
    @cvar stage_name_peak_calling: C{bsf.Stage.name} for the peak calling stage
    @type stage_name_peak_calling: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode | None
    @ivar colour_default: Default UCSC Genome Browser Track Hub RGB colour
    @type colour_default: str | None
    @ivar colour_dict: Python C{dict} of
        Python C{str} factor name key and
        Python C{str} RGB colour value data
    @type colour_dict: dict[str, str] | None
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode | None
    @ivar genome_sizes_path: Reference genome (chromosome) sizes file path
    @type genome_sizes_path: str | unicode | None
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: None | str | unicode
    """

    name = 'ChIP-seq Analysis'
    prefix = 'chipseq'

    stage_name_alignment = '_'.join((prefix, 'alignment'))
    stage_name_peak_calling = '_'.join((prefix, 'peak_calling'))
    stage_name_diff_bind = '_'.join((prefix, 'diff_bind'))

    @classmethod
    def get_prefix_chipseq_alignment(cls, paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param paired_reads_name: PairedReads name
        @type paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_alignment, paired_reads_name))

    @classmethod
    def get_prefix_chipseq_peak_calling(cls, t_paired_reads_name, c_paired_reads_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param t_paired_reads_name: Treatment PairedReads name
        @type t_paired_reads_name: str
        @param c_paired_reads_name: Control PairedReads name
        @type c_paired_reads_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return cls.stage_name_peak_calling + '_' + t_paired_reads_name + '__' + c_paired_reads_name

    @classmethod
    def get_prefix_diff_bind(cls, factor_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param factor_name: Factor name
        @type factor_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.stage_name_diff_bind, factor_name))

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
            replicate_grouping=True,
            comparison_path=None,
            genome_fasta_path=None,
            genome_sizes_path=None,
            colour_default=None,
            colour_dict=None,
            factor_default=None,
            classpath_picard=None):
        """Initialise a C{bsf.analyses.chip_seq.ChIPSeq}.

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
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param comparison_path: Comparison file path
        @type comparison_path: str | unicode | None
        @param colour_default: Default UCSC Genome Browser Track Hub RGB colour
        @type colour_default: str | None
        @param colour_dict: Python C{dict} of
            Python C{str} factor name key and
            Python C{str} RGB colour value data
        @type colour_dict: dict[str, str] | None
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode | None
        @param genome_sizes_path: Reference genome (chromosome) sizes file path
        @type genome_sizes_path: str | unicode | None
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: None | str | unicode
        @return:
        @rtype:
        """
        super(ChIPSeq, self).__init__(
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
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = True
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        self.comparison_path = comparison_path
        self.colour_default = colour_default
        self.colour_dict = colour_dict
        self.factor_default = factor_default
        self.genome_fasta_path = genome_fasta_path
        self.genome_sizes_path = genome_sizes_path
        self.classpath_picard = classpath_picard

        self._comparison_dict = dict()
        """ @type _comparison_dict: dict[str, ChIPSeqComparison]"""

        self._factor_dict = dict()
        """ @type _factor_dict: dict[str, list[ChIPSeqComparison]] """

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.chip_seq.ChIPSeq} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(ChIPSeq, self).set_configuration(configuration=configuration, section=section)

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.comparison_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_sizes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_sizes_path = configuration.config_parser.get(section=section, option=option)

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        option = 'colour_default'
        if configuration.config_parser.has_option(section=section, option=option):
            self.colour_default = configuration.config_parser.get(section=section, option=option)

        track_hub_section = '.'.join((section, 'TrackHubColours'))
        if configuration.config_parser.has_section(section=track_hub_section):
            if self.colour_dict is None:
                self.colour_dict = dict()

            for option in configuration.config_parser.options(section=track_hub_section):
                self.colour_dict[option] = configuration.config_parser.get(
                    section=track_hub_section,
                    option=option)

        option = 'factor_default'
        if configuration.config_parser.has_option(section=section, option=option):
            self.factor_default = configuration.config_parser.get(section=section, option=option)

        return

    def get_colour(self, factor):
        """Get a UCSC Genome Browser Track Hub RGB colour by a factor name.

        @param factor: ChIP-seq factor name
        @type factor: str
        @return: Track Hub RGB colour
        @rtype: str
        """
        if factor in self.colour_dict:
            return self.colour_dict[factor]
        else:
            return self.colour_default

    # Taken from ConfigParser.RawConfigParser.getboolean()

    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False}

    def run(self):
        """Run a C{bsf.analyses.chip_seq.ChIPSeq} C{bsf.Analysis}.

        @return:
        @rtype:
        """

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file from disk.

                - Column headers for CASAVA folders:
                    - Treatment/Control ProcessedRunFolder:
                        - CASAVA processed run folder name or
                        - C{bsf.Analysis.input_directory} by default
                    - Treatment/Control Project:
                        - CASAVA Project name or
                        - C{bsf.Analysis.project_name} by default
                    - Treatment/Control Sample:
                        - CASAVA Sample name, no default
                - Column headers for independent samples:
                    - Treatment/Control Group:
                    - Treatment/Control Sample:
                    - Treatment/Control Reads:
                    - Treatment/Control File:
            @return:
            @rtype:
            """
            annotation_sheet = bsf.annotation.AnnotationSheet.from_file_path(file_path=self.comparison_path)

            # Unfortunately, two passes through the comparison sheet are required.
            # In the first one merge all Sample objects that share the name.
            # Merging Sample objects is currently the only way to pool PairedReads objects,
            # which is required for ChIP-Seq experiments.

            sample_dict = dict()
            """ @type sample_dict: dict[str, bsf.ngs.Sample] """

            # First pass, merge Sample objects, if they have the same name.
            for row_dict in annotation_sheet.row_dicts:
                for prefix in ('Control', 'Treatment'):
                    name, sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                    for o_sample in sample_list:
                        if o_sample.name not in sample_dict:
                            sample_dict[o_sample.name] = o_sample
                        else:
                            n_sample = bsf.ngs.Sample.from_samples(
                                sample1=sample_dict[o_sample.name],
                                sample2=o_sample)
                            sample_dict[n_sample.name] = n_sample

            # Second pass, add all Sample objects mentioned in a comparison.
            level_1_dict = dict()
            """ @type level_1_dict: dict[str, dict[str, int]] """

            for row_dict in annotation_sheet.row_dicts:

                c_name, c_sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Control')
                t_name, t_sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Treatment')

                # ChIP-Seq experiments use the order treatment versus control in comparisons.
                comparison_key = '__'.join((t_name, c_name))

                # For a successful comparison, both Python list objects of Sample objects have to be defined.

                if not (len(t_sample_list) and len(c_sample_list)):
                    if self.debug > 1:
                        print(
                            'Redundant comparison line with Treatment:', repr(t_name),
                            'samples:', len(t_sample_list),
                            'and Control:', repr(c_name),
                            'samples:', len(c_sample_list))
                    continue

                # Add all control Sample or SampleGroup objects to the Sample list.

                for c_sample in c_sample_list:
                    if self.debug > 1:
                        print('  Control Sample name:', c_sample.name, 'file_path:', c_sample.file_path)
                    if self.debug > 2:
                        sys.stdout.writelines(c_sample.trace(level=1))
                        # Find the Sample in the unified sample dictionary.
                    if c_sample.name in sample_dict:
                        self.add_sample(sample=sample_dict[c_sample.name])

                # Add all treatment Sample or SampleGroup objects to the Sample list.

                for t_sample in t_sample_list:
                    if self.debug > 1:
                        print('  Treatment Sample name:', t_sample.name, 'file_path:', t_sample.file_path)
                    if self.debug > 2:
                        sys.stdout.writelines(t_sample.trace(level=1))
                    if t_sample.name in sample_dict:
                        self.add_sample(sample=sample_dict[t_sample.name])

                if 'Tissue' in row_dict:
                    tissue = row_dict['Tissue']
                else:
                    tissue = str()

                if 'Factor' in row_dict:
                    factor = row_dict['Factor']
                else:
                    factor = self.factor_default

                if 'Condition' in row_dict:
                    condition = row_dict['Condition']
                else:
                    condition = str()

                if 'Treatment' in row_dict:
                    treatment = row_dict['Treatment']
                else:
                    treatment = str()

                if 'DiffBind' in row_dict:
                    if row_dict['DiffBind'].lower() not in ChIPSeq._boolean_states:
                        raise ValueError("Value in field 'DiffBind' is not a boolean: " + repr(row_dict['DiffBind']))
                    diff_bind = ChIPSeq._boolean_states[row_dict['DiffBind'].lower()]
                else:
                    diff_bind = True

                # Automatically create replicate numbers for the sample annotation sheet required by the
                # DiffBind Bioconductor package.
                # Use a first-level dict with replicate key data and second-level dict value data.
                # The second-level dict stores Treatment Sample key data and int value data.

                if 'Replicate' in row_dict:
                    value = row_dict['Replicate']

                    if value not in level_1_dict:
                        level_1_dict[value] = dict()

                    level_2_dict = level_1_dict[value]
                    level_2_dict[comparison_key] = 0

                self._comparison_dict[comparison_key] = ChIPSeqComparison(
                    c_name=c_name,
                    t_name=t_name,
                    c_samples=c_sample_list,
                    t_samples=t_sample_list,
                    factor=factor,
                    tissue=tissue,
                    condition=condition,
                    treatment=treatment,
                    replicate=0,
                    diff_bind=diff_bind)

            # Sort the comparison keys alphabetically and assign replicate numbers into ChIPSeqComparison objects.

            for key_1, level_2_dict in level_1_dict.iteritems():
                i = 1
                for key_2 in sorted(level_2_dict):
                    if not self._comparison_dict[key_2].diff_bind:
                        continue
                    level_2_dict[key_2] = i
                    self._comparison_dict[key_2].replicate = i
                    i += 1

            return

        def run_create_bwa_jobs():
            """Private function to create BWA alignment jobs.

            @return:
            @rtype:
            """
            # Get the BWA index.

            bwa_genome_db = os.path.join(
                bsf.standards.FilePath.get_resource_genome_index(
                    genome_version=self.genome_version,
                    genome_index='bwa'),
                self.genome_version)

            self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

            for sample in self.sample_list:
                if self.debug > 0:
                    print(self, 'Sample name:', sample.name)
                    sys.stdout.writelines(sample.trace(level=1))

                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                for paired_reads_name in sorted(paired_reads_dict):
                    prefix_read_group = self.get_prefix_chipseq_alignment(paired_reads_name=paired_reads_name)

                    file_path_read_group = FilePathAlignment(prefix=prefix_read_group)

                    self.add_runnable(
                        runnable=bsf.Runnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            file_path_object=file_path_read_group,
                            debug=self.debug))

                    # Step 1: Process per lane.

                    bwa = bsf.executables.BWA(name='bwa_mem', analysis=self)

                    bwa_mem = bwa.sub_command

                    # Set BWA mem options.

                    # Allow as many threads as defined in the corresponding Stage object.
                    bwa_mem.add_option_short(key='t', value=str(stage_alignment.threads))
                    # Append FASTA/Q comment to SAM output.
                    bwa_mem.add_switch_short(key='C')
                    # Mark shorter split hits as secondary (for Picard compatibility).
                    bwa_mem.add_switch_short(key='M')
                    # Output warnings and errors only.
                    bwa_mem.add_option_short(key='v', value='2')

                    # Set BWA arguments.

                    bwa_mem.arguments.append(bwa_genome_db)

                    reads1 = list()
                    reads2 = list()

                    # Propagate the SAM read group information around FASTQ files if required.
                    # Please note that only the first read group can be propagated per
                    # PairedReads object.

                    read_group = str()

                    for paired_reads in paired_reads_dict[paired_reads_name]:
                        if paired_reads.reads_1 is not None:
                            reads1.append(paired_reads.reads_1.file_path)
                        if paired_reads.reads_2 is not None:
                            reads2.append(paired_reads.reads_2.file_path)
                        if not read_group and paired_reads.read_group:
                            read_group = paired_reads.read_group

                    if read_group:
                        bwa_mem.add_option_short(key='R', value=read_group)

                    if len(reads1) and not len(reads2):
                        bwa_mem.arguments.append(','.join(reads1))
                    elif len(reads1) and len(reads2):
                        bwa_mem.arguments.append(','.join(reads1))
                        bwa_mem.arguments.append(','.join(reads2))
                    if len(reads2):
                        warnings.warn('Only second reads, but no first reads have been defined.')

                    # Normally, the bwa object would be pushed onto the Stage list.
                    # Experimentally, use Pickler to serialize the Executable object into a file.

                    pickler_dict_align_lane = {
                        'prefix_read_group': stage_alignment.name,
                        'replicate_key': paired_reads_name,
                        'classpath_picard': self.classpath_picard,
                        'bwa_executable': bwa,
                    }

                    pickler_path = os.path.join(
                        self.genome_directory,
                        stage_alignment.name + '_' + paired_reads_name + '.pkl')
                    pickler_file = open(pickler_path, 'wb')
                    pickler = pickle.Pickler(file=pickler_file, protocol=pickle.HIGHEST_PROTOCOL)
                    pickler.dump(obj=pickler_dict_align_lane)
                    pickler_file.close()

                    # Create a bsf_run_bwa.py job to run the pickled object.

                    run_bwa = stage_alignment.add_executable(
                        executable=bsf.process.Executable(
                            name='_'.join((stage_alignment.name, paired_reads_name)),
                            program='bsf_run_bwa.py'))

                    # Only submit this Executable if the final result file does not exist.
                    if (os.path.exists(os.path.join(self.genome_directory, file_path_read_group.aligned_md5)) and
                            os.path.getsize(
                                os.path.join(self.genome_directory, file_path_read_group.aligned_md5))):
                        run_bwa.submit = False

                    # Set run_bwa options.
                    # Read Command options from configuration sections:
                    # [bsf.process.Executable.bsf_run_bwa.py]
                    self.set_command_configuration(command=run_bwa)
                    run_bwa.add_option_long(key='pickler_path', value=pickler_path)
                    run_bwa.add_option_long(key='debug', value=str(self.debug))

            return

        def run_create_bowtie2_jobs():
            """Create Bowtie2 alignment jobs.

            @return:
            @rtype:
            """

            # Get the Bowtie2 index.

            if self.genome_fasta_path.endswith('.fasta'):
                bowtie2_index = self.genome_fasta_path[:-6]
            elif self.genome_fasta_path.endswith('.fa'):
                bowtie2_index = self.genome_fasta_path[:-3]
            else:
                raise Exception(
                    "The 'genome_fasta' option has to end with '.fa' or '.fasta'. " + self.genome_fasta_path)

            # bowtie2_index = os.path.join(bsf.standards.Default.absolute_genomes(self.genome_version),
            #                              'forBowtie2', self.genome_version)

            self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

            for sample in self.sample_list:
                if self.debug > 0:
                    print(self, 'Sample name:', sample.name)
                    sys.stdout.writelines(sample.trace(level=1))

                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                for paired_reads_name in sorted(paired_reads_dict):
                    prefix_read_group = self.get_prefix_chipseq_alignment(paired_reads_name=paired_reads_name)

                    file_path_read_group = FilePathAlignment(prefix=prefix_read_group)

                    runnable_alignment = self.add_runnable(
                        runnable=bsf.Runnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            # FIXME: The bsf.runnables.bowtie2 module requires a lot of work.
                            # FIXME: For the moment, the code assumes that only FASTQ files will be used.
                            # code_module='bsf.runnables.bowtie2',
                            working_directory=self.genome_directory,
                            file_path_object=file_path_read_group,
                            debug=self.debug))
                    self.set_stage_runnable(
                        stage=stage_alignment,
                        runnable=runnable_alignment)
                    # Set dependencies.

                    runnable_alignment.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMakeDirectory(
                            name='make_directory',
                            directory_path=file_path_read_group.output_directory))

                    runnable_step = runnable_alignment.add_runnable_step(
                        runnable_step=bsf.process.RunnableStep(
                            name='bowtie2',
                            program='bowtie2'))

                    # Set Bowtie2 options.

                    runnable_step.add_option_short(key='x', value=bowtie2_index)
                    runnable_step.add_option_long(key='threads', value=str(stage_alignment.threads))

                    file_path_list_1 = list()
                    file_path_list_2 = list()

                    for paired_reads in paired_reads_dict[paired_reads_name]:
                        if paired_reads.reads_1 is not None:
                            file_path_list_1.append(paired_reads.reads_1.file_path)
                        if paired_reads.reads_2 is not None:
                            file_path_list_2.append(paired_reads.reads_2.file_path)

                    if len(file_path_list_1) and not len(file_path_list_2):
                        runnable_step.add_option_short(key='U', value=','.join(file_path_list_1))
                    elif len(file_path_list_1) and len(file_path_list_2):
                        runnable_step.add_option_short(key='1', value=','.join(file_path_list_1))
                    if len(file_path_list_2):
                        runnable_step.add_option_short(key='2', value=','.join(file_path_list_2))

                    # TODO: The following options are properties of the Sample,
                    # PairedReads and Reads objects.
                    # runnable_step.add_switch_short(key='q')
                    # runnable_step.add_switch_short(key='phred33')

                    # TODO: It would be good to have code that parses the Illumina Run Info and writes
                    # this information into the CSV file.

                    # TODO: Andreas' original implementation on Bowtie1 sets -a,
                    # which may not be required for Bowtie2.
                    # See description of option -k in the bowtie2 manual.
                    # TODO: To avoid repeats we may want to set -a?

                    # Set Bowtie2 arguments.

                    runnable_step.stdout_path = file_path_read_group.aligned_sam

                    # Run Picard CleanSam to convert the aligned SAM file into a cleaned SAM file.

                    runnable_step = runnable_alignment.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepPicard(
                            name='picard_clean_sam',
                            obsolete_file_path_list=[file_path_read_group.aligned_sam],
                            java_temporary_path=runnable_alignment.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx4G',
                            java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                            picard_command='CleanSam'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    runnable_step.add_picard_option(key='INPUT', value=file_path_read_group.aligned_sam)
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.cleaned_sam)
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_alignment.get_relative_temporary_directory_path)
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    runnable_step.add_picard_option(key='QUIET', value='false')
                    runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')

                    # Run Picard SortSam to convert the cleaned SAM file into a coordinate sorted BAM file.

                    runnable_step = runnable_alignment.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepPicard(
                            name='picard_sort_sam',
                            obsolete_file_path_list=[file_path_read_group.cleaned_sam],
                            java_temporary_path=runnable_alignment.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx6G',
                            java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                            picard_command='SortSam'))
                    """ @type runnable_step: bsf.process.RunnableStepPicard """
                    runnable_step.add_picard_option(key='INPUT', value=file_path_read_group.cleaned_sam)
                    runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.aligned_bam)
                    runnable_step.add_picard_option(key='SORT_ORDER', value='coordinate')
                    runnable_step.add_picard_option(
                        key='TMP_DIR',
                        value=runnable_alignment.get_relative_temporary_directory_path)
                    runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                    runnable_step.add_picard_option(key='QUIET', value='false')
                    runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')
                    runnable_step.add_picard_option(key='COMPRESSION_LEVEL', value='5')
                    runnable_step.add_picard_option(key='MAX_RECORDS_IN_RAM', value='4000000')
                    runnable_step.add_picard_option(key='CREATE_INDEX', value='true')
                    runnable_step.add_picard_option(key='CREATE_MD5_FILE', value='true')

                    # Create symbolic links from *.bai to *.bam.bai files.

                    runnable_alignment.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepLink(
                            name='link_bai',
                            source_path=file_path_read_group.aligned_bai_link_source,
                            target_path=file_path_read_group.aligned_bai_link_target))

            return

        def run_create_macs2_jobs():
            """Create MACS2 peak caller jobs.

            @return:
            @rtype:
            """

            for comparison_name in sorted(self._comparison_dict):
                chipseq_comparison = self._comparison_dict[comparison_name]
                factor = chipseq_comparison.factor.upper()
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                    for t_paired_reads_name in sorted(t_paired_reads_dict):
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)
                            for c_paired_reads_name in sorted(c_paired_reads_dict):
                                c_prefix_alignment = self.get_prefix_chipseq_alignment(
                                    paired_reads_name=c_paired_reads_name)
                                t_prefix_alignment = self.get_prefix_chipseq_alignment(
                                    paired_reads_name=t_paired_reads_name)
                                c_file_path_alignment = FilePathAlignment(prefix=c_prefix_alignment)
                                t_file_path_alignment = FilePathAlignment(prefix=t_prefix_alignment)

                                prefix_peak_calling = self.get_prefix_chipseq_peak_calling(
                                    t_paired_reads_name=t_paired_reads_name,
                                    c_paired_reads_name=c_paired_reads_name)

                                file_path_peak_calling = FilePathPeakCalling(prefix=prefix_peak_calling)

                                runnable_peak_calling = self.add_runnable(
                                    runnable=bsf.Runnable(
                                        name=prefix_peak_calling,
                                        code_module='bsf.runnables.generic',
                                        working_directory=self.genome_directory,
                                        file_path_object=file_path_peak_calling,
                                        debug=self.debug))
                                executable_peak_calling = self.set_stage_runnable(
                                    stage=stage_peak_calling,
                                    runnable=runnable_peak_calling)
                                executable_peak_calling.dependencies.append(c_prefix_alignment)
                                executable_peak_calling.dependencies.append(t_prefix_alignment)

                                # Add a RunnableStep to create the output directory.

                                runnable_peak_calling.add_runnable_step(
                                    runnable_step=bsf.process.RunnableStepMakeDirectory(
                                        name='make_directory',
                                        directory_path=file_path_peak_calling.output_directory))

                                # Add a RunnableStep for MACS2 call peak.

                                runnable_step_macs2_call_peak = runnable_peak_calling.add_runnable_step(
                                    runnable_step=bsf.process.RunnableStep(
                                        name='macs2_call_peak',
                                        program='macs2',
                                        sub_command=bsf.process.Command(program='callpeak')))
                                # Read RunnableStep options from configuration sections:
                                # [bsf.analyses.chip_seq.ChIPSeq.macs2_call_peak]
                                # [bsf.analyses.chip_seq.ChIPSeq.macs2_call_peak.callpeak]
                                self.set_runnable_step_configuration(runnable_step=runnable_step_macs2_call_peak)

                                runnable_step_macs2_call_peak.sub_command.add_option_long(
                                    key='treatment',
                                    value=t_file_path_alignment.aligned_bam)

                                runnable_step_macs2_call_peak.sub_command.add_option_long(
                                    key='control',
                                    value=c_file_path_alignment.aligned_bam)

                                # MACS2 can cope with directories specified in the --name option, but
                                # the resulting R script has them set too. Hence the R script has to be started
                                # from the genome_directory. However, the R script needs re-writing anyway, because
                                # it would be better to use the PNG rather than the PDF device for plotting.
                                runnable_step_macs2_call_peak.sub_command.add_option_long(
                                    key='name',
                                    value=file_path_peak_calling.name_prefix)

                                # FIXME: The 'gsize' option has to be specified via the configuration.ini file above.
                                # runnable_step_macs2_call_peak.add_option_long(
                                # key='gsize',
                                # value=self.genome_sizes_path)

                                runnable_step_macs2_call_peak.sub_command.add_switch_long(key='bdg')
                                runnable_step_macs2_call_peak.sub_command.add_switch_long(key='SPMR')

                                if factor == 'H3K4ME1':
                                    pass
                                elif factor == 'H3K4ME2':
                                    pass
                                elif factor == 'H3K4ME3':
                                    pass  # Default settings from above.
                                elif factor == 'H3K9ME3':
                                    pass
                                elif factor == 'H3K27AC':
                                    pass
                                elif factor == 'H3K27ME1':
                                    pass
                                elif factor == 'H3K27ME2':
                                    pass
                                elif factor == 'H3K27ME3':
                                    pass
                                elif factor == 'H3K36ME3':
                                    # Parameter setting for H3K36me3 according to Nature Protocols (2012)
                                    # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D)
                                    runnable_step_macs2_call_peak.sub_command.add_switch_long(key='nomodel')
                                    # The shiftsize option is no longer supported in MACS 2.1.0
                                    # mc2.add_option_long(key='shiftsize', value='73')
                                    runnable_step_macs2_call_peak.sub_command.add_option_long(
                                        key='pvalue',
                                        value='1e-3')
                                elif factor == 'H3K56AC':
                                    pass
                                elif factor == 'H4K16AC':
                                    pass
                                elif factor == 'OTHER':
                                    pass
                                else:
                                    warnings.warn(
                                        'Unable to set MACS2 parameters for unknown factor ' + repr(factor) + '.\n' +
                                        'Please use default factor ' + repr(self.factor_default) +
                                        ' or adjust Python code if necessary.',
                                        UserWarning)

                                # Add a RunnableStep to compare BedGraph files.

                                runnable_step_macs2_bdg_cmp = runnable_peak_calling.add_runnable_step(
                                    runnable_step=bsf.process.RunnableStep(
                                        name='macs2_bdg_cmp',
                                        program='macs2',
                                        sub_command=bsf.process.Command(program='bdgcmp')))
                                # Read RunnableStep options from configuration sections:
                                # [bsf.analyses.chip_seq.ChIPSeq.macs2_bdg_cmp]
                                # [bsf.analyses.chip_seq.ChIPSeq.macs2_bdg_cmp.bdgcmp]
                                self.set_runnable_step_configuration(runnable_step=runnable_step_macs2_bdg_cmp)

                                runnable_step_macs2_bdg_cmp.sub_command.add_option_long(
                                    key='tfile',
                                    value=file_path_peak_calling.treatment_bdg)

                                runnable_step_macs2_bdg_cmp.sub_command.add_option_long(
                                    key='cfile',
                                    value=file_path_peak_calling.control_bdg)

                                # Sequencing depth for treatment and control. Aim for setting the --SPMR parameter for
                                # runnable_step_macs2_call_peak to get the track normalised.
                                # --tdepth:
                                # --cdepth:
                                # --pseudocount

                                runnable_step_macs2_bdg_cmp.sub_command.add_option_long(
                                    key='ofile',
                                    value=file_path_peak_calling.comparison_bdg)

                                # --method defaults to ppois i.e. Poisson Pvalue (-log10(pvalue), which yields data
                                # on a logarithmic scale.

                                # runnable_step_macs2_bdg_cmp.sub_command.add_option_long(key='method', value='FE')

                                # Add a RunnableStep to process MACS2 output.

                                process_macs2 = runnable_peak_calling.add_runnable_step(
                                    runnable_step=bsf.process.RunnableStep(
                                        name='process_macs2',
                                        program='bsf_chipseq_process_macs2.bash'))

                                process_macs2.arguments.append(prefix_peak_calling)
                                process_macs2.arguments.append(self.genome_sizes_path)

                                if os.path.exists(file_path_peak_calling.narrow_peaks_bb):
                                    executable_peak_calling.submit = False

            return

        def run_create_diff_bind_jobs():
            """Create Bioconductor DiffBind jobs.

            @return:
            @rtype:
            """

            # Reorganise the ChIPSeqComparison objects by factor.

            for chipseq_comparison in self._comparison_dict.itervalues():
                if not chipseq_comparison.diff_bind:
                    continue
                if chipseq_comparison.factor not in self._factor_dict:
                    self._factor_dict[chipseq_comparison.factor] = list()
                self._factor_dict[chipseq_comparison.factor].append(chipseq_comparison)

            for factor_name in sorted(self._factor_dict):
                if self.debug > 0:
                    print('chipseq factor_name:', factor_name)

                factor_list = self._factor_dict[factor_name]
                factor_list.sort(cmp=lambda x, y: cmp(x, y))

                # No comparison for less than two items.

                if len(factor_list) < 2:
                    continue

                # Create a directory per factor.

                prefix_diff_bind = self.get_prefix_diff_bind(factor_name=factor_name)

                file_path_diff_bind = FilePathDiffBind(prefix=prefix_diff_bind)

                job_dependencies = list()

                # Create a new ChIPSeq DiffBind Sheet per factor.

                dbs = ChIPSeqDiffBindSheet(
                    file_path=os.path.join(
                        self.genome_directory,
                        file_path_diff_bind.sample_annotation_sheet))

                if self.debug > 0:
                    print('ChIPSeqDiffBindSheet file_path:', dbs.file_path)

                for chipseq_comparison in factor_list:
                    if not chipseq_comparison.diff_bind:
                        continue
                    for t_sample in chipseq_comparison.t_samples:
                        t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                        for t_paired_reads_name in sorted(t_paired_reads_dict):
                            for c_sample in chipseq_comparison.c_samples:
                                c_paired_reads_dict = c_sample.get_all_paired_reads(
                                    replicate_grouping=self.replicate_grouping)
                                for c_paired_reads_name in sorted(c_paired_reads_dict):
                                    # Get prefixes and FilePath objects for treatment and control.
                                    c_prefix_alignment = self.get_prefix_chipseq_alignment(
                                        paired_reads_name=c_paired_reads_name)
                                    t_prefix_alignment = self.get_prefix_chipseq_alignment(
                                        paired_reads_name=t_paired_reads_name)
                                    c_file_path_alignment = FilePathAlignment(prefix=c_prefix_alignment)
                                    t_file_path_alignment = FilePathAlignment(prefix=t_prefix_alignment)

                                    # Get prefix and FilePath object for the peak calls.
                                    prefix_peak_calling = self.get_prefix_chipseq_peak_calling(
                                        t_paired_reads_name=t_paired_reads_name,
                                        c_paired_reads_name=c_paired_reads_name)
                                    file_path_peak_calling = FilePathPeakCalling(prefix=prefix_peak_calling)

                                    job_dependencies.append(prefix_peak_calling)

                                    row_dict = dict()
                                    """ @type row_dict: dict[str, str | unicode] """

                                    row_dict['SampleID'] = t_paired_reads_name
                                    row_dict['Tissue'] = chipseq_comparison.tissue
                                    row_dict['Factor'] = chipseq_comparison.factor
                                    row_dict['Condition'] = chipseq_comparison.condition
                                    row_dict['Treatment'] = chipseq_comparison.treatment
                                    row_dict['Replicate'] = str(chipseq_comparison.replicate)
                                    row_dict['bamReads'] = t_file_path_alignment.aligned_bam
                                    row_dict['bamControl'] = c_file_path_alignment.aligned_bam
                                    row_dict['ControlID'] = c_paired_reads_name
                                    row_dict['Peaks'] = file_path_peak_calling.peaks_xls
                                    row_dict['PeakCaller'] = 'macs'
                                    row_dict['PeakFormat'] = 'macs'
                                    # row_dict['ScoreCol'] = str()
                                    # row_dict['LowerBetter'] = str()
                                    # row_dict['Counts'] = str()

                                    dbs.row_dicts.append(row_dict)

                dbs.to_file_path()

                # Create the DiffBind job.

                diff_bind = stage_diff_bind.add_executable(
                    executable=bsf.process.Executable(
                        name=prefix_diff_bind,
                        program='bsf_chipseq_diffbind.R'))
                diff_bind.dependencies.extend(job_dependencies)

                diff_bind.add_option_long(key='factor', value=factor_name)
                diff_bind.add_option_long(key='sample-annotation', value=file_path_diff_bind.sample_annotation_sheet)

                if os.path.exists(os.path.join(self.genome_directory, file_path_diff_bind.correlation_read_counts_png)):
                    diff_bind.submit = False

            return

        # Start of the run() method body.

        super(ChIPSeq, self).run()

        # ChIPSeq requires a genome version.

        if not self.genome_version:
            raise Exception('A ' + self.name + " requires a 'genome_version' configuration option.")

        if not self.comparison_path:
            raise Exception('A ' + self.name + " requires a 'cmp_file' configuration option.")
        self.comparison_path = self.configuration.get_absolute_path(file_path=self.comparison_path)

        if not self.colour_default:
            self.colour_default = '0,0,0'

        if self.colour_dict is None:
            raise Exception('A ' + self.name + " requires a TrackHubColours configuration section.")

        if not self.factor_default:
            self.factor_default = 'OTHER'

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        if not self.genome_fasta_path:
            self.genome_fasta_path = bsf.standards.FilePath.get_resource_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        # TODO: for the genome sizes path, the *.fa.fai file could be used.
        if self.genome_sizes_path:
            self.genome_sizes_path = self.configuration.get_absolute_path(file_path=self.genome_sizes_path)

        if not self.classpath_picard:
            self.classpath_picard = bsf.standards.JavaClassPath.get_picard()
            if not self.classpath_picard:
                raise Exception('A ' + self.name + " requires a 'classpath_picard' configuration option.")

        stage_alignment = self.get_stage(name=self.stage_name_alignment)
        stage_peak_calling = self.get_stage(name=self.stage_name_peak_calling)
        stage_diff_bind = self.get_stage(name=self.stage_name_diff_bind)

        run_read_comparisons()
        # run_create_bwa_jobs()
        run_create_bowtie2_jobs()
        # self._create_macs14_jobs()
        run_create_macs2_jobs()
        run_create_diff_bind_jobs()

        return

    def _create_macs14_jobs(self):
        """Create MACS14 peak caller jobs.

        @return:
        @rtype:
        """

        stage_macs14 = self.get_stage(name='macs14')
        stage_process_macs14 = self.get_stage(name='process_macs14')

        for comparison_name in sorted(self._comparison_dict):
            chipseq_comparison = self._comparison_dict[comparison_name]
            factor = chipseq_comparison.factor.upper()
            for t_sample in chipseq_comparison.t_samples:
                t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                for t_paired_reads_name in sorted(t_paired_reads_dict):
                    for c_sample in chipseq_comparison.c_samples:
                        c_paired_reads_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                        for c_paired_reads_name in sorted(c_paired_reads_dict):
                            macs14 = stage_macs14.add_executable(
                                executable=bsf.process.Executable(
                                    name='chipseq_macs14_' + t_paired_reads_name + '__' + c_paired_reads_name,
                                    program='macs14'))

                            macs14.dependencies.append('chipseq_sam2bam_' + t_paired_reads_name)
                            macs14.dependencies.append('chipseq_sam2bam_' + c_paired_reads_name)

                            process_macs14 = stage_process_macs14.add_executable(
                                executable=bsf.process.Executable(
                                    name='chipseq_process_macs14_' + t_paired_reads_name + '__' + c_paired_reads_name,
                                    program='bsf_chipseq_process_macs14.bash'))
                            process_macs14.dependencies.append(macs14.name)

                            # Read Command options from configuration sections:
                            # [bsf.process.Executable.macs14]
                            self.set_command_configuration(command=macs14)

                            # Set macs14 options.
                            macs14.add_option_long(
                                key='treatment',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_' + t_paired_reads_name,
                                    t_paired_reads_name + '.bam'))
                            macs14.add_option_long(
                                key='control',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_' + c_paired_reads_name,
                                    c_paired_reads_name + '.bam'))

                            # TODO: Experimentally prepend a chipseq_macs14 directory
                            # MACS14 can hopefully also cope with directories specified in the --name option, but
                            # the resulting R script has them set too. Hence the R script has to be started
                            # from the genome_directory. However, the R script needs re-writing anyway, because
                            # it would be better to use the PNG rather than the PDF device for plotting.
                            prefix = 'chipseq_macs14_' + t_paired_reads_name + '__' + c_paired_reads_name

                            replicate_directory = os.path.join(self.genome_directory, prefix)

                            try:
                                os.makedirs(replicate_directory)
                            except OSError as exc:  # Python > 2.5
                                if exc.errno == errno.EEXIST and os.path.isdir(replicate_directory):
                                    pass
                                else:
                                    raise

                            macs14.add_option_long(key='name', value=os.path.join('.', prefix, prefix))

                            # This option has to be specified via the configuration.ini file.
                            # macs14.add_option_long(key='gsize', value=self.genome_sizes_path)
                            macs14.add_switch_long(key='single-profile')
                            macs14.add_switch_long(key='call-subpeaks')
                            macs14.add_switch_long(key='wig')

                            if factor == 'H3K4ME1':
                                pass
                            elif factor == 'H3K4ME2':
                                pass
                            elif factor == 'H3K4ME3':
                                pass  # Default settings from above.
                            elif factor == 'H3K9ME3':
                                pass
                            elif factor == 'H3K27AC':
                                pass
                            elif factor == 'H3K27ME1':
                                pass
                            elif factor == 'H3K27ME2':
                                pass
                            elif factor == 'H3K27ME3':
                                pass
                            elif factor == 'H3K36ME3':
                                # Parameter setting for H3K36me3 according to Nature Protocols (2012)
                                # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D)
                                macs14.add_switch_long(key='nomodel')
                                macs14.add_option_long(key='shiftsize', value='73')
                                macs14.add_option_long(key='pvalue', value='1e-3')
                            elif factor == 'H3K56AC':
                                pass
                            elif factor == 'H4K16AC':
                                pass
                            elif factor == 'OTHER':
                                pass
                            else:
                                warnings.warn(
                                    'Unable to set MACS14 parameters for unknown factor ' + repr(factor) + '.\n' +
                                    'Please use default factor ' + repr(self.factor_default) +
                                    ' or adjust Python code if necessary.',
                                    UserWarning)

                            # Set macs14 arguments.

                            # Set process_macs14 options.
                            # Set process_macs14 arguments.

                            # Specify the output path as in the macs14 --name option.
                            process_macs14.arguments.append(os.path.join('.', prefix, prefix))
                            process_macs14.arguments.append(self.genome_sizes_path)

        return

    def _report_macs14(self):
        """Create a C{bsf.analyses.chip_seq.ChIPSeq} report in HTML format and a UCSC Genome Browser Track Hub.

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

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str | unicode] """

            str_list.append('<h1 id="chip_seq_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.append('<p>\n')
            str_list.append('Next-Generation Sequencing reads are aligned with the short read aligner\n')
            str_list.append('<strong>')
            str_list.append('<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>')
            str_list.append('</strong>,\n')
            str_list.append('before peaks are called with ')
            str_list.append('<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 1.4</a>\n')
            str_list.append('on a treatment and control sample pair.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Construct an automatic UCSC Track Hub link.

            str_list.append('<p id="ucsc_track_hub">\n')
            str_list.append('View Bowtie2 <strong>read alignment</strong> tracks for each sample\n')
            str_list.append('in their genomic context via the project-specific\n')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('.')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="peak_calling_table">\n')
            str_list.append('<thead>\n')

            str_list.append('<tr>\n')
            str_list.append('<th>Peaks</th>\n')
            str_list.append('<th>Negative Peaks</th>\n')
            str_list.append('<th>R Model</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in sorted(self._comparison_dict):
                chipseq_comparison = self._comparison_dict[comparison_name]
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping)
                    for t_paired_reads_name in sorted(t_paired_reads_dict):
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)
                            for c_paired_reads_name in sorted(c_paired_reads_dict):
                                # prefix = 'chipseq_macs14_' + t_paired_reads_name + '__' + c_paired_reads_name
                                str_list.append('<tr>\n')
                                # if treatment and absolute:
                                # str_list.append('<td><strong>' + t_paired_reads_name + '</strong></td>\n')
                                # if not treatment and absolute:
                                # str_list.append('<td><strong>' + c_paired_reads_name + '</strong></td>\n')

                                # Peaks
                                str_list.append('<td>')
                                str_list.append('<a href="' + t_paired_reads_name + '__' + c_paired_reads_name +
                                                '_peaks.xls">')
                                str_list.append('Peaks ' + t_paired_reads_name + ' versus ' + c_paired_reads_name)
                                str_list.append('</a>')
                                str_list.append('</td>\n')
                                # Negative Peaks
                                str_list.append('<td>')
                                str_list.append('<a href="' + t_paired_reads_name + '__' + c_paired_reads_name +
                                                '_negative_peaks.xls">')
                                str_list.append('Negative peaks ' + t_paired_reads_name + ' versus ' +
                                                c_paired_reads_name)
                                str_list.append('</a>')
                                str_list.append('</td>\n')
                                # R Model
                                str_list.append('<td>')
                                str_list.append('<a href="' + t_paired_reads_name + '__' + c_paired_reads_name +
                                                '_model.r">')
                                str_list.append('R model')
                                str_list.append('</a>')
                                str_list.append('</td>\n')

                                str_list.append('</tr>\n')
                                str_list.append('\n')

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

            for comparison_name in sorted(self._comparison_dict):
                chipseq_comparison = self._comparison_dict[comparison_name]
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                    for t_paired_reads_name in sorted(t_paired_reads_dict):
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)
                            for c_paired_reads_name in sorted(c_paired_reads_dict):
                                # prefix = 'chipseq_macs14_' + t_paired_reads_name + '__' + c_paired_reads_name
                                prefix = t_paired_reads_name + '__' + c_paired_reads_name

                                # Add UCSC trackDB entries for each treatment/control and absolute/normalised pair.

                                for treatment in [True, False]:
                                    if treatment:
                                        state = 'treat'
                                    else:
                                        state = 'control'

                                    for absolute in [True, False]:
                                        if absolute:
                                            scaling = 'absolute'
                                        else:
                                            scaling = 'normalised'

                                        #
                                        # Add a UCSC trackDB entry for each bigWig file
                                        #
                                        # Common settings
                                        str_list.append('track ChIP_' + '_'.join((prefix, state, scaling)) + '\n')
                                        str_list.append('type bigWig\n')
                                        str_list.append('shortLabel ChIP_' + '_'.join((prefix, state, scaling)) + '\n')
                                        str_list.append('longLabel ' + scaling.capitalize() +
                                                        ' ChIP-Seq read counts for ' + state + ' of ' +
                                                        t_paired_reads_name + ' versus ' + c_paired_reads_name + '\n')
                                        str_list.append('bigDataUrl ')
                                        str_list.append('/'.join((
                                            prefix + '_MACS_wiggle', state,
                                            '_'.join((prefix, state, 'afterfiting', 'all.bw')))) + '\n')
                                        # str_list.append('html ...\n'
                                        if treatment and not absolute:
                                            str_list.append('visibility full\n')
                                        else:
                                            str_list.append('visibility hide\n')

                                        # Common optional settings
                                        str_list.append(
                                            'color ' + self.get_colour(factor=chipseq_comparison.factor) + '\n')

                                        # bigWig - Signal graphing track settings
                                        str_list.append('graphTypeDefault bar\n')
                                        str_list.append('maxHeightPixels 100:60:20\n')
                                        str_list.append('smoothingWindow off\n')
                                        if absolute:
                                            # Track with absolute scaling.
                                            str_list.append('autoScale on\n')
                                        else:
                                            # Track with relative scaling.
                                            str_list.append('autoScale off\n')
                                            str_list.append('viewLimits 0:40\n')

                                        str_list.append('\n')

                                #
                                # Add a UCSC trackDB entry for each bigBed peaks file
                                #
                                # Common settings
                                str_list.append('track Peaks_' + prefix + '\n')
                                str_list.append('type bigBed\n')
                                str_list.append('shortLabel Peaks_' + prefix + '\n')
                                str_list.append('longLabel ChIP-Seq peaks for ' + t_paired_reads_name +
                                                ' versus ' + c_paired_reads_name + '\n')
                                str_list.append('bigDataUrl ' + prefix + '_peaks.bb\n')
                                # str_list.append('html ...\n')
                                str_list.append('visibility pack\n')

                                # Common optional settings
                                str_list.append('color ' + self.get_colour(factor=chipseq_comparison.factor) + '\n')

                                str_list.append('\n')

            # Add UCSC trackDB entries for each Bowtie2 BAM file.

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                for paired_reads_name in sorted(paired_reads_dict):
                    #
                    # Add a UCSC trackDB entry.
                    #
                    # Common settings
                    str_list.append('track Alignment_' + paired_reads_name + '\n')
                    str_list.append('type bam\n')
                    str_list.append('shortLabel Alignment_' + paired_reads_name + '\n')
                    str_list.append('longLabel Bowtie2 alignment of ' + paired_reads_name + '\n')
                    str_list.append('bigDataUrl ' + paired_reads_name + '.bam\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')

                    # Common optional settings
                    # str_list.append('color ' +
                    #                 bsf.standards.Defaults.web.chipseq_colours[comparison[2].upper()] + '\n')

                    # bam - Compressed Sequence Alignment track settings.

                    str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return

    def report(self):
        """Create a C{bsf.analyses.chip_seq.ChIPSeq} report in HTML format and a UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # contrast_field_names = ['', 'Group1', 'Members1', 'Group2', 'Members2', 'DB.edgeR']

        def report_html():
            """Private function to create a HTML report.

            @return:
            @rtype:
            """
            # Create a symbolic link containing the project name and a UUID.
            link_path = self.create_public_project_link()

            # Write a HTML document.

            str_list = list()
            """ @type str_list: list[str | unicode] """

            str_list.append('<h1 id="chip_seq_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
            str_list.append('\n')

            str_list.append('<p>\n')
            str_list.append('Next-Generation Sequencing reads are aligned with the short read aligner\n')
            str_list.append('<strong>')
            str_list.append('<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>')
            str_list.append('</strong>,\n')
            str_list.append('before peaks are called with ')
            str_list.append('<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 2</a>\n')
            str_list.append('on an enrichment and background sample pair.\n')
            str_list.append('</p>\n')
            str_list.append('\n')

            # Construct an automatic UCSC Track Hub link.

            str_list.append('<p id="ucsc_track_hub">\n')
            str_list.append('View Bowtie2 <strong>read alignment</strong> tracks for each sample\n')
            str_list.append('in their genomic context via the project-specific\n')
            str_list.extend(self.ucsc_hub_html_anchor(link_path=link_path))
            str_list.append('.')
            str_list.append('</p>\n')
            str_list.append('\n')

            str_list.append('<table id="peak_calling_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Comparison</th>\n')
            str_list.append('<th>Peaks</th>\n')
            str_list.append('<th>R Model</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for comparison_name in sorted(self._comparison_dict):
                chipseq_comparison = self._comparison_dict[comparison_name]
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping)
                    for t_paired_reads_name in sorted(t_paired_reads_dict):
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)
                            for c_paired_reads_name in sorted(c_paired_reads_dict):
                                runnable_peak_calling = self.runnable_dict[self.get_prefix_chipseq_peak_calling(
                                    t_paired_reads_name=t_paired_reads_name,
                                    c_paired_reads_name=c_paired_reads_name)]
                                file_path_peak_calling = runnable_peak_calling.file_path_object
                                """ @type file_path_peak_calling: FilePathPeakCalling """

                                str_list.append('<tr>\n')
                                # Comparison
                                str_list.append('<td>')
                                str_list.append('<strong>' + t_paired_reads_name + '</strong>')
                                str_list.append(' versus ')
                                str_list.append('<strong>' + c_paired_reads_name + '</strong>')
                                str_list.append('</td>\n')
                                # Peaks
                                str_list.append('<td><a href="' + file_path_peak_calling.peaks_xls + '">')
                                str_list.append('<abbr title="Tab-Separated Value">XLS</abbr>')
                                str_list.append('</a></td>\n')
                                # R Model
                                str_list.append('<td><a href="' + file_path_peak_calling.model_r + '">')
                                str_list.append('R model')
                                str_list.append('</a></td>\n')
                                str_list.append('</tr>\n')
                                str_list.append('\n')

            str_list.append('</tbody>\n')
            str_list.append('</table>\n')
            str_list.append('\n')

            # Differential binding analysis.

            str_list.append('<h2 id="differential_binding">Differential Binding Analysis</h2>\n')
            str_list.append('\n')

            str_list.append('<table id="differential_binding_table">\n')
            str_list.append('<thead>\n')
            str_list.append('<tr>\n')
            str_list.append('<th>Factor and Contrast</th>\n')
            str_list.append('<th>Correlation Peak Caller</th>\n')
            str_list.append('<th>Correlation Peak Counts</th>\n')
            str_list.append('<th>Correlation Analysis</th>\n')
            str_list.append('<th>MA Plot</th>\n')
            str_list.append('<th>Scatter Plot</th>\n')
            str_list.append('<th>PCA Plot</th>\n')
            str_list.append('<th>Box Plot</th>\n')
            str_list.append('<th>Differntially Bound Sites</th>\n')
            str_list.append('</tr>\n')
            str_list.append('</thead>\n')
            str_list.append('<tbody>\n')

            for factor_name in sorted(self._factor_dict):
                prefix_diff_bind = self.get_prefix_diff_bind(factor_name=factor_name)
                # runnable_diff_bind = self.runnable_dict[prefix_diff_bind]
                # file_path_diff_bind = runnable_diff_bind.file_path_object
                # """ @type file_path_diff_bind: FilePathDiffBind """
                # TODO: For the moment, DiffBind is run via an Executable, not a Runnable.
                file_path_diff_bind = FilePathDiffBind(prefix=prefix_diff_bind)

                str_list.append('<tr>\n')

                # ChIP Factor
                str_list.append('<td><strong>' + factor_name + '</strong></td>\n')
                # Correlation heat map of peak caller scores
                str_list.append('<td>')
                str_list.append('<a href="' + file_path_diff_bind.correlation_peak_caller_score_png + '">')
                str_list.append('<img')
                str_list.append(' alt="DiffBind correlation analysis for factor ' + factor_name + '"')
                str_list.append(' src="' + file_path_diff_bind.correlation_peak_caller_score_png + '"')
                str_list.append(' height="80" width="80">')
                str_list.append('</a>')
                str_list.append('</td>\n')

                # Correlation heat map of counts
                str_list.append('<td>')
                str_list.append('<a href="' + file_path_diff_bind.correlation_read_counts_png + '">')
                str_list.append('<img')
                str_list.append(' alt="DiffBind correlation analysis for factor ' + factor_name + '"')
                str_list.append(' src="' + file_path_diff_bind.correlation_read_counts_png + '"')
                str_list.append(' height="80" width="80">')
                str_list.append('</a>')
                str_list.append('</td>\n')

                # Correlation heat map of differential binding analysis
                str_list.append('<td>')
                str_list.append('<a href="' + file_path_diff_bind.correlation_analysis_png + '">')
                str_list.append('<img')
                str_list.append(' alt="DiffBind correlation analysis for factor ' + factor_name + '"')
                str_list.append(' src="' + file_path_diff_bind.correlation_analysis_png + '"')
                str_list.append(' height="80" width="80">')
                str_list.append('</a>')
                str_list.append('</td>\n')

                str_list.append('<td></td>\n')  # MA Plot
                str_list.append('<td></td>\n')  # Scatter Plot

                # Score-based PCA Plot
                str_list.append('<td>')
                str_list.append('<a href="' + file_path_diff_bind.pca_png + '">')
                str_list.append('<img')
                str_list.append(' alt="Score-based PCA plot for factor ' + factor_name + '"')
                str_list.append(' src="' + file_path_diff_bind.pca_png + '"')
                str_list.append(' height="80" width="80">')
                str_list.append('</a>')
                str_list.append('</td>\n')

                str_list.append('<td></td>\n')  # Box Plot
                str_list.append('<td></td>\n')  # DiffBind Report

                str_list.append('</tr>\n')

                # Read the file of contrasts ...

                file_path = os.path.join(self.genome_directory, file_path_diff_bind.contrasts_csv)

                if not os.path.exists(path=file_path):
                    warnings.warn(
                        'Contrasts table does not exist: ' + file_path,
                        UserWarning)
                    continue

                annotation_sheet = bsf.annotation.AnnotationSheet.from_file_path(
                    file_path=file_path,
                    file_type='excel')

                for row_dict in annotation_sheet.row_dicts:
                    file_path_diff_bind_contrast = FilePathDiffBindContrast(
                        prefix=prefix_diff_bind,
                        group_1=row_dict['Group1'],
                        group_2=row_dict['Group2'])
                    # suffix = '__'.join((row_dict['Group1'], row_dict['Group2']))

                    str_list.append('<tr>\n')

                    str_list.append('<td>' + row_dict['Group1'] + ' vs ' + row_dict['Group2'] + '</td>\n')
                    str_list.append('<td></td>\n')  # Correlation heat map of peak caller scores
                    str_list.append('<td></td>\n')  # Correlation heat map of counts
                    str_list.append('<td></td>\n')  # Correlation heat map of differential binding analysis

                    # MA Plot
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_diff_bind_contrast.ma_plot_png + '">')
                    str_list.append('<img')
                    str_list.append(' alt="DiffBind MA plot for factor ' + factor_name + '"')
                    str_list.append(' src="' + file_path_diff_bind_contrast.ma_plot_png + '"')
                    str_list.append(' height="80" width="80">')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Scatter Plot
                    str_list.append('<td>')
                    str_list.append('<a href="' + file_path_diff_bind_contrast.scatter_plot_png + '">')
                    str_list.append('<img')
                    str_list.append(' alt="DiffBind scatter plot for factor ' + factor_name + '"')
                    str_list.append(' src="' + file_path_diff_bind_contrast.scatter_plot_png + '"')
                    str_list.append(' height="80" width="80">')
                    str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Principal Component Analysis Plot (optional)
                    str_list.append('<td>')
                    if os.path.exists(os.path.join(self.genome_directory, file_path_diff_bind_contrast.pca_plot_png)):
                        str_list.append('<a href="' + file_path_diff_bind_contrast.pca_plot_png + '">')
                        str_list.append('<img')
                        str_list.append(' alt="DiffBind PCA plot for factor ' + factor_name + '"')
                        str_list.append(' src="' + file_path_diff_bind_contrast.pca_plot_png + '"')
                        str_list.append(' height="80" width="80">')
                        str_list.append('</a>')
                    str_list.append('</td>\n')

                    # Box Plot (optional)
                    str_list.append('<td>')
                    if os.path.exists(os.path.join(self.genome_directory, file_path_diff_bind_contrast.box_plot_png)):
                        str_list.append('<a href="' + file_path_diff_bind_contrast.box_plot_png + '">')
                        str_list.append('<img')
                        str_list.append(' alt="DiffBind Box plot for factor ' + factor_name + '"')
                        str_list.append(' src="' + file_path_diff_bind_contrast.box_plot_png + '"')
                        str_list.append(' height="80" width="80">')
                        str_list.append('</a>')
                    str_list.append('</td>\n')

                    # DiffBind report (optional)
                    str_list.append('<td>')
                    if os.path.exists(os.path.join(self.genome_directory, file_path_diff_bind_contrast.dba_report_csv)):
                        str_list.append('<a href="' + file_path_diff_bind_contrast.dba_report_csv + '">')
                        str_list.append('<abbr title="Comma-Separated Value">CSV</abbr>')
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

            @return:
            @rtype:
            """
            str_list = list()
            """ @type str_list: list[str | unicode] """

            # Group via UCSC super tracks.

            str_list.append('track Alignment\n')
            str_list.append('shortLabel QC Alignment\n')
            str_list.append('longLabel ChIP Read Alignment\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group alignment\n')
            str_list.append('\n')

            str_list.append('track Background\n')
            str_list.append('shortLabel QC Background\n')
            str_list.append('longLabel ChIP Background Signal\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group background\n')
            str_list.append('\n')

            str_list.append('track Enrichment\n')
            str_list.append('shortLabel QC Enrichment\n')
            str_list.append('longLabel ChIP Enrichment Signal\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group enrichment\n')
            str_list.append('\n')

            str_list.append('track Comparison\n')
            str_list.append('shortLabel ChIP Intensity\n')
            str_list.append('longLabel ChIP Intensity\n')
            str_list.append('visibility full\n')
            str_list.append('superTrack on\n')
            str_list.append('group comparison\n')
            str_list.append('\n')

            str_list.append('track Peaks\n')
            str_list.append('shortLabel ChIP Peaks\n')
            str_list.append('longLabel ChIP Peaks\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group peaks\n')
            str_list.append('\n')

            str_list.append('track Summits\n')
            str_list.append('shortLabel ChIP Summits\n')
            str_list.append('longLabel ChIP Summits\n')
            str_list.append('visibility hide\n')
            str_list.append('superTrack on\n')
            str_list.append('group summits\n')
            str_list.append('\n')

            # Group via UCSC composite tracks.

            composite_groups = dict()
            """ @type composite_groups: dict[str, bool] """

            for comparison_name in sorted(self._comparison_dict):
                chipseq_comparison = self._comparison_dict[comparison_name]
                factor_name = chipseq_comparison.factor.upper()
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                    for t_paired_reads_name in sorted(t_paired_reads_dict):
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)
                            for c_paired_reads_name in sorted(c_paired_reads_dict):
                                runnable_peak_calling = self.runnable_dict[self.get_prefix_chipseq_peak_calling(
                                    t_paired_reads_name=t_paired_reads_name,
                                    c_paired_reads_name=c_paired_reads_name)]
                                file_path_peak_calling = runnable_peak_calling.file_path_object
                                """ @type file_path_peak_calling: FilePathPeakCalling """

                                prefix_short = '__'.join((t_paired_reads_name, c_paired_reads_name))
                                prefix_long = t_paired_reads_name + ' versus ' + c_paired_reads_name

                                # Add UCSC trackDB entries for each treatment/control and absolute/normalised pair.
                                # NAME_control_lambda.bw
                                # NAME_treat_pileup.bw
                                # NAME_bdgcmp.bw

                                #
                                # Add a factor-specific "background" composite track once
                                #
                                composite_group = '_'.join((factor_name, 'background'))
                                if composite_group not in composite_groups:
                                    composite_groups[composite_group] = True
                                    str_list.append('track ' + composite_group + '\n')
                                    str_list.append('type bigWig\n')
                                    str_list.append('shortLabel ' + composite_group + '\n')
                                    str_list.append("longLabel ChIP background signal for factor '" + factor_name +
                                                    "'\n")
                                    str_list.append('visibility dense\n')
                                    str_list.append('compositeTrack on\n')
                                    str_list.append('parent Background\n')
                                    str_list.append('allButtonPair on\n')
                                    str_list.append('centerLabelsDense on\n')
                                    str_list.append('\n')

                                #
                                # Add a UCSC trackDB entry for each NAME_control_lambda.bw file.
                                #
                                # Common settings
                                str_list.append('track ' + prefix_short + '_background\n')
                                # TODO: The bigWig type must declare the expected signal range.
                                # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                                str_list.append('type bigWig\n')
                                str_list.append('shortLabel ' + prefix_short + '_bkg\n')
                                str_list.append('longLabel ChIP background signal ' + prefix_long + '\n')
                                str_list.append('bigDataUrl ' + file_path_peak_calling.control_bw + '\n')
                                # str_list.append('html ...\n')
                                str_list.append('visibility dense\n')

                                # Common optional settings
                                str_list.append('color ' + self.get_colour(factor=factor_name) + '\n')

                                # bigWig - Signal graphing track settings
                                str_list.append('alwaysZero off\n')
                                str_list.append('autoScale off\n')
                                str_list.append('graphTypeDefault bar\n')
                                str_list.append('maxHeightPixels 100:60:20\n')
                                # str_list.append('maxWindowToQuery 10000000\n')
                                str_list.append('smoothingWindow 5\n')
                                # str_list.append('transformFunc NONE\n')
                                str_list.append('viewLimits 0:15\n')
                                str_list.append('viewLimitsMax 0:40\n')
                                str_list.append('windowingFunction maximum\n')
                                # str_list.append('yLineMark <#>\n')
                                # str_list.append('yLineOnOff on \n')
                                # str_list.append('gridDefault on\n')

                                # Composite track settings
                                str_list.append('parent ' + composite_group + ' on\n')
                                str_list.append('centerLabelsDense off\n')
                                str_list.append('\n')

                                #
                                # Add a factor-specific "enrichment" composite track once
                                #
                                composite_group = '_'.join((factor_name, 'enrichment'))
                                if composite_group not in composite_groups:
                                    composite_groups[composite_group] = True
                                    str_list.append('track ' + composite_group + '\n')
                                    str_list.append('type bigWig\n')
                                    str_list.append('shortLabel ' + factor_name + '_enrichment\n')
                                    str_list.append("longLabel ChIP enrichment signal for factor '" + factor_name +
                                                    "'\n")
                                    str_list.append('visibility dense\n')
                                    str_list.append('compositeTrack on\n')
                                    str_list.append('parent Enrichment\n')
                                    str_list.append('allButtonPair on\n')
                                    str_list.append('centerLabelsDense on\n')
                                    str_list.append('\n')

                                #
                                # Add a UCSC trackDB entry for each NAME_treat_pileup.bw file.
                                #
                                # Common settings
                                str_list.append('track ' + prefix_short + '_enrichment\n')
                                # TODO: The bigWig type must declare the expected signal range.
                                str_list.append('type bigWig\n')
                                str_list.append('shortLabel ' + prefix_short + '_enr\n')
                                str_list.append('longLabel ChIP enrichment signal ' + prefix_long + '\n')
                                str_list.append('bigDataUrl ' + file_path_peak_calling.treatment_bw + '\n')
                                # str_list.append('html ...\n')
                                str_list.append('visibility dense\n')

                                # Common optional settings
                                str_list.append('color ' + self.get_colour(factor=factor_name) + '\n')

                                # bigWig - Signal graphing track settings
                                str_list.append('alwaysZero off\n')
                                str_list.append('autoScale off\n')
                                str_list.append('graphTypeDefault bar\n')
                                str_list.append('maxHeightPixels 100:60:20\n')
                                # str_list.append('maxWindowToQuery 10000000\n')
                                str_list.append('smoothingWindow 5\n')
                                # str_list.append('transformFunc NONE\n')
                                str_list.append('viewLimits 0:15\n')
                                str_list.append('viewLimitsMax 0:40\n')
                                str_list.append('windowingFunction maximum\n')
                                # str_list.append('yLineMark <#>\n')
                                # str_list.append('yLineOnOff on \n')
                                # str_list.append('gridDefault on\n')

                                # Composite track settings
                                str_list.append('parent ' + composite_group + ' on\n')
                                str_list.append('centerLabelsDense off\n')
                                str_list.append('\n')

                                #
                                # Add a factor-specific "intensity" composite track once
                                #
                                composite_group = '_'.join((factor_name, 'intensity'))
                                if composite_group not in composite_groups:
                                    composite_groups[composite_group] = True
                                    str_list.append('track ' + composite_group + '\n')
                                    str_list.append('type bigWig\n')
                                    str_list.append('shortLabel ' + composite_group + '\n')
                                    str_list.append("longLabel ChIP intensity for factor '" + factor_name + "'\n")
                                    str_list.append('visibility full\n')
                                    str_list.append('compositeTrack on\n')
                                    str_list.append('parent Comparison\n')
                                    str_list.append('allButtonPair on\n')
                                    str_list.append('centerLabelsDense on\n')
                                    str_list.append('\n')

                                #
                                # Add a UCSC trackDB entry for each NAME_bdgcmp.bw file.
                                #
                                # Common settings
                                str_list.append('track ' + prefix_short + '_intensity\n')
                                # TODO: The bigWig type must declare the expected signal range.
                                str_list.append('type bigWig\n')
                                str_list.append('shortLabel ' + prefix_short + '_int\n')
                                str_list.append('longLabel ChIP intensity ' + prefix_long + '\n')
                                str_list.append('bigDataUrl ' + file_path_peak_calling.comparison_bw + '\n')
                                # str_list.append('html ...\n')
                                str_list.append('visibility full\n')

                                # Common optional settings
                                str_list.append('color ' + self.get_colour(factor=factor_name) + '\n')

                                # bigWig - Signal graphing track settings
                                str_list.append('alwaysZero off\n')
                                str_list.append('autoScale off\n')
                                str_list.append('graphTypeDefault bar\n')
                                str_list.append('maxHeightPixels 100:60:20\n')
                                # str_list.append('maxWindowToQuery 10000000\n')
                                str_list.append('smoothingWindow 5\n')
                                # str_list.append('transformFunc NONE\n')
                                str_list.append('viewLimits 0:15\n')
                                str_list.append('viewLimitsMax 0:40\n')
                                str_list.append('windowingFunction maximum\n')
                                # str_list.append('yLineMark <#>\n')
                                # str_list.append('yLineOnOff on \n')
                                # str_list.append('gridDefault on\n')

                                # Composite track settings
                                str_list.append('parent ' + composite_group + ' on\n')
                                str_list.append('centerLabelsDense off\n')
                                str_list.append('\n')

                                #
                                # Add a UCSC trackDB entry for each NAME_peaks.bb file.
                                #
                                # Common settings
                                str_list.append('track ' + prefix_short + '_peaks\n')
                                str_list.append('type bigBed\n')
                                str_list.append('shortLabel ' + prefix_short + '_peaks\n')
                                str_list.append('longLabel ChIP peaks ' + prefix_long + '\n')
                                str_list.append('bigDataUrl ' + file_path_peak_calling.narrow_peaks_bb + '\n')
                                # str_list.append('html ...\n')
                                str_list.append('visibility pack\n')

                                # Common optional settings
                                str_list.append('color ' + self.get_colour(factor=factor_name) + '\n')

                                # bigBed - Item or region track settings.

                                # Supertrack settings
                                str_list.append('parent Peaks\n')
                                str_list.append('\n')

                                #
                                # Add a UCSC trackDB entry for each NAME_summits.bb file.
                                #
                                # Common settings
                                str_list.append('track ' + prefix_short + '_summits\n')
                                str_list.append('type bigBed\n')
                                str_list.append('shortLabel ' + prefix_short + '_summits\n')
                                str_list.append('longLabel ChIP summits ' + prefix_long + '\n')
                                str_list.append('bigDataUrl ' + file_path_peak_calling.summits_bb + '\n')
                                # str_list.append('html ...\n')
                                str_list.append('visibility pack\n')

                                # Common optional settings
                                str_list.append('color ' + self.get_colour(factor=factor_name) + '\n')

                                # Supertrack settings
                                str_list.append('parent Summits\n')
                                str_list.append('\n')

            # Add UCSC trackDB entries for each Bowtie2 BAM file.

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)
                for paired_reads_name in sorted(paired_reads_dict):
                    prefix_alignment = self.get_prefix_chipseq_alignment(paired_reads_name=paired_reads_name)
                    runnable_alignment = self.runnable_dict[prefix_alignment]
                    file_path_alignment = runnable_alignment.file_path_object
                    """ @type file_path_alignment: FilePathAlignment """
                    #
                    # Add a UCSC trackDB entry for each NAME.bam file.
                    #
                    # Common settings
                    str_list.append('track ' + paired_reads_name + '_alignment\n')
                    str_list.append('type bam\n')
                    str_list.append('shortLabel ' + paired_reads_name + '_alignment\n')
                    str_list.append('longLabel ' + paired_reads_name + ' ChIP read alignment\n')
                    str_list.append('bigDataUrl ' + file_path_alignment.aligned_bam + '\n')
                    # str_list.append('html ...\n')
                    str_list.append('visibility hide\n')

                    # Common optional settings
                    # str_list.append('color ' +
                    #                 bsf.standards.Defaults.web.chipseq_colours[comparison[2].upper()] + '\n')

                    # bam - Compressed Sequence Alignment track settings

                    # Supertrack settings
                    str_list.append('parent Alignment\n')
                    str_list.append('\n')

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return
