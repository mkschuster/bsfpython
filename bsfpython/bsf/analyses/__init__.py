"""bsf.analyses

A package of classes and methods supporting NGS-specific analyses such as ChIP-Seq or RNA-Seq.
"""

#
# Copyright 2013 - 2016 Michael K. Schuster
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
import re
import warnings

from bsf import Analysis, defaults
from bsf.annotation import AnnotationSheet
from bsf.ngs import Collection, ProcessedRunFolder, Sample
from bsf.executables import Bowtie2, Macs14, Macs2Bdgcmp, Macs2Callpeak, FastQC
from bsf.process import Executable
from bsf.standards import Default


class ChIPSeqComparison(object):
    """ChIP-Seq comparison annotation sheet.

    Attributes:
    @ivar c_name: Control name
    @type c_name: str
    @ivar t_name: Treatment name
    @type t_name: str
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
        """Initialise a C{bsf.analyses.ChIPSeqComparison}.

        @param c_name: Control name
        @type c_name: str
        @param t_name: Treatment name
        @type t_name: str
        @param c_samples: Python C{list} of control C{bsf.ngs.Sample} objects
        @type c_samples: list[bsf.ngs.Sample]
        @param t_samples: Python C{list} of treatment C{bsf.ngs.Sample} objects
        @type t_samples:list[bsf.ngs.Sample]
        @param factor: ChIP factor
        @type factor: str
        @param tissue: Tissue
        @type tissue: str
        @param condition: Condition
        @type condition: str
        @param treatment: Treatment
        @type treatment: str
        @param replicate: replicate number
        @type replicate: int
        @param diff_bind: Run the DiffBind analysis
        @type diff_bind: bool
        @return:
        @rtype:
        """

        super(ChIPSeqComparison, self).__init__()

        # Condition', 'Treatment', 'Replicate',
        # 'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat'

        if c_name is None:
            self.c_name = str()
        else:
            self.c_name = c_name

        if t_name is None:
            self.t_name = str()
        else:
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


class ChIPSeqDiffBindSheet(AnnotationSheet):
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
            AnnotationSheet.check_alphanumeric,
        ],
        'Tissue': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Factor': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Condition': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Treatment': [
            AnnotationSheet.check_alphanumeric,
        ],
        'Replicate': [
            AnnotationSheet.check_numeric,
        ],
        'ControlID': [
            AnnotationSheet.check_alphanumeric,
        ],
        'PeakCaller': [
            AnnotationSheet.check_alphanumeric,
        ],
        'PeakFormat': [
            AnnotationSheet.check_alphanumeric,
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

    def to_file_path(self):
        """Write a C{bsf.analyses.ChIPSeqDiffBindSheet} to a file.

        @return:
        @rtype:
        """

        # Override the method from the super-class to automatically sort before writing to a file.

        self.sort()
        super(ChIPSeqDiffBindSheet, self).to_file_path()

        return


class ChIPSeq(Analysis):
    """The C{bsf.analyses.ChIPSeq} class represents the logic to run a ChIP-Seq-specific C{bsf.Analysis}.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar cmp_file: Comparison file
    @type cmp_file: str | unicode
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode
    @ivar genome_sizes_path: Reference genome (chromosome) sizes file path
    @type genome_sizes_path: str | unicode
    """

    name = 'ChIP-seq Analysis'
    prefix = 'chipseq'

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
            replicate_grouping=True,
            cmp_file=None,
            genome_fasta_path=None,
            genome_sizes_path=None):
        """Initialise a C{bsf.analyses.ChIPSeq}.

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
        @param comparisons: Python C{dict} of Python C{list} objects of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, list[bsf.ngs.Sample]]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param genome_sizes_path: Reference genome (chromosome) sizes file path
        @type genome_sizes_path: str | unicode
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
            comparisons=comparisons,
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = True
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

        if genome_sizes_path is None:
            self.genome_sizes_path = str()
        else:
            self.genome_sizes_path = genome_sizes_path

        self._factor_dict = dict()
        """ @type _factor_dict: dict[str, list[ChIPSeqComparison]] """

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.ChIPSeq} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(ChIPSeq, self).set_configuration(configuration=configuration, section=section)

        # Read a comparison file.

        option = 'replicate_grouping'
        if configuration.config_parser.has_option(section=section, option=option):
            self.replicate_grouping = configuration.config_parser.getboolean(section=section, option=option)

        option = 'cmp_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cmp_file = configuration.config_parser.get(section=section, option=option)

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'genome_sizes'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_sizes_path = configuration.config_parser.get(section=section, option=option)

        return

    # Taken from ConfigParser.RawConfigParser.getboolean()

    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False}

    def _read_comparisons(self, cmp_file):
        """Read a C{bsf.annotation.AnnotationSheet} CSV file from disk.

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
        @param cmp_file: Comparison file path
        @type cmp_file: str | unicode
        @return:
        @rtype:
        """

        if self.debug > 1:
            print '{!r} method _read_comparisons:'.format(self)

        annotation_sheet = AnnotationSheet.from_file_path(file_path=cmp_file)

        # Unfortunately, two passes through the comparison sheet are required.
        # In the first one merge all Sample objects that share the name.
        # Merging Sample objects is currently the only way to pool PairedReads objects,
        # which is required for ChIP-Seq experiments.

        sample_dict = dict()
        """ @type sample_dict: dict[str, bsf.ngs.Sample] """

        # First pass, merge Sample objects, if they have the same name.
        for row_dict in annotation_sheet.row_dicts:
            for prefix in ('Control', 'Treatment'):
                name, samples = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                for o_sample in samples:
                    if o_sample.name not in sample_dict:
                        sample_dict[o_sample.name] = o_sample
                    else:
                        n_sample = Sample.from_samples(sample1=sample_dict[o_sample.name], sample2=o_sample)
                        sample_dict[n_sample.name] = n_sample

        # Second pass, add all Sample objects mentioned in a comparison.
        level_1_dict = dict()
        """ @type level_1_dict: dict[str, dict[str, int]] """

        for row_dict in annotation_sheet.row_dicts:

            c_name, c_samples = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Control')
            t_name, t_samples = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Treatment')

            # ChIP-Seq experiments use the order treatment versus control in comparisons.
            comparison_key = '{}__{}'.format(t_name, c_name)

            # For a successful comparison, both Python list objects of Sample objects have to be defined.

            if not (len(t_samples) and len(c_samples)):
                if self.debug > 1:
                    print 'Redundant comparison line with Treatment {!r} samples {} and Control {!r} samples {}'. \
                        format(t_name, len(t_samples), c_name, len(c_samples))
                continue

            # Add all control Sample or SampleGroup objects to the Sample list.

            for c_sample in c_samples:
                if self.debug > 1:
                    print '  Control Sample name: {!r} file_path: {!r}'.format(c_sample.name, c_sample.file_path)
                    print c_sample.trace(1)
                    # Find the Sample in the unified sample dictionary.
                if c_sample.name in sample_dict:
                    self.add_sample(sample=sample_dict[c_sample.name])

            # Add all treatment Sample or SampleGroup objects to the Sample list.

            for t_sample in t_samples:
                if self.debug > 1:
                    print '  Treatment Sample name: {!r} file_path: {!r}'.format(t_sample.name, t_sample.file_path)
                    print t_sample.trace(1)
                if t_sample.name in sample_dict:
                    self.add_sample(sample=sample_dict[t_sample.name])

            if 'Tissue' in row_dict:
                tissue = row_dict['Tissue']
            else:
                tissue = str()

            if 'Factor' in row_dict:
                factor = row_dict['Factor']
            else:
                factor = defaults.web.chipseq_default_factor

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
                    raise ValueError('Value in field {!r} is not a boolean: {!r}'.format(
                        'DiffBind',
                        row_dict['DiffBind']))
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

            self.comparisons[comparison_key] = ChIPSeqComparison(
                c_name=c_name,
                t_name=t_name,
                c_samples=c_samples,
                t_samples=t_samples,
                factor=factor,
                tissue=tissue,
                condition=condition,
                treatment=treatment,
                replicate=0,
                diff_bind=diff_bind)

        # Sort the comparison keys alphabetically and assign replicate numbers into ChIPSeqComparison objects.

        for key_1 in level_1_dict.keys():
            level_2_dict = level_1_dict[key_1]

            key_2_list = level_2_dict.keys()
            key_2_list.sort(cmp=lambda x, y: cmp(x, y))

            i = 1
            for key_2 in key_2_list:
                if not self.comparisons[key_2].diff_bind:
                    continue
                level_2_dict[key_2] = i
                self.comparisons[key_2].replicate = i
                i += 1

        return

    def run(self):
        """Run a C{bsf.analyses.ChIPSeq} C{bsf.Analysis}.

        @return:
        @rtype:
        """

        super(ChIPSeq, self).run()

        # ChIPSeq requires a genome version.

        if not self.genome_version:
            raise Exception('A ChIPSeq analysis requires a genome_version configuration option.')

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard directory paths.

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file) and not os.path.exists(path=self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        if not self.genome_fasta_path:
            self.genome_fasta_path = Default.absolute_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if self.genome_sizes_path:
            self.genome_sizes_path = os.path.expanduser(self.genome_sizes_path)
            self.genome_sizes_path = os.path.expandvars(self.genome_sizes_path)

        # self._create_bwa_jobs()
        self._create_bowtie2_jobs()
        # self._create_macs14_jobs()
        self._create_macs2_jobs()
        self._create_diffbind_jobs()

        return

    def _create_bowtie2_jobs(self):
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
                "The genome_fasta option has to end with '.fa' or '.fasta'. {!r}".format(self.genome_fasta_path))

        # bowtie2_index = os.path.join(Default.absolute_genomes(self.genome_version),
        #                              'forBowtie2', self.genome_version)

        stage_bowtie2 = self.get_stage(name='bowtie2')

        # Use the bsf_sam2bam.sh script to convert aligned SAM into
        # aligned, sorted, indexed BAM files.

        stage_sam2bam = self.get_stage(name='sam2bam')

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                bowtie2 = stage_bowtie2.add_executable(
                    executable=Bowtie2(
                        name='chipseq_bowtie2_{}'.format(paired_reads_name),
                        analysis=self))

                sam2bam = stage_sam2bam.add_executable(
                    executable=Executable(
                        name='chipseq_sam2bam_{}'.format(paired_reads_name),
                        program='bsf_sam2bam.sh'))
                sam2bam.dependencies.append(bowtie2.name)

                self.set_command_configuration(command=sam2bam)

                # Set Bowtie2 options.

                bowtie2.add_option_short(key='x', value=bowtie2_index)

                reads1 = list()
                reads2 = list()

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if paired_reads.reads_1:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads2.append(paired_reads.reads_2.file_path)

                if len(reads1) and not len(reads2):
                    bowtie2.add_option_short(key='U', value=','.join(reads1))
                elif len(reads1) and len(reads2):
                    bowtie2.add_option_short(key='1', value=','.join(reads1))
                if len(reads2):
                    bowtie2.add_option_short(key='2', value=','.join(reads2))

                # TODO: The following options are properties of the Sample,
                # PairedReads and Reads objects.
                # bowtie2.add_switch_short(key='q')
                # bowtie2.add_switch_short(key='phred33')

                # TODO: It would be good to have code that parses the Illumina Run Info and writes
                # this information into the CSV file.

                # TODO: Andreas' original implementation on Bowtie1 sets -a,
                # which may not be required for Bowtie2.
                # See description of option -k in the bowtie2 manual.
                # TODO: To avoid repeats we may want to set -a?

                bowtie2.add_option_long(key='threads', value=str(stage_bowtie2.threads))

                # Put all sample-specific information into a sub-directory.

                replicate_directory = os.path.join(self.genome_directory,
                                                   'chipseq_bowtie2_{}'.format(paired_reads_name))

                try:
                    os.makedirs(replicate_directory)
                except OSError as exc:  # Python > 2.5
                    if exc.errno == errno.EEXIST and os.path.isdir(replicate_directory):
                        pass
                    else:
                        raise

                # Set Bowtie2 arguments.

                # Store the absolute file_path of the alignment file.

                bowtie2.stdout_path = os.path.join(self.genome_directory,
                                                   replicate_directory,
                                                   paired_reads_name + '.sam')

                # Set bsf_sam2bam.sh options.

                sam2bam.arguments.append(os.path.join(self.genome_directory,
                                                      replicate_directory,
                                                      paired_reads_name))

                if (os.path.exists('{}.bam.bai'.format(paired_reads_name)) and
                        os.path.getsize('{}.bam.bai'.format(paired_reads_name))):
                    bowtie2.submit = False
                    sam2bam.submit = False

        return

    def _create_macs14_jobs(self):
        """Create MACS14 peak caller jobs.
        @return:
        @rtype:
        """

        stage_macs14 = self.get_stage(name='macs14')

        stage_process_macs14 = self.get_stage(name='process_macs14')

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)
            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:
                t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                t_paired_reads_name_list = t_paired_reads_dict.keys()
                t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for t_paired_reads_name in t_paired_reads_name_list:
                    for c_sample in chipseq_comparison.c_samples:
                        c_paired_reads_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_paired_reads_name_list = c_paired_reads_dict.keys()
                        c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                        for c_paired_reads_name in c_paired_reads_name_list:
                            macs14 = stage_macs14.add_executable(
                                executable=Macs14(
                                    name='chipseq_macs14_{}__{}'.format(t_paired_reads_name, c_paired_reads_name),
                                    analysis=self))

                            macs14.dependencies.append('chipseq_sam2bam_' + t_paired_reads_name)
                            macs14.dependencies.append('chipseq_sam2bam_' + c_paired_reads_name)

                            process_macs14 = stage_process_macs14.add_executable(
                                executable=Executable(
                                    name='chipseq_process_macs14_{}__{}'.format(
                                        t_paired_reads_name,
                                        c_paired_reads_name),
                                    program='bsf_chipseq_process_macs14.sh'))
                            process_macs14.dependencies.append(macs14.name)

                            self.set_command_configuration(command=process_macs14)

                            # Set macs14 options.

                            macs14.add_option_long(
                                key='treatment',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(t_paired_reads_name),
                                    '{}.bam'.format(t_paired_reads_name)))
                            macs14.add_option_long(
                                key='control',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(c_paired_reads_name),
                                    '{}.bam'.format(c_paired_reads_name)))

                            # TODO: Experimentally prepend a chipseq_macs14 directory
                            # MACS14 can hopefully also cope with directories specified in the --name option, but
                            # the resulting R script has them set too. Hence the R script has to be started
                            # from the genome_directory. However, the R script needs re-writing anyway, because
                            # it would be better to use the PNG rather than the PDF device for plotting.
                            prefix = 'chipseq_macs14_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)

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
                                    'Unable to set MACS14 parameters for unknown factor {!r}.\n'
                                    'Please use default factor {!r} or adjust Python code if necessary.'.format(
                                        factor, defaults.web.chipseq_default_factor),
                                    UserWarning)

                            # Set macs14 arguments.

                            # Set process_macs14 options.
                            # Set process_macs14 arguments.

                            # Specify the output path as in the macs14 --name option.
                            process_macs14.arguments.append(os.path.join('.', prefix, prefix))
                            process_macs14.arguments.append(self.genome_sizes_path)

        return

    def _create_macs2_jobs(self):
        """Create MACS2 peak caller jobs.

        @return:
        @rtype:
        """

        stage_macs2_callpeak = self.get_stage(name='macs2_callpeak')
        stage_macs2_bdgcmp = self.get_stage(name='macs2_bdgcmp')
        stage_process_macs2 = self.get_stage(name='process_macs2')

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)
            factor = chipseq_comparison.factor.upper()
            for t_sample in chipseq_comparison.t_samples:
                t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                t_paired_reads_name_list = t_paired_reads_dict.keys()
                t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for t_paired_reads_name in t_paired_reads_name_list:
                    for c_sample in chipseq_comparison.c_samples:
                        c_paired_reads_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_paired_reads_name_list = c_paired_reads_dict.keys()
                        c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                        for c_paired_reads_name in c_paired_reads_name_list:
                            macs2_callpeak = stage_macs2_callpeak.add_executable(
                                executable=Macs2Callpeak(
                                    name='chipseq_macs2_callpeak_{}__{}'.format(
                                        t_paired_reads_name,
                                        c_paired_reads_name),
                                    analysis=self))

                            macs2_callpeak.dependencies.append('chipseq_sam2bam_' + t_paired_reads_name)
                            macs2_callpeak.dependencies.append('chipseq_sam2bam_' + c_paired_reads_name)

                            macs2_bdgcmp = stage_macs2_bdgcmp.add_executable(
                                executable=Macs2Bdgcmp(
                                    name='chipseq_macs2_bdgcmp_{}__{}'.format(
                                        t_paired_reads_name,
                                        c_paired_reads_name),
                                    analysis=self))

                            macs2_bdgcmp.dependencies.append(macs2_callpeak.name)

                            process_macs2 = stage_process_macs2.add_executable(
                                executable=Executable(
                                    name='chipseq_process_macs2_{}__{}'.format(
                                        t_paired_reads_name,
                                        c_paired_reads_name),
                                    program='bsf_chipseq_process_macs2.sh'))
                            process_macs2.dependencies.append(macs2_bdgcmp.name)

                            self.set_command_configuration(command=process_macs2)

                            # Set MACS2 callpeak sub-command options.

                            mc2 = macs2_callpeak.sub_command

                            mc2.add_option_long(
                                key='treatment',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(t_paired_reads_name),
                                    '{}.bam'.format(t_paired_reads_name)))

                            mc2.add_option_long(
                                key='control',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(c_paired_reads_name),
                                    '{}.bam'.format(c_paired_reads_name)))

                            # TODO: Experimentally prepend a chipseq_macs2 directory.
                            # MACS2 can cope with directories specified in the --name option, but
                            # the resulting R script has them set too. Hence the R script has to be started
                            # from the genome_directory. However, the R script needs re-writing anyway, because
                            # it would be better to use the PNG rather than the PDF device for plotting.
                            prefix = 'chipseq_macs2_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)

                            replicate_directory = os.path.join(self.genome_directory, prefix)

                            try:
                                os.makedirs(replicate_directory)
                            except OSError as exc:  # Python > 2.5
                                if exc.errno == errno.EEXIST and os.path.isdir(replicate_directory):
                                    pass
                                else:
                                    raise

                            mc2.add_option_long(key='name', value=os.path.join(prefix, prefix))

                            # The 'gsize' option has to be specified via the configuration.ini file.
                            # macs2_callpeak.add_option_long(key='gsize', value=self.genome_sizes_path

                            mc2.add_switch_long(key='bdg')
                            mc2.add_switch_long(key='SPMR')

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
                                mc2.add_switch_long(key='nomodel')
                                # The shiftsize option is no longer supported in MACS 2.1.0
                                # mc2.add_option_long(key='shiftsize', value='73')
                                mc2.add_option_long(key='pvalue', value='1e-3')
                            elif factor == 'H3K56AC':
                                pass
                            elif factor == 'H4K16AC':
                                pass
                            elif factor == 'OTHER':
                                pass
                            else:
                                warnings.warn(
                                    'Unable to set MACS2 parameters for unknown factor {!r}.\n'
                                    'Please use default factor {!r} or adjust Python code if necessary.'.format(
                                        factor, defaults.web.chipseq_default_factor),
                                    UserWarning)

                            # Set macs2_callpeak sub-command arguments.

                            # None to set.

                            # Set macs2_bdgcmp sub-command options.

                            mb2 = macs2_bdgcmp.sub_command

                            mb2.add_option_long(
                                key='tfile',
                                value=os.path.join(prefix, '{}_treat_pileup.bdg'.format(prefix)))

                            mb2.add_option_long(
                                key='cfile',
                                value=os.path.join(prefix, '{}_control_lambda.bdg'.format(prefix)))

                            # Sequencing depth for treatment and control. Aim for setting the --SPMR parameter for
                            # macs2_callpeak to get the track normalised.
                            # --tdepth:
                            # --cdepth:
                            # --pseudocount

                            mb2.add_option_long(
                                key='ofile',
                                value=os.path.join(prefix, '{}_bdgcmp.bdg'.format(prefix)))

                            # --method defaults to ppois i.e. Poisson Pvalue (-log10(pvalue), which yields data
                            # on a logarithmic scale.

                            # mb2.add_option_long(key='--method', value='FE')

                            # Set macs2_bdgcmp arguments.

                            # None to set.

                            # Set process_macs2 options.

                            # None to set.

                            # Set process_macs2 arguments.

                            process_macs2.arguments.append('{}__{}'.format(t_paired_reads_name, c_paired_reads_name))
                            process_macs2.arguments.append(self.genome_sizes_path)

                            if os.path.exists(os.path.join(prefix, '{}_peaks.bb'.format(prefix))):
                                macs2_callpeak.submit = False
                                macs2_bdgcmp.submit = False
                                process_macs2.submit = False

        return

    def _create_diffbind_jobs(self):
        """Create Bioconductor DiffBind jobs.

        @return:
        @rtype:
        """

        stage_diffbind = self.get_stage(name='diffbind')

        # Reorganise the comparisons by factor.

        keys = self.comparisons.keys()
        # keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)
            if not chipseq_comparison.diff_bind:
                continue
            if chipseq_comparison.factor not in self._factor_dict:
                self._factor_dict[chipseq_comparison.factor] = list()
            self._factor_dict[chipseq_comparison.factor].append(chipseq_comparison)

        keys = self._factor_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            factor_list = self._factor_dict[key]
            factor_list.sort(cmp=lambda x, y: cmp(x, y))

            # No comparison for less than two items.

            if len(factor_list) < 2:
                continue

            # Create a directory per factor.

            prefix = 'chipseq_diffbind_{}'.format(key)

            factor_directory = os.path.join(self.genome_directory, prefix)

            try:
                os.makedirs(factor_directory)
            except OSError as exc:  # Python > 2.5
                if exc.errno == errno.EEXIST and os.path.isdir(factor_directory):
                    pass
                else:
                    raise

            job_dependencies = list()

            # Create a new ChIPSeq DiffBind Sheet per factor.

            file_path = os.path.join(factor_directory, 'chipseq_diffbind_{}_samples.csv'.format(key))

            dbs = ChIPSeqDiffBindSheet(file_path=file_path)

            for chipseq_comparison in factor_list:
                assert isinstance(chipseq_comparison, ChIPSeqComparison)
                if not chipseq_comparison.diff_bind:
                    continue

                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                    t_paired_reads_name_list = t_paired_reads_dict.keys()
                    t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                    for t_paired_reads_name in t_paired_reads_name_list:
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)

                            c_paired_reads_name_list = c_paired_reads_dict.keys()
                            c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                            for c_paired_reads_name in c_paired_reads_name_list:
                                row_dict = dict()
                                """ @type row_dict: dict[str, str | unicode] """

                                row_dict['SampleID'] = t_paired_reads_name
                                row_dict['Tissue'] = chipseq_comparison.tissue
                                row_dict['Factor'] = chipseq_comparison.factor
                                row_dict['Condition'] = chipseq_comparison.condition
                                row_dict['Treatment'] = chipseq_comparison.treatment
                                row_dict['Replicate'] = chipseq_comparison.replicate
                                row_dict['bamReads'] = os.path.join(self.genome_directory,
                                                                    'chipseq_bowtie2_{}'.format(t_paired_reads_name),
                                                                    '{}.bam'.format(t_paired_reads_name))
                                row_dict['bamControl'] = os.path.join(self.genome_directory,
                                                                      'chipseq_bowtie2_{}'.format(c_paired_reads_name),
                                                                      '{}.bam'.format(c_paired_reads_name))
                                row_dict['ControlID'] = c_paired_reads_name
                                row_dict['Peaks'] = os.path.join(self.genome_directory,
                                                                 'chipseq_macs2_{}__{}'.
                                                                 format(t_paired_reads_name, c_paired_reads_name),
                                                                 'chipseq_macs2_{}__{}_peaks.xls'.
                                                                 format(t_paired_reads_name, c_paired_reads_name))
                                row_dict['PeakCaller'] = 'macs'
                                row_dict['PeakFormat'] = 'macs'
                                # row_dict['ScoreCol'] = str()
                                # row_dict['LowerBetter'] = str()
                                # row_dict['Counts'] = str()

                                # TODO: Remove once the code works.
                                # ## sas.csv_writer_next(row_dict=row_dict)
                                dbs.row_dicts.append(row_dict)

                                job_dependency = 'chipseq_process_macs2_{}__{}'.format(
                                    t_paired_reads_name,
                                    c_paired_reads_name)
                                job_dependencies.append(job_dependency)

            # TODO: Remove once the code works.
            # ## sas.csv_writer_close()
            dbs.to_file_path()

            # Create the DiffBind job.

            diffbind = stage_diffbind.add_executable(
                executable=Executable(
                    name='chipseq_diffbind_{}'.format(key),
                    program='bsf_chipseq_diffbind.R'))
            diffbind.dependencies.extend(job_dependencies)

            # Add diffbind options.
            self.set_command_configuration(command=diffbind)
            diffbind.add_option_long(key='factor', value=key)
            # diffbind.add_option_long(key='work_directory', value=factor_directory)
            diffbind.add_option_long(key='genome-directory', value=self.genome_directory)
            diffbind.add_option_long(key='sample-annotation', value=file_path)

            if os.path.exists(os.path.join(
                    factor_directory,
                    'chipseq_diffbind_{}_correlation_read_counts.png'.format(key))):
                diffbind.submit = False

        return

    def _report_macs14(self):
        """Create a ChIPSeq report in HTML format and a UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        hub_list = list()
        """ @type hub_list: list[str | unicode] """

        # Write a HTML document.

        report_list = list()
        """ @type report_list: list[str | unicode] """

        report_list += '<h1 id="chip_seq_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
        report_list += '\n'

        report_list += '<p>\n'
        report_list += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
        report_list += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></strong>,\n'
        report_list += 'before peaks are called with '
        report_list += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 1.4</a>\n'
        report_list += 'on a treatment and control sample pair.\n'
        report_list += '</p>\n'
        report_list += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        """ @type options_dict: dict[str, str] """
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        report_list += '<p id="ucsc_track_hub">\n'
        report_list += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
        report_list += 'in their genomic context via the project-specific\n'
        report_list += 'UCSC Genome Browser Track Hub <a href="{}">{}</a>.\n'. \
            format(self.ucsc_track_url(options_dict=options_dict),
                   self.project_name)
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<table id="peak_calling_table">\n'
        report_list += '\n'

        report_list += '<tr>\n'
        report_list += '<th>Peaks</th>\n'
        report_list += '<th>Negative Peaks</th>\n'
        report_list += '<th>R Model</th>\n'
        report_list += '</tr>\n'
        report_list += '\n'

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)
            for t_sample in chipseq_comparison.t_samples:
                t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                t_paired_reads_name_list = t_paired_reads_dict.keys()
                t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for t_paired_reads_name in t_paired_reads_name_list:
                    for c_sample in chipseq_comparison.c_samples:
                        c_paired_reads_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_paired_reads_name_list = c_paired_reads_dict.keys()
                        c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                        for c_paired_reads_name in c_paired_reads_name_list:
                            # prefix = 'chipseq_macs14_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)

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

                                    # Add a UCSC trackDB entry.

                                    # Common trackDb settings.

                                    hub_list += 'track ChIP_{}__{}_{}_{}\n'. \
                                        format(t_paired_reads_name, c_paired_reads_name, state, scaling)
                                    hub_list += 'type bigWig\n'
                                    hub_list += 'shortLabel ChIP_{}__{}_{}_{}\n'. \
                                        format(t_paired_reads_name, c_paired_reads_name, state, scaling)
                                    hub_list += 'longLabel {} ChIP-Seq read counts for {} of {} versus {}\n'. \
                                        format(scaling.capitalize(), state, t_paired_reads_name, c_paired_reads_name)
                                    hub_list += 'bigDataUrl '
                                    hub_list += '{}__{}_MACS_wiggle/{}/{}__{}_{}_afterfiting_all.bw\n'. \
                                        format(t_paired_reads_name, c_paired_reads_name, state,
                                               t_paired_reads_name, c_paired_reads_name, state)
                                    if treatment and not absolute:
                                        hub_list += 'visibility full\n'
                                    else:
                                        hub_list += 'visibility hide\n'
                                        # 'html' is missing from the common settings.

                                    # Common optional settings.

                                    hub_list += 'color {}\n'. \
                                        format(defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                                    # bigWig - Signal graphing track settings.

                                    hub_list += 'graphTypeDefault bar\n'
                                    hub_list += 'maxHeightPixels 100:60:20\n'
                                    hub_list += 'smoothingWindow off\n'

                                    if absolute:
                                        # Track with absolute scaling.
                                        hub_list += 'autoScale on\n'
                                    else:
                                        # Track with relative scaling.
                                        hub_list += 'autoScale off\n'
                                        hub_list += 'viewLimits 0:40\n'

                                    hub_list += '\n'

                            # Add a UCSC trackDB entry for each bigBed peaks file.

                            hub_list += 'track Peaks_{}__{}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'type bigBed\n'
                            hub_list += 'shortLabel Peaks_{}__{}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'longLabel ChIP-Seq peaks for {} versus {}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'bigDataUrl {}__{}_peaks.bb\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'visibility pack\n'
                            # 'html' is missing from the common settings.

                            # Common optional settings.
                            hub_list += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                            hub_list += '\n'

                            # Add web content.

                            report_list += '<tr>\n'

                            # if treatment and absolute:
                            # output += '<td><strong>{}</strong></td>\n'.format(t_paired_reads_name)
                            # if not treatment and absolute:
                            # output += '<td><strong>{}</strong></td>\n'.format(c_paired_reads_name)

                            report_list += '<td><a href="{}__{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                format(t_paired_reads_name, c_paired_reads_name,
                                       t_paired_reads_name, c_paired_reads_name)

                            report_list += '<td>'
                            report_list += '<a href="{}__{}_negative_peaks.xls">Negative peaks {} versus {}</a>'. \
                                format(t_paired_reads_name, c_paired_reads_name,
                                       t_paired_reads_name, c_paired_reads_name)
                            report_list += '</td>\n'

                            report_list += '<td><a href="{}__{}_model.r">R model</a></td>\n'. \
                                format(t_paired_reads_name, c_paired_reads_name,
                                       t_paired_reads_name, c_paired_reads_name)

                            report_list += '</tr>\n'
                            report_list += '\n'

        report_list += '</table>\n'
        report_list += '\n'

        # Add UCSC trackDB entries for each Bowtie2 BAM file.

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                # Add a UCSC trackDB entry.

                # Common trackDb settings.

                hub_list += 'track Alignment_{}\n'.format(paired_reads_name)
                hub_list += 'type bam\n'
                hub_list += 'shortLabel Alignment_{}\n'.format(paired_reads_name)
                hub_list += 'longLabel Bowtie2 alignment of {}\n'. \
                    format(paired_reads_name)
                hub_list += 'bigDataUrl {}.bam\n'. \
                    format(paired_reads_name)
                hub_list += 'visibility hide\n'
                # TODO: The 'html' option is missing.

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                hub_list += '\n'

                # TODO: Including FastQC results would require rearranging the above HTML code.
                # output += '{}__{}_fastqc/fastqc_report.html'. \
                # format(t_paired_reads_name, c_paired_reads_name)

        self.report_to_file(content=report_list)
        self.ucsc_hub_to_file(content=hub_list)

        return

    def report(self):
        """Create a ChIPSeq report in HTML format and a UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # contrast_field_names = ["", "Group1", "Members1", "Group2", "Members2", "DB.edgeR"]

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        hub_list = list()
        """ @type hub_list: list[str | unicode] """

        # Write a HTML document.

        report_list = list()
        """ @type report_list: list[str | unicode] """

        report_list += '<h1 id="chip_seq_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
        report_list += '\n'

        report_list += '<p>\n'
        report_list += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
        report_list += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></strong>,\n'
        report_list += 'before peaks are called with '
        report_list += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 2</a>\n'
        report_list += 'on an enrichment and background sample pair.\n'
        report_list += '</p>\n'
        report_list += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        """ @type options_dict: dict[str, str] """
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        report_list += '<p id="ucsc_track_hub">\n'
        report_list += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
        report_list += 'in their genomic context via the project-specific\n'
        report_list += 'UCSC Genome Browser Track Hub <a href="{}">{}</a>.\n'. \
            format(self.ucsc_track_url(options_dict=options_dict),
                   self.project_name)
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<table id="peak_calling_table">\n'
        report_list += '\n'

        report_list += '<tr>\n'
        report_list += '<th>Peaks</th>\n'
        report_list += '<th>R Model</th>\n'
        report_list += '</tr>\n'
        report_list += '\n'

        # Group via UCSC super tracks.

        hub_list += 'track Alignment\n'
        hub_list += 'shortLabel QC Alignment\n'
        hub_list += 'longLabel ChIP Read Alignment\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group alignment\n'
        hub_list += '\n'

        hub_list += 'track Background\n'
        hub_list += 'shortLabel QC Background\n'
        hub_list += 'longLabel ChIP Background Signal\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group background\n'
        hub_list += '\n'

        hub_list += 'track Enrichment\n'
        hub_list += 'shortLabel QC Enrichment\n'
        hub_list += 'longLabel ChIP Enrichment Signal\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group enrichment\n'
        hub_list += '\n'

        hub_list += 'track Comparison\n'
        hub_list += 'shortLabel ChIP Intensity\n'
        hub_list += 'longLabel ChIP Intensity\n'
        hub_list += 'visibility full\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group comparison\n'
        hub_list += '\n'

        hub_list += 'track Peaks\n'
        hub_list += 'shortLabel ChIP Peaks\n'
        hub_list += 'longLabel ChIP Peaks\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group peaks\n'
        hub_list += '\n'

        hub_list += 'track Summits\n'
        hub_list += 'shortLabel ChIP Summits\n'
        hub_list += 'longLabel ChIP Summits\n'
        hub_list += 'visibility hide\n'
        hub_list += 'superTrack on\n'
        hub_list += 'group summits\n'
        hub_list += '\n'

        # Group via UCSC composite tracks.

        composite_groups = dict()
        """ @type composite_groups: dict[str, int] """

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)
            factor = chipseq_comparison.factor.upper()
            for t_sample in chipseq_comparison.t_samples:
                t_paired_reads_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                t_paired_reads_name_list = t_paired_reads_dict.keys()
                t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for t_paired_reads_name in t_paired_reads_name_list:
                    for c_sample in chipseq_comparison.c_samples:
                        c_paired_reads_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_paired_reads_name_list = c_paired_reads_dict.keys()
                        c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                        for c_paired_reads_name in c_paired_reads_name_list:
                            prefix = 'chipseq_macs2_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)

                            # Add UCSC trackDB entries for each treatment/control and absolute/normalised pair.
                            # NAME_control_lambda.bw
                            # NAME_treat_pileup.bw
                            # NAME_bdgcmp.bw

                            #
                            # Add a UCSC trackDB entry for each NAME_control_lambda.bw file.
                            #

                            composite_group = '{}_background'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                hub_list += 'track {}\n'.format(composite_group)
                                hub_list += 'type bigWig\n'
                                hub_list += 'shortLabel {}\n'.format(composite_group)
                                hub_list += 'longLabel ChIP background signal for factor {!r}\n'. \
                                    format(factor)
                                hub_list += 'visibility dense\n'
                                hub_list += 'compositeTrack on\n'
                                hub_list += 'parent Background\n'
                                hub_list += 'allButtonPair on\n'
                                hub_list += 'centerLabelsDense on\n'
                                hub_list += '\n'

                            # Common trackDb settings.

                            hub_list += 'track {}__{}_background\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            # TODO: The bigWig type must declare the expected signal range.
                            # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                            hub_list += 'type bigWig\n'
                            hub_list += 'shortLabel {}__{}_bkg\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'longLabel ChIP background signal {} versus {}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'bigDataUrl {}/{}_control_lambda.bw\n'. \
                                format(prefix, prefix)
                            hub_list += 'visibility dense\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            hub_list += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            hub_list += 'alwaysZero off\n'
                            hub_list += 'autoScale off\n'
                            hub_list += 'graphTypeDefault bar\n'
                            hub_list += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            hub_list += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            hub_list += 'viewLimits 0:15\n'
                            hub_list += 'viewLimitsMax 0:40\n'
                            hub_list += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            hub_list += 'parent {} on\n'.format(composite_group)
                            hub_list += 'centerLabelsDense off\n'
                            hub_list += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_treat_pileup.bw file.
                            #

                            composite_group = '{}_enrichment'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                hub_list += 'track {}\n'.format(composite_group)
                                hub_list += 'type bigWig\n'
                                hub_list += 'shortLabel {}_enrichment\n'.format(factor)
                                hub_list += 'longLabel ChIP enrichment signal for factor {!r}\n'. \
                                    format(factor)
                                hub_list += 'visibility dense\n'
                                hub_list += 'compositeTrack on\n'
                                hub_list += 'parent Enrichment\n'
                                hub_list += 'allButtonPair on\n'
                                hub_list += 'centerLabelsDense on\n'
                                hub_list += '\n'

                            # Common trackDb settings.

                            hub_list += 'track {}__{}_enrichment\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            # TODO: The bigWig type must declare the expected signal range.
                            hub_list += 'type bigWig\n'
                            hub_list += 'shortLabel {}__{}_enr\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'longLabel ChIP enrichment signal {} versus {}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'bigDataUrl {}/{}_treat_pileup.bw\n'. \
                                format(prefix, prefix)
                            hub_list += 'visibility dense\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            hub_list += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            hub_list += 'alwaysZero off\n'
                            hub_list += 'autoScale off\n'
                            hub_list += 'graphTypeDefault bar\n'
                            hub_list += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            hub_list += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            hub_list += 'viewLimits 0:15\n'
                            hub_list += 'viewLimitsMax 0:40\n'
                            hub_list += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            hub_list += 'parent {} on\n'.format(composite_group)
                            hub_list += 'centerLabelsDense off\n'
                            hub_list += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_bdgcmp.bw file.
                            #

                            composite_group = '{}_intensity'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                hub_list += 'track {}\n'.format(composite_group)
                                hub_list += 'type bigWig\n'
                                hub_list += 'shortLabel {}\n'.format(composite_group)
                                hub_list += 'longLabel ChIP intensity for factor {!r}\n'. \
                                    format(factor)
                                hub_list += 'visibility full\n'
                                hub_list += 'compositeTrack on\n'
                                hub_list += 'parent Comparison\n'
                                hub_list += 'allButtonPair on\n'
                                hub_list += 'centerLabelsDense on\n'
                                hub_list += '\n'

                            # Common trackDb settings.

                            hub_list += 'track {}__{}_intensity\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            # TODO: The bigWig type must declare the expected signal range.
                            hub_list += 'type bigWig\n'
                            hub_list += 'shortLabel {}__{}_int\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'longLabel ChIP intensity {} versus {}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'bigDataUrl {}/{}_bdgcmp.bw\n'. \
                                format(prefix, prefix)
                            hub_list += 'visibility full\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            hub_list += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            hub_list += 'alwaysZero off\n'
                            hub_list += 'autoScale off\n'
                            hub_list += 'graphTypeDefault bar\n'
                            hub_list += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            hub_list += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            hub_list += 'viewLimits 0:15\n'
                            hub_list += 'viewLimitsMax 0:40\n'
                            hub_list += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            hub_list += 'parent {} on\n'.format(composite_group)
                            hub_list += 'centerLabelsDense off\n'
                            hub_list += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_peaks.bb file.
                            #

                            # Common trackDb settings.

                            hub_list += 'track {}__{}_peaks\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'type bigBed\n'
                            hub_list += 'shortLabel {}__{}_peaks\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'longLabel ChIP peaks {} versus {}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'bigDataUrl {}/{}_peaks.bb\n'. \
                                format(prefix, prefix)
                            hub_list += 'visibility pack\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            hub_list += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigBed - Item or region track settings.

                            # Supertrack settings.

                            hub_list += 'parent Peaks\n'
                            hub_list += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_summits.bb file.
                            #

                            # Common trackDb settings.

                            hub_list += 'track {}__{}_summits\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'type bigBed\n'
                            hub_list += 'shortLabel {}__{}_summits\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'longLabel ChIP summits {} versus {}\n'. \
                                format(t_paired_reads_name, c_paired_reads_name)
                            hub_list += 'bigDataUrl {}/{}_summits.bb\n'. \
                                format(prefix, prefix)
                            hub_list += 'visibility pack\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            hub_list += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # Supertrack settings.

                            hub_list += 'parent Summits\n'
                            hub_list += '\n'

                            #
                            # Add web content.
                            #

                            report_list += '<tr>\n'

                            report_list += '<td><a href="{}/{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                format(prefix, prefix, t_paired_reads_name, c_paired_reads_name)

                            report_list += '<td><a href="{}/{}_model.r">R model</a></td>\n'. \
                                format(prefix, prefix)

                            report_list += '</tr>\n'
                            report_list += '\n'

        report_list += '</table>\n'
        report_list += '\n'

        # Add UCSC trackDB entries for each Bowtie2 BAM file.

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                #
                # Add a UCSC trackDB entry for each NAME.bam file.
                #

                # Common trackDb settings.

                hub_list += 'track {}_alignment\n'.format(paired_reads_name)
                hub_list += 'type bam\n'
                hub_list += 'shortLabel {}_alignment\n'.format(paired_reads_name)
                hub_list += 'longLabel {} ChIP read alignment\n'. \
                    format(paired_reads_name)
                hub_list += 'bigDataUrl chipseq_bowtie2_{}/{}.bam\n'. \
                    format(paired_reads_name, paired_reads_name)
                hub_list += 'visibility hide\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                # Supertrack settings.

                hub_list += 'parent Alignment\n'
                hub_list += '\n'

        # Differential binding analysis.

        report_list += '<h2 id="differential_binding">Differential Binding Analysis</h2>\n'
        report_list += '\n'

        report_list += '<table id="differential_binding_table">\n'
        report_list += '\n'

        report_list += '<tr>\n'
        report_list += '<th>Factor and Contrast</th>\n'
        report_list += '<th>Correlation Peak Caller</th>\n'
        report_list += '<th>Correlation Peak Counts</th>\n'
        report_list += '<th>Correlation Analysis</th>\n'
        report_list += '<th>MA Plot</th>\n'
        report_list += '<th>Scatter Plot</th>\n'
        report_list += '<th>PCA Plot</th>\n'
        report_list += '<th>Box Plot</th>\n'
        report_list += '<th>DiffBind Report</th>\n'
        report_list += '</tr>\n'

        keys = self._factor_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            prefix = 'chipseq_diffbind_{}'.format(key)

            report_list += '<tr>\n'

            # ChIP Factor

            report_list += '<td><strong>{}</strong></td>\n'.format(key)

            # Correlation heat map of peak caller scores.

            report_list += '<td>'
            report_list += '<a href="{}/{}_correlation_peak_caller_score.png">'. \
                format(prefix, prefix)
            report_list += '<img alt="DiffBind correlation analysis for factor {}" ' \
                           'src="{}/{}_correlation_peak_caller_score.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            report_list += '</a>'
            report_list += '</td>\n'

            # Correlation heat map of counts.

            report_list += '<td>'
            report_list += '<a href="{}/{}_correlation_read_counts.png">'. \
                format(prefix, prefix)
            report_list += '<img alt="DiffBind correlation analysis for factor {}" ' \
                           'src="{}/{}_correlation_read_counts.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            report_list += '</a>'
            report_list += '</td>\n'

            # Correlation heat map of differential binding analysis.

            report_list += '<td>'
            report_list += '<a href="{}/{}_correlation_analysis.png">'. \
                format(prefix, prefix)
            report_list += '<img alt="DiffBind correlation analysis for factor {}" ' \
                           'src="{}/{}_correlation_analysis.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            report_list += '</a>'
            report_list += '</td>\n'

            report_list += '<td></td>\n'  # MA Plot
            report_list += '<td></td>\n'  # Scatter Plot
            report_list += '<td></td>\n'  # PCA Plot
            report_list += '<td></td>\n'  # Box Plot
            report_list += '<td></td>\n'  # DiffBin Report

            report_list += '</tr>\n'

            # Read the file of contrasts ...

            file_path = os.path.join(self.genome_directory, prefix, '{}_contrasts.csv'.format(prefix))

            if not os.path.exists(path=file_path):
                warnings.warn(
                    'File {!r} does not exist.'.format(file_path),
                    UserWarning)
                continue

            annotation_sheet = AnnotationSheet.from_file_path(file_path=file_path, file_type='excel')

            for row_dict in annotation_sheet.row_dicts:
                suffix = '{}__{}'.format(row_dict['Group1'], row_dict['Group2'])

                report_list += '<tr>\n'

                report_list += '<td>{}</td>\n'.format(suffix)
                report_list += '<td></td>\n'  # Correlation heat map of peak caller scores.
                report_list += '<td></td>\n'  # Correlation heat map of counts.
                report_list += '<td></td>\n'  # Correlation heat map of differential binding analysis.

                # MA Plot

                report_list += '<td>'
                report_list += '<a href="{}/{}_ma_plot_{}.png">'.format(prefix, prefix, suffix)
                report_list += '<img alt="DiffBind MA plot for factor {}" ' \
                               'src="{}/{}_ma_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                report_list += '</a>'
                report_list += '</td>\n'

                # Scatter Plot

                report_list += '<td>'
                report_list += '<a href="{}/{}_scatter_plot_{}.png">'.format(prefix, prefix, suffix)
                report_list += '<img alt="DiffBind scatter plot for factor {}" ' \
                               'src="{}/{}_scatter_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                report_list += '</a>'
                report_list += '</td>\n'

                # Principal Component Analysis Plot

                report_list += '<td>'
                report_list += '<a href="{}/{}_pca_plot_{}.png">'.format(prefix, prefix, suffix)
                report_list += '<img alt="DiffBind PCA plot for factor {}" ' \
                               'src="{}/{}_pca_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                report_list += '</a>'
                report_list += '</td>\n'

                # Box Plot

                report_list += '<td>'
                report_list += '<a href="{}/{}_box_plot_{}.png">'.format(prefix, prefix, suffix)
                report_list += '<img alt="DiffBind Box plot for factor {}" ' \
                               'src="{}/{}_box_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                report_list += '</a>'
                report_list += '</td>\n'

                # DiffBind report

                report_list += '<td><a href="{}/DBA_{}_report_{}.csv">DBA_{}_report_{}</a></td>\n'. \
                    format(prefix, key, suffix, key, suffix)

                report_list += '</tr>\n'

        report_list += '\n'
        report_list += '</table>\n'
        report_list += '\n'

        self.report_to_file(content=report_list)
        self.ucsc_hub_to_file(content=hub_list)

        return


class RunFastQC(Analysis):
    """BSF FastQC-specific Quality Assessment C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar cmp_file: Comparison file
    @type cmp_file: str | unicode
    """

    name = 'FastQC Analysis'
    prefix = 'fastqc'

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
            cmp_file=None):
        """Initialise a C{bsf.analyses.RunFastQC}.

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
        @param stage_list: Python C{list} of BSF C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param comparisons: Python C{dict} of Python C{tuple} objects of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, list[bsf.ngs.Sample]]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @return:
        @rtype:
        """

        super(RunFastQC, self).__init__(
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

        if cmp_file is None:
            self.cmp_file = str()
        else:
            self.cmp_file = cmp_file

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.RunFastQC} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(RunFastQC, self).set_configuration(configuration=configuration, section=section)

        # Read a comparison file.

        # option = 'cmp_file'
        # if configuration.config_parser.has_option(section=section, option=option):
        # self.cmp_file = configuration.config_parser.get(section=section, option=option)
        # Use the sample annotation sheet instead of a separate comparison file.
        option = 'sas_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cmp_file = configuration.config_parser.get(section=section, option=option)

        return

    def _read_comparisons(self, cmp_file):
        """Read a C{bsf.ngs.AnnotationSheet} CSV file from disk.

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
                - Treatment/Control Sample:
                - Treatment/Control Reads:
                - Treatment/Control File:
        @param cmp_file: Comparisons file path
        @type cmp_file: str | unicode
        @return:
        @rtype:
        """

        annotation_sheet = AnnotationSheet.from_file_path(file_path=cmp_file)

        for row_dict in annotation_sheet.row_dicts:
            self.add_sample(sample=self.collection.get_sample_from_row_dict(row_dict=row_dict))

        return

    def run(self):
        """Run a C{bsf.analyses.RunFastQC} C{bsf.Analysis}.

        @return:
        @rtype:
        """

        super(RunFastQC, self).run()

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of BSF Sample objects by the BSF Sample name.
        # This cannot be done in the super-class, because BSF Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_fastqc_jobs()

        return

    def _create_fastqc_jobs(self):
        """Create FASTQC processes.

        @return:
        @rtype:
        """

        # Read configuration options.

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        # Always check each BSF PairedReads object separately.
        replicate_grouping = False

        # TODO: The configuration of projects or samples is still not ideal.
        # Remove both 'samples' and 'projects' options.
        # Ideally, the sample annotation sheet could cope with wildcards,
        # empty cells or even missing columns to register all hierarchical
        # objects underneath.

        if config_parser.has_option(section=config_section, option='projects'):

            # A configuration option 'projects' has been set ...

            prf = ProcessedRunFolder.from_file_path(file_path=self.input_directory, file_type='Automatic')

            projects = config_parser.get(section=config_section, option='projects')
            components = projects.split(',')
            for component in components:
                if component in prf.project_dict:
                    project = prf.project_dict[component]
                    for sample in project.get_all_samples():
                        self.add_sample(sample=sample)
                else:
                    warnings.warn(
                        "Could not find projects component {!r} in RunFolder {!r}.\n"
                        "Check projects configuration value {!r}.".format(component, prf.name, projects),
                        UserWarning)

        elif config_parser.has_option(section=config_section, option='samples'):

            # A configuration option 'samples' has been set ...

            sample = config_parser.get(section=config_section, option='samples')

            if sample == 'auto':

                # Automatically discover sample information from a CASAVA
                # Processed Run Folder and a CASAVA Project name.

                prf = ProcessedRunFolder.from_file_path(file_path=self.input_directory, file_type='Automatic')
                project = prf.project_dict[self.project_name]
                for sample in project.get_all_samples():
                    self.add_sample(sample=sample)

        else:

            # Since the BSF Collection object has complete BSF ProcessedRunFolder objects registered,
            # the sample annotation sheet needs re-reading.

            # TODO: This is now done in set_configuration.
            # annotation_sheet = AnnotationSheet(file_path=sas_file)
            #
            # for row_dict in annotation_sheet.row_dicts:
            # self.add_sample(sample=self.collection.get_sample_from_row_dict(row_dict=row_dict))

            pass

        stage_fastqc = self.get_stage(name='fastqc')

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                fastqc = stage_fastqc.add_executable(
                    executable=FastQC(
                        name='fastqc_{}'.format(paired_reads_name),
                        analysis=self))

                # Set FastQC options.

                fastqc.add_option_long(key='outdir', value=self.project_directory)

                if 'casava' not in fastqc.options and sample.file_type == 'CASAVA':
                    fastqc.add_switch_long(key='casava')

                fastqc.add_option_long(key='threads', value=str(stage_fastqc.threads))

                # Set FastQC arguments.

                reads1 = list()
                reads2 = list()

                for paired_reads in paired_reads_dict[paired_reads_name]:

                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads_1:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads2.append(paired_reads.reads_2.file_path)

                fastqc.arguments.append(' '.join(reads1 + reads2))

        return

    def report(self):
        """Create a C{RunFastQC} report in HTML format.

        @return:
        @rtype:
        """

        # config_parser = self.configuration.config_parser
        # config_section = self.configuration.section_from_instance(self)

        # Get further information.

        # Always check each BSF PairedReads object separately.
        replicate_grouping = False

        # Create a symbolic link containing the project name and a UUID.
        # link_path = self.create_public_project_link(sub_directory=sub_directory)
        # link_name = os.path.basename(link_path.rstrip('/'))

        # Write a HTML document.

        report_list = str()
        """ @type report_list: list[str | unicode] """

        report_list += '<h1 id="fastqc_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
        report_list += '\n'

        report_list += '<p>\n'
        report_list += '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a>\n'
        report_list += 'is a quality control tool for high throughput sequence data.\n'
        report_list += '</p>\n'
        report_list += '\n'

        report_list += '<table id="fastqc_table">\n'
        report_list += '\n'

        report_list += '<tr>\n'
        report_list += '<th>Sample</th>\n'
        report_list += '<th>FastQC Report</th>\n'
        report_list += '<th>Summary</th>\n'
        report_list += '<th>ZIP archive</th>\n'
        report_list += '</tr>\n'
        report_list += '\n'

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping, full=True)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                report_list += '<tr>\n'
                report_list += '<td>{}</td>\n'.format(paired_reads_name)
                report_list += '<td><a href="{}_fastqc/fastqc_report.html"><strong>Report</strong></a></td>\n'. \
                    format(paired_reads_name)
                report_list += '<td><a href="{}_fastqc/summary.txt">Summary</a></td>\n'.format(paired_reads_name)
                report_list += '<td><a href="{}_fastqc.zip">ZIP archive</a></td>\n'.format(paired_reads_name)
                report_list += '</tr>\n'
                report_list += '\n'

        report_list += '</table>\n'
        report_list += '\n'

        self.report_to_file(content=report_list)

        return


class RunBamToFastq(Analysis):
    """BAM or SAM to FASTQ converter sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    """

    name = 'BAM To FastQ'
    prefix = 'bam_to_fastq'

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
            sample_list=None):
        """Initialise a RunBamToFastq.

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
        @param collection: BSF Collection
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{list} objects of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, list[bsf.ngs.Sample]]
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @return:
        @rtype:
        """

        super(RunBamToFastq, self).__init__(
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

        # Nothing else to do for this sub-class ...

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a BSF RunBamToFastq via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(RunBamToFastq, self).set_configuration(configuration=configuration, section=section)

        # Nothing else to do for this sub-class ...

        return

    def run(self):
        """Run a RunBamToFastq C{bsf.Analysis}.

        @return:
        @rtype:
        """

        super(RunBamToFastq, self).run()

        self._convert_bam_to_fastq()

        return

    def _convert_bam_to_fastq(self):
        """Private method to convert all Reads objects in BAM or SAM format into FASTQ format.

        @return:
        @rtype:
        """

        # config_parser = self.configuration.config_parser
        # config_section = self.configuration.section_from_instance(self)

        default = Default.get_global_default()

        # Replicates have to be un-grouped, always!
        # replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')
        replicate_grouping = False

        stage_sam_to_fastq = self.get_stage(name='sam_to_fastq')

        for sample in self.collection.get_all_samples():
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    # In a BSF Paired Reads object, the SAM or BAM file could potentially
                    # occur as reads1 or reads2 instance variable.

                    if paired_reads.reads_1:
                        file_name = paired_reads.reads_1.file_path
                        file_name = file_name.rstrip('/ ')
                        file_name = os.path.basename(file_name)

                        # TODO: The matching part to remove the .bam could be achieved with Bash parameter expansion.
                        match = re.search(pattern=r'(.*)\.bam$', string=file_name)
                        if match:
                            sam_to_fastq = stage_sam_to_fastq.add_executable(
                                executable=Executable(
                                    name='picard_sam_to_fastq_{}_1'.format(paired_reads_name),
                                    program='bsf_bam2fastq.sh'))
                            self.set_command_configuration(command=sam_to_fastq)
                            sam_to_fastq.arguments.append(paired_reads.reads_1.file_path)
                            sam_to_fastq.arguments.append(os.path.join(default.classpath_picard, 'picard.jar'))
                            sam_to_fastq.arguments.append(os.path.join(self.genome_directory, match.group(1)))

                    if paired_reads.reads_2:
                        file_name = paired_reads.reads_2.file_path
                        file_name = file_name.rstrip('/ ')
                        file_name = os.path.basename(file_name)

                        match = re.search(pattern=r'(.*)\.bam$', string=file_name)
                        if match:
                            sam_to_fastq = stage_sam_to_fastq.add_executable(
                                executable=Executable(
                                    name='picard_sam_to_fastq_{}_2'.format(paired_reads_name),
                                    program='bsf_bam2fastq.sh'))
                            self.set_command_configuration(command=sam_to_fastq)
                            sam_to_fastq.arguments.append(paired_reads.reads_2.file_path)
                            sam_to_fastq.arguments.append(os.path.join(default.classpath_picard, 'picard.jar'))
                            sam_to_fastq.arguments.append(os.path.join(self.genome_directory, match.group(1)))

        return
