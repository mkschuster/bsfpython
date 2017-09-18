"""bsf.analyses.chip_seq

A package of classes and methods supporting ChIP-Seq analyses.
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
from pickle import Pickler, HIGHEST_PROTOCOL
import warnings

from bsf import Analysis, defaults, FilePath, Runnable
from bsf.annotation import AnnotationSheet
from bsf.ngs import Sample
from bsf.executables import BWA, Macs14
from bsf.process import Command, Executable, RunnableStep, RunnableStepPicard
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
        """Initialise a C{bsf.analyses.chip_seq.ChIPSeqComparison}.

        @param c_name: Control name
        @type c_name: str
        @param t_name: Treatment name
        @type t_name: str
        @param c_samples: Python C{list} of control C{bsf.ngs.Sample} objects
        @type c_samples: list[bsf.ngs.Sample]
        @param t_samples: Python C{list} of treatment C{bsf.ngs.Sample} objects
        @type t_samples: list[bsf.ngs.Sample]
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


class FilePathChIPSeq(FilePath):
    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.chip_seq.FilePathChIPSeq} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathChIPSeq, self).__init__(prefix=prefix)

        self.replicate_directory = prefix
        self.aligned_sam = os.path.join(prefix, '_'.join((prefix, 'aligned.sam')))
        self.cleaned_sam = os.path.join(prefix, '_'.join((prefix, 'cleaned.sam')))
        self.sorted_bam = os.path.join(prefix, '_'.join((prefix, 'sorted.bam')))
        self.aligned_bam = os.path.join(prefix, '_'.join((prefix, '.bam')))
        self.aligned_bai = os.path.join(prefix, '_'.join((prefix, '.bai')))
        self.aligned_md5 = os.path.join(prefix, '_'.join((prefix, '.bam.md5')))

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
    """The C{bsf.analyses.chip_seq.ChIPSeq} class represents the logic to run a ChIP-Seq-specific C{bsf.Analysis}.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar comparison_path: Comparison file path
    @type comparison_path: str | unicode
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode
    @ivar genome_sizes_path: Reference genome (chromosome) sizes file path
    @type genome_sizes_path: str | unicode
    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
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
            sample_list=None,
            replicate_grouping=True,
            comparison_path=None,
            genome_fasta_path=None,
            genome_sizes_path=None,
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
        @type comparison_path: str | unicode
        @param genome_fasta_path: Reference genome sequence FASTA file path
        @type genome_fasta_path: str | unicode
        @param genome_sizes_path: Reference genome (chromosome) sizes file path
        @type genome_sizes_path: str | unicode
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
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

        if comparison_path is None:
            self.comparison_path = str()
        else:
            self.comparison_path = comparison_path

        if genome_fasta_path is None:
            self.genome_fasta_path = str()
        else:
            self.genome_fasta_path = genome_fasta_path

        if genome_sizes_path is None:
            self.genome_sizes_path = str()
        else:
            self.genome_sizes_path = genome_sizes_path

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
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

        # Read a comparison file.

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

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

        return

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
            annotation_sheet = AnnotationSheet.from_file_path(file_path=self.comparison_path)

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
                            n_sample = Sample.from_samples(sample1=sample_dict[o_sample.name], sample2=o_sample)
                            sample_dict[n_sample.name] = n_sample

            # Second pass, add all Sample objects mentioned in a comparison.
            level_1_dict = dict()
            """ @type level_1_dict: dict[str, dict[str, int]] """

            for row_dict in annotation_sheet.row_dicts:

                c_name, c_sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Control')
                t_name, t_sample_list = self.collection.get_samples_from_row_dict(row_dict=row_dict, prefix='Treatment')

                # ChIP-Seq experiments use the order treatment versus control in comparisons.
                comparison_key = '{}__{}'.format(t_name, c_name)

                # For a successful comparison, both Python list objects of Sample objects have to be defined.

                if not (len(t_sample_list) and len(c_sample_list)):
                    if self.debug > 1:
                        print 'Redundant comparison line with Treatment {!r} samples {} and Control {!r} samples {}'. \
                            format(t_name, len(t_sample_list), c_name, len(c_sample_list))
                    continue

                # Add all control Sample or SampleGroup objects to the Sample list.

                for c_sample in c_sample_list:
                    if self.debug > 1:
                        print '  Control Sample name:', c_sample.name, 'file_path:', c_sample.file_path
                    if self.debug > 2:
                        print c_sample.trace(1)
                        # Find the Sample in the unified sample dictionary.
                    if c_sample.name in sample_dict:
                        self.add_sample(sample=sample_dict[c_sample.name])

                # Add all treatment Sample or SampleGroup objects to the Sample list.

                for t_sample in t_sample_list:
                    if self.debug > 1:
                        print '  Treatment Sample name:', t_sample.name, 'file_path:', t_sample.file_path
                    if self.debug > 2:
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

            for key_1 in level_1_dict.keys():
                level_2_dict = level_1_dict[key_1]

                key_2_list = level_2_dict.keys()
                key_2_list.sort(cmp=lambda x, y: cmp(x, y))

                i = 1
                for key_2 in key_2_list:
                    if not self._comparison_dict[key_2].diff_bind:
                        continue
                    level_2_dict[key_2] = i
                    self._comparison_dict[key_2].replicate = i
                    i += 1

            return

        # Start of the run() method body.

        # Get global defaults.

        default = Default.get_global_default()

        super(ChIPSeq, self).run()

        # ChIPSeq requires a genome version.

        if not self.genome_version:
            raise Exception('A ChIPSeq analysis requires a genome_version configuration option.')

        self.comparison_path = Default.get_absolute_path(file_path=self.comparison_path)

        run_read_comparisons()

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        if not self.genome_fasta_path:
            self.genome_fasta_path = Default.absolute_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if self.genome_sizes_path:
            self.genome_sizes_path = Default.get_absolute_path(file_path=self.genome_sizes_path)

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        # self._create_bwa_jobs()
        self._create_bowtie2_jobs()
        # self._create_macs14_jobs()
        self._create_macs2_jobs()
        self._create_diffbind_jobs()

        return

    def _create_bwa_jobs(self):
        """Create BWA alignment jobs.

        @return:
        @rtype:
        """

        # Get the BWA index.

        # TODO: The BWA index directory needs to be configurable.
        bwa_genome_db = os.path.join(Default.absolute_genome_resource(self.genome_version),
                                     'forBWA_0.7.6a', self.genome_version)

        stage_alignment = self.get_stage(name='chipseq_alignment')

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        for sample in self.sample_list:
            if self.debug > 0:
                print self, 'Sample name:', sample.name
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                prefix_read_group = '_'.join((stage_alignment.name, paired_reads_name))

                file_path_read_group = FilePathChIPSeq(prefix=prefix_read_group)

                self.add_runnable(
                    runnable=Runnable(
                        name=prefix_read_group,
                        code_module='bsf.runnables.chipseq_alignment',
                        working_directory=self.genome_directory,
                        file_path_object=file_path_read_group,
                        debug=self.debug))

                # Step 1: Process per lane.

                bwa = BWA(name='bwa_mem', analysis=self)

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
                    if paired_reads.reads_1:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
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
                    '{}_{}.pkl'.format(stage_alignment.name, paired_reads_name))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_align_lane)
                pickler_file.close()

                # Create a bsf_run_bwa.py job to run the pickled object.

                run_bwa = stage_alignment.add_executable(
                    executable=Executable(
                        name='_'.join((stage_alignment.name, paired_reads_name)),
                        program='bsf_run_bwa.py'))

                # Only submit this Executable if the final result file does not exist.
                if (os.path.exists(os.path.join(self.genome_directory, file_path_read_group.aligned_md5)) and
                        os.path.getsize(
                            os.path.join(self.genome_directory, file_path_read_group.aligned_md5))):
                    run_bwa.submit = False

                    # Set run_bwa options.

                self.set_command_configuration(command=run_bwa)
                run_bwa.add_option_long(key='pickler_path', value=pickler_path)
                run_bwa.add_option_long(key='debug', value=str(self.debug))

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

        stage_alignment = self.get_stage(name='chipseq_alignment')

        # Use the bsf_sam2bam.sh script to convert aligned SAM into
        # aligned, sorted, indexed BAM files.

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        for sample in self.sample_list:
            if self.debug > 0:
                print self, 'Sample name:', sample.name
                print sample.trace(1)

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                prefix = '_'.join((stage_alignment.name, paired_reads_name))

                file_path_read_group = FilePathChIPSeq(prefix=prefix)

                runnable_alignment = self.add_runnable(
                    runnable=Runnable(
                        name=prefix,
                        code_module='bsf.runnables.bowtie2',
                        working_directory=self.genome_directory,
                        file_path_object=file_path_read_group))
                self.set_stage_runnable(
                    stage=stage_alignment,
                    runnable=runnable_alignment)
                # Set dependencies.

                runnable_step = runnable_alignment.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bowtie2',
                        program='bowtie2'))

                # Set Bowtie2 options.

                runnable_step.add_option_short(key='x', value=bowtie2_index)
                runnable_step.add_option_long(key='threads', value=str(stage_alignment.threads))

                reads1 = list()
                reads2 = list()

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if paired_reads.reads_1:
                        reads1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads2.append(paired_reads.reads_2.file_path)

                if len(reads1) and not len(reads2):
                    runnable_step.add_option_short(key='U', value=','.join(reads1))
                elif len(reads1) and len(reads2):
                    runnable_step.add_option_short(key='1', value=','.join(reads1))
                if len(reads2):
                    runnable_step.add_option_short(key='2', value=','.join(reads2))

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
                    runnable_step=RunnableStepPicard(
                        name='picard_clean_sam',
                        java_temporary_path=runnable_alignment.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx4G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='CleanSam'))
                """ @type runnable_step: RunnableStepPicard """
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
                    runnable_step=RunnableStepPicard(
                        name='picard_sort_sam',
                        java_temporary_path=runnable_alignment.get_relative_temporary_directory_path,
                        java_heap_maximum='Xmx6G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='SortSam'))
                """ @type runnable_step: RunnableStepPicard """
                runnable_step.add_picard_option(key='INPUT', value=file_path_read_group.cleaned_sam)
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_read_group.sorted_bam)
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

                # TODO: Create the Executable for the Runnable.

        return

    def _create_macs14_jobs(self):
        """Create MACS14 peak caller jobs.

        @return:
        @rtype:
        """

        stage_macs14 = self.get_stage(name='macs14')
        stage_process_macs14 = self.get_stage(name='process_macs14')

        comparison_name_list = self._comparison_dict.keys()
        comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_name in comparison_name_list:
            chipseq_comparison = self._comparison_dict[comparison_name]
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

        stage_peak_calling = self.get_stage(name='chipseq_peakcalling')

        comparison_name_list = self._comparison_dict.keys()
        comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_name in comparison_name_list:
            chipseq_comparison = self._comparison_dict[comparison_name]
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
                            prefix = '{}_{}__{}'.format(
                                stage_peak_calling.name,
                                t_paired_reads_name,
                                c_paired_reads_name)

                            runnable_peak_calling = self.add_runnable(
                                runnable=Runnable(
                                    name=prefix,
                                    code_module='bsf.runnables.chip_seq_peak_calling',
                                    working_directory=self.genome_directory))

                            macs2_call_peak = runnable_peak_calling.add_runnable_step(
                                runnable_step=RunnableStep(
                                    name='macs2_call_peak',
                                    program='macs2',
                                    sub_command=Command(program='callpeak')))
                            # TODO: Handle the dependencies
                            macs2_call_peak.dependencies.append('chipseq_sam2bam_' + t_paired_reads_name)
                            macs2_call_peak.dependencies.append('chipseq_sam2bam_' + c_paired_reads_name)

                            macs2_bdg_cmp = runnable_peak_calling.add_runnable_step(
                                runnable_step=RunnableStep(
                                    name='chipseq_macs2_bdg_cmp',
                                    program='macs2',
                                    sub_command=Command(program='bdgcmp')
                                ))

                            # TODO: Handle the dependencies
                            macs2_bdg_cmp.dependencies.append(macs2_call_peak.name)

                            process_macs2 = runnable_peak_calling.add_runnable_step(
                                runnable_step=RunnableStep(
                                    name='chipseq_process_macs2',
                                    program='bsf_chipseq_process_macs2.sh'))
                            process_macs2.dependencies.append(macs2_bdg_cmp.name)
                            self.set_command_configuration(command=process_macs2)

                            # Set MACS2 callpeak sub-command options.

                            mc2 = macs2_call_peak.sub_command

                            mc2.add_option_long(
                                key='treatment',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_' + t_paired_reads_name,
                                    t_paired_reads_name + '.bam'))

                            mc2.add_option_long(
                                key='control',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_' + c_paired_reads_name,
                                    c_paired_reads_name + '.bam'))

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
                            # macs2_call_peak.add_option_long(key='gsize', value=self.genome_sizes_path

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

                            # Set macs2_call_peak sub-command arguments.

                            # None to set.

                            # Set macs2_bdg_cmp sub-command options.

                            mb2 = macs2_bdg_cmp.sub_command

                            mb2.add_option_long(
                                key='tfile',
                                value=os.path.join(prefix, prefix + '_treat_pileup.bdg'))

                            mb2.add_option_long(
                                key='cfile',
                                value=os.path.join(prefix, prefix + '_control_lambda.bdg'))

                            # Sequencing depth for treatment and control. Aim for setting the --SPMR parameter for
                            # macs2_call_peak to get the track normalised.
                            # --tdepth:
                            # --cdepth:
                            # --pseudocount

                            mb2.add_option_long(
                                key='ofile',
                                value=os.path.join(prefix, prefix + '_bdgcmp.bdg'))

                            # --method defaults to ppois i.e. Poisson Pvalue (-log10(pvalue), which yields data
                            # on a logarithmic scale.

                            # mb2.add_option_long(key='method', value='FE')

                            # Set macs2_bdg_cmp arguments.

                            # None to set.

                            # Set process_macs2 options.

                            # None to set.

                            # Set process_macs2 arguments.

                            process_macs2.arguments.append('{}__{}'.format(t_paired_reads_name, c_paired_reads_name))
                            process_macs2.arguments.append(self.genome_sizes_path)

                            if os.path.exists(os.path.join(prefix, prefix + '_peaks.bb')):
                                macs2_call_peak.submit = False
                                macs2_bdg_cmp.submit = False
                                process_macs2.submit = False

        return

    def _create_diffbind_jobs(self):
        """Create Bioconductor DiffBind jobs.

        @return:
        @rtype:
        """

        stage_diffbind = self.get_stage(name='diffbind')

        # Reorganise the ChIPSeqComparison objects by factor.

        comparison_name_list = self._comparison_dict.keys()
        # comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_name in comparison_name_list:
            chipseq_comparison = self._comparison_dict[comparison_name]
            if not chipseq_comparison.diff_bind:
                continue
            if chipseq_comparison.factor not in self._factor_dict:
                self._factor_dict[chipseq_comparison.factor] = list()
            self._factor_dict[chipseq_comparison.factor].append(chipseq_comparison)

        comparison_name_list = self._factor_dict.keys()
        comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_name in comparison_name_list:
            factor_list = self._factor_dict[comparison_name]
            factor_list.sort(cmp=lambda x, y: cmp(x, y))

            # No comparison for less than two items.

            if len(factor_list) < 2:
                continue

            # Create a directory per factor.

            prefix = 'chipseq_diffbind_' + comparison_name

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

            file_path = os.path.join(factor_directory, 'chipseq_diffbind_{}_samples.csv'.format(comparison_name))

            dbs = ChIPSeqDiffBindSheet(file_path=file_path)

            for chipseq_comparison in factor_list:
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
                                row_dict['bamReads'] = os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_' + t_paired_reads_name,
                                    t_paired_reads_name + '.bam')
                                row_dict['bamControl'] = os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_' + c_paired_reads_name,
                                    c_paired_reads_name + '.bam')
                                row_dict['ControlID'] = c_paired_reads_name
                                row_dict['Peaks'] = os.path.join(
                                    self.genome_directory,
                                    'chipseq_macs2_{}__{}'.format(t_paired_reads_name, c_paired_reads_name),
                                    'chipseq_macs2_{}__{}_peaks.xls'.format(t_paired_reads_name, c_paired_reads_name))
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
                    name='chipseq_diffbind_' + comparison_name,
                    program='bsf_chipseq_diffbind.R'))
            diffbind.dependencies.extend(job_dependencies)

            # Add diffbind options.
            self.set_command_configuration(command=diffbind)
            diffbind.add_option_long(key='factor', value=comparison_name)
            # diffbind.add_option_long(key='work_directory', value=factor_directory)
            diffbind.add_option_long(key='genome-directory', value=self.genome_directory)
            diffbind.add_option_long(key='sample-annotation', value=file_path)

            if os.path.exists(os.path.join(
                    factor_directory,
                    'chipseq_diffbind_{}_correlation_read_counts.png'.format(comparison_name))):
                diffbind.submit = False

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

            str_list += '<h1 id="chip_seq_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
            str_list += '\n'

            str_list += '<p>\n'
            str_list += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
            str_list += '<strong>'
            str_list += '<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>'
            str_list += '</strong>,\n'
            str_list += 'before peaks are called with '
            str_list += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 1.4</a>\n'
            str_list += 'on a treatment and control sample pair.\n'
            str_list += '</p>\n'
            str_list += '\n'

            # Construct an automatic UCSC Track Hub link.

            str_list += '<p id="ucsc_track_hub">\n'
            str_list += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
            str_list += 'in their genomic context via the project-specific\n'
            str_list += self.ucsc_hub_html_anchor(link_path=link_path)
            str_list += '.'
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<table id="peak_calling_table">\n'
            str_list += '<thead>\n'

            str_list += '<tr>\n'
            str_list += '<th>Peaks</th>\n'
            str_list += '<th>Negative Peaks</th>\n'
            str_list += '<th>R Model</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for comparison_name in comparison_name_list:
                chipseq_comparison = self._comparison_dict[comparison_name]
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping)

                    t_paired_reads_name_list = t_paired_reads_dict.keys()
                    t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                    for t_paired_reads_name in t_paired_reads_name_list:
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)

                            c_paired_reads_name_list = c_paired_reads_dict.keys()
                            c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                            for c_paired_reads_name in c_paired_reads_name_list:
                                # prefix = 'chipseq_macs14_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)
                                str_list += '<tr>\n'
                                # if treatment and absolute:
                                # str_list += '<td><strong>{}</strong></td>\n'.format(t_paired_reads_name)
                                # if not treatment and absolute:
                                # str_list += '<td><strong>{}</strong></td>\n'.format(c_paired_reads_name)

                                # Peaks
                                str_list += '<td><a href="{}__{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name,
                                           t_paired_reads_name, c_paired_reads_name)
                                # Negative Peaks
                                str_list += '<td>'
                                str_list += '<a href="{}__{}_negative_peaks.xls">Negative peaks {} versus {}</a>'. \
                                    format(t_paired_reads_name, c_paired_reads_name,
                                           t_paired_reads_name, c_paired_reads_name)
                                str_list += '</td>\n'
                                # R Model
                                str_list += '<td><a href="{}__{}_model.r">R model</a></td>\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name,
                                           t_paired_reads_name, c_paired_reads_name)

                                str_list += '</tr>\n'
                                str_list += '\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            self.report_to_file(content=str_list)

            return

        def report_hub():
            """Private function to create a UCSC Track Hub.

            @return:
            @rtype:
            """

            str_list = list()
            """ @type str_list: list[str | unicode] """

            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for comparison_name in comparison_name_list:
                chipseq_comparison = self._comparison_dict[comparison_name]
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

                                        #
                                        # Add a UCSC trackDB entry for each bigWig file
                                        #
                                        # Common settings
                                        str_list += 'track ChIP_{}__{}_{}_{}\n'. \
                                            format(t_paired_reads_name, c_paired_reads_name, state, scaling)
                                        str_list += 'type bigWig\n'
                                        str_list += 'shortLabel ChIP_{}__{}_{}_{}\n'. \
                                            format(t_paired_reads_name, c_paired_reads_name, state, scaling)
                                        str_list += 'longLabel {} ChIP-Seq read counts for {} of {} versus {}\n'. \
                                            format(scaling.capitalize(), state, t_paired_reads_name,
                                                   c_paired_reads_name)
                                        str_list += 'bigDataUrl '
                                        str_list += '{}__{}_MACS_wiggle/{}/{}__{}_{}_afterfiting_all.bw\n'. \
                                            format(t_paired_reads_name, c_paired_reads_name, state,
                                                   t_paired_reads_name, c_paired_reads_name, state)
                                        # str_list += 'html ...\n'
                                        if treatment and not absolute:
                                            str_list += 'visibility full\n'
                                        else:
                                            str_list += 'visibility hide\n'

                                        # Common optional settings
                                        str_list += 'color {}\n'. \
                                            format(defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                                        # bigWig - Signal graphing track settings
                                        str_list += 'graphTypeDefault bar\n'
                                        str_list += 'maxHeightPixels 100:60:20\n'
                                        str_list += 'smoothingWindow off\n'
                                        if absolute:
                                            # Track with absolute scaling.
                                            str_list += 'autoScale on\n'
                                        else:
                                            # Track with relative scaling.
                                            str_list += 'autoScale off\n'
                                            str_list += 'viewLimits 0:40\n'

                                        str_list += '\n'

                                #
                                # Add a UCSC trackDB entry for each bigBed peaks file
                                #
                                # Common settings
                                str_list += 'track Peaks_{}__{}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'type bigBed\n'
                                str_list += 'shortLabel Peaks_{}__{}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'longLabel ChIP-Seq peaks for {} versus {}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'bigDataUrl {}__{}_peaks.bb\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                # str_list += 'html ...\n'
                                str_list += 'visibility pack\n'

                                # Common optional settings
                                str_list += 'color {}\n'. \
                                    format(defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                                str_list += '\n'

            # Add UCSC trackDB entries for each Bowtie2 BAM file.

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                paired_reads_name_list = paired_reads_dict.keys()
                paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for paired_reads_name in paired_reads_name_list:
                    #
                    # Add a UCSC trackDB entry.
                    #
                    # Common settings
                    str_list += 'track Alignment_{}\n'.format(paired_reads_name)
                    str_list += 'type bam\n'
                    str_list += 'shortLabel Alignment_{}\n'.format(paired_reads_name)
                    str_list += 'longLabel Bowtie2 alignment of {}\n'. \
                        format(paired_reads_name)
                    str_list += 'bigDataUrl {}.bam\n'. \
                        format(paired_reads_name)
                    # str_list += 'html ...\n'
                    str_list += 'visibility hide\n'

                    # Common optional settings
                    # str_list += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                    # bam - Compressed Sequence Alignment track settings.

                    str_list += '\n'

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

        # contrast_field_names = ["", "Group1", "Members1", "Group2", "Members2", "DB.edgeR"]

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

            str_list += '<h1 id="chip_seq_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
            str_list += '\n'

            str_list += '<p>\n'
            str_list += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
            str_list += '<strong>'
            str_list += '<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>'
            str_list += '</strong>,\n'
            str_list += 'before peaks are called with '
            str_list += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 2</a>\n'
            str_list += 'on an enrichment and background sample pair.\n'
            str_list += '</p>\n'
            str_list += '\n'

            # Construct an automatic UCSC Track Hub link.

            str_list += '<p id="ucsc_track_hub">\n'
            str_list += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
            str_list += 'in their genomic context via the project-specific\n'
            str_list += self.ucsc_hub_html_anchor(link_path=link_path)
            str_list += '.'
            str_list += '</p>\n'
            str_list += '\n'

            str_list += '<table id="peak_calling_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Peaks</th>\n'
            str_list += '<th>R Model</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for comparison_name in comparison_name_list:
                chipseq_comparison = self._comparison_dict[comparison_name]
                for t_sample in chipseq_comparison.t_samples:
                    t_paired_reads_dict = t_sample.get_all_paired_reads(
                        replicate_grouping=self.replicate_grouping)

                    t_paired_reads_name_list = t_paired_reads_dict.keys()
                    t_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                    for t_paired_reads_name in t_paired_reads_name_list:
                        for c_sample in chipseq_comparison.c_samples:
                            c_paired_reads_dict = c_sample.get_all_paired_reads(
                                replicate_grouping=self.replicate_grouping)

                            c_paired_reads_name_list = c_paired_reads_dict.keys()
                            c_paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                            for c_paired_reads_name in c_paired_reads_name_list:
                                prefix = 'chipseq_macs2_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)

                                str_list += '<tr>\n'
                                # Peaks
                                str_list += '<td><a href="{}/{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                    format(prefix, prefix, t_paired_reads_name, c_paired_reads_name)
                                # R Model
                                str_list += '<td><a href="{}/{}_model.r">R model</a></td>\n'. \
                                    format(prefix, prefix)
                                str_list += '</tr>\n'
                                str_list += '\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

            # Differential binding analysis.

            str_list += '<h2 id="differential_binding">Differential Binding Analysis</h2>\n'
            str_list += '\n'

            str_list += '<table id="differential_binding_table">\n'
            str_list += '<thead>\n'
            str_list += '<tr>\n'
            str_list += '<th>Factor and Contrast</th>\n'
            str_list += '<th>Correlation Peak Caller</th>\n'
            str_list += '<th>Correlation Peak Counts</th>\n'
            str_list += '<th>Correlation Analysis</th>\n'
            str_list += '<th>MA Plot</th>\n'
            str_list += '<th>Scatter Plot</th>\n'
            str_list += '<th>PCA Plot</th>\n'
            str_list += '<th>Box Plot</th>\n'
            str_list += '<th>DiffBind Report</th>\n'
            str_list += '</tr>\n'
            str_list += '</thead>\n'
            str_list += '<tbody>\n'

            factor_name_list = self._factor_dict.keys()
            factor_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for factor_name in factor_name_list:
                prefix = 'chipseq_diffbind_{}'.format(factor_name)

                str_list += '<tr>\n'

                # ChIP Factor
                str_list += '<td><strong>{}</strong></td>\n'.format(factor_name)
                # Correlation heat map of peak caller scores
                str_list += '<td>'
                str_list += '<a href="{}/{}_correlation_peak_caller_score.png">'. \
                    format(prefix, prefix)
                str_list += '<img alt="DiffBind correlation analysis for factor {}" ' \
                            'src="{}/{}_correlation_peak_caller_score.png" height="80" width="80">'. \
                    format(factor_name, prefix, prefix)
                str_list += '</a>'
                str_list += '</td>\n'

                # Correlation heat map of counts
                str_list += '<td>'
                str_list += '<a href="{}/{}_correlation_read_counts.png">'. \
                    format(prefix, prefix)
                str_list += '<img alt="DiffBind correlation analysis for factor {}" ' \
                            'src="{}/{}_correlation_read_counts.png" height="80" width="80">'. \
                    format(factor_name, prefix, prefix)
                str_list += '</a>'
                str_list += '</td>\n'

                # Correlation heat map of differential binding analysis
                str_list += '<td>'
                str_list += '<a href="{}/{}_correlation_analysis.png">'. \
                    format(prefix, prefix)
                str_list += '<img alt="DiffBind correlation analysis for factor {}" ' \
                            'src="{}/{}_correlation_analysis.png" height="80" width="80">'. \
                    format(factor_name, prefix, prefix)
                str_list += '</a>'
                str_list += '</td>\n'

                str_list += '<td></td>\n'  # MA Plot
                str_list += '<td></td>\n'  # Scatter Plot
                str_list += '<td></td>\n'  # PCA Plot
                str_list += '<td></td>\n'  # Box Plot
                str_list += '<td></td>\n'  # DiffBin Report

                str_list += '</tr>\n'

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

                    str_list += '<tr>\n'

                    str_list += '<td>{}</td>\n'.format(suffix)
                    str_list += '<td></td>\n'  # Correlation heat map of peak caller scores
                    str_list += '<td></td>\n'  # Correlation heat map of counts
                    str_list += '<td></td>\n'  # Correlation heat map of differential binding analysis

                    # MA Plot
                    str_list += '<td>'
                    str_list += '<a href="{}/{}_ma_plot_{}.png">'.format(prefix, prefix, suffix)
                    str_list += '<img alt="DiffBind MA plot for factor {}" ' \
                                'src="{}/{}_ma_plot_{}.png" height="80" width="80">'. \
                        format(factor_name, prefix, prefix, suffix)
                    str_list += '</a>'
                    str_list += '</td>\n'

                    # Scatter Plot
                    str_list += '<td>'
                    str_list += '<a href="{}/{}_scatter_plot_{}.png">'.format(prefix, prefix, suffix)
                    str_list += '<img alt="DiffBind scatter plot for factor {}" ' \
                                'src="{}/{}_scatter_plot_{}.png" height="80" width="80">'. \
                        format(factor_name, prefix, prefix, suffix)
                    str_list += '</a>'
                    str_list += '</td>\n'

                    # Principal Component Analysis Plot
                    str_list += '<td>'
                    str_list += '<a href="{}/{}_pca_plot_{}.png">'.format(prefix, prefix, suffix)
                    str_list += '<img alt="DiffBind PCA plot for factor {}" ' \
                                'src="{}/{}_pca_plot_{}.png" height="80" width="80">'. \
                        format(factor_name, prefix, prefix, suffix)
                    str_list += '</a>'
                    str_list += '</td>\n'

                    # Box Plot
                    str_list += '<td>'
                    str_list += '<a href="{}/{}_box_plot_{}.png">'.format(prefix, prefix, suffix)
                    str_list += '<img alt="DiffBind Box plot for factor {}" ' \
                                'src="{}/{}_box_plot_{}.png" height="80" width="80">'. \
                        format(factor_name, prefix, prefix, suffix)
                    str_list += '</a>'
                    str_list += '</td>\n'

                    # DiffBind report
                    str_list += '<td><a href="{}/DBA_{}_report_{}.csv">DBA_{}_report_{}</a></td>\n'. \
                        format(prefix, factor_name, suffix, factor_name, suffix)

                    str_list += '</tr>\n'

            str_list += '</tbody>\n'
            str_list += '</table>\n'
            str_list += '\n'

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

            str_list += 'track Alignment\n'
            str_list += 'shortLabel QC Alignment\n'
            str_list += 'longLabel ChIP Read Alignment\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group alignment\n'
            str_list += '\n'

            str_list += 'track Background\n'
            str_list += 'shortLabel QC Background\n'
            str_list += 'longLabel ChIP Background Signal\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group background\n'
            str_list += '\n'

            str_list += 'track Enrichment\n'
            str_list += 'shortLabel QC Enrichment\n'
            str_list += 'longLabel ChIP Enrichment Signal\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group enrichment\n'
            str_list += '\n'

            str_list += 'track Comparison\n'
            str_list += 'shortLabel ChIP Intensity\n'
            str_list += 'longLabel ChIP Intensity\n'
            str_list += 'visibility full\n'
            str_list += 'superTrack on\n'
            str_list += 'group comparison\n'
            str_list += '\n'

            str_list += 'track Peaks\n'
            str_list += 'shortLabel ChIP Peaks\n'
            str_list += 'longLabel ChIP Peaks\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group peaks\n'
            str_list += '\n'

            str_list += 'track Summits\n'
            str_list += 'shortLabel ChIP Summits\n'
            str_list += 'longLabel ChIP Summits\n'
            str_list += 'visibility hide\n'
            str_list += 'superTrack on\n'
            str_list += 'group summits\n'
            str_list += '\n'

            # Group via UCSC composite tracks.

            composite_groups = dict()
            """ @type composite_groups: dict[str, bool] """
            comparison_name_list = self._comparison_dict.keys()
            comparison_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for comparison_name in comparison_name_list:
                chipseq_comparison = self._comparison_dict[comparison_name]
                factor = chipseq_comparison.factor.upper()
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
                                prefix = 'chipseq_macs2_{}__{}'.format(t_paired_reads_name, c_paired_reads_name)

                                # Add UCSC trackDB entries for each treatment/control and absolute/normalised pair.
                                # NAME_control_lambda.bw
                                # NAME_treat_pileup.bw
                                # NAME_bdgcmp.bw

                                #
                                # Add a factor-specific "background" composite track once
                                #
                                composite_group = '{}_background'.format(factor)
                                if composite_group not in composite_groups:
                                    composite_groups[composite_group] = True
                                    str_list += 'track {}\n'.format(composite_group)
                                    str_list += 'type bigWig\n'
                                    str_list += 'shortLabel {}\n'.format(composite_group)
                                    str_list += 'longLabel ChIP background signal for factor {!r}\n'. \
                                        format(factor)
                                    str_list += 'visibility dense\n'
                                    str_list += 'compositeTrack on\n'
                                    str_list += 'parent Background\n'
                                    str_list += 'allButtonPair on\n'
                                    str_list += 'centerLabelsDense on\n'
                                    str_list += '\n'

                                #
                                # Add a UCSC trackDB entry for each NAME_control_lambda.bw file.
                                #
                                # Common settings
                                str_list += 'track {}__{}_background\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                # TODO: The bigWig type must declare the expected signal range.
                                # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                                str_list += 'type bigWig\n'
                                str_list += 'shortLabel {}__{}_bkg\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'longLabel ChIP background signal {} versus {}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'bigDataUrl {}/{}_control_lambda.bw\n'. \
                                    format(prefix, prefix)
                                # str_list += 'html {}\n'.format()
                                str_list += 'visibility dense\n'

                                # Common optional settings
                                str_list += 'color {}\n'. \
                                    format(defaults.web.get_chipseq_colour(factor=factor))

                                # bigWig - Signal graphing track settings
                                str_list += 'alwaysZero off\n'
                                str_list += 'autoScale off\n'
                                str_list += 'graphTypeDefault bar\n'
                                str_list += 'maxHeightPixels 100:60:20\n'
                                # str_list += 'maxWindowToQuery 10000000\n'
                                str_list += 'smoothingWindow 5\n'
                                # str_list += 'transformFunc NONE\n'
                                str_list += 'viewLimits 0:15\n'
                                str_list += 'viewLimitsMax 0:40\n'
                                str_list += 'windowingFunction maximum\n'
                                # str_list += 'yLineMark <#>\n'
                                # str_list += 'yLineOnOff on \n'
                                # str_list += 'gridDefault on\n'

                                # Composite track settings
                                str_list += 'parent {} on\n'.format(composite_group)
                                str_list += 'centerLabelsDense off\n'
                                str_list += '\n'

                                #
                                # Add a factor-specific "enrichment" composite track once
                                #
                                composite_group = '{}_enrichment'.format(factor)
                                if composite_group not in composite_groups:
                                    composite_groups[composite_group] = True
                                    str_list += 'track {}\n'.format(composite_group)
                                    str_list += 'type bigWig\n'
                                    str_list += 'shortLabel {}_enrichment\n'.format(factor)
                                    str_list += 'longLabel ChIP enrichment signal for factor {!r}\n'. \
                                        format(factor)
                                    str_list += 'visibility dense\n'
                                    str_list += 'compositeTrack on\n'
                                    str_list += 'parent Enrichment\n'
                                    str_list += 'allButtonPair on\n'
                                    str_list += 'centerLabelsDense on\n'
                                    str_list += '\n'

                                #
                                # Add a UCSC trackDB entry for each NAME_treat_pileup.bw file.
                                #
                                # Common settings
                                str_list += 'track {}__{}_enrichment\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                # TODO: The bigWig type must declare the expected signal range.
                                str_list += 'type bigWig\n'
                                str_list += 'shortLabel {}__{}_enr\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'longLabel ChIP enrichment signal {} versus {}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'bigDataUrl {}/{}_treat_pileup.bw\n'. \
                                    format(prefix, prefix)
                                # str_list += 'html {}\n'.format()
                                str_list += 'visibility dense\n'

                                # Common optional settings
                                str_list += 'color {}\n'. \
                                    format(defaults.web.get_chipseq_colour(factor=factor))

                                # bigWig - Signal graphing track settings
                                str_list += 'alwaysZero off\n'
                                str_list += 'autoScale off\n'
                                str_list += 'graphTypeDefault bar\n'
                                str_list += 'maxHeightPixels 100:60:20\n'
                                # str_list += 'maxWindowToQuery 10000000\n'
                                str_list += 'smoothingWindow 5\n'
                                # str_list += 'transformFunc NONE\n'
                                str_list += 'viewLimits 0:15\n'
                                str_list += 'viewLimitsMax 0:40\n'
                                str_list += 'windowingFunction maximum\n'
                                # str_list += 'yLineMark <#>\n'
                                # str_list += 'yLineOnOff on \n'
                                # str_list += 'gridDefault on\n'

                                # Composite track settings
                                str_list += 'parent {} on\n'.format(composite_group)
                                str_list += 'centerLabelsDense off\n'
                                str_list += '\n'

                                #
                                # Add a factor-specific "intensity" composite track once
                                #
                                composite_group = '{}_intensity'.format(factor)
                                if composite_group not in composite_groups:
                                    composite_groups[composite_group] = True
                                    str_list += 'track {}\n'.format(composite_group)
                                    str_list += 'type bigWig\n'
                                    str_list += 'shortLabel {}\n'.format(composite_group)
                                    str_list += 'longLabel ChIP intensity for factor {!r}\n'. \
                                        format(factor)
                                    str_list += 'visibility full\n'
                                    str_list += 'compositeTrack on\n'
                                    str_list += 'parent Comparison\n'
                                    str_list += 'allButtonPair on\n'
                                    str_list += 'centerLabelsDense on\n'
                                    str_list += '\n'

                                #
                                # Add a UCSC trackDB entry for each NAME_bdgcmp.bw file.
                                #
                                # Common settings
                                str_list += 'track {}__{}_intensity\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                # TODO: The bigWig type must declare the expected signal range.
                                str_list += 'type bigWig\n'
                                str_list += 'shortLabel {}__{}_int\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'longLabel ChIP intensity {} versus {}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'bigDataUrl {}/{}_bdgcmp.bw\n'. \
                                    format(prefix, prefix)
                                # str_list += 'html {}\n'.format()
                                str_list += 'visibility full\n'

                                # Common optional settings
                                str_list += 'color {}\n'. \
                                    format(defaults.web.get_chipseq_colour(factor=factor))

                                # bigWig - Signal graphing track settings
                                str_list += 'alwaysZero off\n'
                                str_list += 'autoScale off\n'
                                str_list += 'graphTypeDefault bar\n'
                                str_list += 'maxHeightPixels 100:60:20\n'
                                # str_list += 'maxWindowToQuery 10000000\n'
                                str_list += 'smoothingWindow 5\n'
                                # str_list += 'transformFunc NONE\n'
                                str_list += 'viewLimits 0:15\n'
                                str_list += 'viewLimitsMax 0:40\n'
                                str_list += 'windowingFunction maximum\n'
                                # str_list += 'yLineMark <#>\n'
                                # str_list += 'yLineOnOff on \n'
                                # str_list += 'gridDefault on\n'

                                # Composite track settings
                                str_list += 'parent {} on\n'.format(composite_group)
                                str_list += 'centerLabelsDense off\n'
                                str_list += '\n'

                                #
                                # Add a UCSC trackDB entry for each NAME_peaks.bb file.
                                #
                                # Common settings
                                str_list += 'track {}__{}_peaks\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'type bigBed\n'
                                str_list += 'shortLabel {}__{}_peaks\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'longLabel ChIP peaks {} versus {}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'bigDataUrl {}/{}_peaks.bb\n'. \
                                    format(prefix, prefix)
                                # str_list += 'html {}\n'.format()
                                str_list += 'visibility pack\n'

                                # Common optional settings
                                str_list += 'color {}\n'. \
                                    format(defaults.web.get_chipseq_colour(factor=factor))

                                # bigBed - Item or region track settings.

                                # Supertrack settings
                                str_list += 'parent Peaks\n'
                                str_list += '\n'

                                #
                                # Add a UCSC trackDB entry for each NAME_summits.bb file.
                                #
                                # Common settings
                                str_list += 'track {}__{}_summits\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'type bigBed\n'
                                str_list += 'shortLabel {}__{}_summits\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'longLabel ChIP summits {} versus {}\n'. \
                                    format(t_paired_reads_name, c_paired_reads_name)
                                str_list += 'bigDataUrl {}/{}_summits.bb\n'. \
                                    format(prefix, prefix)
                                # str_list += 'html {}\n'.format()
                                str_list += 'visibility pack\n'

                                # Common optional settings
                                str_list += 'color {}\n'. \
                                    format(defaults.web.get_chipseq_colour(factor=factor))

                                # Supertrack settings
                                str_list += 'parent Summits\n'
                                str_list += '\n'

            # Add UCSC trackDB entries for each Bowtie2 BAM file.

            for sample in self.sample_list:
                paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                paired_reads_name_list = paired_reads_dict.keys()
                paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

                for paired_reads_name in paired_reads_name_list:
                    #
                    # Add a UCSC trackDB entry for each NAME.bam file.
                    #
                    # Common settings
                    str_list += 'track {}_alignment\n'.format(paired_reads_name)
                    str_list += 'type bam\n'
                    str_list += 'shortLabel {}_alignment\n'.format(paired_reads_name)
                    str_list += 'longLabel {} ChIP read alignment\n'. \
                        format(paired_reads_name)
                    str_list += 'bigDataUrl chipseq_bowtie2_{}/{}.bam\n'. \
                        format(paired_reads_name, paired_reads_name)
                    # str_list += 'html {}\n'.format()
                    str_list += 'visibility hide\n'

                    # Common optional settings
                    # str_list += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                    # bam - Compressed Sequence Alignment track settings

                    # Supertrack settings
                    str_list += 'parent Alignment\n'
                    str_list += '\n'

            self.ucsc_hub_to_file(content=str_list)

            return

        report_html()
        report_hub()

        return
