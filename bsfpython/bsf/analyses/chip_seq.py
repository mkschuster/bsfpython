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

from bsf import Analysis, defaults, Runnable
from bsf.annotation import AnnotationSheet, ChIPSeqDiffBindSheet
from bsf.data import Sample
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
    @ivar c_samples: Python C{list} of control C{bsf.data.Sample} objects
    @type c_samples: list[bsf.data.Sample]
    @ivar t_samples: Python C{list} of treatment C{bsf.data.Sample} objects
    @type t_samples: list[bsf.data.Sample]
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
        """Initialise a C{bsf.analyses.chip_seq.ChIPSeqComparison} object.

        @param c_name: Control name
        @type c_name: str
        @param t_name: Treatment name
        @type t_name: str
        @param c_samples: Python C{list} of control C{bsf.data.Sample} objects
        @type c_samples: list[bsf.data.Sample]
        @param t_samples: Python C{list} of treatment C{bsf.data.Sample} objects
        @type t_samples: list[bsf.data.Sample]
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
            self.replicate = int(x=0)
        else:
            self.replicate = replicate

        if diff_bind is None:
            self.diff_bind = True
        else:
            self.diff_bind = diff_bind

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
    @ivar cmp_file: Comparison file
    @type cmp_file: str | unicode
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
            drms_list=None,
            collection=None,
            comparisons=None,
            samples=None,
            replicate_grouping=True,
            cmp_file=None,
            genome_fasta_path=None,
            genome_sizes_path=None,
            classpath_picard=None):
        """Initialise a C{bsf.analyses.chip_seq.ChIPSeq} object.

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
        @param drms_list: Python C{list} of C{DRMS} objects
        @type drms_list: list[DRMS]
        @param collection: C{bsf.data.Collection}
        @type collection: bsf.data.Collection
        @param comparisons: Python C{dict} of Python C{list} objects of C{bsf.data.Sample} objects
        @type comparisons: dict[str, list[bsf.data.Sample]]
        @param samples: Python C{list} of C{bsf.data.Sample} objects
        @type samples: list[bsf.data.Sample]
        @param replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
        @type replicate_grouping: bool
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
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
            drms_list=drms_list,
            collection=collection,
            comparisons=comparisons,
            samples=samples)

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

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.chip_seq.ChIPSeq} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a
        configuration option remain unchanged.
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

        # Get the Picard tools Java Archive (JAR) class path directory.

        option = 'classpath_picard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_picard = configuration.config_parser.get(section=section, option=option)

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
        level1_dict = dict()

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

                if value in level1_dict:
                    level2_dict = level1_dict[value]
                else:
                    level2_dict = dict()
                    level1_dict[value] = level2_dict

                level2_dict[comparison_key] = 0

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

        # Sort the comparison keys alphabetically and assign replicate numbers into ChiPSeqComparison objects.

        for key1 in level1_dict.keys():

            level2_dict = level1_dict[key1]

            keys2 = level2_dict.keys()
            keys2.sort(cmp=lambda x, y: cmp(x, y))

            i = 1
            for key2 in keys2:
                if not self.comparisons[key2].diff_bind:
                    continue
                level2_dict[key2] = i
                self.comparisons[key2].replicate = i
                i += 1

        return

    def run(self):
        """Run this C{bsf.analyses.chip_seq.ChIPSeq} C{bsf.Analysis}.

        @return:
        @rtype:
        """

        # Get global defaults.

        default = Default.get_global_default()

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

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        if not self.genome_fasta_path:
            self.genome_fasta_path = Default.absolute_genome_fasta(
                genome_version=self.genome_version,
                genome_index='bowtie2')

        if self.genome_sizes_path:
            self.genome_sizes_path = os.path.expanduser(self.genome_sizes_path)
            self.genome_sizes_path = os.path.expandvars(self.genome_sizes_path)

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
        bwa_genome_db = os.path.join(Default.absolute_genomes(self.genome_version),
                                     'forBWA_0.7.6a', self.genome_version)

        alignment_drms = self.get_drms(name='chipseq_alignment')

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
            # Python str key and Python list of Python list objects
            # of PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                prefix = '_'.join((alignment_drms.name, replicate_key))

                file_path_dict = {
                    'temporary_directory': '_'.join((prefix, 'temporary')),
                    'replicate_directory': prefix,
                    'aligned_sam': os.path.join(prefix, '_'.join((prefix, 'aligned.sam'))),
                    'cleaned_sam': os.path.join(prefix, '_'.join((prefix, 'cleaned.sam'))),
                    'sorted_bam': os.path.join(prefix, '_'.join((prefix, 'sorted.bam'))),
                }

                self.add_runnable(
                    runnable=Runnable(
                        name=prefix,
                        code_module='bsf.runnables.chipseq_alignment',
                        working_directory=self.genome_directory,
                        file_path_dict=file_path_dict,
                        debug=self.debug))

                # Step 1: Process per lane.

                bwa = BWA(name='bwa_mem', analysis=self)

                bwa_mem = bwa.sub_command

                # Set BWA mem options.

                # Allow as many threads as defined in the corresponding DRMS object.
                bwa_mem.add_option_short(key='t', value=str(alignment_drms.threads))
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

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)
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

                file_path_chipseq_alignment = {
                    # TODO: The name for the aligned BAM is constructed by the bsf_run_bwa.py script.
                    # It is currently based on the alignment_drms.name and replicate_key.
                    # The script should also be changed to pre-set all file names beforehand.
                    'aligned_bam': '{}_{}.bam'.format(alignment_drms.name, replicate_key),
                    'aligned_bai': '{}_{}.bai'.format(alignment_drms.name, replicate_key),
                    'aligned_md5': '{}_{}.bam.md5'.format(alignment_drms.name, replicate_key),
                }

                # Normally, the bwa object would be pushed onto the drms list.
                # Experimentally, use Pickler to serialize the Executable object into a file.

                pickler_dict_align_lane = {
                    'prefix': alignment_drms.name,
                    'replicate_key': replicate_key,
                    'classpath_picard': self.classpath_picard,
                    'bwa_executable': bwa,
                }

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(alignment_drms.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_align_lane)
                pickler_file.close()

                # Create a bsf_run_bwa.py job to run the pickled object.

                run_bwa = alignment_drms.add_executable(
                    executable=Executable(
                        name='_'.join((alignment_drms.name, replicate_key)),
                        program='bsf_run_bwa.py'))

                # Only submit this Executable if the final result file does not exist.
                if (os.path.exists(os.path.join(self.genome_directory, file_path_chipseq_alignment['aligned_md5'])) and
                        os.path.getsize(
                            os.path.join(self.genome_directory, file_path_chipseq_alignment['aligned_md5']))):
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

        alignment_drms = self.get_drms(name='chipseq_alignment')

        # Initialise the Distributed Resource Management System objects for Bowtie2.

        # bowtie2_drms = self.add_drms(
        #     drms=DRMS.from_analysis(
        #         name='chipseq_bowtie2',
        #         working_directory=self.genome_directory,
        #         analysis=self))

        # Use the bsf_sam2bam.sh script to convert aligned SAM into
        # aligned, sorted, indexed BAM files.

        # sam2bam_drms = self.add_drms(
        #     drms=DRMS.from_analysis(
        #         name='sam2bam',
        #         working_directory=self.genome_directory,
        #         analysis=self))

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
            # Python str key and Python list of Python list objects
            # of PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                prefix = '_'.join((alignment_drms.name, replicate_key))

                file_path_dict = {
                    'temporary_directory': '_'.join((prefix, 'temporary')),
                    'replicate_directory': prefix,
                    'aligned_sam': os.path.join(prefix, '_'.join((prefix, 'aligned.sam'))),
                    'cleaned_sam': os.path.join(prefix, '_'.join((prefix, 'cleaned.sam'))),
                    'sorted_bam': os.path.join(prefix, '_'.join((prefix, '.bam'))),
                }

                runnable_alignment = self.add_runnable(
                    runnable=Runnable(
                        name=prefix,
                        code_module='bsf.runnables.bowtie2',
                        working_directory=self.genome_directory,
                        file_path_dict=file_path_dict))
                self.set_drms_runnable(
                    drms=alignment_drms,
                    runnable=runnable_alignment)
                # Set dependencies.

                runnable_step = runnable_alignment.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bowtie2',
                        program='bowtie2'))

                # Set Bowtie2 options.

                runnable_step.add_option_short(key='x', value=bowtie2_index)
                runnable_step.add_option_long(key='threads', value=str(alignment_drms.threads))

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

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

                runnable_step.stdout_path = file_path_dict['aligned_sam']

                # Run Picard CleanSam to convert the aligned SAM file into a cleaned SAM file.

                runnable_step = runnable_alignment.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='picard_clean_sam',
                        java_temporary_path=file_path_dict['temporary_directory'],
                        java_heap_maximum='Xmx4G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='CleanSam'))
                assert isinstance(runnable_step, RunnableStepPicard)
                runnable_step.add_picard_option(key='INPUT', value=file_path_dict['aligned_sam'])
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict['cleaned_sam'])
                runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict['temporary_directory'])
                runnable_step.add_picard_option(key='VERBOSITY', value='WARNING')
                runnable_step.add_picard_option(key='QUIET', value='false')
                runnable_step.add_picard_option(key='VALIDATION_STRINGENCY', value='STRICT')

                # Run Picard SortSam to convert the cleaned SAM file into a coordinate sorted BAM file.

                runnable_step = runnable_alignment.add_runnable_step(
                    runnable_step=RunnableStepPicard(
                        name='picard_sort_sam',
                        java_temporary_path=file_path_dict['temporary_directory'],
                        java_heap_maximum='Xmx6G',
                        java_jar_path=os.path.join(self.classpath_picard, 'picard.jar'),
                        picard_command='SortSam'))
                assert isinstance(runnable_step, RunnableStepPicard)
                runnable_step.add_picard_option(key='INPUT', value=file_path_dict['cleaned_sam'])
                runnable_step.add_picard_option(key='OUTPUT', value=file_path_dict['sorted_bam'])
                runnable_step.add_picard_option(key='SORT_ORDER', value='coordinate')
                runnable_step.add_picard_option(key='TMP_DIR', value=file_path_dict['temporary_directory'])
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
        """Create Macs14 peak caller jobs.

        @return:
        @rtype:
        """

        macs14_drms = self.get_drms(name='macs14')
        process_macs14_drms = self.get_drms(name='process_macs14')

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)

            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
                # Python str key and Python list of Python list objects
                # of PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_replicate_keys = c_replicate_dict.keys()
                        c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                        for c_replicate_key in c_replicate_keys:

                            macs14 = macs14_drms.add_executable(
                                executable=Macs14(
                                    name='chipseq_macs14_{}__{}'.format(t_replicate_key, c_replicate_key),
                                    analysis=self))

                            macs14.dependencies.append('chipseq_sam2bam_' + t_replicate_key)
                            macs14.dependencies.append('chipseq_sam2bam_' + c_replicate_key)

                            process_macs14 = process_macs14_drms.add_executable(
                                executable=Executable(
                                    name='chipseq_process_macs14_{}__{}'.format(t_replicate_key, c_replicate_key),
                                    program='bsf_chipseq_process_macs14.sh'))
                            process_macs14.dependencies.append(macs14.name)
                            self.set_command_configuration(command=process_macs14)

                            # Set macs14 options.

                            macs14.add_option_long(
                                key='treatment',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(t_replicate_key),
                                    '{}.aligned.sorted.bam'.format(t_replicate_key)))
                            macs14.add_option_long(
                                key='control',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(c_replicate_key),
                                    '{}.aligned.sorted.bam'.format(c_replicate_key)))

                            # TODO: Experimentally prepend a chipseq_macs14 directory
                            # MACS14 can hopefully also cope with directories specified in the --name option, but
                            # the resulting R script has them set too. Hence the R script has to be started
                            # from the genome_directory. However, the R script needs re-writing anyway, because
                            # it would be better to use the PNG rather than the PDF device for plotting.
                            prefix = 'chipseq_macs14_{}__{}'.format(t_replicate_key, c_replicate_key)

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
        """Create Macs2 peak caller jobs.

        @return:
        @rtype:
        """

        peakcalling_drms = self.get_drms(name='chipseq_peakcalling')
        # macs2_callpeak_drms = self.get_drms(name='macs2_callpeak')
        # macs2_bdgcmp_drms = self.get_drms(name='macs2_bdgcmp')
        # process_macs2_drms = self.get_drms(name='process_macs2')

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)

            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
                # Python str key and Python list of Python list objects
                # of PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_replicate_keys = c_replicate_dict.keys()
                        c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                        for c_replicate_key in c_replicate_keys:

                            prefix = '{}_{}__{}'.format(peakcalling_drms.name, t_replicate_key, c_replicate_key)

                            file_path_dict = {
                                'temporary_directory': '_'.join((prefix, 'temporary')),
                                'replicate_directory': prefix,
                            }

                            runnable = self.add_runnable(
                                runnable=Runnable(
                                    name=prefix,
                                    code_module='bsf.runnables.chip_seq_peak_calling',
                                    working_directory=self.genome_directory,
                                    file_path_dict=file_path_dict))

                            macs2_call_peak = runnable.add_runnable_step(
                                runnable_step=RunnableStep(
                                    name='macs2_call_peak',
                                    program='macs2',
                                    sub_command=Command(program='callpeak')))
                            # TODO: Handle the dependencies
                            macs2_call_peak.dependencies.append('chipseq_sam2bam_' + t_replicate_key)
                            macs2_call_peak.dependencies.append('chipseq_sam2bam_' + c_replicate_key)

                            macs2_bdg_cmp = runnable.add_runnable_step(
                                runnable_step=RunnableStep(
                                    name='chipseq_macs2_bdg_cmp',
                                    program='macs2',
                                    sub_command=Command(program='bdgcmp')
                                ))

                            # TODO: Handle the dependencies
                            macs2_bdg_cmp.dependencies.append(macs2_call_peak.name)

                            process_macs2 = runnable.add_runnable_step(
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
                                    'chipseq_bowtie2_{}'.format(t_replicate_key),
                                    '{}.aligned.sorted.bam'.format(t_replicate_key)))

                            mc2.add_option_long(
                                key='control',
                                value=os.path.join(
                                    self.genome_directory,
                                    'chipseq_bowtie2_{}'.format(c_replicate_key),
                                    '{}.aligned.sorted.bam'.format(c_replicate_key)))

                            # TODO: Experimentally prepend a chipseq_macs2 directory.
                            # MACS2 can cope with directories specified in the --name option, but
                            # the resulting R script has them set too. Hence the R script has to be started
                            # from the genome_directory. However, the R script needs re-writing anyway, because
                            # it would be better to use the PNG rather than the PDF device for plotting.
                            prefix = 'chipseq_macs2_{}__{}'.format(t_replicate_key, c_replicate_key)

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

                            mb2 = macs2_bdg_cmp.sub_command

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

                            process_macs2.arguments.append('{}__{}'.format(t_replicate_key, c_replicate_key))
                            process_macs2.arguments.append(self.genome_sizes_path)

        return

    def _create_diffbind_jobs(self):
        """Create Bioconductor DiffBind jobs.

        @return:
        @rtype:
        """

        diffbind_drms = self.get_drms(name='diffbind')

        # Reorganise the comparisons by factor.

        self._factor_dict = dict()

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

                for t_sample in chipseq_comparison.t_samples:

                    t_replicate_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                    # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
                    # Python str key and Python list of Python list objects
                    # of PairedReads objects.

                    t_replicate_keys = t_replicate_dict.keys()
                    t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                    for t_replicate_key in t_replicate_keys:

                        for c_sample in chipseq_comparison.c_samples:

                            c_replicate_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                            c_replicate_keys = c_replicate_dict.keys()
                            c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                            for c_replicate_key in c_replicate_keys:
                                row_dict = dict()

                                row_dict['SampleID'] = t_replicate_key
                                row_dict['Tissue'] = chipseq_comparison.tissue
                                row_dict['Factor'] = chipseq_comparison.factor
                                row_dict['Condition'] = chipseq_comparison.condition
                                row_dict['Treatment'] = chipseq_comparison.treatment
                                row_dict['Replicate'] = chipseq_comparison.replicate
                                row_dict['bamReads'] = os.path.join(self.genome_directory,
                                                                    'chipseq_bowtie2_{}'.format(t_replicate_key),
                                                                    '{}.aligned.sorted.bam'.format(t_replicate_key))
                                row_dict['bamControl'] = os.path.join(self.genome_directory,
                                                                      'chipseq_bowtie2_{}'.format(c_replicate_key),
                                                                      '{}.aligned.sorted.bam'.format(c_replicate_key))
                                row_dict['ControlID'] = c_replicate_key
                                row_dict['Peaks'] = os.path.join(self.genome_directory,
                                                                 'chipseq_macs2_{}__{}'.
                                                                 format(t_replicate_key, c_replicate_key),
                                                                 'chipseq_macs2_{}__{}_peaks.xls'.
                                                                 format(t_replicate_key, c_replicate_key))
                                row_dict['PeakCaller'] = 'macs'
                                row_dict['PeakFormat'] = 'macs'
                                # row_dict['ScoreCol'] = str()
                                # row_dict['LowerBetter'] = str()
                                # row_dict['Counts'] = str()

                                # TODO: Remove once the code works.
                                # ## sas.csv_writer_next(row_dict=row_dict)
                                dbs.row_dicts.append(row_dict)

                                job_dependency = 'chipseq_process_macs2_{}__{}'.format(t_replicate_key, c_replicate_key)
                                job_dependencies.append(job_dependency)

            # TODO: Remove once the code works.
            # ## sas.csv_writer_close()
            dbs.to_file_path()

            # Create the DiffBind job.

            diffbind = diffbind_drms.add_executable(
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

        return

    def _report_macs14(self):
        """Create a C{bsf.analyses.chip_seq.ChIPSeq} report in HTML format and a UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        output_hub = str()

        # Write a HTML document.

        output_html = str()

        output_html += '<h1 id="chip_seq_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
        output_html += '\n'

        output_html += '<p>\n'
        output_html += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
        output_html += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></strong>,\n'
        output_html += 'before peaks are called with '
        output_html += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 1.4</a>\n'
        output_html += 'on a treatment and control sample pair.\n'
        output_html += '</p>\n'
        output_html += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output_html += '<p id="ucsc_track_hub">\n'
        output_html += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
        output_html += 'in their genomic context via the project-specific\n'
        output_html += 'UCSC Genome Browser Track Hub <a href="{}">{}</a>.\n'. \
            format(self.ucsc_track_url(options_dict=options_dict),
                   self.project_name)
        output_html += '</p>\n'
        output_html += '\n'

        output_html += '<table id="peak_calling_table">\n'
        output_html += '\n'

        output_html += '<tr>\n'
        output_html += '<th>Peaks</th>\n'
        output_html += '<th>Negative Peaks</th>\n'
        output_html += '<th>R Model</th>\n'
        output_html += '</tr>\n'
        output_html += '\n'

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
                # Python str key and Python list of Python list objects
                # of PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_replicate_keys = c_replicate_dict.keys()
                        c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                        for c_replicate_key in c_replicate_keys:

                            # prefix = 'chipseq_macs14_{}__{}'.format(t_replicate_key, c_replicate_key)

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

                                    output_hub += 'track ChIP_{}__{}_{}_{}\n'. \
                                        format(t_replicate_key, c_replicate_key, state, scaling)
                                    output_hub += 'type bigWig\n'
                                    output_hub += 'shortLabel ChIP_{}__{}_{}_{}\n'. \
                                        format(t_replicate_key, c_replicate_key, state, scaling)
                                    output_hub += 'longLabel {} ChIP-Seq read counts for {} of {} versus {}\n'. \
                                        format(scaling.capitalize(), state, t_replicate_key, c_replicate_key)
                                    output_hub += 'bigDataUrl '
                                    output_hub += '{}__{}_MACS_wiggle/{}/{}__{}_{}_afterfiting_all.bw\n'. \
                                        format(t_replicate_key, c_replicate_key, state,
                                               t_replicate_key, c_replicate_key, state)
                                    if treatment and not absolute:
                                        output_hub += 'visibility full\n'
                                    else:
                                        output_hub += 'visibility hide\n'
                                        # 'html' is missing from the common settings.

                                    # Common optional settings.

                                    output_hub += 'color {}\n'. \
                                        format(defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                                    # bigWig - Signal graphing track settings.

                                    output_hub += 'graphTypeDefault bar\n'
                                    output_hub += 'maxHeightPixels 100:60:20\n'
                                    output_hub += 'smoothingWindow off\n'

                                    if absolute:
                                        # Track with absolute scaling.
                                        output_hub += 'autoScale on\n'
                                    else:
                                        # Track with relative scaling.
                                        output_hub += 'autoScale off\n'
                                        output_hub += 'viewLimits 0:40\n'

                                    output_hub += '\n'

                            # Add a UCSC trackDB entry for each bigBed peaks file.

                            output_hub += 'track Peaks_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'type bigBed\n'
                            output_hub += 'shortLabel Peaks_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'longLabel ChIP-Seq peaks for {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'bigDataUrl {}__{}_peaks.bb\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'visibility pack\n'
                            # 'html' is missing from the common settings.

                            # Common optional settings.
                            output_hub += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                            output_hub += '\n'

                            # Add web content.

                            output_html += '<tr>\n'

                            # if treatment and absolute:
                            # output += '<td><strong>{}</strong></td>\n'.format(t_replicate_key)
                            # if not treatment and absolute:
                            # output += '<td><strong>{}</strong></td>\n'.format(c_replicate_key)

                            output_html += '<td><a href="{}__{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                format(t_replicate_key, c_replicate_key,
                                       t_replicate_key, c_replicate_key)

                            output_html += '<td>'
                            output_html += '<a href="{}__{}_negative_peaks.xls">Negative peaks {} versus {}</a>'. \
                                format(t_replicate_key, c_replicate_key,
                                       t_replicate_key, c_replicate_key)
                            output_html += '</td>\n'

                            output_html += '<td><a href="{}__{}_model.r">R model</a></td>\n'. \
                                format(t_replicate_key, c_replicate_key,
                                       t_replicate_key, c_replicate_key)

                            output_html += '</tr>\n'
                            output_html += '\n'

        output_html += '</table>\n'
        output_html += '\n'

        # Add UCSC trackDB entries for each Bowtie2 BAM file.

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
            # Python str key and Python list of Python list objects
            # of PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                # Add a UCSC trackDB entry.

                # Common trackDb settings.

                output_hub += 'track Alignment_{}\n'.format(replicate_key)
                output_hub += 'type bam\n'
                output_hub += 'shortLabel Alignment_{}\n'.format(replicate_key)
                output_hub += 'longLabel Bowtie2 alignment of {}\n'. \
                    format(replicate_key)
                output_hub += 'bigDataUrl {}.aligned.sorted.bam\n'. \
                    format(replicate_key)
                output_hub += 'visibility hide\n'
                # TODO: The 'html' option is missing.

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                output_hub += '\n'

                # TODO: Including FastQC results would require rearranging the above HTML code.
                # output += '{}__{}_fastqc/fastqc_report.html'. \
                # format(t_replicate_key, c_replicate_key)

        self.report_to_file(content=output_html)
        self.ucsc_hub_to_file(content=output_hub)

        return

    def report(self):
        """Create a C{bsf.analyses.chip_seq.ChIPSeq} report in HTML format and a UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # contrast_field_names = ["", "Group1", "Members1", "Group2", "Members2", "DB.edgeR"]

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        output_hub = str()

        # Write a HTML document.

        output_html = str()

        output_html += '<h1 id="chip_seq_analysis">{} {}</h1>\n'.format(self.project_name, self.name)
        output_html += '\n'

        output_html += '<p>\n'
        output_html += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
        output_html += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></strong>,\n'
        output_html += 'before peaks are called with '
        output_html += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 2</a>\n'
        output_html += 'on an enrichment and background sample pair.\n'
        output_html += '</p>\n'
        output_html += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output_html += '<p id="ucsc_track_hub">\n'
        output_html += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
        output_html += 'in their genomic context via the project-specific\n'
        output_html += 'UCSC Genome Browser Track Hub <a href="{}">{}</a>.\n'. \
            format(self.ucsc_track_url(options_dict=options_dict),
                   self.project_name)
        output_html += '</p>\n'
        output_html += '\n'

        output_html += '<table id="peak_calling_table">\n'
        output_html += '\n'

        output_html += '<tr>\n'
        output_html += '<th>Peaks</th>\n'
        output_html += '<th>R Model</th>\n'
        output_html += '</tr>\n'
        output_html += '\n'

        # Group via UCSC super tracks.

        output_hub += 'track Alignment\n'
        output_hub += 'shortLabel QC Alignment\n'
        output_hub += 'longLabel ChIP Read Alignment\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group alignment\n'
        output_hub += '\n'

        output_hub += 'track Background\n'
        output_hub += 'shortLabel QC Background\n'
        output_hub += 'longLabel ChIP Background Signal\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group background\n'
        output_hub += '\n'

        output_hub += 'track Enrichment\n'
        output_hub += 'shortLabel QC Enrichment\n'
        output_hub += 'longLabel ChIP Enrichment Signal\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group enrichment\n'
        output_hub += '\n'

        output_hub += 'track Comparison\n'
        output_hub += 'shortLabel ChIP Intensity\n'
        output_hub += 'longLabel ChIP Intensity\n'
        output_hub += 'visibility full\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group comparison\n'
        output_hub += '\n'

        output_hub += 'track Peaks\n'
        output_hub += 'shortLabel ChIP Peaks\n'
        output_hub += 'longLabel ChIP Peaks\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group peaks\n'
        output_hub += '\n'

        output_hub += 'track Summits\n'
        output_hub += 'shortLabel ChIP Summits\n'
        output_hub += 'longLabel ChIP Summits\n'
        output_hub += 'visibility hide\n'
        output_hub += 'superTrack on\n'
        output_hub += 'group summits\n'
        output_hub += '\n'

        # Group via UCSC composite tracks.

        composite_groups = dict()

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]

            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
                # Python str key and Python list of Python list objects
                # of PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                        c_replicate_keys = c_replicate_dict.keys()
                        c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                        for c_replicate_key in c_replicate_keys:

                            prefix = 'chipseq_macs2_{}__{}'.format(t_replicate_key, c_replicate_key)

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
                                output_hub += 'track {}\n'.format(composite_group)
                                output_hub += 'type bigWig\n'
                                output_hub += 'shortLabel {}\n'.format(composite_group)
                                output_hub += 'longLabel ChIP background signal for factor {!r}\n'. \
                                    format(factor)
                                output_hub += 'visibility dense\n'
                                output_hub += 'compositeTrack on\n'
                                output_hub += 'parent Background\n'
                                output_hub += 'allButtonPair on\n'
                                output_hub += 'centerLabelsDense on\n'
                                output_hub += '\n'

                            # Common trackDb settings.

                            output_hub += 'track {}__{}_background\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                            output_hub += 'type bigWig\n'
                            output_hub += 'shortLabel {}__{}_bkg\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'longLabel ChIP background signal {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'bigDataUrl {}/{}_control_lambda.bw\n'. \
                                format(prefix, prefix)
                            output_hub += 'visibility dense\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            output_hub += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            output_hub += 'alwaysZero off\n'
                            output_hub += 'autoScale off\n'
                            output_hub += 'graphTypeDefault bar\n'
                            output_hub += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            output_hub += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            output_hub += 'viewLimits 0:15\n'
                            output_hub += 'viewLimitsMax 0:40\n'
                            output_hub += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            output_hub += 'parent {} on\n'.format(composite_group)
                            output_hub += 'centerLabelsDense off\n'
                            output_hub += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_treat_pileup.bw file.
                            #

                            composite_group = '{}_enrichment'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                output_hub += 'track {}\n'.format(composite_group)
                                output_hub += 'type bigWig\n'
                                output_hub += 'shortLabel {}_enrichment\n'.format(factor)
                                output_hub += 'longLabel ChIP enrichment signal for factor {!r}\n'. \
                                    format(factor)
                                output_hub += 'visibility dense\n'
                                output_hub += 'compositeTrack on\n'
                                output_hub += 'parent Enrichment\n'
                                output_hub += 'allButtonPair on\n'
                                output_hub += 'centerLabelsDense on\n'
                                output_hub += '\n'

                            # Common trackDb settings.

                            output_hub += 'track {}__{}_enrichment\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            output_hub += 'type bigWig\n'
                            output_hub += 'shortLabel {}__{}_enr\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'longLabel ChIP enrichment signal {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'bigDataUrl {}/{}_treat_pileup.bw\n'. \
                                format(prefix, prefix)
                            output_hub += 'visibility dense\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            output_hub += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            output_hub += 'alwaysZero off\n'
                            output_hub += 'autoScale off\n'
                            output_hub += 'graphTypeDefault bar\n'
                            output_hub += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            output_hub += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            output_hub += 'viewLimits 0:15\n'
                            output_hub += 'viewLimitsMax 0:40\n'
                            output_hub += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            output_hub += 'parent {} on\n'.format(composite_group)
                            output_hub += 'centerLabelsDense off\n'
                            output_hub += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_bdgcmp.bw file.
                            #

                            composite_group = '{}_intensity'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                output_hub += 'track {}\n'.format(composite_group)
                                output_hub += 'type bigWig\n'
                                output_hub += 'shortLabel {}\n'.format(composite_group)
                                output_hub += 'longLabel ChIP intensity for factor {!r}\n'. \
                                    format(factor)
                                output_hub += 'visibility full\n'
                                output_hub += 'compositeTrack on\n'
                                output_hub += 'parent Comparison\n'
                                output_hub += 'allButtonPair on\n'
                                output_hub += 'centerLabelsDense on\n'
                                output_hub += '\n'

                            # Common trackDb settings.

                            output_hub += 'track {}__{}_intensity\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            output_hub += 'type bigWig\n'
                            output_hub += 'shortLabel {}__{}_int\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'longLabel ChIP intensity {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'bigDataUrl {}/{}_bdgcmp.bw\n'. \
                                format(prefix, prefix)
                            output_hub += 'visibility full\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            output_hub += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            output_hub += 'alwaysZero off\n'
                            output_hub += 'autoScale off\n'
                            output_hub += 'graphTypeDefault bar\n'
                            output_hub += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            output_hub += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            output_hub += 'viewLimits 0:15\n'
                            output_hub += 'viewLimitsMax 0:40\n'
                            output_hub += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            output_hub += 'parent {} on\n'.format(composite_group)
                            output_hub += 'centerLabelsDense off\n'
                            output_hub += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_peaks.bb file.
                            #

                            # Common trackDb settings.

                            output_hub += 'track {}__{}_peaks\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'type bigBed\n'
                            output_hub += 'shortLabel {}__{}_peaks\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'longLabel ChIP peaks {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'bigDataUrl {}/{}_peaks.bb\n'. \
                                format(prefix, prefix)
                            output_hub += 'visibility pack\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            output_hub += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # bigBed - Item or region track settings.

                            # Supertrack settings.

                            output_hub += 'parent Peaks\n'
                            output_hub += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_summits.bb file.
                            #

                            # Common trackDb settings.

                            output_hub += 'track {}__{}_summits\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'type bigBed\n'
                            output_hub += 'shortLabel {}__{}_summits\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'longLabel ChIP summits {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            output_hub += 'bigDataUrl {}/{}_summits.bb\n'. \
                                format(prefix, prefix)
                            output_hub += 'visibility pack\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            output_hub += 'color {}\n'. \
                                format(defaults.web.get_chipseq_colour(factor=factor))

                            # Supertrack settings.

                            output_hub += 'parent Summits\n'
                            output_hub += '\n'

                            #
                            # Add web content.
                            #

                            output_html += '<tr>\n'

                            output_html += '<td><a href="{}/{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                format(prefix, prefix, t_replicate_key, c_replicate_key)

                            output_html += '<td><a href="{}/{}_model.r">R model</a></td>\n'. \
                                format(prefix, prefix)

                            output_html += '</tr>\n'
                            output_html += '\n'

        output_html += '</table>\n'
        output_html += '\n'

        # Add UCSC trackDB entries for each Bowtie2 BAM file.

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # bsf.data.Sample.get_all_paired_reads() returns a Python dict of
            # Python str key and Python list of Python list objects
            # of PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                #
                # Add a UCSC trackDB entry for each NAME.aligned.sorted.bam file.
                #

                # Common trackDb settings.

                output_hub += 'track {}_alignment\n'.format(replicate_key)
                output_hub += 'type bam\n'
                output_hub += 'shortLabel {}_alignment\n'.format(replicate_key)
                output_hub += 'longLabel {} ChIP read alignment\n'. \
                    format(replicate_key)
                output_hub += 'bigDataUrl chipseq_bowtie2_{}/{}.aligned.sorted.bam\n'. \
                    format(replicate_key, replicate_key)
                output_hub += 'visibility hide\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                # Supertrack settings.

                output_hub += 'parent Alignment\n'
                output_hub += '\n'

        # Differential binding analysis.

        output_html += '<h2 id="differential_binding">Differential Binding Analysis</h2>\n'
        output_html += '\n'

        output_html += '<table id="differential_binding_table">\n'
        output_html += '\n'

        output_html += '<tr>\n'
        output_html += '<th>Factor and Contrast</th>\n'
        output_html += '<th>Correlation Peak Caller</th>\n'
        output_html += '<th>Correlation Peak Counts</th>\n'
        output_html += '<th>Correlation Analysis</th>\n'
        output_html += '<th>MA Plot</th>\n'
        output_html += '<th>Scatter Plot</th>\n'
        output_html += '<th>PCA Plot</th>\n'
        output_html += '<th>Box Plot</th>\n'
        output_html += '<th>DiffBind Report</th>\n'
        output_html += '</tr>\n'

        keys = self._factor_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            prefix = 'chipseq_diffbind_{}'.format(key)

            output_html += '<tr>\n'

            # ChIP Factor

            output_html += '<td><strong>{}</strong></td>\n'.format(key)

            # Correlation heat map of peak caller scores.

            output_html += '<td>'
            output_html += '<a href="{}/{}_correlation_peak_caller_score.png">'. \
                format(prefix, prefix)
            output_html += '<img alt="DiffBind correlation analysis for factor {}" ' \
                           'src="{}/{}_correlation_peak_caller_score.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            output_html += '</a>'
            output_html += '</td>\n'

            # Correlation heat map of counts.

            output_html += '<td>'
            output_html += '<a href="{}/{}_correlation_read_counts.png">'. \
                format(prefix, prefix)
            output_html += '<img alt="DiffBind correlation analysis for factor {}" ' \
                           'src="{}/{}_correlation_read_counts.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            output_html += '</a>'
            output_html += '</td>\n'

            # Correlation heat map of differential binding analysis.

            output_html += '<td>'
            output_html += '<a href="{}/{}_correlation_analysis.png">'. \
                format(prefix, prefix)
            output_html += '<img alt="DiffBind correlation analysis for factor {}" ' \
                           'src="{}/{}_correlation_analysis.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            output_html += '</a>'
            output_html += '</td>\n'

            output_html += '<td></td>\n'  # MA Plot
            output_html += '<td></td>\n'  # Scatter Plot
            output_html += '<td></td>\n'  # PCA Plot
            output_html += '<td></td>\n'  # Box Plot
            output_html += '<td></td>\n'  # DiffBin Report

            output_html += '</tr>\n'

            # Read the file of contrasts ...

            file_path = os.path.join(self.genome_directory, prefix, '{}_contrasts.csv'.format(prefix))

            if not os.path.exists(path=file_path):
                warnings.warn(
                    'File {!r} does not exist.'.format(file_path),
                    UserWarning)
                continue

            sas = AnnotationSheet.from_file_path(file_path=file_path)

            for row_dict in sas.row_dicts:
                suffix = '{}__{}'.format(row_dict['Group1'], row_dict['Group2'])

                output_html += '<tr>\n'

                output_html += '<td>{}</td>\n'.format(suffix)
                output_html += '<td></td>\n'  # Correlation heat map of peak caller scores.
                output_html += '<td></td>\n'  # Correlation heat map of counts.
                output_html += '<td></td>\n'  # Correlation heat map of differential binding analysis.

                # MA Plot

                output_html += '<td>'
                output_html += '<a href="{}/{}_ma_plot_{}.png">'.format(prefix, prefix, suffix)
                output_html += '<img alt="DiffBind MA plot for factor {}" ' \
                               'src="{}/{}_ma_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output_html += '</a>'
                output_html += '</td>\n'

                # Scatter Plot

                output_html += '<td>'
                output_html += '<a href="{}/{}_scatter_plot_{}.png">'.format(prefix, prefix, suffix)
                output_html += '<img alt="DiffBind scatter plot for factor {}" ' \
                               'src="{}/{}_scatter_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output_html += '</a>'
                output_html += '</td>\n'

                # Principal Component Analysis Plot

                output_html += '<td>'
                output_html += '<a href="{}/{}_pca_plot_{}.png">'.format(prefix, prefix, suffix)
                output_html += '<img alt="DiffBind PCA plot for factor {}" ' \
                               'src="{}/{}_pca_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output_html += '</a>'
                output_html += '</td>\n'

                # Box Plot

                output_html += '<td>'
                output_html += '<a href="{}/{}_box_plot_{}.png">'.format(prefix, prefix, suffix)
                output_html += '<img alt="DiffBind Box plot for factor {}" ' \
                               'src="{}/{}_box_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output_html += '</a>'
                output_html += '</td>\n'

                # DiffBind report

                output_html += '<td><a href="{}/DBA_{}_report_{}.csv">DBA_{}_report_{}</a></td>\n'. \
                    format(prefix, key, suffix, key, suffix)

                output_html += '</tr>\n'

        output_html += '\n'
        output_html += '</table>\n'
        output_html += '\n'

        self.report_to_file(content=output_html)
        self.ucsc_hub_to_file(content=output_hub)

        return
