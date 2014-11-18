"""Bio.BSF.Analysis

A package of classes and methods supporting NGS-specific analyses such as ChIP-Seq or RNA-Seq.
"""

#
# Copyright 2013 Michael K. Schuster
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
import string
import warnings

from Bio.BSF import Analysis, Configuration, Default, Defaults, DRMS, Executable
from Bio.BSF.annotation import AnnotationSheet, SampleAnnotationSheet
from Bio.BSF.Data import Collection, ProcessedRunFolder, Sample
from Bio.BSF.Executables import Bowtie2, Macs14, Macs2Bdgcmp, Macs2Callpeak, FastQC


class ChIPSeqComparison(object):
    def __init__(self, c_name, t_name, c_samples, t_samples,
                 factor, tissue=None, condition=None, treatment=None, replicate=None):
        """Initialise a C{ChIPSeqComparison} object.

        @param c_name: Control name
        @type c_name: str
        @param t_name: Treatment name
        @type t_name: str
        @param c_samples: Python C{list} of control C{Sample} objects
        @type c_samples: list
        @param t_samples: Python C{list} of treatment C{Sample} objects
        @type t_samples:list
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
        @return: Nothing
        @rtype: None
        """

        # Condition', 'Treatment', 'Replicate',
        # 'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat'

        if c_name:
            self.c_name = c_name
        else:
            self.c_name = str()

        if t_name:
            self.t_name = t_name
        else:
            self.t_name = str()

        if c_samples:
            self.c_samples = c_samples
        else:
            self.c_samples = list()

        if t_samples:
            self.t_samples = t_samples
        else:
            self.t_samples = list()

        if factor:
            self.factor = factor
        else:
            self.factor = str()

        if tissue:
            self.tissue = tissue
        else:
            self.tissue = str()

        if condition:
            self.condition = condition
        else:
            self.condition = str()

        if treatment:
            self.treatment = treatment
        else:
            self.treatment = str()

        if replicate:
            self.replicate = replicate
        else:
            self.replicate = int(x=0)


class ChIPSeqDiffBindSheet(AnnotationSheet):
    """ChIP-Seq Bioconductor DiffBind annotation sheet class.

    Attributes:
    @cvar _file_type: File type (i.e. I{excel} or I{excel-tab} defined in the csv module Dialect class)
    @type _file_type: str
    @cvar _header_line: Header line exists
    @type _header_line: bool
    @cvar _field_names: Python C{list} of Python C{str} (field name) objects
    @type _field_names: list
    @cvar _test_methods: Python C{dict} of Python C{str} (field name) key data and
        Python C{list} of Python C{function} value data
    @type _test_methods: dict
    """

    _file_type = 'excel'

    _header_line = True

    _field_names = [
        'SampleID', 'Tissue', 'Factor', 'Condition', 'Treatment', 'Replicate',
        'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat'
    ]

    _test_methods = dict(
        SampleID=[
            AnnotationSheet.check_alphanumeric
        ],
        Tissue=[
            AnnotationSheet.check_alphanumeric
        ],
        Factor=[
            AnnotationSheet.check_alphanumeric
        ],
        Condition=[
            AnnotationSheet.check_alphanumeric
        ],
        Treatment=[
            AnnotationSheet.check_alphanumeric
        ],
        Replicate=[
            AnnotationSheet.check_numeric
        ],
        ControlID=[
            AnnotationSheet.check_alphanumeric
        ],
        PeakCaller=[
            AnnotationSheet.check_alphanumeric
        ],
        PeakFormat=[
            AnnotationSheet.check_alphanumeric
        ]
    )

    def sort(self):
        """Sort by columns I{Tissue}, I{Factor}, I{Condition}, I{Treatment} and I{Replicate}.
        """
        self.row_dicts.sort(
            cmp=lambda x, y:
            cmp(x['Tissue'], y['Tissue']) or
            cmp(x['Factor'], y['Factor']) or
            cmp(x['Condition'], y['Condition']) or
            cmp(x['Treatment'], y['Treatment']) or
            cmp(int(x['Replicate']), int(y['Replicate'])))

    def write_to_file(self):
        """Write a C{ChIPSeqDiffBindSheet} to a file.
        """

        # Override the method from the super-class to automatically sort before writing to a file.

        self.sort()
        super(ChIPSeqDiffBindSheet, self).write_to_file()


class ChIPSeq(Analysis):
    """ChIP-Seq-specific C{Analysis} sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new C{ChIPSeq} object from a UNIX-style configuration file via the C{Configuration} class.

        @param config_file: UNIX-style configuration file
        @type config_file: str | unicode
        @return: C{ChIPSeq}
        @rtype: ChIPSeq
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new C{ChIPSeq} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{ChIPSeq}
        @rtype: ChIPSeq
        """

        assert isinstance(configuration, Configuration)

        chipseq = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.ChIPSeq" section specifies defaults for this BSF Analysis.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        chipseq.set_Configuration(chipseq.configuration, section=section)

        return chipseq

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 cmp_file=None):
        """Initialise a C{ChIPSeq} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{Analysis}-wide project directory,
            normally under the C{Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{Analysis}-wide genome directory,
            normally under the C{Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param drms_list: Python C{list} of C{DRMS} objects
        @type drms_list: list
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{list} objects of C{Sample} objects
        @type comparisons: dict
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        @param cmp_file: Comparison file
        @type cmp_file: str | unicode
        @return: Nothing
        @rtype: None
        """

        super(ChIPSeq, self).__init__(configuration=configuration,
                                      project_name=project_name, genome_version=genome_version,
                                      input_directory=input_directory, output_directory=output_directory,
                                      project_directory=project_directory, genome_directory=genome_directory,
                                      e_mail=e_mail, debug=debug, drms_list=drms_list,
                                      collection=collection, comparisons=comparisons, samples=samples)

        # Sub-class specific ...

        if cmp_file:
            self.cmp_file = cmp_file
        else:
            self.cmp_file = str()

    def set_Configuration(self, configuration, section):
        """Set instance variables of a C{ChIPSeq} object via a section of a
        C{Configuration} object.

        Instance variables without a
        configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        super(ChIPSeq, self).set_Configuration(configuration=configuration, section=section)

        # Read a comparison file.

        if configuration.config_parser.has_option(section=section, option='cmp_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')

    def _read_comparisons(self, cmp_file):
        """Read a C{SampleAnnotationSheet} CSV file from disk.

        Column headers for CASAVA folders:
            Treatment/Control ProcessedRunFolder:
                CASAVA processed run folder name or
                C{Analysis.input_directory} by default
            Treatment/Control Project:
                CASAVA Project name or
                C{Analysis.project_name} by default
            Treatment/Control Sample:
                CASAVA Sample name, no default
        Column headers for independent samples:
            Treatment/Control Sample
            Treatment/Control File
            Treatment/Control Group
        @param cmp_file: Comparison file path
        @type cmp_file: str | unicode
        @return: Nothing
        @rtype: None
        """

        if self.debug > 1:
            print '{!r} method _read_comparisons:'.format(self)

        sas = SampleAnnotationSheet.read_from_file(file_path=cmp_file)

        # Unfortunately, two passes through the comparison sheet are required.
        # In the first one merge all BSF Sample objects that share the name.
        # Merging BSF Sample objects is currently the only way to pool BSF PairedReads objects,
        # which is required for ChIP-Seq experiments.

        sample_dict = dict()

        # First pass, merge BSF Sample objects, if they have the same name.
        for row_dict in sas.row_dicts:
            for prefix in ('Control', 'Treatment'):
                name, samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                for o_sample in samples:
                    if o_sample.name not in sample_dict:
                        sample_dict[o_sample.name] = o_sample
                    else:
                        n_sample = Sample.from_Samples(sample1=sample_dict[o_sample.name], sample2=o_sample)
                        sample_dict[n_sample.name] = n_sample

        # Second pass, add all BSF Sample objects mentioned in a comparison.
        level1_dict = dict()

        for row_dict in sas.row_dicts:

            c_name, c_samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix='Control')
            t_name, t_samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix='Treatment')

            # ChIP-Seq experiments use the order treatment versus control in comparisons.
            comparison_key = '{}__{}'.format(t_name, c_name)

            # For a successful comparison, both Python list objects of BSF Sample objects have to be defined.

            if not (len(t_samples) and len(c_samples)):
                if self.debug > 1:
                    print 'Redundant comparison line with Treatment {!r} samples {} and Control {!r} samples {}'. \
                        format(t_name, len(t_samples), c_name, len(c_samples))
                continue

            # Add all control BSF Sample or BSF SampleGroup objects to the Sample list.

            for c_sample in c_samples:
                if self.debug > 1:
                    print '  Control Sample name: {!r} file_path:{!r}'.format(c_sample.name, c_sample.file_path)
                    print c_sample.trace(1)
                    # Find the BSF Sample in the unified sample dictionary.
                if c_sample.name in sample_dict:
                    self.add_Sample(sample=sample_dict[c_sample.name])

            # Add all treatment BSF Sample or BSF SampleGroup objects to the Sample list.

            for t_sample in t_samples:
                if self.debug > 1:
                    print '  Treatment Sample name: {!r} file_path:{!r}'.format(t_sample.name, t_sample.file_path)
                    print t_sample.trace(1)
                if t_sample.name in sample_dict:
                    self.add_Sample(sample=sample_dict[t_sample.name])

            if 'Tissue' in row_dict:
                tissue = row_dict['Tissue']
            else:
                tissue = str()

            if 'Factor' in row_dict:
                factor = row_dict['Factor']
            else:
                factor = Defaults.web.chipseq_default_factor

            if 'Condition' in row_dict:
                condition = row_dict['Condition']
            else:
                condition = str()

            if 'Treatment' in row_dict:
                treatment = row_dict['Treatment']
            else:
                treatment = str()

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

            self.comparisons[comparison_key] = ChIPSeqComparison(c_name=c_name,
                                                                 t_name=t_name,
                                                                 c_samples=c_samples,
                                                                 t_samples=t_samples,
                                                                 factor=factor,
                                                                 tissue=tissue,
                                                                 condition=condition,
                                                                 treatment=treatment,
                                                                 replicate=0)

        # Sort the comparison keys alphabetically and assign replicate numbers into ChiPSeqComparison objects.

        for key1 in level1_dict.keys():

            level2_dict = level1_dict[key1]

            keys2 = level2_dict.keys()
            keys2.sort(cmp=lambda x, y: cmp(x, y))

            i = 1
            for key2 in keys2:
                level2_dict[key2] = i
                self.comparisons[key2].replicate = i
                i += 1

    def run(self):
        """Run this C{ChIPSeq} C{Analysis}.
        """

        super(ChIPSeq, self).run()

        # ChIPSeq requires a genome version.

        if not self.genome_version:
            raise Exception('A ChIPSeq analysis requires a genome_version configuration option.')

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file) and not os.path.exists(path=self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of BSF Sample objects by the BSF Sample name.
        # This cannot be done in the super-class, because BSF Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_bowtie2_jobs()
        # self._create_Macs14_jobs()
        self._create_macs2_jobs()
        self._create_diffbind_jobs()

    def _create_bowtie2_jobs(self):
        """Create Bowtie2 alignment jobs.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        # Get the Bowtie2 index.

        bowtie2_index = os.path.join(Default.absolute_genomes(self.genome_version),
                                     'forBowtie2', self.genome_version)

        # Initialise the Distributed Resource Management System objects for Bowtie2.

        bowtie2_drms = DRMS.from_Analysis(name='bowtie2',
                                          work_directory=self.genome_directory,
                                          analysis=self)
        self.drms_list.append(bowtie2_drms)

        # Use the bsf_sam2bam.sh script to convert aligned SAM into
        # aligned, sorted, indexed BAM files.

        sam2bam_drms = DRMS.from_Analysis(name='sam2bam',
                                          work_directory=self.genome_directory,
                                          analysis=self)
        self.drms_list.append(sam2bam_drms)

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                bowtie2 = Bowtie2(name='chipseq_bowtie2_{}'.format(replicate_key),
                                  analysis=self)
                bowtie2_drms.add_Executable(bowtie2)

                sam2bam = Executable.from_Analysis(name='chipseq_sam2bam_{}'.format(replicate_key),
                                                   program='bsf_sam2bam.sh',
                                                   analysis=self)
                sam2bam_drms.add_Executable(sam2bam)

                sam2bam.dependencies.append(bowtie2.name)

                # Set Bowtie2 options.

                bowtie2.add_OptionShort(key='x', value=bowtie2_index)

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

                if len(reads1) and not len(reads2):
                    bowtie2.add_OptionShort(key='U', value=string.join(words=reads1, sep=','))
                elif len(reads1) and len(reads2):
                    bowtie2.add_OptionShort(key='1', value=string.join(words=reads1, sep=','))
                if len(reads2):
                    bowtie2.add_OptionShort(key='2', value=string.join(words=reads2, sep=','))

                # TODO: The following options are properties of the Sample,
                # PairedReads and Reads objects.
                # bowtie2.add_SwitchShort(key='q')
                # bowtie2.add_SwitchShort(key='phred33')

                # TODO: It would be good to have code that parses the Illumina Run Info and writes
                # this information into the CSV file.

                # TODO: Andreas' original implementation on Bowtie1 sets -a,
                # which may not be required for Bowtie2.
                # See description of option -k in the bowtie2 manual.
                # TODO: To avoid repeats we may want to set -a?

                bowtie2.add_OptionLong(key='threads', value=str(bowtie2_drms.threads))

                # Put all sample-specific information into a sub-directory.

                replicate_directory = os.path.join(self.genome_directory,
                                                   'chipseq_bowtie2_{}'.format(replicate_key))

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
                                                   replicate_key + '.aligned.sam')

                # Set bsf_sam2bam.sh options.

                sam2bam.arguments.append(os.path.join(self.genome_directory,
                                                      replicate_directory,
                                                      replicate_key))

    def _create_macs14_jobs(self):
        """Create MACS14 peak caller jobs.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        genome_sizes = config_parser.get(section=config_section, option='genome_sizes')
        genome_sizes = os.path.expanduser(genome_sizes)
        genome_sizes = os.path.expandvars(genome_sizes)

        macs14_drms = DRMS.from_Analysis(name='macs14',
                                         work_directory=self.genome_directory,
                                         analysis=self)

        self.drms_list.append(macs14_drms)

        process_macs14_drms = DRMS.from_Analysis(name='process_macs14',
                                                 work_directory=self.genome_directory,
                                                 analysis=self)

        self.drms_list.append(process_macs14_drms)

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)

            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                # Python str key and Python list of Python list objects
                # of Bio.BSF.Data.PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                        c_replicate_keys = c_replicate_dict.keys()
                        c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                        for c_replicate_key in c_replicate_keys:

                            macs14 = Macs14(name='chipseq_macs14_{}__{}'.format(t_replicate_key, c_replicate_key),
                                            analysis=self)

                            macs14.dependencies.append('chipseq_sam2bam_' + t_replicate_key)
                            macs14.dependencies.append('chipseq_sam2bam_' + c_replicate_key)

                            macs14_drms.add_Executable(macs14)

                            process_macs14 = Executable.from_Analysis(
                                name='chipseq_process_macs14_{}__{}'.format(t_replicate_key, c_replicate_key),
                                program='bsf_chipseq_process_macs14.sh',
                                analysis=self)

                            process_macs14.dependencies.append(macs14.name)

                            process_macs14_drms.add_Executable(process_macs14)

                            # Set macs14 options.

                            macs14.add_OptionLong(key='treatment',
                                                  value=os.path.join(self.genome_directory,
                                                                     'chipseq_bowtie2_{}'.format(t_replicate_key),
                                                                     '{}.aligned.sorted.bam'.format(t_replicate_key)))
                            macs14.add_OptionLong(key='control',
                                                  value=os.path.join(self.genome_directory,
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

                            macs14.add_OptionLong(key='name', value=os.path.join('.', prefix, prefix))

                            # This option has to be specified via the configuration.ini file.
                            # macs14.add_OptionLong(key='gsize', value=genome_sizes)
                            macs14.add_SwitchLong(key='single-profile')
                            macs14.add_SwitchLong(key='call-subpeaks')
                            macs14.add_SwitchLong(key='wig')

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
                                macs14.add_SwitchLong(key='nomodel')
                                macs14.add_OptionLong(key='shiftsize', value='73')
                                macs14.add_OptionLong(key='pvalue', value='1e-3')
                            elif factor == 'H3K56AC':
                                pass
                            elif factor == 'H4K16AC':
                                pass
                            elif factor == 'OTHER':
                                pass
                            else:
                                warnings.warn(
                                    'Unable to set MACS14 parameters for unknown factor {!r}.\n'
                                    'Please use default factor {!r} or adjust Python code if necessary.'.
                                    format(factor, Defaults.web.chipseq_default_factor),
                                    UserWarning)

                            # Set macs14 arguments.

                            # Set process_macs14 options.
                            # Set process_macs14 arguments.

                            # Specify the output path as in the macs14 --name option.
                            process_macs14.arguments.append(os.path.join('.', prefix, prefix))
                            process_macs14.arguments.append(genome_sizes)

    def _create_macs2_jobs(self):
        """Create MACS2 peak caller jobs.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        genome_sizes = config_parser.get(section=config_section, option='genome_sizes')
        genome_sizes = os.path.expanduser(genome_sizes)
        genome_sizes = os.path.expandvars(genome_sizes)

        macs2_callpeak_drms = DRMS.from_Analysis(name='macs2_callpeak',
                                                 work_directory=self.genome_directory,
                                                 analysis=self)
        self.drms_list.append(macs2_callpeak_drms)

        macs2_bdgcmp_drms = DRMS.from_Analysis(name='macs2_bdgcmp',
                                               work_directory=self.genome_directory,
                                               analysis=self)
        self.drms_list.append(macs2_bdgcmp_drms)

        process_macs2_drms = DRMS.from_Analysis(name='process_macs2',
                                                work_directory=self.genome_directory,
                                                analysis=self)
        self.drms_list.append(process_macs2_drms)

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)

            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                # Python str key and Python list of Python list objects
                # of Bio.BSF.Data.PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                        c_replicate_keys = c_replicate_dict.keys()
                        c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                        for c_replicate_key in c_replicate_keys:

                            macs2_callpeak = Macs2Callpeak(
                                name='chipseq_macs2_callpeak_{}__{}'.format(t_replicate_key, c_replicate_key),
                                analysis=self)
                            macs2_callpeak_drms.add_Executable(macs2_callpeak)
                            macs2_callpeak.dependencies.append('chipseq_sam2bam_' + t_replicate_key)
                            macs2_callpeak.dependencies.append('chipseq_sam2bam_' + c_replicate_key)

                            macs2_bdgcmp = Macs2Bdgcmp(
                                name='chipseq_macs2_bdgcmp_{}__{}'.format(t_replicate_key, c_replicate_key),
                                analysis=self)
                            macs2_bdgcmp_drms.add_Executable(macs2_bdgcmp)
                            macs2_bdgcmp.dependencies.append(macs2_callpeak.name)

                            process_macs2 = Executable.from_Analysis(
                                name='chipseq_process_macs2_{}__{}'.format(t_replicate_key, c_replicate_key),
                                program='bsf_chipseq_process_macs2.sh',
                                analysis=self)
                            process_macs2_drms.add_Executable(process_macs2)
                            process_macs2.dependencies.append(macs2_bdgcmp.name)

                            # Set MACS2 callpeak sub-command options.

                            mc2 = macs2_callpeak.sub_command

                            mc2.add_OptionLong(key='treatment',
                                               value=os.path.join(self.genome_directory,
                                                                  'chipseq_bowtie2_{}'.format(t_replicate_key),
                                                                  '{}.aligned.sorted.bam'.format(t_replicate_key)))

                            mc2.add_OptionLong(key='control',
                                               value=os.path.join(self.genome_directory,
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

                            mc2.add_OptionLong(key='name', value=os.path.join(prefix, prefix))

                            # The 'gsize' option has to be specified via the configuration.ini file.
                            # macs2_callpeak.add_OptionLong(key='gsize', value=genome_sizes

                            mc2.add_SwitchLong(key='bdg')
                            mc2.add_SwitchLong(key='SPMR')

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
                                mc2.add_SwitchLong(key='nomodel')
                                mc2.add_OptionLong(key='shiftsize', value='73')
                                mc2.add_OptionLong(key='pvalue', value='1e-3')
                            elif factor == 'H3K56AC':
                                pass
                            elif factor == 'H4K16AC':
                                pass
                            elif factor == 'OTHER':
                                pass
                            else:
                                warnings.warn(
                                    'Unable to set MACS2 parameters for unknown factor {!r}.\n'
                                    'Please use default factor {!r} or adjust Python code if necessary.'.
                                    format(factor, Defaults.web.chipseq_default_factor),
                                    UserWarning)

                            # Set macs2_callpeak sub-command arguments.

                            # None to set.

                            # Set macs2_bdgcmp sub-command options.

                            mb2 = macs2_bdgcmp.sub_command

                            mb2.add_OptionLong(key='tfile',
                                               value=os.path.join(prefix,
                                                                  '{}_treat_pileup.bdg'.format(prefix)))

                            mb2.add_OptionLong(key='cfile',
                                               value=os.path.join(prefix,
                                                                  '{}_control_lambda.bdg'.format(prefix)))

                            # Sequencing depth for treatment and control. Aim for setting the --SPMR parameter for
                            # macs2_callpeak to get the track normalised.
                            # --tdepth:
                            # --cdepth:
                            # --pseudocount

                            mb2.add_OptionLong(key='output',
                                               value=os.path.join(prefix,
                                                                  '{}_bdgcmp.bdg'.format(prefix)))

                            # --method defaults to ppois i.e. Poisson Pvalue (-log10(pvalue), which yields data
                            # on a logarithmic scale.

                            # mb2.add_OptionLong(key='--method', value='FE')

                            # Set macs2_bdgcmp arguments.

                            # None to set.

                            # Set process_macs2 options.

                            # None to set.

                            # Set process_macs2 arguments.

                            process_macs2.arguments.append('{}__{}'.format(t_replicate_key, c_replicate_key))
                            process_macs2.arguments.append(genome_sizes)

    def _create_diffbind_jobs(self):
        """Create Bioconductor DiffBind jobs.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        diffbind_drms = DRMS.from_Analysis(name='diffbind',
                                           work_directory=self.genome_directory,
                                           analysis=self)
        self.drms_list.append(diffbind_drms)

        # Reorganise the comparisons by factor.

        self._factor_dict = dict()

        keys = self.comparisons.keys()
        # keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]
            assert isinstance(chipseq_comparison, ChIPSeqComparison)

            if chipseq_comparison.factor in self._factor_dict:
                factor_list = self._factor_dict[chipseq_comparison.factor]
            else:
                factor_list = list()
                self._factor_dict[chipseq_comparison.factor] = factor_list

            factor_list.append(chipseq_comparison)

        keys = self._factor_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

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

            # Create a new BSF ChIPSeq DiffBind Sheet per factor.

            file_path = os.path.join(factor_directory, 'chipseq_diffbind_{}_samples.csv'.format(key))

            dbs = ChIPSeqDiffBindSheet(file_path=file_path)

            factor_list = self._factor_dict[key]
            factor_list.sort(cmp=lambda x, y: cmp(x, y))

            for chipseq_comparison in factor_list:

                for t_sample in chipseq_comparison.t_samples:

                    t_replicate_dict = t_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                    # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                    # Python str key and Python list of Python list objects
                    # of Bio.BSF.Data.PairedReads objects.

                    t_replicate_keys = t_replicate_dict.keys()
                    t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                    for t_replicate_key in t_replicate_keys:

                        for c_sample in chipseq_comparison.c_samples:

                            c_replicate_dict = c_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

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
            dbs.write_to_file()

            # Create the DiffBind job.

            diffbind = Executable.from_Analysis(name='chipseq_diffbind_{}'.format(key),
                                                program='bsf_chipseq_diffbind.R',
                                                analysis=self)
            diffbind_drms.add_Executable(diffbind)
            diffbind.dependencies.extend(job_dependencies)

            # Add diffbind options.

            diffbind.add_OptionLong(key='factor', value=key)
            # diffbind.add_OptionLong(key='work_directory', value=factor_directory)
            diffbind.add_OptionLong(key='genome_directory', value=self.genome_directory)
            diffbind.add_OptionLong(key='sample_annotation', value=file_path)

    def _report_macs14(self):
        """Create a BSF ChIPSeq report in HTML format and a UCSC Genome Browser Track Hub.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        # ucsc_location = config_parser.get(section=config_section, option='ucsc_location')
        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        track_output = str()

        # Write a HTML document.

        output = str()

        output += Defaults.web.html_header(title='{} ChIP-Seq Analysis'.format(self.project_name))
        output += '<body>\n'
        output += '\n'

        output += '<h1>{} ChIP-Seq Analysis</h1>\n'.format(self.project_name)
        output += '\n'

        output += '<p>\n'
        output += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
        output += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></strong>,\n'
        output += 'before peaks are called with '
        output += '<a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 1.4</a>\n'
        output += 'on a treatment and control sample pair.\n'
        output += '</p>\n'
        output += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output += '<p>\n'
        output += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
        output += 'in their genomic context via the project-specific\n'
        output += 'UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.\n'. \
            format(Defaults.web.ucsc_track_url(options_dict=options_dict, host_name=default.ucsc_host_name),
                   self.project_name)
        output += '</p>\n'
        output += '\n'

        output += '<table>\n'
        output += '\n'

        output += '<tr>\n'
        output += '<th>Peaks</th>\n'
        output += '<th>Negative Peaks</th>\n'
        output += '<th>R Model</th>\n'
        output += '</tr>\n'
        output += '\n'

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                # Python str key and Python list of Python list objects
                # of Bio.BSF.Data.PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

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

                                    track_output += 'track ChIP_{}__{}_{}_{}\n'. \
                                        format(t_replicate_key, c_replicate_key, state, scaling)
                                    track_output += 'type bigWig\n'
                                    track_output += 'shortLabel ChIP_{}__{}_{}_{}\n'. \
                                        format(t_replicate_key, c_replicate_key, state, scaling)
                                    track_output += 'longLabel {} ChIP-Seq read counts for {} of {} versus {}\n'. \
                                        format(scaling.capitalize(), state, t_replicate_key, c_replicate_key)
                                    track_output += 'bigDataUrl '
                                    track_output += './{}__{}_MACS_wiggle/{}/{}__{}_{}_afterfiting_all.bw\n'. \
                                        format(t_replicate_key, c_replicate_key, state,
                                               t_replicate_key, c_replicate_key, state)
                                    if treatment and not absolute:
                                        track_output += 'visibility full\n'
                                    else:
                                        track_output += 'visibility hide\n'
                                        # 'html' is missing from the common settings.

                                    # Common optional settings.

                                    track_output += 'color {}\n'. \
                                        format(Defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                                    # bigWig - Signal graphing track settings.

                                    track_output += 'graphTypeDefault bar\n'
                                    track_output += 'maxHeightPixels 100:60:20\n'
                                    track_output += 'smoothingWindow off\n'

                                    if absolute:
                                        # Track with absolute scaling.
                                        track_output += 'autoScale on\n'
                                    else:
                                        # Track with relative scaling.
                                        track_output += 'autoScale off\n'
                                        track_output += 'viewLimits 0:40\n'

                                    track_output += '\n'

                            # Add a UCSC trackDB entry for each bigBed peaks file.

                            track_output += 'track Peaks_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'type bigBed\n'
                            track_output += 'shortLabel Peaks_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP-Seq peaks for {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'bigDataUrl ./{}__{}_peaks.bb\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'visibility pack\n'
                            # 'html' is missing from the common settings.

                            # Common optional settings.
                            track_output += 'color {}\n'. \
                                format(Defaults.web.get_chipseq_colour(factor=chipseq_comparison.factor))

                            track_output += '\n'

                            # Add web content.

                            output += '<tr>\n'

                            # if treatment and absolute:
                            # output += '<td><strong>{}</strong></td>\n'.format(t_replicate_key)
                            # if not treatment and absolute:
                            # output += '<td><strong>{}</strong></td>\n'.format(c_replicate_key)

                            output += '<td><a href="./{}__{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                format(t_replicate_key, c_replicate_key,
                                       t_replicate_key, c_replicate_key)

                            output += '<td>'
                            output += '<a href="./{}__{}_negative_peaks.xls">Negative peaks {} versus {}</a>'. \
                                format(t_replicate_key, c_replicate_key,
                                       t_replicate_key, c_replicate_key)
                            output += '</td>\n'

                            output += '<td><a href="./{}__{}_model.r">R model</a></td>\n'. \
                                format(t_replicate_key, c_replicate_key,
                                       t_replicate_key, c_replicate_key)

                            output += '</tr>\n'
                            output += '\n'

        output += '</table>\n'
        output += '\n'

        # Add UCSC trackDB entries for each Bowtie2 BAM file.

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                # Add a UCSC trackDB entry.

                # Common trackDb settings.

                track_output += 'track Alignment_{}\n'.format(replicate_key)
                track_output += 'type bam\n'
                track_output += 'shortLabel Alignment_{}\n'.format(replicate_key)
                track_output += 'longLabel Bowtie2 alignment of {}\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./{}.aligned.sorted.bam\n'. \
                    format(replicate_key)
                track_output += 'visibility hide\n'
                # TODO: The 'html' option is missing.

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                track_output += '\n'

                # TODO: Including FastQC results would require rearranging the above HTML code.
                # output += '{}__{}_fastqc/fastqc_report.html'. \
                # format(t_replicate_key, c_replicate_key)

        output += '</body>\n'
        output += Defaults.web.html_footer()

        file_path = os.path.join(self.genome_directory, 'chipseq_report.html')

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        # Create the UCSC Genome Browser Track Hub.

        self.ucsc_hub_write_hub(prefix='chipseq')
        self.ucsc_hub_write_genomes(prefix='chipseq')
        self.ucsc_hub_write_tracks(output=track_output, prefix='chipseq')

    def report(self):
        """Create a BSF ChIPSeq report in HTML format and a UCSC Genome Browser Track Hub.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        # ucsc_location = config_parser.get(section=config_section, option='ucsc_location')
        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        contrast_field_names = ["", "Group1", "Members1", "Group2", "Members2", "DB.edgeR"]

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        track_output = str()

        # Write a HTML document.

        output = str()

        output += Defaults.web.html_header(title='{} ChIP-Seq Analysis'.format(self.project_name))
        output += '<body>\n'
        output += '\n'

        output += '<h1>{} ChIP-Seq Analysis</h1>\n'.format(self.project_name)
        output += '\n'

        output += '<p>\n'
        output += 'Next-Generation Sequencing reads are aligned with the short read aligner\n'
        output += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a></strong>,\n'
        output += 'before peaks are called with <a href="http://liulab.dfci.harvard.edu/MACS/index.html">MACS 2</a>\n'
        output += 'on an enrichment and background sample pair.\n'
        output += '</p>\n'
        output += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output += '<p>\n'
        output += 'View Bowtie2 <strong>read alignment</strong> tracks for each sample\n'
        output += 'in their genomic context via the project-specific\n'
        output += 'UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.\n'. \
            format(Defaults.web.ucsc_track_url(options_dict=options_dict, host_name=default.ucsc_host_name),
                   self.project_name)
        output += '</p>\n'
        output += '\n'

        output += '<table>\n'
        output += '\n'

        output += '<tr>\n'
        output += '<th>Peaks</th>\n'
        output += '<th>R Model</th>\n'
        output += '</tr>\n'
        output += '\n'

        # Group via UCSC super tracks.

        track_output += 'track Alignment\n'
        track_output += 'shortLabel QC Alignment\n'
        track_output += 'longLabel ChIP Read Alignment\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignment\n'
        track_output += '\n'

        track_output += 'track Background\n'
        track_output += 'shortLabel QC Background\n'
        track_output += 'longLabel ChIP Background Signal\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group background\n'
        track_output += '\n'

        track_output += 'track Enrichment\n'
        track_output += 'shortLabel QC Enrichment\n'
        track_output += 'longLabel ChIP Enrichment Signal\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group enrichment\n'
        track_output += '\n'

        track_output += 'track Comparison\n'
        track_output += 'shortLabel ChIP Intensity\n'
        track_output += 'longLabel ChIP Intensity\n'
        track_output += 'visibility full\n'
        track_output += 'superTrack on\n'
        track_output += 'group comparison\n'
        track_output += '\n'

        track_output += 'track Peaks\n'
        track_output += 'shortLabel ChIP Peaks\n'
        track_output += 'longLabel ChIP Peaks\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group peaks\n'
        track_output += '\n'

        track_output += 'track Summits\n'
        track_output += 'shortLabel ChIP Summits\n'
        track_output += 'longLabel ChIP Summits\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group summits\n'
        track_output += '\n'

        # Group via UCSC composite tracks.

        composite_groups = dict()

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            chipseq_comparison = self.comparisons[key]

            factor = chipseq_comparison.factor.upper()

            for t_sample in chipseq_comparison.t_samples:

                t_replicate_dict = t_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                # Python str key and Python list of Python list objects
                # of Bio.BSF.Data.PairedReads objects.

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:

                    for c_sample in chipseq_comparison.c_samples:

                        c_replicate_dict = c_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

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
                                track_output += 'track {}\n'.format(composite_group)
                                track_output += 'type bigWig\n'
                                track_output += 'shortLabel {}\n'.format(composite_group)
                                track_output += 'longLabel ChIP background signal for factor {!r}\n'. \
                                    format(factor)
                                track_output += 'visibility dense\n'
                                track_output += 'compositeTrack on\n'
                                track_output += 'parent Background\n'
                                track_output += 'allButtonPair on\n'
                                track_output += 'centerLabelsDense on\n'
                                track_output += '\n'

                            # Common trackDb settings.

                            track_output += 'track {}__{}_background\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                            track_output += 'type bigWig\n'
                            track_output += 'shortLabel {}__{}_bkg\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP background signal {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'bigDataUrl ./{}/{}_control_lambda.bw\n'. \
                                format(prefix, prefix)
                            track_output += 'visibility dense\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            track_output += 'color {}\n'. \
                                format(Defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            track_output += 'alwaysZero off\n'
                            track_output += 'autoScale off\n'
                            track_output += 'graphTypeDefault bar\n'
                            track_output += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            track_output += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            track_output += 'viewLimits 0:15\n'
                            track_output += 'viewLimitsMax 0:40\n'
                            track_output += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            track_output += 'parent {} on\n'.format(composite_group)
                            track_output += 'centerLabelsDense off\n'
                            track_output += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_treat_pileup.bw file.
                            #

                            composite_group = '{}_enrichment'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                track_output += 'track {}\n'.format(composite_group)
                                track_output += 'type bigWig\n'
                                track_output += 'shortLabel {}_enrichment\n'.format(factor)
                                track_output += 'longLabel ChIP enrichment signal for factor {!r}\n'. \
                                    format(factor)
                                track_output += 'visibility dense\n'
                                track_output += 'compositeTrack on\n'
                                track_output += 'parent Enrichment\n'
                                track_output += 'allButtonPair on\n'
                                track_output += 'centerLabelsDense on\n'
                                track_output += '\n'

                            # Common trackDb settings.

                            track_output += 'track {}__{}_enrichment\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            track_output += 'type bigWig\n'
                            track_output += 'shortLabel {}__{}_enr\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP enrichment signal {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'bigDataUrl ./{}/{}_treat_pileup.bw\n'. \
                                format(prefix, prefix)
                            track_output += 'visibility dense\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            track_output += 'color {}\n'. \
                                format(Defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            track_output += 'alwaysZero off\n'
                            track_output += 'autoScale off\n'
                            track_output += 'graphTypeDefault bar\n'
                            track_output += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            track_output += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            track_output += 'viewLimits 0:15\n'
                            track_output += 'viewLimitsMax 0:40\n'
                            track_output += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            track_output += 'parent {} on\n'.format(composite_group)
                            track_output += 'centerLabelsDense off\n'
                            track_output += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_bdgcmp.bw file.
                            #

                            composite_group = '{}_intensity'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                track_output += 'track {}\n'.format(composite_group)
                                track_output += 'type bigWig\n'
                                track_output += 'shortLabel {}\n'.format(composite_group)
                                track_output += 'longLabel ChIP intensity for factor {!r}\n'. \
                                    format(factor)
                                track_output += 'visibility full\n'
                                track_output += 'compositeTrack on\n'
                                track_output += 'parent Comparison\n'
                                track_output += 'allButtonPair on\n'
                                track_output += 'centerLabelsDense on\n'
                                track_output += '\n'

                            # Common trackDb settings.

                            track_output += 'track {}__{}_intensity\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            track_output += 'type bigWig\n'
                            track_output += 'shortLabel {}__{}_int\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP intensity {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'bigDataUrl ./{}/{}_bdgcmp.bw\n'. \
                                format(prefix, prefix)
                            track_output += 'visibility full\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            track_output += 'color {}\n'. \
                                format(Defaults.web.get_chipseq_colour(factor=factor))

                            # bigWig - Signal graphing track settings.

                            track_output += 'alwaysZero off\n'
                            track_output += 'autoScale off\n'
                            track_output += 'graphTypeDefault bar\n'
                            track_output += 'maxHeightPixels 100:60:20\n'
                            # track_output += 'maxWindowToQuery 10000000\n'
                            track_output += 'smoothingWindow 5\n'
                            # track_output += 'transformFunc NONE\n'
                            track_output += 'viewLimits 0:15\n'
                            track_output += 'viewLimitsMax 0:40\n'
                            track_output += 'windowingFunction maximum\n'
                            # track_output += 'yLineMark <#>\n'
                            # track_output += 'yLineOnOff on \n'
                            # track_output += 'gridDefault on\n'

                            # Composite track settings.

                            track_output += 'parent {} on\n'.format(composite_group)
                            track_output += 'centerLabelsDense off\n'
                            track_output += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_peaks.bb file.
                            #

                            # Common trackDb settings.

                            track_output += 'track {}__{}_peaks\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'type bigBed\n'
                            track_output += 'shortLabel {}__{}_peaks\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP peaks {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'bigDataUrl ./{}/{}_peaks.bb\n'. \
                                format(prefix, prefix)
                            track_output += 'visibility pack\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            track_output += 'color {}\n'. \
                                format(Defaults.web.get_chipseq_colour(factor=factor))

                            # bigBed - Item or region track settings.

                            # Supertrack settings.

                            track_output += 'parent Peaks\n'
                            track_output += '\n'

                            #
                            # Add a UCSC trackDB entry for each NAME_summits.bb file.
                            #

                            # Common trackDb settings.

                            track_output += 'track {}__{}_summits\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'type bigBed\n'
                            track_output += 'shortLabel {}__{}_summits\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP summits {} versus {}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'bigDataUrl ./{}/{}_summits.bb\n'. \
                                format(prefix, prefix)
                            track_output += 'visibility pack\n'
                            # track_output += 'html {}\n'.format()

                            # Common optional settings.

                            track_output += 'color {}\n'. \
                                format(Defaults.web.get_chipseq_colour(factor=factor))

                            # Supertrack settings.

                            track_output += 'parent Summits\n'
                            track_output += '\n'

                            #
                            # Add web content.
                            #

                            output += '<tr>\n'

                            output += '<td><a href="./{}/{}_peaks.xls">Peaks {} versus {}</a></td>\n'. \
                                format(prefix, prefix, t_replicate_key, c_replicate_key)

                            output += '<td><a href="./{}/{}_model.r">R model</a></td>\n'. \
                                format(prefix, prefix)

                            output += '</tr>\n'
                            output += '\n'

        output += '</table>\n'
        output += '\n'

        # Add UCSC trackDB entries for each Bowtie2 BAM file.

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                #
                # Add a UCSC trackDB entry for each NAME.aligned.sorted.bam file.
                #

                # Common trackDb settings.

                track_output += 'track {}_alignment\n'.format(replicate_key)
                track_output += 'type bam\n'
                track_output += 'shortLabel {}_alignment\n'.format(replicate_key)
                track_output += 'longLabel {} ChIP read alignment\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./chipseq_bowtie2_{}/{}.aligned.sorted.bam\n'. \
                    format(replicate_key, replicate_key)
                track_output += 'visibility hide\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                # Supertrack settings.

                track_output += 'parent Alignment\n'
                track_output += '\n'

        # Differential binding analysis.

        output += '<h2>Differential Binding Analysis</h2>\n'
        output += '\n'

        output += '<table>\n'
        output += '\n'

        output += '<tr>\n'
        output += '<th>Factor and Contrast</th>\n'
        output += '<th>Correlation Peak Caller</th>\n'
        output += '<th>Correlation Peak Counts</th>\n'
        output += '<th>Correlation Analysis</th>\n'
        output += '<th>MA Plot</th>\n'
        output += '<th>Scatter Plot</th>\n'
        output += '<th>PCA Plot</th>\n'
        output += '<th>Box Plot</th>\n'
        output += '<th>DiffBind Report</th>\n'
        output += '</tr>\n'

        keys = self._factor_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            prefix = 'chipseq_diffbind_{}'.format(key)

            output += '<tr>\n'

            # ChIP Factor

            output += '<td><strong>{}</strong></td>\n'.format(key)

            # Correlation heat map of peak caller scores.

            output += '<td>'
            output += '<a href="./{}/{}_correlation_peak_caller_score.png">'. \
                format(prefix, prefix)
            output += '<img alt="DiffBind correlation analysis for factor {}" ' \
                      'src="./{}/{}_correlation_peak_caller_score.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            output += '</a>'
            output += '</td>\n'

            # Correlation heat map of counts.

            output += '<td>'
            output += '<a href="./{}/{}_correlation_read_counts.png">'. \
                format(prefix, prefix)
            output += '<img alt="DiffBind correlation analysis for factor {}" ' \
                      'src="./{}/{}_correlation_read_counts.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            output += '</a>'
            output += '</td>\n'

            # Correlation heat map of differential binding analysis.

            output += '<td>'
            output += '<a href="./{}/{}_correlation_analysis.png">'. \
                format(prefix, prefix)
            output += '<img alt="DiffBind correlation analysis for factor {}" ' \
                      'src="./{}/{}_correlation_analysis.png" height="80" width="80">'. \
                format(key, prefix, prefix)
            output += '</a>'
            output += '</td>\n'

            output += '<td></td>\n'  # MA Plot
            output += '<td></td>\n'  # Scatter Plot
            output += '<td></td>\n'  # PCA Plot
            output += '<td></td>\n'  # Box Plot
            output += '<td></td>\n'  # DiffBin Report

            output += '</tr>\n'

            # Read the file of contrasts ...

            file_path = os.path.join(self.genome_directory, prefix, '{}_contrasts.csv'.format(prefix))

            if not os.path.exists(path=file_path):
                warnings.warn(
                    'File {!r} does not exist.'.format(file_path),
                    UserWarning)
                continue

            sas = SampleAnnotationSheet(file_path=file_path, field_names=contrast_field_names)
            sas.csv_reader_open()

            for row_dict in sas._csv_reader_object:
                suffix = '{}__{}'.format(row_dict['Group1'], row_dict['Group2'])

                output += '<tr>\n'

                output += '<td>{}</td>\n'.format(suffix)
                output += '<td></td>\n'  # Correlation heat map of peak caller scores.
                output += '<td></td>\n'  # Correlation heat map of counts.
                output += '<td></td>\n'  # Correlation heat map of differential binding analysis.

                # MA Plot

                output += '<td>'
                output += '<a href="./{}/{}_ma_plot_{}.png">'.format(prefix, prefix, suffix)
                output += '<img alt="DiffBind MA plot for factor {}" ' \
                          'src="./{}/{}_ma_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output += '</a>'
                output += '</td>\n'

                # Scatter Plot

                output += '<td>'
                output += '<a href="./{}/{}_scatter_plot_{}.png">'.format(prefix, prefix, suffix)
                output += '<img alt="DiffBind scatter plot for factor {}" ' \
                          'src="./{}/{}_scatter_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output += '</a>'
                output += '</td>\n'

                # Principal Component Analysis Plot

                output += '<td>'
                output += '<a href="./{}/{}_pca_plot_{}.png">'.format(prefix, prefix, suffix)
                output += '<img alt="DiffBind PCA plot for factor {}" ' \
                          'src="./{}/{}_pca_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output += '</a>'
                output += '</td>\n'

                # Box Plot

                output += '<td>'
                output += '<a href="./{}/{}_box_plot_{}.png">'.format(prefix, prefix, suffix)
                output += '<img alt="DiffBind Box plot for factor {}" ' \
                          'src="./{}/{}_box_plot_{}.png" height="80" width="80">'. \
                    format(key, prefix, prefix, suffix)
                output += '</a>'
                output += '</td>\n'

                # DiffBind report

                output += '<td><a href="./{}/DBA_{}_report_{}.csv">DBA_{}_report_{}</a></td>\n'. \
                    format(prefix, key, suffix, key, suffix)

                output += '</tr>\n'

            sas.csv_reader_close()

        output += '\n'
        output += '</table>\n'
        output += '\n'

        output += '</body>\n'
        output += Defaults.web.html_footer()

        file_path = os.path.join(self.genome_directory, 'chipseq_report.html')

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        # Create the UCSC Genome Browser Track Hub.

        self.ucsc_hub_write_hub(prefix='chipseq')
        self.ucsc_hub_write_genomes(prefix='chipseq')
        self.ucsc_hub_write_tracks(output=track_output, prefix='chipseq')


class RunFastQC(Analysis):
    """BSF FastQC-specific Quality Assessment Analysis sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new C{RunFastQC} object from a UNIX-style configuration file via the C{Configuration} class.

        @param config_file: UNIX-style configuration file
        @type config_file: str | unicode
        @return: C{RunFastQC}
        @rtype: RunFastQC
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new C{RunFastQC} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{RunFastQC}
        @rtype: RunFastQC
        """

        assert isinstance(configuration, Configuration)

        run_fastqc = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.RunFastQC" section specifies defaults for this BSF Analysis.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        run_fastqc.set_Configuration(run_fastqc.configuration, section=section)

        return run_fastqc

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):
        """Initialise a C{RunFastQC} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{Analysis}-wide project directory,
            normally under the C{Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{Analysis}-wide genome directory,
            normally under the C{Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param drms_list: Python C{list} of BSF C{DRMS} objects
        @type drms_list: list
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Python C{tuple} objects of C{Sample} objects
        @type comparisons: dict
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        """

        super(RunFastQC, self).__init__(configuration=configuration,
                                        project_name=project_name, genome_version=genome_version,
                                        input_directory=input_directory, output_directory=output_directory,
                                        project_directory=project_directory, genome_directory=genome_directory,
                                        e_mail=e_mail, debug=debug, drms_list=drms_list,
                                        collection=collection, comparisons=comparisons, samples=samples)

        # Nothing else to do for this sub-class ...

    def set_Configuration(self, configuration, section):
        """Set instance variables of a C{RunFastQC} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        super(RunFastQC, self).set_Configuration(configuration=configuration, section=section)

        # Read a comparison file.

        # if configuration.config_parser.has_option(section=section, option='cmp_file'):
        # self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')
        # Use the sample annotation sheet instead of a separate comparison file.
        if configuration.config_parser.has_option(section=section, option='sas_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='sas_file')

    def _read_comparisons(self, cmp_file):
        """Read a C{SampleAnnotationSheet} CSV file from disk.

        Column headers for CASAVA folders:
            Treatment/Control ProcessedRunFolder:
                CASAVA processed run folder name or
                C{Analysis.input_directory} by default
            Treatment/Control Project:
                CASAVA Project name or
                C{Analysis.project_name} by default
            Treatment/Control Sample:
                CASAVA Sample name, no default
        Column headers for independent samples:
            Treatment/Control Sample:
            Treatment/Control File:
        @param cmp_file: Comparisons file path
        @type cmp_file: str | unicode
        @return: Nothing
        @rtype: None
        """

        sas = SampleAnnotationSheet(file_path=cmp_file)
        sas.csv_reader_open()

        for row_dict in sas._csv_reader_object:
            self.add_Sample(sample=self.collection.get_Sample_from_row_dict(row_dict=row_dict))

        sas.csv_reader_close()

    def run(self):
        """Run this C{RunFastQC} C{Analysis}.
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

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_FastQC_jobs()

    def _create_FastQC_jobs(self):
        """Initialise a C{RunFastQC} object.
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
                if component in prf.projects:
                    project = prf.projects[component]
                    for sample in project.get_all_Samples():
                        self.add_Sample(sample=sample)
                else:
                    warnings.warn(
                        "Could not find projects component {!r} in RunFolder {!r}.\n"
                        "Check projects configuration value {!r}.".
                        format(component, prf.name, projects),
                        UserWarning)

        elif config_parser.has_option(section=config_section, option='samples'):

            # A configuration option 'samples' has been set ...

            sample = config_parser.get(section=config_section, option='samples')

            if sample == 'auto':

                # Automatically discover sample information from a CASAVA
                # Processed Run Folder and a CASAVA Project name.

                prf = ProcessedRunFolder.from_file_path(file_path=self.input_directory, file_type='Automatic')
                project = prf.projects[self.project_name]
                for sample in project.get_all_Samples():
                    self.add_Sample(sample=sample)

        else:

            # Since the BSF Collection object has complete BSF ProcessedRunFolder objects registered,
            # the sample annotation sheet needs re-reading.

            # TODO: This is now done in set_Configuration.
            # sas = SampleAnnotationSheet(file_path=sas_file)
            # sas.open_csv()
            #
            # for row_dict in sas._csv_reader:
            # self.add_sample(sample=self.collection.get_Sample_from_row_dict(row_dict=row_dict))
            #
            # sas.close_csv()

            pass

        fastqc_drms = DRMS.from_Analysis(name='fastqc',
                                         work_directory=self.project_directory,
                                         analysis=self)

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                fastqc = FastQC(name='fastqc_{}'.format(replicate_key),
                                analysis=self)

                fastqc_drms.add_Executable(fastqc)

                # Set FastQC options.

                fastqc.add_OptionLong(key='outdir', value=self.project_directory)

                if not 'casava' in fastqc.options and sample.file_type == 'CASAVA':
                    fastqc.add_SwitchLong(key='casava')

                fastqc.add_OptionLong(key='threads', value=str(fastqc_drms.threads))

                # Set FastQC arguments.

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:

                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

                fastqc.arguments.append(string.join(words=reads1 + reads2, sep=' '))

        self.drms_list.append(fastqc_drms)

    def report(self):
        """Create C{RunFastQC} report in HTML format.
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        # Get further information.

        # Always check each BSF PairedReads object separately.
        replicate_grouping = False

        # Let the user specify the report sub-directory or default to the projects directory.
        if config_parser.has_option(section=config_section, option='sub_directory'):
            sub_directory = config_parser.get(section=config_section, option='sub_directory')
        else:
            sub_directory = str()

        if not sub_directory:
            sub_directory = Default.get_global_default().url_relative_projects

        # Create a symbolic link containing the project name and a UUID.
        link_path = self.create_public_project_link(sub_directory=sub_directory)
        link_name = os.path.basename(link_path.rstrip('/'))

        # Write a HTML document.

        output = str()

        output += Defaults.web.html_header(title='{} FastQC Analysis'.format(self.project_name))
        output += '<body>\n'
        output += '\n'

        output += '<h1>{} FastQC Analysis</h1>\n'.format(self.project_name)
        output += '\n'

        output += '<p>\n'
        output += '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a>\n'
        output += 'is a quality control tool for high throughput sequence data.\n'
        output += '</p>\n'
        output += '\n'

        output += '<table>\n'
        output += '\n'

        output += '<tr>\n'
        output += '<th>Sample</th>\n'
        output += '<th>FastQC Report</th>\n'
        output += '<th>Summary</th>\n'
        output += '<th>ZIP archive</th>\n'
        output += '</tr>\n'
        output += '\n'

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.
            # Since FastQC is run on each replicate this needs the full name.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping, full=True)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                output += '<tr>\n'
                output += '<td>{}</td>\n'.format(replicate_key)
                output += '<td><a href="./{}_fastqc/fastqc_report.html"><strong>Report</strong></a></td>\n'. \
                    format(replicate_key)
                output += '<td><a href="./{}_fastqc/summary.txt">Summary</a></td>\n'.format(replicate_key)
                output += '<td><a href="./{}_fastqc.zip">ZIP archive</a></td>\n'.format(replicate_key)
                output += '</tr>\n'
                output += '\n'

        output += '</table>\n'
        output += '\n'

        output += '</body>\n'
        output += Defaults.web.html_footer()

        file_path = os.path.join(self.project_directory, 'fastqc_report.html')

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()


class RunBamToFastq(Analysis):
    """BSF BAM or SAM to FASTQ converter sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new BSF RunBamToFastq object from a UNIX-style configuration file via the BSF Configuration class.

        @param config_file: UNIX-style configuration file
        @type config_file: str | unicode
        @return: BSF RunBamToFastq
        @rtype: RunBamToFastq
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new BSF RunBamToFastq object from a BSF Configuration object.

        @param configuration: BSF Configuration
        @type configuration: Configuration
        @return: BSF RunBamToFastq
        @rtype: RunBamToFastq
        """

        assert isinstance(configuration, Configuration)

        run_bam_to_fastq = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.RunBamToFastq" section specifies defaults for this BSF Analysis.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        run_bam_to_fastq.set_Configuration(run_bam_to_fastq.configuration, section=section)

        return run_bam_to_fastq

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):
        """Initialise a Bio.BSF.Analysis.RunBamToFastq object.

        @param configuration: BSF Configuration
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: BSF Analysis-wide input directory
        @type input_directory: str
        @param output_directory: BSF Analysis-wide output directory
        @type output_directory: str
        @param project_directory: BSF Analysis-wide project directory,
            normally under the BSF Analysis-wide output directory
        @type project_directory: str
        @param genome_directory: BSF Analysis-wide genome directory,
            normally under the BSF Analysis-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param drms_list: Python list of BSF DRMS objects
        @type drms_list: list
        @param collection: BSF Collection
        @type collection: Collection
        @param comparisons: Python dict of Python list objects of BSF Sample objects
        @type comparisons: dict
        @param samples: Python list of BSF Sample objects
        @type samples: list
        @return: Nothing
        @rtype: None
        """

        super(RunBamToFastq, self).__init__(configuration=configuration,
                                            project_name=project_name, genome_version=genome_version,
                                            input_directory=input_directory, output_directory=output_directory,
                                            project_directory=project_directory, genome_directory=genome_directory,
                                            e_mail=e_mail, debug=debug, drms_list=drms_list,
                                            collection=collection, comparisons=comparisons, samples=samples)

        # Nothing else to do for this sub-class ...

    def set_Configuration(self, configuration, section):
        """Set instance variables of a BSF RunBamToFastq object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: BSF Configuration
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        super(RunBamToFastq, self).set_Configuration(configuration=configuration, section=section)

        # Nothing else to do for this sub-class ...

    def run(self):
        """Run this BSF RunBamToFastq analysis.

        @return: Nothing
        @rtype: None
        """

        super(RunBamToFastq, self).run()

        self._convert_BamToFastq()

    def _convert_BamToFastq(self):
        """Private method to convert all Reads objects in BAM or SAM format into FASTQ format.

        @return: Nothing
        @rtype: None
        """

        # config_parser = self.configuration.config_parser
        # config_section = self.configuration.section_from_instance(self)

        default = Default.get_global_default()

        # Replicates have to be un-grouped, always!
        # replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')
        replicate_grouping = False

        sam_to_fastq_drms = DRMS.from_Analysis(name='sam_to_fastq',
                                               work_directory=self.genome_directory,
                                               analysis=self)

        self.drms_list.append(sam_to_fastq_drms)

        for sample in self.collection.get_all_Samples():

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                for paired_reads in replicate_dict[replicate_key]:

                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    # In a BSF Paired Reads object, the SAM or BAM file could potentially
                    # occur as reads1 or reads2 instance variable.

                    if paired_reads.reads1:
                        file_name = str(paired_reads.reads1.file_path)
                        file_name = file_name.rstrip('/ ')
                        file_name = os.path.basename(file_name)

                        # TODO: The matching part to remove the .bam could be achieved with Bash parameter expansion.
                        match = re.search(pattern=r'(.*)\.bam$', string=file_name)
                        if match:
                            sam_to_fastq = Executable.from_Analysis(
                                name='picard_sam_to_fastq_{}_1'.format(replicate_key),
                                program='bsf_bam2fastq.sh',
                                analysis=self)

                            sam_to_fastq_drms.add_Executable(sam_to_fastq)

                            sam_to_fastq.arguments.append(paired_reads.reads1.file_path)
                            sam_to_fastq.arguments.append(os.path.join(default.classpath_picard,
                                                                       'SamToFastq.jar'))
                            sam_to_fastq.arguments.append(os.path.join(self.genome_directory,
                                                                       match.group(1)))

                    if paired_reads.reads2:
                        file_name = str(paired_reads.reads2.file_path)
                        file_name = file_name.rstrip('/ ')
                        file_name = os.path.basename(file_name)

                        match = re.search(pattern=r'(.*)\.bam$', string=file_name)
                        if match:
                            sam_to_fastq = Executable.from_Analysis(
                                name='picard_sam_to_fastq_{}_2'.format(replicate_key),
                                program='bsf_bam2fastq.sh',
                                analysis=self)

                            sam_to_fastq_drms.add_Executable(sam_to_fastq)

                            sam_to_fastq.arguments.append(paired_reads.reads2.file_path)
                            sam_to_fastq.arguments.append(os.path.join(default.classpath_picard,
                                                                       'SamToFastq.jar'))
                            sam_to_fastq.arguments.append(os.path.join(self.genome_directory,
                                                                       match.group(1)))
