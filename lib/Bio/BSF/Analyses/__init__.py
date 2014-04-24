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
from pickle import Pickler, HIGHEST_PROTOCOL
import re
import string
import warnings

from Bio.BSF import Analysis, Configuration, Default, Defaults, DRMS, Executable
from Bio.BSF.Data import Collection, ProcessedRunFolder, Sample, AnnotationSheet, SampleAnnotationSheet
from Bio.BSF.Executables import Bowtie2, BWA, Macs14, Macs2Bdgcmp, Macs2Callpeak, Cuffdiff, Cufflinks, Cuffmerge, \
    TopHat, FastQC


class ChIPSeqComparison(object):
    def __init__(self, c_name, t_name, c_samples, t_samples,
                 factor, tissue=None, condition=None, treatment=None, replicate=None):

        """Initialise a BSF ChIPSeqComparison object.
        :param self: BSF ChIPSeqComparison object
        :type self: ChIPSeqComparison
        :param c_name: Control name
        :type c_name: str
        :param t_name: Treatment name
        :type t_name: str
        :param c_samples: Python list of control BSF Sample objects
        :type c_samples: list
        :param t_samples: Python list of treatment BSF Sample objects
        :type t_samples:list
        :param factor: ChIP factor
        :type factor: str
        :param tissue: Tissue
        :type tissue: str
        :param condition: Condition
        :type condition: str
        :param treatment: Treatment
        :type treatment: str
        :param replicate: replicate number
        :type replicate: int
        :return Nothing
        :rtype: None
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
    field_names = ['SampleID', 'Tissue', 'Factor', 'Condition', 'Treatment', 'Replicate',
                   'bamReads', 'bamControl', 'ControlID', 'Peaks', 'PeakCaller', 'PeakFormat']

    test_methods = dict(SampleID=AnnotationSheet.check_alphanumeric,
                        Tissue=AnnotationSheet.check_alphanumeric,
                        Factor=AnnotationSheet.check_alphanumeric,
                        Condition=AnnotationSheet.check_alphanumeric,
                        Treatment=AnnotationSheet.check_alphanumeric,
                        Replicate=AnnotationSheet.check_numeric,
                        ControlID=AnnotationSheet.check_alphanumeric,
                        PeakCaller=AnnotationSheet.check_alphanumeric,
                        PeakFormat=AnnotationSheet.check_alphanumeric)

    def __init__(self, file_path=None, file_type=None, name=None, row_dicts=None):
        """Initialise a BSF ChIPSeq DiffBind Sheet object.

        :param self: BSF ChIPSeq DiffBind Sheet
        :type self: ChIPSeqDiffBindSheet
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. ...)
        :type file_type: str
        :param name: Name
        :type name: str
        :param row_dicts: Python list of Python dict objects
        :type row_dicts: list
        :return: Nothing
        :rtype: None
        """

        super(ChIPSeqDiffBindSheet, self).__init__(file_path=file_path,
                                                   file_type=file_type,
                                                   name=name,
                                                   field_names=ChIPSeqDiffBindSheet.field_names,
                                                   row_dicts=row_dicts)

    def sort(self):
        self.row_dicts.sort(cmp=lambda x, y:
                            cmp(x['Tissue'], y['Tissue']) or
                            cmp(x['Factor'], y['Factor']) or
                            cmp(x['Condition'], y['Condition']) or
                            cmp(x['Treatment'], y['Treatment']) or
                            cmp(int(x['Replicate']), int(y['Replicate'])))

    def write_to_file(self):
        """Write a BSF ChIPSeq DiffBind Sheet to a file.

        :param self: BSF ChIPSeq DiffBind Sheet
        :type self: ChIPSeqDiffBindSheet
        :return: Nothing
        :rtype: None
        """

        # Override the method from the super-class to automatically sort before writing to a file.

        self.sort()
        super(ChIPSeqDiffBindSheet, self).write_to_file()


class ChIPSeq(Analysis):
    """BSF ChIP-Seq-specific BSF Analysis sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new BSF ChIPSeq object from a UNIX-style configuration file via the BSF Configuration class.

        :param cls: Class
        :type cls: ChIPSeq
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF ChIPSeq
        :rtype: ChIPSeq
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new BSF ChIPSeq object from a BSF Configuration object.

        :param cls: Class
        :type cls: ChIPSeq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF ChIPSeq
        :rtype: ChIPSeq
        """

        assert isinstance(configuration, Configuration)

        chipseq = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.ChIPSeq" section specifies defaults for this BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        chipseq.set_Configuration(chipseq.configuration, section=section)

        return chipseq

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 cmp_file=None):

        """Initialise a Bio.BSF.Analysis.ChIPSeq object.

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: BSF Analysis-wide input directory
        :type input_directory: str
        :param output_directory: BSF Analysis-wide output directory
        :type output_directory: str
        :param project_directory: BSF Analysis-wide project directory,
         normally under the BSF Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: BSF Analysis-wide genome directory,
         normally under the BSF Analysis-wide project directory
        :type genome_directory: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of BSF DRMS objects
        :type drms_list: list
        :param collection: BSF Collection
        :type collection: Collection
        :param comparisons: Python dict of Python list objects of BSF Sample objects
        :type comparisons: dict
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :param cmp_file: Comparison file
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
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

        """Set instance variables of a BSF ChIPSeq object via a section of a
        BSF Configuration object.

        Instance variables without a
        configuration option remain unchanged.
        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(ChIPSeq, self).set_Configuration(configuration=configuration, section=section)

        # Read a comparison file.

        if configuration.config_parser.has_option(section=section, option='cmp_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')

    def _read_comparisons(self, cmp_file):

        """Read a BSF SampleAnnotationSheet CSV file from disk.

        Column headers for CASAVA folders:
          Treatment/Control ProcessedRunFolder:
            CASAVA processed run folder name or
            Bio.BSF.Analysis.input_directory by default.
          Treatment/Control Project:
            CASAVA Project name or
            Bio.BSF.Analysis.project_name by default.
          Treatment/Control Sample:
            CASAVA Sample name, no default.
        Column headers for independent samples:
          Treatment/Control Sample:
          Treatment/Control File:
          Treatment/Control Group:
        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :param cmp_file: Comparison file path
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        if self.debug > 1:
            print '{!r} method _read_comparisons:'.format(self)

        sas = SampleAnnotationSheet(file_path=cmp_file)

        # Unfortunately, two passes through the comparison sheet are required.
        # In the first one merge all BSF Sample objects that share the name.
        # Merging BSF Sample objects is currently the only way to pool BSF PairedReads objects,
        # which is required for ChIP-Seq experiments.

        sample_dict = dict()

        # First pass, merge BSF Sample objects, if they have the same name.
        sas.csv_reader_open()
        for row_dict in sas._csv_reader:
            for prefix in ('Control', 'Treatment'):
                name, samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                for o_sample in samples:
                    if o_sample.name not in sample_dict:
                        sample_dict[o_sample.name] = o_sample
                    else:
                        n_sample = Sample.from_Samples(sample1=sample_dict[o_sample.name], sample2=o_sample)
                        sample_dict[n_sample.name] = n_sample
        sas.csv_reader_close()

        # Second pass, add all BSF Sample objects mentioned in a comparison.
        level1_dict = dict()

        sas.csv_reader_open()
        for row_dict in sas._csv_reader:

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

        sas.csv_reader_close()

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

        """Run this BSF ChIPSeq analysis.

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
        """

        super(ChIPSeq, self).run()

        # ChIPSeq requires a genome version.

        if not self.genome_version:
            message = 'A ChIPSeq analysis requires a genome_version configuration option.'
            raise Exception(message)

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of BSF Sample objects by the BSF Sample name.
        # This cannot be done in the super-class, because BSF Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_Bowtie2_jobs()
        # self._create_Macs14_jobs()
        self._create_Macs2_jobs()
        self._create_diffbind_jobs()

    def _create_Bowtie2_jobs(self):

        """Create Bowtie2 alignment jobs.

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
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
                    bowtie2.add_OptionShort(key='U', value=string.join(reads1, ','))
                elif len(reads1) and len(reads2):
                    bowtie2.add_OptionShort(key='1', value=string.join(reads1, ','))
                if len(reads2):
                    bowtie2.add_OptionShort(key='2', value=string.join(reads2, ','))

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

    def _create_Macs14_jobs(self):

        """Create MACS14 peak caller jobs.

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
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
                                message = 'Unable to set MACS14 parameters for unknown factor {!r}. ' \
                                          'Please use default factor {!r} or adjust Python code if necessary.'. \
                                    format(factor, Defaults.web.chipseq_default_factor)
                                # raise Exception(message)
                                warnings.warn(message, UserWarning)

                            # Set macs14 arguments.

                            # Set process_macs14 options.
                            # Set process_macs14 arguments.

                            process_macs14.arguments.append(macs14.options['name'].value)
                            process_macs14.arguments.append(genome_sizes)

    def _create_Macs2_jobs(self):

        """Create MACS2 peak caller jobs.

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
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
                                message = 'Unable to set MACS2 parameters for unknown factor {!r}. ' \
                                          'Please use default factor {!r} or adjust Python code if necessary.'. \
                                    format(factor, Defaults.web.chipseq_default_factor)
                                # raise Exception(message)
                                warnings.warn(message, UserWarning)

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

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
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
                                ### sas.csv_writer_next(row_dict=row_dict)
                                dbs.row_dicts.append(row_dict)

                                job_dependency = 'chipseq_process_macs2_{}__{}'.format(t_replicate_key, c_replicate_key)
                                job_dependencies.append(job_dependency)

            # TODO: Remove once the code works.
            ### sas.csv_writer_close()
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

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
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
        output += 'View Bowtie2 <strong>read alignments</strong> tracks for each sample\n'
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
                            #     output += '<td><strong>{}</strong></td>\n'.format(t_replicate_key)
                            # if not treatment and absolute:
                            #     output += '<td><strong>{}</strong></td>\n'.format(c_replicate_key)

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

                track_output += 'track Alignments_{}\n'.format(replicate_key)
                track_output += 'type bam\n'
                track_output += 'shortLabel Alignments_{}\n'.format(replicate_key)
                track_output += 'longLabel Bowtie2 alignments of {}\n'. \
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
                #     format(t_replicate_key, c_replicate_key)

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

        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :return: Nothing
        :rtype: None
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
        output += 'on a treatment and control sample pair.\n'
        output += '</p>\n'
        output += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/chipseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output += '<p>\n'
        output += 'View Bowtie2 <strong>read alignments</strong> tracks for each sample\n'
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

        track_output += 'track Alignments\n'
        track_output += 'shortLabel Alignments\n'
        track_output += 'longLabel Bowtie2 alignment tracks\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Control\n'
        track_output += 'shortLabel Control\n'
        track_output += 'longLabel Local lambda control\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group control\n'
        track_output += '\n'

        track_output += 'track Treatment\n'
        track_output += 'shortLabel Treatment\n'
        track_output += 'longLabel Treatment\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group treatment\n'
        track_output += '\n'

        track_output += 'track Comparison\n'
        track_output += 'shortLabel Comparison\n'
        track_output += 'longLabel Comparison of Treatment with Control\n'
        track_output += 'visibility full\n'
        track_output += 'superTrack on\n'
        track_output += 'group comparison\n'
        track_output += '\n'

        track_output += 'track Peaks\n'
        track_output += 'shortLabel Peaks\n'
        track_output += 'longLabel Peaks\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group peaks\n'
        track_output += '\n'

        track_output += 'track Summits\n'
        track_output += 'shortLabel Summits\n'
        track_output += 'longLabel Summits\n'
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

                            composite_group = 'ChIP_{}_control'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                track_output += 'track {}\n'.format(composite_group)
                                track_output += 'type bigWig\n'
                                track_output += 'shortLabel ChIP_{}_control\n'.format(factor)
                                track_output += 'longLabel ChIP-Seq control track for factor {!r}\n'. \
                                    format(factor)
                                track_output += 'visibility dense\n'
                                track_output += 'compositeTrack on\n'
                                track_output += 'parent Control\n'
                                track_output += 'allButtonPair on\n'
                                track_output += 'centerLabelsDense on\n'
                                track_output += '\n'

                            # Common trackDb settings.

                            track_output += 'track ChIP_{}__{}_control\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                            track_output += 'type bigWig\n'
                            track_output += 'shortLabel C_{}__{}_ctl\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP-Seq control lambda of {} versus {}\n'. \
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

                            composite_group = 'ChIP_{}_treatment'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                track_output += 'track {}\n'.format(composite_group)
                                track_output += 'type bigWig\n'
                                track_output += 'shortLabel ChIP_{}_treatment\n'.format(factor)
                                track_output += 'longLabel ChIP-Seq treatment track for factor {!r}\n'. \
                                    format(factor)
                                track_output += 'visibility dense\n'
                                track_output += 'compositeTrack on\n'
                                track_output += 'parent Treatment\n'
                                track_output += 'allButtonPair on\n'
                                track_output += 'centerLabelsDense on\n'
                                track_output += '\n'

                            # Common trackDb settings.

                            track_output += 'track ChIP_{}__{}_treatment\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            track_output += 'type bigWig\n'
                            track_output += 'shortLabel C_{}__{}_trt\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP-Seq pileup for treatment of {} versus {}\n'. \
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

                            composite_group = 'ChIP_{}_comparison'.format(factor)
                            if composite_group not in composite_groups:
                                composite_groups[composite_group] = 1
                                track_output += 'track {}\n'.format(composite_group)
                                track_output += 'type bigWig\n'
                                track_output += 'shortLabel ChIP_{}_comparison\n'.format(factor)
                                track_output += 'longLabel ChIP-Seq comparison track for factor {!r}\n'. \
                                    format(factor)
                                track_output += 'visibility full\n'
                                track_output += 'compositeTrack on\n'
                                track_output += 'parent Comparison\n'
                                track_output += 'allButtonPair on\n'
                                track_output += 'centerLabelsDense on\n'
                                track_output += '\n'

                            # Common trackDb settings.

                            track_output += 'track ChIP_{}__{}_comparison\n'. \
                                format(t_replicate_key, c_replicate_key)
                            # TODO: The bigWig type must declare the expected signal range.
                            track_output += 'type bigWig\n'
                            track_output += 'shortLabel C_{}__{}_cmp\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP-Seq comparison for treatment of {} versus {}\n'. \
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

                            track_output += 'track Peaks_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'type bigBed\n'
                            track_output += 'shortLabel Peaks_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP-Seq peaks for {} versus {}\n'. \
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

                            track_output += 'track Summits_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'type bigBed\n'
                            track_output += 'shortLabel Summits_{}__{}\n'. \
                                format(t_replicate_key, c_replicate_key)
                            track_output += 'longLabel ChIP-Seq summits for {} versus {}\n'. \
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

                track_output += 'track Alignments_{}\n'.format(replicate_key)
                track_output += 'type bam\n'
                track_output += 'shortLabel Alignments_{}\n'.format(replicate_key)
                track_output += 'longLabel Bowtie2 alignments of {}\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./chipseq_bowtie2_{}/{}.aligned.sorted.bam\n'. \
                    format(replicate_key, replicate_key)
                track_output += 'visibility hide\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                # track_output += 'color {}\n'.format(Defaults.web.chipseq_colours[comparison[2].upper()])

                # bam - Compressed Sequence Alignment track settings.

                # Supertrack settings.

                track_output += 'parent Alignments\n'
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
                message = 'File {!r} does not exist.'.format(file_path)
                warnings.warn(message, UserWarning)
                continue

            sas = SampleAnnotationSheet(file_path=file_path, field_names=contrast_field_names)
            sas.csv_reader_open()

            for row_dict in sas._csv_reader:
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


class RNASeq(Analysis):
    """BSF RNA-Seq-specific BSF Analysis sub-class.

    Attributes:
    :ivar cmp_file: Comparison file
    :type cmp_file: str, unicode
    """

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new BSF RNASeq object from a UNIX-style configuration file via the BSF Configuration class.

        :param cls: Class
        :type cls: RNASeq
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF RNASeq
        :rtype: RNASeq
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new BSF RNASeq object from a BSF Configuration object.

        :param cls: Class
        :type cls: RNASeq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF RNASeq
        :rtype: RNASeq
        """

        assert isinstance(configuration, Configuration)

        rnaseq = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.RNASeq" section specifies defaults for this BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        rnaseq.set_Configuration(rnaseq.configuration, section=section)

        return rnaseq

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 cmp_file=None):

        """Initialise a Bio.BSF.Analysis.RNASeq object.

        :param self: BSF RNASeq
        :type self: RNASeq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: BSF Analysis-wide input directory
        :type input_directory: str
        :param output_directory: BSF Analysis-wide output directory
        :type output_directory: str
        :param project_directory: BSF Analysis-wide project directory,
        normally under the BSF Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: BSF Analysis-wide genome directory,
        normally under the BSF Analysis-wide project directory
        :type genome_directory: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of BSF DRMS objects
        :type drms_list: list
        :param collection: BSF Collection
        :type collection: Collection
        :param comparisons: Python dict of Python tuple objects of BSF Sample objects
        :type comparisons: dict
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :param cmp_file: Comparison file
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        super(RNASeq, self).__init__(configuration=configuration,
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

        """Set instance variables of a BSF RNASeq object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF RNASeq
        :type self: RNASeq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(RNASeq, self).set_Configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        if configuration.config_parser.has_option(section=section, option='cmp_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')

    def _read_comparisons(self, cmp_file):

        """Read a BSF SampleAnnotationSheet CSV file from disk.

        Column headers for CASAVA folders:
          Treatment/Control ProcessedRunFolder:
            CASAVA processed run folder name or
            Bio.BSF.Analysis input_directory by default.
          Treatment/Control Project:
            CASAVA Project name or
            Bio.BSF.Analysis project_name by default.
          Treatment/Control Sample:
            CASAVA Sample name, no default.
        Column headers for independent samples:
          Treatment/Control Sample:
          Treatment/Control File:
        :param self: BSF ChIPSeq
        :type self: ChIPSeq
        :param cmp_file: Comparisons file path
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        if self.debug > 1:
            print '{!r} method _read_comparisons:'.format(self)

        sas = SampleAnnotationSheet(file_path=cmp_file)
        sas.csv_reader_open()

        for row_dict in sas._csv_reader:

            # In addition to defining samples allow the definition of groups.
            # If the row dictionary has a 'Group' key then the BSF Sample in the same row gets added to the group.
            # So 'ProcessedRunFolder',Project','Sample','Group' defines the groups, while ...
            # 'Control Group','Treatment Group' defines a comparison, as does ...
            # 'Control Group','Treatment ProcessedRunFolder','Treatment Project','Treatment Sample'
            # Define groups of BSF sample objects to compare these groups rather than
            # individual BSF Sample objects.

            # This can be one or more BSF Sample objects depending on 'Group' or 'Sample' column entries.
            c_name, c_samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix='Control')
            t_name, t_samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix='Treatment')

            # TODO: This should be extended to cope with a time-series...

            # TODO: There are two ways to interpret a group. RNA-Seq needs the pool variant of merging.

            if 0:

                # Run each BSF Sample in the first (control) group against
                # each BSF Sample in the second (treatment) group.

                for c_sample in c_samples:

                    if self.debug > 1:
                        print '  Control Sample name: {}'.format(c_sample.name)
                        # print c_sample.trace(1)

                    self.add_Sample(sample=c_sample)

                    for t_sample in t_samples:

                        if self.debug > 1:
                            print '    Treatment Sample name: {}'.format(t_sample.name)
                            # print t_sample.trace(2)

                        self.add_Sample(sample=t_sample)

                        key = '{}__{}'.format(c_sample.name, t_sample.name)
                        self.comparisons[key] = ([c_sample], [t_sample])

            else:

                # Run a pool of all BSF Sample objects in the first group against
                # a pool of all BSF Sample objects in the second group.

                # For a successful comparison, both Python list objects of BSF Sample objects have to be defined.

                if not (len(t_samples) and len(c_samples)):
                    if self.debug > 1:
                        print 'Comparison line with {} Treatment samples and {} Control samples'. \
                            format(len(t_samples), len(c_samples))
                    continue

                # Add all control BSF Sample objects to the Sample dictionary.

                for c_sample in c_samples:

                    if self.debug > 1:
                        print '  Control Sample name: {!r} file_path:{!r}'.format(c_sample.name, c_sample.file_path)
                        # print c_sample.trace(1)

                    self.add_Sample(sample=c_sample)

                # Add all treatment BSF Sample objects to the Sample dictionary.

                for t_sample in t_samples:

                    if self.debug > 1:
                        print '  Treatment Sample name: {!r} file_path:{!r}'.format(t_sample.name, t_sample.file_path)
                        # print t_sample.trace(1)

                    self.add_Sample(sample=t_sample)

                # TODO: Use a RNASeqComparison object here?
                key = '{}__{}'.format(c_name, t_name)
                self.comparisons[key] = (c_samples, t_samples)

        sas.csv_reader_close()

    def run(self):

        """Run this BSF RNASeq analysis.

        :param self: BSF RNASeq analysis
        :type self: RNASeq
        :return: Nothing
        :rtype: None
        """

        super(RNASeq, self).run()

        # RNASeq requires a genome version.

        if not self.genome_version:
            message = 'A RNASeq analysis requires a genome_version configuration option.'
            raise Exception(message)

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of BSF Sample objects by the BSF Sample name.
        # This cannot be done in the super-class, because BSF Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_TopHat_Cufflinks_jobs()
        self._create_Cuffmerge_Cuffdiff_jobs()

    def _create_TopHat_Cufflinks_jobs(self):

        """Create TopHat aligner jobs and Cufflinks transcript assembler jobs.

        :param self: BSF RNASeq
        :type self: RNASeq
        :return: Nothing
        :rtype: None
        """

        if self.debug > 1:
            print '{!r} method _create_TopHat_Cufflinks_jobs:'.format(self)

        # Read configuration options.

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')
        transcriptome = config_parser.get(section=config_section, option='transcriptome')
        # TODO: These really are properties of the Reads, PairedReads or Sample objects rather than an Analysis.
        insert_size = config_parser.getint(section=config_section, option='insert_size')
        read_length = config_parser.getint(section=config_section, option='read_length')

        genome_sizes = config_parser.get(section=config_section, option='genome_sizes')
        genome_sizes = os.path.expanduser(genome_sizes)
        genome_sizes = os.path.expandvars(genome_sizes)

        if config_parser.has_option(section=config_section, option='tophat_hold'):
            tophat_hold = config_parser.getboolean(section=config_section, option='tophat_hold')
        else:
            tophat_hold = False

        mate_inner_dist = insert_size - 2 * read_length

        genome_fasta = Default.absolute_genome_fasta(self.genome_version, 'bowtie2')

        # Get the Bowtie2 index

        bowtie2_index = os.path.join(Default.absolute_genomes(self.genome_version),
                                     'forBowtie2', self.genome_version)

        # Check if transcriptome is an absolute path and
        # prepend the annotation default if not.

        if not os.path.isabs(transcriptome):
            aga = Default.absolute_genome_annotation
            transcriptome = os.path.join(aga(self.genome_version), transcriptome)

        # Initialise the Distributed Resource Management System objects for
        # TopHat and Cufflinks Executable objects.

        tophat_drms = DRMS.from_Analysis(name='tophat',
                                         work_directory=self.genome_directory,
                                         analysis=self)
        self.drms_list.append(tophat_drms)

        process_tophat_drms = DRMS.from_Analysis(name='process_tophat',
                                                 work_directory=self.genome_directory,
                                                 analysis=self)
        self.drms_list.append(process_tophat_drms)

        cufflinks_drms = DRMS.from_Analysis(name='cufflinks',
                                            work_directory=self.genome_directory,
                                            analysis=self)
        self.drms_list.append(cufflinks_drms)

        process_cufflinks_drms = DRMS.from_Analysis(name='process_cufflinks',
                                                    work_directory=self.genome_directory,
                                                    analysis=self)
        self.drms_list.append(process_cufflinks_drms)

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

                # Create a new TopHat Executable.

                tophat = TopHat(name='rnaseq_tophat_{}'.format(replicate_key),
                                analysis=self)
                tophat.hold = tophat_hold

                tophat_drms.add_Executable(tophat)

                # Set TopHat options.

                tophat.add_OptionLong(key='GTF', value=transcriptome)
                tophat.add_OptionLong(key='output-dir',
                                      value=os.path.join(self.genome_directory, tophat.name))
                tophat.add_OptionLong(key='num-threads', value=str(tophat_drms.threads))
                # TODO: These really are properties of the Reads, PairedReads or Sample objects.
                tophat.add_OptionLong(key='mate-inner-dist', value=str(mate_inner_dist))
                if config_parser.has_option(section=config_section, option='mate-std-dev'):
                    tophat.add_OptionLong(key='mate-std-dev',
                                          value=config_parser.getint(section=config_section,
                                                                     option='mate-std-dev'))

                # Set TopHat arguments.

                tophat.arguments.append(bowtie2_index)

                # Set TopHat arguments for reads1 and reads2.

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:

                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

                tophat.arguments.append(string.join(reads1, ','))
                tophat.arguments.append(string.join(reads2, ','))

                # Create a new process_tophat Executable

                process_tophat = Executable.from_Analysis(name='rnaseq_process_tophat_{}'.format(replicate_key),
                                                          program='bsf_rnaseq_process_tophat2.sh',
                                                          analysis=self)

                process_tophat.dependencies.append(tophat.name)

                process_tophat_drms.add_Executable(process_tophat)

                # Set process_tophat options.

                # Set process_tophat arguments.

                process_tophat.arguments.append(os.path.join(self.genome_directory, tophat.name))
                process_tophat.arguments.append(genome_sizes)

                # Create a new Cufflinks Executable.

                cufflinks = Cufflinks(name='rnaseq_cufflinks_{}'.format(replicate_key), analysis=self)

                cufflinks.dependencies.append(tophat.name)

                cufflinks_drms.add_Executable(cufflinks)

                # Set Cufflinks options.

                cufflinks.add_OptionLong(key='output-dir', value=os.path.join(self.genome_directory, cufflinks.name))
                cufflinks.add_OptionLong(key='num-threads', value=str(cufflinks_drms.threads))

                # Cufflinks has a GTF option, in which case it will not assemble
                # novel transcripts and a GTF-guide option in which case it will
                # assemble novel transcripts.
                # TODO: The 'GTF-guide' versus 'GTF' options may have to be configurable ...
                cufflinks.add_OptionLong(key='GTF-guide', value=transcriptome)

                # TODO: The 'mask-file' options would be good to implement.
                # Annotated mitochondrial transcripts, rRNAs and other abundant
                # transcripts should be excluded to make abundance estimates more robust.

                cufflinks.add_OptionLong(key='frag-bias-correct', value=genome_fasta)

                # TODO: The 'multi-read-correct' option may have to be configurable.
                cufflinks.add_SwitchLong(key='multi-read-correct')

                # TODO: Explore other Advanced Abundance Estimation Options

                # Set Cufflinks arguments.

                aligned_reads = os.path.join(tophat.options['output-dir'].value, 'accepted_hits.bam')
                cufflinks.arguments.append(aligned_reads)

                # Create a new process_cufflinks Executable

                process_cufflinks = Executable.from_Analysis(name='rnaseq_process_cufflinks_{}'.format(replicate_key),
                                                             program='bsf_rnaseq_process_cufflinks.R',
                                                             analysis=self)

                process_cufflinks.dependencies.append(cufflinks.name)

                process_cufflinks_drms.add_Executable(process_cufflinks)

                # Set ProcessCufflinks options.

                # TODO: The BioMart data set needs to be configured somewhere ...
                # It is specific for the genome assembly.
                # For the moment it is set in the rnaseq_config.ini file.

                # TODO: For the moment, the species-specific BioMart_data_set option needs specifying
                # in the configuration file.
                # process_cufflinks.add_OptionLong(key='biomart_data_set', value='mmusculus_gene_ensembl)
                process_cufflinks.add_OptionLong(key='sample', value=replicate_key)
                process_cufflinks.add_OptionLong(key='genome_directory', value=self.genome_directory)

    def _create_Cuffmerge_Cuffdiff_jobs(self):

        """Create Cuffdiff differential expression jobs.

        :param self: BSF RNASeq
        :type self: RNASeq
        :return: Nothing
        :rtype: None
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')
        transcriptome = config_parser.get(section=config_section, option='transcriptome')

        # Check if transcriptome is an absolute path and
        # prepend the annotation default if not.

        if not os.path.isabs(transcriptome):
            aga = Default.absolute_genome_annotation
            transcriptome = os.path.join(aga(self.genome_version), transcriptome)

        genome_fasta = Default.absolute_genome_fasta(self.genome_version, 'bowtie2')

        cuffmerge_drms = DRMS.from_Analysis(name='cuffmerge',
                                            work_directory=self.genome_directory,
                                            analysis=self)

        self.drms_list.append(cuffmerge_drms)

        cuffdiff_drms = DRMS.from_Analysis(name='cuffdiff',
                                           work_directory=self.genome_directory,
                                           analysis=self)

        self.drms_list.append(cuffdiff_drms)

        process_cuffdiff_drms = DRMS.from_Analysis(name='process_cuffdiff',
                                                   work_directory=self.genome_directory,
                                                   analysis=self)

        self.drms_list.append(process_cuffdiff_drms)

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            c_samples, t_samples = self.comparisons[key]

            # Create a new Cuffmerge Executable.

            cuffmerge = Cuffmerge(name='rnaseq_cuffmerge_{}'.format(key),
                                  analysis=self)

            cuffmerge_drms.add_Executable(cuffmerge)

            # Set Cuffmerge options.

            cuffmerge.add_OptionLong(key='output-dir',
                                     value=os.path.join(self.genome_directory, cuffmerge.name))

            cuffmerge.add_OptionLong(key='num-threads', value=str(cuffmerge_drms.threads))

            cuffmerge.add_OptionLong(key='ref-gtf', value=transcriptome)

            cuffmerge.add_OptionLong(key='ref-sequence',
                                     value=Default.absolute_genome_fasta(self.genome_version, 'bowtie2'))

            # Set Cuffmerge arguments.

            # Create an assembly manifest file to merge all replicates of each BSF Sample ...

            assembly_name = cuffmerge.name + '_assembly.txt'
            assembly_path = os.path.join(self.genome_directory, assembly_name)
            assembly_file = open(name=assembly_path, mode='w')

            # Create a new Cuffdiff Executable.

            cuffdiff = Cuffdiff(name='rnaseq_cuffdiff_{}'.format(key), analysis=self)

            cuffdiff.dependencies.append(cuffmerge.name)

            cuffdiff_drms.add_Executable(cuffdiff)

            # Set Cuffdiff options.

            cuffdiff.add_OptionLong(key='output-dir',
                                    value=os.path.join(self.genome_directory, cuffdiff.name))

            cuffdiff.add_OptionLong(key='num-threads', value=str(cuffdiff_drms.threads))

            (c_name, t_name) = key.split('__')
            cuffdiff.add_OptionLong(key='labels', value='{},{}'.format(c_name, t_name))

            cuffdiff.add_OptionLong(key='frag-bias-correct', value=genome_fasta)

            c_alignments = list()
            t_alignments = list()

            for c_sample in c_samples:

                # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                # Python str key and Python list of Python list objects
                # of Bio.BSF.Data.PairedReads objects.

                c_replicate_dict = c_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                c_replicate_keys = c_replicate_dict.keys()
                c_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for c_replicate_key in c_replicate_keys:
                    # Add the Cufflinks assembled transcripts to the Cuffmerge manifest.
                    transcripts_path = os.path.join(self.genome_directory,
                                                    'rnaseq_cufflinks_{}'.format(c_replicate_key),
                                                    'transcripts.gtf')
                    assembly_file.write(transcripts_path + '\n')

                    # Wait for each TopHat and Cufflinks replicate to finish,
                    # before Cuffmerge can run.

                    cuffmerge.dependencies.append('rnaseq_cufflinks_{}'.format(c_replicate_key))

                    # Add the control BAM file to Cuffdiff ...

                    c_alignments.append(os.path.join(self.genome_directory,
                                                     'rnaseq_tophat_{}'.format(c_replicate_key),
                                                     'accepted_hits.bam'))

            for t_sample in t_samples:

                # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                # Python str key and Python list of Python list objects
                # of Bio.BSF.Data.PairedReads objects.

                t_replicate_dict = t_sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                t_replicate_keys = t_replicate_dict.keys()
                t_replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                for t_replicate_key in t_replicate_keys:
                    # Add the Cufflinks assembled transcripts to the Cuffmerge manifest.
                    transcripts_path = os.path.join(self.genome_directory,
                                                    'rnaseq_cufflinks_{}'.format(t_replicate_key),
                                                    'transcripts.gtf')
                    assembly_file.write(transcripts_path + '\n')

                    # Wait for each TopHat and Cufflinks replicate to finish,
                    # before Cuffmerge can run.

                    cuffmerge.dependencies.append('rnaseq_cufflinks_{}'.format(t_replicate_key))

                    # Add the treatment BAM file to Cuffdiff ...

                    t_alignments.append(os.path.join(self.genome_directory,
                                                     'rnaseq_tophat_{}'.format(t_replicate_key),
                                                     'accepted_hits.bam'))

            # Set Cuffmerge arguments.

            assembly_file.close()

            # Add the assembly manifest file as Cuffmerge argument.
            cuffmerge.arguments.append(assembly_path)

            # Add the Cuffmerge merged assembly as Cuffdiff output.
            cuffdiff.arguments.append(os.path.join(cuffmerge.options['output-dir'].value,
                                                   'merged.gtf'))

            # Add the aligned BAM files as Cuffdiff argument.
            cuffdiff.arguments.append(string.join(c_alignments, ','))
            cuffdiff.arguments.append(string.join(t_alignments, ','))

            # Create a new Process Cuffdiff Executable.

            process_cuffdiff = Executable.from_Analysis(name='rnaseq_process_cuffdiff_{}'.format(key),
                                                        program='bsf_rnaseq_process_cuffdiff.R',
                                                        analysis=self)

            process_cuffdiff.dependencies.append(cuffdiff.name)

            process_cuffdiff_drms.add_Executable(process_cuffdiff)

            # Set Process Cuffdiff options.

            process_cuffdiff.add_OptionLong(key='comparison', value=key)
            process_cuffdiff.add_OptionLong(key='genome_directory', value=self.genome_directory)

            # Set Process Cuffdiff arguments.

            # None.

    def report(self):

        """Create a BSF RNASeq report in HTML format and a UCSC Genome Browser Track Hub.

        :param self: BSF RNASeq
        :type self: RNASeq
        :return: Nothing
        :rtype: None
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        if config_parser.has_option(section=config_section, option='ucsc_location'):
            ucsc_location = config_parser.get(section=config_section, option='ucsc_location')
        else:
            ucsc_location = str()

        # Create a symbolic link containing the project name and a UUID.
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        # This code only needs the public URL.

        track_output = str()

        # Write a HTML document.

        output = str()

        output += Defaults.web.html_header(title='{} RNA-Seq Analysis'.format(self.project_name))
        output += '<body>\n'
        output += '\n'

        output += '<h1>{} RNA-Seq Analysis</h1>\n'.format(self.project_name)
        output += '\n'

        # TopHat and Cufflinks table.

        output += '<h2>Transcriptome Browsing</h2>\n'
        output += '\n'

        output += '<h3>Read Alignments</h3>\n'
        output += '\n'

        output += '<p>\n'
        output += '<strong><a href="http://tophat.cbcb.umd.edu/manual.html">TopHat</a></strong> '
        output += 'aligns RNA-Seq reads to a genome in order to identify '
        output += 'exon-exon splice junctions. It is built on the ultrafast\n'
        output += 'short read mapping program\n'
        output += '<strong><a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie 2</a></strong>.\n'
        output += '</p>\n'
        output += '\n'

        # Construct an automatic UCSC Track Hub link.

        options_dict = dict()
        options_dict['db'] = self.genome_version
        options_dict['hubUrl'] = '{}/{}/rnaseq_hub.txt'. \
            format(Default.url_absolute_projects(), link_name)

        output += '<p>\n'
        output += 'View TopHat <strong>read alignments</strong> tracks for each sample\n'
        output += 'in their genomic context via the project-specific\n'
        output += 'UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.\n'. \
            format(Defaults.web.ucsc_track_url(options_dict=options_dict, host_name=default.ucsc_host_name),
                   self.project_name)
        output += '</p>\n'
        output += '\n'

        output += '<h3>Splice Junctions, Insertions and Deletions</h3>\n'
        output += '\n'

        output += '<p>\n'
        output += 'TopHat reports <strong>splice junctions</strong> on the basis of RNA-Seq\n'
        output += 'read alignments in UCSC BED track format.\n'
        output += 'Each junction consists of two connected BED blocks,\n'
        output += 'where each block is as long as the maximal overhang\n'
        output += 'of any read spanning the junction. The score is\n'
        output += 'the number of alignments spanning the junction.\n'
        output += 'UCSC BED tracks of <strong>insertions</strong> and\n'
        output += '<strong>deletions</strong> are also reported by TopHat.\n'
        output += '</p>\n'

        output += '<p>\n'
        output += 'View the corresponding TopHat tracks for junctions, deletions and insertions\n'
        output += 'for each sample in their genomic context via the project-specific\n'
        output += 'UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.\n'. \
            format(Defaults.web.ucsc_track_url(options_dict=options_dict, host_name=default.ucsc_host_name),
                   self.project_name)
        output += '</p>\n'
        output += '\n'

        #output += '<p>\n'
        #output += 'Follow the links below to attach\n'
        #output += 'Tophat junction, deletion and insertion annotation to the\n'
        #output += 'UCSC Genome Browser. Since each file needs transferring to\n'
        #output += 'the UCSC site, subsequent pages will take some time to load.\n'
        #output += '</p>\n'

        output += '<h2>Gene Expression Profiles</h2>\n'
        output += '\n'

        output += '<p>\n'
        output += '<strong><a href="http://cufflinks.cbcb.umd.edu/howitworks.html">Cufflinks</a></strong>\n'
        output += 'assembles aligned RNA-Seq reads into transcripts,\n'
        output += 'estimates their abundances, and tests for differential\n'
        output += 'expression and regulation transcriptome-wide.\n'
        output += 'It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of\n'
        output += 'transcripts. Cufflinks then estimates the relative abundances of these transcripts based\n'
        output += 'on how many reads support each one, taking into account biases in library preparation protocols.\n'
        output += '</p>\n'

        output += '<p>\n'
        output += 'The Cufflinks <strong>assembled transcripts</strong> can be attached to the \n'
        output += 'UCSC Genome Browser, by following the "Transcript Assembly" links\n'
        output += 'below.\n'
        output += 'The isoforms.fpkm_tracking and genes.fpkm_tracking files\n'
        output += 'contain the estimated isoform or gene expression values in the generic\n'
        output += '<a href="http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format">FPKM Tracking format</a>.\n'
        output += '</p>\n'

        output += '<p>\n'
        output += 'Please see a more detailed description of Cufflinks\n'
        output += '<a href="http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output">output</a>.\n'
        output += '</p>\n'

        output += '<table>\n'
        output += '\n'
        output += '<tr>\n'
        output += '<th>Sample</th>\n'
        output += '<th>Assembled Transcripts</th>\n'
        output += '<th>Gene FPKM</th>\n'
        output += '<th>Transcript FPKM</th>\n'
        output += '<th>Genes (Symbols)</th>\n'
        output += '<th>Isoforms (Symbols)</th>\n'
        output += '</tr>\n'
        output += '\n'

        # Group via UCSC super tracks.

        track_output += 'track Alignments\n'
        track_output += 'shortLabel Alignments\n'
        track_output += 'longLabel Tophat2 alignments tracks\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Coverage\n'
        track_output += 'shortLabel Coverage\n'
        track_output += 'longLabel TopHat alignment coverage\n'
        track_output += 'visibility full\n'  # full or show?
        track_output += 'superTrack on\n'
        track_output += 'group coverage\n'
        track_output += '\n'

        track_output += 'track Deletions\n'
        track_output += 'shortLabel Deletions\n'
        track_output += 'longLabel Tophat2 deletions tracks\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Insertions\n'
        track_output += 'shortLabel Insertions\n'
        track_output += 'longLabel Tophat2 insertions tracks\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Junctions\n'
        track_output += 'shortLabel Junctions\n'
        track_output += 'longLabel Tophat2 exon junctions tracks\n'
        track_output += 'visibility show\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

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

                # TopHat produces accepted_hits.bam, deletions.bb,
                # insertions.bb and junctions.bb files.

                #
                # Add a trackDB entry for each accepted_hits.bam file.
                #

                # Common trackDb settings.

                track_output += 'track alignments_{}\n'. \
                    format(replicate_key)
                track_output += 'type bam\n'
                track_output += 'shortLabel alignments_{}\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} RNA-Seq read alignments\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./rnaseq_tophat_{}/accepted_hits.bam\n'. \
                    format(replicate_key)
                track_output += 'visibility dense\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                track_output += 'color {}\n'. \
                    format('0,0,0')

                # Compressed Sequence Alignment track settings.

                # None so far.

                # Composite track settings.

                track_output += 'parent Alignments\n'
                track_output += '\n'

                #
                # Add a trackDB entry for each accepted_hits.bw file.
                #

                # Common trackDB settings.

                track_output += 'track coverage_{}\n'. \
                    format(replicate_key)
                # TODO: The bigWig type must declare the expected signal range.
                # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                track_output += 'type bigWig\n'
                track_output += 'shortLabel coverage_{}\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} RNA-Seq alignment coverage\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./rnaseq_tophat_{}/accepted_hits.bw\n'. \
                    format(replicate_key)
                track_output += 'visibility full\n'
                # track_output += 'html {}\n'.format()

                # Common optional settings.

                track_output += 'color {}\n'. \
                    format('0,0,0')

                # bigWig - Signal graphing track settings.

                track_output += 'alwaysZero on\n'
                track_output += 'autoScale on\n'
                track_output += 'graphTypeDefault bar\n'
                track_output += 'maxHeightPixels 100:60:20\n'
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

                track_output += 'parent Coverage\n'
                track_output += 'centerLabelsDense off\n'
                track_output += '\n'

                #
                # Add a trackDB entry for each deletions.bb file.
                #

                track_output += 'track deletions_{}\n'. \
                    format(replicate_key)
                track_output += 'type bigBed\n'
                track_output += 'shortLabel deletions_{}\n'. \
                    format(replicate_key)
                track_output += 'longLabel RNA-Seq deletions for {}\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./rnaseq_tophat_{}/deletions.bb\n'. \
                    format(replicate_key)
                track_output += 'visibility hide\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                track_output += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                track_output += 'parent Deletions\n'
                track_output += '\n'

                # Insertions

                track_output += 'track insertions_{}\n'. \
                    format(replicate_key)
                track_output += 'type bigBed\n'
                track_output += 'shortLabel insertions_{}\n'. \
                    format(replicate_key)
                track_output += 'longLabel RNA-Seq insertions for {}\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./rnaseq_tophat_{}/insertions.bb\n'. \
                    format(replicate_key)
                track_output += 'visibility hide\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                track_output += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                track_output += 'parent Insertions\n'
                track_output += '\n'

                # Junctions

                track_output += 'track junctions_{}\n'. \
                    format(replicate_key)
                track_output += 'type bigBed\n'
                track_output += 'shortLabel junctions_{}\n'. \
                    format(replicate_key)
                track_output += 'longLabel RNA-Seq junctions for {}\n'. \
                    format(replicate_key)
                track_output += 'bigDataUrl ./rnaseq_tophat_{}/junctions.bb\n'. \
                    format(replicate_key)
                track_output += 'visibility pack\n'
                # 'html' is missing from the common settings.

                # Common optional settings.

                track_output += 'color {}\n'. \
                    format('0,0,0')

                # Composite track settings.

                track_output += 'parent Junctions\n'
                track_output += '\n'

                # Cufflinks produces genes.fpkm_tracking, isoforms.fpkm_tracking,
                # skipped.gtf and transcripts.gtf.

                # UCSC options dictionary.

                options_dict = dict()
                options_dict['db'] = self.genome_version
                if ucsc_location:
                    options_dict['position'] = ucsc_location

                # UCSC track dictionary.

                track_dict = dict()
                track_dict['name'] = 'transcripts_{}'. \
                    format(replicate_key)
                track_dict['description'] = '"Cufflinks transcript assembly of {}"'. \
                    format(replicate_key)
                track_dict['track_type'] = 'gtf'
                track_dict['visibility'] = '3'  # Visibility pack
                track_dict['color'] = '0,0,0'
                track_dict['db'] = self.genome_version

                prefix = 'rnaseq_cufflinks_{}'.format(replicate_key)

                output += '<tr>\n'
                output += '<td>{}</td>\n'.format(replicate_key)
                output += '<td><a href="{}">Transcript Assembly</a></td>\n'. \
                    format(Defaults.web.ucsc_track_url(options_dict=options_dict,
                                                       track_dict=track_dict,
                                                       host_name=default.ucsc_host_name))
                output += '<td><a href="./{}/genes.fpkm_tracking">Genes FPKM</a></td>\n'. \
                    format(prefix)
                output += '<td><a href="./{}/isoforms.fpkm_tracking">Isoforms FPKM</a></td>\n'. \
                    format(prefix)
                # Add files from bsf_process_cufflinks.R
                output += '<td><a href="./{}/{}_genes_fpkm_tracking.txt">Genes (Symbols)</a></td>\n'. \
                    format(prefix, prefix)
                output += '<td><a href="./{}/{}_isoforms_fpkm_tracking.txt">Isoforms (Symbols)</a></td>\n'. \
                    format(prefix, prefix)
                output += '</tr>\n'
                output += '\n'

        output += '</table>\n'
        output += '\n'

        # Cuffmerge produces merged.gtf.
        # GTF files would need converting into BED and bigBed to be compatible with UCSC Track Hubs.

        # Cuffdiff produces cds_exp.diff, gene_exp.diff, isoform_exp.diff
        # promoters.diff, splicing.diff and tss_group_exp.diff amongst many others.

        output += '<h2>Differential Expression</h2>\n'
        output += '\n'

        output += '<p>\n'
        output += '<strong><a href="http://cufflinks.cbcb.umd.edu/howitworks.html#diff">Cuffdiff</a></strong>\n'
        output += 'finds significant changes in transcript\n'
        output += 'expression, splicing, and promoter use.'
        output += '</p>\n'
        output += '\n'

        output += '<table>\n'
        output += '\n'

        output += '<tr>\n'
        output += '<th>Comparison</th>\n'
        output += '<th>Coding Sequences</th>\n'
        output += '<th>Genes</th>\n'
        output += '<th>Isoforms</th>\n'
        output += '<th>Promoters</th>\n'
        output += '<th>Splicing</th>\n'
        output += '<th>Transcription Start Sites</th>\n'
        output += '<th>Gene FPKM Replicates</th>\n'
        output += '<th>Isoform FPKM Replicates</th>\n'
        output += '</tr>\n'

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            prefix = 'rnaseq_process_cuffdiff_{}'.format(key)

            output += '<tr>\n'
            output += '<td>{}</td>\n'. \
                format(key)

            # output += '<td><a href="./rnaseq_cuffdiff_{}/cds_exp.diff">Coding Sequences</a></td>\n'. \
            #     format(key)
            # output += '<td><a href="./rnaseq_cuffdiff_{}/gene_exp.diff"><strong>Genes</strong></a></td>\n'. \
            #     format(key)
            # output += '<td><a href="./rnaseq_cuffdiff_{}/isoform_exp.diff">Isoforms</a></td>\n'. \
            #     format(key)
            # output += '<td><a href="./rnaseq_cuffdiff_{}/promoters.diff">Promoters</a></td>\n'. \
            #     format(key)
            # output += '<td><a href="./rnaseq_cuffdiff_{}/splicing.diff">Splicing</a></td>\n'. \
            #     format(key)
            # output += '<td><a href="./rnaseq_cuffdiff_{}/tss_group_exp.diff">Transcription Start Sites</a></td>\n'. \
            #     format(key)

            # Link to comparison-specific symbolic links in the directory after cummRbund processing.

            output += '<td><a href="./{}/{}_cds_exp_diff.txt">Coding Sequences</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_gene_exp_diff.txt"><strong>Genes</strong></a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_isoform_exp_diff.txt">Isoforms</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_promoters_diff.txt">Promoters</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_splicing_diff.txt">Splicing</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_tss_group_exp_diff.txt">Transcription Start Sites</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_genes_fpkm_replicates.txt">Gene FPKM Replicates</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_isoforms_fpkm_replicates.txt">Isoform FPKM Replicates</a></td>\n'. \
                format(prefix, prefix)

            output += '</tr>\n'

        output += '</table>\n'
        output += '\n'

        # Show cummeRbund quality plots.

        output += '<h2>Quality Plots</h2>\n'
        output += '\n'

        output += '<p>\n'
        output += '</p>\n'
        output += '\n'

        output += '<table>\n'

        output += '<tr>\n'
        output += '<th>Comparison</th>\n'
        output += '<th>Dispersion Plot - Genes</th>\n'
        output += '<th>Dispersion Plot - Isoforms</th>\n'
        output += '<th>Squared Coefficient of Variation - Genes</th>\n'
        output += '<th>Squared Coefficient of Variation - Isoforms</th>\n'
        output += '<th>Density Plot without Replicates - Genes</th>\n'
        output += '<th>Density Plot with Replicates - Genes</th>\n'
        output += '<th>Density Plot without Replicates - Isoforms</th>\n'
        output += '<th>Density Plot with Replicates - Isoforms</th>\n'
        output += '<th>Box Plot without Replicates - Genes</th>\n'
        output += '<th>Box Plot with Replicates - Genes</th>\n'
        # output += '<th>Box Plot without Replicates - Isoforms</th>\n'
        # output += '<th>Box Plot with Replicates - Isoforms</th>\n'
        output += '<th>Scatter Matrix Plot - Genes</th>\n'
        output += '<th>Scatter Matrix Plot - Isoforms</th>\n'
        output += '<th>Volcano Matrix Plot - Genes</th>\n'
        output += '<th>Multidimensional Scaling Plot - Genes</th>\n'
        output += '<th>Principal Component Analysis Plot - Genes</th>\n'
        output += '</tr>\n'

        for key in keys:

            prefix = 'rnaseq_process_cuffdiff_{}'.format(key)

            output += '<tr>\n'
            output += '<td>{}</td>\n'.format(key)

            # Dispersion Plots for Genes and Isoforms

            img_source = './{}/{}_genes_dispersion.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Dispersion Plot - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            img_source = './{}/{}_isoforms_dispersion.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Dispersion Plot - Isoforms - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            # Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms

            img_source = './{}/{}_genes_scv.png'.format(prefix, prefix)
            if os.path.exists(path=os.path.join(self.genome_directory, img_source)):
                output += '<td><a href="{}">'.format(img_source)
                output += '<img alt="Squared Coefficient of Variation (SCV) - Genes - {}" ' \
                          'src="{}" height="80" width="80" />'.format(key, img_source)
                output += '</a></td>\n'
            else:
                output += '<td></td>\n'

            img_source = './{}/{}_isoforms_scv.png'.format(prefix, prefix)
            if os.path.exists(path=os.path.join(self.genome_directory, img_source)):
                output += '<td><a href="{}">'.format(img_source)
                output += '<img alt="Squared Coefficient of Variation (SCV) - Isoforms - {}" ' \
                          'src="{}" height="80" width="80" />'.format(key, img_source)
                output += '</a></td>\n'
            else:
                output += '<td></td>\n'

            # Density Plots for Genes without and with Replicates

            img_source = './{}/{}_genes_density_wo_replicates.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Density Plot without Replicates - Genes- {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            img_source = './{}/{}_genes_density_w_replicates.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Density Plot with Replicates - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            # Density Plots for Isoforms without and with Replicates

            img_source = './{}/{}_isoforms_density_wo_replicates.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Density Plot without Replicates - Isoforms - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            img_source = './{}/{}_isoforms_density_w_replicates.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Density Plot with Replicates - Isoforms - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            # Box Plots for Genes without and with Replicates

            img_source = './{}/{}_genes_box_wo_replicates.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Box Plot without Replicates - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            img_source = './{}/{}_genes_box_w_replicates.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Box Plot with Replicates - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            # Box Plots for Isoforms with and without Replicates

            if 0:
                img_source = './{}/{}_isoforms_box_wo_replicates.png'.format(prefix, prefix)
                output += '<td><a href="{}">'.format(img_source)
                output += '<img alt="Box Plot without Replicates - Isoforms - {}" ' \
                          'src="{}" height="80" width="80" />'.format(key, img_source)
                output += '</a></td>\n'

                img_source = './{}/{}_isoforms_box_w_replicates.png'.format(prefix, prefix)
                output += '<td><a href="{}">'.format(img_source)
                output += '<img alt="Box Plot with Replicates - Isoforms - {}" ' \
                          'src="{}" height="80" width="80" />'.format(key, img_source)
                output += '</a></td>\n'

            # Scatter Matrix Plot for Genes and Isoforms

            img_source = './{}/{}_genes_scatter_matrix.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Scatter Matrix Plot - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            img_source = './{}/{}_isoforms_scatter_matrix.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Scatter Matrix Plot - Isoforms - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            # Volcano Matrix Plot for Genes

            img_source = './{}/{}_genes_volcano_matrix.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Volcano Matrix Plot - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            # Multidimensional Scaling Plot for Genes

            img_source = './{}/{}_genes_mds.png'.format(prefix, prefix)
            if os.path.exists(path=os.path.join(self.genome_directory, img_source)):
                output += '<td><a href="{}">'.format(img_source)
                output += '<img alt="Multidimensional Scaling Plot - Genes - {}" ' \
                          'src="{}" height="80" width="80" />'.format(key, img_source)
                output += '</a></td>\n'
            else:
                output += '<td></td>\n'

            # Principal Component Analysis Plot for Genes

            img_source = './{}/{}_genes_pca.png'.format(prefix, prefix)
            output += '<td><a href="{}">'.format(img_source)
            output += '<img alt="Principal Component Analysis Plot - Genes - {}" ' \
                      'src="{}" height="80" width="80" />'.format(key, img_source)
            output += '</a></td>\n'

            output += '</tr>\n'

        output += '</table>\n'

        output += '</body>\n'
        output += Defaults.web.html_footer()

        file_path = os.path.join(self.genome_directory, 'rnaseq_report.html')

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        # Create the UCSC Genome Browser Track Hub.

        self.ucsc_hub_write_hub(prefix='rnaseq')
        self.ucsc_hub_write_genomes(prefix='rnaseq')
        self.ucsc_hub_write_tracks(output=track_output, prefix='rnaseq')


class RunFastQC(Analysis):
    """BSF FastQC-specific Quality Assessment Analysis sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new BSF RunFastQC object from a UNIX-style configuration file via the BSF Configuration class.

        :param cls: Class
        :type cls: RunFastQC
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF RunFastQC
        :rtype: RunFastQC
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new BSF RunFastQC object from a BSF Configuration object.

        :param cls: Class
        :type cls: RunFastQC
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF RunFastQC
        :rtype: RunFastQC
        """

        assert isinstance(configuration, Configuration)

        run_fastqc = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.RunFastQC" section specifies defaults for this BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        run_fastqc.set_Configuration(run_fastqc.configuration, section=section)

        return run_fastqc

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):

        """Initialise a Bio.BSF.Analysis.RunFastQC object.

        :param self: BSF RunFastQC
        :type self: RunFastQC
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: BSF Analysis-wide input directory
        :type input_directory: str
        :param output_directory: BSF Analysis-wide output directory
        :type output_directory: str
        :param project_directory: BSF Analysis-wide project directory,
        normally under the BSF Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: BSF Analysis-wide genome directory,
        normally under the BSF Analysis-wide project directory
        :type genome_directory: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of BSF DRMS objects
        :type drms_list: list
        :param collection: BSF Collection
        :type collection: Collection
        :param comparisons: Python dict of Python tuple objects of BSF Sample objects
        :type comparisons: dict
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :return: Nothing
        :rtype: None
        """

        super(RunFastQC, self).__init__(configuration=configuration,
                                        project_name=project_name, genome_version=genome_version,
                                        input_directory=input_directory, output_directory=output_directory,
                                        project_directory=project_directory, genome_directory=genome_directory,
                                        e_mail=e_mail, debug=debug, drms_list=drms_list,
                                        collection=collection, comparisons=comparisons, samples=samples)

        # Nothing else to do for this sub-class ...

    def set_Configuration(self, configuration, section):

        """Set instance variables of a BSF RunFastQC object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF RunFastQC
        :type self: RunFastQC
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(RunFastQC, self).set_Configuration(configuration=configuration, section=section)

        # Read a comparison file.

        # if configuration.config_parser.has_option(section=section, option='cmp_file'):
        #     self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')
        # Use the sample annotation sheet instead of a separate comparison file.
        if configuration.config_parser.has_option(section=section, option='sas_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='sas_file')

    def _read_comparisons(self, cmp_file):

        """Read a BSF SampleAnnotationSheet CSV file from disk.

        Column headers for CASAVA folders:
          Treatment/Control ProcessedRunFolder:
            CASAVA processed run folder name or
            Bio.BSF.Analysis input_directory by default.
          Treatment/Control Project:
            CASAVA Project name or
            Bio.BSF.Analysis project_name by default.
          Treatment/Control Sample:
            CASAVA Sample name, no default.
        Column headers for independent samples:
          Treatment/Control Sample:
          Treatment/Control File:
        :param self: BSF RunFastQC
        :type self: RunFastQC
        :param cmp_file: Comparisons file path
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        sas = SampleAnnotationSheet(file_path=cmp_file)
        sas.csv_reader_open()

        for row_dict in sas._csv_reader:
            self.add_Sample(sample=self.collection.get_Sample_from_row_dict(row_dict=row_dict))

        sas.csv_reader_close()

    def run(self):

        """Run this BSF RunFastQC analysis.

        :param self: BSF RunFastQC analysis
        :type self: RunFastQC
        :return: Nothing
        :rtype: None
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

        """Initialise a BSF RunFastQC object.

        :param self: BSF RunFastQC
        :type self: RunFastQC
        :return: Nothing
        :rtype: None
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
                    message = "Could not find projects component {} in RunFolder {}.\n" \
                              "Check projects configuration value {!r}.". \
                        format(component, prf.name, projects)
                    warnings.warn(message, UserWarning)

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
            #     self.add_sample(sample=self.collection.get_Sample_from_row_dict(row_dict=row_dict))
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

                fastqc.arguments.append(string.join(reads1 + reads2, ' '))

        self.drms_list.append(fastqc_drms)

    def report(self):

        """Create a RunFastQC report in HTML format.

        :param self: BSF RunFastQC
        :type self: RunFastQC
        :return: Nothing
        :rtype: None
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

        :param cls: Class
        :type cls: RunBamToFastq
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF RunBamToFastq
        :rtype: RunBamToFastq
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new BSF RunBamToFastq object from a BSF Configuration object.

        :param cls: Class
        :type cls: RunBamToFastq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF RunBamToFastq
        :rtype: RunBamToFastq
        """

        assert isinstance(configuration, Configuration)

        run_bam_to_fastq = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.RunBamToFastq" section specifies defaults for this BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        run_bam_to_fastq.set_Configuration(run_bam_to_fastq.configuration, section=section)

        return run_bam_to_fastq

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):

        """Initialise a Bio.BSF.Analysis.RunBamToFastq object.

        :param self: BSF RunBamToFastq
        :type self: RunBamToFastq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: BSF Analysis-wide input directory
        :type input_directory: str
        :param output_directory: BSF Analysis-wide output directory
        :type output_directory: str
        :param project_directory: BSF Analysis-wide project directory,
        normally under the BSF Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: BSF Analysis-wide genome directory,
        normally under the BSF Analysis-wide project directory
        :type genome_directory: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of BSF DRMS objects
        :type drms_list: list
        :param collection: BSF Collection
        :type collection: Collection
        :param comparisons: Python dict of Python list objects of BSF Sample objects
        :type comparisons: dict
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :return: Nothing
        :rtype: None
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
        :param self: BSF RunBamToFastq
        :type self: RunBamToFastq
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(RunBamToFastq, self).set_Configuration(configuration=configuration, section=section)

        # Nothing else to do for this sub-class ...

    def run(self):

        """Run this BSF RunBamToFastq analysis.

        :param self: BSF RunBamToFastq analysis
        :type self: RunBamToFastq
        :return: Nothing
        :rtype: None
        """

        super(RunBamToFastq, self).run()

        self._convert_BamToFastq()

    def _convert_BamToFastq(self):

        """Private method to convert all Reads objects in BAM or SAM format into FASTQ format.

        :param self: BSF RunBamToFastq
        :type self: RunBamToFastq
        :return: Nothing
        :rtype: None
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


class VariantCalling(Analysis):
    """BSF Variant Calling-specific BSF Analysis sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new BSF VariantCalling object from a UNIX-style configuration file via the BSF Configuration class.

        :param cls: Class
        :type cls: VariantCalling
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF VariantCalling
        :rtype: VariantCalling
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new BSF VariantCalling object from a BSF Configuration object.

        :param cls: Class
        :type cls: VariantCalling
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF VariantCalling
        :rtype: VariantCalling
        """

        assert isinstance(configuration, Configuration)

        variant_calling = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.VariantCalling" section specifies defaults for this BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        variant_calling.set_Configuration(variant_calling.configuration, section=section)

        return variant_calling

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 cmp_file=None):
        """Initialise a Bio.BSF.Analysis.VariantCalling object.

        :param self: BSF VariantCalling
        :type self: VariantCalling
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: BSF Analysis-wide input directory
        :type input_directory: str
        :param output_directory: BSF Analysis-wide output directory
        :type output_directory: str
        :param project_directory: BSF Analysis-wide project directory,
         normally under the BSF Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: BSF Analysis-wide genome directory,
         normally under the BSF Analysis-wide project directory
        :type genome_directory: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of BSF DRMS objects
        :type drms_list: list
        :param collection: BSF Collection
        :type collection: Collection
        :param comparisons: Python dict of Python list objects of BSF Sample objects
        :type comparisons: dict
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :param cmp_file: Comparison file
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        super(VariantCalling, self).__init__(configuration=configuration,
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
        """Set instance variables of a BSF VariantCalling object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF VariantCalling
        :type self: VariantCalling
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(VariantCalling, self).set_Configuration(configuration=configuration, section=section)

        # Read a comparison file.

        # if configuration.config_parser.has_option(section=section, option='cmp_file'):
        #     self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')
        # Use the sample annotation sheet instead of a separate comparison file.
        if configuration.config_parser.has_option(section=section, option='sas_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='sas_file')

    def _read_comparisons(self, cmp_file):
        """Read a BSF SampleAnnotationSheet CSV file from disk.

        Column headers for CASAVA folders:
          Treatment/Control ProcessedRunFolder:
            CASAVA processed run folder name or
            Bio.BSF.Analysis input_directory by default.
          Treatment/Control Project:
            CASAVA Project name or
            Bio.BSF.Analysis project_name by default.
          Treatment/Control Sample:
            CASAVA Sample name, no default.
        Column headers for independent samples:
          Treatment/Control Sample:
          Treatment/Control File:
        :param self: BSF RunFastQC
        :type self: RunFastQC
        :param cmp_file: Comparisons file path
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        sas = SampleAnnotationSheet(file_path=cmp_file)
        sas.csv_reader_open()

        for row_dict in sas._csv_reader:
            self.add_Sample(sample=self.collection.get_Sample_from_row_dict(row_dict=row_dict))

        sas.csv_reader_close()

    def run(self):
        """Run this BSF Variant Calling analysis.

        :param self: BSF VariantCalling
        :type self: VariantCalling
        :return: Nothing
        :rtype: None
        """

        super(VariantCalling, self).run()

        # Variant Calling requires a genome version.

        if not self.genome_version:
            message = 'A Variant Calling analysis requires a genome_version configuration option.'
            raise Exception(message)

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        # Real comparisons would be required for somatic mutation calling.
        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of BSF Sample objects by the BSF Sample name.
        # This cannot be done in the super-class, because BSF Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_alignment_jobs()

    def _create_alignment_jobs(self):
        """Create BWA alignment and post-processing jobs.

        The GATK documentation recommends mapping, de-duplicating, re-aligning and re-calibrating
        reads per sample and per lane.
        :param self: BSF Analysis
        :type self: Analysis
        :return: Nothing
        :rtype: None
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        bwa_genome_db = config_parser.get(section=config_section, option='bwa_genome_db')

        classpath_gatk = str()
        if config_parser.has_option(section=config_section, option='classpath_gatk'):
            classpath_gatk = config_parser.get(section=config_section, option='classpath_gatk')

        classpath_picard = str()
        if config_parser.has_option(section=config_section, option='classpath_picard'):
            classpath_picard = config_parser.get(section=config_section, option='classpath_picard')

        known_sites_realignment = list()
        if config_parser.has_option(section=config_section, option='known_sites_realignment'):
            temporary_str = config_parser.get(section=config_section, option='known_sites_realignment')
            known_sites_realignment.extend(temporary_str.split(','))

        # TODO: Prepend the GATK bundle default path.

        known_sites_recalibration = list()
        if config_parser.has_option(section=config_section, option='known_sites_recalibration'):
            temporary_str = config_parser.get(section=config_section, option='known_sites_recalibration')
            known_sites_recalibration.extend(temporary_str.split(','))

        # TODO: Prepend the GATK bundle default path.

        known_sites_discovery = list()
        if config_parser.has_option(section=config_section, option='known_sites_discovery'):
            temporary_str = config_parser.get(section=config_section, option='known_sites_discovery')
            known_sites_discovery.extend(temporary_str.split(','))

        # TODO: Prepend the GATK bundle default path.

        # path_known_insertions_deletions_list = list()
        # if config_parser.has_option(section=config_section, option='path_known_insertions_deletions'):
        #     path_known_insertions_deletions_str = config_parser.get(
        #         section=config_section,
        #         option='path_known_insertions_deletions')
        #     path_known_insertions_deletions_list.extend(path_known_insertions_deletions_str.split(','))
        #
        # for path_known_insertions_deletions in path_known_insertions_deletions_list:
        #     if os.path.isabs(path_known_insertions_deletions):
        #         # TODO: Prepend the GATK bundle default path.
        #         pass
        #
        # path_known_sites_list = list()
        # if config_parser.has_option(section=config_section, option='path_known_sites'):
        #     path_known_sites_str = config_parser.get(
        #         section=config_section,
        #         option='path_known_sites')
        #     path_known_sites_list.extend(path_known_sites_str.split(','))
        #
        # for path_known_sites in path_known_sites_list:
        #     if os.path.isabs(path_known_sites):
        #         # TODO: Prepend the GATK bundle default path.
        #         pass

        # Initialise the Distributed Management System objects for the run_bwa script.

        vc_align_lane_drms = DRMS.from_Analysis(name='variant_calling_align_lane',
                                                work_directory=self.genome_directory,
                                                analysis=self)
        self.drms_list.append(vc_align_lane_drms)

        # Initialise the Distributed Resource Management System object for the
        # bsf_run_variant_calling_processing_lane.py script.

        run_vc_process_lane_drms = DRMS.from_Analysis(name='variant_calling_process_lane',
                                                      work_directory=self.genome_directory,
                                                      analysis=self)
        self.drms_list.append(run_vc_process_lane_drms)

        # Initialise the Distributed Resource Management System object for the
        # bsf_run_variant_calling_processing_sample.py script.

        vc_process_sample_drms = DRMS.from_Analysis(name='variant_calling_process_sample',
                                                    work_directory=self.genome_directory,
                                                    analysis=self)
        self.drms_list.append(vc_process_sample_drms)

        # Initialise the Distributed Resource Management System object for the
        # bsf_run_variant_calling_process_cohort.py script.

        vc_process_cohort_drms = DRMS.from_Analysis(name='variant_calling_process_cohort',
                                                    work_directory=self.genome_directory,
                                                    analysis=self)
        self.drms_list.append(vc_process_cohort_drms)

        cohort_name = 'cohort'
        vc_process_cohort_dependencies = list()
        vc_process_cohort_replicates = list()

        for sample in self.samples:

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # Go back to the BSF Project to get its name and set is as the RG library name.
            # This is obsolete if the read group line can be set from the initial BAM file.
            # project = sample.weak_reference_project()
            # if project is None:
            #     library_name = 'Default'
            # else:
            #     library_name = project.name

            vc_process_sample_dependencies = list()
            vc_process_sample_replicates = list()

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                # TODO: Not sure such a file_path dictionary works.
                # file_path_dict = dict()
                # file_path_dict['aligned_bam'] = '{}{}.bam'.format(prefix, replicate_key)
                # file_path_dict['aligned_bai'] = '{}{}.bai'.format(prefix, replicate_key)
                # file_path_dict['aligned_md5'] = '{}{}.bam.md5'.format(prefix, replicate_key)
                # file_path_dict['duplicates_marked_bam'] = '{}{}_duplicates_marked.bam'.format(prefix, replicate_key)
                # file_path_dict['path_duplicates_marked_bai'] = '{}{}_duplicates_marked.bai'.format(prefix, replicate_key)
                # file_path_dict['duplicates_marked_md5'] = '{}{}_duplicates_marked.bam.md5'.format(prefix, replicate_key)
                # file_path_dict['realigner_targets'] = '{}{}_realigner.intervals'.format(prefix, replicate_key)
                # file_path_dict['realigned_bam'] = '{}{}_realigned.bam'.format(prefix, replicate_key)
                # file_path_dict['realigned_bai'] = '{}{}_realigned.bai'.format(prefix, replicate_key)
                # file_path_dict['realigned_md5'] = '{}{}_realigned.bam.md5'.format(prefix, replicate_key)
                # file_path_dict['recalibration_table'] = '{}{}_recalibration.table'.format(prefix, replicate_key)
                # file_path_dict['recalibrated_bam'] = '{}{}_recalibrated.bam'.format(prefix, replicate_key)
                # file_path_dict['recalibrated_bai'] = '{}{}_recalibrated.bai'.format(prefix, replicate_key)
                # file_path_dict['recalibrated_md5'] = '{}{}_recalibrated.bam.md5'.format(prefix, replicate_key)
                # file_path_dict['alignment_summary_metrics'] = '{}{}_alignment_summary_metrics.csv'.format(prefix, replicate_key)

                bwa = BWA(name='variant_calling_bwa_{}'.format(replicate_key), analysis=self)
                # Instead of adding the BWA Executable to the DRMS, it gets serialised into the pickler_file.
                # bwa_drms.add_Executable(bwa)

                bwa_mem = bwa.sub_command

                # Set BWA mem options.

                # bwa_mem.add_OptionShort(key='t', value=str(vc_align_lane_drms.threads))
                bwa_mem.add_OptionShort(key='t', value=str(1))  # TODO: For the moment, use only one thread.
                bwa_mem.add_SwitchShort(key='C')  # Append FASTA/Q comment to SAM output.
                bwa_mem.add_SwitchShort(key='M')  # Mark  shorter split hits as secondary (for Picard compatibility).
                # The @RG line, including the platform unit (PU) gets taken from the original BAM file.
                # bwa_mem.add_OptionShort(
                #     key='R',
                #     value="@RG\\tID:{}\\tSM:{}\\tLB:{}\\tPL:ILLUMINA". \
                #           format(replicate_key, sample.name, library_name))
                bwa_mem.add_OptionShort(key='v', value='2')  # Output warnings and errors only.

                # Set BWA arguments.

                bwa_mem.arguments.append(bwa_genome_db)

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:
                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

                if len(reads1) and not len(reads2):
                    bwa_mem.arguments.append(string.join(reads1, ','))
                elif len(reads1) and len(reads2):
                    bwa_mem.arguments.append(string.join(reads1, ','))
                    bwa_mem.arguments.append(string.join(reads2, ','))
                if len(reads2):
                    warning = 'Only second reads, but no first reads have been defined.'
                    warnings.warn(warning)

                # Normally, the bwa object would be pushed onto the drms list.
                # Experimentally, use Pickler to serialize the Executable object into a file.

                pickler_dict_align_lane = dict()
                pickler_dict_align_lane['prefix'] = vc_align_lane_drms.name
                pickler_dict_align_lane['replicate_key'] = replicate_key
                pickler_dict_align_lane['classpath_gatk'] = classpath_gatk
                pickler_dict_align_lane['classpath_picard'] = classpath_picard
                pickler_dict_align_lane['bwa_executable'] = bwa

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(vc_align_lane_drms.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_align_lane)
                pickler_file.close()

                # Create a bsf_run_bwa.py job to run the pickled object.

                run_bwa = Executable.from_Analysis(
                    name='{}_{}'.format(vc_align_lane_drms.name, replicate_key),
                    program='bsf_run_bwa.py',
                    analysis=self)
                vc_align_lane_drms.add_Executable(executable=run_bwa)

                # Set run_bwa options.

                run_bwa.add_OptionLong(key='pickler_path', value=pickler_path)
                run_bwa.add_OptionLong(key='debug', value=str(self.debug))

                # Create a bsf_run_variant_calling_process_lane.py job.

                pickler_dict_process_lane = dict()
                pickler_dict_process_lane['prefix'] = run_vc_process_lane_drms.name
                pickler_dict_process_lane['replicate_key'] = replicate_key
                pickler_dict_process_lane['classpath_gatk'] = classpath_gatk
                pickler_dict_process_lane['classpath_picard'] = classpath_picard
                # The 'path_replicate' gets defined in the bsf_run_variant_calling_align_lane.py script.
                pickler_dict_process_lane['path_replicate'] = '{}_{}.bam'. \
                    format(vc_align_lane_drms.name, replicate_key)
                pickler_dict_process_lane['path_reference_sequence'] = bwa_genome_db
                pickler_dict_process_lane['known_sites_realignment'] = known_sites_realignment
                pickler_dict_process_lane['known_sites_recalibration'] = known_sites_recalibration

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(run_vc_process_lane_drms.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_process_lane)
                pickler_file.close()

                vc_process_lane = Executable.from_Analysis(
                    name='{}_{}'.format(run_vc_process_lane_drms.name, replicate_key),
                    program='bsf_run_variant_calling_process_lane.py',
                    analysis=self)
                vc_process_lane.dependencies.append(run_bwa.name)
                run_vc_process_lane_drms.add_Executable(vc_process_lane)

                # Set variant_calling_run_process_lane options.

                vc_process_lane.add_OptionLong(key='pickler_path', value=pickler_path)
                vc_process_lane.add_OptionLong(key='debug', value=str(self.debug))

                # Set dependencies for the next stage.
                vc_process_sample_dependencies.append(vc_process_lane.name)
                # Add the result of the bsf_run_variant_calling_process_lane.py script.
                vc_process_sample_replicates.append(
                    '{}_{}_recalibrated.bam'.format(run_vc_process_lane_drms.name, replicate_key))

            # Finally, write the pickler file for processing per sample.

            pickler_dict_process_sample = dict()
            pickler_dict_process_sample['prefix'] = vc_process_sample_drms.name
            pickler_dict_process_sample['sample_key'] = sample.name
            pickler_dict_process_sample['replicate_file_list'] = vc_process_sample_replicates
            pickler_dict_process_sample['path_reference_sequence'] = bwa_genome_db
            pickler_dict_process_sample['known_sites_realignment'] = known_sites_realignment
            pickler_dict_process_sample['known_sites_recalibration'] = known_sites_recalibration
            pickler_dict_process_sample['known_sites_discovery'] = known_sites_discovery
            pickler_dict_process_sample['classpath_gatk'] = classpath_gatk
            pickler_dict_process_sample['classpath_picard'] = classpath_picard

            pickler_path = os.path.join(
                self.genome_directory,
                '{}_{}.pkl'.format(vc_process_sample_drms.name, sample.name))
            pickler_file = open(pickler_path, 'wb')
            pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
            pickler.dump(obj=pickler_dict_process_sample)
            pickler_file.close()

            vc_process_sample = Executable.from_Analysis(
                name='{}_{}'.format(vc_process_sample_drms.name, sample.name),
                program='bsf_run_variant_calling_process_sample.py',
                analysis=self)
            vc_process_sample.dependencies.extend(vc_process_sample_dependencies)
            vc_process_sample_drms.add_Executable(vc_process_sample)

            # Set variant_calling_run_process_sample options.

            vc_process_sample.add_OptionLong(key='pickler_path', value=pickler_path)
            vc_process_sample.add_OptionLong(key='debug', value=str(self.debug))

            # Set dependencies for the next stage.
            vc_process_cohort_dependencies.append(vc_process_sample.name)
            # Add the result of the bsf_run_variant_calling_process_sample.py script.
            vc_process_cohort_replicates.append(
                '{}_{}_raw_variants.vcf'.format(vc_process_cohort_drms.name, sample.name))

        # Finally, write the Pickler dict file for processing the cohort.

        pickler_dict_process_cohort = dict()
        pickler_dict_process_cohort['prefix'] = vc_process_sample_drms.name
        pickler_dict_process_cohort['cohort_key'] = cohort_name
        pickler_dict_process_cohort['replicate_file_list'] = vc_process_cohort_replicates
        pickler_dict_process_cohort['path_reference_sequence'] = bwa_genome_db
        pickler_dict_process_cohort['known_sites_realignment'] = known_sites_realignment
        pickler_dict_process_cohort['known_sites_recalibration'] = known_sites_recalibration
        pickler_dict_process_cohort['known_sites_discovery'] = known_sites_discovery
        pickler_dict_process_cohort['classpath_gatk'] = classpath_gatk
        pickler_dict_process_cohort['classpath_picard'] = classpath_picard

        pickler_path = os.path.join(
            self.genome_directory,
            '{}_{}.pkl'.format(vc_process_cohort_drms.name, cohort_name))
        pickler_file = open(pickler_path, 'wb')
        pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
        pickler.dump(obj=pickler_dict_process_cohort)
        pickler_file.close()

        vc_process_cohort = Executable.from_Analysis(
            name='{}_{}'.format(vc_process_cohort_drms.name, cohort_name),
            program='bsf_run_variant_calling_process_cohort.py',
            analysis=self)
        vc_process_cohort.dependencies.extend(vc_process_cohort_dependencies)
        vc_process_cohort_drms.add_Executable(vc_process_cohort)

        # Set variant_calling_run_process_sample options.

        vc_process_cohort.add_OptionLong(key='pickler_path', value=pickler_path)
        vc_process_cohort.add_OptionLong(key='debug', value=str(self.debug))
