"""Bio.BSF.Analysis.RNASeq

A package of classes and methods supporting RNA-Seq analyses.
"""

#
# Copyright 2014 Michael K. Schuster
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

import os.path
from pickle import Pickler, HIGHEST_PROTOCOL
import string

from Bio.BSF import Analysis, Configuration, Default, Defaults, DRMS, Executable
from Bio.BSF.Data import SampleAnnotationSheet
from Bio.BSF.Executables import Cuffdiff, Cufflinks, Cuffmerge, TopHat


class Tuxedo(Analysis):
    """Tuxedo RNASeq Analysis sub-class.

    Attributes:
    :ivar cmp_file: Comparison file
    :type cmp_file: str, unicode
    """

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new Tuxedo object from a UNIX-style configuration file via the Configuration class.

        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: Tuxedo
        :rtype: Tuxedo
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new Tuxedo object from a Configuration object.

        :param configuration: Configuration
        :type configuration: Configuration
        :return: Tuxedo
        :rtype: Tuxedo
        """

        assert isinstance(configuration, Configuration)

        rnaseq = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.RNASeq.Tuxedo" section specifies defaults for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        rnaseq.set_Configuration(rnaseq.configuration, section=section)

        return rnaseq

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None,
                 cmp_file=None):

        """Initialise a Bio.BSF.Analysis.RNASeq.Tuxedo object.

        :param configuration: Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: Analysis-wide input directory
        :type input_directory: str
        :param output_directory: Analysis-wide output directory
        :type output_directory: str
        :param project_directory: Analysis-wide project directory,
        normally under the Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: Analysis-wide genome directory,
        normally under the Analysis-wide project directory
        :type genome_directory: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of DRMS objects
        :type drms_list: list
        :param collection: Collection
        :type collection: Collection
        :param comparisons: Python dict of Python tuple objects of Sample objects
        :type comparisons: dict
        :param samples: Python list of Sample objects
        :type samples: list
        :param cmp_file: Comparison file
        :type cmp_file: str, unicode
        :return: Nothing
        :rtype: None
        """

        super(Tuxedo, self).__init__(configuration=configuration,
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

        """Set instance variables of a Tuxedo object via a section of a Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param configuration: Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(Tuxedo, self).set_Configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        if configuration.config_parser.has_option(section=section, option='cmp_file'):
            self.cmp_file = configuration.config_parser.get(section=section, option='cmp_file')

    def _read_comparisons(self, cmp_file):

        """Read a SampleAnnotationSheet CSV file specifying comparisons from disk.

        Column headers for CASAVA folders:
          Treatment/Control/Point N ProcessedRunFolder:
            CASAVA processed run folder name or
            Analysis input_directory by default.
          Treatment/Control/Point N Project:
            CASAVA Project name or
            Analysis project_name by default.
          Treatment/Control/Point N Sample:
            CASAVA Sample name, no default.
        Column headers for independent samples:
          Treatment/Control/Point N Sample:
          Treatment/Control/Point N File:
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

            # In addition to defining samples, allow also the definition of groups in comparison files.
            # If the row dictionary has a 'Group' key, then the Sample in the same row gets added to the group.
            # So, 'ProcessedRunFolder','Project','Sample','Group' defines the groups, while ...
            # 'Control Group','Treatment Group' defines a comparison, as does ...
            # 'Control Group','Treatment ProcessedRunFolder','Treatment Project','Treatment Sample'

            # Get Sample objects for classical 'Control' and 'Treatment' keys,
            # before looking up a series of 'Point N' keys.
            i = -2
            key = str()
            comparison_groups = list()
            while 1:
                i += 1
                if i == -1:
                    prefix = 'Control'
                elif i == 0:
                    prefix = 'Treatment'
                else:
                    prefix = 'Point {}'.format(i)
                # Get Sample objects for 'Point N' keys for as long as they are defined.
                # The Bio.BSF.Data.Collection.get_Sample_from_row_dict method can return one or more Sample objects,
                # depending on 'Group' or 'Sample' column entries.
                # In RNA-Seq experiments, entire pools of Sample objects (replicates) are compared with each other.
                group_name, group_samples = self.collection.get_Samples_from_row_dict(row_dict=row_dict, prefix=prefix)
                if group_name and len(group_samples):
                    key += group_name
                    key += '__'
                    comparison_groups.append((group_name, group_samples))
                    # Also expand each Python list of Sample objects to get all those Sample objects
                    # that this Analysis needs considering.
                    for sample in group_samples:
                        if self.debug > 1:
                            print '  {} Sample name: {!r} file_path:{!r}'.format(prefix, sample.name, sample.file_path)
                            # print sample.trace(1)
                        self.add_Sample(sample=sample)
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

            # Truncate the last '__' separator off the key string.

            self.comparisons[key[:-2]] = comparison_groups

        sas.csv_reader_close()

    def run(self):

        """Run this Tuxedo analysis.

        :return: Nothing
        :rtype: None
        """

        super(Tuxedo, self).run()

        # Tuxedo requires a genome version.

        if not self.genome_version:
            raise Exception('A Tuxedo analysis requires a genome_version configuration option.')

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
        # This cannot be done in the super-class, because Sample objects are only put into the
        # Analysis.samples list by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        self._create_TopHat_Cufflinks_jobs()
        self._create_Cuffmerge_Cuffdiff_jobs()

    def _create_TopHat_Cufflinks_jobs(self):

        """Create TopHat aligner jobs and Cufflinks transcript assembler jobs.

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

        if config_parser.has_option(section=config_section, option='genome_fasta'):
            genome_fasta = config_parser.get(section=config_section, option='genome_fasta')
        else:
            genome_fasta = Default.absolute_genome_fasta(self.genome_version, 'bowtie2')

        # Get the Bowtie2 index

        if config_parser.has_option(section=config_section, option='bowtie2_index'):
            bowtie2_index = config_parser.get(section=config_section, option='bowtie2_index')
        else:
            bowtie2_index = os.path.join(Default.absolute_genomes(self.genome_version),
                                         'forBowtie2', self.genome_version)

        if config_parser.has_option(section=config_section, option='novel_transcripts'):
            novel_transcripts = config_parser.getboolean(section=config_section, option='novel_transcripts')
        else:
            novel_transcripts = True

        # Check if transcriptome is an absolute path and
        # prepend the annotation default if not.

        if not os.path.isabs(transcriptome):
            aga = Default.absolute_genome_annotation
            transcriptome = os.path.join(aga(self.genome_version), transcriptome)

        # Initialise the Distributed Resource Management System (DRMS) objects for
        # TopHat and Cufflinks Executable objects.

        run_tophat_drms = DRMS.from_Analysis(
            name='rnaseq_run_tophat',
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(run_tophat_drms)

        process_tophat_drms = DRMS.from_Analysis(
            name='rnaseq_process_tophat',
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(process_tophat_drms)

        run_cufflinks_drms = DRMS.from_Analysis(
            name='rnaseq_run_cufflinks',
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(run_cufflinks_drms)

        process_cufflinks_drms = DRMS.from_Analysis(
            name='rnaseq_process_cufflinks',
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

                # Create a new rnaseq_tophat Executable, which is run via the rnaseq_run_tophat Executable below.

                tophat = TopHat(
                    name=string.join(words=('rnaseq_tophat', replicate_key), sep='_'),
                    analysis=self)
                tophat.hold = tophat_hold

                # Set rnaseq_tophat options.

                tophat.add_OptionLong(key='GTF',
                                      value=transcriptome)
                tophat.add_OptionLong(key='output-dir',
                                      value=os.path.join(self.genome_directory, tophat.name))
                tophat.add_OptionLong(key='num-threads',
                                      value=str(run_tophat_drms.threads))
                # TODO: These really are properties of the Reads, PairedReads or Sample objects.
                tophat.add_OptionLong(key='mate-inner-dist',
                                      value=str(mate_inner_dist))
                if config_parser.has_option(section=config_section, option='mate-std-dev'):
                    tophat.add_OptionLong(key='mate-std-dev',
                                          value=config_parser.getint(section=config_section,
                                                                     option='mate-std-dev'))

                # Set rnaseq_tophat arguments.

                tophat.arguments.append(bowtie2_index)

                # Set rnaseq_tophat arguments for reads1 and reads2.

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:

                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

                tophat.arguments.append(string.join(words=reads1, sep=','))
                tophat.arguments.append(string.join(words=reads2, sep=','))

                # Create a new rnaseq_run_tophat Executable.

                pickler_dict_run_tophat = dict()
                pickler_dict_run_tophat['prefix'] = run_tophat_drms.name
                pickler_dict_run_tophat['replicate_key'] = replicate_key
                pickler_dict_run_tophat['tophat_executable'] = tophat

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(run_tophat_drms.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_run_tophat)
                pickler_file.close()

                run_tophat = Executable.from_Analysis(
                    name=string.join(words=(run_tophat_drms.name, replicate_key), sep='_'),
                    program='bsf_run_rnaseq_tophat.py',
                    analysis=self)
                run_tophat_drms.add_Executable(run_tophat)
                # The rnaseq_run_tophat Executable requires no DRMS.dependencies.

                # Set rnaseq_run_tophat options.

                run_tophat.add_OptionLong(key='pickler_path', value=pickler_path)
                run_tophat.add_OptionLong(key='debug', value=str(self.debug))

                # Create a new rnaseq_process_tophat Executable

                process_tophat = Executable.from_Analysis(
                    name=string.join(words=(process_tophat_drms.name, replicate_key), sep='_'),
                    program='bsf_rnaseq_process_tophat2.sh',
                    analysis=self)
                process_tophat_drms.add_Executable(process_tophat)
                process_tophat.dependencies.append(run_tophat.name)

                # Set rnaseq_process_tophat options.

                # Set rnaseq_process_tophat arguments.

                process_tophat.arguments.append(os.path.join(self.genome_directory, tophat.name))
                process_tophat.arguments.append(genome_sizes)

                # Create a new rnaseq_cufflinks Executable, which is run via the rnaseq_run_cufflinks Executable below.

                cufflinks = Cufflinks(
                    name=string.join(words=('rnaseq_cufflinks', replicate_key), sep='_'),
                    analysis=self)

                # Set rnaseq_cufflinks options.

                cufflinks.add_OptionLong(key='output-dir',
                                         value=os.path.join(self.genome_directory, cufflinks.name))
                cufflinks.add_OptionLong(key='num-threads',
                                         value=str(run_cufflinks_drms.threads))

                # Cufflinks has a GTF option, in which case it will not assemble
                # novel transcripts and a GTF-guide option in which case it will
                # assemble novel transcripts.

                if novel_transcripts:
                    cufflinks.add_OptionLong(key='GTF-guide',
                                             value=transcriptome)
                else:
                    cufflinks.add_OptionLong(key='GTF',
                                             value=transcriptome)

                # TODO: The 'mask-file' options would be good to implement.
                # Annotated mitochondrial transcripts, rRNAs and other abundant
                # transcripts should be excluded to make abundance estimates more robust.

                cufflinks.add_OptionLong(key='frag-bias-correct',
                                         value=genome_fasta)

                # TODO: The 'multi-read-correct' option may have to be configurable.
                cufflinks.add_SwitchLong(key='multi-read-correct')

                # TODO: Explore other Advanced Abundance Estimation Options

                # Set rnaseq_cufflinks arguments.

                cufflinks.arguments.append(os.path.join(self.genome_directory, tophat.name, 'accepted_hits.bam'))

                # Create a new rnaseq_run_cufflinks Executable.

                pickler_dict_run_cufflinks = dict()
                pickler_dict_run_cufflinks['prefix'] = run_cufflinks_drms.name
                pickler_dict_run_cufflinks['replicate_key'] = replicate_key
                pickler_dict_run_cufflinks['cufflinks_executable'] = cufflinks

                pickler_path = os.path.join(
                    self.genome_directory,
                    '{}_{}.pkl'.format(run_cufflinks_drms.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_run_cufflinks)
                pickler_file.close()

                run_cufflinks = Executable.from_Analysis(
                    name=string.join(words=(run_cufflinks_drms.name, replicate_key), sep='_'),
                    program='bsf_run_rnaseq_cufflinks.py',
                    analysis=self)
                run_cufflinks_drms.add_Executable(run_cufflinks)
                run_cufflinks.dependencies.append(run_tophat.name)

                # Set rnaseq_run_cufflinks options.

                run_cufflinks.add_OptionLong(key='pickler_path', value=pickler_path)
                run_cufflinks.add_OptionLong(key='debug', value=str(self.debug))

                # Create a new rnaseq_process_cufflinks Executable.

                process_cufflinks = Executable.from_Analysis(
                    name=string.join(words=(process_cufflinks_drms.name, replicate_key), sep='_'),
                    program='bsf_rnaseq_process_cufflinks.R',
                    analysis=self)
                process_cufflinks_drms.add_Executable(process_cufflinks)
                process_cufflinks.dependencies.append(run_cufflinks.name)

                # Set rnaseq_process_cufflinks options.

                # TODO: The BioMart data set needs to be configured somewhere ...
                # It is specific for the genome assembly.
                # For the moment it is set in the rnaseq_config.ini file.

                # TODO: For the moment, the species-specific BioMart_data_set option needs specifying
                # in the configuration file.
                # process_cufflinks.add_OptionLong(key='biomart_data_set', value='mmusculus_gene_ensembl)
                process_cufflinks.add_OptionLong(key='sample',
                                                 value=replicate_key)
                process_cufflinks.add_OptionLong(key='genome_directory',
                                                 value=self.genome_directory)

    def _create_Cuffmerge_Cuffdiff_jobs(self):

        """Create Cuffdiff differential expression jobs.

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

        # Initialise the Distributed Resource Management System (DRMS) objects for
        # Cuffmerge and Cuffdiff Executable objects.

        cuffmerge_drms = DRMS.from_Analysis(
            name='rnaseq_cuffmerge',
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(cuffmerge_drms)

        run_cuffdiff_drms = DRMS.from_Analysis(
            name='rnaseq_run_cuffdiff',
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(run_cuffdiff_drms)

        process_cuffdiff_drms = DRMS.from_Analysis(
            name='rnaseq_process_cuffdiff',
            work_directory=self.genome_directory,
            analysis=self)
        self.drms_list.append(process_cuffdiff_drms)

        keys = self.comparisons.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            # Create a new rnaseq_cuffmerge Executable.

            cuffmerge = Cuffmerge(
                name=string.join(words=(cuffmerge_drms.name, key), sep='_'),
                analysis=self)
            cuffmerge_drms.add_Executable(cuffmerge)
            # DRMS.dependencies for rnaseq_cuffmerge are added for each rnaseq_run_cufflinks replicate below.

            # Set rnaseq_cuffmerge options.

            cuffmerge.add_OptionLong(key='output-dir',
                                     value=os.path.join(self.genome_directory, cuffmerge.name))
            cuffmerge.add_OptionLong(key='num-threads',
                                     value=str(cuffmerge_drms.threads))
            cuffmerge.add_OptionLong(key='ref-gtf',
                                     value=transcriptome)
            cuffmerge.add_OptionLong(key='ref-sequence',
                                     value=Default.absolute_genome_fasta(self.genome_version, 'bowtie2'))

            # Set rnaseq_cuffmerge arguments.

            # Create an assembly manifest file to merge all replicates of each Sample object ...

            assembly_name = string.join(words=(cuffmerge.name, 'assembly.txt'), sep='_')
            assembly_path = os.path.join(self.genome_directory, assembly_name)
            assembly_file = open(name=assembly_path, mode='w')

            # Create a new rnaseq_cuffdiff Executable.

            cuffdiff = Cuffdiff(
                name=string.join(words=('rnaseq_cuffdiff', key), sep='_'),
                analysis=self)
            # The rnaseq_cuffdiff Executable is run via the rnaseq_run_cuffdiff Executable below.

            # Set rnaseq_cuffdiff options.

            cuffdiff.add_OptionLong(key='output-dir',
                                    value=os.path.join(self.genome_directory, cuffdiff.name))
            cuffdiff.add_OptionLong(key='num-threads',
                                    value=str(run_cuffdiff_drms.threads))
            cuffdiff.add_OptionLong(key='frag-bias-correct',
                                    value=genome_fasta)

            # Process rnaseq_cuffmerge and rnaseq_cuffdiff arguments in parallel.

            # Keep a Python list of Python lists of TopHat aligned BAM files per comparison group.
            cuffdiff_alignments = list()
            cuffdiff_labels = list()

            for group_name, group_samples in self.comparisons[key]:

                cuffdiff_labels.append(group_name)
                alignments_list = list()
                cuffdiff_alignments.append(alignments_list)

                for sample in group_samples:

                    # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
                    # Python str key and Python list of Python list objects
                    # of Bio.BSF.Data.PairedReads objects.

                    replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

                    replicate_keys = replicate_dict.keys()
                    replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                    for replicate_key in replicate_keys:
                        # Add the Cufflinks assembled transcripts to the Cuffmerge manifest.
                        transcripts_path = os.path.join(
                            self.genome_directory,
                            string.join(words=('rnaseq_cufflinks', replicate_key), sep='_'),
                            'transcripts.gtf')
                        assembly_file.write(transcripts_path + '\n')

                        # Wait for each TopHat and Cufflinks replicate to finish,
                        # before Cuffmerge can run.

                        cuffmerge.dependencies.append(
                            string.join(words=('rnaseq_run_cufflinks', replicate_key), sep='_'))

                        # Add the TopHat accepted hits BAM file to Cuffdiff ...

                        alignments_list.append(os.path.join(
                            self.genome_directory,
                            string.join(words=('rnaseq_tophat', replicate_key), sep='_'),
                            'accepted_hits.bam'))

            assembly_file.close()

            # Add the assembly manifest file as Cuffmerge argument.
            cuffmerge.arguments.append(assembly_path)

            cuffdiff.add_OptionLong(key='labels', value=string.join(words=cuffdiff_labels, sep=','))

            # Add the Cuffmerge merged assembly as Cuffdiff output.
            cuffdiff.arguments.append(os.path.join(self.genome_directory, cuffmerge.name, 'merged.gtf'))

            # Add the TopHat aligned BAM files per point as Cuffdiff arguments.
            for alignments_list in cuffdiff_alignments:
                cuffdiff.arguments.append(string.join(words=alignments_list, sep=','))

            # Create a new rnaseq_run_cuffdiff Executable.

            pickler_dict_run_cuffdiff = dict()
            pickler_dict_run_cuffdiff['prefix'] = run_cuffdiff_drms.name
            pickler_dict_run_cuffdiff['replicate_key'] = key
            pickler_dict_run_cuffdiff['cuffdiff_executable'] = cuffdiff

            pickler_path = os.path.join(
                self.genome_directory,
                '{}_{}.pkl'.format(run_cuffdiff_drms.name, key))
            pickler_file = open(pickler_path, 'wb')
            pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
            pickler.dump(obj=pickler_dict_run_cuffdiff)
            pickler_file.close()

            run_cuffdiff = Executable.from_Analysis(
                name=string.join(words=(run_cuffdiff_drms.name, key), sep='_'),
                program='bsf_run_rnaseq_cuffdiff.py',
                analysis=self)
            run_cuffdiff_drms.add_Executable(run_cuffdiff)
            run_cuffdiff.dependencies.append(cuffmerge.name)

            # Set rnaseq_run_cuffdiff options.

            run_cuffdiff.add_OptionLong(key='pickler_path', value=pickler_path)
            run_cuffdiff.add_OptionLong(key='debug', value=str(self.debug))

            # Create a new rnaseq_process_cuffdiff Executable.

            process_cuffdiff = Executable.from_Analysis(
                name=string.join(words=(process_cuffdiff_drms.name, key), sep='_'),
                program='bsf_rnaseq_process_cuffdiff.R',
                analysis=self)
            process_cuffdiff_drms.add_Executable(process_cuffdiff)
            process_cuffdiff.dependencies.append(run_cuffdiff.name)

            # Set rnaseq_process_cuffdiff options.

            process_cuffdiff.add_OptionLong(key='comparison',
                                            value=key)
            process_cuffdiff.add_OptionLong(key='genome_directory',
                                            value=self.genome_directory)

            # Set rnaseq_process_cuffdiff arguments.

            # None so far.

    def report(self):

        """Create a Tuxedo report in HTML format and a UCSC Genome Browser Track Hub.

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
        output += 'exon-exon splice junctions. It is built on the ultra fast\n'
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
        output += '<th>Aligned BAM file</th>\n'
        output += '<th>Aligned BAI file</th>\n'
        output += '</tr>\n'
        output += '\n'

        # Group via UCSC super tracks.

        track_output += 'track Alignments\n'
        track_output += 'shortLabel Alignments\n'
        track_output += 'longLabel TopHat RNA-Seq read alignments\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Coverage\n'
        track_output += 'shortLabel Coverage\n'
        track_output += 'longLabel TopHat RNA-Seq alignment coverage\n'
        track_output += 'visibility full\n'
        track_output += 'superTrack on\n'
        track_output += 'group coverage\n'
        track_output += '\n'

        track_output += 'track Deletions\n'
        track_output += 'shortLabel Deletions\n'
        track_output += 'longLabel TopHat RNA-Seq deletions\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Insertions\n'
        track_output += 'shortLabel Insertions\n'
        track_output += 'longLabel TopHat RNA-Seq insertions\n'
        track_output += 'visibility hide\n'
        track_output += 'superTrack on\n'
        track_output += 'group alignments\n'
        track_output += '\n'

        track_output += 'track Junctions\n'
        track_output += 'shortLabel Junctions\n'
        track_output += 'longLabel TopHat RNA-Seq splice junctions\n'
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

                track_output += 'track {}_alignments\n'. \
                    format(replicate_key)
                track_output += 'type bam\n'
                track_output += 'shortLabel {}_alignments\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} TopHat RNA-Seq read alignments\n'. \
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

                track_output += 'track {}_coverage\n'. \
                    format(replicate_key)
                # TODO: The bigWig type must declare the expected signal range.
                # The signal range of a bigWig file would be available via the UCSC tool bigWigInfo.
                track_output += 'type bigWig\n'
                track_output += 'shortLabel {}_coverage\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} TopHat RNA-Seq alignment coverage\n'. \
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

                track_output += 'track {}_deletions\n'. \
                    format(replicate_key)
                track_output += 'type bigBed\n'
                track_output += 'shortLabel {}_deletions\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} TopHat RNA-Seq deletions\n'. \
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
                track_output += 'shortLabel {}_insertions\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} TopHat RNA-Seq insertions\n'. \
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

                track_output += 'track {}_junctions\n'. \
                    format(replicate_key)
                track_output += 'type bigBed\n'
                track_output += 'shortLabel {}_junctions\n'. \
                    format(replicate_key)
                track_output += 'longLabel {} TopHat RNA-Seq splice junctions\n'. \
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
                track_dict['name'] = '{}_transcripts'. \
                    format(replicate_key)
                track_dict['description'] = '"{} Cufflinks RNA-Seq transcript assembly"'. \
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
                output += '<td><a href="./{}/rnaseq_tophat_{}_accepted_hits.bam">Aligned BAM</td>\n'. \
                    format(prefix, replicate_key)
                output += '<td><a href="./{}/rnaseq_tophat_{}_accepted_hits.bam.bai">Aligned BAI</td>\n'. \
                    format(prefix, replicate_key)
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
        output += '<th>Replicates</th>\n'
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

            # Link to comparison-specific symbolic links in the directory after cummeRbund processing.

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
            output += '<td><a href="./{}/{}_replicates.txt">Replicates</a></td>\n'. \
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
