"""bsf.analyses.rna_seq

A package of classes and methods supporting RNA-Seq analyses.
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

import os.path
from pickle import Pickler, HIGHEST_PROTOCOL
import re

from bsf import Analysis, defaults, DRMS, Runnable
from bsf.annotation import AnnotationSheet, TuxedoSamplePairSheet
from bsf.data import PairedReads, Sample
from bsf.executables import TopHat
from bsf.process import Executable, RunnableStep
from bsf.standards import Configuration, Default


class Tuxedo(Analysis):
    """Tuxedo RNASeq Analysis sub-class.

    Attributes:
    @cvar drms_name_run_tophat: C{DRMS.name} for the run Tophat C{Analysis} stage
    @type drms_name_run_tophat: str
    @cvar drms_name_process_tophat: C{DRMS.name} for the process Tophat C{Analysis} stage
    @type drms_name_process_tophat: str
    @cvar drms_name_run_cufflinks: C{DRMS.name} for the run Cufflinks C{Analysis} stage
    @type drms_name_run_cufflinks: str
    @cvar drms_name_process_cufflinks: C{DRMS.name} for the process Cufflinks C{Analysis} stage
    @type drms_name_process_cufflinks: str
    @cvar drms_name_run_cuffmerge: C{DRMS.name} for the Cuffmerge C{Analysis} stage
    @type drms_name_run_cuffmerge: str
    @cvar drms_name_run_cuffquant: C{DRMS.name} for the Cuffquant C{Analysis} stage
    @type drms_name_run_cuffquant: str
    @cvar drms_name_run_cuffnorm: C{DRMS.name} for the Cuffnorm C{Analysis} stage
    @type drms_name_run_cuffnorm: str
    @cvar drms_name_run_cuffdiff: C{DRMS.name} for the Cuffdiff C{Analysis} stage
    @type drms_name_run_cuffdiff: str
    @cvar drms_name_process_cuffdiff: C{DRMS.name} for the process Cuffdiff C{Analysis} stage
    @type drms_name_process_cuffdiff: str
    @ivar replicate_grouping: Group all replicates into a single Tophat and Cufflinks process
    @type replicate_grouping: bool
    @ivar cmp_file: Comparison file
    @type cmp_file: str | unicode
    @ivar genome_fasta_path: Reference genome sequence FASTA file path
    @type genome_fasta_path: str | unicode
    @ivar transcriptome_gtf_path: Reference transcriptome GTF file path
    @type transcriptome_gtf_path: str | unicode
    @ivar transcriptome_index_path: Tophat transcriptome index path
    @type transcriptome_index_path: str | unicode
    @ivar mask_gtf_path: GTF file path to mask transcripts
    @type mask_gtf_path: str | unicode
    @ivar multi_read_correction: Apply multi-read correction
    @type multi_read_correction: bool
    @ivar library_type: Library type
        Cuffquant and Cuffdiff I{fr-unstranded} (default), I{fr-firststrand} or I{fr-secondstrand}
    @type library_type: str
    @ivar no_length_correction: Do not correct for transcript lengths as in 3-prime sequencing
    @type no_length_correction: bool
    """

    # Replicate stage
    drms_name_run_tophat = 'rnaseq_run_tophat'
    drms_name_process_tophat = 'rnaseq_process_tophat'
    drms_name_run_cufflinks = 'rnaseq_run_cufflinks'
    drms_name_process_cufflinks = 'rnaseq_process_cufflinks'

    # Comparison stage
    drms_name_run_cuffmerge = 'rnaseq_cuffmerge'
    drms_name_run_cuffquant = 'rnaseq_cuffquant'
    drms_name_run_cuffnorm = 'rnaseq_cuffnorm'
    drms_name_run_cuffdiff = 'rnaseq_cuffdiff'
    drms_name_process_cuffdiff = 'rnaseq_process_cuffdiff'

    @classmethod
    def get_prefix_rnaseq_run_tophat(cls, replicate_key):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param replicate_key: Replicate key
        @type replicate_key: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return '_'.join((cls.drms_name_run_tophat, replicate_key))

    @classmethod
    def get_prefix_rnaseq_run_cuffmerge(cls, comparison_key):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param comparison_key: Comparison key
        @type comparison_key: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return '_'.join((cls.drms_name_run_cuffmerge, comparison_key))

    @classmethod
    def get_prefix_rnaseq_run_cuffquant(cls, comparison_key):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param comparison_key: Comparison key
        @type comparison_key: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return '_'.join((cls.drms_name_run_cuffquant, comparison_key))

    @classmethod
    def get_prefix_rnaseq_run_cuffnorm(cls, comparison_key):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param comparison_key: Comparison key
        @type comparison_key: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return '_'.join((cls.drms_name_run_cuffnorm, comparison_key))

    @classmethod
    def get_prefix_rnaseq_run_cuffdiff(cls, comparison_key):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param comparison_key: Comparison key
        @type comparison_key: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return '_'.join((cls.drms_name_run_cuffdiff, comparison_key))

    @classmethod
    def get_prefix_rnaseq_process_cuffdiff(cls, comparison_key):
        """Get a Python C{str} for setting C{Executable.dependencies} in other C{Analysis} objects

        @param comparison_key: Comparison key
        @type comparison_key: str
        @return: The dependency string for an C{Executable} of this C{Analysis}
        @rtype: str
        """
        return '_'.join((cls.drms_name_process_cuffdiff, comparison_key))

    def __init__(self, configuration=None, project_name=None, genome_version=None, input_directory=None,
                 output_directory=None, project_directory=None, genome_directory=None, e_mail=None, debug=0,
                 drms_list=None, collection=None, comparisons=None, samples=None, replicate_grouping=False,
                 cmp_file=None, genome_fasta_path=None, transcriptome_gtf_path=None, transcriptome_index_path=None,
                 mask_gtf_path=None, multi_read_correction=None, library_type=None, no_length_correction=False):
        """Initialise a C{Tuxedo} object.

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
        @param comparisons: Python C{dict} of Python C{str} (comparison name) key objects and
            Python C{tuple} value objects of C{Sample.name} and Python C{list} of C{Sample} objects
        @type comparisons: dict[str, (Sample.name, list[Sample])]
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
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
                drms_list=drms_list,
                collection=collection,
                comparisons=comparisons,
                samples=samples)

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
        """Set instance variables of a C{Tuxedo} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
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

        option = 'genome_fasta'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_fasta_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_gtf'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_gtf_path = configuration.config_parser.get(section=section, option=option)

        option = 'transcriptome_index'
        if configuration.config_parser.has_option(section=section, option=option):
            self.transcriptome_index_path = configuration.config_parser.get(section=section, option=option)

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

    def _calculate_comparisons_from_samples(self):
        """Calculate the comparisons on the basis of the sample annotation sheet alone.
        An all-against-all sample comparison gets configured.
        @return:
        @rtype:
        """

        # Without a comparison file path, simply add all Sample objects from the Collection.
        self.samples.extend(self.collection.get_all_samples())

        # Create a global comparison by adding all samples under their sample name as group name.

        comparison_groups = list()
        for sample in self.samples:
            assert isinstance(sample, Sample)
            # Add a tuple of group (i.e. sample) name and a Python list of the Sample object.
            comparison_groups.append((sample.name, [sample]))

        # Sort the list of comparison groups by Sample.name.
        comparison_groups.sort(cmp=lambda x, y: cmp(x[0], y[0]))

        self.comparisons['global'] = comparison_groups

        return

    def _calculate_comparisons_from_groups(self):
        """Calculate the comparisons on the basis of the sample annotation sheet alone.
        An all-against-all group comparison gets configured.

        @return:
        @rtype:
        """

        # Without a comparison file path, simply add all Sample objects from the Collection.
        self.samples.extend(self.collection.get_all_samples())

        # Create a global comparison by adding all sample groups.

        comparison_groups = list()
        for key in self.collection.sample_groups.keys():
            comparison_groups.append((key, self.collection.sample_groups[key]))

        # Sort the list of comparison groups by group name.
        comparison_groups.sort(cmp=lambda x, y: cmp(x[0], y[0]))

        self.comparisons['global'] = comparison_groups

        return

    def _read_comparisons(self, cmp_file):
        """Read an C{AnnotationSheet} CSV file specifying comparisons from disk.

        All C{Sample} objects referenced in a comparison are added from the C{Collection} to the C{Analysis} object.

            - Column headers for CASAVA folders:
                - Treatment/Control/Point N ProcessedRunFolder:
                    - CASAVA processed run folder name or
                    - C{Analysis.input_directory} by default
                - Treatment/Control/Point N Project:
                    - CASAVA Project name or
                    - C{Analysis.project_name} by default
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
            assert isinstance(row_dict, dict)
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
                assert isinstance(group_name, str)
                assert isinstance(group_samples, list)
                if group_name and len(group_samples):
                    key += group_name
                    key += '__'
                    comparison_groups.append((group_name, group_samples))
                    # Also expand each Python list of Sample objects to get all those Sample objects
                    # that this Analysis needs considering.
                    for sample in group_samples:
                        assert isinstance(sample, Sample)
                        if self.debug > 1:
                            print '  {} Sample name: {!r} file_path:{!r}'.format(prefix, sample.name, sample.file_path)
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

            if 'Comparison Name' in row_dict:
                # For ridiculously large comparisons involving loads of groups or samples a comparison name can be
                # explicitly specified. Any non-word characters get replaced by underscore characters.
                key = str(row_dict['Comparison Name'])
                key = re.sub(pattern=regular_expression, repl='_', string=key)
            else:
                # Truncate the last '__' separator off the key string.
                key = key.rstrip('_')

            self.comparisons[key] = comparison_groups

        return

    def run(self):
        """Run this C{Tuxedo} analysis.
        @return:
        @rtype:
        """

        super(Tuxedo, self).run()

        # Tuxedo requires a genome version.

        if not self.genome_version:
            raise Exception('A Tuxedo analysis requires a genome_version configuration option.')

        if self.cmp_file:
            # A comparison file path has been provided.
            if self.cmp_file == '*groups*':
                # The special file name *groups* creates the comparisons on the basis of an
                # all-against-all group comparison.
                self._calculate_comparisons_from_groups()
            elif self.cmp_file == '*samples*':
                # The special file name *samples* creates the comparisons on the basis of an
                # all-against-all sample comparison.
                self._calculate_comparisons_from_samples()
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
            self.samples.extend(self.collection.get_all_samples())

        # Experimentally, sort the Python list of Sample objects by the Sample name.
        # This cannot be done in the super-class, because Sample objects are only put into the
        # Analysis.samples list by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Define the reference genome FASTA file path.
        # If it does not exist, construct it from defaults.

        if not self.genome_fasta_path:
            self.genome_fasta_path = Default.absolute_genome_fasta(
                    genome_version=self.genome_version,
                    genome_index='bowtie2')

        # Define the reference transcriptome GTF file path.
        # Check if transcriptome_gtf is an absolute path and
        # prepend the annotation default if not.

        if self.transcriptome_index_path:
            self.transcriptome_index_path = Default.get_absolute_path(
                    file_path=self.transcriptome_index_path,
                    default_path=Default.absolute_genomes(genome_version=self.genome_version))

            if not os.path.isdir(self.transcriptome_index_path):
                raise Exception("Reference transcriptome index directory {!r} does not exist.".
                                format(self.transcriptome_index_path))

            transcriptome_index = os.path.basename(self.transcriptome_index_path)

            # The tophat --transcript-index process puts a GFF file into the index directory
            # that really is a GTF file. A symbolic link to a GTF file is needed to make the
            # process cuffdiff script work.

            self.transcriptome_gtf_path = os.path.join(
                    self.transcriptome_index_path,
                    '.'.join((transcriptome_index, 'gtf')))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        elif self.transcriptome_gtf_path:
            self.transcriptome_gtf_path = Default.get_absolute_path(
                    file_path=self.transcriptome_gtf_path,
                    default_path=Default.absolute_genome_annotation(genome_version=self.genome_version))

            if not os.path.exists(self.transcriptome_gtf_path):
                raise Exception("Reference transcriptome GTF file {!r} does not exist.".
                                format(self.transcriptome_gtf_path))

        if not self.transcriptome_gtf_path:
            raise Exception("Reference transcriptome GTF file not defined.")

        self._create_tophat_cufflinks_jobs()
        self._create_cuffmerge_cuffdiff_jobs()

        return

    def _create_tophat_cufflinks_jobs(self):
        """Create a TopHat aligner process and a Cufflinks transcript assembler process
        for each sample or replicate.
        @return:
        @rtype:
        """

        if self.debug > 1:
            print '{!r} method _create_TopHat_Cufflinks_jobs:'.format(self)

        # Read configuration options.

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        # TODO: Move the ConfigParser code.
        genome_sizes = config_parser.get(section=config_section, option='genome_sizes')
        genome_sizes = os.path.expanduser(genome_sizes)
        genome_sizes = os.path.expandvars(genome_sizes)

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

        # Initialise the Distributed Resource Management System (DRMS) objects for
        # TopHat and Cufflinks Executable objects.

        drms_run_tophat = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_run_tophat,
                working_directory=self.genome_directory,
                analysis=self))

        drms_process_tophat = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_process_tophat,
                working_directory=self.genome_directory,
                analysis=self))

        drms_run_cufflinks = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_run_cufflinks,
                working_directory=self.genome_directory,
                analysis=self))

        drms_process_cufflinks = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_process_cufflinks,
                working_directory=self.genome_directory,
                analysis=self))

        process_cufflinks_dependencies = list()

        for sample in self.samples:
            assert isinstance(sample, Sample)
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # bsf.data.Sample.get_all_paired_reads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of bsf.data.PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                assert isinstance(replicate_key, str)

                # Create a Tophat Runnable per replicate (replicate_key).

                # TODO: Activate the new code once the bsf_run_rnaseq_tophat.py script has been retired.

                # prefix_run_tophat = '_'.join((drms_run_tophat.name, replicate_key))

                file_path_dict_tophat = dict(
                        # The output directory deviates from the prefix_run_tophat that itself is based on
                        # drms_run_tophat.name. Both rnaseq_run_tophat and rnaseq_process_tophat processes should
                        # use the same rnaseq_tophat directory.
                        # output_directory=prefix_run_tophat,
                        output_directory='_'.join(('rnaseq_tophat', replicate_key))
                )

                # runnable_run_tophat = self.add_runnable(
                #         runnable=Runnable(
                #                 name=self.get_prefix_rnaseq_run_tophat(replicate_key=replicate_key),
                #                 code_module='bsf.runnables.generic',
                #                 working_directory=self.genome_directory,
                #                 file_path_dict=file_path_dict_tophat,
                #                 debug=self.debug))
                # executable_run_tophat = self.set_drms_runnable(
                #         drms=drms_run_tophat,
                #         runnable=runnable_run_tophat)
                #
                # Create a new Tophat RunnableStep.

                # tophat = runnable_run_tophat.add_runnable_step(
                #     runnable_step=RunnableStep(
                #     name='tophat2',
                #     program='tophat2'))

                tophat = TopHat(
                        name='_'.join(('rnaseq_tophat', replicate_key)),
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
                        value=file_path_dict_tophat['output_directory'])
                tophat.add_option_long(
                        key='num-threads',
                        value=str(drms_run_tophat.threads))
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

                tophat.arguments.append(bowtie2_index)

                # Set rnaseq_tophat arguments for reads1 and reads2.

                reads1 = list()
                reads2 = list()

                for paired_reads in replicate_dict[replicate_key]:
                    assert isinstance(paired_reads, PairedReads)
                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    if paired_reads.reads1:
                        reads1.append(paired_reads.reads1.file_path)
                    if paired_reads.reads2:
                        reads2.append(paired_reads.reads2.file_path)

                tophat.arguments.append(','.join(reads1))
                tophat.arguments.append(','.join(reads2))

                # Create a new rnaseq_run_tophat Executable.
                # TODO: The following code block is required as long as the bsf_run_rnaseq_tophat.py script
                # has not been retired.

                pickler_dict_run_tophat = dict()
                pickler_dict_run_tophat['prefix'] = drms_run_tophat.name
                pickler_dict_run_tophat['replicate_key'] = replicate_key
                pickler_dict_run_tophat['tophat_executable'] = tophat

                pickler_path = os.path.join(
                        self.genome_directory,
                        '{}_{}.pkl'.format(drms_run_tophat.name, replicate_key))
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_run_tophat)
                pickler_file.close()

                run_tophat = drms_run_tophat.add_executable(
                        executable=Executable(
                                name='_'.join((drms_run_tophat.name, replicate_key)),
                                program='bsf_run_rnaseq_tophat.py'))

                # The rnaseq_run_tophat Executable requires no DRMS.dependencies.

                # Set rnaseq_run_tophat options.
                self.set_command_configuration(command=run_tophat)
                run_tophat.add_option_long(key='pickler_path', value=pickler_path)
                run_tophat.add_option_long(key='debug', value=str(self.debug))

                # TODO: End of code block.

                # Create a new rnaseq_process_tophat Executable.

                process_tophat = drms_process_tophat.add_executable(
                        executable=Executable(
                                name='_'.join((drms_process_tophat.name, replicate_key)),
                                program='bsf_rnaseq_process_tophat2.sh'))
                process_tophat.dependencies.append(run_tophat.name)

                # Set rnaseq_process_tophat options.

                # Set rnaseq_process_tophat arguments.
                self.set_command_configuration(command=process_tophat)
                process_tophat.arguments.append(file_path_dict_tophat['output_directory'])
                process_tophat.arguments.append(genome_sizes)

                # TODO: Switch from an external Bash script to a set of Runnable and RunnableStep objects.
                # Since the Bash script includes Perl code to reset the BED score field to 0, rather than
                # re-scale it properly, it would be good to write a new bsf.runnables.process_tophat module
                # to implement this in Python code.

                # Create a process_tophat Runnable.

                # file_path_dict_process_tophat = dict(
                #
                # )
                #
                # runnable_process_tophat = self.add_runnable(runnable=Runnable(
                #     name='rnaseq_process_tophat',
                #     code_module='bsf.runnables.generic',
                #     working_directory=self.genome_directory,
                #     file_path_dict=file_path_dict_process_tophat,
                #     debug=self.debug))

                # Create a Cufflinks Runnable per replicate (replicate_key).

                prefix_run_cufflinks = '_'.join((drms_run_cufflinks.name, replicate_key))

                file_path_dict_cufflinks = dict(
                        # The output directory deviates from the prefix_run_cufflinks that itself is based on
                        # drms_run_cufflinks.name. Both rnaseq_run_cufflinks and rnaseq_process_cufflinks processes
                        # should use the same rnaseq_cufflinks directory.
                        # output_directory=prefix_run_cufflinks,
                        output_directory='_'.join(('rnaseq_cufflinks', replicate_key)),
                        tophat_accepted_hits=os.path.join(file_path_dict_tophat['output_directory'],
                                                          'accepted_hits.bam')
                )

                runnable_run_cufflinks = self.add_runnable(
                        runnable=Runnable(
                                name=prefix_run_cufflinks,
                                code_module='bsf.runnables.generic',
                                working_directory=self.genome_directory,
                                file_path_dict=file_path_dict_cufflinks,
                                debug=self.debug))

                # Create an Executable for running the Cufflinks Runnable.

                executable_run_cufflinks = self.set_drms_runnable(
                        drms=drms_run_cufflinks,
                        runnable=runnable_run_cufflinks)
                executable_run_cufflinks.dependencies.append(run_tophat.name)

                # Create a new Cufflinks RunnableStep.

                cufflinks = runnable_run_cufflinks.add_runnable_step(
                        runnable_step=RunnableStep(
                                name='cufflinks',
                                program='cufflinks'))

                # Create a new rnaseq_cufflinks Executable, which is run via the rnaseq_run_cufflinks Executable below.

                # Set Cufflinks options.

                # Cufflinks has a GTF option, in which case it will not assemble
                # novel transcripts and a GTF-guide option in which case it will
                # assemble novel transcripts.

                if novel_transcripts:
                    cufflinks.add_option_long(
                            key='GTF-guide',
                            value=self.transcriptome_gtf_path)
                else:
                    cufflinks.add_option_long(
                            key='GTF',
                            value=self.transcriptome_gtf_path)

                if self.mask_gtf_path:
                    cufflinks.add_option_long(
                            key='mask-file',
                            value=self.mask_gtf_path)

                cufflinks.add_option_long(
                        key='frag-bias-correct',
                        value=self.genome_fasta_path)

                if self.multi_read_correction:
                    cufflinks.add_switch_long(
                            key='multi-read-correct')

                if self.library_type:
                    cufflinks.add_option_long(
                            key='library-type',
                            value=self.library_type)

                # Cufflinks has a --library-norm-method option, but only one option (classic-fpkm) seems supported.

                if self.no_length_correction:
                    cufflinks.add_switch_long(
                            key='no-length-correction')

                cufflinks.add_option_long(
                        key='output-dir',
                        value=file_path_dict_cufflinks['output_directory'])

                cufflinks.add_option_long(
                        key='num-threads',
                        value=str(drms_run_cufflinks.threads))

                cufflinks.add_switch_long(
                        key='quiet')

                cufflinks.add_switch_long(
                        key='no-update-check')

                # Set Cufflinks arguments.

                cufflinks.arguments.append(file_path_dict_cufflinks['tophat_accepted_hits'])

                # Add the run_cufflinks dependency for process_cufflinks.
                process_cufflinks_dependencies.append(executable_run_cufflinks.name)

        # Create one process_cufflinks Executable to process all sub-directories.

        process_cufflinks = drms_process_cufflinks.add_executable(
                executable=Executable(
                        name=drms_process_cufflinks.name,
                        program='bsf_rnaseq_process_cufflinks.R'))
        process_cufflinks.dependencies.extend(process_cufflinks_dependencies)

        # Set process_cufflinks options.
        self.set_command_configuration(command=process_cufflinks)
        # None so far.

        return

    def _create_cuffmerge_cuffdiff_jobs(self):
        """Create a Cuffmerge and a Cuffdiff process for each comparison.
        @return:
        @rtype:
        """

        run_cuffquant_before_cuffdiff = False

        # Initialise the Distributed Resource Management System (DRMS) objects for
        # Cuffmerge and Cuffdiff Executable objects.

        drms_run_cuffmerge = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_run_cuffmerge,
                working_directory=self.genome_directory,
                analysis=self))

        drms_run_cuffquant = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_run_cuffquant,
                working_directory=self.genome_directory,
                analysis=self))

        drms_run_cuffnorm = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_run_cuffnorm,
                working_directory=self.genome_directory,
                analysis=self))

        drms_run_cuffdiff = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_run_cuffdiff,
                working_directory=self.genome_directory,
                analysis=self))

        drms_process_cuffdiff = self.add_drms(drms=DRMS.from_analysis(
                name=self.drms_name_process_cuffdiff,
                working_directory=self.genome_directory,
                analysis=self))

        comparison_keys = self.comparisons.keys()
        comparison_keys.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_key in comparison_keys:
            assert isinstance(comparison_key, str)

            # TODO: Should the comparison prefix also include the project name or number?
            prefix_cuffmerge = self.get_prefix_rnaseq_run_cuffmerge(comparison_key=comparison_key)

            file_path_dict_cuffmerge = dict(
                    output_directory=prefix_cuffmerge,
                    merged_gtf=os.path.join(prefix_cuffmerge, 'merged.gtf'),
                    assembly_txt='_'.join((prefix_cuffmerge, 'assembly.txt')))

            runnable_run_cuffmerge = self.add_runnable(
                    runnable=Runnable(
                            name=prefix_cuffmerge,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            file_path_dict=file_path_dict_cuffmerge,
                            debug=self.debug))
            executable_run_cuffmerge = self.set_drms_runnable(
                    drms=drms_run_cuffmerge,
                    runnable=runnable_run_cuffmerge)

            # Create a new Cuffmerge RunnableStep.

            cuffmerge = runnable_run_cuffmerge.add_runnable_step(
                    runnable_step=RunnableStep(
                            name='cuffmerge',
                            program='cuffmerge'))

            # Set rnaseq_cuffmerge options.

            cuffmerge.add_option_long(
                    key='output-dir',
                    value=file_path_dict_cuffmerge['output_directory'])
            cuffmerge.add_option_long(
                    key='num-threads',
                    value=str(drms_run_cuffmerge.threads))
            cuffmerge.add_option_long(
                    key='ref-gtf',
                    value=self.transcriptome_gtf_path)
            cuffmerge.add_option_long(
                    key='ref-sequence',
                    value=self.genome_fasta_path)

            # Set rnaseq_cuffmerge arguments.

            # Create an assembly manifest file to merge all replicates of each Sample object.
            # This file requires an absolute path, because the working directory is not set at the stage of
            # job submission.

            assembly_path = os.path.join(self.genome_directory, '_'.join((prefix_cuffmerge, 'assembly.txt')))
            assembly_file = open(name=assembly_path, mode='w')

            # Process rnaseq_cuffmerge and rnaseq_cuffdiff arguments in parallel.

            # Create a Python list of Python list objects of Cuffquant abundances per comparison group.
            cuffdiff_cuffnorm_abundances = list()
            cuffdiff_cuffnorm_alignments = list()
            cuffdiff_cuffnorm_dependencies = list()
            cuffdiff_cuffnorm_labels = list()

            for group_name, group_samples in self.comparisons[comparison_key]:
                assert isinstance(group_name, str)
                assert isinstance(group_samples, list)
                cuffdiff_cuffnorm_labels.append(group_name)
                per_group_abundances_list = list()
                cuffdiff_cuffnorm_abundances.append(per_group_abundances_list)
                per_group_alignments_list = list()
                cuffdiff_cuffnorm_alignments.append(per_group_alignments_list)

                for sample in group_samples:
                    assert isinstance(sample, Sample)
                    # bsf.data.Sample.get_all_paired_reads returns a Python dict of
                    # Python str comparison_key and Python list of Python list objects
                    # of bsf.data.PairedReads objects.

                    replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

                    replicate_keys = replicate_dict.keys()
                    replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

                    for replicate_key in replicate_keys:
                        assert isinstance(comparison_key, str)

                        # Add the Cufflinks assembled transcripts to the Cuffmerge manifest.

                        transcripts_path = os.path.join(
                                self.genome_directory,
                                '_'.join(('rnaseq_cufflinks', replicate_key)),
                                'transcripts.gtf')
                        assembly_file.write(transcripts_path + '\n')

                        # Wait for each Cufflinks replicate to finish, before Cuffmerge can run.

                        executable_run_cuffmerge.dependencies.append(
                                '_'.join((self.drms_name_run_cufflinks, replicate_key)))

                        # Create a Cuffquant Runnable per comparison (comparison_key) and replicate (replicate_key)
                        # on the basis of the Cuffmerge GTF file.

                        prefix_cuffquant = '_'.join((drms_run_cuffquant.name, comparison_key, replicate_key))

                        file_path_dict_cuffquant = dict(
                                output_directory=prefix_cuffquant,
                                abundances=os.path.join(prefix_cuffquant, 'abundances.cxb'),
                                merged_gtf=file_path_dict_cuffmerge['merged_gtf'],
                                tophat_directory='_'.join(('rnaseq_tophat', replicate_key)),
                                tophat_accepted_hits=os.path.join(
                                        '_'.join(('rnaseq_tophat', replicate_key)),
                                        'accepted_hits.bam'))

                        runnable_run_cuffquant = self.add_runnable(
                                runnable=Runnable(
                                        name=prefix_cuffquant,
                                        code_module='bsf.runnables.generic',
                                        working_directory=self.genome_directory,
                                        file_path_dict=file_path_dict_cuffquant,
                                        debug=self.debug))
                        executable_run_cuffquant = self.set_drms_runnable(
                                drms=drms_run_cuffquant,
                                runnable=runnable_run_cuffquant)

                        # Each Cuffquant process depends on Cuffmerge.

                        executable_run_cuffquant.dependencies.append(executable_run_cuffmerge.name)

                        # Create a new cuffquant RunnableStep.

                        cuffquant = runnable_run_cuffquant.add_runnable_step(
                                runnable_step=RunnableStep(
                                        name='cuffquant',
                                        program='cuffquant'))

                        # Set Cuffquant options.

                        cuffquant.add_option_long(
                                key='output-dir',
                                value=file_path_dict_cuffquant['output_directory'])
                        cuffquant.add_option_long(
                                key='num-threads',
                                value=str(drms_run_cuffquant.threads))
                        if self.mask_gtf_path:
                            cuffquant.add_option_long(
                                    key='mask-file',
                                    value=self.mask_gtf_path)
                        cuffquant.add_option_long(
                                key='frag-bias-correct',
                                value=self.genome_fasta_path)
                        if self.multi_read_correction:
                            cuffquant.add_switch_long(
                                    key='multi-read-correct')
                        if self.library_type:
                            cuffquant.add_option_long(
                                    key='library-type',
                                    value=self.library_type)
                        if self.no_length_correction:
                            cuffquant.add_switch_long(
                                    key='no-length-correction')
                        cuffquant.add_switch_long(
                                key='quiet')
                        cuffquant.add_switch_long(
                                key='no-update-check')

                        # Set Cuffquant arguments.
                        # Add the Cuffmerge GTF file and the TopHat BAM file as Cuffquant arguments.

                        cuffquant.arguments.append(file_path_dict_cuffquant['merged_gtf'])
                        cuffquant.arguments.append(file_path_dict_cuffquant['tophat_accepted_hits'])

                        # Add the Cuffquant abundances file to the Cuffdiff list.

                        per_group_abundances_list.append(file_path_dict_cuffquant['abundances'])

                        # Add the TopHat BAM file to the Cuffdiff alignments list.

                        per_group_alignments_list.append(file_path_dict_cuffquant['tophat_accepted_hits'])

                        # Add the Cuffquant Runnable process name to the Cuffdiff and Cuffnorm dependencies list.

                        cuffdiff_cuffnorm_dependencies.append(executable_run_cuffquant.name)

            assembly_file.close()

            # Add the assembly manifest file as Cuffmerge argument.
            cuffmerge.arguments.append(assembly_path)

            # Create a Cuffnorm Runnable per comparison.

            prefix_cuffnorm = '_'.join((drms_run_cuffnorm.name, comparison_key))

            file_path_dict_cuffnorm = dict(
                    output_directory=prefix_cuffnorm,
                    merged_gtf=file_path_dict_cuffmerge['merged_gtf'])

            runnable_run_cuffnorm = self.add_runnable(
                    runnable=Runnable(
                            name=prefix_cuffnorm,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            file_path_dict=file_path_dict_cuffnorm,
                            debug=self.debug))
            executable_run_cuffnorm = self.set_drms_runnable(
                    drms=drms_run_cuffnorm,
                    runnable=runnable_run_cuffnorm)

            executable_run_cuffnorm.dependencies.extend(cuffdiff_cuffnorm_dependencies)

            # Create a new Cuffnorm RunnableStep.

            cuffnorm = runnable_run_cuffnorm.add_runnable_step(
                    runnable_step=RunnableStep(
                            name='cuffnorm',
                            program='cuffnorm'))

            # Set Cuffnorm options.

            cuffnorm.add_option_long(
                    key='output-dir',
                    value=file_path_dict_cuffnorm['output_directory'])
            cuffnorm.add_option_long(
                    key='labels',
                    value=','.join(cuffdiff_cuffnorm_labels))
            cuffnorm.add_option_long(
                    key='num-threads',
                    value=str(drms_run_cuffnorm.threads))
            if self.library_type:
                cuffnorm.add_option_long(
                        key='library-type',
                        value=self.library_type)
            cuffnorm.add_switch_long(
                    key='quiet')
            cuffnorm.add_switch_long(
                    key='no-update-check')

            # Add the Cuffmerge GTF file as first Cuffnorm argument.
            cuffnorm.arguments.append(file_path_dict_cuffnorm['merged_gtf'])

            # Add the Cuffquant abundances files per point as Cuffnorm arguments.
            for per_group_abundances_list in cuffdiff_cuffnorm_abundances:
                assert isinstance(per_group_abundances_list, list)
                cuffnorm.arguments.append(','.join(per_group_abundances_list))

            # Create a Cuffdiff Runnable per comparison.

            prefix_cuffdiff = '_'.join((drms_run_cuffdiff.name, comparison_key))

            file_path_dict_cuffdiff = dict(
                    output_directory=prefix_cuffdiff,
                    merged_gtf=file_path_dict_cuffmerge['merged_gtf'])

            runnable_run_cuffdiff = self.add_runnable(
                    runnable=Runnable(
                            name=prefix_cuffdiff,
                            code_module='bsf.runnables.generic',
                            working_directory=self.genome_directory,
                            file_path_dict=file_path_dict_cuffdiff,
                            debug=self.debug))
            executable_run_cuffdiff = self.set_drms_runnable(
                    drms=drms_run_cuffdiff,
                    runnable=runnable_run_cuffdiff)

            if run_cuffquant_before_cuffdiff:
                # Add all executable_run_cuffquant dependencies to the executable_run_cuffdiff process.
                executable_run_cuffdiff.dependencies.extend(cuffdiff_cuffnorm_dependencies)
            else:
                # Add the executable_run_cuffmerge dependency to the executable_run_cuffdiff process.
                executable_run_cuffdiff.dependencies.append(executable_run_cuffmerge.name)

            # Create a new Cuffdiff RunnableStep.

            cuffdiff = runnable_run_cuffdiff.add_runnable_step(
                    runnable_step=RunnableStep(
                            name='cuffdiff',
                            program='cuffdiff'))

            # Set Cuffdiff options.

            cuffdiff.add_option_long(
                    key='output-dir',
                    value=file_path_dict_cuffdiff['output_directory'])
            cuffdiff.add_option_long(
                    key='labels',
                    value=','.join(cuffdiff_cuffnorm_labels))
            cuffdiff.add_option_long(
                    key='num-threads',
                    value=str(drms_run_cuffdiff.threads))
            if self.mask_gtf_path:
                cuffdiff.add_option_long(
                        key='mask-file',
                        value=self.mask_gtf_path)
            cuffdiff.add_option_long(
                    key='frag-bias-correct',
                    value=self.genome_fasta_path)
            if self.multi_read_correction:
                cuffdiff.add_switch_long(
                        key='multi-read-correct')
            if self.library_type:
                cuffdiff.add_option_long(
                        key='library-type',
                        value=self.library_type)
            if self.no_length_correction:
                cuffdiff.add_switch_long(
                        key='no-length-correction')
            cuffdiff.add_switch_long(
                    key='quiet')
            cuffdiff.add_switch_long(
                    key='no-update-check')

            # Add the Cuffmerge GTF file as first Cuffdiff argument.
            cuffdiff.arguments.append(file_path_dict_cuffdiff['merged_gtf'])

            # Cuffdiff seems to have a problem with Cuffquant abundances files in that the isoforms.count_tracking
            # files show ridiculously low numbers such as 1e-319 for some splice variants. Usually, other splice
            # variants in the same cluster seem fine.
            if run_cuffquant_before_cuffdiff:
                # Add the Cuffquant abundances files per comparison group as Cuffdiff arguments.
                for per_group_abundances_list in cuffdiff_cuffnorm_abundances:
                    assert isinstance(per_group_abundances_list, list)
                    cuffdiff.arguments.append(','.join(per_group_abundances_list))
            else:
                # Add the TopHat BAM files per comparison group as Cuffdiff arguments.
                for per_group_alignments_list in cuffdiff_cuffnorm_alignments:
                    assert isinstance(per_group_alignments_list, list)
                    cuffdiff.arguments.append(','.join(per_group_alignments_list))

            # Create a new rnaseq_process_cuffdiff Executable.

            process_cuffdiff = drms_process_cuffdiff.add_executable(
                    executable=Executable(
                            name='_'.join((drms_process_cuffdiff.name, comparison_key)),
                            program='bsf_rnaseq_process_cuffdiff.R'))
            process_cuffdiff.dependencies.append(executable_run_cuffdiff.name)

            # Set rnaseq_process_cuffdiff options.
            self.set_command_configuration(command=process_cuffdiff)
            process_cuffdiff.add_option_long(
                    key='comparison-name',
                    value=comparison_key)
            process_cuffdiff.add_option_long(
                    key='gtf-file',
                    value=file_path_dict_cuffdiff['merged_gtf'])
            process_cuffdiff.add_option_long(
                    key='genome-version',
                    value=self.genome_version)

            # Set rnaseq_process_cuffdiff arguments.

            # None so far.

        return

    def report(self):
        """Create a C{Tuxedo} report in HTML format and a UCSC Genome Browser Track Hub.
        @return:
        @rtype:
        """

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

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

        output += defaults.web.html_header(title='{} RNA-Seq Analysis'.format(self.project_name))
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
        # http://tophat.cbcb.umd.edu/manual.html
        output += '<strong><a href="http://ccb.jhu.edu/software/tophat/index.shtml">TopHat</a></strong> '
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
        output += 'UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.\n'.format(
                defaults.web.ucsc_track_url(
                        options_dict=options_dict,
                        host_name=default.ucsc_host_name),
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
        output += 'UCSC Genome Browser Track Hub <a href="{}" target="UCSC">{}</a>.\n'.format(
                defaults.web.ucsc_track_url(
                        options_dict=options_dict,
                        host_name=default.ucsc_host_name),
                self.project_name)
        output += '</p>\n'
        output += '\n'

        # output += '<p>\n'
        # output += 'Follow the links below to attach\n'
        # output += 'Tophat junction, deletion and insertion annotation to the\n'
        # output += 'UCSC Genome Browser. Since each file needs transferring to\n'
        # output += 'the UCSC site, subsequent pages will take some time to load.\n'
        # output += '</p>\n'

        output += '<h2>Gene Expression Profiles</h2>\n'
        output += '\n'

        output += '<p>\n'
        # http://cufflinks.cbcb.umd.edu/howitworks.html
        output += '<strong><a href="http://cole-trapnell-lab.github.io/cufflinks/">Cufflinks</a></strong>\n'
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
        # http://cufflinks.cbcb.umd.edu/manual.html#fpkm_tracking_format
        output += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                  'fpkm-tracking-format">FPKM Tracking format</a>.\n'
        output += 'The isofroms.count_tracking and genes.count_tracking files\n'
        output += 'contain the scaled isofrom or gene count values in the generic\n'
        output += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                  'count-tracking-format">Count Tracking format</a>.\n'
        output += '</p>\n'

        output += '<p>\n'
        output += 'Please see a more detailed description of\n'
        # http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_output
        output += '<a href="http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#' \
                  'output-formats-used-in-the-cufflinks-suite">Cufflinks output</a>.\n'
        output += '</p>\n'

        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr>\n'
        output += '<th>Sample</th>\n'
        output += '<th>Assembled Transcripts</th>\n'
        output += '<th>Gene FPKM</th>\n'
        output += '<th>Transcript FPKM</th>\n'
        output += '<th>Genes (Symbols)</th>\n'
        output += '<th>Isoforms (Symbols)</th>\n'
        output += '<th>Aligned BAM file</th>\n'
        output += '<th>Aligned BAI file</th>\n'
        output += '<th>Unaligned BAM file</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

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
            assert isinstance(sample, Sample)
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

            # bsf.data.Sample.get_all_paired_reads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of bsf.data.PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                assert isinstance(replicate_key, str)
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
                options_dict['hgt.customText'] = '{}/{}/{}/rnaseq_cufflinks_{}/transcripts.gtf'.format(
                        Default.url_absolute_projects(), link_name, self.genome_version, replicate_key)
                if ucsc_location:
                    options_dict['position'] = ucsc_location

                # UCSC track dictionary.

                track_dict = dict()
                track_dict['name'] = '{}_transcripts'. \
                    format(replicate_key)
                track_dict['description'] = '"{} Cufflinks RNA-Seq transcript assembly"'. \
                    format(replicate_key)
                track_dict['track_type'] = 'gtf'
                track_dict['visibility'] = 'squish'
                track_dict['color'] = '0,0,0'
                track_dict['db'] = self.genome_version

                prefix = 'rnaseq_cufflinks_{}'.format(replicate_key)

                output += '<tr>\n'
                output += '<td>{}</td>\n'.format(replicate_key)
                output += '<td><a href="{}">Transcript Assembly</a></td>\n'.format(
                        defaults.web.ucsc_track_url(
                                options_dict=options_dict,
                                # It does not seem as if a GTF file could be attached via a track line.
                                # The URL has to be supplied via the hgt.customText option.
                                # To get around this problem the GTF file would have to be converted into BigBED format,
                                # yet convincing tools to do this do not seem to be available.
                                # track_dict=track_dict,
                                host_name=default.ucsc_host_name))
                output += '<td><a href="./{}/genes.fpkm_tracking">Genes FPKM</a></td>\n'. \
                    format(prefix)
                output += '<td><a href="./{}/isoforms.fpkm_tracking">Isoforms FPKM</a></td>\n'. \
                    format(prefix)
                # Add files from bsf_process_cufflinks.R
                output += '<td><a href="./{}/{}_genes_fpkm_tracking.tsv">Genes (Symbols)</a></td>\n'. \
                    format(prefix, prefix)
                output += '<td><a href="./{}/{}_isoforms_fpkm_tracking.tsv">Isoforms (Symbols)</a></td>\n'. \
                    format(prefix, prefix)
                output += '<td><a href="./{}/rnaseq_tophat_{}_accepted_hits.bam">Aligned BAM</a></td>\n'. \
                    format(prefix, replicate_key)
                output += '<td><a href="./{}/rnaseq_tophat_{}_accepted_hits.bam.bai">Aligned BAI</a></td>\n'. \
                    format(prefix, replicate_key)
                output += '<td><a href="./{}/rnaseq_tophat_{}_unaligned.bam">Unaligned BAM</a></td>\n'. \
                    format(prefix, replicate_key)
                output += '</tr>\n'

        output += '</tbody>\n'
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

        output += '<h3>All Genes</h3>\n'

        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr>\n'
        output += '<th>Comparison</th>\n'
        output += '<th>Samples</th>\n'
        output += '<th>Replicates</th>\n'
        output += '<th>Coding Sequences</th>\n'
        output += '<th>Genes</th>\n'
        output += '<th>Isoforms</th>\n'
        output += '<th>Promoters</th>\n'
        output += '<th>Splicing</th>\n'
        output += '<th>Transcription Start Sites</th>\n'
        output += '<th>Gene FPKM Replicates</th>\n'
        output += '<th>Gene Count Replicates</th>\n'
        output += '<th>Isoform FPKM Replicates</th>\n'
        output += '<th>Isoform Count Replicates</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        comparison_keys = self.comparisons.keys()
        comparison_keys.sort(cmp=lambda x, y: cmp(x, y))

        for comparison_key in comparison_keys:
            assert isinstance(comparison_key, str)
            prefix = 'rnaseq_process_cuffdiff_{}'.format(comparison_key)

            output += '<tr>\n'
            output += '<td>{}</td>\n'. \
                format(comparison_key)

            # Link to comparison-specific symbolic links in the directory after cummeRbund processing.

            output += '<td><a href="./{}/{}_samples.tsv">Samples</a></td>'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_replicates.tsv">Replicates</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_cds_exp_diff.tsv">Coding Sequences</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_genes_exp_diff.tsv"><strong>Genes</strong></a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_isoforms_exp_diff.tsv">Isoforms</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_promoters_diff.tsv">Promoters</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_splicing_diff.tsv">Splicing</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_tss_group_exp_diff.tsv">Transcription Start Sites</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_genes_fpkm_replicates.tsv">Gene FPKM Replicates</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_genes_counts_replicates.tsv">Gene Count Replicates</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_isoforms_fpkm_replicates.tsv">Isoform FPKM Replicates</a></td>\n'. \
                format(prefix, prefix)
            output += '<td><a href="./{}/{}_isoforms_counts_replicates.tsv">Isoform Count Replicates</a></td>\n'. \
                format(prefix, prefix)

            output += '</tr>\n'

            # Read sample pair information if available.

            sample_pair_path = os.path.join(
                    self.genome_directory,
                    prefix,
                    '_'.join((prefix, 'sample_pairs.tsv')))

            if os.path.exists(sample_pair_path):

                sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                for row_dict in sample_pair_sheet.row_dicts:
                    assert isinstance(row_dict, dict)
                    output += '<tr>\n'
                    output += '<td></td>'  # Comparison
                    output += '<td colspan="3"><strong>{}</strong> versus <strong>{}</strong></td>\n'. \
                        format(row_dict['V1'], row_dict['V2'])  # Sample
                    output += '<td><a href="./{}/{}_{}_{}_genes_diff.tsv"><strong>Genes</strong></a></td>\n'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '<td><a href="./{}/{}_{}_{}_isoforms_diff.tsv">Isoforms</a></td>\n'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '<td colspan="5"></td>'
                    output += '</tr>\n'

        output += '</tbody>\n'
        output += '</table>\n'
        output += '\n'

        output += '<h3>Significant Genes</h3>\n'

        output += '<table>\n'
        output += '<thead>\n'
        output += '<tr>\n'
        output += '<th>Comparison</th>\n'
        output += '<th>Genes</th>\n'
        output += '<th>Isoforms</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        for comparison_key in comparison_keys:
            assert isinstance(comparison_key, str)
            prefix = 'rnaseq_process_cuffdiff_{}'.format(comparison_key)

            output += '<tr>\n'
            output += '<td>{}</td>\n'.format(comparison_key)

            output += '<td><a href="./{}/{}_genes_significance_matrix.pdf">'.format(prefix, prefix)
            output += '<img alt="Significance Matrix Plot - Genes - {}" ' \
                      'src="./{}/{}_genes_significance_matrix.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_isoforms_significance_matrix.pdf">'.format(prefix, prefix)
            output += '<img alt="Significance Matrix Plot - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_significance_matrix.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '</tr>\n'

        output += '</tbody>\n'
        output += '</table>\n'

        # Show cummeRbund quality plots.

        output += '<h2>Quality Plots</h2>\n'
        output += '\n'

        output += '<p>\n'
        output += '</p>\n'
        output += '\n'

        output += '<table>\n'
        output += '<thead>\n'
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
        output += '<th>Box Plot without Replicates - Isoforms</th>\n'
        output += '<th>Box Plot with Replicates - Isoforms</th>\n'
        output += '<th>Scatter Matrix Plot - Genes</th>\n'
        output += '<th>Scatter Matrix Plot - Isoforms</th>\n'
        output += '<th>Dendrogram Plot</th>\n'
        output += '<th>Volcano Matrix Plot - Genes</th>\n'
        output += '<th>Multidimensional Scaling Plot - Genes</th>\n'
        output += '<th>Principal Component Analysis Plot - Genes</th>\n'
        output += '</tr>\n'
        output += '</thead>\n'
        output += '<tbody>\n'

        for comparison_key in comparison_keys:
            assert isinstance(comparison_key, str)
            prefix = 'rnaseq_process_cuffdiff_{}'.format(comparison_key)

            output += '<tr>\n'
            output += '<td>{}</td>\n'.format(comparison_key)

            # Dispersion Plots for Genes and Isoforms

            output += '<td><a href="./{}/{}_genes_dispersion.pdf">'.format(prefix, prefix)
            output += '<img alt="Dispersion Plot - Genes - {}" ' \
                      'src="./{}/{}_genes_dispersion.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_isoforms_dispersion.pdf">'.format(prefix, prefix)
            output += '<img alt="Dispersion Plot - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_dispersion.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Squared Coefficient of Variation (SCV) Plots for Genes and Isoforms

            if os.path.exists(
                    path=os.path.join(self.genome_directory, './{}/{}_genes_scv.png'.format(prefix, prefix))):
                output += '<td><a href="./{}/{}_genes_scv.pdf">'.format(prefix, prefix)
                output += '<img alt="Squared Coefficient of Variation (SCV) - Genes - {}" ' \
                          'src="./{}/{}_genes_scv.png" height="80" width="80" />'. \
                    format(comparison_key, prefix, prefix)
                output += '</a></td>\n'
            else:
                output += '<td></td>\n'

            if os.path.exists(
                    path=os.path.join(self.genome_directory, './{}/{}_isoforms_scv.png'.format(prefix, prefix))):
                output += '<td><a href="./{}/{}_isoforms_scv.pdf">'.format(prefix, prefix)
                output += '<img alt="Squared Coefficient of Variation (SCV) - Isoforms - {}" ' \
                          'src="./{}/{}_isoforms_scv.png" height="80" width="80" />'. \
                    format(comparison_key, prefix, prefix)
                output += '</a></td>\n'
            else:
                output += '<td></td>\n'

            # Density Plots for Genes without and with Replicates

            output += '<td><a href="./{}/{}_genes_density_wo_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Density Plot without Replicates - Genes- {}" ' \
                      'src="./{}/{}_genes_density_wo_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_genes_density_w_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Density Plot with Replicates - Genes - {}" ' \
                      'src="./{}/{}_genes_density_w_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Density Plots for Isoforms without and with Replicates

            output += '<td><a href="./{}/{}_isoforms_density_wo_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Density Plot without Replicates - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_density_wo_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_isoforms_density_w_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Density Plot with Replicates - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_density_w_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Box Plots for Genes without and with Replicates

            output += '<td><a href="./{}/{}_genes_box_wo_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Box Plot without Replicates - Genes - {}" ' \
                      'src="./{}/{}_genes_box_wo_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_genes_box_w_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Box Plot with Replicates - Genes - {}" ' \
                      'src="./{}/{}_genes_box_w_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Box Plots for Isoforms with and without Replicates

            output += '<td><a href="./{}/{}_isoforms_box_wo_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Box Plot without Replicates - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_box_wo_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_isoforms_box_w_replicates.pdf">'.format(prefix, prefix)
            output += '<img alt="Box Plot with Replicates - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_box_w_replicates.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Scatter Matrix Plot for Genes and Isoforms

            output += '<td><a href="./{}/{}_genes_scatter_matrix.pdf">'.format(prefix, prefix)
            output += '<img alt="Scatter Matrix Plot - Genes - {}" ' \
                      'src="./{}/{}_genes_scatter_matrix.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '<td><a href="./{}/{}_isoforms_scatter_matrix.pdf">'.format(prefix, prefix)
            output += '<img alt="Scatter Matrix Plot - Isoforms - {}" ' \
                      'src="./{}/{}_isoforms_scatter_matrix.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Dendrogram Plot for Genes

            output += '<td><a href="./{}/{}_genes_dendrogram.pdf">'.format(prefix, prefix)
            output += '<img alt="Dendrogram Plot - Genes - {}" ' \
                      'src="./{}/{}_genes_dendrogram.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Volcano Matrix Plot for Genes

            output += '<td><a href="./{}/{}_genes_volcano_matrix.pdf">'.format(prefix, prefix)
            output += '<img alt="Volcano Matrix Plot - Genes - {}" ' \
                      'src="./{}/{}_genes_volcano_matrix.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            # Multidimensional Scaling Plot for Genes

            if os.path.exists(
                    path=os.path.join(self.genome_directory, './{}/{}_genes_mds.png'.format(prefix, prefix))):
                output += '<td><a href="./{}/{}_genes_mds.pdf">'.format(prefix, prefix)
                output += '<img alt="Multidimensional Scaling Plot - Genes - {}" ' \
                          'src="./{}/{}_genes_mds.png" height="80" width="80" />'. \
                    format(comparison_key, prefix, prefix)
                output += '</a></td>\n'
            else:
                output += '<td></td>\n'

            # Principal Component Analysis Plot for Genes

            output += '<td><a href="./{}/{}_genes_pca.pdf">'.format(prefix, prefix)
            output += '<img alt="Principal Component Analysis Plot - Genes - {}" ' \
                      'src="./{}/{}_genes_pca.png" height="80" width="80" />'. \
                format(comparison_key, prefix, prefix)
            output += '</a></td>\n'

            output += '</tr>\n'

            # Read sample pair information if available.

            sample_pair_path = os.path.join(self.genome_directory, prefix, '_'.join((prefix, 'sample_pairs.tsv')))

            if os.path.exists(sample_pair_path):

                sample_pair_sheet = TuxedoSamplePairSheet.from_file_path(file_path=sample_pair_path)

                for row_dict in sample_pair_sheet.row_dicts:
                    assert isinstance(row_dict, dict)
                    output += '<tr>\n'

                    output += '<td></td>\n'
                    output += '<td colspan="10"><strong>{}</strong> versus <strong>{}</strong></td>\n'. \
                        format(row_dict['V1'], row_dict['V2'])

                    output += '<td><a href="./{}/{}_{}_{}_genes_scatter.pdf">'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '<img alt="Scatter Plot on genes {} versus {}" ' \
                              'src="./{}/{}_{}_{}_genes_scatter.png" height="80" width="80" />'. \
                        format(row_dict['V1'], row_dict['V2'], prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '</a></td>\n'

                    output += '<td></td>\n'

                    output += '<td><a href="./{}/{}_{}_{}_maplot.pdf">'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '<img alt="M vs A Plot on genes {} versus {}" ' \
                              'src="./{}/{}_{}_{}_maplot.png" height="80" width="80" />'. \
                        format(row_dict['V1'], row_dict['V2'], prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '</a></td>\n'

                    output += '<td><a href="./{}/{}_{}_{}_genes_volcano.pdf">'. \
                        format(prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '<img alt="Volcano Plot on genes {} versus {}" ' \
                              'src="./{}/{}_{}_{}_genes_volcano.png" height="80" width="80" />'. \
                        format(row_dict['V1'], row_dict['V2'], prefix, prefix, row_dict['V1'], row_dict['V2'])
                    output += '</a></td>\n'

                    output += '<td colspan="2"></td>\n'

                    output += '</tr>\n'

        output += '</tbody>\n'
        output += '</table>\n'

        output += '</body>\n'
        output += defaults.web.html_footer()

        file_path = os.path.join(self.genome_directory, 'rnaseq_report.html')

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        # Create the UCSC Genome Browser Track Hub.

        self.ucsc_hub_write_hub(prefix='rnaseq')
        self.ucsc_hub_write_genomes(prefix='rnaseq')
        self.ucsc_hub_write_tracks(output=track_output, prefix='rnaseq')

        return
