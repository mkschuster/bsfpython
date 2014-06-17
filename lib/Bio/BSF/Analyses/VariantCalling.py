"""Bio.BSF.Analysis.VariantCalling

A package of classes and methods supporting variant calling analyses.
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
import warnings

from Bio.BSF import Analysis, Command, Configuration, Default, DRMS, Executable
from Bio.BSF.Data import SampleAnnotationSheet
from Bio.BSF.Executables import BWA


class VariantCallingGATK(Analysis):
    """BSF VariantCallingGATK sub-class.

    Attributes:
    None
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new BSF VariantCallingGATK object from a UNIX-style configuration file via the
        BSF Configuration class.

        :param cls: Class
        :type cls: VariantCallingGATK
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF VariantCallingGATK
        :rtype: VariantCallingGATK
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new BSF VariantCallingGATK object from a BSF Configuration object.

        :param cls: Class
        :type cls: VariantCallingGATK
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF VariantCallingGATK
        :rtype: VariantCallingGATK
        """

        assert isinstance(configuration, Configuration)

        variant_calling = cls(configuration=configuration)

        # A "Bio.BSF.Analysis.VariantCalling.VariantCallingGATK" section specifies defaults for this BSF Analysis.

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
        """Initialise a Bio.BSF.Analysis.VariantCalling.VariantCallingGATK object.

        :param self: BSF VariantCallingGATK
        :type self: VariantCallingGATK
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

        super(VariantCallingGATK, self).__init__(
            configuration=configuration,
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
        """Set instance variables of a BSF VariantCallingGATK object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF VariantCallingGATK
        :type self: VariantCallingGATK
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(VariantCallingGATK, self).set_Configuration(configuration=configuration, section=section)

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

    def _read_vqsr_configuration(self, variation_type=None, gatk_bundle_version=None):
        """Private method to read variant quality score recalibration (VQSR) configuration information.

        :param variation_type: Variation type 'indel' or 'snp'
        :type variation_type: str
        :return: :rtype: Python dict of Python str (resource name) and Python dict values
        :rtype: dict
        """
        if variation_type not in ('indel', 'snp'):
            message = "Variation type has to be 'indel' or 'snp', not {!r}.".format(variation_type)
            warnings.warn(message)

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        vqsr_dict = dict()

        resource_option = "vqsr_resources_{}".format(variation_type)
        if config_parser.has_option(section=config_section, option=resource_option):
            temporary_str = config_parser.get(section=config_section, option=resource_option)
            for resource in string.split(s=temporary_str, sep=','):
                resource_section = "vqsr_{}_{}".format(variation_type, resource)
                if config_parser.has_section(section=resource_section):
                    if resource in vqsr_dict:
                        resource_dict = vqsr_dict[resource]
                    else:
                        resource_dict = dict()
                        vqsr_dict[resource] = resource_dict
                    if config_parser.has_option(section=resource_section, option='known'):
                        resource_dict['known'] = config_parser.get(section=resource_section, option='known')
                    if config_parser.has_option(section=resource_section, option='training'):
                        resource_dict['training'] = config_parser.get(section=resource_section, option='training')
                    if config_parser.has_option(section=resource_section, option='truth'):
                        resource_dict['truth'] = config_parser.get(section=resource_section, option='truth')
                    if config_parser.has_option(section=resource_section, option='prior'):
                        resource_dict['prior'] = config_parser.get(section=resource_section, option='prior')
                    if config_parser.has_option(section=resource_section, option='file_path'):
                        temporary_str = str(config_parser.get(section=resource_section, option='file_path'))
                        if not os.path.isabs(temporary_str):
                            temporary_str = os.path.join(
                                Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                                             genome_version=self.genome_version),
                                temporary_str)
                        resource_dict['file_path'] = temporary_str
                else:
                    message = "Missing configuration section '{}' declared in option resource_section '{}'.". \
                        format(config_section, temporary_str)
                    warnings.warn(message)

        return vqsr_dict

    def run(self):
        """Run this BSF VariantCallingGATK analysis.

        :param self: BSF VariantCallingGATK
        :type self: VariantCallingGATK
        :return: Nothing
        :rtype: None
        """

        super(VariantCallingGATK, self).run()

        # VariantCallingGATK requires a genome version.

        if not self.genome_version:
            message = "A VariantCallingGATK analysis requires a genome_version configuration option."
            raise Exception(message)

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        self.cmp_file = os.path.expanduser(path=self.cmp_file)
        self.cmp_file = os.path.expandvars(path=self.cmp_file)

        if not os.path.isabs(self.cmp_file) and not os.path.exists(self.cmp_file):
            self.cmp_file = os.path.join(self.project_directory, self.cmp_file)

        # Real comparisons would be required for somatic mutation calling.
        self._read_comparisons(cmp_file=self.cmp_file)

        # Experimentally, sort the Python list of BSF Sample objects by the BSF Sample name.
        # This cannot be done in the super-class, because BSF Samples are only put into the Analysis.samples list
        # by the _read_comparisons method.

        self.samples.sort(cmp=lambda x, y: cmp(x.name, y.name))

        # Get global defaults.

        default = Default.get_global_default()

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')

        if config_parser.has_option(section=config_section, option='gatk_bundle_version'):
            gatk_bundle_version = config_parser.get(section=config_section, option='gatk_bundle_version')
        else:
            message = "A VariantCallingGATK analysis requires a gatk_bundle_version configuration option."
            raise Exception(message)

        if config_parser.has_option(section=config_section, option='classpath_gatk'):
            classpath_gatk = config_parser.get(section=config_section, option='classpath_gatk')
        else:
            classpath_gatk = default.classpath_gatk

        if config_parser.has_option(section=config_section, option='classpath_picard'):
            classpath_picard = config_parser.get(section=config_section, option='classpath_picard')
        else:
            classpath_picard = default.classpath_picard

        if config_parser.has_option(section=config_section, option='classpath_snpeff'):
            classpath_snpeff = config_parser.get(section=config_section, option='classpath_snpeff')
        else:
            classpath_snpeff = default.classpath_snpeff

        if config_parser.has_option(section=config_section, option='skip_mark_duplicates'):
            skip_mark_duplicates = config_parser.getboolean(section=config_section, option='skip_mark_duplicates')
        else:
            skip_mark_duplicates = False

        bwa_genome_db = str(config_parser.get(section=config_section, option='bwa_genome_db'))
        if not os.path.isabs(bwa_genome_db):
            bwa_genome_db = os.path.join(
                Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                             genome_version=self.genome_version),
                bwa_genome_db
            )

        snpeff_genome_version = config_parser.get(section=config_section, option='snpeff_genome_version')

        include_intervals_list = list()
        if config_parser.has_option(section=config_section, option='include_intervals'):
            temporary_str = config_parser.get(section=config_section, option='include_intervals')
            include_intervals_list.extend(temporary_str.split(','))

        exclude_intervals_list = list()
        if config_parser.has_option(section=config_section, option='exclude_intervals'):
            temporary_str = config_parser.get(section=config_section, option='exclude_intervals')
            exclude_intervals_list.extend(temporary_str.split(','))

        # Comma-separated list of VCF files with known variant sites for the
        # GATK RealignerTargetCreator and IndelRealigner steps.
        known_sites_realignment = list()
        if config_parser.has_option(section=config_section, option='known_sites_realignment'):
            temporary_str = config_parser.get(section=config_section, option='known_sites_realignment')
            for file_path in temporary_str.split(','):
                if not os.path.isabs(file_path):
                    file_path = os.path.join(Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                                                          genome_version=self.genome_version),
                                             file_path)
                known_sites_realignment.append(file_path)

        # Comma-separated list of VCF files with known variant sites for the
        # GATK BaseRecalibrator and PrintReads steps.
        known_sites_recalibration = list()
        if config_parser.has_option(section=config_section, option='known_sites_recalibration'):
            temporary_str = config_parser.get(section=config_section, option='known_sites_recalibration')
            for file_path in temporary_str.split(','):
                if not os.path.isabs(file_path):
                    file_path = os.path.join(Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                                                          genome_version=self.genome_version),
                                             file_path)
                known_sites_recalibration.append(file_path)

        # Single VCF file of known sites for the
        # GATK HaplotypeCaller and GenotypeGVCFs steps.
        known_sites_discovery = str()
        if config_parser.has_option(section=config_section, option='known_sites_discovery'):
            known_sites_discovery = str(config_parser.get(section=config_section, option='known_sites_discovery'))
            if not os.path.isabs(known_sites_discovery):
                known_sites_discovery = os.path.join(
                    Default.absolute_gatk_bundle(gatk_bundle_version=gatk_bundle_version,
                                                 genome_version=self.genome_version),
                    known_sites_discovery)

        if config_parser.has_option(section=config_section, option='truth_sensitivity_filter_level_snp'):
            truth_sensitivity_filter_level_snp = config_parser.get(section=config_section,
                                                                   option='truth_sensitivity_filter_level_snp')
        else:
            truth_sensitivity_filter_level_snp = str()

        if config_parser.has_option(section=config_section, option='truth_sensitivity_filter_level_indel'):
            truth_sensitivity_filter_level_indel = config_parser.get(section=config_section,
                                                                     option='truth_sensitivity_filter_level_indel')
        else:
            truth_sensitivity_filter_level_indel = str()

        vqsr_annotations_indel_list = list()
        if config_parser.has_option(section=config_section, option='vqsr_annotations_indel'):
            temporary_str = config_parser.get(section=config_section, option='vqsr_annotations_indel')
            for annotation in string.split(s=temporary_str, sep=','):
                vqsr_annotations_indel_list.append(annotation)

        vqsr_annotations_snp_list = list()
        if config_parser.has_option(section=config_section, option='vqsr_annotations_snp'):
            temporary_str = config_parser.get(section=config_section, option='vqsr_annotations_snp')
            for annotation in string.split(s=temporary_str, sep=','):
                vqsr_annotations_snp_list.append(annotation)

        vqsr_resources_indel_dict = self._read_vqsr_configuration(variation_type='indel',
                                                                  gatk_bundle_version=gatk_bundle_version)
        vqsr_resources_snp_dict = self._read_vqsr_configuration(variation_type='snp',
                                                                gatk_bundle_version=gatk_bundle_version)

        # Initialise the Distributed Management System objects for the run_bwa script.

        vc_align_lane_drms = DRMS.from_Analysis(name='variant_calling_align_lane',
                                                work_directory=self.genome_directory,
                                                analysis=self)
        self.drms_list.append(vc_align_lane_drms)

        # Initialise the Distributed Resource Management System object for the
        # bsf_run_variant_calling_process_lane.py script.

        vc_process_lane_drms = DRMS.from_Analysis(name='variant_calling_process_lane',
                                                  work_directory=self.genome_directory,
                                                  analysis=self)
        self.drms_list.append(vc_process_lane_drms)

        # Initialise the Distributed Resource Management System object for the
        # bsf_run_variant_calling_process_sample.py script.

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

            vc_process_sample_dependencies = list()
            vc_process_sample_replicates = list()

            # Bio.BSF.Data.Sample.get_all_PairedReads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of Bio.BSF.Data.PairedReads objects.

            replicate_dict = sample.get_all_PairedReads(replicate_grouping=replicate_grouping)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:

                # Step 1: Process per lane.

                bwa = BWA(name='variant_calling_bwa_{}'.format(replicate_key), analysis=self)
                # Instead of adding the BWA Executable to the DRMS, it gets serialised into the pickler_file.
                # bwa_drms.add_Executable(bwa)

                bwa_mem = bwa.sub_command

                # Set BWA mem options.

                # Allow as many threads as defined in the corresponding DRMS object.
                bwa_mem.add_OptionShort(key='t', value=str(vc_align_lane_drms.threads))
                # Append FASTA/Q comment to SAM output.
                bwa_mem.add_SwitchShort(key='C')
                # Mark  shorter split hits as secondary (for Picard compatibility).
                bwa_mem.add_SwitchShort(key='M')
                # Output warnings and errors only.
                bwa_mem.add_OptionShort(key='v', value='2')

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

                prefix_lane = '{}_{}'.format(vc_process_lane_drms.name, replicate_key)

                # Lane-specific file paths

                file_path_lane = dict(
                    temporary_lane=prefix_lane + '_temporary',
                    # TODO: The name for the aligned BAM is constructed by the bsf_run_bwa.py script.
                    # It is currently based on the vc_align_lane_drms.name and replicate_key.
                    # The script should also be changed to pre-set all file names beforehand.
                    aligned_bam='{}_{}.bam'.format(vc_align_lane_drms.name, replicate_key),
                    aligned_bai='{}_{}.bai'.format(vc_align_lane_drms.name, replicate_key),
                    aligned_md5='{}_{}.bam.md5'.format(vc_align_lane_drms.name, replicate_key),
                    duplicates_marked_bam=prefix_lane + '_duplicates_marked.bam',
                    duplicates_marked_bai=prefix_lane + '_duplicates_marked.bai',
                    duplicates_marked_md5=prefix_lane + '_duplicates_marked.bam.md5',
                    duplicate_metrics=prefix_lane + '_duplicate_metrics.csv',
                    realigner_targets=prefix_lane + '_realigner.intervals',
                    realigned_bam=prefix_lane + '_realigned.bam',
                    realigned_bai=prefix_lane + '_realigned.bai',
                    recalibration_table_pre=prefix_lane + '_recalibration_pre.table',
                    recalibration_table_post=prefix_lane + '_recalibration_post.table',
                    recalibration_plot=prefix_lane + '_recalibration_report.pdf',
                    recalibrated_bam=prefix_lane + '_recalibrated.bam',
                    recalibrated_bai=prefix_lane + '_recalibrated.bai',
                    alignment_summary_metrics=prefix_lane + '_alignment_summary_metrics.csv')

                # Lane-specific pickler_dict

                pickler_dict_process_lane = dict(
                    file_path_dict=file_path_lane,
                    prefix=vc_process_lane_drms.name,
                    replicate_key=replicate_key)

                # Run the Picard MarkDuplicates step, unless configured to skip it.

                if not skip_mark_duplicates:
                    java_process = Executable(name='picard_mark_duplicates',
                                              program='java',
                                              sub_command=Command(command=str()))
                    java_process.add_SwitchShort(key='d64')
                    java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'MarkDuplicates.jar'))
                    java_process.add_SwitchShort(key='Xmx6G')
                    java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                    sub_command = java_process.sub_command
                    sub_command.add_OptionPair(key='INPUT', value=file_path_lane['aligned_bam'])
                    sub_command.add_OptionPair(key='OUTPUT', value=file_path_lane['duplicates_marked_bam'])
                    sub_command.add_OptionPair(key='METRICS_FILE', value=file_path_lane['duplicate_metrics'])
                    # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
                    # as otherwise, optical duplicates could not be detected. This is a consequence of using
                    # Illumina2bam rather than Picard ExtractIlluminaBarcodes, IlluminaBasecallsToFastq and
                    # IlluminaBasecallsToSam.
                    # See BioStar post: http://www.biostars.org/p/12538/
                    # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                    sub_command.add_OptionPair(key='READ_NAME_REGEX',
                                               value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*')
                    sub_command.add_OptionPair(key='TMP_DIR', value=file_path_lane['temporary_lane'])
                    sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
                    sub_command.add_OptionPair(key='QUIET', value='false')
                    sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
                    sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
                    sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
                    sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
                    sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

                    pickler_dict_process_lane[java_process.name] = java_process

                # Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

                java_process = Executable(name='gatk_realigner_target_creator',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionLong(key='analysis_type', value='RealignerTargetCreator')
                sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
                for interval in exclude_intervals_list:
                    sub_command.add_OptionLong(key='excludeIntervals', value=interval)
                for interval in include_intervals_list:
                    sub_command.add_OptionLong(key='intervals', value=interval)
                for file_path in known_sites_realignment:
                    sub_command.add_OptionLong(key='known', value=file_path)
                if skip_mark_duplicates:
                    sub_command.add_OptionLong(key='input_file', value=file_path_lane['aligned_bam'])
                else:
                    sub_command.add_OptionLong(key='input_file', value=file_path_lane['duplicates_marked_bam'])
                sub_command.add_OptionLong(key='out', value=file_path_lane['realigner_targets'])

                pickler_dict_process_lane[java_process.name] = java_process

                # Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

                java_process = Executable(name='gatk_indel_realigner',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionLong(key='analysis_type', value='IndelRealigner')
                sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
                for interval in exclude_intervals_list:
                    sub_command.add_OptionLong(key='excludeIntervals', value=interval)
                for interval in include_intervals_list:
                    sub_command.add_OptionLong(key='intervals', value=interval)
                for file_path in known_sites_realignment:
                    sub_command.add_OptionLong(key='knownAlleles', value=file_path)
                if skip_mark_duplicates:
                    sub_command.add_OptionLong(key='input_file', value=file_path_lane['aligned_bam'])
                else:
                    sub_command.add_OptionLong(key='input_file', value=file_path_lane['duplicates_marked_bam'])
                sub_command.add_OptionLong(key='targetIntervals', value=file_path_lane['realigner_targets'])
                sub_command.add_OptionLong(key='out', value=file_path_lane['realigned_bam'])

                pickler_dict_process_lane[java_process.name] = java_process

                # Run the GATK BaseRecalibrator step as a first-pass walker for the GATK PrintReads step.

                java_process = Executable(name='gatk_base_recalibrator_pre',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionLong(key='analysis_type', value='BaseRecalibrator')
                sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
                for interval in exclude_intervals_list:
                    sub_command.add_OptionLong(key='excludeIntervals', value=interval)
                for interval in include_intervals_list:
                    sub_command.add_OptionLong(key='intervals', value=interval)
                for file_path in known_sites_recalibration:
                    sub_command.add_OptionLong(key='knownSites', value=file_path)
                sub_command.add_OptionLong(key='input_file', value=file_path_lane['realigned_bam'])
                sub_command.add_OptionLong(key='out', value=file_path_lane['recalibration_table_pre'])

                pickler_dict_process_lane[java_process.name] = java_process

                # Run the GATK BaseRecalibrator on-the-fly recalibration step to generate plots.

                java_process = Executable(name='gatk_base_recalibrator_post',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionLong(key='analysis_type', value='BaseRecalibrator')
                sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
                for interval in exclude_intervals_list:
                    sub_command.add_OptionLong(key='excludeIntervals', value=interval)
                for interval in include_intervals_list:
                    sub_command.add_OptionLong(key='intervals', value=interval)
                for file_path in known_sites_recalibration:
                    sub_command.add_OptionLong(key='knownSites', value=file_path)
                sub_command.add_OptionLong(key='BQSR', value=file_path_lane['recalibration_table_pre'])
                sub_command.add_OptionLong(key='input_file', value=file_path_lane['realigned_bam'])
                sub_command.add_OptionLong(key='out', value=file_path_lane['recalibration_table_post'])

                pickler_dict_process_lane[java_process.name] = java_process

                # Run the GATK AnalyzeCovariates step to create a recalibration plot.

                java_process = Executable(name='gatk_analyze_covariates',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionLong(key='analysis_type', value='AnalyzeCovariates')
                sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
                for interval in exclude_intervals_list:
                    sub_command.add_OptionLong(key='excludeIntervals', value=interval)
                for interval in include_intervals_list:
                    sub_command.add_OptionLong(key='intervals', value=interval)
                sub_command.add_OptionLong(key='afterReportFile', value=file_path_lane['recalibration_table_post'])
                sub_command.add_OptionLong(key='beforeReportFile', value=file_path_lane['recalibration_table_pre'])
                sub_command.add_OptionLong(key='plotsReportFile', value=file_path_lane['recalibration_plot'])
                # sub_command.add_OptionLong(key='logging_level', value='DEBUG')

                pickler_dict_process_lane[java_process.name] = java_process

                # Run the GATK PrintReads step as second-pass walker after the BaseRecalibrator step.

                java_process = Executable(name='gatk_print_reads',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionLong(key='analysis_type', value='PrintReads')
                sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
                for interval in exclude_intervals_list:
                    sub_command.add_OptionLong(key='excludeIntervals', value=interval)
                for interval in include_intervals_list:
                    sub_command.add_OptionLong(key='intervals', value=interval)
                sub_command.add_OptionLong(key='input_file', value=file_path_lane['realigned_bam'])
                sub_command.add_OptionLong(key='BQSR', value=file_path_lane['recalibration_table_pre'])
                sub_command.add_OptionLong(key='out', value=file_path_lane['recalibrated_bam'])

                pickler_dict_process_lane[java_process.name] = java_process

                # Run the Picard CollectAlignmentSummaryMetrics step.

                java_process = Executable(name='picard_collect_alignment_summary_metrics',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard,
                                                                           'CollectAlignmentSummaryMetrics.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_lane['temporary_lane'])

                sub_command = java_process.sub_command
                sub_command.add_OptionPair(key='INPUT', value=file_path_lane['recalibrated_bam'])
                sub_command.add_OptionPair(key='OUTPUT', value=file_path_lane['alignment_summary_metrics'])
                sub_command.add_OptionPair(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS')
                sub_command.add_OptionPair(key='REFERENCE_SEQUENCE', value=bwa_genome_db)
                sub_command.add_OptionPair(key='TMP_DIR', value=file_path_lane['temporary_lane'])
                sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
                sub_command.add_OptionPair(key='QUIET', value='false')
                sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
                sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
                sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
                sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
                sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

                pickler_dict_process_lane[java_process.name] = java_process

                # Write the Pickler dict file for processing the lane.

                pickler_path = os.path.join(self.genome_directory, prefix_lane + '.pkl')
                pickler_file = open(pickler_path, 'wb')
                pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
                pickler.dump(obj=pickler_dict_process_lane)
                pickler_file.close()

                # Create a BSF Executable for processing the lane.

                vc_process_lane = Executable.from_Analysis(
                    name=prefix_lane,
                    program='bsf_run_variant_calling_process_lane.py',
                    analysis=self)
                vc_process_lane_drms.add_Executable(vc_process_lane)

                vc_process_lane.dependencies.append(run_bwa.name)

                # Set variant_calling_run_process_lane options.

                vc_process_lane.add_OptionLong(key='pickler_path', value=pickler_path)
                vc_process_lane.add_OptionLong(key='debug', value=str(self.debug))

                # Set dependencies for the next stage.
                vc_process_sample_dependencies.append(vc_process_lane.name)
                # Add the result of the bsf_run_variant_calling_process_lane.py script.
                vc_process_sample_replicates.append(file_path_lane['recalibrated_bam'])

            # Step 2: Process per sample.
            #
            #   Picard MergeSamFiles
            #   Picard MarkDuplicates
            #   GATK RealignerTargetCreator
            #   GATK IndelRealigner
            #   Picard CollectAlignmentSummaryMetrics
            #   GATK HaplotypeCaller

            prefix_sample = '{}_{}'.format(vc_process_sample_drms.name, sample.name)

            file_path_sample = dict(
                temporary_sample=prefix_sample + '_temporary',
                merged_bam=prefix_sample + '_merged.bam',
                merged_bai=prefix_sample + '_merged.bai',
                merged_md5=prefix_sample + '_merged.bam.md5',
                duplicates_marked_bam=prefix_sample + '_duplicates_marked.bam',
                duplicates_marked_bai=prefix_sample + '_duplicates_marked.bai',
                duplicates_marked_md5=prefix_sample + '_duplicates_marked.bam.md5',
                duplicate_metrics=prefix_sample + '_duplicate_metrics.csv',
                realigner_targets=prefix_sample + '_realigner.intervals',
                realigned_bam=prefix_sample + '_realigned.bam',
                realigned_bai=prefix_sample + '_realigned.bai',
                alignment_summary_metrics=prefix_sample + '_alignment_summary_metrics.csv',
                raw_variants=prefix_sample + '_raw_variants.vcf')

            # Sample-specific pickler_dict

            pickler_dict_process_sample = dict(
                file_path_dict=file_path_sample,
                prefix=vc_process_sample_drms.name,
                sample_key=sample.name)

            # Run the Picard MergeSamFiles step.

            java_process = Executable(name='picard_merge_sam_files',
                                      program='java',
                                      sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'MergeSamFiles.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_sample['temporary_sample'])

            sub_command = java_process.sub_command
            for file_path in vc_process_sample_replicates:
                sub_command.add_OptionPair(key='INPUT', value=file_path)
            sub_command.add_OptionPair(key='OUTPUT', value=file_path_sample['merged_bam'])
            sub_command.add_OptionPair(key='COMMENT', value='Merged from the following files:')
            for file_path in vc_process_sample_replicates:
                sub_command.add_OptionPair(key='COMMENT', value=file_path)
            sub_command.add_OptionPair(key='TMP_DIR', value=file_path_sample['temporary_sample'])
            sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
            sub_command.add_OptionPair(key='QUIET', value='false')
            sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
            sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
            sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
            sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
            sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

            pickler_dict_process_sample[java_process.name] = java_process

            # Run the Picard MarkDuplicates step, unless configured to skip it.
            # Optical duplicates should already have been flagged in the lane-specific processing step.

            if not skip_mark_duplicates:
                java_process = Executable(name='picard_mark_duplicates',
                                          program='java',
                                          sub_command=Command(command=str()))
                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'MarkDuplicates.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_sample['temporary_sample'])

                sub_command = java_process.sub_command
                sub_command.add_OptionPair(key='INPUT', value=file_path_sample['merged_bam'])
                sub_command.add_OptionPair(key='OUTPUT', value=file_path_sample['duplicates_marked_bam'])
                sub_command.add_OptionPair(key='METRICS_FILE', value=file_path_sample['duplicate_metrics'])
                # Since read names typically contain a dash and an underscore, the READ_NAME_REGEX needs adjusting,
                # as otherwise optical duplicates could not be detected. This is a consequence of using Illumina2bam
                # rather than Picard ExtractIlluminaBarcodes, IlluminaBasecallsToFastq and IlluminaBasecallsToSam.
                # See BioStar post: http://www.biostars.org/p/12538/
                # Default:  [a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                # Adjusted: [a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*
                sub_command.add_OptionPair(key='READ_NAME_REGEX',
                                           value='[a-zA-Z0-9_-]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*')
                sub_command.add_OptionPair(key='TMP_DIR', value=file_path_sample['temporary_sample'])
                sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
                sub_command.add_OptionPair(key='QUIET', value='false')
                sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
                sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
                sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
                sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
                sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

                pickler_dict_process_sample[java_process.name] = java_process

            # Run the GATK RealignerTargetCreator step as the first-pass walker for the IndelRealigner step.

            java_process = Executable(name='gatk_realigner_target_creator',
                                      program='java',
                                      sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_sample['temporary_sample'])

            sub_command = java_process.sub_command
            sub_command.add_OptionLong(key='analysis_type', value='RealignerTargetCreator')
            sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
            for interval in exclude_intervals_list:
                sub_command.add_OptionLong(key='excludeIntervals', value=interval)
            for interval in include_intervals_list:
                sub_command.add_OptionLong(key='intervals', value=interval)
            for file_path in known_sites_realignment:
                sub_command.add_OptionLong(key='known', value=file_path)
            if skip_mark_duplicates:
                sub_command.add_OptionLong(key='input_file', value=file_path_sample['merged_bam'])
            else:
                sub_command.add_OptionLong(key='input_file', value=file_path_sample['duplicates_marked_bam'])
            sub_command.add_OptionLong(key='out', value=file_path_sample['realigner_targets'])

            pickler_dict_process_sample[java_process.name] = java_process

            # Run the GATK IndelRealigner step as a second-pass walker after the GATK RealignerTargetCreator step.

            java_process = Executable(name='gatk_indel_realigner',
                                      program='java',
                                      sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_sample['temporary_sample'])

            sub_command = java_process.sub_command
            sub_command.add_OptionLong(key='analysis_type', value='IndelRealigner')
            sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
            for interval in exclude_intervals_list:
                sub_command.add_OptionLong(key='excludeIntervals', value=interval)
            for interval in include_intervals_list:
                sub_command.add_OptionLong(key='intervals', value=interval)
            for file_path in known_sites_realignment:
                sub_command.add_OptionLong(key='knownAlleles', value=file_path)
            if skip_mark_duplicates:
                sub_command.add_OptionLong(key='input_file', value=file_path_sample['merged_bam'])
            else:
                sub_command.add_OptionLong(key='input_file', value=file_path_sample['duplicates_marked_bam'])
            sub_command.add_OptionLong(key='targetIntervals', value=file_path_sample['realigner_targets'])
            sub_command.add_OptionLong(key='out', value=file_path_sample['realigned_bam'])
            # For debugging only.
            # sub_command.add_OptionLong(key='logging_level', value='DEBUG')

            pickler_dict_process_sample[java_process.name] = java_process

            # Run the Picard CollectAlignmentSummaryMetrics step.

            java_process = Executable(name='picard_collect_alignment_summary_metrics',
                                      program='java',
                                      sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard,
                                                                       'CollectAlignmentSummaryMetrics.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_sample['temporary_sample'])

            sub_command = java_process.sub_command
            sub_command.add_OptionPair(key='INPUT', value=file_path_sample['realigned_bam'])
            sub_command.add_OptionPair(key='OUTPUT', value=file_path_sample['alignment_summary_metrics'])
            sub_command.add_OptionPair(key='METRIC_ACCUMULATION_LEVEL', value='ALL_READS')
            sub_command.add_OptionPair(key='REFERENCE_SEQUENCE', value=bwa_genome_db)
            sub_command.add_OptionPair(key='TMP_DIR', value=file_path_sample['temporary_sample'])
            sub_command.add_OptionPair(key='VERBOSITY', value='WARNING')
            sub_command.add_OptionPair(key='QUIET', value='false')
            sub_command.add_OptionPair(key='VALIDATION_STRINGENCY', value='STRICT')
            sub_command.add_OptionPair(key='COMPRESSION_LEVEL', value='5')
            sub_command.add_OptionPair(key='MAX_RECORDS_IN_RAM', value='4000000')
            sub_command.add_OptionPair(key='CREATE_INDEX', value='true')
            sub_command.add_OptionPair(key='CREATE_MD5_FILE', value='true')

            pickler_dict_process_sample[java_process.name] = java_process

            # Run the GATK HaplotypeCaller per sample.

            java_process = Executable(name='gatk_haplotype_caller', program='java', sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_sample['temporary_sample'])

            sub_command = java_process.sub_command
            sub_command.add_OptionLong(key='analysis_type', value='HaplotypeCaller')
            sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
            for interval in exclude_intervals_list:
                sub_command.add_OptionLong(key='excludeIntervals', value=interval)
            for interval in include_intervals_list:
                sub_command.add_OptionLong(key='intervals', value=interval)
            sub_command.add_OptionLong(key='genotyping_mode', value='DISCOVERY')
            sub_command.add_OptionLong(key='standard_min_confidence_threshold_for_emitting', value='10')
            sub_command.add_OptionLong(key='standard_min_confidence_threshold_for_calling', value='30')
            sub_command.add_OptionLong(key='emitRefConfidence', value='GVCF')
            if known_sites_discovery:
                sub_command.add_OptionLong(key='dbsnp', value=known_sites_discovery)
            sub_command.add_OptionLong(key='input_file', value=file_path_sample['realigned_bam'])
            sub_command.add_OptionLong(key='out', value=file_path_sample['raw_variants'])
            # Parameter to pass to the VCF/BCF IndexCreator
            sub_command.add_OptionLong(key='variant_index_type', value='LINEAR')
            sub_command.add_OptionLong(key='variant_index_parameter', value='128000')

            pickler_dict_process_sample[java_process.name] = java_process

            # Write the Pickler dict file for processing the sample.

            pickler_path = os.path.join(self.genome_directory, prefix_sample + '.pkl')
            pickler_file = open(pickler_path, 'wb')
            pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
            pickler.dump(obj=pickler_dict_process_sample)
            pickler_file.close()

            # Create a BSF Executable for processing the sample.

            vc_process_sample = Executable.from_Analysis(
                name=prefix_sample,
                program='bsf_run_variant_calling_process_sample.py',
                analysis=self)
            vc_process_sample_drms.add_Executable(vc_process_sample)

            vc_process_sample.dependencies.extend(vc_process_sample_dependencies)

            # Set variant_calling_run_process_sample options.

            vc_process_sample.add_OptionLong(key='pickler_path', value=pickler_path)
            vc_process_sample.add_OptionLong(key='debug', value=str(self.debug))

            # Record dependencies for the next stage.
            vc_process_cohort_dependencies.append(vc_process_sample.name)
            # Add the result of the bsf_run_variant_calling_process_sample.py script.
            vc_process_cohort_replicates.append(file_path_sample['raw_variants'])

        # Step 3: Process per cohort.
        #
        #   GATK CombineGVCFs
        #   GATK GenotypeGVCFs
        #   GATK VariantRecalibrator for SNPs
        #   GATK VariantRecalibrator for INDELs
        #   GATK ApplyRecalibration for SNPs
        #   GATK ApplyRecalibration for INDELs

        prefix_cohort = '{}_{}'.format(vc_process_cohort_drms.name, cohort_name)

        file_path_cohort = dict(
            temporary_cohort=prefix_cohort + '_temporary',
            gvcf_combined=prefix_cohort + '_gvcf_combined.vcf',
            gvcf_genotyped_raw=prefix_cohort + '_gvcf_genotyped_raw_snp_raw_indel.vcf',
            gvcf_recalibrated_snp_raw_indel=prefix_cohort + '_gvcf_recalibrated_snp_raw_indel.vcf',
            gvcf_recalibrated_snp_recalibrated_indel=prefix_cohort + '_gvcf_recalibrated_snp_recalibrated_indel.vcf',
            snpeff_annotated=prefix_cohort + '_snpeff.vcf',
            snpeff_stats=prefix_cohort + '_snpeff_summary.html',
            gvcf_annotated=prefix_cohort + '_gvcf_annotated.vcf',
            recalibration_indel=prefix_cohort + '_recalibration_indel.recal',
            recalibration_snp=prefix_cohort + '_recalibration_snp.recal',
            tranches_indel=prefix_cohort + '_recalibration_indel.tranches',
            tranches_snp=prefix_cohort + '_recalibration_snp.tranches',
            plots_indel=prefix_cohort + '_recalibration_indel.R',
            plots_snp=prefix_cohort + '_recalibration_snp.R')

        # Cohort-specific pickler_dict

        pickler_dict_process_cohort = dict(
            file_path_dict=file_path_cohort,
            prefix=vc_process_cohort_drms.name,
            cohort_key=cohort_name)

        # Run the GATK CombineGVCFs step. It is only required for hierarchically merging samples before GenotypeGVCFs,
        # if too many samples need processing.

        java_process = Executable(name='gatk_combine_gvcfs',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx8G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='CombineGVCFs')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        for file_path in vc_process_cohort_replicates:
            sub_command.add_OptionLong(key='variant', value=file_path)
        sub_command.add_OptionLong(key='out', value=file_path_cohort['gvcf_combined'])

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the GATK GenotypeGVCFs step.

        java_process = Executable(name='gatk_genotype_gvcfs',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx8G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='GenotypeGVCFs')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        if known_sites_discovery:
            sub_command.add_OptionLong(key='dbsnp', value=known_sites_discovery)
        sub_command.add_OptionLong(key='variant', value=file_path_cohort['gvcf_combined'])
        sub_command.add_OptionLong(key='out', value=file_path_cohort['gvcf_genotyped_raw'])

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the GATK VariantRecalibrator for SNPs.

        java_process = Executable(name='gatk_variant_recalibrator_snp',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx8G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='VariantRecalibrator')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        sub_command.add_OptionLong(key='mode', value='SNP')
        for resource in vqsr_resources_snp_dict.keys():
            resource_option = 'resource:{},known={},training={},truth={},prior={}'. \
                format(resource,
                       vqsr_resources_snp_dict[resource]['known'],
                       vqsr_resources_snp_dict[resource]['training'],
                       vqsr_resources_snp_dict[resource]['truth'],
                       vqsr_resources_snp_dict[resource]['prior'])
            sub_command.add_OptionLong(key=resource_option, value=vqsr_resources_snp_dict[resource]['file_path'])
        for annotation in vqsr_annotations_snp_list:
            sub_command.add_OptionLong(key='use_annotation', value=annotation)
        sub_command.add_OptionLong(key='input', value=file_path_cohort['gvcf_genotyped_raw'])
        sub_command.add_OptionLong(key='recal_file', value=file_path_cohort['recalibration_snp'])
        sub_command.add_OptionLong(key='tranches_file', value=file_path_cohort['tranches_snp'])
        sub_command.add_OptionLong(key='rscript_file', value=file_path_cohort['plots_snp'])

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the GATK VariantRecalibrator for INDELs.

        java_process = Executable(name='gatk_variant_recalibrator_indel',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx8G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='VariantRecalibrator')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        sub_command.add_OptionLong(key='mode', value='INDEL')
        for resource in vqsr_resources_indel_dict.keys():
            resource_option = 'resource:{},known={},training={},truth={},prior={}'. \
                format(resource,
                       vqsr_resources_indel_dict[resource]['known'],
                       vqsr_resources_indel_dict[resource]['training'],
                       vqsr_resources_indel_dict[resource]['truth'],
                       vqsr_resources_indel_dict[resource]['prior'])
            sub_command.add_OptionLong(key=resource_option, value=vqsr_resources_indel_dict[resource]['file_path'])
        for annotation in vqsr_annotations_indel_list:
            sub_command.add_OptionLong(key='use_annotation', value=annotation)
        sub_command.add_OptionLong(key='maxGaussians', value='4')  # TODO: Would be good to have this configurable.
        sub_command.add_OptionLong(key='input', value=file_path_cohort['gvcf_genotyped_raw'])
        sub_command.add_OptionLong(key='recal_file', value=file_path_cohort['recalibration_indel'])
        sub_command.add_OptionLong(key='tranches_file', value=file_path_cohort['tranches_indel'])
        sub_command.add_OptionLong(key='rscript_file', value=file_path_cohort['plots_indel'])

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the GATK ApplyRecalibration step for SNPs.

        java_process = Executable(name='gatk_apply_recalibration_snp',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx8G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='ApplyRecalibration')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        sub_command.add_OptionLong(key='mode', value='SNP')
        sub_command.add_OptionLong(key='input', value=file_path_cohort['gvcf_genotyped_raw'])
        sub_command.add_OptionLong(key='recal_file', value=file_path_cohort['recalibration_snp'])
        sub_command.add_OptionLong(key='tranches_file', value=file_path_cohort['tranches_snp'])
        sub_command.add_OptionLong(key='out', value=file_path_cohort['gvcf_recalibrated_snp_raw_indel'])
        # TODO: The lodCutoff (VQSLOD score) filter is not applied for the moment.
        if truth_sensitivity_filter_level_snp:
            sub_command.add_OptionLong(key='ts_filter_level', value=truth_sensitivity_filter_level_snp)

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the GATK ApplyRecalibration step for INDELs.

        java_process = Executable(name='gatk_apply_recalibration_indel',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx8G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='ApplyRecalibration')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        sub_command.add_OptionLong(key='mode', value='INDEL')
        sub_command.add_OptionLong(key='input', value=file_path_cohort['gvcf_recalibrated_snp_raw_indel'])
        sub_command.add_OptionLong(key='recal_file', value=file_path_cohort['recalibration_indel'])
        sub_command.add_OptionLong(key='tranches_file', value=file_path_cohort['tranches_indel'])
        sub_command.add_OptionLong(key='out', value=file_path_cohort['gvcf_recalibrated_snp_recalibrated_indel'])
        # TODO: The lodCutoff (VQSLOD score) filter is not applied for the moment.
        if truth_sensitivity_filter_level_indel:
            sub_command.add_OptionLong(key='ts_filter_level', value=truth_sensitivity_filter_level_indel)

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the snpEff tool for functional variant annotation.

        java_process = Executable(name='snpeff',
                                  program='java',
                                  sub_command=Command(command='eff'))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_snpeff, 'snpEff.jar'))
        java_process.add_SwitchShort(key='Xmx6G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])
        java_process.stdout_path = file_path_cohort['snpeff_annotated']

        sub_command = java_process.sub_command
        sub_command.add_SwitchShort(key='download')
        sub_command.add_OptionShort(key='o', value='gatk')
        sub_command.add_OptionShort(key='stats', value=file_path_cohort['snpeff_stats'])
        sub_command.add_OptionShort(key='config', value=os.path.join(classpath_snpeff, 'snpEff.config'))

        sub_command.arguments.append(snpeff_genome_version)
        sub_command.arguments.append(file_path_cohort['gvcf_recalibrated_snp_recalibrated_indel'])

        pickler_dict_process_cohort[java_process.name] = java_process

        # Run the GATK VariantAnnotator

        java_process = Executable(name='gatk_variant_annotator',
                                  program='java',
                                  sub_command=Command(command=str()))
        java_process.add_SwitchShort(key='d64')
        java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
        java_process.add_SwitchShort(key='Xmx6G')
        java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

        sub_command = java_process.sub_command
        sub_command.add_OptionLong(key='analysis_type', value='VariantAnnotator')
        sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
        for interval in exclude_intervals_list:
            sub_command.add_OptionLong(key='excludeIntervals', value=interval)
        for interval in include_intervals_list:
            sub_command.add_OptionLong(key='intervals', value=interval)
        if known_sites_discovery:
            sub_command.add_OptionLong(key='dbsnp', value=known_sites_discovery)
        sub_command.add_OptionLong(key='variant', value=file_path_cohort['gvcf_recalibrated_snp_recalibrated_indel'])
        sub_command.add_OptionLong(key='annotation', value='SnpEff')
        sub_command.add_OptionLong(key='snpEffFile', value=file_path_cohort['snpeff_annotated'])
        sub_command.add_OptionLong(key='out', value=file_path_cohort['gvcf_annotated'])

        pickler_dict_process_cohort[java_process.name] = java_process

        # Re-process the cohort by sample.

        pickler_dict_process_cohort['sample_names'] = list()

        for sample in self.samples:

            pickler_dict_process_cohort['sample_names'].append(sample.name)

            # Run the GATK SelectVariants step to split multi-sample VCF files into one per sample.

            file_path_cohort['sample_vcf_' + sample.name] = prefix_cohort + '_sample_{}.vcf'.format(sample.name)

            java_process = Executable(name='gatk_select_variants_sample_' + sample.name,
                                      program='java',
                                      sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

            sub_command = java_process.sub_command
            sub_command.add_OptionLong(key='analysis_type', value='SelectVariants')
            sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
            for interval in exclude_intervals_list:
                sub_command.add_OptionLong(key='excludeIntervals', value=interval)
            for interval in include_intervals_list:
                sub_command.add_OptionLong(key='intervals', value=interval)

            sub_command.add_OptionLong(key='variant', value=file_path_cohort['gvcf_annotated'])
            sub_command.add_OptionLong(key='out', value=file_path_cohort['sample_vcf_' + sample.name])
            sub_command.add_OptionLong(key='sample_name', value=sample.name)
            sub_command.add_SwitchLong(key='excludeNonVariants')

            pickler_dict_process_cohort[java_process.name] = java_process

            # Run the GATK VariantsToTable step.

            file_path_cohort['sample_csv_' + sample.name] = prefix_cohort + '_sample_{}.csv'.format(sample.name)

            java_process = Executable(name='gatk_variants_to_table_sample_' + sample.name,
                                      program='java',
                                      sub_command=Command(command=str()))
            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_gatk, 'GenomeAnalysisTK.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_cohort['temporary_cohort'])

            sub_command = java_process.sub_command
            sub_command.add_OptionLong(key='analysis_type', value='VariantsToTable')
            sub_command.add_OptionLong(key='reference_sequence', value=bwa_genome_db)
            for interval in exclude_intervals_list:
                sub_command.add_OptionLong(key='excludeIntervals', value=interval)
            for interval in include_intervals_list:
                sub_command.add_OptionLong(key='intervals', value=interval)

            sub_command.add_OptionLong(key='variant', value=file_path_cohort['sample_vcf_' + sample.name])
            sub_command.add_OptionLong(key='out', value=file_path_cohort['sample_csv_' + sample.name])
            sub_command.add_SwitchLong(key='allowMissingData')
            sub_command.add_SwitchLong(key='showFiltered')
            sub_command.add_OptionLong(key='fields', value='CHROM')
            sub_command.add_OptionLong(key='fields', value='POS')
            sub_command.add_OptionLong(key='fields', value='ID')
            sub_command.add_OptionLong(key='fields', value='REF')
            sub_command.add_OptionLong(key='fields', value='ALT')
            sub_command.add_OptionLong(key='fields', value='QUAL')
            sub_command.add_OptionLong(key='fields', value='FILTER')
            sub_command.add_OptionLong(key='fields', value='AF')
            sub_command.add_OptionLong(key='fields', value='VQSLOD')
            sub_command.add_OptionLong(key='genotypeFields', value='GT')
            sub_command.add_OptionLong(key='genotypeFields', value='AD')
            sub_command.add_OptionLong(key='genotypeFields', value='DP')
            sub_command.add_OptionLong(key='genotypeFields', value='GQ')
            sub_command.add_OptionLong(key='genotypeFields', value='PL')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_EFFECT')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_IMPACT')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_FUNCTIONAL_CLASS')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_CODON_CHANGE')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_AMINO_ACID_CHANGE')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_GENE_NAME')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_GENE_BIOTYPE')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_TRANSCRIPT_ID')
            sub_command.add_OptionLong(key='fields', value='SNPEFF_EXON_ID')

            pickler_dict_process_cohort[java_process.name] = java_process

        # Write the Pickler dict file for processing the cohort.

        pickler_path = os.path.join(self.genome_directory, prefix_cohort + '.pkl')
        pickler_file = open(pickler_path, 'wb')
        pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
        pickler.dump(obj=pickler_dict_process_cohort)
        pickler_file.close()

        # Create a BSF Executable for processing the cohort.

        vc_process_cohort = Executable.from_Analysis(
            name=prefix_cohort,
            program='bsf_run_variant_calling_process_cohort.py',
            analysis=self)
        vc_process_cohort_drms.add_Executable(vc_process_cohort)

        vc_process_cohort.dependencies.extend(vc_process_cohort_dependencies)

        # Set variant_calling_run_process_cohort options.

        vc_process_cohort.add_OptionLong(key='pickler_path', value=pickler_path)
        vc_process_cohort.add_OptionLong(key='debug', value=str(self.debug))
