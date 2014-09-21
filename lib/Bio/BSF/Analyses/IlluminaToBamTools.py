"""Bio.BSF.Analyses.IlluminaToBamTools

A package of classes and methods supporting analyses of the Illumina2Bam-Tools package.
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


import csv
import os
import string

from Bio.BSF import Analysis, Command, Configuration, Default, DRMS, Executable, Runnable
from Bio.BSF.Data import SampleAnnotationSheet
from Bio.BSF.Illumina import RunFolder


class IlluminaToBam(Analysis):
    """IlluminaToBam Analysis sub-class to convert Illumina BCL to a BAM or SAM file.

    Attributes:
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new IlluminaToBam object from a UNIX-style configuration file via the Configuration class.

        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: IlluminaToBam
        :rtype: IlluminaToBam
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new IlluminaToBam object from a Configuration object.

        :param configuration: Configuration
        :type configuration: Configuration
        :return: IlluminaToBam
        :rtype: IlluminaToBam
        """

        assert isinstance(configuration, Configuration)

        itb = cls(configuration=configuration)

        # A "Bio.BSF.Analyses.IlluminaToBamTools.IlluminaToBam" section specifies defaults
        # for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        itb.set_Configuration(itb.configuration, section=section)

        return itb

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):
        """Initialise a Bio.BSF.Analyses.Illumina2BamTools.IlluminaToBam object.

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
        """

        super(IlluminaToBam, self).__init__(
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

    def set_Configuration(self, configuration, section):

        """Set instance variables of an IlluminaToBam object via a section of a Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param configuration: Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(IlluminaToBam, self).set_Configuration(configuration=configuration, section=section)

        # Sub-class specific ...

    def run(self):

        """Run this IlluminaToBam analysis.

        Convert an Illumina flow-cell into lane-specific archive BAM files.

        To convert an Illumina flow cell, Illumina2bam is run first, setting the SAM Read Group (@RG)
        library name (LB) and sample name (SM) to 'flow-cell identifier.lane'.
        The resulting archive BAM file is then sorted by query name with Picard SortSam.
        Depends on the bsf_illumina2bam.sh Bash script to run this efficiently.
        """

        # Read configuration options.

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        default = Default.get_global_default()

        # Get the Illumina Run Folder

        illumina_run_folder = config_parser.get(section=config_section, option='illumina_run_folder')
        illumina_run_folder = os.path.expanduser(path=illumina_run_folder)
        illumina_run_folder = os.path.expandvars(path=illumina_run_folder)

        if not os.path.isabs(illumina_run_folder):
            illumina_run_folder = os.path.join(Default.absolute_runs_illumina(), illumina_run_folder)

        irf = RunFolder.from_file_path(file_path=illumina_run_folder)

        if not os.path.exists(path=os.path.join(illumina_run_folder, 'RTAComplete.txt')):
            raise Exception('Run folder {!r} is not complete.'.format(illumina_run_folder))

        # The experiment name (e.g. BSF_0000) is used as the prefix for archive BAM files.
        # Read it from the configuration file or from the
        # Run Parameters of the Illumina Run Folder.

        if config_parser.has_option(section=config_section, option='experiment_name'):
            experiment_name = config_parser.get(section=config_section, option='experiment_name')
        else:
            experiment_name = str()

        if not experiment_name:
            experiment_name = irf.run_parameters.get_experiment_name()

        # The project name is a concatenation of the experiment name and the Illumina flow cell identifier.
        # In case it has not been specified in the configuration file, read it from the
        # Run Information of the Illumina Run Folder.

        if not self.project_name:
            self.project_name = string.join(words=(experiment_name, irf.run_information.flow_cell), sep='_')

        # In contrast to normal analyses that use the Default.absolute_samples() as input directory,
        # this one needs to use the Default.absolute_runs_illumina() directory.
        # Make sure an absolute path is defined, before the super-class method Bio.BSF.Analysis.run() is called.

        self.input_directory = os.path.expanduser(path=self.input_directory)
        self.input_directory = os.path.expandvars(path=self.input_directory)

        if not os.path.isabs(self.input_directory):
            self.input_directory = os.path.join(Default.absolute_runs_illumina(), self.input_directory)

        # In contrast to normal analyses that use the Default.absolute_projects() output directory,
        # this one needs to use the Default.absolute_sequences() directory.
        # Make sure an absolute path is defined, before the super-class method Bio.BSF.Analysis.run() is called.

        self.output_directory = os.path.expanduser(path=self.output_directory)
        self.output_directory = os.path.expandvars(path=self.output_directory)

        if not os.path.isabs(self.output_directory):
            self.output_directory = os.path.join(Default.absolute_sequences(), self.output_directory)

        # Get sequencing centre information.

        if config_parser.has_option(section=config_section, option='sequencing_centre'):
            sequencing_centre = config_parser.get(section=config_section, option='sequencing_centre')
        else:
            sequencing_centre = str()

        if not sequencing_centre:
            sequencing_centre = default.operator_sequencing_centre

        # Get the Illumina2bam Java class path.

        if config_parser.has_option(section=config_section, option='classpath_illumina2bam'):
            classpath_illumina2bam = config_parser.get(section=config_section, option='classpath_illumina2bam')
        else:
            classpath_illumina2bam = str()

        if not classpath_illumina2bam:
            classpath_illumina2bam = default.classpath_illumina2bam

        # Get the Picard Java class path.

        if config_parser.has_option(section=config_section, option='classpath_picard'):
            classpath_picard = config_parser.get(section=config_section, option='classpath_picard')
        else:
            classpath_picard = str()
        if not classpath_picard:
            classpath_picard = default.classpath_picard

        # At this stage the configuration is set. Now call run on the super-class to create the necessary directories.

        super(IlluminaToBam, self).run()

        itb_drms = DRMS.from_Analysis(
            name='illumina_to_bam',
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(itb_drms)

        for lane in range(0 + 1, irf.run_information.flow_cell_layout.lane_count + 1):

            lane_str = str(lane)
            prefix = string.join(words=(itb_drms.name, self.project_name, lane_str), sep='_')

            file_path_dict = dict(
                # In contrast to other code modules,
                # sequence archive file names do not use the 'illumina_to_bam' prefix.
                temporary_directory=string.join((self.project_name, lane_str, 'temporary'), sep='_'),
                sorted_bam=string.join((self.project_name, lane_str, 'sorted.bam'), sep='_'),
                sorted_md5=string.join((self.project_name, lane_str, 'sorted.bam.md5'), sep='_'),
                unsorted_bam=string.join((self.project_name, lane_str, 'unsorted.bam'), sep='_'),
                unsorted_md5=string.join((self.project_name, lane_str, 'unsorted.bam.md5'), sep='_')
            )

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            runnable = Runnable(
                name=prefix,
                code_module='Bio.BSF.Runnables.IlluminaToBam',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict)
            self.add_runnable(runnable=runnable)

            # Run Illumina2Bam-Tools Illumina2bam

            java_process = Executable(name='illumina_to_bam', program='java', sub_command=Command(command=str()))
            runnable.add_executable(executable=java_process)

            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_illumina2bam, 'Illumina2bam.jar'))
            java_process.add_SwitchShort(key='Xmx4G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_dict['temporary_directory'])

            sub_command = java_process.sub_command

            # RUN_FOLDER
            sub_command.add_OptionPair(
                key='INTENSITY_DIR',
                value=os.path.join(illumina_run_folder, 'Data', 'Intensities'))
            # BASECALLS_DIR
            sub_command.add_OptionPair(
                key='LANE',
                value=lane_str)
            sub_command.add_OptionPair(
                key='OUTPUT',
                value=file_path_dict['unsorted_bam'])
            sub_command.add_OptionPair(
                key='GENERATE_SECONDARY_BASE_CALLS',
                value='false')
            sub_command.add_OptionPair(
                key='PF_FILTER',
                value='false')
            sub_command.add_OptionPair(
                key='READ_GROUP_ID',
                value=string.join((irf.run_information.flow_cell, lane_str), sep='_'))
            # SAMPLE_ALIAS
            sub_command.add_OptionPair(
                key='LIBRARY_NAME',
                value=string.join((irf.run_information.flow_cell, lane_str), sep='_'))
            # STUDY_NAME
            # PLATFORM_UNIT
            # RUN_START_DATE
            sub_command.add_OptionPair(
                key='SEQUENCING_CENTER',
                value=sequencing_centre)
            # PLATFORM
            # FIRST_TILE
            # TILE_LIMIT
            # BARCODE_SEQUENCE_TAG_NAME
            # BARCODE_QUALITY_TAG_NAME
            # SECOND_BARCODE_SEQUENCE_TAG_NAME
            # SECOND_BARCODE_QUALITY_TAG_NAME
            # FIRST_CYCLE
            # FINAL_CYCLE
            # FIRST_INDEX_CYCLE
            # FINAL_INDEX_CYCLE
            sub_command.add_OptionPair(
                key='TMP_DIR',
                value=file_path_dict['temporary_directory'])
            sub_command.add_OptionPair(
                key='VERBOSITY',
                value='WARNING')
            # QUIET
            # VALIDATION_STRINGENCY
            # COMPRESSION_LEVEL
            sub_command.add_OptionPair(
                key='MAX_RECORDS_IN_RAM',
                value='4000000')
            sub_command.add_OptionPair(
                key='CREATE_INDEX',
                value='false')
            sub_command.add_OptionPair(
                key='CREATE_MD5_FILE',
                value='true')
            # OPTIONS_FILE

            # Run Picard SortSam

            java_process = Executable(
                name='picard_sort_sam',
                program='java',
                sub_command=Command(command=str()))
            runnable.add_executable(executable=java_process)

            java_process.add_SwitchShort(key='d64')
            java_process.add_OptionShort(key='jar', value=os.path.join(classpath_picard, 'SortSam.jar'))
            java_process.add_SwitchShort(key='Xmx6G')
            java_process.add_OptionPair(key='-Djava.io.tmpdir', value=file_path_dict['temporary_directory'])

            sub_command = java_process.sub_command

            sub_command.add_OptionPair(
                key='INPUT',
                value=file_path_dict['unsorted_bam'])
            sub_command.add_OptionPair(
                key='OUTPUT',
                value=file_path_dict['sorted_bam'])
            sub_command.add_OptionPair(
                key='SORT_ORDER',
                value='queryname')
            sub_command.add_OptionPair(
                key='TMP_DIR',
                value=file_path_dict['temporary_directory'])
            sub_command.add_OptionPair(
                key='VERBOSITY',
                value='WARNING')
            # QUIET
            # VALIDATION_STRINGENCY
            # COMPRESSION_LEVEL
            sub_command.add_OptionPair(
                key='MAX_RECORDS_IN_RAM',
                value='4000000')
            sub_command.add_OptionPair(
                key='CREATE_INDEX',
                value='false')
            sub_command.add_OptionPair(
                key='CREATE_MD5_FILE',
                value='true')
            # OPTIONS_FILE

            # Submit the corresponding BSF Executable for the BSF Runner job into the DRMS.
            # Should the Runnable object have dependencies just like the Executable class already has?
            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            itb = Executable.from_analysis_runnable(analysis=self, runnable_name=runnable.name)
            itb_drms.add_Executable(executable=itb)

            # Only submit this Executable if the final result file does not exist.
            if (os.path.exists(
                    os.path.join(self.project_directory, file_path_dict['sorted_md5']))
                and os.path.getsize(
                    os.path.join(self.project_directory, file_path_dict['sorted_md5']))):
                itb.submit = False


class BamIndexDecoder(Analysis):
    """BamIndexDecoder Analysis sub-class to decode sequence archive BAM files into sample-specific BAM files.

    Attributes:
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new BamIndexDecoder object from a UNIX-style configuration file via the Configuration class.

        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BamIndexDecoder
        :rtype: BamIndexDecoder
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new BamIndexDecoder object from a Configuration object.

        :param configuration: Configuration
        :type configuration: Configuration
        :return: BamIndexDecoder
        :rtype: BamIndexDecoder
        """

        assert isinstance(configuration, Configuration)

        itb = cls(configuration=configuration)

        # A "Bio.BSF.Analyses.IlluminaToBamTools.BamIndexDecoder" section specifies defaults
        # for this Analysis sub-class.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        itb.set_Configuration(itb.configuration, section=section)

        return itb

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):
        """Initialise a Bio.BSF.Analyses.Illumina2BamTools.BamIndexDecoder object.

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
        """

        super(BamIndexDecoder, self).__init__(
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

    def set_Configuration(self, configuration, section):

        """Set instance variables of a BamIndexDecoder object via a section of a Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param configuration: Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        super(BamIndexDecoder, self).set_Configuration(configuration=configuration, section=section)

        # Sub-class specific ...

    def run(self):

        """Decode an archive BAM file produced with Illumina2bam into sample-specific BAM files.

        :return: Nothing
        :rtype: None
        """

        # The standard BSF Python *comma-separated* value sample sheet needs to be transformed into
        # a Picard *tab-separated* value sample sheet.
        # lane, barcode_sequence_1, barcode_sequence_2, sample_name, library_name
        # barcode_sequence, barcode_name, library_name, sample_name, description

        # Read configuration options.

        config_parser = self.configuration.config_parser
        config_section = self.configuration.section_from_instance(self)

        default = Default.get_global_default()

        # Load from the configuration file and override with the default if necessary.

        # In contrast to normal analyses that use the Default.absolute_samples() as input directory,
        # this one needs to use the Default.absolute_sequences() directory.
        # Make sure an absolute path is defined, before the super-class method Bio.BSF.Analysis.run() is called.

        self.input_directory = os.path.expanduser(path=self.input_directory)
        self.input_directory = os.path.expandvars(path=self.input_directory)

        if not os.path.isabs(self.input_directory):
            self.input_directory = os.path.join(Default.absolute_sequences(), self.input_directory)

        # In contrast to normal analyses that use the Defaults.absolute_projects() output directory,
        # this one needs to use the Default.absolute_samples() directory.
        # Make sure an absolute path is defined, before the super-class method Bio.BSF.Analysis.run() is called.

        self.output_directory = os.path.expanduser(path=self.output_directory)
        self.output_directory = os.path.expandvars(path=self.output_directory)

        if not os.path.isabs(self.output_directory):
            self.output_directory = os.path.join(Default.absolute_samples(), self.output_directory)

        # Get the library annotation sheet.
        # The library annotation sheet is deliberately not passed in via sas_file,
        # as the Analysis.run() method reads that option into a BSF Collection object.

        if config_parser.has_option(section=config_section, option='library_file'):
            library_file = config_parser.get(section=config_section, option='library_file')
        else:
            library_file = '{}_libraries.csv'.format(self.project_name)

        library_file = os.path.expanduser(path=library_file)
        library_file = os.path.expandvars(path=library_file)

        if not os.path.exists(path=library_file):
            raise Exception(
                'A library_file option is missing from configuration file section {!r}.'.
                format(config_section))

        # Get the Illumina2bam Java class path.

        if config_parser.has_option(section=config_section, option='classpath_illumina2bam'):
            classpath_illumina2bam = config_parser.get(section=config_section, option='classpath_illumina2bam')
        else:
            classpath_illumina2bam = str()

        if not classpath_illumina2bam:
            classpath_illumina2bam = default.classpath_illumina2bam

        # Get the Picard Java class path.

        if config_parser.has_option(section=config_section, option='classpath_picard'):
            classpath_picard = config_parser.get(section=config_section, option='classpath_picard')
        else:
            classpath_picard = str()
        if not classpath_picard:
            classpath_picard = default.classpath_picard

        # At this stage the configuration is set. Now call run on the super-class to create the necessary directories.

        super(BamIndexDecoder, self).run()

        bid_drms = DRMS.from_Analysis(
            name='bam_index_decoder',
            work_directory=self.project_directory,
            analysis=self)
        self.drms_list.append(bid_drms)

        index_by_lane = dict()

        field_names_1 = ['lane', 'barcode_sequence_1', 'barcode_sequence_2', 'sample_name', 'library_name']
        field_names_2 = ['barcode_sequence', 'barcode_name', 'library_name', 'sample_name', 'description']
        field_names_3 = ['ProcessedRunFolder', 'Project', 'Sample', 'Reads1', 'File1']

        sas = SampleAnnotationSheet(file_path=library_file, field_names=field_names_1)

        sas.csv_reader_open()

        for row_dict in sas._csv_reader:

            if row_dict['lane'] in index_by_lane:
                lane_list = index_by_lane[row_dict['lane']]
            else:
                lane_list = list()
                index_by_lane[row_dict['lane']] = lane_list

            lane_list.append(row_dict)

        sas.csv_reader_close()

        sample_csv_name = os.path.join(self.project_directory,
                                       '{}_samples.csv'.format(self.project_name))
        sample_csv_file = open(name=sample_csv_name, mode='wb')
        sample_csv_writer = csv.DictWriter(f=sample_csv_file, fieldnames=field_names_3)
        sample_csv_writer.writeheader()

        keys = index_by_lane.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            # The key represents the lane number as a Python str.

            prefix = string.join(words=(bid_drms.name, self.project_name, key), sep='_')

            file_path_dict = dict(
                temporary_directory=string.join(words=(self.project_name, key, 'temporary'), sep='_'),
                samples_directory=string.join(words=(self.project_name, key, 'samples'), sep='_'),
                barcode=string.join(words=(self.project_name, key, 'barcode.tsv'), sep='_'),
                metrics=string.join(words=(self.project_name, key, 'metrics.tsv'), sep='_'),
                input=os.path.join(self.input_directory,
                                   self.project_name,
                                   string.join(words=(self.project_name, key, 'sorted.bam'), sep='_'))
            )

            require_decoding = 0
            file_handle_barcode = open(name=os.path.join(self.project_directory, file_path_dict['barcode']), mode='w')
            file_handle_barcode.write(string.join(words=field_names_2, sep='\t') + '\n')

            for row_dict in index_by_lane[key]:

                if len(row_dict['barcode_sequence_1']) or len(row_dict['barcode_sequence_2']):
                    require_decoding = 1

                # Write the lane-specific tab-delimited Picard barcode file.
                file_handle_barcode.write(
                    string.join(words=(row_dict['barcode_sequence_1'] + row_dict['barcode_sequence_2'],
                                       row_dict['sample_name'],
                                       row_dict['library_name'],
                                       row_dict['sample_name'],
                                       ''), sep='\t') + '\n')

                # Write the flow-cell-specific sample annotation sheet.
                sample_dict = dict(
                    ProcessedRunFolder=self.project_name,
                    Project=row_dict['library_name'],
                    Sample=row_dict['sample_name'],
                    Reads1=string.join(words=(self.project_name, key, row_dict['sample_name']), sep='_'),
                    File1=os.path.join(self.genome_directory, file_path_dict['samples_directory'],
                                       '{}_{}#{}.bam'.format(self.project_name, key, row_dict['sample_name'])))

                sample_csv_writer.writerow(rowdict=sample_dict)

            file_handle_barcode.close()

            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            runnable = Runnable(
                name=prefix,
                code_module='Bio.BSF.Runnables.BamIndexDecoder',
                working_directory=self.project_directory,
                file_path_dict=file_path_dict)
            self.add_runnable(runnable=runnable)

            if require_decoding:

                # Run the BamIndexDecoder if there is at least one line containing a barcode sequence.

                java_process = Executable(
                    name='bam_index_decoder',
                    program='java',
                    sub_command=Command(command=str()))
                runnable.add_executable(executable=java_process)

                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(
                    key='jar',
                    value=os.path.join(classpath_illumina2bam, 'BamIndexDecoder.jar'))
                java_process.add_SwitchShort(key='Xmx4G')
                java_process.add_OptionPair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict['temporary_directory'])

                sub_command = java_process.sub_command

                sub_command.add_OptionPair(
                    key='INPUT',
                    value=file_path_dict['input'])
                # OUTPUT
                sub_command.add_OptionPair(
                    key='OUTPUT_DIR',
                    value=file_path_dict['samples_directory'])
                sub_command.add_OptionPair(
                    key='OUTPUT_PREFIX',
                    value=string.join(words=(self.project_name, key), sep='_'))
                sub_command.add_OptionPair(
                    key='OUTPUT_FORMAT',
                    value='bam')
                # BARCODE_TAG_NAME
                # BARCODE_QUALITY_TAG_NAME
                # BARCODE
                sub_command.add_OptionPair(
                    key='BARCODE_FILE',
                    value=file_path_dict['barcode'])
                sub_command.add_OptionPair(
                    key='METRICS_FILE',
                    value=file_path_dict['metrics'])
                # MAX_MISMATCHES
                # MIN_MISMATCH_DELTA
                # MAX_NO_CALLS
                # CONVERT_LOW_QUALITY_TO_NO_CALL
                # MAX_LOW_QUALITY_TO_CONVERT
                sub_command.add_OptionPair(
                    key='TMP_DIR',
                    value=file_path_dict['temporary_directory'])
                sub_command.add_OptionPair(
                    key='VERBOSITY',
                    value='WARNING')
                # QUIET
                # VALIDATION_STRINGENCY
                # COMPRESSION_LEVEL
                # MAX_RECORDS_IN_RAM
                sub_command.add_OptionPair(
                    key='CREATE_INDEX',
                    value='false')
                sub_command.add_OptionPair(
                    key='CREATE_MD5_FILE',
                    value='true')
                # OPTIONS_FILE

            else:

                # Run Picard CollectAlignmentSummaryMetrics if there is no line containing a barcode sequence.

                java_process = Executable(
                    name='picard_collect_alignment_summary_metrics',
                    program='java',
                    sub_command=Command(command=str()))
                runnable.add_executable(executable=java_process)

                java_process.add_SwitchShort(key='d64')
                java_process.add_OptionShort(
                    key='jar',
                    value=os.path.join(classpath_picard, 'CollectAlignmentSummaryMetrics.jar'))
                java_process.add_SwitchShort(key='Xmx6G')
                java_process.add_OptionPair(
                    key='-Djava.io.tmpdir',
                    value=file_path_dict['temporary_directory'])

                sub_command = java_process.sub_command

                # MAX_INSERT_SIZE
                # ADAPTER_SEQUENCE
                sub_command.add_OptionPair(
                    key='METRIC_ACCUMULATION_LEVEL',
                    value='READ_GROUP')
                # IS_BISULFITE_SEQUENCED
                sub_command.add_OptionPair(
                    key='INPUT',
                    value=file_path_dict['input'])
                sub_command.add_OptionPair(
                    key='OUTPUT',
                    value=file_path_dict['metrics'])
                # REFERENCE_SEQUENCE
                # ASSUME_SORTED
                # STOP_AFTER
                sub_command.add_OptionPair(
                    key='TMP_DIR',
                    value=file_path_dict['temporary_directory'])
                sub_command.add_OptionPair(
                    key='VERBOSITY',
                    value='WARNING')
                sub_command.add_OptionPair(
                    key='QUIET',
                    value='false')
                sub_command.add_OptionPair(
                    key='VALIDATION_STRINGENCY',
                    value='STRICT')
                sub_command.add_OptionPair(
                    key='COMPRESSION_LEVEL',
                    value='5')
                sub_command.add_OptionPair(
                    key='MAX_RECORDS_IN_RAM',
                    value='4000000')
                sub_command.add_OptionPair(
                    key='CREATE_INDEX',
                    value='true')
                sub_command.add_OptionPair(
                    key='CREATE_MD5_FILE',
                    value='true')
                # OPTIONS_FILE

                # Add a symbolic link to the BSF Sequence Archive file.
                file_path_dict['link_name'] = string.join(words=(self.project_name, key, '.bam'), sep='_')

            # Submit the corresponding BSF Executable for the BSF Runner job into the DRMS.
            # NOTE: The Runnable.name has to match the Executable.name that gets submitted via the DRMS.
            bid = Executable.from_analysis_runnable(analysis=self, runnable_name=runnable.name)
            bid_drms.add_Executable(executable=bid)

            # Since the dependency is managed in the project database of the IlluminaToBam analysis,
            # which is stored in the BSF sequence archive directory, the process identifier cannot be resolved.
            bid.dependencies.append(string.join(words=('illumina_to_bam', self.project_name, key), sep='_'))

            # Only submit this Executable if the final result file does not exist.
            if (os.path.exists(
                    os.path.join(self.project_directory, file_path_dict['metrics']))
                and os.path.getsize(
                    os.path.join(self.project_directory, file_path_dict['metrics']))):
                bid.submit = False

        sample_csv_file.close()
