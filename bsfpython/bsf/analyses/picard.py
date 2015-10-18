"""bsf.analyses.picard

A package of classes and methods modelling Picard analyses data files and data directories.
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


import csv
import errno
import os.path
import warnings
import weakref

from bsf import Analysis, Command, Configuration, Default, DRMS, Executable, Runnable, RunnableStep,\
    RunnableStepMakeDirectory
from bsf.annotation import LibraryAnnotationSheet
from bsf.data import Reads, PairedReads, Sample
from bsf.illumina import RunFolder

import pysam


def _process_row_dict(barcode_dict, row_dict, prefix=None):
    """Private function to read fields from a Python row_dict object, index by the 'lane' field in the barcode_dict.

    @param barcode_dict: A Python dict of 'lane' key data and Python list objects of lane_dict value data
    @type barcode_dict: dict
    @param row_dict: Python row_dict object
    @type row_dict: dict
    @param prefix: Optional prefix
        (e.g. '[Control] lane', ...)
    @type prefix: str
    @return:
    @rtype:
    """

    sample_dict = dict()

    if not prefix:
        prefix = str()

    # Mandatory key 'lane'.
    key1 = str(prefix + ' lane').lstrip(' ')
    sample_dict['lane'] = row_dict[key1]

    # Mandatory key 'barcode_sequence_1'.
    key1 = str(prefix + ' barcode_sequence_1').lstrip(' ')
    sample_dict['barcode_sequence_1'] = row_dict[key1]

    # Optional key 'barcode_sequence_2'.
    key1 = str(prefix + ' barcode_sequence_2').lstrip(' ')
    if key1 in row_dict and row_dict[key1]:
        sample_dict['barcode_sequence_2'] = row_dict[key1]
    else:
        sample_dict['barcode_sequence_2'] = str()

    # Mandatory key 'library_name'.
    key1 = str(prefix + ' library_name').lstrip(' ')
    sample_dict['library_name'] = row_dict[key1]

    # Mandatory key 'sample_name'.
    key1 = str(prefix + ' sample_name').lstrip(' ')
    sample_dict['sample_name'] = row_dict[key1]

    if sample_dict['lane'] in barcode_dict:
        lane_list = barcode_dict[sample_dict['lane']]
        assert isinstance(lane_list, list)
    else:
        lane_list = list()
        barcode_dict[sample_dict['lane']] = lane_list

    lane_list.append(sample_dict)

    return


def extract_illumina_barcodes(config_path):
    """Convert an Illumina Run Folder into BAM files.

    @param config_path: Configuration file
    @type config_path: str | unicode
    @return:
    @rtype:
    """

    default = Default.get_global_default()

    analysis = Analysis.from_config_file_path(config_path=config_path)

    cp = analysis.configuration.config_parser
    section = '.'.join((__name__, analysis.__name__))

    if cp.has_option(section=section, option='max_mismatches'):
        max_mismatches = cp.get(section=section, option='max_mismatches')
    else:
        max_mismatches = str()

    if cp.has_option(section=section, option='min_base_quality'):
        min_base_quality = cp.get(section=section, option='min_base_quality')
    else:
        min_base_quality = str()

    # Get the Illumina Run Folder

    illumina_run_folder = cp.get(section=section, option='illumina_run_folder')
    illumina_run_folder = os.path.expanduser(path=illumina_run_folder)
    illumina_run_folder = os.path.expandvars(path=illumina_run_folder)
    if not os.path.isabs(illumina_run_folder):
        os.path.join(Default.absolute_runs_illumina(), illumina_run_folder)

    irf = RunFolder.from_file_path(file_path=illumina_run_folder)

    base_calls_directory = irf.get_base_calls_directory

    # Read the barcodes file ...

    barcode_path = cp.get(section=section, option='barcode_file')
    barcode_path = os.path.expanduser(barcode_path)
    barcode_path = os.path.expandvars(barcode_path)
    # TODO: Prepend path defaults.

    barcode_dict = dict()

    library_annotation_sheet = LibraryAnnotationSheet.from_file_path(file_path=barcode_path)

    for row_dict in library_annotation_sheet.row_dicts:
        assert isinstance(row_dict, dict)
        _process_row_dict(row_dict=row_dict, prefix=analysis.sas_prefix, barcode_dict=barcode_dict)

    # Picard ExtractIlluminaBarcodes

    eib_drms = analysis.add_drms(drms=DRMS.from_analysis(
        name='ExtractIlluminaBarcodes',
        working_directory=analysis.genome_directory,
        analysis=analysis))

    # Picard IlluminaBasecallsToSam

    ibs_drms = analysis.add_drms(drms=DRMS.from_analysis(
        name='IlluminaBasecallsToSam',
        working_directory=analysis.genome_directory,
        analysis=analysis))

    # For each lane in the barcode_dict ...

    keys = barcode_dict.keys()
    keys.sort(cmp=lambda x, y: cmp(x, y))

    for key in keys:
        assert isinstance(key, str)

        # Make a directory L001 - L008 for each lane.
        # TODO: Check if lane is numeric or a string.

        lane_path = os.path.join(analysis.project_directory, 'L{:03d}'.format(int(key)))

        try:
            os.makedirs(lane_path)
        except OSError as exc:  # Python > 2.5
            if exc.errno == errno.EEXIST and os.path.isdir(lane_path):
                pass
            else:
                raise

        barcodes_path = os.path.join(lane_path, '{}_L{:03d}_input_barcodes.txt'.format(irf.flow_cell, int(key)))
        metrics_path = os.path.join(lane_path, '{}_L{:03d}_output_metrics.txt'.format(irf.flow_cell, int(key)))
        library_path = os.path.join(lane_path, '{}_L{:03d}_library.txt'.format(irf.flow_cell, int(key)))

        lane_list = barcode_dict[key]
        assert isinstance(lane_list, list)

        # Check whether all barcodes are of the same length in a particular lane.

        # TODO: Generalise this to any number of barcodes.
        bc1_length = 0
        bc2_length = 9

        for lane_dict in lane_list:
            assert isinstance(lane_dict, dict)

            if lane_dict['barcode_sequence_1'] == 'NoIndex' or not lane_dict['barcode_sequence_1']:
                bc_length = -1
            else:
                bc_length = len(lane_dict['barcode_sequence_1'])

            if not bc1_length:
                bc1_length = bc_length
            else:
                if bc1_length != bc_length:
                    # Barcode lengths do not match ...
                    warnings.warn(
                        'The length {} of barcode 1 {!r} does not match the length ({}) of other barcodes.'.
                        format(bc_length, lane_dict['barcode_sequence_1'], bc1_length),
                        UserWarning)

            if lane_dict['barcode_sequence_2'] == 'NoIndex' or not lane_dict['barcode_sequence_2']:
                bc_length = -1
            else:
                bc_length = len(lane_dict['barcode_sequence_2'])

            if not bc2_length:
                bc2_length = bc_length
            else:
                if bc2_length != bc_length:
                    # Barcode lengths do not match ...
                    warnings.warn(
                        'The length {} of barcode 2 {!r} does not match the length ({}) of other barcodes.'.
                        format(bc_length, lane_dict['barcode_sequence_2'], bc2_length),
                        UserWarning)

        # TODO: Get the read structure from the IRF and the bc_lengths above ...

        for iread in irf.run_information.reads:
            if iread.number == 1:
                if iread.cycles != bc1_length:
                    # TODO: Finish this!
                    warnings.warn()

        # Create a BARCODE_FILE file for Picard ExtractIlluminaBarcodes.
        field_names = ['barcode_sequence_1', 'barcode_sequence_2', 'barcode_name', 'library_name']
        barcodes_file = open(barcodes_path, 'w')
        csv_writer = csv.DictWriter(f=barcodes_file, fieldnames=field_names, dialect=csv.excel_tab)
        csv_writer.writeheader()
        for lane_dict in lane_list:
            # Create a new row dict to adjust column names to the ones required by ExtractIlluminaBarcodes.
            row_dict = dict()
            row_dict['barcode_sequence_1'] = lane_dict['barcode_sequence_1']
            row_dict['barcode_sequence_2'] = lane_dict['barcode_sequence_2']
            row_dict['barcode_name'] = lane_dict['sample_name']
            row_dict['library_name'] = lane_dict['library_name']
            csv_writer.writerow(rowdict=row_dict)
        barcodes_file.close()

        # Create a LIBRARY_PARAMS file for Picard IlluminaBasecallsToSam.
        field_names = ['OUTPUT', 'SAMPLE_ALIAS', 'LIBRARY_NAME', 'BARCODE_1', 'BARCODE_2']
        library_file = open(library_path, 'w')
        csv_writer = csv.DictWriter(f=library_file, fieldnames=field_names, dialect=csv.excel_tab)
        csv_writer.writeheader()
        for lane_dict in lane_list:
            # Create a new row_dict to adjust column names to the ones required by IlluminaBasecallsToSam.
            row_dict = dict()
            row_dict['OUTPUT'] = '{}_{}_L{:03d}_unmapped.bam'.format(lane_dict['sample_name'], irf.flow_cell, int(key))
            row_dict['SAMPLE_ALIAS'] = lane_dict['sample_name']
            row_dict['LIBRARY_NAME'] = lane_dict['library_name']
            row_dict['BARCODE_1'] = lane_dict['barcode_sequence_1']
            row_dict['BARCODE_2'] = lane_dict['barcode_sequence_2']
            csv_writer.writerow(rowdict=row_dict)
        library_file.close()

        # Picard ExtractIlluminaBarcodes

        eib = eib_drms.add_executable(executable=Executable.from_analysis(
            name='eib',
            program='java',
            analysis=analysis))

        # Set Java options for Picard ExtractIlluminaBarcodes.

        eib.add_switch_short(key='Xmx6G')  # TODO: Make this configurable somewhere ...

        if default.classpath_picard:
            eib.add_option_short(key='cp', value=default.classpath_picard)

        eib.add_option_short(key='jar', value='ExtractIlluminaBarcodes.jar')

        # These Picard 'options' can be just arguments ...

        eib.arguments.append('BASECALLS_DIR={}'.format(base_calls_directory))
        # OUTPUT_DIR for s_l_t_barcode.txt files not set. Defaults to BASECALLS_DIR.
        eib.arguments.append('LANE={}'.format(int(key)))
        eib.arguments.append('READ_STRUCTURE={}'.format())  # TODO: This would also require barcode lengths.
        # BARCODE not set, since BARCODE_FILE gets used.
        eib.arguments.append('BARCODE_FILE={}'.format(barcodes_path))
        eib.arguments.append('METRICS_FILE={}'.format(metrics_path))
        if max_mismatches:
            eib.arguments.append('MAX_MISMATCHES={}').format()
            # MIN_MISMATCH_DELTA
        # MAX_NO_CALLS
        if min_base_quality:
            eib.arguments.append('MINIMUM_BASE_QUALITY={}'.format(min_base_quality))
            # MINIMUM_QUALITY
        # COMPRESS_OUTPUTS for s_l_t_barcode.txt files
        eib.arguments.append('NUM_PROCESSORS={}'.format(int(eib_drms.threads)))
        eib.arguments.append('COMPRESS_OUTPUTS=TRUE')

        # Picard IlluminaBasecallsToSam

        ibs = ibs_drms.add_executable(executable=Executable.from_analysis(
            name='ibs',
            program='java',
            analysis=analysis))

        # Set Java options for Picard IlluminaBasecallsToSam.

        ibs.add_switch_short(key='Xmx6G')  # TODO: Make this configurable somewhere ...

        if default.classpath_picard:
            ibs.add_option_short(key='cp', value=default.classpath_picard)

        ibs.add_option_short(key='jar', value='IlluminaBasecallsToSam.jar')

        # These Picard 'options' can be just arguments ...

        ibs.arguments.append('BASECALLS_DIR={}'.format(base_calls_directory))
        ibs.arguments.append('LANE={}'.format(int(key)))
        ibs.arguments.append('RUN_BARCODE={}'.format())  # TODO:
        ibs.arguments.append('READ_GROUP_ID={}'.format())  # TODO:
        ibs.arguments.append('SEQUENCING_CENTER={}'.format())  # TODO: Get from Defaults?
        ibs.arguments.append('RUN_START_DATE={}'.format())  # TODO:
        ibs.arguments.append('PLATFORM={}'.format())  # TODO: Defaults to illumina
        ibs.arguments.append('READ_STRUCTURE={}'.format())  # TODO
        ibs.arguments.append('LIBRARY_PARAMS={}'.format(library_path))
        ibs.arguments.append('ADAPTERS_TO_CHECK={}'.format())  # TODO
        ibs.arguments.append('NUM_PROCESSORS={}'.format(int(ibs_drms.threads)))

    return


class SamToFastq(Analysis):
    """The C{SamToFastq} class represents the logic to run the Picard SamToFastq analysis.

    Attributes:

    @ivar classpath_picard: Picard tools Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar include_non_pf_reads: Include non-pass filer reads
    @type include_non_pf_reads: bool
    """

    def __init__(self, configuration=None, project_name=None, genome_version=None, input_directory=None,
                 output_directory=None, project_directory=None, genome_directory=None, e_mail=None, debug=0,
                 drms_list=None, collection=None, comparisons=None, samples=None, classpath_picard=None,
                 include_non_pf_reads=False):
        """Initialise a C{SamToFastq} object.

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
        @param comparisons: Python C{dict} of Python C{tuple} objects of C{Sample} objects
        @type comparisons: dict
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list
        @param classpath_picard: Picard tools Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param include_non_pf_reads: Include non-pass filer reads
        @type include_non_pf_reads: bool
        @return:
        @rtype:
        """

        super(SamToFastq, self).__init__(
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

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if include_non_pf_reads is None:
            self.include_non_pass_filter_reads = False
        else:
            assert isinstance(include_non_pf_reads, bool)
            self.include_non_pass_filter_reads = include_non_pf_reads

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{SamToFastq} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(SamToFastq, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the Picard tools Java Archive (JAR) class path directory.

        if configuration.config_parser.has_option(section=section, option='classpath_picard'):
            self.classpath_picard = configuration.config_parser.get(
                section=section,
                option='classpath_picard')

        if configuration.config_parser.has_option(section=section, option='include_non_pass_filter_reads'):
            self.include_non_pass_filter_reads = configuration.config_parser.getboolean(
                section=section,
                option='include_non_pass_filter_reads')

        return

    def _read_comparisons(self):

        self.samples.extend(self.collection.get_all_samples())

        return

    def run(self):
        """Run the C{SamToFastq} C{Analysis} to convert a I{BAM} or I{SAM} file into I{FASTQ} files.

        This method changes the C{Collection} object of this C{Analysis} to update with FASTQ file paths.
        @return:
        @rtype:
        """

        super(SamToFastq, self).run()

        default = Default.get_global_default()

        # Get the Picard tools Java Archive (JAR) class path directory.

        if not self.classpath_picard:
            self.classpath_picard = default.classpath_picard

        self._read_comparisons()

        prune_sas_dependencies = list()

        # Picard SamToFastq

        drms_picard_stf = self.add_drms(drms=DRMS.from_analysis(
            name='picard_sam_to_fastq',
            working_directory=self.project_directory,
            analysis=self))

        for sample in self.samples:
            assert isinstance(sample, Sample)

            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(level=1)

            # bsf.data.Sample.get_all_paired_reads returns a Python dict of
            # Python str key and Python list of Python list objects
            # of bsf.data.PairedReads objects.

            replicate_dict = sample.get_all_paired_reads(replicate_grouping=False)

            replicate_keys = replicate_dict.keys()
            replicate_keys.sort(cmp=lambda x, y: cmp(x, y))

            for replicate_key in replicate_keys:
                assert isinstance(replicate_key, str)

                for paired_reads in replicate_dict[replicate_key]:
                    assert isinstance(paired_reads, PairedReads)

                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    # Apply some sanity checks.

                    if paired_reads.reads2 and not paired_reads.reads1:
                        raise Exception('PairedReads object with reads1 but no reads2 object.', UserWarning)

                    reads = paired_reads.reads1
                    if reads.file_path.endswith('.bam'):
                        bam_file_path = reads.file_path
                        prefix_picard_stf = '_'.join((drms_picard_stf.name, replicate_key))

                        file_path_dict_picard_stf = dict(
                            temporary_directory='_'.join((prefix_picard_stf, 'temporary')),
                            output_directory=os.path.join(self.project_directory, prefix_picard_stf),
                        )

                        # Get the SAM header of a BAM file to extract the read group (@RG), amongst other things.

                        # Open the BAM file, while not checking sequence (@SQ) entries.
                        # De-multiplexed, unaligned BAM files have no reference sequence dictionary.
                        # The check_sq option exists in the calignment code, yet, does not seem to be part of the
                        # function interface.

                        alignment_file = pysam.AlignmentFile(reads.file_path, 'rb', check_sq=False)

                        for read_group in alignment_file.header['RG']:
                            platform_unit = read_group['PU'].replace('#', '_')
                            read_group_list = ['@RG']
                            read_group_list.extend(map(lambda x: '{}:{}'.format(x, read_group[x]), read_group.keys()))
                            if read_group == alignment_file.header['RG'][0]:
                                # For the first read group, modify the PairedReads object in place.
                                paired_reads.read_group = '\\t'.join(read_group_list)
                                paired_reads.reads1.name = platform_unit + '_1'
                                paired_reads.reads1.file_path = os.path.join(
                                    file_path_dict_picard_stf['output_directory'],
                                    platform_unit + '_1.fastq')
                                paired_reads.reads2 = Reads(
                                    file_path=os.path.join(
                                        file_path_dict_picard_stf['output_directory'],
                                        platform_unit + '_2.fastq'),
                                    file_type=paired_reads.reads1.file_type,
                                    name=platform_unit + '_2',
                                    lane=paired_reads.reads1.lane,
                                    read=paired_reads.reads1.read,
                                    chunk=paired_reads.reads1.chunk,
                                    weak_reference_paired_reads=weakref.ref(paired_reads))
                            else:
                                # For further read groups, create new PairedReads objects.
                                reads1 = Reads(
                                    file_path=os.path.join(
                                        file_path_dict_picard_stf['output_directory'],
                                        platform_unit + '_1.fastq'),
                                    file_type=paired_reads.reads1.file_type,
                                    name=platform_unit + '_1',
                                    lane=paired_reads.reads1.lane,
                                    read='R1',
                                    chunk=paired_reads.reads1.chunk)
                                reads2 = Reads(
                                    file_path=os.path.join(
                                        file_path_dict_picard_stf['output_directory'],
                                        platform_unit + '_2.fastq'),
                                    file_type=paired_reads.reads2.file_type,
                                    name=platform_unit + '_2',
                                    lane=paired_reads.reads1.lane,
                                    read='R2',
                                    chunk=paired_reads.reads1.chunk)
                                new_paired_reads = PairedReads(
                                    reads1=reads1,
                                    reads2=reads2,
                                    read_group='\\t'.join(read_group_list))

                                reads1.weak_reference_paired_reads = weakref.ref(new_paired_reads)
                                reads2.weak_reference_paired_reads = weakref.ref(new_paired_reads)
                                sample.add_paired_reads(paired_reads=new_paired_reads)

                        alignment_file.close()

                        # Create a Runnable for running the Picard SamToFastq analysis.

                        runnable_picard_stf = self.add_runnable(runnable=Runnable(
                            name=prefix_picard_stf,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory,
                            file_path_dict=file_path_dict_picard_stf))

                        # Create an Executable for running the Picard SamToFastq Runnable.

                        drms_picard_stf.add_executable(
                            executable=Executable.from_analysis_runnable(
                                analysis=self,
                                runnable_name=runnable_picard_stf.name))

                        # Record the Executable.name for the prune_sas dependency.

                        prune_sas_dependencies.append(runnable_picard_stf.name)

                        # Create a new RunnableStepMakeDirectory in preparation of the Picard program.

                        runnable_picard_stf.add_runnable_step(
                            runnable_step=RunnableStepMakeDirectory(
                                name='mkdir',
                                directory_path=file_path_dict_picard_stf['output_directory']))

                        # Create a new RunnableStep for the Picard SamToFastq program.

                        java_process = runnable_picard_stf.add_runnable_step(runnable_step=RunnableStep(
                            name='sam_to_fastq',
                            program='java',
                            sub_command=Command()))

                        java_process.add_switch_short(
                            key='d64')
                        java_process.add_option_short(
                            key='jar',
                            value=os.path.join(self.classpath_picard, 'SamToFastq.jar'))
                        java_process.add_switch_short(
                            key='Xmx2G')
                        java_process.add_option_pair(
                            key='-Djava.io.tmpdir',
                            value=file_path_dict_picard_stf['temporary_directory'])

                        # Set Picard SamToFastq options.

                        sub_command = java_process.sub_command

                        sub_command.add_option_pair(
                            key='INPUT',
                            value=bam_file_path)
                        sub_command.add_option_pair(
                            key='OUTPUT_PER_RG',
                            value='true')
                        sub_command.add_option_pair(
                            key='OUTPUT_DIR',
                            value=file_path_dict_picard_stf['output_directory'])
                        # RE_REVERSE
                        # INTERLEAVE
                        if self.include_non_pass_filter_reads:
                            sub_command.add_option_pair(
                                key='INCLUDE_NON_PF_READS',
                                value='true')
                        else:
                            sub_command.add_option_pair(
                                key='INCLUDE_NON_PF_READS',
                                value='false')
                        # CLIPPING_ATTRIBUTE
                        # CLIPPING_ACTION
                        # READ1_TRIM
                        # READ1_MAX_BASES_TO_WRITE
                        # READ2_TRIM
                        # READ2_MAX_BASES_TO_WRITE
                        # INCLUDE_NON_PRIMARY_ALIGNMENTS
                        sub_command.add_option_pair(
                            key='TMP_DIR',
                            value=file_path_dict_picard_stf['temporary_directory'])
                        # VERBOSITY defaults to 'INFO'.
                        sub_command.add_option_pair(
                            key='VERBOSITY',
                            value='WARNING')
                        # QUIET defaults to 'false'.
                        sub_command.add_option_pair(
                            key='QUIET',
                            value='false')
                        # VALIDATION_STRINGENCY defaults to 'STRICT'.
                        sub_command.add_option_pair(
                            key='VALIDATION_STRINGENCY',
                            value='STRICT')
                        # COMPRESSION_LEVEL defaults to '5'.
                        # MAX_RECORDS_IN_RAM defaults to '500000'.
                        # CREATE_INDEX defaults to 'false'.
                        # CREATE_MD5_FILE defaults to 'false'.
                        # OPTIONS_FILE

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'sam_to_fastq_original.csv'))),
            name='_'.join((self.project_name, 'sam_to_fastq')))

        annotation_sheet.write_to_file()

        # Create a Runnable for pruning the sample annotation sheet.

        prefix_prune_sas = '_'.join((drms_picard_stf.name, self.project_name))

        file_path_dict_prune_sas = dict(
            temporary_directory='_'.join((prefix_prune_sas, 'temporary')),
            output_directory=prefix_prune_sas,
        )

        runnable_prune_sas = self.add_runnable(runnable=Runnable(
            name=prefix_prune_sas,
            code_module='bsf.runnables.picard_sam_to_fastq_sample_sheet',
            working_directory=self.project_directory,
            file_path_dict=file_path_dict_prune_sas))

        # Create an Executable for running the Runnable for pruning the sample annotation sheet.

        executable_prune_sas = drms_picard_stf.add_executable(
            executable=Executable.from_analysis_runnable(
                analysis=self,
                runnable_name=runnable_prune_sas.name))

        executable_prune_sas.dependencies.extend(prune_sas_dependencies)

        # Create a new RunnableStep.

        prune_sas = runnable_prune_sas.add_runnable_step(
            runnable_step=RunnableStep(
                name='prune_sample_annotation_sheet',
            ))

        prune_sas.add_option_long(key='sas_path', value=annotation_sheet.file_path)

        return annotation_sheet
