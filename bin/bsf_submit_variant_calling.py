#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
#
#  Biomedical Sequencing Facility (BSF), part of the genomics core facility
#  of the Research Center for Molecular Medicine (CeMM) of the
#  Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
#  This file is part of BSF Python.
#
#  BSF Python is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BSF Python is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.
#
#
#  BSF Python script to drive the Variant Calling analysis pipeline.
#
import logging
import os
from argparse import ArgumentParser

from bsf.analyses.variant_calling import VariantCallingGATK, FilePathProcessSample, FilePathProcessReadGroup

argument_parser = ArgumentParser(
    description=VariantCallingGATK.name + ' driver script.')

argument_parser.add_argument(
    '--dry-run',
    action='store_false',
    default=True,
    dest='drms_submit',
    help='dry run',
    required=False)

argument_parser.add_argument(
    '--logging-level',
    choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
    default='INFO',
    dest='logging_level',
    help='Logging level [INFO]',
    required=False)

argument_parser.add_argument(
    '--stage',
    help='limit job submission to a particular Analysis stage',
    required=False)

argument_parser.add_argument(
    '--find-missing',
    action='store_true',
    dest='find_missing',
    help='identify missing files',
    required=False)

argument_parser.add_argument(
    'configuration',
    help='configuration (*.ini) file path')

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

analysis = VariantCallingGATK.from_config_file_path(config_path=name_space.configuration)

analysis.run()
analysis.check_state()

if name_space.find_missing:
    # Work through the BSF Sample table to see, which files are missing.

    for sample in analysis.sample_list:
        if sample.is_excluded():
            print("Excluded sample:", sample.name)
            continue

        prefix_sample = VariantCallingGATK.get_prefix_process_sample(sample_name=sample.name)
        file_path_sample = FilePathProcessSample(prefix=prefix_sample)

        sample_missing = False
        for file_name_sample in (
                file_path_sample.realigned_bam,
                file_path_sample.realigned_bai,
                file_path_sample.realigned_md5,
                file_path_sample.realigned_bam_bai,
                file_path_sample.raw_variants_gvcf_vcf,
                file_path_sample.raw_variants_gvcf_tbi,
                file_path_sample.alignment_summary_metrics,
                file_path_sample.duplicate_metrics
        ):
            if not os.path.exists(os.path.join(analysis.genome_directory, file_name_sample)):
                sample_missing = True
                print("  Sample file missing:", file_name_sample)

        if sample_missing:
            read_group_missing = False
            for paired_reads in sample.paired_reads_list:
                paired_reads_name = paired_reads.get_name()
                prefix_process_lane = VariantCallingGATK.get_prefix_process_lane(paired_reads_name=paired_reads_name)
                file_path_read_group = FilePathProcessReadGroup(prefix=prefix_process_lane)
                for file_name_read_group in (
                        file_path_read_group.recalibrated_bam,
                        file_path_read_group.recalibrated_bai,
                        file_path_read_group.recalibrated_md5
                ):
                    if not os.path.exists(os.path.join(analysis.genome_directory, file_name_read_group)):
                        read_group_missing = True
                        print("    ReadGroup file missing:", file_name_read_group)

            if read_group_missing:
                print("Read group(s) missing.")
            else:
                print("Read group(s) complete!")

        if sample_missing:
            print("Missing sample:", sample.name)
        else:
            print("Complete sample:", sample.name)
else:
    analysis.submit(name=name_space.stage, drms_submit=name_space.drms_submit)

    print(analysis.name)
    print('Project name:      ', analysis.project_name)
    print('Genome version:    ', analysis.genome_version)
    print('Input directory:   ', analysis.input_directory)
    print('Project directory: ', analysis.project_directory)
    print('Genome directory:  ', analysis.genome_directory)
    if name_space.stage == 'report':
        print('Report URL:        ', analysis.get_html_report_url())
