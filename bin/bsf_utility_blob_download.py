#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
#  Copyright 2013 - 2021 Michael K. Schuster
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
#  BSF Python script to download from the Microsoft Azure Storage Blob Service.
#
import errno
import os
from argparse import ArgumentParser

from bsf.cloud import get_azure_blob_service_client, azure_block_blob_download
from bsf.standards import Configuration, StandardFilePath as StandardsFilePath

argument_parser = ArgumentParser(
    description='Microsoft Azure Storage Blob Service download script.')

argument_parser.add_argument(
    '--debug',
    default=0,
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--account',
    help='Microsoft Azure Storage Account name',
    required=True,
    type=str)

argument_parser.add_argument(
    '--container',
    help='Microsoft Azure Storage Blob Service container name',
    required=True,
    type=str)

argument_parser.add_argument(
    '--sequence-items',
    dest='sequence_items',
    help='Sequence item names (i.e. experiment_flowcell_lane)',
    required=False,
    type=str)

argument_parser.add_argument(
    '--irf-items',
    dest='irf_items',
    help='Illumina Run Folder item names (i.e. 200623_K00288_0200_BHHCTCBBXY.tar.gz)',
    required=False,
    type=str)

argument_parser.add_argument(
    '--threads',
    default=1,
    help='Maximum number of concurrent download threads',
    required=False,
    type=int)

name_space = argument_parser.parse_args()

azure_blob_service_client = get_azure_blob_service_client(account_name=name_space.account)

# For Sequence items i.e. unaligned BAM and MD5 files.
if name_space.sequence_items:
    sequences_directory = StandardsFilePath.get_sequences()
    for sequence_item in Configuration.list_from_csv(csv_string=name_space.sequence_items):
        print('Sequence item:', sequence_item)
        # Get the experiment directory name by removing the lane suffix.
        experiment_name = '_'.join(sequence_item.split('_')[:-1])
        print('Experiment name:', experiment_name)

        if not os.path.exists(sequences_directory):
            raise Exception('Sequences directory does not exist: ' + repr(sequences_directory))

        experiment_directory = os.path.join(sequences_directory, experiment_name)
        if not os.path.exists(experiment_directory):
            try:
                os.makedirs(experiment_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        for suffix in ('.bam', '.bam.md5'):
            file_path = os.path.join(experiment_directory, sequence_item + suffix)
            if os.path.exists(file_path):
                print('File exists already:', file_path)
            else:
                blob_properties = azure_block_blob_download(
                    azure_blob_service_client=azure_blob_service_client,
                    container=name_space.container,
                    blob=experiment_name + '/' + sequence_item + suffix,
                    file_path=file_path,
                    max_concurrency=name_space.threads)

                print('Azure Blob name:', blob_properties.name)
                print('Azure Blob size:', blob_properties.size)
                print('Azure Blob ETag:', blob_properties.etag)
                print('Azure Blob Last Modified:', blob_properties.last_modified.isoformat())

# For Illumina Run Folders.
if name_space.irf_items:
    run_directory = StandardsFilePath.get_illumina_run(absolute=True)
    for irf_item in Configuration.list_from_csv(csv_string=name_space.irf_items):
        print('IRF item:', irf_item)
        for suffix in ('', '.md5'):
            file_path = os.path.join(run_directory, irf_item + suffix)
            if os.path.exists(file_path):
                print('File exists already:', file_path)
            else:
                blob_properties = azure_block_blob_download(
                    azure_blob_service_client=azure_blob_service_client,
                    container=name_space.container,
                    blob=irf_item + suffix,
                    file_path=file_path,
                    max_concurrency=name_space.threads)

                print('Azure Blob name:', blob_properties.name)
                print('Azure Blob size:', blob_properties.size)
                print('Azure Blob ETag:', blob_properties.etag)
                print('Azure Blob Last Modified:', blob_properties.last_modified.isoformat())
