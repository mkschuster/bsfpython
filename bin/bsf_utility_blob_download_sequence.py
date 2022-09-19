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
#  BSF Python script to download from the Microsoft Azure Storage Blob Service.
#
import errno
import logging
import os
from argparse import ArgumentParser

from bsf.cloud import get_azure_blob_service_client, azure_block_blob_download, azure_container_exists
from bsf.standards import StandardFilePath as StandardsFilePath

argument_parser = ArgumentParser(
    description='Microsoft Azure Storage Blob Service sequence download script.')

argument_parser.add_argument(
    '--account',
    help='Microsoft Azure Storage Account name',
    required=True)

argument_parser.add_argument(
    '--container',
    help='Microsoft Azure Storage Blob Service container name',
    required=True)

argument_parser.add_argument(
    '--logging-level',
    choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
    default='INFO',
    dest='logging_level',
    help='Logging level [INFO]',
    required=False)

argument_parser.add_argument(
    '--threads',
    default=1,
    help='maximum number of concurrent download threads [1]',
    required=False,
    type=int)

argument_parser.add_argument(
    '--sequence-path',
    default=StandardsFilePath.get_sequences(),
    dest='sequence_path',
    help=f'sequence directory path for storing sequence items [{StandardsFilePath.get_sequences()}]',
    required=False)

argument_parser.add_argument(
    'sequence_item',
    dest='sequence_items',
    help='sequence item names (i.e., experiment_flowcell_lane)',
    nargs='+')

name_space = argument_parser.parse_args()

if not os.path.isdir(name_space.sequence_path):
    raise Exception(f'Sequences directory {name_space.sequence_path!r} does not exist.')

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

azure_blob_service_client = get_azure_blob_service_client(account_name=name_space.account)

if not azure_container_exists(
        azure_blob_service_client=azure_blob_service_client,
        container=name_space.container):
    raise Exception(f'Azure Blob Container {name_space.container!r} does not exist.')

for sequence_item in name_space.sequence_item:
    logging.debug('Sequence item: %r', sequence_item)

    # Get the experiment directory name by removing the lane suffix.
    experiment_name = '_'.join(sequence_item.split('_')[:-1])
    logging.debug('Experiment name: %r', experiment_name)

    experiment_directory = os.path.join(name_space.sequence_path, experiment_name)

    if not os.path.exists(experiment_directory):
        try:
            os.makedirs(experiment_directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

    for suffix in ('.bam', '.bam.md5'):
        file_path = os.path.join(experiment_directory, sequence_item + suffix)

        if os.path.exists(file_path):
            logging.info('File exists already: %r', file_path)
            continue

        blob_properties = azure_block_blob_download(
            azure_blob_service_client=azure_blob_service_client,
            container=name_space.container,
            blob=experiment_name + '/' + sequence_item + suffix,
            file_path=file_path,
            max_concurrency=name_space.threads)

        logging.info('Azure Blob name: %r', blob_properties.name)
        logging.info('Azure Blob size: %r', blob_properties.size)
        logging.info('Azure Blob ETag: %r', blob_properties.etag)
        logging.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())
