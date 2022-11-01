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
import logging
import os
from argparse import ArgumentParser

from bsf.cloud import get_azure_blob_service_client, azure_block_blob_download, azure_container_exists

argument_parser = ArgumentParser(
    description='Microsoft Azure Storage Blob Service download script.')

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
    '--blob-path',
    default='.',
    dest='blob_path',
    help='directory path for storing downloaded blobs [.]',
    required=False)

argument_parser.add_argument(
    'blob_item',
    help='blob item (e.g., archive_file.tar.gz)',
    nargs='+')

name_space = argument_parser.parse_args()

if not os.path.isdir(name_space.blob_path):
    raise Exception(f'The blob directory path {name_space.blob_path!r} does not exist.')

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

azure_blob_service_client = get_azure_blob_service_client(account_name=name_space.account)

if not azure_container_exists(
        azure_blob_service_client=azure_blob_service_client,
        container=name_space.container):
    raise Exception(f'The Azure Blob Container {name_space.container!r} does not exist.')

for blob_item in name_space.blob_item:
    logging.debug('Blob item: %r', blob_item)

    for suffix in ('', '.md5'):
        file_path = os.path.join(name_space.blob_path, blob_item + suffix)

        if os.path.exists(file_path):
            logging.info('File exists already: %r', file_path)
            continue

        blob_properties = azure_block_blob_download(
            azure_blob_service_client=azure_blob_service_client,
            container=name_space.container,
            blob=blob_item + suffix,
            file_path=file_path,
            max_concurrency=name_space.threads)

        logging.info('Azure Blob name: %r', blob_properties.name)
        logging.info('Azure Blob size: %r', blob_properties.size)
        logging.info('Azure Blob ETag: %r', blob_properties.etag)
        logging.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())
