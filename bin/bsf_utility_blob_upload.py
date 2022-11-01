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
#  BSF Python script to upload to the Microsoft Azure Storage Blob Service.
#
import logging
import os
from argparse import ArgumentParser

from bsf.cloud import get_azure_blob_service_client, azure_block_blob_upload, azure_container_exists

argument_parser = ArgumentParser(
    description='Microsoft Azure Storage Blob Service upload script.')

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
    '--standard-blob-tier',
    dest='standard_blob_tier',
    help='Microsoft Azure Standard Blob Tier name (i.e., Archive, Cool or Hot)',
    required=False)

argument_parser.add_argument(
    '--threads',
    default=1,
    help='maximum number of concurrent download threads',
    required=False,
    type=int)

argument_parser.add_argument(
    '--retain-path',
    action='store_true',
    dest='retain_path',
    help='retain the local path in the blob path',
    required=False)

argument_parser.add_argument(
    'file_path',
    help='file path (e.g., archive_file.tar.gz)',
    nargs='+')

name_space = argument_parser.parse_args()

if name_space.logging_level:
    logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
    logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

    logging.basicConfig(level=name_space.logging_level)

azure_blob_service_client = get_azure_blob_service_client(account_name=name_space.account)

if not azure_container_exists(
        azure_blob_service_client=azure_blob_service_client,
        container=name_space.container):
    raise Exception(f'The Azure Blob Container {name_space.container!r} does not exist.')

for file_path in name_space.file_path:
    if not os.path.isfile(file_path):
        raise Exception(f'The file path {file_path!r} does not exist.')

    if name_space.retain_path:
        # The local path needs rewriting into a URL schema.
        path_list = []
        drive_str, path_str = os.path.splitdrive(file_path)
        while 1:
            path_str, folder_str = os.path.split(path_str)

            if folder_str != '':
                path_list.append(folder_str)
            elif path_str != '':
                path_list.append(path_str)
            else:
                break

        path_list.reverse()

        blob_path = '/'.join(path_list)
    else:
        blob_path = None

    logging.debug('Uploading local path: %r blob path: %r', file_path, blob_path)

    blob_properties = azure_block_blob_upload(
        file_path=file_path,
        azure_blob_service_client=azure_blob_service_client,
        container=name_space.container,
        blob=blob_path,
        standard_blob_tier=name_space.standard_blob_tier,
        max_concurrency=name_space.threads)

    logging.info('Azure Blob name: %r', blob_properties.name)
    logging.info('Azure Blob size: %r', blob_properties.size)
    logging.info('Azure Blob ETag: %r', blob_properties.etag)
    logging.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())
