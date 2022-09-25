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
import os
from argparse import ArgumentParser

from bsf.cloud import get_azure_blob_service_client, azure_block_blob_upload, azure_container_exists

argument_parser = ArgumentParser(
    description='Microsoft Azure Storage Blob Service upload script.')

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
    '--standard-blob-tier',
    dest='standard_blob_tier',
    help='Microsoft Azure Standard Blob Tier name (i.e., Archive, Cool or Hot)',
    required=False,
    type=str)

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
    help='retain the local path in the blob path')

argument_parser.add_argument(
    'file_path',
    help='file path (e.g., archive_file.tar.gz)',
    nargs='+')

name_space = argument_parser.parse_args()

azure_blob_service_client = get_azure_blob_service_client(account_name=name_space.account)

if not azure_container_exists(
        azure_blob_service_client=azure_blob_service_client,
        container=name_space.container):
    raise Exception(f'Azure Blob Container {name_space.container!r} does not exist.')

for file_path in name_space.file_path:
    if not os.path.isfile(file_path):
        raise Exception('File path ' + repr(file_path) + ' does not exist.')

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

    print('Uploading local path:', repr(file_path), 'blob path:', repr(blob_path))

    if name_space.debug == 0:
        blob_properties = azure_block_blob_upload(
            file_path=file_path,
            azure_blob_service_client=azure_blob_service_client,
            container=name_space.container,
            blob=blob_path,
            standard_blob_tier=name_space.standard_blob_tier,
            max_concurrency=name_space.threads)

        print('  Azure Blob name:', blob_properties.name)
        print('  Azure Blob size:', blob_properties.size)
        print('  Azure Blob ETag:', blob_properties.etag)
        print('  Azure Blob Last Modified:', blob_properties.last_modified.isoformat())
