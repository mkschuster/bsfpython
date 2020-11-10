#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
#
# BSF Python script to upload to the Microsoft Azure Storage Blob Service.
#
#
# Copyright 2013 - 2019 Michael K. Schuster
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
#
import os
from argparse import ArgumentParser

from bsf.cloud import get_azure_blob_service_client, azure_block_blob_upload

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
    'files',
    help='file paths',
    nargs='+')

name_space = argument_parser.parse_args()

azure_blob_service_client = get_azure_blob_service_client(account_name=name_space.account)

for file_path in name_space.files:
    if not os.path.isfile(path=file_path):
        raise Exception('File path ' + repr(file_path) + ' does not exist.')

    print('Uploading blob file path:', file_path)

    blob_properties = azure_block_blob_upload(
        file_path=file_path,
        azure_blob_service_client=azure_blob_service_client,
        container=name_space.container)

    print('  Azure Blob name:', blob_properties.name)
    print('  Azure Blob size:', blob_properties.size)
    print('  Azure Blob ETag:', blob_properties.etag)
    print('  Azure Blob Last Modified:', blob_properties.last_modified.isoformat())
