# -*- coding: utf-8 -*-
"""Cloud services module.

"""
#  Copyright 2013 - 2019 Michael K. Schuster
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

import json
import os

import azure.storage.blob
import azure.storage.blob.models

import bsf.standards


def get_azure_secrets_dict(account_name):
    """Get Microsoft Azure secrets from a JSON configuration file.

    @param account_name: Microsoft Azure account name
    @type account_name: str
    @return: Python C{dict} with Microsoft Azure secrets i.e. I{azure_name} and I{azure_key}
    @rtype: dict | None
    """
    file_path = bsf.standards.Secrets.get_azure_file_path()

    if file_path and os.path.exists(file_path):
        with open(file=file_path, mode='rt') as io_object:
            secrets_dict = json.load(io_object)

        if account_name not in secrets_dict:
            raise Exception('Account name ' + repr(account_name) + ' not in secrets file ' +
                            repr(bsf.standards.Secrets.get_azure_file_path()))
        account_dict = secrets_dict[account_name]

        if 'azure_name' not in account_dict:
            raise Exception('The Azure secrets JSON file requires an "azure_name" dict entry.')

        if 'azure_key' not in account_dict:
            raise Exception('The Azure secrets JSON file requires an "azure_key" dict entry.')

        return account_dict
    else:
        raise Exception('Could not get an Azure secrets JSON file.')


def get_azure_block_blob_service(account_name):
    """Get a Microsoft Azure C{BlockBlobService} object.

    Automatically sets the C{azure.storage.blob.BlockBlobService.MAX_BLOCK_SIZE} class variable
    to support 100 MiB blocks.
    @param account_name: Microsoft Azure account name
    @type account_name: str
    @return: Microsoft Azure C{BlockBlobService} object
    @rtype: azure.storage.blob.BlockBlobService
    """
    secrets_dict = get_azure_secrets_dict(account_name=account_name)

    azure.storage.blob.BlockBlobService.MAX_BLOCK_SIZE = 100 * 1024 * 1024

    return azure.storage.blob.BlockBlobService(
        account_name=secrets_dict['azure_name'],
        account_key=secrets_dict['azure_key'])


def azure_block_blob_upload(block_blob_service, container_name, file_path, blob_name=None):
    """Upload a block blob into a Microsoft Azure I{BlockBlobService} container.

    @param block_blob_service: Microsoft Azure C{BlockBlobService} object
    @type block_blob_service: azure.storage.blob.BlockBlobService
    @param container_name: Container name
    @type container_name: str
    @param file_path: Local file path
    @type file_path: str
    @param blob_name: The Blob name
    @type blob_name: str | None
    @return: A C{azure.storage.blob.models.ResourceProperties} object
    @rtype: azure.storage.blob.models.ResourceProperties
    """
    if not block_blob_service.exists(container_name=container_name):
        raise Exception('Container ' + container_name + ' does not exists.')

    if not blob_name:
        blob_name = os.path.basename(file_path)

    return block_blob_service.create_blob_from_path(
        container_name=container_name,
        blob_name=blob_name,
        file_path=file_path)


def azure_block_blob_download(block_blob_service, container_name, blob_name, file_path=None):
    """Download a block blob from a Microsoft Azure I{BlockBlobService} container.

    @param block_blob_service: Microsoft Azure C{BlockBlobService} object
    @type block_blob_service: azure.storage.blob.BlockBlobService
    @param container_name: Container name
    @type container_name: str
    @param blob_name: The Blob name
    @type blob_name: str | None
    @param file_path: Local file path
    @type file_path: str
    @return: A C{azure.storage.blob.models.Blob} object
    @rtype: azure.storage.blob.models.Blob
    """
    if not block_blob_service.exists(container_name=container_name):
        raise Exception('Container ' + container_name + ' does not exists.')

    if not file_path:
        file_path = os.path.basename(blob_name)

    return block_blob_service.get_blob_to_path(
        container_name=container_name,
        blob_name=blob_name,
        file_path=file_path)
