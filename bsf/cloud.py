# -*- coding: utf-8 -*-
"""Cloud services module.

"""
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
import json
import os
from typing import Dict, IO, Union

from azure.storage.blob import BlobClient, BlobServiceClient, ContainerClient, ContainerProperties, BlobProperties

from bsf.standards import Secrets


def get_azure_secrets_dict(account_name: str) -> Dict[str, str]:
    """Get Microsoft Azure secrets from a JSON configuration file.

    @param account_name: I{Microsoft Azure Storage Account} name
    @type account_name: str
    @return: Python C{dict} with I{Microsoft Azure Storage Account} secrets
    @rtype: dict[str, str]
    """
    file_path = Secrets.get_azure_file_path()

    if file_path and os.path.exists(file_path):
        with open(file=file_path, mode='rt') as text_io:
            secrets_dict: Dict[str, Dict[str, str]] = json.load(fp=text_io)

        if account_name not in secrets_dict:
            raise Exception('Microsoft Azure Storage Account name ' + repr(account_name) + ' not in secrets file ' +
                            repr(file_path))

        return secrets_dict[account_name]
    else:
        raise Exception('Could not get a Microsoft Azure secrets JSON file.')


def get_azure_blob_name(blob: Union[BlobProperties, str]) -> str:
    """Get an Azure blob name from a C{azure.storage.blob.BlobProperties} or a Python C{str} object.

    @param blob: C{azure.storage.blob.BlobProperties} or Python C{str} object
    @type blob: BlobProperties | str
    @return: Blob name
    @rtype: str
    """
    if isinstance(blob, BlobProperties):
        return blob.name
    elif isinstance(blob, str):
        return blob
    else:
        raise Exception("The 'blob' option has to be of type 'BlobProperties' or 'str'.")


def get_container_name(container: Union[ContainerProperties, str]) -> str:
    """Get a container name from a C{azure.storage.blob.ContainerProperties} or a Python C{str} object.

    @param container: C{azure.storage.blob.ContainerProperties} or Python C{str} object
    @type container: ContainerProperties | str
    @return: Container name
    @rtype: str
    """
    if isinstance(container, ContainerProperties):
        return container.name
    elif isinstance(container, str):
        return container
    else:
        raise Exception("The 'container' option has to be of type 'ContainerProperties' or 'str'.")


def get_azure_blob_service_client(
        account_name: str,
        max_block_size: int = 100 * 1024 * 1024,
        retries_total: int = 10,
        logging_enable: bool = False) -> BlobServiceClient:
    """Get an C{azure.storage.blob.BlobServiceClient} object.

    For uploading large files, the default block size of 4 * 1024 * 1024 (4 MiB) is not big enough,
    since the number of blocks in a blob is limited to 50,000, which allows for 200,000 MiB
    or 195 GiB blobs. Changing the maximum block size to 100 MiB allows for 4882 GiB or 4 TiB blobs.
    @param account_name: I{Microsoft Azure Storage Account} name
    @type account_name: str
    @param max_block_size: Maximum block size defaults to 100 MiB, which is als the maximum Azure supports
    @type max_block_size: int
    @param retries_total: Number of retries
    @type retries_total: int
    @param logging_enable: Enable logging
    @type logging_enable: bool
    @return: C{azure.storage.blob.BlobServiceClient} object
    @rtype: BlobServiceClient
    """
    secrets_dict = get_azure_secrets_dict(account_name=account_name)

    if 'azure_connection' not in secrets_dict:
        raise Exception('The Azure secrets JSON file requires an "azure_connection" dict entry.')

    return BlobServiceClient.from_connection_string(
        conn_str=secrets_dict['azure_connection'],
        max_block_size=max_block_size,
        logging_enable=logging_enable,
        retries_total=retries_total)


def get_azure_container_client(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str]) -> ContainerClient:
    """Get an C{azure.storage.blob.ContainerClient} object.

    @param azure_blob_service_client: C{azure.storage.blob.BlobServiceClient} object
    @type azure_blob_service_client: BlobServiceClient
    @param container: C{azure.storage.blob.ContainerProperties} object or Python C{str} container name
    @type container: ContainerProperties | str
    @return: C{azure.storage.blob.ContainerClient} object
    @rtype: ContainerClient
    """
    return azure_blob_service_client.get_container_client(container=container)


def get_azure_blob_client(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str],
        blob: Union[BlobProperties, str]) -> BlobClient:
    """Get an C{azure.storage.blob.BlobClient} object.

    @param azure_blob_service_client: C{azure.storage.blob.BlobServiceClient} object
    @type azure_blob_service_client: BlobServiceClient
    @param container: C{azure.storage.blob.ContainerProperties} object or Python C{str} container name
    @type container: ContainerProperties | str
    @param blob: C{azure.storage.blob.BlobProperties} or Python C{str} blob name
    @type blob: BlobProperties | str
    @return: C{azure.storage.blob.BlobClient} object
    @rtype: BlobClient
    """
    return azure_blob_service_client.get_blob_client(container=container, blob=blob)


def azure_container_exists(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str]) -> bool:
    """Check if a container exists in a I{Microsoft Azure Storage Account}.

    @param azure_blob_service_client: C{azure.storage.blob.BlobServiceClient} object
    @type azure_blob_service_client: BlobServiceClient
    @param container: C{azure.storage.blob.ContainerProperties} object or Python C{str} container name
    @type container: ContainerProperties | str
    @return: C{True} if the container exists, C{False} otherwise
    @rtype: bool
    """
    container_name = get_container_name(container=container)

    result = False
    for container_properties in azure_blob_service_client.list_containers(name_starts_with=container_name):
        if container_properties.name == container_name:
            result = True
            break

    return result


def azure_block_blob_upload(
        file_path: str,
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str],
        blob: Union[BlobProperties, str, None] = None,
        max_concurrency: int = 4,
        logging_enable: bool = False) -> BlobProperties:
    """Upload a block blob into a I{Microsoft Azure Storage Account} I{Container}.

    @param file_path: Local file path
    @type file_path: str
    @param azure_blob_service_client: C{azure.storage.blob.BlobServiceClient} object
    @type azure_blob_service_client: BlobServiceClient
    @param container: C{azure.storage.blob.ContainerProperties} object or Python C{str} container name
    @type container: ContainerProperties | str
    @param blob: C{azure.storage.blob.BlobProperties} or Python C{str} blob name
    @type blob: BlobProperties | str | None
    @param max_concurrency: Maximum number of network connections
    @type max_concurrency: int
    @param logging_enable: Enable logging via the Python C{logger} module
    @type logging_enable: bool
    @return: C{azure.storage.blob.BlobProperties}
    @rtype: BlobProperties
    """
    # Test if the container exists.
    if not azure_container_exists(azure_blob_service_client=azure_blob_service_client, container=container):
        raise Exception(
            'Container ' + repr(get_container_name(container=container)) + ' does not exist in this account.')

    if blob is None:
        blob = os.path.basename(file_path)

    azure_blob_client = azure_blob_service_client.get_blob_client(container=container, blob=blob)

    with open(file=file_path, mode='rb') as binary_io:
        # The blob_type option defaults to azure.storage.blob.BlobType.BlockBlob
        # The retries_total option cannot be used in the upload_blob() call.
        azure_blob_client.upload_blob(
            data=binary_io,
            max_concurrency=max_concurrency,
            logging_enable=logging_enable)

    return azure_blob_client.get_blob_properties()


def azure_block_blob_download(
        file_path: str,
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str],
        blob: Union[BlobProperties, str],
        max_concurrency: int = 4,
        logging_enable: bool = False) -> BlobProperties:
    """Download a block blob from a I{Microsoft Azure Storage Account} I{Container}.

    @param file_path: Local file path
    @type file_path: str
    @param azure_blob_service_client: C{azure.storage.blob.BlobServiceClient} object
    @type azure_blob_service_client: BlobServiceClient
    @param container: C{azure.storage.blob.ContainerProperties} object or Python C{str} container name
    @type container: ContainerProperties | str
    @param blob: C{azure.storage.blob.BlobProperties} or Python C{str} blob name
    @type blob: BlobProperties | str
    @param max_concurrency: Maximum number of network connections
    @type max_concurrency: int
    @param logging_enable: Enable logging via the Python C{logger} module
    @type logging_enable: bool
    @return: C{azure.storage.blob.BlobProperties}
    @rtype: BlobProperties
    """
    # Test if the container exists.
    if not azure_container_exists(azure_blob_service_client=azure_blob_service_client, container=container):
        raise Exception(
            'Container ' + repr(get_container_name(container=container)) + ' does not exist in this account.')

    if not file_path:
        # The Azure Storage Blob Service always uses URL-compliant slash characters as path separators.
        file_path = get_azure_blob_name(blob=blob).split('/')[-1:]

    azure_blob_client = azure_blob_service_client.get_blob_client(container=container, blob=blob)

    with open(file=file_path, mode='wb') as binary_io:
        azure_blob_client.download_blob(
            max_concurrency=max_concurrency,
            logging_enable=logging_enable).readinto(stream=binary_io)

    return azure_blob_client.get_blob_properties()
