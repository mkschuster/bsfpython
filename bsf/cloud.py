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
"""The :py:mod:`bsf.cloud` module provides functions to interact with cloud services.
"""
import json
import os
from typing import IO, Union

from azure.storage.blob import BlobServiceClient, \
    ContainerClient, ContainerProperties, \
    BlobClient, BlobProperties, \
    StandardBlobTier

from bsf.standards import Secrets


def get_azure_secrets_dict(account_name: str) -> dict[str, str]:
    """Get Microsoft Azure secrets from a JSON configuration file.

    :param account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str
    :return: A Python :py:class:`dict` object with :literal:`Microsoft Azure Storage Account` secrets.
    :rtype: dict[str, str]
    """
    file_path = Secrets.get_azure_file_path()

    if file_path and os.path.exists(file_path):
        with open(file=file_path, mode='rt') as input_text_io:
            secrets_dict: dict[str, dict[str, str]] = json.load(fp=input_text_io)

        if account_name not in secrets_dict:
            raise Exception(
                f'The Microsoft Azure Storage Account name {account_name!r} '
                f'is missing from secrets file {file_path!r}.')

        return secrets_dict[account_name]
    else:
        raise Exception('Could not get a Microsoft Azure secrets JSON file.')


def get_azure_blob_name(blob: Union[BlobProperties, str]) -> str:
    """Get an Azure blob name from a :py:class:`azure.storage.blob.BlobProperties` or a Python :py:class:`str` object.

    :param blob: A :py:class:`azure.storage.blob.BlobProperties` or Python :py:class:`str` object.
    :type blob: BlobProperties | str
    :return: A blob name.
    :rtype: str
    """
    if isinstance(blob, BlobProperties):
        return blob.name
    elif isinstance(blob, str):
        return blob
    else:
        raise Exception("The 'blob' option has to be of type 'BlobProperties' or 'str'.")


def get_azure_container_name(container: Union[ContainerProperties, str]) -> str:
    """Get a container name from a :py:class:`azure.storage.blob.ContainerProperties` or a
    Python :py:class:`str` object.

    :param container: A :py:class:`azure.storage.blob.ContainerProperties` or Python :py:class:`str` object.
    :type container: ContainerProperties | str
    :return: A container name.
    :rtype: str
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
    """Get an :py:class:`azure.storage.blob.BlobServiceClient` object.

    For uploading large files, the default block size of 4 * 1024 * 1024 (4 MiB) is not big enough,
    since the number of blocks in a blob is limited to 50,000, which allows for 200,000 MiB
    or 195 GiB blobs. Changing the maximum block size to 100 MiB allows for 4882 GiB or 4 TiB blobs.

    :param account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str
    :param max_block_size: A maximum block size, which defaults to 100 MiB, which is also the maximum Azure supports.
    :type max_block_size: int
    :param retries_total: A number of retries.
    :type retries_total: int
    :param logging_enable: Enable logging via the Python :py:class:`logging.Logger` class.
    :type logging_enable: bool
    :return: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :rtype: BlobServiceClient
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
    """Get an :py:class:`azure.storage.blob.ContainerClient` object.

    :param azure_blob_service_client: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :type azure_blob_service_client: BlobServiceClient
    :param container: A :py:class:`azure.storage.blob.ContainerProperties` object or
        Python :py:class:`str` (container name) object.
    :type container: ContainerProperties | str
    :return: A :py:class:`azure.storage.blob.ContainerClient` object.
    :rtype: ContainerClient
    """
    return azure_blob_service_client.get_container_client(container=container)


def get_azure_blob_client(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str],
        blob: Union[BlobProperties, str]) -> BlobClient:
    """Get an :py:class:`azure.storage.blob.BlobClient` object.

    :param azure_blob_service_client: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :type azure_blob_service_client: BlobServiceClient
    :param container: A :py:class:`azure.storage.blob.ContainerProperties` object or
        Python :py:class:`str` (container name) object.
    :type container: ContainerProperties | str
    :param blob: A :py:class:`azure.storage.blob.BlobProperties` or Python :py:class:`str` (blob name) object.
    :type blob: BlobProperties | str
    :return: A :py:class:`azure.storage.blob.BlobClient` object.
    :rtype: BlobClient
    """
    return azure_blob_service_client.get_blob_client(container=container, blob=blob)


def azure_container_exists(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str]) -> bool:
    """Check if a container exists in a :literal:`Microsoft Azure Storage Account`.

    :param azure_blob_service_client: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :type azure_blob_service_client: BlobServiceClient
    :param container: A :py:class:`azure.storage.blob.ContainerProperties` object or
        Python :py:class:`str` (container name) object.
    :type container: ContainerProperties | str
    :return: :py:const:`True` if the container exists, :py:const:`False` otherwise.
    :rtype: bool
    """
    container_name = get_azure_container_name(container=container)

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
        standard_blob_tier: Union[StandardBlobTier, str, None] = None,
        max_concurrency: int = 4,
        logging_enable: bool = False) -> BlobProperties:
    """Upload a block blob into a :literal:`Microsoft Azure Storage Account` :literal:`Container`.

    :param file_path: A local file path.
    :type file_path: str
    :param azure_blob_service_client: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :type azure_blob_service_client: BlobServiceClient
    :param container: A :py:class:`azure.storage.blob.ContainerProperties` object or
        Python :py:class:`str` (container name) object.
    :type container: ContainerProperties | str
    :param blob: A :py:class:`azure.storage.blob.BlobProperties` object or
        Python :py:class:`str` (blob name) object.
    :type blob: BlobProperties | str | None
    :param standard_blob_tier: A :py:class:`azure.storage.blob.StandardBlobTier` object or
        Python :py:class:`str` standard blob tier enumerated value
        (i.e., :literal:`Archive`, :literal:`Cool`, :literal:`Hot`).
    :type standard_blob_tier: StandardBlobTier | str | None
    :param max_concurrency: Maximum number of network connections.
    :type max_concurrency: int
    :param logging_enable: Enable logging via the Python :py:class:`logging.Logger` class.
    :type logging_enable: bool
    :return: A :py:class:`azure.storage.blob.BlobProperties` object.
    :rtype: BlobProperties
    """
    # Test if the container exists.
    if not azure_container_exists(azure_blob_service_client=azure_blob_service_client, container=container):
        raise Exception(
            f'An Azure Blob Container {get_azure_container_name(container=container)!r} '
            f'does not exist in this account.')

    if blob is None:
        blob = os.path.basename(file_path)

    azure_blob_client = azure_blob_service_client.get_blob_client(container=container, blob=blob)

    # If the standard_blob_tier option was passed in as str object,
    # it needs converting into a StandardBlobTier object unless it was an empty str.
    if isinstance(standard_blob_tier, str):
        if standard_blob_tier:
            standard_blob_tier = StandardBlobTier(value=standard_blob_tier.title())
        else:
            standard_blob_tier = None

    with open(file=file_path, mode='rb') as input_binary_io:
        # The blob_type option defaults to azure.storage.blob.BlobType.BlockBlob
        # The retries_total option cannot be used in the upload_blob() call.
        azure_blob_client.upload_blob(
            data=input_binary_io,
            standard_blob_tier=standard_blob_tier,
            max_concurrency=max_concurrency,
            logging_enable=logging_enable)

    return azure_blob_client.get_blob_properties()


def azure_block_blob_download_io(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str],
        blob: Union[BlobProperties, str],
        file_io: IO,
        max_concurrency: int = 4,
        logging_enable: bool = False) -> BlobProperties:
    """Download a block blob from a :literal:`Microsoft Azure Storage Account` :literal:`Container`.

    :param azure_blob_service_client: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :type azure_blob_service_client: BlobServiceClient
    :param container: A :py:class:`azure.storage.blob.ContainerProperties` object or
        Python :py:class:`str` (container name) object.
    :type container: ContainerProperties | str
    :param blob: A :py:class:`azure.storage.blob.BlobProperties` object or
        Python :py:class:`str` (blob name) object.
    :type blob: BlobProperties | str
    :param file_io: A Python :py:class:`IO` object.
    :type file_io: IO
    :param max_concurrency: A maximum number of network connections.
    :type max_concurrency: int
    :param logging_enable: Enable logging via the Python :py:class:`logging.Logger` class.
    :type logging_enable: bool
    :return: A :py:class:`azure.storage.blob.BlobProperties` object.
    :rtype: BlobProperties
    """
    # Test if the container exists.
    if not azure_container_exists(azure_blob_service_client=azure_blob_service_client, container=container):
        raise Exception(
            f'An Azure Blob Container {get_azure_container_name(container=container)!r} '
            f'does not exist in this account.')

    azure_blob_client = azure_blob_service_client.get_blob_client(container=container, blob=blob)

    azure_blob_client.download_blob(
        max_concurrency=max_concurrency,
        logging_enable=logging_enable).readinto(stream=file_io)

    return azure_blob_client.get_blob_properties()


def azure_block_blob_download(
        azure_blob_service_client: BlobServiceClient,
        container: Union[ContainerProperties, str],
        blob: Union[BlobProperties, str],
        file_path: str = None,
        max_concurrency: int = 4,
        logging_enable: bool = False) -> BlobProperties:
    """Download a block blob from a :literal:`Microsoft Azure Storage Account` :literal:`Container`.

    :param azure_blob_service_client: A :py:class:`azure.storage.blob.BlobServiceClient` object.
    :type azure_blob_service_client: BlobServiceClient
    :param container: A :py:class:`azure.storage.blob.ContainerProperties` object or
        Python :py:class:`str` (container name) object.
    :type container: ContainerProperties | str
    :param blob: A :py:class:`azure.storage.blob.BlobProperties` object or
        Python :py:class:`str` (blob name) object.
    :type blob: BlobProperties | str
    :param file_path: A local file path.
    :type file_path: str | None
    :param max_concurrency: A maximum number of network connections.
    :type max_concurrency: int
    :param logging_enable: Enable logging via the Python :py:class:`logging.Logger` class.
    :type logging_enable: bool
    :return: A :py:class:`azure.storage.blob.BlobProperties` object.
    :rtype: BlobProperties
    """
    # Test if the container exists.
    if not azure_container_exists(azure_blob_service_client=azure_blob_service_client, container=container):
        raise Exception(
            f'An Azure Blob Container {get_azure_container_name(container=container)!r} '
            f'does not exist in this account.')

    if not file_path:
        # The Azure Storage Blob Service always uses URL-compliant slash characters as path separators.
        file_path = get_azure_blob_name(blob=blob).split('/')[-1:]

    azure_blob_client = azure_blob_service_client.get_blob_client(container=container, blob=blob)

    with open(file=file_path, mode='wb') as output_binary_io:
        azure_blob_client.download_blob(
            max_concurrency=max_concurrency,
            logging_enable=logging_enable).readinto(stream=output_binary_io)

    return azure_blob_client.get_blob_properties()
