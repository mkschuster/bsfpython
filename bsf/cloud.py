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
import errno
import json
import logging
import os
from argparse import ArgumentParser
from typing import IO, Optional, Union

from azure.storage.blob import BlobServiceClient, \
    ContainerClient, ContainerProperties, \
    BlobClient, BlobProperties, BlobBlock, \
    StandardBlobTier

from bsf.standards import Secrets, StandardFilePath


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


def console_blob_download(
        account_name: Optional[str] = None,
        container_name: Optional[str] = None,
        blob_item_list: Optional[list[str]] = None,
        blob_path: Optional[str] = None,
        max_concurrency: Optional[int] = None) -> int:
    """Console function to download from the :emphasis:`Microsoft Azure Storage Blob Service`.

    :param account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str | None
    :param container_name: A Python :py:class:`str` (container name) object.
    :type container_name: str | None
    :param blob_item_list: A Python :py:class:`list` object of Python :py:Class:`str` (blob name) objects.
    :type blob_item_list: list[str] | None
    :param blob_path: Directory path for downloaded blobs.
    :type blob_path: str | None
    :param max_concurrency: Maximum number of network connections.
    :type max_concurrency: int | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not os.path.isdir(blob_path):
        raise Exception(f'The blob directory path {blob_path!r} does not exist.')

    azure_blob_service_client = get_azure_blob_service_client(account_name=account_name)

    if not azure_container_exists(
            azure_blob_service_client=azure_blob_service_client,
            container=container_name):
        raise Exception(f'The Azure Blob Container {container_name!r} does not exist.')

    for blob_item in blob_item_list:
        logging.debug('Blob item: %r', blob_item)

        for suffix in ('', '.md5'):
            file_path = os.path.join(blob_path, blob_item + suffix)

            if os.path.exists(file_path):
                logging.info('File exists already: %r', file_path)
                continue

            blob_properties = azure_block_blob_download(
                azure_blob_service_client=azure_blob_service_client,
                container=container_name,
                blob=blob_item + suffix,
                file_path=file_path,
                max_concurrency=max_concurrency)

            logging.info('Azure Blob name: %r', blob_properties.name)
            logging.info('Azure Blob size: %r', blob_properties.size)
            logging.info('Azure Blob ETag: %r', blob_properties.etag)
            logging.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())

    return 0


def entry_point_blob_download() -> int:
    """Console entry point to download from the :emphasis:`Microsoft Azure Storage Blob Service`.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Microsoft Azure Storage Blob Service download script.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--account',
        required=True,
        help='Microsoft Azure Storage Account name')

    argument_parser.add_argument(
        '--container',
        required=True,
        help='Microsoft Azure Storage Blob Service container name')

    argument_parser.add_argument(
        '--threads',
        default=1,
        help='maximum number of concurrent download threads [1]',
        type=int)

    argument_parser.add_argument(
        '--blob-path',
        default='.',
        help='directory path for storing downloaded blobs [.]')

    argument_parser.add_argument(
        'blob_item',
        nargs='+',
        help='one or more blob items (e.g., archive_file.tar.gz)',
        dest='blob_item_list')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_blob_download(
        account_name=name_space.account,
        container_name=name_space.container,
        blob_item_list=name_space.blob_item_list,
        blob_path=name_space.blob_path,
        max_concurrency=name_space.threads)


def console_blob_download_sequence(
        account_name: Optional[str] = None,
        container_name: Optional[str] = None,
        sequence_item_list: Optional[list[str]] = None,
        sequence_path: Optional[str] = None,
        max_concurrency: Optional[int] = None) -> int:
    """Console function to download sequences in BAM format from the
    :emphasis:`Microsoft Azure Storage Blob Service`.

    :param account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str | None
    :param container_name: A Python :py:class:`str` (container name) object.
    :type container_name: str | None
    :param sequence_item_list: A Python :py:class:`list` object of Python :py:class:`str` (sequence item) objects.
    :type sequence_item_list: list[str] | None
    :param sequence_path: Directory path for downloaded blobs.
    :type sequence_path: str | None
    :param max_concurrency: Maximum number of network connections.
    :type max_concurrency: int | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    if not os.path.isdir(sequence_path):
        raise Exception(f'The sequences directory {sequence_path!r} does not exist.')

    azure_blob_service_client = get_azure_blob_service_client(account_name=account_name)

    if not azure_container_exists(
            azure_blob_service_client=azure_blob_service_client,
            container=container_name):
        raise Exception(f'The Azure Blob Container {container_name!r} does not exist.')

    for sequence_item in sequence_item_list:
        logging.debug('Sequence item: %r', sequence_item)

        # Get the experiment directory name by removing the lane suffix.
        experiment_name = '_'.join(sequence_item.split('_')[:-1])
        logging.debug('Experiment name: %r', experiment_name)

        experiment_directory = os.path.join(sequence_path, experiment_name)

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
                container=container_name,
                blob=experiment_name + '/' + sequence_item + suffix,
                file_path=file_path,
                max_concurrency=max_concurrency)

            logging.info('Azure Blob name: %r', blob_properties.name)
            logging.info('Azure Blob size: %r', blob_properties.size)
            logging.info('Azure Blob ETag: %r', blob_properties.etag)
            logging.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())

    return 0


def entry_point_blob_download_sequence() -> int:
    """Console entry point to download sequences in BAM format from the
    :emphasis:`Microsoft Azure Storage Blob Service`.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Microsoft Azure Storage Blob Service sequence download script.')

    argument_parser.add_argument(
        '--account',
        required=True,
        help='Microsoft Azure Storage Account name')

    argument_parser.add_argument(
        '--container',
        required=True,
        help='Microsoft Azure Storage Blob Service container name')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--threads',
        default=1,
        type=int,
        help='maximum number of concurrent download threads [1]')

    argument_parser.add_argument(
        '--sequence-path',
        default=StandardFilePath.get_sequences(),
        help=f'sequence directory path for storing sequence items [{StandardFilePath.get_sequences()}]')

    argument_parser.add_argument(
        'sequence_item',
        nargs='+',
        help='one or more sequence item names (i.e., experiment_flowcell_lane)',
        dest='sequence_item_list')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_blob_download_sequence(
        account_name=name_space.account,
        container_name=name_space.container,
        sequence_item_list=name_space.sequence_item_list,
        sequence_path=name_space.sequence_path,
        max_concurrency=name_space.threads)


def console_blob_upload(
        account_name: Optional[str] = None,
        container_name: Optional[str] = None,
        standard_blob_tier: Optional[str] = None,
        max_concurrency: Optional[int] = None,
        retain_path: Optional[bool] = None,
        file_path_list: Optional[list[str]] = None) -> int:
    """Console function to upload to the :emphasis:`Microsoft Azure Storage Blob Service`.

    :param account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str | None
    :param container_name: A Python :py:class:`str` (container name) object.
    :type container_name: str | None
    :param standard_blob_tier: A standard blob tier.
    :type standard_blob_tier: str | None
    :param file_path_list: A Python :py:class:`list` object of Python :py:class:`str` (file path) objects.
    :type file_path_list: list[str] | None
    :param max_concurrency: Maximum number of network connections.
    :type max_concurrency: int | None
    :param retain_path: Request retaining the path.
    :type retain_path: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    azure_blob_service_client = get_azure_blob_service_client(account_name=account_name)

    if not azure_container_exists(
            azure_blob_service_client=azure_blob_service_client,
            container=container_name):
        raise Exception(f'The Azure Blob Container {container_name!r} does not exist.')

    for file_path in file_path_list:
        if not os.path.isfile(file_path):
            raise Exception(f'The file path {file_path!r} does not exist.')

        if retain_path:
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
            container=container_name,
            blob=blob_path,
            standard_blob_tier=standard_blob_tier,
            max_concurrency=max_concurrency)

        logging.info('Azure Blob name: %r', blob_properties.name)
        logging.info('Azure Blob size: %r', blob_properties.size)
        logging.info('Azure Blob ETag: %r', blob_properties.etag)
        logging.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())

    return 0


def entry_point_blob_upload() -> int:
    """Console entry point to upload to the :emphasis:`Microsoft Azure Storage Blob Service`.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Microsoft Azure Storage Blob Service upload script.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    argument_parser.add_argument(
        '--account',
        required=True,
        help='Microsoft Azure Storage Account name')

    argument_parser.add_argument(
        '--container',
        required=True,
        help='Microsoft Azure Storage Blob Service container name')

    argument_parser.add_argument(
        '--standard-blob-tier',
        choices=('Archive', 'Cool', 'Hot'),
        help='Microsoft Azure Standard Blob Tier name (i.e., Archive, Cool or Hot)')

    argument_parser.add_argument(
        '--threads',
        default=1,
        type=int,
        help='maximum number of concurrent download threads')

    argument_parser.add_argument(
        '--retain-path',
        action='store_true',
        help='retain the local path in the blob path')

    argument_parser.add_argument(
        'file_path',
        nargs='+',
        help='one or more file paths (e.g., archive_file.tar.gz)',
        dest='file_path_list')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_blob_upload(
        account_name=name_space.account,
        container_name=name_space.container,
        standard_blob_tier=name_space.standard_blob_tier,
        max_concurrency=name_space.threads,
        retain_path=name_space.retain_path,
        file_path_list=name_space.file_path_list)


def console_blob_list(
        account_name: Optional[str] = None,
        container_name: Optional[str] = None,
        blob_prefix: Optional[str] = None,
        blob_name: Optional[str] = None) -> int:
    """Console function to list objects of the :emphasis:`Microsoft Azure Storage Blob Service`.

    - List of :emphasis:`blobs` of a :emphasis:`container`
    - List of :emphasis:`blocks` of a :emphasis:`blob`

    :param account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str | None
    :param container_name: A Python :py:class:`str` (container name) object.
    :type container_name: str | None
    :param blob_prefix: A blob prefix.
    :type blob_prefix: str | None
    :param blob_name: A blob name.
    :type blob_name: str | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    azure_blob_service_client = get_azure_blob_service_client(account_name=account_name)

    if not azure_container_exists(
            azure_blob_service_client=azure_blob_service_client,
            container=container_name):
        raise Exception(f'The Azure Blob Container {container_name!r} does not exists.')

    if blob_name:
        # If a blob name was provided, list all blocks of the blob.

        azure_blob_client = azure_blob_service_client.get_blob_client(container=container_name, blob=blob_name)

        committed_block_list: list[BlobBlock]
        uncommitted_block_list: list[BlobBlock]

        (committed_block_list, uncommitted_block_list) = azure_blob_client.get_block_list(block_list_type='uncommitted')

        print('List of blocks in blob:', blob_name)

        counter = 0
        print('  Committed Blocks:')
        for blob_block in committed_block_list:
            counter += 1
            print('    Blob Block id:', blob_block.id,
                  'size:', blob_block.size,
                  'state:', blob_block.state,
                  'counter:', counter)

        counter = 0
        print('  Uncommitted Blocks:')
        for blob_block in uncommitted_block_list:
            counter += 1
            print('    Blob Block id:', blob_block.id,
                  'size:', blob_block.size,
                  'state:', blob_block.state,
                  'counter:', counter)
    else:
        # If a particular blob name was not provided, list all blobs in the container.

        azure_container_client = azure_blob_service_client.get_container_client(container=container_name)

        print('List of blobs in container:', container_name)

        blob_properties: BlobProperties
        for blob_properties in azure_container_client.list_blobs(
                name_starts_with=blob_prefix,
                include=['snapshots', 'metadata', 'uncommittedblobs', 'copy', 'deleted']):
            print('  Blob name:', blob_properties.name,
                  'etag:', blob_properties.etag,
                  'size:', blob_properties.size,
                  'type:', blob_properties.blob_type,
                  'archive_status:', blob_properties.archive_status)

            if blob_properties.etag is None and not blob_properties.deleted:
                answer = input('Delete this uncommitted blob? [y/N] ')

                if answer == 'y':
                    azure_container_client.delete_blob(blob=blob_properties)

    return 0


def entry_point_blob_list() -> int:
    """Console entry point to list objects of the :emphasis:`Microsoft Azure Storage Blob Service`.

    - List of :emphasis:`blobs` of a :emphasis:`container`
    - List of :emphasis:`blocks` of a :emphasis:`blob`

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Microsoft Azure Storage Blob Service blob or block listing.')

    argument_parser.add_argument(
        '--account',
        required=True,
        help='Microsoft Azure Storage Account name')

    argument_parser.add_argument(
        '--container',
        required=True,
        help='Microsoft Azure Storage Blob Service container name')

    argument_parser.add_argument(
        '--blob-prefix',
        help='blob name prefix')

    argument_parser.add_argument(
        '--blob-name',
        help='blob name, if given list blocks for this blob')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='logging level [WARNING]')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    return console_blob_list(
        account_name=name_space.account,
        container_name=name_space.container,
        blob_prefix=name_space.blob_prefix,
        blob_name=name_space.blob_name)
