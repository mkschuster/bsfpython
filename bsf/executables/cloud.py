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
"""The :py:mod:`bsf.executables.cloud` module provides classes and functions to transfer files to and from the
:literal:`Microsoft Azure Storage Blob Service`.
"""
import logging
import sys
from argparse import ArgumentParser
from subprocess import Popen
from typing import Optional

from bsf.argument import Argument
from bsf.cloud import get_azure_blob_service_client, azure_block_blob_download, azure_block_blob_upload
from bsf.connector import Connector
from bsf.process import Command, RunnableStep

module_logger = logging.getLogger(name=__name__)


class RunnableStepAzureBlockBlob(RunnableStep):
    """The :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlob` class supports the
     :literal:`Microsoft Azure Storage Blob Service`.

    :ivar account_name: A :literal:`Microsoft Azure Storage Account` name.
    :type account_name: str
    :ivar container_name: A :literal:`Microsoft Azure Storage Blob Service` container name.
    :type container_name: str
    :ivar source_path: A source (local) file path.
    :type source_path: str | None
    :ivar target_path: A target (blob) file path.
    :type target_path: str | None
    :ivar max_concurrency: A maximum number of :literal:`Microsoft Azure Storage Blob Service` network connections.
    :type max_concurrency: int
    :ivar logging_enable: Enable :literal:`Microsoft Azure Storage Blob Service` logging via the
        Python :py:class:`logging.Logger` class.
    :type logging_enable: bool
    """

    def __init__(
            self,
            name: Optional[str],
            program: Optional[str] = None,
            options: Optional[dict[str, list[Argument]]] = None,
            arguments: Optional[list[str]] = None,
            sub_command: Optional[Command] = None,
            stdin: Optional[Connector] = None,
            stdout: Optional[Connector] = None,
            stderr: Optional[Connector] = None,
            dependencies: Optional[list[str]] = None,
            hold: Optional[bool] = None,
            submit: bool = True,
            process_identifier: Optional[str] = None,
            process_name: Optional[str] = None,
            sub_process: Optional[Popen] = None,
            obsolete_file_path_list: Optional[list[str]] = None,
            account_name: str = None,
            container_name: str = None,
            source_path: Optional[str] = None,
            target_path: Optional[str] = None,
            max_concurrency: int = 4,
            logging_enable: bool = False) -> None:
        """Initialise a :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlob` object.

        :param name: A name.
        :type name: str | None
        :param program: A program.
        :type program: str | None
        :param options: A Python :py:class:`dict` object of
            Python :py:class:`str` (:py:attr:`bsf.argument.Argument.key`) key and
            Python :py:class:`list` value objects of :py:class:`bsf.argument.Argument` objects.
        :type options: dict[str, list[Argument]] | None
        :param arguments: A Python :py:class:`list` object of Python :py:class:`str` (program argument) objects.
        :type arguments: list[str] | None
        :param sub_command: A subordinate :py:class:`bsf.process.Command` object.
        :type sub_command: Command | None
        :param stdin: A standard input :literal:`STDIN` :py:class:`bsf.connector.Connector` object.
        :type stdin: Connector | None
        :param stdout: A standard output :literal:`STDOUT` :py:class:`bsf.connector.Connector` object.
        :type stdout: Connector | None
        :param stderr: A standard error :literal:`STDERR` :py:class:`bsf.connector.Connector` object.
        :type stderr: Connector | None
        :param dependencies: A Python :py:class:`list` object of
            Python :py:class:`str` (:py:attr:`bsf.process.Executable.name`) objects
            in the context of :py:class:`bsf.analysis.Stage` dependencies.
        :type dependencies: list[str] | None
        :param hold: Request a hold on job scheduling.
        :type hold: bool | None
        :param submit: Request the submission via the :py:meth:`bsf.analysis.Stage.submit` method.
        :type submit: bool
        :param process_identifier: A process identifier.
        :type process_identifier: str | None
        :param process_name: A process name.
        :type process_name: str | None
        :param sub_process: A :py:class:`subprocess.Popen` object.
        :type sub_process: Popen | None
        :param obsolete_file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (file path) objects
            that can be removed after successfully completing the :py:meth:`bsf.process.RunnableStep.run` method.
        :type obsolete_file_path_list: list[str] | None
        :param account_name: A :literal:`Microsoft Azure Storage Account` name.
        :type account_name: str
        :param container_name: A :literal:`Microsoft Azure Storage Blob Service` container name.
        :type container_name: str
        :param source_path: A source (local) file path.
        :type source_path: str | None
        :param target_path: A target (blob) file path.
        :type target_path: str | None
        :param max_concurrency: A maximum number of
            :literal:`Microsoft Azure Storage Blob Service` network connections.
        :type max_concurrency: int
        :param logging_enable: Enable :literal:`Microsoft Azure Storage Blob Service` logging via the
            Python :py:class:`logging.Logger` class.
        :type logging_enable: bool
        """
        super(RunnableStepAzureBlockBlob, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.account_name = account_name
        self.container_name = container_name
        self.source_path = source_path
        self.target_path = target_path
        self.max_concurrency = max_concurrency
        self.logging_enable = logging_enable

        return

    def __repr__(self) -> str:
        return \
            f'{super(RunnableStepAzureBlockBlob, self).__repr__()[:-1]}, ' \
            f'account_name={self.account_name!r}, ' \
            f'container_name={self.container_name!r}, ' \
            f'source_path={self.source_path!r}, ' \
            f'target_path={self.target_path!r}, ' \
            f'max_concurrency={self.max_concurrency!r}, ' \
            f'logging_enable={self.logging_enable!r})'


class RunnableStepAzureBlockBlobUpload(RunnableStepAzureBlockBlob):
    """The :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlobUpload` class uploads local file paths to
    block blobs.

    :ivar standard_blob_tier: A :literal:`Microsoft Azure Storage Blob Service`
        Python :py:class:`str` standard blob tier enumerated value
        (i.e., :literal:`Archive`, :literal:`Cool`, :literal:`Hot`).
    :type standard_blob_tier: str | None
    """

    def __init__(
            self,
            name: Optional[str],
            program: Optional[str] = None,
            options: Optional[dict[str, list[Argument]]] = None,
            arguments: Optional[list[str]] = None,
            sub_command: Optional[Command] = None,
            stdin: Optional[Connector] = None,
            stdout: Optional[Connector] = None,
            stderr: Optional[Connector] = None,
            dependencies: Optional[list[str]] = None,
            hold: Optional[bool] = None,
            submit: bool = True,
            process_identifier: Optional[str] = None,
            process_name: Optional[str] = None,
            sub_process: Optional[Popen] = None,
            obsolete_file_path_list: Optional[list[str]] = None,
            account_name: str = None,
            container_name: str = None,
            source_path: Optional[str] = None,
            target_path: Optional[str] = None,
            max_concurrency: int = 4,
            logging_enable: bool = False,
            standard_blob_tier: Optional[str] = None) -> None:
        """Initialise a :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlob` object.

        :param name: A name.
        :type name: str | None
        :param program: A program.
        :type program: str | None
        :param options: A Python :py:class:`dict` object of
            Python :py:class:`str` (:py:attr:`bsf.argument.Argument.key`) key and
            Python :py:class:`list` value objects of :py:class:`bsf.argument.Argument` objects.
        :type options: dict[str, list[Argument]] | None
        :param arguments: A Python :py:class:`list` object of Python :py:class:`str` (program argument) objects.
        :type arguments: list[str] | None
        :param sub_command: A subordinate :py:class:`bsf.process.Command` object.
        :type sub_command: Command | None
        :param stdin: A standard input :literal:`STDIN` :py:class:`bsf.connector.Connector` object.
        :type stdin: Connector | None
        :param stdout: A standard output :literal:`STDOUT` :py:class:`bsf.connector.Connector` object.
        :type stdout: Connector | None
        :param stderr: A standard error :literal:`STDERR` :py:class:`bsf.connector.Connector` object.
        :type stderr: Connector | None
        :param dependencies: A Python :py:class:`list` object of
            Python :py:class:`str` (:py:attr:`bsf.process.Executable.name`) objects
            in the context of :py:class:`bsf.analysis.Stage` dependencies.
        :type dependencies: list[str] | None
        :param hold: Request a hold on job scheduling.
        :type hold: bool | None
        :param submit: Request the submission via the :py:meth:`bsf.analysis.Stage.submit` method.
        :type submit: bool
        :param process_identifier: A process identifier.
        :type process_identifier: str | None
        :param process_name: A process name.
        :type process_name: str | None
        :param sub_process: A :py:class:`subprocess.Popen` object.
        :type sub_process: Popen | None
        :param obsolete_file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (file path) objects
            that can be removed after successfully completing the :py:meth:`bsf.process.RunnableStep.run` method.
        :type obsolete_file_path_list: list[str] | None
        :param account_name: A :literal:`Microsoft Azure Storage Account` name.
        :type account_name: str
        :param container_name: A :literal:`Microsoft Azure Storage Blob Service` container name.
        :type container_name: str
        :param source_path: A source (local) file path.
        :type source_path: str | None
        :param target_path: A target (blob) file path.
        :type target_path: str | None
        :param max_concurrency: The maximum number of
            :literal:`Microsoft Azure Storage Blob Service` network connections.
        :type max_concurrency: int
        :param logging_enable: Enable :literal:`Microsoft Azure Storage Blob Service` logging via the
            Python :py:class:`logging.Logger` class.
        :type logging_enable: bool
        :param standard_blob_tier: A :literal:`Microsoft Azure Storage Blob Service`
            Python :py:class:`str` standard blob tier enumerated value
            (i.e., :literal:`Archive`, :literal:`Cool`, :literal:`Hot`).
        :type standard_blob_tier: str | None
        """
        super(RunnableStepAzureBlockBlobUpload, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list,
            account_name=account_name,
            container_name=container_name,
            source_path=source_path,
            target_path=target_path,
            max_concurrency=max_concurrency,
            logging_enable=logging_enable
        )

        self.standard_blob_tier = standard_blob_tier

        return

    def __repr__(self) -> str:
        return \
            f'{super(RunnableStepAzureBlockBlob, self).__repr__()[:-1]}, ' \
            f'account_name={self.account_name!r}, ' \
            f'container_name={self.container_name!r}, ' \
            f'source_path={self.source_path!r}, ' \
            f'target_path={self.target_path!r}, ' \
            f'max_concurrency={self.max_concurrency!r}, ' \
            f'logging_enable={self.logging_enable!r}, ' \
            f'standard_blob_tier={self.standard_blob_tier!r})'

    def run(self) -> None:
        """Run a :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlobUpload` object.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (exception) objects.
        :rtype: list[str] | None
        """
        blob_properties = azure_block_blob_upload(
            file_path=self.source_path,
            azure_blob_service_client=get_azure_blob_service_client(account_name=self.account_name),
            container=self.container_name,
            blob=self.target_path,
            standard_blob_tier=self.standard_blob_tier,
            max_concurrency=self.max_concurrency,
            logging_enable=self.logging_enable)

        module_logger.info('Azure Blob name: %r', blob_properties.name)
        module_logger.info('Azure Blob size: %r', blob_properties.size)
        module_logger.info('Azure Blob ETag: %r', blob_properties.etag)
        module_logger.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())

        return None


class RunnableStepAzureBlockBlobDownload(RunnableStepAzureBlockBlob):
    """The :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlobDownload` class downloads block blobs to
    local file paths.
    """

    def run(self) -> None:
        """Run a :py:class:`bsf.executables.cloud.RunnableStepAzureBlockBlobDownload` object.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (exception) objects.
        :rtype: list[str] | None
        """
        blob_properties = azure_block_blob_download(
            azure_blob_service_client=get_azure_blob_service_client(account_name=self.account_name),
            container=self.container_name,
            blob=self.target_path,
            file_path=self.source_path,
            max_concurrency=self.max_concurrency,
            logging_enable=self.logging_enable)

        module_logger.info('Azure Blob name: %r', blob_properties.name)
        module_logger.info('Azure Blob size: %r', blob_properties.size)
        module_logger.info('Azure Blob ETag: %r', blob_properties.etag)
        module_logger.info('Azure Blob Last Modified: %r', blob_properties.last_modified.isoformat())

        return None


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description='Module console script.')

    argument_parser.add_argument(
        '--logging-level',
        default='WARNING',
        choices=('CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'),
        help='Logging level [WARNING]')

    argument_parser.add_argument(
        '--action',
        choices=('upload', 'download'),
        required=True,
        help='Action (i.e. upload, download)')

    argument_parser.add_argument(
        '--account-name',
        required=True,
        help='Account name')

    argument_parser.add_argument(
        '--container-name',
        required=True,
        help='Container name')

    argument_parser.add_argument(
        '--source-path',
        required=True,
        help='Source (local) file path')

    argument_parser.add_argument(
        '--target-path',
        required=True,
        help='Target (blob) file path')

    argument_parser.add_argument(
        '--maximum-concurrency',
        default=4,
        type=int,
        help='Maximum number of network connections [4]')

    argument_parser.add_argument(
        '--enable-logging',
        action='store_true',
        help='Enable logging [False]',
        dest='logging_enable')

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    if name_space.action == 'upload':
        runnable_step = RunnableStepAzureBlockBlobUpload(
            name='block_blob_upload',
            account_name=name_space.account_name,
            container_name=name_space.container_name,
            source_path=name_space.source_path,
            target_path=name_space.target_path,
            max_concurrency=name_space.maximum_concurrency,
            logging_enable=name_space.logging_enable)

        runnable_step.run()

    if name_space.action == 'download':
        runnable_step = RunnableStepAzureBlockBlobDownload(
            name='block_blob_download',
            account_name=name_space.account_name,
            container_name=name_space.container_name,
            source_path=name_space.source_path,
            target_path=name_space.target_path,
            max_concurrency=name_space.maximum_concurrency,
            logging_enable=name_space.logging_enable)

        runnable_step.run()

    return 0


if __name__ == '__main__':
    sys.exit(main())
