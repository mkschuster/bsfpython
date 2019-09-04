#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""Microsoft Azure cloud module

A package of classes and functions to transfer files to and from the
Microsoft Azure cloud service.
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

from __future__ import print_function

import argparse

import bsf.cloud
import bsf.process


class RunnableStepAzureBlockBlob(bsf.process.RunnableStep):
    """The C{bsf.executables.cloud.RunnableStepAzureBlockBlob} class supports the Azure C{BlockBlobService}.

    Attributes:
    @ivar account_name: Microsoft Azure storage account name
    @type account_name: str
    @ivar container_name: Microsoft Azure storage container name
    @type container_name: str
    @ivar source_path: Source (local) file path
    @type source_path: str | None
    @ivar target_path: Target (blob) file path
    @type target_path: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            account_name=None,
            container_name=None,
            source_path=None,
            target_path=None):
        """Initialise a C{bsf.executables.cloud.RunnableStepAzureBlockBlob} object.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} or C{unicode} (program argument) objects
        @type arguments: list[str | unicode] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: bsf.connector.Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: bsf.connector.Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: bsf.connector.Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: subprocess.Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str | unicode] | None
        @param account_name: Microsoft Azure storage account name
        @type account_name: str
        @param container_name: Microsoft Azure storage container name
        @type container_name: str
        @param source_path: Source (local) file path
        @type source_path: str | None
        @param target_path: Target (blob) file path
        @type target_path: str | None
        @return:
        @rtype:
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

        return


class RunnableStepAzureBlockBlobUpload(RunnableStepAzureBlockBlob):
    """The C{bsf.executables.cloud.RunnableStepAzureBlockBlobUpload} class uploads file paths as block blobs.

    Attributes:
    """

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.executables.cloud.RunnableStepAzureBlockBlobUpload}.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return 0 on success
        @rtype: int
        """
        resource_properties = bsf.cloud.azure_block_blob_upload(
            block_blob_service=bsf.cloud.get_azure_block_blob_service(account_name=self.account_name),
            container_name=self.container_name,
            file_path=self.source_path,
            blob_name=self.target_path)

        print('Azure BlockBlob name:', self.target_path)
        print('Azure BlockBlob ETag:', resource_properties.etag)
        print('Azure BlockBlob Last Modified:', resource_properties.last_modified.isoformat())

        return 0


class RunnableStepAzureBlockBlobDownload(RunnableStepAzureBlockBlob):
    """The C{bsf.executables.cloud.RunnableStepAzureBlockBlobDownload} class downloads block blobs to file paths.

    Attributes:
    """

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.executables.cloud.RunnableStepAzureBlockBlobDownload}.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return 0 on success
        @rtype: int
        """
        block_blob = bsf.cloud.azure_block_blob_download(
            block_blob_service=bsf.cloud.get_azure_block_blob_service(account_name=self.account_name),
            container_name=self.container_name,
            file_path=self.source_path,
            blob_name=self.target_path)

        print('Azure BlockBlob name:', block_blob.name)
        print('Azure BlockBlob ETag:', block_blob.properties.etag)
        print('Azure BlockBlob Last Modified:', block_blob.properties.last_modified.isoformat())

        return 0


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--debug',
        help='Debug level',
        required=False,
        type=int)

    argument_parser.add_argument(
        '--action',
        dest='action',
        help='Action (i.e. upload, download)',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--account-name',
        dest='account_name',
        help='Account name',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--container-name',
        dest='container_name',
        help='Container name',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--source-path',
        dest='source_path',
        help='Source (local) file path',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--target-path',
        dest='target_path',
        help='Target (blob) file path',
        required=True,
        type=str)

    name_space = argument_parser.parse_args()

    if name_space.action == 'upload':
        runnable_step = RunnableStepAzureBlockBlobUpload(
            name='blob_upload',
            account_name=name_space.account_name,
            container_name=name_space.container_name,
            source_path=name_space.source_path,
            target_path=name_space.target_path)

        runnable_step.run(debug=name_space.debug)

    if name_space.action == 'download':
        runnable_step = RunnableStepAzureBlockBlobDownload(
            name='blob_upload',
            account_name=name_space.account_name,
            container_name=name_space.container_name,
            source_path=name_space.source_path,
            target_path=name_space.target_path)

        runnable_step.run(debug=name_space.debug)