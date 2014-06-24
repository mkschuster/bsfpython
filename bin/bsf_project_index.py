#! /usr/bin/env python
#
# BSF Python script to create a project index HTML document.
#
#
# Copyright 2013 Michael K. Schuster
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

import argparse
import inspect
import os.path
import re
from stat import *
import string

from Bio.BSF import Default
import Bio.BSF.Defaults.web as web


def scan_directory(report_dict, directory_root, directory_path=None):

    """
    Scan a directory recursively for *_report.html files and add them to a Python dict of
    directory_path key data and Python list of report_type value data.
    :param report_dict: Python dict of directory_path key data and Python list of report_type value data
    :type report_dict: dict
    :param directory_root: Directory root
    :type directory_root: str, unicode
    :param directory_path: Directory path
    :type directory_path: str, unicode
    :return: Nothing
    :rtype: None
    """

    if not directory_path:
        directory_path = '.'

    directory_absolute = os.path.join(directory_root, directory_path)

    for file_name in os.listdir(directory_absolute):
        file_path = os.path.join(directory_absolute, file_name)
        mode = os.stat(file_path).st_mode
        if S_ISDIR(mode):
            scan_directory(report_dict=report_dict,
                           directory_root=directory_root,
                           directory_path=os.path.join(directory_path, file_name))
        elif S_ISREG(mode):
            match = re.search(pattern=r'^(.*)_report.html$', string=file_name)
            if match:
                report_type = match.group(1)
                if not directory_path in report_dict:
                    report_dict[directory_path] = list()
                report_dict[directory_path].append(report_type)


def scan_projects(project_name):

    """
    Scan a public_html project directory for full project directory names to match the given prefix.
    Resolve something like BSA_0001 to BSA_0001_ASP14_40d0aa70bf854ac99599caf5a52a9aa3
    :param project_name: Project name or prefix
    :type project_name: str, unicode
    """

    directory_path = os.path.join(Default.absolute_public_html(), 'projects')

    for file_name in os.listdir(directory_path):
        match = re.search(pattern=r'^{}'.format(project_name), string=file_name)
        if match:
            return file_name


parser = argparse.ArgumentParser(description='Create index.html documents in BSF public_html/project folders.')

parser.add_argument('--project', required=True, help='Project identifier')

args = parser.parse_args()

project_name = str(args.project)
project_directory = str(args.project)

if not os.path.isabs(project_directory):

    # If the project is not an absolute path,
    # prepend the absolute public html directory.
    # TODO: This does not deal with sub-directories i.e. public_html/projects correctly.

    project_directory = os.path.join(Default.absolute_public_html(), 'projects', project_directory)

    if not os.path.isdir(project_directory):

        # If the absolute directory path still does not exist,
        # scan the projects directory for an entry that begins with the project name.

        project_name = scan_projects(project_name=args.project)

        if not project_name:
            raise Exception('Cannot locate project directory for project {!r}.'.format(args.project))

        project_directory = os.path.join(Default.absolute_public_html(), 'projects', project_name)


report_dict = dict()

scan_directory(report_dict=report_dict, directory_root=project_directory, directory_path='.')

# Remove the random string from the project_directory name if it exists
# BSA_0001_ASP14_40d0aa70bf854ac99599caf5a52a9aa3

components = project_name.split('_')
if re.search(pattern=r'^[0-9a-f]{32,32}$', string=components[-1]):
    project_name = string.join(components[:-1], '_')

# Now evaluate the report_dictionary.

# TODO: Open an index.html file in the project-specific directory.
# TODO: Check if it already exists and eventually ask for confirmation to over-write ...

output = str()

output += web.html_header(title='Project {} Overview'.format(project_name),
                          source=inspect.getfile(inspect.currentframe()))

keys = report_dict.keys()
keys.sort(cmp=lambda x, y: cmp(x, y))

for key in keys:

    # TODO: Unfortunately, fastqc_report.html documents exist both in the top-level directory,
    # as well as directly from FastQC in sample-specific sub-directories.

    report_list = report_dict[key]

    output += '<p>{}</p>\n'.format(key)
    output += '<ul>\n'

    for report in report_list:
        output += '<li>{}</li>\n'.format(report)

    output += '</ul>\n'

output += web.html_footer()

print output
