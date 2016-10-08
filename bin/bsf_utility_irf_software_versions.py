#! /usr/bin/env python
#
# BSF Python utility script to report software versions for Illumina Run Folders.
#
#
# Copyright 2013 - 2016 Michael K. Schuster
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

from argparse import ArgumentParser
import os
import re
import xml

from bsf.illumina import RunFolder
from bsf.annotation import AnnotationSheet


argument_parser = ArgumentParser(
    description='List software version of Illumina Run Folders.')

argument_parser.add_argument(
    '--debug',
    help='debug level',
    required=False,
    type=int)

argument_parser.add_argument(
    '--directory',
    help='directory of illumina run folders',
    required=True,
    type=str)

argument_parser.add_argument(
    '--output-file',
    dest='output_file',
    help='output (*.csv) file path',
    required=True,
    type=str)

argument_parser.add_argument(
    '--ascending',
    action='store_true',
    help='sort flow cells in ascending order rather than in descending by default')

name_space = argument_parser.parse_args()

file_name_list = os.listdir(name_space.directory)
file_name_list.sort(cmp=lambda x, y: cmp(x, y))

if not name_space.ascending:
    file_name_list.reverse()

field_names = ['run_identifier', 'application_name', 'application_version', 'rta_version']
annotation_sheet = AnnotationSheet(file_path=name_space.output_file, header=True, field_names=field_names)

# Illumina Run Folders obey a pattern and additionally directories with just Sequence Analysis Viewer (SAV)
# information should also be allowed.

irf_pattern = re.compile(pattern=r'^[0-9]{6,6}_.*(?:_sav)?$')

for file_name in file_name_list:
    if name_space.debug:
        print 'File name: {!r}'.format(file_name)
    # Process just entries that obey the Illumina Run Folder pattern.
    match = re.search(pattern=irf_pattern, string=file_name)
    if not match:
        print 'No match: {!r}'.format(file_name)
        continue
    file_path = os.path.join(name_space.directory, file_name)
    if not os.path.exists(os.path.join(file_path, 'runParameters.xml')):
        print 'Directory {!r} not an Illumina Run Folder'.format(file_name)
        continue
    # Temporarily catch IOError and xml.etree.ElementTree.ParseError exceptions
    # that result from a broken FhGFS file system.
    try:
        irf = RunFolder.from_file_path(file_path=file_path)
    except IOError as exception:
        if name_space.debug:
            print "\t".join((file_name, '?', '?', '?'))
        continue
    except xml.etree.ElementTree.ParseError:
        if name_space.debug:
            print "\t".join((file_name, '?', '?', '?'))
        continue
    except:
        print 'Exception in run folder {!r}'.format(file_name)
        raise

    root_node = irf.run_parameters.element_tree.getroot()
    annotation_sheet.row_dicts.append({
        'run_identifier': irf.run_parameters.get_run_identifier,
        'application_name': irf.run_parameters.get_application_name,
        'application_version': irf.run_parameters.get_application_version,
        'rta_version': irf.run_parameters.get_real_time_analysis_version,
    })

annotation_sheet.to_file_path()
