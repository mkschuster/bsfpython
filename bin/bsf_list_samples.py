#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-
#
# BSF Python script to list a hierarchy of BSF ProcessedRunFolder, BSF Project and
# BSF Sample objects as CSV file.
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
import re
from argparse import ArgumentParser
from csv import DictWriter

from bsf.ngs import ProcessedRunFolder
from bsf.standards import Configuration, FilePath as StandardsFilePath

argument_parser = ArgumentParser(description='List projects and samples.')

# Only a --full flag, with out a value.
argument_parser.add_argument(
    '--full',
    const=True,
    default=False,
    dest='full',
    help='Full listing with file paths',
    nargs='?',
    required=False)

argument_parser.add_argument(
    '--input',
    dest='input_directory',
    help='Input directory',
    required=True)

argument_parser.add_argument(
    '--output',
    dest='output_file',
    help='Output (CSV) file',
    required=True)

name_space = argument_parser.parse_args()

# Open the output file, create a csv.DictWriter and write a header line.

# For Python2.7, the open() function has to use the binary 'b' flag.
# For Python3, the open() function has to use newline=''.
csv_file = open(file=name_space.output_file, mode='w', newline='')

csv_fields = ['ProcessedRunFolder', 'Project', 'Sample']

if name_space.full:
    csv_fields.extend(['File1', 'Reads1', 'File2', 'Reads2'])

csv_writer = DictWriter(f=csv_file, fieldnames=csv_fields)

csv_writer.writeheader()

# Assemble the input directory.

input_directory = Configuration.get_absolute_path(
    file_path=name_space.input_directory,
    default_path=StandardsFilePath.get_sequences(absolute=True))

prf = ProcessedRunFolder.from_file_path(file_path=input_directory, file_type='Automatic')
for project_name in sorted(prf.project_dict):
    project = prf.project_dict[project_name]
    for sample_name in sorted(project.sample_dict):
        sample = project.sample_dict[sample_name]

        row_dict = {'ProcessedRunFolder': prf.name, 'Project': project.name, 'Sample': sample.name}

        if name_space.full:
            for paired_reads in sample.paired_reads_list:
                if paired_reads.reads_1 is not None:
                    row_dict['File1'] = paired_reads.reads_1.file_path
                    # row_dict['Reads1'] = paired_reads.reads1.name
                    # Deduce the Reads.name from the base name without file extensions.
                    file_name = os.path.basename(paired_reads.reads_1.file_path.rstrip('/ '))
                    # Remove any file-extensions like .bam, .fastq.gz, ...
                    match = re.search(pattern=r'^([^.]*)', string=file_name)
                    if match:
                        file_name = match.group(1)
                    row_dict['Reads1'] = file_name
                else:
                    row_dict['File1'] = str()
                    row_dict['Reads1'] = str()
                if paired_reads.reads_2 is not None:
                    row_dict['File2'] = paired_reads.reads_2.file_path
                    # row_dict['Reads2'] = paired_reads.reads2.name
                    # Deduce the Reads.name from the base name without file extensions.
                    file_name = os.path.basename(paired_reads.reads_2.file_path.rstrip('/ '))
                    # Remove any file-extensions like .bam, .fastq.gz, ...
                    match = re.search(pattern=r'^([^.]*)', string=file_name)
                    if match:
                        file_name = match.group(1)
                    row_dict['Reads2'] = file_name
                else:
                    row_dict['File2'] = str()
                    row_dict['Reads2'] = str()
                # A line for each PairedReads replicate.
                csv_writer.writerow(rowdict=row_dict)
        else:
            # A line for each Sample, containing PairedReads replicates.
            csv_writer.writerow(rowdict=row_dict)

csv_file.close()
