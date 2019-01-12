#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# BSF Python script to extract sequences form Exonerate alignments.
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

from __future__ import print_function

import re

import Bio.SeqIO
import Bio.SeqRecord
import numpy

import bsf.exonerate

vulgar_pattern = re.compile(pattern='^vulgar: (.*)')
identifier_pattern = re.compile(pattern='^\\w+\\|(\\w+)\\|([^ ]+)')

# sample_name = 'Plasmid'
sample_name = 'cDNA'

record_dict = Bio.SeqIO.index_db(
    index_filename='sequences_{}_pass_hq_2d.idx'.format(sample_name),
    filenames=['sequences_{}_pass_hq_2d.fasta'.format(sample_name)],
    format='fasta')

print('Number of query sequence records:', len(record_dict))

# TODO: This could possibly run Exonerate as a sub process.
seq_io_dict = dict()

# TODO: To analyse these pairwise alignments it would be good to create a profile from all pairwise alignments.
# First run: Establish the maximum target (reference) region that appears in the alignments.
# Create a two-dimensional array (reference position vs base or gap) and use the minimum target position to
# scale the array to 0.
# Second run: Parse each individual position in the alignment, compare to reference and fill in matrix.
# Evaluate the matrix and ideally visualise the results.

vulgar_list = bsf.exonerate.parse_alignment_file(
    file_path='alignments_{}_pass_hq_2d_sequence_db.txt'.format(sample_name))

# Process the list of VULGAR objects.

min_target_vulgar = min(vulgar_list, key=lambda x: x.t_start_natural())
max_target_vulgar = max(vulgar_list, key=lambda x: x.t_end_natural())

min_target_position = min_target_vulgar.t_start_natural()
max_target_position = max_target_vulgar.t_end_natural()

array_length = max_target_position - min_target_position + 1

# Use these values to initialise a two-dimensional array i.e. sequence length and A, C, G, T and - letters.
profile_array = numpy.zeros(shape=(array_length, 5), dtype=numpy.int32)

print('Profile Array shape:', profile_array.shape)

for vulgar in vulgar_list:
    # Fetch the corresponding SeqRecord object.
    q_record = record_dict[vulgar.q_name]
    """ @type q_record: Bio.SeqRecord.SeqRecord """

    # Check if the sequence was in forward or reverse orientation and reverse complement if necessary.
    orientation = int(vulgar.q_strand + '1') * int(vulgar.t_strand + '1')
    if orientation < 0:
        # BioPython type ambiguity.
        # In method Bio.SeqRecord.SeqRecord.reverse_complement(), the id parameter can be of type basestring or bool
        # at the same time.
        q_reverse_complement = q_record.reverse_complement(
            id=q_record.id + '_rc',
            name=q_record.name + '_rc',
            description=True)
        q_record = q_reverse_complement

    identifier_match = re.match(pattern=identifier_pattern, string=vulgar.t_name)
    if identifier_match:
        identifier = identifier_match.group(1)
    else:
        identifier = vulgar.t_name

    if identifier in seq_io_dict:
        seq_io = seq_io_dict[identifier]
    else:
        seq_io = open("test_aligned_{}_to_{}.fasta".format(sample_name, identifier), 'wt')
        seq_io_dict[identifier] = seq_io

    Bio.SeqIO.write(sequences=q_record, handle=seq_io, format='fasta')

    # TODO: Populate the matrix at this stage.

    print('Query:', vulgar.q_name)
    for exonerate_tuple in vulgar.triplet_list:
        print('Operation:', ' '.join(exonerate_tuple))

# Clean up stage: close all Bio.SeqIO file handles.

for seq_io in seq_io_dict.itervalues():
    seq_io.close()
