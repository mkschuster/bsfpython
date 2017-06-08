#! /usr/bin/env python
#
# BSF Python script to test library functions.
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


import math
import os

import pysam


genome_path = '/data/prod/ngs_resources/gatk_bundle/2.8/b37/human_g1k_v37_decoy.fasta'


def _partition_indices(tile_list, chunk_size):
    return [tile_list[offs:offs + chunk_size] for offs in range(0, len(tile_list), chunk_size)]


def _read_sequence_dict(fasta_path, tiles=0, width=0):

    interval_list = list()
    total_length = 0  # int

    dict_path = os.path.splitext(fasta_path)[0] + '.dict'
    if not os.path.exists(dict_path):
        raise Exception("Picard sequence dictionary {!r} does not exist.".format(dict_path))

    alignment_file = pysam.AlignmentFile(dict_path, 'r')
    # Summarise sequence lengths to get the total length.
    for sq_entry in alignment_file.header['SQ']:
        assert isinstance(sq_entry, dict)
        # print repr(sq_entry)
        total_length += sq_entry['LN']  # int

    if tiles:
        tile_length = float(total_length) / float(tiles)  # float
    elif width:
        tile_length = float(width)  # float
    else:
        # The intervals are just the natural sequence regions.
        # Thus the start coordinate is always 1 and the end coordinate is the sequence length (@SQ SL).
        # 1 2 3 ... 7 8 9
        # start = 1
        # end = 9
        # length = end - start + 1 = 9 - 1 + 1 = 9
        for sq_entry in alignment_file.header['SQ']:
            assert isinstance(sq_entry, dict)
            interval_list.append((sq_entry['SN'], 1, sq_entry['LN']))

        return interval_list

    print "Genome length: {:d} number of tiles: {:d} tile_length: {:,.3f}".format(total_length, tiles, tile_length)

    current_length = 0.0  # float
    sq_list = list()
    for sq_entry in alignment_file.header['SQ']:
        assert isinstance(sq_entry, dict)
        # print "SQ SN:{!r} LN:{:,d}".format(sq_entry['SN'], sq_entry['LN'])
        sq_start = 0.0  # float
        sq_length = float(sq_entry['LN'])
        while sq_start < sq_length:  # float
            # The sequence end is the minimum of the sequence start plus remaining tile length or the sequence length.
            sq_end = min(sq_start + tile_length - current_length, sq_length)
            sq_list.append((sq_entry['SN'], int(math.floor(sq_start + 1.0)), int(math.floor(sq_end))))
            current_length += sq_end - sq_start
            sq_start = sq_end

            if math.floor(current_length) >= math.floor(tile_length):
                # If a tile is complete, append the sequence list to the interval list and reset both list and length.
                interval_list.append(sq_list)
                sq_list = []
                current_length = 0.0

    if len(sq_list):
        interval_list.append(sq_list)

    # For gathering ...

    _partition_indices(tile_list=range(0, len(interval_list)), chunk_size=25)

    return interval_list


intervals = _read_sequence_dict(fasta_path=genome_path, tiles=100)
print "Interval list length {:d}".format(len(intervals))
for index in range(0, len(intervals)):
    print "{:d} {!r}".format(index, intervals[index])
    interval_length = 0
    for sq_region in intervals[index]:
        interval_length += sq_region[2] - sq_region[1] + 1
        print "   {!r} length {:,d}".format(sq_region, sq_region[2] - sq_region[1] + 1)

    print "  Length: {:,d}".format(interval_length)
