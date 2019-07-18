# -*- coding: utf-8 -*-
"""Intervals module

A package of classes and methods modelling (genome) intervals.
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
import math
import os

import pysam


class Interval(object):
    """Intervals with sequence region name, start and end coordinates.

    Attributes:
    @ivar name: Sequence region name
    @type name: str | unicode
    @ivar start: Start coordinate
    @type start: int
    @ivar end: End coordinate
    @type end: int
    """

    def __init__(self, name, start, end):
        """Initialise an C{bsf.intervals.Interval} object.

        @param name: Sequence region name
        @type name: str | unicode
        @param start: Start coordinate
        @type start: int
        @param end: End coordinate
        @type end: int
        """
        self.name = name
        self.start = start
        self.end = end

        return

    def __len__(self):
        """Calculate the length.

        @return: Length
        @rtype: int
        """
        return self.end - self.start + 1

    def __str__(self):
        return 'Interval(name={!r}, start={!r}, end={!r})'.format(self.name, self.start, self.end)
        # return ':'.join((self.name, str(self.start), str(self.end)))

    def to_gatk_interval(self):
        """Convert into a GATK interval string (i.e. name:start-end).

        @return: GATK interval string
        @rtype: str | unicode
        """
        return self.name + ':' + str(self.start) + '-' + str(self.end)


class Container(object):
    """Container for C{bsf.intervals.Interval} objects that keeps a running sum.

    Attributes:
    @ivar interval_list: Python C{list} of C{Interval} objects
    @type interval_list: list[bsf.intervals.Interval]
    @ivar sum: Sum
    @type sum: int
    """

    def __init__(self):
        """Initialise a C{bsf.intervals.Container} object.

        @return:
        @rtype:
        """
        self.interval_list = list()

        self.sum = 0

        return

    def __len__(self):
        """Calculate the length of a C{Container}, which is the number of C{Interval} objects.

        @return: Length
        @rtype: int
        """
        return len(self.interval_list)

    def append(self, interval):
        """Append an C{bsf.intervals.Interval} object.

        @param interval: C{bsf.intervals.Interval}
        @type interval: bsf.intervals.Interval
        @return:
        @rtype:
        """
        self.interval_list.append(interval)
        self.sum += len(interval)

        return

    def __str__(self):
        """ Printable representation """
        return 'Container(sum={:d}, interval_list={!r})'.format(self.sum, self.interval_list)


def get_interval_tiles(interval_path=None, tile_number=None, tile_width=None, natural=None, packed=None):
    """Create interval tiles on the basis of a Picard ScatterIntervalsByNs interval list.

    Partition a list into sub-lists whose sums do not exceed a maximum using a First Fit Decreasing algorithm.
    @see: https://en.wikipedia.org/wiki/Bin_packing_problem#First-fit_algorithm
    @see: https://stackoverflow.com/questions/7392143/python-implementations-of-packing-algorithm
    @param interval_path: Picard ScatterIntervalsByNs interval list file path
    @type interval_path: str | unicode | None
    @param tile_number: Number of tiles
    @type tile_number: int | None
    @param tile_width: Width of tile
    @type tile_width: int | None
    @param natural: Tile on (natural) sequence regions (i.e. SAM @SQ entries)
    @type natural: bool | None
    @param packed: Pack C{Interval} objects into C{Container} objects
    @type packed: bool | None
    @return: Python C{list} of C{bsf.intervals.Container} objects
    @rtype: list[bsf.intervals.Container]
    """
    if not os.path.exists(interval_path):
        raise Exception('Interval file ' + repr(interval_path) + ' does not exists.')

    container_list = list()
    """ @type container_list: list[Container] """

    interval_list = list()
    """ @type interval_list: list[Interval] """

    total_length = 0
    """ @type total_length: int """

    with open(file=interval_path, mode='rt') as input_file:
        for line_str in input_file:
            if line_str.startswith('@'):
                continue
            field_list = line_str.strip().split('\t')
            if field_list[4] == 'ACGTmer':
                interval = Interval(name=field_list[0], start=int(field_list[1]), end=int(field_list[2]))
                interval_list.append(interval)
                total_length += len(interval)

    if natural:
        # If neither width nor number was requested, return a Container with all Interval objects.
        container = Container()
        container_list.append(container)

        for interval in interval_list:
            container.append(interval=interval)

        return container_list

    if tile_number is not None and tile_number > 0:
        tile_length = float(total_length) / float(tile_number)
    elif tile_width is not None and tile_width > 0:
        tile_length = float(tile_width)
    else:
        # Do not tile at all, if neither a number of tiles nor a tile width was provided.
        # Return a list of a single Container with an empty Interval i.e. sequence region, start and end coordinates.
        container = Container()
        container.append(interval=Interval(name='', start=0, end=0))
        container_list.append(container)

        return container_list

    if packed:
        # Pack the interval list using a First Fit Decreasing algorithm.
        # Sort interval tuples by length in descending order.
        for interval in sorted(interval_list, key=lambda item: len(item), reverse=True):
            # Try to fit the item into a bin.
            for container in container_list:
                if container.sum + len(interval) <= tile_length:
                    container.append(interval=interval)
                    break
            else:
                container = Container()
                container_list.append(container)
                container.append(interval=interval)
    else:
        container = Container()
        for interval in interval_list:
            if len(container) > 0 and container.sum + len(interval) > tile_length:
                container_list.append(container)
                container = Container()

            container.append(interval=interval)

        container_list.append(container)

    return container_list


def get_genome_tiles(dictionary_path, tile_number=None, tile_width=None, natural=None):
    """Create genome tiles on the basis of a Picard CreateSequenceDictionary sequence dictionary.

    The tiles are created on the basis of a Picard sequence dictionary accompanying the genome FASTA file.
    Sequence regions can serve as natural tiles, alternatively a number of tiles or a tile width can be
    requested. If neither tiles nor width are requested or both are 0, no tiling is attempted and a
    single C{bsf.intervals.Container} with an empty C{bsf.intervals.Interval} is put onto the Python C{list}.
    @param dictionary_path: Picard CreateSequenceDictionary sequence dictionary file path
    @type dictionary_path: str | unicode
    @param tile_number: Number of tiles
    @type tile_number: int | None
    @param tile_width: Width of tile
    @type tile_width: int | None
    @param natural: Tile on (natural) sequence regions (i.e. SAM @SQ entries)
    @type natural: bool | None
    @return: Python C{list} of C{bsf.intervals.Container} objects
    @rtype: list[bsf.intervals.Container]
    """
    if not os.path.exists(dictionary_path):
        raise Exception('Picard sequence dictionary ' + repr(dictionary_path) + ' does not exist.')

    container_list = list()
    """ @type container_list: list[Container] """

    total_length = 0
    """ @type total_length: int """

    alignment_file = pysam.AlignmentFile(dictionary_path, 'rt')

    # Summarise sequence lengths to get the total length.
    for sq_entry in alignment_file.header['SQ']:
        """ @type sq_entry: dict """
        total_length += int(sq_entry['LN'])

    if natural:
        # The intervals are just the natural sequence regions.
        # Thus the start coordinate is always 1 and the end coordinate is the sequence length (@SQ SL).
        # 1 2 3 ... 7 8 9
        # start = 1
        # end = 9
        # length = end - start + 1 = 9 - 1 + 1 = 9
        for sq_entry in alignment_file.header['SQ']:
            container = Container()
            container.append(Interval(name=str(sq_entry['SN']), start=1, end=int(sq_entry['LN'])))
            container_list.append(container)

        return container_list

    if tile_number is not None and tile_number > 0:
        tile_length = float(total_length) / float(tile_number)
    elif tile_width is not None and tile_width > 0:
        tile_length = float(tile_width)
    else:
        # Do not tile at all, if neither a number of tiles nor a tile width was provided.
        # Return a list of a single Container with an empty Interval i.e. sequence region, start and end coordinates.
        container = Container()
        container.append(Interval(name='', start=0, end=0))
        container_list.append(container)

        return container_list

    current_length = 0.0
    """ @type current_length: float """
    container = Container()
    for sq_entry in alignment_file.header['SQ']:
        sq_start = 0.0
        """ @type sq_start: float """
        sq_length = float(sq_entry['LN'])
        """ @type sq_length: float """
        while sq_start < sq_length:  # float
            # The sequence end is the minimum of the sequence start plus remaining tile length or
            # the sequence length.
            sq_end = min(sq_start + tile_length - current_length, sq_length)
            """ @type seq_end: float """
            container.append(
                Interval(
                    name=str(sq_entry['SN']),
                    start=int(math.floor(sq_start + 1.0)),
                    end=int(math.floor(sq_end))))
            current_length += sq_end - sq_start
            sq_start = sq_end

            if math.floor(current_length) >= math.floor(tile_length):
                # If a tile is complete, append the Container to the container list and reset both
                # Container and length.
                container_list.append(container)
                container = Container()
                current_length = 0.0

    if container.sum:
        container_list.append(container)

    return container_list


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--dictionary-path',
        dest='dictionary_path',
        help='Picard CreateSequenceDictionary sequence dictionary file path',
        required=False,
        type=str)

    argument_parser.add_argument(
        '--interval-path',
        dest='interval_path',
        help='Picard ScatterIntervalsByNs interval list file path',
        required=False,
        type=str)

    argument_parser.add_argument(
        '--tile-number',
        default=100,
        dest='tile_number',
        help='Number of tiles',
        required=False,
        type=int)

    argument_parser.add_argument(
        '--packed',
        action='store_true',
        help='Pack Interval objects into Container objects')

    name_space = argument_parser.parse_args()

    if name_space.interval_path:
        print('Interval tiles:')
        _container_list = get_interval_tiles(
            interval_path=name_space.interval_path,
            tile_number=name_space.tile_number,
            packed=name_space.packed)
        for _container in _container_list:
            print('  Container(sum={:,})'.format(_container.sum))
            for _interval in _container.interval_list:
                print('    ' + str(_interval))
        print('Number of Container objects:', len(_container_list))

    if name_space.dictionary_path:
        print('Genome tiles:')
        _container_list = get_genome_tiles(
            dictionary_path=name_space.dictionary_path,
            tile_number=name_space.tile_number)
        for _container in _container_list:
            print('  Container(sum={:,})'.format(_container.sum))
            for _interval in _container.interval_list:
                print('    ' + str(_interval))
        print('Number of Container objects:', len(_container_list))
