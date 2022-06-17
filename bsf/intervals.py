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
"""The :py:mod:`bsf.intervals` module provides classes modelling (genome) intervals.
"""
import math
import os
from argparse import ArgumentParser
from typing import Dict, List

from pysam import AlignmentFile


class Interval(object):
    """Intervals with sequence region name, start and end coordinates.

    The :py:class:`bsf.intervals.Interval` class reflects Picard-style interval lists and its coordinates are 1-based.

    :ivar name: Sequence region name
    :type name: str
    :ivar start: Start coordinate
    :type start: int
    :ivar end: End coordinate
    :type end: int
    """

    def __init__(self, name, start, end):
        """Initialise a :py:class:`bsf.intervals.Interval` object.

        :param name: Sequence region name
        :type name: str
        :param start: Start coordinate
        :type start: int
        :param end: End coordinate
        :type end: int
        """
        self.name = name
        self.start = start
        self.end = end

        return

    def __bool__(self):
        """Test for truth.

        Since Picard-style interval lists are 1-based, all components need to be defined.

        :return: :py:const:`True` if start and end are defined, :py:const:`False` otherwise.
        :rtype: bool
        """
        return bool(self.name) and bool(self.start) and bool(self.end)

    def __len__(self):
        """Calculate the length.

        :return: The length.
        :rtype: int
        """
        return self.end - self.start + 1

    def __str__(self):
        """Convert into a string representation

        :return: A string representation.
        :rtype: str
        """
        return f'Interval(name={self.name!r}, start={self.start!r}, end={self.end!r})'
        # return ':'.join((self.name, str(self.start), str(self.end)))

    def to_gatk_interval(self):
        """Convert into a GATK interval string (i.e., name:start-end).

        :return: A GATK interval string.
        :rtype: str
        """
        return self.name + ':' + str(self.start) + '-' + str(self.end)


class Container(object):
    """Container for :py:class:`bsf.intervals.Interval` class that keeps a running sum.

    :ivar interval_list: Python :py:class:`list` object of :py:class:`bsf.intervals.Interval` objects
    :type interval_list: list[Interval]
    :ivar sum: Sum
    :type sum: int
    """

    def __init__(self):
        """Initialise a :py:class:`bsf.intervals.Container` object.
        """
        self.interval_list = list()

        self.sum = 0

        return

    def __len__(self):
        """Calculate the length of a :py:class:`bsf.intervals.Container`,
        which is the number of :py:class:`bsf.intervals.Interval` objects.

        :return: The length.
        :rtype: int
        """
        return len(self.interval_list)

    def append(self, interval):
        """Append an :py:class:`bsf.intervals.Interval` object.

        :param interval: A :py:class:`bsf.intervals.Interval` object
        :type interval: Interval
        """
        self.interval_list.append(interval)
        self.sum += len(interval)

        return

    def __str__(self):
        """Get a string representation.

        :return: A string representation.
        :rtype: str
        """
        return 'Container(sum={:d}, interval_list={!r})'.format(self.sum, self.interval_list)


def get_interval_tiles(interval_path=None, tile_number=None, tile_width=None, acgt=None, natural=None, packed=None):
    """Create tiles on the basis of a Picard-style interval list.

    The Picard :literal:`ScatterIntervalsByNs` interval list is a special case where only :literal:`ACGTmer`
    intervals need considering.

    Partition a list into sub-lists, which sums do not exceed a maximum using a :emphasis:`First Fit Decreasing`
    algorithm.

    See also:
        - `First-Fit algorithm <https://en.wikipedia.org/wiki/Bin_packing_problem#First-fit_algorithm>`_
        - `Python implementations of packing algorithm <https://stackoverflow.com/questions/7392143/python-implementations-of-packing-algorithm>`_

    :param interval_path: A Picard :literal:`ScatterIntervalsByNs` interval list file path.
    :type interval_path: str | None
    :param tile_number: Number of tiles
    :type tile_number: int | None
    :param tile_width: Width of tile
    :type tile_width: int | None
    :param acgt: Picard :literal:`ScatterIntervalsByNs` mode reading just :literal:`ACGTmer` intervals
    :type acgt: bool | None
    :param natural: Tile on (natural) sequence regions (i.e., SAM :literal:`@SQ` entries)
    :type natural: bool | None
    :param packed: Pack :py:class:`bsf.intervals.Interval` objects into :py:class:`bsf.intervals.Container` objects
    :type packed: bool | None
    :return: A Python :py:class:`list` object of :py:class:`bsf.intervals.Container` objects.
    :rtype: list[Container]
    """
    if not os.path.exists(interval_path):
        raise Exception('Interval file ' + repr(interval_path) + ' does not exists.')

    container_list: List[Container] = list()
    interval_list: List[Interval] = list()
    total_length: int = 0

    with open(file=interval_path, mode='rt') as input_file:
        for line_str in input_file:
            if line_str.startswith('@'):
                continue
            field_list = line_str.strip().split('\t')
            if acgt and field_list[4] != 'ACGTmer':
                continue
            interval = Interval(name=field_list[0], start=int(field_list[1]), end=int(field_list[2]))
            interval_list.append(interval)
            total_length += len(interval)

    if natural:
        # If neither width nor number was requested, return one Container with all Interval objects.
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
        # Return a list of a single Container with an empty Interval i.e., sequence region, start and end coordinates.
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
    """Create tiles on the basis of a Picard :literal:`CreateSequenceDictionary` sequence dictionary.

    The tiles are created on the basis of a Picard sequence dictionary accompanying the genome FASTA file.
    Sequence regions can serve as natural tiles, alternatively a number of tiles or a tile width can be
    requested. If neither tiles nor width are requested or both are 0, no tiling is attempted and a
    single :py:class:`bsf.intervals.Container` object with an empty :py:class:`bsf.intervals.Interval`
    is put onto the Python :py:class:`list` object.

    :param dictionary_path: Picard :literal:`CreateSequenceDictionary` sequence dictionary file path
    :type dictionary_path: str
    :param tile_number: Number of tiles
    :type tile_number: int | None
    :param tile_width: Width of tile
    :type tile_width: int | None
    :param natural: Tile on (natural) sequence regions (i.e., SAM :literal:`@SQ` entries)
    :type natural: bool | None
    :return: A Python :py:class:`list` object of :py:class:`bsf.intervals.Container` objects.
    :rtype: list[Container]
    """
    if not os.path.exists(dictionary_path):
        raise Exception('Picard sequence dictionary ' + repr(dictionary_path) + ' does not exist.')

    sq_entry: Dict
    container_list: List[Container] = list()
    total_length: int = 0

    alignment_file = AlignmentFile(dictionary_path, 'rt')

    # Summarise sequence lengths to get the total length.
    for sq_entry in alignment_file.header['SQ']:
        total_length += int(sq_entry['LN'])

    if natural:
        # The intervals are just the natural sequence regions.
        # Thus, the start coordinate is always 1 and the end coordinate is the sequence length (@SQ LN).
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
        # Return a list of a single Container with an empty Interval i.e., sequence region, start and end coordinates.
        container = Container()
        container.append(Interval(name='', start=0, end=0))
        container_list.append(container)

        return container_list

    current_length: float = 0.0
    container = Container()
    for sq_entry in alignment_file.header['SQ']:
        sq_start: float = 0.0
        sq_length: float = float(sq_entry['LN'])
        while sq_start < sq_length:  # float
            # The sequence end is the minimum of the sequence start plus remaining tile length or
            # the sequence length.
            sq_end: float = min(sq_start + tile_length - current_length, sq_length)
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
    argument_parser = ArgumentParser(
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
        help='Picard interval list file path',
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
        '--acgt',
        action='store_true',
        help='Consider only ACGTmer intervals defined by Picard ScatterIntervalsByNs')

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
            acgt=name_space.acgt,
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
