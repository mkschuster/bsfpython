# -*- coding: utf-8 -*-
"""Exonerate module.

A package of classes and methods supporting Guy Slater's Exonerate alignment tool.
"""
#  Copyright 2013 - 2021 Michael K. Schuster
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
import re

vulgar_pattern = re.compile(pattern='^vulgar: (.*)')


class VULGAR(object):
    """This class models the Verbose Useful Labelled Gapped Alignment Report (C{bsf.exonerate.VULGAR}) as defined by
    Guy Slater's Exonerate generic alignment tool.
    https://www.ebi.ac.uk/~guy/exonerate/

    Attributes:
    @ivar q_name: Query name
    @type q_name: str | None
    @ivar q_start: Query start
    @type q_start: str | None
    @ivar q_end: Query end
    @type q_end: str | None
    @ivar q_strand: Query strand
    @type q_strand: str | None
    @ivar t_name: Target name
    @type t_name: str | None
    @ivar t_start: Target start
    @type t_start: str | None
    @ivar t_end: Target end
    @type t_end: str | None
    @ivar t_strand: Target strand
    @type t_strand: str | None
    @ivar score: Score
    @type score: str | None
    @ivar triplet_list: List of (I{operation}, I{query_length}, I{target_length}) tuples
    @type triplet_list: list[(str, str, str)] | None
    """

    @classmethod
    def from_vulgar_str(cls, vulgar_str=None):
        """Create a new Verbose Useful Labelled Gapped Alignment Report (C{bsf.exonerate.VULGAR}) object
        from a Python C{str} (VULGAR) object.

        @param vulgar_str: VULGAR string
        @type vulgar_str: str
        @return: C{bsf.exonerate.VULGAR} object
        @rtype: VULGAR
        """

        self = cls()

        (
            self.q_name, self.q_start, self.q_end, self.q_strand,
            self.t_name, self.t_start, self.t_end, self.t_strand,
            self.score, triplet_str
        ) = vulgar_str.split(' ', 9)

        triplet_list = triplet_str.split(' ')
        # Further split into triplets by means of a list comprehension.
        self.triplet_list = [triplet_list[i: i + 3] for i in range(0, len(triplet_list), 3)]

        return self

    def __init__(
            self,
            q_name=None, q_start=None, q_end=None, q_strand=None,
            t_name=None, t_start=None, t_end=None, t_strand=None,
            score=None, triplet_list=None):
        """Initialise a new C{bsf.exonerate.VULGAR} object.

        @param q_name: Query name
        @type q_name: str | None
        @param q_start: Query start
        @type q_start: str | None
        @param q_end: Query end
        @type q_end: str | None
        @param q_strand: Query strand
        @type q_strand: str | None
        @param t_name: Target name
        @type t_name: str | None
        @param t_start: Target start
        @type t_start: str | None
        @param t_end: Target end
        @type t_end: str | None
        @param t_strand: Target strand
        @type t_strand: str | None
        @param score: Score
        @type score: str | None
        @param triplet_list: List of (I{operation}, I{query_length}, I{target_length}) tuples
        @type triplet_list: list[(str, str, str)] | None
        """

        super(VULGAR, self).__init__()

        self.q_name = q_name
        self.q_start = q_start
        self.q_end = q_end
        self.q_strand = q_strand
        self.t_name = t_name
        self.t_start = t_start
        self.t_end = t_end
        self.t_strand = t_strand
        self.score = score
        self.triplet_list = triplet_list

        return

    def t_start_natural(self):
        """Return the natural target start coordinate.

        @return: Natural start coordinate
        @rtype: int
        """
        if self.t_strand == '+':
            return int(self.t_start)
        elif self.t_strand == '-':
            return int(self.t_end)
        else:
            # If not defined return target start
            return int(self.t_start)

    def t_end_natural(self):
        """Return the natural end coordinate.

        @return: Natural end coordinate
        @rtype: int
        """
        if self.t_strand == '+':
            return int(self.t_end)
        elif self.t_strand == '-':
            return int(self.t_start)
        else:
            # If not defined return target end
            return int(self.t_end)


def parse_alignment_file(file_path):
    """Parse an Exonerate alignment file.

    @param file_path: Alignment file path
    @type file_path: str
    @return: Python C{list} of C{bsf.exonerate.VULGAR} objects
    @rtype: list[VULGAR]
    """
    vulgar_list = list()
    """ @type vulgar_list: list[VULGAR] """

    with open(file=file_path, mode='rt') as input_file:
        for line_str in input_file:
            vulgar_match = re.search(pattern=vulgar_pattern, string=line_str)
            if vulgar_match:
                vulgar = VULGAR.from_vulgar_str(vulgar_str=vulgar_match.group(1))
                vulgar_list.append(vulgar)

    return vulgar_list
