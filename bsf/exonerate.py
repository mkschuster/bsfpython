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
"""The :py:mod:`bsf.exonerate` module provides classes supporting Guy Slater's Exonerate alignment tool.
"""
import re
from typing import Optional, TypeVar

VULGARType = TypeVar(name='VULGARType', bound='VULGAR')


class VULGAR(object):
    """The :py:class:`bsf.exonerate.VULGAR` class models the Verbose Useful Labelled Gapped Alignment Report (VULGAR)
    as defined by Guy Slater's `Exonerate <https://www.ebi.ac.uk/~guy/exonerate/>`_ generic alignment tool.

    :ivar q_name: A :literal:`query` name.
    :type q_name: str | None
    :ivar q_start: A :literal:`query` start.
    :type q_start: str | None
    :ivar q_end: A :literal:`query` end.
    :type q_end: str | None
    :ivar q_strand: A :literal:`query` strand.
    :type q_strand: str | None
    :ivar t_name: A :literal:`target` name.
    :type t_name: str | None
    :ivar t_start: A :literal:`target` start.
    :type t_start: str | None
    :ivar t_end: A :literal:`target` end.
    :type t_end: str | None
    :ivar t_strand: A :literal:`target` strand.
    :type t_strand: str | None
    :ivar score: A score.
    :type score: str | None
    :ivar triplet_list: A Python :py:class:`list` object of
        (:literal:`operation`, :literal:`query_length`, :literal:`target_length`)
        Python :py:class:`tuple` objects.
    :type triplet_list: list[(str, str, str)] | None
    """

    @classmethod
    def from_vulgar_str(cls, vulgar_str: str) -> VULGARType:
        """Create a new Verbose Useful Labelled Gapped Alignment Report (:py:class:`bsf.exonerate.VULGAR`) object
        from a Python :py:class:`str` (VULGAR) object.

        :param vulgar_str: A VULGAR string.
        :type vulgar_str: str
        :return: A :py:class:`bsf.exonerate.VULGAR` object.
        :rtype: VULGAR
        """
        vulgar = cls()

        (
            vulgar.q_name, vulgar.q_start, vulgar.q_end, vulgar.q_strand,
            vulgar.t_name, vulgar.t_start, vulgar.t_end, vulgar.t_strand,
            vulgar.score, triplet_str
        ) = vulgar_str.split(' ', 9)

        triplet_list = triplet_str.split(' ')
        # Further, split into triplets by means of a list comprehension.
        vulgar.triplet_list = [triplet_list[i: i + 3] for i in range(0, len(triplet_list), 3)]

        return vulgar

    def __init__(
            self,
            q_name: Optional[str] = None,
            q_start: Optional[str] = None,
            q_end: Optional[str] = None,
            q_strand: Optional[str] = None,
            t_name: Optional[str] = None,
            t_start: Optional[str] = None,
            t_end: Optional[str] = None,
            t_strand: Optional[str] = None,
            score: Optional[str] = None,
            triplet_list: Optional[list[tuple[str, str, str]]] = None) -> None:
        """Initialise a new :py:class:`bsf.exonerate.VULGAR` object.

        :param q_name: A :literal:`query` name.
        :type q_name: str | None
        :param q_start: A :literal:`query` start.
        :type q_start: str | None
        :param q_end: A :literal:`query` end.
        :type q_end: str | None
        :param q_strand: A :literal:`query` strand.
        :type q_strand: str | None
        :param t_name: A :literal:`target` name.
        :type t_name: str | None
        :param t_start: A :literal:`target` start.
        :type t_start: str | None
        :param t_end: A :literal:`target` end.
        :type t_end: str | None
        :param t_strand: A :literal:`target` strand.
        :type t_strand: str | None
        :param score: A score.
        :type score: str | None
        :param triplet_list: A Python :py:class:`list` object of
            (:literal:`operation`, :literal:`query_length`, :literal:`target_length`)
            Python :py:class:`tuple` objects.
        :type triplet_list: list[(str, str, str)] | None
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

    def t_start_natural(self) -> int:
        """Return a natural target start coordinate.

        :return: A natural start coordinate.
        :rtype: int
        """
        if self.t_strand == '+':
            return int(self.t_start)
        elif self.t_strand == '-':
            return int(self.t_end)
        else:
            # If not defined return target start
            return int(self.t_start)

    def t_end_natural(self) -> int:
        """Return a natural end coordinate.

        :return: A natural end coordinate.
        :rtype: int
        """
        if self.t_strand == '+':
            return int(self.t_end)
        elif self.t_strand == '-':
            return int(self.t_start)
        else:
            # If not defined return target end
            return int(self.t_end)


def parse_alignment_file(file_path: str) -> list[VULGAR]:
    """Parse an alignment file.

    :param file_path: Alignment file path.
    :type file_path: str
    :return: A Python :py:class:`list` object of :py:class:`bsf.exonerate.VULGAR` objects.
    :rtype: list[VULGAR]
    """
    vulgar_pattern = re.compile(pattern=r'^vulgar: (.*)')

    vulgar_list: list[VULGAR] = list()

    with open(file=file_path, mode='rt') as input_text_io:
        for line_str in input_text_io:
            vulgar_match = re.search(pattern=vulgar_pattern, string=line_str)
            if vulgar_match:
                vulgar = VULGAR.from_vulgar_str(vulgar_str=vulgar_match.group(1))
                vulgar_list.append(vulgar)

    return vulgar_list
