# -*- coding: utf-8 -*-
"""Executables module

A package of classes and methods supporting executable programs and scripts.
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

import bsf.process


class BWA(bsf.process.Executable):
    """Burrows-Wheeler Aligner (C{bsf.executables.BWA}) class.

    Reference: http://bio-bwa.sourceforge.net/
    Usage: bwa mem db_prefix reads.fq [mates.fq]
    """

    def __init__(self, name, analysis):
        """Initialise a C{bsf.executables.BWA} object.

        @param name: Name
        @type name: str
        @param analysis: C{bsf.analysis.Analysis}
        @type analysis: bsf.analysis.Analysis
        @return:
        @rtype:
        """
        super(BWA, self).__init__(name=name, program='bwa', sub_command=bsf.process.Command(program='mem'))

        # The options have to be set for the 'mem' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_configuration(configuration=analysis.configuration, section=section)

        # Set default BWA mem options.

        # None for the moment.


class TopHat(bsf.process.Executable):
    """C{bsf.executables.TopHat} RNA-Seq aligner class.

    Reference: http://tophat.cbcb.umd.edu/manual.html
    Usage: tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]
    Arguments:
    <ebwt_base> Base name of the index to be searched.
    <reads1_1[,...,readsN_1]>
    <[reads1_2,...readsN_2]>
    """

    def __init__(self, name, analysis):
        """Initialise a C{bsf.executables.TopHat} object.

        @param name: Name
        @type name: str
        @param analysis: C{bsf.analysis.Analysis}
        @type analysis: bsf.analysis.Analysis
        @return:
        @rtype:
        """
        super(TopHat, self).__init__(name=name, program='tophat2')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default TopHat options.

        # None for the moment.
