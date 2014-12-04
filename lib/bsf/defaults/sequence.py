"""bsf.defaults.sequence

A package to centralise sequence information.
"""

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


# Illumina TruSeq DNA v1/v2/LT and RNA v1/v2/LT and ChIP Sample Prep Kits.

illumina_truseq_dna_rna_indices = \
    [
        'NNNNNN',
        'ATCACG',  # 1
        'CGATGT',  # 2
        'TTAGGC',  # 3
        'TGACCA',  # 4
        'ACAGTG',  # 5
        'GCCAAT',  # 6
        'CAGATC',  # 7
        'ACTTGA',  # 8
        'GATCAG',  # 9
        'TAGCTT',  # 10
        'GGCTAC',  # 11
        'CTTGTA',  # 12
        'AGTCAA',  # 13
        'AGTTCC',  # 14
        'ATGTCA',  # 15
        'CCGTCC',  # 16
        'CTAGAG',  # 17 Reserved, but reverse complement of small RNA kit.
        'GTCCGC',  # 18
        'GTGAAA',  # 19
        'GTGGCC',  # 20
        'GTTTCG',  # 21
        'CGTACG',  # 22
        'GAGTGG',  # 23
        'GGTAGC',  # 24 Reserved, but reverse complement of small RNA kit.
        'ACTGAT',  # 25
        'ATGAGC',  # 26 Reserved, but reverse complement of small RNA kit.
        'ATTCCT'   # 27
    ]


# Illumina oligonucleotide sequences for TruSeq Small RNA Sample Prep Kits.

illumina_truseq_small_rna_indices = \
    [
        'NNNNNN',
        'CGTGAT',  # RPI1
        'ACATCG',  # RPI2
        'GCCTAA',  # RPI3
        'TGGTCA',  # RPI4
        'CACTGT',  # RPI5
        'ATTGGC',  # RPI6
        'GATCTG',  # RPI7
        'TCAAGT',  # RPI8
        'CTGATC',  # RPI9
        'AAGCTA',  # RPI10
        'GTAGCC',  # RPI11
        'TACAAG',  # RPI12
        'TTGACT',  # RPI13
        'GGAACT',  # RPI14
        'TGACAT',  # RPI15
        'GGACGG',  # RPI16
        'CTCTAC',  # RPI17
        'GCGGAC',  # RPI18
        'TTTCAC',  # RPI19
        'GGCCAC',  # RPI20
        'CGAAAC',  # RPI21
        'CGTACG',  # RPI22
        'CCACTC',  # RPI23
        'GCTACC',  # RPI24
        'ATCAGT',  # RPI25
        'GCTCAT',  # RPI26
        'AGGAAT',  # RPI27
        'CTTTTG',  # RPI28
        'TAGTTG',  # RPI29
        'CCGGTG',  # RPI30
        'ATCGTG',  # RPI31
        'TGAGTG',  # RPI32
        'CGCCTG',  # RPI33
        'GCCATG',  # RPI34
        'AAAATG',  # RPI35
        'TGTTGG',  # RPI36
        'ATTCCG',  # RPI37
        'AGCTAG',  # RPI38
        'GTATAG',  # RPI39
        'TCTGAG',  # RPI40
        'GTCGTC',  # RPI41
        'CGATTA',  # RPI42
        'GCTGTA',  # RPI43
        'ATTATA',  # RPI44
        'GAATGA',  # RPI45
        'TCGGGA',  # RPI46
        'CTTCGA',  # RPI47
        'TGCCGA'   # RPI48
    ]
