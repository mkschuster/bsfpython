#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# BSF GNU Bourne-Again (Bash) script to convert TopHat2 output files into
# binary formats suitable for UCSC Genome Browser track hubs.
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

if [[ -z "${LANG}" ]];
then
    declare -x LANG='C';
else
    LANG='C';
fi

if test "$#" -lt '3';
then
    echo "Error: $(basename "${0}")  Too few arguments." 1>&2 \
    || exit 1;
    echo "Usage: $(basename "${0}")  <directory> <genome_fasta> <genome_sizes>" 1>&2 \
    || exit 1;
    exit 1;
fi

declare directory="${1}";
declare genome_fasta="${2}";
declare genome_sizes="${3}";

declare -x TMPDIR="{$directory}";

# Function to convert UCSC BED to BigBed files.

function process_bed () {

    declare directory="${1}";
    declare prefix="${2}";
    declare genome_sizes="${3}";

    if [[ -f "${directory}/${prefix}.bb" ]];
    then
        return 0;
    fi

    # The UCSC bedToBigBed utility requires that track lines are stripped out. Sigh!

    # The BED format specification requires that the fifth field
    # is an integer score between 0 and 1000.
    # MACS does not seem to know or care about that.
    # Use a Perl one-liner to correct this.
    # The bedToBigBed does neither understand a 'stdin' option,
    # nor is it capable of reading from standard input.
    # Convert via a temporary *.txt file.

    # TODO: The score would need a complete re-scaling between 0 and 1000. Set to 0 in the meantime.

    grep --extended-regexp --invert-match '^track|^browser' \
        "${directory}/${prefix}.bed" \
        | perl -e 'while (<>) { chomp; my @fields = split q{ }; $fields[4] = 0; print join(qq{\t}, @fields), qq{\n}; }' \
        > "${directory}/${prefix}.txt" \
        || exit 1;

    sort -k1,1 -k2,2n "${directory}/${prefix}.txt" \
        > "${directory}/${prefix}_sorted.txt" \
        || exit 1;

    bedToBigBed \
        "${directory}/${prefix}_sorted.txt" \
        "${genome_sizes}" \
        "${directory}/${prefix}.bb" \
        || exit 1;

    rm "${directory}/${prefix}.txt" || exit 1;
    rm "${directory}/${prefix}_sorted.txt" || exit 1;

    return 0;

}

declare -a prefixes;

prefixes[0]='deletions';
prefixes[1]='insertions';
prefixes[2]='junctions';

for prefix in "${prefixes[@]}";
do
    process_bed "${directory}" "${prefix}" "${genome_sizes}";
done

# Re-sort the accepted_hits.bam file with the newer samtools.
samtools sort \
    -l 9 \
    -m 2G \
    -o "${directory}/accepted_hits_sorted.bam" \
    --reference "${genome_fasta}" \
    --threads 3 \
    "${directory}/accepted_hits.bam" \
    || exit 1;

mv "${directory}/accepted_hits_sorted.bam" "${directory}/accepted_hits.bam" || exit 1;

# Index the sorted BAM file accepted_hits.bam.

samtools index \
    -@ 3 \
    "${directory}/accepted_hits.bam" \
    || exit 1;

# Use RseQC bam.wig.py to create a UCSC Wig file.

# TODO: This should set the --strand rule option.

bam2wig.py \
    --input-file "${directory}/accepted_hits.bam" \
    --chromSize "${genome_sizes}" \
    --out-prefix "${directory}/accepted_hits" \
    || exit 1;

function process_wig {

    declare directory="${1}";
    declare prefix="${2}";
    declare genome_sizes="${3}";

    if [[ -f "${directory}/accepted_hits${prefix}.wig" ]];
    then
        wigToBigWig \
            -clip \
            "${directory}/accepted_hits${prefix}.wig" \
            "${genome_sizes}" \
            "${directory}/accepted_hits${prefix}.bw" \
            || exit 1;

        rm "${directory}/accepted_hits${prefix}.wig" \
            || exit 1;
    fi

    return 0;

}

# Process by RseQC names for unstranded versus stranded data.

prefixes[0]='';
prefixes[1]='_Forward';
prefixes[2]='_Reverse';

for prefix in "${prefixes[@]}";
do
    process_wig "${directory}" "${prefix}" "${genome_sizes}";
done

echo "All done." || exit 1;
