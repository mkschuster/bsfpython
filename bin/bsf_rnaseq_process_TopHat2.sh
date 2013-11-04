#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to convert TopHat2 output files into
# binary formats suitable for UCSC Genome Browser track hubs.
#
#
# Copyright 2013 Michael Schuster
# CeMM - Research Center for Molecular Medicine of the Austrian Academy of Sciences
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

if [ -z "${LANG}" ]; then
    declare -x LANG='C'
else
    LANG='C'
fi

if test "$#" -lt '2'; then
    echo "Error: bsf_rnaseq_process_tophat2.sh Too few arguments." 1>&2 \
    || exit 1
    echo "Usage: bsf_rnaseq_process_tophat2.sh <directory> <chromosome_sizes>" 1>&2 \
    || exit 1
    exit 1
fi

declare directory="${1}"
declare chromosome_sizes="${2}"

# Function to convert UCSC BED to BigBed files.

function process_bed () {

    declare directory="${1}"
    declare chromosome_sizes="${2}"
    declare prefix="${3}"

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
        || exit 1

    sort -k1,1 -k2,2n "${directory}/${prefix}.txt" \
        > "${directory}/${prefix}_sorted.txt" \
        || exit 1

    bedToBigBed \
        "${directory}/${prefix}_sorted.txt" \
        "${chromosome_sizes}" \
        "${directory}/${prefix}.bb" \
        || exit 1

    rm "${directory}/${prefix}.txt" || exit 1
    rm "${directory}/${prefix}_sorted.txt" || exit 1

    return 0

}

declare -a prefixes

prefixes[0]='deletions'
prefixes[1]='insertions'
prefixes[2]='junctions'

for prefix in "${prefixes[@]}"; do
    process_bed "${directory}" "${chromosome_sizes}" "${prefix}"
done

# Index the sorted BAM file accepted_hits.bam.

samtools index "${directory}/accepted_hits.bam" || exit 1

exit 0
