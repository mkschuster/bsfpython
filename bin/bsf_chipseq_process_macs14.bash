#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to convert MACS14 BED files into
# BigBed format and WIG files into BigWig format.
#
#
# Copyright 2013 - 2018 Michael K. Schuster
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

if [ -z "${LANG}" ]; then
    declare -x LANG='C'
else
    LANG='C'
fi

if test "$#" -lt '2'; then
    echo "Error: $(basename ${0}) Too few arguments." 1>&2 \
    || exit 1
    echo "Usage: $(basename ${0}) <prefix> <chromosome_sizes>" 1>&2 \
    || exit 1
    exit 1
fi

declare prefix="${1}"
declare chromosome_sizes="${2}"

# The UCSC wigToBigWig utility requires that browser and track lines are stripped out.
# Sigh!

# NAME_treat_afterfiting_all.wig.gz

echo "$(date) wigToBigWig: ${prefix}_treat_afterfiting_all.wig.gz" || exit 1
echo "$(date) wigToBigWig: ${prefix}_treat_afterfiting_all.wig.gz" 1>&2 || exit 1

zgrep --extended-regexp --invert-match '^track|^browser' \
    "${prefix}_MACS_wiggle/treat/${prefix}_treat_afterfiting_all.wig.gz" \
    | wigToBigWig -clip \
    stdin \
    "${chromosome_sizes}" \
    "${prefix}_MACS_wiggle/treat/${prefix}_treat_afterfiting_all.bw" \
    || exit 1

# NAME_control_afterfiting_all.wig.gz

echo "$(date) wigToBigWig: ${prefix}_control_afterfiting_all.wig.gz" || exit 1
echo "$(date) wigToBigWig: ${prefix}_control_afterfiting_all.wig.gz" 1>&2 || exit 1

zgrep --extended-regexp --invert-match '^track|^browser' \
    "${prefix}_MACS_wiggle/control/${prefix}_control_afterfiting_all.wig.gz" \
    | wigToBigWig -clip \
    stdin \
    "${chromosome_sizes}" \
    "${prefix}_MACS_wiggle/control/${prefix}_control_afterfiting_all.bw" \
    || exit 1

# The UCSC bedToBigBed utility requires that track lines are stripped out. Sigh!

# MACs seems to generate BED files that are sorted already,
# so that the following line is not needed.
# sort -k1,1 -k2,2n unsorted.bed > sorted.bed
#
# The BED format specification requires that the fifth field
# is an integer score between 0 and 1000.
# MACS does not seem to know or care about that.
# Use a Perl one-liner to correct this.
# The bedToBigBed does neither understand a 'stdin' option,
# nor is it capable of reading from standard input.
# Convert via a temporary *.txt file.

# TODO: The score would need a complete re-scaling between 0 and 1000. Set to 0 in the meantime.

# NAME_peaks.bed

echo "$(date) bedToBigBed: ${prefix}_peaks.bed" || exit 1
echo "$(date) bedToBigBed: ${prefix}_peaks.bed" 1>&2 || exit 1

grep --extended-regexp --invert-match '^track|^browser' \
    "${prefix}_peaks.bed" \
    | perl -e 'while (<>) { chomp; my @fields = split q{ }; $fields[4] = 0; print join(qq{\t}, @fields), qq{\n}; }' \
    > "${prefix}_peaks.txt" \
    || exit 1

bedToBigBed \
    "${prefix}_peaks.txt" \
    "${chromosome_sizes}" \
    "${prefix}_peaks.bb" \
    || exit 1

rm "${prefix}_peaks.txt" || exit 1

# NAME_summits.bed

echo "$(date) bedToBigBed: ${prefix}_summits.bed" || exit 1
echo "$(date) bedToBigBed: ${prefix}_summits.bed" 1>&2 || exit 1

grep --extended-regexp --invert-match '^track|^browser' \
    "${prefix}_summits.bed" \
    | perl -e 'while (<>) { chomp; my @fields = split q{ }; $fields[4] = 0; print join(qq{\t}, @fields), qq{\n}; }' \
    > "${prefix}_summits.txt" \
    || exit 1

bedToBigBed \
    "${prefix}_summits.txt" \
    "${chromosome_sizes}" \
    "${prefix}_summits.bb" \
    || exit 1

rm "${prefix}_summits.txt" || exit 1

exit 0
