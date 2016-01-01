#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to convert MACS2 bedGraph files into
# BigWig format and BED files into the BigBED format.
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

if [ -z "${LANG}" ]; then
    declare -x LANG='C'
else
    LANG='C'
fi

if test "$#" -lt '2'; then
    echo "Error: bsf_chipseq_process_macs2.sh Too few arguments." 1>&2 \
    || exit 1
    echo "Usage: bsf_chipseq_process_macs2.sh <prefix> <chromosome_sizes>" 1>&2 \
    || exit 1
    exit 1
fi

declare prefix="chipseq_macs2_${1}"
declare chromosome_sizes="${2}"

# Create a temporary directory for sorting as the default one may not be big enough.

declare temporary_directory="${prefix}/temporary"
mkdir -p "${temporary_directory}" || exit 1;
declare -x TMPDIR="${temporary_directory}";

# Convert the following MACS2 bedGraph files into the BigWig format.
# prefix_control_lambda.bdg
# prefix_treat_pileup.bdg
# prefix_treat_pvalue.bdg  Obsolete since MACS 2.0.10
# prefix_treat_qvalue.bdg  Obsolete since MACS 2.0.10

# The UCSC Genome Browser recommends the following sort command.
# sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph

# https://gist.github.com/taoliu/2469050

# TODO: Activate the deletion of the rather big *.bdg files after debugging ...

# NAME_control_lambda.bdg

if [ -f "./${prefix}/${prefix}_control_lambda.bdg" ]; then

    echo "$(date) bedGraphToBigWig: ${prefix}_control_lambda.bdg" || exit 1
    echo "$(date) bedGraphToBigWig: ${prefix}_control_lambda.bdg" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_control_lambda.bdg" \
        | sort -k1,1 -k2,2n \
        | slopBed -i - -g "${chromosome_sizes}" -b 0 \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_clipped.bdg" \
        || exit 1

    bedGraphToBigWig \
        "./${prefix}/${prefix}_clipped.bdg" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_control_lambda.bw" \
        || exit 1

    rm "./${prefix}/${prefix}_clipped.bdg" || exit 1
    mv "${prefix}/${prefix}_control_lambda.bdg" "${prefix}/${prefix}_control_lambda.bdg.bak" || exit 1
#   rm "${prefix}/${prefix}_control_lambda.bdg" || exit 1

fi

if [ -f "./${prefix}/${prefix}_treat_pileup.bdg" ]; then

    echo "$(date) bedGraphToBigWig: ${prefix}_treat_pileup.bdg" || exit 1
    echo "$(date) bedGraphToBigWig: ${prefix}_treat_pileup.bdg" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_treat_pileup.bdg" \
        | sort -k1,1 -k2,2n \
        | slopBed -i - -g "${chromosome_sizes}" -b 0 \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_clipped.bdg" \
        || exit 1

    bedGraphToBigWig \
        "./${prefix}/${prefix}_clipped.bdg" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_treat_pileup.bw" \
        || exit 1

    rm "./${prefix}/${prefix}_clipped.bdg" || exit 1
    mv "./${prefix}/${prefix}_treat_pileup.bdg" "./${prefix}/${prefix}_treat_pileup.bdg.bak" || exit 1
#   rm "./${prefix}/${prefix}_treat_pileup.bdg" || exit 1

fi

# NAME_bdgcmp.bdg

if [ -f "./${prefix}/${prefix}_bdgcmp.bdg" ]; then

    echo "$(date) bedGraphToBigWig: ${prefix}_bdgcmp.bdg" || exit 1
    echo "$(date) bedGraphToBigWig: ${prefix}_bdgcmp.bdg" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_bdgcmp.bdg" \
        | sort -k1,1 -k2,2n \
        | slopBed -i - -g "${chromosome_sizes}" -b 0 \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_clipped.bdg" \
        || exit 1

    bedGraphToBigWig \
        "./${prefix}/${prefix}_clipped.bdg" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_bdgcmp.bw" \
        || exit 1

    rm "./${prefix}/${prefix}_clipped.bdg" || exit 1
    mv "./${prefix}/${prefix}_bdgcmp.bdg" "./${prefix}/${prefix}_bdgcmp.bdg.bak" || exit 1
#   rm "./${prefix}/${prefix}_bdgcmp.bdg" || exit 1

fi

# NAME_treat_pvalue.bdg
# This file is obsolete since version 2.0.10.

if [ -f "./${prefix}/${prefix}_treat_pvalue.bdg" ]; then

    echo "$(date) bedGraphToBigWig: ${prefix}_treat_pvalue.bdg" || exit 1
    echo "$(date) bedGraphToBigWig: ${prefix}_treat_pvalue.bdg" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_treat_pvalue.bdg" \
        | sort -k1,1 -k2,2n \
        | slopBed -i - -g "${chromosome_sizes}" -b 0 \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_clipped.bdg" \
        || exit 1

    bedGraphToBigWig \
        "./${prefix}/${prefix}_clipped.bdg" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_treat_pvalue.bw" \
        || exit 1

    rm "./${prefix}/${prefix}_clipped.bdg" || exit 1
    rm "./${prefix}/${prefix}_treat_pvalue.bdg" || exit 1

fi

# NAME_treat_qvalue.bdg
# This file is obsolete since version 2.0.10.

if [ -f "./${prefix}/${prefix}_treat_qvalue.bdg" ]; then

    echo "$(date) bedGraphToBigWig: ${prefix}_treat_qvalue.bdg" || exit 1
    echo "$(date) bedGraphToBigWig: ${prefix}_treat_qvalue.bdg" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_treat_qvalue.bdg" \
        | sort -k1,1 -k2,2n \
        | slopBed -i - -g "${chromosome_sizes}" -b 0 \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_clipped.bdg" \
        || exit 1

    bedGraphToBigWig \
        "./${prefix}/${prefix}_clipped.bdg" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_treat_qvalue.bw" \
        || exit 1

    rm "./${prefix}/${prefix}_clipped.bdg" || exit 1
    rm "./${prefix}/${prefix}_treat_qvalue.bdg" || exit 1

fi

# Convert the following BED files into the BigBED format.
# prefix_peaks.bed
# prefix_summits.bed

# The UCSC bedToBigBed utility requires that track lines are stripped out. Sigh!
# The UCSC Genome Browser recommends the following sort command.
# sort -k1,1 -k2,2n unsorted.bed > sorted.bed

# TODO: The score would need a complete re-scaling between 0 and 1000. Set to 0 in the meantime.
# The BED format specification requires that the fifth field
# is an integer score between 0 and 1000.
# MACS does not seem to know or care about that.
# Use a Perl one-liner to correct this.
# The bedToBigBed does neither understand a 'stdin' option,
# nor is it capable of reading from standard input.
# Convert via a temporary *.txt file.


# NAME_peaks.bed
# This file seems obsolete in version 2.1.0.

if [ -f "./${prefix}/${prefix}_peaks.bed" ]; then

    echo "$(date) bedToBigBed: ${prefix}_peaks.bed" || exit 1
    echo "$(date) bedToBigBed: ${prefix}_peaks.bed" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_peaks.bed" \
        | sort -k1,1 -k2,2n \
        | perl -e 'while (<>) { chomp; my @fields = split q{ }; $fields[4] = 0; print join(qq{\t}, @fields), qq{\n}; }' \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_peaks.txt" \
        || exit 1

    bedToBigBed \
        "./${prefix}/${prefix}_peaks.txt" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_peaks.bb" \
        || exit 1

    rm "./${prefix}/${prefix}_peaks.txt" || exit 1
#   rm "./${prefix}/${prefix}_peaks.bed" || exit 1

fi


# NAME_peaks.narrowPeak

if [ -f "./${prefix}/${prefix}_peaks.narrowPeak" ]; then

    echo "$(date) bedToBigBed: ${prefix}_peaks.narrowPeak" || exit 1
    echo "$(date) bedToBigBed: ${prefix}_peaks.narrowPeak" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_peaks.narrowPeak" \
        | sort -k1,1 -k2,2n \
        | perl -e 'while (<>) { chomp; my @fields = split q{ }; $fields[4] = 0; print join(qq{\t}, @fields), qq{\n}; }' \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_peaks.txt" \
        || exit 1

    bedToBigBed -type=bed6+4 \
        "./${prefix}/${prefix}_peaks.txt" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_peaks.bb" \
        || exit 1

    rm "./${prefix}/${prefix}_peaks.txt" || exit 1
#   rm "./${prefix}/${prefix}_peaks.bed" || exit 1

fi


# NAME_summits.bed

if [ -f "./${prefix}/${prefix}_summits.bed" ]; then

    echo "$(date) bedToBigBed: ${prefix}_summits.bed" || exit 1
    echo "$(date) bedToBigBed: ${prefix}_summits.bed" 1>&2 || exit 1

    grep --extended-regexp --invert-match '^track|^browser' \
        "./${prefix}/${prefix}_summits.bed" \
        | sort -k1,1 -k2,2n \
        | perl -e 'while (<>) { chomp; my @fields = split q{ }; $fields[4] = 0; print join(qq{\t}, @fields), qq{\n}; }' \
        | bedClip 'stdin' "${chromosome_sizes}" "./${prefix}/${prefix}_summits.txt" \
        || exit 1

    bedToBigBed \
        "./${prefix}/${prefix}_summits.txt" \
        "${chromosome_sizes}" \
        "./${prefix}/${prefix}_summits.bb" \
        || exit 1

    rm "./${prefix}/${prefix}_summits.txt" || exit 1
#   rm "./${prefix}/${prefix}_summits.bed" || exit 1

fi

rm -R "${temporary_directory}" || exit 1;

exit 0
