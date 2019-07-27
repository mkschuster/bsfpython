#!/usr/bin/env bash
# -*- coding: utf-8 -*-
#
# BSF GNU Bourne-Again (Bash) script to convert MACS2 bedGraph files into
# BigWig format and BED files into the BigBED format.
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
    declare -x LANG='C'
else
    LANG='C'
fi

if test "$#" -lt '2';
then
    echo "Error: $(basename "${0}") Too few arguments." 1>&2 \
    || exit 1
    echo "Usage: $(basename "${0}") <prefix> <chromosome_sizes>" 1>&2 \
    || exit 1
    exit 1
fi

declare prefix="${1}"
declare chromosome_sizes="${2}"

# Create a temporary directory for sorting as the default one may not be big enough.

if [[ -d "${prefix}_temporary" ]];
then
    # Rely on the temporary directory created by the bsf.process.Runnable.
    declare -x TMPDIR="${prefix}_temporary";
else
    # Create a new one.
    declare temporary_directory="${prefix}/temporary";
    mkdir -p "${temporary_directory}" || exit 1;
    declare -x TMPDIR="${temporary_directory}";
fi

# Convert the following MACS2 bedGraph files into the BigWig format.
#
# callpeak:
#   prefix_control_lambda.bdg
#   prefix_treat_pileup.bdg
#   prefix_treat_pvalue.bdg  Obsolete since MACS 2.0.10
#   prefix_treat_qvalue.bdg  Obsolete since MACS 2.0.10
#
# bdgcmp:
#   prefix_ppois.bdg
#   prefix_subtract.bdg
#   prefix_logFE.bdg

# The UCSC Genome Browser recommends the following sort command.
# sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph

# https://gist.github.com/taoliu/2469050

declare -a suffixes=('control_lambda' 'treat_pileup' 'ppois' 'subtract' 'logFE');

for suffix in "${suffixes[@]}";
do
    if [[ -f "${prefix}/${prefix}_${suffix}.bdg" && ! -s "${prefix}/${prefix}_${suffix}.bw" ]];
    then
        echo "$(date) bedGraphToBigWig: ${prefix}_${suffix}.bdg" || exit 1;

        grep --extended-regexp --invert-match '^track|^browser' \
            "${prefix}/${prefix}_${suffix}.bdg" \
            | sort -k 1,1 -k 2,2n \
            | slopBed -i - -g "${chromosome_sizes}" -b 0 \
            | bedClip 'stdin' "${chromosome_sizes}" "${prefix}/${prefix}_${suffix}_clipped.bdg" \
            || exit 1;

        if [[ -s "${prefix}/${prefix}_${suffix}_clipped.bdg" ]];
        then
            # Run bedGraphToBigWig only, if the file size is greater than zero.
            bedGraphToBigWig \
                "${prefix}/${prefix}_${suffix}_clipped.bdg" \
                "${chromosome_sizes}" \
                "${prefix}/${prefix}_${suffix}.bw" \
                || exit 1;

            bigWigInfo "${prefix}/${prefix}_${suffix}.bw" > "${prefix}/${prefix}_${suffix}_bwi.txt" \
                || exit 1;
        fi

        rm "${prefix}/${prefix}_${suffix}_clipped.bdg" || exit 1;
        rm "${prefix}/${prefix}_${suffix}.bdg" || exit 1;
    fi
done

# Convert the following BED files into the BigBED format.
#
# prefix_peaks.bed         Obsolete since MACS 2.1.0
# prefix_summits.bed
# prefix_peaks.narrowPeak  Linked to prefix_narrow_peaks.bed symbolically.

if [[ -f "${prefix}/${prefix}_peaks.narrowPeak" && ! -h "${prefix}/${prefix}_narrow_peaks.bed" ]];
then
    cd "${prefix}" || exit 1;
    ln -s \
        "${prefix}_peaks.narrowPeak" \
        "${prefix}_narrow_peaks.bed" \
        || exit 1;
    ln -s \
        "${prefix}_peaks.xls" \
        "${prefix}_peaks.tsv" \
        || exit 1;
    cd '-' || exit 1;
fi

# Convert peak BED files into bigBed format.
# UCSC bigNarrowPeak AutoSql file taken from:
# https://genome-source.gi.ucsc.edu/gitlist/kent.git/blob/master/src/hg/lib/bigNarrowPeak.as

cat > "${prefix}/${prefix}_ucsc_big_narrow_peak.as" << EOF
table bigNarrowPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string  chrom;       "Reference sequence chromosome or scaffold"
    uint    chromStart;  "Start position in chromosome"
    uint    chromEnd;    "End position in chromosome"
    string  name;        "Name given to a region (preferably unique). Use . if no name is assigned"
    uint    score;       "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1] strand;      "+ or - or . for unknown"
    float   signalValue; "Measurement of average enrichment for the region"
    float   pValue;      "Statistical significance of signal value (-log10). Set to -1 if not used."
    float   qValue;      "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
    int     peak;        "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)
EOF

cat > "${prefix}/${prefix}_ucsc_big_summits.as" << EOF
table bigPeakSummit
"BED4+1 Peak summits of signal enrichment based on pooled, normalized (interpreted) data."
(
    string  chrom;       "Reference sequence chromosome or scaffold"
    uint    chromStart;  "Start position in chromosome"
    uint    chromEnd;    "End position in chromosome"
    string  name;        "Name given to a region (preferably unique). Use . if no name is assigned"
    float   pValue;      "Statistical significance of signal value (-log10). Set to -1 if not used."
)
EOF

# The BED format specification requires that the fifth field is an integer score between 0 and 1000,
# which is largely met by MACS2. Use a Perl one-liner to cut scores at 1000.

declare -a suffixes=('peaks' 'summits' 'narrow_peaks');

for suffix in "${suffixes[@]}";
do
    if [[ -f "${prefix}/${prefix}_${suffix}.bed" && ! -s "${prefix}/${prefix}_${suffix}.bb" ]];
    then
        echo "$(date) bedToBigBed: ${prefix}_${suffix}.bed" || exit 1;

        # The UCSC bedToBigBed utility requires that track lines are stripped out. Sigh!
        grep --extended-regexp --invert-match '^track|^browser' \
            "${prefix}/${prefix}_${suffix}.bed" \
            | sort -k 1,1 -k 2,2n \
            | perl -a -n -e 'if ($F[4] > 1000) { $F[4] = 1000 }; print join(qq{\t}, @F), qq{\n};' \
            | bedClip 'stdin' "${chromosome_sizes}" "${prefix}/${prefix}_${suffix}_clipped.bed" \
            || exit 1;

        if [[ -s "${prefix}/${prefix}_${suffix}_clipped.bed" ]];
        then
            # Run bedToBigBed only, if the file size is greater than zero.
            declare cl='';
            cl+='bedToBigBed';
            case ${suffix} in
                ( narrow_peaks )
                    cl+=" -as=${prefix}/${prefix}_ucsc_big_narrow_peak.as";
                    cl+=' -type=bed6+4';
                    ;;
                ( summits )
                    cl+=" -as=${prefix}/${prefix}_ucsc_big_summits.as";
                    cl+=' -type=bed4+1';
                    ;;
            esac
            cl+=" ${prefix}/${prefix}_${suffix}_clipped.bed";
            cl+=" ${chromosome_sizes}";
            cl+=" ${prefix}/${prefix}_${suffix}.bb";
            eval "${cl}" || exit 1;

            bigBedInfo "${prefix}/${prefix}_${suffix}.bb" > "${prefix}/${prefix}_${suffix}_bbi.txt" \
                || exit 1;
        fi

        rm "${prefix}/${prefix}_${suffix}_clipped.bed" || exit 1;
    fi
done

rm "${prefix}/${prefix}_ucsc_big_narrow_peak.as" || exit 1;
rm "${prefix}/${prefix}_ucsc_big_summits.as" || exit 1;

if [[ -f "${prefix}/${prefix}_model.r" ]];
then
    Rscript "${prefix}/${prefix}_model.r";
fi

if [[ -f "${prefix}/${prefix}_model.pdf" ]];
then
    convert "${prefix}/${prefix}_model.pdf" "${prefix}/${prefix}_model.png";
fi

if [[ -d "${temporary_directory}" ]];
then
    rm -R "${temporary_directory}" || exit 1;
fi
