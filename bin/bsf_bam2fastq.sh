#! /bin/bash
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
#
#
# BSF GNU Bourne-Again (Bash) script to convert an (unmapped) BAM file
# into FASTQ files via Picard SamToFastq and compress them with GNU Gzip.
#
# Usage: bsf_bam2fastq.sh <bam_file> [prefix]
#   bam_file: The path to the BAM file ot be converted
#   prefix: An optional prefix for the FASTQ files, _R1_001.fastq and
#           _R2_001.fastq get appended.
#           If not specified, built from the bam_file parameter by removing
#           .bam and [._]unmapped are removed from the end.
#   output_directory: The directory in which to place the FASTQ files.
#                     Defaults to the directory name of the bam_file parameter.

if [ -z "${LANG}" ]; then
    declare -x LANG='C'
else
    LANG='C'
fi

if test "$#" -lt '1'; then
    echo "Error: bsf_bam2fastq.sh Too few arguments." 1>&2 \
    || exit 1
    echo "Usage: bsf_bam2fastq.sh <bam_file> [output_directory]" 1>&2 \
    || exit 1
    exit 1
fi

declare bam_file="${1}"

if test -n "${2}"; then
    declare directory="${2}"
else
    declare directory="$(dirname "${bam_file}")"
fi

# Get the base name of the bam_file parameter, remove the .bam extension
# and subsequently a potential [._]unmapped part.

declare prefix="$(basename "${bam_file}")"
prefix="${prefix%.bam}"
prefix="${prefix%[._]unmapped}"

declare read1="${directory}/${prefix}_R1_001.fastq"
declare read2="${directory}/${prefix}_R2_001.fastq"

# TODO: Alternatively, pass the Picard directory location into this script.
if test -z "${NGS_PICARD}"; then
    echo "Environment variable NGS_PICARD has to be declared."
    exit 1
fi

if test -f "${read1}.gz" -o -f "${read2}.gz"; then
    echo "BAM file ${bam_file} has already been converted."
else
    java -Xmx2g -jar "${NGS_PICARD}/SamToFastq.jar" \
        INPUT="${bam_file}" \
        FASTQ="${read1}" \
        SECOND_END_FASTQ="${read2}" \
        || exit 1
fi

if test -f "${read1}"; then
    # A regular file exists ...
    if test -s "${read1}"; then
        # ... and has a size grater than zero.
        gzip --best "${read1}" || exit 1
    else
        # ... and is empty.
        rm "${read1}" || exit 1
    fi
fi

if test -f "${read2}"; then
    # A regular file exists ...
    if test -s "${read2}"; then
        # ... and has a size grater than zero.
        gzip --best "${read2}" || exit 1
    else
        # ... and is empty.
        rm "${read2}" || exit 1
    fi
fi

exit 0
