#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to decode an
# Illumina2bam archive BAM file into sample-specific BAM files.
#
# Usage: bsf_bam_index_decoder.sh
#   prefix                  # e.g. BSF_0052_D2C55ACXX_1
#   bam_index_decoder_jar
#
# Copyright 2013 Michael Schuster
# CeMM - Research Center for Molecular Medicine of the
# Austrian Academy of Sciences
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

function bsf_error() {

    cat << EOF_ERROR

 #######
 #        #####   #####    ####   #####
 #        #    #  #    #  #    #  #    #
 #        #    #  #    #  #    #  #    #
 #####    #    #  #    #  #    #  #    #
 #        #####   #####   #    #  #####
 #        #   #   #   #   #    #  #   #
 #######  #    #  #    #   ####   #    #

EOF_ERROR

    exit 1
}

if [ -z "${LANG}" ]; then
    declare -x LANG='C'
else
    LANG='C'
fi

if [ "$#" -ne '2' ]; then
    echo "Error: ${0} Too few arguments." 1>&2  \
 || bsf_error
    echo "Usage: ${0} prefix illumina2bam_jar" 1>&2  \
 || bsf_error
    exit 1
fi

declare prefix="${1}"
declare bam_index_decoder_jar="${2}"

mkdir -p "${prefix}_temporary" || bsf_error
mkdir -p "${prefix}_samples" || bsf_error

# Work with sorted or unsorted archive BAM files, where sorted take precedence.

declare input_file=''
if [ -f "${prefix}_sorted.bam" ]; then
    input_file="${prefix}_sorted.bam"
elif [ -f "${prefix}_unsorted.bam" ]; then
    input_file="${prefix}_unsorted.bam"
fi

java  \
 -Xmx4G  \
 -jar "${bam_index_decoder_jar}"  \
 INPUT="${input_file}"  \
 OUTPUT_DIR="${prefix}_samples"  \
 OUTPUT_PREFIX="${prefix}"  \
 OUTPUT_FORMAT='bam'  \
 BARCODE_FILE="${prefix}_barcode.csv"  \
 METRICS_FILE="${prefix}_metrics.txt"  \
 TMP_DIR="${prefix}_temporary" \
 CREATE_MD5_FILE='true'  \
 VERBOSITY='WARNING'  \
 || bsf_error

rm -R "${prefix}_temporary" || bsf_error
