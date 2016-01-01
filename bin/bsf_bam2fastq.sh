#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to convert an (unmapped) BAM or SAM file
# into FASTQ files via Picard SamToFastq before compressing with GNU Gzip.
#
# Usage: bsf_bam2fastq.sh input_file sam_to_fastq_jar [output_directory]
#   input_file:        The path to the BAM or SAM file to be converted
#   sam_to_fastq_jar:  File path to Java Archive (JAR) file (SamToFastq.jar)
#   output_directory:  The directory in which to place the FASTQ files.
#                      Defaults to the directory name of the input_file parameter.
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

function bsf_error () {

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

if test "$#" -lt '2'; then
    echo "Error: Too few arguments." 1>&2 \
    || bsf_error
    echo "Usage: $(basename "${0}") bam_file sam_to_fastq_jar [output_directory]" 1>&2 \
    || bsf_error
    exit 1
fi

declare input_file="${1}"
declare picard_jar="${2}"

if test -n "${3}"; then
    declare directory="${3}"
else
    declare directory="$(dirname "${input_file}")"
fi

# Get the base name of the input_file parameter, remove the .bam or .sam extension
# and subsequently a potential [._]unmapped part.

declare prefix="$(basename "${input_file}")"
prefix="${prefix%.[bs]am}"
prefix="${prefix%[._]unmapped}"

# Construct file names for reads one and two in the style of
# the Illumina CASAVA software package.

declare read1="${directory}/${prefix}_R1_001.fastq"
declare read2="${directory}/${prefix}_R2_001.fastq"

if test -f "${read1}.gz" -o -f "${read2}.gz"; then
    echo "BAM or SAM file ${input_file} has already been converted." \
    || bsf_error
else
    java -Xmx2g -jar "${picard_jar}" \
        'SamToFastq' \
        INPUT="${input_file}" \
        FASTQ="${read1}" \
        SECOND_END_FASTQ="${read2}" \
        || bsf_error
fi

if test -f "${read1}"; then
    # A regular file exists ...
    if test -s "${read1}"; then
        # ... and has a size grater than zero.
        gzip --best "${read1}" \
        || bsf_error
    else
        # ... and is empty.
        rm "${read1}" \
        || bsf_error
    fi
fi

if test -f "${read2}"; then
    # A regular file exists ...
    if test -s "${read2}"; then
        # ... and has a size grater than zero.
        gzip --best "${read2}" \
        || bsf_error
    else
        # ... and is empty.
        rm "${read2}" \
        || bsf_error
    fi
fi

exit 0
