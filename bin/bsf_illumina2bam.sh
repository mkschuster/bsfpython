#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to convert a lane of an
# Illumina flow-cell into an archive BAM file sorted by query name.
#
# Usage: bsf_illumina2bam.sh
#   run_folder       # e.g. 130926_SN181_0391_AD28APACXX
#   lane             # e.g. 1
#   center           # e.g. BSF
#   experiment       # e.g. BSF_0000
#   illumina2bam_jar # file path to Java Archive (JAR) file
#   sort_sam_jar     # file path to Java Archive (JAR) file
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

if [ "$#" -ne '6' ]; then
    echo "Error: ${0} Too few arguments." 1>&2 \
    || bsf_error
    echo "Usage: ${0} run_folder lane sequencing_center experiment illumina2bam_jar sort_sam_jar" 1>&2 \
    || bsf_error
    exit 1
fi

declare run_folder="${1%/}"
declare lane="${2}"
declare center="${3}"
declare experiment="${4}"
declare illumina2bam_jar="${5}"
declare sort_sam_jar="${6}"
declare barcode="${run_folder##*_[A-B]}"
declare prefix="${experiment}_${barcode}_${lane}"

if [ ! -d "${run_folder}" ]; then
    echo "Error: ${0} Illumina Run Folder ${run_folder} does not exist." 1>&2 \
    || bsf_error
    exit 1
fi

if [ ! -f "${run_folder}/RTAComplete.txt" ]; then
    echo "Error: ${0} Illumina Run Folder not complete." 1>&2 \
    || bsf_error
    exit 1
fi

if [ ! -f "${illumina2bam_jar}" ]; then
    echo "Error: ${0} Illumina2bam Java Archive (JAR) file does not exist." 1>&2 \
    || bsf_error
    exit 1
fi

if [ ! -f "${sort_sam_jar}" ]; then
    echo "Error: ${0} Picard SortSam Java Archive (JAR) file does not exist." 1>&2 \
    || bsf_error
    exit 1
fi

mkdir -p "${prefix}_temporary" || bsf_error

# If the unsorted.bam file and the corresponding md5 file are present,
# the Illumina2bam stage does not need re-running.

if [ -f "${prefix}_unsorted.bam" ] && [ -f "${prefix}_unsorted.bam.md5" ]; then
    echo "Skipping Illumina2bam step that has already run successfully." \
    || bsf_error
    echo "Files '${prefix}_unsorted.bam' and '${prefix}_unsorted.bam.md5' are present." \
    ||bsf_error
else
    java  \
    -Xmx4G  \
    -jar "${illumina2bam_jar}"  \
    INTENSITY_DIR="${run_folder}/Data/Intensities"  \
    LANE="${lane}"  \
    OUTPUT="${prefix}_unsorted.bam"  \
    GENERATE_SECONDARY_BASE_CALLS='false'  \
    PF_FILTER='false'  \
    READ_GROUP_ID="${barcode}_${lane}"  \
    LIBRARY_NAME="${barcode}_${lane}"  \
    SEQUENCING_CENTER="${center}"  \
    TMP_DIR="${prefix}_temporary"  \
    VERBOSITY='WARNING'  \
    MAX_RECORDS_IN_RAM='4000000'  \
    CREATE_INDEX='false'  \
    CREATE_MD5_FILE='true'  \
    || bsf_error
fi

if [ -f "${prefix}_sorted.bam" ] && [ -f "${prefix}_sorted.bam.md5" ]; then
    echo "Skipping SortSam step that has already run successfully." \
    || bsf_error
    echo "Files '${prefix}_sorted.bam' and '${prefix}_sorted.bam.md5' are present." \
    || bsf_error
else
    java  \
    -Xmx4G  \
    -jar "${sort_sam_jar}"  \
    INPUT="${prefix}_unsorted.bam"  \
    OUTPUT="${prefix}_sorted.bam"  \
    SORT_ORDER='queryname'  \
    TMP_DIR="${prefix}_temporary"  \
    VERBOSITY='WARNING'  \
    MAX_RECORDS_IN_RAM='4000000'  \
    CREATE_INDEX='false'  \
    CREATE_MD5_FILE='true'  \
    || bsf_error
fi

rm "${prefix}_unsorted.bam" || bsf_error
rm "${prefix}_unsorted.bam.md5" || bsf_error
rm -R "${prefix}_temporary" || bsf_error
