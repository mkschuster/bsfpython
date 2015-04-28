#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to unpack an Illumina Run Folder
# from a directory of tape archive files. This script depends upon
# GNU tar functionality and will skip and not overwrite older files.
#
#  1. Check if the Illumina Run Folder does already exist.
#  2. Extract the Illumina Run Folder.
#  3. Extract the Intensities directory.
#  4. Extract the Lane directories.
#  5. All done.
#
# Usage: bsf_run_irf_restore.bash archive_directory [output_directory] [force]
#
#   archive_directory: Path to the archive directory containing
#                      tape archive files.
#
#   output_directory: Path to the output directory, in which the
#                     Illumina Run Folder will be recreated. If not specified,
#                     defaults to the current working directory.
#   force: Forces restoring over an already existing Illumina Run Folder.
#          Older files will not be overwritten.
#
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

# Function definitions.

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

if [ "$#" -eq '0' ]; then
    echo "Error: Too few arguments." 1>&2 \
        || bsf_error
    echo "Usage: $(basename ${0}) archive_directory [output_directory] [force]" 1>&2 \
        || bsf_error
    exit 1
fi

declare force='FALSE'
declare archive_directory="${1%/}"
shift

if [ ! -d "${archive_directory}" ]; then
    echo "Error: An archive directory '${archive_directory}' does not exist." 1>&2 \
        || bsf_error
    exit 1
fi

declare archive_name="$(basename "${archive_directory}")"
declare illumina_run_name="${archive_name%_archive}"

# Check if 'force' was set as the second argument.

if [ "${1}" = 'force' ]; then
    force='TRUE'
    shift
fi

# Check if an output directory was set as the second argument.

if [ -n "${1}" ]; then
    declare output_directory="${2}"
    shift
    # The output directory has to exist.
    if [ ! -d "${output_directory}" ]; then
        echo "Error: The output directory '${output_directory}' does not exist." 1>&2 \
            || bsf_error
        exit 1
    fi
else
    declare output_directory="${PWD}"
fi

# Check if force was set as the third argument.

if [ "${1}" = 'force' ]; then
    force='TRUE'
    shift
fi

exec 3>>"${output_directory}/${illumina_run_name}_restore.txt" \
    || bsf_error

# 1. Check if the Illumina Run Folder does already exist.

if [ -d "${output_directory}/${illumina_run_name}" ]; then
    if [ "${force}" = 'TRUE' ]; then
        echo "Info: Forced restoring of an existing Illumina Run Folder ..." 1>&2 \
	        || bsf_error
        echo "INFO: Forced restoring of an existing Illumina Run Folder ..." 1>&3 \
	        || bsf_error
    else
        echo "Error: An Illumina Run Folder already exists at this location." 1>&2 \
            || bsf_error
        echo "Error: An Illumina Run Folder already exists at this location." 1>&3 \
            || bsf_error
        echo "       '${output_directory}/${illumina_run_name}'" 1>&2 \
            || bsf_error
        echo "       '${output_directory}/${illumina_run_name}'" 1>&3 \
            || bsf_error
        echo "       Use the 'force' parameter to start restoring regardless." 1>&2 \
	        || bsf_error
        echo "       Use the 'force' parameter to start restoring regardless." 1>&3 \
	        || bsf_error
        exit 1
    fi
fi

echo "Illumina Run Folder: ${illumina_run_name}" 1>&3 \
    ||bsf_error

cd "${output_directory}" \
    || bsf_error

# 2. Extract the Illumina Run Folder.

echo "$(date): Started extracting the Illumina Run Folder" 1>&3 \
    || bsf_error
tar \
    --extract \
    --file "${archive_directory}/${illumina_run_name}_Folder.tar" \
    || bsf_error
echo "$(date): Finished extracting the Illumina Run Folder" 1>&3 \
    || bsf_error

# 3. Extract the Intensities directory.

echo "$(date): Started extracting the Intensities directory" 1>&3 \
    || bsf_error
tar \
    --extract \
    --file "${archive_directory}/${illumina_run_name}_Intensities.tar" \
    || bsf_error
echo "$(date): Finished extracting the Intensities directory" 1>&3 \
    || bsf_error

# 4. Extract the Lane directories.

for i in 1 2 3 4 5 6 7 8; do

    declare lane_name=$(printf 'L%03u' "${i}")
    declare lane_path="${archive_directory}/${illumina_run_name}_${lane_name}.tar"

    if [ -f "${lane_path}" ]; then
        echo "$(date): Started extracting the Lanes directory ${lane_name}" 1>&3 \
            || bsf_error
        tar \
            --extract \
            --file "${lane_path}" \
            || bsf_error
        echo "$(date): Finished extracting the Lanes directory ${lane_name}" 1>&3 \
            || bsf_error
    else
        echo "$(date): No archive file for Lanes directory ${lane_name}" 1>&3 \
            || bsf_error
    fi

    unset lane_name
    unset lane_path

done

# 5. All done.

echo "$(date): Finished everything ..." 1>&3 \
    || bsf_error

exit 0
