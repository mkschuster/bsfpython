#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to pack an Illumina Run Folder
# for tape archiving.
#
#  1. Check if the Illumina Run has finished by testing for an
#     RTAComplete.txt file.
#  2. Check if an archive process is already running by testing for an
#     archive directory.
#  3. Create an archive directory.
#  4. Size the native Illumina Run Folder via the du utility.
#  5. Reset the file permissions for all directories via the find utility.
#  6. Reset the file permissions for all regular files via the find utility.
#  7. Run the GNU tar utility over each Data/Intensities/L00[1-8] directory,
#     before deleting the directory.
#  8. Run the GNU tar utility over the remaining Data/Intensities directory,
#     before deleting the Intensities directory.
#  9. Run the GNU tar utility over the remaining Illumina Run folder.
# 10. Record the archive file sizes via the ls utility.
#
# Usage: bsf_archive_run_folder.sh illumina_run_folder [output_folder] [force]
#
#   illumina_run_folder: Path to the Illumina Run Folder
#
#   output_folder: Path to the output folder, in which an archive folder
#                  will be created. If not specified, defaults to the current
#                  working directory.
#
#   force: Forces archiving of incomplete run folders
#          (RTAComplete.txt is missing) or a restart of archiving after an
#          archive folder has already been created.
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
    echo "Error: ${0} Too few arguments." 1>&2 \
        || bsf_error
    echo "Usage: ${0} bsf_irf_archive.sh illumina_run_folder [output_folder] [force]" 1>&2 \
        || bsf_error
    exit 1
fi

declare illumina_run_folder="${1%/}"
shift

if [ ! -d "${illumina_run_folder}" ]; then
    echo "Error: ${0} '${illumina_run_folder}' does not exist." 1>&2 \
        || bsf_error
    exit 1
fi

declare illumina_run_name="$(basename "${illumina_run_folder}")"

# Check if 'force' was set as the second argument.

if [ "${1}" = 'force' ]; then
    declare force='TRUE'
    shift
fi

# Check if an archive directory was set as the second argument.

if [ -n "${1}" ]; then
    declare output_folder="${2}"
    shift
    mkdir -p "${output_folder}" \
        || bsf_error
else
    declare output_folder="${PWD}"
fi

# Check if force was set as the third argument.
if [ "${1}" = 'force' ]; then
    declare force='TRUE'
    shift
else
    declare force='FALSE'
fi

exec 3>>"${output_folder}/${illumina_run_name}_archive.txt" \
    || bsf_error

# 0. Check whether Picard ExtractIlluminaBarcodes has written any
# s_<lane>_<tile>_barcode.txt files into the BaseCalls directory.
# Keeping them is rather pointless and they should be removed.
# http://picard.sourceforge.net/command-line-overview.shtml#ExtractIlluminaBarcodes

if [ -n "$(find "${illumina_run_folder}/Data/Intensities/BaseCalls" -name '*_barcode.txt' -print -quit)" ]; then
    echo "Illumina Run Folder contains Picard ExtractIlluminaBarcodes files." 1>&2 \
	    || bsf_error
    echo "Illumina Run Folder contains Picard ExtractIlluminaBarcodes files." 1>&3 \
	    || bsf_error
    echo "Delete those files before restarting the packing process." 1>&2 \
	    || bsf_error
    echo "Delete those files before restarting the packing process." 1>&3 \
	    || bsf_error
    exit 1
fi

# 1. Check whether the RTAComplete.txt file exists in the Illumina Run Folder
# to prevent archiving and deleting an active folder.
# Otherwise require force to start archiving.

if [ ! -f "${illumina_run_folder}/RTAComplete.txt" ]; then
    if [ "${force}" = 'TRUE' ]; then
	    echo "Archive Illumina Run Folder forced incomplete folder ..." \
	        || bsf_error
	    echo "Archive Illumina Run Folder forced incomplete folder ..." 1>&3 \
	        || bsf_error
    else
	    echo "Archive Illumina Run Folder got an incomplete folder." 1>&2 \
	        || bsf_error
	    echo "The RTAComplete.txt file is missing - Exiting." 1>&2 \
	        || bsf_error
	    echo "Use the 'force' parameter to start archiving regardless." 1>&2 \
	        || bsf_error
	    exit 1
    fi
fi

# 2. Use the archive_folder to check, whether an archive process is
# already running. If so, stop here, unless a 'force' parameter has
# been set, in which case resume a archiving process that failed.

declare archive_folder="${output_folder}/${illumina_run_name}_archive"

if [ -d "${archive_folder}" ]; then
    if [ "${force}" = 'TRUE' ]; then
	    echo "Archive Illumina Run Folder forced a re-start ..." \
	        || bsf_error
	    echo "Archive Illumina Run Folder forced a re-start ..." 1>&3 \
	        || bsf_error
    else
	    echo "Archive Illumina Run Folder process already running. Exiting." \
	        || bsf_error
	    echo "Use the 'force' parameter to start archiving regardless." \
	        || bsf_error
	    exit 1
    fi
else
    # 3. Create a archive directory.
    mkdir -p "${archive_folder}" \
        || bsf_error
    # 4. Size the native Illumina Run Folder.
    # The size is only meaningful in the first run.
    echo "Archive Illumina Run Folder: ${illumina_run_name}" 1>&3 \
        || bsf_error
    echo "Size: $(du -k -s "${illumina_run_folder}")" 1>&3 \
        || bsf_error
fi

# 5. Reset all file permissions for directories.

echo "$(date): Started resetting directory access permissions" 1>&3 \
    || bsf_error
find "${illumina_run_folder}" -type d -execdir chmod u=rwx,go=rx {} \+ \
    || bsf_error
echo "$(date): Finished resetting directory permissions" 1>&3 \
    || bsf_error

# 6. Reset all file permissions for regular files.

echo "$(date): Started resetting file access permissions" 1>&3 \
    || bsf_error
find "${illumina_run_folder}" -type f -execdir chmod u=rw,go=r {} \+ \
    || bsf_error
echo "$(date): Finished resetting file permissions" 1>&3 \
    || bsf_error

# Start packing into GNU tape archive files.

declare archive_prefix="${archive_folder}/${illumina_run_name}"

# 7. Run the GNU tar utility over each L00[1-8] directory (8 * 11 % = 88%),
#    before removing it.
#    Each lane directory contains one cluster location file *.clocs per tile and
#    one sub-directory per cycle containing one cluster intensity files *.cif per tile.

declare intensities_folder="${illumina_run_folder}/Data/Intensities"

# PyCharm suggests replacing the expansion with its result.
# for i in {1..8}; do

for i in 1 2 3 4 5 6 7 8; do

    declare lane_name=$(printf 'L%03u' "${i}")
    declare lane_folder="${intensities_folder}/${lane_name}"
    declare state_file="${archive_prefix}_${lane_name}_removing.txt"

    if [ -d "${lane_folder}" ] && [ ! -e "${state_file}" ]; then
	    # The Lane folder is still there and has not already been
	    # deleted partially. The tar process can re-start regardless.
	    echo "$(date): Started archiving ${intensities_folder}/${lane_name}/" 1>&3 \
	        || bsf_error
	    tar -c -f "${archive_prefix}_${lane_name}.tar" \
	        "${intensities_folder}/${lane_name}/" \
	        || bsf_error
	    echo "$(date): Finished archiving ${intensities_folder}/${lane_name}/" 1>&3 \
	        || bsf_error
    fi

    if [ -d "${lane_folder}" ]; then
	    # The Lane folder is still there, so delete it. Touch the
	    # state file before and remove it after successful completion.
    	echo "$(date): Started removing ${intensities_folder}/${lane_name}/" 1>&3 \
	        || bsf_error
	    touch "${state_file}" \
	        || bsf_error
	    rm -f -R "${intensities_folder}/${lane_name}" \
	        || bsf_error
	    echo "$(date): Finished removing ${intensities_folder}/${lane_name}/" 1>&3 \
	        || bsf_error
	    rm "${state_file}" \
	        || bsf_error
    fi

done

unset lane_name
unset lane_folder
unset state_file

# 8. Run the GNU tar utility over the remaining Intensities folder (11 %),
#    before removing it.
#    Aside the lane directories that have been packed and removed before,
#    the Intensities folder contains:
#    * The BaseCalls directory
#    * The Offsets directory
#    * A BaseCalls config.xml file
#    * The RTAConfiguration.xml file

declare state_file="${archive_prefix}_Intensities_removing.txt"

if [ -d "${intensities_folder}" ] && [ ! -e "${state_file}" ]; then
    # The Intensities folder is still there and has not already been
    # deleted partially. The tar process can re-start regardless.
    echo "$(date): Started archiving ${intensities_folder}/" 1>&3 \
	    || bsf_error
    tar -c -f "${archive_prefix}_Intensities.tar" \
	"${intensities_folder}/" \
	    || bsf_error
    echo "$(date): Finished archiving ${intensities_folder}/" 1>&3 \
	    || bsf_error
fi

if [ -d "${intensities_folder}" ]; then
    # The Intensities folder is still there, so delete it. Touch the
    # state file before and remove it after successful completion.
    echo "$(date): Started removing ${intensities_folder}/" 1>&3 \
	    || bsf_error
    touch "${state_file}" \
	    || bsf_error
    rm -f -R "${intensities_folder}" \
	    || bsf_error
    echo "$(date): Finished removing ${intensities_folder}/" 1>&3 \
	    || bsf_error
    rm "${state_file}" \
	    || bsf_error
fi

unset state_file
unset intensities_folder

# 9. Archive, but do not remove the remaining Illumina Run Folder. (<1%)

echo "$(date): Started archiving ${illumina_run_folder}/" 1>&3 \
    || bsf_error
tar -c -f "${archive_prefix}_Folder.tar" \
    "${illumina_run_folder}/" \
    || bsf_error
echo "$(date): Finished archiving ${illumina_run_folder}/" 1>&3 \
    || bsf_error

# 10. Record archive folder file sizes.

echo "$(date): Archive folder file sizes:" 1>&3 \
    || bsf_error

ls -la "${archive_folder}" 1>&3 \
    || bsf_error

# 11. All done ...

echo "$(date): Finished everything ..." 1>&3 \
    || bsf_error

exit 0
