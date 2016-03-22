#! /bin/bash
#
# BSF GNU Bourne-Again (Bash) script to convert an aligned SAM file
# into an aligned sorted and indexed BAM file.
# The samtools application has to be in the PATH.
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

if test "$#" -lt '1'; then
    echo "Error: bsf_sam2bam.sh Too few arguments." 1>&2 \
    || exit 1
    echo "Usage: bsf_sam2bam.sh <prefix>" 1>&2 \
    || exit 1
    exit 1
fi

declare prefix="${1}"

# Convert the SAM (-S) input file into an uncompressed (-u)
# BAM (-b) output stream and sort the BAM output stream.
# Warning: The samtools sort command attaches a *.bam to the prefix.

samtools view -b -S -u "${prefix}.sam" | \
    samtools sort -o "${prefix}.bam" -T "$(dirname ${prefix})" - || exit 1

# Index the aligned and sorted BAM file.

samtools index "${prefix}.bam" "${prefix}.bai" || exit 1

rm "${prefix}.sam" || exit 1

exit 0
