#! /usr/bin/env python
#
# BSF Python script to list multiplexing indices (barcodes).
#
#
# Copyright 2012-2013 Andreas Schoenegger
# Copyright 2013 Michael K. Schuster
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

from collections import OrderedDict
import platform
import os
import sys


# TODO: dual indices not supported
# TODO: GNU Zip FASTQ files could be read directly: http://docs.python.org/2/library/gzip.html
# TODO: Use the argparse module for command-line options.

if len(sys.argv) < 2:
    print "input filename is missing!"
    exit(1)

input_filename = sys.argv[1]
file_type = os.path.splitext(input_filename)[-1]     # supported formats: .fastq/.sam/.bam

log_after_x_processed_reads = 1000000


def parse_sam_file(input_filename):

    """
    Parse a SAM file, processing barcodes in the 12th column.
    :param input_filename: File path
    :type input_filename: str, unicode
    :return: Python dict of barcode key and count value data
    :rtype: dict
    """

    # TODO: breaks if barcode not in 12th column or tag has different name (expecting BC:Z:....)

    print "Processing SAM file (expecting barcode in the 12th column of the SAM file as a tag, " \
          "like BC:Z:GGAACC (Picard - IlluminaToSam creates these unmapped BAM files)"

    input_file = open(input_filename)
    barcodes = dict()

    i = 0
    for line in input_file:
        # ignore SAM header (starts with @)
        if not line.startswith("@"):
            # split each line at tabs, extract 12 column, split barcode tag
            bc = line.split("\t")[11].split(":")[2]

            # store barcode in a dict
            if not barcodes.has_key(bc):
                barcodes[bc] = 0

            barcodes[bc] += 1

        # log message
        i += 1
        if i % log_after_x_processed_reads == 0:
            print "Read " + str(i) + " lines."

    input_file.close()

    return barcodes


# converts a BAM file to SAM, runs parse_sam_file on the temporary SAM file
def parse_bam_file(input_filename):

    """
    Parse a BAM file.
    Internally converts a BA;M file to a SAM file and subsequently runs parse_sam_file over it.
    :param input_filename: File path
    :type input_filename: str, unicode
    :return: Python dict of barcode key and count value data
    :rtype: dict
    """

    if platform.system() == "Windows":
        raise Exception("BAM files not supported on Windows")

    print "Processing BAM file, running samtools view <bamfile> > <tmp_samfile> (samtools required in PATH!)"

    tmp_samfile = os.path.splitext(input_filename)[0]+ ".tmp.sam"

    exit_code = os.system("samtools view -h " + input_filename + " > " + tmp_samfile)

    if not exit_code == 0:
        raise Exception("Error executing samtools. Check if samtools is in your PATH (e.g. which samtools)")

    barcodes = parse_sam_file(tmp_samfile)

    os.remove(tmp_samfile)

    return barcodes


# parses a FASTQ file, processes barcode in header line of each FASTQ entry
# TODO: breaks if barcode has /1 or /2 at the end (I think FASTQ files created by Picard from PE BAM files have these)
def parse_fastq_file(input_filename):

    print "Processing FASTQ file (expecting barcode in the header)"

    input_file = open(input_filename)
    barcodes = dict()

    i = 0
    for line in input_file:
        if i % 4 == 0:
            bc = line.rstrip().split(":").pop()
            if not barcodes.has_key(bc):
                barcodes[bc] = 0

            barcodes[bc] += 1

        i += 1
        if i % log_after_x_processed_reads == 0:
            print "Read " + str(i) + " lines."

    input_file.close()
    return barcodes


if file_type == ".fastq":
    barcodes = parse_fastq_file(input_filename)

elif file_type == ".sam":
    barcodes = parse_sam_file(input_filename)

elif file_type == ".bam":
    barcodes = parse_bam_file(input_filename)

else:
    raise Exception()


# write report
outfile = open("unmatched_barcode_report.csv", "w")
# print dict sorted by highest number of occurence
for bc in OrderedDict(sorted(barcodes.items(), reverse=True, key=lambda t: t[1])):
    message = bc + ";" + str(barcodes[bc])

    outfile.write(message + "\n")
    print message

# Check again for Illumina Barcodes:
outfile.close()

print "--------"

illumina_adapters = ["ATCACG","CGATGT","TTAGGC","TGACCA","ACAGTG","GCCAAT",
                     "CAGATC","ACTTGA","GATCAG","TAGCTT","GGCTAC","CTTGTA",
                     "AGTCAA","AGTTCC","ATGTCA","CCGTCC","GTCCGC","GTGAAA",
                     "GTGGCC","GTTTCG","CGTACG","GAGTGG","ACTGAT","ATTCCT"]

outfile = open("unmatched_barcode_report.illumina.csv", "w")

for adapter in illumina_adapters:
    count=0
    if barcodes.has_key(adapter):
        count = barcodes[adapter]

    message = adapter + ";" + str(count)

    outfile.write(message + "\n")
    print message

outfile.close()
