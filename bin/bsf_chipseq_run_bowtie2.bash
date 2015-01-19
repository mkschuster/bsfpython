#! /bin/bash

# Bowtie2 is a Perl script, but requires bowtie2-align to be in the same directory.
# Since SLURM copies the submitted script into its spool directory, this is not the case.

bowtie2 ${@} || exit 1;
