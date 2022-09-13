#!/usr/bin/env bash
#
# Copyright 2013 - 2022 Michael K. Schuster
#
# BSF Bash script to compress a RNA-seq DESeq2 project directory via Zip.

declare directory_path

directory_path="${1%%/}"

zip \
  --recurse-paths \
  --exclude '..' '*.bash' '*.pkl' '*_completed.txt' '*bsfpython_*_jobs.db' '*/.htaccess' '*/analysis/*' '*/geo/*' '*/scratch/*' \
  --include '*rnaseq_deseq_*' \
  --suffixes '.bam:.bai:.bb:.bw:.rds:.RData:.gz:.zip:.Z' \
  --verbose \
  "${directory_path}_rnaseq_deseq.zip" \
  "${directory_path}" ||
  exit 1

echo 'All done.' || exit 1
