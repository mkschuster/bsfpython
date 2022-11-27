#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
#
#  Biomedical Sequencing Facility (BSF), part of the genomics core facility
#  of the Research Center for Molecular Medicine (CeMM) of the
#  Austrian Academy of Sciences and the Medical University of Vienna (MUW).
#
#
#  This file is part of BSF Python.
#
#  BSF Python is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BSF Python is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with BSF Python.  If not, see <http://www.gnu.org/licenses/>.
#
#
#  BSF Python script to install BSF Python.
#
import setuptools

with open(file='README.md', mode='r') as text_io:
    long_description = text_io.read()

# https://docs.python.org/3/distutils/setupscript.html
# https://setuptools.pypa.io/en/latest/references/keywords.html
setuptools.setup(
    # https://packaging.python.org/specifications/core-metadata/
    # metadata_version='2.1',
    # Additional Meta-Data
    name='bsfpython-MSchuster',
    version='2022.11.27',
    author='Michael Schuster',
    author_email='mschuster@cemm.oeaw.ac.at',
    maintainer='Michael Schuster',
    maintainer_email='mschuster@cemm.oeaw.ac.at',
    url='https://github.com/mkschuster/bsfpython',
    description='BSF Python Library',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # download_url='',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    # platforms=[],
    # keywords=[],
    # license='',

    packages=setuptools.find_packages(),
    python_requires='>=3.10',
    install_requires=[
        'python-dateutil',
        'azure-storage-blob',
        'biopython',
        'pysam',
        'SQLAlchemy',
        'mysql-connector-python',
        'openpyxl',
    ],
    scripts=[
        'bin/bsf_runner.py',

        # Analysis submission scripts.
        'bin/bsf_submit_bam_index_decoder.py',
        'bin/bsf_submit_bowtie2.py',
        'bin/bsf_submit_bwa_mem.py',
        'bin/bsf_submit_chipseq.py',
        'bin/bsf_submit_ega_cryptor.py',
        'bin/bsf_submit_fastqc.py',
        'bin/bsf_submit_hisat2.py',
        'bin/bsf_submit_illumina_to_bam.py',
        'bin/bsf_submit_irf_archive.py',
        'bin/bsf_submit_irf_processor.py',
        'bin/bsf_submit_irf_restore.py',
        'bin/bsf_submit_kallisto.py',
        'bin/bsf_submit_picard_collect_hiseq_x_pf_fail_metrics.py',
        'bin/bsf_submit_picard_downsample_sam.py',
        'bin/bsf_submit_picard_extract_irf.py',
        'bin/bsf_submit_picard_illumina_demultiplex_sam.py',
        'bin/bsf_submit_picard_illumina_multiplex_sam.py',
        'bin/bsf_submit_picard_sam_to_fastq.py',
        'bin/bsf_submit_rnaseq_deseq.py',
        'bin/bsf_submit_rnaseq_tuxedo.py',
        'bin/bsf_submit_star.py',
        'bin/bsf_submit_tophat2.py',
        'bin/bsf_submit_trimmomatic.py',
        'bin/bsf_submit_variant_calling.py',

        # Utility scripts.
        'bin/bsf_utility_blob_download.py',
        'bin/bsf_utility_blob_download_sequence.py',
        'bin/bsf_utility_blob_upload.py',
        'bin/bsf_utility_check_md5sums.py',
        'bin/bsf_utility_convert_bed_to_interval_list.py',
        'bin/bsf_utility_convert_cadd.py',
        'bin/bsf_utility_convert_sequence_dict_to_chromosome_sizes.py',
        'bin/bsf_utility_convert_vcf_to_23andMe.py',
        'bin/bsf_utility_count_samples_per_lane.py',
        'bin/bsf_utility_meta_data.py',
        'bin/bsf_utility_extract_exonerate_sequences.py',
        'bin/bsf_utility_fix_sample_links.py',
        'bin/bsf_utility_irf_software_versions.py',
        'bin/bsf_utility_irf_summary.py',
        'bin/bsf_utility_library_reverse_complement.py',
        'bin/bsf_utility_ltfs_ltfscp.py',
        'bin/bsf_utility_plot_insert_size.py',
        'bin/bsf_utility_sam_bc_qt_to_fastq.py',
        'bin/bsf_utility_trace_pickler_file.py',
        'bin/bsf_utility_update_ltfs_archive.py',
        'bin/bsf_utility_update_md5sums.py',
        'bin/bsf_utility_uuid_append.py',
        'bin/bsf_utility_validate_library_annotation.py',
        'bin/bsf_utility_zip_bowtie2.bash',
        'bin/bsf_utility_zip_rnaseq_deseq.bash',
        'bin/bsf_utility_zip_star.bash',

        # Legacy scripts.
        'bin/bsf_chipseq_process_macs14.bash',
        'bin/bsf_chipseq_process_macs2.bash',
        'bin/bsf_rnaseq_process_tophat2.bash',
        'bin/bsf_run_bwa.py',
        'bin/bsf_run_rnaseq_tophat.py',
    ]
)
