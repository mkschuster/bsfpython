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
"""The :py:mod:`setup` module is a script to install BSF Python.
"""

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
    version='2023.6.22',
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
        # Analysis submission scripts.
        'bin/bsf_submit_irf_processor.py',

        # Utility scripts.
        'bin/bsf_utility_convert_bed_to_interval_list.py',
        'bin/bsf_utility_convert_cadd.py',
        'bin/bsf_utility_convert_sequence_dict_to_chromosome_sizes.py',
        'bin/bsf_utility_convert_vcf_to_23andMe.py',
        'bin/bsf_utility_count_samples_per_lane.py',
        'bin/bsf_utility_extract_exonerate_sequences.py',
        'bin/bsf_utility_fix_sample_links.py',
        'bin/bsf_utility_plot_insert_size.py',
        'bin/bsf_utility_sam_bc_qt_to_fastq.py',
        'bin/bsf_utility_uuid_append.py',
        'bin/bsf_utility_zip_bowtie2.bash',
        'bin/bsf_utility_zip_rnaseq_deseq.bash',
        'bin/bsf_utility_zip_star.bash',

        # Legacy scripts.
        'bin/bsf_chipseq_process_macs14.bash',
        'bin/bsf_chipseq_process_macs2.bash',
        'bin/bsf_rnaseq_process_tophat2.bash',
        'bin/bsf_run_bwa.py',
        'bin/bsf_run_rnaseq_tophat.py',
    ],
    entry_points={
        'console_scripts': [
            'bsf_run_runnable = bsf.procedure:Runnable.entry_point_run',
            'bsf_submit_bam_index_decoder = bsf.analyses.illumina_to_bam_tools:BamIndexDecoder.entry_point_submit',
            'bsf_submit_illumina_to_bam = bsf.analyses.illumina_to_bam_tools:IlluminaToBam.entry_point_submit',
            'bsf_submit_bowtie1 = bsf.analyses.bowtie:Bowtie1.entry_point_submit',
            'bsf_submit_bowtie2 = bsf.analyses.bowtie:Bowtie2.entry_point_submit',
            'bsf_submit_bwa_mem = bsf.analyses.bwa:MaximalExactMatches.entry_point_submit',
            'bsf_submit_chipseq = bsf.analyses.chipseq:ChIPSeq.entry_point_submit',
            'bsf_submit_ega_cryptor = bsf.analyses.ega:EGACryptor.entry_point_submit',
            'bsf_submit_fastqc = bsf.analyses.fastqc:FastQC.entry_point_submit',
            'bsf_submit_hisat2 = bsf.analyses.hisat:Hisat2.entry_point_submit',
            'bsf_submit_irf_archive = bsf.analyses.illumina_run_folder:IlluminaRunFolderArchive.entry_point_submit',
            'bsf_submit_irf_restore = bsf.analyses.illumina_run_folder:IlluminaRunFolderRestore.entry_point_submit',
            'bsf_submit_kallisto = bsf.analyses.kallisto:Kallisto.entry_point_submit',
            'bsf_submit_picard_collect_hiseq_x_pf_fail_metrics = '
            'bsf.analyses.picard:CollectHiSeqXPfFailMetrics.entry_point_submit',
            'bsf_submit_picard_downsample_sam = bsf.analyses.picard:DownsampleSam.entry_point_submit',
            'bsf_submit_picard_extract_irf = bsf.analyses.picard:ExtractIlluminaRunFolder.entry_point_submit',
            'bsf_submit_picard_illumina_demultiplex_sam = '
            'bsf.analyses.picard:IlluminaDemultiplexSam.entry_point_submit',
            'bsf_submit_picard_illumina_multiplex_sam = '
            'bsf.analyses.picard:IlluminaMultiplexSam.entry_point_submit',
            'bsf_submit_variant_calling = bsf.analyses.variant_calling:VariantCallingGATK.entry_point_submit',
            'bsf_submit_picard_sam_to_fastq = bsf.analyses.picard:SamToFastq.entry_point_submit',
            'bsf_submit_rnaseq_deseq = bsf.analyses.rnaseq:DESeq.entry_point_submit',
            'bsf_submit_rnaseq_tuxedo = bsf.analyses.rnaseq:Tuxedo.entry_point_submit',
            'bsf_submit_star = bsf.analyses.star:Star.entry_point_submit',
            'bsf_submit_tophat2 = bsf.analyses.tophat:Tophat2.entry_point_submit',
            'bsf_submit_trimmomatic = bsf.analyses.trimmomatic:Trimmomatic.entry_point_submit',
            'bsf_utility_blob_download = bsf.cloud:entry_point_blob_download',
            'bsf_utility_blob_download_sequence = bsf.cloud:entry_point_blob_download_sequence',
            'bsf_utility_blob_upload = bsf.cloud:entry_point_blob_upload',
            'bsf_utility_blob_list = bsf.cloud:entry_point_blob_list',
            'bsf_utility_irf_software_versions = bsf.illumina:RunFolder.entry_point_software_versions',
            'bsf_utility_irf_summary = bsf.illumina:RunFolder.entrypoint_summary',
            'bsf_utility_ltfs_archive_update = bsf.ltfs:entry_point_archive_update',
            'bsf_utility_ltfs_drive_archive = bsf.ltfs:entry_point_drive_archive',
            'bsf_utility_ltfs_drive_mount = bsf.ltfs:entry_point_drive_mount',
            'bsf_utility_ltfs_drive_restore = bsf.ltfs:entry_point_drive_restore',
            'bsf_utility_ltfs_drive_size = bsf.ltfs:entry_point_drive_size',
            'bsf_utility_ltfs_library_archive = bsf.ltfs:entry_point_library_archive',
            'bsf_utility_ltfs_library_mount = bsf.ltfs:entry_point_library_mount',
            'bsf_utility_ltfs_library_multi_restore = bsf.ltfs:entrypoint_library_multi_restore',
            'bsf_utility_ltfs_library_restore = bsf.ltfs:entry_point_library_restore',
            'bsf_utility_ltfs_library_size = bsf.ltfs:entry_point_library_size',
            'bsf_utility_ltfs_ltfscp = bsf.ltfs:entry_point_ltfscp',
            'bsf_utility_md5sum_check = bsf.md5sum:MD5SumArchive.entry_point_check',
            'bsf_utility_md5sum_update = bsf.md5sum:MD5SumArchive.entry_point_update',
            'bsf_utility_reset_demultiplexing = bsf.analyses.picard:IlluminaDemultiplexSam.entry_point_reset',
            'bsf_utility_library_annotation_validate = '
            'bsf.analyses.illumina_to_bam_tools:LibraryAnnotationSheet.entry_point_validate',
            'bsf_utility_library_annotation_reverse = '
            'bsf.analyses.illumina_to_bam_tools:LibraryAnnotationSheet.entry_point_reverse',
            'bsf_utility_runnable_trace = bsf.procedure:Runnable.entry_point_trace',
            'bsf_utility_variant_calling_find_missing = '
            'bsf.analyses.variant_calling:VariantCallingGATK.entry_point_find_missing',
        ],
    })
