#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
#
# BSF Python script to create boxplots of insert sizes of given bam files.
#
# Copyright 2017 - 2019 Bekir Erguner
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
#

from __future__ import print_function

import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns

matplotlib.use('Agg')

parser = argparse.ArgumentParser(description="Plot insert sizes from bam files")
parser.add_argument('--output',
                    help='Output .png file',
                    type=str)
parser.add_argument('--input',
                    help='Input bam file(s)',
                    type=str,
                    nargs='+')
parser.add_argument('--count',
                    help='Number of reads to be used creating the plots',
                    type=int,
                    default=1000000)
args = parser.parse_args()

data_to_plot = []
sample_names = []
for i in args.input:
    print(i)
    sample_name = os.path.basename(i).replace('.bam', '')
    sample_names.append(sample_name)
    bam = pysam.AlignmentFile(i, 'rb')
    count = 0
    insert_sizes = []
    for read in bam:
        if count > args.count:
            break
        if read.is_paired and read.is_read1 and read.is_proper_pair and abs(read.template_length) < 1000:
            count += 1
            insert_sizes.append(abs(read.template_length))
    data_to_plot.append(insert_sizes)
    bam.close()

sns.set(style="whitegrid")

# Set up dataframe
my_data = {}
for i in range(len(sample_names)):
    my_data[sample_names[i]] = data_to_plot[i]
df = pd.DataFrame(data=my_data)

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=((len(sample_names) * 2) / 3, 6))

# Draw a violin plot with a narrower bandwidth than the default
sns.violinplot(data=df, palette="Set3", bw=.2, cut=1, linewidth=1)

# Finalize the figure
ax.set(ylim=(1, 1000))

sns.despine(left=True, bottom=True)
plt.xticks(rotation=90)
# Save the figure
f.savefig(args.output, bbox_inches='tight')
