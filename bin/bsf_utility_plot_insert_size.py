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
#  BSF Python script to create box plots of insert sizes of given bam files.
#
import os
from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns

matplotlib.use('Agg')

argument_parser = ArgumentParser(
    description="Plot insert sizes from bam files")

argument_parser.add_argument(
    '--output',
    help='Output .png file',
    type=str)

argument_parser.add_argument(
    '--input',
    help='Input bam file(s)',
    type=str,
    nargs='+')

argument_parser.add_argument(
    '--count',
    help='Number of reads to be used creating the plots',
    type=int,
    default=1000000)

name_space = argument_parser.parse_args()

data_to_plot = []
sample_names = []
for input_path in name_space.input:
    print(input_path)

    sample_name = os.path.basename(input_path).replace('.bam', '')
    sample_names.append(sample_name)

    alignment_file = pysam.AlignmentFile(filename=input_path, mode='rb')

    count = 0
    insert_sizes = []
    for aligned_segment in alignment_file.fetch():
        if count > name_space.count:
            break

        if aligned_segment.is_paired and aligned_segment.is_read1 and aligned_segment.is_proper_pair and \
                abs(aligned_segment.template_length) < 1000:
            count += 1
            insert_sizes.append(abs(aligned_segment.template_length))

    data_to_plot.append(insert_sizes)

    alignment_file.close()

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
f.savefig(name_space.output, bbox_inches='tight')
