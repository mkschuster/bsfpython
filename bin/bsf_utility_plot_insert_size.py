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
"""The :py:mod:`bin.bsf_utility_plot_insert_size` module is a script to
create box plots of insert sizes of given :emphasis:`BAM` files.
"""

import os
import sys
from argparse import ArgumentParser

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns


def run(
        input_path_list: list[str],
        output_path: str,
        read_count: int) -> int:
    """Run function.

    :param input_path_list: A Python :py:class:`list` object of Python :py:class:`str` (input bam file path) objects.
    :type input_path_list: list[str]
    :param output_path: An output PNG file path.
    :type output_path: str
    :param read_count: A number of reads to be used creating the plots.
    :type read_count: bool | None
    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    matplotlib.use(backend='Agg')

    plot_data_list = []
    sample_name_list = []
    for input_path in input_path_list:
        print(input_path)

        sample_name = os.path.basename(input_path).replace('.bam', '')
        sample_name_list.append(sample_name)

        alignment_file = pysam.AlignmentFile(filename=input_path, mode='rb')

        count = 0
        insert_size_list = []
        for aligned_segment in alignment_file.fetch():
            if count > read_count:
                break

            if aligned_segment.is_paired and aligned_segment.is_read1 and aligned_segment.is_proper_pair and \
                    abs(aligned_segment.template_length) < 1000:
                count += 1
                insert_size_list.append(abs(aligned_segment.template_length))

        plot_data_list.append(insert_size_list)

        alignment_file.close()

    sns.set(style="whitegrid")

    # Set up dataframe
    my_data = {}
    for i in range(len(sample_name_list)):
        my_data[sample_name_list[i]] = plot_data_list[i]

    df = pd.DataFrame(data=my_data)

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=((len(sample_name_list) * 2) / 3, 6))

    # Draw a violin plot with a narrower bandwidth than the default
    sns.violinplot(data=df, palette="Set3", bw=.2, cut=1, linewidth=1)

    # Finalize the figure
    ax.set(ylim=(1, 1000))

    sns.despine(left=True, bottom=True)
    plt.xticks(rotation=90)
    # Save the figure
    fig.savefig(output_path, bbox_inches='tight')

    return 0


def main() -> int:
    """Main function.

    :return: A :py:class:`SystemExit` status value.
    :rtype: int
    """
    argument_parser = ArgumentParser(
        description="Plot insert sizes from bam files")

    argument_parser.add_argument(
        '--output-path',
        help='output PNG file path')

    argument_parser.add_argument(
        '--input-path',
        nargs='+',
        help='input bam file path(s)',
        dest='input_path_list')

    argument_parser.add_argument(
        '--read-count',
        default=1000000,
        type=int,
        help='number of reads to be used creating the plots')

    name_space = argument_parser.parse_args()

    return run(
        input_path_list=name_space.input_path_list,
        output_path=name_space.output_path,
        read_count=name_space.read_count)


if __name__ == '__main__':
    sys.exit(main())
