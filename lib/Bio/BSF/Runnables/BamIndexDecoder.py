"""Bio.BSF.Runnables.IlluminaToBam

A package of classes and methods to run BamIndexDecoder processes.
"""

#
# Copyright 2013 - 2014 Michael K. Schuster
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


import errno
import os
import shutil

from Bio.BSF import Runnable


def run_bam_index_decoder(runnable):
    """Run the 'bam_index_decoder' Executable defined in the BSF Runnable.

    :param runnable: BSF Runnable
    :type runnable: Runnable
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(path=runnable.file_path_dict['metrics']) \
            and os.path.getsize(filename=runnable.file_path_dict['metrics']):
        return

    runnable.run_executable(name='bam_index_decoder')


def run_picard_collect_alignment_summary_metrics(runnable):
    """Run the 'picard_collect_alignment_summary_metrics' Executable defined in the BSF Runnable.

    :param runnable: BSF Runnable
    :type runnable: Runnable
    :return: Nothing
    :rtype: None
    """

    if os.path.exists(path=runnable.file_path_dict['metrics']) \
            and os.path.getsize(filename=runnable.file_path_dict['metrics']):
        return

    runnable.run_executable(name='picard_collect_alignment_summary_metrics')

    # Add a relative symbolic link to the original BSF sequence archive file.

    os.symlink(
        os.path.relpath(runnable.file_path_dict['input'], runnable.file_path_dict['samples_directory']),
        runnable.file_path_dict['link'])


def run(runnable):
    """Run the the BSF Runnable.

    :param runnable: BSF Runnable
    :type runnable: Runnable
    :return: Nothing
    :rtype: None
    """

    path_temporary = runnable.file_path_dict['temporary_directory']

    if not os.path.isdir(path_temporary):
        try:
            os.makedirs(path_temporary)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    path_samples = runnable.file_path_dict['samples_directory']

    if not os.path.isdir(path_samples):
        try:
            os.makedirs(path_samples)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    #  Run all Executable objects of this Runnable.

    if 'bam_index_decoder' in runnable.executable_dict:
        run_bam_index_decoder(runnable=runnable)
    elif 'picard_collect_alignment_summary_metrics' in runnable.executable_dict:
        run_picard_collect_alignment_summary_metrics(runnable=runnable)
    else:
        raise Exception("Runnable.executable_dict with neither 'bam_index_decoder' "
                        "nor 'picard_collect_alignment_summary_metrics' key.")

    # Remove the temporary directory and everything within it.

    shutil.rmtree(path=path_temporary, ignore_errors=False)

    # Job done.
