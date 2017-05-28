"""bsf.runnables.picard_sam_to_fastq_sample_sheet

A package of classes and methods to clean a sample sheet for the Picard SamToFastq analysis.
"""

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


import os

from bsf.argument import OptionLong
from bsf.ngs import Collection, ProcessedRunFolder, Project, Sample, PairedReads
from bsf.process import RunnableStep


def run(runnable):
    """Run the the C{bsf.Runnable}.

    @param runnable: C{bsf.Runnable}
    @type runnable: bsf.Runnable
    @return:
    @rtype:
    """

    # If the Runnable status file exists, there is nothing to do and
    # this Runnable should not have been submitted in the first place.

    if os.path.exists(runnable.runnable_status_file_path(success=True)):
        return

    # Do the work.

    runnable_step = runnable.runnable_step_list[0]
    assert isinstance(runnable_step, RunnableStep)
    argument = runnable_step.options['sas_path'][0]
    assert isinstance(argument, OptionLong)
    old_file_path = argument.value
    if old_file_path.endswith('_original.csv'):
        new_file_path = old_file_path[:-13] + '_samples.csv'
    else:
        new_file_path = old_file_path

    collection = Collection.from_sas_path(
        file_path='',
        file_type='',
        name='picard_sam_to_fastq',
        sas_path=argument.value)

    for prf in collection.processed_run_folder_dict.itervalues():
        assert isinstance(prf, ProcessedRunFolder)
        for project in prf.project_dict.itervalues():
            assert isinstance(project, Project)
            for sample in project.sample_dict.itervalues():
                assert isinstance(sample, Sample)
                new_paired_reads_list = list()
                for paired_reads in sample.paired_reads_list:
                    assert isinstance(paired_reads, PairedReads)
                    paired_reads_keep = False
                    if paired_reads.reads_1 is not None:
                        if os.path.exists(paired_reads.reads_1.file_path):
                            if os.path.getsize(paired_reads.reads_1.file_path):
                                paired_reads_keep = True
                            else:
                                os.remove(paired_reads.reads_1.file_path)
                        else:
                            # The PairedReads object does not have a meaningful Reads object in reads1.
                            paired_reads.reads_1 = None
                    if paired_reads.reads_2 is not None:
                        if os.path.exists(paired_reads.reads_2.file_path):
                            if os.path.getsize(paired_reads.reads_2.file_path):
                                paired_reads_keep = True
                            else:
                                os.remove(paired_reads.reads_2.file_path)
                        else:
                            # The PairedReads object does not have a meaningful Reads object in reads2.
                            paired_reads.reads_2 = None
                    if paired_reads.reads_1 is None and paired_reads.reads_2 is None:
                        paired_reads_keep = False
                    if paired_reads_keep:
                        new_paired_reads_list.append(paired_reads)
                sample.paired_reads_list = new_paired_reads_list
                # The Sample object could have lost all its PairedReads objects.
                # The PairedReads objects may no longer have the correct weakref to their Sample.

    collection.to_sas_path(name='picard_sam_to_fastq', file_path=new_file_path)

    runnable_step.remove_obsolete_file_paths()

    # Upon success, create a Runnable-specific status file that indicates completion for the whole Runnable.

    runnable.runnable_status_file_remove()
    runnable.runnable_status_file_create(success=True)

    return
