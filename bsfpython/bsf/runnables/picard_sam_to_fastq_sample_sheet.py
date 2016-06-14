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

from bsf import Runnable
from bsf.argument import OptionLong
from bsf.data import Collection, ProcessedRunFolder, Project, Sample, PairedReads
from bsf.process import RunnableStep


# TODO: The following methods are a copied from the bsf.runnables.generic module.
# _runnable_step_remove_obsolete_file_paths
# _runnable_step_status_file_path
# _runnable_step_status_file_create
# _runnable_step_status_file_remove
# It would be good to define them as methods of the Runnable or RunnableStep method.


def _runnable_step_remove_obsolete_file_paths(runnable_step):
    """Remove the list of file path objects that the RunnableStep declared to be obsolete.

    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @return: Nothing
    @rtype: None
    """
    if runnable_step is None:
        return

    for file_path in runnable_step.obsolete_file_path_list:
        assert isinstance(file_path, (str, unicode))
        if os.path.exists(file_path):
            os.remove(file_path)


def _runnable_step_status_file_path(runnable, runnable_step):
    """Get the status file path for a C{RunnableStep} of a C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @return: Status file path
    @rtype: str
    """

    return '_'.join((runnable.name, runnable_step.name, 'completed.txt'))


def _runnable_step_status_file_create(runnable, runnable_step):
    """Create an empty status file for a C{RunnableStep} of a C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @return: Nothing
    @rtype: None
    """
    if runnable_step is None:
        return

    status_path = _runnable_step_status_file_path(runnable=runnable, runnable_step=runnable_step)
    open(status_path, 'w').close()


def _runnable_step_status_file_remove(runnable, runnable_step):
    """Remove the status file for a C{RunnableStep} of a C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @param runnable_step: C{RunnableStep}
    @type runnable_step: RunnableStep
    @return: Nothing
    @rtype: None
    """

    if runnable_step is None:
        return

    status_path = _runnable_step_status_file_path(runnable=runnable, runnable_step=runnable_step)
    if os.path.exists(status_path):
        os.remove(status_path)


def run(runnable):
    """Run the the C{Runnable}.

    @param runnable: C{Runnable}
    @type runnable: Runnable
    @return: Nothing
    @rtype: None
    """

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

    for prf in collection.processed_run_folders.itervalues():
        assert isinstance(prf, ProcessedRunFolder)
        for project in prf.projects.itervalues():
            assert isinstance(project, Project)
            for sample in project.samples.itervalues():
                assert isinstance(sample, Sample)
                new_paired_reads_list = list()
                for paired_reads in sample.paired_reads_list:
                    assert isinstance(paired_reads, PairedReads)
                    paired_reads_keep = False
                    if paired_reads.reads1 is not None:
                        if os.path.exists(paired_reads.reads1.file_path):
                            if os.path.getsize(paired_reads.reads1.file_path):
                                paired_reads_keep = True
                            else:
                                os.remove(paired_reads.reads1.file_path)
                        else:
                            # The PairedReads object does not have a meaningful Reads object in reads1.
                            paired_reads.reads1 = None
                    if paired_reads.reads2 is not None:
                        if os.path.exists(paired_reads.reads2.file_path):
                            if os.path.getsize(paired_reads.reads2.file_path):
                                paired_reads_keep = True
                            else:
                                os.remove(paired_reads.reads2.file_path)
                        else:
                            # The PairedReads object does not have a meaningful Reads object in reads2.
                            paired_reads.reads2 = None
                    if paired_reads.reads1 is None and paired_reads.reads2 is None:
                        paired_reads_keep = False
                    if paired_reads_keep:
                        new_paired_reads_list.append(paired_reads)
                sample.paired_reads_list = new_paired_reads_list
                # The Sample object could have lost all its PairedReads objects.
                # The PairedReads objects may no longer have the correct weakref to their Sample.

    collection.to_sas_path(name='picard_sam_to_fastq', file_path=new_file_path)

    return
