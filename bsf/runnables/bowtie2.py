# -*- coding: utf-8 -*-
"""Bowtie2 Runnable module

A package of classes and methods to run the Bowtie2 short read aligner.
"""
#  Copyright 2013 - 2019 Michael K. Schuster
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
from __future__ import print_function

import errno
import os
import re
import shutil
import sys

import bsf
import bsf.analyses.aligner
import bsf.argument
import bsf.process
import bsf.standards


class ReadGroupContainer(object):
    """The C{bsf.runnables.bowtie2.ReadGroupContainer} models a SAM or BAM read group (@RG) line.

    Attributes:
    @ivar rg_string: SAM Read Group (@RG)
    @type rg_string: str
    @ivar fastq_1_path: FASTQ R1 file path
    @type fastq_1_path: str | unicode
    @ivar fastq_2_path: FASTQ R2 file path
    @type fastq_2_path: str | unicode
    """

    def __init__(
            self,
            rg_string=None,
            fastq_1_path=None,
            fastq_2_path=None):
        """Initialise a C{bsf.runnables.bowtie2.ReadGroupContainer} object.

        @param rg_string: SAM Read Group (@RG)
        @type rg_string: str
        @param fastq_1_path: FASTQ R1 file path
        @type fastq_1_path: str | unicode
        @param fastq_2_path: FASTQ R2 file path
        @type fastq_2_path: str | unicode
        @return:
        @rtype:
        """

        super(ReadGroupContainer, self).__init__()

        if rg_string:
            self.rg_string = rg_string
        else:
            self.rg_string = str()

        if fastq_1_path:
            self.fastq_1_path = fastq_1_path
        else:
            self.fastq_1_path = str()

        if fastq_2_path:
            self.fastq_2_path = fastq_2_path
        else:
            self.fastq_2_path = str()

    def get_rg_element(self, tag):

        if len(tag) != 2:
            raise Exception('SAM element tag ' + repr(tag) + ' is not exactly two characters long.')

        tag_string = tag + ':'
        for element in self.rg_string.split():
            if element[:3] == tag_string:
                return element[4:]

    def get_rg_id(self):

        return self.get_rg_element(tag='ID')

    def get_rg_pu(self):

        return self.get_rg_element(tag='PU')

    def get_rg_elements(self):
        element_list = list()
        """ @type element_list: list[str] """
        for element in self.rg_string.split():
            if element[:3] != 'ID:':
                element_list.append(element)

        return element_list


def run_picard_sam_to_fastq(runnable, bam_file_path):
    """Run Picard SamToFastq on a BAM file and convert into a FASTQ file pair per read group (@RG).

    Expand a BAM file into a pair of FASTQ files per SAM read group.
    @param runnable: C{bsf.procedure.ConsecutiveRunnable}
    @type runnable: bsf.procedure.ConsecutiveRunnable
    @param bam_file_path: BAM file path
    @type bam_file_path: str | unicode
    @return:
    @rtype:
    """
    # Propagate SAM header lines @PG and @RG around FASTQ files into the final BAM file. Sigh!

    bam_file_name = os.path.basename(bam_file_path)
    sam_file_path = os.path.join(runnable.temporary_directory_path(absolute=False), bam_file_name[:-4] + '.sam')

    samtools = bsf.process.Executable(
        name='samtools_view',
        program='samtools',
        sub_command=bsf.process.Command(program='view'),
        stdout_path=sam_file_path)

    samtools_view = samtools.sub_command
    samtools_view.add_switch_short(key='H')
    samtools_view.arguments.append(bam_file_path)

    child_return_code = samtools.run()
    if child_return_code:
        raise Exception(
            'Could not complete the ' + repr(samtools.name) + ' step on the BAM file for the replicate.')

    sam_pg_list = list()
    """ @type sam_pg_list: list[str] """
    sam_rg_list = list()
    """ @type sam_rg_list: list[str] """

    with open(sam_file_path, 'rt') as output_file:
        for line_str in output_file:
            if line_str.startswith('@PG'):
                sam_pg_list.append(line_str.rstrip())
            if line_str.startswith('@RG'):
                sam_rg_list.append(line_str.rstrip())

    os.remove(sam_file_path)

    # At this stage, the SAM @PG and @RG lines are stored internally.
    # Now run Picard SamToFastq to convert.

    java_process = bsf.process.Executable(name='sam_to_fastq', program='java', sub_command=bsf.process.Command())
    java_process.add_switch_short(key='d64')
    java_process.add_switch_short(key='server')
    java_process.add_switch_short(key='Xmx2G')
    java_process.add_option_pair(key='-Djava.io.tmpdir', value=runnable.temporary_directory_path(absolute=False))

    picard_process = java_process.sub_command
    picard_process.add_option_short(key='jar', value=os.path.join(
        bsf.standards.JavaClassPath.get_picard(), 'picard.jar'))
    picard_process.sub_command = bsf.process.Command(program='SamToFastq')

    sam_to_fastq = picard_process.sub_command
    sam_to_fastq.add_option_pair(key='INPUT', value=bam_file_path)
    sam_to_fastq.add_option_pair(key='OUTPUT_PER_RG', value='true')
    sam_to_fastq.add_option_pair(key='OUTPUT_DIR', value=runnable.temporary_directory_path(absolute=False))
    sam_to_fastq.add_option_pair(key='INCLUDE_NON_PF_READS', value='false')  # TODO: Make this configurable.
    sam_to_fastq.add_option_pair(key='TMP_DIR', value=runnable.temporary_directory_path(absolute=False))
    sam_to_fastq.add_option_pair(key='VERBOSITY', value='WARNING')
    sam_to_fastq.add_option_pair(key='QUIET', value='false')
    sam_to_fastq.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')

    child_return_code = java_process.run()
    if child_return_code:
        raise Exception('Could not complete the ' + repr(java_process.name) + ' step.')

    rgc_list = list()
    """ @type rgc_list: list[ReadGroupContainer] """
    regular_expression = re.compile(pattern='\\W')
    # Expect a FASTQ file pair for each read group.
    for sam_rg in sam_rg_list:
        rgc = ReadGroupContainer(rg_string=sam_rg)
        rgc_list.append(rgc)
        for element in sam_rg.split():
            if element[:3] == 'PU:':
                # Picard builds file names from the PU string, but replaces all non-alphanumeric characters.
                sam_pu = re.sub(pattern=regular_expression, repl='_', string=element[3:])
                sam_pu_r1_path = os.path.join(runnable.temporary_directory_path(absolute=False), sam_pu + '_1.fastq')
                sam_pu_r2_path = os.path.join(runnable.temporary_directory_path(absolute=False), sam_pu + '_2.fastq')

                if os.path.exists(sam_pu_r1_path) and os.path.getsize(sam_pu_r1_path):
                    rgc.fastq_1_path = sam_pu_r1_path
                if os.path.exists(sam_pu_r2_path) and os.path.getsize(sam_pu_r2_path):
                    rgc.fastq_2_path = sam_pu_r2_path

                if not rgc.fastq_1_path:
                    raise Exception(
                        'SAM Read Group ' + repr(sam_rg) + ' without R1 FASTQ file ' + repr(sam_pu_r1_path) + '.')

    return rgc_list


def run_bowtie2(runnable):
    """Run Bowtie2.

    @param runnable: C{bsf.procedure.ConsecutiveRunnable}
    @type runnable: bsf.procedure.ConsecutiveRunnable
    @return:
    @rtype:
    """

    # Put all sample-specific information into a sub-directory.

    file_path_read_group = runnable.file_path_object
    """ @type file_path_read_group: bsf.analyses.aligner.FilePathAlign """

    if not os.path.isdir(file_path_read_group.output_directory):
        try:
            os.makedirs(file_path_read_group.output_directory)
        except OSError as exception:  # Python > 2.5
            if exception.errno != errno.EEXIST:
                raise

    if os.path.exists(path=file_path_read_group.aligned_sam) \
            and os.path.getsize(filename=file_path_read_group.aligned_sam):
        return

    # FIXME: The ConsecutiveRunnable object now only has a list of RunnableStep objects.
    # FIXME: The bsf.executables.Bowtie2 object is also deprecated.
    # bowtie2 = runnable.executable_dict['bowtie2']
    # The following bowtie2 definition is only a placeholder.
    # bowtie2 = bsf.executables.Bowtie2(name='bowtie2', analysis=Analysis())

    # FIXME: Run bowtie via a RunnableStep to avoid SLURM copying it into the spool directory.
    # Bowtie2 is a Perl script, but requires bowtie2-align to be in the same directory.
    # Since SLURM copies the submitted script into its spool directory, this is not the case.
    bowtie2 = bsf.process.Executable(name='bowtie2', program='bowtie2')

    # aligned_sam = bowtie2.stdout_path

    # TODO: For the moment, convert only files set in the bowtie2 -U option.
    # Pop the original list of -U options off the Bowtie2 Executable object.
    if 'U' in bowtie2.options:
        option_list_u = list(bowtie2.options['U'])
    else:
        option_list_u = list()

    if '1' in bowtie2.options:
        option_list_1 = list(bowtie2.options['1'])
    else:
        option_list_1 = list()

    if '2' in bowtie2.options:
        option_list_2 = list(bowtie2.options['2'])
    else:
        option_list_2 = list()

    rgc_list = list()
    """ @type rgc_list: list[ReadGroupContainer] """
    fastq_list = list()
    """ @type fastq_list: list[str | unicode] """

    # Expand all BAM files into ReadGroupContainers representing FASTQ files and

    for option in option_list_u:
        """ @type option: bsf.argument.OptionShort """
        assert isinstance(option, bsf.argument.OptionShort)
        for file_path in option.value.split(','):
            if file_path.endswith('.bam'):
                rgc_list.extend(run_picard_sam_to_fastq(runnable=runnable, bam_file_path=file_path))
            else:
                fastq_list.append(file_path)

    # First, process all ReadGroupContainer objects.

    sam_file_path_list = list()
    """ @type sam_file_path_list: list[str | unicode] """

    for rgc in rgc_list:
        # Clear Bowtie2 options.
        bowtie2.options.pop('U', None)
        bowtie2.options.pop('1', None)
        bowtie2.options.pop('2', None)

        if rgc.fastq_1_path and rgc.fastq_2_path:
            bowtie2.add_option_short(key='1', value=rgc.fastq_1_path)
            bowtie2.add_option_short(key='2', value=rgc.fastq_2_path)
        elif rgc.fastq_1_path:
            bowtie2.add_option_short(key='U', value=rgc.fastq_1_path)
        elif rgc.fastq_2_path:
            raise Exception('Received a FASTQ R2 file ' + repr(rgc.fastq_2_path) + ', but no FASTQ R1 file.')
        else:
            raise Exception('Received neither a FASTQ R1 nor a FASTQ R2 file.')

        # Set the STDOUT path by removing '_1.fastq' from the temporary R1 FASTQ file.
        sam_file_path = rgc.fastq_1_path[:-8] + '.sam'
        bowtie2.stdout_path = sam_file_path
        sam_file_path_list.append(sam_file_path)

        sys.stdout.writelines(bowtie2.trace(level=1))

        child_return_code = bowtie2.run()
        if child_return_code:
            raise Exception('Could not complete the ' + repr(bowtie2.name) + ' step.')

        # Remove the temporary FASTQ files.

        if os.path.exists(rgc.fastq_1_path):
            os.remove(rgc.fastq_1_path)
        if os.path.exists(rgc.fastq_2_path):
            os.remove(rgc.fastq_2_path)

    # Second, process the remaining FASTQ files.

    if len(fastq_list):

        bowtie2.options.pop('U', None)
        bowtie2.options.pop('1', None)
        bowtie2.options.pop('2', None)

        bowtie2.add_option_short(key='U', value=','.join(fastq_list))

        sam_file_path = os.path.join(
            runnable.temporary_directory_path(absolute=False),
            'unpaired_fastq_files.sam')
        bowtie2.stdout_path = sam_file_path
        sam_file_path_list.append(sam_file_path)

        sys.stdout.writelines(bowtie2.trace(level=1))

        child_return_code = bowtie2.run()
        if child_return_code:
            raise Exception('Could not complete the ' + repr(bowtie2.name) + ' step.')

            # Do not delete the non-temporary FASTQ files.

    if len(option_list_1):

        bowtie2.options.pop('U', None)
        bowtie2.options.pop('1', None)
        bowtie2.options.pop('2', None)

        bowtie2.add_option_short(key='1', value=','.join(option_list_1))

        if len(option_list_2):
            bowtie2.add_option_short(key='2', value=','.join(option_list_2))

        sam_file_path = os.path.join(
            runnable.temporary_directory_path(absolute=False),
            'paired_fastq_files.sam')
        bowtie2.stdout_path = sam_file_path
        sam_file_path_list.append(sam_file_path)

        sys.stdout.writelines(bowtie2.trace(level=1))

        child_return_code = bowtie2.run()
        if child_return_code:
            raise Exception('Could not complete the ' + bowtie2.name + ' step.')

            # Do not delete the non-temporary FASTQ files.

    # TODO: At this stage, a list of SAM files exists, which needs merging.
    if len(sam_file_path_list) > 1:
        # SAM files need merging.

        java_process = bsf.process.Executable(name='sam_to_fastq', program='java', sub_command=bsf.process.Command())
        java_process.add_switch_short(key='d64')
        java_process.add_switch_short(key='server')
        java_process.add_switch_short(key='Xmx2G')
        java_process.add_option_pair(key='-Djava.io.tmpdir',
                                     value=runnable.temporary_directory_path(absolute=False))

        picard_process = java_process.sub_command
        picard_process.add_option_short(key='jar', value=os.path.join(
            bsf.standards.JavaClassPath.get_picard(), 'picard.jar'))
        picard_process.sub_command = bsf.process.Command(program='SamToFastq')

        sam_to_fastq = picard_process.sub_command
        # sam_to_fastq.add_option_pair(key='INPUT', value=bam_file_path)  # TODO:
        sam_to_fastq.add_option_pair(key='OUTPUT_PER_RG', value='true')
        sam_to_fastq.add_option_pair(key='OUTPUT_DIR', value=runnable.temporary_directory_path(absolute=False))
        sam_to_fastq.add_option_pair(key='INCLUDE_NON_PF_READS', value='false')  # TODO: Make this configurable.
        sam_to_fastq.add_option_pair(key='TMP_DIR', value=runnable.temporary_directory_path(absolute=False))
        sam_to_fastq.add_option_pair(key='VERBOSITY', value='WARNING')
        sam_to_fastq.add_option_pair(key='QUIET', value='false')
        sam_to_fastq.add_option_pair(key='VALIDATION_STRINGENCY', value='STRICT')


def run(runnable):
    """Run the C{bsf.procedure.ConsecutiveRunnable}.

    @param runnable: C{bsf.procedure.ConsecutiveRunnable}
    @type runnable: bsf.procedure.ConsecutiveRunnable
    @return:
    @rtype:
    """

    if not os.path.isdir(runnable.temporary_directory_path(absolute=False)):
        try:
            os.makedirs(runnable.temporary_directory_path(absolute=False))
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    file_path_read_group = runnable.file_path_object
    """ @type file_path_read_group: bsf.analyses.aligner.FilePathAlign """

    if not os.path.isdir(file_path_read_group.output_directory):
        try:
            os.makedirs(file_path_read_group.output_directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Run all Executable objects of this ConsecutiveRunnable.

    run_bowtie2(runnable=runnable)

    # Remove the temporary directory and everything within it.

    shutil.rmtree(path=runnable.temporary_directory_path(absolute=False), ignore_errors=False)

    # Job done.
