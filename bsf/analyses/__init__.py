"""bsf.analyses

A package of classes and methods supporting NGS-specific analyses such as ChIP-Seq or RNA-Seq.
"""

#
# Copyright 2013 - 2018 Michael K. Schuster
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


from __future__ import print_function

import os
import re
import sys

import bsf
import bsf.ngs
import bsf.process
import bsf.standards


class RunBamToFastq(bsf.Analysis):
    """BAM or SAM to FASTQ converter sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    """

    name = 'BAM To FastQ'
    prefix = 'bam_to_fastq'

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python C{str} for setting C{bsf.process.Executable.dependencies} in other C{bsf.Analysis} objects

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: The dependency string for a C{bsf.process.Executable} of this C{bsf.Analysis}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    def __init__(
            self,
            configuration=None,
            project_name=None,
            genome_version=None,
            input_directory=None,
            output_directory=None,
            project_directory=None,
            genome_directory=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None):
        """Initialise a RunBamToFastq.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: BSF Collection
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @return:
        @rtype:
        """

        super(RunBamToFastq, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            e_mail=e_mail,
            debug=debug,
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list)

        # Nothing else to do for this sub-class ...

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a BSF RunBamToFastq via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(RunBamToFastq, self).set_configuration(configuration=configuration, section=section)

        # Nothing else to do for this sub-class ...

        return

    def run(self):
        """Run a RunBamToFastq C{bsf.Analysis}.

        @return:
        @rtype:
        """

        super(RunBamToFastq, self).run()

        # Replicates have to be un-grouped, always!
        # replicate_grouping = config_parser.getboolean(section=config_section, option='replicate_grouping')
        replicate_grouping = False

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())

        for sample in self.collection.get_all_samples():
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping)
            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    # In a BSF Paired Reads object, the SAM or BAM file could potentially
                    # occur as reads1 or reads2 instance variable.

                    if paired_reads.reads_1 is not None:
                        file_name = paired_reads.reads_1.file_path
                        file_name = file_name.rstrip('/ ')
                        file_name = os.path.basename(file_name)

                        # TODO: The matching part to remove the .bam could be achieved with Bash parameter expansion.
                        match = re.search(pattern=r'(.*)\.bam$', string=file_name)
                        if match:
                            executable_read_group = stage_read_group.add_executable(
                                executable=bsf.process.Executable(
                                    name='picard_sam_to_fastq_{}_1'.format(paired_reads_name),
                                    program='bsf_bam2fastq.sh'))
                            self.set_command_configuration(command=executable_read_group)
                            executable_read_group.arguments.append(paired_reads.reads_1.file_path)
                            executable_read_group.arguments.append(
                                os.path.join(bsf.standards.JavaClassPath.get_picard(), 'picard.jar'))
                            executable_read_group.arguments.append(os.path.join(self.genome_directory, match.group(1)))

                    if paired_reads.reads_2 is not None:
                        file_name = paired_reads.reads_2.file_path
                        file_name = file_name.rstrip('/ ')
                        file_name = os.path.basename(file_name)

                        match = re.search(pattern=r'(.*)\.bam$', string=file_name)
                        if match:
                            executable_read_group = stage_read_group.add_executable(
                                executable=bsf.process.Executable(
                                    name='picard_sam_to_fastq_{}_2'.format(paired_reads_name),
                                    program='bsf_bam2fastq.sh'))
                            self.set_command_configuration(command=executable_read_group)
                            executable_read_group.arguments.append(paired_reads.reads_2.file_path)
                            executable_read_group.arguments.append(
                                os.path.join(bsf.standards.JavaClassPath.get_picard(), 'picard.jar'))
                            executable_read_group.arguments.append(os.path.join(self.genome_directory, match.group(1)))

        return
