# -*- coding: utf-8 -*-
"""FastQC Analysis module

A package of classes and methods supporting the FastQC tool.
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
import os
import sys
import urllib.parse

import bsf
import bsf.annotation
import bsf.executables
import bsf.ngs
import bsf.procedure


class FilePathFastQCReadGroup(bsf.procedure.FilePath):
    """The C{bsf.analyses.fastqc.FilePathFastQCReadGroup} models read group-specific FastQC file paths.

    Attributes:
    @ivar archive: GNU Zip compressed archive
    @type archive: str | unicode
    @ivar report: FastQC HTML report
    @type report: str | unicode
    """

    def __init__(self, prefix, file_prefix):
        """Initialise a C{bsf.analyses.fastqc.FilePathFastQCReadGroup} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @param file_prefix: File prefix
        @type file_prefix: str | unicode
        @return:
        @rtype:
        """
        super(FilePathFastQCReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        self.archive = os.path.join(prefix, file_prefix + '_fastqc.zip')
        self.report = os.path.join(prefix, file_prefix + '_fastqc.html')

        return


class FastQC(bsf.Analysis):
    """BSF FastQC-specific Quality Assessment C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    """

    name = 'FastQC Analysis'
    prefix = 'fastqc'

    @classmethod
    def get_stage_name_read_group(cls):
        """Get a particular C{bsf.Stage.name}.

        @return: C{bsf.Stage.name}
        @rtype: str
        """
        return '_'.join((cls.prefix, 'read_group'))

    @classmethod
    def get_prefix_read_group(cls, read_group_name):
        """Get a Python C{str} prefix representing a C{bsf.procedure.Runnable}.

        @param read_group_name: Read group name
        @type read_group_name: str
        @return: Python C{str} prefix representing a C{bsf.procedure.Runnable}
        @rtype: str
        """
        return '_'.join((cls.get_stage_name_read_group(), read_group_name))

    @classmethod
    def get_file_path_read_group(cls, read_group_name, file_path):
        """Get a C{FilePathFastQCReadGroup} object.

        @param read_group_name: Read group name
        @type read_group_name: str
        @param file_path: File path (*.bam, *.fastq, *. fastq.gz, ...)
        @type file_path: str
        @return: C{FilePathFastQCReadGroup} object
        @rtype: FilePathFastQCReadGroup
        """
        def get_file_prefix(_file_path):
            """Private function to isolate a file prefix from '*.bam', '*.fastq', '*.fastq.gz', ... file paths.

            @param _file_path: File path
            @type _file_path: str
            @return: File prefix
            @rtype: str | unicode
            """
            root_path, extension = os.path.splitext(os.path.basename(_file_path))

            if extension:
                return get_file_prefix(_file_path=root_path)
            else:
                return root_path

        return FilePathFastQCReadGroup(
            prefix=cls.get_prefix_read_group(read_group_name=read_group_name),
            file_prefix=get_file_prefix(_file_path=file_path))

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
        """Initialise a C{bsf.analyses.fastqc.FastQC}.

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
        @param stage_list: Python C{list} of BSF C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @return:
        @rtype:
        """

        super(FastQC, self).__init__(
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

        # Sub-class specific ...

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.fastqc.FastQC} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(FastQC, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        return

    def run(self):
        """Run a C{bsf.analyses.fastqc.FastQC} C{bsf.Analysis}.

        @return:
        @rtype:
        """

        # Always check each BSF PairedReads object separately.
        replicate_grouping = False

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        super(FastQC, self).run()

        run_read_comparisons()

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())

        # Sort the Python list of Sample objects by Sample.name.

        self.sample_list.sort(key=lambda item: item.name)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            # The FastQC analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping, exclude=False)
            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                    file_path_read_group = self.get_file_path_read_group(
                        read_group_name=paired_reads_name,
                        file_path=paired_reads.reads_1.file_path)

                    # Create a Runnable and an Executable for running the FastQC analysis.

                    runnable_read_group = self.add_runnable_consecutive(
                        runnable=bsf.procedure.ConsecutiveRunnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory))
                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                    # Create a new RunnableStepMakeDirectory in preparation of the FastQC program.

                    runnable_step = bsf.process.RunnableStepMakeDirectory(
                            name='mkdir',
                            directory_path=file_path_read_group.output_directory)
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    runnable_step = bsf.process.RunnableStep(
                            name=prefix_read_group,
                            program='fastqc')
                    runnable_read_group.add_runnable_step(runnable_step=runnable_step)

                    runnable_step.add_switch_long(key='quiet')
                    runnable_step.add_option_long(key='outdir', value=file_path_read_group.output_directory)
                    runnable_step.add_option_long(key='threads', value=str(stage_read_group.threads))

                    if sample.file_type == 'CASAVA' and 'casava' not in runnable_step.options:
                        runnable_step.add_switch_long(key='casava')

                    # Add the Reads file paths.

                    if paired_reads.reads_1 is not None:
                        runnable_step.arguments.append(paired_reads.reads_1.file_path)

                    if paired_reads.reads_2 is not None:
                        runnable_step.arguments.append(paired_reads.reads_2.file_path)

        return

    def report(self):
        """Create a C{bsf.analyses.fastqc.FastQC} report in HTML format.

        @return:
        @rtype:
        """
        # Always check each BSF PairedReads object separately.
        replicate_grouping = False

        # Create a symbolic link containing the project name and a UUID.
        self.create_public_project_link()

        # Write a HTML document.

        report_list = list()
        """ @type report_list: list[str | unicode] """

        report_list.append('<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n')
        report_list.append('\n')

        report_list.append('<p>\n')
        report_list.append('<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a>\n')
        report_list.append('is a quality control tool for high throughput sequence data.\n')
        report_list.append('</p>\n')
        report_list.append('\n')

        report_list.append('<table id="fastq_report_table">\n')
        report_list.append('\n')

        report_list.append('<thead>\n')
        report_list.append('<tr>\n')
        report_list.append('<th>Sample</th>\n')
        report_list.append('<th>FastQC Report</th>\n')
        report_list.append('<th>FastQC Archive</th>\n')
        report_list.append('</tr>\n')
        report_list.append('</thead>\n')
        report_list.append('<tbody>\n')

        for sample in self.sample_list:
            # The FastQC analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(
                replicate_grouping=replicate_grouping,
                exclude=False,
                full=True)

            if not paired_reads_dict:
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue

            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    file_path_read_group = self.get_file_path_read_group(
                        read_group_name=paired_reads_name,
                        file_path=paired_reads.reads_1.file_path)

                    report_list.append('<tr>\n')
                    report_list.append('<td>')
                    report_list.append(paired_reads_name)
                    report_list.append('</td>\n')
                    report_list.append('<td>')
                    report_list.append('<a href="' + urllib.parse.quote(string=file_path_read_group.report) + '">')
                    report_list.append('<strong><abbr title="Hypertext Markup Language">HTML</abbr></strong>')
                    report_list.append('</a>')
                    report_list.append('</td>\n')
                    report_list.append('<td>')
                    report_list.append('<a href="' + urllib.parse.quote(string=file_path_read_group.archive) + '">')
                    report_list.append('<abbr title="PKWARE, Inc. ZIP">ZIP</abbr>')
                    report_list.append('</a>')
                    report_list.append('</td>\n')
                    report_list.append('</tr>\n')
                    report_list.append('\n')

        report_list.append('</tbody>\n')
        report_list.append('</table>\n')
        report_list.append('\n')

        self.report_to_file(content=report_list)

        return
