"""bsf.analyses.ega

A package of classes and methods supporting the EGA Cryptor tool.
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


import os
import sys

import bsf
import bsf.process


class FilePathEGACryptorReadGroup(bsf.FilePath):
    """The C{bsf.analyses.ega.FilePathEGACryptorReadGroup} models read group-specific EGA Cryptor file paths.

    Attributes:
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.ega.FilePathEGACryptorReadGroup} object.

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathEGACryptorReadGroup, self).__init__(prefix=prefix)

        self.output_directory = prefix

        return


class EGACryptor(bsf.Analysis):
    """BSF EGA Cryptor C{bsf.Analysis} sub-class.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    """

    name = 'EGA Cryptor Analysis'
    prefix = 'ega_cryptor'

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
            sample_list=None,
            classpath_ega_cryptor=None):
        """Initialise a C{bsf.analyses.ega.EGACryptor}.

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
        @param classpath_ega_cryptor: EGA Cryptor tool Java Archive (JAR) class path directory
        @type classpath_ega_cryptor: str | unicode | None
        @return:
        @rtype:
        """
        super(EGACryptor, self).__init__(
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

        self.classpath_ega_cryptor = classpath_ega_cryptor

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.ega.EGACryptor} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(EGACryptor, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        # Get the EGA Cryptor tool Java Archive (JAR) class path directory.

        option = 'classpath_ega_cryptor'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_ega_cryptor = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a C{bsf.analyses.ega.EGACryptor} C{bsf.Analysis}.

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

        super(EGACryptor, self).run()

        run_read_comparisons()

        if not self.classpath_ega_cryptor:
            raise Exception('A ' + self.name + " analysis requires a 'classpath_ega_cryptor' configuration option.")

        stage_read_group = self.get_stage(name=self.get_stage_name_read_group())

        self.sample_list.sort(cmp=lambda x, y: cmp(x.name, y.name))

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', sample.name)
                sys.stdout.writelines(sample.trace(level=1))

            # The EGACryptor analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=replicate_grouping, exclude=False)
            for paired_reads_name in sorted(paired_reads_dict):
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print(self, 'PairedReads name:', paired_reads.get_name())

                    prefix_read_group = self.get_prefix_read_group(read_group_name=paired_reads_name)

                    file_path_read_group = FilePathEGACryptorReadGroup(
                        prefix=prefix_read_group)

                    # Create a Runnable and an Executable for running the EGA Cryptor analysis.

                    runnable_read_group = self.add_runnable(
                        runnable=bsf.Runnable(
                            name=prefix_read_group,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory,
                            file_path_object=file_path_read_group))
                    self.set_stage_runnable(stage=stage_read_group, runnable=runnable_read_group)

                    # Create a new RunnableStepMakeDirectory in preparation of the EGA Cryptor program.

                    runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepMakeDirectory(
                            name='mkdir',
                            directory_path=file_path_read_group.output_directory))

                    link_target_path = os.path.join(
                        file_path_read_group.output_directory,
                        os.path.basename(paired_reads.reads_1.file_path))
                    runnable_read_group.add_runnable_step(runnable_step=bsf.process.RunnableStepLink(
                        name='link',
                        source_path=paired_reads.reads_1.file_path,
                        target_path=link_target_path))

                    # Create a RunnableStep to run the Java-based EGA Cryptor.

                    runnable_step_read_group = runnable_read_group.add_runnable_step(
                        runnable_step=bsf.process.RunnableStepJava(
                            name='ega_cryptor',
                            java_temporary_path=runnable_read_group.temporary_directory_path(absolute=False),
                            java_heap_maximum='Xmx4G',
                            java_jar_path=self.classpath_ega_cryptor))

                    # Use a sequence of sub-Command objects to separate options that have to appear
                    # in a particular order. Sigh!

                    runnable_step_read_group.sub_command.sub_command = bsf.process.Command()
                    sub_command = runnable_step_read_group.sub_command.sub_command
                    sub_command.add_switch_short(key='p')

                    sub_command.sub_command = bsf.process.Command()
                    sub_command = sub_command.sub_command
                    sub_command.add_switch_short(key='fm')

                    sub_command.sub_command = bsf.process.Command()
                    sub_command = sub_command.sub_command
                    sub_command.add_switch_short(key='file')

                    sub_command.arguments.append(link_target_path)

        return
