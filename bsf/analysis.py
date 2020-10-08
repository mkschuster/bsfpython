# -*- coding: utf-8 -*-
"""Analysis module.

A package of classes and methods modelling analyses.
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

import datetime
import errno
import getpass
import html
import importlib
import inspect
import os
import sys
import urllib.parse
import uuid
import warnings

from bsf.ngs import Collection, Sample
from bsf.procedure import Runnable, ConcurrentRunnable, ConsecutiveRunnable
from bsf.process import Command, Executable, RunnableStep
from bsf.standards import Configuration, StandardFilePath, Operator, Genome, Transcriptome, UCSC, URL


class Analysis(object):
    """The C{bsf.analysis.Analysis} class represents a high-level analysis.

    It consists of one or more C{bsf.analysis.Stage} objects that may run one or more
    C{bsf.process.Executable} or C{bsf.process.RunnableStep} objects (programs).

    Attributes:
    @cvar name: C{bsf.analysis.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.analysis.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar ucsc_name_hub: UCSC Genome Browser Track Hub "hub" file name
    @type ucsc_name_hub: str
    @cvar ucsc_name_genomes: UCSC Genome Browser Track Hub "genomes" file name
    @type ucsc_name_genomes: str
    @cvar ucsc_name_tracks: UCSC Genome Browser Track Hub "tracks" file name
    @type ucsc_name_tracks: str
    @ivar configuration: C{bsf.standards.Configuration}
    @type configuration: Configuration
    @ivar debug: Debug level
    @type debug: int
    @ivar project_name: Project name (arbitrary)
    @type project_name: str | None
    @ivar genome_version: Genome version (e.g. hg19, mm10, GRCh37, GRCm38, ...)
    @type genome_version: str | None
    @ivar cache_directory: C{bsf.analysis.Analysis}-wide cache directory
    @type cache_directory: str | None
    @ivar input_directory: Input directory
    @type input_directory: str | None
    @ivar output_directory: Output directory, user-specified including a genome version sub-directory
    @type output_directory: str | None
    @ivar project_directory: Project-specific directory
    @type project_directory: str | None
    @ivar genome_directory: Genome-specific directory
    @type genome_directory: str | None
    @ivar sas_file: Sample Annotation Sheet (SAS) file path
    @type sas_file: str | None
    @ivar sas_prefix: A prefix to columns in a Sample Annotation Sheet
        (e.g. Control Sample, Treatment Sample, ...)
    @type sas_prefix: str | None
    @ivar e_mail: e-Mail address for a UCSC Genome Browser Track Hub
    @type e_mail: str | None
    @ivar stage_list: Python C{list} of C{bsf.analysis.Stage} objects
    @type stage_list: list[Stage]
    @ivar runnable_dict: Python C{dict} of Python C{str} (C{bsf.procedure.Runnable.name}) key data and
        C{bsf.procedure.Runnable} value data
    @type runnable_dict: dict[Runnable.name, Runnable]
    @ivar collection: C{bsf.ngs.Collection}
    @type collection: Collection
    @ivar sample_list: Python C{list} of C{bsf.ngs.Sample} objects
    @type sample_list: list[Sample]
    """

    name = 'Analysis'
    prefix = 'analysis'

    ucsc_name_hub = 'hub.txt'
    ucsc_name_genomes = 'genomes.txt'
    ucsc_name_tracks = 'tracks.txt'

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{bsf.analysis.Analysis} from a UNIX-style configuration file path.

        The configuration file in C{bsf.standards.Configuration.global_file_path} is read as default,
        before the project-specific one gets read, if it is not the same file.

        @param config_path: UNIX-style configuration file path
        @type config_path: str
        @return: C{bsf.analysis.Analysis}
        @rtype: Analysis
        """
        return cls.from_configuration(
            configuration=Configuration.from_file_path_list(
                file_path_list=[Configuration.global_file_path, config_path]))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{bsf.analysis.Analysis} from a C{bsf.standards.Configuration}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @return: C{bsf.analysis.Analysis}
        @rtype: Analysis
        """
        analysis = cls(configuration=configuration)

        # A "module.class" configuration section specifies defaults for this Analysis or sub-class
        # i.e. "bsf.analysis.Analysis" or "bsf.analyses.*", respectively.

        analysis.set_configuration(
            configuration=analysis.configuration,
            section=Configuration.section_from_instance(instance=analysis))

        return analysis

    def __init__(
            self,
            configuration=None,
            project_name=None,
            genome_version=None,
            cache_directory=None,
            input_directory=None,
            output_directory=None,
            project_directory=None,
            genome_directory=None,
            sas_file=None,
            sas_prefix=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            runnable_dict=None,
            collection=None,
            sample_list=None):
        """Initialise a C{bsf.analysis.Analysis}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration | None
        @param project_name: Project name
        @type project_name: str | None
        @param genome_version: Genome version
        @type genome_version: str | None
        @param cache_directory: C{bsf.analysis.Analysis}-wide cache directory
        @type cache_directory: str | None
        @param input_directory: C{bsf.analysis.Analysis}-wide input directory
        @type input_directory: str | None
        @param output_directory: C{bsf.analysis.Analysis}-wide output directory
        @type output_directory: str | None
        @param project_directory: C{bsf.analysis.Analysis}-wide project directory,
            normally under the C{bsf.analysis.Analysis}-wide output directory
        @type project_directory: str | None
        @param genome_directory: C{bsf.analysis.Analysis}-wide genome directory,
            normally under the C{bsf.analysis.Analysis}-wide project directory
        @type genome_directory: str | None
        @param sas_file: Sample Annotation Sheet (SAS) file path
        @type sas_file: str | None
        @param sas_prefix: A prefix to columns in a Sample Annotation Sheet
            (e.g. Control Sample, Treatment Sample, ...)
        @type sas_prefix: str | None
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str | None
        @param debug: Integer debugging level
        @type debug: int | None
        @param stage_list: Python C{list} of C{bsf.analysis.Stage} objects
        @type stage_list: list[Stage] | None
        @param runnable_dict: Python C{dict} of Python C{str} (C{bsf.procedure.Runnable.name}) and
            C{bsf.procedure.Runnable} value data
        @type runnable_dict: dict[Runnable.name, Runnable] | None
        @param collection: C{bsf.ngs.Collection}
        @type collection: Collection | None
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[Sample] | None
        """
        super(Analysis, self).__init__()

        if configuration is None:
            self.configuration = Configuration()
        else:
            self.configuration = configuration

        self.project_name = project_name
        self.genome_version = genome_version
        self.cache_directory = cache_directory
        self.input_directory = input_directory
        self.output_directory = output_directory
        self.project_directory = project_directory
        self.genome_directory = genome_directory
        self.sas_file = sas_file
        self.sas_prefix = sas_prefix
        self.e_mail = e_mail
        self.debug = debug

        if stage_list is None:
            self.stage_list = list()
        else:
            self.stage_list = stage_list

        if runnable_dict is None:
            self.runnable_dict = dict()
        else:
            self.runnable_dict = runnable_dict

        if collection is None:
            self.collection = Collection()
        else:
            self.collection = collection

        if sample_list is None:
            self.sample_list = list()
        else:
            self.sample_list = sample_list

        return

    def trace(self, level):
        """Trace a C{bsf.analysis.Analysis}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  project_name: {!r}\n'.format(indent, self.project_name))
        str_list.append('{}  genome_version: {!r}\n'.format(indent, self.genome_version))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  input_directory: {!r}\n'.format(indent, self.input_directory))
        str_list.append('{}  output_directory: {!r}\n'.format(indent, self.output_directory))
        str_list.append('{}  genome_directory: {!r}\n'.format(indent, self.genome_directory))
        str_list.append('{}  sas_file: {!r}\n'.format(indent, self.sas_file))
        str_list.append('{}  sas_prefix: {!r}\n'.format(indent, self.sas_prefix))
        str_list.append('{}  e_mail: {!r}\n'.format(indent, self.e_mail))
        str_list.append('{}  debug: {!r}\n'.format(indent, self.debug))
        str_list.append('{}  stage_list: {!r}\n'.format(indent, self.stage_list))
        str_list.append('{}  runnable_dict: {!r}\n'.format(indent, self.runnable_dict))
        str_list.append('{}  collection: {!r}\n'.format(indent, self.collection))
        str_list.append('{}  sample_list: {!r}\n'.format(indent, self.sample_list))

        str_list.append('{}  Python dict of Runnable objects:\n'.format(indent))
        for runnable_name in sorted(self.runnable_dict):
            str_list.append('{}    Key: {!r} Runnable: {!r}\n'.format(
                indent, runnable_name, self.runnable_dict[runnable_name]))
            runnable = self.runnable_dict[runnable_name]
            str_list.extend(runnable.trace(level=level + 2))

        str_list.append('{}  Python List of Sample objects:\n'.format(indent))
        for sample in self.sample_list:
            str_list.append('{}    Sample name: {!r} file_path: {!r}\n'.format(indent, sample.name, sample.file_path))

        if self.collection:
            str_list.extend(self.collection.trace(level=level + 1))

        return str_list

    def add_stage(self, stage):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.analysis.Stage}.

        If the C{bsf.analysis.Stage} exists already in the C{bsf.analysis.Analysis.stage_list} the method returns the
        already existing C{bsf.analysis.Stage}.

        @param stage: C{bsf.analysis.Stage}
        @type stage: Stage
        @return: C{bsf.analysis.Stage}
        @rtype: Stage
        """
        if stage not in self.stage_list:
            self.stage_list.append(stage)

        return stage

    def add_runnable(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.procedure.Runnable}.

        @param runnable: C{bsf.procedure.Runnable}
        @type runnable: Runnable
        @return: C{bsf.procedure.Runnable}
        @rtype: Runnable
        @raise Exception: The C{bsf.procedure.Runnable.name} already exists in the
            C{bsf.analysis.Analysis.runnable_dict}
        """
        if runnable.name in self.runnable_dict:
            raise Exception('A Runnable with name ' + repr(runnable.name) +
                            ' already exists in Analysis ' + repr(self.project_name) + '.')
        else:
            self.runnable_dict[runnable.name] = runnable

        return runnable

    def add_runnable_concurrent(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.procedure.ConcurrentRunnable}.

        @param runnable: C{bsf.procedure.ConcurrentRunnable}
        @type runnable: ConcurrentRunnable
        @return: C{bsf.procedure.ConcurrentRunnable}
        @rtype: ConcurrentRunnable
        @raise Exception: The C{bsf.procedure.Runnable.name} already exists in the
            C{bsf.analysis.Analysis.runnable_dict}
        """
        return self.add_runnable(runnable=runnable)

    def add_runnable_consecutive(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.procedure.ConsecutiveRunnable}.

        @param runnable: C{bsf.procedure.ConsecutiveRunnable}
        @type runnable: ConsecutiveRunnable
        @return: C{bsf.procedure.ConsecutiveRunnable}
        @rtype: ConsecutiveRunnable
        @raise Exception: The C{bsf.procedure.Runnable.name} already exists in the
            C{bsf.analysis.Analysis.runnable_dict}
        """
        return self.add_runnable(runnable=runnable)

    def add_sample(self, sample):
        """Add a C{bsf.ngs.Sample} to the Python C{list} of C{bsf.ngs.Sample} objects.

        If the C{bsf.ngs.Sample} already exists in the C{bsf.analysis.Analysis}, the method just returns.
        The check is based on the Python 'in' comparison operator and in lack of a specific
        __cmp__ method, relies on object identity (i.e. address).

        @param sample: C{bsf.ngs.Sample}
        @type sample: Sample
        """
        if sample not in self.sample_list:
            self.sample_list.append(sample)

        return

    def get_annotation_file(self, prefix_list, suffix):
        """Get a project and genome-specific annotation file.

        Based on the project name, a list of file name prefixes and one file name suffix,
        the file name that exists in the file system will be returned.
        @param prefix_list: Python C{list} of Python C{str} (prefix) objects
        @type prefix_list: list[str] | None
        @param suffix: File name suffix
        @type suffix: str
        @return: File name or None
        @rtype: str | None
        """
        if prefix_list is None:
            prefix_list = [self.prefix]

        for prefix in prefix_list:
            # Preferentially test with the genome version.
            file_name = '_'.join((self.project_name, self.genome_version, prefix, suffix))
            if self.debug > 0:
                print('Checking annotation sheet: ', repr(file_name))
            if os.path.exists(file_name):
                return file_name

            # Fall-back test without the genome version.
            file_name = '_'.join((self.project_name, prefix, suffix))
            if self.debug > 0:
                print('Checking annotation sheet: ', repr(file_name))
            if os.path.exists(file_name):
                return file_name

        return

    def get_stage(self, name):
        """Get a C{bsf.analysis.Stage} from a C{bsf.analysis.Analysis}.

        If the C{bsf.analysis.Stage} does not exist, it is created and initialised via the
        C{bsf.standards.Configuration} in C{bsf.analysis.Analysis.configuration}.
        Reads from configuration file sections
        I{[bsf.analysis.Stage]}
        I{[bsf.analysis.Analysis.Stage]} or I{[bsf.analyses.*.Stage]}
        I{[bsf.analysis.Analysis.Stage.name]} or I{[bsf.analyses.*.Stage.name]}

        @param name: Name
        @type name: str
        @return: C{bsf.analysis.Stage}
        @rtype: Stage
        """
        # Check if a Stage with this the name already exists and if so, return it.
        for stage in self.stage_list:
            if stage.name == name:
                return stage

        # Initialise a new Stage and add it to the Python list of Stage objects.

        stage = Stage(name=name, working_directory=self.genome_directory)
        self.stage_list.append(stage)

        # A "bsf.analysis.Stage" section specifies defaults for all Stage objects of an Analysis.

        section = Configuration.section_from_instance(instance=stage)
        stage.set_configuration(configuration=self.configuration, section=section)

        if self.debug > 1:
            print('Stage configuration section:', repr(section))

        # A "bsf.analysis.Analysis.Stage" or "bsf.analyses.*.Stage" pseudo-class section specifies
        # Analysis-specific or sub-class-specific options for the Stage, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=self), 'Stage'))
        stage.set_configuration(configuration=self.configuration, section=section)

        if self.debug > 1:
            print('Stage configuration section:', repr(section))

        # A "bsf.analysis.Analysis.Stage.name" or "bsf.analyses.*.Stage.name" section specifies defaults
        # for a particular Stage of an Analysis or sub-class, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=self), 'Stage', stage.name))
        stage.set_configuration(configuration=self.configuration, section=section)

        if self.debug > 1:
            print('Stage configuration section:', repr(section))

        return stage

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analysis.Analysis} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @raise Exception: The specified section does not exist
        """
        if not configuration.config_parser.has_section(section=section):
            raise Exception(
                'Section ' + repr(section) +
                ' not defined in Configuration files:\n' + repr(configuration.file_path_list))

        # The configuration section is available.

        option = 'debug'
        if configuration.config_parser.has_option(section=section, option=option):
            self.debug = configuration.config_parser.getint(section=section, option=option)

        option = 'project_name'
        if configuration.config_parser.has_option(section=section, option=option):
            self.project_name = configuration.config_parser.get(section=section, option=option)

        option = 'cache_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.cache_directory = configuration.config_parser.get(section=section, option=option)

        option = 'input_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.input_directory = configuration.config_parser.get(section=section, option=option)

        option = 'output_directory'
        if configuration.config_parser.has_option(section=section, option=option):
            self.output_directory = configuration.config_parser.get(section=section, option=option)

        option = 'genome_version'
        if configuration.config_parser.has_option(section=section, option=option):
            self.genome_version = configuration.config_parser.get(section=section, option=option)

        option = 'sas_file'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sas_file = configuration.config_parser.get(section=section, option=option)

        option = 'sas_prefix'
        if configuration.config_parser.has_option(section=section, option=option):
            self.sas_prefix = configuration.config_parser.get(section=section, option=option)

        option = 'e_mail'
        if configuration.config_parser.has_option(section=section, option=option):
            self.e_mail = configuration.config_parser.get(section=section, option=option)

        return

    def set_command_configuration(self, command):
        """Set default C{bsf.argument.Argument} objects for a C{bsf.process.Command}.

        @param command: C{bsf.process.Command}
        @type command: Command
        """
        # TODO: Phase out this method.
        # Once all accessory scripts are converted to Runnable and bsf.process.RunnableStep objects,
        # this method becomes redundant.
        section = Configuration.section_from_instance(instance=command)

        # For plain bsf.process.Executable objects append the value of the
        # bsf.process.Executable.program to make the configuration section more meaningful.

        if section == 'bsf.process.Executable':
            section += '.'
            section += command.program

        if self.debug > 1:
            print('Command configuration section:', repr(section))

        command.set_configuration(configuration=self.configuration, section=section)

        return

    def set_runnable_step_configuration(self, runnable_step):
        """Set default C{bsf.argument.Argument} objects for a C{bsf.process.RunnableStep}.

        This method reads configuration section(s)
        "Analysis.__class__.__name__"."RunnableStep.name"[."Command.name"]*
        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: RunnableStep
        """

        def _set_configuration(command, section):
            """Recursively set default C{bsf.argument.Argument} objects for a C{bsf.process.RunnableStep}.

            The method sets the default C{bsf.argument.Argument} objects for a C{bsf.process.RunnableStep},
            as well as for its contained sub C{bsf.process.Command} objects.
            @param command: C{bsf.process.Command}
            @type command: Command
            @param section: Configuration section
            @type section: str
            """
            if command.name:
                section += '.' + command.name
            else:
                section += '.' + command.program

            if self.debug > 1:
                print('RunnableStep configuration section:', repr(section))

            command.set_configuration(configuration=self.configuration, section=section)

            if command.sub_command is not None:
                _set_configuration(command=command.sub_command, section=section)

            return

        # Initially, the configuration section prefix is based on the Analysis class name and the Analysis Stage.name.

        _set_configuration(command=runnable_step, section=self.configuration.section_from_instance(instance=self))

        return

    def set_stage_runnable(self, stage, runnable):
        """Create a C{bsf.process.Executable} to assign a C{bsf.procedure.Runnable} to a C{bsf.analysis.Stage}.

        In case the file in C{bsf.procedure.Runnable.get_relative_status_path} exists already,
        C{bsf.process.Executable.submit} will be set to C{False}.

        @param stage: C{bsf.analysis.Stage}
        @type stage: Stage
        @param runnable: C{bsf.procedure.Runnable}
        @type runnable: Runnable
        @return: C{bsf.process.Executable}
        @rtype: Executable
        @raise Exception: A C{bsf.procedure.Runnable.name} does not exist in C{bsf.analysis.Analysis.runnable_dict}
        @raise Exception: A C{bsf.analysis.Stage} does not exist in C{bsf.analysis.Analysis.stage_list}
        """
        if stage not in self.stage_list:
            raise Exception('A Stage with name ' + repr(stage.name) +
                            ' does not exist in the Analysis with name ' + repr(self.project_name) + '.')

        if runnable.name not in self.runnable_dict:
            raise Exception('A Runnable with name ' + repr(runnable.name) +
                            ' does not exist in the Analysis with name ' + repr(self.project_name) + '.')

        if stage.template_script:
            if os.path.isabs(stage.template_script):
                template_script = stage.template_script
            else:
                template_script = os.path.join(StandardFilePath.get_template_scripts(), stage.template_script)

            executable = Executable(name=runnable.name, program=template_script)
            executable.arguments.append(runnable.pickler_path)
        else:
            executable = Executable(name=runnable.name, program=Runnable.runner_script)
            executable.add_option_long(key='pickler-path', value=runnable.pickler_path)

        # Only submit the bsf.process.Executable if the status file does not exist already.
        if os.path.exists(runnable.runnable_status_file_path(success=True, absolute=True)):
            executable.submit = False

        stage.add_executable(executable=executable)

        return executable

    def run(self):
        """Run a C{bsf.analysis.Analysis}.

        @raise Exception: An C{bsf.analysis.Analysis.project_name} has not been defined
        """
        if self.debug is None:
            self.debug = 0

        if not self.project_name:
            raise Exception('A ' + self.name + " requires a 'project_name' configuration option.")

        # Some analyses such as FastQC do not require a genome_version,
        # nor a genome_version-specific output directory.
        # Also, add the e-mail address for UCSC track hubs into the genome subclass.

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard directory paths.

        self.cache_directory = self.configuration.get_absolute_path(
            file_path=self.cache_directory,
            default_path=StandardFilePath.get_cache())

        self.input_directory = self.configuration.get_absolute_path(
            file_path=self.input_directory,
            default_path=StandardFilePath.get_samples(absolute=True))

        self.output_directory = self.configuration.get_absolute_path(
            file_path=self.output_directory,
            default_path=StandardFilePath.get_projects(absolute=True))

        # As a safety measure, to prevent creation of rogue directory paths, the output_directory has to exist.

        if not os.path.isdir(self.output_directory):
            raise Exception('The Analysis output_directory ' + repr(self.output_directory) + ' does not exist.')

        # Define project_directory and genome_directory instance variables.
        # If a genome_version option is present, append
        # it to the project_directory instance variable.
        # This allows analyses run against more than one directory and
        # simplifies UCSC Genome Browser track hub creation.

        # FIXME: The Analysis(project_directory) option is ignored.
        # Change the project_directory instance variable into a private (i.e. _project_directory) instance variable.
        # Then, access via @property get_project_directory.
        self.project_directory = os.path.join(self.output_directory, self.project_name)

        # FIXME: The Analysis(genome_directory) option is ignored.
        # Change the genome_directory instance variable into a private (i.e. _genome_directory) instance variable.
        # Then, access via @property get_genome_directory.
        if self.genome_version:
            self.genome_directory = os.path.join(self.project_directory, self.genome_version)
        else:
            self.genome_directory = self.project_directory

        if not os.path.isdir(self.genome_directory):
            try:
                os.makedirs(self.genome_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        if not self.e_mail:
            self.e_mail = Operator.get_e_mail()
            if not self.e_mail:
                raise Exception('A ' + self.name + " requires an 'e_mail' configuration option.")

        if self.sas_file:
            # Populate a bsf.ngs.Collection from a SampleAnnotationSheet.
            self.sas_file = self.configuration.get_absolute_path(file_path=self.sas_file)

            # NOTE: Do no longer search for configuration files and sample annotation sheets in the project directory.
            # if not os.path.isabs(self.sas_file) and not os.path.exists(self.sas_file):
            #     self.sas_file = os.path.join(self.project_directory, self.sas_file)

            self.collection = Collection.from_sas_path(
                file_path=self.input_directory,
                file_type='Automatic',
                name=self.project_name,
                sas_path=self.sas_file,
                sas_prefix=self.sas_prefix)

            if self.debug > 1:
                print(self, 'Collection name:', repr(self.collection.name))
                sys.stdout.writelines(self.collection.trace(level=1))
        else:
            # Create an empty bsf.ngs.Collection.
            self.collection = Collection()

        return

    def report(self):
        """Create a C{bsf.analysis.Analysis} report.

        The method must be implemented in a sub-class.
        """
        warnings.warn(
            "The 'report' method must be implemented in the sub-class.",
            UserWarning)

        return

    # @staticmethod
    # def escape_url(url):
    #     url = url.replace(' ', '%20')
    #     url = url.replace('&', '%26')
    #     url = url.replace('+', '%2B')
    #     url = url.replace(';', '%3B')
    #     url = url.replace('=', '%3D')
    #     url = url.replace('?', '%3F')
    #     return url
    #
    # @staticmethod
    # def escape_html(html):
    #     html = html.replace('&', '&amp;')
    #     html = html.replace('<', '&lt;')
    #     html = html.replace('>', '&gt;')
    #     return html
    #
    # @staticmethod
    # def escape_space(text):
    #     return text.replace(' ', '%20')

    @staticmethod
    def get_html_anchor(prefix, suffix, text):
        """Create a HTML anchor element with a relative reference path.

        <a href="prefix/prefix_suffix">text</a>
        @param prefix: Prefix
        @type prefix: str
        @param suffix: Suffix
        @type suffix: str
        @param text: Link text
        @type text: str
        @return: HTML anchor element
        @rtype: str
        """
        return '<a href="' + prefix + '/' + prefix + '_' + suffix + '">' + text + '</a>'

    @staticmethod
    def get_html_image(prefix, suffix, text, height=None, width=None):
        """Create a HTML image element with a relative source path.

        <img alt="text" src="prefix/prefix_suffix" height="80" width="80" />
        @param prefix: Prefix
        @type prefix: str
        @param suffix: Suffix
        @type suffix: str
        @param text: Alternative text
        @type text: str
        @param height: Image height attribute
        @type height: str
        @param width: Image width attribute
        @type width: str
        @return: HTML image element
        @rtype: str
        """
        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('<img')
        str_list.append('alt="' + text + '"')
        str_list.append('src="' + prefix + '/' + prefix + '_' + suffix + '"')

        if height:  # not None and not empty
            str_list.append('height="' + height + '"')

        if width:  # not None and not empty
            str_list.append('width="' + width + '"')

        str_list.append('/>')

        return ' '.join(str_list)

    @staticmethod
    def get_html_genome(genome_version=None):
        """Get a genome description HTML paragraph.

        @param genome_version: Genome version
        @type genome_version: str | None
        @return: Genome description
        @rtype: list[str]
        """
        str_list = list()
        """ @type str_list: list[str] """

        if genome_version is None:
            return str_list

        species = Genome.get_species(genome_version=genome_version)
        # Without species information, the description is not useful.
        if species is None:
            warnings.warn("No species information for genome version '" + genome_version +
                          "' in the central configuration file.")
            return str_list

        str_list.append('<p>')
        str_list.append('<strong>Genome:</strong> ')
        str_list.append('<i>' + species + '</i> ')
        str_list.append('genome assembly ')
        str_list.append(genome_version)

        date = Genome.get_date(genome_version=genome_version)
        if date is not None:
            str_list.append(' (' + date + ')')

        str_list.append('</p>\n')

        return str_list

    @staticmethod
    def get_html_transcriptome(transcriptome_version=None):
        """Get a transcriptome description HTML paragraph.

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: str | None
        @return: Transcriptome description
        @rtype: list[str]
        """
        str_list = list()
        """ @type str_list: list[str] """

        if transcriptome_version is None:
            return str_list

        species = Transcriptome.get_species(transcriptome_version=transcriptome_version)
        # Without species information, the description is not useful.
        if species is None:
            warnings.warn("No species information for transcriptome version '" + transcriptome_version +
                          "' in the central configuration file.")
            return str_list

        str_list.append('<p>')
        str_list.append('<strong>Transcriptome:</strong> ')
        str_list.append('<i>' + species + '</i> ')
        str_list.append('transcriptome ')
        str_list.append(transcriptome_version)

        date = Transcriptome.get_date(transcriptome_version=transcriptome_version)
        if date is not None:
            str_list.append(' (' + date + ')')

        str_list.append('</p>\n')

        return str_list

    def get_html_header(
            self,
            strict=True,
            creator=None,
            source=None,
            title=None):
        """Get the header section of a XHTML 1.0 document.

        @param strict: XHTML 1.0 Strict or XHTML 1.0 Transitional Document Type Declaration,
            defaults to XHTML 1.0 Strict
        @type strict: bool
        @param creator: Dublin Core DC.Creator meta field value,
            defaults to the USER environment variable
        @type creator: str
        @param source: Dublin Core DC.Source meta field value,
            defaults to the Python script file path
        @type source: str
        @param title: Title element value,
            defaults to a concatenation of
            C{bsf.analysis.Analysis.project_name} and
            C{bsf.analysis.Analysis.report_name}
        @type title: str
        @return: XHTML 1.0 header section as Python C{list} of Python C{str} objects
        @rtype: list[str]
        """
        if creator is None or not creator:
            creator = getpass.getuser()
            # The getpass.getuser method just relies on environment variables,
            # but at least works under Unix and Windows.

        if source is None or not source:
            source = inspect.getfile(inspect.currentframe())

        if title is None or not title:
            title = ' '.join((self.project_name, self.name))

        str_list = list()
        """ @type str_list: list[str] """

        if strict:
            str_list.append('<!DOCTYPE html PUBLIC ')
            str_list.append('"-//W3C//DTD XHTML 1.0 Strict//EN" ')
            str_list.append('"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n')
        else:
            str_list.append('<!DOCTYPE html PUBLIC ')
            str_list.append('"-//W3C//DTD XHTML 1.0 Transitional//EN" ')
            str_list.append('"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')

        str_list.append('\n')
        str_list.append('<html xmlns="http://www.w3.org/1999/xhtml">\n')
        str_list.append('<head>\n')
        str_list.append('<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />\n')
        # str_list.append('<link rel="stylesheet" href="/' +
        #                 urllib.parse.quote(string=URL.get_relative_projects()) +
        #                 '/bsfpython.css" type="text/css" />\n')
        str_list.append('<link rel="schema.DC" href="http://purl.org/DC/elements/1.0/" />\n')
        str_list.append('<meta name="DC.Creator" content="' + html.escape(s=creator, quote=True) + '" />\n')
        str_list.append('<meta name="DC.Date" content="' + datetime.datetime.now().isoformat() + '" />\n')
        str_list.append('<meta name="DC.Source" content="' + html.escape(s=source, quote=True) + '" />\n')
        str_list.append('<meta name="DC.Title" content="' + html.escape(s=title, quote=True) + '" />\n')
        str_list.append('<style type="text/css">\n')
        str_list.append('  .left    { text-align: left; }\n')
        str_list.append('  .right   { text-align: right; }\n')
        str_list.append('  .center  { text-align: center; }\n')
        str_list.append('  .justify { text-align: justify; }\n')
        str_list.append('  .start   { text-align: start; }\n')
        str_list.append('  .end     { text-align: end; }\n')
        str_list.append('  body { font-family: sans-serif; }\n')
        str_list.append('  h1 { color: #40B9D4; }\n')
        str_list.append('  a { color: #40B9D4; text-decoration: none; }\n')
        str_list.append('  a:hover, a:focus { color: #2a6496; text-decoration: underline; }\n')
        str_list.append('</style>\n')
        str_list.append('<title>' + html.escape(s=title, quote=True) + '</title>\n')
        str_list.append('</head>\n')
        str_list.append('\n')
        str_list.append('<body>\n')
        str_list.append('\n')

        return str_list

    def get_html_footer(
            self,
            contact=None,
            institution=None,
            url_protocol=None,
            url_host_name=None,
            title=None):
        """Get the footer section of a XHTML 1.0 document.

        @param contact: Institution contact e-mail address,
            defaults to C{bsf.standards.Operator.get_contact()}
        @type contact: str
        @param institution: Institution name to be inserted into 'This report was generated by ...',
            defaults to C{bsf.standards.Operator.get_institution()}
        @type institution: str
        @param url_protocol: The protocol section of the institution URL (i.e. http, https, ...),
            defaults to C{bsf.standards.URL.get_protocol()}
        @type url_protocol: str
        @param url_host_name: The host name section of the institution URL (e.g. biomedical-sequencing.at),
            defaults to C{bsf.standards.URL.get_host_name()}
        @type url_host_name: str
        @param title: Title element value,
            defaults to a concatenation of
            C{bsf.analysis.Analysis.project_name} and
            C{bsf.analysis.Analysis.report_name}
        @type title: str
        @return: XHTML 1.0 footer section as Python C{list} of Python C{str} objects
        @rtype: list[str]
        """
        if not contact:
            contact = Operator.get_contact()

        if not institution:
            institution = Operator.get_institution()

        if not url_protocol:
            url_protocol = URL.get_protocol()

        if not url_host_name:
            url_host_name = URL.get_host_name()

        if not title:
            title = ' '.join((self.project_name, self.name))

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('<hr class="footer" />\n')
        str_list.append('<p class="footer">\n')
        if institution:
            str_list.append('This report was generated by ' + html.escape(s=institution, quote=True) + '\n')

        # Method bsf.standards.URL.get_absolute_base() would also return the default URL.
        if url_protocol:
            str_list.append('<a href="' + url_protocol + '://' +
                            url_host_name + '/">' + html.escape(s=url_host_name, quote=True) + '</a>.\n')
        else:
            str_list.append('<a href="//' +
                            url_host_name + '/">' + html.escape(s=url_host_name, quote=True) + '</a>.\n')

        str_list.append('<br class="footer" />\n')
        if contact:
            # The e-mail address outside of an URL still needs HTML quoting.
            # After URL quoting, nothing critical needing HTML escaping should be left.
            str_list.append(
                'Contact: <a href="mailto:' + urllib.parse.quote(string=contact) +
                '?subject=' + urllib.parse.quote(string=title) + '">')
            str_list.append(html.escape(s=contact, quote=True))
            str_list.append('</a>')
        str_list.append('</p>\n')
        str_list.append('</body>\n')
        str_list.append('</html>\n')
        str_list.append('\n')

        return str_list

    def get_html_report(
            self,
            content,
            strict=True,
            title=None,
            creator=None,
            source=None,
            contact=None,
            institution=None,
            url_protocol=None,
            url_host_name=None):
        """Get a report as a XHTML 1.0 document.

        The method automatically concatenates the
        XHTML header C{bsf.analysis.Analysis.get_html_header}, the XHTML content and the
        XHTML footer C{bsf.analysis.Analysis.get_html_footer} before returning the report.
        @param content: XHTML 1.0 content
        @type content: list[str]
        @param strict: XHTML 1.0 Strict or XHTML 1.0 Transitional Document Type Declaration,
            defaults to XHTML 1.0 Strict
        @type strict: bool
        @param creator: Dublin Core DC.Creator meta field value,
            defaults to the USER environment variable
        @type creator: str
        @param source: Dublin Core DC.Source meta field value,
            defaults to the Python script file path
        @type source: str
        @param title: Title element value,
            defaults to a concatenation of
            C{bsf.analysis.Analysis.project_name} and
            C{bsf.analysis.Analysis.report_name}
        @type title: str
        @param contact: Institution contact e-mail address,
            defaults to C{bsf.standards.Operator.get_contact()}
        @type contact: str
        @param institution: Institution name to be inserted into 'This report was generated by ...',
            defaults to C{bsf.standards.Operator.get_institution()}
        @type institution: str
        @param url_protocol: The protocol section of the institution URL (i.e. http, https, ...),
            defaults to C{bsf.standards.URL.get_protocol()}
        @type url_protocol: str
        @param url_host_name: The host name section of the institution URL (e.g. biomedical-sequencing.at),
            defaults to C{bsf.standards.URL.get_host_name()}
        @type url_host_name: str
        @return: XHTML 1.0 report as Python C{list} of Python C{str} objects
        @rtype: list[str]
        """
        str_list = list()
        """ @type str_list: list[str] """

        str_list.extend(
            self.get_html_header(
                strict=strict,
                creator=creator,
                source=source,
                title=title))
        str_list.extend(content)
        str_list.extend(
            self.get_html_footer(
                contact=contact,
                institution=institution,
                url_protocol=url_protocol,
                url_host_name=url_host_name,
                title=title))

        return str_list

    def report_to_file(
            self,
            content,
            prefix=None,
            strict=True,
            title=None,
            creator=None,
            source=None,
            contact=None,
            institution=None,
            url_protocol=None,
            url_host_name=None):
        """Write a XHTML 1.0 report I{prefix_report.html} file into the C{bsf.analysis.Analysis.genome_directory}.

        The method automatically concatenates the
        XHTML header C{bsf.analysis.Analysis.get_html_header}, the XHTML content and the
        XHTML footer C{bsf.analysis.Analysis.get_html_footer} before writing the file.
        @param content: XHTML 1.0 content
        @type content: list[str]
        @param prefix: A file name prefix (e.g. chipseq, rnaseq, ...), defaults to C{bsf.analysis.Analysis.prefix}
        @type prefix: str
        @param strict: XHTML 1.0 Strict or XHTML 1.0 Transitional Document Type Declaration,
            defaults to XHTML 1.0 Strict
        @type strict: bool
        @param creator: Dublin Core DC.Creator meta field value,
            defaults to the USER environment variable
        @type creator: str
        @param source: Dublin Core DC.Source meta field value,
            defaults to the Python script file path
        @type source: str
        @param title: Title element value,
            defaults to a concatenation of
            C{bsf.analysis.Analysis.project_name} and
            C{bsf.analysis.Analysis.report_name}
        @type title: str
        @param contact: Institution contact e-mail address,
            defaults to C{bsf.standards.Operator.get_contact()}
        @type contact: str
        @param institution: Institution name to be inserted into 'This report was generated by ...',
            defaults to C{bsf.standards.Operator.get_institution()}
        @type institution: str
        @param url_protocol: The protocol section of the institution URL (i.e. http, https, ...),
            defaults to C{bsf.standards.URL.get_protocol()}
        @type url_protocol: str
        @param url_host_name: The host name section of the institution URL (e.g. biomedical-sequencing.at),
            defaults to C{bsf.standards.URL.get_host_name()}
        @type url_host_name: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        with open(
                file=os.path.join(self.genome_directory, '_'.join((prefix, 'report.html'))),
                mode='wt') as output_file:
            output_file.writelines(
                self.get_html_report(
                    content=content,
                    strict=strict,
                    title=title,
                    creator=creator,
                    source=source,
                    contact=contact,
                    institution=institution,
                    url_protocol=url_protocol,
                    url_host_name=url_host_name))

        return

    def create_project_genome_directory(self):
        """Check for and create a project or genome directory if necessary.

        C{bsf.analysis.Analysis.project_directory}
        C{bsf.analysis.Analysis.genome_directory}
        @raise Exception: Output (genome) directory does not exist
        """
        if not os.path.isdir(self.genome_directory):
            answer = input(
                'Output (genome) directory ' + repr(self.genome_directory) + ' does not exist.\n' +
                'Create? [Y/n] ')

            if not answer or answer == 'Y' or answer == 'y':
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.makedirs(self.genome_directory)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
            else:
                raise Exception(
                    'Output (genome) directory ' + repr(self.genome_directory) + ' does not exist.')

        return

    def create_public_project_link(self, sub_directory=None):
        """Create a symbolic link from the public HTML directory to the project directory if not already there.

        The link will be placed in the specified sub directory and contain
        the project name followed by a 128 bit hexadecimal UUID string.
        If not specified, the sub directory defaults to the value of C{bsf.standards.URL.get_relative_projects()}.

        @param sub_directory: C{bsf.analysis.Analysis}-specific directory
        @type sub_directory: str
        @return: Symbolic link to the project directory
        @rtype: str
        @raise Exception: Public HTML path does not exist
        """
        if sub_directory is None:
            sub_directory = URL.get_relative_projects()

        html_path = os.path.join(StandardFilePath.get_public_html(absolute=True), sub_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the html_path directory has to exist.

        if not os.path.isdir(html_path):
            raise Exception(
                'The public HTML directory path ' + repr(html_path) + ' does not exist.\n' +
                'Please check the optional sub-directory name ' + repr(sub_directory) + '.')

        # The target_path consists of the absolute public_html directory,
        # the analysis-specific sub-directory, the project name and a 128 bit hexadecimal UUID string.

        target_path = os.path.join(html_path, '_'.join((self.project_name, uuid.uuid4().hex)))

        # The final_path holds the final symbolic link target path.
        # It can be the target_path assembled above or point to an already existing one.

        final_path = target_path

        for file_name in os.listdir(html_path):
            file_path = os.path.join(html_path, file_name)
            if os.path.islink(file_path):
                source_path = os.readlink(file_path)
                if not os.path.isabs(source_path):
                    source_path = os.path.join(html_path, source_path)
                source_path = os.path.normpath(source_path)
                if not os.path.exists(source_path):
                    # Both paths for os.path.samefile have to exist.
                    # Check for dangling symbolic links.
                    warnings.warn(
                        'Dangling symbolic link ' + repr(source_path) + ' to ' + repr(file_path),
                        UserWarning)
                    continue
                if os.path.samefile(source_path, self.project_directory):
                    # Reset the final_path to the already existing file_path.
                    final_path = file_path
                    # Do not break out here to discover all dangling symbolic links.

        if final_path != target_path:
            # A symbolic link already exists.
            # Ask the user to re-create the symbolic link.
            answer = input(
                'Public HTML link ' + repr(self.project_directory) + ' to ' + repr(final_path) + ' exists.\n' +
                'Re-create? [y/N] ')

            if not answer or answer.upper() == 'N':
                print('Public HTML link ' + repr(self.project_directory) + ' to ' + repr(final_path) + ' not reset.')
            else:
                try:
                    os.remove(final_path)
                except OSError as exception:
                    if exception.errno != errno.ENOENT:
                        raise
                try:
                    os.symlink(os.path.relpath(self.project_directory, html_path), target_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
                final_path = target_path
        else:
            # A symbolic link does not exist.
            # Ask the user to create a symbolic link.
            answer = input(
                'Public HTML link ' + repr(self.project_directory) + ' to ' + repr(final_path) + ' does not exist.\n' +
                'Create? [Y/n] ')

            if not answer or answer.upper() == 'Y':
                try:
                    os.symlink(os.path.relpath(self.project_directory, html_path), final_path)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
            else:
                print('Public HTML link ' + repr(self.project_directory) + ' to ' + repr(final_path) + ' not set.')

        return final_path

    def ucsc_track_url(
            self,
            options_dict=None,
            browser_dict=None,
            track_dict=None,
            ucsc_protocol=None,
            ucsc_host_name=None):
        """Return a URL to automatically attach a UCSC Genome Browser track.

        @param options_dict: Python C{dict} of Python C{str} URL option key value pairs
        @type options_dict: dict[str, str]
        @param browser_dict: Python C{dict} of Python C{str} browser line key value pairs
        @type browser_dict: dict[str, str]
        @param track_dict: Python C{dict} of Python C{str} track line (hgct_customText) key value pairs
        @type track_dict: dict[str, str]
        @param ucsc_protocol: UCSC Genome Browser URL protocol (i.e. http, https, ...)
            defaults to the value of C{bsf.standards.UCSC.get_protocol()}
        @type ucsc_protocol: str
        @param ucsc_host_name: UCSC Genome Browser URL host name,
            defaults to the value of C{bsf.standards.UCSC.get_host_name()}
        @type ucsc_host_name: str
        @return: A URL to attach a track to the UCSC Genome Browser
        @rtype: str
        """
        if options_dict is None:
            options_dict = dict()

        if 'db' not in options_dict:
            options_dict['db'] = Genome.resolve_ucsc_alias(genome_version=self.genome_version)

        # UCSC "browser" configuration dictionary.

        if browser_dict:
            pass

        # UCSC "track" configuration dictionary.

        if track_dict:
            options_dict['hgct_customText'] = 'track'
            for key in sorted(track_dict):
                options_dict['hgct_customText'] += ' ' + key + '=' + track_dict[key]

        # UCSC protocol

        if not ucsc_protocol:
            ucsc_protocol = UCSC.get_protocol()

        # UCSC host name

        if not ucsc_host_name:
            ucsc_host_name = UCSC.get_host_name()

        # Strip leading colons to support protocol-independent URLs.
        return html.escape(
            s='{}://{}/cgi-bin/hgTracks?{}'.format(
                urllib.parse.quote(string=ucsc_protocol),
                urllib.parse.quote(string=ucsc_host_name),
                urllib.parse.urlencode(query=options_dict)).lstrip(':'),
            quote=True)

    def ucsc_hub_url(self, link_path, options_dict=None):
        """Return a URL to automatically attach a UCSC Genome Browser Track Hub.

        @param link_path: Symbolic link path in the public HTML directory including project name and a UUID
        @type link_path: str
        @param options_dict: Python C{dict} of Python C{str} URL option key value pairs
        @type options_dict: dict[str, str]
        @return: A URL to automatically attach a UCSC Genome Browser Track Hub
        @rtype: str
        """
        if options_dict is None:
            options_dict = dict()

        if 'hubUrl' not in options_dict:
            # The track hub URL requires the link name, i.e. the link path base name, to be inserted.
            link_name = os.path.basename(link_path.rstrip('/'))
            options_dict['hubUrl'] = '/'.join((
                URL.get_absolute_projects(),
                link_name,
                '_'.join((self.prefix, self.ucsc_name_hub))))

        return self.ucsc_track_url(options_dict=options_dict)

    def ucsc_hub_html_anchor(self, link_path):
        """Return a XHTML 1.0 anchor element to automatically attach a UCSC Genome Browser Track Hub.

        @param link_path: Symbolic link path in the public HTML directory including project name and a UUID
        @type link_path: str
        @see: C{bsf.analysis.Analysis.create_public_project_link}
        @return: XHTML 1.0 anchor element as Python C{list} of Python C{str} objects
        @rtype: list[str]
        """
        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('UCSC Genome Browser Track Hub ')
        str_list.append('<a href="' + self.ucsc_hub_url(link_path=link_path) + '">' + self.project_name + '</a>')

        return str_list

    def ucsc_hub_write_hub(self, prefix=None):
        """Write a UCSC Track Hub I{prefix_hub.txt} file into the C{bsf.analysis.Analysis.project_directory}.

        The C{bsf.analysis.Analysis.project_directory} is one level above the C{bsf.analysis.Analysis.genome_directory}.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('hub ' + '_'.join((self.project_name, prefix)) + '\n')
        str_list.append('shortLabel ' + '_'.join((self.project_name, prefix)) + '\n')
        str_list.append('longLabel Project ' + '_'.join((self.project_name, prefix)) + '\n')
        str_list.append('genomesFile ' + '_'.join((prefix, self.ucsc_name_genomes)) + '\n')
        str_list.append('email ' + self.e_mail + '\n')

        with open(
                file=os.path.join(self.project_directory, '_'.join((prefix, self.ucsc_name_hub))),
                mode='wt') as output_file:
            output_file.writelines(str_list)

        return

    def ucsc_hub_write_genomes(self, prefix=None):
        """Write a UCSC Track Hub I{prefix_genomes.txt} file into the C{bsf.analysis.Analysis.project_directory}.

        The C{bsf.analysis.Analysis.project_directory} is one level above the C{bsf.analysis.Analysis.genome_directory}.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        # If the file exists, read it first to retain any other genome assembly entries.
        genome_version_dict = dict()
        """ @type genome_version_dict: dict[str, str] """

        file_path = os.path.join(self.project_directory, '_'.join((prefix, self.ucsc_name_genomes)))

        if os.path.exists(file_path):
            genome_version = None
            with open(file=file_path, mode='rt') as input_file:
                for line_str in input_file:
                    line_str = line_str.strip()
                    if not line_str:
                        continue
                    field_list = line_str.split()
                    if len(field_list) != 2:
                        warnings.warn('Malformed line ' + repr(line_str) + ' in UCSC genomes file ' +
                                      repr(file_path) + '.\n' +
                                      'Expected exactly two components after line splitting.')
                    if field_list[0] == 'genome':
                        if genome_version is not None:
                            warnings.warn('Malformed line ' + repr(line_str) + ' in UCSC genomes file ' +
                                          repr(file_path) + '.\n' +
                                          "Got more than one 'genomes' lines in succession.")
                        genome_version = field_list[1]
                    if field_list[0] == 'trackDb':
                        if genome_version is None:
                            warnings.warn('Malformed line ' + repr(line_str) + ' in UCSC genomes file ' +
                                          repr(file_path) + '.\n' +
                                          "Got a 'trackDb' line without a preceding genomes line.")
                        else:
                            genome_version_dict[genome_version] = field_list[1]
                            genome_version = None

        # Resolve an eventual alias for the UCSC genome assembly name in "genome_version/prefix_tracks.txt".
        genome_version_dict[Genome.resolve_ucsc_alias(genome_version=self.genome_version)] = \
            '/'.join((self.genome_version, '_'.join((prefix, self.ucsc_name_tracks))))

        str_list = list()
        """ @type str_list: list[str] """

        for genome_version in sorted(genome_version_dict):
            str_list.append('genome ' + genome_version + '\n')
            str_list.append('trackDb ' + genome_version_dict[genome_version] + '\n')
            str_list.append('\n')

        with open(file=file_path, mode='wt') as output_file:
            output_file.writelines(str_list)

        return

    def ucsc_hub_write_tracks(self, content, prefix=None):
        """Write a UCSC Track Hub I{prefix_tracks.txt} file into the C{bsf.analysis.Analysis.genome_directory}.

        @param content: Content
        @type content: list[str]
        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        with open(
                file=os.path.join(self.genome_directory, '_'.join((prefix, self.ucsc_name_tracks))),
                mode='wt') as output_file:
            output_file.writelines(content)

        return

    def ucsc_hub_to_file(self, content, prefix=None):
        """Write UCSC Genome Browser Track Hub files to disk.

        The method writes a I{prefix_hub.txt} and a I{prefix_genomes.txt} file into the
        C{bsf.analysis.Analysis.project_directory}, above the C{bsf.analysis.Analysis.genome_directory}, as well as a
        I{prefix_tracks.txt} file into the C{bsf.analysis.Analysis.genome_directory}.

        @param content: Content of the track database file
        @type content: list[str]
        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...), defaults to C{bsf.analysis.Analysis.prefix}
        @type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        self.ucsc_hub_write_hub(prefix=prefix)
        self.ucsc_hub_write_genomes(prefix=prefix)
        self.ucsc_hub_write_tracks(content=content, prefix=prefix)

        return

    def check_state(self):
        """Check the state of each C{bsf.analysis.Stage}.
        """
        for stage in self.stage_list:
            stage.check_state(debug=self.debug)

        return

    def submit(self, name=None):
        """Submit each C{bsf.analysis.Stage}.

        Submits each C{bsf.process.Executable} of either all C{bsf.analysis.Stage} objects or a named one and pickles
        each C{bsf.procedure.Runnable}.

        @param name: Only submit C{bsf.process.Executable} objects linked to C{bsf.analysis.Stage.name}
        @type name: str
        """
        # Pickle all bsf.procedure.Runnable objects.

        for runnable in self.runnable_dict.values():
            runnable.to_pickler_path()

        # Submit all bsf.process.Executable objects of all bsf.analysis.Stage objects.

        submit = 0

        for stage in self.stage_list:
            if name:
                if name == stage.name:
                    submit += 1
                else:
                    continue
            stage.submit(debug=self.debug)

            if self.debug:
                print(repr(stage))
                sys.stdout.writelines(stage.trace(level=1))

        if name:
            if name == 'report':
                self.report()
            elif not submit:
                name_list = [stage.name for stage in self.stage_list]
                name_list.append('report')
                print('Valid Analysis Stage names are:', repr(name_list))

        return


class Stage(object):
    """The C{bsf.analysis.Stage} class represents a stage of a C{bsf.analysis.Analysis}.

    A C{bsf.analysis.Stage} represents C{bsf.process.Executable} or C{bsf.procedure.Runnable} objects that share
    similar resource requirements of a I{Distributed Resource Management System} (I{DRMS}).

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar working_directory: Working directory path
    @type working_directory: str
    @ivar implementation: Implementation (e.g. I{sge}, I{slurm}, ...)
    @type implementation: str
    @ivar memory_free_mem: Memory limit (free physical)
    @type memory_free_mem: str | None
    @ivar memory_free_swap: Memory limit (free swap)
    @type memory_free_swap: str | None
    @ivar memory_free_virtual: Memory limit (free virtual)
    @type memory_free_virtual: str | None
    @ivar memory_limit_hard: Memory limit (hard)
    @type memory_limit_hard: str | None
    @ivar memory_limit_soft: Memory limit (soft)
    @type memory_limit_soft: str | None
    @ivar node_list_exclude: List of nodes to exclude
    @type node_list_exclude: list[str] | None
    @ivar node_list_include: List of nodes to include
    @type node_list_include: list[str] | None
    @ivar time_limit: Time limit
    @type time_limit: str | None
    @ivar parallel_environment: Parallel environment
    @type parallel_environment: str | None
    @ivar queue: Queue
    @type queue: str | None
    @ivar threads: Number of threads
    @type threads: int
    @ivar hold: Hold on job scheduling
    @type hold: str | None
    @ivar is_script: C{bsf.process.Executable} objects represent shell scripts,
        or alternatively binary programs
    @type is_script: bool
    @ivar executable_list: Python C{list} of C{bsf.process.Executable} objects
    @type executable_list: list[Executable]
    """

    def __init__(
            self,
            name,
            working_directory,
            implementation=None,
            memory_free_mem=None,
            memory_free_swap=None,
            memory_free_virtual=None,
            memory_limit_hard=None,
            memory_limit_soft=None,
            node_list_exclude=None,
            node_list_include=None,
            time_limit=None,
            parallel_environment=None,
            queue=None,
            threads=1,
            hold=None,
            is_script=False,
            template_script=None,
            executable_list=None):
        """Initialise a C{bsf.analysis.Stage}.

        @param name: Name
        @type name: str
        @param working_directory: Working directory
        @type working_directory: str
        @param implementation: Implementation (e.g. I{sge}, I{slurm}, ...)
        @type implementation: str
        @param memory_free_mem: Memory limit (free physical)
        @type memory_free_mem: str | None
        @param memory_free_swap: Memory limit (free swap)
        @type memory_free_swap: str | None
        @param memory_free_virtual: Memory limit (free virtual)
        @type memory_free_virtual: str | None
        @param memory_limit_hard: Memory limit (hard)
        @type memory_limit_hard: str | None
        @param memory_limit_soft: Memory limit (soft)
        @type memory_limit_soft: str | None
        @param node_list_exclude: List of nodes to exclude
        @type node_list_exclude: list[str] | None
        @param node_list_include: List of nodes to include
        @type node_list_include: list[str] | None
        @param time_limit: Time limit
        @type time_limit: str | None
        @param parallel_environment: Parallel environment
        @type parallel_environment: str | None
        @param queue: Queue
        @type queue: str | None
        @param threads: Number of threads
        @type threads: int
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param is_script: C{bsf.process.Executable} objects represent shell scripts,
            or alternatively binary programs
        @type is_script: bool
        @param template_script: Template script for submission.
        @type template_script: str
        @param executable_list: Python C{list} of C{bsf.process.Executable} objects
        @type executable_list: list[Executable]
        """
        super(Stage, self).__init__()

        if name is None:
            self.name = str()
        else:
            self.name = name

        if working_directory is None:
            self.working_directory = str()
        else:
            self.working_directory = working_directory

        if implementation is None:
            self.implementation = str()
        else:
            self.implementation = implementation

        self.memory_free_mem = memory_free_mem
        self.memory_free_swap = memory_free_swap
        self.memory_free_virtual = memory_free_virtual
        self.memory_limit_hard = memory_limit_hard
        self.memory_limit_soft = memory_limit_soft
        self.node_list_exclude = node_list_exclude
        self.node_list_include = node_list_include
        self.time_limit = time_limit
        self.parallel_environment = parallel_environment
        self.queue = queue

        if threads is None:
            self.threads = 1
        else:
            assert isinstance(threads, int)
            self.threads = threads

        self.hold = hold
        # FIXME: Does not seem to be used. Remove!

        if is_script is None:
            self.is_script = False
        else:
            assert isinstance(is_script, bool)
            self.is_script = is_script

        self.template_script = template_script

        if executable_list is None:
            self.executable_list = list()
        else:
            self.executable_list = executable_list

        return

    def trace(self, level):
        """Trace a C{bsf.analysis.Stage}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  name:                 {!r}\n'.format(indent, self.name))
        str_list.append('{}  working_directory:    {!r}\n'.format(indent, self.working_directory))
        str_list.append('{}  implementation:       {!r}\n'.format(indent, self.implementation))
        str_list.append('{}  memory_free_mem:      {!r}\n'.format(indent, self.memory_free_mem))
        str_list.append('{}  memory_free_swap:     {!r}\n'.format(indent, self.memory_free_swap))
        str_list.append('{}  memory_free_virtual:  {!r}\n'.format(indent, self.memory_free_virtual))
        str_list.append('{}  memory_limit_hard:    {!r}\n'.format(indent, self.memory_limit_hard))
        str_list.append('{}  memory_limit_soft:    {!r}\n'.format(indent, self.memory_limit_soft))
        str_list.append('{}  node_list_exclude:    {!r}\n'.format(indent, self.node_list_exclude))
        str_list.append('{}  node_list_include:    {!r}\n'.format(indent, self.node_list_include))
        str_list.append('{}  time_limit:           {!r}\n'.format(indent, self.time_limit))
        str_list.append('{}  queue:                {!r}\n'.format(indent, self.queue))
        str_list.append('{}  parallel_environment: {!r}\n'.format(indent, self.parallel_environment))
        str_list.append('{}  threads:              {!r}\n'.format(indent, self.threads))
        str_list.append('{}  hold:                 {!r}\n'.format(indent, self.hold))
        str_list.append('{}  is_script:            {!r}\n'.format(indent, self.is_script))
        str_list.append('{}  template_script:      {!r}\n'.format(indent, self.template_script))

        str_list.append('{}  executable_list:\n'.format(indent))

        for executable in self.executable_list:
            str_list.extend(executable.trace(level=level + 2))

        return str_list

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analysis.Stage} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """
        if not configuration.config_parser.has_section(section=section):
            raise Exception(
                'Section ' + repr(section) + ' not defined in Configuration files:\n' +
                repr(configuration.file_path_list))

        # The configuration section is available.

        option = 'hold'
        if configuration.config_parser.has_option(section=section, option=option):
            self.hold = configuration.config_parser.get(section=section, option=option)

        option = 'implementation'
        if configuration.config_parser.has_option(section=section, option=option):
            self.implementation = configuration.config_parser.get(section=section, option=option)

        option = 'is_script'
        if configuration.config_parser.has_option(section=section, option=option):
            self.is_script = configuration.config_parser.getboolean(section=section, option=option)

        option = 'memory_free_mem'
        if configuration.config_parser.has_option(section=section, option=option):
            self.memory_free_mem = configuration.config_parser.get(section=section, option=option)

        option = 'memory_free_swap'
        if configuration.config_parser.has_option(section=section, option=option):
            self.memory_free_swap = configuration.config_parser.get(section=section, option=option)

        option = 'memory_free_virtual'
        if configuration.config_parser.has_option(section=section, option=option):
            self.memory_free_virtual = configuration.config_parser.get(section=section, option=option)

        option = 'memory_hard'
        if configuration.config_parser.has_option(section=section, option=option):
            self.memory_limit_hard = configuration.config_parser.get(section=section, option=option)

        option = 'memory_soft'
        if configuration.config_parser.has_option(section=section, option=option):
            self.memory_limit_soft = configuration.config_parser.get(section=section, option=option)

        option = 'node_list_exclude'
        if configuration.config_parser.has_option(section=section, option=option):
            self.node_list_exclude = configuration.get_list_from_csv(section=section, option=option)

        option = 'node_list_include'
        if configuration.config_parser.has_option(section=section, option=option):
            self.node_list_include = configuration.get_list_from_csv(section=section, option=option)

        option = 'time_limit'
        if configuration.config_parser.has_option(section=section, option=option):
            self.time_limit = configuration.config_parser.get(section=section, option=option)

        option = 'parallel_environment'
        if configuration.config_parser.has_option(section=section, option=option):
            self.parallel_environment = configuration.config_parser.get(section=section, option=option)

        option = 'queue'
        if configuration.config_parser.has_option(section=section, option=option):
            self.queue = configuration.config_parser.get(section=section, option=option)

        option = 'threads'
        if configuration.config_parser.has_option(section=section, option=option):
            self.threads = configuration.config_parser.getint(section=section, option=option)

        option = 'template_script'
        if configuration.config_parser.has_option(section=section, option=option):
            self.template_script = configuration.config_parser.get(section=section, option=option)

        return

    def add_executable(self, executable):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.process.Executable}.

        @param executable: C{bsf.process.Executable}
        @type executable: Executable
        @return: C{bsf.process.Executable}
        @rtype: Executable
        """
        self.executable_list.append(executable)

        return executable

    def check_state(self, debug=0):
        """Check the state of each C{bsf.process.Executable}.

        @param debug: Debug level
        @type debug: int
        """
        # Dynamically import the module specific for the configured DRMS implementation.

        python_module = importlib.import_module(name='.'.join(('bsf', 'drms', self.implementation)))
        python_module.check_state(stage=self, debug=debug)

        return

    def submit(self, debug=0):
        """Submit a command line for each C{bsf.process.Executable}.

        @param debug: Debug level
        @type debug: int
        """
        # Dynamically import the module specific for the configured DRMS implementation.

        python_module = importlib.import_module(name='.'.join(('bsf', 'drms', self.implementation)))
        python_module.submit(stage=self, debug=debug)

        return
