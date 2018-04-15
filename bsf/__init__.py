"""bsf

A package of classes and methods specific to the Biomedical Sequencing Facility (BSF).
Reference: http://www.biomedical-sequencing.at/
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

import cgi
import datetime
import errno
import getpass
import importlib
import inspect
import os
import pickle
import urllib
import uuid
import warnings

from bsf.ngs import Collection, Sample
from bsf.process import Command, Executable, RunnableStep
from bsf.standards import Configuration, Default, Operator


def _comma_separated_to_list(value_string):
    return filter(lambda x: x != '', map(lambda x: x.strip(), value_string.split(',')))


class Analysis(object):
    """The C{bsf.Analysis} class represents a high-level analysis.

    It consists of one or more C{bsf.Stage} objects that may run one or more
    C{bsf.process.Executable} or C{bsf.process.Runnable} objects (programs).

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar ucsc_name_hub: UCSC Genome Browser Track Hub "hub" file name
    @type ucsc_name_hub: str
    @cvar ucsc_name_genomes: UCSC Genome Browser Track Hub "genomes" file name
    @type ucsc_name_genomes: str
    @cvar ucsc_name_tracks: UCSC Genome Browser Track Hub "tracks" file name
    @type ucsc_name_tracks: str
    @ivar configuration: C{bsf.standards.Configuration}
    @type configuration: bsf.standards.Configuration
    @ivar debug: Debug level
    @type debug: int
    @ivar project_name: Project name (arbitrary)
    @type project_name: str
    @ivar genome_version: Genome version (e.g. hg19, mm10, GRCh37, GRCm38, ...)
    @type genome_version: str
    @ivar input_directory: Input directory
    @type input_directory: str | unicode
    @ivar output_directory: Output directory, user-specified including a genome version sub-directory
    @type output_directory: str | unicode
    @ivar project_directory: Project-specific directory
    @type project_directory: str | unicode
    @ivar genome_directory: Genome-specific directory
    @type genome_directory: str | unicode
    @ivar stage_list: Python C{list} of C{bsf.Stage} objects
    @type stage_list: list[bsf.Stage]
    @ivar runnable_dict: Python C{dict} of Python C{str} (C{bsf.Runnable.name}) key data and C{bsf.Runnable} value data
    @type runnable_dict: dict[bsf.Runnable.name, bsf.Runnable]
    @ivar collection: C{bsf.ngs.Collection}
    @type collection: bsf.ngs.Collection
    @ivar sample_list: Python C{list} of C{bsf.ngs.Sample} objects
    @type sample_list: list[bsf.ngs.Sample]
    """

    name = 'Analysis'
    prefix = 'analysis'

    ucsc_name_hub = 'hub.txt'
    ucsc_name_genomes = 'genomes.txt'
    ucsc_name_tracks = 'tracks.txt'

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{bsf.Analysis} from a UNIX-style configuration file path.

        The configuration file in C{bsf.standards.Configuration.global_file_path} is read as default,
        before the project-specific one gets read, if it is not the same file.

        @param config_path: UNIX-style configuration file path
        @type config_path: str | unicode
        @return: C{bsf.Analysis}
        @rtype: bsf.Analysis
        """

        return cls.from_configuration(
            configuration=Configuration.from_file_path_list(
                file_path_list=[Configuration.global_file_path, config_path]))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{bsf.Analysis} from a C{bsf.standards.Configuration}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @return: C{bsf.Analysis}
        @rtype: bsf.Analysis
        """

        assert isinstance(configuration, Configuration)

        analysis = cls(configuration=configuration)

        # A "module.class" configuration section specifies defaults for this Analysis or sub-class
        # i.e. "bsf.Analysis" or "bsf.analyses.*", respectively.

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
        """Initialise a C{bsf.Analysis}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param cache_directory: C{bsf.Analysis}-wide cache directory
        @type cache_directory: str
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
        @param sas_file: Sample Annotation Sheet (SAS) file path
        @type sas_file: str | unicode
        @param sas_prefix: A prefix to columns in a Sample Annotation Sheet
            (e.g. Control Sample, Treatment Sample, ...)
        @type sas_prefix: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param runnable_dict: Python C{dict} of Python C{str} (C{bsf.Runnable.name}) and C{bsf.Runnable} value data
        @type runnable_dict: dict[bsf.Runnable.name, bsf.Runnable]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @return:
        @rtype:
        """

        super(Analysis, self).__init__()

        if configuration is None:
            self.configuration = Configuration()
        else:
            assert isinstance(configuration, Configuration)
            self.configuration = configuration

        if project_name is None:
            self.project_name = str()
        else:
            self.project_name = project_name

        if genome_version is None:
            self.genome_version = str()
        else:
            self.genome_version = genome_version

        if cache_directory is None:
            self.cache_directory = str()
        else:
            self.cache_directory = cache_directory

        if input_directory is None:
            self.input_directory = str()
        else:
            self.input_directory = input_directory

        if output_directory is None:
            self.output_directory = str()
        else:
            self.output_directory = output_directory

        if project_directory is None:
            self.project_directory = str()
        else:
            self.project_directory = project_directory

        if genome_directory is None:
            self.genome_directory = str()
        else:
            self.genome_directory = genome_directory

        if sas_file is None:
            self.sas_file = str()
        else:
            self.sas_file = sas_file

        if sas_prefix is None:
            self.sas_prefix = str()
        else:
            self.sas_prefix = sas_prefix

        if e_mail is None:
            self.e_mail = str()
        else:
            self.e_mail = e_mail

        if debug is None:
            self.debug = int(x=0)
        else:
            assert isinstance(debug, int)
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
            assert isinstance(collection, Collection)
            self.collection = collection

        if sample_list is None:
            self.sample_list = list()
        else:
            self.sample_list = sample_list

        return

    def trace(self, level):
        """Trace a C{bsf.Analysis}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  project_name: {!r}\n'.format(indent, self.project_name)
        output += '{}  genome_version: {!r}\n'.format(indent, self.genome_version)
        output += '{}  cache_directory: {!r}\n'.format(indent, self.cache_directory)
        output += '{}  input_directory: {!r}\n'.format(indent, self.input_directory)
        output += '{}  output_directory: {!r}\n'.format(indent, self.output_directory)
        output += '{}  genome_directory: {!r}\n'.format(indent, self.genome_directory)
        output += '{}  sas_file: {!r}\n'.format(indent, self.sas_file)
        output += '{}  sas_prefix: {!r}\n'.format(indent, self.sas_prefix)
        output += '{}  e_mail: {!r}\n'.format(indent, self.e_mail)
        output += '{}  debug: {!r}\n'.format(indent, self.debug)
        output += '{}  stage_list: {!r}\n'.format(indent, self.stage_list)
        output += '{}  runnable_dict: {!r}\n'.format(indent, self.runnable_dict)
        output += '{}  collection: {!r}\n'.format(indent, self.collection)
        output += '{}  sample_list: {!r}\n'.format(indent, self.sample_list)

        output += '{}  Python dict of Runnable objects:\n'.format(indent)
        runnable_name_list = self.runnable_dict.keys()
        runnable_name_list.sort(cmp=lambda x, y: cmp(x, y))
        for runnable_name in runnable_name_list:
            output += '{}    Key: {!r} Runnable: {!r}\n'.format(
                indent, runnable_name, self.runnable_dict[runnable_name])
            runnable = self.runnable_dict[runnable_name]
            output += runnable.trace(level=level + 2)

        output += '{}  Python List of Sample objects:\n'.format(indent)
        for sample in self.sample_list:
            output += '{}    Sample name: {!r} file_path: {!r}\n'.format(indent, sample.name, sample.file_path)

        if self.collection:
            output += self.collection.trace(level + 1)

        return output

    def add_stage(self, stage):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.Stage}.

        If the C{bsf.Stage} exists already in the C{bsf.Analysis.stage_list} the method returns the
        already existing C{bsf.Stage}.

        @param stage: C{bsf.Stage}
        @type stage: bsf.Stage
        @return: C{bsf.Stage}
        @rtype: bsf.Stage
        """
        assert isinstance(stage, Stage)

        if stage not in self.stage_list:
            self.stage_list.append(stage)

        return stage

    def add_runnable(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.Runnable}.

        @param runnable: C{bsf.Runnable}
        @type runnable: bsf.Runnable
        @return: C{bsf.Runnable}
        @rtype: bsf.Runnable
        @raise Exception: The C{bsf.Runnable.name} already exists in the C{bsf.Analysis.runnable_dict}
        """
        assert isinstance(runnable, Runnable)

        if runnable.name in self.runnable_dict:
            raise Exception('A Runnable with name ' + repr(runnable.name) +
                            ' already exists in Analysis ' + repr(self.project_name) + '.')
        else:
            self.runnable_dict[runnable.name] = runnable

        return runnable

    def add_sample(self, sample):
        """Add a C{bsf.ngs.Sample} to the Python C{list} of C{bsf.ngs.Sample} objects.

        If the C{bsf.ngs.Sample} already exists in the C{bsf.Analysis}, the method just returns.
        The check is based on the Python 'in' comparison operator and in lack of a specific
        __cmp__ method, relies on object identity (i.e. address).

        @param sample: C{bsf.ngs.Sample}
        @type sample: bsf.ngs.Sample
        @return:
        @rtype:
        """
        assert isinstance(sample, Sample)

        if sample not in self.sample_list:
            self.sample_list.append(sample)

        return

    def get_stage(self, name):
        """Get a C{bsf.Stage} from a C{bsf.Analysis}.

        If the C{bsf.Stage} does not exist, it is created and initialised via the
        C{bsf.standards.Configuration} in C{bsf.Analysis.configuration}.
        Reads from configuration file sections
        I{[bsf.Stage]}
        I{[bsf.Analysis.Stage]} or I{[bsf.analyses.*.Stage]}
        I{[bsf.Analysis.Stage.name]} or I{[bsf.analyses.*.Stage.name]}

        @param name: Name
        @type name: str
        @return: C{bsf.Stage}
        @rtype: bsf.Stage
        """
        # Check if a Stage with this the name already exists and if so, return it.
        for stage in self.stage_list:
            if stage.name == name:
                return stage

        # Initialise a new Stage and add it to the Python list of Stage objects.

        stage = Stage(name=name, working_directory=self.genome_directory)
        self.stage_list.append(stage)

        # A "bsf.Stage" section specifies defaults for all Stage objects of an Analysis.

        section = Configuration.section_from_instance(instance=stage)
        stage.set_configuration(configuration=self.configuration, section=section)

        if self.debug > 1:
            print('Stage configuration section:', repr(section))

        # A "bsf.Analysis.Stage" or "bsf.analyses.*.Stage" pseudo-class section specifies
        # Analysis-specific or sub-class-specific options for the Stage, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=self), 'Stage'))
        stage.set_configuration(configuration=self.configuration, section=section)

        if self.debug > 1:
            print('Stage configuration section:', repr(section))

        # A "bsf.Analysis.Stage.name" or "bsf.analyses.*.Stage.name" section specifies defaults
        # for a particular Stage of an Analysis or sub-class, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=self), 'Stage', stage.name))
        stage.set_configuration(configuration=self.configuration, section=section)

        if self.debug > 1:
            print('Stage configuration section:', repr(section))

        return stage

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.Analysis} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @raise Exception: The specified section does not exist
        @return:
        @rtype:
        """
        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

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
        @type command: bsf.process.Command
        @return:
        @rtype:
        """
        # TODO: Phase out this method.
        # Once all accessory scripts are converted to Runnable and RunnableStep objects this method becomes redundant.
        assert isinstance(command, Command)

        section = Configuration.section_from_instance(instance=command)

        # For plain Executable objects append the value of the
        # Executable.program to make the configuration section more meaningful.

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
        @type runnable_step: bsf.process.RunnableStep
        @return:
        @rtype:
        """

        def _set_configuration(command, section):
            """Recursively set default C{bsf.argument.Argument} objects for a C{bsf.process.RunnableStep}.

            The method sets the default C{bsf.argument.Argument} objects for a C{bsf.process.RunnableStep},
            as well as for its contained sub C{bsf.process.Command} objects.
            @param command: C{bsf.process.Command}
            @type command: bsf.process.Command
            @param section: Configuration section
            @type section: str
            @return:
            @rtype:
            """
            assert isinstance(command, Command)
            assert isinstance(prefix, str)

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

        assert isinstance(runnable_step, RunnableStep)

        # Initially, the configuration section prefix is based on the Analysis class name and the Analysis Stage.name.

        prefix = Configuration.section_from_instance(instance=self)

        _set_configuration(command=runnable_step, section=prefix)

        return

    def set_stage_runnable(self, stage, runnable):
        """Create a C{bsf.process.Executable} to assign a C{bsf.Runnable} to a C{bsf.Stage}.

        In case the file in C{bsf.Runnable.get_relative_status_path} exists already,
        C{bsf.process.Executable.submit} will be set to C{False}.

        @param stage: C{bsf.Stage}
        @type stage: bsf.Stage
        @param runnable: C{bsf.Runnable}
        @type runnable: bsf.Runnable
        @return: C{bsf.process.Executable}
        @rtype: bsf.process.Executable
        @raise Exception: A C{bsf.Runnable.name} does not exist in C{bsf.Analysis.runnable_dict}
        @raise Exception: A C{bsf.Stage} does not exist in C{bsf.Analysis.stage_list}
        """
        assert isinstance(stage, Stage)
        assert isinstance(runnable, Runnable)

        if stage not in self.stage_list:
            raise Exception('A Stage with name ' + repr(stage.name) +
                            ' does not exist in the Analysis with name ' + repr(self.project_name) + '.')

        if runnable.name not in self.runnable_dict:
            raise Exception('A Runnable with name ' + repr(runnable.name) +
                            ' does not exist in the Analysis with name ' + repr(self.project_name) + '.')

        executable = Executable(name=runnable.name, program=Runnable.runner_script)
        executable.add_option_long(key='pickler-path', value=runnable.pickler_path)

        # Only submit the Executable if the status file does not exist already.
        if os.path.exists(runnable.absolute_runnable_status_file_path(success=True)):
            executable.submit = False

        stage.add_executable(executable=executable)

        return executable

    def run(self):
        """Run a C{bsf.Analysis}.

        @raise Exception: An C{bsf.Analysis.project_name} has not been defined
        @return:
        @rtype:
        """

        if not self.project_name:
            raise Exception('An Analysis project_name has not been defined.')

        # Some analyses such as FastQC do not require a genome_version,
        # nor a genome_version-specific output directory.
        # Also, add the e-mail address for UCSC track hubs into the genome subclass.

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard directory paths.

        self.cache_directory = self.configuration.get_absolute_path(
            file_path=self.cache_directory,
            default_path=Default.absolute_cache())

        self.input_directory = self.configuration.get_absolute_path(
            file_path=self.input_directory,
            default_path=Default.absolute_samples())

        self.output_directory = self.configuration.get_absolute_path(
            file_path=self.output_directory,
            default_path=Default.absolute_projects())

        # As a safety measure, to prevent creation of rogue directory paths, the output_directory has to exist.

        if not os.path.isdir(self.output_directory):
            raise Exception('The Analysis output_directory ' + repr(self.output_directory) + ' does not exist.')

        # Define project_directory and genome_directory instance variables.
        # If a genome_version option is present, append
        # it to the project_directory instance variable.
        # This allows analyses run against more than one directory and
        # simplifies UCSC Genome Browser track hub creation.

        self.project_directory = os.path.join(self.output_directory, self.project_name)

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
                raise Exception("An 'Analysis' requires an 'e_mail' configuration option.")

        if self.sas_file:
            # Populate a Collection from a SampleAnnotationSheet.
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
                print(self.collection.trace(1))
        else:
            # Create an empty Collection.
            self.collection = Collection()

        return

    def report(self):
        """Create a C{bsf.Analysis} report.

        The method must be implemented in a sub-class.

        @return:
        @rtype:
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
            defaults to a concatenation of C{bsf.Analysis.project_name} and C{bsf.Analysis.report_name}
        @type title: str
        @return: XHTML 1.0 header section as Python C{list} of Python C{str} or C{unicode} objects
        @rtype: list[str | unicode]
        """

        default = Default.get_global_default()

        if creator is None or not creator:
            creator = getpass.getuser()
            # The getpass.getuser method just relies on environment variables,
            # but at least works under Unix and Windows.

        if source is None or not source:
            source = inspect.getfile(inspect.currentframe())

        if title is None or not title:
            title = ' '.join((self.project_name, self.name))

        str_list = list()
        """ @type str_list: list[str | unicode] """

        if strict:
            str_list += '<!DOCTYPE html PUBLIC '
            str_list += '"-//W3C//DTD XHTML 1.0 Strict//EN" '
            str_list += '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n'
        else:
            str_list += '<!DOCTYPE html PUBLIC '
            str_list += '"-//W3C//DTD XHTML 1.0 Transitional//EN" '
            str_list += '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n'

        str_list += '\n'
        str_list += '<html xmlns="http://www.w3.org/1999/xhtml">\n'
        str_list += '<head>\n'
        str_list += '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />\n'
        str_list += '<link rel="stylesheet" href="/' + urllib.quote(s=default.url_relative_projects) + \
                    '/bsfpython.css" type="text/css" />\n'
        str_list += '<link rel="schema.DC" href="http://purl.org/DC/elements/1.0/" />\n'
        str_list += '<meta name="DC.Creator" content="' + cgi.escape(s=creator, quote=True) + '" />\n'
        str_list += '<meta name="DC.Date" content="' + datetime.datetime.now().isoformat() + '" />\n'
        str_list += '<meta name="DC.Source" content="' + cgi.escape(s=source, quote=True) + '" />\n'
        str_list += '<meta name="DC.Title" content="' + cgi.escape(s=title, quote=True) + '" />\n'
        str_list += '<title>' + cgi.escape(s=title, quote=True) + '</title>\n'
        str_list += '</head>\n'
        str_list += '\n'
        str_list += '<body>\n'
        str_list += '\n'

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
            defaults to C{bsf.standards.Default.operator_contact}
        @type contact: str
        @param institution: Institution name to be inserted into 'This report was generated by ...',
            defaults to C{bsf.standards.Default.operator_institution}
        @type institution: str
        @param url_protocol: The protocol section of the institution URL (i.e. http, https, ...),
            defaults to C{bsf.standards.Default.url_protocol}
        @type url_protocol: str
        @param url_host_name: The host name section of the institution URL (e.g. biomedical-sequencing.at),
            defaults to C{bsf.standards.Default.url_host_name}
        @type url_host_name: str
        @param title: Title element value,
            defaults to a concatenation of C{bsf.Analysis.project_name} and C{bsf.Analysis.report_name}
        @type title: str
        @return: XHTML 1.0 footer section as Python C{list} of Python C{str} or C{unicode} objects
        @rtype: list[str | unicode]
        """

        default = Default.get_global_default()

        if not contact:
            contact = Operator.get_contact()

        if not institution:
            institution = Operator.get_institution()

        if not url_protocol:
            url_protocol = default.url_protocol

        if not url_host_name:
            url_host_name = default.url_host_name

        if not title:
            title = ' '.join((self.project_name, self.name))

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list += '<hr class="footer" />\n'
        str_list += '<p class="footer">\n'
        if institution:
            str_list += 'This report was generated by ' + cgi.escape(s=institution, quote=True) + '\n'

        # Method bsf.standards.Default.url_absolute_base() would also return the default URL.
        if url_protocol:
            str_list += '<a href="' + url_protocol + '://' + \
                        url_host_name + '/">' + cgi.escape(s=url_host_name, quote=True) + '</a>.\n'
        else:
            str_list += '<a href="//' + \
                        url_host_name + '/">' + cgi.escape(s=url_host_name, quote=True) + '</a>.\n'

        str_list += '<br class="footer" />\n'
        if contact:
            # The e-mail address outside of an URL still needs HTML quoting.
            # After URL quoting, nothing critical needing HTML escaping should be left.
            str_list += 'Contact: <a href="mailto:' + urllib.quote(s=contact) + \
                        '?subject=' + urllib.quote(s=title) + '">'
            str_list += cgi.escape(s=contact, quote=True)
            str_list += '</a>'
        str_list += '</p>\n'
        str_list += '</body>\n'
        str_list += '</html>\n'
        str_list += '\n'

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

        The method automatically concatenates the XHTML header C{bsf.Analysis.get_html_header}, the XHTML content and
        the XHTML footer C{bsf.Analysis.get_html_footer} before returning the report.

        @param content: XHTML 1.0 content
        @type content: list[str | unicode]
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
            defaults to a concatenation of C{bsf.Analysis.project_name} and C{bsf.Analysis.report_name}
        @type title: str
        @param contact: Institution contact e-mail address,
            defaults to C{bsf.standards.Default.operator_contact}
        @type contact: str
        @param institution: Institution name to be inserted into 'This report was generated by ...',
            defaults to C{bsf.standards.Default.operator_institution}
        @type institution: str
        @param url_protocol: The protocol section of the institution URL (i.e. http, https, ...),
            defaults to C{bsf.standards.Default.url_protocol}
        @type url_protocol: str
        @param url_host_name: The host name section of the institution URL (e.g. biomedical-sequencing.at),
            defaults to C{bsf.standards.Default.url_host_name}
        @type url_host_name: str
        @return: XHTML 1.0 report as Python C{list} of Python C{str} or C{unicode} objects
        @rtype: list[str | unicode]
        """
        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list += self.get_html_header(
            strict=strict,
            creator=creator,
            source=source,
            title=title)
        str_list += content
        str_list += self.get_html_footer(
            contact=contact,
            institution=institution,
            url_protocol=url_protocol,
            url_host_name=url_host_name,
            title=title)

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
        """Write a XHTML 1.0 report I{prefix_report.html} file into the C{bsf.Analysis.genome_directory}.

        The method automatically concatenates the XHTML header C{bsf.Analysis.get_html_header}, the XHTML content and
        the XHTML footer C{bsf.Analysis.get_html_footer} before writing the file.

        @param content: XHTML 1.0 content
        @type content: list[str | unicode]
        @param prefix: A file name prefix (e.g. chipseq, rnaseq, ...), defaults to C{bsf.Analysis.prefix}
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
            defaults to a concatenation of C{bsf.Analysis.project_name} and C{bsf.Analysis.report_name}
        @type title: str
        @param contact: Institution contact e-mail address,
            defaults to C{bsf.standards.Default.operator_contact}
        @type contact: str
        @param institution: Institution name to be inserted into 'This report was generated by ...',
            defaults to C{bsf.standards.Default.operator_institution}
        @type institution: str
        @param url_protocol: The protocol section of the institution URL (i.e. http, https, ...),
            defaults to C{bsf.standards.Default.url_protocol}
        @type url_protocol: str
        @param url_host_name: The host name section of the institution URL (e.g. biomedical-sequencing.at),
            defaults to C{bsf.standards.Default.url_host_name}
        @type url_host_name: str
        @return:
        @rtype:
        """

        if prefix is None or not prefix:
            prefix = self.prefix

        file_handle = open(os.path.join(self.genome_directory, '_'.join((prefix, 'report.html'))), 'w')
        file_handle.writelines(
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
        file_handle.close()

        return

    def create_project_genome_directory(self):
        """Check for and create a C{bsf.Analysis.project_directory} or C{bsf.Analysis.genome_directory} if necessary.

        @return:
        @rtype:
        @raise Exception: Output (genome) directory does not exist
        """

        if not os.path.isdir(self.genome_directory):
            answer = raw_input(
                'Output (genome) directory ' + repr(self.genome_directory) + ' does not exist.\n'
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
        If not specified, the sub directory defaults to the Default.url_relative_projects instance variable.

        @param sub_directory: C{bsf.Analysis}-specific directory
        @type sub_directory: str
        @return: Symbolic link to the project directory
        @rtype: str
        @raise Exception: Public HTML path does not exist
        """
        default = Default.get_global_default()

        if sub_directory is None:
            sub_directory = default.url_relative_projects

        html_path = os.path.join(default.absolute_public_html(), sub_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the html_path directory has to exist.

        if not os.path.isdir(html_path):
            raise Exception(
                'The public HTML directory path ' + repr(html_path) + ' does not exist.\n'
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
                if not os.path.exists(path=source_path):
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
            answer = raw_input(
                'Public HTML link ' + repr(self.project_directory) + ' to ' + repr(final_path) + ' exists.\n'
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
            answer = raw_input(
                'Public HTML link ' + repr(self.project_directory) + ' to ' + repr(final_path) + ' does not exist.\n'
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
        @type ucsc_protocol: str
        @param ucsc_host_name: UCSC Genome Browser URL host name,
            defaults to C{bsf.standards.Default.ucsc_host_name}
        @type ucsc_host_name: str
        @return: A URL to attach a track to the UCSC Genome Browser
        @rtype: str
        """

        default = Default.get_global_default()

        if options_dict is None:
            options_dict = dict()

        if 'db' not in options_dict:
            options_dict['db'] = default.genome_alias_ucsc(genome_version=self.genome_version)

        # UCSC "browser" configuration dictionary.

        if browser_dict:
            pass

        # UCSC "track" configuration dictionary.

        if track_dict:
            key_list = track_dict.keys()
            key_list.sort(cmp=lambda x, y: cmp(x, y))

            options_dict['hgct_customText'] = 'track'
            for key in key_list:
                options_dict['hgct_customText'] += ' ' + key + '=' + track_dict[key]

        # UCSC protocol

        if not ucsc_protocol:
            ucsc_protocol = default.ucsc_protocol

        # UCSC host name

        if not ucsc_host_name:
            ucsc_host_name = default.ucsc_host_name

        # Strip leading colons to support protocol-independent URLs.
        return cgi.escape(
            s='{}://{}/cgi-bin/hgTracks?{}'.format(
                urllib.quote(s=ucsc_protocol),
                urllib.quote(s=ucsc_host_name),
                urllib.urlencode(query=options_dict)).lstrip(':'),
            quote=True)

    def ucsc_hub_url(self, link_path, options_dict=None):
        """Return a URL to automatically attach a UCSC Genome Browser Track Hub.

        @param link_path: Symbolic link path in the public HTML directory including project name and a UUID
        @type link_path: str | unicode
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
                Default.url_absolute_projects(),
                link_name,
                '_'.join((self.prefix, self.ucsc_name_hub))))

        return self.ucsc_track_url(options_dict=options_dict)

    def ucsc_hub_html_anchor(self, link_path):
        """Return a XHTML 1.0 anchor element to automatically attach a UCSC Genome Browser Track Hub.

        @param link_path: Symbolic link path in the public HTML directory including project name and a UUID
        @type link_path: str | unicode
        @see: C{bsf.Analysis.create_public_project_link}
        @return: XHTML 1.0 anchor element as Python C{list} of Python C{str} or C{unicode} objects
        @rtype: list[str | unicode]
        """
        str_list = list()
        str_list += 'UCSC Genome Browser Track Hub '
        str_list += '<a href="' + self.ucsc_hub_url(link_path=link_path) + '">' + self.project_name + '</a>'

        return str_list

    def ucsc_hub_write_hub(self, prefix=None):
        """Write a UCSC Track Hub I{prefix_hub.txt} file into the C{bsf.Analysis.project_directory}.

        The C{bsf.Analysis.project_directory} is one level above the C{bsf.Analysis.genome_directory}.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        @return:
        @rtype:
        """

        if prefix is None or not prefix:
            prefix = self.prefix

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list += 'hub ' + '_'.join((self.project_name, prefix)) + '\n'
        str_list += 'shortLabel ' + '_'.join((self.project_name, prefix)) + '\n'
        str_list += 'longLabel Project ' + '_'.join((self.project_name, prefix)) + '\n'
        str_list += 'genomesFile ' + '_'.join((prefix, self.ucsc_name_genomes)) + '\n'
        str_list += 'email ' + self.e_mail + '\n'

        file_handle = open(os.path.join(self.project_directory, '_'.join((prefix, self.ucsc_name_hub))), 'w')
        file_handle.writelines(str_list)
        file_handle.close()

        return

    def ucsc_hub_write_genomes(self, prefix=None):
        """Write a UCSC Track Hub I{prefix_genomes.txt} file into the C{bsf.Analysis.project_directory}.

        The C{bsf.Analysis.project_directory} is one level above the C{bsf.Analysis.genome_directory}.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        @return:
        @rtype:
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        file_path = os.path.join(self.project_directory, '_'.join((prefix, self.ucsc_name_genomes)))

        # If the file exists, read it first to retain any other genome assembly entries.
        genome_version_dict = dict()
        """ @type genome_version_dict: dict[str, str] """
        if os.path.exists(file_path):
            genome_version = None
            file_handle = open(file_path, 'r')
            for line in file_handle:
                line = line.strip()
                if not line:
                    continue
                line_list = line.split()
                if len(line_list) != 2:
                    warnings.warn('Malformed line ' + repr(line) + ' in UCSC genomes file ' +
                                  repr(file_path) + '.\n' +
                                  'Expected exactly two components after line splitting.')
                if line_list[0] == 'genome':
                    if genome_version is not None:
                        warnings.warn('Malformed line ' + repr(line) + ' in UCSC genomes file ' +
                                      repr(file_path) + '.\n' +
                                      "Got more than one 'genomes' lines in succession.")
                    genome_version = line_list[1]
                if line_list[0] == 'trackDb':
                    if genome_version is None:
                        warnings.warn('Malformed line ' + repr(line) + ' in UCSC genomes file ' +
                                      repr(file_path) + '.\n' +
                                      "Got a 'trackDb' line without a preceding genomes line.")
                    else:
                        genome_version_dict[genome_version] = line_list[1]
                        genome_version = None
            file_handle.close()

        # Resolve an eventual alias for the UCSC genome assembly name in "genome_version/prefix_tracks.txt".
        genome_version_dict[Default.genome_alias_ucsc(genome_version=self.genome_version)] = \
            '/'.join((self.genome_version, '_'.join((prefix, self.ucsc_name_tracks))))

        str_list = list()
        """ @type str_list: list[str | unicode] """

        genome_version_list = genome_version_dict.keys()
        genome_version_list.sort(cmp=lambda x, y: cmp(x, y))
        for genome_version in genome_version_list:
            str_list += 'genome ' + genome_version + '\n'
            str_list += 'trackDb ' + genome_version_dict[genome_version] + '\n'
            str_list += '\n'

        file_handle = open(file_path, 'w')
        file_handle.writelines(str_list)
        file_handle.close()

        return

    def ucsc_hub_write_tracks(self, content, prefix=None):
        """Write a UCSC Track Hub I{prefix_tracks.txt} file into the C{bsf.Analysis.genome_directory}.

        @param content: Content
        @type content: list[str | unicode]
        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        @return:
        @rtype:
        """

        if prefix is None or not prefix:
            prefix = self.prefix

        file_handle = open(os.path.join(self.genome_directory, '_'.join((prefix, self.ucsc_name_tracks))), 'w')
        file_handle.writelines(content)
        file_handle.close()

        return

    def ucsc_hub_to_file(self, content, prefix=None):
        """Write UCSC Genome Browser Track Hub files to disk.

        The method writes a I{prefix_hub.txt} and a I{prefix_genomes.txt} file into the
        C{bsf.Analysis.project_directory}, above the C{bsf.Analysis.genome_directory}, as well as a
        I{prefix_tracks.txt} file into the C{bsf.Analysis.genome_directory}.

        @param content: Content of the track database file
        @type content: list[str | unicode]
        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...), defaults to C{bsf.Analysis.prefix}
        @type prefix: str
        @return:
        @rtype:
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        self.ucsc_hub_write_hub(prefix=prefix)
        self.ucsc_hub_write_genomes(prefix=prefix)
        self.ucsc_hub_write_tracks(content=content, prefix=prefix)

        return

    def check_state(self):
        """Check the state of each C{bsf.Stage}.

        @return:
        @rtype:
        """
        for stage in self.stage_list:
            stage.check_state(debug=self.debug)

        return

    def submit(self, name=None):
        """Submit each C{bsf.Stage}.

        Submits each C{bsf.process.Executable} of either all C{bsf.Stage} objects or a named one and pickles
        each C{bsf.Runnable}.

        @param name: Only submit C{bsf.process.Executable} objects linked to C{bsf.Stage.name}
        @type name: bsf.Stage.name
        @return:
        @rtype:
        """

        # Pickle all Runnable objects.

        for runnable in self.runnable_dict.values():
            runnable.to_pickler_path()

        # Submit all Executable objects of all Stage objects.

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
                print(stage.trace(1))

        if name:
            if name == 'report':
                self.report()
            elif not submit:
                name_list = [stage.name for stage in self.stage_list]
                name_list.append('report')
                print('Valid Analysis Stage names are:', repr(name_list))

        return


class Stage(object):
    """The C{bsf.Stage} class represents a stage of a C{bsf.Analysis}.

    A C{bsf.Stage} represents C{bsf.process.Executable} or C{bsf.Runnable} objects that share
    similar resource requirements of a I{Distributed Resource Management System} (I{DRMS}).

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar working_directory: Working directory path
    @type working_directory: str
    @ivar implementation: Implementation (e.g. I{sge}, I{slurm}, ...)
    @type implementation: str
    @ivar memory_free_mem: Memory limit (free physical)
    @type memory_free_mem: str
    @ivar memory_free_swap: Memory limit (free swap)
    @type memory_free_swap: str
    @ivar memory_free_virtual: Memory limit (free virtual)
    @type memory_free_virtual: str
    @ivar memory_limit_hard: Memory limit (hard)
    @type memory_limit_hard: str
    @ivar memory_limit_soft: Memory limit (soft)
    @type memory_limit_soft: str
    @ivar node_list_exclude: List of nodes to exclude
    @type node_list_exclude: list[str]
    @ivar node_list_include: List of nodes to include
    @type node_list_include: list[str]
    @ivar time_limit: Time limit
    @type time_limit: str
    @ivar parallel_environment: Parallel environment
    @type parallel_environment: str
    @ivar queue: Queue
    @type queue: str
    @ivar threads: Number of threads
    @type threads: int
    @ivar hold: Hold on job scheduling
    @type hold: str
    @ivar is_script: C{bsf.process.Executable} objects represent shell scripts,
        or alternatively binary programs
    @type is_script: bool
    @ivar executable_list: Python C{list} of C{bsf.process.Executable} objects
    @type executable_list: list[bsf.process.Executable]
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
            executable_list=None):
        """Initialise a C{bsf.Stage}.

        @param name: Name
        @type name: str
        @param working_directory: Working directory
        @type working_directory: str
        @param implementation: Implementation (e.g. I{sge}, I{slurm}, ...)
        @type implementation: str
        @param memory_free_mem: Memory limit (free physical)
        @type memory_free_mem: str
        @param memory_free_swap: Memory limit (free swap)
        @type memory_free_swap: str
        @param memory_free_virtual: Memory limit (free virtual)
        @type memory_free_virtual: str
        @param memory_limit_hard: Memory limit (hard)
        @type memory_limit_hard: str
        @param memory_limit_soft: Memory limit (soft)
        @type memory_limit_soft: str
        @param node_list_exclude: List of nodes to exclude
        @type node_list_exclude: list[str]
        @param node_list_include: List of nodes to include
        @type node_list_include: list[str]
        @param time_limit: Time limit
        @type time_limit: str
        @param parallel_environment: Parallel environment
        @type parallel_environment: str
        @param queue: Queue
        @type queue: str
        @param threads: Number of threads
        @type threads: int
        @param hold: Hold on job scheduling
        @type hold: str
        @param is_script: C{bsf.process.Executable} objects represent shell scripts,
            or alternatively binary programs
        @type is_script: bool
        @param executable_list: Python C{list} of C{bsf.process.Executable} objects
        @type executable_list: list[bsf.process.Executable]
        @return:
        @rtype:
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

        if memory_free_mem is None:
            self.memory_free_mem = str()
        else:
            self.memory_free_mem = memory_free_mem

        if memory_free_swap is None:
            self.memory_free_swap = str()
        else:
            self.memory_free_swap = memory_free_swap

        if memory_free_virtual is None:
            self.memory_free_virtual = str()
        else:
            self.memory_free_virtual = memory_free_virtual

        if memory_limit_hard is None:
            self.memory_limit_hard = str()
        else:
            self.memory_limit_hard = memory_limit_hard

        if memory_limit_soft is None:
            self.memory_limit_soft = str()
        else:
            self.memory_limit_soft = memory_limit_soft

        if node_list_exclude is None:
            self.node_list_exclude = list()
        else:
            self.node_list_exclude = node_list_exclude

        if node_list_include is None:
            self.node_list_include = list()
        else:
            self.node_list_include = node_list_include

        if time_limit is None:
            self.time_limit = str()
        else:
            self.time_limit = time_limit

        if parallel_environment is None:
            self.parallel_environment = str()
        else:
            self.parallel_environment = parallel_environment

        if queue is None:
            self.queue = str()
        else:
            self.queue = queue

        if threads is None:
            self.threads = int(x=1)
        else:
            assert isinstance(threads, int)
            self.threads = threads

        if hold is None:
            self.hold = str()
        else:
            self.hold = hold

        if is_script is None:
            self.is_script = False
        else:
            assert isinstance(is_script, bool)
            self.is_script = is_script

        if executable_list is None:
            self.executable_list = list()
        else:
            self.executable_list = executable_list

        return

    def trace(self, level):
        """Trace a C{bsf.Stage}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  name:                 {!r}\n'.format(indent, self.name)
        output += '{}  working_directory:    {!r}\n'.format(indent, self.working_directory)
        output += '{}  implementation:       {!r}\n'.format(indent, self.implementation)
        output += '{}  memory_free_mem:      {!r}\n'.format(indent, self.memory_free_mem)
        output += '{}  memory_free_swap:     {!r}\n'.format(indent, self.memory_free_swap)
        output += '{}  memory_free_virtual:  {!r}\n'.format(indent, self.memory_free_virtual)
        output += '{}  memory_limit_hard:    {!r}\n'.format(indent, self.memory_limit_hard)
        output += '{}  memory_limit_soft:    {!r}\n'.format(indent, self.memory_limit_soft)
        output += '{}  node_list_exclude:    {!r}\n'.format(indent, self.node_list_exclude)
        output += '{}  node_list_include:    {!r}\n'.format(indent, self.node_list_include)
        output += '{}  time_limit:           {!r}\n'.format(indent, self.time_limit)
        output += '{}  queue:                {!r}\n'.format(indent, self.queue)
        output += '{}  parallel_environment: {!r}\n'.format(indent, self.parallel_environment)
        output += '{}  threads:              {!r}\n'.format(indent, self.threads)
        output += '{}  hold:                 {!r}\n'.format(indent, self.hold)
        output += '{}  is_script:            {!r}\n'.format(indent, self.is_script)

        output += '{}  executable_list:\n'.format(indent)

        for executable in self.executable_list:
            output += executable.trace(level=level + 2)

        return output

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.Stage} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

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
            self.node_list_exclude = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    configuration.config_parser.get(section=section, option=option).split(',')))

        option = 'node_list_include'
        if configuration.config_parser.has_option(section=section, option=option):
            self.node_list_include = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    configuration.config_parser.get(section=section, option=option).split(',')))

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
            self.threads = configuration.config_parser.get(section=section, option=option)

        return

    def add_executable(self, executable):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.process.Executable}.

        @param executable: C{bsf.process.Executable}
        @type executable: bsf.process.Executable
        @return: C{bsf.process.Executable}
        @rtype: bsf.process.Executable
        """

        assert isinstance(executable, Executable)

        self.executable_list.append(executable)

        return executable

    def check_state(self, debug=0):
        """Check the state of each C{bsf.process.Executable}.

        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        # Dynamically import the module specific for the configured DRMS implementation.

        python_module = importlib.import_module(name='.'.join((__name__, 'drms', self.implementation)))
        python_module.check_state(stage=self, debug=debug)

        return

    def submit(self, debug=0):
        """Submit a command line for each C{bsf.process.Executable}.

        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        # Dynamically import the module specific for the configured DRMS implementation.

        python_module = importlib.import_module(name='.'.join((__name__, 'drms', self.implementation)))
        python_module.submit(stage=self, debug=debug)

        return


class FilePath(object):
    """The C{bsf.FilePath} class represents formalised file path information for the C{bsf.Runnable} class.

    Each C{bsf.Runnable} class is expected to define its corresponding C{bsf.FilePath} sub-class.
    Attributes:
    @ivar prefix: File path prefix
    @type prefix: str | unicode
    #ivar temporary_directory: Temporary directory path
    #type temporary_directory: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.FilePath}.

        @param prefix: File path prefix
        @type prefix: str | unicode
        @return:
        @rtype:
        """

        self.prefix = prefix
        # self.temporary_directory = prefix + '_temporary'

        return


class Runnable(object):
    """The C{bsf.Runnable} class represents one or more C{bsf.process.Executable} objects for the I{Runner} script.

    A C{bsf.Runnable} holds all information to run one or more C{bsf.process.Executable} objects through the
    C{bsf.Runnable.runner_script}. It can be thought of a GNU Bash script that executes as set of
    C{bsf.process.RunnableStep} objects reflecting commands of a GNU Bash script.

    Attributes:
    @cvar runner_script: Name of the I{Runner} script
    @type runner_script: str | unicode
    @ivar name: Name
    @type name: str
    @ivar code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
        C{bsf.process.Executable} objects via the C{bsf.Runnable.runner_script}.
    @type code_module: str
    @ivar cache_directory: Cache directory
    @type cache_directory: str | unicode
    @ivar cache_path_dict: Python C{dict} of Python C{str} (name) key and
        Python C{str} (file_path) value data of files that will be copied into the C{bsf.Runnable.cache_directory}
    @type cache_path_dict: dict[str, str | unicode]
    @ivar file_path_object: C{bsf.FilePath}
    @type file_path_object: bsf.FilePath
    @ivar runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
    @type runnable_step_list: list[bsf.process.RunnableStep]
    @ivar working_directory: Working directory to write C{pickle.Pickler} files
    @type working_directory: str | unicode
    @ivar debug: Debug level
    @type debug: int
    """

    runner_script = 'bsf_runner.py'

    def __init__(
            self,
            name,
            code_module,
            working_directory,
            cache_directory=None,
            cache_path_dict=None,
            file_path_object=None,
            runnable_step_list=None,
            debug=0):
        """Initialise a C{bsf.Runnable}.

        @param name: Name
        @type name: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{bsf.process.Executable} objects via the C{bsf.Runnable.runner_script}
        @type code_module: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str | unicode
        @param cache_directory: Cache directory
        @type cache_directory: str | unicode
        @param cache_path_dict: Python C{dict} of Python C{str} (name) key and
            Python C{str} (file_path) value data of files that will be copied into the C{bsf.Runnable.cache_directory}
        @type cache_path_dict: dict[str, str | unicode]
        @param file_path_object: C{bsf.FilePath}
        @type file_path_object: bsf.FilePath
        @param runnable_step_list: Python C{list} of C{bsf.process.RunnableStep} objects
        @type runnable_step_list: list[bsf.process.RunnableStep]
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        super(Runnable, self).__init__()

        self.name = name  # Can be None.
        self.code_module = code_module  # Can be None.
        self.working_directory = working_directory  # Can be None.

        if cache_directory is None:
            self.cache_directory = str()
        else:
            self.cache_directory = cache_directory

        if cache_path_dict is None:
            self.cache_path_dict = dict()
        else:
            self.cache_path_dict = cache_path_dict

        if file_path_object is None:
            self.file_path_object = FilePath(prefix='default_file_path')
        else:
            self.file_path_object = file_path_object

        if runnable_step_list is None:
            self.runnable_step_list = list()
        else:
            self.runnable_step_list = runnable_step_list

        if debug is None:
            self.debug = int(x=0)
        else:
            assert isinstance(debug, int)
            self.debug = debug

        return

    def trace(self, level=1):
        """Trace a C{bsf.Runnable}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  name: {!r}\n'.format(indent, self.name)
        output += '{}  code_module: {!r}\n'.format(indent, self.code_module)
        output += '{}  working_directory: {!r}\n'.format(indent, self.working_directory)
        output += '{}  cache_directory: {!r}\n'.format(indent, self.cache_directory)
        output += '{}  cache_path_dict: {!r}\n'.format(indent, self.cache_path_dict)
        output += '{}  file_path_object: {!r}\n'.format(indent, self.file_path_object)
        output += '{}  runnable_step_list: {!r}\n'.format(indent, self.runnable_step_list)
        output += '{}  debug: {!r}\n'.format(indent, self.debug)

        output += '{}  Python dict of Python str (cache path) objects:\n'.format(indent)
        key_list = self.cache_path_dict.keys()
        key_list.sort(cmp=lambda x, y: cmp(x, y))
        for key in key_list:
            output += '{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.cache_path_dict[key])

        output += '{}  Python list of RunnableStep objects:\n'.format(indent)
        for runnable_step in self.runnable_step_list:
            output += runnable_step.trace(level=level + 1)

        return output

    def add_runnable_step(self, runnable_step=None):
        """Convenience method to facilitate initialising, adding and returning a C{bsf.process.RunnableStep}.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep
        @return: C{bsf.process.RunnableStep}
        @rtype: bsf.process.RunnableStep
        """

        if runnable_step is None:
            return

        assert isinstance(runnable_step, RunnableStep)

        self.runnable_step_list.append(runnable_step)

        return runnable_step

    @property
    def pickler_path(self):
        """Get the Python C{pickle.Pickler} file path.

        @return: Python C{pickle.Pickler} file path
        @rtype: str | unicode
        """

        return os.path.join(self.working_directory, '.'.join((self.name, 'pkl')))

    def to_pickler_path(self):
        """Write this C{bsf.Runnable} as a Python C{pickle.Pickler} file into the working directory.

        @return:
        @rtype:
        """

        pickler_file = open(self.pickler_path, 'wb')
        pickler = pickle.Pickler(file=pickler_file, protocol=pickle.HIGHEST_PROTOCOL)
        pickler.dump(obj=self)
        pickler_file.close()

        return

    @classmethod
    def from_pickler_path(cls, file_path):
        """Create a C{bsf.Runnable} from a Python C{pickle.Pickler} file via Python C{pickle.Unpickler}.

        @param file_path: File path to a Python C{pickle.Pickler} file
        @type file_path: str | unicode
        @return: C{bsf.Runnable}
        @rtype: bsf.Runnable
        """

        pickler_file = open(file_path, 'rb')
        unpickler = pickle.Unpickler(file=pickler_file)
        runnable = unpickler.load()
        """ @type runnable: bsf.Runnable """
        pickler_file.close()

        # Did the Unpickler really return a Runnable object?
        assert isinstance(runnable, Runnable)

        return runnable

    @property
    def get_relative_cache_directory_path(self):
        """Get the relative cache directory path of a C{bsf.Runnable}.

        @return: Relative cache directory path (i.e. C{bsf.Runnable.name}_cache)
        @rtype: str
        """

        return '_'.join((self.name, 'cache'))

    @property
    def get_relative_temporary_directory_path(self):
        """Get the relative temporary directory path of a C{bsf.Runnable}.

        @return: Relative temporary directory path (i.e. C{bsf.Runnable.name}_temporary)
        @rtype: str
        """

        return '_'.join((self.name, 'temporary'))

    @property
    def get_absolute_cache_directory_path(self):
        """Get the absolute cache directory path including the C{bsf.Runnable.cache_directory}.

        If C{bsf.Runnable.cache_directory} is not defined, C{bsf.Runnable.working_directory} will be prepended.
        Since the relative cache directory path includes the C{bsf.Runnable.name},
        the directory is C{bsf.Runnable}-specific.
        (i.e. C{bsf.Runnable.cache_directory}/C{bsf.Runnable.name}_cache or
        C{bsf.Runnable.working_directory}/C{bsf.Runnable.name}_cache)

        @return: Absolute cache directory path
        @rtype: str
        """

        if self.cache_directory:
            return Configuration.get_absolute_path(
                file_path=self.get_relative_cache_directory_path,
                default_path=self.cache_directory)
        else:
            return Configuration.get_absolute_path(
                file_path=self.get_relative_cache_directory_path,
                default_path=self.working_directory)

    def get_absolute_cache_file_path(self, file_path):
        """Get the absolute cache file path for a file path.

        @param file_path: Default file path
        @type file_path: str | unicode
        @return: Absolute cache file path
        @rtype: str
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        return os.path.join(self.get_absolute_cache_directory_path, file_name)

    @property
    def get_absolute_temporary_directory_path(self):
        """Get the absolute temporary directory path including the C{bsf.Runnable.working_directory}.

        @return: Absolute temporary directory path
            (i.e. C{bsf.Runnable.working_directory}/C{bsf.Runnable.name}_temporary)
        @rtype: str
        """

        return os.path.join(self.working_directory, self.get_relative_temporary_directory_path)

    def runnable_status_file_path(self, success=True):
        """Get the status file path for a C{bsf.Runnable}.

        @param success: Successful completion
        @type success: bool
        @return: Status file path
        @rtype: str
        """
        if success:
            return '_'.join((self.name, 'completed.txt'))
        else:
            return '_'.join((self.name, 'failed.txt'))

    def runnable_status_file_create(self, success=True):
        """Create an empty status file for a C{bsf.Runnable}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @param success: Successful completion
        @type success: bool
        @return:
        @rtype:
        """
        status_path = self.runnable_status_file_path(success=success)
        open(status_path, 'w').close()

        return

    def runnable_status_file_remove(self):
        """Remove the status file for a C{bsf.Runnable}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @return:
        @rtype:
        """
        # Automatically remove both status files, successful or not.

        status_path = self.runnable_status_file_path(success=True)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        status_path = self.runnable_status_file_path(success=False)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        return

    def runnable_step_status_file_path(self, runnable_step, success=True):
        """Get the status file path for a C{bsf.process.RunnableStep} of a C{bsf.Runnable}.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep
        @param success: Successful completion
        @type success: bool
        @return: Status file path
        @rtype: str
        """
        assert isinstance(runnable_step, RunnableStep)

        if success:
            return '_'.join((self.name, runnable_step.name, 'completed.txt'))
        else:
            return '_'.join((self.name, runnable_step.name, 'failed.txt'))

    def runnable_step_status_file_create(self, runnable_step, success=True):
        """Create an empty status file for a C{bsf.process.RunnableStep} of a C{bsf.Runnable}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep
        @param success: Successful completion
        @type success: bool
        @return:
        @rtype:
        """
        assert isinstance(runnable_step, RunnableStep)

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=success)
        open(status_path, 'w').close()

        return

    def runnable_step_status_file_remove(self, runnable_step):
        """Remove the status file for a C{bsf.process.RunnableStep} of a C{bsf.Runnable}.

        This method is mainly used by C{bsf.runnable.generic} and related modules.

        @param runnable_step: C{bsf.process.RunnableStep}
        @type runnable_step: bsf.process.RunnableStep
        @return:
        @rtype:
        """
        assert isinstance(runnable_step, RunnableStep)

        if runnable_step is None:
            return

        # Automatically remove both status files, successful or not.

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=True)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        status_path = self.runnable_step_status_file_path(runnable_step=runnable_step, success=False)
        try:
            os.remove(status_path)
        except OSError as exception:
            if exception.errno != errno.ENOENT:
                raise

        return

    def absolute_runnable_status_file_path(self, success=True):
        """Get the absolute status file path including the C{bsf.Runnable.working_directory}.

        @return: Absolute status file path
            C{bsf.Runnable.working_directory}/C{bsf.Runnable.name}_completed.txt or
            C{bsf.Runnable.working_directory}/C{bsf.Runnable.name}_failed.txt
        @rtype: str
        """
        return os.path.join(self.working_directory, self.runnable_status_file_path(success=success))
