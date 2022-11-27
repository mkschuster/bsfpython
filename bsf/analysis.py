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
"""The :py:mod:`bsf.analysis` module provides classes modelling analyses.

    - The :py:class:`bsf.analysis.Analysis` class represents a high-level analysis that involves one or more stages.
    - The :py:class:`bsf.analysis.Stage` class represents processes that require similar computational resources.
"""
import datetime
import errno
import getpass
import html
import importlib
import inspect
import logging
import os
import urllib.parse
import uuid
from typing import Dict, List, Optional

from bsf.ngs import Collection, Sample
from bsf.procedure import Runnable, ConcurrentRunnable, ConsecutiveRunnable
from bsf.process import Command, Executable, RunnableStep
from bsf.standards import Configuration, StandardFilePath, Operator, Genome, Transcriptome, UCSC, URL

module_logger = logging.getLogger(name=__name__)


class Stage(object):
    """The :py:class:`bsf.analysis.Stage` class represents a stage of a :py:class:`bsf.analysis.Analysis` object.

    A :py:class:`bsf.analysis.Stage` represents :py:class:`bsf.process.Executable` or
    :py:class:`bsf.procedure.Runnable` objects that share similar resource requirements of a
    :literal:`Distributed Resource Management System` (:literal:`DRMS`).

    :ivar name: A name
    :type name: str
    :ivar working_directory: A working directory path.
    :type working_directory: str
    :ivar implementation: An implementation (e.g., :literal:`sge`, :literal:`slurm`, ...).
    :type implementation: str
    :ivar memory_free_mem: A memory limit (free physical).
    :type memory_free_mem: str | None
    :ivar memory_free_swap: A memory limit (free swap).
    :type memory_free_swap: str | None
    :ivar memory_free_virtual: A memory limit (free virtual).
    :type memory_free_virtual: str | None
    :ivar memory_limit_hard: A memory limit (hard).
    :type memory_limit_hard: str | None
    :ivar memory_limit_soft: A memory limit (soft).
    :type memory_limit_soft: str | None
    :ivar node_list_exclude: A list of nodes to exclude.
    :type node_list_exclude: list[str] | None
    :ivar node_list_include: A list of nodes to include.
    :type node_list_include: list[str] | None
    :ivar time_limit: A time limit.
    :type time_limit: str | None
    :ivar parallel_environment: A parallel environment.
    :type parallel_environment: str | None
    :ivar queue: A queue name.
    :type queue: str | None
    :ivar reservation: A reservation name.
    :type reservation: str | None
    :ivar threads: A number of threads.
    :type threads: int
    :ivar hold: Request a hold on job scheduling.
    :type hold: str | None
    :ivar is_script: The :py:class:`bsf.process.Executable` objects represent shell scripts,
        or alternatively binary programs.
    :type is_script: bool
    :ivar executable_list: A Python :py:class:`list` object of :py:class:`bsf.process.Executable` objects.
    :type executable_list: list[Executable]
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
            reservation=None,
            threads=1,
            hold=None,
            is_script=False,
            template_script=None,
            executable_list=None):
        """Initialise a :py:class:`bsf.analysis.Stage` object.

        :param name: An name.
        :type name: str
        :param working_directory: A working directory.
        :type working_directory: str
        :param implementation: An implementation (e.g., :literal:`sge`, :literal:`slurm`, ...).
        :type implementation: str
        :param memory_free_mem: A memory limit (free physical).
        :type memory_free_mem: str | None
        :param memory_free_swap: A memory limit (free swap).
        :type memory_free_swap: str | None
        :param memory_free_virtual: A memory limit (free virtual).
        :type memory_free_virtual: str | None
        :param memory_limit_hard: A memory limit (hard).
        :type memory_limit_hard: str | None
        :param memory_limit_soft: A memory limit (soft).
        :type memory_limit_soft: str | None
        :param node_list_exclude: A list of nodes to exclude.
        :type node_list_exclude: list[str] | None
        :param node_list_include: A list of nodes to include.
        :type node_list_include: list[str] | None
        :param time_limit: A time limit.
        :type time_limit: str | None
        :param parallel_environment: A parallel environment.
        :type parallel_environment: str | None
        :param queue: A queue name.
        :type queue: str | None
        :param reservation: A reservation name.
        :type reservation: str | None
        :param threads: A number of threads.
        :type threads: int
        :param hold: Request a hold on job scheduling.
        :type hold: bool | None
        :param is_script: The :py:class:`bsf.process.Executable` objects represent shell scripts,
            or alternatively binary programs.
        :type is_script: bool
        :param template_script: Template script for submission.
        :type template_script: str
        :param executable_list: A Python :py:class:`list` object of :py:class:`bsf.process.Executable` objects.
        :type executable_list: list[Executable]
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
        self.reservation = reservation

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

    def __repr__(self):
        return \
            f'Stage(' \
            f'name={self.name!r}, ' \
            f'working_directory={self.working_directory!r}, ' \
            f'implementation={self.implementation!r}, ' \
            f'memory_free_mem={self.memory_free_mem!r}, ' \
            f'memory_free_swap={self.memory_free_swap!r}, ' \
            f'memory_free_virtual={self.memory_free_virtual!r}, ' \
            f'memory_limit_hard={self.memory_limit_hard!r}, ' \
            f'memory_limit_soft={self.memory_limit_soft!r}, ' \
            f'node_list_exclude={self.node_list_exclude!r}, ' \
            f'node_list_include={self.node_list_include!r}, ' \
            f'time_limit={self.time_limit!r}, ' \
            f'parallel_environment={self.parallel_environment!r}, ' \
            f'queue={self.queue!r}, ' \
            f'reservation{self.reservation!r}, ' \
            f'threads={self.threads!r}, ' \
            f'hold={self.hold!r}, ' \
            f'is_script={self.is_script!r}, ' \
            f'template_script={self.template_script!r}, ' \
            f'executable_list={self.executable_list!r})'

    def trace(self, level):
        """Trace a :py:class:`bsf.analysis.Stage` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

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
        str_list.append('{}  reservation:          {!r}\n'.format(indent, self.reservation))
        str_list.append('{}  threads:              {!r}\n'.format(indent, self.threads))
        str_list.append('{}  hold:                 {!r}\n'.format(indent, self.hold))
        str_list.append('{}  is_script:            {!r}\n'.format(indent, self.is_script))
        str_list.append('{}  template_script:      {!r}\n'.format(indent, self.template_script))

        str_list.append('{}  executable_list:\n'.format(indent))

        for executable in self.executable_list:
            str_list.extend(executable.trace(level=level + 2))

        return str_list

    def set_configuration(self, configuration, section):
        """Set instance variables of a :py:class:`bsf.analysis.Stage` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        """
        if not configuration.config_parser.has_section(section=section):
            raise Exception(
                f'A section {section!r} is not defined in Configuration files:\n'
                f'{configuration.file_path_list!r}')

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

        option = 'reservation'
        if configuration.config_parser.has_option(section=section, option=option):
            self.reservation = configuration.config_parser.get(section=section, option=option)

        option = 'threads'
        if configuration.config_parser.has_option(section=section, option=option):
            self.threads = configuration.config_parser.getint(section=section, option=option)

        option = 'template_script'
        if configuration.config_parser.has_option(section=section, option=option):
            self.template_script = configuration.config_parser.get(section=section, option=option)

        return

    def add_executable(self, executable):
        """Convenience method to facilitate initialising, adding and returning a
        :py:class:`bsf.process.Executable` object.

        :param executable: A :py:class:`bsf.process.Executable` object.
        :type executable: Executable
        :return: A :py:class:`bsf.process.Executable` object.
        :rtype: Executable
        """
        self.executable_list.append(executable)

        return executable

    def check_state(self):
        """Check the state of each :py:class:`bsf.process.Executable` object.
        """
        # Dynamically import the module specific for the configured DRMS implementation.

        module_type = importlib.import_module(name='.'.join(('bsf', 'drms', self.implementation)))

        check_state_function = getattr(module_type, 'check_state')
        check_state_function(stage=self)

        return

    def submit(self, drms_submit=None):
        """Submit a command line for each :py:class:`bsf.process.Executable` object.

        :param drms_submit: Submit to the :emphasis:`Distributed Resource Management System` (DRMS).
        :type drms_submit: bool | None
        """
        # Dynamically import the module specific for the configured DRMS implementation.

        module_type = importlib.import_module(name='.'.join(('bsf', 'drms', self.implementation)))

        submit_function = getattr(module_type, 'submit')
        submit_function(stage=self, drms_submit=drms_submit)

        return


class Analysis(object):
    """The :py:class:`bsf.analysis.Analysis` class represents a high-level analysis.

    It consists of one or more :py:class:`bsf.analysis.Stage` objects that may run one or more
    :py:class:`bsf.process.Executable` or :py:class:`bsf.process.RunnableStep` objects (programs).

    :cvar name: :py:attr:`bsf.analysis.Analysis.name` that should be overridden by subclasses.
    :type name: str
    :cvar prefix: :py:attr:`bsf.analysis.Analysis.prefix` that should be overridden by subclasses.
    :type prefix: str
    :cvar ucsc_name_hub: UCSC Genome Browser Track Hub :literal:`hub` file name.
    :type ucsc_name_hub: str
    :cvar ucsc_name_genomes: UCSC Genome Browser Track Hub :literal:`genomes` file name.
    :type ucsc_name_genomes: str
    :cvar ucsc_name_tracks: UCSC Genome Browser Track Hub :literal:`tracks` file name.
    :type ucsc_name_tracks: str
    :ivar configuration: A :py:class:`bsf.standards.Configuration` object.
    :type configuration: Configuration
    :ivar project_name: A project name.
    :type project_name: str | None
    :ivar genome_version: A genome assembly version.
    :type genome_version: str | None
    :ivar cache_directory: A cache directory path.
    :type cache_directory: str | None
    :ivar input_directory: An input directory path.
    :type input_directory: str | None
    :ivar output_directory: An output directory path.
    :type output_directory: str | None
    :ivar project_directory: A project directory path, normally under the output directory path.
    :type project_directory: str | None
    :ivar genome_directory: A genome directory path, normally under the project directory path.
    :type genome_directory: str | None
    :ivar report_style_path: Report :literal:`CSS` file path.
    :type report_style_path: str | None
    :ivar report_header_path: Report header :literal:`XHTML 1.0` file path.
    :type report_header_path: str | None
    :ivar report_footer_path: Report footer :literal:`XHTML 1.0` file path.
    :type report_footer_path: str | None
    :ivar sas_file: A Sample Annotation Sheet (SAS) file path.
    :type sas_file: str | None
    :ivar sas_prefix: A prefix to columns in a Sample Annotation Sheet
        (e.g., :literal:`[Control] Sample`, :literal:`[Treatment] Sample`, ...).
    :type sas_prefix: str | None
    :ivar e_mail: An e-mail address for a UCSC Genome Browser Track Hub.
    :type e_mail: str | None
    :ivar stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
    :type stage_list: list[Stage]
    :ivar runnable_dict: A Python :py:class:`dict` object of
        Python :py:class:`str` (:py:attr:`bsf.procedure.Runnable.name`) key objects and
        :py:class:`bsf.procedure.Runnable` value objects.
    :type runnable_dict: dict[Runnable.name, Runnable]
    :ivar collection: A :py:class:`bsf.ngs.Collection` object.
    :type collection: Collection
    :ivar sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
    :type sample_list: list[Sample]
    """

    name = 'Analysis'
    prefix = 'analysis'

    ucsc_name_hub = 'hub.txt'
    ucsc_name_genomes = 'genomes.txt'
    ucsc_name_tracks = 'tracks.txt'

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new :py:class:`bsf.analysis.Analysis` object from a UNIX-style configuration file path.

        The configuration file in :py:attr:`bsf.standards.Configuration.global_file_path` is read as default,
        before the project-specific one gets read, if it is not the same file.

        :param config_path: A UNIX-style configuration file path.
        :type config_path: str
        :return: A :py:class:`bsf.analysis.Analysis` object.
        :rtype: Analysis
        """
        return cls.from_configuration(
            configuration=Configuration.from_file_path_list(
                file_path_list=[Configuration.get_global_file_path(), config_path]))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new :py:class:`bsf.analysis.Analysis` object from a :py:class:`bsf.standards.Configuration` object.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :return: A :py:class:`bsf.analysis.Analysis` object.
        :rtype: Analysis
        """
        analysis = cls(configuration=configuration)

        # A "module.class" configuration section specifies defaults for this Analysis or sub-class
        # (i.e., "bsf.analysis.Analysis" or "bsf.analyses.*", respectively).

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
            report_style_path=None,
            report_header_path=None,
            report_footer_path=None,
            sas_file=None,
            sas_prefix=None,
            e_mail=None,
            stage_list=None,
            runnable_dict=None,
            collection=None,
            sample_list=None):
        """Initialise a :py:class:`bsf.analysis.Analysis` object.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration | None
        :param project_name: A project name.
        :type project_name: str | None
        :param genome_version: A genome assembly version.
        :type genome_version: str | None
        :param cache_directory: A cache directory path.
        :type cache_directory: str | None
        :param input_directory: An input directory path.
        :type input_directory: str | None
        :param output_directory: An output directory path.
        :type output_directory: str | None
        :param project_directory: A project directory path, normally under the output directory path.
        :type project_directory: str | None
        :param genome_directory: A genome directory path, normally under the project directory path.
        :type genome_directory: str | None
        :param report_style_path: Report :literal:`CSS` file path.
        :type report_style_path: str | None
        :param report_header_path: Report header :literal:`XHTML 1.0` file path.
        :type report_header_path: str | None
        :param report_footer_path: Report footer :literal:`XHTML 1.0` file path.
        :type report_footer_path: str | None
        :param sas_file: A Sample Annotation Sheet (SAS) file path.
        :type sas_file: str | None
        :param sas_prefix: A prefix to columns in a Sample Annotation Sheet
            (e.g., :literal:`[Control] Sample`, :literal:`[Treatment] Sample`, ...).
        :type sas_prefix: str | None
        :param e_mail: An e-mail address for a UCSC Genome Browser Track Hub.
        :type e_mail: str | None
        :param stage_list: A Python :py:class:`list` object of :py:class:`bsf.analysis.Stage` objects.
        :type stage_list: list[Stage] | None
        :param runnable_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (:py:attr:`bsf.procedure.Runnable.name`) key objects and
            :py:class:`bsf.procedure.Runnable` value objects.
        :type runnable_dict: dict[Runnable.name, Runnable] | None
        :param collection: A :py:class:`bsf.ngs.Collection` object.
        :type collection: Collection | None
        :param sample_list: A Python :py:class:`list` object of :py:class:`bsf.ngs.Sample` objects.
        :type sample_list: list[Sample] | None
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
        self.report_style_path = report_style_path
        self.report_header_path = report_header_path
        self.report_footer_path = report_footer_path
        self.sas_file = sas_file
        self.sas_prefix = sas_prefix
        self.e_mail = e_mail

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

        self._public_project_link_path: Optional[str] = None

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.analysis.Analysis` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  project_name: {!r}\n'.format(indent, self.project_name))
        str_list.append('{}  genome_version: {!r}\n'.format(indent, self.genome_version))
        str_list.append('{}  cache_directory: {!r}\n'.format(indent, self.cache_directory))
        str_list.append('{}  input_directory: {!r}\n'.format(indent, self.input_directory))
        str_list.append('{}  output_directory: {!r}\n'.format(indent, self.output_directory))
        str_list.append('{}  genome_directory: {!r}\n'.format(indent, self.genome_directory))
        str_list.append('{}  report_style_path: {!r}\n'.format(indent, self.report_style_path))
        str_list.append('{}  report_header_path: {!r}\n'.format(indent, self.report_header_path))
        str_list.append('{}  report_footer_path: {!r}\n'.format(indent, self.report_footer_path))
        str_list.append('{}  sas_file: {!r}\n'.format(indent, self.sas_file))
        str_list.append('{}  sas_prefix: {!r}\n'.format(indent, self.sas_prefix))
        str_list.append('{}  e_mail: {!r}\n'.format(indent, self.e_mail))
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
        """Convenience method to facilitate initialising, adding and returning a
        :py:class:`bsf.analysis.Stage` object.

        If the :py:class:`bsf.analysis.Stage` exists already in the
        :py:attr:`bsf.analysis.Analysis.stage_list` attribute,
        the method returns the already existing :py:class:`bsf.analysis.Stage` object.

        :param stage: A :py:class:`bsf.analysis.Stage` object.
        :type stage: Stage
        :return: A :py:class:`bsf.analysis.Stage` object.
        :rtype: Stage
        """
        if stage not in self.stage_list:
            self.stage_list.append(stage)

        return stage

    def add_runnable(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a
        :py:class:`bsf.procedure.Runnable` object.

        :param runnable: A :py:class:`bsf.procedure.Runnable` object.
        :type runnable: Runnable
        :return: A :py:class:`bsf.procedure.Runnable` object.
        :rtype: Runnable
        :raise Exception: The :py:attr:`bsf.procedure.Runnable.name` already exists in the
            :py:attr:`bsf.analysis.Analysis.runnable_dict`.
        """
        if runnable.name in self.runnable_dict:
            raise Exception(f'A Runnable with name {runnable.name!r} already exists in '
                            f'Analysis {self.project_name!r}.')
        else:
            self.runnable_dict[runnable.name] = runnable

        return runnable

    def add_runnable_concurrent(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a
        :py:class:`bsf.procedure.ConcurrentRunnable` object.

        :param runnable: A :py:class:`bsf.procedure.ConcurrentRunnable` object.
        :type runnable: ConcurrentRunnable
        :return: A :py:class:`bsf.procedure.ConcurrentRunnable` object.
        :rtype: ConcurrentRunnable
        :raise Exception: The :py:attr:`bsf.procedure.Runnable.name` already exists in the
            :py:attr:`bsf.analysis.Analysis.runnable_dict`.
        """
        return self.add_runnable(runnable=runnable)

    def add_runnable_consecutive(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a
        :py:class:`bsf.procedure.ConsecutiveRunnable` object.

        :param runnable: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :type runnable: ConsecutiveRunnable
        :return: A :py:class:`bsf.procedure.ConsecutiveRunnable` object.
        :rtype: ConsecutiveRunnable
        :raise Exception: The :py:attr:`bsf.procedure.Runnable.name` already exists in the
            :py:attr:`bsf.analysis.Analysis.runnable_dict`.
        """
        return self.add_runnable(runnable=runnable)

    def add_sample(self, sample):
        """Add a :py:class:`bsf.ngs.Sample` object to the Python :py:class:`list` of :py:class:`bsf.ngs.Sample` objects.

        If the :py:class:`bsf.ngs.Sample` object already exists in the :py:class:`bsf.analysis.Analysis` object,
        the method just returns.
        The check is based on the Python :literal:`in` comparison operator and in lack of a specific
        :py:meth:`__cmp__` method, relies on object identity (i.e., address).

        :param sample: A :py:class:`bsf.ngs.Sample` object.
        :type sample: Sample
        """
        if sample not in self.sample_list:
            self.sample_list.append(sample)

        return

    def get_annotation_file(self, prefix_list, suffix):
        """Get a project and genome-specific annotation file.

        Based on the project name, a list of file name prefixes and one file name suffix,
        the file name that exists in the file system will be returned.

        :param prefix_list: A Python :py:class:`list` object of Python :py:class:`str` (prefix) objects.
        :type prefix_list: list[str] | None
        :param suffix: A file name suffix.
        :type suffix: str
        :return: A file name.
        :rtype: str | None
        """
        if prefix_list is None:
            prefix_list = [self.prefix]

        for prefix in prefix_list:
            # Preferentially test with the genome version.
            file_name = '_'.join((self.project_name, self.genome_version, prefix, suffix))

            module_logger.debug('Checking for annotation sheet: %r', file_name)

            if os.path.exists(file_name):
                return file_name

            # Fall-back test without the genome version.
            file_name = '_'.join((self.project_name, prefix, suffix))

            module_logger.debug('Checking annotation sheet: %r', file_name)

            if os.path.exists(file_name):
                return file_name

        return

    def get_stage(self, name):
        """Get a :py:class:`bsf.analysis.Stage` object from a :py:class:`bsf.analysis.Analysis` object.

        If the :py:class:`bsf.analysis.Stage` object does not exist, it is created and initialised via the
        :py:class:`bsf.standards.Configuration` object in :py:attr:`bsf.analysis.Analysis.configuration`.

        Reads from configuration file sections
        :literal:`[bsf.analysis.Stage]`
        :literal:`[bsf.analysis.Analysis.Stage]` or :literal:`[bsf.analyses.*.Stage]`
        :literal:`[bsf.analysis.Analysis.Stage.name]` or :literal:`[bsf.analyses.*.Stage.name]`.

        :param name: A name.
        :type name: str
        :return: A :py:class:`bsf.analysis.Stage` object.
        :rtype: Stage
        """
        # Check if a Stage with this name already exists and if so, return it.
        for stage in self.stage_list:
            if stage.name == name:
                return stage

        # Initialise a new Stage and add it to the Python list of Stage objects.

        stage = Stage(name=name, working_directory=self.genome_directory)
        self.stage_list.append(stage)

        # A "bsf.analysis.Stage" section specifies defaults for all Stage objects of an Analysis.

        section = Configuration.section_from_instance(instance=stage)
        stage.set_configuration(configuration=self.configuration, section=section)

        module_logger.log(logging.DEBUG - 1, 'Stage configuration section: %r', section)

        # A "bsf.analysis.Analysis.Stage" or "bsf.analyses.*.Stage" pseudo-class section specifies
        # Analysis-specific or sub-class-specific options for the Stage, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=self), 'Stage'))
        stage.set_configuration(configuration=self.configuration, section=section)

        module_logger.log(logging.DEBUG - 1, 'Stage configuration section: %r', section)

        # A "bsf.analysis.Analysis.Stage.name" or "bsf.analyses.*.Stage.name" section specifies defaults
        # for a particular Stage of an Analysis or subclass, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=self), 'Stage', stage.name))
        stage.set_configuration(configuration=self.configuration, section=section)

        module_logger.log(logging.DEBUG - 1, 'Stage configuration section: %r', section)

        return stage

    def set_configuration(self, configuration, section):
        """Set instance variables of a :py:class:`bsf.analysis.Analysis` object
        via a section of a :py:class:`bsf.standards.Configuration` object.

        Instance variables without a configuration option remain unchanged.

        :param configuration: A :py:class:`bsf.standards.Configuration` object.
        :type configuration: Configuration
        :param section: A configuration file section.
        :type section: str
        :raise Exception: If the specified section does not exist.
        """
        if not configuration.config_parser.has_section(section=section):
            raise Exception(f'Section {section!r} is not defined in Configuration files:\n'
                            f'{configuration.file_path_list!r}')

        # The configuration section is available.

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

        option = 'report_style_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.report_style_path = configuration.config_parser.get(section=section, option=option)

        option = 'report_header_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.report_header_path = configuration.config_parser.get(section=section, option=option)

        option = 'report_footer_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.report_footer_path = configuration.config_parser.get(section=section, option=option)

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
        """Set default :py:class:`bsf.argument.Argument` objects for a :py:class:`bsf.process.Command` object.

        :param command: A :py:class:`bsf.process.Command` object.
        :type command: Command
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

        module_logger.log(logging.DEBUG - 1, 'Command configuration section: %r', section)

        command.set_configuration(configuration=self.configuration, section=section)

        return

    def set_runnable_step_configuration(self, runnable_step, tag=None):
        """Set default :py:class:`bsf.argument.Argument` objects for a :py:class:`bsf.process.RunnableStep` object.

        This method reads configuration section(s)
        :literal:`"Analysis.__class__.__name__"."RunnableStep.name"[."Command.name"]*`

        :param runnable_step: A :py:class:`bsf.process.RunnableStep` object.
        :type runnable_step: RunnableStep
        :param tag: An optional tag appended to the last sub-:py:class:`bsf.process.Command` of this
            :py:class:`.bsf.process.RunnableStep` object.
            This is useful to configure ChIP-seq peak calling by the ChIP factor (e.g., H3K27me3).
        :type tag: str | None
        """

        def _set_configuration(command, section):
            """Recursively set default :py:class:`bsf.argument.Argument` objects for a
            :py:class:`bsf.process.RunnableStep` object.

            The method sets the default :py:class:`bsf.argument.Argument` objects for a
            :py:class:`bsf.process.RunnableStep` object,
            as well as for its contained sub :py:class:`bsf.process.Command` objects.

            :param command: A :py:class:`bsf.process.Command` object.
            :type command: Command
            :param section: A Configuration section.
            :type section: str
            """
            if command.name:
                section += '.' + command.name
            else:
                section += '.' + command.program

            module_logger.log(logging.DEBUG - 1, 'RunnableStep configuration section: %r', section)

            command.set_configuration(configuration=self.configuration, section=section)

            if command.sub_command is None:
                # If no more sub-Command is available, check for a final tag and configure it directly.
                if tag:
                    section += '.' + tag
                    command.set_configuration(configuration=self.configuration, section=section)
            else:
                # If a sub-Command is available, recurse.
                _set_configuration(command=command.sub_command, section=section)

            return

        # Initially, the configuration section prefix is based on the Analysis class name and the Analysis Stage.name.

        _set_configuration(command=runnable_step, section=self.configuration.section_from_instance(instance=self))

        return

    def set_stage_runnable(self, stage, runnable):
        """Create a :py:class:`bsf.process.Executable` object to assign a
        :py:class:`bsf.procedure.Runnable` object to a :py:class:`bsf.analysis.Stage` object.

        In case the file in :py:meth:`bsf.procedure.Runnable.get_relative_status_path` exists already,
        :py:attr:`bsf.process.Executable.submit` will be set to :py:class:`False`.

        :param stage: A :py:class:`bsf.analysis.Stage` object.
        :type stage: Stage
        :param runnable: A :py:class:`bsf.procedure.Runnable` object.
        :type runnable: Runnable
        :return: A :py:class:`bsf.process.Executable` object.
        :rtype: Executable
        :raise Exception: A :py:attr:`bsf.procedure.Runnable.name` object does not exist in the
            :py:attr:`bsf.analysis.Analysis.runnable_dict` object.
        :raise Exception: A :py:class:`bsf.analysis.Stage` object does not exist in the
            :py:attr:`bsf.analysis.Analysis.stage_list` attribute.
        """
        if stage not in self.stage_list:
            raise Exception(f'A Stage.name {stage.name!r} does not exist in Analysis {self.project_name!r}.')

        if runnable.name not in self.runnable_dict:
            raise Exception(f'A Runnable.name {runnable.name!r} does not exist in Analysis {self.project_name!r}.')

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
        """Run a :py:class:`bsf.analysis.Analysis` object.

        :raise Exception: An :py:attr:`bsf.analysis.Analysis.project_name` has not been defined
        """
        if not self.project_name:
            raise Exception(f"A {self.name!s} requires a 'project_name' configuration option.")

        # Some analyses such as FastQC do not require a genome_version,
        # nor a genome_version-specific output directory.
        # Also, add the e-mail address for UCSC track hubs into the genome subclass.

        # Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        # expand any environment variables (i.e., on UNIX ${NAME} or $NAME)
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
            raise Exception(f"The 'output_directory' {self.output_directory!r} does not exist.")

        # Define project_directory and genome_directory instance variables.
        # If a genome_version option is present, append
        # it to the project_directory instance variable.
        # This allows analyses run against more than one directory and
        # simplifies UCSC Genome Browser track hub creation.

        self.project_directory = self.configuration.get_absolute_path(
            file_path=self.project_directory,
            default_path=os.path.join(self.output_directory, self.project_name))

        if self.genome_version:
            self.genome_directory = self.configuration.get_absolute_path(
                file_path=self.genome_directory,
                default_path=os.path.join(self.project_directory, self.genome_version))
        else:
            self.genome_directory = self.project_directory

        if not os.path.isdir(self.genome_directory):
            try:
                os.makedirs(self.genome_directory)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise exception

        if self.report_style_path:
            self.report_style_path = self.configuration.get_absolute_path(
                file_path=self.report_style_path,
                default_path=StandardFilePath.get_template_documents(absolute=True))

        if self.report_header_path:
            self.report_header_path = self.configuration.get_absolute_path(
                file_path=self.report_header_path,
                default_path=StandardFilePath.get_template_documents(absolute=True))

        if self.report_footer_path:
            self.report_footer_path = self.configuration.get_absolute_path(
                file_path=self.report_footer_path,
                default_path=StandardFilePath.get_template_documents(absolute=True))

        if not self.e_mail:
            self.e_mail = Operator.get_e_mail()
            if not self.e_mail:
                raise Exception(f"A {self.name!s} requires an 'e_mail' configuration option.")

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

            module_logger.log(logging.DEBUG - 1, 'Collection name: %r', self.collection.name)
            module_logger.log(logging.DEBUG - 2, 'Collection: %r', self.collection)
        else:
            # Create an empty bsf.ngs.Collection.
            self.collection = Collection()

        return

    def report(self):
        """Create a :literal:`XHTML 1.0` report.

        The method must be fully implemented in a subclass.
        """
        module_logger.warning(
            "The 'report' method must be implemented in the %r sub-class.",
            self.__class__.__name__)

        return

    @staticmethod
    def get_html_anchor(prefix, suffix, text):
        """Get a :literal:`XHTML 1.0` :literal:`anchor` element with a relative reference path.

        :literal:`<a href="prefix/prefix_suffix">text</a>`

        :param prefix: A prefix.
        :type prefix: str
        :param suffix: A suffix.
        :type suffix: str
        :param text: A link text.
        :type text: str
        :return: A :literal:`XHTML 1.0` :literal:`anchor` element.
        :rtype: str
        """
        return '<a href="' + prefix + '/' + prefix + '_' + suffix + '">' + text + '</a>'

    @staticmethod
    def get_html_image(prefix, suffix, text, height=None, width=None):
        """Get a :literal:`XHTML 1.0` :literal:`img` element with a relative source path.

        :literal:`<img alt="text" src="prefix/prefix_suffix" height="80" width="80" />`

        :param prefix: A prefix.
        :type prefix: str
        :param suffix: A suffix.
        :type suffix: str
        :param text: An alternative text.
        :type text: str
        :param height: An image height attribute.
        :type height: str
        :param width: An image width attribute.
        :type width: str
        :return: A :literal:`XHTML 1.0` :literal:`img` element.
        :rtype: str
        """
        str_list: List[str] = list()

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
        """Get a genome description :literal:`XHTML 1.0` paragraph.

        :param genome_version: A genome assembly version.
        :type genome_version: str | None
        :return: A genome description :literal:`XHTML 1.0` paragraph.
        :rtype: list[str]
        """
        str_list: List[str] = list()

        if genome_version is None:
            return str_list

        description = Genome.get_description(genome_version=genome_version)
        if description:
            str_list.append('<p>')
            str_list.append('<strong>Genome:</strong> ')
            str_list.append(description)
            str_list.append('</p>\n')

            return str_list

        species = Genome.get_species(genome_version=genome_version)
        # Without species information, the description is not useful.
        if species is None:
            module_logger.warning(
                'No species information for genome version %r in the central configuration file.',
                genome_version)

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
        """Get a transcriptome description :literal:`XHTML 1.0` paragraph.

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str | None
        :return: A transcriptome description :literal:`XHTML 1.0` paragraph.
        :rtype: list[str]
        """
        str_list: List[str] = list()

        if transcriptome_version is None:
            return str_list

        description = Transcriptome.get_description(transcriptome_version=transcriptome_version)
        if description:
            str_list.append('<p>')
            str_list.append('<strong>Transcriptome:</strong> ')
            str_list.append(description)
            str_list.append('</p>\n')

            return str_list

        species = Transcriptome.get_species(transcriptome_version=transcriptome_version)
        # Without species information, the description is not useful.
        if species is None:
            module_logger.warning(
                'No species information for transcriptome version %r in the central configuration file.',
                transcriptome_version)

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
        """Get a header section of an :literal:`XHTML 1.0` document.

        :param strict: Either a :literal:`XHTML 1.0 Strict` or a :literal:`XHTML 1.0 Transitional`
            Document Type Declaration, defaults to :literal:`XHTML 1.0 Strict`.
        :type strict: bool
        :param creator: A Dublin Core :literal:`DC.Creator` meta field value,
            defaults to the :literal:`USER` environment variable.
        :type creator: str
        :param source: A Dublin Core :literal:`DC.Source` meta field value,
            defaults to the Python script file path.
        :type source: str
        :param title: A title element value,
            defaults to a concatenation of the
            :py:attr:`bsf.analysis.Analysis.project_name` and
            :py:attr:`bsf.analysis.Analysis.report_name` attributes.
        :type title: str
        :return: A :literal:`XHTML 1.0` header section as Python :py:class:`list` of Python :py:class:`str` objects.
        :rtype: list[str]
        """
        if creator is None or not creator:
            creator = getpass.getuser()
            # The getpass.getuser method just relies on environment variables,
            # but at least works under Unix and Windows.

        if source is None or not source:
            source = inspect.getfile(inspect.currentframe())

        if title is None or not title:
            title = ' '.join((self.project_name, self.name))

        str_list: List[str] = list()

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
        str_list.append('<link rel="schema.DC" href="http://purl.org/DC/elements/1.0/" />\n')
        str_list.append('<meta name="DC.Creator" content="' + html.escape(s=creator, quote=True) + '" />\n')
        str_list.append('<meta name="DC.Date" content="' + datetime.datetime.now().isoformat() + '" />\n')
        str_list.append('<meta name="DC.Source" content="' + html.escape(s=source, quote=True) + '" />\n')
        str_list.append('<meta name="DC.Title" content="' + html.escape(s=title, quote=True) + '" />\n')
        str_list.append('<style type="text/css">\n')

        # If available, import a custom CSS document.
        if self.report_style_path:
            with open(file=self.report_style_path, mode='rt') as input_text_io:
                str_list.extend(input_text_io)
        else:
            str_list.append('  .left    { text-align: left; }\n')
            str_list.append('  .right   { text-align: right; }\n')
            str_list.append('  .center  { text-align: center; }\n')
            str_list.append('  .justify { text-align: justify; }\n')
            str_list.append('  .start   { text-align: start; }\n')
            str_list.append('  .end     { text-align: end; }\n')
            str_list.append('  body { font-family: sans-serif; }\n')
            str_list.append('  h1 { color: #40B9D4; }\n')
            str_list.append('  a { color: #40B9D4; text-decoration: none; }\n')
            str_list.append('  a:hover, a:focus { color: #2A6496; text-decoration: underline; }\n')

        str_list.append('</style>\n')
        str_list.append('<title>' + html.escape(s=title, quote=True) + '</title>\n')
        str_list.append('</head>\n')
        str_list.append('\n')
        str_list.append('<body>\n')
        str_list.append('\n')

        # Include additional lines from an optional, Analysis-specific report header template.
        if self.report_header_path:
            with open(file=self.report_header_path, mode='rt') as input_text_io:
                str_list.extend(input_text_io)

        return str_list

    def get_html_footer(
            self,
            contact=None,
            institution=None,
            url_protocol=None,
            url_host_name=None,
            title=None):
        """Get a footer section of an :literal:`XHTML 1.0` document.

        :param contact: An institution contact e-mail address,
            defaults to the value returned by the :py:meth:`bsf.standards.Operator.get_contact` method.
        :type contact: str
        :param institution: An institution name to be inserted into 'This report was generated by ...',
            defaults to the value returned by the :py:meth:`bsf.standards.Operator.get_institution` method.
        :type institution: str
        :param url_protocol: A protocol section of the institution URL (i.e., http, https, ...),
            defaults to the value returned by the :py:meth:`bsf.standards.URL.get_protocol` method.
        :type url_protocol: str
        :param url_host_name: A host name section of the institution URL (e.g., biomedical-sequencing.at),
            defaults to the value returned by :py:meth:`bsf.standards.URL.get_host_name` method.
        :type url_host_name: str
        :param title: A title element value,
            defaults to a concatenation of the
            :py:attr:`bsf.analysis.Analysis.project_name` and
            :py:attr:`bsf.analysis.Analysis.report_name` attributes.
        :type title: str
        :return: A :literal:`XHTML 1.0` footer section as Python :py:class:`list` of Python :py:class:`str` objects.
        :rtype: list[str]
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

        str_list: List[str] = list()

        # Include additional lines from an optional, Analysis-specific report footer template.
        if self.report_footer_path:
            with open(file=self.report_footer_path, mode='rt') as input_text_io:
                str_list.extend(input_text_io)

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
            # The e-mail address outside a URL still needs HTML quoting.
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
        """Get a report as an :literal:`XHTML 1.0` document.

        The method automatically concatenates the
        :literal:`XHTML 1.0` header :py:meth:`bsf.analysis.Analysis.get_html_header`, the
        :literal:`XHTML 1.0` content and the
        :literal:`XHTML 1.0` footer :py:meth:`bsf.analysis.Analysis.get_html_footer` before returning the report.

        :param content: A :literal:`XHTML 1.0` content.
        :type content: list[str]
        :param strict: Either a :literal:`XHTML 1.0 Strict` or a :literal:`XHTML 1.0 Transitional`
            Document Type Declaration, defaults to :literal:`XHTML 1.0 Strict`.
        :type strict: bool
        :param creator: A Dublin Core :literal:`DC.Creator` meta field value,
            defaults to the :literal:`USER` environment variable.
        :type creator: str
        :param source: A Dublin Core :literal:`DC.Source` meta field value,
            defaults to the Python script file path.
        :type source: str
        :param title: A title element value,
            defaults to a concatenation of the
            :py:attr:`bsf.analysis.Analysis.project_name` and
            :py:attr:`bsf.analysis.Analysis.report_name` attributes.
        :type title: str
        :param contact: An institution contact e-mail address,
            defaults to the value returned by the :py:meth:`bsf.standards.Operator.get_contact` method.
        :type contact: str
        :param institution: An institution name to be inserted into 'This report was generated by ...',
            defaults to the value returned by the :py:meth:`bsf.standards.Operator.get_institution` method.
        :type institution: str
        :param url_protocol: A protocol section of the institution URL (i.e., http, https, ...),
            defaults to the value returned by the :py:meth:`bsf.standards.URL.get_protocol` method.
        :type url_protocol: str
        :param url_host_name: A host name section of the institution URL (e.g., biomedical-sequencing.at),
            defaults to the value returned by the :py:meth:`bsf.standards.URL.get_host_name` attribute.
        :type url_host_name: str
        :return: A :literal:`XHTML 1.0` report as Python :py:class:`list` of Python :py:class:`str` objects.
        :rtype: list[str]
        """
        str_list: List[str] = list()

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
        """Write an :literal:`XHTML 1.0` report :literal:`prefix_report.html` file into the
        :py:attr:`bsf.analysis.Analysis.genome_directory`.

        The method automatically concatenates the
        :literal:`XHTML 1.0` header :py:meth:`bsf.analysis.Analysis.get_html_header`, the
        :literal:`XHTML 1.0` content and the
        :literal:`XHTML 1.0` footer :py:meth:`bsf.analysis.Analysis.get_html_footer` before writing the file.

        :param content: A :literal:`XHTML 1.0` content.
        :type content: list[str]
        :param prefix: A file name prefix (e.g., chipseq, rnaseq, ...),
            defaults to the :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix: str
        :param strict: Either a :literal:`XHTML 1.0 Strict` or a :literal:`XHTML 1.0 Transitional`
            Document Type Declaration, defaults to :literal:`XHTML 1.0 Strict`.
        :type strict: bool
        :param creator: A Dublin Core :literal:`DC.Creator` meta field value,
            defaults to the :literal:`USER` environment variable.
        :type creator: str
        :param source: A Dublin Core :literal:`DC.Source` meta field value,
            defaults to the Python script file path.
        :type source: str
        :param title: A title element value,
            defaults to a concatenation of the
            :py:attr:`bsf.analysis.Analysis.project_name` and
            :py:attr:`bsf.analysis.Analysis.report_name` attributes.
        :type title: str
        :param contact: An institution contact e-mail address,
            defaults to the value returned by the :py:meth:`bsf.standards.Operator.get_contact` method.
        :type contact: str
        :param institution: An institution name to be inserted into 'This report was generated by ...',
            defaults to the value returned by the :py:meth:`bsf.standards.Operator.get_institution` method.
        :type institution: str
        :param url_protocol: A protocol section of the institution URL (i.e., http, https, ...),
            defaults to the value returned by the :py:meth:`bsf.standards.URL.get_protocol` method.
        :type url_protocol: str
        :param url_host_name: A host name section of the institution URL (e.g., biomedical-sequencing.at),
            defaults to the value returned by the :py:class:`bsf.standards.URL.get_host_name` method.
        :type url_host_name: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        with open(
                file=os.path.join(self.genome_directory, '_'.join((prefix, 'report.html'))),
                mode='wt') as output_text_io:
            output_text_io.writelines(
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

    def create_public_project_link(self, sub_directory=None, reset=None):
        """Create a symbolic link from the public HTML directory to the project directory if not already there.

        The link will be placed in the specified subdirectory under :py:meth:`StandardFilePath.get_public_html` and
        contain the project name followed by a 128 bit hexadecimal UUID string.
        If not specified, the subdirectory defaults to the value of :py:meth:`bsf.standards.URL.get_relative_projects`.

        All symbolic links are checked to identify an already existing one and
        any dangling links encountered are reported.

        :param sub_directory: A :py:class:`bsf.analysis.Analysis`-specific directory.
        :type sub_directory: str
        :param reset: Replace an existing symbolic link with a new one.
        :type reset: bool | None
        :return: A symbolic link to the project directory.
        :rtype: str
        :raise Exception: Public HTML path does not exist
        """
        # Cache the results since reading directories and checking symbolic links is quite involved.
        # The public project path consists of the absolute public_html directory,
        # the analysis-specific subdirectory, the project name and a 128 bit hexadecimal UUID string.

        if self._public_project_link_path and not reset:
            return self._public_project_link_path

        if sub_directory is None:
            sub_directory = URL.get_relative_projects()

        html_path = os.path.join(StandardFilePath.get_public_html(absolute=True), sub_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the html_path directory has to exist.

        if not os.path.isdir(html_path):
            raise Exception(
                f'The public HTML directory path {html_path!r} does not exist.\n'
                f'Please check the optional sub-directory name {sub_directory!r}.')

        # Check all symbolic links to find one already pointing to the source path and
        # to identify any danging symbolic links.

        for target_name in os.listdir(html_path):
            target_path = os.path.join(html_path, target_name)
            if os.path.islink(target_path):
                source_path = os.readlink(target_path)
                if not os.path.isabs(source_path):
                    source_path = os.path.join(html_path, source_path)
                source_path = os.path.normpath(source_path)

                if not os.path.exists(source_path):
                    # For os.path.samefile both paths have to exist.
                    # Check the source path and report dangling symbolic links for any analysis project.
                    module_logger.warning('Dangling symbolic link %r to %r.', target_path, source_path)
                    continue

                if os.path.samefile(source_path, self.project_directory):
                    # Set the target path to the already existing file_path, but report duplicate symbolic links.
                    # Do not break out here to discover all dangling symbolic links.
                    if self._public_project_link_path:
                        module_logger.warning(
                            'More than one symbolic link points to this project directory. primary: %r secondary: %r',
                            target_path,
                            self._public_project_link_path)
                    else:
                        self._public_project_link_path = target_path

        if reset:
            try:
                os.remove(self._public_project_link_path)
            except OSError as exception:
                if exception.errno != errno.ENOENT:
                    raise exception

            self._public_project_link_path = None

        if not self._public_project_link_path:
            self._public_project_link_path = os.path.join(html_path, '_'.join((self.project_name, uuid.uuid4().hex)))

            try:
                os.symlink(os.path.relpath(self.project_directory, html_path), self._public_project_link_path, True)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise exception

        return self._public_project_link_path

    def get_html_report_url(self, sub_directory=None, prefix=None):
        """Get an HTML report URL.

        :param sub_directory: A :py:class:`bsf.analysis.Analysis`-specific directory.
        :type sub_directory: str
        :param prefix: A file name prefix (e.g., chipseq, rnaseq, ...),
            defaults to the :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix: str
        :return:
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        link_name = os.path.basename(self.create_public_project_link(sub_directory=sub_directory))

        if self.genome_version:
            return f'{URL.get_absolute_projects()}/{link_name}/{self.genome_version}/{prefix}_report.html'
        else:
            return f'{URL.get_absolute_projects()}/{link_name}/{prefix}_report.html'

    def ucsc_track_url(
            self,
            options_dict=None,
            browser_dict=None,
            track_dict=None,
            ucsc_protocol=None,
            ucsc_host_name=None):
        """Return a URL to automatically attach a UCSC Genome Browser track.

        :param options_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (URL option key value pair) objects.
        :type options_dict: dict[str, str]
        :param browser_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (browser line key value pair) objects.
        :type browser_dict: dict[str, str]
        :param track_dict: A Python :py:class:`dict` object of
            Python :py:class:`str` (track line hgct_customText) key value pair objects.
        :type track_dict: dict[str, str]
        :param ucsc_protocol: A UCSC Genome Browser URL protocol (i.e., http, https, ...)
            defaults to the value returned by the :py:meth:`bsf.standards.UCSC.get_protocol` method.
        :type ucsc_protocol: str
        :param ucsc_host_name: A UCSC Genome Browser URL host name,
            defaults to the value returned by the :py:meth:`bsf.standards.UCSC.get_host_name` method.
        :type ucsc_host_name: str
        :return: A URL to attach a track to the UCSC Genome Browser.
        :rtype: str
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

        :param link_path: A symbolic link path in the public HTML directory including project name and a UUID.
        :type link_path: str
        :param options_dict: A Python :py:class:`dict` of
            Python :py:class:`str` (URL option key) key and
            Python :py:class:`str` (URL option value) value objects.
        :type options_dict: dict[str, str]
        :return: A URL to automatically attach a UCSC Genome Browser Track Hub.
        :rtype: str
        """
        if options_dict is None:
            options_dict = dict()

        if 'hubUrl' not in options_dict:
            # The track hub URL requires the link name (i.e., the link path base name) to be inserted.
            link_name = os.path.basename(link_path.rstrip('/'))
            options_dict['hubUrl'] = '/'.join((
                URL.get_absolute_projects(),
                link_name,
                '_'.join((self.prefix, self.ucsc_name_hub))))

        return self.ucsc_track_url(options_dict=options_dict)

    def ucsc_hub_html_anchor(self, link_path):
        """Return an :literal:`XHTML 1.0` anchor element to automatically attach a UCSC Genome Browser Track Hub.

        See also the :py:meth:`bsf.analysis.Analysis.create_public_project_link` method.

        :param link_path: A symbolic link path in the public HTML directory including project name and a UUID.
        :type link_path: str
        :return: A :literal:`XHTML 1.0` :literal:`anchor` element as
            Python :py:class:`list` of Python :py:class:`str` objects.
        :rtype: list[str]
        """
        str_list: List[str] = list()

        str_list.append('UCSC Genome Browser Track Hub ')
        str_list.append('<a href="' + self.ucsc_hub_url(link_path=link_path) + '">' + self.project_name + '</a>')

        return str_list

    def ucsc_hub_write_hub(self, prefix=None):
        """Write a UCSC Track Hub :literal:`prefix_hub.txt` file into the
        :py:attr:`bsf.analysis.Analysis.project_directory`.

        The :py:attr:`bsf.analysis.Analysis.project_directory` is one level above
        the :py:attr:`bsf.analysis.Analysis.genome_directory`.

        :param prefix: A hub prefix (e.g., chipseq, rnaseq, ...).
        :type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        str_list: List[str] = list()

        str_list.append('hub ' + '_'.join((self.project_name, prefix)) + '\n')
        str_list.append('shortLabel ' + '_'.join((self.project_name, prefix)) + '\n')
        str_list.append('longLabel Project ' + '_'.join((self.project_name, prefix)) + '\n')
        str_list.append('genomesFile ' + '_'.join((prefix, self.ucsc_name_genomes)) + '\n')
        str_list.append('email ' + self.e_mail + '\n')

        with open(
                file=os.path.join(self.project_directory, '_'.join((prefix, self.ucsc_name_hub))),
                mode='wt') as output_text_io:
            output_text_io.writelines(str_list)

        return

    def ucsc_hub_write_genomes(self, prefix=None):
        """Write a UCSC Track Hub :literal:`prefix_genomes.txt` file into the
        :py:attr:`bsf.analysis.Analysis.project_directory`.

        The :py:attr:`bsf.analysis.Analysis.project_directory` is one level above
        the :py:attr:`bsf.analysis.Analysis.genome_directory`.

        :param prefix: A hub prefix (e.g., chipseq, rnaseq, ...).
        :type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        # If the file exists, read it first to retain any other genome assembly entries.
        genome_version_dict: Dict[str, str] = dict()

        file_path = os.path.join(self.project_directory, '_'.join((prefix, self.ucsc_name_genomes)))

        if os.path.exists(file_path):
            genome_version = None
            with open(file=file_path, mode='rt') as input_text_io:
                for line_str in input_text_io:
                    line_str = line_str.strip()
                    if not line_str:
                        continue
                    field_list = line_str.split()

                    if len(field_list) != 2:
                        module_logger.warning(
                            'Malformed line %r in UCSC genomes file %r. '
                            'Expected exactly two components after line splitting.', line_str, file_path)

                    if field_list[0] == 'genome':
                        if genome_version is not None:
                            module_logger.warning(
                                "Malformed line %r in UCSC genomes file %r. "
                                "Got more than one 'genomes' lines in succession.", line_str, file_path)
                        genome_version = field_list[1]

                    if field_list[0] == 'trackDb':
                        if genome_version is None:
                            module_logger.warning(
                                "Malformed line %r in UCSC genomes file %r. "
                                "Got a 'trackDb' line without a preceding genomes line.", line_str, file_path)
                        else:
                            genome_version_dict[genome_version] = field_list[1]
                            genome_version = None

        # Resolve an eventual alias for the UCSC genome assembly name in "genome_version/prefix_tracks.txt".
        genome_version_dict[Genome.resolve_ucsc_alias(genome_version=self.genome_version)] = \
            '/'.join((self.genome_version, '_'.join((prefix, self.ucsc_name_tracks))))

        str_list: List[str] = list()

        for genome_version in sorted(genome_version_dict):
            str_list.append('genome ' + genome_version + '\n')
            str_list.append('trackDb ' + genome_version_dict[genome_version] + '\n')
            str_list.append('\n')

        with open(file=file_path, mode='wt') as output_text_io:
            output_text_io.writelines(str_list)

        return

    def ucsc_hub_write_tracks(self, content, prefix=None):
        """Write a UCSC Track Hub :literal:`prefix_tracks.txt` file into the
        :py:attr:`bsf.analysis.Analysis.genome_directory`.

        :param content: Content for a UCSC Track Hub.
        :type content: list[str]
        :param prefix: A hub prefix (e.g., chipseq, rnaseq, ...).
        :type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        with open(
                file=os.path.join(self.genome_directory, '_'.join((prefix, self.ucsc_name_tracks))),
                mode='wt') as output_text_io:
            output_text_io.writelines(content)

        return

    def ucsc_hub_to_file(self, content, prefix=None):
        """Write UCSC Genome Browser Track Hub files to disk.

        The method writes :literal:`prefix_hub.txt` and :literal:`prefix_genomes.txt` files into
        the :py:attr:`bsf.analysis.Analysis.project_directory}, above
        the :py:attr:`bsf.analysis.Analysis.genome_directory`, as well as a
        :literal:`prefix_tracks.txt` file into
        the :py:attr:`bsf.analysis.Analysis.genome_directory`.

        :param content: Content for the track database file.
        :type content: list[str]
        :param prefix: A hub prefix (e.g., chipseq, rnaseq, ...),
            defaults to the value of the :py:attr:`bsf.analysis.Analysis.prefix` attribute.
        :type prefix: str
        """
        if prefix is None or not prefix:
            prefix = self.prefix

        self.ucsc_hub_write_hub(prefix=prefix)
        self.ucsc_hub_write_genomes(prefix=prefix)
        self.ucsc_hub_write_tracks(content=content, prefix=prefix)

        return

    def ucsc_hub_bigwig_info_signal_range(self, file_path):
        """Read the bigWig signal range from a bigWig information file and return a
        UCSC Track Hub :literal:`type bigWig` line.

        :param file_path: A UCSC bigWigInfo file path.
        :type file_path: str
        :return: A UCSC Track Hub :literal:`type bigWig` line with optional signal range values.
        :rtype: str
        """
        if os.path.isabs(file_path):
            file_path_absolute = file_path
        else:
            file_path_absolute = os.path.join(self.genome_directory, file_path)

        if os.path.exists(file_path_absolute):
            minimum = None
            maximum = None

            with open(file=file_path_absolute, mode='rt') as input_text_io:
                for line in input_text_io:
                    if line.startswith('min:'):
                        line_list = line.split()
                        minimum = line_list[1].strip()
                    if line.startswith('max:'):
                        line_list = line.split()
                        maximum = line_list[1].strip()

            return 'type bigWig ' + minimum + ' ' + maximum + '\n'
        else:
            return 'type bigWig\n'

    def check_state(self):
        """Check the state of each :py:class:`bsf.analysis.Stage` object.
        """
        for stage in self.stage_list:
            stage.check_state()

        return

    def submit(self, name=None, drms_submit=None):
        """Submit each :py:class:`bsf.analysis.Stage` object.

        Submits each :py:class:`bsf.process.Executable` of either all :py:class:`bsf.analysis.Stage` objects or
        a named one and pickles each :py:class:`bsf.procedure.Runnable` object.

        :param name: Only submit :py:class:`bsf.process.Executable` objects linked to this
            :py:attr:`bsf.analysis.Stage.name` attribute.
        :type name: str
        :param drms_submit: Submit to the :emphasis:`Distributed Resource Management System` (DRMS).
        :type drms_submit: bool | None
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
            stage.submit(drms_submit=drms_submit)

            module_logger.debug('Stage.name: %r', stage.name)
            module_logger.log(logging.DEBUG - 2, 'Stage: %r', stage)

        if name:
            if name == 'report':
                self.report()
            elif not submit:
                name_list = [stage.name for stage in self.stage_list]
                name_list.append('report')
                module_logger.warning('Valid Analysis Stage names are: %r', name_list)

        return
