"""Bio.BSF

A package of classes and methods specific to the Biomedical Sequencing Facility (BSF).
Reference: http://www.biomedical-sequencing.at/
"""

#
# Copyright 2013 Michael K. Schuster
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


from ConfigParser import SafeConfigParser
import datetime
import errno
import importlib
import os
import re
from stat import *
import string
from subprocess import PIPE, Popen
import sys
from threading import Lock, Thread
import uuid
import warnings

from Bio.BSF import Defaults
from Bio.BSF.Data import Collection, Sample, SampleGroup
from Bio.BSF.Argument import *


class Analysis(object):
    """BSF Analysis class.

    The BSF Analysis class represents a high-level analysis that may use one or more
    algorithms (BSF Executable).

    Attributes:
    :ivar configuration: BSF Configuration
    :type configuration: Configuration
    :ivar debug: Debug level
    :type debug: int
    :ivar project_name: Project name (arbitrary)
    :type project_name: str
    :ivar genome_version: Genome version (e.g. hg19, mm10, GRCh37, GRCm38, ...)
    :type genome_version: str
    :ivar input_directory: Input directory
    :type input_directory: str, unicode
    :ivar output_directory: Output directory, user-specified including a genome version sub-directory
    :type output_directory: str, unicode
    :ivar project_directory: Project-specific directory
    :type project_directory: str, unicode
    :ivar genome_directory: Genome-specific directory
    :type genome_directory: str, unicode
    :ivar drms_list: Python list of BSF DRMS objects
    :type drms_list: list
    :ivar collection: BSF Collection
    :type collection: Collection
    :ivar comparisons: Python dict of comparisons
    :type comparisons: dict
    :ivar samples: Python list of BSF Sample objects
    :type samples: list
    """

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new BSF Analysis object from a UNIX-style configuration file via the BSF Configuration class.

        :param cls: Class
        :type cls: Analysis
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF Analysis
        :rtype: Analysis
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new BSF Analysis object from a BSF Configuration object.

        :param cls: Class
        :type cls: Analysis
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF Analysis
        :rtype: Analysis
        """

        assert isinstance(configuration, Configuration)

        # Set a minimal set of global defaults.

        default = Default.get_global_default()

        analysis = cls(configuration=configuration, e_mail=default.operator_e_mail)

        # A "Bio.BSF.Analysis.*" section specifies defaults for this BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        analysis.set_Configuration(analysis.configuration, section=section)

        return analysis

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 sas_file=None, sas_prefix=None, e_mail=None, debug=0, drms_list=None,
                 collection=None, comparisons=None, samples=None):

        """Initialise a Bio.BSF.Analysis object.

        :param self: BSF Analysis
        :type self: Analysis
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param project_name: Project name
        :type project_name: str
        :param genome_version: Genome version
        :type genome_version: str
        :param input_directory: BSF Analysis-wide input directory
        :type input_directory: str
        :param output_directory: BSF Analysis-wide output directory
        :type output_directory: str
        :param project_directory: BSF Analysis-wide project directory,
        normally under the BSF Analysis-wide output directory
        :type project_directory: str
        :param genome_directory: BSF Analysis-wide genome directory,
        normally under the BSF Analysis-wide project directory
        :type genome_directory: str
        :param sas_file: Sample Annotation Sheet (SAS) file path
        :type sas_file: str, unicode
        :param sas_prefix: A prefix to columns in a Sample Annotation Sheet
        (e.g. Control Sample, Treatment Sample, ...)
        :type sas_prefix: str
        :param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        :type e_mail: str
        :param debug: Integer debugging level
        :type debug: int
        :param drms_list: Python list of BSF DRMS objects
        :type drms_list: list
        :param collection: BSF Collection
        :type collection: Collection
        :param comparisons: Python dict of Analysis-specific objects
        (i.e. Python tuple for RNA-Seq and ChIPSeqComparison for ChIPSeq)
        :type comparisons: dict
        :param samples: Python list of BSF Sample objects
        :type samples: list
        :return: Nothing
        :rtype: None
        """

        if configuration:
            self.configuration = configuration
        else:
            self.configuration = Configuration()

        if project_name:
            self.project_name = project_name
        else:
            self.project_name = str()

        if genome_version:
            self.genome_version = genome_version
        else:
            self.genome_version = str()

        if input_directory:
            self.input_directory = input_directory
        else:
            self.input_directory = str()

        if output_directory:
            self.output_directory = output_directory
        else:
            self.output_directory = str()

        if project_directory:
            self.project_directory = project_directory
        else:
            self.project_directory = str()

        if genome_directory:
            self.genome_directory = genome_directory
        else:
            self.genome_directory = str()

        if sas_file:
            self.sas_file = sas_file
        else:
            self.sas_file = str()

        if sas_prefix:
            self.sas_prefix = sas_prefix
        else:
            self.sas_prefix = str()

        if e_mail:
            self.e_mail = e_mail
        else:
            self.e_mail = str()

        self.debug = debug

        if drms_list:
            self.drms_list = drms_list
        else:
            self.drms_list = list()

        self.collection = collection

        if comparisons:
            self.comparisons = comparisons
        else:
            self.comparisons = dict()

        if samples:
            self.samples = samples
        else:
            self.samples = list()

    def trace(self, level):

        """Trace a BSF Analysis object.

        :param self: BSF Analysis
        :type self: Analysis
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  project_name:     {!r}\n'.format(indent, self.project_name)
        output += '{}  genome_version:   {!r}\n'.format(indent, self.genome_version)
        output += '{}  input_directory:  {!r}\n'.format(indent, self.input_directory)
        output += '{}  output_directory: {!r}\n'.format(indent, self.output_directory)
        output += '{}  genome_directory: {!r}\n'.format(indent, self.genome_directory)
        output += '{}  sas_file:         {!r}\n'.format(indent, self.sas_file)
        output += '{}  sas_prefix        {!r}\n'.format(indent, self.sas_prefix)
        output += '{}  e_mail:           {!r}\n'.format(indent, self.e_mail)
        output += '{}  debug:            {!r}\n'.format(indent, self.debug)
        output += '{}  drms_list         {!r}\n'.format(indent, self.drms_list)
        output += '{}  collection        {!r}\n'.format(indent, self.collection)
        output += '{}  comparisons       {!r}\n'.format(indent, self.comparisons)
        output += '{}  samples           {!r}\n'.format(indent, self.samples)

        output += '{}  Python List of BSF Sample objects:'.format(indent)
        for sample in self.samples:
            output += '{}    BSF Sample name: {!r} file_path: {!r}'.format(indent, sample.name, sample.file_path)

        if self.collection:
            output += self.collection.trace(level + 1)

        return output

    def add_Sample(self, sample):

        """Add a BSF Sample object to the Python list of BSF Sample objects if it does not already exist.

        The check is based on the >Python 'in' comparison operator and in lack of a specific
        Bio.BSF.Data.__cmp__ method, relies on object identity (i.e. address).
        :param self: BSF Analysis
        :type self: Analysis
        :param sample: BSF Sample
        :type sample: Sample
        :return: Nothing
        :rtype: None
        """

        if not sample:
            return

        assert isinstance(sample, Sample) or isinstance(sample, SampleGroup)

        if sample not in self.samples:
            self.samples.append(sample)

    def set_Configuration(self, configuration, section):

        """Set instance variables of a BSF Analysis object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF Analysis
        :type self: Analysis
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        assert isinstance(configuration, Configuration)

        if not configuration.config_parser.has_section(section=section):
            message = 'Section {!r} not defined in BSF Configuration file {!r}.'. \
                format(section, configuration.config_file)
            warnings.warn(message, UserWarning)

            return

        # The configuration section is available.

        if configuration.config_parser.has_option(section=section, option='debug'):
            self.debug = configuration.config_parser.getint(section=section, option='debug')

        if configuration.config_parser.has_option(section=section, option='project_name'):
            self.project_name = configuration.config_parser.get(section=section, option='project_name')

        if configuration.config_parser.has_option(section=section, option='input_directory'):
            self.input_directory = configuration.config_parser.get(section=section, option='input_directory')

        if configuration.config_parser.has_option(section=section, option='output_directory'):
            self.output_directory = configuration.config_parser.get(section=section, option='output_directory')

        if configuration.config_parser.has_option(section=section, option='genome_version'):
            self.genome_version = configuration.config_parser.get(section=section, option='genome_version')

        if configuration.config_parser.has_option(section=section, option='sas_file'):
            self.sas_file = configuration.config_parser.get(section=section, option='sas_file')

        if configuration.config_parser.has_option(section=section, option='sas_prefix'):
            self.sas_prefix = configuration.config_parser.get(section=section, option='sas_prefix')

        if configuration.config_parser.has_option(section=section, option='e_mail'):
            self.e_mail = configuration.config_parser.get(section=section, option='e_mail')

    def run(self):

        """Run the BSF Analysis.

        :param self: BSF Analysis
        :type self: Analysis
        :return: Nothing
        :rtype: None
        """

        if not self.project_name:
            message = 'A BSF Analysis project_name has not been defined.'
            raise Exception(message)

        # Some analyses such as FastQC do not require a genome_version,
        # nor a genome_version-specific output directory.
        # Also, add the e-mail address for UCSC track hubs into the genome subclass.

        # Expand an eventual user part i.e. on UNIX ~ or ~user and
        # expand any environment variables i.e. on UNIX ${NAME} or $NAME
        # Check if an absolute path has been provided, if not,
        # automatically prepend standard BSF directory paths.

        self.input_directory = os.path.expanduser(path=self.input_directory)
        self.input_directory = os.path.expandvars(path=self.input_directory)

        if not os.path.isabs(self.input_directory):
            self.input_directory = os.path.join(Default.absolute_samples(), self.input_directory)

        self.output_directory = os.path.expanduser(path=self.output_directory)
        self.output_directory = os.path.expandvars(path=self.output_directory)

        if not os.path.isabs(self.output_directory):
            self.output_directory = os.path.join(Default.absolute_projects(), self.output_directory)

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

        if self.sas_file:

            # Populate a BSF Collection from a Sample Annotation Sheet.

            self.sas_file = os.path.expanduser(path=self.sas_file)
            self.sas_file = os.path.expandvars(path=self.sas_file)

            if not os.path.isabs(self.sas_file):
                self.sas_file = os.path.join(self.project_directory, self.sas_file)

            self.collection = Collection.from_sas(file_path=self.input_directory,
                                                  file_type='Automatic',
                                                  name=self.project_name,
                                                  sas_file=self.sas_file,
                                                  sas_prefix=self.sas_prefix)

            if self.debug:
                print '{!r} Collection name: {}'.format(self, self.collection.name)
                print self.collection.trace(1)

        else:

            # Create an empty BSF Collection.

            self.collection = Collection()

        self.create_project_genome_directory()

    def report(self):

        """Create a BSF Analysis report.

        :param self: BSF Analysis
        :type self: Analysis
        :return: Nothing
        :rtype: None
        """

        message = "The 'report' method must be implemented in the sub-class."
        warnings.warn(message, UserWarning)

    def create_project_genome_directory(self):

        """Check and create a BSF Analysis project_directory or genome_directory if necessary.

        :parm self: BSF Analysis
        :type self: Analysis
        :return: Nothing
        :rtype: None
        """

        if not os.path.isdir(self.genome_directory):
            message = 'Output (genome) directory {!r} does not exist. Create? [Y/n] '.format(self.genome_directory)
            answer = raw_input(message)

            if not answer or answer == 'Y' or answer == 'y':
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.makedirs(self.genome_directory)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
            else:
                message = 'Output (genome) directory {!r} does not exist.'.format(self.genome_directory)
                raise Exception(message)

    def create_public_project_link(self, sub_directory=None):

        """Create a symbolic link from the web directory to the project directory if not already there.

        The link will be placed in the sub directory and contain
        the project name followed by a 128 bit hexadecimal UUID string.
        :param self: BSF Analysis
        :type self: Analysis
        :param sub_directory: BSF Analysis-specific directory
        :type sub_directory: str
        :return: Symbolic link to the project directory
        :rtype: str
        """

        # The html_path consists of the absolute public_html directory and
        # the analysis-specific sub-directory.

        html_path = Default.absolute_public_html()

        if sub_directory:
            html_path = os.path.join(html_path, sub_directory)

        # Do not automatically create "new" paths based on mis-spelt sub_directory names.

        if not os.path.isdir(html_path):
            message = "Public HTML path {} does not exist. Check sub-directory {} name.". \
                format(html_path, sub_directory)
            raise Exception(message)

        # The link_name consists of the absolute public_html directory,
        # the analysis-specific sub-directory, the project name and a 128 bit hexadecimal UUID string.

        link_name = os.path.join(html_path, string.join([self.project_name, uuid.uuid4().hex], '_'))

        # While checking for already existing symbolic links,
        # the path_name holds the complete path for each link in the sub-directory.

        path_name = str()

        # The link_final holds the final symbolic link. It can be the link_name assembled above or
        # point to an already existing one.

        link_final = link_name

        link_exists = False

        for file_name in os.listdir(html_path):
            path_name = os.path.join(html_path, file_name)
            mode = os.lstat(path_name).st_mode
            if S_ISLNK(mode):
                target_name = os.readlink(path_name)
                if not os.path.isabs(target_name):
                    target_name = os.path.join(os.path.dirname(html_path), target_name)
                if not os.path.exists(path=target_name):
                    # Both paths for os.path.samefile have to exist.
                    # Check for dangling symbolic links.
                    message = 'Dangling symbolic link {} to {}'.format(path_name, target_name)
                    warnings.warn(message, UserWarning)
                    continue
                if os.path.samefile(target_name, self.project_directory):
                    link_exists = True
                    link_final = path_name  # Reset link_final to the already existing path_name.
                    break

        if link_exists:
            # Ask the user to re-create the symbolic link.
            message = 'Public HTML link {!r} to {!r} does exist. Re-create? [y/N] '. \
                format(path_name, self.project_directory)
            answer = raw_input(message)

            if not answer or answer == 'N' or answer == 'n':
                print 'Public HTML link {!r} to {!r} not reset.'. \
                    format(path_name, self.project_directory)
            else:
                try:
                    os.remove(path_name)
                except OSError as exception:
                    if exception.errno != errno.ENOENT:
                        raise
                    # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.symlink(self.project_directory, link_name)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
        else:
            # Ask the user to create a symbolic link.
            message = 'Public HTML link {!r} to {!r} does not exist. Create? [Y/n] '. \
                format(link_name, self.project_directory)
            answer = raw_input(message)

            if not answer or answer == 'Y' or answer == 'y':
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.symlink(self.project_directory, link_name)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
            else:
                print 'Public HTML link {!r} to {!r} not set.'. \
                    format(link_name, self.project_directory)

        return link_final

    def ucsc_hub_write_hub(self, prefix=None):

        """Write a UCSC Track Hub hub.txt file into the project directory, above the genome directory.

        :param self: BSF Analysis
        :type self: Analysis
        :param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        :type prefix: str
        :return: Nothing
        :rtype: None
        """

        if prefix:
            file_name = '{}_hub.txt'.format(prefix)
        else:
            file_name = 'hub.txt'

        output = str()

        if prefix:
            output += 'hub {}_{}\n'.format(self.project_name, prefix)
            output += 'shortLabel {}_{}\n'.format(self.project_name, prefix)
            output += 'longLabel Project {}_{}\n'.format(self.project_name, prefix)
            output += 'genomesFile {}_genomes.txt\n'.format(prefix)
        else:
            output += 'hub {}\n'.format(self.project_name)
            output += 'shortLabel {}\n'.format(self.project_name)
            output += 'longLabel Project {}\n'.format(self.project_name)
            output += 'genomesFile genomes.txt\n'

        output += 'email {}\n'.format(self.e_mail)

        # The [prefix_]hub.txt goes into the project directory above the genome directory.
        file_path = os.path.join(self.project_directory, file_name)

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

    def ucsc_hub_write_genomes(self, prefix=None):

        """Write a UCSC Track Hub genomes.txt file into the project directory, above the genome directory.

        :param self: BSF Analysis
        :type self: Analysis
        :param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        :type prefix: str
        :return: Nothing
        :rtype: None
        """

        if prefix:
            file_name = '{}_genomes.txt'.format(prefix)
        else:
            file_name = 'genomes.txt'

        output = str()

        output += 'genome {}\n'.format(self.genome_version)
        if prefix:
            output += 'trackDb {}/{}_trackDB.txt\n'.format(self.genome_version, prefix)
        else:
            output += 'trackDb {}/trackDB.txt\n'.format(self.genome_version)

        # The [prefix_]genomes.txt goes into the project directory above the genome directory.
        file_path = os.path.join(self.project_directory, file_name)

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

    def ucsc_hub_write_tracks(self, output, prefix=None):

        """Write a UCSC Track Hub trackDB.txt file into the genome directory.

        :param self: BSF Analysis
        :type self: Analysis
        :param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        :type prefix: str
        :return: Nothing
        :rtype: None
        """

        if prefix:
            file_name = '{}_trackDB.txt'.format(prefix)
        else:
            file_name = 'trackDB.txt'

        # The [prefix_]trackDB.txt goes into the genome directory under the project directory.
        file_path = os.path.join(self.genome_directory, file_name)

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()


class Configuration(object):
    """BSF Configuration class.

    The BSF Configuration class represents a UNIX-style initialisation (*.ini) file and
    an associated Python SafeConfigParser object to parse the file.

    Attributes:
    :ivar config_file: Configuration file path
    :type config_file: str, unicode
    :ivar config_parser: Python SafeConfigParser
    :type config_parser: SafeConfigParser
    """

    @staticmethod
    def section_from_instance(instance):

        """Get a configuration section from a Python instance.

        :param instance: A Python instance (or object)
        :type instance: object
        :return: Configuration section string
        :rtype: str
        """

        match = re.search(pattern=r'^<([^ ]+)\s+', string=repr(instance))

        if match:
            return match.group(1)
        else:
            return str()

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new BSF Configuration object based on a configuration file path.

        :param cls: Class
        :type cls: Configuration
        :param config_file: Configuration file path.
        Both, user and variable expansion gets applied.
        :type config_file: str, unicode
        :return: BSF Configuration
        :rtype: Configuration
        """

        assert isinstance(config_file, basestring)

        config_file = os.path.expanduser(path=config_file)
        config_file = os.path.expandvars(path=config_file)

        # Since ConfigParser options are used as command line options,
        # they have to be case sensitive.
        # Hence, override optionxform() with str().

        config_parser = SafeConfigParser()
        config_parser.optionxform = str

        configuration = cls(config_file=config_file, config_parser=config_parser)

        files = configuration.config_parser.read(configuration.config_file)

        if len(files) == 0:
            message = 'Could not find configuration file {!r}.'.format(configuration.config_file)
            raise Exception(message)

        return configuration

    def __init__(self, config_file=None, config_parser=None):

        """Initialise a BSF Configuration object.

        :param self: BSF Configuration
        :type self: Configuration
        :param config_file: Configuration file path
        :type config_file: str, unicode
        :param config_parser: Python SafeConfigParser
        :type config_parser: SafeConfigParser
        :return: Nothing
        :rtype: None
        """

        if config_file:
            self.config_file = config_file
        else:
            self.config_file = str()

        if config_parser:
            self.config_parser = config_parser
        else:
            self.config_parser = SafeConfigParser()

    def trace(self, level):

        """Trace a BSF Configuration object.

        :param self: BSF Configuration
        :type self: Configuration
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  config_file:   {!r}\n'.format(indent, self.config_file)
        output += '{}  config_parser: {!r}\n'.format(indent, self.config_parser)

        return output

    def get_expanded_directory(self, config_section, config_option):

        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e. on UNIX ~ or ~user and
        any environment variables i.e. on UNIX ${NAME} or $NAME.
        :param self: BSF Configuration
        :type self: Configuration
        :param config_section: Configuration section string
        :type config_section: str
        :param config_option: Configuration option string
        :type config_option: str
        :return: Expanded directory
        :rtype: str
        """

        directory = self.config_parser.get(config_section, config_option)
        directory = os.path.expanduser(path=directory)
        directory = os.path.expandvars(path=directory)

        return directory


class Default(object):
    """The BSF Default class specifies the application or library default configuration.

    Attributes:
    :cvar global_default: Global BSF Default
    :type global_default: Default
    :ivar classpath_picard: Picard Java Archive (JAR) class path directory
    :type classpath_picard: str, unicode
    :ivar classpath_illumina2bam: Illumina2bam Java Archive (JAR) class path directory
    :type classpath_illumina2bam: str, unicode
    :ivar directory_home: Home directory for all data
    :type directory_home: str, unicode
    :ivar directory_runs_illumina: Sub-directory for Illumina runs
    :type directory_runs_illumina: str, unicode
    :ivar directory_sequences: Sub-directory for sequences
    :type directory_sequences: str, unicode
    :ivar directory_samples: Sub-directory for processed samples
    :type directory_samples: str, unicode
    :ivar directory_projects: Sub-directory for processed projects
    :type directory_projects: str, unicode
    :ivar directory_public_html: Sub-directory for the web server
    :type directory_public_html: str, unicode
    :ivar directory_genomes: Directory for genomes and their annotation
    :type directory_genomes: str, unicode
    :ivar directory_annotations: Sub-directory for genome annotations
    :type directory_annotations: str, unicode
    :ivar indices: Python dict of program name key and index directory name value data
    :type indices: dict
    :ivar drms_implementation: DRMS implementation (e.g. Bash, SGE)
    :type drms_implementation: str
    :ivar drms_maximum_threads: DRMS maximum threads
    :type drms_maximum_threads: str
    :ivar drms_memory_limit_hard: DRMS memory limit hard
    :type drms_memory_limit_hard: str
    :ivar drms_memory_limit_soft: DRMS memory limit soft
    :type drms_memory_limit_soft: str
    :ivar operator_e_mail: Operator e-mail
    :type operator_e_mail: str
    :ivar ucsc_host_name: UCSC Genome Browser host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...)
    :type ucsc_host_name: str
    :ivar url_protocol: URL protocol (i.e. HTTP)
    :type url_protocol: str
    :ivar url_host_name: URL host name
    :type url_host_name:str
    :ivar url_relative_projects: Sub-directory for analysis projects
    :type url_relative_projects: str
    :ivar url_relative_chip_seq: Sub-directory for ChIPSeq experiments
    :type url_relative_chip_seq:str
    :ivar url_relative_dna_seq: Sub-directory for general DNA sequencing
    :type url_relative_dna_seq: str
    :ivar url_relative_rna_seq: Sub-directory for RNA-Seq experiments
    :type url_relative_rna_seq: str
    """

    global_default = None

    @staticmethod
    def get_global_default():

        """Get the global BSF Default configuration and initialise it, if not already done so.

        :return: BSF Default
        :rtype: Default
        """

        if not Default.global_default:
            Default.global_default = Default.from_global_file()

        return Default.global_default

    @classmethod
    def from_global_file(cls):

        """Create a new BSF Default object from the global default configuration file.

        The default configuration is based on the file $HOME/.bsfpython.ini in the user's home directory.
        :param cls: Class
        :type cls: Default
        :return: BSF Default
        :rtype: Default
        """

        config_file = '~/.bsfpython.ini'
        config_file = os.path.expanduser(path=config_file)

        return cls.from_config_file(config_file=config_file)

    @classmethod
    def from_config_file(cls, config_file):

        """Create a new BSF Default object from a UNIX-style configuration file.

        :param cls: Class
        :type cls: Default
        :param config_file: UNIX-style configuration file
        :type config_file: str, unicode
        :return: BSF Default
        :rtype: Default
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):

        """Create a new BSF Default objects from a BSF Configuration object.

        :param cls: Class
        :type cls: Default
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :return: BSF Default
        :rtype: Default
        """

        assert isinstance(configuration, Configuration)

        default = cls()

        default.set_Configuration(configuration=configuration)

        return default

    def __init__(self, classpath_picard=None, classpath_illumina2bam=None, directory_home=None,
                 directory_runs_illumina=None, directory_sequences=None, directory_samples=None,
                 directory_projects=None, directory_public_html=None, directory_genomes=None,
                 directory_annotations=None, indices=None, drms_implementation=None, drms_maximum_threads=None,
                 drms_memory_limit_hard=None, drms_memory_limit_soft=None, operator_e_mail=None,
                 operator_sequencing_centre=None, ucsc_host_name=None, url_protocol=None, url_host_name=None,
                 url_relative_projects=None, url_relative_chip_seq=None, url_relative_dna_seq=None,
                 url_relative_rna_seq=None):

        """Initialise a BSF Default object.

        :param self: BSF Default
        :type self: Default
        :param classpath_picard: Picard Java Archive (JAR) class path directory
        :type classpath_picard: str, unicode
        :param classpath_illumina2bam: Illumina2bam Java Archive (JAR) class path directory
        :type classpath_illumina2bam: str, unicode
        :param directory_home: Home directory for all data
        :type directory_home: str, unicode
        :param directory_runs_illumina: Sub-directory for Illumina runs
        :type directory_runs_illumina: str, unicode
        :param directory_sequences: Sub-directory for sequences
        :type directory_sequences: str, unicode
        :param directory_samples: Sub-directory for processed samples
        :type directory_samples: str, unicode
        :param directory_projects: Sub-directory for processed projects
        :type directory_projects: str, unicode
        :param directory_public_html: Sub-directory for the web server
        :type directory_public_html: str, unicode
        :param directory_genomes: Directory for genomes and their annotation
        :type directory_genomes: str, unicode
        :param directory_annotations: Sub-directory for genome annotations
        :type directory_annotations: str, unicode
        :param indices: Python dict of program name key and index directory name value data
        :type indices: dict
        :param drms_implementation: DRMS implementation (e.g. Bash, SGE)
        :type drms_implementation: str
        :param drms_maximum_threads: DRMS maximum threads
        :type drms_maximum_threads: str
        :param drms_memory_limit_hard: DRMS memory limit hard
        :type drms_memory_limit_hard: str
        :param drms_memory_limit_soft: DRMS memory limit soft
        :type drms_memory_limit_soft: str
        :param operator_e_mail: Operator e-mail
        :type operator_e_mail: str
        :param operator_sequencing_centre: BAM sequencing centre code
        :type operator_sequencing_centre: str
        :param ucsc_host_name: UCSC Genome Browser host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...)
        :type ucsc_host_name: str
        :param url_protocol: URL protocol (i.e. HTTP)
        :type url_protocol: str
        :param url_host_name: URL host name
        :type url_host_name:str
        :param url_relative_projects: Sub-directory for analysis projects
        :type url_relative_projects: str
        :param url_relative_chip_seq: Sub-directory for ChIPSeq experiments
        :type url_relative_chip_seq:str
        :param url_relative_dna_seq: Sub-directory for general DNA sequencing
        :type url_relative_dna_seq: str
        :param url_relative_rna_seq: Sub-directory for RNA-Seq experiments
        :type url_relative_rna_seq: str
        :return: Nothing
        :rtype: None
        """

        # Set Java class path information.

        if classpath_illumina2bam:
            self.classpath_illumina2bam = classpath_illumina2bam
        else:
            self.classpath_illumina2bam = str()

        if classpath_picard:
            self.classpath_picard = classpath_picard
        else:
            self.classpath_picard = str()

        # Set directory information.

        if directory_home:
            self.directory_home = directory_home
        else:
            self.directory_home = str()

        if directory_runs_illumina:
            self.directory_runs_illumina = directory_runs_illumina
        else:
            self.directory_runs_illumina = str()

        if directory_sequences:
            self.directory_sequences = directory_sequences
        else:
            self.directory_sequences = str()

        if directory_samples:
            self.directory_samples = directory_samples
        else:
            self.directory_samples = str()

        if directory_projects:
            self.directory_projects = directory_projects
        else:
            self.directory_projects = str()

        if directory_public_html:
            self.directory_public_html = directory_public_html
        else:
            self.directory_public_html = str()

        if directory_genomes:
            self.directory_genomes = directory_genomes
        else:
            self.directory_genomes = str()

        if directory_annotations:
            self.directory_annotations = directory_annotations
        else:
            self.directory_annotations = str()

        # Set index information.

        if indices:
            self.indices = indices
        else:
            self.indices = dict()

        # Set DRMS information.

        if drms_implementation:
            self.drms_implementation = drms_implementation
        else:
            self.drms_implementation = str()

        if drms_maximum_threads:
            self.drms_maximum_threads = drms_maximum_threads
        else:
            self.drms_maximum_threads = str()

        if drms_memory_limit_hard:
            self.drms_memory_limit_hard = drms_memory_limit_hard
        else:
            self.drms_memory_limit_hard = str()

        if drms_memory_limit_soft:
            self.drms_memory_limit_soft = drms_memory_limit_soft
        else:
            self.drms_memory_limit_soft = str()

        # Set operator information.

        if operator_e_mail:
            self.operator_e_mail = operator_e_mail
        else:
            self.operator_e_mail = str()

        if operator_sequencing_centre:
            self.operator_sequencing_centre = operator_sequencing_centre
        else:
            self.operator_sequencing_centre = str()

        # Set UCSC Genome Browser information.

        if ucsc_host_name:
            self.ucsc_host_name = ucsc_host_name
        else:
            self.ucsc_host_name = str()

        # Set URL information.

        if url_protocol:
            self.url_protocol = url_protocol
        else:
            self.url_protocol = str()

        if url_host_name:
            self.url_host_name = url_host_name
        else:
            self.url_host_name = str()

        if url_relative_projects:
            self.url_relative_projects = url_relative_projects
        else:
            self.url_relative_projects = str()

        if url_relative_chip_seq:
            self.url_relative_chip_seq = url_relative_chip_seq
        else:
            self.url_relative_chip_seq = str()

        if url_relative_dna_seq:
            self.url_relative_dna_seq = url_relative_dna_seq
        else:
            self.url_relative_dna_seq = str()

        if url_relative_rna_seq:
            self.url_relative_rna_seq = url_relative_rna_seq
        else:
            self.url_relative_rna_seq = str()

    def set_Configuration(self, configuration):

        """Set instance variables of a BSF Default object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF Default
        :type self: Default
        :param configuration: BSF Configuration
        :type configuration: Configuration
        """

        assert isinstance(configuration, Configuration)

        cp = configuration.config_parser

        # Reading configuration cannot be done via a single Python dict,
        # because each option really needs defining.

        section = 'classpath'

        self.classpath_illumina2bam = cp.get(section=section, option='illumina2bam')
        self.classpath_picard = cp.get(section=section, option='picard')

        section = 'directories'

        self.directory_home = cp.get(section=section, option='home')
        self.directory_runs_illumina = cp.get(section=section, option='runs_illumina')
        self.directory_sequences = cp.get(section=section, option='sequences')
        self.directory_samples = cp.get(section=section, option='samples')
        self.directory_projects = cp.get(section=section, option='projects')
        self.directory_public_html = cp.get(section=section, option='public_html')
        self.directory_genomes = cp.get(section=section, option='genomes')
        self.directory_annotations = cp.get(section=section, option='annotations')

        section = 'indices'

        for option in cp.options(section=section):
            self.indices[option] = cp.get(section=section, option=option)

        section = 'drms'

        self.drms_implementation = cp.get(section=section, option='implementation')
        self.drms_maximum_threads = cp.get(section=section, option='maximum_threads')
        self.drms_memory_limit_hard = cp.get(section=section, option='memory_limit_hard')
        self.drms_memory_limit_soft = cp.get(section=section, option='memory_limit_soft')
        # self.drms_memory_free_mem = cp.get(section=section, option='memory_free_mem')
        # self.drms_memory_free_swap = cp.get(section=section, option='memory_free_swap')
        # self.drms_memory_free_virtual = cp.get(section=section, option='memory_free_virtual')
        self.drms_parallel_environment = cp.get(section=section, option='parallel_environment')
        self.drms_queue = cp.get(section=section, option='queue')

        section = 'operator'

        self.operator_e_mail = cp.get(section=section, option='e_mail')
        self.operator_sequencing_centre = cp.get(section=section, option='sequencing_centre')

        section = 'ucsc'

        self.ucsc_host_name = cp.get(section=section, option='ucsc_host_name')

        section = 'url'

        self.url_protocol = cp.get(section=section, option='protocol')
        self.url_host_name = cp.get(section=section, option='host_name')
        self.url_relative_projects = cp.get(section=section, option='relative_projects')
        # TODO: With the new public_html/projects directory, these should become obsolete.
        self.url_relative_chip_seq = cp.get(section=section, option='relative_chip_seq')
        self.url_relative_dna_seq = cp.get(section=section, option='relative_dna_seq')
        self.url_relative_rna_seq = cp.get(section=section, option='relative_rna_seq')

    @staticmethod
    def absolute_home():

        """
        Get the absolute directory path for the home directory.

        :return: Absolute path to the home directory
        :rtype; str, unicode
        """

        default = Default.get_global_default()

        return default.directory_home

    @staticmethod
    def absolute_runs_illumina():

        """Get the absolute directory path for Illumina runs.

        :return: Absolute path to the Illumina runs directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_runs_illumina):
            return default.directory_runs_illumina
        else:
            return os.path.join(default.directory_home, default.directory_runs_illumina)

    @staticmethod
    def absolute_sequences():

        """Get the absolute directory path for processed lanes.

        :return: Absolute path to the processed lanes directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_sequences):
            return default.directory_sequences
        else:
            return os.path.join(default.directory_home, default.directory_sequences)

    @staticmethod
    def absolute_projects():

        """Get the absolute directory path for projects.

        :return: Absolute path to the projects directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_projects):
            return default.directory_projects
        else:
            return os.path.join(default.directory_home, default.directory_projects)

    @staticmethod
    def absolute_samples():

        """Get the absolute directory path for processed samples.

        :return: Absolute path to the processed samples directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_samples):
            return default.directory_samples
        else:
            return os.path.join(default.directory_home, default.directory_samples)

    @staticmethod
    def absolute_public_html():

        """Get the absolute directory path for public HTML documents.

        :return: Absolute path to the public HTML directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_public_html):
            return default.directory_public_html
        else:
            return os.path.join(default.directory_home, default.directory_public_html)

    @staticmethod
    def absolute_genomes(genome_version):

        """Get the absolute directory path for genomes.

        :param genome_version: The genome version (e.g. mm10, ...)
        :type genome_version: str
        :return: Absolute path to the genomes directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if genome_version:
            return os.path.join(default.directory_genomes, genome_version)
        else:
            return default.directory_genomes

    @staticmethod
    def absolute_genome_annotation(genome_version):

        """Get the absolute directory path for genome annotation.

        :param genome_version: The genome version (e.g. mm10, ...)
        :type genome_version: str
        :return: Absolute path to the genome annotation directory
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        return os.path.join(default.directory_genomes, genome_version, default.directory_annotations)

    @staticmethod
    def absolute_genome_fasta(genome_version, genome_index):

        """Get the absolute file path to a genome in FASTA format.

        :param genome_version: Genome version (e.g. mm10, ...)
        :type genome_version: str
        :param genome_index: Genome index (e.g. bowtie2, ...)
        :type genome_index: str
        :return: Absolute path to the genome FASTA file
        :rtype: str, unicode
        """

        default = Default.get_global_default()

        if not genome_index in default.indices:
            message = 'Unknown genome index name {}'.format(genome_index)
            raise Exception(message)

        return os.path.join(Default.absolute_genomes(genome_version),
                            default.indices[genome_index], genome_version + '.fa')

    @staticmethod
    def url_absolute_base():

        """Return the absolute URL to the web site.

        :return: URL string
        :rtype: str
        """

        default = Default.get_global_default()

        return '{}://{}'.format(default.url_protocol, default.url_host_name)

    @staticmethod
    def url_absolute_projects():

        """Return the absolute URL to the analysis projects directory.

        :return: URL string
        :rtype: str
        """

        default = Default.get_global_default()

        return '{}/{}'.format(default.url_absolute_base(), default.url_relative_projects)

    @staticmethod
    def url_absolute_chip_seq():

        """Return the absolute URL to ChIP-Seq experiments.

        :return: URL string
        :rtype: str
        """

        # TODO: With the new public_html/projects directory, this should become obsolete.

        default = Default.get_global_default()

        return '{}/{}'.format(default.url_absolute_base(), default.url_relative_chip_seq)

    @staticmethod
    def url_absolute_dna_seq():

        """Return the absolute URL to DNA sequencing experiments.

        :return: URL string
        :rtype: str
        """

        # TODO: With the new public_html/projects directory, this should become obsolete.

        default = Default.get_global_default()

        return '{}/{}'.format(default.url_absolute_base(), default.url_relative_dna_seq)

    @staticmethod
    def url_absolute_rna_seq():

        """Return the absolute URL to RNA-Seq experiments.

        :return: URL string
        :rtype: str
        """

        # TODO: With the new public_html/projects directory, this should become obsolete.

        default = Default.get_global_default()

        return '{}/{}'.format(default.url_absolute_base(), default.url_relative_rna_seq)


class DRMS(object):
    """BSF Distributed Resource Management System (DRMS) class.

    The BSF DRMS class represents a Distributed Resource Management System or
    batch job scheduler.

    Attributes:
    :ivar name: Name
    :type name: str
    :ivar work_directory: Work directory path
    :type work_directory: str
    :ivar implementation: Implementation (e.g. SGE, ...)
    :type implementation: str
    :ivar memory_free_mem: Memory limit (free)
    :type memory_free_mem: str
    :ivar memory_limit_hard: Memory limit (hard)
    :type memory_limit_hard: str
    :ivar memory_limit_soft: Memory limit (soft)
    :type memory_limit_soft: str
    :ivar parallel_environment: Parallel environment
    :type parallel_environment: str
    :ivar queue: Queue
    :type queue: str
    :ivar threads: Number of threads
    :type threads: int
    :ivar hold: Hold on job scheduling
    :type hold: str
    :ivar is_script: BSF Executable objects represent shell scripts,
    or alternatively binary programs
    :type is_script: bool
    :ivar executables: Python list of BSF Executable objects
    :type executables: list
    """

    @classmethod
    def from_Analysis(cls, name, work_directory, analysis):

        """Create a BSF DRMS object from a BSF Analysis object.

        :param cls: Class
        :type cls: DRMS
        :param name: Name
        :type name: str
        :param work_directory: Work directory
        :type work_directory: str
        :param analysis: BSF Analysis
        :type analysis: Analysis
        :return: BSF DRMS object
        :rtype: DRMS
        """

        assert isinstance(analysis, Analysis)

        # Set a minimal set of global defaults.

        default = Default.get_global_default()

        drms = cls(name=name, work_directory=work_directory,
                   implementation=default.drms_implementation,
                   memory_limit_hard=default.drms_memory_limit_hard,
                   memory_limit_soft=default.drms_memory_limit_soft,
                   parallel_environment=default.drms_parallel_environment,
                   queue=default.drms_queue)

        # A "Bio.BSF.DRMS" section specifies defaults for all BSF DRMS objects of a BSF Analysis.

        section = '{}.{}'.format(__name__, cls.__name__)
        drms.set_Configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 0:
            print 'DRMS configuration section: {}'.format(section)

        # A "Bio.BSF.Analysis.*.DRMS" pseudo-class section specifies
        # BSF Analysis-specific options for the BSF DRMS.

        section = '{}.DRMS'.format(analysis.configuration.section_from_instance(analysis))
        drms.set_Configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 0:
            print 'DRMS configuration section: {}'.format(section)

        # A "Bio.BSF.Analysis.*.DRMS.name" section specifies defaults
        # for a particular BSF DRMS objects of a BSF Analysis.

        section = '{}.DRMS.{}'.format(Configuration.section_from_instance(analysis), drms.name)
        drms.set_Configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 0:
            print 'DRMS configuration section: {}'.format(section)

        return drms

    @classmethod
    def from_Configuration(cls, name, work_directory, configuration, section):

        """Create a BSF DRMS object from a BSF Configuration object.

        :param cls: Class
        :type cls: DRMS
        :param name: Name
        :type name: str
        :param work_directory: Work directory
        :type work_directory: str
        :param configuration: BSF Configuration object
        :type configuration: Configuration
        :param section: Configuration section string
        :type section: str
        :return: BSF DRMS object
        :rtype: DRMS
        """

        assert isinstance(configuration, Configuration)

        # Set a minimal set of global defaults.

        default = Default.get_global_default()

        drms = cls(name=name, work_directory=work_directory,
                   implementation=default.drms_implementation,
                   memory_limit_hard=default.drms_memory_limit_hard,
                   memory_limit_soft=default.drms_memory_limit_soft,
                   parallel_environment=default.drms_parallel_environment,
                   queue=default.drms_queue)

        drms.set_Configuration(configuration=configuration, section=section)

        return drms

    def __init__(self, name, work_directory,
                 implementation=None,
                 memory_free_mem=None,
                 memory_free_swap=None,
                 memory_free_virtual=None,
                 memory_limit_hard=None,
                 memory_limit_soft=None,
                 parallel_environment=None,
                 queue=None,
                 threads=1,
                 hold=None,
                 is_script=False,
                 executables=None):

        """Initialise a BSF DRMS object.

        :param memory_limit_soft:
        :param self: BSF DRMS
        :type self: DRMS
        :param name: Name
        :type name: str
        :param work_directory: Work directory
        :type work_directory: str
        :param implementation: Implementation (e.g. SGE, ...)
        :type implementation: str
        :param memory_free_mem: Memory limit (free physical)
        :type memory_free_mem: str
        :param memory_free_swap: Memory limit (free swap)
        :type memory_free_swap: str
        :param memory_free_virtual: Memory limit (free virtual)
        :type memory_free_virtual: str
        :param memory_limit_hard: Memory limit (hard)
        :type memory_limit_hard: str
        :param memory_limit_soft: Memory limit (soft)
        :type memory_limit_soft: str
        :param parallel_environment: Parallel environment
        :type parallel_environment: str
        :param queue: Queue
        :type queue: str
        :param threads: Number of threads
        :type threads: int
        :param hold: Hold on job scheduling
        :type hold: str
        :param is_script: BSF Executable objects represent shell scripts,
        or alternatively binary programs
        :type is_script: bool
        :param executables: Python list of BSF Executable objects
        :type executables: list
        :return: Nothing
        :rtype: None
        """

        if name:
            self.name = name
        else:
            self.name = str()

        if work_directory:
            self.work_directory = work_directory
        else:
            self.work_directory = str()

        if implementation:
            self.implementation = implementation
        else:
            self.implementation = str()

        if memory_free_mem:
            self.memory_free_mem = memory_free_mem
        else:
            self.memory_free_mem = str()

        if memory_free_swap:
            self.memory_free_swap = memory_free_swap
        else:
            self.memory_free_swap = str()

        if memory_free_virtual:
            self.memory_free_virtual = memory_free_virtual
        else:
            self.memory_free_virtual = str()

        if memory_limit_hard:
            self.memory_limit_hard = memory_limit_hard
        else:
            self.memory_limit_hard = str()

        if memory_limit_soft:
            self.memory_limit_soft = memory_limit_soft
        else:
            self.memory_limit_soft = str()

        if parallel_environment:
            self.parallel_environment = parallel_environment
        else:
            self.parallel_environment = str()

        if queue:
            self.queue = queue
        else:
            self.queue = str()

        self.threads = threads

        if hold:
            self.hold = hold
        else:
            self.hold = str()

        self.is_script = is_script

        if executables:
            self.executables = executables
        else:
            self.executables = list()

    def trace(self, level):

        """Trace a BSF DRMS object.

        :param self: BSF DRMS
        :type self: DRMS
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  name:                 {!r}\n'. \
            format(indent, self.name)
        output += '{}  work_directory:       {!r}\n'. \
            format(indent, self.work_directory)
        output += '{}  implementation:       {!r}\n'. \
            format(indent, self.implementation)
        output += '{}  memory_free_mem:      {!r}\n'. \
            format(indent, self.memory_free_mem)
        output += '{}  memory_free_swap:     {!r}\n'. \
            format(indent, self.memory_free_swap)
        output += '{}  memory_free_virtual:  {!r}\n'. \
            format(indent, self.memory_free_virtual)
        output += '{}  memory_limit_hard:    {!r}\n'. \
            format(indent, self.memory_limit_hard)
        output += '{}  memory_limit_soft:    {!r}\n'. \
            format(indent, self.memory_limit_soft)
        output += '{}  queue:                {!r}\n'. \
            format(indent, self.queue)
        output += '{}  parallel_environment: {!r}\n'. \
            format(indent, self.parallel_environment)
        output += '{}  threads:              {!r}\n'. \
            format(indent, self.threads)
        output += '{}  hold:                 {!r}\n'. \
            format(indent, self.hold)
        output += '{}  is_script:            {!r}\n'. \
            format(indent, self.is_script)

        output += '{}  executables:\n'.format(indent)

        for executable in self.executables:
            output += executable.trace(level=level + 2)

        return output

    def set_Configuration(self, configuration, section):

        """Set instance variables of a BSF DRMS object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF DRMS
        :type self: DRMS
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section
        :type section: str
        """

        assert isinstance(configuration, Configuration)

        if not configuration.config_parser.has_section(section=section):
            message = 'Section {!r} not defined in BSF Configuration file {!r}.'. \
                format(section, configuration.config_file)
            warnings.warn(message, UserWarning)

            return

        # The configuration section is available.

        if configuration.config_parser.has_option(section=section, option='hold'):
            self.hold = configuration.config_parser.get(section=section,
                                                        option='hold')

        if configuration.config_parser.has_option(section=section, option='implementation'):
            self.implementation = configuration.config_parser.get(section=section,
                                                                  option='implementation')

        if configuration.config_parser.has_option(section=section, option='is_script'):
            self.is_script = configuration.config_parser.getboolean(section=section,
                                                                    option='is_script')

        if configuration.config_parser.has_option(section=section, option='memory_free_mem'):
            self.memory_free_mem = configuration.config_parser.get(section=section,
                                                                   option='memory_free_mem')

        if configuration.config_parser.has_option(section=section, option='memory_free_swap'):
            self.memory_free_swap = configuration.config_parser.get(section=section,
                                                                    option='memory_free_swap')

        if configuration.config_parser.has_option(section=section, option='memory_free_virtual'):
            self.memory_free_virtual = configuration.config_parser.get(section=section,
                                                                       option='memory_free_virtual')

        if configuration.config_parser.has_option(section=section, option='memory_hard'):
            self.memory_limit_hard = configuration.config_parser.get(section=section,
                                                                     option='memory_hard')

        if configuration.config_parser.has_option(section=section, option='memory_soft'):
            self.memory_limit_soft = configuration.config_parser.get(section=section,
                                                                     option='memory_soft')

        if configuration.config_parser.has_option(section=section, option='parallel_environment'):
            self.parallel_environment = configuration.config_parser.get(section=section,
                                                                        option='parallel_environment')

        if configuration.config_parser.has_option(section=section, option='queue'):
            self.queue = configuration.config_parser.get(section=section,
                                                         option='queue')

        if configuration.config_parser.has_option(section=section, option='threads'):
            self.threads = configuration.config_parser.get(section=section,
                                                           option='threads')

    def add_Executable(self, executable):

        """Add a BSF Executable object.

        :param self: BSF DRMS
        :type self: DRMS
        :param executable: BSF Executable
        :type executable: Executable
        :return: Nothing
        :rtype: None
        """

        assert isinstance(executable, Executable)

        self.executables.append(executable)

    def submit(self, debug=0):

        """Submit a command line for each BSF Executable object.

        :param self: BSF DRMS
        :type self: DRMS
        :param debug: Debug level
        :type debug: int
        :return: Nothing
        :rtype: None
        """

        # Dynamically import the module specific for the configured DRMS implementation.

        module = importlib.import_module(__name__ + '.DRMS.' + self.implementation)

        module.submit(self, debug=debug)


class Command(object):
    """BSF Command class.

    The BSF Command class represents one (subordinate) command,
    its options and arguments and possibly another subordinate command.

    Attributes:
    :ivar command: Command or program
    :type command: str
    :ivar options: Python dict of option keys and values
    :type options: dict
    :ivar arguments: Python list of arguments
    :type arguments: list
    :ivar sub_command: Subordinate BSF Command
    :type sub_command: Command
    """

    def __init__(self, command, options=None, arguments=None, sub_command=None):

        """Initialise a BSF Command object.

        :param command: Command
        :type command: str
        :param options: Python dict of program option and value pairs
        :type options: dict
        :param arguments: Python list of program arguments
        :type arguments: list
        :param sub_command: Subordinate BSF Command object
        :type sub_command: Command
        :return: Nothing
        :rtype: None
        """

        self.command = command

        if options:
            self.options = options
        else:
            self.options = dict()

        if arguments:
            self.arguments = arguments
        else:
            self.arguments = list()

        self.sub_command = sub_command

    def trace(self, level):

        """Trace a BSF Command object.

        :param self: BSF Command
        :type self: Command
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  command:            {!r}\n'. \
            format(indent, self.command)

        # List all options

        output += '{}  options:\n'.format(indent)

        for key in self.options.keys():
            output += '{}    key: {!r} value: {!r}\n'.format(indent, key, self.options[key])
            output += self.options[key].trace(level=level + 2)

        # List all arguments

        output += '{}  arguments:\n'.format(indent)

        i = 0
        for argument in self.arguments:
            output += '{}    {:2d}: {!r}\n'.format(indent, i, argument)
            i += 1

        if self.sub_command:
            output += self.sub_command.trace(level=level + 1)

        return output

    def add_Argument(self, argument, override):

        """Add a Bio.BSF.Argument.Switch, Bio.BSF.Argument.Option or one of its sub-classes.

        The sub-classes are
        Bio.BSF.Argument.SwitchLong,
        Bio.BSF.Argument.SwitchShort,
        Bio.BSF.Argument.Option,
        Bio.BSF.Argument.OptionLong,
        Bio.BSF.Argument.OptionShort or
        Bio.BSF.Argument.OptionPair.
        :param self: Bio.BSF.Command
        :type self: Command
        :param argument: Bio.BSF.Argument or sub-class thereof
        (Bio.BSF.Argument.Option or Bio.BSF.Argument.Switch)
        :type argument: Argument
        :param override: Override existing switch or option without warning.
        :type override: bool
        :return: Nothing
        :rtype: None
        """

        assert isinstance(argument, Argument)

        if not override and argument.key in self.options:
            message = 'Overwriting a Bio.BSF.Argument.Switch or Bio.BSF.Argument.Option ' \
                      'with key {!r} that exits already in Command {!r}.'. \
                format(argument.key, self.command)
            warnings.warn(message, UserWarning)

        self.options[argument.key] = argument

    def add_SwitchLong(self, key, override=False):

        """Initialise a Bio.BSF.Argument.SwitchLong object and add it to this Bio.BSF.Command object.

        :param self: BSF Command
        :type self: Command
        :param key: Key
        :type key: str
        :return: Nothing
        :rtype: None
        """

        self.add_Argument(argument=SwitchLong(key=key), override=override)

    def add_SwitchShort(self, key, override=False):

        """Initialise a Bio.BSF.Argument.SwitchShort object and add it to this Bio.BSF.Command object.

        :param self: BSF Command
        :type self: Command
        :param key: Key
        :type key: str
        :return: Nothing
        :rtype: None
        """

        self.add_Argument(argument=SwitchShort(key=key), override=override)

    def add_OptionLong(self, key, value, override=False):

        """Initialise a Bio.BSF.Argument.OptionLong object and add it to this Bio.BSF.Command object.

        :param self: BSF Command
        :type self: Command
        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Nothing
        :rtype: None
        """

        self.add_Argument(argument=OptionLong(key=key, value=value), override=override)

    def add_OptionShort(self, key, value, override=False):

        """Initialise a Bio.BSF.Argument.OptionShort object and add it to this Bio.BSF.Command object.

        :param self: BSF Command
        :type self: Command
        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Nothing
        :rtype: None
        """

        self.add_Argument(argument=OptionShort(key=key, value=value), override=override)

    def add_OptionPair(self, key, value, override=False):

        """Initialise a Bio.BSF.Argument.OptionPair object and add it to this Bio.BSF.Command object.

        :param self: BSF Command
        :type self: Command
        :param key: Key
        :type key: str
        :param value: Value
        :type value: str
        :return: Nothing
        :rtype: None
        """

        self.add_Argument(argument=OptionPair(key=key, value=value), override=override)

    def set_Configuration(self, configuration, section):

        """Set instance variables of a BSF Command object via a section of a BSF Configuration object.

        Instance variables without a configuration option remain unchanged.
        :param self: BSF Command
        :type self: Command
        :param configuration: BSF Configuration
        :type configuration: Configuration
        :param section: Configuration file section, defaults to instance class
        :type section: str
        """

        assert isinstance(configuration, Configuration)

        if not configuration.config_parser.has_section(section=section):
            message = 'Section {!r} not defined in BSF Configuration file {!r}.'. \
                format(section, configuration.config_file)
            warnings.warn(message, UserWarning)

            return

        # The configuration section is available.

        for option in configuration.config_parser.options(section=section):
            argument = Argument.from_key_value(key=option,
                                               value=configuration.config_parser.get(section=section,
                                                                                     option=option))

            self.add_Argument(argument=argument, override=False)

    def command_list(self):

        """Assemble the command line from program, options and arguments.

        :param self: BSF Command
        :type self: Command
        :return: Python list of program, options, switches and arguments
        :rtype: list
        """

        command_line = list()

        if self.command:
            command_line.append(self.command)

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            argument = self.options[key]

            if isinstance(argument, SwitchLong):
                command_line.append('--{}'.format(argument.key))
            elif isinstance(argument, SwitchShort):
                command_line.append('-{}'.format(argument.key))
            elif isinstance(argument, OptionLong):
                command_line.append('--{}'.format(argument.key))
                if argument.value:
                    command_line.append(str(argument.value))
            elif isinstance(argument, OptionShort):
                command_line.append('-{}'.format(argument.key))
                if argument.value:
                    command_line.append(str(argument.value))
            elif isinstance(argument, OptionPair):
                command_line.append('{}={}'.format(argument.key, argument.value))
            else:
                message = 'Unexpected object {!r} in Bio.BSF.Command.options dict'.format(argument)
                warnings.warn(message, UserWarning)

        # Add all arguments.

        for argument in self.arguments:
            command_line.append(str(argument))

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line.extend(self.sub_command.command_list())

        return command_line

    def command_str(self):

        """Assemble the command line from program, options, switches and arguments.

        :param self: BSF Command
        :type self: Command
        :return: A Python str of program, options, switches and arguments
        :rtype: str
        """

        command_line = str()

        if self.command:
            command_line += self.command

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:

            argument = self.options[key]

            if isinstance(argument, SwitchLong):
                command_line += ' --{}'.format(argument.key)
            elif isinstance(argument, SwitchShort):
                command_line += ' -{}'.format(argument.key)
            elif isinstance(argument, OptionLong):
                command_line += ' --{} {}'.format(argument.key, argument.value)
            elif isinstance(argument, OptionShort):
                command_line += ' -{} {}'.format(argument.key, argument.value)
            elif isinstance(argument, OptionPair):
                command_line += ' {}={}'.format(argument.key, argument.value)
            else:
                message = 'Unexpected object {!r} in Bio.BSF.Command.options dict'.format(argument)
                warnings.warn(message, UserWarning)

        # Add all arguments.

        for argument in self.arguments:
            command_line += ' {}'.format(argument)

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line += ' ' + self.sub_command.command_str()

        return command_line


class Executable(Command):
    """BSF Executable class.

    The BSF Executable class represents an executable program,
    its options and arguments.

    Attributes:
    :ivar name: Name in the context of a Bio.BSF.DRMS dependency
    :type name: str
    :ivar program: Program (executable or full file path)
    :type program: str
    :ivar options: Python dict of option keys and values
    :type options: dict
    :ivar arguments: Python list of arguments
    :type arguments: list
    :ivar sub_command: Subordinate BSF Command
    :type sub_command: Command
    :ivar stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
    :type stdout_path: str, unicode
    :ivar stderr_path: Standard error (STDERR) redirection in Bash (2>word)
    :type stderr_path: str, unicode
    :ivar dependencies: Python list of BSF Executable name strings in the
                        context of BSF DRMS dependencies
    :type dependencies: list
    :ivar hold: Hold on job scheduling
    :type hold: str
    :ivar process_identifier: Process identifier
    :type process_identifier: str
    :ivar process_name: Process name
    :type process_name: str
    """

    @classmethod
    def from_Analysis(cls, name, program, analysis):

        """Create a BSF Executable object from a BSF Analysis object.

        :param cls: Class
        :type cls: Executable
        :param name: Name
        :type name: str
        :param program: Program
        :type program: str
        :param analysis: BSF Analysis
        :type analysis: Analysis
        :return: BSF Executable object
        :rtype: Executable
        """

        assert isinstance(analysis, Analysis)

        # Initialise a BSF Executable object with default values.

        executable = cls(name=name, program=program)

        section = Configuration.section_from_instance(executable)

        # For plain Bio.BSF.Executable objects append the value of the
        # Bio.BSF.Executable.command to make this more meaningful.

        if section == 'Bio.BSF.Executable':
            section += '.{}'.format(executable.command)

        if analysis.debug > 0:
            print 'Executable configuration section: {}'.format(section)

        executable.set_Configuration(configuration=analysis.configuration, section=section)

        return executable

    @classmethod
    def from_Configuration(cls, name, program, configuration, section):

        """Create a BSF Executable object from a BSF Configuration object.

        :param cls: Class
        :type cls: Executable
        :param name: Name
        :type name: str
        :param program: Program
        :type program: str
        :param configuration: BSF Configuration object
        :type configuration: Configuration
        :param section: Configuration section string
        :type section: str
        :return: BSF Executable object
        :rtype: Executable
        """

        assert isinstance(configuration, Configuration)

        executable = cls(name=name, program=program)

        executable.set_Configuration(configuration=configuration, section=section)

        return executable

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 process_identifier=None, process_name=None):

        """Initialise a BSF Executable object.

        :param name:
        :param self: BSF Executable
        :type self: Executable
        :param program: Program
        :type program: str
        :param options: Python dict of program option and value pairs
        :type options: dict
        :param arguments: Python list of program arguments
        :type arguments: list
        :param sub_command: Subordinate BSF Command
        :type sub_command: Command
        :param stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
        :type stdout_path: str, unicode
        :param stderr_path: Standard error (STDERR) redirection in Bash (2>word)
        :type stderr_path: str, unicode
        :param dependencies: Python list of BSF Executable
         name strings in the context of BSF DRMS dependencies
        :type dependencies: list
        :param hold: Hold on job scheduling
        :type hold: str
        :param process_identifier: Process identifier
        :type process_identifier: str
        :param process_name: Process name
        :type process_name: str
        :return: Nothing
        :rtype: None
        """

        self.name = name

        super(Executable, self).__init__(command=program, options=options, arguments=arguments,
                                         sub_command=sub_command)

        if stderr_path:
            self.stderr_path = stderr_path
        else:
            self.stderr_path = str()

        if stdout_path:
            self.stdout_path = stdout_path
        else:
            self.stdout_path = str()

        if dependencies:
            self.dependencies = dependencies
        else:
            self.dependencies = list()

        if hold:
            self.hold = hold
        else:
            self.hold = str()

        if process_identifier:
            self.process_identifier = process_identifier
        else:
            self.process_identifier = str()

        if process_name:
            self.process_name = process_name
        else:
            self.process_name = str()

    def trace(self, level):

        """Trace a BSF Executable object.

        :param self: BSF Executable
        :type self: Executable
        :param level: Indentation level
        :type level: int
        :return: Trace information
        :rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  name:               {!r}\n'. \
            format(indent, self.name)
        output += '{}  stdout:             {!r}\n'. \
            format(indent, self.stdout_path)
        output += '{}  stderr:             {!r}\n'. \
            format(indent, self.stderr_path)
        output += '{}  hold:               {!r}\n'. \
            format(indent, self.hold)
        output += '{}  process_identifier: {!r}\n'. \
            format(indent, self.process_identifier)
        output += '{}  process_name:       {!r}\n'. \
            format(indent, self.process_name)

        # List all dependencies.

        output += '{}  dependencies:\n'.format(indent)

        i = 0
        for dependency in self.dependencies:
            output += '{}    {:2d} {!r}\n'.format(indent, i, dependency)
            i += 1

        # Trace the Command super-class.

        output += super(Executable, self).trace(level=level + 1)

        return output

    def command_list(self):

        """Assemble the command line from program, options and arguments.

        :param self: BSF Executable
        :type self: Executable
        :return: Python list of program, options and arguments
        :rtype: list
        """

        command = list()

        command.extend(super(Executable, self).command_list())

        # The stdout_path and stderr_path gets appended in specific modules.

        return command

    def command_str(self):

        """Assemble the command line from program, options, switches and arguments.

        :param self: BSF Executable
        :type self: Executable
        :return: A Python str of program, options, switches and arguments
        :rtype: str
        """

        command = str()

        command += super(Executable, self).command_str()

        # The stdout_path and stderr_path gets appended in specific modules.

        return command


class Runnable(object):
    """BSF Runnable class.

    The BSF Runnable class represents a sub process.

    Attributes:
    """

    @staticmethod
    def process_stream(file_type, file_handle, thread_lock, file_path=None, debug=0):
        """BSF Runnable function to process STDOUT or STDERR from the child process as a thread.

        :param file_type: File handle type STDOUT or STDERR
        :type file_type: str
        :param file_handle: The STDOUT or STDERR file handle
        :type file_handle: file
        :param thread_lock: A Python threading.Lock object
        :type thread_lock: threading.Lock
        :param file_path: STDOUT file path
        :type file_path: str, unicode
        :param debug: Debug level
        :type debug: int
        :return: Nothing
        :rtype: None
        """

        if file_type not in ('STDOUT', 'STDERR'):
            message = 'The file_type has to be either STDOUT or STDERR.'
            raise Exception(message)

        thread_lock.acquire(True)
        if debug > 0:
            print '[{}] Started BSF Runner {} processor in module {}.'. \
                format(datetime.datetime.now().isoformat(), file_type, __name__)
        output_file = None
        if file_path:
            output_file = open(file_path, 'w')
            if debug > 0:
                print "[{}] Opened {} file '{}'.". \
                    format(datetime.datetime.now().isoformat(), file_type, file_path)
        thread_lock.release()

        for line in file_handle:
            thread_lock.acquire(True)
            if output_file:
                output_file.write(line)
            else:
                print '[{}] {}: {}'.format(datetime.datetime.now().isoformat(), file_type, line.rstrip())
            thread_lock.release()

        thread_lock.acquire(True)
        if debug > 0:
            print '[{}] Received EOF on {} pipe.'.format(datetime.datetime.now().isoformat(), file_type)
        if output_file:
            output_file.close()
            if debug > 0:
                print "[{}] Closed {} file '{}'.". \
                    format(datetime.datetime.now().isoformat(), file_type, file_path)
        thread_lock.release()

    @staticmethod
    def process_stdout(stdout_handle, thread_lock, stdout_path=None, debug=0):
        """BSF Runnable function to process STDOUT from the child process as a thread.

        :param stdout_handle: The STDOUT file handle
        :type stdout_handle: file
        :param thread_lock: A Python threading.Lock object
        :type thread_lock: threading.Lock
        :param stdout_path: STDOUT file path
        :type stdout_path: str, unicode
        :param debug: Debug level
        :type debug: int
        :return: Nothing
        :rtype: None
        """

        return Runnable.process_stream(file_type='STDOUT', file_handle=stdout_handle,
                                       thread_lock=thread_lock, file_path=stdout_path,
                                       debug=debug)

    @staticmethod
    def process_stderr(stderr_handle, thread_lock, stderr_path=None, debug=0):
        """BSF Runnable function to process STDERR from the child process as a thread.

        :param stderr_handle: The STDERR file handle
        :type stderr_handle: file
        :param thread_lock: A Python threading.Lock object
        :type thread_lock: threading.Lock
        :param stderr_path: STDOUT file path
        :type stderr_path: str, unicode
        :param debug: Debug level
        :type debug: int
        :return: Nothing
        :rtype: None
        """

        return Runnable.process_stream(file_type='STDERR', file_handle=stderr_handle,
                                       thread_lock=thread_lock, file_path=stderr_path,
                                       debug=debug)

    @staticmethod
    def run(executable, max_loop_counter=3, max_thread_joins=10, thread_join_timeout=10, debug=1):
        """BSF Runnable function to run a BSF Executable object as Python subprocess.

        :param executable: BSF Executable
        :type executable: Executable
        :param max_loop_counter: Maximum number of retries
        :type max_loop_counter: int
        :param max_thread_joins: Maximum number of attempts to join the output threads
        :type max_thread_joins: int
        :param thread_join_timeout: Timeout for each attempt to join the output threads
        :type thread_join_timeout: int
        :param debug: Debug level
        :type debug: int
        :return: Return value of the child in the Python subprocess.
        Negative values indicate that the child received a signal.
        :rtype: int
        """

        on_posix = 'posix' in sys.builtin_module_names

        loop_counter = 0
        child_return_code = 0

        while loop_counter < max_loop_counter:

            child_process = Popen(args=executable.command_list(),
                                  bufsize=-1,
                                  stdin=PIPE,
                                  stdout=PIPE,
                                  stderr=PIPE,
                                  shell=False,
                                  close_fds=on_posix)

            # Two threads, thread_out and thread_err reading STDOUT and STDERR, respectively,
            # should make sure that buffers are not filling up.

            thread_lock = Lock()

            thread_out = Thread(target=Runnable.process_stdout,
                                kwargs={'stdout_handle': child_process.stdout,
                                        'thread_lock': thread_lock,
                                        'stdout_path': executable.stdout_path,
                                        'debug': debug})
            thread_out.daemon = True  # Thread dies with the program.
            thread_out.start()

            thread_err = Thread(target=Runnable.process_stderr,
                                kwargs={'stderr_handle': child_process.stderr,
                                        'thread_lock': thread_lock,
                                        'stderr_path': executable.stderr_path,
                                        'debug': debug})
            thread_err.daemon = True  # Thread dies with the program.
            thread_err.start()

            # Wait for the child process to finish.

            child_return_code = child_process.wait()

            thread_join_counter = 0

            while thread_out.is_alive() and thread_join_counter < max_thread_joins:
                thread_lock.acquire(True)
                if debug > 0:
                    print '[{}] Waiting for STDOUT processor to finish.'. \
                        format(datetime.datetime.now().isoformat())
                thread_lock.release()

                thread_out.join(timeout=thread_join_timeout)
                thread_join_counter += 1

            thread_join_counter = 0

            while thread_err.is_alive() and thread_join_counter < max_thread_joins:
                thread_lock.acquire(True)
                if debug > 0:
                    print '[{}] Waiting for STDERR processor to finish.'. \
                        format(datetime.datetime.now().isoformat())
                thread_lock.release()

                thread_err.join(timeout=thread_join_timeout)
                thread_join_counter += 1

            if child_return_code > 0:
                if debug > 0:
                    print '[{}] Child process {} failed with exit code {}'. \
                        format(datetime.datetime.now().isoformat(), executable.name, +child_return_code)
                loop_counter += 1
            elif child_return_code < 0:
                if debug > 0:
                    print '[{}] Child process {} received signal {}.'. \
                        format(datetime.datetime.now().isoformat(), executable.name, -child_return_code)
            else:
                if debug > 0:
                    print '[{}] Child process {} completed successfully {}.'. \
                        format(datetime.datetime.now().isoformat(), executable.name, +child_return_code)
                break

        else:
            if debug > 0:
                print "[{}] BSF Runnable '{}' exceeded the maximum re-run counter {}." \
                    .format(datetime.datetime.now().isoformat(), executable.name, max_loop_counter)

        return child_return_code

    @staticmethod
    def evaluate_return_code(executable, return_code):
        """Evaluate a return code from the run method.
        :param executable: BSF Executable
        :type executable: Executable
        :param return_code: Return code
        :type return_code: int
        :return: Nothing
        :rtype: None
        """

        if return_code > 0:
            print "[{}] Child process '{}' failed with return code {}". \
                format(datetime.datetime.now().isoformat(), executable.name, +return_code)
        elif return_code < 0:
            print "[{}] Child process '{}' received signal {}.". \
                format(datetime.datetime.now().isoformat(), executable.name, -return_code)
        else:
            print "[{}] Child process '{}' completed with return code {}.". \
                format(datetime.datetime.now().isoformat(), executable.name, +return_code)
