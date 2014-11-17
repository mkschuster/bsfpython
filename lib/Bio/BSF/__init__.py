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
from pickle import Pickler, Unpickler, HIGHEST_PROTOCOL
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
    """The Analysis class represents a high-level analysis that may use one or more
    programs (Executable objects).

    Attributes:
    @ivar configuration: Configuration
    @type configuration: Configuration
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
    @ivar drms_list: Python list of DRMS objects
    @type drms_list: list
    @ivar runnable_dict: Python dict of Python str (Runnable.name) key data and Runnable value data
    @type runnable_dict: dict
    @ivar collection: Collection
    @type collection: Collection
    @ivar comparisons: Python dict of comparisons
    @type comparisons: dict
    @ivar samples: Python list of Sample objects
    @type samples: list
    """

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new Analysis object from a UNIX-style configuration file via the Configuration class.

        @param config_file: UNIX-style configuration file
        @type config_file: str | unicode
        @return: Analysis
        @rtype: Analysis
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new Analysis object from a Configuration object.

        @param configuration: Configuration
        @type configuration: Configuration
        @return: Analysis
        @rtype: Analysis
        """

        assert isinstance(configuration, Configuration)

        # Set a minimal set of global defaults.

        default = Default.get_global_default()

        analysis = cls(configuration=configuration, e_mail=default.operator_e_mail)

        # A "Bio.BSF.Analysis.*" section specifies defaults for this Analysis.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        analysis.set_Configuration(analysis.configuration, section=section)

        return analysis

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 sas_file=None, sas_prefix=None, e_mail=None, debug=0, drms_list=None,
                 runnable_dict=None, collection=None, comparisons=None, samples=None):
        """Initialise an Analysis object.

        @param configuration: Configuration
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: Analysis-wide input directory
        @type input_directory: str
        @param output_directory: Analysis-wide output directory
        @type output_directory: str
        @param project_directory: Analysis-wide project directory,
            normally under the Analysis-wide output directory
        @type project_directory: str
        @param genome_directory: Analysis-wide genome directory,
            normally under the Analysis-wide project directory
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
        @param drms_list: Python list of DRMS objects
        @type drms_list: list
        @param runnable_dict: Python dict of Python str (Runnable.name) and Runnable value data
        @type runnable_dict: dict
        @param collection: Collection
        @type collection: Collection
        @param comparisons: Python dict of Analysis-specific objects
            (i.e. Python tuple for RNA-Seq and ChIPSeqComparison for ChIPSeq)
        @type comparisons: dict
        @param samples: Python list of Sample objects
        @type samples: list
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

        if runnable_dict:
            self.runnable_dict = runnable_dict
        else:
            self.runnable_dict = dict()

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
        """Trace an Analysis object.

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
        output += '{}  input_directory: {!r}\n'.format(indent, self.input_directory)
        output += '{}  output_directory: {!r}\n'.format(indent, self.output_directory)
        output += '{}  genome_directory: {!r}\n'.format(indent, self.genome_directory)
        output += '{}  sas_file: {!r}\n'.format(indent, self.sas_file)
        output += '{}  sas_prefix: {!r}\n'.format(indent, self.sas_prefix)
        output += '{}  e_mail: {!r}\n'.format(indent, self.e_mail)
        output += '{}  debug: {!r}\n'.format(indent, self.debug)
        output += '{}  drms_list: {!r}\n'.format(indent, self.drms_list)
        output += '{}  runnable_dict: {!r}\n'.format(indent, self.runnable_dict)
        output += '{}  collection: {!r}\n'.format(indent, self.collection)
        output += '{}  comparisons: {!r}\n'.format(indent, self.comparisons)
        output += '{}  samples: {!r}\n'.format(indent, self.samples)

        output += '{}  Python List of Sample objects:'.format(indent)
        for sample in self.samples:
            output += '{}    Sample name: {!r} file_path: {!r}'.format(indent, sample.name, sample.file_path)

        if self.collection:
            output += self.collection.trace(level + 1)

        return output

    def add_runnable(self, runnable):
        """Add a Runnable.

        @param runnable: Runnable
        @type runnable: Runnable
        @raise Exception: A Runnable.name already exists in the Analysis
        """

        assert isinstance(runnable, Runnable)

        if runnable.name in self.runnable_dict:
            raise Exception("A Runnable object with name {!r} already exists in Analysis {!r}".
                            format(runnable.name, self.project_name))
        else:
            self.runnable_dict[runnable.name] = runnable

    def add_Sample(self, sample):
        """Add a Sample object to the Python list of Sample objects if it does not already exist.

        The check is based on the Python 'in' comparison operator and in lack of a specific
        __cmp__ method, relies on object identity (i.e. address).
        @param sample: Sample
        @type sample: Sample
        """

        if not sample:
            return

        assert isinstance(sample, Sample) or isinstance(sample, SampleGroup)

        if sample not in self.samples:
            self.samples.append(sample)

    def set_Configuration(self, configuration, section):
        """Set instance variables of an Analysis object via a section of a Configuration object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: Configuration
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @raise Exception: The specified section does not exist
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            raise Exception(
                'Section {!r} not defined in Configuration file {!r}.'.
                format(section, configuration.config_file))

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
        """Run the Analysis.

        @raise Exception: An Analysis.project_name has not been defined
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

        self.input_directory = os.path.expanduser(path=self.input_directory)
        self.input_directory = os.path.expandvars(path=self.input_directory)

        if not os.path.isabs(self.input_directory):
            self.input_directory = os.path.join(Default.absolute_samples(), self.input_directory)

        self.output_directory = os.path.expanduser(path=self.output_directory)
        self.output_directory = os.path.expandvars(path=self.output_directory)

        if not os.path.isabs(self.output_directory):
            self.output_directory = os.path.join(Default.absolute_projects(), self.output_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the output_directory has to exist.

        if not os.path.isdir(self.output_directory):
            raise Exception('The Analysis output_directory {!r} does not exist.'.format(self.output_directory))

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

        if self.sas_file:

            # Populate a Collection from a SampleAnnotationSheet.

            self.sas_file = os.path.expanduser(path=self.sas_file)
            self.sas_file = os.path.expandvars(path=self.sas_file)

            if not os.path.isabs(self.sas_file) and not os.path.exists(self.sas_file):
                self.sas_file = os.path.join(self.project_directory, self.sas_file)

            self.collection = Collection.from_sas(file_path=self.input_directory,
                                                  file_type='Automatic',
                                                  name=self.project_name,
                                                  sas_file=self.sas_file,
                                                  sas_prefix=self.sas_prefix)

            if self.debug > 1:
                print '{!r} Collection name: {!r}'.format(self, self.collection.name)
                print self.collection.trace(1)

        else:

            # Create an empty Collection.

            self.collection = Collection()

    def report(self):
        """Create an Analysis report.
        """

        warnings.warn(
            "The 'report' method must be implemented in the sub-class.",
            UserWarning)

    def create_project_genome_directory(self):
        """Check and create an Analysis project_directory or genome_directory if necessary.

        @raise Exception: Output (genome) directory does not exist
        """

        if not os.path.isdir(self.genome_directory):
            answer = raw_input(
                "Output (genome) directory {!r} does not exist.\n"
                'Create? [Y/n] '.format(self.genome_directory))

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
                    'Output (genome) directory {!r} does not exist.'.format(self.genome_directory))

    def create_public_project_link(self, sub_directory=None):
        """Create a symbolic link from the web directory to the project directory if not already there.

        The link will be placed in the sub directory and contain
        the project name followed by a 128 bit hexadecimal UUID string.
        @param sub_directory: Analysis-specific directory
        @type sub_directory: str
        @return: Symbolic link to the project directory
        @rtype: str
        @raise Exception: Public HTML path does not exist
        """

        # The html_path consists of the absolute public_html directory and
        # the analysis-specific sub-directory.

        html_path = Default.absolute_public_html()

        if sub_directory:
            html_path = os.path.join(html_path, sub_directory)

        # As a safety measure, to prevent creation of rogue directory paths, the html_path directory has to exist.

        if not os.path.isdir(html_path):
            raise Exception(
                "The public HTML path {!r} does not exist.\n"
                "Please check the optional sub-directory name {!r}.".format(html_path, sub_directory))

        # The link_name consists of the absolute public_html directory,
        # the analysis-specific sub-directory, the project name and a 128 bit hexadecimal UUID string.

        link_name = os.path.join(html_path, string.join(words=(self.project_name, uuid.uuid4().hex), sep='_'))

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
                    target_name = os.path.join(html_path, target_name)
                if not os.path.exists(path=target_name):
                    # Both paths for os.path.samefile have to exist.
                    # Check for dangling symbolic links.
                    warnings.warn(
                        'Dangling symbolic link {!r} to {!r}'.format(path_name, target_name),
                        UserWarning)
                    continue
                if os.path.samefile(target_name, self.project_directory):
                    link_exists = True
                    link_final = path_name  # Reset link_final to the already existing path_name.
                    break

        if link_exists:
            # Ask the user to re-create the symbolic link.
            answer = raw_input(
                "Public HTML link {!r} to {!r} does exist.\n"
                "Re-create? [y/N] ".format(path_name, self.project_directory))

            if not answer or answer == 'N' or answer == 'n':
                print 'Public HTML link {!r} to {!r} not reset.'. \
                    format(path_name, self.project_directory)
            else:
                try:
                    os.remove(path_name)
                except OSError as exception:
                    # In principle, a race condition could occur as the directory
                    # could have been created after its existence has been checked.
                    if exception.errno != errno.ENOENT:
                        raise
                try:
                    os.symlink(os.path.relpath(self.project_directory, html_path), link_name)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
        else:
            # Ask the user to create a symbolic link.
            answer = raw_input(
                'Public HTML link {!r} to {!r} does not exist.\n'
                'Create? [Y/n] '.format(link_name, self.project_directory))

            if not answer or answer == 'Y' or answer == 'y':
                # In principle, a race condition could occur as the directory
                # could have been created after its existence has been checked.
                try:
                    os.symlink(os.path.relpath(self.project_directory, html_path), link_name)
                except OSError as exception:
                    if exception.errno != errno.EEXIST:
                        raise
            else:
                print 'Public HTML link {!r} to {!r} not set.'. \
                    format(link_name, self.project_directory)

        return link_final

    def ucsc_hub_write_hub(self, prefix=None):
        """Write a UCSC Track Hub hub.txt file into the project directory, above the genome directory.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
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

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
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

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
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

    def submit(self, drms_name=None):
        """Submit each DRMS object and pickle each Runnable object.

        @param drms_name: Only submit Executables linked to DRMS name
        @type drms_name: str
        """

        # Pickle all Runnable objects.

        for key in self.runnable_dict.keys():
            self.runnable_dict[key].to_pickler_file()

        # Submit all Executable objects of all Distributed Resource Management System objects.

        submit = 0

        for drms in self.drms_list:

            if drms_name:
                if drms_name == drms.name:
                    submit += 1
                else:
                    continue

            drms.submit(debug=self.debug)

            if self.debug:
                print repr(drms)
                print drms.trace(1)

        if drms_name:
            if drms_name == 'report':
                self.report()
            elif not submit:
                name_list = [drms.name for drms in self.drms_list]
                name_list.append('report')
                print 'Valid Analysis DRMS names are: {!r}'.format(name_list)


class Configuration(object):
    """Configuration class representing a UNIX-style initialisation (*.ini) file and
    an associated Python SafeConfigParser object to parse the file.

    Attributes:
    @ivar config_file: Configuration file path
    @type config_file: str | unicode
    @ivar config_parser: Python SafeConfigParser
    @type config_parser: SafeConfigParser
    """

    @staticmethod
    def section_from_instance(instance):
        """Get a configuration section from a Python instance.

        @param instance: A Python instance (or object)
        @type instance: object
        @return: Configuration section string
        @rtype: str
        """

        match = re.search(pattern=r'^<([^ ]+)\s+', string=repr(instance))

        if match:
            return match.group(1)
        else:
            return str()

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new Configuration object based on a configuration file path.

        Both, user and variable expansion gets applied to the file path.
        @param config_file: Configuration file path
        @type config_file: str | unicode
        @return: Configuration
        @rtype: Configuration
        @raise Exception: Configuration file does not exist
        """

        assert isinstance(config_file, (str, unicode))

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
            raise Exception(
                'Configuration file {!r} does not exist.'.format(configuration.config_file))

        return configuration

    def __init__(self, config_file=None, config_parser=None):
        """Initialise a Configuration object.

        @param config_file: Configuration file path
        @type config_file: str | unicode
        @param config_parser: Python SafeConfigParser
        @type config_parser: SafeConfigParser
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
        """Trace a Configuration object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
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
        @param config_section: Configuration section string
        @type config_section: str
        @param config_option: Configuration option string
        @type config_option: str
        @return: Expanded directory
        @rtype: str
        """

        directory = self.config_parser.get(config_section, config_option)
        directory = os.path.expanduser(path=directory)
        directory = os.path.expandvars(path=directory)

        return directory


class Default(object):
    """The Default class specifies the application or library default configuration.

    Attributes:
    @cvar global_default: Global Default
    @type global_default: Default
    @cvar global_file_path: Global configuration file
    @type global_file_path: str | unicode
    @ivar classpath_picard: Picard Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar classpath_illumina2bam: Illumina2bam Java Archive (JAR) class path directory
    @type classpath_illumina2bam: str | unicode
    @ivar classpath_snpeff: snpEff Java Archive (JAR) class path directory
    @type classpath_snpeff: str | unicode
    @ivar directory_home: Home directory for all data
    @type directory_home: str | unicode
    @ivar directory_runs_illumina: Sub-directory for Illumina runs
    @type directory_runs_illumina: str | unicode
    @ivar directory_sequences: Sub-directory for sequences
    @type directory_sequences: str | unicode
    @ivar directory_samples: Sub-directory for processed samples
    @type directory_samples: str | unicode
    @ivar directory_projects: Sub-directory for processed projects
    @type directory_projects: str | unicode
    @ivar directory_public_html: Sub-directory for the web server
    @type directory_public_html: str | unicode
    @ivar directory_genomes: Directory for genomes and their annotation
    @type directory_genomes: str | unicode
    @ivar directory_annotations: Sub-directory for genome annotations
    @type directory_annotations: str | unicode
    @ivar directory_gatk_bundle: Sub-directory for GATK bundle data
    @type directory_gatk_bundle: str | unicode
    @ivar directory_snpeff_data: snpEff database directory
    @type directory_snpeff_data: str | unicode
    @ivar indices: Python dict of program name key and index directory name value data
    @type indices: dict
    @ivar drms_implementation: DRMS implementation (e.g. Bash, SGE)
    @type drms_implementation: str
    @ivar drms_maximum_threads: DRMS maximum threads
    @type drms_maximum_threads: str
    @ivar drms_memory_limit_hard: DRMS memory limit hard
    @type drms_memory_limit_hard: str
    @ivar drms_memory_limit_soft: DRMS memory limit soft
    @type drms_memory_limit_soft: str
    @ivar drms_time_limit: DRMS time limit
    @type drms_time_limit: str
    @ivar drms_parallel_environment: DRMS parallel environment
    @type drms_parallel_environment: str
    @ivar drms_queue: DRMS queue
    @type drms_queue: str
    @ivar operator_e_mail: Operator e-mail
    @type operator_e_mail: str
    @ivar ucsc_host_name: UCSC Genome Browser host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...)
    @type ucsc_host_name: str
    @ivar url_protocol: URL protocol (i.e. HTTP)
    @type url_protocol: str
    @ivar url_host_name: URL host name
    @type url_host_name:str
    @ivar url_relative_projects: Sub-directory for analysis projects
    @type url_relative_projects: str
    """

    global_default = None

    global_file_path = '~/.bsfpython.ini'
    global_file_path = os.path.expanduser(path=global_file_path)
    global_file_path = os.path.expandvars(path=global_file_path)

    @staticmethod
    def get_global_default():
        """Get the global Default configuration and initialise it, if not already done so.

        @return: Default
        @rtype: Default
        """

        if not Default.global_default:
            Default.global_default = Default.from_global_file()

        return Default.global_default

    @classmethod
    def from_global_file(cls):
        """Create a new Default object from the global default configuration file.

        The default configuration is based on the file $HOME/.bsfpython.ini in the user's home directory.
        @return: Default
        @rtype: Default
        """

        return cls.from_config_file(config_file=Default.global_file_path)

    @classmethod
    def from_config_file(cls, config_file):
        """Create a new Default object from a UNIX-style configuration file.

        @param config_file: UNIX-style configuration file
        @type config_file: str | unicode
        @return: Default
        @rtype: Default
        """

        return cls.from_Configuration(configuration=Configuration.from_config_file(config_file=config_file))

    @classmethod
    def from_Configuration(cls, configuration):
        """Create a new Default objects from a Configuration object.

        @param configuration: Configuration
        @type configuration: Configuration
        @return: Default
        @rtype: Default
        """

        assert isinstance(configuration, Configuration)

        default = cls()

        default.set_Configuration(configuration=configuration)

        return default

    def __init__(self, classpath_gatk=None, classpath_illumina2bam=None, classpath_picard=None, classpath_snpeff=None,
                 directory_home=None, directory_runs_illumina=None, directory_sequences=None, directory_samples=None,
                 directory_projects=None, directory_public_html=None, directory_genomes=None,
                 directory_annotations=None, directory_gatk_bundle=None, directory_snpeff_data=None,
                 indices=None, drms_implementation=None,
                 drms_maximum_threads=None, drms_memory_limit_hard=None, drms_memory_limit_soft=None,
                 drms_time_limit=None, drms_parallel_environment=None, drms_queue=None,
                 operator_e_mail=None, operator_sequencing_centre=None, ucsc_host_name=None, url_protocol=None,
                 url_host_name=None, url_relative_projects=None):
        """Initialise a Default object.

        @param classpath_gatk: Genome Analysis Toolkit Java Archive (JAR) class path directory
        @type classpath_gatk: str | unicode
        @param classpath_illumina2bam: Illumina2bam Java Archive (JAR) class path directory
        @type classpath_illumina2bam: str | unicode
        @param classpath_picard: Picard Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param classpath_snpeff: snpEff Java Archive (JAR) class path directory
        @type classpath_snpeff: str | unicode
        @param directory_home: Home directory for all data
        @type directory_home: str | unicode
        @param directory_runs_illumina: Sub-directory for Illumina runs
        @type directory_runs_illumina: str | unicode
        @param directory_sequences: Sub-directory for sequences
        @type directory_sequences: str | unicode
        @param directory_samples: Sub-directory for processed samples
        @type directory_samples: str | unicode
        @param directory_projects: Sub-directory for processed projects
        @type directory_projects: str | unicode
        @param directory_public_html: Sub-directory for the web server
        @type directory_public_html: str | unicode
        @param directory_genomes: Directory for genomes and their annotation
        @type directory_genomes: str | unicode
        @param directory_annotations: Sub-directory for genome annotations
        @type directory_annotations: str | unicode
        @param directory_gatk_bundle: Sub-directory for GATK bundle data
        @type directory_gatk_bundle: str | unicode
        @param directory_snpeff_data: snpEff database directory
        @type directory_snpeff_data: str | unicode
        @param indices: Python dict of program name key and index directory name value data
        @type indices: dict
        @param drms_implementation: DRMS implementation (e.g. Bash, SGE)
        @type drms_implementation: str
        @param drms_maximum_threads: DRMS maximum threads
        @type drms_maximum_threads: str
        @param drms_memory_limit_hard: DRMS memory limit hard
        @type drms_memory_limit_hard: str
        @param drms_memory_limit_soft: DRMS memory limit soft
        @type drms_memory_limit_soft: str
        @param drms_time_limit: DRMS time limit
        @type drms_time_limit: str
        @param drms_parallel_environment: DRMS parallel environment
        @type drms_parallel_environment: str
        @param drms_queue: DRMS queue
        @type drms_queue: str
        @param operator_e_mail: Operator e-mail
        @type operator_e_mail: str
        @param operator_sequencing_centre: BAM sequencing centre code
        @type operator_sequencing_centre: str
        @param ucsc_host_name: UCSC Genome Browser host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...)
        @type ucsc_host_name: str
        @param url_protocol: URL protocol (i.e. HTTP)
        @type url_protocol: str
        @param url_host_name: URL host name
        @type url_host_name:str
        @param url_relative_projects: Sub-directory for analysis projects
        @type url_relative_projects: str
        """

        # Set Java class path information.

        if classpath_gatk:
            self.classpath_gatk = classpath_gatk
        else:
            self.classpath_gatk = str()

        if classpath_illumina2bam:
            self.classpath_illumina2bam = classpath_illumina2bam
        else:
            self.classpath_illumina2bam = str()

        if classpath_picard:
            self.classpath_picard = classpath_picard
        else:
            self.classpath_picard = str()

        if classpath_snpeff:
            self.classpath_snpeff = classpath_snpeff
        else:
            self.classpath_snpeff = str()

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

        if directory_gatk_bundle:
            self.directory_gatk_bundle = directory_gatk_bundle
        else:
            self.directory_gatk_bundle = str()

        if directory_snpeff_data:
            self.directory_snpeff_data = directory_snpeff_data
        else:
            self.directory_snpeff_data = str()

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

        if drms_time_limit:
            self.drms_time_limit = drms_time_limit
        else:
            self.drms_time_limit = str()

        if drms_parallel_environment:
            self.drms_parallel_environment = drms_parallel_environment
        else:
            self.drms_parallel_environment = str()

        if drms_queue:
            self.drms_queue = drms_queue
        else:
            self.drms_queue = str()

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

    def set_Configuration(self, configuration):
        """Set instance variables of a Default object via a section of a Configuration object.

        For each instance variable a configuration option has to be present.
        @param configuration: Configuration
        @type configuration: Configuration
        """

        assert isinstance(configuration, Configuration)

        cp = configuration.config_parser

        # Reading configuration cannot be done via a single Python dict,
        # because each option really needs defining.

        section = 'classpath'

        self.classpath_gatk = cp.get(section=section, option='gatk')
        self.classpath_illumina2bam = cp.get(section=section, option='illumina2bam')
        self.classpath_picard = cp.get(section=section, option='picard')
        self.classpath_snpeff = cp.get(section=section, option='snpeff')

        section = 'directories'

        self.directory_home = cp.get(section=section, option='home')
        self.directory_runs_illumina = cp.get(section=section, option='runs_illumina')
        self.directory_sequences = cp.get(section=section, option='sequences')
        self.directory_samples = cp.get(section=section, option='samples')
        self.directory_projects = cp.get(section=section, option='projects')
        self.directory_public_html = cp.get(section=section, option='public_html')
        self.directory_genomes = cp.get(section=section, option='genomes')
        self.directory_annotations = cp.get(section=section, option='annotations')
        self.directory_gatk_bundle = cp.get(section=section, option='gatk_bundle')
        self.directory_snpeff_data = cp.get(section=section, option='snpeff_data')

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
        self.drms_time_limit = cp.get(section=section, option='time_limit')
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

    @staticmethod
    def absolute_home():
        """
        Get the absolute directory path for the home directory.

        @return: Absolute path to the home directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        return default.directory_home

    @staticmethod
    def absolute_runs_illumina():
        """Get the absolute directory path for Illumina runs.

        @return: Absolute path to the Illumina runs directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_runs_illumina):
            return default.directory_runs_illumina
        else:
            return os.path.join(default.directory_home, default.directory_runs_illumina)

    @staticmethod
    def absolute_sequences():
        """Get the absolute directory path for processed lanes.

        @return: Absolute path to the processed lanes directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_sequences):
            return default.directory_sequences
        else:
            return os.path.join(default.directory_home, default.directory_sequences)

    @staticmethod
    def absolute_projects():
        """Get the absolute directory path for projects.

        @return: Absolute path to the projects directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_projects):
            return default.directory_projects
        else:
            return os.path.join(default.directory_home, default.directory_projects)

    @staticmethod
    def absolute_samples():
        """Get the absolute directory path for processed samples.

        @return: Absolute path to the processed samples directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_samples):
            return default.directory_samples
        else:
            return os.path.join(default.directory_home, default.directory_samples)

    @staticmethod
    def absolute_public_html():
        """Get the absolute directory path for public HTML documents.

        @return: Absolute path to the public HTML directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_public_html):
            return default.directory_public_html
        else:
            return os.path.join(default.directory_home, default.directory_public_html)

    @staticmethod
    def absolute_gatk_bundle(gatk_bundle_version, genome_version):
        """Get the absolute directory path for the Genome Analysis Toolkit bundle.

        @param gatk_bundle_version: The GATK bundle version
        @type gatk_bundle_version: str
        @param genome_version: The genome version (e.g. b37, ...)
        @type genome_version: str
        @return Absolute path to the GATK bundle directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        file_path = str(default.directory_gatk_bundle)

        if gatk_bundle_version:
            file_path = os.path.join(file_path, gatk_bundle_version)

        if genome_version:
            file_path = os.path.join(file_path, genome_version)

        return file_path

    @staticmethod
    def absolute_genomes(genome_version):
        """Get the absolute directory path for genomes.

        @param genome_version: The genome version (e.g. mm10, ...)
        @type genome_version: str
        @return: Absolute path to the genomes directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if genome_version:
            return os.path.join(default.directory_genomes, genome_version)
        else:
            return default.directory_genomes

    @staticmethod
    def absolute_genome_annotation(genome_version):
        """Get the absolute directory path for genome annotation.

        @param genome_version: The genome version (e.g. mm10, ...)
        @type genome_version: str
        @return: Absolute path to the genome annotation directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        return os.path.join(default.directory_genomes, genome_version, default.directory_annotations)

    @staticmethod
    def absolute_genome_fasta(genome_version, genome_index):
        """Get the absolute file path to a genome in FASTA format.

        @param genome_version: Genome version (e.g. mm10, ...)
        @type genome_version: str
        @param genome_index: Genome index (e.g. bowtie2, ...)
        @type genome_index: str
        @return: Absolute path to the genome FASTA file
        @rtype: str | unicode
        @raise Exception: Unknown genome index name
        """

        default = Default.get_global_default()

        if not genome_index in default.indices:
            raise Exception(
                'Unknown genome index name {!r}.'.format(genome_index))

        return os.path.join(Default.absolute_genomes(genome_version),
                            default.indices[genome_index],
                            genome_version + '.fa')

    @staticmethod
    def url_absolute_base():
        """Return the absolute URL to the web site.

        @return: URL string
        @rtype: str
        """

        default = Default.get_global_default()

        return '{}://{}'.format(default.url_protocol, default.url_host_name)

    @staticmethod
    def url_absolute_projects():
        """Return the absolute URL to the analysis projects directory.

        @return: URL string
        @rtype: str
        """

        default = Default.get_global_default()

        return string.join(words=(default.url_absolute_base(), default.url_relative_projects), sep='/')


class DRMS(object):
    """The Distributed Resource Management System (DRMS) class represents a Distributed Resource Management System or
    batch job scheduler.

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar work_directory: Work directory path
    @type work_directory: str
    @ivar implementation: Implementation (e.g. SGE, ...)
    @type implementation: str
    @ivar memory_free_mem: Memory limit (free)
    @type memory_free_mem: str
    @ivar memory_limit_hard: Memory limit (hard)
    @type memory_limit_hard: str
    @ivar memory_limit_soft: Memory limit (soft)
    @type memory_limit_soft: str
    @ivar parallel_environment: Parallel environment
    @type parallel_environment: str
    @ivar queue: Queue
    @type queue: str
    @ivar threads: Number of threads
    @type threads: int
    @ivar hold: Hold on job scheduling
    @type hold: str
    @ivar is_script: Executable objects represent shell scripts,
        or alternatively binary programs
    @type is_script: bool
    @ivar executables: Python list of Executable objects
    @type executables: list
    """

    @classmethod
    def from_Analysis(cls, name, work_directory, analysis):
        """Create a DRMS object from an Analysis object.

        @param name: Name
        @type name: str
        @param work_directory: Work directory
        @type work_directory: str
        @param analysis: Analysis
        @type analysis: Analysis
        @return: DRMS object
        @rtype: DRMS
        """

        assert isinstance(analysis, Analysis)

        drms = cls(name=name, work_directory=work_directory)

        # Set a minimal set of global defaults.

        drms.set_Default(default=Default.get_global_default())

        # A "Bio.BSF.DRMS" section specifies defaults for all DRMS objects of an Analysis.

        section = string.join(words=(__name__, cls.__name__), sep='.')
        drms.set_Configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 1:
            print 'DRMS configuration section: {!r}.'.format(section)

        # A "Bio.BSF.Analysis.*.DRMS" pseudo-class section specifies
        # Analysis-specific options for the DRMS.

        section = string.join(words=(analysis.configuration.section_from_instance(analysis), 'DRMS'), sep='.')
        drms.set_Configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 1:
            print 'DRMS configuration section: {!r}.'.format(section)

        # A "Bio.BSF.Analysis.*.DRMS.name" section specifies defaults
        # for a particular DRMS objects of an Analysis.

        section = string.join(words=(Configuration.section_from_instance(analysis), 'DRMS', drms.name), sep='.')
        drms.set_Configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 1:
            print 'DRMS configuration section: {!r}.'.format(section)

        return drms

    @classmethod
    def from_Configuration(cls, name, work_directory, configuration, section):
        """Create a DRMS object from a Configuration object.

        @param name: Name
        @type name: str
        @param work_directory: Work directory
        @type work_directory: str
        @param configuration: Configuration object
        @type configuration: Configuration
        @param section: Configuration section string
        @type section: str
        @return: DRMS object
        @rtype: DRMS
        """

        assert isinstance(configuration, Configuration)

        drms = cls(name=name, work_directory=work_directory)

        # Set a minimal set of global defaults before setting the Configuration.

        drms.set_Default(default=Default.get_global_default())
        drms.set_Configuration(configuration=configuration, section=section)

        return drms

    def __init__(self, name, work_directory,
                 implementation=None,
                 memory_free_mem=None,
                 memory_free_swap=None,
                 memory_free_virtual=None,
                 memory_limit_hard=None,
                 memory_limit_soft=None,
                 time_limit=None,
                 parallel_environment=None,
                 queue=None,
                 threads=1,
                 hold=None,
                 is_script=False,
                 executables=None):
        """Initialise a DRMS object.

        @param name: Name
        @type name: str
        @param work_directory: Work directory
        @type work_directory: str
        @param implementation: Implementation (e.g. SGE, ...)
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
        @param is_script: Executable objects represent shell scripts,
            or alternatively binary programs
        @type is_script: bool
        @param executables: Python list of Executable objects
        @type executables: list
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

        if time_limit:
            self.time_limit = time_limit
        else:
            self.time_limit = str()

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
        """Trace a DRMS object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
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
        output += '{}  time_limit:           {!r}\n'. \
            format(indent, self.time_limit)
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
        """Set instance variables of a DRMS object via a section of a Configuration object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: Configuration
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            raise Exception(
                'Section {!r} not defined in Configuration file {!r}.'.
                format(section, configuration.config_file))

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

        if configuration.config_parser.has_option(section=section, option='time_limit'):
            self.time_limit = configuration.config_parser.get(section=section,
                                                              option='time_limit')

        if configuration.config_parser.has_option(section=section, option='parallel_environment'):
            self.parallel_environment = configuration.config_parser.get(section=section,
                                                                        option='parallel_environment')

        if configuration.config_parser.has_option(section=section, option='queue'):
            self.queue = configuration.config_parser.get(section=section,
                                                         option='queue')

        if configuration.config_parser.has_option(section=section, option='threads'):
            self.threads = configuration.config_parser.get(section=section,
                                                           option='threads')

    def set_Default(self, default):
        """Set instance variables of a DRMS object via a Default object.

        @param default: Default
        @type default: Default
        """

        assert isinstance(default, Default)

        self.implementation = default.drms_implementation
        # is_script
        # memory_free_mem
        # memory_free_swap
        # memory_free_virtual
        self.memory_limit_hard = default.drms_memory_limit_hard
        self.memory_limit_soft = default.drms_memory_limit_soft
        self.time_limit = default.drms_time_limit
        self.parallel_environment = default.drms_parallel_environment
        self.queue = default.drms_queue
        # threads

    def add_Executable(self, executable):
        """Add a Executable object.

        @param executable: Executable
        @type executable: Executable
        """

        assert isinstance(executable, Executable)

        self.executables.append(executable)

    def submit(self, debug=0):
        """Submit a command line for each Executable object.

        @param debug: Debug level
        @type debug: int
        """

        # Dynamically import the module specific for the configured DRMS implementation.

        module = importlib.import_module(string.join(words=(__name__, 'DRMS', self.implementation), sep='.'))

        module.submit(drms=self, debug=debug)


class Command(object):
    """Command class representing one (subordinate) command,
    its options and arguments and possibly another subordinate command.

    Attributes:
    @ivar command: Command or program
    @type command: str
    @ivar options: Python dict of option keys and values
    @type options: dict
    @ivar arguments: Python list of arguments
    @type arguments: list
    @ivar sub_command: Subordinate Command
    @type sub_command: Command
    """

    def __init__(self, command, options=None, arguments=None, sub_command=None):
        """Initialise a Command object.

        @param command: Command
        @type command: str
        @param options: Python dict of program option and value pairs
        @type options: dict
        @param arguments: Python list of program arguments
        @type arguments: list
        @param sub_command: Subordinate Command object
        @type sub_command: Command
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
        """Trace a Command object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  command:            {!r}\n'. \
            format(indent, self.command)

        # List all options

        output += '{}  options:\n'.format(indent)

        for key in self.options.keys():
            output += '{}    key: {!r} Argument objects:\n'.format(indent, key)
            for argument in self.options[key]:
                output += argument.trace(level=level + 2)

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
        """Add an Argument or one of its sub-classes.

        @param argument: Argument
        @type argument: Argument
        @param override: Override existing Argument without warning
        @type override: bool
        """

        assert isinstance(argument, Argument)
        assert isinstance(override, bool)

        if not override and argument.key in self.options:
            warnings.warn(
                'Adding an Argument with key {!r} that exits already in Command {!r}.'.
                format(argument.key, self.command),
                UserWarning)

        if argument.key in self.options:
            arguments_list = self.options[argument.key]
        else:
            arguments_list = list()
            self.options[argument.key] = arguments_list

        arguments_list.append(argument)

    def add_SwitchLong(self, key, override=False):
        """Initialise and add a SwitchLong object.

        @param key: Key
        @type key: str
        @param override: Override existing Argument without warning
        @type override: bool
        """

        self.add_Argument(argument=SwitchLong(key=key), override=override)

    def add_SwitchShort(self, key, override=False):
        """Initialise and add a SwitchShort object.

        @param key: Key
        @type key: str
        @param override: Override existing Argument without warning
        @type override: bool
        """

        self.add_Argument(argument=SwitchShort(key=key), override=override)

    def add_OptionLong(self, key, value, override=False):
        """Initialise and add an OptionLong object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing Argument without warning
        @type override: bool
        """

        self.add_Argument(argument=OptionLong(key=key, value=value), override=override)

    def add_OptionShort(self, key, value, override=False):
        """Initialise and add an OptionShort object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing Argument without warning
        @type override: bool
        """

        self.add_Argument(argument=OptionShort(key=key, value=value), override=override)

    def add_OptionPair(self, key, value, override=False):
        """Initialise and add an OptionPair object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing Argument without warning
        @type override: bool
        """

        self.add_Argument(argument=OptionPair(key=key, value=value), override=override)

    def set_argument(self, argument, override):
        """Set an Argument or one of its sub-classes.

        @param argument: Argument
        @type argument: Argument
        @param override: Override existing Argument without warning
        @type override: bool
        """
        assert isinstance(argument, Argument)
        assert isinstance(override, bool)

        if not override and argument.key in self.options:
            warnings.warn(
                'Setting an Argument with key {!r} that exits already in Command {!r}.'.
                format(argument.key, self.command),
                UserWarning)

        self.options[argument.key] = [argument]

    def set_switch_long(self, key, override=False):
        """Initialise and set a SwitchLong object.

        @param key: Key
        @type key: str
        @param override: Override existing Argument without warning
        @type override: bool
        """
        self.set_argument(argument=SwitchLong(key=key), override=override)

    def set_switch_short(self, key, override=False):
        """Initialise and set a SwitchShort object.

        @param key: Key
        @type key: str
        @param override: Override existing Argument without warning
        @type override: bool
        """
        self.set_argument(argument=SwitchShort(key=key), override=override)

    def set_option_long(self, key, value, override=False):
        """Initialise and set an OptionLong object.

        @param key: Key
        @type key: str
        @param value: Value
        @param override: Override existing Argument without warning
        @type override: bool
        @type value: str | unicode
        """
        self.set_argument(argument=OptionLong(key=key, value=value), override=override)

    def set_option_short(self, key, value, override=False):
        """Initialise and set an OptionShort object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing Argument without warning
        @type override: bool
        """
        self.set_argument(argument=OptionShort(key=key, value=value), override=override)

    def set_option_pair(self, key, value, override=False):
        """Initialise and set an OptionPair object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing Argument without warning
        @type override: bool
        """
        self.set_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_Configuration(self, configuration, section):
        """Set instance variables of a Command object via a section of a Configuration object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: Configuration
        @type configuration: Configuration
        @param section: Configuration file section, defaults to instance class
        @type section: str
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            warnings.warn(
                'Section {!r} not defined in Configuration file {!r}.'.
                format(section, configuration.config_file),
                UserWarning)

            return

        # The configuration section is available.

        for option in configuration.config_parser.options(section=section):
            argument = Argument.from_key_value(key=option,
                                               value=configuration.config_parser.get(section=section,
                                                                                     option=option))

            self.add_Argument(argument=argument, override=False)

    def command_list(self):
        """Assemble the command line from program, options and arguments.

        @return: Python list of program, options, switches and arguments
        @rtype: list
        """

        command_line = list()

        if self.command:
            command_line.append(self.command)

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            options_list = self.options[key]
            for argument in options_list:
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
                    warnings.warn(
                        'Unexpected object {!r} in Command.options dict.'.
                        format(argument),
                        UserWarning)

        # Add all arguments.

        for argument in self.arguments:
            command_line.append(str(argument))

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line.extend(self.sub_command.command_list())

        return command_line

    def command_str(self):
        """Assemble the command line from program, options, switches and arguments.

        @return: A Python str of program, options, switches and arguments
        @rtype: str
        """

        command_line = str()

        if self.command:
            command_line += self.command

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            options_list = self.options[key]
            for argument in options_list:
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
                    warnings.warn(
                        'Unexpected object {!r} in Command.options dict.'.
                        format(argument),
                        UserWarning)

        # Add all arguments.

        for argument in self.arguments:
            command_line += ' '
            command_line += argument

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line += ' '
            command_line += self.sub_command.command_str()

        return command_line


class Executable(Command):
    """The Executable class represents an executable program,
    its options and arguments.

    Attributes:
    @ivar name: Name in the context of a DRMS dependency
    @type name: str
    @ivar program: Program (executable or full file path)
    @type program: str
    @ivar options: Python dict of option keys and values
    @type options: dict
    @ivar arguments: Python list of arguments
    @type arguments: list
    @ivar sub_command: Subordinate Command
    @type sub_command: Command
    @ivar stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
    @type stdout_path: str | unicode
    @ivar stderr_path: Standard error (STDERR) redirection in Bash (2>word)
    @type stderr_path: str | unicode
    @ivar dependencies: Python list of Executable name strings in the
        context of DRMS dependencies
    @type dependencies: list
    @ivar hold: Hold on job scheduling
    @type hold: str
    @ivar submit: Submit the Executable into the DRMS
    @type submit: bool
    @ivar process_identifier: Process identifier
    @type process_identifier: str
    @ivar process_name: Process name
    @type process_name: str
    """

    @classmethod
    def from_Analysis(cls, name, program, analysis):
        """Create an Executable object from an Analysis object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param analysis: Analysis
        @type analysis: Analysis
        @return: Executable object
        @rtype: Executable
        """

        assert isinstance(analysis, Analysis)

        # Initialise an Executable object with default values.

        executable = cls(name=name, program=program)

        section = Configuration.section_from_instance(executable)

        # For plain Executable objects append the value of the
        # Executable.command to make this more meaningful.

        if section == 'Bio.BSF.Executable':
            section += '.'
            section += executable.command

        if analysis.debug > 1:
            print 'Executable configuration section: {!r}.'.format(section)

        executable.set_Configuration(configuration=analysis.configuration, section=section)

        return executable

    @classmethod
    def from_analysis_runnable(cls, analysis, runnable_name):
        """Create an Executable to submit a Runnable into a DRMS.

        @param analysis: Analysis
        @type analysis: Analysis
        @param runnable_name: Runnable name
        @type runnable_name: str
        @return: Executable
        @rtype: Executable
        @raise Exception: A Runnable.name does not exist in Analysis.name
        """

        assert isinstance(analysis, Analysis)

        if not runnable_name in analysis.runnable_dict:
            raise Exception("A Runnable object with name {!r} does not exist in the Analysis object with name {!r}.".
                            format(runnable_name, analysis.project_name))

        runnable = analysis.runnable_dict[runnable_name]
        executable = cls(name=runnable.name, program=Runnable.runner_script)
        executable.set_Configuration(configuration=analysis.configuration, section=runnable.code_module)
        executable.add_OptionLong(key='pickler_path', value=runnable.pickler_path)
        # executable.add_OptionLong(key='runnable_name', value=runnable.name)
        # executable.add_OptionLong(key='debug', value=str(analysis.debug))

        return executable

    @classmethod
    def from_Configuration(cls, name, program, configuration, section):
        """Create an Executable object from a Configuration object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param configuration: Configuration
        @type configuration: Configuration
        @param section: Configuration section string
        @type section: str
        @return: Executable
        @rtype: Executable
        """

        assert isinstance(configuration, Configuration)

        executable = cls(name=name, program=program)

        executable.set_Configuration(configuration=configuration, section=section)

        return executable

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None):
        """Initialise an Executable object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options: Python dict of program option and value pairs
        @type options: dict
        @param arguments: Python list of program arguments
        @type arguments: list
        @param sub_command: Subordinate Command
        @type sub_command: Command
        @param stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (STDERR) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python list of Executable
            name strings in the context of DRMS dependencies
        @type dependencies: list
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the Executable into the DRMS
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
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

        self.submit = submit

        if process_identifier:
            self.process_identifier = process_identifier
        else:
            self.process_identifier = str()

        if process_name:
            self.process_name = process_name
        else:
            self.process_name = str()

    def trace(self, level):
        """Trace an Executable object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
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
        output += '{}  submit:             {!r}\n'. \
            format(indent, self.submit)
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

        @return: Python list of program, options and arguments
        @rtype: list
        """

        command = list()

        command.extend(super(Executable, self).command_list())

        # The stdout_path and stderr_path gets appended in specific modules.

        return command

    def command_str(self):
        """Assemble the command line from program, options, switches and arguments.

        @return: A Python str of program, options, switches and arguments
        @rtype: str
        """

        command = str()

        command += super(Executable, self).command_str()

        # The stdout_path and stderr_path gets appended in specific modules.

        return command


class Runnable(object):
    """The Runnable class holds all information to run one or more Executable objects through the
    Runner script.

    Attributes:
    @cvar runner_script: Name of the Runner script
    @type runner_script: str | unicode
    @ivar name: Name
    @type name: str
    @ivar code_module: The name of a module, usually in Runnables that implements the logic required to run
        Executable objects via the Runner script.
    @type code_module: str
    @ivar executable_dict: Python dict of Python str (Executable.name) key data and Executable value data
    @type executable_dict: dict
    @ivar file_path_dict: Python dict of Python str (name) key data and Python str (file_path) value data
    @type file_path_dict: dict
    @ivar working_directory: Working directory to write Pickler files
    @type working_directory: str | unicode
    """

    runner_script = 'bsf_runner.py'

    @staticmethod
    def process_stream(file_type, file_handle, thread_lock, file_path=None, debug=0):
        """Runnable function to process STDOUT or STDERR from the child process as a thread.

        @param file_type: File handle type STDOUT or STDERR
        @type file_type: str
        @param file_handle: The STDOUT or STDERR file handle
        @type file_handle: file
        @param thread_lock: A Python threading.Lock object
        @type thread_lock: thread.lock
        @param file_path: STDOUT file path
        @type file_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @raise Exception: The file_type has to be either STDOUT or STDERR
        """

        if file_type not in ('STDOUT', 'STDERR'):
            raise Exception('The file_type has to be either STDOUT or STDERR.')

        thread_lock.acquire(True)
        if debug > 0:
            print '[{}] Started Runner {} processor in module {}.'. \
                format(datetime.datetime.now().isoformat(), file_type, __name__)
        output_file = None
        if file_path:
            output_file = open(file_path, 'w')
            if debug > 0:
                print '[{}] Opened {} file {!r}.'. \
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
                print '[{}] Closed {} file {!r}.'. \
                    format(datetime.datetime.now().isoformat(), file_type, file_path)
        thread_lock.release()

    @staticmethod
    def process_stdout(stdout_handle, thread_lock, stdout_path=None, debug=0):
        """Runnable function to process STDOUT from the child process as a thread.

        @param stdout_handle: The STDOUT file handle
        @type stdout_handle: file
        @param thread_lock: A Python threading.Lock object
        @type thread_lock: thread.lock
        @param stdout_path: STDOUT file path
        @type stdout_path: str | unicode
        @param debug: Debug level
        @type debug: int
        """

        return Runnable.process_stream(file_type='STDOUT', file_handle=stdout_handle,
                                       thread_lock=thread_lock, file_path=stdout_path,
                                       debug=debug)

    @staticmethod
    def process_stderr(stderr_handle, thread_lock, stderr_path=None, debug=0):
        """Runnable function to process STDERR from the child process as a thread.

        @param stderr_handle: The STDERR file handle
        @type stderr_handle: file
        @param thread_lock: A Python threading.Lock object
        @type thread_lock: thread.lock
        @param stderr_path: STDOUT file path
        @type stderr_path: str | unicode
        @param debug: Debug level
        @type debug: int
        """

        return Runnable.process_stream(file_type='STDERR', file_handle=stderr_handle,
                                       thread_lock=thread_lock, file_path=stderr_path,
                                       debug=debug)

    @staticmethod
    def run(executable, max_loop_counter=1, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Runnable function to run an Executable object as Python subprocess.

        @param executable: Executable
        @type executable: Executable
        @param max_loop_counter: Maximum number of retries
        @type max_loop_counter: int
        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """

        on_posix = 'posix' in sys.builtin_module_names

        loop_counter = 0
        child_return_code = 0

        while loop_counter < max_loop_counter:

            child_process = Popen(args=executable.command_list(),
                                  bufsize=0,
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
                    print '[{}] Child process {!r} failed with exit code {}'. \
                        format(datetime.datetime.now().isoformat(), executable.name, +child_return_code)
                loop_counter += 1
            elif child_return_code < 0:
                if debug > 0:
                    print '[{}] Child process {!r} received signal {}.'. \
                        format(datetime.datetime.now().isoformat(), executable.name, -child_return_code)
            else:
                if debug > 0:
                    print '[{}] Child process {!r} completed successfully {}.'. \
                        format(datetime.datetime.now().isoformat(), executable.name, +child_return_code)
                break

        else:
            if debug > 0:
                print '[{}] Runnable {!r} exceeded the maximum re-run counter {}.' \
                    .format(datetime.datetime.now().isoformat(), executable.name, max_loop_counter)

        return child_return_code

    @staticmethod
    def evaluate_return_code(executable, return_code):
        """Evaluate a return code from the run method.

        @param executable: Executable
        @type executable: Executable
        @param return_code: Return code
        @type return_code: int
        """

        if return_code > 0:
            print '[{}] Child process {!r} failed with return code {}'. \
                format(datetime.datetime.now().isoformat(), executable.name, +return_code)
        elif return_code < 0:
            print '[{}] Child process {!r} received signal {}.'. \
                format(datetime.datetime.now().isoformat(), executable.name, -return_code)
        else:
            print '[{}] Child process {!r} completed with return code {}.'. \
                format(datetime.datetime.now().isoformat(), executable.name, +return_code)

    def __init__(self, name, code_module, working_directory, file_path_dict=None, executable_dict=None):
        """Initialise a Runnable object.

        @param name: Name
        @type name: str
        @param code_module: The Runnables module that implements the logic for this Runnable
        @type code_module: str
        @param working_directory: Working directory for writing a Python Pickler file
        @type working_directory: str | unicode
        @param file_path_dict: Python dict of Python str (name) key data and Python str (file_path) value data
        @type file_path_dict: dict
        @param executable_dict: Python dict of Python str (Executable.name) key data and Executable value data
        @type executable_dict: dict
        """

        self.name = name
        self.code_module = code_module
        self.working_directory = working_directory

        if file_path_dict:
            self.file_path_dict = file_path_dict
        else:
            self.file_path_dict = dict()

        if executable_dict:
            self.executable_dict = executable_dict
        else:
            self.executable_dict = dict()

    def add_executable(self, executable):
        """Add an Executable

        @param executable: Executable
        @type executable: Executable
        @raise Exception: An Executable.name already exists in the Runnable object
        """

        if not executable:
            return

        if executable.name in self.executable_dict:
            raise Exception("An Executable object with name {!r} already exists in Runnable object {!r}.".
                            format(executable.name, self.name))
        else:
            self.executable_dict[executable.name] = executable

    def run_executable(self, name):
        """Run an Executable defined in a Runnable.

        @param name: Executable name
        @type name: str
        @raise Exception: Child process failed with return code or received a signal
        """

        executable = self.executable_dict[name]
        child_return_code = Runnable.run(executable=executable)

        if child_return_code > 0:
            raise Exception('[{}] Child process {!r} failed with return code {}'.
                            format(datetime.datetime.now().isoformat(), executable.name, +child_return_code))
        elif child_return_code < 0:
            raise Exception('[{}] Child process {!r} received signal {}.'.
                            format(datetime.datetime.now().isoformat(), executable.name, -child_return_code))

    @property
    def pickler_path(self):
        """Get the Python Pickler file path.

        @return: Python Pickler file path
        @rtype: str
        """

        return os.path.join(self.working_directory, string.join(words=(self.name, 'pkl'), sep='.'))

    def to_pickler_file(self):
        """Write this object as a Python Pickler file into the working directory.
        """

        pickler_file = open(self.pickler_path, 'wb')
        pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
        pickler.dump(obj=self)
        pickler_file.close()

    @classmethod
    def from_picker_file(cls, file_path):
        """Create a Runnable object from a Python Pickler file via Python Unpickler.

        @param file_path: File path to a Picker file
        @type file_path: str | unicode
        @return: Runnable
        @rtype: Runnable
        """

        pickler_file = open(file_path, 'rb')
        unpickler = Unpickler(file=pickler_file)
        runnable = unpickler.load()
        pickler_file.close()

        assert isinstance(runnable, Runnable)

        return runnable
