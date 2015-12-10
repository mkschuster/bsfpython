"""bsf

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
import shutil
from stat import *
from subprocess import PIPE, Popen
import sys
from threading import Lock, Thread
import time
import uuid
import warnings

from bsf import defaults
from bsf.data import Collection, Sample
from bsf.argument import *


class Analysis(object):
    """The C{Analysis} class represents a high-level analysis that may run one or more
    C{Executable} objects (programs).

    Attributes:
    @ivar configuration: C{Configuration}
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
    @ivar drms_list: Python C{list} of C{DRMS} objects
    @type drms_list: list[DRMS]
    @ivar runnable_dict: Python C{dict} of Python C{str} (C{Runnable.name}) key data and C{Runnable} value data
    @type runnable_dict: dict[Runnable.name, Runnable]
    @ivar collection: C{Collection}
    @type collection: Collection
    @ivar comparisons: Python C{dict} of comparisons
    @type comparisons: dict[str, any]
    @ivar samples: Python C{list} of C{Sample} objects
    @type samples: list[Sample]
    """

    @classmethod
    def from_config_file_path(cls, config_path):
        """Create a new C{Analysis} object from a UNIX-style configuration file path via the C{Configuration} class.

        @param config_path: UNIX-style configuration file path
        @type config_path: str | unicode
        @return: C{Analysis}
        @rtype: Analysis
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{Analysis} object from a C{Configuration} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @return: C{Analysis}
        @rtype: Analysis
        """

        assert isinstance(configuration, Configuration)

        # Set a minimal set of global defaults.

        default = Default.get_global_default()

        analysis = cls(configuration=configuration, e_mail=default.operator_e_mail)

        # A "module.class" configuration section specifies defaults for this Analysis or sub-class
        # i.e. "bsf.Analysis" or "bsf.analyses.*", respectively.

        analysis.set_configuration(
            configuration=analysis.configuration,
            section=Configuration.section_from_instance(instance=analysis))

        return analysis

    def __init__(self, configuration=None,
                 project_name=None, genome_version=None,
                 input_directory=None, output_directory=None,
                 project_directory=None, genome_directory=None,
                 sas_file=None, sas_prefix=None, e_mail=None, debug=0, drms_list=None,
                 runnable_dict=None, collection=None, comparisons=None, samples=None):
        """Initialise an C{Analysis} object.

        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{Analysis}-wide project directory,
            normally under the C{Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{Analysis}-wide genome directory,
            normally under the C{Analysis}-wide project directory
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
        @param drms_list: Python C{list} of C{DRMS} objects
        @type drms_list: list[DRMS]
        @param runnable_dict: Python C{dict} of Python C{str} (C{Runnable.name}) and C{Runnable} value data
        @type runnable_dict: dict[Runnable.name, Runnable]
        @param collection: C{Collection}
        @type collection: Collection
        @param comparisons: Python C{dict} of Analysis-specific objects
            (i.e. Python tuple for RNA-Seq and ChIPSeqComparison for ChIPSeq)
        @type comparisons: dict[str, Any]
        @param samples: Python C{list} of C{Sample} objects
        @type samples: list[Sample]
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

        if drms_list is None:
            self.drms_list = list()
        else:
            self.drms_list = drms_list

        if runnable_dict is None:
            self.runnable_dict = dict()
        else:
            self.runnable_dict = runnable_dict

        if collection is None:
            self.collection = Collection()
        else:
            assert isinstance(collection, Collection)
            self.collection = collection

        if comparisons is None:
            self.comparisons = dict()
        else:
            self.comparisons = comparisons

        if samples is None:
            self.samples = list()
        else:
            self.samples = samples

        return

    def trace(self, level):
        """Trace an C{Analysis} object.

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

        output += '{}  Python dict of Runnable objects:\n'.format(indent)
        keys = self.runnable_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))
        for key in keys:
            assert isinstance(key, str)
            output += '{}    Key: {!r} Runnable: {!r}\n'.format(indent, key, self.runnable_dict[key])
            runnable = self.runnable_dict[key]
            assert isinstance(runnable, Runnable)
            output += runnable.trace(level=level + 2)

        output += '{}  Python List of Sample objects:\n'.format(indent)
        for sample in self.samples:
            assert isinstance(sample, Sample)
            output += '{}    Sample name: {!r} file_path: {!r}\n'.format(indent, sample.name, sample.file_path)

        if self.collection:
            output += self.collection.trace(level + 1)

        return output

    def add_drms(self, drms):
        """Convenience method to facilitate initialising, adding and returning a C{DRMS} object.

        @param drms: C{DRMS}
        @type drms: DRMS
        @return: C{DRMS}
        @rtype: DRMS
        """
        assert isinstance(drms, DRMS)

        self.drms_list.append(drms)

        return drms

    def add_runnable(self, runnable):
        """Convenience method to facilitate initialising, adding and returning a C{Runnable}.

        @param runnable: C{Runnable}
        @type runnable: Runnable
        @return: C{Runnable}
        @rtype: Runnable
        @raise Exception: The C{Runnable.name} already exists in the C{Analysis}
        """

        assert isinstance(runnable, Runnable)

        if runnable.name in self.runnable_dict:
            raise Exception("A Runnable object with name {!r} already exists in Analysis {!r}".
                            format(runnable.name, self.project_name))
        else:
            self.runnable_dict[runnable.name] = runnable

        return runnable

    def add_sample(self, sample):
        """Add a C{Sample} object to the Python C{list} of C{Sample} objects if it does not already exist.

        The check is based on the Python 'in' comparison operator and in lack of a specific
        __cmp__ method, relies on object identity (i.e. address).
        @param sample: C{Sample}
        @type sample: Sample
        @return:
        @rtype:
        """

        assert isinstance(sample, Sample)

        if sample not in self.samples:
            self.samples.append(sample)

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of an C{Analysis} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: Configuration
        @type configuration: Configuration
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
                'Section {!r} not defined in Configuration file {!r}.'.
                format(section, configuration.config_path))

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

        return

    def run(self):
        """Run the C{Analysis}.

        @raise Exception: An C{Analysis.project_name} has not been defined
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

        self.input_directory = Default.get_absolute_path(
            file_path=self.input_directory,
            default_path=Default.absolute_samples())

        self.output_directory = Default.get_absolute_path(
            file_path=self.output_directory,
            default_path=Default.absolute_projects())

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

            self.collection = Collection.from_sas_path(
                file_path=self.input_directory,
                file_type='Automatic',
                name=self.project_name,
                sas_path=self.sas_file,
                sas_prefix=self.sas_prefix)

            if self.debug > 1:
                print '{!r} Collection name: {!r}'.format(self, self.collection.name)
                print self.collection.trace(1)
        else:
            # Create an empty Collection.
            self.collection = Collection()

        return

    def report(self):
        """Create an C{Analysis} report.

        The method must be implemented in a sub-class.
        @return:
        @rtype:
        """

        warnings.warn(
            "The 'report' method must be implemented in the sub-class.",
            UserWarning)

        return

    def create_project_genome_directory(self):
        """Check and create an C{Analysis.project_directory} or C{Analysis.genome_directory} if necessary.

        @return:
        @rtype:
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

        return

    def create_public_project_link(self, sub_directory=None):
        """Create a symbolic link from the web directory to the project directory if not already there.

        The link will be placed in the sub directory and contain
        the project name followed by a 128 bit hexadecimal UUID string.
        @param sub_directory: C{Analysis}-specific directory
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

        link_name = os.path.join(html_path, '_'.join((self.project_name, uuid.uuid4().hex)))

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
        """Write a UCSC Track Hub I{hub.txt} file into the C{Analysis.project_directory},
        above the C{Analysis.genome_directory}.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        @return:
        @rtype:
        """

        output = str()

        if prefix is None or not prefix:
            file_name = 'hub.txt'
            output += 'hub {}\n'.format(self.project_name)
            output += 'shortLabel {}\n'.format(self.project_name)
            output += 'longLabel Project {}\n'.format(self.project_name)
            output += 'genomesFile genomes.txt\n'
        else:
            file_name = '{}_hub.txt'.format(prefix)
            output += 'hub {}_{}\n'.format(self.project_name, prefix)
            output += 'shortLabel {}_{}\n'.format(self.project_name, prefix)
            output += 'longLabel Project {}_{}\n'.format(self.project_name, prefix)
            output += 'genomesFile {}_genomes.txt\n'.format(prefix)

        output += 'email {}\n'.format(self.e_mail)

        # The [prefix_]hub.txt goes into the project directory above the genome directory.
        file_path = os.path.join(self.project_directory, file_name)

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        return

    def ucsc_hub_write_genomes(self, prefix=None):
        """Write a UCSC Track Hub I{genomes.txt} file into the C{Analysis.project_directory},
        above the C{Analysis.genome_directory}.

        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        @return:
        @rtype:
        """

        output = str()

        output += 'genome {}\n'.format(self.genome_version)
        if prefix is None or not prefix:
            file_name = 'genomes.txt'
            output += 'trackDb {}/trackDB.txt\n'.format(self.genome_version)
        else:
            file_name = '{}_genomes.txt'.format(prefix)
            output += 'trackDb {}/{}_trackDB.txt\n'.format(self.genome_version, prefix)

        # The [prefix_]genomes.txt goes into the project directory above the genome directory.
        file_path = os.path.join(self.project_directory, file_name)

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        return

    def ucsc_hub_write_tracks(self, output, prefix=None):
        """Write a UCSC Track Hub I{trackDB.txt} file into the C{Analysis.genome_directory}.

        @param output: Content
        @type output: str
        @param prefix: A hub prefix (e.g. chipseq, rnaseq, ...)
        @type prefix: str
        @return:
        @rtype:
        """

        if prefix is None or not prefix:
            file_name = 'trackDB.txt'
        else:
            file_name = '{}_trackDB.txt'.format(prefix)

        # The [prefix_]trackDB.txt goes into the genome directory under the project directory.
        file_path = os.path.join(self.genome_directory, file_name)

        file_handle = open(name=file_path, mode='w')
        file_handle.write(output)
        file_handle.close()

        return

    def submit(self, drms_name=None):
        """Submit each C{DRMS} object and pickle each C{Runnable} object.

        @param drms_name: Only submit C{Executable} objects linked to C{DRMS.name}
        @type drms_name: str
        @return:
        @rtype:
        """

        # Pickle all Runnable objects.

        for key in self.runnable_dict.keys():
            assert isinstance(key, str)
            self.runnable_dict[key].to_pickler_path()

        # Submit all Executable objects of all Distributed Resource Management System objects.

        submit = 0

        for drms in self.drms_list:
            assert isinstance(drms, DRMS)
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

        return


class Configuration(object):
    """The C{Configuration} class represents a UNIX-style initialisation (*.ini) file and
    an associated Python C{SafeConfigParser} object to parse the file.

    Attributes:
    @ivar config_path: C{Configuration} file path
    @type config_path: str | unicode
    @ivar config_parser: Python C{SafeConfigParser}
    @type config_parser: SafeConfigParser
    """

    @staticmethod
    def section_from_instance(instance):
        """Get a configuration section string composed of the Python module and Python class name from a
        Python instance.

        @param instance: A Python instance (i.e. object)
        @type instance: object
        @return: Configuration section string
        @rtype: str
        """

        # For Python "type" instances the "__name__" instance variable provides the name of the type, while for
        # Python "object" instances the "__class__" variable provides the Python "type" object.

        if isinstance(instance, type):
            return '.'.join((instance.__module__, instance.__name__))
        else:
            return '.'.join((instance.__module__, instance.__class__.__name__))

    @classmethod
    def from_config_path(cls, config_path):
        """Create a new C{Configuration} object based on a configuration file path.

        Both, user and variable expansion gets applied to the file path.
        @param config_path: Configuration file path
        @type config_path: str | unicode
        @return: C{Configuration}
        @rtype: Configuration
        @raise Exception: Configuration file path does not exist
        """

        assert isinstance(config_path, (str, unicode))

        config_path = os.path.expanduser(path=config_path)
        config_path = os.path.expandvars(path=config_path)

        # Since ConfigParser options are used as command line options,
        # they have to be case sensitive.
        # Hence, override optionxform() with str().

        config_parser = SafeConfigParser()
        config_parser.optionxform = str

        configuration = cls(config_path=config_path, config_parser=config_parser)

        files = configuration.config_parser.read(configuration.config_path)

        if len(files) == 0:
            raise Exception(
                'Configuration file {!r} does not exist.'.format(configuration.config_path))

        return configuration

    def __init__(self, config_path=None, config_parser=None):
        """Initialise a C{Configuration} object.

        @param config_path: Configuration file path
        @type config_path: str | unicode
        @param config_parser: Python C{SafeConfigParser}
        @type config_parser: SafeConfigParser
        @return:
        @rtype:
        """

        super(Configuration, self).__init__()

        if config_path is None:
            self.config_path = str()
        else:
            self.config_path = config_path

        if config_parser is None:
            self.config_parser = SafeConfigParser()
        else:
            assert isinstance(config_parser, SafeConfigParser)
            self.config_parser = config_parser

        return

    def trace(self, level):
        """Trace a C{Configuration} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  config_path:   {!r}\n'.format(indent, self.config_path)
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
    """The C{Default} class specifies the application or library default configuration.

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
    @ivar indices: Python C{dict} of program name key and index directory name value data
    @type indices: dict[str, str]
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
    def get_absolute_path(file_path, default_path=None):
        """Return an absolute file path.

        Expand an eventual user part i.e. on UNIX ~ or ~user and
        expand any environment variables i.e. on UNIX ${NAME} or $NAME
        Check if an absolute path has been provided, if not,
        automatically prepend default directory paths.
        Finally, normalise the path.

        @param file_path: File path
        @type file_path: str | unicode
        @param default_path: Default absolute path
        @type default_path: str | unicode
        @return: Absolute path
        @rtype: str | unicode
        """

        absolute_path = os.path.expanduser(path=file_path)
        absolute_path = os.path.expandvars(path=absolute_path)

        if default_path and not os.path.isabs(absolute_path):
            absolute_path = os.path.join(default_path, absolute_path)

        return os.path.normpath(absolute_path)

    @staticmethod
    def get_global_default():
        """Get the global Default configuration and initialise it, if not already done so.

        @return: Default
        @rtype: Default
        """

        if Default.global_default is None:
            Default.global_default = Default.from_global_file_path()

        return Default.global_default

    @classmethod
    def from_global_file_path(cls):
        """Create a new Default object from the global default configuration file.

        The default configuration is based on the file $HOME/.bsfpython.ini in the user's home directory.
        @return: Default
        @rtype: Default
        """

        return cls.from_config_path(config_path=Default.global_file_path)

    @classmethod
    def from_config_path(cls, config_path):
        """Create a new Default object from a UNIX-style configuration file.

        @param config_path: UNIX-style configuration file path
        @type config_path: str | unicode
        @return: Default
        @rtype: Default
        """

        return cls.from_configuration(configuration=Configuration.from_config_path(config_path=config_path))

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new Default objects from a Configuration object.

        @param configuration: Configuration
        @type configuration: Configuration
        @return: Default
        @rtype: Default
        """

        assert isinstance(configuration, Configuration)

        default = cls()

        default.set_configuration(configuration=configuration)

        return default

    def __init__(self, classpath_gatk=None, classpath_illumina2bam=None, classpath_picard=None, classpath_snpeff=None,
                 directory_home=None, directory_runs_illumina=None, directory_sequences=None, directory_samples=None,
                 directory_projects=None, directory_public_html=None, directory_genomes=None,
                 directory_annotations=None, directory_gatk_bundle=None, directory_intervals=None,
                 directory_snpeff_data=None,
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
        @param directory_intervals: Directory for interval list files
        @type directory_intervals: str | unicode
        @param directory_snpeff_data: snpEff database directory
        @type directory_snpeff_data: str | unicode
        @param indices: Python C{dict} of program name key and index directory name value data
        @type indices: dict[str, str]
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
        @return:
        @rtype:
        """

        super(Default, self).__init__()

        # Set Java class path information.

        if classpath_gatk is None:
            self.classpath_gatk = str()
        else:
            self.classpath_gatk = classpath_gatk

        if classpath_illumina2bam is None:
            self.classpath_illumina2bam = str()
        else:
            self.classpath_illumina2bam = classpath_illumina2bam

        if classpath_picard is None:
            self.classpath_picard = str()
        else:
            self.classpath_picard = classpath_picard

        if classpath_snpeff is None:
            self.classpath_snpeff = str()
        else:
            self.classpath_snpeff = classpath_snpeff

        # Set directory information.

        if directory_home is None:
            self.directory_home = str()
        else:
            self.directory_home = directory_home

        if directory_runs_illumina is None:
            self.directory_runs_illumina = str()
        else:
            self.directory_runs_illumina = directory_runs_illumina

        if directory_sequences is None:
            self.directory_sequences = str()
        else:
            self.directory_sequences = directory_sequences

        if directory_samples is None:
            self.directory_samples = str()
        else:
            self.directory_samples = directory_samples

        if directory_projects is None:
            self.directory_projects = str()
        else:
            self.directory_projects = directory_projects

        if directory_public_html is None:
            self.directory_public_html = str()
        else:
            self.directory_public_html = directory_public_html

        if directory_genomes is None:
            self.directory_genomes = str()
        else:
            self.directory_genomes = directory_genomes

        if directory_annotations is None:
            self.directory_annotations = str()
        else:
            self.directory_annotations = directory_annotations

        if directory_gatk_bundle is None:
            self.directory_gatk_bundle = str()
        else:
            self.directory_gatk_bundle = directory_gatk_bundle

        if directory_intervals is None:
            self.directory_intervals = str()
        else:
            self.directory_intervals = directory_intervals

        if directory_snpeff_data is None:
            self.directory_snpeff_data = str()
        else:
            self.directory_snpeff_data = directory_snpeff_data

        # Set index information.

        if indices is None:
            self.indices = dict()
        else:
            self.indices = indices

        # Set DRMS information.

        if drms_implementation is None:
            self.drms_implementation = str()
        else:
            self.drms_implementation = drms_implementation

        if drms_maximum_threads is None:
            self.drms_maximum_threads = str()
        else:
            self.drms_maximum_threads = drms_maximum_threads

        if drms_memory_limit_hard is None:
            self.drms_memory_limit_hard = str()
        else:
            self.drms_memory_limit_hard = drms_memory_limit_hard

        if drms_memory_limit_soft is None:
            self.drms_memory_limit_soft = str()
        else:
            self.drms_memory_limit_soft = drms_memory_limit_soft

        if drms_time_limit is None:
            self.drms_time_limit = str()
        else:
            self.drms_time_limit = drms_time_limit

        if drms_parallel_environment is None:
            self.drms_parallel_environment = str()
        else:
            self.drms_parallel_environment = drms_parallel_environment

        if drms_queue is None:
            self.drms_queue = str()
        else:
            self.drms_queue = drms_queue

        # Set operator information.

        if operator_e_mail is None:
            self.operator_e_mail = str()
        else:
            self.operator_e_mail = operator_e_mail

        if operator_sequencing_centre is None:
            self.operator_sequencing_centre = str()
        else:
            self.operator_sequencing_centre = operator_sequencing_centre

        # Set UCSC Genome Browser information.

        if ucsc_host_name is None:
            self.ucsc_host_name = str()
        else:
            self.ucsc_host_name = ucsc_host_name

        # Set URL information.

        if url_protocol is None:
            self.url_protocol = str()
        else:
            self.url_protocol = url_protocol

        if url_host_name is None:
            self.url_host_name = str()
        else:
            self.url_host_name = url_host_name

        if url_relative_projects is None:
            self.url_relative_projects = str()
        else:
            self.url_relative_projects = url_relative_projects

        return

    def set_configuration(self, configuration):
        """Set instance variables of a Default object via a section of a Configuration object.

        For each instance variable a configuration option has to be present.
        @param configuration: Configuration
        @type configuration: Configuration
        @return:
        @rtype:
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
        self.directory_intervals = cp.get(section=section, option='intervals')
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

        return

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
        @return: Absolute path to the GATK bundle directory
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
    def absolute_intervals():
        """Get the absolute directory path for interval list files.

        @return: Absolute path to the interval list directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_intervals):
            return default.directory_intervals
        else:
            return os.path.join(default.directory_home, default.directory_intervals)

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

        if genome_index not in default.indices:
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

        return '/'.join((default.url_absolute_base(), default.url_relative_projects))


class DRMS(object):
    """The I{Distributed Resource Management System} (C{DRMS}) class represents a
    I{Distributed Resource Management System} or batch job scheduler.

    Attributes:
    @ivar name: Name
    @type name: str
    @ivar working_directory: Working directory path
    @type working_directory: str
    @ivar implementation: Implementation (e.g. I{sge}, I{slurm}, ...)
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
    @ivar is_script: C{Executable} objects represent shell scripts,
        or alternatively binary programs
    @type is_script: bool
    @ivar executables: Python C{list} of C{Executable} objects
    @type executables: list[Executable]
    """

    @classmethod
    def from_analysis(cls, name, working_directory, analysis):
        """Create a C{DRMS} object from an C{Analysis} object.

        @param name: Name
        @type name: str
        @param working_directory: Working directory
        @type working_directory: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        @return: C{DRMS} object
        @rtype: DRMS
        """

        assert isinstance(analysis, Analysis)

        drms = cls(name=name, working_directory=working_directory)

        # Set a minimal set of global defaults.

        drms.set_default(default=Default.get_global_default())

        # A "bsf.DRMS" section specifies defaults for all DRMS objects of an Analysis.

        section = Configuration.section_from_instance(instance=drms)
        drms.set_configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 1:
            print 'DRMS configuration section: {!r}.'.format(section)

        # A "bsf.Analysis.DRMS" or "bsf.analyses.*.DRMS" pseudo-class section specifies
        # Analysis-specific or sub-class-specific options for the DRMS, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=analysis), 'DRMS'))
        drms.set_configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 1:
            print 'DRMS configuration section: {!r}.'.format(section)

        # A "bsf.Analysis.DRMS.name" or "bsf.analyses.*.DRMS.name" section specifies defaults
        # for a particular DRMS object of an Analysis or sub-class, respectively.

        section = '.'.join((Configuration.section_from_instance(instance=analysis), 'DRMS', drms.name))
        drms.set_configuration(configuration=analysis.configuration, section=section)

        if analysis.debug > 1:
            print 'DRMS configuration section: {!r}.'.format(section)

        return drms

    @classmethod
    def from_configuration(cls, name, work_directory, configuration, section):
        """Create a C{DRMS} object from a C{Configuration} object.

        @param name: Name
        @type name: str
        @param work_directory: Work directory
        @type work_directory: str
        @param configuration: C{Configuration} object
        @type configuration: Configuration
        @param section: Configuration section string
        @type section: str
        @return: C{DRMS} object
        @rtype: DRMS
        """

        assert isinstance(configuration, Configuration)

        drms = cls(name=name, work_directory=work_directory)

        # Set a minimal set of global defaults before setting the Configuration.

        drms.set_default(default=Default.get_global_default())
        drms.set_configuration(configuration=configuration, section=section)

        return drms

    def __init__(self, name, working_directory, implementation=None, memory_free_mem=None, memory_free_swap=None,
                 memory_free_virtual=None, memory_limit_hard=None, memory_limit_soft=None, time_limit=None,
                 parallel_environment=None, queue=None, threads=1, hold=None, is_script=False, executables=None):
        """Initialise a C{DRMS} object.

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
        @param is_script: C{Executable} objects represent shell scripts,
            or alternatively binary programs
        @type is_script: bool
        @param executables: Python C{list} of C{Executable} objects
        @type executables: list[Executable]
        @return:
        @rtype:
        """

        super(DRMS, self).__init__()

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

        if executables is None:
            self.executables = list()
        else:
            self.executables = executables

        return

    def trace(self, level):
        """Trace a C{DRMS} object.

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
        output += '{}  working_directory:       {!r}\n'. \
            format(indent, self.working_directory)
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
            assert isinstance(executable, Executable)
            output += executable.trace(level=level + 2)

        return output

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{DRMS} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            raise Exception(
                'Section {!r} not defined in Configuration file {!r}.'.
                format(section, configuration.config_path))

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

        return

    def set_default(self, default):
        """Set instance variables of a C{DRMS} object via a C{Default} object.

        @param default: C{Default} object
        @type default: Default
        @return:
        @rtype:
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

        return

    def add_executable(self, executable):
        """Convenience method to facilitate initialising, adding and returning an C{Executable}.

        @param executable: C{Executable}
        @type executable: Executable
        @return: C{Executable}
        @rtype: Executable
        """

        assert isinstance(executable, Executable)

        self.executables.append(executable)

        return executable

    def submit(self, debug=0):
        """Submit a command line for each C{Executable} object.

        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        # Dynamically import the module specific for the configured DRMS implementation.

        module = importlib.import_module('.'.join((__name__, 'drms', self.implementation)))

        module.submit(drms=self, debug=debug)

        return


class Command(object):
    """C{Command} class representing one program, its options and arguments and possibly
    another subordinate C{Command}.

    Attributes:
    @ivar program: Program
    @type program: str
    @ivar options: Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
        C{Argument} objects
    @type options: dict[Argument.key, list[Argument]]
    @ivar arguments: Python C{list} of Python C{str} (program argument) objects
    @type arguments: list[str]
    @ivar sub_command: Subordinate C{Command}
    @type sub_command: Command
    """

    def __init__(self, program=None, options=None, arguments=None, sub_command=None):
        """Initialise a C{Command} object.

        @param program: Program
        @type program: str
        @param options: Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str]
        @param sub_command: Subordinate C{Command} object
        @type sub_command: Command
        @return:
        @rtype:
        """

        super(Command, self).__init__()

        self.program = program  # Can be None.

        if options is None:
            self.options = dict()
        else:
            self.options = options

        if arguments is None:
            self.arguments = list()
        else:
            self.arguments = arguments

        self.sub_command = sub_command  # Can be None.

        return

    def trace(self, level):
        """Trace a C{Command} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  program:            {!r}\n'. \
            format(indent, self.program)

        # List all options

        output += '{}  options:\n'.format(indent)

        for key in self.options.keys():
            assert isinstance(key, str)
            output += '{}    key: {!r} Argument objects:\n'.format(indent, key)
            for argument in self.options[key]:
                assert isinstance(argument, Argument)
                output += argument.trace(level=level + 2)

        # List all arguments

        output += '{}  arguments:\n'.format(indent)

        i = 0
        for argument in self.arguments:
            assert isinstance(argument, str)
            output += '{}    {:2d}: {!r}\n'.format(indent, i, argument)
            i += 1

        if self.sub_command:
            output += self.sub_command.trace(level=level + 1)

        return output

    def add_argument(self, argument, override):
        """Add an C{Argument} or one of its sub-classes.

        @param argument: C{Argument}
        @type argument: Argument
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        assert isinstance(argument, Argument)
        assert isinstance(override, bool)

        if not override and argument.key in self.options:
            warnings.warn(
                'Adding an Argument with key {!r} that exits already in Command.program {!r}.'.
                format(argument.key, self.program),
                UserWarning)

        if argument.key in self.options:
            arguments_list = self.options[argument.key]
            assert isinstance(arguments_list, list)
        else:
            arguments_list = list()
            self.options[argument.key] = arguments_list

        arguments_list.append(argument)

        return

    def add_switch_long(self, key, override=False):
        """Initialise and add a C{SwitchLong} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=SwitchLong(key=key), override=override)

    def add_switch_short(self, key, override=False):
        """Initialise and add a C{SwitchShort} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=SwitchShort(key=key), override=override)

    def add_option_long(self, key, value, override=False):
        """Initialise and add an C{OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=OptionLong(key=key, value=value), override=override)

    def add_option_short(self, key, value, override=False):
        """Initialise and add an C{OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=OptionShort(key=key, value=value), override=override)

    def add_option_pair(self, key, value, override=False):
        """Initialise and add an C{OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.add_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_argument(self, argument, override):
        """Set an C{Argument} or one of its sub-classes.

        @param argument: C{Argument}
        @type argument: Argument
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        assert isinstance(argument, Argument)
        assert isinstance(override, bool)

        if not override and argument.key in self.options:
            warnings.warn(
                'Setting an Argument with key {!r} that exits already in Command.program {!r}.'.
                format(argument.key, self.program),
                UserWarning)

        self.options[argument.key] = [argument]

        return

    def set_switch_long(self, key, override=False):
        """Initialise and set a C{SwitchLong} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=SwitchLong(key=key), override=override)

    def set_switch_short(self, key, override=False):
        """Initialise and set a C{SwitchShort} object.

        @param key: Key
        @type key: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=SwitchShort(key=key), override=override)

    def set_option_long(self, key, value, override=False):
        """Initialise and set an C{OptionLong} object.

        @param key: Key
        @type key: str
        @param value: Value
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @type value: str | unicode
        @return:
        @rtype:
        """
        return self.set_argument(argument=OptionLong(key=key, value=value), override=override)

    def set_option_short(self, key, value, override=False):
        """Initialise and set an C{OptionShort} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=OptionShort(key=key, value=value), override=override)

    def set_option_pair(self, key, value, override=False):
        """Initialise and set an C{OptionPair} object.

        @param key: Key
        @type key: str
        @param value: Value
        @type value: str | unicode
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """
        return self.set_argument(argument=OptionPair(key=key, value=value), override=override)

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{Command} object via a section of a C{Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration file section, defaults to instance class
        @type section: str
        @return:
        @rtype:
        """

        assert isinstance(configuration, Configuration)
        assert isinstance(section, str)

        if not configuration.config_parser.has_section(section=section):
            warnings.warn(
                'Section {!r} not defined in Configuration file {!r}.'.
                format(section, configuration.config_path),
                UserWarning)

            return

        # The configuration section is available.

        for option in configuration.config_parser.options(section=section):
            self.add_argument(
                argument=Argument.from_key_value(
                    key=option,
                    value=configuration.config_parser.get(
                        section=section,
                        option=option)),
                override=False)

        return

    def command_list(self):
        """Assemble the command line from program, options and arguments.

        @return: Python C{list} of program, options, switches and arguments
        @rtype: list[str | unicode]
        """

        command_line = list()

        if self.program:
            command_line.append(self.program)

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            assert isinstance(key, str)
            options_list = self.options[key]
            assert isinstance(options_list, list)
            for argument in options_list:
                assert isinstance(argument, Argument)
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
            assert isinstance(argument, str)
            command_line.append(str(argument))

        # Expand a subordinate command, if defined.

        if self.sub_command:
            command_line.extend(self.sub_command.command_list())

        return command_line

    def command_str(self):
        """Assemble the command line from program, options, switches and arguments.

        @return: A Python C{str} of program, options, switches and arguments
        @rtype: str
        """

        command_line = str()

        if self.program:
            command_line += self.program

        # Add all options and switches in alphabetical order.

        keys = self.options.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))

        for key in keys:
            assert isinstance(key, str)
            options_list = self.options[key]
            assert isinstance(options_list, list)
            for argument in options_list:
                assert isinstance(argument, Argument)
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
            assert isinstance(argument, str)
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
    @ivar options: Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
        C{Argument} objects
    @type options: dict[Argument.key, list[Argument]]
    @ivar arguments: Python C{list} of Python C{str} or C{unicode} (argument) objects
    @type arguments: list[str | unicode]
    @ivar sub_command: Subordinate Command
    @type sub_command: Command
    @ivar stdout_path: Standard output (STDOUT) redirection in Bash (1>word)
    @type stdout_path: str | unicode
    @ivar stderr_path: Standard error (STDERR) redirection in Bash (2>word)
    @type stderr_path: str | unicode
    @ivar dependencies: Python C{list} of C{Executable.name} properties in the
        context of C{DRMS} dependencies
    @type dependencies: list[Executable.name]
    @ivar hold: Hold on job scheduling
    @type hold: str
    @ivar submit: Submit the Executable into the DRMS
    @type submit: bool
    @ivar maximum_attempts: Maximum number of attempts to run this C{Executable}
    @type maximum_attempts: int
    @ivar process_identifier: Process identifier
    @type process_identifier: str
    @ivar process_name: Process name
    @type process_name: str
    """

    @staticmethod
    def process_stream(file_type, file_handle, thread_lock, file_path=None, debug=0):
        """C{Executable} function to process I{STDOUT} or I{STDERR} from the child process as a thread.

        @param file_type: File handle type I{STDOUT} or I{STDERR}
        @type file_type: str
        @param file_handle: The I{STDOUT} or I{STDERR} file handle
        @type file_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: thread.lock
        @param file_path: I{STDOUT} file path
        @type file_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        @raise Exception: The file_type has to be either I{STDOUT} or I{STDERR}
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

        return

    @staticmethod
    def process_stdout(stdout_handle, thread_lock, stdout_path=None, debug=0):
        """C{Executable} function to process I{STDOUT} from the child process as a thread.

        @param stdout_handle: The I{STDOUT} file handle
        @type stdout_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: thread.lock
        @param stdout_path: I{STDOUT} file path
        @type stdout_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        return Executable.process_stream(
            file_type='STDOUT',
            file_handle=stdout_handle,
            thread_lock=thread_lock,
            file_path=stdout_path,
            debug=debug)

    @staticmethod
    def process_stderr(stderr_handle, thread_lock, stderr_path=None, debug=0):
        """C{Executable} function to process I{STDERR} from the child process as a thread.

        @param stderr_handle: The I{STDERR} file handle
        @type stderr_handle: file
        @param thread_lock: A Python C{threading.Lock} object
        @type thread_lock: thread.lock
        @param stderr_path: I{STDERR} file path
        @type stderr_path: str | unicode
        @param debug: Debug level
        @type debug: int
        @return:
        @rtype:
        """

        return Executable.process_stream(
            file_type='STDERR',
            file_handle=stderr_handle,
            thread_lock=thread_lock,
            file_path=stderr_path,
            debug=debug)

    @classmethod
    def from_analysis(cls, name, program, analysis):
        """Create an C{Executable} object from an C{Analysis} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        @return: C{Executable}
        @rtype: Executable
        """

        assert isinstance(analysis, Analysis)

        # Initialise an Executable object with default values.

        executable = cls(name=name, program=program)

        section = Configuration.section_from_instance(executable)

        # For plain Executable objects append the value of the
        # Executable.program to make this more meaningful.

        if section == 'bsf.Executable':
            section += '.'
            section += executable.program

        if analysis.debug > 1:
            print 'Executable configuration section: {!r}.'.format(section)

        executable.set_configuration(configuration=analysis.configuration, section=section)

        return executable

    @classmethod
    def from_analysis_runnable(cls, analysis, runnable_name):
        """Create an C{Executable} to submit a C{Runnable} into a C{DRMS}.

        In case C{Runnable.get_relative_status_path} exists already, C{Executable.submit} will be set to C{False}.
        @param analysis: C{Analysis}
        @type analysis: Analysis
        @param runnable_name: C{Runnable.name}
        @type runnable_name: str
        @return: Executable
        @rtype: Executable
        @raise Exception: A C{Runnable.name} does not exist in C{Analysis}
        """

        assert isinstance(analysis, Analysis)

        if runnable_name not in analysis.runnable_dict:
            raise Exception("A Runnable object with name {!r} does not exist in the Analysis object with name {!r}.".
                            format(runnable_name, analysis.project_name))

        runnable = analysis.runnable_dict[runnable_name]
        assert isinstance(runnable, Runnable)
        executable = cls(name=runnable.name, program=Runnable.runner_script)
        # TODO: Read configuration files for RunnableStep objects rather than Runnable objects.
        # Since bsf.Runnable.code_module objects such as 'bsf.runnables.generic' can be very generic,
        # it makes no sense to read standard configuration options from a Configuration object.
        # It would be better to read standard configuration options for RunnableStep objects.
        # executable.set_configuration(configuration=analysis.configuration, section=runnable.code_module)
        executable.add_option_long(key='pickler-path', value=runnable.pickler_path)

        # Only submit the Executable if the status file does not exist already.
        if os.path.exists(runnable.get_absolute_status_path):
            executable.submit = False

        return executable

    @classmethod
    def from_configuration(cls, name, program, configuration, section):
        """Create an C{Executable} object from a C{Configuration} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param configuration: C{Configuration}
        @type configuration: Configuration
        @param section: Configuration section string
        @type section: str
        @return: C{Executable}
        @rtype: Executable
        """

        assert isinstance(configuration, Configuration)

        executable = cls(name=name, program=program)

        executable.set_configuration(configuration=configuration, section=section)

        return executable

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, maximum_attempts=1, process_identifier=None, process_name=None):
        """Initialise an C{Executable} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param maximum_attempts: Maximum number of attempts to run this C{Executable}
        @type maximum_attempts: int
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @return:
        @rtype:
        """

        super(Executable, self).__init__(
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command)

        self.name = name  # Can be None.

        if stderr_path is None:
            self.stderr_path = str()
        else:
            self.stderr_path = stderr_path

        if stdout_path is None:
            self.stdout_path = str()
        else:
            self.stdout_path = stdout_path

        if dependencies is None:
            self.dependencies = list()
        else:
            self.dependencies = dependencies

        if hold is None:
            self.hold = str()
        else:
            self.hold = hold

        if submit is None:
            self.submit = True
        else:
            assert isinstance(submit, bool)
            self.submit = submit

        if maximum_attempts is None:
            self.maximum_attempts = int(x=1)
        else:
            self.maximum_attempts = maximum_attempts

        if process_identifier is None:
            self.process_identifier = str()
        else:
            self.process_identifier = process_identifier

        if process_name is None:
            self.process_name = str()
        else:
            self.process_name = process_name

        return

    def trace(self, level):
        """Trace an C{Executable} object.

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
        output += '{}  maximum_attempts: {!r}\n'. \
            format(indent, self.maximum_attempts)
        output += '{}  process_identifier: {!r}\n'. \
            format(indent, self.process_identifier)
        output += '{}  process_name:       {!r}\n'. \
            format(indent, self.process_name)

        # List all dependencies.

        output += '{}  dependencies:\n'.format(indent)

        i = 0
        for dependency in self.dependencies:
            assert isinstance(dependency, str)
            output += '{}    {:2d} {!r}\n'.format(indent, i, dependency)
            i += 1

        # Trace the Command super-class.

        output += super(Executable, self).trace(level=level + 1)

        return output

    def command_list(self):
        """Assemble the command line from program, options and arguments.

        @return: Python C{list} of program, options and arguments
        @rtype: list[str | unicode]
        """

        command = list()

        command.extend(super(Executable, self).command_list())

        # The stdout_path and stderr_path gets appended in specific modules.

        return command

    def command_str(self):
        """Assemble the command line from program, options, switches and arguments.

        @return: A Python C{str} of program, options, switches and arguments
        @rtype: str
        """

        command = str()

        command += super(Executable, self).command_str()

        # The stdout_path and stderr_path gets appended in specific modules.

        return command

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run an C{Executable} object via the Python C{subprocess.Popen} class.

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

        child_return_code = 0
        attempt_counter = 0

        while attempt_counter < self.maximum_attempts:

            child_process = Popen(
                args=self.command_list(),
                bufsize=0,
                stdin=PIPE,
                stdout=PIPE,
                stderr=PIPE,
                shell=False,
                close_fds=on_posix)

            # Two threads, thread_out and thread_err reading STDOUT and STDERR, respectively,
            # should make sure that buffers are not filling up.

            thread_lock = Lock()

            thread_out = Thread(
                target=Executable.process_stdout,
                kwargs=dict(
                    stdout_handle=child_process.stdout,
                    thread_lock=thread_lock,
                    stdout_path=self.stdout_path,
                    debug=debug))
            thread_out.daemon = True  # Thread dies with the program.
            thread_out.start()

            thread_err = Thread(
                target=Executable.process_stderr,
                kwargs=dict(
                    stderr_handle=child_process.stderr,
                    thread_lock=thread_lock,
                    stderr_path=self.stderr_path,
                    debug=debug))
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
                        format(datetime.datetime.now().isoformat(), self.name, +child_return_code)
                attempt_counter += 1
            elif child_return_code < 0:
                if debug > 0:
                    print '[{}] Child process {!r} received signal {}.'. \
                        format(datetime.datetime.now().isoformat(), self.name, -child_return_code)
            else:
                if debug > 0:
                    print '[{}] Child process {!r} completed successfully {}.'. \
                        format(datetime.datetime.now().isoformat(), self.name, +child_return_code)
                break

        else:
            if debug > 0:
                print '[{}] Runnable {!r} exceeded the maximum retry counter {}.' \
                    .format(datetime.datetime.now().isoformat(), self.name, self.maximum_attempts)

        return child_return_code

    def evaluate_return_code(self, return_code):
        """Evaluate a return code from the run method.

        @param return_code: Return code
        @type return_code: int
        @return:
        @rtype:
        """

        if return_code > 0:
            print '[{}] Child process {!r} failed with return code {}'. \
                format(datetime.datetime.now().isoformat(), self.name, +return_code)
        elif return_code < 0:
            print '[{}] Child process {!r} received signal {}.'. \
                format(datetime.datetime.now().isoformat(), self.name, -return_code)
        else:
            print '[{}] Child process {!r} completed with return code {}.'. \
                format(datetime.datetime.now().isoformat(), self.name, +return_code)

        return


class RunnableStep(Executable):
    """The C{RunnableStep} represents a step in a C{Runnable} class.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{RunnableStep}
    @type obsolete_file_path_list: list[str | unicode]
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None):
        """Initialise a C{RunnableStep} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @return:
        @rtype:
        """

        super(RunnableStep, self).__init__(
            name=name, program=program, options=options, arguments=arguments,  sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold, submit=submit,
            process_identifier=process_identifier, process_name=process_name)

        if obsolete_file_path_list is None:
            self.obsolete_file_path_list = list()
        else:
            self.obsolete_file_path_list = obsolete_file_path_list

        return

    def trace(self, level=1):
        """Trace a C{RunnableStep} object.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  obsolete_file_path_list: {!r}\n'.format(indent, self.obsolete_file_path_list)
        output += super(RunnableStep, self).trace(level=level + 1)

        return output


class RunnableStepJava(RunnableStep):
    """The C{RunnableStepJava} class represents a C{RunnableStep} with all peculiarities of Java programs.

    Attributes:
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None,
                 java_temporary_path=None, java_heap_maximum=None, java_jar_path=None):
        """Create a C{RunnableStep} for a Java program.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @return: C{RunnableStep}
        @rtype: RunnableStep
        """

        super(RunnableStepJava, self).__init__(
            name=name,
            program=program, options=options, arguments=arguments, sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold,
            submit=submit, process_identifier=process_identifier, process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        # JavaVM command
        if self.program is None:
            self.program = 'java'

        # JavaVM options
        if 'd64' not in self.options:
            self.add_switch_short(key='d64')

        if 'server' not in self.options:
            self.add_switch_short(key='server')

        if java_heap_maximum and java_heap_maximum not in self.options:
            self.add_switch_short(key=java_heap_maximum)

        if '-Djava.io.tmpdir' not in self.options:
            self.add_option_pair(key='-Djava.io.tmpdir', value=java_temporary_path)

        if self.sub_command is None:
            # The Picard command line interface is a bit broken, as the -jar option needs to come last,
            # just before the Picard command. GATK does this slightly better with the --analysis_type option.
            # To be on the safe side, an empty sub command is required to separate the -jar option from the
            # other JavaVM options.
            self.sub_command = Command()
            if java_jar_path is not None:
                self.sub_command.add_option_short(key='jar', value=java_jar_path)

        return


class RunnableStepPicard(RunnableStepJava):
    """The C{RunnableStepPicard} class represents a C{RunnableStepJava} specific to Picard tools.

    Attributes:
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None,
                 java_temporary_path=None, java_heap_maximum=None, java_jar_path=None,
                 picard_classpath=None, picard_command=None):
        """Create a C{RunnableStep} for a Picard algorithm.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate Command
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param java_temporary_path: Temporary directory path for the Java Virtual Machine
        @type java_temporary_path: str | unicode
        @param java_heap_maximum: Java heap maximum size (-Xmx option)
        @type java_heap_maximum: str
        @param java_jar_path: Java archive file path
        @type java_jar_path: str | unicode
        @param picard_classpath: Picard class path
        @type picard_classpath: str | unicode
        @param picard_command: Picard command
        @type picard_command: str
        @return: C{RunnableStep}
        @rtype: RunnableStep
        """

        super(RunnableStepPicard, self).__init__(
            name=name,
            program=program, options=options, arguments=arguments, sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold,
            submit=submit, process_identifier=process_identifier, process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list,
            java_temporary_path=java_temporary_path, java_heap_maximum=java_heap_maximum, java_jar_path=java_jar_path)

        # Set the Picard classpath and the Picard Java archive.
        if 'jar' not in self.sub_command.options:
            self.sub_command.add_option_short(key='jar', value=os.path.join(picard_classpath, 'picard.jar'))

        # The Picard algorithm is then another sub-command.
        if self.sub_command.sub_command is None:
            self.sub_command.sub_command = Command(program=picard_command)

        return

    def add_picard_option(self, key, value, override=False):
        """Add an option to the Picard command.

        @param key: Option key
        @type key: str
        @param value: Option value
        @type value: str
        @param override: Override existing C{Argument} without warning
        @type override: bool
        @return:
        @rtype:
        """

        return self.sub_command.sub_command.add_option_pair(key=key, value=value, override=override)


class RunnableStepLink(RunnableStep):
    """The C{RunnableStepLink} represents a step in a C{Runnable} class.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{RunnableStep}
    @type obsolete_file_path_list: list[str | unicode]
    @ivar source_path: Source path
    @type source_path: str | unicode
    @ivar target_path: Target path
    @type target_path: str | unicode
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None, source_path=None, target_path=None):
        """Initialise a C{RunnableStepLink} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param source_path: Source path
        @type source_path: str | unicode
        @param target_path: Target path
        @type target_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepLink, self).__init__(
            name=name, program=program, options=options, arguments=arguments,  sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold, submit=submit,
            process_identifier=process_identifier, process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if source_path is None:
            self.source_path = str()
        else:
            self.source_path = source_path

        if target_path is None:
            self.target_path = str()
        else:
            self.target_path = target_path

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

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

        if self.source_path and self.target_path and not os.path.exists(self.target_path):
            try:
                os.symlink(self.source_path, self.target_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return 0


class RunnableStepMakeDirectory(RunnableStep):
    """The C{RunnableStepMakeDirectory} represents a step in a C{Runnable} class.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{RunnableStep}
    @type obsolete_file_path_list: list[str | unicode]
    @ivar directory_path: Directory path
    @type directory_path: str | unicode
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None, directory_path=None):
        """Initialise a C{RunnableStepMakeDirectory} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param directory_path: Directory path
        @type directory_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepMakeDirectory, self).__init__(
            name=name, program=program, options=options, arguments=arguments,  sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold, submit=submit,
            process_identifier=process_identifier, process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if directory_path is None:
            self.directory_path = str()
        else:
            self.directory_path = directory_path

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

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

        if self.directory_path and not os.path.isdir(self.directory_path):
            try:
                os.makedirs(self.directory_path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

        return 0


class RunnableStepMove(RunnableStep):
    """The C{RunnableStepMove} represents a step in a C{Runnable} class.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{RunnableStep}
    @type obsolete_file_path_list: list[str | unicode]
    @ivar source_path: Source path
    @type source_path: str | unicode
    @ivar target_path: Target path
    @type target_path: str | unicode
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None, source_path=None, target_path=None):
        """Initialise a C{RunnableStepMove} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param source_path: Source path
        @type source_path: str | unicode
        @param target_path: Target path
        @type target_path: str | unicode
        @return:
        @rtype:
        """

        super(RunnableStepMove, self).__init__(
            name=name, program=program, options=options, arguments=arguments,  sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold, submit=submit,
            process_identifier=process_identifier, process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if source_path is None:
            self.source_path = str()
        else:
            self.source_path = source_path

        if target_path is None:
            self.target_path = str()
        else:
            self.target_path = target_path

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

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

        if self.source_path and self.target_path:
            # os.rename(self.source_path, self.target_path)
            shutil.move(src=self.source_path, dst=self.target_path)

        return 0


class RunnableStepSleep(RunnableStep):
    """The C{RunnableStepSleep} represents a step in a C{Runnable} class.

    Attributes:
    @ivar obsolete_file_path_list: Python C{list} of file paths that can be removed
        after successfully completing this C{RunnableStep}
    @type obsolete_file_path_list: list[str | unicode]
    @ivar sleep_time: Sleep time in seconds
    @type sleep_time: float
    """

    def __init__(self, name,
                 program=None, options=None, arguments=None, sub_command=None,
                 stdout_path=None, stderr_path=None, dependencies=None, hold=None,
                 submit=True, process_identifier=None, process_name=None,
                 obsolete_file_path_list=None, sleep_time=None):
        """Initialise a C{RunnableStepSleep} object.

        @param name: Name
        @type name: str
        @param program: Program
        @type program: str
        @param options:  Python C{dict} of Python C{str} (C{Argument.key}) key and Python C{list} value objects of
            C{Argument} objects
        @type options: dict[Argument.key, list[Argument]]
        @param arguments: Python C{list} of program arguments
        @type arguments: list[str | unicode]
        @param sub_command: Subordinate C{Command}
        @type sub_command: Command
        @param stdout_path: Standard output (I{STDOUT}) redirection in Bash (1>word)
        @type stdout_path: str | unicode
        @param stderr_path: Standard error (I{STDERR}) redirection in Bash (2>word)
        @type stderr_path: str | unicode
        @param dependencies: Python C{list} of C{Executable.name}
            properties in the context of C{DRMS} dependencies
        @type dependencies: list[Executable.name]
        @param hold: Hold on job scheduling
        @type hold: str
        @param submit: Submit the C{Executable} into the C{DRMS}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str
        @param process_name: Process name
        @type process_name: str
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{RunnableStep}
        @type obsolete_file_path_list: list[str | unicode]
        @param sleep_time: Sleep time in seconds
        @type sleep_time: float
        @return:
        @rtype:
        """

        super(RunnableStepSleep, self).__init__(
            name=name, program=program, options=options, arguments=arguments,  sub_command=sub_command,
            stdout_path=stdout_path, stderr_path=stderr_path, dependencies=dependencies, hold=hold, submit=submit,
            process_identifier=process_identifier, process_name=process_name,
            obsolete_file_path_list=obsolete_file_path_list)

        if sleep_time is None:
            self.sleep_time = float()
        else:
            assert isinstance(sleep_time, float)
            self.sleep_time = sleep_time

        return

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{RunnableStepMakeDirectory} object.

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

        time.sleep(self.sleep_time)

        return 0


class Runnable(object):
    """The C{Runnable} class holds all information to run one or more C{Executable} objects through the
    I{Runner} script. It can be thought of a script that executes runnable steps.

    Attributes:
    @cvar runner_script: Name of the I{Runner} script
    @type runner_script: str | unicode
    @ivar name: Name
    @type name: str
    @ivar code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
        C{Executable} objects via the I{Runner} script.
    @type code_module: str
    @ivar executable_dict: Python C{dict} of Python C{str} (C{Executable.name}) key data and C{Executable} value data
    @type executable_dict: dict[Executable.name, Executable]
    @ivar file_path_dict: Python C{dict} of Python C{str} (name) key data and Python C{str} (file_path) value data
    @type file_path_dict: dict[str, str | unicode]
    @ivar runnable_step_list: Python C{list} of C{RunnableStep} objects
    @type runnable_step_list: list[RunnableStep]
    @ivar working_directory: Working directory to write C{pickle.Pickler} files
    @type working_directory: str | unicode
    @ivar debug: Debug level
    @type debug: int
    """

    runner_script = 'bsf_runner.py'

    def __init__(self, name, code_module, working_directory, file_path_dict=None, executable_dict=None,
                 runnable_step_list=None, debug=0):
        """Initialise a C{Runnable} object.

        @param name: Name
        @type name: str
        @param code_module: The name of a module, usually in C{bsf.runnables} that implements the logic required to run
            C{Executable} objects via the I{Runner} script
        @type code_module: str
        @param working_directory: Working directory for writing a Python C{pickle.Pickler} file
        @type working_directory: str | unicode
        @param file_path_dict: Python C{dict} of Python C{str} (name) key data and
            Python C{str} (file_path) value data
        @type file_path_dict: dict[Executable.name, Executable]
        @param executable_dict: Python C{dict} of Python C{str} (C{Executable.name}) key data and
            C{Executable} value data
        @type executable_dict: dict[str, str | unicode]
        @param runnable_step_list: Python C{list} of C{RunnableStep} objects
        @type runnable_step_list: list[RunnableStep]
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        super(Runnable, self).__init__()

        self.name = name  # Can be None.
        self.code_module = code_module  # Can be None.
        self.working_directory = working_directory  # Can be None.

        if file_path_dict is None:
            self.file_path_dict = dict()
        else:
            self.file_path_dict = file_path_dict

        if executable_dict is None:
            self.executable_dict = dict()
        else:
            self.executable_dict = executable_dict

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
        """Trace a C{Runnable} object.

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
        output += '{}  file_path_dict: {!r}\n'.format(indent, self.file_path_dict)
        output += '{}  executable_dict: {!r}\n'.format(indent, self.executable_dict)
        output += '{}  runnable_step_list: {!r}\n'.format(indent, self.runnable_step_list)
        output += '{}  debug: {!r}\n'.format(indent, self.debug)

        output += '{}  Python dict of Python str (file path) objects:\n'.format(indent)
        keys = self.file_path_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))
        for key in keys:
            assert isinstance(key, str)
            output += '{}    Key: {!r} file_path: {!r}\n'.format(indent, key, self.file_path_dict[key])

        output += '{}  Python dict of Executable objects:\n'.format(indent)
        keys = self.executable_dict.keys()
        keys.sort(cmp=lambda x, y: cmp(x, y))
        for key in keys:
            assert isinstance(key, str)
            output += '{}    Key: {!r} Executable: {!r}\n'.format(indent, key, self.executable_dict[key])
            executable = self.executable_dict[key]
            assert isinstance(executable, Executable)
            output += executable.trace(level=level + 2)

        output += '{}  Python list of RunnableStep objects:\n'.format(indent)
        for runnable_step in self.runnable_step_list:
            assert isinstance(runnable_step, RunnableStep)
            output += runnable_step.trace(level=level + 1)

        return output

    def add_executable(self, executable):
        """Add an C{Executable}.

        @param executable: C{Executable}
        @type executable: Executable
        @return: C{Executable}
        @rtype: Executable
        @raise Exception: An C{Executable.name} already exists in the C{Runnable}
        """

        assert isinstance(executable, Executable)

        if executable.name in self.executable_dict:
            raise Exception("An Executable object with name {!r} already exists in Runnable object {!r}.".
                            format(executable.name, self.name))
        else:
            self.executable_dict[executable.name] = executable

        return executable

    def add_runnable_step(self, runnable_step=None):
        """Convenience method to facilitate initialising, adding and retuning a C{RunnableStep}.

        @param runnable_step: C{RunnableStep}
        @type runnable_step: RunnableStep
        @return: C{RunnableStep}
        @rtype: RunnableStep
        """

        if runnable_step is None:
            return

        assert isinstance(runnable_step, RunnableStep)

        self.runnable_step_list.append(runnable_step)

        return runnable_step

    def run_executable(self, name):
        """Run an C{Executable} defined in the C{Runnable} object.

        @param name: C{Executable.name}
        @type name: str
        @return:
        @rtype:
        @raise Exception: Child process failed with return code or received a signal
        """

        executable = self.executable_dict[name]
        assert isinstance(executable, Executable)
        child_return_code = executable.run()

        if child_return_code > 0:
            raise Exception('[{}] Child process {!r} failed with return code {}'.
                            format(datetime.datetime.now().isoformat(), executable.name, +child_return_code))
        elif child_return_code < 0:
            raise Exception('[{}] Child process {!r} received signal {}.'.
                            format(datetime.datetime.now().isoformat(), executable.name, -child_return_code))
        else:
            return

    @property
    def pickler_path(self):
        """Get the Python C{pickle.Pickler} file path.

        @return: Python C{pickle.Pickler} file path
        @rtype: str | unicode
        """

        return os.path.join(self.working_directory, '.'.join((self.name, 'pkl')))

    def to_pickler_path(self):
        """Write this C{Runnable} object as a Python C{pickle.Pickler} file into the working directory.
        @return:
        @rtype:
        """

        pickler_file = open(self.pickler_path, 'wb')
        pickler = Pickler(file=pickler_file, protocol=HIGHEST_PROTOCOL)
        pickler.dump(obj=self)
        pickler_file.close()

        return

    @classmethod
    def from_pickler_path(cls, file_path):
        """Create a C{Runnable} object from a Python C{pickle.Pickler} file via Python C{pickle.Unpickler}.

        @param file_path: File path to a Python C{pickle.Pickler} file
        @type file_path: str | unicode
        @return: C{Runnable}
        @rtype: Runnable
        """

        pickler_file = open(file_path, 'rb')
        unpickler = Unpickler(file=pickler_file)
        runnable = unpickler.load()
        pickler_file.close()

        assert isinstance(runnable, Runnable)

        return runnable

    @property
    def get_relative_status_path(self):
        """Get the relative status file path indicating successful completion of this C{Runnable}.

        @return: Relative status file path
        @rtype: str
        """

        return '_'.join((self.name, 'completed.txt'))

    @property
    def get_relative_temporary_directory_path(self):
        """Get the relative temporary directory path for this C{Runnable}.

        @return: Relative temporary directory path
        @rtype: str
        """

        return '_'.join((self.name, 'temporary'))

    @property
    def get_absolute_status_path(self):
        """Get the absolute status file path including the C{Runnable.working_directory}.

        @return: Absolute status file path
        @rtype: str
        """

        return os.path.join(self.working_directory, self.get_relative_status_path)

    @property
    def get_absolute_temporary_directory_path(self):
        """Get the absolute temporary directory path including the  C{Runnable.working_directory}.

        @return: Absolute temporary directory path
        @rtype: str
        """

        return os.path.join(self.working_directory, self.get_relative_temporary_directory_path)
