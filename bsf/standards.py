"""bsf.standards

A package of classes and methods modelling configuration and default information.
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


import ConfigParser
import os


class Configuration(object):
    """The C{bsf.standards.Configuration} class represents one or more UNIX-style initialisation (*.ini) files.

    A C{bsf.standards.Configuration} has an associated Python C{ConfigParser.SafeConfigParser} to parse the file(s).

    Attributes:
    @cvar global_configuration: Global C{bsf.standards.Configuration}
    @type global_configuration: bsf.standards.Configuration
    @cvar global_file_path: Global configuration file
    @type global_file_path: str | unicode
    @ivar file_path_list: C{bsf.standards.Configuration} file path
    @type file_path_list: list[str | unicode]
    @ivar config_parser: Python C{ConfigParser.SafeConfigParser}
    @type config_parser: ConfigParser.SafeConfigParser
    """

    global_configuration = None

    global_file_path = '~/.bsfpython.ini'

    @staticmethod
    def get_global_configuration():
        """Get the global C{bsf.standards.Configuration} and initialise it, if not already done so.

        @return: C{bsf.standards.Configuration}
        @rtype: bsf.standards.Configuration
        """

        if Configuration.global_configuration is None:
            Configuration.global_configuration = Configuration.from_file_path_list(
                file_path_list=[Configuration.global_file_path])

        return Configuration.global_configuration

    @staticmethod
    def get_absolute_path(file_path=None, default_path=None):
        """Return an absolute file path.

        Expand an eventual user part (i.e. on UNIX ~ or ~user) and
        expand any environment variables (i.e. on UNIX ${NAME} or $NAME).
        Check if an absolute path has been provided, if not,
        automatically prepend default directory paths,
        which again has user part and environment variables expanded.
        Finally, normalise the file path.

        @param file_path: File path
        @type file_path: str | unicode
        @param default_path: Default absolute path
        @type default_path: str | unicode
        @return: Absolute path
        @rtype: None | str | unicode
        """
        if file_path:
            file_path = os.path.expanduser(file_path)
            file_path = os.path.expandvars(file_path)

            if os.path.isabs(file_path):
                return os.path.normpath(file_path)

        if default_path:
            default_path = os.path.expanduser(default_path)
            default_path = os.path.expandvars(default_path)

            if file_path:
                return os.path.normpath(os.path.join(default_path, file_path))
            else:
                return os.path.normpath(default_path)

        return file_path

    @staticmethod
    def section_from_instance(instance):
        """Get a configuration section Python C{str} composed of the Python module and Python class name.

        @param instance: A Python instance (i.e. object)
        @type instance: object
        @return: Configuration file section string
        @rtype: str
        """

        # For Python "type" instances the "__name__" instance variable provides the name of the type, while for
        # Python "object" instances the "__class__" variable provides the Python "type" object.

        if isinstance(instance, type):
            return '.'.join((instance.__module__, instance.__name__))
        else:
            return '.'.join((instance.__module__, instance.__class__.__name__))

    @classmethod
    def from_file_path_list(cls, file_path_list):
        """Create a C{bsf.standards.Configuration} object based on a Python C{list} of Python C{str} file paths.

        Both, user and variable expansion gets applied to each file path.
        Identical files are read only once.

        @param file_path_list: Python C{list} of Python C{str} or C{unicode} configuration file path objects
        @type file_path_list: list[str | unicode]
        @return: C{bsf.standards.Configuration}
        @rtype: bsf.standards.Configuration
        @raise Exception: Configuration file path does not exist
        """
        # Expand each file_path for user and variable names.
        temporary_list = map(lambda x: cls.get_absolute_path(file_path=x), file_path_list)
        """ @type temporary_list: list[str | unicode] """

        file_path_list = list()
        for temporary_path in temporary_list:
            file_exists = False
            for file_path in file_path_list:
                if os.path.samefile(temporary_path, file_path):
                    file_exists = True
                    break
            if not file_exists:
                file_path_list.append(temporary_path)

        # Since ConfigParser options are used as command line options,
        # they have to be case sensitive.
        # Hence, override optionxform() with str().

        config_parser = ConfigParser.SafeConfigParser()
        config_parser.optionxform = str

        configuration = cls(file_path_list=file_path_list, config_parser=config_parser)

        configuration._config_path_list = configuration.config_parser.read(filenames=configuration.file_path_list)

        if len(configuration._config_path_list) == 0:
            raise Exception(
                'None of the configuration files exists:\n' + repr(configuration.file_path_list))

        return configuration

    def __init__(self, file_path_list=None, config_parser=None):
        """Initialise a C{bsf.standards.Configuration}.

        @param file_path_list: Python C{list} of Python C{str} or C{unicode} configuration file path objects
        @type file_path_list: list[str | unicode]
        @param config_parser: Python C{ConfigParser.SafeConfigParser}
        @type config_parser: ConfigParser.SafeConfigParser
        @return:
        @rtype:
        """

        super(Configuration, self).__init__()

        if file_path_list is None:
            self.file_path_list = list()
        else:
            assert isinstance(file_path_list, list)
            self.file_path_list = file_path_list

        if config_parser is None:
            self.config_parser = ConfigParser.SafeConfigParser()
        else:
            assert isinstance(config_parser, ConfigParser.SafeConfigParser)
            self.config_parser = config_parser

        self._config_path_list = None
        """ @type _config_path_list: list[str | unicode] """

        return

    def trace(self, level):
        """Trace a C{bsf.standards.Configuration}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: list[str | unicode]
        """
        indent = '  ' * level

        str_list = list()
        """ @type str_list: list[str | unicode] """

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  file_path_list: {!r}\n'.format(indent, self.file_path_list))
        str_list.append('{}  config_parser:  {!r}\n'.format(indent, self.config_parser))

        return str_list

    def get_expanded_directory(self, section, option):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e. on UNIX ~ or ~user and
        any environment variables i.e. on UNIX ${NAME} or $NAME.
        @param section: Configuration file section string
        @type section: str | unicode
        @param option: Configuration file option string
        @type option: str | unicode
        @return: Expanded directory
        @rtype: None | str | unicode
        """
        return self.get_absolute_path(file_path=self.config_parser.get(section=section, option=option))


class Default(object):
    """The C{bsf.standards.Default} class represents the application or library default configuration.

    Attributes:
    @cvar global_default: Global C{bsf.standards.Default}
    @type global_default: bsf.standards.Default
    @ivar directory_cache: Local cache directory on the compute node (e.g. /dev/shm)
    @type directory_cache: str | unicode
    @ivar directory_home: Home directory for all data
    @type directory_home: str | unicode
    @ivar directory_illumina_run: Sub-directory for Illumina Run Folders
    @type directory_illumina_run: str | unicode
    @ivar directory_illumina_sav: Sub-directory for Illumina Sequence Analysis Viewer (SAV) Folders
    @type directory_illumina_sav: str | unicode
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
    @ivar directory_transcriptomes: Sub-directory for transcriptomes
    @type directory_transcriptomes: str | unicode
    @ivar directory_gatk_bundle: Sub-directory for GATK bundle data
    @type directory_gatk_bundle: str | unicode
    @ivar directory_intervals: Directory for interval list files
    @type directory_intervals: str | unicode
    @ivar directory_cosmic: Sub-directory for COSMIC data
    @type directory_cosmic: str | unicode
    @ivar directory_snpeff_data: snpEff database directory
    @type directory_snpeff_data: str | unicode
    @ivar indices: Python C{dict} of program name key and index directory name value data
    @type indices: dict[str, str]
    @ivar genome_aliases_ucsc_dict: Alias of genome assembly names for the UCSC Genome Browser
    @type genome_aliases_ucsc_dict: dict[str, str]
    """

    global_default = None

    @staticmethod
    def get_global_default():
        """Get the global C{bsf.standards.Default} configuration and initialise it, if not already done so.

        @return: C{bsf.standards.Default}
        @rtype: bsf.standards.Default
        """

        if Default.global_default is None:
            Default.global_default = Default.from_configuration(
                configuration=Configuration.get_global_configuration())

        return Default.global_default

    @classmethod
    def from_configuration(cls, configuration):
        """Create a new C{bsf.standards.Default} object from a C{bsf.standards.Configuration} object.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @return: C{bsf.standards.Default}
        @rtype: bsf.standards.Default
        """
        assert isinstance(configuration, Configuration)

        default = cls()
        default.set_configuration(configuration=configuration)

        return default

    def __init__(
            self,
            directory_cache=None,
            directory_home=None,
            directory_illumina_run=None,
            directory_illumina_sav=None,
            directory_sequences=None,
            directory_samples=None,
            directory_projects=None,
            directory_public_html=None,
            directory_genomes=None,
            directory_transcriptomes=None,
            directory_gatk_bundle=None,
            directory_intervals=None,
            directory_cosmic=None,
            directory_snpeff_data=None,
            indices=None,
            genome_aliases_ucsc_dict=None):
        """Initialise a C{bsf.standards.Default} object.

        @param directory_cache: Local cache directory on the compute node (e.g. /dev/shm)
        @type directory_cache: str | unicode
        @param directory_home: Home directory for all data
        @type directory_home: str | unicode
        @param directory_illumina_run: Sub-directory for Illumina Run Folders
        @type directory_illumina_run: str | unicode
        @param directory_illumina_sav: Sub-directory for Illumina Sequence Analysis Viewer (SAV) Folders
        @type directory_illumina_sav: str | unicode
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
        @param directory_transcriptomes: Sub-directory for transcriptomes
        @type directory_transcriptomes: str | unicode
        @param directory_gatk_bundle: Sub-directory for GATK bundle data
        @type directory_gatk_bundle: str | unicode
        @param directory_intervals: Directory for interval list files
        @type directory_intervals: str | unicode
        @param directory_cosmic: Sub-directory for COSMIC data
        @type directory_cosmic: str | unicode
        @param directory_snpeff_data: snpEff database directory
        @type directory_snpeff_data: str | unicode
        @param indices: Python C{dict} of program name key and index directory name value data
        @type indices: dict[str, str]
        @param genome_aliases_ucsc_dict: Alias of genome assembly names for the UCSC Genome Browser
        @type genome_aliases_ucsc_dict: dict[str, str]
        @return:
        @rtype:
        """

        super(Default, self).__init__()

        # Set directory information.

        if directory_cache is None:
            self.directory_cache = str()
        else:
            self.directory_cache = directory_cache

        if directory_home is None:
            self.directory_home = str()
        else:
            self.directory_home = directory_home

        if directory_illumina_run is None:
            self.directory_illumina_run = str()
        else:
            self.directory_illumina_run = directory_illumina_run

        if directory_illumina_sav is None:
            self.directory_illumina_sav = str()
        else:
            self.directory_illumina_sav = directory_illumina_sav

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

        if directory_transcriptomes is None:
            self.directory_transcriptomes = str()
        else:
            self.directory_transcriptomes = directory_transcriptomes

        if directory_gatk_bundle is None:
            self.directory_gatk_bundle = str()
        else:
            self.directory_gatk_bundle = directory_gatk_bundle

        if directory_intervals is None:
            self.directory_intervals = str()
        else:
            self.directory_intervals = directory_intervals

        if directory_cosmic is None:
            self.directory_cosmic = str()
        else:
            self.directory_cosmic = directory_cosmic

        if directory_snpeff_data is None:
            self.directory_snpeff_data = str()
        else:
            self.directory_snpeff_data = directory_snpeff_data

        # Set index information.

        if indices is None:
            self.indices = dict()
        else:
            self.indices = indices

        # Set Genome Aliases for the UCSC Genome Browser.

        if genome_aliases_ucsc_dict is None:
            self.genome_aliases_ucsc_dict = dict()
        else:
            self.genome_aliases_ucsc_dict = genome_aliases_ucsc_dict

        return

    def set_configuration(self, configuration):
        """Set instance variables of a C{bsf.standards.Default} object via a C{bsf.standards.Configuration} section.

        For each instance variable a configuration option has to be present.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @return:
        @rtype:
        """
        assert isinstance(configuration, Configuration)

        cp = configuration.config_parser

        # Reading configuration cannot be done via a single Python dict,
        # because each option really needs defining.

        section = 'directories'

        self.directory_cache = cp.get(section=section, option='cache')
        self.directory_home = cp.get(section=section, option='home')
        self.directory_illumina_run = cp.get(section=section, option='illumina_run')
        self.directory_illumina_sav = cp.get(section=section, option='illumina_sav')
        self.directory_sequences = cp.get(section=section, option='sequences')
        self.directory_samples = cp.get(section=section, option='samples')
        self.directory_projects = cp.get(section=section, option='projects')
        self.directory_public_html = cp.get(section=section, option='public_html')
        self.directory_genomes = cp.get(section=section, option='genomes')
        self.directory_transcriptomes = cp.get(section=section, option='transcriptomes')
        self.directory_gatk_bundle = cp.get(section=section, option='gatk_bundle')
        self.directory_intervals = cp.get(section=section, option='intervals')
        self.directory_cosmic = cp.get(section=section, option='cosmic')
        self.directory_snpeff_data = cp.get(section=section, option='snpeff_data')

        section = 'indices'

        for option in cp.options(section=section):
            self.indices[option] = cp.get(section=section, option=option)

        section = 'genome_aliases_ucsc'

        if cp.has_section(section=section):
            for option in cp.options(section=section):
                self.genome_aliases_ucsc_dict[option] = cp.get(section=section, option=option)

        return

    @staticmethod
    def absolute_cache():
        """
        Get the absolute directory path for the cache directory.

        @return: Absolute path to the cache directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        return default.directory_cache

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
    def absolute_illumina_run():
        """Get the absolute directory path for Illumina Run Folders.

        @return: Absolute path to the Illumina Run Folder directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_illumina_run):
            return default.directory_illumina_run
        else:
            return os.path.join(default.directory_home, default.directory_illumina_run)

    @staticmethod
    def absolute_illumina_sav():
        """Get the absolute directory path for Illumina Sequence Analysis Viewer (SAV) Folders.

        @return: Absolute path to the Illumina Sequence Analysis Viewer (SAV) Folders directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_illumina_sav):
            return default.directory_illumina_sav
        else:
            return os.path.join(default.directory_home, default.directory_illumina_sav)

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
        """Get the absolute directory path for the Genome Analysis Toolkit (GATK) bundle.

        @param gatk_bundle_version: The GATK bundle version
        @type gatk_bundle_version: str
        @param genome_version: The genome version (e.g. b37, ...)
        @type genome_version: str
        @return: Absolute path to the GATK bundle directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        file_path = default.directory_gatk_bundle

        if gatk_bundle_version:
            file_path = os.path.join(file_path, gatk_bundle_version)

        if genome_version:
            file_path = os.path.join(file_path, genome_version)

        return file_path

    @staticmethod
    def absolute_genome_resource(genome_version):
        """Get the absolute path to the genome resource directory.

        @param genome_version: The genome version (e.g. mm10, ...)
        @type genome_version: str
        @return: Absolute path to the genome resource directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if genome_version:
            return os.path.join(default.directory_genomes, genome_version)
        else:
            return default.directory_genomes

    @staticmethod
    def absolute_genome_index(genome_version, genome_index):
        """Get the absolute file path to a genome index.

        @param genome_version: Genome version (e.g. mm10, ...)
        @type genome_version: str
        @param genome_index: Genome index (e.g. bowtie2, ...)
        @type genome_index: str
        @return: Absolute path to the genome index
        @rtype: str | unicode
        @raise Exception: Unknown genome index name
        """

        default = Default.get_global_default()

        if genome_index not in default.indices:
            raise Exception('Unknown genome index name ' + repr(genome_index) + '.')

        return os.path.join(
            Default.absolute_genome_resource(genome_version=genome_version),
            default.indices[genome_index],
            genome_version)

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

        return Default.absolute_genome_index(genome_version=genome_version, genome_index=genome_index) + '.fa'

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
    def absolute_cosmic():
        """Get the absolute directory path for COSMIC data files.

        @return: Absolute path to the COSMIC directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if os.path.isabs(default.directory_cosmic):
            return default.directory_cosmic
        else:
            return os.path.join(default.directory_home, default.directory_cosmic)

    @staticmethod
    def absolute_transcriptome_resource(transcriptome_version):
        """Get the absolute path to the transcriptome resource directory.

        @param transcriptome_version: The transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: str
        @return: Absolute path to the transcriptome resource directory
        @rtype: str | unicode
        """

        default = Default.get_global_default()

        if transcriptome_version:
            return os.path.join(default.directory_transcriptomes, transcriptome_version)
        else:
            return default.directory_transcriptomes

    @staticmethod
    def absolute_transcriptome_index(transcriptome_version, transcriptome_index):
        """Get the absolute file path to a transcriptome index.

        @param transcriptome_version: Transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: str
        @param transcriptome_index: Transcriptome index (e.g. tophat, ...)
        @type transcriptome_index: str
        @return: Absolute path to the transcriptome index
        @rtype: str | unicode
        @raise Exception: Unknown transcriptome index name
        """

        default = Default.get_global_default()

        if transcriptome_index not in default.indices:
            raise Exception('Unknown transcriptome index name ' + repr(transcriptome_index) + '.')

        return os.path.join(
            Default.absolute_transcriptome_resource(transcriptome_version=transcriptome_version),
            default.indices[transcriptome_index],
            transcriptome_version)

    @staticmethod
    def absolute_transcriptome_gtf(transcriptome_version, transcriptome_index):
        """Get the absolute file path to a transcriptome in GTF format.
        
        @param transcriptome_version: Transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: str
        @param transcriptome_index: Transcriptome index (e.g. tophat, ...)
        @type transcriptome_index: str
        @return: Absolute path to the transcriptome GTF file
        @rtype: str | unicode
        @raise Exception: Unknown transcriptome index name
        """

        return Default.absolute_transcriptome_index(
            transcriptome_version=transcriptome_version,
            transcriptome_index=transcriptome_index) + '.gtf'

    @staticmethod
    def genome_alias_ucsc(genome_version):
        """Resolve a genome version to an eventual, UCSC-specific alias.

        If an alias has not been defined the original genome version will be returned.
        @param genome_version: Genome version
        @type genome_version: str
        @return: UCSC genome version
        @rtype: str
        """
        default = Default.get_global_default()

        if default.genome_aliases_ucsc_dict is not None and genome_version in default.genome_aliases_ucsc_dict:
            return default.genome_aliases_ucsc_dict[genome_version]
        else:
            return genome_version


class InitialisationBase(object):
    """The C{bsf.standards.InitialisationBase} class is the base class for global configuration defaults.

    The defaults are read from the a particular section of the global configuration file.
    Attributes:
    @cvar section: C{SafeConfigParser} section
    @type section: str | unicode
    """
    section = None

    @classmethod
    def get(cls, option=None):
        """Get the value for a configuration option in the section defined by the sub-class section class variable.

        This method is a re-implementation of the C{SafeConfigParser.get()} method that returns C{None}
        upon non-existing sections or options.
        @param option: Configuration option
        @type option: None | str | unicode
        @return: Configuration value
        @rtype: None | str | unicode
        """
        if not cls.section:
            return

        if not option:
            return

        if Configuration.get_global_configuration().config_parser.has_option(
                section=cls.section,
                option=option):
            return Configuration.get_global_configuration().config_parser.get(
                section=cls.section,
                option=option)
        else:
            return


class JavaClassPath(InitialisationBase):
    """The C{bsf.standards.JavaClassPath} class models Java class path defaults.

    The defaults are read from the [classpath] section of the global configuration file.
    Attributes:
    @cvar section: C{SafeConfigParser} section for Java class path defaults
    @type section: str | unicode
    """
    section = 'classpath'

    @classmethod
    def get_gatk(cls):
        """Get the GATK Java class path.

        @return: GATK Java class path
        @rtype: None | str | unicode
        """
        return cls.get(option='gatk')

    @classmethod
    def get_illumina2bam(cls):
        """Get the Illumina2bam tools Java class path.

        @return: Illumina2bam tools Java class path
        @rtype: None | str | unicode
        """
        return cls.get(option='illumina2bam')

    @classmethod
    def get_picard(cls):
        """Get the Picard tools Java class path.

        @return: Picard tools Java class path
        @rtype: None | str | unicode
        """
        return cls.get(option='picard')

    @classmethod
    def get_snpeff(cls):
        """Get the snpEff Java class path.

        @return: snpEff Java class path
        @rtype: None | str | unicode
        """
        return cls.get(option='snpeff')

    @classmethod
    def get_trimmomatic(cls):
        """Get the Trimmomatic Java class path.

        @return: Trimmomatic Java class path
        @rtype: None | str | unicode
        """
        return cls.get(option='trimmomatic')

    @classmethod
    def get_vcf_filter(cls):
        """Get the VCF.Filter Java class path.

        @return: VCF.Filter Java class path
        @rtype: None | str | unicode
        """
        return cls.get(option='vcf_filter')


class EnsemblVEP(object):
    """The C{bsf.standards.EnsemblVEP} class models Ensembl Variant Effect Predictor (VEP) defaults.

    The defaults are read form the [ensembl_vep_{genome_version}] section of the global configuration file.
    Attributes:
    @cvar section: C{SafeConfigParser} section for the Ensembl Variant Effect Predictor (VEP)
    @type section: str | unicode
    """
    section = 'ensembl_vep'

    @classmethod
    def get_section(cls, genome_version=None):
        """Get the Configuration section from a genome assembly version.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Configuration section
        @rtype: str
        """
        if genome_version:
            return '_'.join((cls.section, genome_version))
        else:
            return cls.section

    @classmethod
    def get(cls, option=None, genome_version=None):
        """Get the value for a configuration option.

        @param option: Configuration option
        @type option: None | str | unicode
        @param genome_version: Genome assembly version
        @type genome_version: None | str | unicode
        @return: Configuration value
        @rtype: None | str | unicode
        """
        if not option:
            return

        section = cls.get_section(genome_version=genome_version)
        if Configuration.get_global_configuration().config_parser.has_option(section=section, option=option):
            return Configuration.get_global_configuration().config_parser.get(section=section, option=option)
        else:
            return

    @classmethod
    def get_expanded_directory(cls, option=None, genome_version=None):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e. on UNIX ~ or ~user and
        any environment variables i.e. on UNIX ${NAME} or $NAME.
        @param option: Configuration option
        @type option: None | str | unicode
        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Expanded directory
        @rtype: None | str | unicode
        """
        if not option:
            return

        section = cls.get_section(genome_version=genome_version)
        if Configuration.get_global_configuration().config_parser.has_option(section=section, option=option):
            return Configuration.get_global_configuration().get_expanded_directory(section=section, option=option)
        else:
            return

    @classmethod
    def get_directory_cache(cls, genome_version=None):
        """Get the cache directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Cache directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_cache', genome_version=genome_version)

    @classmethod
    def get_directory_fasta(cls, genome_version=None):
        """Get the FASTA directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: FASTA directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_fasta', genome_version=genome_version)

    @classmethod
    def get_directory_plugin(cls, genome_version=None):
        """Get the plug-ins directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Plug-ins directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_plugin', genome_version=genome_version)

    @classmethod
    def get_directory_source(cls, genome_version=None):
        """Get the source directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Source directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_source', genome_version=genome_version)

    @classmethod
    def get_name_assembly(cls, genome_version=None):
        """Get the genome assembly name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Genome assembly name
        @rtype: None | str | unicode
        """
        return cls.get(option='name_assembly', genome_version=genome_version)

    @classmethod
    def get_name_species(cls, genome_version=None):
        """Get the scientific species name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Scientific species name
        @rtype: None | str | unicode
        """
        return cls.get(option='name_species', genome_version=genome_version)

    @classmethod
    def get_sql_user(cls, genome_version=None):
        """Get the SQL database user name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database user name
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_user', genome_version=genome_version)

    @classmethod
    def get_sql_pass(cls, genome_version=None):
        """Get the SQL database password.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database password
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_pass', genome_version=genome_version)

    @classmethod
    def get_sql_host(cls, genome_version=None):
        """Get the SQL database host name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database host name
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_host', genome_version=genome_version)

    @classmethod
    def get_sql_port(cls, genome_version=None):
        """Get the SQL database TCP/IP port number.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database TCP/IP port number
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_port', genome_version=genome_version)


class Operator(InitialisationBase):
    """The C{bsf.standards.Operator} class models operator defaults.

    The defaults are read form the [operator] section of the global configuration file.
    Attributes:
    @cvar section: C{SafeConfigParser} section for the operator
    @type section: str | unicode
    """
    section = 'operator'

    @classmethod
    def get_contact(cls):
        """Get the operator contact information.

        @return: Operator contact information
        @rtype: None | str | unicode
        """
        return cls.get(option='contact')

    @classmethod
    def get_e_mail(cls):
        """Get the operator e-mail information.

        @return: Operator e-mail information
        @rtype: None | str | unicode
        """
        return cls.get(option='e_mail')

    @classmethod
    def get_institution(cls):
        """Get the operator institution information.

        @return: Operator institution information
        @rtype: None | str | unicode
        """
        return cls.get(option='institution')

    @classmethod
    def get_sequencing_centre(cls):
        """Get the operator sequencing centre information.

        @return: Operator sequencing centre information
        @rtype: None | str | unicode
        """
        return cls.get(option='sequencing_centre')


class UCSC(InitialisationBase):
    """The C{bsf.standards.UCSC} class models UCSC Genome Browser uniform resource locator (URL) defaults.

    The defaults are read form the [ucsc] section of the global configuration file.
    Attributes:
    @cvar section: C{SafeConfigParser} section for the operator
    @type section: str | unicode
    """

    section = 'ucsc'

    @classmethod
    def get_protocol(cls):
        """Get the UCSC Genome Browser URL protocol (i.e. 'http' or 'https').

        @return: Protocol
        @rtype: str
        """
        return cls.get(option='protocol')

    @classmethod
    def get_host_name(cls):
        """Get the UCSC Genome Browser URL host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...).

        @return: Host name
        @rtype: str
        """
        return cls.get(option='host_name')


class URL(InitialisationBase):
    """The C{bsf.standards.URL} class models web server uniform resource locator (URL) defaults.

    The defaults are read form the [url] section of the global configuration file.
    Attributes:
    @cvar section: C{SafeConfigParser} section for the operator
    @type section: str | unicode
    """

    section = 'url'

    @classmethod
    def get_protocol(cls):
        """Get the web server URL protocol (i.e. 'http' or 'https').

        @return: Protocol
        @rtype: str
        """
        return cls.get(option='protocol')

    @classmethod
    def get_host_name(cls):
        """Get the web server host name.

        @return: Host name
        @rtype: str
        """
        return cls.get(option='host_name')

    @classmethod
    def get_relative_dna(cls):
        """Get the relative URL to the DNA directory.

        @return: Relative 'DNA' URL path
        @rtype: str
        """
        return cls.get(option='relative_dna')

    @classmethod
    def get_relative_projects(cls):
        """Get the relative URL to the analysis projects directory.

        @return: Relative 'projects' URL path
        @rtype: str
        """
        return cls.get(option='relative_projects')

    @classmethod
    def get_absolute_base(cls):
        """Get the absolute URL to the base directory.

        @return: URL string
        @rtype: str
        """
        if cls.get_protocol():
            return cls.get_protocol() + '://' + cls.get_host_name()
        else:
            return '//' + cls.get_host_name()

    @classmethod
    def get_absolute_dna(cls):
        """Get the absolute URL to the DNA directory.

        @return: URL string
        @rtype: str
        """
        return '/'.join((cls.get_absolute_base(), cls.get_relative_dna()))

    @classmethod
    def get_absolute_projects(cls):
        """Get the absolute URL to the analysis projects directory.

        @return: URL string
        @rtype: str
        """
        return '/'.join((cls.get_absolute_base(), cls.get_relative_projects()))
