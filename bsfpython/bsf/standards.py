"""bsf.standards

A package of classes and methods modelling configuration and default information.
"""

#
# Copyright 2013 - 2016 Michael K. Schuster
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
import os


class Configuration(object):
    """The C{bsf.standards.Configuration} class represents one or more UNIX-style initialisation (*.ini) files.

    A C{bsf.standards.Configuration} has an associated Python C{ConfigParser.SafeConfigParser} to parse the file(s).

    Attributes:
    @ivar file_path_list: C{bsf.standards.Configuration} file path
    @type file_path_list: list[str | unicode]
    @ivar config_parser: Python C{ConfigParser.SafeConfigParser}
    @type config_parser: ConfigParser.SafeConfigParser
    """

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

        @param file_path_list: Python C{list} of Python C{str} or C{unicode} configuration file path objects
        @type file_path_list: list[str | unicode]
        @return: C{bsf.standards.Configuration}
        @rtype: bsf.standards.Configuration
        @raise Exception: Configuration file path does not exist
        """
        assert isinstance(file_path_list, list)

        # Expand each file_path for user and variable names.
        expanded_list = list()
        """ @type expanded_list: list[str | unicode] """

        for file_path in file_path_list:
            file_path = os.path.expanduser(path=file_path)
            file_path = os.path.expandvars(path=file_path)
            file_path = os.path.normpath(path=file_path)
            expanded_list.append(file_path)

        # Since ConfigParser options are used as command line options,
        # they have to be case sensitive.
        # Hence, override optionxform() with str().

        config_parser = SafeConfigParser()
        config_parser.optionxform = str

        configuration = cls(file_path_list=expanded_list, config_parser=config_parser)

        configuration._config_path_list = configuration.config_parser.read(filenames=configuration.file_path_list)

        if len(configuration._config_path_list) == 0:
            raise Exception(
                'None of the configuration files exists: {!r}'.format(configuration.file_path_list))

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
            self.config_parser = SafeConfigParser()
        else:
            assert isinstance(config_parser, SafeConfigParser)
            self.config_parser = config_parser

        self._config_path_list = None
        """ @type _config_path_list: list[str | unicode] """

        return

    def trace(self, level):
        """Trace a C{bsf.standards.Configuration}.

        @param level: Indentation level
        @type level: int
        @return: Trace information
        @rtype: str
        """

        indent = '  ' * level
        output = str()
        output += '{}{!r}\n'.format(indent, self)
        output += '{}  file_path_list: {!r}\n'.format(indent, self.file_path_list)
        output += '{}  config_parser:  {!r}\n'.format(indent, self.config_parser)

        return output

    def get_expanded_directory(self, config_section, config_option):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e. on UNIX ~ or ~user and
        any environment variables i.e. on UNIX ${NAME} or $NAME.
        @param config_section: Configuration file section string
        @type config_section: str
        @param config_option: Configuration file option string
        @type config_option: str
        @return: Expanded directory
        @rtype: str
        """

        directory = self.config_parser.get(config_section, config_option)
        directory = os.path.expanduser(path=directory)
        directory = os.path.expandvars(path=directory)
        directory = os.path.normpath(path=directory)

        return directory


class Default(object):
    """The C{bsf.standards.Default} class represents the application or library default configuration.

    Attributes:
    @cvar global_default: Global C{bsf.standards.Default}
    @type global_default: bsf.standards.Default
    @cvar global_file_path: Global configuration file
    @type global_file_path: str | unicode
    @ivar classpath_gatk: Genome Analysis Toolkit Java Archive (JAR) class path directory
    @type classpath_gatk: str | unicode
    @ivar classpath_illumina2bam: Illumina2bam Java Archive (JAR) class path directory
    @type classpath_illumina2bam: str | unicode
    @ivar classpath_picard: Picard Java Archive (JAR) class path directory
    @type classpath_picard: str | unicode
    @ivar classpath_snpeff: snpEff Java Archive (JAR) class path directory
    @type classpath_snpeff: str | unicode
    @ivar classpath_vcf_filter: VCF.Filter Java Archive (JAR) class path directory
    @type classpath_vcf_filter: str | unicode
    @ivar directory_cache: Local cache directory on the compute node (e.g. /dev/shm)
    @type directory_cache: str | unicode
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
    @ivar directory_transcriptomes: Sub-directory for transcriptomes
    @type directory_transcriptomes: str | unicode
    @ivar directory_gatk_bundle: Sub-directory for GATK bundle data
    @type directory_gatk_bundle: str | unicode
    @ivar directory_intervals: Directory for interval list files
    @type directory_intervals: str | unicode
    @ivar directory_snpeff_data: snpEff database directory
    @type directory_snpeff_data: str | unicode
    @ivar directory_vep_cache: Ensembl Variant Effect Predictor (VEP) cache directory
    @type directory_vep_cache: str | unicode
    @ivar directory_vep_plugins: Ensembl Variant Effect Predictor (VEP) plug-ins directory
    @type directory_vep_plugins: str | unicode
    @ivar directory_vep_src: Ensembl Variant Effect Predictor (VEP) source directory
    @type directory_vep_src: str | unicode
    @ivar indices: Python C{dict} of program name key and index directory name value data
    @type indices: dict[str, str]
    @ivar operator_contact: Contact e-mail address
    @type operator_contact: str
    @ivar operator_e_mail: Operator e-mail
    @type operator_e_mail: str
    @ivar operator_institution: Institution name
    @type operator_institution: str
    @ivar operator_sequencing_centre: BAM sequencing centre code
    @type operator_sequencing_centre: str
    @ivar genome_aliases_ucsc_dict: Alias of genome assembly names for the UCSC Genome Browser
    @type genome_aliases_ucsc_dict: dict[str, str]
    @ivar ucsc_protocol: UCSC Genome Browser URL protocol (i.e. http, https, ...)
    @type ucsc_protocol: str
    @ivar ucsc_host_name: UCSC Genome Browser URL host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...)
    @type ucsc_host_name: str
    @ivar url_protocol: URL protocol (i.e. HTTP)
    @type url_protocol: str
    @ivar url_host_name: URL host name
    @type url_host_name:str
    @ivar url_relative_dna: Sub-directory for DNA sequences
    @type url_relative_dna: str
    @ivar url_relative_projects: Sub-directory for analysis projects
    @type url_relative_projects: str
    """

    global_default = None

    global_file_path = '~/.bsfpython.ini'
    global_file_path = os.path.expanduser(path=global_file_path)
    global_file_path = os.path.expandvars(path=global_file_path)
    global_file_path = os.path.normpath(path=global_file_path)

    @staticmethod
    def get_absolute_path(file_path, default_path=None):
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
        @rtype: str | unicode
        """

        absolute_path = os.path.expanduser(path=file_path)
        absolute_path = os.path.expandvars(path=absolute_path)

        if default_path and not os.path.isabs(absolute_path):
            default_path = os.path.expanduser(path=default_path)
            default_path = os.path.expandvars(path=default_path)
            absolute_path = os.path.join(default_path, absolute_path)

        return os.path.normpath(absolute_path)

    @staticmethod
    def get_global_default():
        """Get the global C{bsf.standards.Default} configuration and initialise it, if not already done so.

        @return: C{bsf.standards.Default}
        @rtype: bsf.standards.Default
        """

        if Default.global_default is None:
            Default.global_default = Default.from_global_file_path()

        return Default.global_default

    @classmethod
    def from_global_file_path(cls):
        """Create a new C{bsf.standards.Default} object from the global default configuration file.

        The default configuration is based on the file $HOME/.bsfpython.ini in the user's home directory.

        @return: C{bsf.standards.Default}
        @rtype: bsf.standards.Default
        """

        return cls.from_config_path(config_path=Default.global_file_path)

    @classmethod
    def from_config_path(cls, config_path):
        """Create a new C{bsf.standards.Default} object from a UNIX-style configuration file.

        @param config_path: UNIX-style configuration file path
        @type config_path: str | unicode
        @return: C{bsf.standards.Default}
        @rtype: bsf.standards.Default
        """

        return cls.from_configuration(configuration=Configuration.from_file_path_list(file_path_list=[config_path]))

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
            classpath_gatk=None,
            classpath_illumina2bam=None,
            classpath_picard=None,
            classpath_snpeff=None,
            classpath_vcf_filter=None,
            directory_cache=None,
            directory_home=None,
            directory_runs_illumina=None,
            directory_sequences=None,
            directory_samples=None,
            directory_projects=None,
            directory_public_html=None,
            directory_genomes=None,
            directory_transcriptomes=None,
            directory_gatk_bundle=None,
            directory_intervals=None,
            directory_snpeff_data=None,
            directory_vep_cache=None,
            directory_vep_plugins=None,
            directory_vep_src=None,
            indices=None,
            operator_contact=None,
            operator_e_mail=None,
            operator_institution=None,
            operator_sequencing_centre=None,
            genome_aliases_ucsc_dict=None,
            ucsc_protocol=None,
            ucsc_host_name=None,
            url_protocol=None,
            url_host_name=None,
            url_relative_dna=None,
            url_relative_projects=None):
        """Initialise a C{bsf.standards.Default}.

        @param classpath_gatk: Genome Analysis Toolkit Java Archive (JAR) class path directory
        @type classpath_gatk: str | unicode
        @param classpath_illumina2bam: Illumina2bam Java Archive (JAR) class path directory
        @type classpath_illumina2bam: str | unicode
        @param classpath_picard: Picard Java Archive (JAR) class path directory
        @type classpath_picard: str | unicode
        @param classpath_snpeff: snpEff Java Archive (JAR) class path directory
        @type classpath_snpeff: str | unicode
        @param classpath_vcf_filter: VCF.Filter Java Archive (JAR) class path directory
        @type classpath_vcf_filter: str | unicode
        @param directory_cache: Local cache directory on the compute node (e.g. /dev/shm)
        @type directory_cache: str | unicode
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
        @param directory_transcriptomes: Sub-directory for transcriptomes
        @type directory_transcriptomes: str | unicode
        @param directory_gatk_bundle: Sub-directory for GATK bundle data
        @type directory_gatk_bundle: str | unicode
        @param directory_intervals: Directory for interval list files
        @type directory_intervals: str | unicode
        @param directory_snpeff_data: snpEff database directory
        @type directory_snpeff_data: str | unicode
        @param directory_vep_cache: Ensembl Variant Effect Predictor (VEP) cache directory
        @type directory_vep_cache: str | unicode
        @param directory_vep_plugins: Ensembl Variant Effect Predictor (VEP) plug-ins directory
        @type directory_vep_plugins: str | unicode
        @param directory_vep_src: Ensembl Variant Effect Predictor (VEP) source directory
        @type directory_vep_src: str | unicode
        @param indices: Python C{dict} of program name key and index directory name value data
        @type indices: dict[str, str]
        @param operator_contact: Contact e-mail address
        @type operator_contact: str
        @param operator_e_mail: Operator e-mail
        @type operator_e_mail: str
        @param operator_institution: Institution name
        @type operator_institution: str
        @param operator_sequencing_centre: BAM sequencing centre code
        @type operator_sequencing_centre: str
        @param genome_aliases_ucsc_dict: Alias of genome assembly names for the UCSC Genome Browser
        @type genome_aliases_ucsc_dict: dict[str, str]
        @param ucsc_protocol: UCSC Genome Browser URL protocol (i.e. http, https, ...)
        @type ucsc_protocol: str
        @param ucsc_host_name: UCSC Genome Browser URL host name (e.g. genome.ucsc.edu, genome-euro.ucsc.edu, ...)
        @type ucsc_host_name: str
        @param url_protocol: URL protocol (i.e. HTTP)
        @type url_protocol: str
        @param url_host_name: URL host name
        @type url_host_name:str
        @param url_relative_dna: Sub-directory for DNA sequences
        @type url_relative_dna: str
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

        if classpath_vcf_filter is None:
            self.classpath_vcf_filter = str()
        else:
            self.classpath_vcf_filter = classpath_vcf_filter

        # Set directory information.

        if directory_cache is None:
            self.directory_cache = str()
        else:
            self.directory_cache = directory_cache

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

        if directory_snpeff_data is None:
            self.directory_snpeff_data = str()
        else:
            self.directory_snpeff_data = directory_snpeff_data

        if directory_vep_cache is None:
            self.directory_vep_cache = str()
        else:
            self.directory_vep_cache = directory_vep_cache

        if directory_vep_plugins is None:
            self.directory_vep_plugins = str()
        else:
            self.directory_vep_plugins = directory_vep_plugins

        if directory_vep_src is None:
            self.directory_vep_src = str()
        else:
            self.directory_vep_src = directory_vep_src

        # Set index information.

        if indices is None:
            self.indices = dict()
        else:
            self.indices = indices

        # Set operator information.

        if operator_contact is None:
            self.operator_contact = str()
        else:
            self.operator_contact = operator_contact

        if operator_e_mail is None:
            self.operator_e_mail = str()
        else:
            self.operator_e_mail = operator_e_mail

        if operator_institution is None:
            self.operator_institution = str()
        else:
            self.operator_institution = operator_institution

        if operator_sequencing_centre is None:
            self.operator_sequencing_centre = str()
        else:
            self.operator_sequencing_centre = operator_sequencing_centre

        # Set Genome Aliases for the UCSC Genome Browser.

        if genome_aliases_ucsc_dict is None:
            self.genome_aliases_ucsc_dict = dict()
        else:
            self.genome_aliases_ucsc_dict = genome_aliases_ucsc_dict

        # Set UCSC Genome Browser information.

        if ucsc_protocol is None:
            self.ucsc_protocol = str()
        else:
            self.ucsc_protocol = ucsc_protocol

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

        if url_relative_dna is None:
            self.url_relative_dna = str()
        else:
            self.url_relative_dna = url_relative_dna

        if url_relative_projects is None:
            self.url_relative_projects = str()
        else:
            self.url_relative_projects = url_relative_projects

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

        section = 'classpath'

        self.classpath_gatk = cp.get(section=section, option='gatk')
        self.classpath_illumina2bam = cp.get(section=section, option='illumina2bam')
        self.classpath_picard = cp.get(section=section, option='picard')
        self.classpath_snpeff = cp.get(section=section, option='snpeff')
        self.classpath_vcf_filter = cp.get(section=section, option='vcf_filter')

        section = 'directories'

        self.directory_cache = cp.get(section=section, option='cache')
        self.directory_home = cp.get(section=section, option='home')
        self.directory_runs_illumina = cp.get(section=section, option='runs_illumina')
        self.directory_sequences = cp.get(section=section, option='sequences')
        self.directory_samples = cp.get(section=section, option='samples')
        self.directory_projects = cp.get(section=section, option='projects')
        self.directory_public_html = cp.get(section=section, option='public_html')
        self.directory_genomes = cp.get(section=section, option='genomes')
        self.directory_transcriptomes = cp.get(section=section, option='transcriptomes')
        self.directory_gatk_bundle = cp.get(section=section, option='gatk_bundle')
        self.directory_intervals = cp.get(section=section, option='intervals')
        self.directory_snpeff_data = cp.get(section=section, option='snpeff_data')
        self.directory_vep_cache = cp.get(section=section, option='vep_cache')
        self.directory_vep_plugins = cp.get(section=section, option='vep_plugins')
        self.directory_vep_src = cp.get(section=section, option='vep_src')

        section = 'indices'

        for option in cp.options(section=section):
            self.indices[option] = cp.get(section=section, option=option)

        section = 'operator'

        self.operator_contact = cp.get(section=section, option='contact')
        self.operator_e_mail = cp.get(section=section, option='e_mail')
        self.operator_institution = cp.get(section=section, option='institution')
        self.operator_sequencing_centre = cp.get(section=section, option='sequencing_centre')

        section = 'genome_aliases_ucsc'

        if cp.has_section(section=section):
            for option in cp.options(section=section):
                self.genome_aliases_ucsc_dict[option] = cp.get(section=section, option=option)

        section = 'ucsc'

        self.ucsc_protocol = cp.get(section=section, option='protocol')
        self.ucsc_host_name = cp.get(section=section, option='host_name')

        section = 'url'

        self.url_protocol = cp.get(section=section, option='protocol')
        self.url_host_name = cp.get(section=section, option='host_name')
        self.url_relative_dna = cp.get(section=section, option='relative_dna')
        self.url_relative_projects = cp.get(section=section, option='relative_projects')

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
            raise Exception(
                'Unknown genome index name {!r}.'.format(genome_index))

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
            raise Exception(
                'Unknown transcriptome index name {!r}.'.format(transcriptome_index))

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
    def url_absolute_base():
        """Return the absolute URL to the web site.

        @return: URL string
        @rtype: str
        """

        default = Default.get_global_default()

        # Strip leading colons to support protocol-independent URLs.
        return '{}://{}'.format(default.url_protocol, default.url_host_name).lstrip(':')

    @staticmethod
    def url_absolute_dna():
        """Return the absolute URL to the DNA directory.

        @return: URL string
        @rtype: str
        """

        default = Default.get_global_default()

        return '/'.join((default.url_absolute_base(), default.url_relative_dna))

    @staticmethod
    def url_absolute_projects():
        """Return the absolute URL to the analysis projects directory.

        @return: URL string
        @rtype: str
        """

        default = Default.get_global_default()

        return '/'.join((default.url_absolute_base(), default.url_relative_projects))

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
