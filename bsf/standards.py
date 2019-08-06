# -*- coding: utf-8 -*-
"""Standards module

A package of classes and methods modelling configuration and default information.
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
import configparser
import os
import stat


class Configuration(object):
    """The C{bsf.standards.Configuration} class represents one or more UNIX-style initialisation (*.ini) files.

    A C{bsf.standards.Configuration} has an associated Python C{configparser.ConfigParser} to parse the file(s).

    Attributes:
    @cvar global_configuration: Global C{bsf.standards.Configuration}
    @type global_configuration: bsf.standards.Configuration
    @cvar global_file_path: Global configuration file
    @type global_file_path: str | unicode
    @ivar file_path_list: C{bsf.standards.Configuration} file path
    @type file_path_list: list[str | unicode]
    @ivar config_parser: Python C{configparser.ConfigParser}
    @type config_parser: configparser.ConfigParser
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
        temporary_list = [cls.get_absolute_path(file_path=x) for x in file_path_list]

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

        config_parser = configparser.ConfigParser()
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
        @param config_parser: Python C{configparser.ConfigParser}
        @type config_parser: configparser.ConfigParser
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
            self.config_parser = configparser.ConfigParser()
        else:
            assert isinstance(config_parser, configparser.ConfigParser)
            self.config_parser = config_parser

        self._config_path_list = None
        """ @type _config_path_list: list[str | unicode] | None """

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

    @staticmethod
    def list_from_csv(csv_string):
        """Convert a comma-separated Python C{str} into a Python C{list} of Python C{str} objects.

        All elements are stripped and only non-empty elements are appended to the list.
        @param csv_string:
        @return:
        """
        return [x.strip() for x in csv_string.split(',') if x.strip()]

    def get_list_from_csv(self, section, option):
        """Convert a comma-separated Python C{str} into a Python C{list} of Python C{str} objects.

        All elements are stripped and only non-empty elements are appended to the list.
        @param section: Configuration file section string
        @type section: str | unicode
        @param option: Configuration file option string
        @type option: str | unicode
        @return: Python C{list} of Python C{str} objects
        @rtype: list[str] | None
        """
        csv_string = self.config_parser.get(section=section, option=option)
        if csv_string:
            return self.list_from_csv(csv_string=csv_string)


class BaseSection(object):
    """The C{bsf.standards.BaseSection} class is the base class for a global configuration section.

    The defaults are read from the [{section}] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = None

    @classmethod
    def get(cls, option=None):
        """Get the value for a configuration option in the section defined by the sub-class section class variable.

        This method is a re-implementation of the C{configparser.ConfigParser.get()} method that returns C{None}
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

    @classmethod
    def getboolean(cls, option=None):
        """Get the value for a configuration option in the section defined by the sub-class section class variable.

        This method is a re-implementation of the C{configparser.ConfigParser.getboolean()} method that returns C{None}
        upon non-existing sections or options.
        @param option: Configuration option
        @type option: None | str | unicode
        @return: Configuration value
        @rtype: None | bool
        """
        if not cls.section:
            return

        if not option:
            return

        if Configuration.get_global_configuration().config_parser.has_option(
                section=cls.section,
                option=option):
            return Configuration.get_global_configuration().config_parser.getboolean(
                section=cls.section,
                option=option)
        else:
            return


class BaseSectionVersion(object):
    """The C{bsf.standards.BaseSectionVersion} class is the base class for a global configuration section and version.

    The defaults are read from the [{section}_{version}] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = None

    @classmethod
    def get_section(cls, version=None):
        """Get the the section defined by the sub-class section class variable and a version.

        @param version: Version
        @type version: None | str
        @return: Configuration section
        @rtype: str
        """
        if version:
            return '_'.join((cls.section, version))
        else:
            return cls.section

    @classmethod
    def get(cls, option=None, version=None):
        """Get the value for a configuration option in the section defined by the sub-class section class variable.

        @param option: Configuration option
        @type option: None | str | unicode
        @param version: Version
        @type version: None | str | unicode
        @return: Configuration value
        @rtype: None | str | unicode
        """
        if not cls.section:
            return

        if not option:
            return

        section = cls.get_section(version=version)
        if Configuration.get_global_configuration().config_parser.has_option(section=section, option=option):
            return Configuration.get_global_configuration().config_parser.get(section=section, option=option)
        else:
            return

    @classmethod
    def get_expanded_directory(cls, option=None, version=None):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e. on UNIX ~ or ~user and
        any environment variables i.e. on UNIX ${NAME} or $NAME.
        @param option: Configuration option
        @type option: None | str | unicode
        @param version: Version
        @type version: None | str
        @return: Expanded directory
        @rtype: None | str | unicode
        """
        if not option:
            return

        section = cls.get_section(version=version)
        if Configuration.get_global_configuration().config_parser.has_option(section=section, option=option):
            return Configuration.get_global_configuration().get_expanded_directory(section=section, option=option)
        else:
            return


class JavaClassPath(BaseSection):
    """The C{bsf.standards.JavaClassPath} class models Java class path defaults.

    The defaults are read from the [classpath] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
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


class EnsemblVEP(BaseSectionVersion):
    """The C{bsf.standards.EnsemblVEP} class models Ensembl Variant Effect Predictor (VEP) defaults.

    The defaults are read from the [ensembl_vep_{genome_version}] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = 'ensembl_vep'

    @classmethod
    def get_directory_cache(cls, genome_version=None):
        """Get the cache directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Cache directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_cache', version=genome_version)

    @classmethod
    def get_directory_fasta(cls, genome_version=None):
        """Get the FASTA directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: FASTA directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_fasta', version=genome_version)

    @classmethod
    def get_directory_plugin(cls, genome_version=None):
        """Get the plug-ins directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Plug-ins directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_plugin', version=genome_version)

    @classmethod
    def get_directory_source(cls, genome_version=None):
        """Get the source directory path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Source directory path
        @rtype: None | str | unicode
        """
        return cls.get_expanded_directory(option='directory_source', version=genome_version)

    @classmethod
    def get_name_assembly(cls, genome_version=None):
        """Get the genome assembly name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Genome assembly name
        @rtype: None | str | unicode
        """
        return cls.get(option='name_assembly', version=genome_version)

    @classmethod
    def get_name_species(cls, genome_version=None):
        """Get the scientific species name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Scientific species name
        @rtype: None | str | unicode
        """
        return cls.get(option='name_species', version=genome_version)

    @classmethod
    def get_sql_user(cls, genome_version=None):
        """Get the SQL database user name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database user name
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_user', version=genome_version)

    @classmethod
    def get_sql_pass(cls, genome_version=None):
        """Get the SQL database password.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database password
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_pass', version=genome_version)

    @classmethod
    def get_sql_host(cls, genome_version=None):
        """Get the SQL database host name.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database host name
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_host', version=genome_version)

    @classmethod
    def get_sql_port(cls, genome_version=None):
        """Get the SQL database TCP/IP port number.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: SQL database TCP/IP port number
        @rtype: None | str | unicode
        """
        return cls.get(option='sql_port', version=genome_version)

    @classmethod
    def get_ofc_path(cls, genome_version=None):
        """Get the output fields configuration (TSV) file path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Output fields configuration (TSV) file path
        @rtype: None | str | unicode
        """
        return cls.get(option='ofc_path', version=genome_version)

    @classmethod
    def get_soc_path(cls, genome_version=None):
        """Get the Sequence Ontology configuration (TSV) file path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Sequence Ontology configuration (TSV) file path
        @rtype: None | str | unicode
        """
        return cls.get(option='soc_path', version=genome_version)

    @classmethod
    def get_refseq_alignments_path(cls, genome_version=None):
        """Get the NCBI RefSeq alignments (BAM) file path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: NCBI RefSeq alignments (BAM) file path
        @rtype: None | str | unicode
        """
        return cls.get(option='refseq_alignments_path', version=genome_version)

    @classmethod
    def get_cadd_path(cls, genome_version=None):
        """Get the Combined Annotation Dependent Depletion (CADD) file path.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Combined Annotation Dependent Depletion (CADD) file path
        @rtype: None | str | unicode
        """
        return cls.get(option='cadd_path', version=genome_version)


class Genome(BaseSectionVersion):
    """The C{bsf.standards.Genome} class models Genome defaults.

    The defaults are read from the [genome_{genome_version}] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = 'genome'

    @classmethod
    def get_date(cls, genome_version=None):
        """Get the release date.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Release date
        @rtype: None | str | unicode
        """
        return cls.get(option='date', version=genome_version)

    @classmethod
    def get_fasta_suffix(cls, genome_version=None):
        """Get the FASTA suffix.

        The suffix could be 'fa' or 'fasta' as in "*.fa" or "*.fasta".
        THe NCBI uses fna for nucleotide FASTA files.
        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Description
        @rtype: None | str | unicode
        """
        return cls.get(option='fasta_suffix', version=genome_version)

    @classmethod
    def get_description(cls, genome_version=None):
        """Get the description.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Description
        @rtype: None | str | unicode
        """
        return cls.get(option='description', version=genome_version)

    @classmethod
    def get_provider(cls, genome_version=None):
        """Get the provider.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Provider
        @rtype: None | str | unicode
        """
        return cls.get(option='provider', version=genome_version)

    @classmethod
    def get_species(cls, genome_version=None):
        """Get the species.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Species
        @rtype: None | str | unicode
        """
        return cls.get(option='species', version=genome_version)

    @classmethod
    def get_ucsc(cls, genome_version=None):
        """Get the UCSC Genome Browser alias.

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: UCSC Genome Browser alias
        @rtype: None | str | unicode
        """
        return cls.get(option='ucsc', version=genome_version)

    @classmethod
    def get_uri(cls, genome_version=None):
        """Get the uniform resource identifier (URI).

        @param genome_version: Genome assembly version
        @type genome_version: None | str
        @return: Uniform resource identifier (URI)
        @rtype: None | str | unicode
        """
        return cls.get(option='uri', version=genome_version)

    @classmethod
    def resolve_ucsc_alias(cls, genome_version):
        """Resolve a genome version to an eventual UCSC Genome Browser-specific alias.

        If an alias has not been defined, the original genome version will be returned.
        @param genome_version: Genome version
        @type genome_version: str
        @return: UCSC Genome Browser assembly version
        @rtype: str
        """
        ucsc_version = cls.get_ucsc(genome_version=genome_version)

        if ucsc_version is None:
            return genome_version
        else:
            return ucsc_version


class Transcriptome(BaseSectionVersion):
    """The C{bsf.standards.Transcriptome} class models Transcriptome defaults.

    The defaults are read from the [transcriptome_{transcriptome_version}] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = 'transcriptome'

    @classmethod
    def get_date(cls, transcriptome_version=None):
        """Get the release date.

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: None | str
        @return: Release date
        @rtype: None | str | unicode
        """
        return cls.get(option='date', version=transcriptome_version)

    @classmethod
    def get_description(cls, transcriptome_version=None):
        """Get the description.

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: None | str
        @return: Description
        @rtype: None | str | unicode
        """
        return cls.get(option='description', version=transcriptome_version)

    @classmethod
    def get_genome(cls, transcriptome_version=None):
        """Get the genome version.

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: None | str
        @return: Genome version
        @rtype: None | str | unicode
        """
        return cls.get(option='genome', version=transcriptome_version)

    @classmethod
    def get_provider(cls, transcriptome_version=None):
        """Get the provider.

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: None | str
        @return: Provider
        @rtype: None | str | unicode
        """
        return cls.get(option='provider', version=transcriptome_version)

    @classmethod
    def get_species(cls, transcriptome_version=None):
        """Get the species.

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: None | str
        @return: Species
        @rtype: None | str | unicode
        """
        return cls.get(option='species', version=transcriptome_version)

    @classmethod
    def get_uri(cls, transcriptome_version=None):
        """Get the uniform resource identifier (URI).

        @param transcriptome_version: Transcriptome version
        @type transcriptome_version: None | str
        @return: Uniform resource identifier (URI)
        @rtype: None | str | unicode
        """
        return cls.get(option='uri', version=transcriptome_version)


class FilePath(BaseSection):
    """The C{bsf.standards.FilePath} class models file path defaults.

    The defaults are read from the [directories] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """

    section = 'directories'

    @classmethod
    def get_cache(cls):
        """Get the (absolute) cache directory path locally on the compute node (e.g. /dev/shm).

        @return: Cache directory path
        @rtype: None | str | unicode
        """
        return cls.get(option='cache')

    @classmethod
    def get_home(cls):
        """Get the (absolute) Home directory path.

        @return: Home directory path
        @rtype: None | str | unicode
        """
        return cls.get(option='home')

    @classmethod
    def _prepend_home(cls, absolute=True, file_path=None):
        """Private class method to prepend the I{home} directory.

        @param absolute: Absolute file path
        @type absolute: bool
        @param file_path: File path
        @type file_path: None | str | unicode
        @return: File path
        @rtype: None | str | unicode
        """
        if file_path is None:
            return

        if absolute and not os.path.isabs(file_path):
            return os.path.join(cls.get_home(), file_path)
        else:
            return file_path

    @classmethod
    def get_illumina_run(cls, absolute=True):
        """Get the Illumina Run Folder (IRF) directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Illumina Run Folder directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='illumina_run'))

    @classmethod
    def get_illumina_sav(cls, absolute=True):
        """Get the Illumina Sequence Analysis Viewer (SAV) directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Illumina Sequence Analysis Viewer directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='illumina_sav'))

    @classmethod
    def get_sequences(cls, absolute=True):
        """Get the Sequences directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Sequences directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='sequences'))

    @classmethod
    def get_samples(cls, absolute=True):
        """Get the Samples directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Samples directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='samples'))

    @classmethod
    def get_projects(cls, absolute=True):
        """Get the Analysis Projects directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Analysis projects directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='projects'))

    @classmethod
    def get_public_html(cls, absolute=True):
        """Get the Web Server directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Web server directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='public_html'))

    @classmethod
    def get_resource(cls, absolute=True):
        """Get the Resources directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Resources directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_home(absolute=absolute, file_path=cls.get(option='resources'))

    @classmethod
    def _prepend_resource(cls, absolute=True, file_path=None):
        """Private class method to prepend the I{resource} directory.

        @param absolute: Absolute file path
        @type absolute: bool
        @param file_path: File path
        @type file_path: None | str | unicode
        @return: File path
        @rtype: None | str | unicode
        """
        if file_path is None:
            return

        if absolute and not os.path.isabs(file_path):
            return os.path.join(cls.get_resource(absolute=absolute), file_path)
        else:
            return file_path

    @classmethod
    def get_resource_genome(cls, genome_version=None, absolute=True):
        """Get the Genome resource directory path.

        @param genome_version: The genome version (e.g. mm10, ...)
        @type genome_version: None | str
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Genome directory path
        @rtype: None | str | unicode
        """
        file_path = cls.get(option='genomes')

        if genome_version:
            file_path = os.path.join(file_path, genome_version)

        return cls._prepend_resource(absolute=absolute, file_path=file_path)

    @classmethod
    def get_resource_genome_index(cls, genome_version, genome_index=None, absolute=True):
        """Get a genome index resource directory path.

        In case the genome_index is not specified, the resource genome directory will be returned.
        @param genome_version: Genome version (e.g. mm10, ...)
        @type genome_version: str
        @param genome_index: Genome index (e.g. bowtie2, ...)
        @type genome_index: str | None
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Genome index resource directory path
        @rtype: None | str | unicode
        @raise Exception: Unknown genome index name
        """
        if genome_index is None:
            return cls.get_resource_genome(genome_version=genome_version, absolute=absolute)

        index_directory = Index.get(option=genome_index)

        if index_directory is None:
            raise Exception("Unknown genome index name '" + repr(genome_index) + "'.")
        else:
            return os.path.join(
                cls.get_resource_genome(genome_version=genome_version, absolute=absolute),
                index_directory)

    @classmethod
    def get_resource_genome_fasta(cls, genome_version, genome_index=None, absolute=True):
        """Get a genome FASTA resource file path.

        @param genome_version: Genome version (e.g. mm10, ...)
        @type genome_version: str
        @param genome_index: Genome index (e.g. bowtie2, ...)
        @type genome_index: str | None
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Genome FASTA resource file path
        @rtype: str | unicode
        @raise Exception: Unknown genome index name
        """
        fasta_suffix = Genome.get_fasta_suffix(genome_version=genome_version)
        if not fasta_suffix:
            fasta_suffix = 'fa'

        return os.path.join(
            cls.get_resource_genome_index(
                genome_version=genome_version,
                genome_index=genome_index,
                absolute=absolute),
            '.'.join((genome_version, fasta_suffix)))

    @classmethod
    def get_resource_genome_fasta_index(cls, genome_version, genome_index=None, absolute=True):
        """Get a genome FASTA index (*.fai) resource file path.

        @param genome_version: Genome version (e.g. mm10, ...)
        @type genome_version: str
        @param genome_index: Genome index (e.g. bowtie2, ...)
        @type genome_index: str | None
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Genome FASTA resource file path
        @rtype: str | unicode
        """
        return cls.get_resource_genome_fasta(
            genome_version=genome_version,
            genome_index=genome_index,
            absolute=absolute) + '.fai'

    @classmethod
    def get_resource_transcriptome(cls, transcriptome_version=None, absolute=True):
        """Get the transcriptome resource directory path.

        @param transcriptome_version: The transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: None | str
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Transcriptome resource directory path
        @rtype: None | str | unicode
        """
        file_path = cls.get(option='transcriptomes')

        if transcriptome_version:
            file_path = os.path.join(file_path, transcriptome_version)

        return cls._prepend_resource(absolute=absolute, file_path=file_path)

    @classmethod
    def get_resource_transcriptome_index(cls, transcriptome_version, transcriptome_index, absolute=True):
        """Get a transcriptome index resource directory path.

        @param transcriptome_version: Transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: str
        @param transcriptome_index: Transcriptome index (e.g. star, tophat, ...)
        @type transcriptome_index: str
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Transcriptome index resource directory path
        @rtype: None | str | unicode
        @raise Exception: Unknown transcriptome index name
        """
        index_directory = Index.get(option=transcriptome_index)

        if index_directory is None:
            raise Exception("Unknown transcriptome index name '" + repr(transcriptome_index) + "'.")
        else:
            return os.path.join(
                cls.get_resource_transcriptome(transcriptome_version=transcriptome_version, absolute=absolute),
                index_directory)

    @classmethod
    def get_resource_transcriptome_gtf(cls, transcriptome_version, transcriptome_index, basic=True, absolute=True):
        """Get a transcriptome GTF resource file path.

        @param transcriptome_version: Transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: str
        @param transcriptome_index: Transcriptome index (e.g. star, tophat, ...)
        @type transcriptome_index: str
        @param basic: Basic transcriptome
        @type basic: bool
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Transcriptome GTF resource file path
        @rtype: None | str | unicode
        @raise Exception: Unknown transcriptome index name
        """
        if basic:
            file_name = transcriptome_version + '.basic.gtf'
        else:
            file_name = transcriptome_version + '.gtf'

        return os.path.join(
            cls.get_resource_transcriptome_index(
                transcriptome_version=transcriptome_version,
                transcriptome_index=transcriptome_index,
                absolute=absolute),
            file_name)

    @classmethod
    def get_resource_transcriptome_txdb(cls, transcriptome_version, transcriptome_index, basic=True, absolute=True):
        """Get a transcriptome TxDb resource file path.

        @param transcriptome_version: Transcriptome version (e.g. mm10_e87, ...)
        @type transcriptome_version: str
        @param transcriptome_index: Transcriptome index (e.g. star, tophat, ...)
        @type transcriptome_index: str
        @param basic: Basic transcriptome
        @type basic: bool
        @param absolute: Absolute file path
        @type absolute: bool
        @return: Transcriptome TxDb resource file path
        @rtype: None | str | unicode
        @raise Exception: Unknown transcriptome index name
        """
        if basic:
            file_name = transcriptome_version + '.basic.sqlite'
        else:
            file_name = transcriptome_version + '.sqlite'

        return os.path.join(
            cls.get_resource_transcriptome_index(
                transcriptome_version=transcriptome_version,
                transcriptome_index=transcriptome_index,
                absolute=absolute),
            file_name)

    @classmethod
    def get_resource_gatk_bundle(cls, gatk_bundle_version=None, genome_version=None, absolute=True):
        """Get the GATK Bundle resource directory path.

        @param gatk_bundle_version: The GATK bundle version
        @type gatk_bundle_version: None | str
        @param genome_version: The genome version (e.g. b37, ...)
        @type genome_version: None | str
        @param absolute: Absolute file path
        @type absolute: bool
        @return: GATK Bundle resource directory path
        @rtype: None | str | unicode
        """
        file_path = cls.get(option='gatk_bundle')

        if gatk_bundle_version:
            file_path = os.path.join(file_path, gatk_bundle_version)

        if genome_version:
            file_path = os.path.join(file_path, genome_version)

        return cls._prepend_resource(absolute=absolute, file_path=file_path)

    @classmethod
    def get_resource_intervals(cls, absolute=True):
        """Get the Target Intervals resource directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Target Intervals resource directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_resource(absolute=absolute, file_path=cls.get(option='intervals'))

    @classmethod
    def get_resource_cadd(cls, absolute=True):
        """Get the Combined Annotation Dependent Depletion (CADD) resource directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Combined Annotation Dependent Depletion (CADD) resource directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_resource(absolute=absolute, file_path=cls.get(option='cadd'))

    @classmethod
    def get_resource_cosmic(cls, absolute=True):
        """Get the Catalogue Of Somatic Mutations In Cancer (COSMIC) resource directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: Catalogue Of Somatic Mutations In Cancer (COSMIC) resource directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_resource(absolute=absolute, file_path=cls.get(option='cosmic'))

    @classmethod
    def get_resource_snpeff_data(cls, absolute=True):
        """Get the snpEff Data resource directory path.

        @param absolute: Absolute file path
        @type absolute: bool
        @return: snpEff Data resource directory path
        @rtype: None | str | unicode
        """
        return cls._prepend_resource(absolute=absolute, file_path=cls.get(option='snpeff_data'))


class Index(BaseSection):
    """The C{bsf.standards.Index} class models genome or transcriptome index directory defaults.

    The defaults are read from the [indices] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = 'indices'


class Operator(BaseSection):
    """The C{bsf.standards.Operator} class models operator defaults.

    The defaults are read from the [operator] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
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


class UCSC(BaseSection):
    """The C{bsf.standards.UCSC} class models UCSC Genome Browser uniform resource locator (URL) defaults.

    The defaults are read from the [ucsc] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
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


class URL(BaseSection):
    """The C{bsf.standards.URL} class models web server uniform resource locator (URL) defaults.

    The defaults are read from the [url] section of the global configuration file.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
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


class VendorQualityFilter(BaseSection):
    """The C{bsf.standards.VendorQualityFilter} class models (SAM) Vendor Quality Filter defaults.

    The defaults are read from the [VendorQualityFilter] section of the global configuration file.
    For each flow cell type a boolean specifies whether vendor quality filtering should be applied or not.
    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """

    section = 'vendor_quality_filter'

    @classmethod
    def get_vendor_quality_filter(cls, flow_cell_type=None):
        """Get the vendor quality filter setting for a particular flow cell type.

        The flow cell type is accessible via property C{bsf.illumina.RunParameters.get_flow_cell_type},
        practically directly via C{bsf.illumina.RunFolder.run_parameters.get_flow_cell_type}.
        @param flow_cell_type: FLow cell (chemistry) type
        @type flow_cell_type: str | unicode | None
        @return: Vendor quality filer setting
        @rtype: bool
        """
        if not flow_cell_type:
            return

        # To avoid emitting cryptic error messages, check for the presence of the flow cell type before
        # and inform the user about an eventual problem with a hopefully more meaningful message.

        if not Configuration.get_global_configuration().config_parser.has_option(
                section=cls.section,
                option=flow_cell_type):
            raise Exception('Flow cell type ' + repr(flow_cell_type) +
                            ' is not defined in the ' + repr(cls.section) +
                            ' section of the standard configuration file ' +
                            repr(Configuration.global_file_path) + '\n')

        return cls.getboolean(option=flow_cell_type)


class Secrets(BaseSection):
    """The C{bsf.standards.Secrets} class models file paths to configuration files with secrets.

    Attributes:
    @cvar section: C{configparser.ConfigParser} section
    @type section: str | unicode
    """
    section = 'secrets'

    user_mask = stat.S_IRWXG | stat.S_IRWXO

    @classmethod
    def get_azure_file_path(cls):
        """Get the configuration file path with Microsoft Azure secrets.

        Also checks that the file is only readable by the user and not accessible for group and other.
        @return: File path to configuration file with secrets.
        @rtype: str | None
        """
        file_path = os.path.normpath(os.path.expandvars(os.path.expanduser(cls.get(option='azure_file_path'))))

        if not file_path:
            return None

        file_stat_result = os.stat(path=file_path, follow_symlinks=True)

        if file_stat_result.st_mode & cls.user_mask:
            raise Exception(
                'Secrets configuration file {} has file mode {:#0o}, '
                'but should obey user mask {:#0o}.'.format(
                    file_path, file_stat_result.st_mode, cls.user_mask))

        return file_path
