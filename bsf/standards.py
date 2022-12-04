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
"""The :py:mod:`bsf.standards` module provides classes modelling configuration and default information.
"""
import os
import re
import stat
from configparser import ConfigParser
from typing import Dict, List, Optional
from xml.etree.ElementTree import ElementTree


def get_irf_path(name):
    """Return the absolute file path for an Illumina Run Folder (IRF) name.

    This function first checks for existence in :py:meth:`bsf.standards.StandardFilePath.get_illumina_run`, before
    checking in :py:meth:`bsf.standards.StandardFilePath.get_illumina_sav`.

    :param name: An Illumina Run Folder (IRF) name.
    :type name: str
    :return: An absolute file path.
    :rtype: str | None
    """
    # Check the Illumina Run Folder directory.
    file_path = Configuration.get_absolute_path(
        file_path=name,
        default_path=StandardFilePath.get_illumina_run(absolute=True))
    if os.path.exists(file_path):
        return file_path

    # Check the Illumina Sequence Analysis Viewer directory.
    file_path = Configuration.get_absolute_path(
        file_path=name,
        default_path=StandardFilePath.get_illumina_sav(absolute=True))
    if os.path.exists(file_path):
        return file_path

    # To avoid clashes, Illumina Run Folders get a '_sav' suffix upon archiving,
    # so the name as such should not exist in the above
    # Illumina Sequence Analysis Viewer directory.
    # Now, append the '_sav' suffix customary for SAV folders and check once more.
    file_path += '_sav'
    if os.path.exists(file_path):
        return file_path

    return


class SafeFileName(object):
    """The :py:class:`bsf.standards.SafeFileName` class represents a regular expression pattern
    to make file names safe.
    """
    _re_pattern: Optional[re.Pattern] = None

    @classmethod
    def get_safe_file_name(cls, file_name):
        """Get a safe file name by replacing special characters with underscore characters.

        This function is modelled after the `HTSJDK <http://samtools.github.io/htsjdk/>`_
        :literal:`makeFileNameSafe()` method of the
        :literal:`htsjdk.samtools.util.IOUtil` class and uses the following pattern.

        :literal:`[\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]`

        :param file_name: File name
        :type file_name: str
        :return: Safe file name
        :rtype: str
        """
        # https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/util/IOUtil.java#L779
        if not cls._re_pattern:
            cls._re_pattern = re.compile(pattern="[\\s!\"#$%&'()*/:;<=>?@\\[\\]\\\\^`{|}~]")

        return re.sub(pattern=cls._re_pattern, repl='_', string=file_name)


class Configuration(object):
    """The :py:class:`bsf.standards.Configuration` class represents one or more UNIX-style initialisation
    (:literal:`*.ini`) files.

    A :py:class:`bsf.standards.Configuration` object has an associated Python :py:class:`configparser.ConfigParser`
    object to parse the file(s).

    :ivar file_path_list: A Python :py:class:`list` object of
        Python :py:class:`str` (configuration file path) objects.
    :type file_path_list: list[str]
    :ivar config_parser: A Python :py:class:`configparser.ConfigParser` object.
    :type config_parser: ConfigParser
    """

    _global_configuration = None
    _global_environment = 'BSF_PYTHON_INI'
    _global_file_path = '~/.bsfpython.ini'

    @classmethod
    def get_global_file_path(cls):
        """Get the global UNIX-style initialisation (:literal:`*.ini`) file path.

        The global configuration file path is based on the value of
        environment variable :literal:`BSF_PYTHON_INI` if defined,
        and defaults to :literal:`~/.bsfpython.ini` otherwise.

        :return: Global INI configuration file path
        :rtype: str
        """
        if cls._global_environment in os.environ and os.environ[cls._global_environment]:
            file_path = os.environ[cls._global_environment]
        else:
            file_path = cls._global_file_path

        file_path = os.path.expanduser(file_path)
        file_path = os.path.expandvars(file_path)
        file_path = os.path.normpath(file_path)

        return file_path

    @classmethod
    def get_global_configuration(cls):
        """Get a global :py:class:`bsf.standards.Configuration` object and initialise it, if not already done so.

        :return: A global :py:class:`bsf.standards.Configuration` object.
        :rtype: Configuration
        """
        if cls._global_configuration is None:
            cls._global_configuration = cls.from_file_path_list(
                file_path_list=[cls.get_global_file_path()])

        return cls._global_configuration

    @staticmethod
    def get_absolute_path(file_path=None, default_path=None):
        """Return an absolute file path.

        Expand an eventual user part (i.e., on UNIX ~ or ~user) and
        expand any environment variables (i.e., on UNIX ${NAME} or $NAME).
        Check if an absolute path has been provided, if not,
        automatically prepend default directory paths,
        which again has user part and environment variables expanded.
        Finally, normalise the file path.

        :param file_path: A file path.
        :type file_path: str | None
        :param default_path: A default absolute path.
        :type default_path: str | None
        :return: An absolute file path.
        :rtype: str | None
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
        """Get A configuration section Python :py:class:`str` object composed of the
        Python module and Python class name of an object instance.

        :param instance: A Python instance (i.e., object).
        :type instance: object
        :return: A configuration file section string.
        :rtype: str
        """
        # For Python "type" instances the "__name__" instance variable provides the name of the type, while for
        # Python "object" instances the "__class__" variable provides the Python "type" object.

        if isinstance(instance, type):
            return '.'.join((instance.__module__, instance.__name__))
        else:
            return '.'.join((instance.__module__, instance.__class__.__name__))

    @classmethod
    def from_file_path_list(cls, file_path_list):
        """Create a :py:class:`bsf.standards.Configuration` object based on a
        Python :py:class:`list` object of Python :py:class:`str` (file path) objects.

        Both, user and variable expansion gets applied to each file path.
        Identical files are read only once.

        :param file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (configuration file path) objects.
        :type file_path_list: list[str]
        :return: A :py:class:`bsf.standards.Configuration` object.
        :rtype: Configuration
        :raise Exception: If the configuration file path does not exist.
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
        # they have to be case-sensitive.
        # Hence, override method optionxform() with str().

        config_parser = ConfigParser()
        config_parser.optionxform = str

        configuration = cls(file_path_list=file_path_list, config_parser=config_parser)

        configuration._config_path_list = configuration.config_parser.read(filenames=configuration.file_path_list)

        if len(configuration._config_path_list) == 0:
            raise Exception(f'None of the configuration files exist:\n'
                            f'{configuration.file_path_list!r}')

        return configuration

    def __init__(self, file_path_list=None, config_parser=None):
        """Initialise a :py:class:`bsf.standards.Configuration` object.

        :param file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (configuration file path) objects.
        :type file_path_list: list[str] | None
        :param config_parser: A Python :py:class:`configparser.ConfigParser` object.
        :type config_parser: ConfigParser | None
        """
        super(Configuration, self).__init__()

        if file_path_list is None:
            self.file_path_list = list()
        else:
            self.file_path_list = file_path_list

        if config_parser is None:
            self.config_parser = ConfigParser()
        else:
            self.config_parser = config_parser

        self._config_path_list: Optional[List[str]] = None

        return

    def trace(self, level):
        """Trace a :py:class:`bsf.standards.Configuration` object.

        :param level: Indentation level
        :type level: int
        :return: Trace information.
        :rtype: list[str]
        """
        indent = '  ' * level

        str_list: List[str] = list()

        str_list.append('{}{!r}\n'.format(indent, self))
        str_list.append('{}  file_path_list: {!r}\n'.format(indent, self.file_path_list))
        str_list.append('{}  config_parser:  {!r}\n'.format(indent, self.config_parser))

        return str_list

    def get_expanded_directory(self, section, option):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e., on UNIX ~ or ~user and
        any environment variables i.e., on UNIX ${NAME} or $NAME.

        :param section: A configuration section string.
        :type section: str
        :param option: A configuration option string.
        :type option: str
        :return: An expanded directory path.
        :rtype: str | None
        """
        return self.get_absolute_path(file_path=self.config_parser.get(section=section, option=option))

    @staticmethod
    def list_from_csv(csv_string):
        """Convert a comma-separated Python :py:class:`str` object into a
        Python :py:class:`list` object of Python :py:class:`str` objects.

        All elements are stripped and only non-empty elements are appended to the list.

        :param csv_string: A Python :py:class:`str` object of comma-separated values.
        :type csv_string: str
        :return: A Python :py:class:`list` object of stripped Python :py:class:`str` objects.
        :rtype: list[str]
        """
        return [x.strip() for x in csv_string.split(',') if x.strip()]

    def get_list_from_csv(self, section, option):
        """Convert a comma-separated Python :py:class:`str` object into a
        Python :py:class:`list` object of Python :py:class:`str` objects.

        All elements are stripped and only non-empty elements are appended to the list.

        :param section: A configuration section string.
        :type section: str
        :param option: A configuration option string.
        :type option: str
        :return: A Python :py:class:`list` object of Python :py:class:`str` objects.
        :rtype: list[str] | None
        """
        csv_string = self.config_parser.get(section=section, option=option)
        if csv_string:
            return self.list_from_csv(csv_string=csv_string)


class BaseSection(object):
    """The :py:class:`bsf.standards.BaseSection` class is the base class for a global configuration section.

    The defaults are read from the :literal:`[{section}]` section of the global configuration file.
    """
    _section = None

    @classmethod
    def get(cls, option):
        """Get a value for a configuration option in a configuration section defined by a
        :py:attr:`bsf.standards.BaseSection.section` class variable.

        This method is a re-implementation of the :py:meth:`configparser.ConfigParser.get` method
        that returns :literal:`None` upon non-existing sections or options.

        :param option: A configuration option.
        :type option: str
        :return: A configuration option value.
        :rtype: str | None
        """
        if not cls._section:
            return

        if not option:
            return

        if Configuration.get_global_configuration().config_parser.has_option(
                section=cls._section,
                option=option):
            return Configuration.get_global_configuration().config_parser.get(
                section=cls._section,
                option=option)
        else:
            return

    @classmethod
    def getboolean(cls, option):
        """Get a value for a configuration option in a configuration section defined by a
        :py:attr:`bsf.standards.BaseSection.section` class variable.

        This method is a re-implementation of the :py:meth:`configparser.ConfigParser.getboolean` method
        that returns :literal:`None` upon non-existing sections or options.

        :param option: A configuration option.
        :type option: str
        :return: A configuration option value.
        :rtype: bool | None
        """
        if not cls._section:
            return

        if not option:
            return

        if Configuration.get_global_configuration().config_parser.has_option(
                section=cls._section,
                option=option):
            return Configuration.get_global_configuration().config_parser.getboolean(
                section=cls._section,
                option=option)
        else:
            return

    @classmethod
    def get_expanded_directory(cls, option):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e., on UNIX ~ or ~user and
        any environment variables i.e., on UNIX ${NAME} or $NAME.

        :param option: A configuration option.
        :type option: str
        :return: An expanded directory path.
        :rtype: str | None
        """
        if not cls._section:
            return

        if not option:
            return

        if Configuration.get_global_configuration().config_parser.has_option(section=cls._section, option=option):
            return Configuration.get_global_configuration().get_expanded_directory(section=cls._section, option=option)
        else:
            return


class BaseSectionVersion(object):
    """The :py:class:`bsf.standards.BaseSectionVersion` class is the base class for a global configuration section
    and version.

    The defaults are read from the :literal:`[{section}_{version}]` section of the global configuration file.
    """
    _section = None

    @classmethod
    def get_section(cls, version):
        """Get a configuration section defined by a
        :py:attr:`bsf.standards.BaseSectionVersion.section` class variable and a version.

        :param version: A version.
        :type version: str
        :return: A configuration section.
        :rtype: str
        """
        if version:
            return '_'.join((cls._section, version))
        else:
            return cls._section

    @classmethod
    def get(cls, option, version):
        """Get a value for a configuration option in a configuration section defined by a
        :py:attr:`bsf.standards.BaseSectionVersion.section` class variable.

        :param option: A configuration option.
        :type option: str
        :param version: A version.
        :type version: str
        :return: A configuration value.
        :rtype: str | None
        """
        if not cls._section:
            return

        if not option:
            return

        section = cls.get_section(version=version)
        if Configuration.get_global_configuration().config_parser.has_option(section=section, option=option):
            return Configuration.get_global_configuration().config_parser.get(section=section, option=option)
        else:
            return

    @classmethod
    def get_expanded_directory(cls, option, version):
        """Get configuration information for a directory and expand it.

        The expansion includes an eventual user part i.e., on UNIX ~ or ~user and
        any environment variables i.e., on UNIX ${NAME} or $NAME.

        :param option: A configuration option.
        :type option: str
        :param version: A version.
        :type version: str
        :return: An expanded directory path.
        :rtype: str | None
        """
        if not option:
            return

        section = cls.get_section(version=version)
        if Configuration.get_global_configuration().config_parser.has_option(section=section, option=option):
            return Configuration.get_global_configuration().get_expanded_directory(section=section, option=option)
        else:
            return


class JavaArchive(BaseSection):
    """The :py:class:`bsf.standards.JavaArchive` class models Java Archive (JAR) defaults.

    The defaults are read from the :literal:`[java_archive]` section of the global configuration file.
    """

    _section = 'java_archive'

    @classmethod
    def get_fgbio(cls):
        """Get a Fulcrum Genomics (fgbio) Java Archive (JAR) file path.

        :return: A Fulcrum Genomics (fgbio) Java Archive (JAR) file path.
        :rtype: str | None
        """
        return cls.get(option='fgbio')

    @classmethod
    def get_gatk(cls):
        """Get a GATK Java Archive (JAR) file path.

        :return: A GATK Java Archive (JAR) file path.
        :rtype: str | None
        """
        return cls.get(option='gatk')

    @classmethod
    def get_picard(cls):
        """Get a Picard Java Archive (JAR) file path.

        :return: A Picard Java Archive (JAR) file path.
        :rtype: str | None
        """
        return cls.get(option='picard')

    @classmethod
    def get_snpeff(cls):
        """Get a snpEff Java Archive (JAR) file path.

        :return: A snpEff Java Archive (JAR) file path.
        :rtype: str | None
        """
        return cls.get(option='snpeff')

    @classmethod
    def get_trimmomatic(cls):
        """Get a Trimmomatic Java Archive (JAR) file path.

        :return: A Trimmomatic Java Archive (JAR) file path.
        :rtype: str | None
        """
        return cls.get(option='trimmomatic')

    @classmethod
    def get_vcf_filter(cls):
        """Get a VCF.Filter Java Archive (JAR) file path.

        :return: A VCF.Filter Java Archive (JAR) file path.
        :rtype: str | None
        """
        return cls.get(option='vcf_filter')


class JavaClassPath(BaseSection):
    """The :py:class:`bsf.standards.JavaClassPath` class models Java Class Path directory defaults.

    The defaults are read from the :literal:`[java_classpath]` section of the global configuration file.
    """
    _section = 'java_classpath'

    @classmethod
    def get_illumina2bam(cls):
        """Get an Illumina2bam tools Java Class Path directory.

        :return: An Illumina2bam tools Java Class Path directory.
        :rtype: str | None
        """
        return cls.get(option='illumina2bam')


class EnsemblVEP(BaseSectionVersion):
    """The :py:class:`bsf.standards.EnsemblVEP` class models Ensembl Variant Effect Predictor (VEP) defaults.

    The defaults are read from the :literal:`[ensembl_vep_{genome_version}]` section of the global configuration file.
    """
    _section = 'ensembl_vep'

    @classmethod
    def get_directory_cache(cls, genome_version):
        """Get a cache directory path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A cache directory path.
        :rtype: str | None
        """
        return cls.get_expanded_directory(option='directory_cache', version=genome_version)

    @classmethod
    def get_directory_fasta(cls, genome_version):
        """Get a FASTA directory path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A FASTA directory path.
        :rtype: str | None
        """
        return cls.get_expanded_directory(option='directory_fasta', version=genome_version)

    @classmethod
    def get_directory_plugin(cls, genome_version):
        """Get a plug-ins directory path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A plug-ins directory path.
        :rtype: str | None
        """
        return cls.get_expanded_directory(option='directory_plugin', version=genome_version)

    @classmethod
    def get_directory_source(cls, genome_version):
        """Get a source directory path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A source directory path.
        :rtype: str | None
        """
        return cls.get_expanded_directory(option='directory_source', version=genome_version)

    @classmethod
    def get_name_assembly(cls, genome_version):
        """Get a genome assembly name for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A genome assembly name.
        :rtype: str | None
        """
        return cls.get(option='name_assembly', version=genome_version)

    @classmethod
    def get_name_species(cls, genome_version):
        """Get a scientific species name for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A scientific species name.
        :rtype: str | None
        """
        return cls.get(option='name_species', version=genome_version)

    @classmethod
    def get_sql_user(cls, genome_version):
        """Get an SQL database username for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: An SQL database username.
        :rtype: str | None
        """
        return cls.get(option='sql_user', version=genome_version)

    @classmethod
    def get_sql_pass(cls, genome_version):
        """Get an SQL database password for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: An SQL database password.
        :rtype: str | None
        """
        return cls.get(option='sql_pass', version=genome_version)

    @classmethod
    def get_sql_host(cls, genome_version):
        """Get an SQL database host name for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: An SQL database host name.
        :rtype: str | None
        """
        return cls.get(option='sql_host', version=genome_version)

    @classmethod
    def get_sql_port(cls, genome_version):
        """Get an SQL database TCP/IP port number for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: An SQL database TCP/IP port number.
        :rtype: str | None
        """
        return cls.get(option='sql_port', version=genome_version)

    @classmethod
    def get_ofc_path(cls, genome_version):
        """Get an output fields configuration (TSV) file path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: An output fields configuration (TSV) file path.
        :rtype: str | None
        """
        return cls.get(option='ofc_path', version=genome_version)

    @classmethod
    def get_soc_path(cls, genome_version):
        """Get a Sequence Ontology configuration (TSV) file path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A Sequence Ontology configuration (TSV) file path.
        :rtype: str | None
        """
        return cls.get(option='soc_path', version=genome_version)

    @classmethod
    def get_refseq_alignments_path(cls, genome_version):
        """Get an NCBI RefSeq alignments (BAM) file path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: An NCBI RefSeq alignments (BAM) file path.
        :rtype: str | None
        """
        return cls.get(option='refseq_alignments_path', version=genome_version)

    @classmethod
    def get_cadd_path(cls, genome_version):
        """Get a Combined Annotation Dependent Depletion (CADD) file path for an Ensembl VEP genome version.

        :param genome_version: A VEP genome assembly version.
        :type genome_version: str
        :return: A Combined Annotation Dependent Depletion (CADD) file path.
        :rtype: str | None
        """
        return cls.get(option='cadd_path', version=genome_version)


class Genome(BaseSectionVersion):
    """The :py:class:`bsf.standards.Genome` class models Genome defaults.

    The defaults are read from the :literal:`[genome_{genome_version}]` section of the global configuration file.
    """
    _section = 'genome'

    @classmethod
    def get_black_list(cls, genome_version):
        """Get an (ENCODE) black list file path of problematic regions.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: An (ENCODE) Black list file path.
        :rtype: str | None
        """
        return cls.get(option='black_list', version=genome_version)

    @classmethod
    def get_date(cls, genome_version):
        """Get a release date.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A release date.
        :rtype: str | None
        """
        return cls.get(option='date', version=genome_version)

    @classmethod
    def get_effective_size(cls, genome_version):
        """Get an effective genome size.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: An effective genome size.
        :rtype: str | None
        """
        return cls.get(option='effective_size', version=genome_version)

    @classmethod
    def get_fasta_suffix(cls, genome_version):
        """Get a FASTA suffix.

        The suffix could be :literal:`fa` or :literal:`fasta` as in :literal:`*.fa` or :literal:`*.fasta`.
        The NCBI uses a :literal:`faa` or :literal:`fna` suffix for FASTA files of amino-acids or nucleic-acids,
        respectively.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A FASTA suffix.
        :rtype: str | None
        """
        return cls.get(option='fasta_suffix', version=genome_version)

    @classmethod
    def get_description(cls, genome_version):
        """Get A description.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A description.
        :rtype: str | None
        """
        return cls.get(option='description', version=genome_version)

    @classmethod
    def get_provider(cls, genome_version):
        """Get a provider.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A provider.
        :rtype: str | None
        """
        return cls.get(option='provider', version=genome_version)

    @classmethod
    def get_species(cls, genome_version):
        """Get a species.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A species.
        :rtype: str | None
        """
        return cls.get(option='species', version=genome_version)

    @classmethod
    def get_ucsc(cls, genome_version):
        """Get a UCSC Genome Browser alias.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A UCSC Genome Browser alias.
        :rtype: str | None
        """
        return cls.get(option='ucsc', version=genome_version)

    @classmethod
    def get_uri(cls, genome_version):
        """Get a Uniform Resource Identifier (URI).

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A Uniform Resource Identifier (URI).
        :rtype: str | None
        """
        return cls.get(option='uri', version=genome_version)

    @classmethod
    def resolve_ucsc_alias(cls, genome_version):
        """Resolve a genome version to an eventual UCSC Genome Browser-specific alias.

        If an alias has not been defined, the original genome version will be returned.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A UCSC Genome Browser assembly version.
        :rtype: str
        """
        ucsc_version = cls.get_ucsc(genome_version=genome_version)

        if ucsc_version is None:
            return genome_version
        else:
            return ucsc_version


class SnpEff(BaseSectionVersion):
    """The :py:class:`bsf.standards.SnpEff` class models snpEff defaults.

    The defaults are read from the :literal:`[snpeff_{genome_version}]` section of the global configuration file.
    """
    _section = 'snpeff'

    @classmethod
    def get_genome_version(cls, genome_version):
        """Get A snpEff genome version.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A snpEff genome version.
        :rtype: str | None
        """
        return cls.get(option='genome_version', version=genome_version)


class Transcriptome(BaseSectionVersion):
    """The :py:class:`bsf.standards.Transcriptome` class models Transcriptome defaults.

    The defaults are read from the :literal:`[transcriptome_{transcriptome_version}]` section of the
    global configuration file.
    """
    _section = 'transcriptome'

    @classmethod
    def get_date(cls, transcriptome_version):
        """Get a release date.

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str
        :return: A release date.
        :rtype: str | None
        """
        return cls.get(option='date', version=transcriptome_version)

    @classmethod
    def get_description(cls, transcriptome_version):
        """Get a description.

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str
        :return: A description.
        :rtype: str | None
        """
        return cls.get(option='description', version=transcriptome_version)

    @classmethod
    def get_genome(cls, transcriptome_version):
        """Get a genome version.

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str
        :return: A genome version.
        :rtype: str | None
        """
        return cls.get(option='genome', version=transcriptome_version)

    @classmethod
    def get_provider(cls, transcriptome_version):
        """Get a provider.

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str
        :return: A provider.
        :rtype: str | None
        """
        return cls.get(option='provider', version=transcriptome_version)

    @classmethod
    def get_species(cls, transcriptome_version):
        """Get a species.

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str
        :return: A species.
        :rtype: str | None
        """
        return cls.get(option='species', version=transcriptome_version)

    @classmethod
    def get_uri(cls, transcriptome_version):
        """Get a Uniform Resource Identifier (URI).

        :param transcriptome_version: A transcriptome version.
        :type transcriptome_version: str
        :return: A Uniform Resource Identifier (URI).
        :rtype: str | None
        """
        return cls.get(option='uri', version=transcriptome_version)


class StandardFilePath(BaseSection):
    """The :py:class:`bsf.standards.StandardFilePath` class models file path defaults.

    The defaults are read from the :literal:`[directories]` section of the global configuration file.
    """

    _section = 'directories'

    @classmethod
    def get_cache(cls):
        """Get a cache directory path locally on a compute-node (e.g., /dev/shm).

        :return: A cache directory path.
        :rtype: str | None
        """
        return cls.get_expanded_directory(option='cache')

    @classmethod
    def get_home(cls):
        """Get a home directory path.

        :return: A home directory path.
        :rtype: str | None
        """
        return cls.get_expanded_directory(option='home')

    @classmethod
    def _prepend_home(cls, file_path, absolute=True):
        """Private class method to prepend a :literal:`home` directory path.

        :param file_path: A file path.
        :type file_path: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A :literal:`home` file path.
        :rtype: str | None
        """
        if file_path is None:
            return

        if absolute and not os.path.isabs(file_path):
            return os.path.join(cls.get_home(), file_path)
        else:
            return file_path

    @classmethod
    def get_illumina_run(cls, absolute=True):
        """Get an Illumina Run Folder (IRF) directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: An Illumina Run Folder directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='illumina_run'), absolute=absolute)

    @classmethod
    def get_illumina_sav(cls, absolute=True):
        """Get an Illumina Sequence Analysis Viewer (SAV) directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: An Illumina Sequence Analysis Viewer (SAV) directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='illumina_sav'), absolute=absolute)

    @classmethod
    def get_sequences(cls, absolute=True):
        """Get a sequence directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A sequence directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='sequences'), absolute=absolute)

    @classmethod
    def get_samples(cls, absolute=True):
        """Get a sample directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A sample directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='samples'), absolute=absolute)

    @classmethod
    def get_projects(cls, absolute=True):
        """Get a project directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A project directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='projects'), absolute=absolute)

    @classmethod
    def get_public_html(cls, absolute=True):
        """Get a web server directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A web server directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='public_html'), absolute=absolute)

    @classmethod
    def get_template_documents(cls, absolute=True):
        """Get a template documents directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A template documents directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='template_documents'), absolute=absolute)

    @classmethod
    def get_template_scripts(cls, absolute=True):
        """Get a template script directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A template script directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='template_scripts'), absolute=absolute)

    @classmethod
    def get_resource(cls, absolute=True):
        """Get a resource directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A resource directory path.
        :rtype: str | None
        """
        return cls._prepend_home(file_path=cls.get_expanded_directory(option='resources'), absolute=absolute)

    @classmethod
    def _prepend_resource(cls, file_path, absolute=True):
        """Private class method to prepend a :literal:`resource` directory path.

        :param file_path: A file path.
        :type file_path: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A :literal:`resource` file path.
        :rtype: str | None
        """
        if file_path is None:
            return

        if absolute and not os.path.isabs(file_path):
            return os.path.join(cls.get_resource(absolute=absolute), file_path)
        else:
            return file_path

    @classmethod
    def get_resource_genome(cls, genome_version, absolute=True):
        """Get a genome resource directory path.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A genome resource directory path.
        :rtype: str | None
        """
        file_path = cls.get_expanded_directory(option='genomes')

        if genome_version:
            file_path = os.path.join(file_path, genome_version)

        return cls._prepend_resource(file_path=file_path, absolute=absolute)

    @classmethod
    def get_resource_genome_black_list(cls, genome_version):
        """Get a genome black list resource file path.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :return: A genome black list resource file path.
        :rtype: str | None
        """
        black_list_file_path = Genome.get_black_list(genome_version=genome_version)

        if not black_list_file_path:
            return None

        if os.path.isabs(black_list_file_path):
            return black_list_file_path
        else:
            return os.path.join(
                cls.get_resource_genome(genome_version=genome_version, absolute=True),
                black_list_file_path)

    @classmethod
    def get_resource_genome_index(cls, genome_version, genome_index=None, absolute=True):
        """Get a genome index resource directory path.

        In case the genome_index is not specified, the resource genome directory will be returned.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :param genome_index: A genome index (e.g., bowtie2, ...).
        :type genome_index: str | None
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A genome index resource directory path.
        :rtype: str | None
        :raise Exception: If the genome index name is unknown.
        """
        if genome_index is None:
            return cls.get_resource_genome(genome_version=genome_version, absolute=absolute)

        index_directory = Index.get(option=genome_index)

        if index_directory is None:
            raise Exception(f'Unknown genome index name {genome_index!r}.')
        else:
            return os.path.join(
                cls.get_resource_genome(genome_version=genome_version, absolute=absolute),
                index_directory)

    @classmethod
    def get_resource_genome_fasta(cls, genome_version, genome_index=None, absolute=True):
        """Get a genome FASTA resource file path.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :param genome_index: A genome index (e.g., bowtie2, ...).
        :type genome_index: str | None
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A genome FASTA resource file path.
        :rtype: str
        :raise Exception: If the genome index name is unknown.
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
        """Get a genome FASTA index (:literal:`*.fai`) resource file path.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :param genome_index: A genome index (e.g., bowtie2, ...).
        :type genome_index: str | None
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A genome FASTA index resource file path.
        :rtype: str
        """
        return cls.get_resource_genome_fasta(
            genome_version=genome_version,
            genome_index=genome_index,
            absolute=absolute) + '.fai'

    @classmethod
    def get_resource_transcriptome(cls, transcriptome_version, absolute=True):
        """Get a transcriptome resource directory path.

        :param transcriptome_version: A transcriptome version (e.g., mm10_e87, ...).
        :type transcriptome_version: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A transcriptome resource directory path.
        :rtype: str | None
        """
        file_path = cls.get_expanded_directory(option='transcriptomes')

        if transcriptome_version:
            file_path = os.path.join(file_path, transcriptome_version)

        return cls._prepend_resource(file_path=file_path, absolute=absolute)

    @classmethod
    def get_resource_transcriptome_index(cls, transcriptome_version, transcriptome_index, absolute=True):
        """Get a transcriptome index resource directory path.

        :param transcriptome_version: A transcriptome version (e.g., mm10_e87, ...)
        :type transcriptome_version: str
        :param transcriptome_index: A transcriptome index (e.g., star, tophat, ...).
        :type transcriptome_index: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A transcriptome index resource directory path.
        :rtype: str | None
        :raise Exception: If the transcriptome index name is unknown.
        """
        index_directory = Index.get(option=transcriptome_index)

        if index_directory is None:
            raise Exception(f'Unknown transcriptome index name {transcriptome_index!r}.')
        else:
            return os.path.join(
                cls.get_resource_transcriptome(transcriptome_version=transcriptome_version, absolute=absolute),
                index_directory)

    @classmethod
    def get_resource_transcriptome_gtf(cls, transcriptome_version, transcriptome_index, basic=True, absolute=True):
        """Get a transcriptome GTF resource file path.

        :param transcriptome_version: A transcriptome version (e.g., mm10_e87, ...).
        :type transcriptome_version: str
        :param transcriptome_index: A transcriptome index (e.g., star, tophat, ...).
        :type transcriptome_index: str
        :param basic: Get a basic transcriptome.
        :type basic: bool
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A transcriptome GTF resource file path.
        :rtype: str | None
        :raise Exception: If the transcriptome index name is unknown.
        """
        if basic:
            file_name = transcriptome_version + '_basic.gtf'
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

        :param transcriptome_version: A transcriptome version (e.g., mm10_e87, ...).
        :type transcriptome_version: str
        :param transcriptome_index: A transcriptome index (e.g., star, tophat, ...).
        :type transcriptome_index: str
        :param basic: Get a basic transcriptome.
        :type basic: bool
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A transcriptome TxDb resource file path.
        :rtype: str | None
        :raise Exception: If the transcriptome index name is unknown.
        """
        if basic:
            file_name = transcriptome_version + '_basic.sqlite'
        else:
            file_name = transcriptome_version + '.sqlite'

        return os.path.join(
            cls.get_resource_transcriptome_index(
                transcriptome_version=transcriptome_version,
                transcriptome_index=transcriptome_index,
                absolute=absolute),
            file_name)

    @classmethod
    def get_resource_gatk_bundle(cls, genome_version, gatk_bundle_version, absolute=True):
        """Get a GATK bundle resource directory path.

        :param genome_version: A genome assembly version.
        :type genome_version: str
        :param gatk_bundle_version: A GATK bundle version.
        :type gatk_bundle_version: str
        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A GATK bundle resource directory path.
        :rtype: str | None
        """
        file_path = cls.get_expanded_directory(option='gatk_bundle')

        if genome_version:
            file_path = os.path.join(file_path, genome_version)

        if gatk_bundle_version:
            file_path = os.path.join(file_path, gatk_bundle_version)

        return cls._prepend_resource(file_path=file_path, absolute=absolute)

    @classmethod
    def get_resource_intervals(cls, absolute=True):
        """Get a Target Intervals resource directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A target intervals resource directory path.
        :rtype: str | None
        """
        return cls._prepend_resource(file_path=cls.get_expanded_directory(option='intervals'), absolute=absolute)

    @classmethod
    def get_resource_cadd(cls, absolute=True):
        """Get A Combined Annotation Dependent Depletion (CADD) resource directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A Combined Annotation Dependent Depletion (CADD) resource directory path.
        :rtype: str | None
        """
        return cls._prepend_resource(file_path=cls.get_expanded_directory(option='cadd'), absolute=absolute)

    @classmethod
    def get_resource_cosmic(cls, absolute=True):
        """Get a Catalogue Of Somatic Mutations In Cancer (COSMIC) resource directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A Catalogue Of Somatic Mutations In Cancer (COSMIC) resource directory path.
        :rtype: str | None
        """
        return cls._prepend_resource(file_path=cls.get_expanded_directory(option='cosmic'), absolute=absolute)

    @classmethod
    def get_resource_snpeff_data(cls, absolute=True):
        """Get a snpEff data resource directory path.

        :param absolute: Get an absolute file path.
        :type absolute: bool
        :return: A snpEff data resource directory path.
        :rtype: str | None
        """
        return cls._prepend_resource(file_path=cls.get_expanded_directory(option='snpeff_data'), absolute=absolute)


class Index(BaseSection):
    """The :py:class:`bsf.standards.Index` class models genome or transcriptome index directory defaults.

    The defaults are read from the :literal:`[indices]` section of the global configuration file.
    """
    _section = 'indices'


class Operator(BaseSection):
    """The :py:class:`bsf.standards.Operator` class models the operator's defaults.

    The defaults are read from the :literal:`[operator]` section of the global configuration file.
    """
    _section = 'operator'

    @classmethod
    def get_contact(cls):
        """Get an operator contact information.

        :return: An operator contact information.
        :rtype: str | None
        """
        return cls.get(option='contact')

    @classmethod
    def get_e_mail(cls):
        """Get an operator e-mail information.

        :return: An operator e-mail information.
        :rtype: str | None
        """
        return cls.get(option='e_mail')

    @classmethod
    def get_institution(cls):
        """Get an operator institution information.

        :return: An operator institution information.
        :rtype: str | None
        """
        return cls.get(option='institution')

    @classmethod
    def get_sequencing_centre(cls):
        """Get an operator sequencing centre information.

        :return: An operator sequencing centre information.
        :rtype: str | None
        """
        return cls.get(option='sequencing_centre')


class UCSC(BaseSection):
    """The :py:class:`bsf.standards.UCSC` class models UCSC Genome Browser Uniform Resource Locator (URL) defaults.

    The defaults are read from the :literal:`[ucsc]` section of the global configuration file.
    """

    _section = 'ucsc'

    @classmethod
    def get_protocol(cls):
        """Get a UCSC Genome Browser URL protocol (i.e., 'http' or 'https').

        :return: A UCSC Genome Browser URL protocol.
        :rtype: str | None
        """
        return cls.get(option='protocol')

    @classmethod
    def get_host_name(cls):
        """Get a UCSC Genome Browser URL host name (e.g., genome.ucsc.edu, genome-euro.ucsc.edu, ...).

        :return: A UCSC Genome Browser URL host name.
        :rtype: str | None
        """
        return cls.get(option='host_name')


class URL(BaseSection):
    """The :py:class:`bsf.standards.URL` class models web server Uniform Resource Locator (URL) defaults.

    The defaults are read from the :literal:`[url]` section of the global configuration file.
    """

    _section = 'url'

    @classmethod
    def get_protocol(cls):
        """Get a web server URL protocol (i.e., 'http' or 'https').

        :return: A web server URL protocol.
        :rtype: str | None
        """
        return cls.get(option='protocol')

    @classmethod
    def get_host_name(cls):
        """Get a web server host name.

        :return: A web server URL host name.
        :rtype: str | None
        """
        return cls.get(option='host_name')

    @classmethod
    def get_relative_dna(cls):
        """Get a relative URL to the DNA directory.

        :return: A relative 'DNA' URL path.
        :rtype: str | None
        """
        return cls.get(option='relative_dna')

    @classmethod
    def get_relative_projects(cls):
        """Get a relative URL to the analysis projects directory.

        :return: A relative 'projects' URL path.
        :rtype: str | None
        """
        return cls.get(option='relative_projects')

    @classmethod
    def get_absolute_base(cls):
        """Get an absolute URL to the base directory.

        :return: An absolute URL to the base directory.
        :rtype: str
        """
        if cls.get_protocol():
            return cls.get_protocol() + '://' + cls.get_host_name()
        else:
            return '//' + cls.get_host_name()

    @classmethod
    def get_absolute_dna(cls):
        """Get an absolute URL to the DNA directory.

        :return: An absolute URL to the DNA directory.
        :rtype: str
        """
        return '/'.join((cls.get_absolute_base(), cls.get_relative_dna()))

    @classmethod
    def get_absolute_projects(cls):
        """Get an absolute URL to the analysis projects directory.

        :return: An absolute URL to the analysis projects directory.
        :rtype: str
        """
        return '/'.join((cls.get_absolute_base(), cls.get_relative_projects()))


class VendorQualityFilter(BaseSection):
    """The :py:class:`bsf.standards.VendorQualityFilter` class models (SAM) Vendor Quality Filter defaults.

    The defaults are read from the :literal:`[VendorQualityFilter]` section of the global configuration file.
    For each flow cell type, a boolean value specifies whether vendor quality filtering should be applied or not.
    """

    _section = 'vendor_quality_filter'

    @classmethod
    def get_vendor_quality_filter(cls, flow_cell_type):
        """Get a vendor quality filter setting for a particular flow cell type.

        The flow cell type is accessible via property :py:meth:`bsf.illumina.RunParameters.get_flow_cell_type`,
        practically directly via :py:meth:`bsf.illumina.RunFolder.run_parameters.get_flow_cell_type`.

        :param flow_cell_type: A fLow cell (chemistry) type.
        :type flow_cell_type: str
        :return: A vendor quality filter setting.
        :rtype: bool | None
        """
        if not flow_cell_type:
            return

        # To avoid emitting cryptic error messages, check for the presence of the flow cell type before
        # and inform the user about an eventual problem with a hopefully more meaningful message.

        if not Configuration.get_global_configuration().config_parser.has_option(
                section=cls._section,
                option=flow_cell_type):
            raise Exception(f'The flow cell type {flow_cell_type!r} is not defined in the {cls._section!r} '
                            f'section of the standard configuration file {Configuration.get_global_file_path()!r}.')

        return cls.getboolean(option=flow_cell_type)


class Secrets(BaseSection):
    """The :py:class:`bsf.standards.Secrets` class models file paths to configuration files with secrets.
    """
    _section = 'secrets'

    _user_mask = stat.S_IRWXG | stat.S_IRWXO

    @classmethod
    def get_file_path(cls, option):
        """Get a configuration file path with secrets.

        The configuration path is defined in the specified option under the
        :literal:`[secrets]` configuration section.
        Also checks that the file is only readable by the user and not accessible for group and other.

        :param option: A configuration option.
        :type option: str
        :return: A file path to configuration file with secrets.
        :rtype: str | None
        """
        file_path = cls.get_expanded_directory(option=option)

        if not file_path:
            return None

        path_stat_result = os.stat(path=file_path, follow_symlinks=True)

        if path_stat_result.st_mode & cls._user_mask:
            raise Exception(
                f'The secrets configuration file {file_path!r} has file mode {path_stat_result.st_mode:#0o}, '
                f'but should obey user mask {cls._user_mask:#0o}.')

        return file_path

    @classmethod
    def get_azure_file_path(cls):
        """Get a configuration file path with :literal:`Microsoft Azure` secrets.

        Also checks that the file is only readable by the user and not accessible for group and other.

        :return: A file path to configuration file with secrets.
        :rtype: str | None
        """
        return cls.get_file_path(option='azure_file_path')

    @classmethod
    def get_mysql_file_path(cls):
        """Get a configuration file path with :literal:`Oracle MySQL` secrets.

        Also checks that the file is only readable by the user and not accessible for group and other.

        :return: A file path to configuration file with secrets.
        :rtype: str | None
        """
        return cls.get_file_path(option='mysql_file_path')


class Central(BaseSection):
    """The :py:class:`bsf.standards.Central` class models the central XML configuration document.
    """

    _section = 'central'

    _global_element_tree = None
    _global_environment = 'BSF_PYTHON_XML'
    _global_file_path = '~/.bsfpython.xml'

    @classmethod
    def get_global_file_path(cls):
        """Get the global XML configuration (:literal:`*.xml`) file path.

        The global configuration file path is based on the value of
        environment variable :literal:`BSF_PYTHON_XML` if defined,
        the :literal:`configuration_xml` option in the :literal:`[central]` section of the
        global UNIX-style initialisation (:literal:`*.ini`) file path if defined,
        or defaults to :literal:`~/.bsfpython.xml` otherwise.

        :return: Global XML configuration file path
        :rtype: str
        """
        if cls._global_environment in os.environ and os.environ[cls._global_environment]:
            file_path = os.environ[cls._global_environment]
        elif cls.get(option='configuration_xml'):
            file_path = cls.get(option='configuration_xml')
        else:
            file_path = cls._global_file_path

        file_path = os.path.expanduser(file_path)
        file_path = os.path.expandvars(file_path)
        file_path = os.path.normpath(file_path)

        return file_path

    @classmethod
    def get_element_tree(cls):
        """Get a central :py:class:`xml.etree.ElementTree.ElementTree` object.

        :return: A central :py:class:`xml.etree.ElementTree.ElementTree` object.
        :rtype: ElementTree
        """
        if cls._global_element_tree is None:
            cls._global_element_tree = ElementTree(file=cls.get_global_file_path())

        return cls._global_element_tree


class CentralIndexDirectories(object):
    """The :py:class:`bsf.standards.CentralIndexDirectories` class models index directory names of the
    central XML configuration document.

    :ivar bowtie1: A Bowtie1 index directory name.
    :type bowtie1: str | None
    :ivar bowtie2: A Bowtie2 index directory name.
    :type bowtie2: str | None
    :ivar bwa: A BWA index directory name.
    :type bwa: str | None
    :ivar hisat2: A Hisat2 index directory name.
    :type hisat2: str | None
    :ivar kallisto: A Kallisto index directory name.
    :type kallisto: str | None
    :ivar star: A STAR index directory name.
    :type star: str | None
    :ivar tophat2: A Tophat2 index directory name.
    :type tophat2: str | None
    """

    def __init__(self):
        """Initialise a :py:class:`bsf.standards.CentralIndexDirectories` object.
        """
        element_tree = Central.get_element_tree()

        self.bowtie1 = element_tree.find(path='IndexDirectories/IndexDirectoryBowtie1').text
        self.bowtie2 = element_tree.find(path='IndexDirectories/IndexDirectoryBowtie2').text
        self.bwa = element_tree.find(path='IndexDirectories/IndexDirectoryBWA').text
        # self.hisat1 = element_tree.find(path='IndexDirectories/IndexDirectoryHisat1').text
        self.hisat2 = element_tree.find(path='IndexDirectories/IndexDirectoryHisat2').text
        self.kallisto = element_tree.find(path='IndexDirectories/IndexDirectoryKallisto').text
        self.star = element_tree.find(path='IndexDirectories/IndexDirectorySTAR').text
        # self.tophat1 = element_tree.find(path='IndexDirectories/IndexDirectoryTophat1').text
        self.tophat2 = element_tree.find(path='IndexDirectories/IndexDirectoryTophat2').text


class CentralVendorQualityFilters(object):
    """The :py:class:`bsf.standards.CentralVendorQualityFilters` class models vendor quality filter settings of the
    central XML configuration document.
    """

    # Python dict of (boolean state) Python str objects and Python bool value objects.
    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False
    }

    @classmethod
    def str_to_bool(cls, value):
        """Convert Python :py:class:`str` objects to Python :py:class:`bool` objects.

        :param value: A Python :py:class:`str` object.
        :type value: str
        :return: A Python :py:class:`bool` object.
        :rtype: bool | None
        """
        value = value.lower()
        if value in cls._boolean_states:
            return cls._boolean_states[value]
        else:
            return None

    def __init__(self):
        """Initialise a :py:class:`bsf.standards.CentralVendorQualityFilters` object.
        """
        self.filter_dict: Dict[str, bool] = dict()

        for element in Central.get_element_tree().find(path='VendorQualityFilters'):
            self.filter_dict[element.text] = self.str_to_bool(value=element.get(key='filter'))

        return
