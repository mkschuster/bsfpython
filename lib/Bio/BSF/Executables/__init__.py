"""Bio.BSF.Executable

A package of classes and methods supporting executable programs and scripts.
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


from Bio.BSF import Command, Executable


class Bowtie1(Executable):
    """BSF Bowtie1 short read aligner class.

    Reference: http://bowtie-bio.sourceforge.net/manual.shtml
    """

    def __init__(self, name, analysis):
        """Initialise a Bowtie1 object.

        :param self: BSF Bowtie1 object
        :type self: Bowtie1
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(Bowtie1, self).__init__(name=name, program='bowtie')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Bowtie1 options.


class Bowtie2(Executable):
    """BSF Bowtie2 short read aligner class."""

    def __init__(self, name, analysis):
        """Initialise a Bowtie2 object.

        :param self: BSF Bowtie2 object
        :type self: Bowtie2
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(Bowtie2, self).__init__(name=name, program='bowtie2')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Bowtie2 options.


class BWA07(Executable):
    """BSF Burrows Wheeler Aligner version 0.7 class.

    Reference: http://bio-bwa.sourceforge.net/
    Usage: bwa mem db_prefix reads.fq [mates.fq]
    """

    def __init__(self, name, analysis):
        """Initialise a BWA07 object.

        :param self: BSF BWA07 object
        :type self: BWA07
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        super(BWA07, self).__init__(name=name, program='bwa', sub_command=Command(command='mem'))

        # The options have to be set for the 'mem' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default BWA mem options.

        # None for the moment.


class Cuffdiff(Executable):
    """BSF Cuffdiff differential expression class.

    Reference: http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff
    Usage: cuffdiff [options]* <transcripts.gtf>
    <sample1_replicate1.sam[,...,sample1_replicateM]>
    <sample2_replicate1.sam[,...,sample2_replicateM.sam]>...
    [sampleN.sam_replicate1.sam[,...,sample2_replicateM.sam]]
    """

    def __init__(self, name, analysis):
        """Initialise a Cuffdiff object.

        :param self: BSF Cuffdiff object
        :type self: Cuffdiff
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(Cuffdiff, self).__init__(name=name, program='cuffdiff')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Cuffdiff options.

        if not ('library-type' in self.options and self.options['library-type']):
            self.add_OptionLong(key='library-type', value='fr-unstranded')

        # if not ('num-threads' in self.options and self.options['num-threads']):
        #     self.add_OptionLong(key='num-threads', value='1')

        if not 'quiet' in self.options:
            self.add_SwitchLong(key='quiet')

        if not 'no-update-check' in self.options:
            self.add_SwitchLong(key='no-update-check')


class Cufflinks(Executable):
    """BSF Cufflinks transcript assembler class.

    Reference: http://cufflinks.cbcb.umd.edu/manual.html
    Usage: cufflinks [options]* <aligned_reads.(sam/bam)>
    """

    def __init__(self, name, analysis):
        """Initialise a Cufflinks object.

        :param self: BSF Cufflinks object
        :type self: Cufflinks
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(Cufflinks, self).__init__(name=name, program='cufflinks')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Cufflinks options.

        if not 'library-type' in self.options:
            self.add_OptionLong(key='library-type', value='fr-unstranded')

        if not ('num-threads' in self.options and self.options['num-threads']):
            self.add_OptionLong(key='num-threads', value='1')

        if not 'quiet' in self.options:
            self.add_SwitchLong(key='quiet')

        if not 'no-update-check' in self.options:
            self.add_SwitchLong(key='no-update-check')


class Cuffmerge(Executable):
    """BSF Cuffmerge transcript assembly merge class.

    Reference: http://cufflinks.cbcb.umd.edu/manual.html#cuffmerge
    Usage: cuffmerge [options]* <assembly_GTF_list.txt>
    """

    def __init__(self, name, analysis):
        """Initialise a Cuffmerge object.

        :param self: BSF Cuffmerge object
        :type self: Cuffmerge
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # super(Cuffmerge, self).__init__(name=name, program='cuffmerge')
        # TODO: Experimentally change this so that the new bsf_runner.py script gets used.
        # Although this seems rather successful, so far, the runnable option is not ideal.
        # Maybe the DRMS object should supply a mechanism to execute via the BSF Runner script.

        # assert isinstance(analysis, Analysis)

        super(Cuffmerge, self).__init__(name=name, program='bsf_runner.py')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # TODO: Should this be refactored so that the bsf_runner.py script can be used from a DRMS or Executable object?
        self.add_OptionLong(key='runnable', value='Cuffmerge')

        # Set default Cuffmerge options.

        # if not ('num-threads' in self.options and self.options['num-threads']):
        #     self.add_OptionLong(key='num-threads', value='1')


class TopHat(Executable):
    """BSF TopHat RNA-Seq aligner class.

    Reference: http://tophat.cbcb.umd.edu/manual.html
    Usage: tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]
    Arguments:
    <ebwt_base> Base name of the index to be searched.
    <reads1_1[,...,readsN_1]>
    <[reads1_2,...readsN_2]>
    """

    def __init__(self, name, analysis):
        """Initialise a TopHat object.

        :param self: BSF TopHat object
        :type self: TopHat
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(TopHat, self).__init__(name=name, program='tophat2')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default TopHat options.

        if not ('library-type' in self.options and self.options['library-type']):
            self.add_OptionLong(key='library-type', value='fr-unstranded')

        if not ('num-threads' in self.options and self.options['num-threads']):
            self.add_OptionLong(key='num-threads', value='1')

        if not ('coverage-search' in self.options and self.options['coverage-search']):
            self.add_SwitchLong(key='no-coverage-search')


class Macs14(Executable):
    """BSF Model-based Analysis for ChIP-Seq (MACS) version 1.4 peak-caller class.

    Reference: http://liulab.dfci.harvard.edu/MACS/
    """

    def __init__(self, name, analysis):
        """Initialise a Macs14 object.

        :param self: BSF Macs14 object
        :type self: Macs14
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(Macs14, self).__init__(name=name, program='macs14')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        if not ('gsize' in self.options and self.options['gsize']):
            message = "A 'gsize' option is required in the [{}] configuration section".format(section)
            raise Exception(message)


class Macs2Bdgcmp(Executable):
    """BSF Model-based Analysis for ChIP-Seq (MACS) version 2 bedGraph comparison class.

    Reference: http://liulab.dfci.harvard.edu/MACS/
    """

    def __init__(self, name, analysis):
        """Initialise a MACS2 BedGraph Comparison object.

        :param self: BSF Macs2Bdgmp object
        :type self: Macs2Bdgcmp
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        super(Macs2Bdgcmp, self).__init__(name=name, program='macs2', sub_command=Command(command='bdgcmp'))

        # The options have to be set for the 'bdgcmp' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        # None for the moment.


class Macs2Callpeak(Executable):
    """BSF Model-based Analysis for ChIP-Seq (MACS) version 2 peak-caller class.

    Reference: http://liulab.dfci.harvard.edu/MACS/
    """

    def __init__(self, name, analysis):
        """Initialise a Macs2Callpeak object.

        :param self: BSF Macs2Callpeak object
        :type self: Macs2Callpeak
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(Macs2Callpeak, self).__init__(name=name, program='macs2', sub_command=Command(command='callpeak'))

        # The options have to be set for the 'callpeak' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        if not ('gsize' in self.sub_command.options and self.sub_command.options['gsize']):
            message = "A 'gsize' option is required in the [{}] configuration option".format(section)
            raise Exception(message)


class FastQC(Executable):
    """BSF FastQC quality checker class.

    Reference: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """

    def __init__(self, name, analysis):
        """Initialise a FastQC object.

        :param self: BSF FastQC object
        :type self: FastQC
        :param name: Name
        :type name: str
        :param analysis: Bio.BSF.Analysis object or a sub-class thereof
        :type analysis: Analysis
        :return: Nothing
        :rtype: None
        """

        # assert isinstance(analysis, Analysis)

        super(FastQC, self).__init__(name=name, program='fastqc')

        section = analysis.configuration.section_from_instance(self)
        self.set_Configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        self.add_SwitchLong(key='quiet')
