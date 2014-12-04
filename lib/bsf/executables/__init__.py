"""bsf.executables

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


from bsf import Command, Executable


class Bowtie1(Executable):
    """Bowtie1 short read aligner class.

    Reference: http://bowtie-bio.sourceforge.net/manual.shtml
    """

    def __init__(self, name, analysis):
        """Initialise a C{Bowtie1} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(Bowtie1, self).__init__(name=name, program='bowtie')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Bowtie1 options.


class Bowtie2(Executable):
    """Bowtie2 short read aligner class."""

    def __init__(self, name, analysis):
        """Initialise a C{Bowtie2} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(Bowtie2, self).__init__(name=name, program='bowtie2')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Bowtie2 options.


class BWA(Executable):
    """Burrows-Wheeler Aligner version class.

    Reference: http://bio-bwa.sourceforge.net/
    Usage: bwa mem db_prefix reads.fq [mates.fq]
    """

    def __init__(self, name, analysis):
        """Initialise a C{BWA} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        super(BWA, self).__init__(name=name, program='bwa', sub_command=Command(command='mem'))

        # The options have to be set for the 'mem' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_configuration(configuration=analysis.configuration, section=section)

        # Set default BWA mem options.

        # None for the moment.


class Cuffdiff(Executable):
    """Cuffdiff differential expression class.

    Reference: http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff
    Usage: cuffdiff [options]* <transcripts.gtf>
    <sample1_replicate1.sam[,...,sample1_replicateM]>
    <sample2_replicate1.sam[,...,sample2_replicateM.sam]>...
    [sampleN.sam_replicate1.sam[,...,sample2_replicateM.sam]]
    """

    def __init__(self, name, analysis):
        """Initialise a C{Cuffdiff} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(Cuffdiff, self).__init__(name=name, program='cuffdiff')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Cuffdiff options.

        if not ('library-type' in self.options and self.options['library-type']):
            self.add_option_long(key='library-type', value='fr-unstranded')

        # if not ('num-threads' in self.options and self.options['num-threads']):
        #     self.add_option_long(key='num-threads', value='1')

        if not 'quiet' in self.options:
            self.add_switch_long(key='quiet')

        if not 'no-update-check' in self.options:
            self.add_switch_long(key='no-update-check')


class Cufflinks(Executable):
    """Cufflinks transcript assembler class.

    Reference: http://cufflinks.cbcb.umd.edu/manual.html
    Usage: cufflinks [options]* <aligned_reads.(sam/bam)>
    """

    def __init__(self, name, analysis):
        """Initialise a C{Cufflinks} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(Cufflinks, self).__init__(name=name, program='cufflinks')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Cufflinks options.

        if not 'library-type' in self.options:
            self.add_option_long(key='library-type', value='fr-unstranded')

        if not ('num-threads' in self.options and self.options['num-threads']):
            self.add_option_long(key='num-threads', value='1')

        if not 'quiet' in self.options:
            self.add_switch_long(key='quiet')

        if not 'no-update-check' in self.options:
            self.add_switch_long(key='no-update-check')


class Cuffmerge(Executable):
    """Cuffmerge transcript assembly merge class.

    Reference: http://cufflinks.cbcb.umd.edu/manual.html#cuffmerge
    Usage: cuffmerge [options]* <assembly_GTF_list.txt>
    """

    def __init__(self, name, analysis):
        """Initialise a C{Cuffmerge} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # super(Cuffmerge, self).__init__(name=name, program='cuffmerge')
        # TODO: Experimentally change this so that the new bsf_run_rnaseq_cuffmerge.py script gets used.
        # Although this seems rather successful, so far, the runnable option is not ideal.
        # Maybe the DRMS object should supply a mechanism to execute via the BSF Runner script.

        # assert isinstance(analysis, Analysis)

        super(Cuffmerge, self).__init__(name=name, program='bsf_run_rnaseq_cuffmerge.py')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # TODO: Should this be refactored so that the bsf_run_rnaseq_cuffmerge.py script can be used from a DRMS or
        # Executable object?
        self.add_option_long(key='runnable', value='Cuffmerge')

        # Set default Cuffmerge options.

        # if not ('num-threads' in self.options and self.options['num-threads']):
        #     self.add_option_long(key='num-threads', value='1')


class TopHat(Executable):
    """TopHat RNA-Seq aligner class.

    Reference: http://tophat.cbcb.umd.edu/manual.html
    Usage: tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]
    Arguments:
    <ebwt_base> Base name of the index to be searched.
    <reads1_1[,...,readsN_1]>
    <[reads1_2,...readsN_2]>
    """

    def __init__(self, name, analysis):
        """Initialise a C{TopHat} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(TopHat, self).__init__(name=name, program='tophat2')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default TopHat options.

        if not ('library-type' in self.options and self.options['library-type']):
            self.add_option_long(key='library-type', value='fr-unstranded')

        if not ('num-threads' in self.options and self.options['num-threads']):
            self.add_option_long(key='num-threads', value='1')

        if not ('coverage-search' in self.options and self.options['coverage-search']):
            self.add_switch_long(key='no-coverage-search')


class Macs14(Executable):
    """Model-based Analysis for ChIP-Seq (MACS) version 1.4 peak-caller class.

    Reference: http://liulab.dfci.harvard.edu/MACS/
    """

    def __init__(self, name, analysis):
        """Initialise a C{Macs14} (MACS 1.4) object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(Macs14, self).__init__(name=name, program='macs14')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        if not ('gsize' in self.options and self.options['gsize']):
            raise Exception(
                "A 'gsize' option is required in the {!r} configuration section.".
                format(section))


class Macs2Bdgcmp(Executable):
    """Model-based Analysis for ChIP-Seq (MACS) version 2 bedGraph comparison class.

    Reference: http://liulab.dfci.harvard.edu/MACS/
    """

    def __init__(self, name, analysis):
        """Initialise a C{Macs2Bdgcmp} (MACS2 BedGraph Comparison) object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        super(Macs2Bdgcmp, self).__init__(name=name, program='macs2', sub_command=Command(command='bdgcmp'))

        # The options have to be set for the 'bdgcmp' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        # None for the moment.


class Macs2Callpeak(Executable):
    """Model-based Analysis for ChIP-Seq (MACS) version 2 peak-caller class.

    Reference: http://liulab.dfci.harvard.edu/MACS/
    """

    def __init__(self, name, analysis):
        """Initialise a C{Macs2Callpeak} (MACS 2.0 peak caller) object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(Macs2Callpeak, self).__init__(name=name, program='macs2', sub_command=Command(command='callpeak'))

        # The options have to be set for the 'callpeak' sub-command.
        section = analysis.configuration.section_from_instance(self)
        self.sub_command.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        if not ('gsize' in self.sub_command.options and self.sub_command.options['gsize']):
            raise Exception(
                "A 'gsize' option is required in the {!r} configuration section.".
                format(section))


class FastQC(Executable):
    """FastQC quality checker class.

    Reference: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """

    def __init__(self, name, analysis):
        """Initialise a C{FastQC} object.

        @param name: Name
        @type name: str
        @param analysis: C{Analysis}
        @type analysis: Analysis
        """

        # assert isinstance(analysis, Analysis)

        super(FastQC, self).__init__(name=name, program='fastqc')

        section = analysis.configuration.section_from_instance(self)
        self.set_configuration(configuration=analysis.configuration, section=section)

        # Set default Macs options.

        self.add_switch_long(key='quiet')
