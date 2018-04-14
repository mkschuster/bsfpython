"""bsf.analyses.bowtie

A package of classes and methods supporting Bowtie alignment analyses.
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

from __future__ import print_function

from bsf import Analysis, Runnable
from bsf.process import RunnableStep


class Bowtie1(Analysis):
    """The C{bsf.analyses.bowtie.Bowtie1} class represents the logic to run the Bowtie1 short read aligner.

    Attributes:
    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_bowtie1: C{bsf.Stage.name} for the Bowtie1 alignment C{bsf.Analysis} stage
    @type stage_name_bowtie1: str
    @ivar replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects for processing or run them separately
    @type replicate_grouping: bool
    @ivar bwa_genome_db: Genome sequence file path with BWA index
    @type bwa_genome_db: str | unicode
    """

    name = 'Bowtie1 Analysis'
    prefix = 'bowtie1'

    stage_name_bowtie1 = 'bowtie1'

    def __init__(
            self,
            configuration=None,
            project_name=None,
            genome_version=None,
            input_directory=None,
            output_directory=None,
            project_directory=None,
            genome_directory=None,
            e_mail=None,
            debug=0,
            stage_list=None,
            collection=None,
            sample_list=None,
            replicate_grouping=False,
            bwa_genome_db=None):
        """Initialise a C{bsf.analyses.bowtie.Bowtie1}.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param project_name: Project name
        @type project_name: str
        @param genome_version: Genome version
        @type genome_version: str
        @param input_directory: C{bsf.Analysis}-wide input directory
        @type input_directory: str
        @param output_directory: C{bsf.Analysis}-wide output directory
        @type output_directory: str
        @param project_directory: C{bsf.Analysis}-wide project directory,
            normally under the C{bsf.Analysis}-wide output directory
        @type project_directory: str
        @param genome_directory: C{bsf.Analysis}-wide genome directory,
            normally under the C{bsf.Analysis}-wide project directory
        @type genome_directory: str
        @param e_mail: e-Mail address for a UCSC Genome Browser Track Hub
        @type e_mail: str
        @param debug: Integer debugging level
        @type debug: int
        @param stage_list: Python C{list} of C{bsf.Stage} objects
        @type stage_list: list[bsf.Stage]
        @param collection: C{bsf.ngs.Collection}
        @type collection: bsf.ngs.Collection
        @param sample_list: Python C{list} of C{bsf.ngs.Sample} objects
        @type sample_list: list[bsf.ngs.Sample]
        @param replicate_grouping: Group individual C{bsf.ngs.PairedReads} objects for processing or
            run them separately
        @type replicate_grouping: bool
        @param bwa_genome_db: Genome sequence file path with BWA index
        @type bwa_genome_db: str | unicode
        @return:
        @rtype:
        """

        super(Bowtie1, self).__init__(
            configuration=configuration,
            project_name=project_name,
            genome_version=genome_version,
            input_directory=input_directory,
            output_directory=output_directory,
            project_directory=project_directory,
            genome_directory=genome_directory,
            e_mail=e_mail,
            debug=debug,
            stage_list=stage_list,
            collection=collection,
            sample_list=sample_list)

        # Sub-class specific ...

        if replicate_grouping is None:
            self.replicate_grouping = False
        else:
            assert isinstance(replicate_grouping, bool)
            self.replicate_grouping = replicate_grouping

        if bwa_genome_db is None:
            self.bwa_genome_db = str()
        else:
            self.bwa_genome_db = bwa_genome_db

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.bowtie.Bowtie1} via a C{bsf.standards.Configuration} section.

        Instance variables without a configuration option remain unchanged.

        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """

        super(Bowtie1, self).set_configuration(configuration=configuration, section=section)

        config_parser = configuration.config_parser

        option = 'replicate_grouping'
        if config_parser.has_option(section=section, option=option):
            self.replicate_grouping = config_parser.getboolean(section=section, option=option)

        # Get the genome database.

        option = 'bwa_genome_db'
        if config_parser.has_option(section=section, option=option):
            self.bwa_genome_db = config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run a C{bsf.analyses.bowtie.Bowtie1} analysis.

        @return:
        @rtype:
        """

        # Get global defaults.

        # default = Default.get_global_default()

        super(Bowtie1, self).run()

        stage_bowtie1 = self.get_stage(name=self.stage_name_bowtie1)

        for sample in self.sample_list:
            if self.debug > 0:
                print(self, 'Sample name:', repr(sample.name))
                print(sample.trace(1))

            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=self.replicate_grouping, exclude=True)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                if not len(paired_reads_dict[paired_reads_name]):
                    # Skip replicate keys, which PairedReads objects have all been excluded.
                    continue

                prefix_bowtie1 = '_'.join((stage_bowtie1.name, paired_reads_name))

                # Create a Runnable and Executable for processing each Bowtie1 alignment.

                runnable_bowtie1 = self.add_runnable(
                    runnable=Runnable(
                        name=prefix_bowtie1,
                        code_module='bsf.runnables.generic',
                        working_directory=self.genome_directory,
                        cache_directory=self.cache_directory,
                        debug=self.debug))
                self.set_stage_runnable(
                    stage=stage_bowtie1,
                    runnable=runnable_bowtie1)

                runnable_step = runnable_bowtie1.add_runnable_step(
                    runnable_step=RunnableStep(
                        name='bowtie1',
                        obsolete_file_path_list=[
                        ]))
                """ @type runnable_step: RunnableStep """
                runnable_step.arguments.append(self.bwa_genome_db)

                reads_list_1 = list()
                reads_list_2 = list()

                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if paired_reads.reads_1:
                        reads_list_1.append(paired_reads.reads_1.file_path)
                    if paired_reads.reads_2:
                        reads_list_2.append(paired_reads.reads_2.file_path)

                if len(reads_list_1) and not len(reads_list_2):
                    # For Bowtie1 unpaired reads are an argument, paired come with options -1 <m1> and -2 <m2>.
                    runnable_step.arguments.append(','.join(reads_list_1))
                elif len(reads_list_1) and len(reads_list_2):
                    runnable_step.add_option_short(key='1', value=','.join(reads_list_1))
                if len(reads_list_2):
                    runnable_step.add_option_short(key='2', value=','.join(reads_list_2))

        return
