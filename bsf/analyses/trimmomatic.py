"""bsf.analyses.trimmomatic

A package of classes and methods supporting the Trimmomatic tool.
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


import os

from bsf import Analysis, FilePath, Runnable
from bsf.ngs import Reads, PairedReads
from bsf.process import Command, RunnableStep, RunnableStepJava, RunnableStepMakeDirectory


class FilePathTrimmomatic(FilePath):
    """The C{bsf.analyses.trimmomatic.FilePathTrimmomatic} models files in a sample-specific Trimmomatic directory.

    Attributes:
    @ivar output_directory: Output directory
    @type output_directory: str | unicode
    @ivar trim_log_tsv: Trimmomatic trim log Tab-Separated Value (TSV) file path
    @type trim_log_tsv: str | unicode
    @ivar summary_tsv: Summary Tab-Separated Value (TSV) file path
    @type summary_tsv: str | unicode
    @ivar coverage_png: Coverage Portable Network Graphics (PNG) file path
    @type coverage_png: str | unicode
    @ivar frequency_png: Frequency Portable Network Graphics (PNG) file path
    @type frequency_png: str | unicode
    @ivar surviving_png: Surviving Portable Network Graphics (PNG) file path
    @type surviving_png: str | unicode
    @ivar reads_1p: First Reads paired
    @type reads_1p: str | unicode
    @ivar reads_1u: First Reads unpaired
    @type reads_1u: str | unicode
    @ivar reads_2p: Second Reads paired
    @type reads_2p: str | unicode
    @ivar reads_2u: Second Reads unpaired
    @type reads_2u: str | unicode
    """

    def __init__(self, prefix):
        """Initialise a C{bsf.analyses.trimmomatic.FilePathTrimmomatic} object

        @param prefix: Prefix
        @type prefix: str | unicode
        @return:
        @rtype
        """
        super(FilePathTrimmomatic, self).__init__(prefix=prefix)

        self.output_directory = prefix
        # Automatic GNU Zip-compression of trim log files does not work.
        self.trim_log_tsv = prefix + '_trim_log.tsv'
        self.summary_tsv = prefix + '_summary.tsv'  # Defined by the R script.
        self.coverage_png = prefix + '_coverage.png'  # Defined by the R script.
        self.frequency_png = prefix + '_frequency.png'  # Defined by the R script.
        self.surviving_png = prefix + '_surviving.png'  # Defined by the R script.
        self.reads_1p = ''
        self.reads_1u = ''
        self.reads_2p = ''
        self.reads_2u = ''

        return


class Trimmomatic(Analysis):
    """The C{bsf.analyses.trimmomatic.Trimmomatic} class represents the logic to run the Trimmomatic analysis.

    Attributes:

    @cvar name: C{bsf.Analysis.name} that should be overridden by sub-classes
    @type name: str
    @cvar prefix: C{bsf.Analysis.prefix} that should be overridden by sub-classes
    @type prefix: str
    @cvar stage_name_trimmomatic: C{bsf.Stage.name} for the Trimmomatic stage
    @type stage_name_trimmomatic: str
    @ivar adapter_path: Adapter file path
    @type adapter_path: str | unicode
    @ivar trimming_step_pe_list: Colon-separated Trimmomatic steps for paired-end data
    @type trimming_step_pe_list: list[str | unicode]
    @ivar trimming_step_se_list: Colon-separated Trimmomatic steps for single-end data
    @type trimming_step_se_list: list[str | unicode]
    @ivar classpath_trimmomatic: Trimmomatic tool Java Archive (JAR) class path directory
    @type classpath_trimmomatic: str | unicode
    """

    name = 'Trimmomatic Analysis'
    prefix = 'trimmomatic'

    # stage_name_trimmomatic = '_'.join((prefix, 'trimmomatic'))
    stage_name_trimmomatic = prefix

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
            adapter_path=None,
            trimming_step_pe_list=None,
            trimming_step_se_list=None,
            classpath_trimmomatic=None):
        """Initialise a C{bsf.analyses.trimmomatic.Trimmomatic} object.

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
        @param adapter_path: Adapter file path
        @type adapter_path: str | unicode
        @param trimming_step_pe_list: Colon-separated Trimmomatic steps for paired-end data
        @type trimming_step_pe_list: list[str | unicode]
        @param trimming_step_se_list: Colon-separated Trimmomatic steps for single-end data
        @type trimming_step_se_list: list[str | unicode]
        @param classpath_trimmomatic: Trimmomatic tool Java Archive (JAR) class path directory
        @type classpath_trimmomatic: str | unicode
        @return:
        @rtype:
        """

        super(Trimmomatic, self).__init__(
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

        if adapter_path is None:
            self.adapter_path = str()
        else:
            self.adapter_path = adapter_path

        if trimming_step_pe_list is None:
            self.trimming_step_pe_list = list()
        else:
            self.trimming_step_pe_list = trimming_step_pe_list

        if trimming_step_se_list is None:
            self.trimming_step_se_list = list()
        else:
            self.trimming_step_se_list = trimming_step_se_list

        if classpath_trimmomatic is None:
            self.classpath_trimmomatic = str()
        else:
            self.classpath_trimmomatic = classpath_trimmomatic

        return

    def set_configuration(self, configuration, section):
        """Set instance variables of a C{bsf.analyses.trimmomatic.Trimmomatic} object via a section of a
        C{bsf.standards.Configuration} object.

        Instance variables without a configuration option remain unchanged.
        @param configuration: C{bsf.standards.Configuration}
        @type configuration: bsf.standards.Configuration
        @param section: Configuration file section
        @type section: str
        @return:
        @rtype:
        """
        super(Trimmomatic, self).set_configuration(configuration=configuration, section=section)

        # Sub-class specific ...

        option = 'adapter_path'
        if configuration.config_parser.has_option(section=section, option=option):
            self.adapter_path = configuration.config_parser.get(section=section, option=option)

        # Get the list of default trimming steps for paired-end data mode.

        option = 'trimming_steps_pe'
        if configuration.config_parser.has_option(section=section, option=option):
            self.trimming_step_pe_list = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    configuration.config_parser.get(section=section, option=option).split(',')))

        # Get the list of default trimming steps for single-end data mode.

        option = 'trimming_steps_se'
        if configuration.config_parser.has_option(section=section, option=option):
            self.trimming_step_se_list = filter(
                lambda x: x != '',
                map(
                    lambda x: x.strip(),
                    configuration.config_parser.get(section=section, option=option).split(',')))

        # Get the Trimmomatic tool Java Archive (JAR) class path directory.

        option = 'classpath_trimmomatic'
        if configuration.config_parser.has_option(section=section, option=option):
            self.classpath_trimmomatic = configuration.config_parser.get(section=section, option=option)

        return

    def run(self):
        """Run the C{bsf.analyses.trimmomatic.Trimmomatic} C{bsf.Analysis}.

        This method changes the C{bsf.ngs.Collection} object of this C{bsf.Analysis} to update with FASTQ file paths.
        @return:
        @rtype:
        """

        def run_adjust_illumina_clip_path(trimming_step_list):
            """Private function to adjust the adapter FASTA file path of ILLUMINACLIP trimming steps.

            If the file path is not absolute, prepend the value of the adapter_path
            instance variable.
            @param trimming_step_list: Python C{list} of trimming steps.
            @type trimming_step_list: list[str | unicode]
            @return:
            @rtype:
            """
            for i in range(0, len(trimming_step_list)):
                if trimming_step_list[i].startswith('ILLUMINACLIP'):
                    component_list = trimming_step_list[i].split(':')
                    if not os.path.isabs(component_list[1]):
                        component_list[1] = os.path.join(self.adapter_path, component_list[1])
                    trimming_step_list[i] = ':'.join(component_list)
            return

        def run_read_comparisons():
            """Private function to read a C{bsf.annotation.AnnotationSheet} CSV file specifying comparisons from disk.

            This implementation just adds all C{bsf.ngs.Sample} objects from the
            C{bsf.Analysis.collection} instance variable (i.e. C{bsf.ngs.Collection}) to the
            C{bsf.Analysis.sample_list} instance variable.
            @return:
            @rtype:
            """

            self.sample_list.extend(self.collection.get_all_samples())

            return

        # Start of the run() method body.

        super(Trimmomatic, self).run()

        # default = Default.get_global_default()

        if not os.path.isabs(self.adapter_path):
            self.adapter_path = os.path.join(
                os.path.dirname(self.classpath_trimmomatic),
                'adapters')

        run_adjust_illumina_clip_path(trimming_step_list=self.trimming_step_pe_list)
        run_adjust_illumina_clip_path(trimming_step_list=self.trimming_step_se_list)

        # Get the Trimmomatic tool Java Archive (JAR) class path directory.

        # if not self.classpath_trimmomatic:
        #     self.classpath_trimmomatic = default.classpath_trimmomatic

        run_read_comparisons()

        # Trimmomatic

        stage_trimmomatic = self.get_stage(name=self.stage_name_trimmomatic)

        for sample in self.sample_list:
            if self.debug > 0:
                print self, 'Sample name:', sample.name
                print sample.trace(level=1)

            sample_step_list = list()
            """ @type sample_step_list: list[str | unicode] """
            if 'Trimmomatic Steps' in sample.annotation_dict:
                for trimming_step in sample.annotation_dict['Trimmomatic Steps']:
                    sample_step_list.extend(
                        filter(lambda x: x != '', map(lambda x: x.strip(), trimming_step.split(','))))
                run_adjust_illumina_clip_path(trimming_step_list=sample_step_list)

            # The Trimmomatic analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=False)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print self, 'PairedReads name:', paired_reads.get_name()

                    paired_reads_step_list = list()
                    if 'Trimmomatic Steps' in paired_reads.annotation_dict:
                        for trimming_step in paired_reads.annotation_dict['Trimmomatic Steps']:
                            paired_reads_step_list.extend(
                                filter(lambda x: x != '', map(lambda x: x.strip(), trimming_step.split(','))))
                        run_adjust_illumina_clip_path(trimming_step_list=paired_reads_step_list)

                    # Apply some sanity checks.

                    # Maybe this case should be allowed after Trimmomatic trimming,
                    # where only the second Read survives.
                    if paired_reads.reads_2 and not paired_reads.reads_1:
                        raise Exception('PairedReads object with reads1 but no reads2 object.', UserWarning)

                    prefix_trimmomatic = '_'.join((stage_trimmomatic.name, paired_reads_name))

                    if self.debug > 0:
                        print 'Trimmomatic Prefix:', prefix_trimmomatic

                    file_path_trimmomatic = FilePathTrimmomatic(prefix=prefix_trimmomatic)

                    # Create a Runnable and an Executable for running the Trimmomatic analysis.

                    runnable_trimmomatic = self.add_runnable(
                        runnable=Runnable(
                            name=prefix_trimmomatic,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory,
                            file_path_object=file_path_trimmomatic))
                    self.set_stage_runnable(stage=stage_trimmomatic, runnable=runnable_trimmomatic)

                    # Create a new RunnableStepMakeDirectory in preparation of the Trimmomatic program.

                    runnable_trimmomatic.add_runnable_step(
                        runnable_step=RunnableStepMakeDirectory(
                            name='mkdir',
                            directory_path=file_path_trimmomatic.output_directory))

                    # Create a RunnableStep for the Trimmomatic program.

                    runnable_step_trimmomatic = runnable_trimmomatic.add_runnable_step(
                        runnable_step=RunnableStepJava(
                            name='trimmomatic',
                            java_temporary_path=runnable_trimmomatic.get_relative_temporary_directory_path,
                            java_heap_maximum='Xmx4G',
                            java_jar_path=self.classpath_trimmomatic))
                    """ @type runnable_step_trimmomatic: bsf.process.RunnableStepJava """

                    if paired_reads.reads_2 is None or not paired_reads.reads_2.name:
                        # FIXME: For the moment, PairedReads.reads2 is always defined.
                        runnable_step_trimmomatic.sub_command.sub_command = Command(program='SE')
                    else:
                        runnable_step_trimmomatic.sub_command.sub_command = Command(program='PE')

                    # Add options to the sub command.
                    sub_command = runnable_step_trimmomatic.sub_command.sub_command
                    sub_command.add_option_short(key='trimlog', value=file_path_trimmomatic.trim_log_tsv)

                    if paired_reads.reads_2 is None or not paired_reads.reads_2.name:
                        # FIXME: For the moment, PairedReads.reads2 is always defined.
                        file_path_trimmomatic.reads_1u = os.path.join(
                            file_path_trimmomatic.output_directory,
                            paired_reads.reads_1.name + 'U.fastq.gz')

                        sub_command.arguments.append(paired_reads.reads_1.file_path)
                        sub_command.arguments.append(file_path_trimmomatic.reads_1u)

                        # Update unpaired Reads information.

                        paired_reads.reads_1.name += 'U'
                        paired_reads.reads_1.file_path = os.path.join(
                            self.genome_directory,
                            file_path_trimmomatic.reads_1u)
                    else:
                        file_path_trimmomatic.reads_1p = os.path.join(
                            file_path_trimmomatic.output_directory,
                            paired_reads.reads_1.name + 'P.fastq.gz')
                        file_path_trimmomatic.reads_1u = os.path.join(
                            file_path_trimmomatic.output_directory,
                            paired_reads.reads_1.name + 'U.fastq.gz')
                        file_path_trimmomatic.reads_2p = os.path.join(
                            file_path_trimmomatic.output_directory,
                            paired_reads.reads_2.name + 'P.fastq.gz')
                        file_path_trimmomatic.reads_2u = os.path.join(
                            file_path_trimmomatic.output_directory,
                            paired_reads.reads_2.name + 'U.fastq.gz')

                        sub_command.arguments.append(paired_reads.reads_1.file_path)
                        sub_command.arguments.append(paired_reads.reads_2.file_path)
                        sub_command.arguments.append(file_path_trimmomatic.reads_1p)
                        sub_command.arguments.append(file_path_trimmomatic.reads_1u)
                        sub_command.arguments.append(file_path_trimmomatic.reads_2p)
                        sub_command.arguments.append(file_path_trimmomatic.reads_2u)

                        # Update paired Reads information.

                        paired_reads.reads_1.name += 'P'
                        paired_reads.reads_1.file_path = os.path.join(
                            self.genome_directory,
                            file_path_trimmomatic.reads_1p)
                        paired_reads.reads_2.name += 'P'
                        paired_reads.reads_2.file_path = os.path.join(
                            self.genome_directory,
                            file_path_trimmomatic.reads_2p)

                        # Add unpaired Reads 1 and 2 as separate PairedReads objects to this sample.

                        sample.add_paired_reads(
                            paired_reads=PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=Reads(
                                    name=paired_reads.reads_1.name[:-1] + 'U',
                                    file_path=os.path.join(
                                        self.genome_directory,
                                        file_path_trimmomatic.reads_1u)),
                                exclude=paired_reads.exclude,
                                index_1=paired_reads.index_1,
                                index_2=paired_reads.index_2,
                                read_group=paired_reads.read_group))

                        sample.add_paired_reads(
                            paired_reads=PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=Reads(
                                    name=paired_reads.reads_2.name[:-1] + 'U',
                                    file_path=os.path.join(
                                        self.genome_directory,
                                        file_path_trimmomatic.reads_2u)),
                                exclude=paired_reads.exclude,
                                index_1=paired_reads.index_1,
                                index_2=paired_reads.index_2,
                                read_group=paired_reads.read_group))

                    # Append trimming steps in order of specificity read group, sample and default.

                    if len(paired_reads_step_list):
                        sub_command.arguments.extend(paired_reads_step_list)
                    elif len(sample_step_list):
                        sub_command.arguments.extend(sample_step_list)
                    else:
                        if paired_reads.reads_2 is None or not paired_reads.reads_2.name:
                            # FIXME: For the moment, PairedReads.reads2 is always defined.
                            sub_command.arguments.extend(self.trimming_step_se_list)
                        else:
                            sub_command.arguments.extend(self.trimming_step_pe_list)

                    # Create a new RunnableStep to aggregate the trim log file.

                    runnable_step_trimmomatic_summary = runnable_trimmomatic.add_runnable_step(
                        runnable_step=RunnableStep(
                            name='trimmomatic_summary',
                            program='bsf_trimmomatic_summary.R',
                            obsolete_file_path_list=[
                                file_path_trimmomatic.trim_log_tsv,
                            ]))
                    """ @type runnable_step_trimmomatic_summary: bsf.process.RunnableStep """
                    runnable_step_trimmomatic_summary.add_option_long(
                        key='file_path',
                        value=file_path_trimmomatic.trim_log_tsv)

        # Convert the (modified) Collection object into a SampleAnnotationSheet object and write it to disk.

        annotation_sheet = self.collection.to_sas(
            file_path=os.path.join(
                self.project_directory,
                '_'.join((self.project_name, 'trimmomatic_samples.csv'))),
            name='_'.join((self.project_name, 'trimmomatic')))

        annotation_sheet.to_file_path()

        return annotation_sheet

    def report(self):
        """Create a C{bsf.analyses.trimmomatic.Trimmomatic} report in HTML format and a
        UCSC Genome Browser Track Hub.

        @return:
        @rtype:
        """

        # Create a symbolic link containing the project name and a UUID.
        self.create_public_project_link()

        # Write a HTML document.

        report_list = list()
        """ @type report_list: list[str | unicode] """

        report_list += '<h1 id="' + self.prefix + '_analysis">' + self.project_name + ' ' + self.name + '</h1>\n'
        report_list += '\n'

        report_list += '<h2 id="aliquot_and_sample_level">Aliquot and Sample Level</h2>\n'
        report_list += '\n'
        report_list += '<table id="aliquot_and_sample_table">\n'
        report_list += '<thead>\n'
        report_list += '<tr>\n'
        report_list += '<th>Sample</th>\n'
        report_list += '<th>Aliquot</th>\n'
        report_list += '<th>Coverage</th>\n'
        report_list += '<th>Frequency</th>\n'
        report_list += '<th>Summary</th>\n'
        report_list += '</tr>\n'
        report_list += '</thead>\n'
        report_list += '<tbody>\n'

        for sample in self.sample_list:
            # The Trimmomatic analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=False)

            paired_reads_name_list = paired_reads_dict.keys()
            if not len(paired_reads_name_list):
                # Skip Sample objects, which PairedReads objects have all been excluded.
                continue
            # paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            # This analysis is special in that read group names carry 'P' or 'U' suffices and samples carry additional
            # read groups after trimming that do no longer correspond to the initial Runnable objects. Sigh.
            paired_reads_name_list = dict(map(lambda x: (x[:-1], True), paired_reads_name_list)).keys()
            paired_reads_name_list.sort()

            report_list += '<tr>\n'
            report_list += '<td class="left">' + sample.name + '</td>\n'
            report_list += '<td class="left"></td>\n'  # Aliquot
            report_list += '<td class="center"></td>\n'  # Coverage PNG
            report_list += '<td class="center"></td>\n'  # Frequency PNG
            report_list += '<td class="center"></td>\n'  # Summary TSV
            report_list += '</tr>\n'

            for paired_reads_name in paired_reads_name_list:
                # The second read may still not be there.
                if '_'.join((self.stage_name_trimmomatic, paired_reads_name)) not in self.runnable_dict:
                    continue

                runnable_trimmomatic = self.runnable_dict[
                    '_'.join((self.stage_name_trimmomatic, paired_reads_name))]
                file_path_trimmomatic = runnable_trimmomatic.file_path_object
                """ @type file_path_trimmomatic: FilePathTrimmomatic """

                report_list += '<tr>\n'
                # Sample
                report_list += '<td class="left"></td>\n'
                # Aliquot
                report_list += '<td class="left">' + paired_reads_name + '</td>\n'
                # Coverage
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_trimmomatic.coverage_png + '">'
                report_list += '<img alt="Coverage ' + runnable_trimmomatic.name + '"'
                report_list += ' src="' + file_path_trimmomatic.coverage_png + '"'
                report_list += ' height="100" width="100" />'
                report_list += '</a>'
                report_list += '</td>\n'
                # Frequency
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_trimmomatic.frequency_png + '">PNG</a>'
                report_list += '</td>\n'
                # The frequency plots provide little information that does not necessarily justify
                # adding another set of images onto the HTML report.
                report_list += '<td class="center">'
                report_list += '<a href="' + file_path_trimmomatic.summary_tsv + '">TSV</a>'
                report_list += '</td>\n'
                report_list += '</tr>\n'

        report_list += '</tbody>\n'
        report_list += '</table>\n'
        report_list += '\n'

        self.report_to_file(content=report_list)

        return
