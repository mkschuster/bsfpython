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

from bsf import Analysis, Default, Runnable
from bsf.ngs import Reads, PairedReads
from bsf.process import Command, RunnableStep, RunnableStepJava, RunnableStepMakeDirectory
from bsf.standards import Configuration


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
            comparisons=None,
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
        @param comparisons: Python C{dict} of Python C{tuple} objects of C{bsf.ngs.Sample} objects
        @type comparisons: dict[str, tuple[bsf.ngs.Sample]]
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
            comparisons=comparisons,
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
        assert isinstance(configuration, Configuration)

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

        def _adjust_illumina_clip_path(trimming_step_list):
            """Adjust the adapter FASTA file path of ILLUMINACLIP trimming steps.

            If the file path is not absolute, prepend the value of the adapter_path
            instance variable.
            @param trimming_step_list: Python C{list} of trimming steps.
            @type trimming_step_list: list[str | unicode]
            @return:
            @rtype:
            """
            for i in range(1, len(trimming_step_list)):
                if trimming_step_list[i].startswith('ILLUMINACLIP'):
                    component_list = trimming_step_list[i].split(':')
                    if not os.path.isabs(component_list[2]):
                        component_list[2] = os.path.join(self.adapter_path, component_list[2])
                    trimming_step_list[i] = ':'.join(component_list)
            return

        def _read_comparisons():
            """Read a comparison annotation sheet to select those C{Sample} objects that should be analysed.

            By default, the C{Collection} object can automatically register all C{Sample} objects.
            @return:
            @rtype:
            """
            self.sample_list.extend(self.collection.get_all_samples())
            return

        super(Trimmomatic, self).run()

        # default = Default.get_global_default()

        if not os.path.isabs(self.adapter_path):
            self.adapter_path = os.path.join(
                os.path.dirname(self.classpath_trimmomatic),
                'adapters')

        _adjust_illumina_clip_path(trimming_step_list=self.trimming_step_pe_list)
        _adjust_illumina_clip_path(trimming_step_list=self.trimming_step_se_list)

        # Get the Trimmomatic tool Java Archive (JAR) class path directory.

        # if not self.classpath_trimmomatic:
        #     self.classpath_trimmomatic = default.classpath_trimmomatic

        _read_comparisons()

        # Trimmomatic

        stage_trimmomatic = self.get_stage(name=self.stage_name_trimmomatic)

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(level=1)

            sample_step_list = list()
            if 'Trimmomatic Steps' in sample.annotation_dict:
                for trimming_step in sample.annotation_dict['Trimmomatic Steps']:
                    sample_step_list.extend(
                        filter(lambda x: x != '', map(lambda x: x.strip(), trimming_step.split(','))))
                _adjust_illumina_clip_path(trimming_step_list=sample_step_list)

            # The Trimmomatic analysis does not obey excluded PairedReads objects,
            # more high-level analyses generally do.
            paired_reads_dict = sample.get_all_paired_reads(replicate_grouping=False, exclude=False)

            paired_reads_name_list = paired_reads_dict.keys()
            paired_reads_name_list.sort(cmp=lambda x, y: cmp(x, y))

            for paired_reads_name in paired_reads_name_list:
                for paired_reads in paired_reads_dict[paired_reads_name]:
                    if self.debug > 0:
                        print '{!r} PairedReads name: {}'.format(self, paired_reads.get_name())

                    paired_reads_step_list = list()
                    if 'Trimmomatic Steps' in paired_reads.annotation_dict:
                        for trimming_step in paired_reads.annotation_dict['Trimmomatic Steps']:
                            paired_reads_step_list.extend(
                                filter(lambda x: x != '', map(lambda x: x.strip(), trimming_step.split(','))))
                        _adjust_illumina_clip_path(trimming_step_list=paired_reads_step_list)

                    # Apply some sanity checks.

                    # Maybe this case should be allowed after Trimmomatic trimming,
                    # where only the second Read survives.
                    if paired_reads.reads_2 and not paired_reads.reads_1:
                        raise Exception('PairedReads object with reads1 but no reads2 object.', UserWarning)

                    prefix_trimmomatic = '_'.join((stage_trimmomatic.name, paired_reads_name))

                    if self.debug > 0:
                        print 'Trimmomatic Prefix: {}'.format(prefix_trimmomatic)

                    file_path_dict_trimmomatic = {
                        'temporary_directory': '_'.join((prefix_trimmomatic, 'temporary')),
                        'output_directory': os.path.join(self.project_directory, prefix_trimmomatic),
                        # Automatic GNU Zip-compression of trim log files does not work.
                        'trim_log_tsv': '_'.join((prefix_trimmomatic, 'trim_log.tsv')),
                        'summary_tsv': '_'.join((prefix_trimmomatic, 'summary.tsv')),  # Defined by the R script.
                        'coverage_png': '_'.join((prefix_trimmomatic, 'coverage.png')),  # Defined by the R script.
                        'frequency_png': '_'.join((prefix_trimmomatic, 'frequency.png')),  # Defined by the R script.
                        'surviving_png': '_'.join((prefix_trimmomatic, 'surviving.png')),  # Defined by the R script.
                    }

                    # Create a Runnable and an Executable for running the Trimmomatic analysis.

                    runnable_trimmomatic = self.add_runnable(
                        runnable=Runnable(
                            name=prefix_trimmomatic,
                            code_module='bsf.runnables.generic',
                            working_directory=self.project_directory,
                            file_path_dict=file_path_dict_trimmomatic))
                    self.set_stage_runnable(stage=stage_trimmomatic, runnable=runnable_trimmomatic)

                    # Create a new RunnableStepMakeDirectory in preparation of the Trimmomatic program.

                    runnable_trimmomatic.add_runnable_step(
                        runnable_step=RunnableStepMakeDirectory(
                            name='mkdir',
                            directory_path=file_path_dict_trimmomatic['output_directory']))

                    # Create a RunnableStep for the Trimmomatic program.

                    runnable_step_trimmomatic = runnable_trimmomatic.add_runnable_step(
                        runnable_step=RunnableStepJava(
                            name='trimmomatic',
                            java_temporary_path=file_path_dict_trimmomatic['temporary_directory'],
                            java_heap_maximum='Xmx4G',
                            java_jar_path=self.classpath_trimmomatic))

                    if paired_reads.reads_2 is None:
                        runnable_step_trimmomatic.sub_command.sub_command = Command(program='SE')
                    else:
                        runnable_step_trimmomatic.sub_command.sub_command = Command(program='PE')

                    # Add options to the sub command.
                    sub_command = runnable_step_trimmomatic.sub_command.sub_command
                    sub_command.add_option_short(key='trimlog', value=file_path_dict_trimmomatic['trim_log_tsv'])

                    if paired_reads.reads_2 is None:
                        file_path_1u = os.path.join(
                            file_path_dict_trimmomatic['output_directory'],
                            paired_reads.reads_1.name + 'U.fastq.gz')

                        sub_command.arguments.append(paired_reads.reads_1.file_path)
                        sub_command.arguments.append(file_path_1u)

                        # Update unpaired Reads information.

                        paired_reads.reads_1.name += 'U'
                        paired_reads.reads_1.file_path = file_path_1u
                    else:
                        file_path_1p = os.path.join(
                            file_path_dict_trimmomatic['output_directory'],
                            paired_reads.reads_1.name + 'P.fastq.gz')
                        file_path_1u = os.path.join(
                            file_path_dict_trimmomatic['output_directory'],
                            paired_reads.reads_1.name + 'U.fastq.gz')
                        file_path_2p = os.path.join(
                            file_path_dict_trimmomatic['output_directory'],
                            paired_reads.reads_2.name + 'P.fastq.gz')
                        file_path_2u = os.path.join(
                            file_path_dict_trimmomatic['output_directory'],
                            paired_reads.reads_2.name + 'U.fastq.gz')

                        sub_command.arguments.append(paired_reads.reads_1.file_path)
                        sub_command.arguments.append(paired_reads.reads_2.file_path)
                        sub_command.arguments.append(file_path_1p)
                        sub_command.arguments.append(file_path_1u)
                        sub_command.arguments.append(file_path_2p)
                        sub_command.arguments.append(file_path_2u)

                        # Update paired Reads information.

                        paired_reads.reads_1.name += 'P'
                        paired_reads.reads_1.file_path = file_path_1p
                        paired_reads.reads_2.name += 'P'
                        paired_reads.reads_2.file_path = file_path_2p

                        # Add unpaired Reads 1 and 2 as separate PairedReads objects to this sample.

                        # TODO: This needs to copy Reads and PairedReads objects.
                        sample.add_paired_reads(
                            paired_reads=PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=Reads(name=paired_reads.reads_1.name[:-1] + 'U', file_path=file_path_1u),
                                exclude=paired_reads.exclude,
                                index_1=paired_reads.index_1,
                                index_2=paired_reads.index_2,
                                read_group=paired_reads.read_group))

                        sample.add_paired_reads(
                            paired_reads=PairedReads(
                                annotation_dict=paired_reads.annotation_dict,
                                reads_1=Reads(name=paired_reads.reads_2.name[:-1] + 'U', file_path=file_path_2u),
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
                        if paired_reads.reads_2 is None:
                            sub_command.arguments.extend(self.trimming_step_se_list)
                        else:
                            sub_command.arguments.extend(self.trimming_step_pe_list)

                    # Create a new RunnableStep to aggregate the trim log file.

                    runnable_step_trimmomatic_summary = runnable_trimmomatic.add_runnable_step(
                        runnable_step=RunnableStep(
                            name='trimmomatic_summary',
                            program='bsf_trimmomatic_summary.R',
                            obsolete_file_path_list=[
                                file_path_dict_trimmomatic['trim_log_tsv'],
                            ]))

                    runnable_step_trimmomatic_summary.add_option_long(
                        key='file_path',
                        value=file_path_dict_trimmomatic['trim_log_tsv'])

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
        default = Default.get_global_default()
        link_path = self.create_public_project_link(sub_directory=default.url_relative_projects)
        link_name = os.path.basename(link_path.rstrip('/'))

        if link_name:
            pass  # Just make link_name used.

        # Write a HTML document.

        output_html = str()

        output_html += '<h1 id="{}_analysis">{} {}</h1>\n'.format(self.prefix, self.project_name, self.name)
        output_html += '\n'

        output_html += '<h2 id="aliquot_and_sample_level">Aliquot and Sample Level</h2>\n'
        output_html += '\n'
        output_html += '<table id="aliquot_and_sample_table">\n'
        output_html += '<thead>\n'
        output_html += '<tr>\n'
        output_html += '<th>Sample</th>\n'
        output_html += '<th>Aliquot</th>\n'
        output_html += '<th>Coverage</th>\n'
        output_html += '<th>Frequency</th>\n'
        output_html += '<th>Summary</th>\n'
        output_html += '</tr>\n'
        output_html += '</thead>\n'
        output_html += '<tbody>\n'

        for sample in self.sample_list:
            if self.debug > 0:
                print '{!r} Sample name: {}'.format(self, sample.name)
                print sample.trace(1)

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

            output_html += '<tr>\n'
            output_html += '<td class="left">{}</td>\n'.format(sample.name)
            output_html += '<td class="left"></td>\n'  # Aliquot
            output_html += '<td class="center"></td>\n'  # Coverage PNG
            output_html += '<td class="center"></td>\n'  # Frequency PNG
            output_html += '<td class="center"></td>\n'  # Summary TSV
            output_html += '</tr>\n'

            for paired_reads_name in paired_reads_name_list:
                # The second read may still not be there.
                if '_'.join((self.stage_name_trimmomatic, paired_reads_name)) not in self.runnable_dict:
                    continue

                runnable_trimmomatic = self.runnable_dict[
                    '_'.join((self.stage_name_trimmomatic, paired_reads_name))]
                assert isinstance(runnable_trimmomatic, Runnable)
                file_path_dict_trimmomatic = runnable_trimmomatic.file_path_dict

                output_html += '<tr>\n'
                # Sample
                output_html += '<td class="left"></td>\n'
                # Aliquot
                output_html += '<td class="left">{}</td>\n'.format(paired_reads_name)
                # Coverage
                output_html += '<td class="center">' \
                               '<a href="{}">' \
                               '<img alt="Coverage {}" src="{}" height="100" width="100" />' \
                               '</a>' \
                               '</td>\n'.format(file_path_dict_trimmomatic['coverage_png'],
                                                runnable_trimmomatic.name,
                                                file_path_dict_trimmomatic['coverage_png'])
                # Frequency
                output_html += '<td class="center"><a href="{}">PNG</a></td>\n'.format(
                    file_path_dict_trimmomatic['frequency_png'])
                # The frequency plots provide little information that does not necessarily justify
                # adding another set of images onto the HTML report.
                output_html += '<td class="center"><a href="{}">TSV</a></td>\n'.format(
                    file_path_dict_trimmomatic['summary_tsv'])
                output_html += '</tr>\n'

        output_html += '</tbody>\n'
        output_html += '</table>\n'
        output_html += '\n'

        self.report_to_file(content=output_html)

        return
