"""bsf.illumina

A package of classes and methods modelling data directories and files
specific for Illumina HiSeq systems.
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


import os.path
import string
import warnings
from xml.etree.ElementTree import ElementTree, Element


class RunInformationFlowcellLayout(object):
    """The C{RunInformationFlowcellLayout} class models
    one I{<FlowcellLayout>} XML element in an I{Illumina Run Information} (RunInfo.xml) document.

    Attributes:
    @ivar lane_count: Number of lanes
    @type lane_count: int
    @ivar surface_count: Number of surfaces
    @type surface_count: int
    @ivar swath_count: Number of swaths
    @type swath_count: int
    @ivar tile_count: Number of tiles
    @type tile_count: int
    """

    def __init__(self, lane_count=0, surface_count=0, swath_count=0, tile_count=0):
        """Initialise a C{RunInformationFlowcellLayout} object.

        @param lane_count: Number of lanes
        @type lane_count: int
        @param surface_count: Number of surfaces
        @type surface_count: int
        @param swath_count: Number of swaths
        @type swath_count: int
        @param tile_count: Number of tiles
        @type tile_count: int
        @return :
        @rtype :
        """

        self.lane_count = lane_count
        self.surface_count = surface_count
        self.swath_count = swath_count
        self.tile_count = tile_count


class RunInformationRead(object):
    """The C{RunInformationRead} class models
    one <Read> XML element in an Illumina Run Information (RunInfo.xml) document.

    Attributes:
    @ivar number: Read number
    @type number: int
    @ivar cycles: Cycle number
    @type cycles: int
    @ivar index: Index read
    @type index: bool
    """

    def __init__(self, number=0, cycles=0, index=False):
        """Initialise a C{RunInformationRead} object.

        @param number: Read number
        @type number: int
        @param cycles: Cycle number
        @type cycles: int
        @param index: Index read
        @type index: bool
        @return :
        @rtype :
        """

        self.number = number
        self.cycles = cycles
        self.index = index


class RunInformation(object):
    """The C{RunInformation} class represents an Illumina
    run information XML (RunInfo.xml) document.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        CASAVA: FASTQ file after post-processing with CASAVA.
        External: other data files.
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar run_identifier: Run identifier e.g. 130724_SN815_0089_BC26JBACXX
    @type run_identifier: str
    @ivar run_number: Run number, which may not have to correspond to the run number in the run identifier e.g. 91
    @type run_number: str
    @ivar flow_cell: Illumina flow-cell identifier e.g. C26JBACXX
    @type flow_cell: str
    @ivar instrument: Illumina instrument serial number e.g. SN815
    @type instrument: str
    @ivar date: Date in YYMMDD format e.g. 130724
    @type date: str
    @ivar reads: Python C{list} of C{RunInformationRead} objects
    @type reads: list
    """

    @staticmethod
    def parse_run_identifier(run_identifier):
        """Split an I{Illumina Run Identifier} into its components.

        Splits the Illumina Run Identifier into <Date>_<Instrument>_<Number>_<FCPosition><Flowcell>.
        This method is particularly useful for older version of RunInfo.xml files, that lack
        <Run>/<Date> and <Run>/<Flowcell> elements.
        @param run_identifier: Illumina Run Identifier (e.g. 130724_SN815_0089_BC26JBACXX)
        @type run_identifier: str
        @return : Python C{list} of Python C{str} objects
        @rtype : list
        """

        # Split into <Date>_<Instrument>_<Number>_<FCPosition><Flowcell>
        components = string.split(run_identifier, sep='_')

        if len(components) != 4:
            warnings.warn('Cannot split Illumina Run Identifier {!r} into its components.'.format(run_identifier))
            return

        # Strip leading zeros from the <Number>, split <FCPosition> and >Flowcell> elements.
        components[2] = components[2].lstrip('0')
        components.append(components[3][1:])
        components[3] = components[3][:1]

        return components

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{RunInformation} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return : C{RunInformation}
        @rtype : RunInformation
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        run_info_tree = ElementTree(file=file_path)
        run_element = run_info_tree.find(path='Run')
        if run_element is None:
            raise Exception('Cannot find the <Run> element in the ElementTree of XML file {!r}.'.format(file_path))
        assert isinstance(run_element, Element)

        # Parse meta-information about the Illumina run

        run_identifier = run_element.get(key='Id')  # e.g. 130724_SN815_0089_BC26JBACXX
        run_number = run_element.get(key='Number')  # e.g. 91
        run_identifier_components = RunInformation.parse_run_identifier(run_identifier=run_identifier)

        xml_flow_cell = run_element.find(path='Flowcell')
        if xml_flow_cell is not None:
            flow_cell = xml_flow_cell.text  # e.g. C26JBACXX
        else:
            flow_cell = run_identifier_components[4]

        xml_instrument = run_element.find(path='Instrument')
        if xml_instrument is not None:
            instrument = xml_instrument.text  # e.g. SN815
        else:
            instrument = run_identifier_components[1]

        xml_date = run_element.find(path='Date')
        if xml_date is not None:
            date = xml_date.text  # e.g. 130724
        else:
            date = run_identifier_components[0]

        xml_second_read = run_element.find(path='SecondRead')
        if xml_second_read is not None:
            second_read = int(xml_second_read.get(key='FirstCycle'))
        else:
            second_read = int(0)

        reads = list()
        number = int(1)

        for read_element in run_element.find(path='Reads'):
            assert isinstance(read_element, Element)

            # <ApplicationName>HiSeq Control Software</ApplicationName>
            # <ApplicationVersion>2.0.12.0</ApplicationVersion>
            # <ApplicationVersion>2.2.58</ApplicationVersion>
            #
            # <ApplicationName>MiSeq Control Software</ApplicationName>
            # <ApplicationVersion>2.5.0.5</ApplicationVersion>
            #
            # <Reads>
            #   <Read Number="1" NumCycles="100" IsIndexedRead="N" />
            #   <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
            #   <Read Number="3" NumCycles="100" IsIndexedRead="N" />
            # </Reads>

            if 'NumCycles' in read_element.keys():
                is_index = read_element.get(key='IsIndexedRead')
                assert isinstance(is_index, str)
                if is_index not in ('Y', 'N'):
                    warning = 'Unexpected value <Read IsIndexedRead="{}"> in Read element attribute IsIndexedRead '. \
                        format(read_element.get(key='IsIndexedRead'))
                    warnings.warn(warning)

                reads.append(RunInformationRead(
                    number=int(read_element.get(key='Number')),
                    cycles=int(read_element.get(key='NumCycles')),
                    index=bool(is_index == 'Y')))

            # <ApplicationName>HiSeq Control Software</ApplicationName>
            # <ApplicationVersion>1.1.37</ApplicationVersion>
            #
            # <Reads>
            #   <Read FirstCycle="1" LastCycle="51" />
            #   <Read FirstCycle="52" LastCycle="102" />
            # </Reads>
            # <SecondRead FirstCycle="52" />

            elif 'FirstCycle' in read_element.keys():

                # This is a bit convoluted. If a read after the first has a FirstCycle attribute of less than the
                # FirstCycle attribute of the SecondRead, then it must be the index read. Sigh!

                reads.append(RunInformationRead(
                    number=number,
                    cycles=int(read_element.get(key='LastCycle')) - int(read_element.get(key='FirstCycle')),
                    index=bool(number > 1 and int(read_element.get(key='FirstCycle')) < second_read)))

                number += 1

        reads.sort(key=lambda read: read.number)

        # Warn if there is not at least one non-index read.

        non_index_reads = filter(lambda read: not read.index, reads)
        if len(non_index_reads) == 0:
            warnings.warn(
                'No non-index read in Illumina RunInformation {!r}'.format(file_path),
                UserWarning)

        # Set a paired_end attribute if more than one read without index is defined?

        # Get the Flow-Cell Layout if it exits.

        xml_flow_cell_layout = run_info_tree.find(path='Run/FlowcellLayout')
        assert isinstance(xml_flow_cell_layout, Element)
        if xml_flow_cell_layout is not None:
            flow_cell_layout = RunInformationFlowcellLayout(
                lane_count=int(xml_flow_cell_layout.get(key='LaneCount')),
                surface_count=int(xml_flow_cell_layout.get(key='SurfaceCount')),
                swath_count=int(xml_flow_cell_layout.get(key='SwathCount')),
                tile_count=int(xml_flow_cell_layout.get(key='TileCount')))
        else:
            flow_cell_layout = None

        iri = cls(file_path=file_path, file_type='xml', name=file_name,
                  run_identifier=run_identifier, run_number=run_number,
                  flow_cell=flow_cell, instrument=instrument, date=date,
                  reads=reads, flow_cell_layout=flow_cell_layout)

        return iri

    def __init__(self, file_path=None, file_type=None, name=None,
                 run_identifier=None, run_number=None, flow_cell=None, instrument=None, date=None,
                 reads=None, flow_cell_layout=None):
        """Initialise a C{RunInformation} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str
        @param name: Name
        @type name: str
        @param run_identifier: Run identifier e.g. I{130724_SN815_0089_BC26JBACXX}
        @type run_identifier: str
        @param run_number: Run number, which may not have to correspond to the run number in the run identifier e.g. 91
        @type run_number: str
        @param flow_cell: Illumina flow-cell identifier
        @type flow_cell: str
        @param instrument: Illumina instrument serial number
        @type instrument: str
        @param date: Date in YYMMDD format
        @type date: str
        @param reads: Python C{list} of C{RunInformationRead} objects
        @type reads: list
        @param flow_cell_layout: C{RunInformationFlowcellLayout}
        @type flow_cell_layout: RunInformationFlowcellLayout
        @return :
        @rtype :
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if run_identifier:
            self.run_identifier = run_identifier
        else:
            self.run_identifier = str()

        if run_number:
            self.run_number = run_number
        else:
            self.run_number = str()

        if date:
            self.date = date
        else:
            self.date = str()

        if instrument:
            self.instrument = str()
        else:
            self.instrument = str()

        if flow_cell:
            self.flow_cell = flow_cell
        else:
            self.flow_cell = str()

        if reads:
            self.reads = reads
        else:
            self.reads = list()

        if flow_cell_layout:
            self.flow_cell_layout = flow_cell_layout
        else:
            self.flow_cell_layout = RunInformationFlowcellLayout()

    @property
    def get_cycle_number(self):
        """Get the total number of cycles.

        @return : Number of cycles
        @rtype : int
        """

        cycle_number = 0

        for read in self.reads:
            assert isinstance(read, RunInformationRead)
            cycle_number += read.cycles

        return cycle_number

    @property
    def get_read_number(self):
        """Get the total number of reads.

        @return : Number of reads
        @rtype : int
        """

        return len(self.reads)

    @property
    def get_read_start_list(self):
        """Get a Python C{list} of cycle numbers at the start of each read.

        @return : Python C{list} of starting cycle for each read
        @rtype : list
        """

        cycle_number = 0
        read_start_list = list()

        for read in self.reads:
            assert isinstance(read, RunInformationRead)
            read_start_list.append(cycle_number)
            cycle_number += read.cycles

        return read_start_list

    @property
    def get_picard_read_structure(self):
        """Get the read structure for the Picard C{ExtractIlluminaBarcodes} module.

        Codes:
        I{T} ... Template
        I{B} ... Barcode
        I{S} ... Skip
        @return : Read structure for Picard C{ExtractIlluminaBarcodes}
        @rtype : str
        """

        # TODO: This needs to become smarter to deal with skipped bases, caused by
        # different read and barcode lengths. Skips would have to be inserted.

        read_structure = str()

        for read in self.reads:
            assert isinstance(read, RunInformationRead)
            read_structure += str(read.cycles)
            if read.index:
                read_structure += 'B'
            else:
                read_structure += 'T'

        return read_structure

    @property
    def get_read_structure_list(self):
        """Get the read structure Python C{list}.

        Codes:
        B ... Base
        I ... Index
        @return : Python C{list} of Python C{str} (read structure) objects
        @rtype : list[str]
        """

        read_structure_list = list()

        for read in self.reads:
            assert isinstance(read, RunInformationRead)
            if read.index:
                read_structure_list.append('{}I'.format(read.cycles))
            else:
                read_structure_list.append('{}B'.format(read.cycles))

        return read_structure_list


class RunParameters(object):
    """The C{RunParameters} class models the contents of runParameters.xml
    files inside an Illumina Run Folder.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar element_tree: C{xml.etree.ElementTree}
    @type element_tree: ElementTree
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{RunParameters} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return : C{RunParameters} object
        @rtype : RunParameters
        """

        file_path = os.path.normpath(file_path)
        return cls(file_path=file_path, element_tree=ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{RunParameters} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: XML Element Tree
        @type element_tree: ElementTree
        @return :
        @rtype :
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if element_tree:
            self.element_tree = element_tree
        else:
            self.element_tree = ElementTree()

        self._run_parameters_version = None

    @property
    def get_run_parameters_version(self):
        """Get the run parameters version of a C{RunParameters} object.

        Returns an empty string for I{HiSeq Control Software} or the text representation of the
        I{<RunParameters>/<RunParametersVersion>} element for I{MiSeq Control Software}.
        @return : Run parameters version or an empty string
        @rtype : str
        """

        if self._run_parameters_version is None:
            element = self.element_tree.find(path='RunParametersVersion')
            if element is None:
                self._run_parameters_version = str()
            else:
                self._run_parameters_version = element.text

        return self._run_parameters_version

    @property
    def get_experiment_name(self):
        """Get the experiment name of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<ExperimentName>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<ExperimentName>} element for
        I{MiSeq Control Software}.
        @return : Experiment name e.g I{BSF_0001}
        @rtype : str
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return self.element_tree.find(path='ExperimentName').text
        else:
            # HiSeq run or else
            return self.element_tree.find(path='Setup/ExperimentName').text

    @property
    def get_flow_cell_barcode(self):
        """Get the flow-cell barcode of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Barcode>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<Barcode>} element for
        I{MiSeq Control Software}.
        @return : Flow-cell barcode e.g BSF_0001
        @rtype : str
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return self.element_tree.find(path='Barcode').text
        else:
            # HiSeq run or else
            return self.element_tree.find(path='Setup/Barcode').text

    @property
    def get_flow_cell_type(self):
        """Get the flow cell of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Flowcell>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<Flowcell>} element for
        I{MiSeq Control Software}.
        @return : Flow-cell type
        @rtype : str
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            # The MiSeq has no concept of Flowcell chemistry version, only a <ReagentKitVersion>.
            return str()
        else:
            # HiSeq run or else
            return self.element_tree.find(path='Setup/Flowcell').text

    @property
    def get_position(self):
        """Get the flow-cell position of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<FCPosition>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Since the element does not exist in older I{HiSeq Control Software} versions an empty string may be returned.
        @return : Flow-cell position e.g. A or B
        @rtype : str
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            # The MiSeq has no concept of Flowcell chemistry version, only a <ReagentKitVersion>.
            return str()
        else:
            # HiSeq run or else
            element = self.element_tree.find(path='Setup/FCPosition')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_run_identifier(self):
        """Get the run identifier of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<RunID>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<RunID>} element for
        I{MiSeq Control Software}.
        @return : Run identifier
        @rtype : str
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return self.element_tree.find(path='RunID').text
        else:
            # HiSeq run or else
            return self.element_tree.find(path='Setup/RunID').text

    @property
    def get_read1(self):
        """Get the read 1 cycle number of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Read1>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        @return : Number of cycles in read 1
        @rtype : str
        @deprecated : The more scalable option is getting read information via the
            C{RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return str()
        else:
            # HiSeq run or else
            element = self.element_tree.find(path='Setup/Read1')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_read2(self):
        """Get the read 2 cycle number of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Read2>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Not every run has a read 2 defined so that an empty string may be returned.
        @return : Number of cycles in read 2
        @rtype : str
        @deprecated : The more scalable option is getting read information via the
            C{RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return str()
        else:
            # HiSeq run or else
            element = self.element_tree.find(path='Setup/Read2')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_index_read1(self):
        """Get the index read 1 cycle number of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<IndexRead1>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Older implementations of the I{HiSeq Control Software} have only a
        I{<RunParameters>/<Setup>/<IndexRead>} element.
        @return : Number of cycles in index read 1
        @rtype : str
        @deprecated : The more scalable option is getting read information via the
            C{RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return str()
        else:
            # HiSeq run or else
            element = self.element_tree.find(path='Setup/IndexRead1')
            if element is not None:
                return element.text

            element = self.element_tree.find(path='Setup/IndexRead')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_index_read2(self):
        """Get the index read 2 cycle number of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<IndexRead2>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Older implementations of the I{HiSeq Control Software} have only a
        I{<RunParameters>/<Setup>/<IndexRead>} element.
        @return : Number of cycles in index read 2
        @rtype : str
        @deprecated : The more scalable option is getting read information via the
            C{RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return str()
        else:
            # HiSeq run or else
            element = self.element_tree.find(path='Setup/IndexRead2')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_real_time_analysis_version(self):
        """Get the Real-Time Analysis (RTA) Version of a C{RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<RTAVersion>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<RTAVersion>} element for
        I{MiSeq Control Software}.
        @return : RTA version
        @rtype : str
        """

        if self.get_run_parameters_version == 'MiSeq_1_1':
            # MiSeq 1_1 run
            return self.element_tree.find(path='RTAVersion').text
        else:
            # HiSeq run or else
            return self.element_tree.find(path='Setup/RTAVersion').text

    @property
    def get_application_name(self):
        """Get the application (i.e. I{HiSeq} or I{MiSeq Control Software}) name.

        Returns the text representation of the I{<RunParameters>/<Setup>/<ApplicationName>} element.
        @return : Application name
        @rtype : str
        """

        return self.element_tree.find(path='Setup/ApplicationName').text

    @property
    def get_application_version(self):
        """Get the application (i.e. I{HiSeq} or I{MiSeq Control Software}) version.

        Returns the text representation of the I{<RunParameters>/<Setup>/<ApplicationVersion>} element.
        @return : Application version
        @rtype : str
        """

        return self.element_tree.find(path='Setup/ApplicationVersion').text


class XMLConfiguration(object):
    """The C{XMLConfiguration} class models the contents of XML configuration
    files inside an Illumina Run Folder.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar element_tree: C{xml.etree.ElementTree}
    @type element_tree: ElementTree
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{XMLConfiguration} object from a file path.
        In case the file path does not exist a C{XMLConfiguration} object
        with an empty C{ElementTree} will be returned.

        @param file_path: File path
        @type file_path: str | unicode
        @return : C{XMLConfiguration} object
        @rtype : XMLConfiguration
        """
        file_path = os.path.normpath(file_path)
        if not os.path.isfile(file_path):
            file_path = None
        return cls(file_path=file_path, element_tree=ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{XMLConfiguration} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: XML Element Tree
        @type element_tree: ElementTree
        @return :
        @rtype :
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if element_tree:
            self.element_tree = element_tree
        else:
            self.element_tree = ElementTree()


class AnalysisConfiguration(XMLConfiguration):
    """The C{AnalysisConfiguration} class models Image and Base Call analysis
    XML configuration files inside and Illumina Run Folder.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar element_tree: C{xml.etree.ElementTree}
    @type element_tree: ElementTree
    """

    def __init__(self, file_path=None, element_tree=None):
        """Initialise an C{AnalysisConfiguration} object.
        Read all I{<Tile>} elements from the I{<Run>/<TileSelection>} element from the XML configuration file
        to initialise an internal Python C{dict} of valid tiles.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: XML Element Tree
        @type element_tree: ElementTree
        @return :
        @rtype :
        """

        super(AnalysisConfiguration, self).__init__(file_path=file_path, element_tree=element_tree)

        if self.element_tree.getroot() is None:
            return

        self._lane_tile_dict = dict()

        for lane_element in self.element_tree.find(path='Run/TileSelection'):
            assert isinstance(lane_element, Element)
            lane_index = lane_element.get(key='Index')
            assert isinstance(lane_index, str)
            if lane_index not in self._lane_tile_dict:
                self._lane_tile_dict[lane_index] = dict()
            lane_dict = self._lane_tile_dict[lane_index]
            for tile_element in lane_element.findall(path='Tile'):
                assert isinstance(tile_element, Element)
                lane_dict[tile_element.text] = True

    def has_lane(self, lane):
        """Check if a particular lane is defined in an C{AnalysisConfiguration} object.

        @param lane: Lane index
        @type lane: str
        @return : Boolean value
        @rtype : bool
        """
        if lane in self._lane_tile_dict:
            return True
        else:
            return False

    def has_lane_tile(self, lane, tile):
        """Check if a particular tile is defined in a lane of an C{AnalysisConfiguration} object.

        @param lane: Lane index
        @type lane: str
        @param tile: Tile index
        @type tile: str
        @return : Boolean value
        @rtype : bool | None
        """
        if lane in self._lane_tile_dict:
            lane_dict = self._lane_tile_dict[lane]
            if tile in lane_dict:
                return True
            else:
                return False
        else:
            return None


class ImageAnalysis(AnalysisConfiguration):
    """The C{ImageAnalysis} class models the contents of the I{IRF/Data/Intensities/config.xml}
    XML configuration file inside an Illumina Run Folder.

    Attributes:
    """

    pass


class BaseCallAnalysis(AnalysisConfiguration):
    """The C{BaseCallAnalysis} class models the contents of the I{IRF/Data/Intensities/BaseCalls/config.xml}
    XML configuration file inside an Illumina Run Folder.

    Attributes:
    """

    pass


class RunFolderNotComplete(Exception):
    pass


class RunFolder(object):
    """The C{RunFolder} class represents an Illumina
    Run Folder copied off the instrument.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @type file_type: str
    @ivar name: Name
    @type name: str
    @ivar date: Date in YYMMDD format
    @type date: str
    @ivar instrument: Illumina instrument serial number
    @type instrument: str
    @ivar run: Run serial number
    @type run: str
    @ivar flow_cell: Flow-cell identifier
    @type flow_cell: str
    @ivar run_information: C{RunInformation} object
    @type run_information: RunInformation
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Construct a C{RunFolder} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return : C{RunFolder} object
        @rtype : RunFolder
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        # Illumina Run Folders obey a "YYMMDD_SN000_Run_PositionFCID"
        # schema.

        components = file_name.split('_')

        irf = cls(
            file_path=file_path,
            file_type='Illumina',
            name=string.join(words=components[:], sep='_'),
            date=components[0],
            instrument=components[1],
            run=components[2],
            flow_cell=components[3],
            run_information=RunInformation.from_file_path(file_path=os.path.join(file_path, 'RunInfo.xml')),
            run_parameters=RunParameters.from_file_path(file_path=os.path.join(file_path, 'runParameters.xml')),
            image_analysis=ImageAnalysis.from_file_path(file_path=os.path.join(
                file_path, 'Data', 'Intensities', 'config.xml')),
            base_call_analysis=BaseCallAnalysis.from_file_path(file_path=os.path.join(
                file_path, 'Data', 'Intensities', 'BaseCalls', 'config.xml')))

        return irf

    def __init__(self, file_path=None, file_type=None, name=None,
                 date=None, instrument=None, run=None, flow_cell=None,
                 run_information=None, run_parameters=None,
                 image_analysis=None, base_call_analysis=None):
        """Initialise a C{RunFolder} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str
        @param name: Name
        @type name: str
        @param date: Date in I{YYMMDD} format
        @type date: str
        @param instrument: Illumina instrument serial number (e.g. I{SN181}, I{SN815}, ...)
        @type instrument: str
        @param run: Run serial number
        @type run: str
        @param flow_cell: The position and flow cell identifier
        @type flow_cell: str
        @param run_information: C{RunInformation}
        @type run_information: RunInformation
        @param image_analysis: C{ImageAnalysis}
        @type image_analysis: ImageAnalysis
        @param base_call_analysis: C{BaseCallAnalysis}
        @type base_call_analysis: BaseCallAnalysis
        @return :
        @rtype :
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if file_type:
            self.file_type = file_type
        else:
            self.file_type = str()

        if name:
            self.name = name
        else:
            self.name = str()

        if date:
            self.date = date
        else:
            self.date = str()

        if instrument:
            self.instrument = instrument
        else:
            self.instrument = str()

        if run:
            self.run = run
        else:
            self.run = str()

        if flow_cell:
            self.flow_cell = flow_cell
        else:
            self.flow_cell = str()

        if run_information:
            self.run_information = run_information
        else:
            self.run_information = RunInformation()

        if run_parameters:
            self.run_parameters = run_parameters
        else:
            self.run_parameters = RunParameters()

        if image_analysis:
            self.image_analysis = image_analysis
        else:
            self.image_analysis = ImageAnalysis()

        if base_call_analysis:
            self.base_call_analysis = base_call_analysis
        else:
            self.base_call_analysis = BaseCallAnalysis()

        self._missing_base_call_tiles = dict()
        self._missing_image_analysis_tiles = dict()

    @property
    def get_base_calls_directory(self):
        """Get the base-calls directory in the I{IRF/Data/Intensities/BaseCalls} hierarchy.

        @return : Illumina base-calls directory
        @rtype : str | unicode
        """

        return os.path.join(self.file_path, 'Data', 'Intensities', 'BaseCalls')

    def _check_tiles_base_call(self):
        """Check for missing I{<Tile>} elements in the I{IRF/Data/Intensities/BaseCalls/config.xml}
        configuration file. This method also builds up a Python C{dict} required for method
        C{_is_missing_base_call_tile}.

        @return :
        @rtype :
        """
        fcl = self.run_information.flow_cell_layout

        for lane in range(1, fcl.lane_count + 1):
            if lane not in self._missing_base_call_tiles:
                self._missing_base_call_tiles[lane] = dict()
            lane_dict = self._missing_base_call_tiles[lane]
            if self.base_call_analysis.has_lane(lane='{:1d}'.format(lane)):
                # Lanes are not defined of the config.xml file could not be read.
                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            tile = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if not self.base_call_analysis.has_lane_tile(lane='{:1d}'.format(lane), tile=tile):
                                lane_dict[tile] = True
                                print "Missing BaseCallAnalysis lane '{:1d}' tile '{}'.".\
                                    format(lane, tile)

    def _check_tiles_image_analysis(self):
        """Check for missing I{<Tile>} elements in the I{IRF/Data/Intensities/config.xml}
        configuration file. This method also builds up a Python C{dict} required for method
        C{_is_missing_image_analysis_tile}.

        @return :
        @rtype :
        """
        fcl = self.run_information.flow_cell_layout

        for lane in range(1, fcl.lane_count + 1):
            if lane not in self._missing_image_analysis_tiles:
                self._missing_image_analysis_tiles[lane] = dict()
            lane_dict = self._missing_image_analysis_tiles[lane]
            if self.image_analysis.has_lane(lane='{:1d}'.format(lane)):
                # Lanes are not defined of the config.xml file could not be read.
                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            tile = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if not self.image_analysis.has_lane_tile(lane='{:1d}'.format(lane), tile=tile):
                                lane_dict[tile] = True
                                print "Missing ImageAnalysis lane '{:1d}' tile '{}'.".\
                                    format(lane, tile)

    def _is_missing_base_call_tile(self, lane, tile):
        """Confirm that a particular I{<Tile>} element is missing from the
        I{IRF/Data/Intensities/config.xml} configuration file.

        @param lane: Lane index
        @type lane: int
        @param tile: Tile name
        @type tile: str
        @return : Boolean value
        @rtype : bool
        """
        if lane in self._missing_base_call_tiles:
            lane_dict = self._missing_base_call_tiles[lane]
            if tile in lane_dict:
                return True
            else:
                return False
        else:
            return False

    def _is_missing_image_analysis_tile(self, lane, tile):
        """Confirm that a particular I{<Tile>} element is missing from the
        I{IRF/Data/Intensities/BaseCalls/config.xml} configuration file.

        @param lane: Lane index
        @type lane: int
        @param tile: Tile name
        @type tile: str
        @return : Boolean value
        @rtype : bool
        """
        if lane in self._missing_image_analysis_tiles:
            lane_dict = self._missing_image_analysis_tiles[lane]
            if tile in lane_dict:
                return True
            else:
                return False
        else:
            return False

    @staticmethod
    def _check_files(directory_dict, directory_path, file_list, debug=0):
        """Check a Python C{list} of file names against a Python C{dict} of directory entries.

        @param directory_dict: Python C{dict} of directory entries
        @type directory_dict: dict
        @param directory_path: Directory path
        @type directory_path: str | unicode
        @param file_list: Python C{list} of file names
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """
        if debug:
            print "Processing directory '{}'".format(directory_path)

        for file_name in file_list:
            if file_name in directory_dict:
                del directory_dict[file_name]
            else:
                print "Missing file '{}' in directory path '{}'".format(file_name, directory_path)

    def _check_data_intensities_base_calls_matrix(self, base_calls_dict, base_calls_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/Matrix/} directory.

        @param base_calls_dict: Python C{dict} of I{IRF/Data/Intensities/BaseCalls/} entries
        @type base_calls_dict: dict
        @param base_calls_path: I{IRF/Data/Intensities/BaseCalls/} path
        @type base_calls_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        matrix_name = 'Matrix'
        if matrix_name in base_calls_dict:
            del base_calls_dict[matrix_name]
        else:
            print "Missing matrix directory: '{}'".format(matrix_name)
            return

        matrix_path = os.path.join(base_calls_path, matrix_name)
        matrix_dict = dict(map(lambda x: (x, 1), os.listdir(matrix_path)))

        if debug:
            print "Processing '{}'".format(matrix_name)

        if rta in '2.5.2':
            # For the HiSeq 3000 platform, process IRF/Data/Intensities/BaseCalls/Matrix/L00[1-8] directories.

            for lane in range(1, fcl.lane_count + 1):

                lane_name = 'L{:03d}'.format(lane)
                if lane_name in matrix_dict:
                    del matrix_dict[lane_name]
                else:
                    print "Missing lane directory: 'L{:03d}'".format(lane)
                    continue

                lane_path = os.path.join(matrix_path, lane_name)
                lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

                # Process IRF/Data/Intensities/BaseCalls/Matrix/L00[1-8]/C[0-9]+.1/ directories.

                for cycle in range(1, self.run_information.get_cycle_number + 1):

                    cycle_name = 'C{:d}.1'.format(cycle)
                    if cycle_name in lane_dict:
                        del lane_dict[cycle_name]
                    else:
                        print "Missing cycle directory '{}' '{}'".format(lane_name, cycle_name)
                        continue

                    cycle_path = os.path.join(lane_path, cycle_name)
                    cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                    for surface in range(1, fcl.surface_count + 1):
                        for swath in range(1, fcl.swath_count + 1):
                            for tile in range(1, fcl.tile_count + 1):
                                # Process tile matrix files.
                                # s_1_1101_matrix.txt
                                # s_1_2228_matrix.txt
                                tile_file = 's_{:1d}_{:1d}{:1d}{:02d}_matrix.txt'.format(lane, surface, swath, tile)
                                if tile_file in cycle_dict:
                                    del cycle_dict[tile_file]
                                else:
                                    print "Missing tile matrix file '{}' '{}' '{}'".\
                                        format(lane_name, cycle_name, tile_file)
        else:
            # HiSeq 2000 has a flat list of matrix.txt files.
            for read in range(1, self.run_information.get_read_number + 1):
                # Process read matrix files.
                # s_1_matrix.txt
                # s_2_matrix.txt
                read_file = 's_{:d}_matrix.txt'.format(read)
                if read_file in matrix_dict:
                    del matrix_dict[read_file]
                else:
                    print "Missing read matrix file '{}' '{}'".format(matrix_name, read_file)

                for lane in range(1, fcl.lane_count + 1):
                    # Process lane matrix files.
                    # s_1_1_matrix.txt
                    # s_8_2_matrix.txt
                    lane_file = 's_{:d}_{:d}_matrix.txt'.format(lane, read)
                    if lane_file in matrix_dict:
                        del matrix_dict[lane_file]
                    else:
                        print "Missing lane matrix file '{}' '{}'".format(matrix_name, lane_file)
                    for surface in range(1, fcl.surface_count + 1):
                        for swath in range(1, fcl.swath_count + 1):
                            for tile in range(1, fcl.tile_count + 1):
                                # Not all tiles have to exists especially after catastrophic events during the
                                # cluster generation step.
                                tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                                if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                                    continue
                                # Process tile matrix files.
                                # s_1_1_1101_matrix.txt
                                # s_8_2_2311_matrix.txt
                                tile_file = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_matrix.txt'.\
                                    format(lane, read, surface, swath, tile)
                                if tile_file in matrix_dict:
                                    del matrix_dict[tile_file]
                                else:
                                    print "Missing tile matrix file '{}' '{}'".\
                                        format(matrix_name, tile_file)

        # The matrix_dict should now be empty.

        if len(matrix_dict):
            print "Matrix '{}' with number of entries: {:d}".format(matrix_name, len(matrix_dict))
            entry_names = matrix_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_data_intensities_base_calls_phasing(self, base_calls_dict, base_calls_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/Phasing/} directory.

        @param base_calls_dict: Python C{dict} of I{IRF/Data/intensities/BaseCalls/} entries
        @type base_calls_dict: dict
        @param base_calls_path: I{IRF/Data/intensities/BaseCalls/} path
        @type base_calls_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        phasing_name = 'Phasing'
        if phasing_name in base_calls_dict:
            del base_calls_dict[phasing_name]
        else:
            print "Missing phasing directory: '{}'".format(phasing_name)
            return

        phasing_path = os.path.join(base_calls_path, phasing_name)
        phasing_dict = dict(map(lambda x: (x, 1), os.listdir(phasing_path)))

        if debug:
            print "Processing '{}'".format(phasing_name)

        for ri_read in self.run_information.reads:
            assert isinstance(ri_read, RunInformationRead)
            # Process read phasing files.
            # s_1_phasing.txt
            # s_2_phasing.txt
            read_file = 's_{:d}_phasing.txt'.format(ri_read.number)
            if read_file in phasing_dict:
                del phasing_dict[read_file]
            else:
                print "Missing read phasing file '{}' '{}'".format(phasing_name, read_file)

            for lane in range(1, fcl.lane_count + 1):
                # Process lane phasing files.
                # s_1_1_phasing.txt
                # s_8_3_phasing.txt
                lane_file = 's_{:d}_{:d}_phasing.txt'.format(lane, ri_read.number)
                if lane_file in phasing_dict:
                    del phasing_dict[lane_file]
                else:
                    print "Missing lane phasing file '{}' '{}'".format(phasing_name, lane_file)
                    continue

                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            # Not all tiles have to exists especially after catastrophic events during the
                            # cluster generation step.
                            tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                                continue
                            if rta in ('1.12.4.2', '1.13.48', '1.17.21.3') and not ri_read.index:
                                # Process tile cycle files, which only exist for payload, but not index reads.
                                # s_1_1_1101_cycle.txt
                                # s_8_1_2316_cycle.txt
                                cycle_file = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_cycle.txt'.\
                                    format(lane, ri_read.number, surface, swath, tile)
                                if cycle_file in phasing_dict:
                                    del phasing_dict[cycle_file]
                                else:
                                    print "Missing tile cycle file '{}' '{}'".\
                                        format(phasing_name, cycle_file)

                            # Process tile phasing files.
                            # s_1_1_1101_phasing.txt
                            # s_8_3_2316_phasing.txt
                            tile_file = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_phasing.txt'.\
                                format(lane, ri_read.number, surface, swath, tile)
                            if tile_file in phasing_dict:
                                del phasing_dict[tile_file]
                            else:
                                print "Missing tile phasing file '{}' '{}'".\
                                    format(phasing_name, tile_file)

                            if rta not in ('1.12.4.2', '1.13.48', '1.17.21.3'):
                                # Process the tile empirical phasing files.
                                tile_file = 'EmpiricalPhasingCorrection_{:1d}_{:1d}_{:1d}{:1d}{:02d}.txt'.\
                                    format(lane, ri_read.number, surface, swath, tile)
                                if tile_file in phasing_dict:
                                    del phasing_dict[tile_file]
                                else:
                                    print "Missing tile empirical phasing file '{}' '{}'".\
                                        format(phasing_name, tile_file)

        if rta in '2.5.2':
            # HiSeq 3000 does not have IRF/DataIntensities/BaseCalls/Phasing/s_lane_cycle_phasing.xml files.
            pass
        else:
            for lane in range(1, fcl.lane_count + 1):
                for read_start in self.run_information.get_read_start_list:
                    file_name = 's_{:1d}_{:02d}_phasing.xml'.format(lane, read_start + 2)
                    if file_name in phasing_dict:
                        del phasing_dict[file_name]
                    else:
                        print "Missing phasing XML file '{}'".format(file_name)

        # The phasing_dict should now be empty.

        if len(phasing_dict):
            print "Phasing '{}' with number of entries: {:d}".format(phasing_name, len(phasing_dict))
            entry_names = phasing_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_data_intensities_base_calls(self, intensities_dict, intensities_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/} directory.

        @param intensities_dict: Python C{dict} of I{IRF/Data/Intensities/} entries
        @type intensities_dict: dict
        @param intensities_path: I{IRF/Data/Intensities/} file path
        @type intensities_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        base_calls_name = 'BaseCalls'
        if base_calls_name in intensities_dict:
            del intensities_dict[base_calls_name]
        else:
            print "Missing base calls directory: '{}'".format(base_calls_name)
            return

        base_calls_path = os.path.join(intensities_path, base_calls_name)
        base_calls_dict = dict(map(lambda x: (x, 1), os.listdir(base_calls_path)))

        if debug:
            print "Processing '{}'".format(base_calls_name)

        # Process the IRF/Data/Intensities/BaseCalls/config.xml file.

        if rta in '2.5.2':
            # HiSeq 3000 does not have the IRF/DataIntensities/BaseCalls/config.xml file.
            pass
        else:
            config_name = 'config.xml'
            if config_name in base_calls_dict:
                del base_calls_dict[config_name]
            else:
                print "Missing configuration file: '{}'".format(config_name)

            # Check for completeness of tiles in the base call configuration XML file (config.xml).
            self._check_tiles_base_call()

        # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/ directories.

        for lane in range(1, fcl.lane_count + 1):

            lane_name = 'L{:03d}'.format(lane)
            if lane_name in base_calls_dict:
                del base_calls_dict[lane_name]
            else:
                print "Missing lane directory: 'L{:03d}'".format(lane)
                continue

            lane_path = os.path.join(base_calls_path, lane_name)
            lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

            # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/C[0-9]+.1 directories.

            for cycle in range(1, self.run_information.get_cycle_number + 1):

                cycle_name = 'C{:d}.1'.format(cycle)
                if cycle_name in lane_dict:
                    del lane_dict[cycle_name]
                else:
                    print "Missing cycle directory '{}' '{}'".format(lane_name, cycle_name)
                    continue

                cycle_path = os.path.join(lane_path, cycle_name)
                cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            # Not all tiles have to exists especially after catastrophic events during the
                            # cluster generation step.
                            tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                                continue
                            # Process tile base call (BCL) files.
                            # s_1_1101.bcl
                            # s_1_2316.bcl
                            tile_file = 's_{:1d}_{:1d}{:1d}{:02d}.bcl'.format(lane, surface, swath, tile)
                            if tile_file in cycle_dict:
                                del cycle_dict[tile_file]
                            else:
                                # Base call (BCL) files can also be GNU Zip compressed.
                                # s_1_1101.bcl.gz
                                # s_1_2316.bcl.gz
                                tile_file += '.gz'
                                if tile_file in cycle_dict:
                                    del cycle_dict[tile_file]
                                else:
                                    print "Missing tile bcl file '{}' '{}' '{}'".\
                                        format(lane_name, cycle_name, tile_file)

                            # Process tile stats files.
                            # s_1_1101.stats
                            # s_1_2316.stats
                            if rta in '2.5.2':
                                # HiSeq 3000 does not have stats files.
                                pass
                            else:
                                tile_file = 's_{:1d}_{:1d}{:1d}{:02d}.stats'.format(lane, surface, swath, tile)
                                if tile_file in cycle_dict:
                                    del cycle_dict[tile_file]
                                else:
                                    print "Missing tile stats file '{}' '{}' '{}'".\
                                        format(lane_name, cycle_name, tile_file)

                # The cycle_dict should now be empty.

                if len(cycle_dict):
                    print "Lane '{}' Cycle '{}' with number of entries: {:d}".\
                        format(lane_name, cycle_name, len(cycle_dict))
                    entry_names = cycle_dict.keys()
                    entry_names.sort(cmp=lambda x, y: cmp(x, y))
                    print '  Remaining entries: {!r}'.format(entry_names)

            # Process control and filter files.

            for surface in range(1, fcl.surface_count + 1):
                for swath in range(1, fcl.swath_count + 1):
                    for tile in range(1, fcl.tile_count + 1):
                        # Not all tiles have to exists especially after catastrophic events during the
                        # cluster generation step.
                        tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                        if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                            continue
                        # Process tile control files.
                        # s_1_1101.control
                        # s_1_2316.control
                        if rta in '2.5.2':
                            # HiSeq 3000 does not have control files.
                            pass
                        else:
                            tile_file = 's_{:1d}_{:1d}{:1d}{:02d}.control'.format(lane, surface, swath, tile)
                            if tile_file in lane_dict:
                                del lane_dict[tile_file]
                            else:
                                print "Missing tile control file '{}' '{}'".\
                                    format(lane_name, tile_file)

                        # Process tile filter files.
                        # s_1_1101.filter
                        # s_1_2316.filter
                        tile_file = 's_{:1d}_{:1d}{:1d}{:02d}.filter'.format(lane, surface, swath, tile)
                        if tile_file in lane_dict:
                            del lane_dict[tile_file]
                        else:
                            print "Missing tile filter file '{}' '{}'".\
                                format(lane_name, tile_file)

            # The lane_dict should now be empty.

            if len(lane_dict):
                print "Data/Intensities/BaseCalls/Lane '{}' with number of entries: {:d}".\
                    format(lane_name, len(lane_dict))
                entry_names = lane_dict.keys()
                entry_names.sort(cmp=lambda x, y: cmp(x, y))
                print '  Remaining entries: {!r}'.format(entry_names)

        # Process the IRF/Data/Intensities/BaseCalls/Matrix/ directory.

        self._check_data_intensities_base_calls_matrix(
            base_calls_dict=base_calls_dict,
            base_calls_path=base_calls_path,
            debug=debug)

        # Process the IRF/Data/Intensities/BaseCalls/Phasing/ directory.

        if rta in '2.5.2':
            # HiSeq 3000 does not have a IRF/Data/Intensities/BaseCalls/Phasing/ directory.
            pass
        else:
            self._check_data_intensities_base_calls_phasing(
                base_calls_dict=base_calls_dict,
                base_calls_path=base_calls_path,
                debug=debug)

        # Check the IRF/Data/Intensities/BaseCalls/SampleSheet.csv file.

        if rta in '1.18.54':
            # Only the MiSeq platform has this file.
            file_name = 'SampleSheet.csv'
            if file_name in base_calls_dict:
                del base_calls_dict[file_name]
            else:
                print "Missing file: '{}".format(file_name)

        # The base_calls_dict should now be empty.

        if len(base_calls_dict):
            print "Data/Intensities/BaseCalls '{}' with number of entries: {:d}".\
                format(base_calls_name, len(base_calls_dict))
            entry_names = base_calls_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    @staticmethod
    def _check_data_intensities_offsets(intensities_dict, intensities_path, debug=0):
        """Check the I{IRF/Data/Intensities/Offsets/} directory.

        @param intensities_dict: Python C{dict} of I{IRF/Data/Intensities/} entries
        @type intensities_dict: dict
        @param intensities_path: I{IRF/Data/Intensities/} file path
        @type intensities_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        directory_name = 'Offsets'
        if directory_name in intensities_dict:
            del intensities_dict[directory_name]
        else:
            print "Missing offsets directory: '{}'".format(directory_name)
            return

        directory_path = os.path.join(intensities_path, directory_name)
        directory_dict = dict(map(lambda x: (x, 1), os.listdir(directory_path)))

        if debug:
            print "Processing '{}'".format(directory_name)

        # Check the IRF/Data/Intensities/Offsets/offsets.txt file.

        file_name = 'offsets.txt'
        if file_name in directory_dict:
            del directory_dict[file_name]
        else:
            print "Missing offsets file '{}'".format(file_name)

        # Check the IRF/Data/Intensities/Offsets/SubTileOffsets.txt file.

        file_name = 'SubTileOffsets.txt'
        if file_name in directory_dict:
            del directory_dict[file_name]
        else:
            print "Missing offsets file '{}".format(file_name)

        # The directory_dict should now be empty.

        if len(directory_dict):
            print "Data/Intensities/Offsets '{}' with number of entries: {:d}".\
                format(directory_name, len(directory_dict))
            entry_names = directory_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_data_intensities(self, data_dict, data_path, debug=0):
        """Check the I{IRF/Data/Intensities/} directory.

        @param data_dict: Python C{dict} of I{IRF/Data/} entries
        @type data_dict: dict
        @param data_path: I{IRF/Data/} path
        @type data_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        intensities_name = 'Intensities'
        if intensities_name in data_dict:
            del data_dict[intensities_name]
        else:
            print "Missing intensities directory '{}'".format(intensities_name)
            return

        intensities_path = os.path.join(data_path, intensities_name)
        intensities_dict = dict(map(lambda x: (x, 1), os.listdir(intensities_path)))

        if debug:
            print "Processing '{}'".format(intensities_name)

        # Build a list of cycle numbers that have no error map such as the last cycle of a read and index cycles.
        no_error_cycles = list()
        cycles = 1
        for ri_read in self.run_information.reads:
            assert isinstance(ri_read, RunInformationRead)
            if ri_read.index:
                no_error_cycles.extend(range(cycles, cycles + ri_read.cycles))
            else:
                no_error_cycles.append(cycles + ri_read.cycles - 1)
            cycles += ri_read.cycles

        if debug:
            print "Cycles without errorMap files: {!r}".format(no_error_cycles)

        # Check the IRF/Data/Intensities/BaseCalls/ directory.

        self._check_data_intensities_base_calls(
            intensities_dict=intensities_dict,
            intensities_path=intensities_path,
            debug=debug)

        if rta in '2.5.2':
            # The HiSeq 3000 platform has:
            # s.locs

            file_name = 's.locs'
            if file_name in intensities_dict:
                del intensities_dict[file_name]
            else:
                print "Missing locations file: '{}'".format(file_name)
        else:
            # The HiSeq 2000 platform has:
            # config.xml
            # L00[1-8]
            # Offsets
            # RTAConfiguration.xml

            # Process the IRF/Data/Intensities/config.xml file.

            file_name = 'config.xml'
            if file_name in intensities_dict:
                del intensities_dict[file_name]
            else:
                print "Missing configuration file: '{}'".format(file_name)

            # Check for completeness of tiles in the image analysis configuration XML file (config.xml).

            self._check_tiles_image_analysis()

            # Process IRF/Data/Intensities/L00[1-8]/ directories.

            for lane in range(1, fcl.lane_count + 1):

                lane_name = 'L{:03d}'.format(lane)
                if lane_name in intensities_dict:
                    del intensities_dict[lane_name]
                else:
                    print "Missing lane directory: 'L{:03d}'".format(lane)
                    continue

                lane_path = os.path.join(intensities_path, lane_name)
                lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            # Not all tiles have to exists especially after catastrophic events during the
                            # cluster generation step.
                            tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if self._is_missing_image_analysis_tile(lane=lane, tile=tile_name):
                                continue
                            if rta in '1.18.54':
                                # The MiSeq platform uses locs files. Sigh.
                                # s_1_1101.locs
                                # s_1_2119.locs
                                locs_name = 's_{}_{:d}{:d}{:02d}.locs'.format(lane, surface, swath, tile)
                                if locs_name in lane_dict:
                                    del lane_dict[locs_name]
                                else:
                                    print "Missing locs file '{}' '{}'".format(lane_name, locs_name)
                            else:
                                # s_1_1101.clocs
                                # s_1_2316.clocs
                                clocs_name = 's_{}_{:d}{:d}{:02d}.clocs'.format(lane, surface, swath, tile)
                                if clocs_name in lane_dict:
                                    del lane_dict[clocs_name]
                                else:
                                    print "Missing clocs file '{}' '{}'".format(lane_name, clocs_name)

                # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/C[0-9]+.1 directories.

                if rta not in ('2.5.2', '1.18.64', '1.18.54'):
                    # Not for HiSeq 3000, HiSeq 2500 and MiSeq platforms.
                    for cycle in range(1, self.run_information.get_cycle_number + 1):

                        cycle_name = 'C{:d}.1'.format(cycle)
                        if cycle_name in lane_dict:
                            del lane_dict[cycle_name]
                        else:
                            print "Missing cycle directory '{}' '{}'".format(lane_name, cycle_name)
                            continue

                        cycle_path = os.path.join(lane_path, cycle_name)
                        cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                        for surface in range(1, fcl.surface_count + 1):
                            for swath in range(1, fcl.swath_count + 1):
                                for tile in range(1, fcl.tile_count + 1):
                                    # Not all tiles have to exists especially after catastrophic events during the
                                    # cluster generation step.
                                    tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                                    if self._is_missing_image_analysis_tile(lane=lane, tile=tile_name):
                                        continue
                                    # s_1_1101.cif
                                    # s_1_2316.cif
                                    cif_name = 's_{}_{:d}{:d}{:02d}.cif'.format(lane, surface, swath, tile)
                                    if cif_name in cycle_dict:
                                        del cycle_dict[cif_name]
                                    else:
                                        print "Missing cif file '{}' '{}' '{}'".\
                                            format(lane_name, cycle_name, cif_name)

                                    if rta in ('1.12.4.2', '1.13.48'):
                                        if cycle not in no_error_cycles:
                                            # Process error map files, that do not exist for index cycles.
                                            # s_1_1101.errorMap
                                            # s_1_2308.errorMap
                                            error_name = 's_{}_{:d}{:d}{:02d}.errorMap'.\
                                                format(lane, surface, swath, tile)
                                            if error_name in cycle_dict:
                                                del cycle_dict[error_name]
                                            else:
                                                print "Missing error map file '{}' '{}' '{}'".\
                                                    format(lane_name, cycle_name, error_name)

                                        # Process FWHM map files.
                                        # s_1_1101_T.FWHMMap
                                        # s_1_2308_T.FWHMMap
                                        fwhm_name = 's_{}_{:d}{:d}{:02d}_T.FWHMMap'.\
                                            format(lane, surface, swath, tile)
                                        if fwhm_name in cycle_dict:
                                            del cycle_dict[fwhm_name]
                                        else:
                                            print "Missing FWHM map file '{}' '{}' '{}'".\
                                                format(lane_name, cycle_name, fwhm_name)

                        if len(cycle_dict):
                            print "Data/Intensities/BaseCalls/Lane '{}' Cycle '{}' with number of entries: {:d}".\
                                format(lane_name, cycle_name, len(cycle_dict))
                            entry_names = cycle_dict.keys()
                            entry_names.sort(cmp=lambda x, y: cmp(x, y))
                            print '  Remaining entries: {!r}'.format(entry_names)

                if len(lane_dict):
                    print "Data/Intensities/Lane '{}' with number of entries: {:d}".format(lane_name, len(lane_dict))
                    entry_names = lane_dict.keys()
                    entry_names.sort(cmp=lambda x, y: cmp(x, y))
                    print '  Remaining entries: {!r}'.format(entry_names)

            # Check the IRF/Data/Intensities/Offsets/ directory.

            self._check_data_intensities_offsets(
                intensities_dict=intensities_dict,
                intensities_path=intensities_path,
                debug=debug)

            # Check the IRF/Data/Intensities/RTAConfiguration.xml file.

            file_name = 'RTAConfiguration.xml'
            if file_name in intensities_dict:
                del intensities_dict[file_name]
            else:
                print "Missing Real Time Analysis configuration file '{}'".format(file_name)

        if rta in ('1.12.4.2', '1.13.48'):
            # RTA 1.12.4.2 has position files in addition to clocs files.
            # s_1_1101_pos.txt
            # s_8_2308_pos.txt
            for lane in range(1, fcl.lane_count + 1):
                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            pos_name = 's_{}_{:d}{:d}{:02d}_pos.txt'.format(lane, surface, swath, tile)
                            if pos_name in intensities_dict:
                                del intensities_dict[pos_name]
                            else:
                                print "Missing pos file 'L{:03d}' '{}'".\
                                    format(lane, pos_name)

        # The intensities_dict should now be empty.

        if len(intensities_dict):
            print "Data/Intensities '{}' with number of entries: {:d}".format(intensities_name, len(intensities_dict))
            entry_names = intensities_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_data_tile_status(self, data_dict, data_path, debug=0):
        """Check the I{IRF/Data/TileStatus/} directory.

        @param data_dict: Python C{dict} of I{IRF/Data/} entries
        @type data_dict: dict
        @param data_path: I{IRF/Data/} path
        @type data_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        fcl = self.run_information.flow_cell_layout

        directory_name = 'TileStatus'
        if directory_name in data_dict:
            del data_dict[directory_name]
        else:
            print "Missing tile status directory '{}'".format(directory_name)
            return

        directory_path = os.path.join(data_path, directory_name)
        directory_dict = dict(map(lambda x: (x, 1), os.listdir(directory_path)))

        if debug:
            print "Processing '{}'".format(directory_name)

        for lane in range(1, fcl.lane_count + 1):
            for surface in range(1, fcl.surface_count + 1):
                for swath in range(1, fcl.swath_count + 1):
                    for tile in range(1, fcl.tile_count + 1):
                        tile_prefix = 'TileStatusL{:d}T{:d}{:d}{:02d}'.format(lane, surface, swath, tile)
                        file_name = tile_prefix + '.bin'
                        if file_name in directory_dict:
                            del directory_dict[file_name]
                        else:
                            print "Missing tile status bin file 'L{}' '{}'".format(lane, file_name)

                        file_name = tile_prefix + '.tpl'
                        if file_name in directory_dict:
                            del directory_dict[file_name]
                        else:
                            print "Missing tile status tpl file 'L{}' '{}'".format(lane, file_name)

        # The directory_dict should now be empty.

        if len(data_dict):
            print "Data/TileStatus '{}' with number of entries: {:d}".format(directory_name, len(directory_dict))
            entry_names = directory_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_data(self, folder_dict, folder_path, debug=0):
        """Check the IRF/Data/ directory.

        @param folder_dict: Ptyhon C{dict} of Illumina Run Folder I{IRF/} entries
        @type folder_dict: dict
        @param folder_path: Illumina Run Folder I{IRF/} path
        @type folder_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        rta = self.run_parameters.get_real_time_analysis_version

        data_name = 'Data'
        if data_name in folder_dict:
            del folder_dict[data_name]
        else:
            print "Missing data directory '{}'".format(data_name)
            return

        data_path = os.path.join(folder_path, data_name)
        data_dict = dict(map(lambda x: (x, 1), os.listdir(data_path)))

        if debug:
            print "Processing '{}'".format(data_name)

        # Check the IRF/Data/Intensities/ directory.

        self._check_data_intensities(
            data_path=data_path,
            data_dict=data_dict,
            debug=debug)

        if rta in '2.5.2':
            # HiSeq 3000 does not have a IRF/Data/ImageSize.dat file.
            pass
        else:
            # Check the IRF/Data/ImageSize.dat file.
            file_name = 'ImageSize.dat'
            if file_name in data_dict:
                del data_dict[file_name]
            else:
                print "Missing image size file '{}".format(file_name)

            # Check the IRF/Data/RTALogs directory.
            directory_name = 'RTALogs'
            if directory_name in data_dict:
                del data_dict[directory_name]
            else:
                print "Missing RTA logs directory '{}'".format(directory_name)

        if rta in ('1.12.4.2', '1.13.48'):
            # Check the IRF/Data/reports/ directory.
            # TODO: Check the directory for completeness?
            directory_name = 'reports'
            if directory_name in data_dict:
                del data_dict[directory_name]
            else:
                print "Missing reports directory: '{}'".format(directory_name)

            # Check the IRF/Data/Status_Files/ directory.
            # TODO: Check the directory for completeness?
            directory_name = 'Status_Files'
            if directory_name in data_dict:
                del data_dict[directory_name]
            else:
                print "Missing status directory: '{}'".format(directory_name)

            # Check the IRF/Data/Status.htm file.
            file_name = 'Status.htm'
            if file_name in data_dict:
                del data_dict[file_name]
            else:
                print "Missing status file: '{}'".format(file_name)

        if rta in '1.18.54':
            # Check the IRF/Data/TileStatus/ directory that only exist on the MiSeq platform.
            self._check_data_tile_status(
                data_path=data_path,
                data_dict=data_dict,
                debug=debug)

        # The data_dict should now be empty.

        if len(data_dict):
            print "Data '{}' with number of entries: {:d}".format(data_name, len(data_dict))
            entry_names = data_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_inter_op(self, folder_dict, folder_path, debug=0):
        """Check the I{IRF/InterOp/} directory.

        @param folder_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type folder_dict: dict
        @param folder_path: Illumina Run Folder I{IRF/} path
        @type folder_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        rta = self.run_parameters.get_real_time_analysis_version

        directory_name = 'InterOp'
        if directory_name in folder_dict:
            del folder_dict[directory_name]
        else:
            print "Missing directory: '{}'".format(directory_name)
            return

        directory_path = os.path.join(folder_path, directory_name)
        directory_dict = dict(map(lambda x: (x, 1), os.listdir(directory_path)))

        file_list = [
            # 'ControlMetricsOut.bin',
            'CorrectedIntMetricsOut.bin',
            'ErrorMetricsOut.bin',
            'ExtractionMetricsOut.bin',
            # 'ImageMetricsOut.bin',
            'QMetricsOut.bin',
            'TileMetricsOut.bin',
        ]

        if rta in '2.5.2':
            # HiSeq 3000 platform
            file_list.append('EmpiricalPhasingMetricsOut.bin')
            file_list.append('EventMetricsOut.bin')
            file_list.append('PFGridMetricsOut.bin')
            file_list.append('RegistrationMetricsOut.bin')
        else:
            # other platforms
            file_list.append('ControlMetricsOut.bin')
            if rta not in '1.18.54':
                file_list.append('ImageMetricsOut.bin')

        if rta in '1.18.54':
            file_list.append('IndexMetricsOut.bin')

        self._check_files(
            directory_dict=directory_dict,
            directory_path=directory_path,
            file_list=file_list,
            debug=debug)

        # The directory_dict should now be empty.

        if len(directory_dict):
            print "InterOp '{}' with number of entries: {:d}".format(directory_name, len(directory_dict))
            entry_names = directory_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_periodic_save_rates(self, folder_dict, folder_path, debug=0):
        """Check the I{IRF/PeriodicSaveRates/} directory.

        @param folder_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type folder_dict: dict
        @param folder_path: Illumina Run Folder I{IRF/} path
        @type folder_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        directory_name = 'PeriodicSaveRates'
        if directory_name in folder_dict:
            del folder_dict[directory_name]
        else:
            print "Missing directory: '{}'".format(directory_name)
            return

        directory_path = os.path.join(folder_path, directory_name)
        directory_dict = dict(map(lambda x: (x, 1), os.listdir(directory_path)))

        file_list = [
            'Save All Thumbnails.xml'
        ]

        self._check_files(
            directory_dict=directory_dict,
            directory_path=directory_path,
            file_list=file_list,
            debug=debug)

        # The directory_dict should now be empty.

        if len(directory_dict):
            print "PeriodicSaveRates '{}' with number of entries: {:d}".format(directory_name, len(directory_dict))
            entry_names = directory_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_recipe(self, folder_dict, folder_path, debug=0):
        """Check the I{IRF/Recipe/} directory.

        @param folder_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type folder_dict: dict
        @param folder_path: Illumina Run Folder I{IRF/} path
        @type folder_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        rta = self.run_parameters.get_real_time_analysis_version
        flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.upper()

        directory_name = 'Recipe'
        if directory_name in folder_dict:
            del folder_dict[directory_name]
        else:
            print "Missing directory: '{}'".format(directory_name)
            return

        directory_path = os.path.join(folder_path, directory_name)
        directory_dict = dict(map(lambda x: (x, 1), os.listdir(directory_path)))

        file_list = list()
        if rta in '1.18.54':
            # The MiSeq platform uses the reagent kit barcode.

            file_list.append(
                self.run_parameters.element_tree.find(path='ReagentKitRFIDTag/SerialNumber').text + '.xml')
            file_list.append('RunState.xml')
        else:
            file_list.append(flow_cell_barcode + '.xml')

            if rta not in '2.5.2':
                # The HiSeq 3000 platform does not have a 'FCID_RunState.xml' file.
                file_list.append(flow_cell_barcode + '_RunState.xml')

        self._check_files(
            directory_dict=directory_dict,
            directory_path=directory_path,
            file_list=file_list,
            debug=debug)

        # The directory_dict should now be empty.

        if len(directory_dict):
            print "Recipe '{}' with number of entries: {:d}".format(directory_name, len(directory_dict))
            entry_names = directory_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def _check_thumbnail_images(self, folder_dict, folder_path, debug=0):
        """Check the I{IRF/Thumbnail_Images/} directory.

        @param folder_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type folder_dict: dict
        @param folder_path: Illumina Run Folder I{IRF/} path
        @type folder_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.lower()

        # Helper dict to map surface numbers to abbreviations.

        surface_dict = dict()
        surface_dict[1] = 'bot'
        surface_dict[2] = 'top'

        # Process the IRF/Thumbnail_Images/ directory.

        thumbnail_name = 'Thumbnail_Images'
        if thumbnail_name in folder_dict:
            del folder_dict[thumbnail_name]
        else:
            print "Missing directory: '{}'".format(thumbnail_name)
            return

        thumbnail_path = os.path.join(folder_path, thumbnail_name)
        thumbnail_dict = dict(map(lambda x: (x, 1), os.listdir(thumbnail_path)))

        if debug:
            print "Processing '{}'".format(thumbnail_name)

        for lane in range(1, fcl.lane_count + 1):

            # Process the IRF/Thumbnail_Images/L00[1-8]/ directories.

            lane_name = 'L{:03d}'.format(lane)
            if lane_name in thumbnail_dict:
                del thumbnail_dict[lane_name]
            else:
                print "Missing lane directory: 'L{:03d}'".format(lane)
                continue

            lane_path = os.path.join(thumbnail_path, lane_name)
            lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

            for cycle in range(1, self.run_information.get_cycle_number + 1):

                # Process the IRF/Thumbnail_Images/L00[1-8]/C[0-9]+.1/ directories.

                cycle_name = 'C{:d}.1'.format(cycle)
                if cycle_name in lane_dict:
                    del lane_dict[cycle_name]
                else:
                    print "Missing cycle directory '{}' '{}'".format(lane_name, cycle_name)
                    continue

                cycle_path = os.path.join(lane_path, cycle_name)
                cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for base in ('a', 'c', 'g', 't'):
                            # Process swath image and zprof files.
                            if rta in '1.18.54':
                                # The MiSeq platform does not have swath image and zprof files.
                                pass
                            else:
                                # c6nk1anxx_c001_l1_t001_bot_s1_a.jpg
                                # c6nk1anxx_c001_l1_t001_bot_s1_a.jpg.zprof
                                # c6nk1anxx_c001_l1_t001_top_s3_t.jpg
                                # c6nk1anxx_c001_l1_t001_top_s3_t.jpg.zprof
                                swath_file = '{}_c{:03d}_l{:d}_t{:03d}_{}_s{}_{}.jpg'.\
                                    format(flow_cell_barcode, cycle, lane, 1, surface_dict[surface], swath, base)
                                if swath_file in cycle_dict:
                                    del cycle_dict[swath_file]
                                else:
                                    print "Missing swath image file '{}' '{}' '{}'".\
                                        format(lane_name, cycle_name, swath_file)

                                swath_file += '.zprof'
                                if rta in '2.5.2' and base in ('c', 'g', 't'):
                                    # The HiSeq 3000 platform does not have swath files for bases c, g and t.
                                    pass
                                else:
                                    if swath_file in cycle_dict:
                                        del cycle_dict[swath_file]
                                    else:
                                        print "Missing swath zprof file '{}' '{}' '{}'".\
                                            format(lane_name, cycle_name, swath_file)

                            # Process tile image files.
                            if rta in ('2.5.2', '1.18.54'):
                                # The HiSeq 3000 and MiSeq platforms use lower case bases.
                                pass
                            else:
                                # The HiSeq 2000 platform uses upper case bases.
                                base = base.upper()
                            for tile in range(1, fcl.tile_count + 1):
                                # s_1_1101_A.jpg
                                # s_1_2316_T.jpg
                                tile_file = 's_{:1d}_{:1d}{:1d}{:02d}_{}.jpg'.format(lane, surface, swath, tile, base)
                                if tile_file in cycle_dict:
                                    del cycle_dict[tile_file]
                                else:
                                    print "Missing tile file '{}' '{}' '{}'".\
                                        format(lane_name, cycle_name, tile_file)

                # The cycle_dict should now be empty.

                if len(cycle_dict):
                    print "Thumbnail_Images/Lane '{}' Cycle '{}' with number of entries: {:d}".\
                        format(lane_name, cycle_name, len(cycle_dict))
                    entry_names = cycle_dict.keys()
                    entry_names.sort(cmp=lambda x, y: cmp(x, y))
                    print '  Remaining entries: {!r}'.format(entry_names)

            # The lane_dict should now be empty.

            if len(lane_dict):
                print "Thumbnail_Images/Lane '{}' with number of entries: {:d}".format(lane_name, len(lane_dict))
                entry_names = lane_dict.keys()
                entry_names.sort(cmp=lambda x, y: cmp(x, y))
                print '  Remaining entries: {!r}'.format(entry_names)

        # The thumbnail_dict should now be empty.

        if len(thumbnail_dict):
            print "Thumbnail_Images/Thumbnail_Images '{}' with number of entries: {:d}".\
                format(thumbnail_name, len(thumbnail_dict))
            entry_names = thumbnail_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)

    def check(self, debug=0):
        """Check an Illumina Run Folder regarding its internal directory and file structure and report both,
        missing and additional files.

        @param debug: Integer debugging level
        @type debug: int
        @return :
        @rtype :
        """

        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in (
            '1.12.4.2',
            '1.13.48',
            '1.17.21.3',
            '1.18.54',  # MiSeq
            '1.18.61',
            '1.18.64',
            '2.5.2',
        ):
            raise Exception("Unsupported RTA version: '{}'".format(rta))

        folder_name = os.path.basename(self.file_path)

        folder_path = self.file_path
        folder_dict = dict(map(lambda x: (x, 1), os.listdir(folder_path)))

        # Check the IRF/Data/ directory.

        self._check_data(
            folder_dict=folder_dict,
            folder_path=folder_path,
            debug=debug)

        # Check the IRF/InterOp/ directory.

        self._check_inter_op(
            folder_dict=folder_dict,
            folder_path=folder_path,
            debug=debug)

        # Check the IRF/PeriodicSaveRates/ directory.

        if rta not in '1.18.54':
            # The MiSeq platform does not have this directory.
            self._check_periodic_save_rates(
                folder_dict=folder_dict,
                folder_path=folder_path,
                debug=debug)

        # Check the IRF/Recipe/ directory.

        self._check_recipe(
            folder_dict=folder_dict,
            folder_path=folder_path,
            debug=debug)

        # Check the IRF/Thumbnail_Images/ directory.

        self._check_thumbnail_images(
            folder_dict=folder_dict,
            folder_path=folder_path,
            debug=debug)

        # Check files.

        file_list = [
            'Logs',  # directory
            'RTAComplete.txt',
            'RunInfo.xml',
            'runParameters.xml',
        ]

        if rta in '1.18.54':
            # The MiSeq platform has a sample annotation sheet.
            file_list.append('SampleSheet.csv')
        else:
            # The MiSeq platform does not have this file.
            file_list.append('First_Base_Report.htm')

        if rta in '2.5.2':
            # HiSeq 3000 platform.
            file_list.append('RTAConfiguration.xml')
            file_list.append('RTALogs')  # directory
            file_list.append('SequencingComplete.txt')
            for read in range(1, len(self.run_information.reads) + 1):
                file_list.append('RTARead{:d}Complete.txt'.format(read))
        else:
            file_list.append('Basecalling_Netcopy_complete.txt')
            file_list.append('ImageAnalysis_Netcopy_complete.txt')
            for read in range(1, len(self.run_information.reads) + 1):
                file_list.append('Basecalling_Netcopy_complete_Read{:d}.txt'.format(read))
                file_list.append('ImageAnalysis_Netcopy_complete_Read{:d}.txt'.format(read))

        self._check_files(directory_dict=folder_dict, directory_path=folder_path, file_list=file_list, debug=debug)

        # The thumbnail_dict should now be empty.

        if len(folder_dict):
            print "Folder '{}' with number of entries: {:d}".format(folder_name, len(folder_dict))
            entry_names = folder_dict.keys()
            entry_names.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries: {!r}'.format(entry_names)
