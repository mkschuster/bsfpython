"""Bio.BSF.Illumina

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
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree


class RunInformationFlowcellLayout(object):
    """The C{RunInformationFlowcellLayout} class models
    one <FlowcellLayout> XML element in an Illumina Run Information (RunInfo.xml) document.
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
        """

        self.lane_count = lane_count
        self.surface_count = surface_count
        self.swath_count = swath_count
        self.tile_count = tile_count


class RunInformationRead(object):
    """The C{RunInformationRead} class models
    one <Read> XML element in an Illumina Run Information (RunInfo.xml) document.
    """

    def __init__(self, number=0, cycles=0, index=False):
        """Initialise a C{RunInformationRead} object.

        @param number: Read number
        @type number: int
        @param cycles: Cycle number
        @type cycles: int
        @param index: Index read
        @type index: bool
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
        """Split an Illumina Run Identifier into its components.

        Splits the Illumina Run Identifier into <Date>_<Instrument>_<Number>_<FCPosition><Flowcell>.
        This method is particularly useful for older version of RunInfo.xml files, that lack
        <Run>/<Date> and <Run>/<Flowcell> elements.
        @param run_identifier: Illumina Run Identifier (e.g. 130724_SN815_0089_BC26JBACXX)
        @type run_identifier: str
        @return: Python C{list} of Python C{str} objects
        @rtype: list
        """

        # Split into <Date>_<Instrument>_<Number>_<FCPosition><Flowcell>
        components = string.split(run_identifier, sep='_')

        if len(components) != 4:
            warnings.warn('Cannot split Illumina Run Identifier {!r} into ist components.'.format(run_identifier))
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
        @return: C{RunInformation}
        @rtype: RunInformation
        """

        file_name = os.path.basename(file_path.rstrip('/ '))

        run_info_tree = ET.parse(source=file_path)
        run_info_root = run_info_tree.getroot()

        # Parse meta-information about the Illumina run

        run_identifier = run_info_root.find('Run').attrib['Id']  # e.g. 130724_SN815_0089_BC26JBACXX
        run_number = run_info_root.find('Run').attrib['Number']  # e.g. 91
        run_identifier_components = RunInformation.parse_run_identifier(run_identifier=run_identifier)

        xml_flow_cell = run_info_root.find('Run/Flowcell')
        if xml_flow_cell is not None:
            flow_cell = xml_flow_cell.text  # e.g. C26JBACXX
        else:
            flow_cell = run_identifier_components[4]

        xml_instrument = run_info_root.find('Run/Instrument')
        if xml_instrument is not None:
            instrument = xml_instrument.text  # e.g. SN815
        else:
            instrument = run_identifier_components[1]

        xml_date = run_info_root.find('Run/Date')
        if xml_date is not None:
            date = xml_date.text  # e.g. 130724
        else:
            date = run_identifier_components[0]

        xml_second_read = run_info_root.find('Run/SecondRead')
        if xml_second_read is not None:
            second_read = xml_second_read.attrib['FirstCycle']
        else:
            second_read = 0

        reads = list()
        number = 1

        for xml_read in run_info_root.find('Run/Reads'):

            # <ApplicationName>HiSeq Control Software</ApplicationName>
            # <ApplicationVersion>2.0.12.0</ApplicationVersion>
            #
            # <Reads>
            #   <Read Number="1" NumCycles="100" IsIndexedRead="N" />
            #   <Read Number="2" NumCycles="8" IsIndexedRead="Y" />
            #   <Read Number="3" NumCycles="100" IsIndexedRead="N" />
            # </Reads>

            if 'NumCycles' in xml_read.attrib:
                is_index = bool(0)

                if xml_read.attrib['IsIndexedRead'] == 'Y':
                    is_index = True
                elif xml_read.attrib['IsIndexedRead'] == 'N':
                    is_index = False
                else:
                    warning = 'Unexpected value {} in Read element attribute IsIndexedRead '. \
                        format(xml_read.attrib['IsIndexedRead'])
                    warnings.warn(warning)

                reads.append(RunInformationRead(
                    number=int(xml_read.attrib['Number']),
                    cycles=int(xml_read.attrib['NumCycles']),
                    index=is_index))

            # <ApplicationName>HiSeq Control Software</ApplicationName>
            # <ApplicationVersion>1.1.37</ApplicationVersion>
            #
            # <Reads>
            #   <Read FirstCycle="1" LastCycle="51" />
            #   <Read FirstCycle="52" LastCycle="102" />
            # </Reads>
            # <SecondRead FirstCycle="52" />

            elif 'FirstCycle' in xml_read.attrib:

                # This is a bit convoluted. If a read after the first has a FirstCycle attribute of less than the
                # FirstCycle attribute of the SecondRead, then it must be the index read. Sigh!

                if number > 1 and xml_read.attrib['FirstCycle'] < second_read:
                    is_index = True
                else:
                    is_index = False

                reads.append(RunInformationRead(
                    number=number,
                    cycles=int(xml_read.attrib['LastCycle']) - int(xml_read.attrib['FirstCycle']),
                    index=is_index))

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

        xml_flow_cell_layout = run_info_root.find('Run/FlowcellLayout')
        if xml_flow_cell_layout is not None:
            flow_cell_layout = RunInformationFlowcellLayout(
                lane_count=int(xml_flow_cell_layout.attrib['LaneCount']),
                surface_count=int(xml_flow_cell_layout.attrib['SurfaceCount']),
                swath_count=int(xml_flow_cell_layout.attrib['SwathCount']),
                tile_count=int(xml_flow_cell_layout.attrib['TileCount']))
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
        @param file_type: File type (e.g. CASAVA, External or Automatic)
        @type file_type: str
        @param name: Name
        @type name: str
        @param run_identifier: Run identifier e.g. 130724_SN815_0089_BC26JBACXX
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
    def get_picard_read_structure(self):
        """Get the read structure for the Picard ExtractIlluminaBarcodes module.

        Codes:
        T ... Template
        B ... Barcode
        S ... Skip
        @return: Read structure for Picard ExtractIlluminaBarcodes
        @rtype: str
        """

        # TODO: This needs to become smarter to deal with skipped bases, caused by
        # different read and barcode lengths. Skips would have to be inserted.

        read_structure = str()

        for read in self.reads:
            read_structure += str(read.cycles)
            if read.index:
                read_structure += 'B'
            else:
                read_structure += 'T'

        return read_structure


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
        @return: C{RunParameters} object
        @rtype: RunParameters
        """

        # file_name = os.path.basename(file_path)

        return RunParameters(file_path=file_path, element_tree=ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{RunParameters} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: XML Element Tree
        @type element_tree: ElementTree
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if element_tree:
            self.element_tree = element_tree
        else:
            self.element_tree = ElementTree()

    @property
    def get_experiment_name(self):
        """Get the experiment name of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<ExperimentName>} value.
        @return: Experiment name e.g BSF_0001
        @rtype: str
        """

        return self.element_tree.getroot().find('Setup/ExperimentName').text

    @property
    def get_flow_cell_barcode(self):
        """Get the flow-cell barcode of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<Barcode>} value.
        @return: Flow-cell barcode e.g BSF_0001
        @rtype: str
        """

        return self.element_tree.getroot().find('Setup/Barcode').text

    @property
    def get_flow_cell_type(self):
        """Get the flow cell of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<Flowcell>} value.
        @return: Flow-cell type
        @rtype: str
        """

        return self.element_tree.getroot().find('Setup/Flowcell').text

    @property
    def get_position(self):
        """Get the flow-cell position of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<FCPosition>} value.
        Since the element does not exist in older versions an empty string may be returned.
        @return: Flow-cell position e.g. A or B
        @rtype: str
        """

        element = self.element_tree.getroot().find('Setup/FCPosition')
        if element is not None:
            return element.text
        else:
            return str()

    @property
    def get_run_identifier(self):
        """Get the run identifier of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<RunID>} value.
        @return: Run identifier
        @rtype: str
        """

        return self.element_tree.getroot().find('Setup/RunID').text

    @property
    def get_read1(self):
        """Get the read 1 cycle number of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<Read1>} value.
        @return: Number of cycles in read 1
        @rtype: str
        """

        return self.element_tree.getroot().find('Setup/Read1').text

    @property
    def get_read2(self):
        """Get the read 2 cycle number of a C{RunParameters} object.

        Get the text representation of the C{<Setup>/<Read2>} value.
        @return: Number of cycles in read 2
        @rtype: str
        """

        return self.element_tree.getroot().find('Setup/Read2').text

    @property
    def get_index_read1(self):
        """Get the index read 1 cycle number of a C{RunParameters} object.

        Normally, this corresponds to the text representation of the C{<IndexRead1>} element,
        while older implementations of Illumina runParameters.xml have only an C{<IndexRead>} element.
        @return: Number of cycles in index read 1
        @rtype: str
        """

        element = self.element_tree.getroot().find('Setup/IndexRead1')
        if element is not None:
            return element.text

        element = self.element_tree.getroot().find('Setup/IndexRead')
        if element is not None:
            return element.text
        else:
            return str()

    @property
    def get_index_read2(self):
        """Get the index read 2 cycle number of a C{RunParameters} object.

        Normally, this corresponds to the text representation of the C{<IndexRead2>} element,
        while older implementations of Illumina runParameters.xml have only an C{<IndexRead>} element.
        @return: Number of cycles in index read 2
        @rtype: str
        """

        element = self.element_tree.getroot().find('Setup/IndexRead2')
        if element is not None:
            return element.text
        else:
            return str()


class RunFolder(object):
    """The C{RunFolder} class represents an Illumina
    Run Folder copied off the instrument.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar file_type: File type
        CASAVA: FASTQ file after post-processing with CASAVA
        External: other data files
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
        @return: C{RunFolder} object
        @rtype: RunFolder
        """

        # Since os.path.basename returns an empty string at trailing slashes,
        # use string.rstrip to remove them.

        file_name = os.path.basename(file_path.rstrip('/ '))

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
            run_parameters=RunParameters.from_file_path(file_path=os.path.join(file_path, 'runParameters.xml')))

        return irf

    def __init__(self, file_path=None, file_type=None, name=None,
                 date=None, instrument=None, run=None, flow_cell=None,
                 run_information=None, run_parameters=None):
        """Initialise a C{RunFolder} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param file_type: File type (e.g. CASAVA, External or Automatic)
        @type file_type: str
        @param name: Name
        @type name: str
        @param date: Date in YYMMDD format
        @type date: str
        @param instrument: Illumina instrument serial number (e.g. SN181, SN815, ...)
        @type instrument: str
        @param run: Run serial number
        @type run: str
        @param flow_cell: The position and flow cell identifier
        @type flow_cell: str
        @param run_information: C{RunInformation}
        @type run_information: RunInformation
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

    @property
    def get_base_calls_directory(self):
        """Get the Base-calls directory in the C{RunFolder/Data/Intensities/BaseCalls} hierarchy.

        @return: Illumina Base-calls Directory
        @rtype: str | unicode
        """

        return os.path.join(self.file_path, 'Data', 'Intensities', 'BaseCalls')
