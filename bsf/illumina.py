"""bsf.illumina

A package of classes and methods modelling data directories and files
specific for Illumina HiSeq and MiSeq systems.
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


import datetime
import os
import warnings
import xml.etree.ElementTree

import dateutil.tz

from bsf.annotation import AnnotationSheet


class Adaptors(object):
    """The C{bsf.illumina.Adaptors} class models Illumina Sequencing Adaptor index sequences.

    Attributes:
    @cvar default_class: Default class
    @type default_class: str
    @cvar default_type: Default type
    @type default_type: str
    @ivar adaptor_dict: Adaptor Python C{dict}
    @type adaptor_dict: dict[str, dict[str, dict[str, str]]]
    """
    default_class = 'Default'
    default_type = 'Default'

    def __init__(self, adaptor_dict=None):
        """Initialise a C{bsf.illumina.Adaptors} object.

        @param adaptor_dict: Adaptor Python C{dict}
        @type adaptor_dict: dict[str, dict[str, dict[str, str]]]
        @return:
        @rtype:
        """
        super(Adaptors, self).__init__()

        if adaptor_dict is None:
            self.adaptor_dict = dict()
        else:
            self.adaptor_dict = adaptor_dict

        return

    @classmethod
    def from_file_path(cls, file_path):
        """Instantiate a C{bsf.illumina.Adaptors} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return: C{bsf.illumina.Adaptors}
        @rtype: bsf.illumina.Adaptors
        """
        annotation_sheet = AnnotationSheet.from_file_path(
            file_path=file_path,
            file_type='excel_tab',
            name='Illumina Adaptors')

        adaptors = cls()

        for row_dict in annotation_sheet.row_dicts:
            class_str = row_dict['Class']
            if not class_str:
                class_str = cls.default_class
            if class_str not in adaptors.adaptor_dict:
                adaptors.adaptor_dict[class_str] = dict()
            class_dict = adaptors.adaptor_dict[class_str]

            type_str = row_dict['Type']
            if not type_str:
                type_str = cls.default_type
            if type_str not in class_dict:
                class_dict[type_str] = dict()
            type_dict = class_dict[type_str]

            name_str = row_dict['Name']
            if name_str not in type_dict:
                type_dict[name_str] = row_dict['Sequence']

        return adaptors

    def match(self, sequence, adaptor_class=None, adaptor_type=None, adaptor_name=None):
        """Match an (unknown) adaptor sequence to an adaptor class, type and name.

        @param sequence: Sequence
        @type sequence: str
        @param adaptor_class: Adaptor class
        @type adaptor_class: str
        @param adaptor_type: Adaptor type
        @type adaptor_type: str
        @param adaptor_name: Adaptor name
        @type adaptor_name: str
        @return: Python C{tuple} of Python C{str} (adaptor class), Python C{str} (adaptor type) and
            Python C{str} (adaptor name)
        @rtype: (str, str, str)
        """
        if adaptor_class is not None and adaptor_class:
            adaptor_class_list = [adaptor_class]
        else:
            adaptor_class_list = self.adaptor_dict.keys()

        # Iterate over all adaptor classes.
        for adaptor_class_key in adaptor_class_list:
            class_dict = self.adaptor_dict[adaptor_class_key]
            if adaptor_type is not None and adaptor_type:
                adaptor_type_list = [adaptor_type]
            else:
                adaptor_type_list = class_dict.keys()

            # Iterate over all adaptor types.
            for adaptor_type_key in adaptor_type_list:
                type_dict = class_dict[adaptor_type_key]
                if adaptor_name is not None and adaptor_name:
                    adaptor_name_list = [adaptor_name]
                else:
                    adaptor_name_list = type_dict.keys()

                # Iterate over all adaptor names.
                for adaptor_name_key in adaptor_name_list:
                    if sequence == type_dict[adaptor_name_key]:
                        return adaptor_class_key, adaptor_type_key, adaptor_name_key


class RunInformationFlowcellLayout(object):
    """The C{bsf.illumina.RunInformationFlowcellLayout} class models
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
        """Initialise a C{bsf.illumina.RunInformationFlowcellLayout} object.

        @param lane_count: Number of lanes
        @type lane_count: int
        @param surface_count: Number of surfaces
        @type surface_count: int
        @param swath_count: Number of swaths
        @type swath_count: int
        @param tile_count: Number of tiles
        @type tile_count: int
        @return:
        @rtype:
        """

        super(RunInformationFlowcellLayout, self).__init__()

        self.lane_count = lane_count
        self.surface_count = surface_count
        self.swath_count = swath_count
        self.tile_count = tile_count

        return


class RunInformationRead(object):
    """The C{bsf.illumina.RunInformationRead} class models
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
        """Initialise a C{bsf.illumina.RunInformationRead} object.

        @param number: Read number
        @type number: int
        @param cycles: Cycle number
        @type cycles: int
        @param index: Index read
        @type index: bool
        @return:
        @rtype:
        """

        super(RunInformationRead, self).__init__()

        self.number = number
        self.cycles = cycles
        self.index = index

        return


class RunInformation(object):
    """The C{bsf.illumina.RunInformation} class represents an Illumina
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
    @ivar flow_cell: Illumina flow cell identifier e.g. C26JBACXX
    @type flow_cell: str
    @ivar instrument: Illumina instrument serial number e.g. SN815
    @type instrument: str
    @ivar date: Date in YYMMDD format e.g. 130724
    @type date: str
    @ivar reads: Python C{list} of C{bsf.illumina.RunInformationRead} objects
    @type reads: list[bsf.illumina.RunInformationRead]
    """

    @staticmethod
    def parse_run_identifier(run_identifier):
        """Split an I{bsf.illumina.Illumina Run Identifier} into its components.

        Splits the Illumina Run Identifier into <Date>_<Instrument>_<Number>_<FCPosition><Flowcell>.
        This method is particularly useful for older version of RunInfo.xml files, that lack
        <Run>/<Date> and <Run>/<Flowcell> elements.
        @param run_identifier: Illumina Run Identifier (e.g. 130724_SN815_0089_BC26JBACXX)
        @type run_identifier: str
        @return: Python C{list} of Python C{str} objects
        @rtype: list[str]
        """

        # Split into <Date>_<Instrument>_<Number>_<FCPosition><Flowcell>
        components = run_identifier.split('_')

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
        """Create a C{bsf.illumina.RunInformation} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return: C{bsf.illumina.RunInformation}
        @rtype: bsf.illumina.RunInformation
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        run_info_tree = xml.etree.ElementTree.ElementTree(file=file_path)
        run_element = run_info_tree.find(path='Run')
        """ @type run_element: Element | None """
        if run_element is None:
            raise Exception('Cannot find the <Run> element in the ElementTree of XML file {!r}.'.format(file_path))

        # Parse meta-information about the Illumina run

        run_identifier = run_element.get(key='Id')  # e.g. 130724_SN815_0089_BC26JBACXX
        """ @type run_identifier: str """
        run_number = run_element.get(key='Number')  # e.g. 91
        """ @type run_number: str """
        run_identifier_components = RunInformation.parse_run_identifier(run_identifier=run_identifier)

        xml_flow_cell = run_element.find(path='Flowcell')
        """ @type xml_flow_cell: Element | None """
        if xml_flow_cell is not None:
            flow_cell = xml_flow_cell.text  # e.g. C26JBACXX
            """ @type flow_cell: str """
        else:
            flow_cell = run_identifier_components[4]

        xml_instrument = run_element.find(path='Instrument')
        """ @type xml_instrument: Element | None """
        if xml_instrument is not None:
            instrument = xml_instrument.text  # e.g. SN815
            """ @type instrument: str """
        else:
            instrument = run_identifier_components[1]

        xml_date = run_element.find(path='Date')
        """ @type xml_date: Element | None """
        if xml_date is not None:
            date = xml_date.text  # e.g. 130724
            """ @type date: str """
        else:
            date = run_identifier_components[0]

        xml_second_read = run_element.find(path='SecondRead')
        """ @type xml_second_read: Element | None """
        if xml_second_read is not None:
            second_read = int(xml_second_read.get(key='FirstCycle'))
        else:
            second_read = int(0)

        reads = list()
        """ @type reads: list[bsf.illumina.RunInformationRead] """
        number = int(1)

        for read_element in run_element.find(path='Reads'):
            """ @type read_element: Element | None """
            assert isinstance(read_element, xml.etree.ElementTree.Element)

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
                """ @type is_index: str """
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

        # Get the flow cell layout if it exits.

        xml_flow_cell_layout = run_info_tree.find(path='Run/FlowcellLayout')
        """ @type xml_flow_cell_layout: Element | None """
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

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            run_identifier=None,
            run_number=None,
            flow_cell=None,
            instrument=None,
            date=None,
            reads=None,
            flow_cell_layout=None):
        """Initialise a C{bsf.illumina.RunInformation} object.

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
        @param flow_cell: Illumina flow cell identifier
        @type flow_cell: str
        @param instrument: Illumina instrument serial number
        @type instrument: str
        @param date: Date in YYMMDD format
        @type date: str
        @param reads: Python C{list} of C{bsf.illumina.RunInformationRead} objects
        @type reads: list[bsf.illumina.RunInformationRead]
        @param flow_cell_layout: C{bsf.illumina.RunInformationFlowcellLayout}
        @type flow_cell_layout: bsf.illumina.RunInformationFlowcellLayout
        @return:
        @rtype:
        """

        super(RunInformation, self).__init__()

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

        return

    @property
    def get_cycle_number(self):
        """Get the total number of cycles.

        @return: Number of cycles
        @rtype: int
        """

        cycle_number = 0

        for read in self.reads:
            cycle_number += read.cycles

        return cycle_number

    @property
    def get_iso_date(self):
        """Get the run start date in ISO 8601 format.

        @return: The run start date in ISO 8601 format.
        @rtype: str
        """

        return datetime.datetime(
            int(self.date[0:2]) + 2000,  # year
            int(self.date[2:4]),  # month
            int(self.date[4:6]),  # day
            0,  # hour
            0,  # minute
            0,  # second
            0,  # microsecond
            dateutil.tz.tzlocal()  # tzinfo
        ).isoformat()

    @property
    def get_read_number(self):
        """Get the total number of reads.

        @return: Number of reads
        @rtype: int
        """

        return len(self.reads)

    @property
    def get_read_start_list(self):
        """Get a Python C{list} of cycle numbers at the start of each read.

        @return: Python C{list} of Python C{int} (starting cycle) for each read
        @rtype: list[int]
        """

        cycle_number = 0
        """ @type cycle_number: int """
        read_start_list = list()
        """ @type read_start_list: list[int] """

        for read in self.reads:
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
        @return: Read structure for Picard C{ExtractIlluminaBarcodes}
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

    @property
    def get_read_structure_list(self):
        """Get the read structure Python C{list}.

        Codes:
        B ... Base
        I ... Index
        @return: Python C{list} of Python C{str} (read structure) objects
        @rtype: list[str]
        """

        read_structure_list = list()
        """ @type read_structure_list: list[str] """

        for read in self.reads:
            if read.index:
                read_structure_list.append(str(read.cycles) + 'I')
            else:
                read_structure_list.append(str(read.cycles) + 'B')

        return read_structure_list


class RunParameters(object):
    """The C{bsf.illumina.RunParameters} class models the contents of runParameters.xml
    files inside an Illumina Run Folder.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar element_tree: C{xml.etree.ElementTree.ElementTree}
    @type element_tree: xml.etree.ElementTree.ElementTree
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{bsf.illumina.RunParameters} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return: C{bsf.illumina.RunParameters} object
        @rtype: bsf.illumina.RunParameters
        """

        file_path = os.path.normpath(file_path)
        return cls(file_path=file_path, element_tree=xml.etree.ElementTree.ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{bsf.illumina.RunParameters} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: C{xml.etree.ElementTree.ElementTree}
        @type element_tree: xml.etree.ElementTree.ElementTree
        @return:
        @rtype:
        """

        super(RunParameters, self).__init__()

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if element_tree:
            self.element_tree = element_tree
        else:
            self.element_tree = xml.etree.ElementTree.ElementTree()

        self._run_parameters_version = None
        """ @type _run_parameters_version: str | None """

        return

    @property
    def get_run_parameters_version(self):
        """Get the run parameters version of a C{bsf.illumina.RunParameters} object.

        Returns an empty string for I{HiSeq Control Software} or the text representation of the
        I{<RunParameters>/<RunParametersVersion>} element for I{MiSeq Control Software}.
        @return: Run parameters version or an empty string
        @rtype: str
        """

        if self._run_parameters_version is None:
            element = self.element_tree.find(path='RunParametersVersion')
            if element is None:
                self._run_parameters_version = str()
            else:
                self._run_parameters_version = element.text

        return self._run_parameters_version

    @property
    def get_instrument_type(self):
        """Get the instrument type based on a C{bsf.illumina.RunParameters} object.

        Returns I{MiSeq}, I{NextSeq} or I{HiSeq}
        @return: Instrument type
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1',):
            return 'MiSeq'
        elif self.get_run_parameters_version in ('NextSeq_2_1_0',):
            return 'NextSeq'
        else:
            return 'HiSeq'

    @property
    def get_experiment_name(self):
        """Get the experiment name of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<ExperimentName>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<ExperimentName>} element for
        I{MiSeq Control Software}.
        @return: Experiment name
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq or NextSeq
            return self.element_tree.find(path='ExperimentName').text
        else:
            # HiSeq or else
            return self.element_tree.find(path='Setup/ExperimentName').text

    @property
    def get_flow_cell_barcode(self):
        """Get the flow cell barcode of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Barcode>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<Barcode>} element for
        I{MiSeq Control Software}.
        @return: Flow cell barcode
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1',):
            # MiSeq 1_1
            return self.element_tree.find(path='Barcode').text
        elif self.get_run_parameters_version in ('NextSeq_2_1_0',):
            # NextSeq
            return self.element_tree.find(path='FlowCellSerial').text
        else:
            # HiSeq or else
            return self.element_tree.find(path='Setup/Barcode').text

    @property
    def get_flow_cell_type(self):
        """Get the flow cell chemistry type of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Flowcell>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<ReagentKitVersion>} element for
        I{MiSeq Control Software}.
        @return: Flow cell chemistry type
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1',):
            # MiSeq 1_1
            # The MiSeq has no concept of a <Flowcell> chemistry type, only a <ReagentKitVersion>.
            return self.element_tree.find(path='ReagentKitVersion').text
        elif self.get_run_parameters_version in ('NextSeq_2_1_0',):
            # NextSeq_2_1_0
            # The NextSeq uses <Chemistry>
            return self.element_tree.find(path='Chemistry').text
        else:
            # HiSeq or else
            return self.element_tree.find(path='Setup/Flowcell').text

    @property
    def get_index_type(self):
        """Get the index chemistry type of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Index>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        @return: Index chemistry type
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            # The MiSeq has no concept of a <Index> version, only a <ReagentKitVersion>.
            return str()
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/Index')
            if element is None:
                return str()
            else:
                return element.text

    @property
    def get_pe_type(self):
        """Get the paired-end chemistry type of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Pe>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        @return: Paired-end chemistry type
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            # The MiSeq has no concept of a <Pe> version, only a <ReagentKitVersion>.
            return str()
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/Pe')
            if element is None:
                return str()
            else:
                return element.text

    @property
    def get_sbs_type(self):
        """Get the sequencing-by-synthesis chemistry type of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Sbs>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        @return: Sequencing-by-synthesis chemistry type
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            # The MiSeq has no concept of a <Sbs> chemistry version, only a <ReagentKitVersion>.
            return str()
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/Sbs')
            if element is None:
                return str()
            else:
                return element.text

    @property
    def get_position(self):
        """Get the flow cell position of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<FCPosition>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Since the element does not exist in older I{HiSeq Control Software} versions an empty string may be returned.
        @return: Flow cell position e.g. A or B
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1',):
            # MiSeq 1_1
            # The MiSeq has no concept of a <FCPosition>.
            return str()
        elif self.get_run_parameters_version in ('NextSeq_2_1_0',):
            # The NextSeq has no concept of <FCPosition>, but always uses 'A'.
            return 'A'
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/FCPosition')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_run_identifier(self):
        """Get the run identifier of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<RunID>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<RunID>} element for
        I{MiSeq Control Software}.
        @return: Run identifier
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            return self.element_tree.find(path='RunID').text
        else:
            # HiSeq or else
            return self.element_tree.find(path='Setup/RunID').text

    @property
    def get_read1(self):
        """Get the read 1 cycle number of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Read1>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        @return: Number of cycles in read 1
        @rtype: str
        @deprecated: The more scalable option is getting read information via the
            C{bsf.illumina.RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            return str()
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/Read1')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_read2(self):
        """Get the read 2 cycle number of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<Read2>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Not every run has a read 2 defined so that an empty string may be returned.
        @return: Number of cycles in read 2
        @rtype: str
        @deprecated: The more scalable option is getting read information via the
            C{bsf.illumina.RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            return str()
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/Read2')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_index_read1(self):
        """Get the index read 1 cycle number of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<IndexRead1>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Older implementations of the I{HiSeq Control Software} have only a
        I{<RunParameters>/<Setup>/<IndexRead>} element.
        @return: Number of cycles in index read 1
        @rtype: str
        @deprecated: The more scalable option is getting read information via the
            C{bsf.illumina.RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            return str()
        else:
            # HiSeq or else
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
        """Get the index read 2 cycle number of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<IndexRead2>} element for
        I{HiSeq Control Software} or an empty string for
        I{MiSeq Control Software}.
        Older implementations of the I{HiSeq Control Software} have only a
        I{<RunParameters>/<Setup>/<IndexRead>} element.
        @return: Number of cycles in index read 2
        @rtype: str
        @deprecated: The more scalable option is getting read information via the
            C{bsf.illumina.RunInformation.reads} instance variable.
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            return str()
        else:
            # HiSeq or else
            element = self.element_tree.find(path='Setup/IndexRead2')
            if element is not None:
                return element.text
            else:
                return str()

    @property
    def get_real_time_analysis_version(self):
        """Get the Real-Time Analysis (RTA) Version of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of the I{<RunParameters>/<Setup>/<RTAVersion>} element for
        I{HiSeq Control Software} or the I{<RunParameters>/<RTAVersion>} element for
        I{MiSeq Control Software}.
        @return: RTA version
        @rtype: str
        """

        if self.get_run_parameters_version in ('MiSeq_1_1', 'NextSeq_2_1_0'):
            # MiSeq 1_1
            return self.element_tree.find(path='RTAVersion').text
        else:
            # HiSeq or else
            return self.element_tree.find(path='Setup/RTAVersion').text

    @property
    def get_application_name(self):
        """Get the application (i.e. I{HiSeq} or I{MiSeq Control Software}) name.

        Returns the text representation of the I{<RunParameters>/<Setup>/<ApplicationName>} element.
        @return: Application name
        @rtype: str
        """

        return self.element_tree.find(path='Setup/ApplicationName').text

    @property
    def get_application_version(self):
        """Get the application (i.e. I{HiSeq} or I{MiSeq Control Software}) version.

        Returns the text representation of the I{<RunParameters>/<Setup>/<ApplicationVersion>} element.
        @return: Application version
        @rtype: str
        """

        return self.element_tree.find(path='Setup/ApplicationVersion').text


class XMLConfiguration(object):
    """The C{bsf.illumina.XMLConfiguration} class models the contents of XML configuration
    files inside an Illumina Run Folder.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar element_tree: C{xml.etree.ElementTree.ElementTree}
    @type element_tree: xml.etree.ElementTree.ElementTree
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{bsf.illumina.XMLConfiguration} object from a file path.
        In case the file path does not exist a C{bsf.illumina.XMLConfiguration} object
        with an empty C{xml.etree.ElementTree.ElementTree} will be returned.

        @param file_path: File path
        @type file_path: str | unicode
        @return: C{bsf.illumina.XMLConfiguration} object
        @rtype: bsf.illumina.XMLConfiguration
        """
        file_path = os.path.normpath(file_path)
        if not os.path.isfile(file_path):
            file_path = None
        return cls(file_path=file_path, element_tree=xml.etree.ElementTree.ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{bsf.illumina.XMLConfiguration} object.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: C{xml.etree.ElementTree.ElementTree}
        @type element_tree: xml.etree.ElementTree.ElementTree
        @return:
        @rtype:
        """

        super(XMLConfiguration, self).__init__()

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if element_tree:
            self.element_tree = element_tree
        else:
            self.element_tree = xml.etree.ElementTree.ElementTree()

        return


class AnalysisConfiguration(XMLConfiguration):
    """The C{bsf.illumina.AnalysisConfiguration} class models Image and Base Call analysis
    XML configuration files inside and Illumina Run Folder.

    Attributes:
    @ivar file_path: File path
    @type file_path: str | unicode
    @ivar element_tree: C{xml.etree.ElementTree.ElementTree}
    @type element_tree: xml.etree.ElementTree.ElementTree
    """

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{bsf.illumina.AnalysisConfiguration} object.
        Read all I{<Tile>} elements from the I{<Run>/<TileSelection>} element from the XML configuration file
        to initialise an internal Python C{dict} of valid tiles.

        @param file_path: File path
        @type file_path: str | unicode
        @param element_tree: C{xml.etree.ElementTree.ElementTree}
        @type element_tree: xml.etree.ElementTree.ElementTree
        @return:
        @rtype:
        """

        super(AnalysisConfiguration, self).__init__(file_path=file_path, element_tree=element_tree)

        self._lane_tile_dict = dict()
        """ @type _lane_tile_dict: dict[str, dict[str, bool]] """

        if self.element_tree.getroot() is None:
            return

        for lane_element in self.element_tree.find(path='Run/TileSelection'):
            """ @type lane_element: Element """
            assert isinstance(lane_element, xml.etree.ElementTree.Element)
            lane_index = lane_element.get(key='Index')
            """ @type lane_index: str """
            if lane_index not in self._lane_tile_dict:
                self._lane_tile_dict[lane_index] = dict()
            lane_dict = self._lane_tile_dict[lane_index]
            for tile_element in lane_element.findall(path='Tile'):
                """ @type title_element: Element """
                assert isinstance(tile_element, xml.etree.ElementTree.Element)
                lane_dict[tile_element.text] = True

        return

    def has_lane(self, lane):
        """Check if a particular lane is defined in a C{bsf.illumina.AnalysisConfiguration} object.

        @param lane: Lane index
        @type lane: str
        @return: Boolean value
        @rtype: bool
        """
        return lane in self._lane_tile_dict

    def has_lane_tile(self, lane, tile):
        """Check if a particular tile is defined in a lane of a C{bsf.illumina.AnalysisConfiguration} object.

        @param lane: Lane index
        @type lane: str
        @param tile: Tile index
        @type tile: str
        @return: Boolean value
        @rtype: bool | None
        """
        if lane in self._lane_tile_dict:
            lane_dict = self._lane_tile_dict[lane]
            return tile in lane_dict
        else:
            return None


class ImageAnalysis(AnalysisConfiguration):
    """The C{bsf.illumina.ImageAnalysis} class models the contents of the
    I{IRF/Data/Intensities/config.xml} XML configuration file inside an Illumina Run Folder.

    Attributes:
    """

    pass


class BaseCallAnalysis(AnalysisConfiguration):
    """The C{bsf.illumina.BaseCallAnalysis} class models the contents of the
    I{IRF/Data/Intensities/BaseCalls/config.xml} XML configuration file inside an Illumina Run Folder.

    Attributes:
    """

    pass


class RunFolderNotComplete(Exception):
    pass


class RunFolder(object):
    """The C{bsf.illumina.RunFolder} class represents an Illumina
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
    @ivar flow_cell: Flow cell identifier
    @type flow_cell: str
    @ivar run_information: C{bsf.illumina.RunInformation} object
    @type run_information: bsf.illumina.RunInformation
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Construct a C{bsf.illumina.RunFolder} object from a file path.

        @param file_path: File path
        @type file_path: str | unicode
        @return: C{bsf.illumina.RunFolder}
        @rtype: bsf.illumina.RunFolder
        """

        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        # HiSeq and MiSeq Control Software uses runParameters.xml, while NextSeq Control Software uses a
        # RunParameters.xml file. Why?

        if os.path.exists(os.path.join(file_path, 'runParameters.xml')):
            run_parameters_path = os.path.join(file_path, 'runParameters.xml')
        elif os.path.exists(os.path.join(file_path, 'RunParameters.xml')):
            run_parameters_path = os.path.join(file_path, 'RunParameters.xml')
        else:
            raise Exception('Found neither a runParameters.xml not a RunParameters.xml file.')

        # Illumina Run Folders obey a 'YYMMDD_SN000_Run_PositionFCID' schema.

        components = file_name.split('_')

        irf = cls(
            file_path=file_path,
            file_type='Illumina',
            name='_'.join(components),
            date=components[0],
            instrument=components[1],
            run=components[2],
            flow_cell=components[3],
            run_information=RunInformation.from_file_path(file_path=os.path.join(file_path, 'RunInfo.xml')),
            run_parameters=RunParameters.from_file_path(file_path=run_parameters_path),
            image_analysis=ImageAnalysis.from_file_path(file_path=os.path.join(
                file_path, 'Data', 'Intensities', 'config.xml')),
            base_call_analysis=BaseCallAnalysis.from_file_path(file_path=os.path.join(
                file_path, 'Data', 'Intensities', 'BaseCalls', 'config.xml')))

        return irf

    def __init__(
            self,
            file_path=None,
            file_type=None,
            name=None,
            date=None,
            instrument=None,
            run=None,
            flow_cell=None,
            run_information=None,
            run_parameters=None,
            image_analysis=None,
            base_call_analysis=None):
        """Initialise a C{bsf.illumina.RunFolder} object.

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
        @param run_information: C{bsf.illumina.RunInformation}
        @type run_information: bsf.illumina.RunInformation
        @param image_analysis: C{bsf.illumina.ImageAnalysis}
        @type image_analysis: bsf.illumina.ImageAnalysis
        @param base_call_analysis: C{bsf.illumina.BaseCallAnalysis}
        @type base_call_analysis: bsf.illumina.BaseCallAnalysis
        @return:
        @rtype:
        """

        super(RunFolder, self).__init__()

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
        """ @type _missing_base_call_tiles: dict[int, dict[str, bool]] """

        self._missing_image_analysis_tiles = dict()
        """ @type _missing_image_analysis_tiles: dict[int, dict[str, bool]] """

        return

    @property
    def get_base_calls_directory(self):
        """Get the base-calls directory in the I{IRF/Data/Intensities/BaseCalls} hierarchy.

        @return: Illumina base-calls directory
        @rtype: str | unicode
        """

        return os.path.join(self.file_path, 'Data', 'Intensities', 'BaseCalls')

    def _check_tiles_base_call(self):
        """Check for missing I{<Tile>} elements in the I{IRF/Data/Intensities/BaseCalls/config.xml}
        configuration file. This method also builds up a Python C{dict} required for method
        C{_is_missing_base_call_tile}.

        @return:
        @rtype:
        """
        fcl = self.run_information.flow_cell_layout

        for lane in range(1, fcl.lane_count + 1):
            if lane not in self._missing_base_call_tiles:
                self._missing_base_call_tiles[lane] = dict()
            lane_dict = self._missing_base_call_tiles[lane]
            if self.base_call_analysis.has_lane(lane=str(lane)):
                # Lanes are not defined of the config.xml file could not be read.
                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            tile = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if not self.base_call_analysis.has_lane_tile(lane=str(lane), tile=tile):
                                lane_dict[tile] = True
                                print 'Missing BaseCallAnalysis lane:', lane, 'tile:', tile

        return

    def _check_tiles_image_analysis(self):
        """Check for missing I{<Tile>} elements in the I{IRF/Data/Intensities/config.xml}
        configuration file. This method also builds up a Python C{dict} required for method
        C{_is_missing_image_analysis_tile}.

        @return:
        @rtype:
        """
        fcl = self.run_information.flow_cell_layout

        for lane in range(1, fcl.lane_count + 1):
            if lane not in self._missing_image_analysis_tiles:
                self._missing_image_analysis_tiles[lane] = dict()
            lane_dict = self._missing_image_analysis_tiles[lane]
            if self.image_analysis.has_lane(lane=str(lane)):
                # Lanes are not defined of the config.xml file could not be read.
                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            tile = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if not self.image_analysis.has_lane_tile(lane=str(lane), tile=tile):
                                lane_dict[tile] = True
                                print 'Missing ImageAnalysis lane:', lane, 'tile:', tile

        return

    def _is_missing_base_call_tile(self, lane, tile):
        """Confirm that a particular I{<Tile>} element is missing from the
        I{IRF/Data/Intensities/config.xml} configuration file.

        @param lane: Lane index
        @type lane: int
        @param tile: Tile name
        @type tile: str
        @return: Boolean value
        @rtype: bool
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
        @return: Boolean value
        @rtype: bool
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
    def _check_files(directory_dict, directory_path, file_name_list, debug=0):
        """Check a Python C{list} of file names against a Python C{dict} of directory entries.

        @param directory_dict: Python C{dict} of directory entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: Directory path
        @type directory_path: str | unicode
        @param file_name_list: Python C{list} of file names
        @type file_name_list: list[str | unicode]
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """
        if debug > 0:
            # print 'Processing directory', directory_path
            pass

        for file_name in file_name_list:
            file_path = os.path.join(directory_path, file_name)
            if file_name in directory_dict:
                del directory_dict[file_name]
            else:
                print 'Missing file', file_path

        return

    def _check_data_intensities_base_calls_matrix(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/Matrix/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/Intensities/BaseCalls/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: I{IRF/Data/Intensities/BaseCalls/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        if rta in ('2.4.11', '2.7.3', '2.7.6', '2.7.7'):
            # RTA 2.4.11 (NextSeq) doe not have a IRF/Data/Intensities/BaseCalls/Matrix/ directory.
            # RTA 2.7.3 (HiSeq 3000/4000) does no longer have a IRF/Data/Intensities/BaseCalls/Matrix/ directory.
            return

        _directory_name = 'Matrix'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        if rta in ('2.5.2',):
            # For RTA 2.5.2 (HiSeq 3000/4000) process IRF/Data/Intensities/BaseCalls/Matrix/L00[1-8] directories.
            for lane in range(1, fcl.lane_count + 1):
                lane_name = 'L{:03d}'.format(lane)
                lane_path = os.path.join(_directory_path, lane_name)
                if lane_name in _directory_dict:
                    del _directory_dict[lane_name]
                else:
                    print 'Missing directory', lane_path
                    continue
                lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))
                """ @type lane_dict: dict[str | unicode, int] """

                if debug > 1:
                    print 'Processing directory', lane_path

                # Process IRF/Data/Intensities/BaseCalls/Matrix/L00[1-8]/C[0-9]+.1/ directories.

                for cycle in range(1, self.run_information.get_cycle_number + 1):
                    cycle_name = 'C{:d}.1'.format(cycle)
                    cycle_path = os.path.join(lane_path, cycle_name)
                    if cycle_name in lane_dict:
                        del lane_dict[cycle_name]
                    else:
                        print 'Missing directory', cycle_path
                        continue
                    cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))
                    """ @type cycle_dict: dict[str | unicode, int] """

                    if debug > 2:
                        print 'Processing directory', cycle_path

                    for surface in range(1, fcl.surface_count + 1):
                        for swath in range(1, fcl.swath_count + 1):
                            for tile in range(1, fcl.tile_count + 1):
                                # Process tile matrix files.
                                # s_1_1101_matrix.txt
                                # s_1_2228_matrix.txt
                                _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}_matrix.txt'.format(lane, surface, swath, tile)
                                if _entry_name in cycle_dict:
                                    del cycle_dict[_entry_name]
                                else:
                                    print 'Missing tile matrix file', os.path.join(cycle_path, _entry_name)

                    if len(cycle_dict):
                        print cycle_path, 'with number of entries:', str(len(cycle_dict))
                        entry_name_list = cycle_dict.keys()
                        entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                        print '  Remaining entries:', entry_name_list

                if len(lane_dict):
                    print lane_path, 'with number of entries:', str(len(lane_dict))
                    entry_name_list = lane_dict.keys()
                    entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                    print '  Remaining entries:', entry_name_list
        else:
            # All other platforms have a flat list of matrix.txt files.
            for read in range(1, self.run_information.get_read_number + 1):
                # Process read matrix files.
                # s_1_matrix.txt
                # s_2_matrix.txt
                _entry_name = 's_{:d}_matrix.txt'.format(read)
                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print 'Missing read matrix file', os.path.join(_directory_path, _entry_name)

                for lane in range(1, fcl.lane_count + 1):
                    # Process lane matrix files.
                    # s_1_1_matrix.txt
                    # s_8_2_matrix.txt
                    _entry_name = 's_{:d}_{:d}_matrix.txt'.format(lane, read)
                    if _entry_name in _directory_dict:
                        del _directory_dict[_entry_name]
                    else:
                        print 'Missing lane matrix file', os.path.join(_directory_path, _entry_name)

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
                                _entry_name = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_matrix.txt'. \
                                    format(lane, read, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print 'Missing tile matrix file', os.path.join(_directory_path, _entry_name)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_data_intensities_base_calls_phasing(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/Phasing/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/intensities/BaseCalls/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: I{IRF/Data/intensities/BaseCalls/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        if rta in ('2.4.11', '2.7.3', '2.7.6', '2.7.7'):
            # RTA 2.4.11 (NextSeq) doe not have a IRF/Data/Intensities/BaseCalls/Matrix/ directory.
            # RTA 2.7.3 (HiSeq 3000/4000) does no longer have a IRF/Data/Intensities/BaseCalls/Phasing/ directory.
            return

        _directory_name = 'Phasing'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', os.path.join(directory_path, _directory_name)
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        for ri_read in self.run_information.reads:
            # Process read phasing files.
            # s_1_phasing.txt
            # s_2_phasing.txt
            _entry_name = 's_{:d}_phasing.txt'.format(ri_read.number)
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing file', os.path.join(_directory_path, _entry_name)

            for lane in range(1, fcl.lane_count + 1):
                # Process lane phasing files.
                # s_1_1_phasing.txt
                # s_8_3_phasing.txt
                _entry_name = 's_{:d}_{:d}_phasing.txt'.format(lane, ri_read.number)
                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print 'Missing file', os.path.join(_directory_path, _entry_name)
                    continue

                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for tile in range(1, fcl.tile_count + 1):
                            # Not all tiles have to exists especially after catastrophic events during the
                            # cluster generation step.
                            tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                                continue
                            if rta in ('1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3') and not ri_read.index:
                                # Process tile cycle files, which only exist for payload, but not index reads.
                                # s_1_1_1101_cycle.txt
                                # s_8_1_2316_cycle.txt
                                _entry_name = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_cycle.txt'. \
                                    format(lane, ri_read.number, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print 'Missing file', os.path.join(_directory_path, _entry_name)

                            # Process tile phasing files.
                            # s_1_1_1101_phasing.txt
                            # s_8_3_2316_phasing.txt
                            _entry_name = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_phasing.txt'. \
                                format(lane, ri_read.number, surface, swath, tile)
                            if _entry_name in _directory_dict:
                                del _directory_dict[_entry_name]
                            else:
                                print 'Missing file', os.path.join(_directory_path, _entry_name)

                            if rta not in ('1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3'):
                                # Process the tile empirical phasing files.
                                _entry_name = 'EmpiricalPhasingCorrection_{:1d}_{:1d}_{:1d}{:1d}{:02d}.txt'. \
                                    format(lane, ri_read.number, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print 'Missing file', os.path.join(_directory_path, _entry_name)

        if rta not in ('2.5.2',):
            # RTA 2.5.2 (HiSeq 3000/4000) does not have
            # IRF/DataIntensities/BaseCalls/Phasing/s_{lane}_{cycle}_phasing.xml files.
            for lane in range(1, fcl.lane_count + 1):
                for read_start in self.run_information.get_read_start_list:
                    _entry_name = 's_{:1d}_{:02d}_phasing.xml'.format(lane, read_start + 2)
                    if _entry_name in _directory_dict:
                        del _directory_dict[_entry_name]
                    else:
                        print 'Missing file', os.path.join(_directory_path, _entry_name)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_data_intensities_base_calls(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/Intensities/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: I{IRF/Data/Intensities/} file path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'BaseCalls'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', os.path.join(directory_path, _directory_name)
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        # Process the IRF/Data/Intensities/BaseCalls/config.xml file.

        if rta not in ('2.4.11', '2.5.2', '2.7.3', '2.7.6', '2.7.7'):
            # HiSeq 3000/4000 and NextSeq does not have the IRF/Data/Intensities/BaseCalls/config.xml file.
            _entry_name = 'config.xml'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing file', os.path.join(_directory_path, _entry_name)

            # Check for completeness of tiles in the base call configuration XML file (config.xml).
            self._check_tiles_base_call()

        # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/ directories.

        for lane in range(1, fcl.lane_count + 1):
            lane_name = 'L{:03d}'.format(lane)
            lane_path = os.path.join(_directory_path, lane_name)
            if lane_name in _directory_dict:
                del _directory_dict[lane_name]
            else:
                print 'Missing directory', lane_path
                continue
            lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))
            """ @type lane_dict: dict[str | unicode, int] """

            if debug > 1:
                print 'Processing directory', lane_path

            # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/C[0-9]+.1 directories.

            if rta in ('2.4.11',):
                # For NextSeq
                for cycle in range(1, self.run_information.get_cycle_number + 1):
                    _entry_name = '{:04d}.bcl.bgzf'.format(cycle)
                    if _entry_name in lane_dict:
                        del lane_dict[_entry_name]
                    else:
                        print 'Missing cycle bcl.bgzf file', os.path.join(lane_path, _entry_name)
                    _entry_name = '{:04d}.bcl.bgzf.bci'.format(cycle)
                    if _entry_name in lane_dict:
                        del lane_dict[_entry_name]
                    else:
                        print 'Missing cycle bcl.bgzf.bci file', os.path.join(lane_path, _entry_name)

                _entry_name = 's_{:d}.bci'.format(lane)
                if _entry_name in lane_dict:
                    del lane_dict[_entry_name]
                else:
                    print 'Missing lane bci file', os.path.join(lane_path, _entry_name)

                _entry_name = 's_{:d}.filter'.format(lane)
                if _entry_name in lane_dict:
                    del lane_dict[_entry_name]
                else:
                    print 'Missing lane filter file', os.path.join(lane_path, _entry_name)
            else:
                # For HiSeq and MiSeq
                for cycle in range(1, self.run_information.get_cycle_number + 1):
                    cycle_name = 'C{:d}.1'.format(cycle)
                    cycle_path = os.path.join(lane_path, cycle_name)
                    if cycle_name in lane_dict:
                        del lane_dict[cycle_name]
                    else:
                        print 'Missing directory', cycle_path
                        continue
                    cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))
                    """ @type cycle_dict: dict[str | unicode, int] """

                    if debug > 2:
                        print 'Processing directory', cycle_path

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
                                _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.bcl'.format(lane, surface, swath, tile)
                                if _entry_name in cycle_dict:
                                    del cycle_dict[_entry_name]
                                else:
                                    # Base call (BCL) files can also be GNU Zip compressed.
                                    # s_1_1101.bcl.gz
                                    # s_1_2316.bcl.gz
                                    _entry_name += '.gz'
                                    if _entry_name in cycle_dict:
                                        del cycle_dict[_entry_name]
                                    else:
                                        print 'Missing tile bcl file', os.path.join(cycle_path, _entry_name)

                                # Process tile stats files.
                                # s_1_1101.stats
                                # s_1_2316.stats
                                if rta not in ('2.5.2', '2.7.3', '2.7.6', '2.7.7'):
                                    # HiSeq 3000/4000 does not have stats files.
                                    _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.stats'.format(lane, surface, swath, tile)
                                    if _entry_name in cycle_dict:
                                        del cycle_dict[_entry_name]
                                    else:
                                        print 'Missing tile stats file', os.path.join(cycle_path, _entry_name)

                    if len(cycle_dict):
                        print cycle_path, 'with number of entries:', str(len(cycle_dict))
                        entry_name_list = cycle_dict.keys()
                        entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                        print '  Remaining entries:', entry_name_list

            # Process control and filter files.

            if rta not in ('2.4.11',):
                # Not for NextSeq
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
                            if rta not in ('2.5.2', '2.7.3', '2.7.6', '2.7.7'):
                                # HiSeq 3000/4000 does not have control files.
                                _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.control'.format(lane, surface, swath, tile)
                                if _entry_name in lane_dict:
                                    del lane_dict[_entry_name]
                                else:
                                    print 'Missing tile control file', os.path.join(lane_path, _entry_name)

                            # Process tile filter files.
                            # s_1_1101.filter
                            # s_1_2316.filter
                            _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.filter'.format(lane, surface, swath, tile)
                            if _entry_name in lane_dict:
                                del lane_dict[_entry_name]
                            else:
                                print 'Missing tile filter file', os.path.join(lane_path, _entry_name)

            if len(lane_dict):
                print lane_path, 'with number of entries:', str(len(lane_dict))
                entry_name_list = lane_dict.keys()
                entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                print '  Remaining entries:', entry_name_list

        # Process the IRF/Data/Intensities/BaseCalls/Matrix/ directory.

        self._check_data_intensities_base_calls_matrix(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        # Process the IRF/Data/Intensities/BaseCalls/Phasing/ directory.

        self._check_data_intensities_base_calls_phasing(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        # Check the IRF/Data/Intensities/BaseCalls/SampleSheet.csv file.

        if rta in ('1.18.54',):
            # Only the MiSeq platform has this file.
            _entry_name = 'SampleSheet.csv'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing file', os.path.join(_directory_path, _entry_name)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_data_intensities_offsets(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/Offsets/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/Intensities/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: I{IRF/Data/Intensities/} file path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        rta = self.run_parameters.get_real_time_analysis_version

        if rta in ('2.4.11',):
            return

        _directory_name = 'Offsets'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        # Check the IRF/Data/Intensities/Offsets/offsets.txt file.

        _entry_name = 'offsets.txt'
        if _entry_name in _directory_dict:
            del _directory_dict[_entry_name]
        else:
            print 'Missing file', os.path.join(_directory_path, _entry_name)

        # Check the IRF/Data/Intensities/Offsets/SubTileOffsets.txt file.

        _entry_name = 'SubTileOffsets.txt'
        if _entry_name in _directory_dict:
            del _directory_dict[_entry_name]
        else:
            print 'Missing file', os.path.join(_directory_path, _entry_name)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_data_intensities(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: I{IRF/Data/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'Intensities'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        # Build a list of cycle numbers that have no error map such as the last cycle of a read and index cycles.
        no_error_cycles = list()
        """ @type no_error_cycles: list[int] """

        cycles = 1
        for ri_read in self.run_information.reads:
            if ri_read.index:
                no_error_cycles.extend(range(cycles, cycles + ri_read.cycles))
            else:
                no_error_cycles.append(cycles + ri_read.cycles - 1)
            cycles += ri_read.cycles

        if debug > 0:
            print 'Cycles without errorMap files:', no_error_cycles

        # Check the IRF/Data/Intensities/BaseCalls/ directory.

        self._check_data_intensities_base_calls(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        if rta in ('2.5.2', '2.7.3', '2.7.6', '2.7.7'):
            # The HiSeq 3000/4000 platform has:
            # s.locs

            # Process the IRF/Data/Intensities/s.locs file.

            _entry_name = 's.locs'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing file', os.path.join(_directory_path, _entry_name)
        else:
            # The HiSeq 2000 platform has:
            # config.xml
            # L00[1-8]
            # Offsets
            # RTAConfiguration.xml

            # Process the IRF/Data/Intensities/config.xml file.

            if rta not in ('2.4.11',):
                # Not for NextSeq
                _entry_name = 'config.xml'
                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print 'Missing configuration file', os.path.join(_directory_path, _entry_name)

            # Check for completeness of tiles in the image analysis configuration XML file (config.xml).

            self._check_tiles_image_analysis()

            # Process IRF/Data/Intensities/L00[1-8]/ directories.

            for lane in range(1, fcl.lane_count + 1):
                lane_name = 'L{:03d}'.format(lane)
                lane_path = os.path.join(_directory_path, lane_name)
                if lane_name in _directory_dict:
                    del _directory_dict[lane_name]
                else:
                    print 'Missing directory', lane_path
                    continue
                lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))
                """ @type lane_dict: dict[str | unicode, int] """

                if debug > 1:
                    print 'Processing directory', lane_path

                if rta in ('2.4.11',):
                    _entry_name = 's_{:d}.locs'.format(lane)
                    if _entry_name in lane_dict:
                        del lane_dict[_entry_name]
                    else:
                        print 'Missing file', os.path.join(lane_path, _entry_name)
                else:
                    for surface in range(1, fcl.surface_count + 1):
                        for swath in range(1, fcl.swath_count + 1):
                            for tile in range(1, fcl.tile_count + 1):
                                # Not all tiles have to exists especially after catastrophic events during the
                                # cluster generation step.
                                tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                                if self._is_missing_image_analysis_tile(lane=lane, tile=tile_name):
                                    continue
                                if rta in ('1.18.54',):
                                    # The MiSeq platform uses locs files. Sigh.
                                    # s_1_1101.locs
                                    # s_1_2119.locs
                                    _entry_name = 's_{}_{:d}{:d}{:02d}.locs'.format(lane, surface, swath, tile)
                                    if _entry_name in lane_dict:
                                        del lane_dict[_entry_name]
                                    else:
                                        print 'Missing file', os.path.join(lane_path, _entry_name)
                                else:
                                    # s_1_1101.clocs
                                    # s_1_2316.clocs
                                    _entry_name = 's_{}_{:d}{:d}{:02d}.clocs'.format(lane, surface, swath, tile)
                                    if _entry_name in lane_dict:
                                        del lane_dict[_entry_name]
                                    else:
                                        print 'Missing file', os.path.join(lane_path, _entry_name)

                # Process IRF/Data/Intensities/L00[1-8]/C[0-9]+.1 directories.

                if rta not in ('1.18.54', '1.18.64', '2.4.11'):
                    # Exclude the MiSeq, NextSeq and HiSeq 2500 platforms,
                    # as they do no longer store cycle-specific sub directories with
                    # cluster intensity files (*.cif), error map (*.errorMap) and
                    # full width at half maximum (*.FWHMMap) files.
                    for cycle in range(1, self.run_information.get_cycle_number + 1):
                        cycle_name = 'C{:d}.1'.format(cycle)
                        cycle_path = os.path.join(lane_path, cycle_name)
                        if cycle_name in lane_dict:
                            del lane_dict[cycle_name]
                        else:
                            print 'Missing directory', cycle_path
                            continue
                        cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))
                        """ @type cycle_dict: dict[str | unicode, int] """

                        if debug > 2:
                            print 'Processing directory', cycle_path

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
                                    _entry_name = 's_{}_{:d}{:d}{:02d}.cif'.format(lane, surface, swath, tile)
                                    if _entry_name in cycle_dict:
                                        del cycle_dict[_entry_name]
                                    else:
                                        print 'Missing cif file', os.path.join(cycle_path, _entry_name)

                                    if rta in ('1.12.4', '1.12.4.2', '1.13.48'):
                                        # Older RTA versions store error map (*.errorMap) and
                                        # full width at half maximum map (*.FWHMMap) files.
                                        if cycle not in no_error_cycles:
                                            # Process *.errorMap files, that do not exist for index cycles.
                                            # s_1_1101.errorMap
                                            # s_1_2308.errorMap
                                            _entry_name = 's_{}_{:d}{:d}{:02d}.errorMap'. \
                                                format(lane, surface, swath, tile)
                                            if _entry_name in cycle_dict:
                                                del cycle_dict[_entry_name]
                                            else:
                                                print 'Missing error map file', os.path.join(cycle_path, _entry_name)

                                        # Process *.FWHMMap files.
                                        # s_1_1101_T.FWHMMap
                                        # s_1_2308_T.FWHMMap
                                        _entry_name = 's_{}_{:d}{:d}{:02d}_T.FWHMMap'. \
                                            format(lane, surface, swath, tile)
                                        if _entry_name in cycle_dict:
                                            del cycle_dict[_entry_name]
                                        else:
                                            print 'Missing FWHM map file', os.path.join(cycle_path, _entry_name)

                        if len(cycle_dict):
                            print cycle_path, 'with number of entries:', str(len(cycle_dict))
                            entry_name_list = cycle_dict.keys()
                            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                            print '  Remaining entries:', entry_name_list

                if len(lane_dict):
                    print lane_path, 'with number of entries:', str(len(lane_dict))
                    entry_name_list = lane_dict.keys()
                    entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                    print '  Remaining entries:', entry_name_list

            # Check the IRF/Data/Intensities/Offsets/ directory.

            self._check_data_intensities_offsets(
                directory_dict=_directory_dict,
                directory_path=_directory_path,
                debug=debug)

            # Check the IRF/Data/Intensities/RTAConfiguration.xml file.

            if rta not in ('2.4.11',):
                # Not for NextSeq
                _entry_name = 'RTAConfiguration.xml'
                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print 'Missing Real Time Analysis configuration file', os.path.join(_directory_path, _entry_name)

            if rta in ('1.12.4', '1.12.4.2', '1.13.48'):
                # Older RTA version have position (*_pos.txt) files in addition to cluster location (*.clocs) files.
                # s_1_1101_pos.txt
                # s_8_2308_pos.txt
                for lane in range(1, fcl.lane_count + 1):
                    for surface in range(1, fcl.surface_count + 1):
                        for swath in range(1, fcl.swath_count + 1):
                            for tile in range(1, fcl.tile_count + 1):
                                _entry_name = 's_{}_{:d}{:d}{:02d}_pos.txt'.format(lane, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print 'Missing pos file', os.path.join(_directory_path, _entry_name)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_data_tile_status(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/TileStatus/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: I{IRF/Data/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in ('1.18.54',):
            # Check the IRF/Data/TileStatus/ directory that only exist on the MiSeq platform.
            return

        _directory_name = 'TileStatus'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        for lane in range(1, fcl.lane_count + 1):
            for surface in range(1, fcl.surface_count + 1):
                for swath in range(1, fcl.swath_count + 1):
                    for tile in range(1, fcl.tile_count + 1):
                        tile_prefix = 'TileStatusL{:d}T{:d}{:d}{:02d}'.format(lane, surface, swath, tile)
                        _entry_name = tile_prefix + '.bin'
                        if _entry_name in _directory_dict:
                            del _directory_dict[_entry_name]
                        else:
                            print 'Missing file', os.path.join(_directory_path, _entry_name)

                        _entry_name = tile_prefix + '.tpl'
                        if _entry_name in _directory_dict:
                            del _directory_dict[_entry_name]
                        else:
                            print 'Missing file', os.path.join(_directory_path, _entry_name)

        if len(directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_data(self, directory_dict, directory_path, debug=0):
        """Check the IRF/Data/ directory.

        @param directory_dict: Ptyhon C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'Data'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        # Check the IRF/Data/Intensities/ directory.

        self._check_data_intensities(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        if rta not in ('2.4.11', '2.5.2', '2.7.3', '2.7.6', '2.7.7'):
            # Exclude the HiSeq 3000/4000 and NextSeq platforms.
            # Check the IRF/Data/ImageSize.dat file.
            _entry_name = 'ImageSize.dat'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing file', os.path.join(_directory_path, _entry_name)

            # Check the IRF/Data/RTALogs directory.
            _entry_name = 'RTALogs'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing directory', os.path.join(_directory_path, _entry_name)

        if rta in ('1.12.4', '1.12.4.2', '1.13.48'):
            # Check the IRF/Data/reports/ directory.
            # TODO: Check the directory for completeness?
            _entry_name = 'reports'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing directory', os.path.join(_directory_path, _entry_name)

            # Check the IRF/Data/Status_Files/ directory.
            # TODO: Check the directory for completeness?
            _entry_name = 'Status_Files'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing directory', os.path.join(_directory_path, _entry_name)

            # Check the IRF/Data/Status.htm file.
            _entry_name = 'Status.htm'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print 'Missing file', os.path.join(_directory_path, _entry_name)

        # Process the IRF/Data/TileStatus/ directory.

        self._check_data_tile_status(
            directory_path=_directory_path,
            directory_dict=_directory_dict,
            debug=debug)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_inter_op(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/InterOp/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'InterOp'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        _file_name_list = [
            'CorrectedIntMetricsOut.bin',
            'ErrorMetricsOut.bin',
            'ExtractionMetricsOut.bin',
            'QMetricsOut.bin',
            'TileMetricsOut.bin',
        ]

        if rta in ('1.18.54',):
            # Only on the MiSeq platform.
            _file_name_list.append('IndexMetricsOut.bin')

        if rta in ('2.4.11', '2.5.2', '2.7.3', '2.7.6', '2.7.7'):
            # HiSeq 3000/4000 and NextSeq platforms
            _file_name_list.append('EmpiricalPhasingMetricsOut.bin')
            _file_name_list.append('EventMetricsOut.bin')
            _file_name_list.append('PFGridMetricsOut.bin')
            _file_name_list.append('RegistrationMetricsOut.bin')

        if rta in ('2.7.3', '2.7.6', '2.7.7'):
            # HiSeq 3000/4000 platform, excluding RTA 2.5.2 version
            _file_name_list.append('ColorMatrixMetricsOut.bin')
            _file_name_list.append('FWHMGridMetricsOut.bin')
            _file_name_list.append('StaticRunMetricsOut.bin')

        if rta not in ('1.18.54', '2.4.11', '2.5.2'):
            _file_name_list.append('ImageMetricsOut.bin')

        if rta not in ('2.4.11', '2.5.2', '2.7.3', '2.7.6', '2.7.7'):
            # Other than HiSeq 3000/4000 and NextSeq platforms
            _file_name_list.append('ControlMetricsOut.bin')

        self._check_files(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_periodic_save_rates(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/PeriodicSaveRates/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        rta = self.run_parameters.get_real_time_analysis_version

        if rta in ('2.4.11',):
            return

        _directory_name = 'PeriodicSaveRates'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        _file_name_list = [
            'Save All Thumbnails.xml'
        ]

        self._check_files(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_recipe(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Recipe/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        rta = self.run_parameters.get_real_time_analysis_version
        flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.upper()

        _directory_name = 'Recipe'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        _file_name_list = list()
        """ @type _file_name_list: list[str | unicode] """

        if rta in ('1.18.54',):
            # The MiSeq platform uses the reagent kit barcode.
            _file_name_list.append(
                self.run_parameters.element_tree.find(path='ReagentKitRFIDTag/SerialNumber').text + '.xml')
            _file_name_list.append('RunState.xml')
        else:
            _file_name_list.append(flow_cell_barcode + '.xml')

            if rta not in ('2.5.2', '2.7.3', '2.7.6', '2.7.7'):
                # The HiSeq 3000/4000 platform does not have a 'FCID_RunState.xml' file.
                _file_name_list.append(flow_cell_barcode + '_RunState.xml')

        self._check_files(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def _check_thumbnail_images(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Thumbnail_Images/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str | unicode, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str | unicode
        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        if rta in ('2.4.11',):
            # Not on the NextSeq platform
            return

        flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.lower()

        # Helper dict to map surface numbers to abbreviations.

        surface_dict = dict()
        """ @type surface_dict: dict[int, str] """
        surface_dict[1] = 'bot'
        surface_dict[2] = 'top'

        # Process the IRF/Thumbnail_Images/ directory.

        _directory_name = 'Thumbnail_Images'
        _directory_path = os.path.join(directory_path, _directory_name)
        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print 'Missing directory', _directory_path
            return
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        # Process the IRF/Thumbnail_Images/L00[1-8]/ directories.

        for lane in range(1, fcl.lane_count + 1):
            lane_name = 'L{:03d}'.format(lane)
            if lane_name in _directory_dict:
                del _directory_dict[lane_name]
            else:
                print 'Missing directory', os.path.join(_directory_path, lane_name)
                continue
            lane_path = os.path.join(_directory_path, lane_name)
            lane_dict = dict(map(lambda x: (x, 1), os.listdir(lane_path)))
            """ @type lane_dict: dict[str | unicode, int] """

            if debug > 1:
                print 'Processing directory', lane_path

            # Process the IRF/Thumbnail_Images/L00[1-8]/C[0-9]+.1/ directories.

            for cycle in range(1, self.run_information.get_cycle_number + 1):
                cycle_name = 'C{:d}.1'.format(cycle)
                cycle_path = os.path.join(lane_path, cycle_name)
                if cycle_name in lane_dict:
                    del lane_dict[cycle_name]
                else:
                    print 'Missing directory', cycle_path
                    continue
                cycle_dict = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))
                """ @type cycle_dict: dict[str | unicode, int] """

                if debug > 2:
                    print 'Processing directory', cycle_path

                for surface in range(1, fcl.surface_count + 1):
                    for swath in range(1, fcl.swath_count + 1):
                        for base in ('a', 'c', 'g', 't'):
                            # Process swath image and zprof files.
                            if rta in ('1.18.54',):
                                # The MiSeq platform does not have swath image and zprof files.
                                pass
                            else:
                                # c6nk1anxx_c001_l1_t001_bot_s1_a.jpg
                                # c6nk1anxx_c001_l1_t001_bot_s1_a.jpg.zprof
                                # c6nk1anxx_c001_l1_t001_top_s3_t.jpg
                                # c6nk1anxx_c001_l1_t001_top_s3_t.jpg.zprof
                                _entry_name = '{}_c{:03d}_l{:d}_t{:03d}_{}_s{}_{}.jpg'. \
                                    format(flow_cell_barcode, cycle, lane, 1, surface_dict[surface], swath, base)
                                if _entry_name in cycle_dict:
                                    del cycle_dict[_entry_name]
                                else:
                                    print 'Missing swath image file', os.path.join(cycle_path, _entry_name)

                                _entry_name += '.zprof'
                                if rta in ('2.5.2', '2.7.3', '2.7.6', '2.7.7') and base in ('c', 'g', 't'):
                                    # The HiSeq 3000/4000 platform does not have swath files for bases c, g and t.
                                    pass
                                else:
                                    if _entry_name in cycle_dict:
                                        del cycle_dict[_entry_name]
                                    else:
                                        print 'Missing swath zprof file', os.path.join(cycle_path, _entry_name)

                            # Process tile image files.
                            if rta in ('1.18.54', '2.5.2', '2.7.3', '2.7.6', '2.7.7'):
                                # The HiSeq 3000/4000 and MiSeq platforms use lower case bases.
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
                                    print 'Missing tile file', os.path.join(cycle_path, tile_file)

                if len(cycle_dict):
                    print cycle_path, 'with number of entries:', str(len(cycle_dict))
                    entry_name_list = cycle_dict.keys()
                    entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                    print '  Remaining entries:', entry_name_list

            if len(lane_dict):
                print lane_path, 'with number of entries:', str(len(lane_dict))
                entry_name_list = lane_dict.keys()
                entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
                print '  Remaining entries:', entry_name_list

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return

    def check(self, debug=0):
        """Check an Illumina Run Folder regarding its internal directory and file structure and report both,
        missing and additional files.

        @param debug: Integer debugging level
        @type debug: int
        @return:
        @rtype:
        """

        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in (
                '1.12.4',  # HiSeq Control Software 1.4.5
                '1.12.4.2',  # HiSeq Control Software 1.4.8
                '1.13.48',  # HiSeq Control Software 1.5.15.1
                '1.17.21.3',  # HiSeq Control Software 2.0.12.0
                '1.18.54',  # MiSeq Control Software 2.5.0.5
                '1.18.61',  # HiSeq Control Software 2.2.38
                '1.18.64',  # HiSeq Control Software 2.2.58
                '2.4.11',  # NextSeq Control Software 2.1.0.31
                '2.5.2',  # HiSeq Control Software 3.3.20 (HiSeq 3000/4000)
                '2.7.3',  # HiSeq Control Software 3.3.52 (HiSeq 3000/4000)
                '2.7.6',  # HiSeq Control Software 3.3.76 (HiSeq 3000/4000)
                '2.7.7',  # HiSeq Control Software HD 3.4.0.38 (HiSeq 3000/4000)
        ):
            raise Exception("Unsupported RTA version: '{}'".format(rta))

        # _directory_name = os.path.basename(self.file_path)
        _directory_path = self.file_path
        _directory_dict = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))
        """ @type _directory_dict: dict[str | unicode, int] """

        if debug > 0:
            print 'Processing directory', _directory_path

        # Check the IRF/Data/ directory.

        self._check_data(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        # Check the IRF/InterOp/ directory.

        self._check_inter_op(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        # Check the IRF/PeriodicSaveRates/ directory.

        if rta not in ('1.18.54',):
            # The MiSeq platform does not have this directory.
            self._check_periodic_save_rates(
                directory_dict=_directory_dict,
                directory_path=_directory_path,
                debug=debug)

        # Check the IRF/Recipe/ directory.

        self._check_recipe(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        # Check the IRF/Thumbnail_Images/ directory.

        self._check_thumbnail_images(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        # Check files.

        _file_name_list = [
            'Logs',  # directory
            'RTAComplete.txt',
            'RunInfo.xml',
        ]

        if rta in ('2.4.11',):
            # On the NextSeq platform
            _file_name_list.append('RunParameters.xml')
        else:
            # On all but the NextSeq platform
            _file_name_list.append('runParameters.xml')

        if rta in ('1.18.54',):
            # The MiSeq platform has a sample annotation sheet.
            _file_name_list.append('SampleSheet.csv')

        if rta not in ('1.18.54', '2.4.11'):
            # Not for MiSeq and NextSeq platforms.
            _file_name_list.append('First_Base_Report.htm')

        if rta in ('2.4.11', '2.5.2', '2.7.3', '2.7.6', '2.7.7'):
            # On HiSeq 3000/4000 and NextSeq platforms.
            _file_name_list.append('RTAConfiguration.xml')
            _file_name_list.append('RTALogs')  # directory
            for read in range(1, len(self.run_information.reads) + 1):
                _file_name_list.append('RTARead{:d}Complete.txt'.format(read))
            if rta not in ('2.4.11',):
                # Not on the NextSeq platform
                _file_name_list.append('SequencingComplete.txt')
        else:
            # Other than HiSeq 3000/4000 platforms.
            _file_name_list.append('Basecalling_Netcopy_complete.txt')
            _file_name_list.append('ImageAnalysis_Netcopy_complete.txt')
            for read in range(1, len(self.run_information.reads) + 1):
                _file_name_list.append('Basecalling_Netcopy_complete_Read{:d}.txt'.format(read))
                _file_name_list.append('ImageAnalysis_Netcopy_complete_Read{:d}.txt'.format(read))

        self._check_files(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print _directory_path, 'with number of entries:', str(len(_directory_dict))
            entry_name_list = _directory_dict.keys()
            entry_name_list.sort(cmp=lambda x, y: cmp(x, y))
            print '  Remaining entries:', entry_name_list

        return
