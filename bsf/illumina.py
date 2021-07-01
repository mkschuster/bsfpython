# -*- coding: utf-8 -*-
"""Illumina module.

A package of classes and methods modelling data directories and files
specific for Illumina sequencing systems.
"""
#  Copyright 2013 - 2021 Michael K. Schuster
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
import datetime
import os
import warnings
from typing import Dict, List, Optional, Tuple
from xml.etree.ElementTree import ElementTree

import dateutil.tz
from Bio.Seq import Seq

from bsf.annotation import AnnotationSheet


class Adaptors(object):
    """The C{bsf.illumina.Adaptors} class models Illumina Sequencing Adaptor index sequences.

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

        @param adaptor_dict: Hierarchy of Python C{dict} objects for class, type (e.g. i7, i5),
            name (e.g. D701, D503) and sequence
        @type adaptor_dict: dict[str, dict[str, dict[str, str]]]
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
        @type file_path: str
        @return: C{bsf.illumina.Adaptors}
        @rtype: Adaptors
        """
        annotation_sheet = AnnotationSheet.from_file_path(
            file_path=file_path,
            file_type='excel-tab',
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
        @type adaptor_class: str | None
        @param adaptor_type: Adaptor type
        @type adaptor_type: str | None
        @param adaptor_name: Adaptor name
        @type adaptor_name: str | None
        @return: Python C{tuple} of Python C{str} (adaptor class), Python C{str} (adaptor type),
            Python C{str} (adaptor name) and Python C{bool} (reverse complement)
        @rtype: list[(str, str, str, bool)]
        """
        result_list: List[Tuple[str, str, str, bool]] = list()

        if adaptor_class:  # not None and not empty
            adaptor_class_list = [adaptor_class]
        else:
            adaptor_class_list = [key for key in self.adaptor_dict]

        # Iterate over all adaptor classes.
        for adaptor_class_key in adaptor_class_list:
            class_dict = self.adaptor_dict[adaptor_class_key]
            if adaptor_type:  # not None and not empty
                adaptor_type_list = [adaptor_type]
            else:
                adaptor_type_list = [key for key in class_dict]

            # Iterate over all adaptor types.
            for adaptor_type_key in adaptor_type_list:
                type_dict = class_dict[adaptor_type_key]
                if adaptor_name:  # not None and not empty
                    adaptor_name_list = [adaptor_name]
                else:
                    adaptor_name_list = [key for key in type_dict]

                # Iterate over all adaptor names.
                for adaptor_name_key in adaptor_name_list:
                    bio_seq = Seq(data=type_dict[adaptor_name_key])

                    # Match in forward orientation
                    if sequence == str(bio_seq):
                        result_list.append((adaptor_class_key, adaptor_type_key, adaptor_name_key, False))

                    # Match in reverse orientation
                    if sequence == str(bio_seq.reverse_complement()):
                        result_list.append((adaptor_class_key, adaptor_type_key, adaptor_name_key, True))

        return result_list


class RunInformationFlowcellLayout(object):
    """The C{bsf.illumina.RunInformationFlowcellLayout} class models
    one I{<FlowcellLayout>} XML element in an I{Illumina Run Information} (RunInfo.xml) document.

    @ivar lane_count: Number of lanes
    @type lane_count: int
    @ivar surface_count: Number of surfaces
    @type surface_count: int
    @ivar swath_count: Number of swaths
    @type swath_count: int
    @ivar tile_count: Number of tiles
    @type tile_count: int
    @ivar tile_list: Python C{list} of Python C{str} tile names
    @type tile_list: list[str]
    """

    def __init__(self, lane_count=0, surface_count=0, swath_count=0, tile_count=0, tile_list=None):
        """Initialise a C{bsf.illumina.RunInformationFlowcellLayout} object.

        @param lane_count: Number of lanes
        @type lane_count: int
        @param surface_count: Number of surfaces
        @type surface_count: int
        @param swath_count: Number of swaths
        @type swath_count: int
        @param tile_count: Number of tiles
        @type tile_count: int
        @param tile_list: Python C{list} of Python C{str} tile names
        @type tile_list: list[str]
        """
        super(RunInformationFlowcellLayout, self).__init__()

        self.lane_count = lane_count
        self.surface_count = surface_count
        self.swath_count = swath_count
        self.tile_count = tile_count

        if tile_list is None:
            self.tile_list = list()
        else:
            self.tile_list = tile_list

        return

    def has_tile(self, tile):
        """Test the C{RunInformationFlowcellLayout} for a particular tile.

        @param tile: Flowcell tile
        @type tile: str
        @return: C{True} if the tile is in the C{RunInformationFlowcellLayout}, C{False} otherwise
        @rtype: bool
        """
        if self.tile_list:
            return tile in self.tile_list
        else:
            # The tile list is not populated, return true.
            return True


class RunInformationRead(object):
    """The C{bsf.illumina.RunInformationRead} class models
    one <Read> XML element in an Illumina Run Information (RunInfo.xml) document.

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
        """
        super(RunInformationRead, self).__init__()

        self.number = number
        self.cycles = cycles
        self.index = index

        return


class RunInformation(object):
    """The C{bsf.illumina.RunInformation} class represents an Illumina
    run information XML (RunInfo.xml) document.

    @ivar file_path: File path
    @type file_path: str | None
    @ivar file_type: File type
        CASAVA: FASTQ file after post-processing with CASAVA.
        External: other data files.
    @type file_type: str | None
    @ivar name: Name
    @type name: str | None
    @ivar run_identifier: Run identifier e.g. 130724_SN815_0089_BC26JBACXX
    @type run_identifier: str | None
    @ivar run_number: Run number, which may not have to correspond to the run number in the run identifier e.g. 91
    @type run_number: str | None
    @ivar flow_cell: Illumina flow cell identifier e.g. C26JBACXX
    @type flow_cell: str | None
    @ivar instrument: Illumina instrument serial number e.g. SN815
    @type instrument: str | None
    @ivar date: Date in YYMMDD format e.g. 130724
    @type date: str | None
    @ivar flow_cell_layout: C{bsf.illumina.RunInformationFlowcellLayout}
    @type flow_cell_layout: RunInformationFlowcellLayout | None
    @ivar run_information_read_list: Python C{list} of C{bsf.illumina.RunInformationRead} objects
    @type run_information_read_list: list[RunInformationRead]
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
        component_list = run_identifier.split('_')

        if len(component_list) != 4:
            warnings.warn('Cannot split Illumina Run Identifier ' + repr(run_identifier) + ' into its components.')
            return

        # Strip leading zeros from the <Number>, split <FCPosition> and >Flowcell> elements.
        component_list[2] = component_list[2].lstrip('0')
        component_list.append(component_list[3][1:])
        component_list[3] = component_list[3][:1]

        return component_list

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{bsf.illumina.RunInformation} object from a file path.

        @param file_path: File path
        @type file_path: str
        @return: C{bsf.illumina.RunInformation}
        @rtype: RunInformation
        """
        file_path = os.path.normpath(file_path)
        file_name = os.path.basename(file_path)

        run_info_tree = ElementTree(file=file_path)
        run_element = run_info_tree.find(path='Run')
        if run_element is None:
            raise Exception('Cannot find the <Run> Element in the ElementTree of XML file ' + repr(file_path))

        # Parse meta-information about the Illumina run

        run_identifier = run_element.get(key='Id')  # e.g. 130724_SN815_0089_BC26JBACXX
        run_number = run_element.get(key='Number')  # e.g. 91
        run_identifier_component_list = RunInformation.parse_run_identifier(run_identifier=run_identifier)

        xml_flow_cell = run_element.find(path='Flowcell')
        if xml_flow_cell is None:
            flow_cell = run_identifier_component_list[4]
        else:
            flow_cell = xml_flow_cell.text  # e.g. C26JBACXX

        xml_instrument = run_element.find(path='Instrument')
        if xml_instrument is None:
            instrument = run_identifier_component_list[1]
        else:
            instrument = xml_instrument.text  # e.g. SN815

        xml_date = run_element.find(path='Date')
        if xml_date is None:
            date = run_identifier_component_list[0]
        else:
            date = xml_date.text  # e.g. 130724

        xml_second_read = run_element.find(path='SecondRead')
        if xml_second_read is None:
            second_read = 0
        else:
            second_read = int(xml_second_read.get(key='FirstCycle'))

        run_information_read_list: List[RunInformationRead] = list()
        number = 1

        for read_element in run_element.find(path='Reads').iter(tag='Read'):
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
                if is_index not in ('Y', 'N'):
                    warnings.warn(
                        'Unexpected value <Read IsIndexedRead="' +
                        read_element.get(key='IsIndexedRead') +
                        '"> in Read element attribute IsIndexedRead ',
                        UserWarning)

                run_information_read_list.append(RunInformationRead(
                    number=int(read_element.get(key='Number')),
                    cycles=int(read_element.get(key='NumCycles')),
                    index=is_index == 'Y'))

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

                run_information_read_list.append(RunInformationRead(
                    number=number,
                    cycles=int(read_element.get(key='LastCycle')) - int(read_element.get(key='FirstCycle')),
                    index=bool(number > 1 and int(read_element.get(key='FirstCycle')) < second_read)))

                number += 1

        # Sort by the RunInformationRead.number.

        run_information_read_list.sort(key=lambda _item: _item.number)

        # Warn if there is not at least one non-index read (i.e. RunInformationRead.index).

        non_index_reads_list = [item for item in run_information_read_list if not item.index]
        if len(non_index_reads_list) == 0:
            warnings.warn(
                'No non-index read in Illumina RunInformation: ' + repr(file_path),
                UserWarning)

        # Set a paired_end attribute if more than one read without index is defined?

        # Get the flow cell layout if it exits.

        xml_flow_cell_layout = run_info_tree.find(path='Run/FlowcellLayout')
        if xml_flow_cell_layout is not None:
            xml_tiles = run_info_tree.find(path='Run/FlowcellLayout/TileSet/Tiles')
            tile_name_list: Optional[List[str]] = None
            if xml_tiles is not None:
                tile_name_list = [xml_tile.text for xml_tile in xml_tiles.iter(tag='Tile')]

            flow_cell_layout = RunInformationFlowcellLayout(
                lane_count=int(xml_flow_cell_layout.get(key='LaneCount')),
                surface_count=int(xml_flow_cell_layout.get(key='SurfaceCount')),
                swath_count=int(xml_flow_cell_layout.get(key='SwathCount')),
                tile_count=int(xml_flow_cell_layout.get(key='TileCount')),
                tile_list=tile_name_list)
        else:
            flow_cell_layout = None

        return cls(file_path=file_path,
                   file_type='xml',
                   name=file_name,
                   run_identifier=run_identifier,
                   run_number=run_number,
                   flow_cell=flow_cell,
                   instrument=instrument,
                   date=date,
                   run_information_read_list=run_information_read_list,
                   flow_cell_layout=flow_cell_layout)

    def __init__(self,
                 file_path=None,
                 file_type=None,
                 name=None,
                 run_identifier=None,
                 run_number=None,
                 flow_cell=None,
                 instrument=None,
                 date=None,
                 flow_cell_layout=None,
                 run_information_read_list=None):
        """Initialise a C{bsf.illumina.RunInformation} object.

        @param file_path: File path
        @type file_path: str | None
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str | None
        @param name: Name
        @type name: str | None
        @param run_identifier: Run identifier e.g. I{130724_SN815_0089_BC26JBACXX}
        @type run_identifier: str | None
        @param run_number: Run number, which may not have to correspond to the run number in the run identifier e.g. 91
        @type run_number: str | None
        @param flow_cell: Illumina flow cell identifier
        @type flow_cell: str | None
        @param instrument: Illumina instrument serial number
        @type instrument: str | None
        @param date: Date in YYMMDD format
        @type date: str | None
        @param flow_cell_layout: C{bsf.illumina.RunInformationFlowcellLayout}
        @type flow_cell_layout: RunInformationFlowcellLayout | None
        @param run_information_read_list: Python C{list} of C{bsf.illumina.RunInformationRead} objects
        @type run_information_read_list: list[RunInformationRead]
        """
        super(RunInformation, self).__init__()

        self.file_path = file_path
        self.file_type = file_type
        self.name = name
        self.run_identifier = run_identifier
        self.run_number = run_number
        self.date = date
        self.instrument = instrument
        self.flow_cell = flow_cell
        self.flow_cell_layout = flow_cell_layout

        if run_information_read_list is None:
            self.run_information_read_list = list()
        else:
            self.run_information_read_list = run_information_read_list

        return

    @property
    def get_cycle_number(self):
        """Get the total number of cycles.

        @return: Number of cycles
        @rtype: int
        """
        cycle_number = 0

        for run_information_read in self.run_information_read_list:
            cycle_number += run_information_read.cycles

        return cycle_number

    @property
    def get_iso_date(self):
        """Get the run start date in ISO 8601 format.

        @return: The run start date in ISO 8601 format.
        @rtype: str
        """
        if self.date is None:
            return
        else:
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
        return len(self.run_information_read_list)

    @property
    def get_read_start_list(self):
        """Get a Python C{list} of cycle numbers at the start of each read.

        @return: Python C{list} of Python C{int} (starting cycle) for each read
        @rtype: list[int]
        """
        cycle_number = 1

        read_start_list: List[int] = list()

        for run_information_read in self.run_information_read_list:
            read_start_list.append(cycle_number)
            cycle_number += run_information_read.cycles

        return read_start_list

    @property
    def get_read_end_list(self):
        """Get a Python C{list} of cycle numbers at the end of each read.

        @return: Python C{list} of Python C{int} (ending cycle) for each read
        @rtype: list[int]
        """
        cycle_number = 1

        read_end_list: List[int] = list()

        for run_information_read in self.run_information_read_list:
            cycle_number += run_information_read.cycles
            read_end_list.append(cycle_number - 1)

        return read_end_list

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

        for run_information_read in self.run_information_read_list:
            read_structure += str(run_information_read.cycles)
            if run_information_read.index:
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
        read_structure_list: List[str] = list()

        for run_information_read in self.run_information_read_list:
            if run_information_read.index:
                read_structure_list.append(str(run_information_read.cycles) + 'I')
            else:
                read_structure_list.append(str(run_information_read.cycles) + 'B')

        return read_structure_list

    def has_tile(self, tile):
        """Convenience method to test the C{RunInformationFlowcellLayout} for a particular tile.

        @param tile: Flowcell tile
        @type tile: str
        @return: C{True} if the tile is in the C{RunInformationFlowcellLayout}, C{False} otherwise
        @rtype: bool
        """
        return self.flow_cell_layout.has_tile(tile=tile)


class RunParameters(object):
    """The C{bsf.illumina.RunParameters} class models the contents of runParameters.xml
    files inside an Illumina Run Folder.

    @ivar file_path: File path
    @type file_path: str
    @ivar element_tree: C{xml.etree.ElementTree.ElementTree}
    @type element_tree: ElementTree | None
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{bsf.illumina.RunParameters} object from a file path.

        @param file_path: File path
        @type file_path: str | None
        @return: C{bsf.illumina.RunParameters} object
        @rtype: RunParameters
        """
        file_path = os.path.normpath(file_path)

        if not os.path.isfile(file_path):
            file_path = None

        return cls(file_path=file_path, element_tree=ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{bsf.illumina.RunParameters} object.

        @param file_path: File path
        @type file_path: str | None
        @param element_tree: C{xml.etree.ElementTree.ElementTree}
        @type element_tree: ElementTree | None
        """
        super(RunParameters, self).__init__()

        self.file_path = file_path
        self.element_tree = element_tree

        return

    def xml_paths_to_text(self, xml_paths):
        """Get the text representation of the first element of a tuple of XML paths.

        @param xml_paths: Python C{tuple} of Python C{str} XML path elements
        @type xml_paths: tuple[str]
        @return: Text representation
        @rtype: str
        """
        for xml_path in xml_paths:
            element = self.element_tree.find(path=xml_path)
            if element is not None:
                return element.text
        else:
            return

    @property
    def get_run_parameters_version(self):
        """Get the run parameters version of a C{bsf.illumina.RunParameters} object.

        Returns the text representation of:
            - I{<RunParameters>/<RunParametersVersion>}: I{MiSeq} and I{NextSeq}
        @return: Run parameters version or an empty string
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('RunParametersVersion',))

    @property
    def get_instrument_type(self):
        """Get the instrument type based on a C{bsf.illumina.RunParameters} object.

        Returns I{HiSeq}, I{MiSeq}, I{NextSeq} or I{NovaSeq} depending on C{get_application_name}
        @return: Instrument type
        @rtype: str
        """
        # Simply remove ' Control Software' from the get_application_name property.
        return self.get_application_name[:-17]

    @property
    def get_experiment_name(self):
        """Get the experiment name of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<ExperimentName>}
            - I{MiSeq}:   I{<RunParameters>/<ExperimentName>}
            - I{NextSeq}: I{<RunParameters>/<ExperimentName>}
            - I{NovaSeq}: I{<RunParameters>/<ExperimentName>}
        @return: Experiment name
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/ExperimentName', 'ExperimentName'))

    @property
    def get_flow_cell_barcode(self):
        """Get the flow cell barcode of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<Barcode>}
            - I{MiSeq}:   I{<RunParameters>/<Barcode>}
            - I{NextSeq}: I{<RunParameters>/<FlowCellSerial>}
            - I{NovaSeq}: I{<RunParameters>/<RfidsInfo>/<FlowCellSerialBarcode>}
        @return: Flow cell barcode
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=(
            'Setup/Barcode', 'Barcode', 'FlowCellSerial', 'RfidsInfo/FlowCellSerialBarcode'))

    @property
    def get_flow_cell_type(self):
        """Get the flow cell chemistry type of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<Flowcell>}
            - I{MiSeq}:   I{<RunParameters>/<ReagentKitVersion>}
            - I{NextSeq}: I{<RunParameters>/<Chemistry>}
            - I{NovaSeq}: I{<RunParameters>/<RfidsInfo>/<FlowCellMode>}
        @return: Flow cell chemistry type
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=(
            'Setup/Flowcell', 'ReagentKitVersion', 'Chemistry', 'RfidsInfo/FlowCellMode'))

    @property
    def get_index_type(self):
        """Get the index chemistry type of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<Index>}
            - I{MiSeq}:   I{<RunParameters>/<Setup>/<Index>}
            - I{NextSeq}: I{None}
            - I{NovaSeq}: I{None}
        @return: Index chemistry type
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/Index',))

    @property
    def get_pe_type(self):
        """Get the paired-end chemistry type of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}: I{<RunParameters>/<Setup>/<Pe>}
            - I{MiSeq}: I{<RunParameters>/<Setup>/<Pe>}
            - I{NextSeq}: I{None}
            - I{NovaSeq}: I{None}
        @return: Paired-end chemistry type
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/Pe',))

    @property
    def get_sbs_type(self):
        """Get the sequencing-by-synthesis chemistry type of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<Sbs>}
            - I{MiSeq}:   I{<RunParameters>/<Setup>/<Sbs>}
            - I{NextSeq}: I{None}
            - I{NovaSeq}: I{None}
        @return: Sequencing-by-synthesis chemistry type
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/Sbs',))

    @property
    def get_position(self):
        """Get the flow cell position of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<FCPosition>}
            - I{MiSeq}:   I{None}
            - I{NextSeq}: I{'A'}
            - I{NovaSeq}: I{<RunParameters>/<Side>}

        The I{NextSeq} has no concept of I{<FCPosition>}, but always uses 'A'
        @return: Flow cell position e.g. A or B
        @rtype: str
        """
        if self.get_instrument_type in ('NextSeq',):
            return 'A'
        else:
            return self.xml_paths_to_text(xml_paths=('Setup/FCPosition', 'Side'))

    @property
    def get_run_identifier(self):
        """Get the run identifier of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<RunID>}
            - I{MiSeq}:   I{<RunParameters>/<RunID>}
            - I{NextSeq}: I{<RunParameters>/<RunID>}
            - I{NovaSeq}: I{<RunParameters>/<RunId>}
        @return: Run identifier
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/RunID', 'RunID', 'RunId'))

    @property
    def get_real_time_analysis_version(self):
        """Get the Real-Time Analysis (RTA) Version of a C{bsf.illumina.RunParameters} object.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<RTAVersion>}
            - I{MiSeq}:   I{<RunParameters>/<RTAVersion>}
            - I{NextSeq}: I{<RunParameters>/<RTAVersion>}
            - I{NovaSeq}: I{<RunParameters>/<RtaVersion>}
        @return: RTA version
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/RTAVersion', 'RTAVersion', 'RtaVersion'))

    @property
    def get_application_name(self):
        """Get the application (i.e. I{HiSeq}, I{MiSeq}, I{NextSeq} or I{NovaSeq Control Software}) name.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<ApplicationName>}
            - I{MiSeq}:   I{<RunParameters>/<Setup>/<ApplicationName>}
            - I{NextSeq}: I{<RunParameters>/<Setup>/<ApplicationName>}
            - I{NovaSeq}: I{<RunParameters>/<Application>}
        @return: Application name
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/ApplicationName', 'Application'))

    @property
    def get_application_version(self):
        """Get the application (i.e. I{HiSeq}, I{MiSeq}, I{NextSeq} or I{NovaSeq Control Software}) version.

            - I{HiSeq}:   I{<RunParameters>/<Setup>/<ApplicationVersion>}
            - I{MiSeq}:   I{<RunParameters>/<Setup>/<ApplicationVersion>}
            - I{NextSeq}: I{<RunParameters>/<Setup>/<ApplicationVersion>}
            - I{NovaSeq}: I{<RunParameters>/<ApplicationVersion>}
        @return: Application version
        @rtype: str
        """
        return self.xml_paths_to_text(xml_paths=('Setup/ApplicationVersion', 'ApplicationVersion'))


class XMLConfiguration(object):
    """The C{bsf.illumina.XMLConfiguration} class models the contents of XML configuration
    files inside an Illumina Run Folder.

    @ivar file_path: File path
    @type file_path: str | None
    @ivar element_tree: C{xml.etree.ElementTree.ElementTree}
    @type element_tree: ElementTree | None
    """

    @classmethod
    def from_file_path(cls, file_path):
        """Create a C{bsf.illumina.XMLConfiguration} object from a file path.

        In case the file path does not exist a C{bsf.illumina.XMLConfiguration} object
        with an empty C{xml.etree.ElementTree.ElementTree} will be returned.
        @param file_path: File path
        @type file_path: str
        @return: C{bsf.illumina.XMLConfiguration} object
        @rtype: XMLConfiguration
        """
        file_path = os.path.normpath(file_path)

        if not os.path.isfile(file_path):
            file_path = None

        return cls(file_path=file_path, element_tree=ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{bsf.illumina.XMLConfiguration} object.

        @param file_path: File path
        @type file_path: str | None
        @param element_tree: C{xml.etree.ElementTree.ElementTree}
        @type element_tree: ElementTree | None
        """
        super(XMLConfiguration, self).__init__()

        self.file_path = file_path
        self.element_tree = element_tree

        return


class AnalysisConfiguration(XMLConfiguration):
    """The C{bsf.illumina.AnalysisConfiguration} class models Image and Base Call analysis
    XML configuration files inside an Illumina Run Folder.

    @ivar file_path: File path
    @type file_path: str
    @ivar element_tree: C{xml.etree.ElementTree.ElementTree}
    @type element_tree: ElementTree
    """

    def __init__(self, file_path=None, element_tree=None):
        """Initialise a C{bsf.illumina.AnalysisConfiguration} object.
        Read all I{<Tile>} elements from the I{<Run>/<TileSelection>} element from the XML configuration file
        to initialise an internal Python C{dict} of valid tiles.

        @param file_path: File path
        @type file_path: str
        @param element_tree: C{xml.etree.ElementTree.ElementTree}
        @type element_tree: ElementTree
        """
        super(AnalysisConfiguration, self).__init__(file_path=file_path, element_tree=element_tree)

        self._lane_tile_dict: Dict[str, Dict[str, bool]] = dict()

        if self.element_tree.getroot() is None:
            return

        for lane_element in self.element_tree.find(path='Run/TileSelection').iter(tag='Lane'):
            lane_index = lane_element.get(key='Index')
            if lane_index not in self._lane_tile_dict:
                self._lane_tile_dict[lane_index] = dict()
            lane_dict = self._lane_tile_dict[lane_index]
            for tile_element in lane_element.findall(path='Tile'):
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
            return tile in self._lane_tile_dict[lane]
        else:
            return


class ImageAnalysis(AnalysisConfiguration):
    """The C{bsf.illumina.ImageAnalysis} class models the contents of the
    I{IRF/Data/Intensities/config.xml} XML configuration file inside an Illumina Run Folder.
    """
    pass


class BaseCallAnalysis(AnalysisConfiguration):
    """The C{bsf.illumina.BaseCallAnalysis} class models the contents of the
    I{IRF/Data/Intensities/BaseCalls/config.xml} XML configuration file inside an Illumina Run Folder.
    """
    pass


class RunFolderNotComplete(Exception):
    pass


class RunFolder(object):
    """The C{bsf.illumina.RunFolder} class represents an Illumina Run Folder.

    @cvar rta_dict: Python C{dict} of Python C{str} (RTA version) key objects and
        Python C{str} (RTA description) value objects
    @type rta_dict: dict[str, str]
    @ivar file_path: File path
    @type file_path: str
    @ivar file_type: File type
        I{CASAVA}: FASTQ file after post-processing with CASAVA
        I{External}: other data files
    @type file_type: str
    @ivar date: Date in YYMMDD format
    @type date: str | None
    @ivar instrument: Illumina instrument serial number
    @type instrument: str | None
    @ivar run: Run serial number
    @type run: str | None
    @ivar flow_cell: Flow cell identifier
    @type flow_cell: str | None
    @ivar run_information: C{bsf.illumina.RunInformation}
    @type run_information: RunInformation
    @ivar run_parameters: C{bsf.illumina.RunParameters}
    @type run_parameters: RunParameters
    @ivar image_analysis: C{bsf.illumina.ImageAnalysis}
    @type image_analysis: ImageAnalysis
    @ivar base_call_analysis: C{bsf.illumina.BaseCallAnalysis}
    @type base_call_analysis: BaseCallAnalysis
    """

    rta_dict = {
        '1.12.4': 'HiSeq Control Software 1.4.5 (HiSeq 2000)',
        '1.12.4.2': 'HiSeq Control Software 1.4.8 (HiSeq 2000)',
        '1.13.48': 'HiSeq Control Software 1.5.15.1 (HiSeq 1000/2000) 2012-02-17',
        '1.17.21.3': 'HiSeq Control Software 2.0.12.0 (HiSeq 1000/1500/2000/2500) 2013-06-12',
        '1.18.54': 'MiSeq Control Software 2.5.0.5 (MiSeq) 2014-09-15',
        '1.18.54.4': 'MiSeq Control Software 4.0.0.1769 (MiSeq) 2021-01-26',
        '1.18.61': 'HiSeq Control Software 2.2.38 (HiSeq 1000/1500/2000/2500) 2014-07',
        '1.18.64': 'HiSeq Control Software 2.2.58 (HiSeq 1000/1500/2000/2500) 2014-11-14',
        '1.18.66': 'HiSeq Control Software 2.2.68 (HiSeq 1000/1500/2000/2500) 2015-05-26',
        '2.4.11': 'NextSeq Control Software 2.1.0.31 (NextSeq 500/550) 2017-10-24',
        '2.5.2': 'HiSeq Control Software 3.3.20 (HiSeq 3000/4000) 2015-05-05',
        '2.7.3': 'HiSeq Control Software 3.3.52 (HiSeq 3000/4000)',
        '2.7.6': 'HiSeq Control Software 3.3.76 (HiSeq 3000/4000)',
        '2.7.7': 'HiSeq Control Software HD 3.4.0.38 (HiSeq 3000/4000)',
        'v3.3.3': 'NovaSeq Control Software 1.2.0 (NovaSeq 6000)',
        'v3.4.4': 'NovaSeq Control Software 1.6.0 (NovaSeq 6000)',
    }

    # '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 1000/1500/2000/2500
    # '1.18.54', '1.18.54.4',  # MiSeq
    # '2.4.11',  # NextSeq 500/550
    # '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
    # 'v3.3.3', 'v3.4.4',  # NovaSeq 6000

    @classmethod
    def from_file_path(cls, file_path):
        """Construct a C{bsf.illumina.RunFolder} object from a file path.

        @param file_path: File path
        @type file_path: str
        @return: C{bsf.illumina.RunFolder}
        @rtype: RunFolder
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
            raise Exception('Found neither a runParameters.xml nor a RunParameters.xml file.')

        # Illumina Run Folders obey a 'YYMMDD_SN000_Run_PositionFCID' schema.

        if isinstance(file_name, bytes):
            file_name = file_name.decode()

        component_list = file_name.split('_')

        return cls(
            file_path=file_path,
            file_type='Illumina',
            date=component_list[0],
            instrument=component_list[1],
            run=component_list[2],
            flow_cell=component_list[3],
            run_information=RunInformation.from_file_path(file_path=os.path.join(file_path, 'RunInfo.xml')),
            run_parameters=RunParameters.from_file_path(file_path=run_parameters_path),
            image_analysis=ImageAnalysis.from_file_path(file_path=os.path.join(
                file_path, 'Data', 'Intensities', 'config.xml')),
            base_call_analysis=BaseCallAnalysis.from_file_path(file_path=os.path.join(
                file_path, 'Data', 'Intensities', 'BaseCalls', 'config.xml')))

    def __init__(
            self,
            file_path=None,
            file_type=None,
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
        @type file_path: str
        @param file_type: File type (e.g. I{CASAVA}, I{External} or I{Automatic})
        @type file_type: str
        @param date: Date in I{YYMMDD} format
        @type date: str | None
        @param instrument: Illumina instrument serial number (e.g. I{SN181}, I{SN815}, ...)
        @type instrument: str | None
        @param run: Run serial number
        @type run: str | None
        @param flow_cell: The position and flow cell identifier
        @type flow_cell: str | None
        @param run_information: C{bsf.illumina.RunInformation}
        @type run_information: RunInformation
        @param run_parameters: C{bsf.illumina.RunParameters}
        @type run_parameters: RunParameters
        @param image_analysis: C{bsf.illumina.ImageAnalysis}
        @type image_analysis: ImageAnalysis
        @param base_call_analysis: C{bsf.illumina.BaseCallAnalysis}
        @type base_call_analysis: BaseCallAnalysis
        """
        super(RunFolder, self).__init__()

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        self.file_type = file_type
        self.date = date
        self.instrument = instrument
        self.run = run
        self.flow_cell = flow_cell

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

        self._missing_base_call_tiles: Dict[int, Dict[str, bool]] = dict()

        self._missing_image_analysis_tiles: Dict[int, Dict[str, bool]] = dict()

        return

    @property
    def get_base_calls_directory(self):
        """Get the base-calls directory in the I{IRF/Data/Intensities/BaseCalls} hierarchy.

        @return: Illumina base-calls directory
        @rtype: str
        """
        return os.path.join(self.file_path, 'Data', 'Intensities', 'BaseCalls')

    @property
    def get_name(self):
        """Get the C{bsf.illumina.RunFolder} name.

        @return: Name
        @rtype: str
        """
        return '_'.join((self.date, self.instrument, self.run, self.flow_cell))

    def get_sav_archive_paths(self):
        """Get Illumina Sequence Analysis Viewer (SAV) archive file paths.

        Return paths for files to be extracted from an Illumina Run Folder (IRF) in GNU Tar format,
        which are relevant for viewing in the Illumina Sequence Analysis Viewer (SAV).
        The file paths depend on the Realtime Analysis (RTA) software version.
        @return: Python C{list} object of Python C{str} file path objects
        @rtype: list[str]
        """
        str_list = list()

        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in self.rta_dict:
            raise Exception('Unsupported RTA version: ' + repr(rta))

        # On all instruments.
        str_list.append('Config')
        str_list.append('InterOp')
        str_list.append('RTAComplete.txt')
        str_list.append('RunInfo.xml')

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            str_list.append('Data/Intensities/config.xml')
            str_list.append('Data/Intensities/BaseCalls/config.xml')
            str_list.append('Data/Intensities/RTAConfiguration.xml')

        if rta in (
                '2.4.11',  # NextSeq 500/550
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
        ):
            str_list.append('RTAConfiguration.xml')

        if rta in (
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            str_list.append('RTA3.cfg')

        if rta in (
                '1.18.54.4',  # MiSeq Control Software 4.0.0.1769 (MiSeq)
                '2.4.11',  # NextSeq 500/550
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            str_list.append('RunParameters.xml')
        else:
            str_list.append('runParameters.xml')

        return ['/'.join((self.get_name, item)) for item in str_list]

    def has_compressed_base_calls(self):
        """The Illumina Run Folder has un-compressed base call files.

        @return: True if un-compressed files exist, False otherwise
        @rtype: bool
        """
        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in self.rta_dict:
            raise Exception('Unsupported RTA version: ' + repr(rta))

        if rta in (
                '1.18.54',
                '1.18.54.4',
                # 1.18.54 MiSeq Control Software 2.6.2.1 (MiSeq)
                # 1.18.54.4 MiSeq Control Software 4.0.0.1769 (MiSeq)
                #   IRF/Data/Intensities/BaseCalls/L001/C1.1/s_1_1101.bcl
                #   IRF/Data/Intensities/BaseCalls/L001/C1.1/s_1_1101.stats
                #   Uncompressed base calls, but for consistency they are kept uncompressed.
                '2.4.11',
                # 2.4.11 NextSeq Control Software 2.2.0.4 (NextSeq 500/550)
                #   IRF/Data/Intensities/BaseCalls/L001/*.bcl.bgzf
                #   IRF/Data/Intensities/BaseCalls/L001/*.bcl.bgzf.bci
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',
                # 2.5.2 HiSeq Control Software 3.3.20 (HiSeq 3000/4000)
                # 2.7.3 HiSeq Control Software 3.3.52 (HiSeq 3000/4000)
                # 2.7.6 HiSeq Control Software 3.3.76 (HiSeq 3000/4000)
                # 2.7.7 HiSeq Control Software HD 3.4.0.38 (HiSeq 3000/4000)
                #   IRF/Data/Intensities/BaseCalls/L001/C1.1/s_1_1101.bcl.gz
                'v3.3.3', 'v3.4.4',
                # v3.3.3 NovaSeq Control Software 1.2.0 (NovaSeq 6000)
                # v3.4.4 NovaSeq Control Software 1.6.0 (NovaSeq 6000)
                #   IRF/Data/Intensities/BaseCalls/L001/C1.1/L001_1.cbcl
        ):
            return True
        else:
            return False

    def has_intensities(self):
        """The Illumina Run Folder has intensity files.

        @return: True if intensity files exist, False otherwise
        @rtype: bool
        """
        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in self.rta_dict:
            raise Exception('Unsupported RTA version: ' + repr(rta))

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
        ):
            return True
        else:
            return False

    def _check_tiles_base_call(self):
        """Check for missing I{<Tile>} elements in the I{IRF/Data/Intensities/BaseCalls/config.xml}
        configuration file. This method also builds up a Python C{dict} required for method
        C{_is_missing_base_call_tile}.
        """
        fcl = self.run_information.flow_cell_layout

        for lane in range(0 + 1, fcl.lane_count + 1):
            if lane not in self._missing_base_call_tiles:
                self._missing_base_call_tiles[lane] = dict()
            lane_dict = self._missing_base_call_tiles[lane]
            if self.base_call_analysis.has_lane(lane=str(lane)):
                # Lanes are not defined, if the config.xml file could not be read.
                for surface in range(0 + 1, fcl.surface_count + 1):
                    for swath in range(0 + 1, fcl.swath_count + 1):
                        for tile in range(0 + 1, fcl.tile_count + 1):
                            tile = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if not self.base_call_analysis.has_lane_tile(lane=str(lane), tile=tile):
                                lane_dict[tile] = True
                                print('Missing BaseCallAnalysis lane:', lane, 'tile:', tile)

        return

    def _check_tiles_image_analysis(self):
        """Check for missing I{<Tile>} elements in the XML configuration file.

        Check the I{IRF/Data/Intensities/config.xml} configuration file for missing I{<Tile>} XML elements.
        This method also builds up a Python C{dict} required for method C{_is_missing_image_analysis_tile}.
        """
        fcl = self.run_information.flow_cell_layout

        for lane in range(0 + 1, fcl.lane_count + 1):
            if lane not in self._missing_image_analysis_tiles:
                self._missing_image_analysis_tiles[lane] = dict()
            lane_dict = self._missing_image_analysis_tiles[lane]
            if self.image_analysis.has_lane(lane=str(lane)):
                # Lanes are not defined of the config.xml file could not be read.
                for surface in range(0 + 1, fcl.surface_count + 1):
                    for swath in range(0 + 1, fcl.swath_count + 1):
                        for tile in range(0 + 1, fcl.tile_count + 1):
                            tile = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if not self.image_analysis.has_lane_tile(lane=str(lane), tile=tile):
                                lane_dict[tile] = True
                                print('Missing ImageAnalysis lane:', lane, 'tile:', tile)

        return

    def _is_missing_base_call_tile(self, lane, tile):
        """Confirm that a particular I{<Tile>} element is missing from the XML configuration file.

        Check the I{IRF/Data/Intensities/config.xml} configuration file for missing I{<Tile>} XML elements.
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
        """Confirm that a particular I{<Tile>} element is missing from the XML configuration file.

        Check the I{IRF/Data/Intensities/BaseCalls/config.xml} configuration file for missing I{<Tile>} XML elements.
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
    def _check_file_names(directory_dict, directory_path, file_name_list, debug=0):
        """Check a Python C{list} of file names against a Python C{dict} of directory entries.

        @param directory_dict: Python C{dict} of directory entries
        @type directory_dict: dict[str, int]
        @param directory_path: Directory path
        @type directory_path: str
        @param file_name_list: Python C{list} of file names
        @type file_name_list: list[str]
        @param debug: Integer debugging level
        @type debug: int
        """
        if debug > 0:
            # print('Processing directory', directory_path)
            pass

        for file_name in file_name_list:
            file_path = os.path.join(directory_path, file_name)
            if file_name in directory_dict:
                del directory_dict[file_name]
            else:
                print('Missing file', file_path)

        return

    @staticmethod
    def _check_file_suffixes(directory_dict, directory_path, file_suffix_list, debug=0):
        """Check a Python C{list} of file suffixes against a Python C{dict} of directory entries.

        @param directory_dict: Python C{dict} of directory entries
        @type directory_dict: dict[str, int]
        @param directory_path: Directory path
        @type directory_path: str
        @param file_suffix_list: Python C{list} of file names
        @type file_suffix_list: list[str]
        @param debug: Integer debugging level
        @type debug: int
        """
        if debug > 0:
            # print('Processing directory', directory_path)
            pass

        for file_suffix in file_suffix_list:
            for file_name in directory_dict:
                if file_name.endswith(file_suffix):
                    del directory_dict[file_name]
                    break
            else:
                file_path = os.path.join(directory_path, '*' + file_suffix)
                print('Missing file suffix', file_path)

        return

    def _check_config(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Config/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        rta = self.run_parameters.get_real_time_analysis_version
        # flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.upper()

        _directory_name = 'Config'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        _file_name_list: List[str] = list()
        _file_suffix_list: List[str] = list()

        if rta in (
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            _file_name_list.append('Effective.cfg')
            _file_name_list.append('MiSeqOverride.cfg')
            _file_name_list.append('RTAStart.bat')

        if rta in (
                '1.12.4',  # HiSeq Control Software 1.4.5 (HiSeq 2000)
        ):
            _file_name_list.append('HiSeqControlSoftware.Options.cfg')
            _file_name_list.append('RTAStart.bat')
            _file_name_list.append('Variability_HiSeq.xml')

        if rta in (
                '2.7.7',  # HiSeq Control Software HD 3.4.0.38 (HiSeq 3000/4000)
        ):
            _file_name_list.append('HiSeqControlSoftware.Options.cfg')
            _file_name_list.append('RTAStart.log')
            _file_name_list.append('Variability_HiSeq_E.bin')
            _file_suffix_list.append('_Effective.cfg')
            _file_suffix_list.append('_Override.cfg')
            _file_suffix_list.append('_SortedOverride.cfg')

        if rta in (
                '2.4.11',  # NextSeq 500/550
        ):
            _file_name_list.append('Effective.cfg')
            _file_name_list.append('FirmwareVersions.txt')
            _file_name_list.append('NextSeqCalibration.cfg')
            _file_name_list.append('NextSeqOverride.cfg')

        if rta in (
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('Effective.cfg')
            _file_name_list.append('LaserPowerVariability.xml')
            _file_name_list.append('NovaSeqCalibration.cfg')
            _file_name_list.append('NovaSeqOverride.cfg')
            _file_name_list.append('Options.cfg')

        self._check_file_names(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        self._check_file_suffixes(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_suffix_list=_file_suffix_list,
            debug=debug)

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data_intensities_base_calls_matrix(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/Matrix/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/Intensities/BaseCalls/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: I{IRF/Data/Intensities/BaseCalls/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'Matrix'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        if rta in (
                '2.5.2',  # HiSeq 3000/4000
        ):
            # Process IRF/Data/Intensities/BaseCalls/Matrix/L00[1-8] directories.

            for lane in range(0 + 1, fcl.lane_count + 1):
                lane_name = 'L{:03d}'.format(lane)
                lane_path = os.path.join(_directory_path, lane_name)

                if lane_name in _directory_dict:
                    del _directory_dict[lane_name]
                else:
                    print('Missing directory', lane_path)
                    continue

                lane_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

                if debug > 1:
                    print('Processing directory', lane_path)

                # Process IRF/Data/Intensities/BaseCalls/Matrix/L00[1-8]/C[0-9]+.1/ directories.

                for cycle in range(0 + 1, self.run_information.get_cycle_number + 1):
                    cycle_name = 'C{:d}.1'.format(cycle)
                    cycle_path = os.path.join(lane_path, cycle_name)
                    if cycle_name in lane_dict:
                        del lane_dict[cycle_name]
                    else:
                        print('Missing directory', cycle_path)
                        continue
                    cycle_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                    if debug > 2:
                        print('Processing directory', cycle_path)

                    for surface in range(0 + 1, fcl.surface_count + 1):
                        for swath in range(0 + 1, fcl.swath_count + 1):
                            for tile in range(0 + 1, fcl.tile_count + 1):
                                # Process tile matrix files.
                                # s_1_1101_matrix.txt
                                # s_1_2228_matrix.txt
                                _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}_matrix.txt'.format(lane, surface, swath, tile)
                                if _entry_name in cycle_dict:
                                    del cycle_dict[_entry_name]
                                else:
                                    print('Missing tile matrix file', os.path.join(cycle_path, _entry_name))

                    if len(cycle_dict):
                        print(cycle_path, 'with number of entries:', str(len(cycle_dict)))
                        print('  Remaining entries:', sorted(cycle_dict))

                if len(lane_dict):
                    print(lane_path, 'with number of entries:', str(len(lane_dict)))
                    print('  Remaining entries:', sorted(lane_dict))
        else:
            # All other instruments have a flat list of matrix.txt files.
            for read in range(0 + 1, self.run_information.get_read_number + 1):
                # Process read matrix files.
                # s_1_matrix.txt
                # s_2_matrix.txt
                _entry_name = 's_{:d}_matrix.txt'.format(read)

                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print('Missing read matrix file', os.path.join(_directory_path, _entry_name))

                for lane in range(0 + 1, fcl.lane_count + 1):
                    # Process lane matrix files.
                    # s_1_1_matrix.txt
                    # s_8_2_matrix.txt
                    _entry_name = 's_{:d}_{:d}_matrix.txt'.format(lane, read)

                    if _entry_name in _directory_dict:
                        del _directory_dict[_entry_name]
                    else:
                        print('Missing lane matrix file', os.path.join(_directory_path, _entry_name))

                    for surface in range(0 + 1, fcl.surface_count + 1):
                        for swath in range(0 + 1, fcl.swath_count + 1):
                            for tile in range(0 + 1, fcl.tile_count + 1):
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
                                    print('Missing tile matrix file', os.path.join(_directory_path, _entry_name))

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data_intensities_base_calls_phasing(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/Phasing/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/intensities/BaseCalls/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: I{IRF/Data/intensities/BaseCalls/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'Phasing'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', os.path.join(directory_path, _directory_name))
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        for run_information_read in self.run_information.run_information_read_list:
            # Process read phasing files.
            # s_1_phasing.txt
            # s_2_phasing.txt
            _entry_name = 's_{:d}_phasing.txt'.format(run_information_read.number)

            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing file', os.path.join(_directory_path, _entry_name))

            for lane in range(0 + 1, fcl.lane_count + 1):
                # Process lane phasing files.
                # s_1_1_phasing.txt
                # s_8_3_phasing.txt
                _entry_name = 's_{:d}_{:d}_phasing.txt'.format(lane, run_information_read.number)

                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print('Missing file', os.path.join(_directory_path, _entry_name))
                    continue

                for surface in range(0 + 1, fcl.surface_count + 1):
                    for swath in range(0 + 1, fcl.swath_count + 1):
                        for tile in range(0 + 1, fcl.tile_count + 1):
                            # Not all tiles have to exists especially after catastrophic events during the
                            # cluster generation step.
                            tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                                continue

                            if not run_information_read.index and rta in (
                                    '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3',  # HiSeq 2000
                            ):
                                # Process tile cycle files, which only exist for payload, but not index reads.
                                # s_1_1_1101_cycle.txt
                                # s_8_1_2316_cycle.txt
                                _entry_name = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_cycle.txt'. \
                                    format(lane, run_information_read.number, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print('Missing file', os.path.join(_directory_path, _entry_name))

                            # Process tile phasing files.
                            # s_1_1_1101_phasing.txt
                            # s_8_3_2316_phasing.txt
                            _entry_name = 's_{:1d}_{:1d}_{:1d}{:1d}{:02d}_phasing.txt'. \
                                format(lane, run_information_read.number, surface, swath, tile)
                            if _entry_name in _directory_dict:
                                del _directory_dict[_entry_name]
                            else:
                                print('Missing file', os.path.join(_directory_path, _entry_name))

                            if rta in (
                                    '1.18.61', '1.18.64',  # HiSeq 2000
                                    '1.18.54', '1.18.54.4',  # MiSeq
                                    '2.5.2',  # HiSeq 3000/4000
                            ):
                                # Process the tile empirical phasing files.
                                _entry_name = 'EmpiricalPhasingCorrection_{:1d}_{:1d}_{:1d}{:1d}{:02d}.txt'. \
                                    format(lane, run_information_read.number, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print('Missing file', os.path.join(_directory_path, _entry_name))

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54',  # MiSeq Control Software 2.5.0.5 (MiSeq)
        ):
            # With the exception of RTA 1.18.54.4 (MiSeq) and RTA 2.5.2 (HiSeq 3000/4000),
            # all other HiSeq and MiSeq RTAs have
            # IRF/Data/Intensities/BaseCalls/Phasing/s_{lane}_{cycle}_phasing.xml files.
            for lane in range(0 + 1, fcl.lane_count + 1):
                for read_start in self.run_information.get_read_start_list:
                    _entry_name = 's_{:1d}_{:02d}_phasing.xml'.format(lane, read_start + 1)
                    if _entry_name in _directory_dict:
                        del _directory_dict[_entry_name]
                    else:
                        print('Missing file', os.path.join(_directory_path, _entry_name))

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data_intensities_base_calls(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/BaseCalls/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/Intensities/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: I{IRF/Data/Intensities/} file path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'BaseCalls'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', os.path.join(directory_path, _directory_name))
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        # Process the IRF/Data/Intensities/BaseCalls/config.xml file.

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54', '1.18.54.4',  # Miseq
        ):
            _entry_name = 'config.xml'

            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing file', os.path.join(_directory_path, _entry_name))

            # Check for completeness of tiles in the base call configuration XML file (config.xml).
            self._check_tiles_base_call()

        # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/ directories.

        for lane in range(0 + 1, fcl.lane_count + 1):
            lane_name = 'L{:03d}'.format(lane)
            lane_path = os.path.join(_directory_path, lane_name)

            if lane_name in _directory_dict:
                del _directory_dict[lane_name]
            else:
                print('Missing directory', lane_path)
                continue

            lane_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

            if debug > 1:
                print('Processing directory', lane_path)

            # Process IRF/Data/Intensities/BaseCalls/L00[1-8]/C[0-9]+.1 directories.

            if rta in (
                    '2.4.11',  # NextSeq 500/550
            ):
                for cycle in range(0 + 1, self.run_information.get_cycle_number + 1):
                    _entry_name = '{:04d}.bcl.bgzf'.format(cycle)
                    if _entry_name in lane_dict:
                        del lane_dict[_entry_name]
                    else:
                        print('Missing cycle bcl.bgzf file', os.path.join(lane_path, _entry_name))

                    _entry_name = '{:04d}.bcl.bgzf.bci'.format(cycle)
                    if _entry_name in lane_dict:
                        del lane_dict[_entry_name]
                    else:
                        print('Missing cycle bcl.bgzf.bci file', os.path.join(lane_path, _entry_name))

                _entry_name = 's_{:d}.bci'.format(lane)
                if _entry_name in lane_dict:
                    del lane_dict[_entry_name]
                else:
                    print('Missing lane bci file', os.path.join(lane_path, _entry_name))

                _entry_name = 's_{:d}.filter'.format(lane)
                if _entry_name in lane_dict:
                    del lane_dict[_entry_name]
                else:
                    print('Missing lane filter file', os.path.join(lane_path, _entry_name))
            else:
                # For all other instruments.
                for cycle in range(0 + 1, self.run_information.get_cycle_number + 1):
                    cycle_name = 'C{:d}.1'.format(cycle)
                    cycle_path = os.path.join(lane_path, cycle_name)
                    if cycle_name in lane_dict:
                        del lane_dict[cycle_name]
                    else:
                        print('Missing directory', cycle_path)
                        continue

                    cycle_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                    if debug > 2:
                        print('Processing directory', cycle_path)

                    for surface in range(0 + 1, fcl.surface_count + 1):
                        if rta in (
                                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
                        ):
                            # NovaSeq 6000 has only L001_<surface>.cbcl files.
                            # To make matters even more complicated, NovaSeq SP flow cells only
                            # have base calls for the second surface.
                            if surface == 1 and self.run_parameters.get_flow_cell_type == 'SP':
                                continue

                            _entry_name = '{}_{:d}.cbcl'.format(lane_name, surface)
                            if _entry_name in cycle_dict:
                                del cycle_dict[_entry_name]
                            else:
                                print('Missing cbcl file', os.path.join(cycle_path, _entry_name))
                        else:
                            for swath in range(0 + 1, fcl.swath_count + 1):
                                for tile in range(0 + 1, fcl.tile_count + 1):
                                    # Not all tiles have to exist, especially after catastrophic events during the
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
                                            print('Missing tile bcl file', os.path.join(cycle_path, _entry_name))

                                    # Process tile stats files.
                                    # s_1_1101.stats
                                    # s_1_2316.stats
                                    if rta not in (
                                            '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                                    ):
                                        # HiSeq 3000/4000 does not have stats files.
                                        _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.stats'.format(
                                            lane, surface, swath, tile)
                                        if _entry_name in cycle_dict:
                                            del cycle_dict[_entry_name]
                                        else:
                                            print('Missing tile stats file', os.path.join(cycle_path, _entry_name))

                    if len(cycle_dict):
                        print(cycle_path, 'with number of entries:', str(len(cycle_dict)))
                        print('  Remaining entries:', sorted(cycle_dict))

            # Process control and filter files.

            if rta not in (
                    '2.4.11',  # NextSeq 500/550
            ):
                for surface in range(0 + 1, fcl.surface_count + 1):
                    for swath in range(0 + 1, fcl.swath_count + 1):
                        for tile in range(0 + 1, fcl.tile_count + 1):
                            # Not all tiles have to exist, especially after catastrophic events during the
                            # cluster generation step.
                            tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                            if self._is_missing_base_call_tile(lane=lane, tile=tile_name):
                                continue

                            # Process tile control files.
                            # s_1_1101.control
                            # s_1_2316.control
                            if rta not in (
                                    '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                                    'v3.3.3', 'v3.4.4',  # NovaSeq 6000
                            ):
                                # HiSeq 3000/4000 and NovaSeq 6000 instruments do not have control files.
                                _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.control'.format(lane, surface, swath, tile)
                                if _entry_name in lane_dict:
                                    del lane_dict[_entry_name]
                                else:
                                    print('Missing tile control file', os.path.join(lane_path, _entry_name))

                            # Process tile filter files.
                            # s_1_1101.filter
                            # s_1_2316.filter
                            _tile_name = '{:1d}_{:1d}{:1d}{:02d}'.format(lane, surface, swath, tile)
                            if self.run_information.has_tile(tile=_tile_name):
                                _entry_name = 's_{:1d}_{:1d}{:1d}{:02d}.filter'.format(lane, surface, swath, tile)
                                if _entry_name in lane_dict:
                                    del lane_dict[_entry_name]
                                else:
                                    print('Missing tile filter file', os.path.join(lane_path, _entry_name))

            if len(lane_dict):
                print(lane_path, 'with number of entries:', str(len(lane_dict)))
                print('  Remaining entries:', sorted(lane_dict))

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54', '1.18.54.4',  # MiSeq
                '2.5.2',  # HiSeq 3000/4000
        ):
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
            # Only the MiSeq instrument has this file.
            _entry_name = 'SampleSheet.csv'

            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing file', os.path.join(_directory_path, _entry_name))

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data_intensities_offsets(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/Offsets/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/Intensities/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: I{IRF/Data/Intensities/} file path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        rta = self.run_parameters.get_real_time_analysis_version

        if rta in ('2.4.11',):
            return

        _directory_name = 'Offsets'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        # Check the IRF/Data/Intensities/Offsets/offsets.txt file.

        _entry_name = 'offsets.txt'
        if _entry_name in _directory_dict:
            del _directory_dict[_entry_name]
        else:
            print('Missing file', os.path.join(_directory_path, _entry_name))

        # Check the IRF/Data/Intensities/Offsets/SubTileOffsets.txt file.

        _entry_name = 'SubTileOffsets.txt'
        if _entry_name in _directory_dict:
            del _directory_dict[_entry_name]
        else:
            print('Missing file', os.path.join(_directory_path, _entry_name))

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data_intensities(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/Intensities/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: I{IRF/Data/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'Intensities'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        # Build a list of cycle numbers that have no error map such as the last cycle of a read and index cycles.
        no_error_cycles: List[int] = list()

        cycles = 1
        for run_information_read in self.run_information.run_information_read_list:
            if run_information_read.index:
                no_error_cycles.extend(range(cycles, cycles + run_information_read.cycles))
            else:
                no_error_cycles.append(cycles + run_information_read.cycles - 1)
            cycles += run_information_read.cycles

        if debug > 0:
            print('Cycles without errorMap files:', no_error_cycles)

        # Process the IRF/Data/Intensities/BaseCalls/ directory.

        self._check_data_intensities_base_calls(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        if rta in (
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            # Process the IRF/Data/Intensities/s.locs file.
            _entry_name = 's.locs'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing file', os.path.join(_directory_path, _entry_name))
        else:
            # The HiSeq 2000 instrument has:
            # config.xml
            # L00[1-8]
            # Offsets
            # RTAConfiguration.xml

            if rta not in (
                    '2.4.11',  # NextSeq 500/550
            ):
                # Process the IRF/Data/Intensities/config.xml file.

                _entry_name = 'config.xml'
                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print('Missing configuration file', os.path.join(_directory_path, _entry_name))

            # Check for completeness of tiles in the image analysis configuration XML file (config.xml).

            self._check_tiles_image_analysis()

            # Process IRF/Data/Intensities/L00[1-8]/ directories.

            for lane in range(0 + 1, fcl.lane_count + 1):
                lane_name = 'L{:03d}'.format(lane)
                lane_path = os.path.join(_directory_path, lane_name)

                if lane_name in _directory_dict:
                    del _directory_dict[lane_name]
                else:
                    print('Missing directory', lane_path)
                    continue

                lane_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

                if debug > 1:
                    print('Processing directory', lane_path)

                if rta in (
                        '2.4.11',  # NextSeq 500/550
                ):
                    _entry_name = 's_{:d}.locs'.format(lane)
                    if _entry_name in lane_dict:
                        del lane_dict[_entry_name]
                    else:
                        print('Missing file', os.path.join(lane_path, _entry_name))
                else:
                    for surface in range(0 + 1, fcl.surface_count + 1):
                        for swath in range(0 + 1, fcl.swath_count + 1):
                            for tile in range(0 + 1, fcl.tile_count + 1):
                                # Not all tiles have to exists especially after catastrophic events during the
                                # cluster generation step.
                                tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                                if self._is_missing_image_analysis_tile(lane=lane, tile=tile_name):
                                    continue

                                if rta in (
                                        '1.18.54', '1.18.54.4',  # MiSeq
                                ):
                                    # The MiSeq instrument uses locs files. Sigh.
                                    # s_1_1101.locs
                                    # s_1_2119.locs
                                    _entry_name = 's_{}_{:d}{:d}{:02d}.locs'.format(lane, surface, swath, tile)
                                    if _entry_name in lane_dict:
                                        del lane_dict[_entry_name]
                                    else:
                                        print('Missing file', os.path.join(lane_path, _entry_name))
                                else:
                                    # s_1_1101.clocs
                                    # s_1_2316.clocs
                                    _entry_name = 's_{}_{:d}{:d}{:02d}.clocs'.format(lane, surface, swath, tile)
                                    if _entry_name in lane_dict:
                                        del lane_dict[_entry_name]
                                    else:
                                        print('Missing file', os.path.join(lane_path, _entry_name))

                # Process IRF/Data/Intensities/L00[1-8]/C[0-9]+.1 directories.

                if rta not in (
                        '1.18.54', '1.18.54.4',  # MiSeq
                        '1.18.64',  # HiSeq 1000/1500/2000/2500
                        '2.4.11',  # NextSeq 500/550
                ):
                    # Exclude the MiSeq, NextSeq 500/550 and HiSeq 1000/1500/2000/2500 instruments,
                    # as they do no longer store cycle-specific sub directories with
                    # cluster intensity files (*.cif), error map (*.errorMap) and
                    # full width at half maximum (*.FWHMMap) files.
                    for cycle in range(0 + 1, self.run_information.get_cycle_number + 1):
                        cycle_name = 'C{:d}.1'.format(cycle)
                        cycle_path = os.path.join(lane_path, cycle_name)
                        if cycle_name in lane_dict:
                            del lane_dict[cycle_name]
                        else:
                            print('Missing directory', cycle_path)
                            continue
                        cycle_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                        if debug > 2:
                            print('Processing directory', cycle_path)

                        for surface in range(0 + 1, fcl.surface_count + 1):
                            for swath in range(0 + 1, fcl.swath_count + 1):
                                for tile in range(0 + 1, fcl.tile_count + 1):
                                    # Not all tiles have to exists especially after catastrophic events during the
                                    # cluster generation step.
                                    tile_name = '{:1d}{:1d}{:02d}'.format(surface, swath, tile)
                                    if self._is_missing_image_analysis_tile(lane=lane, tile=tile_name):
                                        continue

                                    # Process *.cif files.
                                    # s_1_1101.cif
                                    # s_1_2316.cif
                                    _entry_name = 's_{}_{:d}{:d}{:02d}.cif'.format(lane, surface, swath, tile)
                                    if _entry_name in cycle_dict:
                                        del cycle_dict[_entry_name]
                                    else:
                                        print('Missing cif file', os.path.join(cycle_path, _entry_name))

                                    if rta in (
                                            '1.12.4', '1.12.4.2', '1.13.48',  # HiSeq 2000
                                    ):
                                        # Older HiSeq 2000 RTA versions store error map (*.errorMap) and
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
                                                print('Missing error map file', os.path.join(cycle_path, _entry_name))

                                        # Process *.FWHMMap files.
                                        # s_1_1101_T.FWHMMap
                                        # s_1_2308_T.FWHMMap
                                        _entry_name = 's_{}_{:d}{:d}{:02d}_T.FWHMMap'. \
                                            format(lane, surface, swath, tile)
                                        if _entry_name in cycle_dict:
                                            del cycle_dict[_entry_name]
                                        else:
                                            print('Missing FWHM map file', os.path.join(cycle_path, _entry_name))

                        if len(cycle_dict):
                            print(cycle_path, 'with number of entries:', str(len(cycle_dict)))
                            print('  Remaining entries:', sorted(cycle_dict))

                if len(lane_dict):
                    print(lane_path, 'with number of entries:', str(len(lane_dict)))
                    print('  Remaining entries:', sorted(lane_dict))

            if rta not in (
                    '2.4.11',  # NextSeq 500/550
            ):
                # Process the IRF/Data/Intensities/Offsets/ directory.
                self._check_data_intensities_offsets(
                    directory_dict=_directory_dict,
                    directory_path=_directory_path,
                    debug=debug)

                # Process the IRF/Data/Intensities/RTAConfiguration.xml file.
                _entry_name = 'RTAConfiguration.xml'
                if _entry_name in _directory_dict:
                    del _directory_dict[_entry_name]
                else:
                    print('Missing Real Time Analysis configuration file', os.path.join(_directory_path, _entry_name))

            if rta in (
                    '1.12.4', '1.12.4.2', '1.13.48',  # HiSeq 2000
            ):
                # Older HiSeq 2000 RTA versions have position (*_pos.txt) files
                # in addition to cluster location (*.clocs) files.
                # s_1_1101_pos.txt
                # s_8_2308_pos.txt
                for lane in range(0 + 1, fcl.lane_count + 1):
                    for surface in range(0 + 1, fcl.surface_count + 1):
                        for swath in range(0 + 1, fcl.swath_count + 1):
                            for tile in range(0 + 1, fcl.tile_count + 1):
                                _entry_name = 's_{}_{:d}{:d}{:02d}_pos.txt'.format(lane, surface, swath, tile)
                                if _entry_name in _directory_dict:
                                    del _directory_dict[_entry_name]
                                else:
                                    print('Missing pos file', os.path.join(_directory_path, _entry_name))

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data_tile_status(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Data/TileStatus/} directory.

        @param directory_dict: Python C{dict} of I{IRF/Data/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: I{IRF/Data/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        fcl = self.run_information.flow_cell_layout

        _directory_name = 'TileStatus'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        for lane in range(0 + 1, fcl.lane_count + 1):
            for surface in range(0 + 1, fcl.surface_count + 1):
                for swath in range(0 + 1, fcl.swath_count + 1):
                    for tile in range(0 + 1, fcl.tile_count + 1):
                        tile_prefix = 'TileStatusL{:d}T{:d}{:d}{:02d}'.format(lane, surface, swath, tile)

                        # Process *.bin files.
                        # TileStatusL001T1101.bin
                        # TileStatusL001T2308.bin
                        _entry_name = tile_prefix + '.bin'
                        if _entry_name in _directory_dict:
                            del _directory_dict[_entry_name]
                        else:
                            print('Missing file', os.path.join(_directory_path, _entry_name))

                        # Process *.tpl files.
                        # TileStatusL001T1101.tpl
                        # TileStatusL001T2308.tpl
                        _entry_name = tile_prefix + '.tpl'
                        if _entry_name in _directory_dict:
                            del _directory_dict[_entry_name]
                        else:
                            print('Missing file', os.path.join(_directory_path, _entry_name))

        if len(directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_data(self, directory_dict, directory_path, debug=0):
        """Check the IRF/Data/ directory.

        @param directory_dict: Ptyhon C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'Data'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        # Check the IRF/Data/Intensities/ directory.

        self._check_data_intensities(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54',  # MiSeq Control Software 2.5.0.5
        ):
            # Only for HiSeq 2000 instruments and MiSeq Control Software 2.5.0.5.
            # Check the IRF/Data/ImageSize.dat file.
            _entry_name = 'ImageSize.dat'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing file', os.path.join(_directory_path, _entry_name))

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            # Only for HiSeq 2000 and MiSeq instruments.
            # Check the IRF/Data/RTALogs directory.
            _entry_name = 'RTALogs'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing directory', os.path.join(_directory_path, _entry_name))

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48',  # HiSeq 2000
        ):
            # HiSeq Control Software 1.4.5, 1.4.8 and 1.5.15.1 (HiSeq 2000)
            # Check the IRF/Data/reports/ directory.
            # NOTE: Should this directory be checked for completeness?
            _entry_name = 'reports'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing directory', os.path.join(_directory_path, _entry_name))

            # Check the IRF/Data/Status_Files/ directory.
            # NOTE: Should this directory be checked for completeness?
            _entry_name = 'Status_Files'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing directory', os.path.join(_directory_path, _entry_name))

            # Check the IRF/Data/Status.htm file.
            _entry_name = 'Status.htm'
            if _entry_name in _directory_dict:
                del _directory_dict[_entry_name]
            else:
                print('Missing file', os.path.join(_directory_path, _entry_name))

        if rta in (
                '1.18.54',  # MiSeq
        ):
            # Process the IRF/Data/TileStatus/ directory.
            self._check_data_tile_status(
                directory_path=_directory_path,
                directory_dict=_directory_dict,
                debug=debug)

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_inter_op(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/InterOp/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        rta = self.run_parameters.get_real_time_analysis_version

        _directory_name = 'InterOp'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        _file_name_list = [
            'CorrectedIntMetricsOut.bin',
            'ErrorMetricsOut.bin',
            'ExtractionMetricsOut.bin',
            'QMetricsOut.bin',
            'TileMetricsOut.bin',
        ]

        if rta in (
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            _file_name_list.append('IndexMetricsOut.bin')

        if rta in (
                '2.4.11',  # NextSeq 500/550
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('EmpiricalPhasingMetricsOut.bin')
            _file_name_list.append('EventMetricsOut.bin')
            _file_name_list.append('PFGridMetricsOut.bin')
            _file_name_list.append('RegistrationMetricsOut.bin')

        if rta in (
                '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
        ):
            # HiSeq 3000/4000 instruments, excluding RTA 2.5.2 version.
            _file_name_list.append('ColorMatrixMetricsOut.bin')
            _file_name_list.append('FWHMGridMetricsOut.bin')
            _file_name_list.append('StaticRunMetricsOut.bin')

        if rta in (
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('AlignmentMetricsOut.bin')
            _file_name_list.append('BasecallingMetricsOut.bin')
            _file_name_list.append('ExtendedTileMetricsOut.bin')
            _file_name_list.append('FWHMGridMetricsOut.bin')
            _file_name_list.append('OpticalModelMetricsOut.bin')
            _file_name_list.append('QMetrics2030Out.bin')
            _file_name_list.append('QMetricsByLaneOut.bin')

            read_start_list = self.run_information.get_read_start_list
            read_end_list = self.run_information.get_read_end_list

            # Qualities are calculated for payload i.e. non-index reads from cycle 25 onwards.
            cycle_number: int = 1
            quality_cycle_list: List[int] = list()
            quality_start_list: List[int] = list()
            quality_cycle_start_read_1 = 0
            for run_information_read in self.run_information.run_information_read_list:
                if not run_information_read.index:
                    # The quality_cycle is the minimum of cycle 25 and the read length minus one.
                    quality_cycle_offset = min(24, run_information_read.cycles - 1)
                    quality_cycle_list.extend(
                        range(
                            cycle_number + quality_cycle_offset,
                            cycle_number + run_information_read.cycles))
                    quality_start_list.append(cycle_number + quality_cycle_offset)
                    if run_information_read.number == 1:
                        # Capture the quality cycle offset for the first non-index read.
                        quality_cycle_start_read_1 = quality_cycle_offset + 1
                cycle_number += run_information_read.cycles

            for cycle in range(0 + 1, self.run_information.get_cycle_number + 1):
                cycle_name = 'C{:d}.1'.format(cycle)
                cycle_path = os.path.join(_directory_path, cycle_name)

                if cycle_name in _directory_dict:
                    del _directory_dict[cycle_name]
                else:
                    print('Missing directory', cycle_path)
                    continue

                cycle_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))
                _cycle_file_name_list = [
                    'BasecallingMetricsOut.bin',
                    'EventMetricsOut.bin',
                    'ExtractionMetricsOut.bin',
                    'ImageMetricsOut.bin',
                    'RegistrationMetricsOut.bin',
                ]

                if cycle in read_start_list:
                    _cycle_file_name_list.append('FWHMGridMetricsOut.bin')

                if cycle in read_end_list:
                    _cycle_file_name_list.append('FWHMGridMetricsOut.bin')
                else:
                    _cycle_file_name_list.append('EmpiricalPhasingMetricsOut.bin')

                if cycle in quality_cycle_list:
                    _cycle_file_name_list.append('AlignmentMetricsOut.bin')
                    _cycle_file_name_list.append('ErrorMetricsOut.bin')

                if cycle in quality_start_list:
                    _cycle_file_name_list.append('TileMetricsOut.bin')

                if cycle == 1:
                    _cycle_file_name_list.append('OpticalModelMetricsOut.bin')

                if cycle == 10:
                    _cycle_file_name_list.append('ExtendedTileMetricsOut.bin')

                if cycle == quality_cycle_start_read_1:
                    # The following files appear at either 25 or the end or read 1, which ever is shorter.
                    _cycle_file_name_list.append('CorrectedIntMetricsOut.bin')
                    _cycle_file_name_list.append('QMetricsOut.bin')
                    _cycle_file_name_list.append('PFGridMetricsOut.bin')

                if cycle > 25:
                    # Irrespective of read 1 length, the following files appear after cycle 25.
                    _cycle_file_name_list.append('CorrectedIntMetricsOut.bin')
                    _cycle_file_name_list.append('QMetricsOut.bin')

                self._check_file_names(
                    directory_dict=cycle_dict,
                    directory_path=cycle_path,
                    file_name_list=_cycle_file_name_list)

                if len(cycle_dict):
                    print(cycle_path, 'with number of entries:', str(len(cycle_dict)))
                    print('  Remaining entries:', sorted(cycle_dict))

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('ImageMetricsOut.bin')

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            _file_name_list.append('ControlMetricsOut.bin')

        self._check_file_names(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_periodic_save_rates(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/PeriodicSaveRates/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        _directory_name = 'PeriodicSaveRates'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        _file_name_list = [
            'Save All Thumbnails.xml'
        ]

        self._check_file_names(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_recipe(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Recipe/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        rta = self.run_parameters.get_real_time_analysis_version
        flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.upper()

        _directory_name = 'Recipe'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        _file_name_list: List[str] = list()

        if rta in (
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            # Process the IRF/Recipe/<ReagentKitRFIDTag/SerialNumber>.xml file.
            _file_name_list.append(
                self.run_parameters.element_tree.find(path='ReagentKitRFIDTag/SerialNumber').text + '.xml')

            # Process the IRF/Recipe/RunState.xml file.
            _file_name_list.append('RunState.xml')
        elif rta in (
                '2.4.11',  # NextSeq 500/550
        ):
            # Process the IRF/Recipe/<ReagentKitSerial>.xml file.
            _file_name_list.append(
                self.run_parameters.element_tree.find(path='ReagentKitSerial').text + '.xml')
        else:
            # Process the IRF/Recipe/<FCID>.xml file.
            _file_name_list.append(flow_cell_barcode + '.xml')

            if rta in (
                    '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                    '1.18.54', '1.18.54.4',  # MiSeq
                    '2.4.11',  # NextSeq 500/550
            ):
                # Process the IRF/Recipe/<FCID>_RunState.xml file.
                _file_name_list.append(flow_cell_barcode + '_RunState.xml')

        self._check_file_names(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def _check_thumbnail_images(self, directory_dict, directory_path, debug=0):
        """Check the I{IRF/Thumbnail_Images/} directory.

        @param directory_dict: Python C{dict} of Illumina Run Folder I{IRF/} entries
        @type directory_dict: dict[str, int]
        @param directory_path: Illumina Run Folder I{IRF/} path
        @type directory_path: str
        @param debug: Integer debugging level
        @type debug: int
        """
        fcl = self.run_information.flow_cell_layout
        rta = self.run_parameters.get_real_time_analysis_version

        if rta in (
                '2.4.11',  # NextSeq 500/550
        ):
            return

        flow_cell_barcode = self.run_parameters.get_flow_cell_barcode.lower()

        read_start_list = self.run_information.get_read_start_list
        # read_end_list = self.run_information.get_read_end_list

        # Helper dict to map surface numbers to abbreviations.

        surface_dict = {1: 'bot', 2: 'top'}

        # Process the IRF/Thumbnail_Images/ directory.

        _directory_name = 'Thumbnail_Images'
        _directory_path = os.path.join(directory_path, _directory_name)

        if _directory_name in directory_dict:
            del directory_dict[_directory_name]
        else:
            print('Missing directory', _directory_path)
            return

        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        # Process the IRF/Thumbnail_Images/L00[1-8]/ directories.

        for lane in range(0 + 1, fcl.lane_count + 1):
            lane_name = 'L{:03d}'.format(lane)

            if lane_name in _directory_dict:
                del _directory_dict[lane_name]
            else:
                print('Missing directory', os.path.join(_directory_path, lane_name))
                continue

            lane_path = os.path.join(_directory_path, lane_name)
            lane_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(lane_path)))

            if debug > 1:
                print('Processing directory', lane_path)

            # Process the IRF/Thumbnail_Images/L00[1-8]/C[0-9]+.1/ directories.

            for cycle in range(0 + 1, self.run_information.get_cycle_number + 1):
                if rta in ('v3.4.4',) and cycle not in read_start_list:
                    # Since RTA v3.4.4, only the first cycles of each read have thumbnail images.
                    continue

                cycle_name = 'C{:d}.1'.format(cycle)
                cycle_path = os.path.join(lane_path, cycle_name)

                if cycle_name in lane_dict:
                    del lane_dict[cycle_name]
                else:
                    print('Missing directory', cycle_path)
                    continue

                cycle_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(cycle_path)))

                if debug > 2:
                    print('Processing directory', cycle_path)

                for surface in range(0 + 1, fcl.surface_count + 1):
                    for swath in range(0 + 1, fcl.swath_count + 1):
                        if rta in (
                                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
                        ):
                            # NovaSeq 6000 instruments have green and red bases.
                            for base in ('green', 'red'):
                                for tile in range(0 + 1, fcl.tile_count + 1):
                                    # FIXME: This should work off the list of tiles configured in RunInfo.xml.
                                    # NovaSeq 6000 only stores thumbnails for tiles that end in 3. Strange.
                                    # s_2_1103_green.png
                                    # s_2_1103_red.png
                                    if (tile - 3) % 10:
                                        continue

                                    # NovaSeq SP flow cells have only particular tiles from the full range defined.
                                    if not self.run_information.has_tile(tile='{:1d}_{:1d}{:1d}{:02d}'.format(
                                            lane, surface, swath, tile)):
                                        continue

                                    tile_file = 's_{:1d}_{:1d}{:1d}{:02d}_{}.png'.format(
                                        lane, surface, swath, tile, base)
                                    if tile_file in cycle_dict:
                                        del cycle_dict[tile_file]
                                    else:
                                        print('Missing tile file', os.path.join(cycle_path, tile_file))
                        else:
                            # NextSeq 500/550 does not have a IRF/Thumbnail_Images directory,
                            # all other instruments have A C G T bases.
                            for base in ('a', 'c', 'g', 't'):
                                # Process swath image and *.jpg.zprof files.
                                if rta in (
                                        '1.18.54', '1.18.54.4',  # MiSeq
                                ):
                                    # The MiSeq instrument does not have swath image and zprof files.
                                    pass
                                else:
                                    # c6nk1anxx_c001_l1_t001_bot_s1_a.jpg
                                    # c6nk1anxx_c001_l1_t001_bot_s1_a.jpg.zprof
                                    # c6nk1anxx_c001_l1_t001_top_s3_t.jpg
                                    # c6nk1anxx_c001_l1_t001_top_s3_t.jpg.zprof
                                    _entry_name = '{}_c{:03d}_l{:d}_t{:03d}_{}_s{}_{}.jpg'.format(
                                        flow_cell_barcode, cycle, lane, 1, surface_dict[surface], swath, base)
                                    if _entry_name in cycle_dict:
                                        del cycle_dict[_entry_name]
                                    else:
                                        print('Missing swath image file', os.path.join(cycle_path, _entry_name))

                                    _entry_name += '.zprof'
                                    if base in ('c', 'g', 't') and rta in (
                                            '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                                    ):
                                        # The HiSeq 3000/4000 instrument does not have swath files for bases c, g and t.
                                        pass
                                    else:
                                        if _entry_name in cycle_dict:
                                            del cycle_dict[_entry_name]
                                        else:
                                            print('Missing swath zprof file', os.path.join(cycle_path, _entry_name))

                                # Process tile image files.
                                if rta in (
                                        '1.18.54', '1.18.54.4',  # MiSeq
                                        '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
                                ):
                                    # The HiSeq 3000/4000 and MiSeq instruments use lower case bases.
                                    pass
                                else:
                                    # The HiSeq 2000 instruments uses upper case bases.
                                    base = base.upper()

                                for tile in range(0 + 1, fcl.tile_count + 1):
                                    # s_1_1101_A.jpg
                                    # s_1_2316_T.jpg
                                    tile_file = 's_{:1d}_{:1d}{:1d}{:02d}_{}.jpg'.format(
                                        lane, surface, swath, tile, base)
                                    if tile_file in cycle_dict:
                                        del cycle_dict[tile_file]
                                    else:
                                        print('Missing tile file', os.path.join(cycle_path, tile_file))

                if len(cycle_dict):
                    print(cycle_path, 'with number of entries:', str(len(cycle_dict)))
                    print('  Remaining entries:', sorted(cycle_dict))

            if len(lane_dict):
                print(lane_path, 'with number of entries:', str(len(lane_dict)))
                print('  Remaining entries:', sorted(lane_dict))

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return

    def check(self, debug=0):
        """Check an Illumina Run Folder with regards to its internal directory and file structure.

        Both, missing and additional files are printed to STDOUT.
        @param debug: Integer debugging level
        @type debug: int
        """
        rta = self.run_parameters.get_real_time_analysis_version

        if rta not in self.rta_dict:
            raise Exception('Unsupported RTA version: ' + repr(rta))

        # _directory_name = os.path.basename(self.file_path)
        _directory_path = self.file_path
        _directory_dict: Dict[str, int] = dict(map(lambda x: (x, 1), os.listdir(_directory_path)))

        if debug > 0:
            print('Processing directory', _directory_path)

        # Check the IRF/Config directory.

        self._check_config(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            debug=debug)

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

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
        ):
            # Check the IRF/PeriodicSaveRates/ directory.

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

        # Check IRF/ files.

        _file_name_list = [
            'Logs',  # directory
            'RTAComplete.txt',
            'RunInfo.xml',
        ]

        if rta in (
                '1.18.54.4',  # MiSeq Control Software 4.0.0.1769 (MiSeq)
                '2.4.11',  # NextSeq 500/550
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('RunParameters.xml')
        else:
            _file_name_list.append('runParameters.xml')

        if rta in (
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            _file_name_list.append('SampleSheet.csv')

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
        ):
            _file_name_list.append('First_Base_Report.htm')

        if rta in (
                '2.4.11',  # NextSeq 500/550
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
        ):
            _file_name_list.append('RTAConfiguration.xml')
            _file_name_list.append('RTALogs')  # directory
            for read_number in range(0 + 1, len(self.run_information.run_information_read_list) + 1):
                _file_name_list.append('RTARead{:d}Complete.txt'.format(read_number))

        if rta in (
                '2.5.2', '2.7.3', '2.7.6', '2.7.7',  # HiSeq 3000/4000
        ):
            _file_name_list.append('SequencingComplete.txt')

        if rta in (
                '1.18.54.4',  # MiSeq Control Software 4.0.0.1769 (MiSeq)
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('CopyComplete.txt')

        if rta in (
                'v3.3.3', 'v3.4.4',  # NovaSeq 6000
        ):
            _file_name_list.append('RTA3.cfg')
            # _file_name_list.append('RunComplete.txt')
            _file_name_list.append('SequenceComplete.txt')

        if rta in (
                '1.12.4', '1.12.4.2', '1.13.48', '1.17.21.3', '1.18.61', '1.18.64',  # HiSeq 2000
                '1.18.54', '1.18.54.4',  # MiSeq
        ):
            _file_name_list.append('Basecalling_Netcopy_complete.txt')
            _file_name_list.append('ImageAnalysis_Netcopy_complete.txt')
            for read_number in range(0 + 1, len(self.run_information.run_information_read_list) + 1):
                _file_name_list.append('Basecalling_Netcopy_complete_Read{:d}.txt'.format(read_number))
                _file_name_list.append('ImageAnalysis_Netcopy_complete_Read{:d}.txt'.format(read_number))

        if rta in (
                '1.18.54.4',  # MiSeq Control Software 4.0.0.1769 (MiSeq)
                '2.4.11',  # NextSeq 500/550
        ):
            _file_name_list.append('RunCompletionStatus.xml')

        if rta in (
                '1.18.54.4',  # MiSeq Control Software 4.0.0.1769 (MiSeq)
        ):
            _file_name_list.extend(
                (
                    'Alignment_1',
                    'CompletedJobInfo.xml',
                    'GenerateFASTQRunStatistics.xml',
                    'InstrumentAnalyticsLogs',
                    'QueuedForAnalysis.txt',
                    'ReportInfo.dat',
                    'RunCheckDetail.txt',
                    'SoftwareVersionsFile.csv',
                )
            )

        self._check_file_names(
            directory_dict=_directory_dict,
            directory_path=_directory_path,
            file_name_list=_file_name_list,
            debug=debug)

        if len(_directory_dict):
            print(_directory_path, 'with number of entries:', str(len(_directory_dict)))
            print('  Remaining entries:', sorted(_directory_dict))

        return
