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


class RunInformationRead(object):

    """BSF Illumina Run Information Reads class.

    The Bio.BSF.Data.Illumina.RunInformationRead class models
    one <Read> element in a BSF Illumina Run Information (RunInfo.xml) document.
    """

    def __init__(self, number=0, cycles=0, index=False):

        """Initialise a Bio.BSF.Data.Illumina.RunInformationRead object.

        :param self: Bio.BSF.Data.Illumina.RunInformationRead object
        :type self: RunInformationRead
        :param number: Read number
        :type number: int
        :param cycles: Cycle number
        :type cycles: int
        :param index: Index read
        :type index: bool
        :return: Nothing
        :rtype: None
        """

        self.number = number
        self.cycles = cycles
        self.index = index


class RunInformation(object):

    """BSF Illumina Run Information class.

    The Bio.BSF.Data.Illumina.RunInformation class represents an Illumina
    run information XML (RunInfo.xml) document.

    Attributes:
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type
                     CASAVA: FASTQ file after post-processing with CASAVA.
                     External: other data files.
    :type file_type: str
    :ivar name: Name
    :type name: str
    :ivar run_identifier: Run identifier e.g. 130724_SN815_0089_BC26JBACXX
    :type run_identifier: str
    :ivar run_number: Run number, which may not have to correspond to the run number in the run identifier e.g. 91
    :type run_number: str
    :ivar flow_cell: Illumina Flow-Cell identifier e.g. C26JBACXX
    :type flow_cell: str
    :ivar instrument: Illumina instrument serial number e.g. SN815
    :type instrument: str
    :ivar date: Date in YYMMDD format e.g. 130724
    :type date: str
    :ivar reads: Python list of Bio.BSF.Data.Illumina.RunInformationRead objects
    :type reads: list
    """

    @classmethod
    def from_file_path(cls, file_path):

        """Create a Bio.BSF.Data.Illumina.RunInformation object from a file path.

        :param cls: Class
        :type cls: RunInformation
        :param file_path: File path
        :type file_path: str, unicode
        :return: BSF RunInformation object
        :rtype: RunInformation
        """

        file_name = os.path.basename(file_path)

        run_info_tree = ET.parse(source=file_path)
        run_info_root = run_info_tree.getroot()

        # Parse meta-information about the Illumina run

        run_identifier = run_info_root.find('Run').attrib['Id']  # e.g. 130724_SN815_0089_BC26JBACXX
        run_number = run_info_root.find('Run').attrib['Number']  # e.g. 91
        flow_cell = run_info_root.find('Run/Flowcell').text  # e.g. C26JBACXX
        instrument = run_info_root.find('Run/Instrument').text  # e.g. SN815
        date = run_info_root.find('Run/Date').text  # e.g. 130724

        xml_list_reads = run_info_tree.find("Run/Reads")

        reads = list()

        for xml_read in xml_list_reads:

            is_index = bool(0)

            if xml_read.attrib['IsIndexedRead'] == 'Y':
                is_index = True
            elif xml_read.attrib['IsIndexedRead'] == 'N':
                is_index = False
            else:
                warning = 'Unexpected value {} in Read element attribute IsIndexedRead '. \
                    format(xml_read.attrib['IsIndexedRead'])
                warnings.warn(warning)

            reads.append(RunInformationRead(number=int(xml_read.attrib['Number']),
                                            cycles=int(xml_read.attrib['NumCycles']),
                                            index=is_index))

        reads.sort(key=lambda read: read.number)

        # Warn if there is not at least one non-index read.

        non_index_reads = filter(lambda read: not read.index, reads)
        if len(non_index_reads) == 0:
            message = 'No non-index read in Illumina RunInformation {}'.format(file_path)
            warnings.warn(message, UserWarning)

        # Set a paired_end attribute if more than one read without index is defined?

        iri = cls(file_path=file_path, file_type='xml', name=file_name,
                  run_identifier=run_identifier, run_number=run_number,
                  flow_cell=flow_cell, instrument=instrument, date=date,
                  reads=reads)

        return iri

    def __init__(self, file_path=None, file_type=None, name=None,
                 run_identifier=None, run_number=None, flow_cell=None, instrument=None, date=None,
                 reads=None):

        """Initialise a Bio.BSF.Data.Illumina.RunInformation object.

        :param self: Bio.BSF.Data.Illumina.RunInformation
        :type self: RunInformation
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External or Automatic)
        :type file_type: str
        :param name: Name
        :type name: str
        :param run_identifier: Run identifier e.g. 130724_SN815_0089_BC26JBACXX
        :type run_identifier: str
        :param run_number: Run number, which may not have to correspond to the run number in the run identifier e.g. 91
        :type run_number: str
        :param flow_cell: Illumina Flow-Cell identifier
        :type flow_cell: str
        :param instrument: Illumina instrument serial number
        :type instrument: str
        :param date: Date in YYMMDD format
        :type date: str
        :param reads: Python list of Bio.BSF.Data.Illumina.RunInformationRead objects
        :type reads: list
        :return: Nothing
        :rtype: None
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

    def picard_read_structure(self):

        """Get the read structure for the Picard ExtractIlluminaBarcodes module.

        Codes:
        T ... Template
        B ... Barcode
        S ... Skip
        :param self: Bio.BSF.Data.Illumina.RunInformation object
        :type self: RunInformation
        :return: Read structure string for Picard ExtractIlluminaBarcodes
        :rtype: str
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

    """BSF Illumina Run Parameters class.

    Attributes:
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar element_tree: xml.etree.ElementTree
    :type element_tree: ElementTree
    """

    @classmethod
    def from_file_path(cls, file_path):

        """Create a Bio.BSF.Data.Illumina.RunParameters object from a file path.

        :param cls: Class
        :type cls: RunParameters
        :param file_path: File path
        :type file_path: str, unicode
        :return: BSF RunParameters object
        :rtype: RunParameters
        """

        # file_name = os.path.basename(file_path)

        return RunParameters(file_path=file_path, element_tree=ElementTree(file=file_path))

    def __init__(self, file_path=None, element_tree=None):

        """Initialise a BSF Illumina Run Parameters object.
        :param self: BSF Illumina Run Parameters
        :type self: RunParameters
        :param file_path: File path
        :type file_path: str, unicode
        :param element_tree: XML Element Tree
        :type element_tree: ElementTree
        :return: Nothing
        :rtype: None
        """

        if file_path:
            self.file_path = file_path
        else:
            self.file_path = str()

        if element_tree:
            self.element_tree = element_tree
        else:
            self.element_tree = ElementTree()

    def get_experiment_name(self):

        """Get the experiment name of a BSF Illumina Run Parameters object.
        :param self: BSF Run Parameters
        :type self: RunParameters
        :return: Experiment name e.g BSF_0001
        :rtype: str
        """

        return self.element_tree.getroot().find('ExperimentName').text

    def get_flow_cell_barcode(self):

        """Get the flow-cell barcode of a BSF Illumina Run Parameters object.
        :param self: BSF Run Parameters
        :type self: RunParameters
        :return: Flow-cell barcode e.g BSF_0001
        :rtype: str
        """

        return self.element_tree.getroot().find('Barcode').text

    def get_flow_cell_type(self):

        """Get the flow cell of a BSF Illumina Run Parameters object.
        :param self: BSF Run Parameters
        :type self: RunParameters
        :return: Flow-cell type
        :rtype: str
        """

        return self.element_tree.getroot().find('Flowcell').text

    def get_position(self):

        """Get the flow-cell position of a BSF Illumina Run Parameters object.
        :param self: BSF Run Parameters
        :type self: RunParameters
        :return: Flow-cell position e.g. A or B
        :rtype: str
        """

        return self.element_tree.getroot().find('FCPosition').text

    def get_run_identifier(self):

        """Get the run identifier of a BSF Illumina Run Parameters object.
        :param self: BSF Run Parameters
        :type self: RunParameters
        :return: Run identifier
        :rtype: str
        """

        return self.element_tree.getroot().find('RunID').text


class RunFolder(object):

    """BSF Illumina Run Folder class.

    The BSF RunFolder class represents an Illumina
    run folder copied off the instrument.

    Attributes:
    :ivar file_path: File path
    :type file_path: str, unicode
    :ivar file_type: File type
                     CASAVA: FASTQ file after post-processing with CASAVA.
                     External: other data files.
    :type file_type: str
    :ivar name: Name
    :type name: str
    :ivar date: Date in YYMMDD format
    :type date: str
    :ivar instrument: Illumina instrument serial number
    :type instrument: str
    :ivar run: Run serial number
    :type run: str
    :ivar flow_cell: Flow-cell identifier.
    :type flow_cell: str
    :ivar run_information: Bio.BSF.Data.Illumina.RunInformation object
    :type run_information: RunInformation
    """

    @classmethod
    def from_file_path(cls, file_path):

        """Construct a Bio.BSF.Data.Illumina.RunFolder object from a file path.

        :param cls: Class
        :type cls: Class
        :param file_path: File path
        :type file_path: str, unicode
        :return: Bio.BSF.Data.Illumina.RunFolder object
        :rtype: RunFolder
        """

        # Since os.path.basename returns an empty string at trailing slashes,
        # use string.rstrip to remove them.

        file_name = os.path.basename(file_path.rstrip('/ '))

        # Illumina Run Folders obey a "YYMMDD_SN000_Run_PositionFCID"
        # schema.

        components = file_name.split('_')

        irf = cls(file_path=file_path, file_type='Illumina', name=string.join(components[:], '_'),
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

        """Initialise a Bio.BSF.Data.Illumina.RunFolder object.

        :param self: Bio.BSF.Data.Illumina.RunFolder
        :type self: RunFolder
        :param file_path: File path
        :type file_path: str, unicode
        :param file_type: File type (e.g. CASAVA, External or Automatic)
        :type file_type: str
        :param name: Name
        :type name: str
        :param date: Date in YYMMDD format
        :type date: str
        :param instrument: Illumina instrument serial number (e.g. SN181, SN815, ...)
        :type instrument: str
        :param run: Run serial number
        :type run: str
        :param flow_cell: The position and flow cell identifier
        :type flow_cell: str
        :param run_information: Bio.BSF.Illumina.RunInformation
        :type run_information: RunInformation
        :return: Nothing
        :rtype: None
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

    def get_base_calls_directory(self):

        """Get the Base-calls directory in the RunFolder/Data/Intensities/BaseCalls hierarchy.

        :param self: Illumina Run Folder
        :type self:RunFolder
        :return: Illumina Base-calls Directory
        :rtype: str, unicode
        """

        return os.path.join(self.file_path, 'Data',  'Intensities', 'BaseCalls')
