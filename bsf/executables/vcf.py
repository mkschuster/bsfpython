#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2013 - 2022 Michael K. Schuster
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
"""The :py:mod:`bsf.executables.vcf` module provides classes and functions to split
Ensembl VEP annotation in the :literal:`CSQ` field of Variant Calling Format (VCF)
`output <http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout>`_ into :literal:`VEP_*` fields.
"""
import logging
from argparse import ArgumentParser
from csv import DictReader
from subprocess import Popen
from typing import Dict, List, Optional, Tuple

from pysam import VariantFile

from bsf.connector import Connector
from bsf.process import Command, Executable, RunnableStep

module_logger = logging.getLogger(name=__name__)


class RunnableStepCsqToVep(RunnableStep):
    """The :py:class:`bsf.executables.vcf.RunnableStepCsqToVep` class expands Ensembl Variant Effect Predictor (VEP)
    :literal:`CSQ INFO` annotation into a set of :literal:`VEP_* INFO` annotation.

    :ivar soc_path: A Sequence Ontology term priority configuration (TSV) file path.
    :type soc_path: str | None
    :ivar ofc_path: An output field configuration (TSV) file path.
    :type ofc_path: str | None
    :ivar vcf_path_old: An old VCF file path.
    :type vcf_path_old: str | None
    :ivar vcf_path_new: A new VCF file path.
    :type vcf_path_new: str | None
    """

    def __init__(
            self,
            name,
            program=None,
            options=None,
            arguments=None,
            sub_command=None,
            stdin=None,
            stdout=None,
            stderr=None,
            dependencies=None,
            hold=None,
            submit=True,
            process_identifier=None,
            process_name=None,
            sub_process=None,
            obsolete_file_path_list=None,
            soc_path=None,
            ofc_path=None,
            vcf_path_old=None,
            vcf_path_new=None):
        """Initialise a :py:class:`bsf.executables.vcf.RunnableStepCsqToVep` object.

        :param name: A name.
        :type name: str | None
        :param program: A program.
        :type program: str | None
        :param options: A Python :py:class:`dict` object of
            Python :py:class:`str` (:py:attr:`bsf.argument.Argument.key`) key and
            Python :py:class:`list` value objects of :py:class:`bsf.argument.Argument` objects.
        :type options: dict[str, list[Argument]] | None
        :param arguments: A Python :py:class:`list` object of Python :py:class:`str` (program argument) objects.
        :type arguments: list[str] | None
        :param sub_command: A subordinate :py:class:`bsf.process.Command` object.
        :type sub_command: Command | None
        :param stdin: A standard input :literal:`STDIN` :py:class:`bsf.connector.Connector` object.
        :type stdin: Connector | None
        :param stdout: A standard output :literal:`STDOUT` :py:class:`bsf.connector.Connector` object.
        :type stdout: Connector | None
        :param stderr: A standard error :literal:`STDERR` :py:class:`bsf.connector.Connector` object.
        :type stderr: Connector | None
        :param dependencies: A Python :py:class:`list` object of
            Python :py:class:`str` (:py:attr:`bsf.process.Executable.name`) objects
            in the context of :py:class:`bsf.analysis.Stage` dependencies.
        :type dependencies: list[str] | None
        :param hold: Request a hold on job scheduling.
        :type hold: bool | None
        :param submit: Request the submission via the :py:meth:`bsf.analysis.Stage.submit` method.
        :type submit: bool
        :param process_identifier: A process identifier.
        :type process_identifier: str | None
        :param process_name: A process name.
        :type process_name: str | None
        :param sub_process: A :py:class:`subprocess.Popen` object.
        :type sub_process: Popen | None
        :param obsolete_file_path_list: A Python :py:class:`list` object of
            Python :py:class:`str` (file path) objects
            that can be removed after successfully completing the :py:meth:`bsf.process.RunnableStep.run` method.
        :type obsolete_file_path_list: list[str] | None
        :param soc_path: A Sequence Ontology term priority configuration (TSV) file path.
        :type soc_path: str | None
        :param ofc_path: An output field configuration (TSV) file path.
        :type ofc_path: str | None
        :param vcf_path_old: An old VCF file path.
        :type vcf_path_old: str | None
        :param vcf_path_new: A new VCF file path.
        :type vcf_path_new: str | None
        """
        super(RunnableStepCsqToVep, self).__init__(
            name=name,
            program=program,
            options=options,
            arguments=arguments,
            sub_command=sub_command,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            dependencies=dependencies,
            hold=hold,
            submit=submit,
            process_identifier=process_identifier,
            process_name=process_name,
            sub_process=sub_process,
            obsolete_file_path_list=obsolete_file_path_list)

        self.soc_path = soc_path
        self.ofc_path = ofc_path
        self.vcf_path_old = vcf_path_old
        self.vcf_path_new = vcf_path_new

        return

    def run(self):
        """Run a :py:class:`bsf.executables.vcf.RunnableStepCsqToVep` object.

        :return: A Python :py:class:`list` object of Python :py:class:`str` (exception) objects.
        :rtype: list[str] | None
        """

        def initialise_allele_list(length):
            """Initialise an empty Python :py:class:`list` object of Python :py:class:`list` object of
            Python :py:class:`str` objects.

            The Python :py:class:`list` object stores one Python :py:class:`list` object per allele and picked feature.

            :param length: A length.
            :type length: int
            :return: A Python :py:class:`list` object of Python :py:class:`list` objects of
                Python :py:class:`str` objects
            :rtype: list[list[str]]
            """
            allele_list = list()

            j = 0
            while j < length:
                allele_list.append(list())
                j += 1

            return allele_list

        def collate_allele_values(allele_list):
            """Collate the list of allele-specific values by selecting only the first list component.

            :param allele_list: A Python :py:class:`list` object of
                Python :py:class:`str` (allele-specific value) objects.
            :type allele_list: list[str]
            :return: A first value.
            :rtype: str | None
            """
            if allele_list:
                return allele_list[0]
            else:
                return None

        def get_consequence_index(consequence):
            """Get the SO index.

            :param consequence: A consequence.
            :type consequence: str
            :return: A consequence index.
            :rtype: int
            """
            # module_logger.debug('get_consequence_index(consequence=%r)', consequence)
            for _so_index, _so_string in enumerate(sequence_ontology_list):
                # module_logger.debug('get_consequence_index: consequence: %r _so_index: %r _so_string: %r',
                #                     consequence, _so_index, _so_string)
                if consequence == _so_string:
                    return _so_index

        # Read the Sequence Ontology (TSV) configuration file.
        # This file provides consequence prioritisation.
        sequence_ontology_list: List[str] = list()

        with open(file=self.soc_path, mode='rt') as input_text_io:
            for row_dict in DictReader(f=input_text_io, dialect='excel-tab'):
                sequence_ontology_list.append(row_dict['SO term'])

        # Load the VEP output field (TSV) configuration file.
        vep_header_dict: Dict[str, Dict[str, str]] = dict()

        with open(file=self.ofc_path, mode='rt') as input_text_io:
            for row_dict in DictReader(f=input_text_io, dialect='excel-tab'):
                # NOTE: Override all types with 'String' to allow multiple values joined by '&' characters.
                row_dict['Type'] = 'String'
                vep_header_dict[row_dict['ID']] = row_dict

        variant_file_old = VariantFile(filename=self.vcf_path_old, mode='r')

        # Copy the header for a new VCF instance.

        variant_header_new = variant_file_old.header.copy()

        if 'CSQ' not in variant_file_old.header.info:
            raise Exception("Cannot convert a VCF file without a 'CSQ' INFO field.")

        # Remove the CSQ info field from the new header.

        # NOTE: Apparently, pysam does not allow deleting an entry from the variant header at the moment.
        # The pysam developers commented out calling the bcf_hdr_sync() function after calling bcf_hdr_remove() in
        # VariantHeaderMetadata.remove_header().
        # https://github.com/pysam-developers/pysam/blob/16ecfee0c12908f5b42ff0542d75666475ca4d1a/pysam/libcbcf.pyx#L1616
        #
        # variant_header_new.info.remove_header('CSQ')

        csq_vmd = variant_file_old.header.info['CSQ']

        # Build a new header with VEP_ INFO fields.
        # Split the description on white space first, the pipe-separated list of VEP field declarations
        # is then in the last block.

        csq_key_list: List[str] = csq_vmd.description.split()[-1].split('|')
        csq_index_allele_num: Optional[int] = None

        for csq_index, csq_key in enumerate(csq_key_list):
            if csq_key == 'ALLELE_NUM':
                csq_index_allele_num = csq_index

            if csq_key in vep_header_dict:
                variant_header_new.info.add(
                    id='VEP_' + csq_key,
                    number=vep_header_dict[csq_key]['Number'],
                    type=vep_header_dict[csq_key]['Type'],
                    description=vep_header_dict[csq_key]['Description'])
                continue

            # Split and expand the 'AFFECTED' and 'TOTAL' fields.
            # INTRON -> INTRON_AFFECTED and INTRON_TOTAL
            # EXON   -> EXON_AFFECTED and EXON_TOTAL
            if csq_key in ('INTRON', 'EXON'):
                for suffix in ('AFFECTED', 'TOTAL'):
                    new_field = '_'.join((csq_key, suffix))
                    variant_header_new.info.add(
                        id='_'.join(('VEP', csq_key, suffix)),
                        number=vep_header_dict[new_field]['Number'],
                        type=vep_header_dict[new_field]['Type'],
                        description=vep_header_dict[new_field]['Description'])
                continue

            # Split and expand the 'start' and 'end' fields.
            # cDNA_position    -> cDNA_position_start and cDNA_position_end
            # CDS_position     -> CDS_position_start and CDS_position_end
            # Protein_position -> Protein_position_start and Protein_position_end
            if csq_key in ('cDNA_position', 'CDS_position', 'Protein_position'):
                for suffix in ('start', 'end'):
                    new_field = '_'.join((csq_key, suffix))
                    variant_header_new.info.add(
                        id='_'.join(('VEP', csq_key, suffix)),
                        number=vep_header_dict[new_field]['Number'],
                        type=vep_header_dict[new_field]['Type'],
                        description=vep_header_dict[new_field]['Description'])
                continue

            module_logger.warning('CSQ VEP field %r missing from field configuration file.', csq_key)

        #     # Print the new header.
        #     for key, vmd in variant_header_new.info.items():
        #         # module_logger.debug('Key: %r VariantMetadata: %r', key, vmd)
        #         # module_logger.debug('VariantMetadata: type: %r, dir: %r', type(vmd), dir(vmd))
        #         # The ID property (vmd.id) is not exposed by pysam.
        #         module_logger.debug(
        #             'Key: %r Value: %r Name: %r Number: %r Type: %r Description: %r',
        #             key, vmd, vmd.name, vmd.number, vmd.type, vmd.description)

        # Open the new VariantFile for writing.

        variant_file_new = VariantFile(filename=self.vcf_path_new, mode='w', header=variant_header_new)

        # Parse the CSQ INFO field of each VariantRecord.

        for variant_record in variant_file_old.fetch():
            variant_record.translate(dst_header=variant_header_new)
            variant_record_info = variant_record.info

            # Initialise a Python dict of VEP key and Python list value data.
            # {
            #   VEP_FIELD 1:
            #     [
            #        # ALLELE 1
            #        [PICK 1 annotation, PICK 2 annotation, ...,],
            #        # ALLELE 2
            #        [PICK 1 annotation, PICK 2 annotation, ...,],
            #        ...,
            #     ],
            #     ...,
            # }
            vep_dict: Dict[str, List[List[str]]] = dict()

            # Total number of alleles including REF and ALT.
            allele_length = len(variant_record.alleles)

            if 'CSQ' not in variant_record_info:
                variant_file_new.write(record=variant_record)
                continue

            # Iterate over all comma-separated allele-transcript blocks.
            csq_record_info: str
            for csq_record_info in variant_record_info['CSQ']:
                csq_value_list = csq_record_info.split('|')

                allele_number = int(csq_value_list[csq_index_allele_num])

                for i in range(0, len(csq_key_list)):
                    if csq_key_list[i] in ('INTRON', 'EXON'):
                        # Splitting of INTRON and EXON keys.
                        if csq_value_list[i]:
                            value_list = csq_value_list[i].split('/')
                            value_affected = value_list[0]
                            # Make sure value_total is defined.
                            if len(value_list) == 2:
                                value_total = value_list[1]
                            else:
                                value_total = value_affected
                        else:
                            value_affected = ''
                            value_total = ''

                        vep_key = '_'.join(('VEP', csq_key_list[i], 'AFFECTED'))
                        if vep_key not in vep_dict:
                            vep_dict[vep_key] = initialise_allele_list(length=allele_length)
                        vep_dict[vep_key][allele_number].append(value_affected)

                        vep_key = '_'.join(('VEP', csq_key_list[i], 'TOTAL'))
                        if vep_key not in vep_dict:
                            vep_dict[vep_key] = initialise_allele_list(length=allele_length)
                        vep_dict[vep_key][allele_number].append(value_total)

                    elif csq_key_list[i] in ('cDNA_position', 'CDS_position', 'Protein_position'):
                        if csq_value_list[i]:
                            value_list = csq_value_list[i].split('-')
                            value_start = value_list[0]
                            # Make sure value end is defined.
                            if len(value_list) == 2:
                                value_end = value_list[1]
                            else:
                                value_end = value_start
                        else:
                            value_start = ''
                            value_end = ''

                        vep_key = '_'.join(('VEP', csq_key_list[i], 'START'))
                        if vep_key not in vep_dict:
                            vep_dict[vep_key] = initialise_allele_list(length=allele_length)
                        vep_dict[vep_key][allele_number].append(value_start)

                        vep_key = '_'.join(('VEP', csq_key_list[i], 'END'))
                        if vep_key not in vep_dict:
                            vep_dict[vep_key] = initialise_allele_list(length=allele_length)
                        vep_dict[vep_key][allele_number].append(value_end)

                    else:
                        # For all other keys.
                        vep_key = '_'.join(('VEP', csq_key_list[i]))
                        if vep_key not in vep_dict:
                            vep_dict[vep_key] = initialise_allele_list(length=allele_length)
                        vep_dict[vep_key][allele_number].append(csq_value_list[i])

            # For each allele, select exactly one index from the list of possible consequences.
            # Thereby, the most severe sequence ontology term that has also the PICK flag gets prioritised.
            # Failing a PICK flag, in case the filter script has removed all, the most severe consequence gets picked.
            priority_list: List[Optional[int]] = list()
            for i in range(0, allele_length):
                # VEP consequences can be '&' separated values (e.g., splice_region_variant&synonymous_variant).
                consequence_tuple_list: List[Tuple[int, str]] = list()
                for vep_consequence_index, vep_consequence in enumerate(vep_dict['VEP_Consequence'][i]):
                    for value in vep_consequence.split('&'):
                        consequence_tuple_list.append((vep_consequence_index, value))

                module_logger.log(logging.DEBUG - 1, 'Consequence list old: %r', consequence_tuple_list)

                consequence_tuple_list.sort(key=lambda item: get_consequence_index(item[1]))

                module_logger.log(logging.DEBUG - 1, 'Consequence list new: %r', consequence_tuple_list)
                if consequence_tuple_list:
                    module_logger.log(
                        logging.DEBUG - 1,
                        'Prioritised index: %r consequence: %r',
                        consequence_tuple_list[0][0],
                        consequence_tuple_list[0][1])

                if consequence_tuple_list:
                    # Search for the entry with the highest consequence priority, which is also picked.
                    for consequence_tuple in consequence_tuple_list:
                        if vep_dict['VEP_PICK'][i][consequence_tuple[0]] == '1':
                            module_logger.log(logging.DEBUG - 1, 'PICK flag set: %r', consequence_tuple[0])
                            priority_list.append(consequence_tuple[0])
                            break
                    else:
                        # Failing an entry with the PICK flag, select the top entry of the priority list.
                        module_logger.log(logging.DEBUG - 1, 'PICK flag not set: %r', consequence_tuple_list[0][0])
                        priority_list.append(consequence_tuple_list[0][0])
                else:
                    # Always append to keep the order of allele indices.
                    priority_list.append(None)

            # Now, prune the VEP dict in order.

            for vep_key in sorted(vep_header_dict):
                new_key = 'VEP_' + vep_key
                if new_key not in vep_dict:
                    continue

                # Prune feature and allele lists by checking for and removing empty entries.
                empty_allele_list = True
                for i in range(0, allele_length):
                    empty_feature_list = True
                    # Apply the prioritised sequence ontology term from above.
                    if priority_list[i] is not None:
                        vep_dict[new_key][i] = [vep_dict[new_key][i][priority_list[i]]]
                    for vep_value in vep_dict[new_key][i]:
                        if vep_value:
                            empty_feature_list = False
                    if empty_feature_list:
                        vep_dict[new_key][i] = []
                    else:
                        empty_allele_list = False
                if empty_allele_list:
                    continue

                # Collate (i.e. only select the first allele-specific value) or completely exclude a field.
                # Join values for each alternative allele, but start at index 1,
                # since 0 represents the reference allele.
                if vep_header_dict[vep_key]['Allele'] == 'collate':
                    # Only select the first value for each alternative allele.
                    variant_record_info[new_key] = \
                        [collate_allele_values(_allele_list) for _allele_list in vep_dict[new_key]][1:]
                elif vep_header_dict[vep_key]['Allele'] == 'exclude':
                    pass
                else:
                    # Join all values for each alternative allele, but start at index 1,
                    # since 0 represents the reference allele.
                    variant_record_info[new_key] = ['|'.join(_allele_list) for _allele_list in vep_dict[new_key]][1:]

            # Remove the original CSQ field from the record.
            del variant_record_info['CSQ']

            variant_file_new.write(record=variant_record)

        variant_file_new.close()
        variant_file_old.close()

        return None


if __name__ == '__main__':
    argument_parser = ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--logging-level',
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'DEBUG1', 'DEBUG2'],
        default='INFO',
        dest='logging_level',
        help='Logging level [INFO]',
        required=False)

    argument_parser.add_argument(
        '--soc-path',
        dest='soc_path',
        help='Sequence Ontology terms (TSV) configuration file path',
        required=True)

    argument_parser.add_argument(
        '--ofc-path',
        dest='ofc_path',
        help='Output fields (TSV) configuration file path',
        required=True)

    argument_parser.add_argument(
        '--old-vcf-path',
        dest='old_vcf_path',
        help='Old (input) VCF file path',
        required=True)

    argument_parser.add_argument(
        '--new-vcf-path',
        dest='new_vcf_path',
        help='New (output) VCF file path',
        required=True)

    name_space = argument_parser.parse_args()

    if name_space.logging_level:
        logging.addLevelName(level=logging.DEBUG - 1, levelName='DEBUG1')
        logging.addLevelName(level=logging.DEBUG - 2, levelName='DEBUG2')

        logging.basicConfig(level=name_space.logging_level)

    runnable_step = RunnableStepCsqToVep(
        name='csq_to_vep',
        soc_path=name_space.soc_path,
        ofc_path=name_space.ofc_path,
        vcf_path_old=name_space.old_vcf_path,
        vcf_path_new=name_space.new_vcf_path)

    runnable_step.run()
