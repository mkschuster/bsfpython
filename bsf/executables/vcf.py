#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""Variant Calling Format (VCF) Executables module.

A package of classes and functions to split Ensembl VEP annotation in the CSQ field into VEP_* fields.
http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
"""
#  Copyright 2013 - 2019 Michael K. Schuster
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

from __future__ import print_function

import argparse
import csv
import warnings

import pysam

import bsf.process


class RunnableStepCsqToVep(bsf.process.RunnableStep):
    """The C{bsf.executables.vcf.RunnableStepCsqToVep} class expands Ensembl Variant Effect Predictor (VEP)
    CSQ INFO annotation into a set of VEP_* INFO annotation.

    Attributes:
    @ivar soc_path: Sequence Ontology term priority configuration (TSV) file path
    @type soc_path: str | None
    @ivar ofc_path: Output field configuration (TSV) file path
    @type ofc_path; str | None
    @ivar vcf_path_old: Old VCF file path
    @type vcf_path_old: str | None
    @ivar vcf_path_new: New VCF file path
    @type vcf_path_new: str | None
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
        """Initialise a C{bsf.executables.vcf.RunnableStepCsqToVep}.

        @param name: Name
        @type name: str | None
        @param program: Program
        @type program: str | None
        @param options: Python C{dict} of Python C{str} (C{bsf.argument.Argument.key}) key and
            Python C{list} value objects of C{bsf.argument.Argument} objects
        @type options: dict[bsf.argument.Argument.key, list[bsf.argument.Argument]] | None
        @param arguments: Python C{list} of Python C{str} (program argument) objects
        @type arguments: list[str] | None
        @param sub_command: Subordinate C{bsf.process.Command}
        @type sub_command: bsf.process.Command | None
        @param stdin: Standard input I{STDIN} C{bsf.connector.Connector}
        @type stdin: bsf.connector.Connector | None
        @param stdout: Standard output I{STDOUT} C{bsf.connector.Connector}
        @type stdout: bsf.connector.Connector | None
        @param stderr: Standard error I{STDERR} C{bsf.connector.Connector}
        @type stderr: bsf.connector.Connector | None
        @param dependencies: Python C{list} of C{bsf.process.Executable.name}
            properties in the context of C{bsf.analysis.Stage} dependencies
        @type dependencies: list[bsf.process.Executable.name] | None
        @param hold: Hold on job scheduling
        @type hold: str | None
        @param submit: Submit the C{bsf.process.Executable} during C{bsf.analysis.Stage.submit}
        @type submit: bool
        @param process_identifier: Process identifier
        @type process_identifier: str | None
        @param process_name: Process name
        @type process_name: str | None
        @param sub_process: C{subprocess.Popen}
        @type sub_process: subprocess.Popen | None
        @param obsolete_file_path_list: Python C{list} of file paths that can be removed
            after successfully completing this C{bsf.process.RunnableStep}
        @type obsolete_file_path_list: list[str] | None
        @param soc_path: Sequence Ontology term priority configuration (TSV) file path
        @type soc_path: str | None
        @param ofc_path: Output field configuration (TSV) file path
        @type ofc_path; str | None
        @param vcf_path_old: Old VCF file path
        @type vcf_path_old: str | None
        @param vcf_path_new: New VCF file path
        @type vcf_path_new: str | None
        @return:
        @rtype:
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

    def run(self, max_thread_joins=10, thread_join_timeout=10, debug=0):
        """Run a C{bsf.executables.vcf.RunnableStepCsqToVep}.

        @param max_thread_joins: Maximum number of attempts to join the output threads
        @type max_thread_joins: int
        @param thread_join_timeout: Timeout for each attempt to join the output threads
        @type thread_join_timeout: int
        @param debug: Debug level
        @type debug: int
        @return: Return value of the child in the Python subprocess,
            negative values indicate that the child received a signal
        @rtype: int
        """

        def initialise_allele_list(length):
            """Initialise an empty Python C{list} of Python C{list} of Python C{str} objects.

            The list stores one Python C{list} object per allele and picked feature.
            @param length: Length
            @return: Python C{list} of Python C{list} of Python C{str} objects
            @rtype: list[list[str]]
            """
            allele_list = list()

            j = 0
            while j < length:
                allele_list.append(list())
                j += 1

            return allele_list

        def collate_allele_values(allele_list):
            """Collate the list of allele-specific values by selecting only the first list component.

            @param allele_list: Python C{list} of allele-specific values
            @type allele_list: list[str]
            @return: First value
            @rtype: str | None
            """

            if allele_list:
                return allele_list[0]
            else:
                return None

        def get_consequence_index(consequence):
            """Get the SO index.

            @param consequence: Consequence
            @type consequence: str
            @return: Consequence index
            @rtype: int
            """
            # print('get_consequence_index:', consequence)
            for _so_index, _so_string in enumerate(sequence_ontology_list):
                # print('get_consequence_index', consequence, '_so_index:', _so_index, '_so_string:', _so_string)
                if consequence == _so_string:
                    return _so_index

        # Read the Sequence Ontology (TSV) configuration file.
        # This file provides consequence prioritisation.
        sequence_ontology_list = list()
        """ @type sequence_ontology_list: list[str] """

        with open(file=self.soc_path, mode='rt') as input_file:
            for row_dict in csv.DictReader(input_file, dialect='excel-tab'):
                """ @type row_dict: dict[str, str] """
                sequence_ontology_list.append(row_dict['SO term'])

        # Load the VEP output field (TSV) configuration file.
        vep_header_dict = dict()
        """ @type vep_header_dict: dict[str, dict[str, str]] """

        with open(file=self.ofc_path, mode='rt') as input_file:
            for row_dict in csv.DictReader(input_file, dialect='excel-tab'):
                """ @type row_dict: dict[str, str] """
                # NOTE: Override all types with 'String' to allow multiple values joined by '&' characters.
                row_dict['Type'] = 'String'
                vep_header_dict[row_dict['ID']] = row_dict

        vf_old = pysam.VariantFile(self.vcf_path_old, 'r')

        # Copy the header for a new VCF instance.

        vh_new = vf_old.header.copy()
        """ @type vh_new: pysam.libcbcf.VariantHeader """

        if 'CSQ' not in vf_old.header.info:
            raise Exception("Cannot convert a VCF file without a 'CSQ' INFO field.")

        # Remove the CSQ info field from the new header.

        # NOTE: Apparently, pysam does not allow deleting an entry from the variant header at the moment.
        # The pysam developers commented out calling the bcf_hdr_sync() function after calling bcf_hdr_remove() in
        # VariantHeaderMetadata.remove_header().
        # https://github.com/pysam-developers/pysam/blob/16ecfee0c12908f5b42ff0542d75666475ca4d1a/pysam/libcbcf.pyx#L1616
        #
        # vh_new.info.remove_header('CSQ')

        csq = vf_old.header.info['CSQ']
        """ @type csq: pysam.libcbcf.VariantMetadata """

        # Build a new header with VEP_ INFO fields.
        # Split the description on white space first, the pipe-separated list of VEP field declarations
        # is then in the last block.

        csq_key_list = csq.description.split()[-1].split('|')
        """ @type csq_key_list: list[str] """
        csq_index_allele_num = None
        """ @type csq_index_allele_num: int | None """

        for index, csq_key in enumerate(csq_key_list):
            if csq_key == 'ALLELE_NUM':
                csq_index_allele_num = index

            if csq_key in vep_header_dict:
                vh_new.info.add(
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
                    vh_new.info.add(
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
                    vh_new.info.add(
                        id='_'.join(('VEP', csq_key, suffix)),
                        number=vep_header_dict[new_field]['Number'],
                        type=vep_header_dict[new_field]['Type'],
                        description=vep_header_dict[new_field]['Description'])
                continue

            warnings.warn('CSQ VEP field ' + repr(csq_key) + ' missing from field configuration file.', UserWarning)

        if debug > 1:
            # Print the new header.
            for key, vmd in vh_new.info.items():
                """ @type key: str """
                """ @type value: pysam.libcbcf.VariantMetadata """
                # print('Key:', repr(key), 'VariantMetadata:', repr(vmd))
                # print(type(vmd), 'dir:', dir(vmd))
                print(
                    'Key:', repr(key),
                    'Value:', repr(vmd),
                    'ID:', repr(vmd.id),
                    'Name:', repr(vmd.name),
                    'Number:', repr(vmd.number),
                    'Type:', repr(vmd.type),
                    'Description:', repr(vmd.description))
            print()

        # Open the new VariantFile for writing.

        vf_new = pysam.VariantFile(self.vcf_path_new, 'w', header=vh_new)

        # Parse the CSQ INFO field of each VariantRecord.

        for vr in vf_old.fetch():
            """ @type vr: pysam.libcbcf.VariantRecord """
            vr.translate(vh_new)
            vri = vr.info
            """ @type vri: pysam.libcbcf.VariantRecordInfo """

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
            vep_dict = dict()
            """ @type vep_dict: dict[str, list[list[str]]] """

            # Total number of alleles including REF and ALT.
            allele_length = len(vr.alleles)

            if 'CSQ' not in vri:
                vf_new.write(vr)
                continue

            # Iterate over all comma-separated allele-transcript blocks.
            for csq_component in vri['CSQ']:
                """ @type csq_component: str """
                csq_value_list = csq_component.split('|')

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
            priority_list = list()
            """ @type priority_list: list[int | None] """
            for i in range(0, allele_length):
                # VEP consequences can be '&' separated values (e.g. splice_region_variant&synonymous_variant).
                consequence_tuple_list = list()
                """ @type consequence_tuple_list: list[(int, str)] """
                for index, vep_consequence in enumerate(vep_dict['VEP_Consequence'][i]):
                    for value in vep_consequence.split('&'):
                        consequence_tuple_list.append((index, value))

                if debug > 1:
                    print('Consequence list old:', repr(consequence_tuple_list))

                consequence_tuple_list.sort(key=lambda item: get_consequence_index(item[1]))

                if debug > 1:
                    print('Consequence list new:', repr(consequence_tuple_list))
                    if consequence_tuple_list:
                        print('Prioritised index:', consequence_tuple_list[0][0],
                              'consequence:', consequence_tuple_list[0][1])
                    print()

                if consequence_tuple_list:
                    # Search for the entry with the highest consequence priority, which is also picked.
                    for consequence_tuple in consequence_tuple_list:
                        if vep_dict['VEP_PICK'][i][consequence_tuple[0]] == '1':
                            if debug > 1:
                                print('PICK flag set:', consequence_tuple[0])
                            priority_list.append(consequence_tuple[0])
                            break
                    else:
                        # Failing an entry with the PICK flag, select the top entry of the priority list.
                        if debug > 1:
                            print('PICK flag not set:', consequence_tuple_list[0][0])
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
                    vri[new_key] = [collate_allele_values(_allele_list) for _allele_list in vep_dict[new_key]][1:]
                elif vep_header_dict[vep_key]['Allele'] == 'exclude':
                    pass
                else:
                    # Join all values for each alternative allele, but start at index 1,
                    # since 0 represents the reference allele.
                    vri[new_key] = ['|'.join(_allele_list) for _allele_list in vep_dict[new_key]][1:]

            # Remove the original CSQ field from the record.
            del vri['CSQ']

            vf_new.write(vr)

        return 0


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(
        description='Module driver script.')

    argument_parser.add_argument(
        '--debug',
        help='Debug level',
        required=False,
        type=int)

    argument_parser.add_argument(
        '--soc-path',
        dest='soc_path',
        help='Sequence Ontology terms (TSV) configuration file path',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--ofc-path',
        dest='ofc_path',
        help='Output fields (TSV) configuration file path',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--old-vcf-path',
        dest='old_vcf_path',
        help='Old (input) VCF file path',
        required=True,
        type=str)

    argument_parser.add_argument(
        '--new-vcf-path',
        dest='new_vcf_path',
        help='New (output) VCF file path',
        required=True,
        type=str)

    name_space = argument_parser.parse_args()

    runnable_step = RunnableStepCsqToVep(
        name='csq_to_vep',
        soc_path=name_space.soc_path,
        ofc_path=name_space.ofc_path,
        vcf_path_old=name_space.old_vcf_path,
        vcf_path_new=name_space.new_vcf_path)

    runnable_step.run(debug=name_space.debug)
