#!/usr/bin/python
# coding=utf-8

import sys

from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio import Seq, SeqIO

from ProtProphXMLParser import ProtXMLParser
from WebLogoGenerator import WebLogoGenerator

"""
Jeroen Merks

Instructions:
 - Install BioPython 1.65: http://biopython.org/wiki/Download

"""


def pretty_print_alignment(prot_seq, peptides):
    """
    # TODO fix alligning overlapping peptides

    :param prot_seq:
    :param peptides:
    :return:
    """
    sorted_peps = sorted(peptides, key=lambda pep_seq: prot_seq.find(pep_seq))

    first_start_pos = prot_seq.find(sorted_peps[0])
    alignment = "" + first_start_pos * " "
    alignment += sorted_peps[0]
    prev_end = len(alignment)

    for peptide in sorted_peps[1:]:
        space_gap = prot_seq.find(peptide) - prev_end
        alignment += (space_gap * " ") + peptide
        prev_end = len(alignment)

    return prot_seq, alignment


def get_prots(min_prob, prot_groups):
    matching_prots = []
    for prot_group in prot_groups:
        for prot in prot_group.get_prots():
            if prot.get_prob() >= min_prob:
                matching_prots.append(prot)
    return matching_prots


def write_readable_results(file, prot, protein_db):
    file.write('Protein ID: %s\n' % prot.get_id())
    file.write('Protein description: %s\n' % prot.get_descr())
    file.write('Probability: %f\n' % prot.get_prob())
    file.write("%s\n%s\n" % (
        pretty_print_alignment(prot.get_seq(protein_db), [pep.get_seq() for pep in prot.get_peptides()])))
    file.write('\n')


def save_readable_prot_infos_min_prob(min_prob, prot_groups, protein_db, readable_prot_file):
    """

    :param readable_prot_file:
    :param protein_db:
    :param prot_groups:
    :param min_prob:
    """

    # Check if probability falls between 0 and 1
    if 1 >= min_prob > 0:
        matching_prots = get_prots(min_prob, prot_groups)
        # Sort matches on probability
        sorted_matching_prot_hits = sorted(matching_prots, key=lambda protein: protein.get_prob(), reverse=True)
    else:
        sys.exit('Min probability must range between 0 and 1.')

    if len(sorted_matching_prot_hits) > 0:
        file = open(readable_prot_file, "w")
        file.write('Nr of proteins matching criteria: %i\n' % len(matching_prots))
        file.write('\n')
        for prot in sorted_matching_prot_hits:
            write_readable_results(file, prot, protein_db)
    else:
        print('None of the proteins met the given critera.')


def find_metap_activity(min_prob, cleavage_loc, motif_range_start, motif_range_end, write_to_fasta, fasta_out,
                        readable_out,
                        prot_groups, protein_db):
    """

    :param min_prob:
    :param cleavage_loc:
    :param motif_range_start:
    :param motif_range_end:
    :param write_to_fasta:
    :param fasta_out:
    :param readable_out:
    :param prot_groups:
    :param protein_db:
    """

    found_metap_activity = False

    # Check if cleavage location is valid
    if cleavage_loc <= 0:
        sys.exit('The cleavage site has to be on a position > 0.')

    # Check if probability falls between 0 and 1
    if 1 >= min_prob > 0:
        matching_prots = get_prots(min_prob, prot_groups)

        # Sort matches on probability
        matching_hits = sorted(matching_prots, key=lambda protgroup: protgroup.get_prob(),
                               reverse=True)
    else:
        sys.exit('Min probability must range between 0 and 1.')

    if len(matching_hits) > 0:

        if write_to_fasta:
            output_handler = open(fasta_out, "w")
        else:
            file = open(readable_out, "w")

        for prot in matching_hits:
            for pept in prot.get_peptides():

                if prot.get_seq(protein_db).find(pept.get_seq()) == cleavage_loc:
                    found_metap_activity = True

                    if write_to_fasta:
                        prot_seq = SeqRecord(
                            Seq.Seq(prot.get_seq(protein_db)[motif_range_start:motif_range_end], generic_protein),
                            id=prot.get_id(), description=prot.get_descr())
                        SeqIO.write(prot_seq, output_handler, "fasta")
                    else:
                        file.write('Nr of proteins matching criteria: %i\n' % len(matching_prots))
                        file.write('\n')
                        write_readable_results(file, prot, protein_db)

        if not found_metap_activity:
            print("No metap-activity has been detected!")

    else:
        print('None of the prot met the given critera.')


def main(prot_prophet_xml, protein_db_name):
    """

    :param prot_prophet_xml:
    :param protein_db_name:
    """

    protein_db = SeqIO.index(protein_db_name, format='fasta')
    parse_results_folder = "parse_results/"
    weblogos_output_folder = 'weblogos_results/'
    readable_out_file = parse_results_folder + 'sorted_prots.txt'
    readable_metap_act_out_file = parse_results_folder + 'sorted_metap_activity_human_readable.txt'
    metap_act_out_fasta_file = parse_results_folder + 'sorted_metap_activity.fasta'

    prot_xml_parser = ProtXMLParser(parse_results_folder, prot_prophet_xml)

    statistics_dict = prot_xml_parser.get_statistics_dict()  # {min_probability: false_positive_error_rate}
    prot_groups = prot_xml_parser.get_prot_groups()
    prots = [prot for prot in prot_groups]

    weblogo_generator = WebLogoGenerator(weblogos_output_folder)

    save_readable_prot_infos_min_prob(statistics_dict[0.004], prots, protein_db, readable_out_file)

    # Test fasta output
    find_metap_activity(min_prob=statistics_dict[0.004], cleavage_loc=2, motif_range_start=0, motif_range_end=5,
                        write_to_fasta=True, fasta_out=metap_act_out_fasta_file,
                        readable_out=readable_metap_act_out_file,
                        prot_groups=prots, protein_db=protein_db)
    # Test human readable output
    find_metap_activity(min_prob=statistics_dict[0.004], cleavage_loc=2, motif_range_start=0, motif_range_end=5,
                        write_to_fasta=False, fasta_out=metap_act_out_fasta_file,
                        readable_out=readable_metap_act_out_file,
                        prot_groups=prot_groups, protein_db=protein_db)

    weblogo_generator.create_weblogo(metap_act_out_fasta_file)


if __name__ == '__main__':
    main("GitHub_test_files/raw.comet.interact.prot.xml",
         'GitHub_test_files/Uniprot_Mt_proteome.fasta')