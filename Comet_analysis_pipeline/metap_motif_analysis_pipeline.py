# !/usr/bin/python
# coding=utf-8
import os
import sys

from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio import SeqIO

from ProteinProphetParser.ProtProphXMLParser import ProtProphXMLParser
from WebLogo33.WebLogoGenerator import WebLogoGenerator


"""
Jeroen Merks

"""

def get_prots_min(min_prob, prot_groups):
    matching_prots = []

    for prot_group in prot_groups:
        for prot in prot_group.get_prots():

            if prot.get_prob() >= min_prob:
                matching_prots.append(prot)

    return matching_prots


def write_human_readable_results(file, prot, protein_db):
    file.write('Protein ID: %s\n' % prot.get_prot_name())
    file.write('Protein description: %s\n' % prot.get_descr())
    file.write('Probability: %f\n' % prot.get_prob())
    file.write('%s\n' % prot.get_seq(protein_db))
    for pep in prot.get_peptides():
        file.write('%s\n' % pep.get_seq())

    # file.write("%s\n%s\n" % (
    #     pretty_print_alignment(prot.get_seq(protein_db), [pep.get_seq() for pep in prot.get_peptides()])))

    file.write('\n')


def save_readable_prot_infos_min_prob(min_prob, prot_groups, protein_db, readable_prot_file):
    """

    :param readable_prot_file:
    :param protein_db:
    :param prot_groups:
    :param min_prob:
    """

    # Check if probability falls between 0 and 1
    if 1 > min_prob or min_prob > 0:
        matching_prots = get_prots_min(min_prob, prot_groups)
        # Sort matches on probability
        # TODO check if this is even necessary, ProteinProphet seems to sort on probability by default?
        sorted_matching_prot_hits = sorted(matching_prots, key=lambda protein: protein.get_prob(), reverse=True)
    else:
        sys.exit('Min probability must range between 0 and 1.')

    if len(sorted_matching_prot_hits) > 0:
        file = open(readable_prot_file, "w")
        file.write('Nr of proteins matching criteria: %i\n' % len(matching_prots))
        file.write('\n')
        for prot in sorted_matching_prot_hits:
            write_human_readable_results(file, prot, protein_db)
        file.close()
    else:
        print('None of the proteins met the given critera.')


def sort_hits(min_prob, prot_groups):
    # Check if probability falls between 0 and 1
    if 1 > min_prob or min_prob > 0:
        matching_prots = get_prots_min(min_prob, prot_groups)

        # Sort matches on probability
        sorten_matching_hits = sorted(matching_prots, key=lambda protgroup: protgroup.get_prob(), reverse=True)
    else:
        sys.exit('Min probability must range between 0 and 1.')

    return matching_prots, sorten_matching_hits


def save_results(write_weblogo_in, motif_range_end, motif_range_start,
                 prot, protein_db, out_file):
    if write_weblogo_in:
        # Write a multiple fasta file with proteis within the given motif range
        prot_seq = SeqRecord(
            Seq(prot.get_seq(protein_db)[motif_range_start:motif_range_end], generic_protein),
            id=prot.get_prot_name(),
            description=prot.get_descr())
        SeqIO.write(prot_seq, out_file, "fasta")

    else:
        write_human_readable_results(out_file, prot, protein_db)


def create_append_out_file(weblogo_out, readable_out, write_to_fasta):
    if write_to_fasta:
        if os.path.exists(weblogo_out):
            os.remove(weblogo_out)
        out_file = open(weblogo_out, "a+")
    else:
        if os.path.exists(readable_out):
            os.remove(readable_out)
        out_file = open(readable_out, "a+")

    return out_file


def create_out_folder(out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    return out_folder


def find_metap_activity(min_prob, cleavage_loc, motif_range_start,
                        motif_range_end, readable_out, prot_groups,
                        protein_db, weblogo_out, write_for_weblogo):
    """

    :param min_prob:
    :param cleavage_loc:
    :param motif_range_start:
    :param motif_range_end:
    :param write_for_weblogo:
    :param weblogo_out:
    :param readable_out:
    :param prot_groups:
    :param protein_db:
    """

    found_metap_activity = False

    # Check if cleavage location is valid
    if cleavage_loc < 0:
        sys.exit('The cleavage site has to be on a position > 0.')

    matching_prots, matching_prot_hits = sort_hits(min_prob, prot_groups)

    if len(matching_prot_hits) > 0:

        out_file = create_append_out_file(weblogo_out, readable_out, write_for_weblogo)

        for prot in matching_prot_hits:
            for pep in prot.get_peptides():

                # Locate the position of the peptide on the protein
                pep_pos = prot.get_seq(protein_db).find(pep.get_seq())

                if pep_pos == cleavage_loc:
                    found_metap_activity = True
                    save_results(write_for_weblogo, motif_range_end,
                                 motif_range_start, prot, protein_db,
                                 out_file)

        if not found_metap_activity:
            print("No metap-activity has been detected!")

    else:
        print('None of the prot met the given critera.')


def choose_prob(error_prob_pairs):
    """

    :param error_prob_pairs:
    :return:
    """

    # TODO make elephant-proof
    print("Pair#\terror rate\tprob")
    for pair_nr, pair in enumerate(error_prob_pairs):
        print("%i:\t\t%f\t%f" % (pair_nr + 1, pair[0], pair[1]))

    print("Which pair of error rate / probability should be used?")
    return error_prob_pairs[int(input(">>> ")) - 1][1]


def run_metap_pipeline(ms_run_code, prot_prophet_xml, protein_db_name):
    """

    :param prot_prophet_xml:
    :param protein_db_name:
    """

    protxml_name = prot_prophet_xml[prot_prophet_xml.rfind("/"):prot_prophet_xml.rfind(".")]

    protein_db = SeqIO.index(protein_db_name, format='fasta')
    parse_results_folder = create_out_folder("protxml_parsing_out/")
    weblogos_output_folder = create_out_folder('weblogo_out/')
    readable_metap_act_out_file = parse_results_folder + protxml_name + ms_run_code + '_human_readable.txt'
    metap_act_weblogo_in = parse_results_folder + protxml_name + ms_run_code + '_weblogo_in.fasta'

    prot_xml_parser = ProtProphXMLParser(parse_results_folder, prot_prophet_xml)
    statistics = prot_xml_parser.get_error_prob_pairs()
    chosen_prob = choose_prob(statistics)
    prots = [prot for prot in prot_xml_parser.get_prot_groups()]
    weblogo_generator = WebLogoGenerator(weblogos_output_folder)

    # Write all hits with a given minimum, for transparency sake
    # save_readable_prot_infos_min_prob(chosen_prob, prots, protein_db, readable_out_file)

    # Fasta output, for WebLogo
    find_metap_activity(min_prob=chosen_prob, cleavage_loc=1,
                        motif_range_start=1, motif_range_end=6,
                        readable_out=readable_metap_act_out_file,
                        prot_groups=prots, protein_db=protein_db,
                        weblogo_out=metap_act_weblogo_in,
                        write_for_weblogo=True)
    # Human readable output
    find_metap_activity(min_prob=chosen_prob, cleavage_loc=1,
                        motif_range_start=1, motif_range_end=6,
                        readable_out=readable_metap_act_out_file,
                        prot_groups=prots, protein_db=protein_db,
                        weblogo_out=metap_act_weblogo_in,
                        write_for_weblogo=False)

    weblogo_generator.create_weblogo(metap_act_weblogo_in)
