# !/usr/bin/python
# coding=utf-8
import os

from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from ProteinProphetParser.ProtProphXMLParser import ProtProphXMLParser
from WebLogo33.WebLogoGenerator import WebLogoGenerator

"""
Jeroen Merks

"""


def create_out_folder(out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    return out_folder


def get_prots_min(min_prob, prot_groups):
    matching_prots = []

    for prot_group in prot_groups:
        for prot in prot_group.get_prots():

            if prot.get_prob() >= min_prob:
                matching_prots.append(prot)

    return matching_prots


def sort_hits(min_prob, prot_groups):
    matching_prots = get_prots_min(min_prob, prot_groups)

    # Sort matches on probability
    sorted_matching_hits = sorted(matching_prots, key=lambda protgroup: protgroup.get_prob(), reverse=True)

    return matching_prots, sorted_matching_hits


def write_readable_result_hit(readable_out_file, prot, protein_db):
    readable_out_file.write('Protein ID: %s\n' % prot.get_prot_name())
    readable_out_file.write('Protein description: %s\n' % prot.get_descr())
    readable_out_file.write('Probability: %f\n' % prot.get_prob())
    readable_out_file.write('%s\n' % prot.get_seq(protein_db))
    for pep in prot.get_peptides():
        readable_out_file.write('%s\n' % pep.get_seq())

    readable_out_file.write('\n')


def save_hit(motif_range_start, motif_range_end, prot, protein_db, weblogo_in_file, readable_out_file):

    # Write a multiple fasta file with proteis within the given motif range
    prot_seq = SeqRecord(
        Seq(prot.get_seq(protein_db)[motif_range_start:motif_range_end], generic_protein),
        id=prot.get_prot_name(),
        description=prot.get_descr())
    SeqIO.write(prot_seq, weblogo_in_file, "fasta")

    write_readable_result_hit(readable_out_file, prot, protein_db)


def find_metap_activity(min_prob, cleavage_loc, motif_range_start,
                        motif_range_end, readable_out_path, prots,
                        protein_db, weblogo_in_path):

    found_metap_activity = False

    weblogo_in_file = open(weblogo_in_path, "a+")
    readable_out_file = open(readable_out_path, "a+")

    matching_prots, matching_prot_hits = sort_hits(min_prob, prots)

    if len(matching_prot_hits) > 0:

        for prot in matching_prot_hits:
            for pep in prot.get_peptides():

                # Locate the position of the peptide on the protein
                pep_pos = prot.get_seq(protein_db).find(pep.get_seq())

                if pep_pos == cleavage_loc:
                    found_metap_activity = True
                    save_hit(motif_range_start, motif_range_end, prot, protein_db, weblogo_in_file, readable_out_file)

        if not found_metap_activity:
            print("No metap-activity has been detected!")

    else:
        print('None of the prot met the given critera.')


def run_metap_pipeline(ms_run_code, prot_prophet_xml, protein_db_path, min_prob, cleavage_loc, motif_range_start,
                       motif_range_end):

    protxml_name = prot_prophet_xml[prot_prophet_xml.rfind("/"):prot_prophet_xml.rfind(".")]

    protein_db = SeqIO.index(protein_db_path, format='fasta')
    parse_results_out_paht = create_out_folder("protxml_parsing_out/")
    weblogos_out_path = create_out_folder('weblogo_out/')
    readable_out_name = parse_results_out_paht + protxml_name + ms_run_code + '_human_readable.txt'
    weblogo_in_name = parse_results_out_paht + protxml_name + ms_run_code + '_weblogo_in.fasta'

    prot_xml_parser = ProtProphXMLParser(parse_results_out_paht, prot_prophet_xml)

    prots = [prot for prot in prot_xml_parser.get_prot_groups()]
    weblogo_generator = WebLogoGenerator(weblogos_out_path)

    print("Parsing prot.xml\n")
    find_metap_activity(min_prob, cleavage_loc, motif_range_start,
                        motif_range_end, readable_out_name, prots,
                        protein_db, weblogo_in_name)

    weblogo_generator.create_weblogo(weblogo_in_name)
