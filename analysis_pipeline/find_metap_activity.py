# !/usr/bin/python
# coding=utf-8
import os
import csv

from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from XinteractParser.XinteractParser import XinteractParser
from WebLogo34.WebLogoGenerator import WebLogoGenerator


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


def filter_and_sort_hits(min_prob, prot_groups):
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


def extract_metap_records_weblogo_in(complete_proteome, metap_db, proteome_without_metap_fasta_out_name, motif_len):
    complete_proteome_db = SeqIO.index(complete_proteome, 'fasta')

    known_ids = [rec.id for rec in SeqIO.parse(complete_proteome, 'fasta')]
    metap_ids = [rec.id for rec in SeqIO.parse(metap_db, 'fasta')]

    out_file = open(proteome_without_metap_fasta_out_name, "w+")

    for known_id in known_ids:
        if known_id not in metap_ids:
            SeqIO.write(complete_proteome_db.get(known_id)[:motif_len], out_file, "fasta")


def cleave_db_motif_len(prot_db, motif_len, proteome_weblogo_in_name):
    complete_proteome_db = SeqIO.parse(prot_db, 'fasta')
    proteome_weblogo_in_file = open(proteome_weblogo_in_name, "w+")

    for rec in complete_proteome_db:
        prot_seq = SeqRecord(
            Seq(str(rec.seq)[:motif_len], generic_protein), id=rec.id, description=rec.description)
        SeqIO.write(prot_seq, proteome_weblogo_in_file, "fasta")

    proteome_weblogo_in_file.close()



def write_maxquant_result_hit(maxquant_csv_writer, prot, protein_db):
    pass


# for pep in prot.get_peptides:
## "Sequence"
# pep.get_seq(protein_db),
## "N-terminal"
# "-",
## "C-terminal"
# "-",
## "Modifications"
# "-",
## "Mass"
# pep.get_neutral_mass(),
## "Mass Fractional Part"
# "-",
## "Protein Groups"
# prot.
## "Proteins"
## "Unique (Groups)"
## "Unique (Proteins)"
## "Acetyl (Protein N-term)"
## "Oxidation (M)"
## "Missed cleavages"
## "Retention time"
## "Calibrated retention time"
## "Charges"
## "PEP"
## "MS/MS scan number"
## "Raw file"
## "Score"
## "Delta score"
## "Intensity"
## "Reverse"
## "Potential contaminant"
## "id"
## "Protein group IDs"
## "Peptide ID"
## "Evidence IDs"
## "MS/MS IDs"
## "Best MS/MS"
## "Oxidation (M) site IDs"
## "MS/MS Count"

# maxquant_csv_writer.writerow(
#     [prot.get_seq(protein_db),
#      "-",
#      "-",
#      "Modifications",
#      "Mass",
#      "Mass Fractional Part",
#      "Protein Groups",
#      "Proteins",
#      "Unique (Groups)",
#      "Unique (Proteins)",
#      "Acetyl (Protein N-term)",
#      "Oxidation (M)",
#      "Missed cleavages",
#      "Retention time",
#      "Calibrated retention time",
#      "Charges",
#      "PEP",
#      "MS/MS scan number",
#      "Raw file",
#      "Score",
#      "Delta score",
#      "Intensity",
#      "Reverse",
#      "Potential contaminant",
#      "id",
#      "Protein group IDs",
#      "Peptide ID",
#      "Evidence IDs",
#      "MS/MS IDs",
#      "Best MS/MS",
#      "Oxidation (M) site IDs",
#      "MS/MS Count"])


def save_hit(motif_len, prot, protein_db, weblogo_in_file, readable_out_file, complete_fasta_out_file,
             maxquant_csv_writer):
    # Write a multiple fasta file with proteis within the given motif range
    prot_seq = SeqRecord(
        Seq(prot.get_seq(protein_db)[:motif_len], generic_protein),
        id=prot.get_prot_name(),
        description=prot.get_descr())
    SeqIO.write(prot_seq, weblogo_in_file, "fasta")

    # Write the complete multiple fasta file with proteis
    prot_seq = SeqRecord(
        Seq(prot.get_seq(protein_db), generic_protein),
        id=prot.get_prot_name(),
        description=prot.get_descr())
    SeqIO.write(prot_seq, complete_fasta_out_file, "fasta")

    # Write the hit in a way a human can easily check it
    write_readable_result_hit(readable_out_file, prot, protein_db)

    # Write a hit in the same way MaxQuant outputs it's results
    write_maxquant_result_hit(maxquant_csv_writer, prot, protein_db)


def find_metap_activity(min_prob, motif_len, prots, protein_db, weblogo_in_path,
                        complete_fasta_out_path,
                        readable_out_path, maxquant_csv_writer):
    weblogo_in_file = open(weblogo_in_path, "a")
    readable_out_file = open(readable_out_path, "a")
    complete_fasta_out_file = open(complete_fasta_out_path, "a")

    matching_prots, matching_prot_hits = filter_and_sort_hits(min_prob, prots)

    if len(matching_prot_hits) > 0:

        for prot in matching_prot_hits:
            for pep in prot.get_peptides():

                # Locate the position of the peptide on the protein
                pep_pos = prot.get_seq(protein_db).find(pep.get_seq())

                if pep_pos == 0:
                    save_hit(motif_len, prot, protein_db, weblogo_in_file, readable_out_file,
                             complete_fasta_out_file,
                             maxquant_csv_writer)

    weblogo_in_file.close()
    readable_out_file.close()
    complete_fasta_out_file.close()


def main(ms_run_code, prot_prophet_xml, protein_db_path, min_prob, motif_len):
    protxml_name = prot_prophet_xml[prot_prophet_xml.rfind("\\"):prot_prophet_xml.find(".")]

    protein_db = SeqIO.index(protein_db_path, format='fasta')
    parse_results_out_paht = create_out_folder("protxml_parsing_out\\")
    weblogos_out_path = create_out_folder('weblogo_out\\')

    metap_weblogo_in_name = parse_results_out_paht + protxml_name + '_metap_weblogo_in.fasta'
    proteome_weblogo_in_name = parse_results_out_paht + protxml_name + '_original_db_weblogo_in.fasta'
    readable_out_name = parse_results_out_paht + protxml_name + '_human_readable_results.txt'
    metap_complete_fasta_out_name = parse_results_out_paht + protxml_name + '_metap_complete_sequences.fasta'
    proteome_without_metap_out_name = parse_results_out_paht + protxml_name + '_original_db_without_metap_weblogo_in.fasta'
    maxquant_out_name = parse_results_out_paht + protxml_name + ms_run_code + '_maxquant.csv'

    maxquant_out_file = open(maxquant_out_name, "w+")
    maquant_csv_writer = csv.writer(maxquant_out_file, quotechar=',', quoting=csv.QUOTE_MINIMAL)

    # Write the first row of the MaxQuant csv file
    # maquant_csv_writer.writerow(["Sequence", "N-terminal", "C-terminal", "Modifications", "Mass",
    #                              "Mass Fractional Part", "Protein Groups", "Proteins", "Unique (Groups)",
    #                              "Unique (Proteins)",
    #                              "Acetyl (Protein N-term)", "Oxidation (M)", "Missed cleavages", "Retention time",
    #                              "Calibrated retention time", "Charges", "PEP", "MS/MS scan number", "Raw file",
    #                              "Score", "Delta score", "Intensity", "Reverse", "Potential contaminant", "id",
    #                              "Protein group IDs", "Peptide ID", "Evidence IDs", "MS/MS IDs", "Best MS/MS",
    #                              "Oxidation (M) site IDs", "MS/MS Count"])

    prot_xml_parser = XinteractParser(parse_results_out_paht, prot_prophet_xml)

    prots = [prot for prot in prot_xml_parser.get_prot_groups()]
    weblogo_generator = WebLogoGenerator(weblogos_out_path)

    find_metap_activity(min_prob, motif_len, prots, protein_db, metap_weblogo_in_name,
                        metap_complete_fasta_out_name, readable_out_name, maquant_csv_writer)

    extract_metap_records_weblogo_in(protein_db_path, metap_weblogo_in_name, proteome_without_metap_out_name, motif_len)

    cleave_db_motif_len(protein_db_path, motif_len, proteome_weblogo_in_name)

    weblogo_generator.create_weblogo(metap_weblogo_in_name)
    weblogo_generator.create_weblogo(proteome_without_metap_out_name)
    weblogo_generator.create_weblogo(proteome_weblogo_in_name)
