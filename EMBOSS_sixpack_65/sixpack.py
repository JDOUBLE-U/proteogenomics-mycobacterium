__author__ = 'Jeroen'

import os

"""
 The program takes the following steps:

 The nucleic acid sequence is read in.
 The required genetic code is read in from the EGC* data files.
 The three forward and three reverse translations are created.
 The name and description are written to the ouput display file.
 Any required regions to be changed to upper case are changed.
 Any required regions to be highlighted in HTML colour tags are changed.
 The reverse sense sequence is placed below the forward sequence.
 The forward translations are placed above the sequences.
 The reverse translation are placed below the sequences.
 The display is written out, split at the ends of lines.
 Any ORFs that are longer than the specified minimum size are written to the output sequence file.
"""

CODON_TABLE = 11  # Bacterial codon table (by NCBI standards)


def choose_min_prot_length():
    print("What is the minimum protein length allowed? "
          "Keep in mind the n-terminus will optionally be cleaved after this step.")
    return int(input(">>> "))


# TODO add option for mandatory M on n-terminus
def choose_mandatory_M():
    print("Make start codons code for Methionine by default? [y/n]")
    if input(">>> ").lower() == "y":
        return " -mstart"
    else:
        return ""


def run_sixpack(genome_file, executable, on_platform):
    out_file = "../" + "Comet_analysis_pipeline/sixpack_out" + genome_file[genome_file.rfind("/"):-len(
        ".fasta")] + "proteome_six_frame.fasta"
    visual_out = out_file[:-len(".fasta")] + ".txt"

    print("EMBOSS Sixpack")
    min_protein_len = choose_min_prot_length()

    # TODO perhaps add the "features format -fformat1" parameter for custom fasta descriptors?
    os.system(executable + " -sequence %s"
                                " -sformat1 fasta "
                                " -snucleotide1"
                                " -supper1 "
                                " -table %i"
                                " -firstorf "
                                " -lastorf "
                                " -outfile %s"
                                " -outseq %s"
                                " -osformat fasta "
                                " -osname %s"
                                " -reverse "
                                " -orfminsize %i"
                                " -number "
                                " -width 120 "
                                " -length 0 "
                                " -margin 10 "
                                " -name "
                                " -description "
                                " -offset 1 "
                                " -nohtml "
                                " -auto"
                                " -verbose"
                                " -mstart"  # Make optional?
              % (genome_file, CODON_TABLE, visual_out, out_file, on_platform, min_protein_len))

    return out_file