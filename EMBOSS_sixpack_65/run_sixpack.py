__author__ = 'Jeroen'

import os

"""
 The program takes the following steps:

 The nucleic acid sequence is read in.
 The three forward and three reverse translations are created.
 The name and description are written to the ouput display file.
 The reverse sense sequence is placed below the forward sequence.
 The forward translations are placed above the sequences.
 The reverse translation are placed below the sequences.
 Any ORFs that are longer than the specified minimum size are written to the output sequence file.
"""


def run_sixpack(genome_file, executable, on_platform, min_prot_length, codon_table):
    print("Running EMBOSS Sixpack")

    fasta_out = "..\\" + "analysis_pipeline\\sixpack_out" + genome_file[genome_file.rfind("\\"):-len(
        ".fasta")] + "_six_frame.fasta"
    # visual_out = fasta_out[:-len(".fasta")] + ".txt"

    # TODO perhaps add the "features format -fformat1" parameter for custom fasta descriptors?
    os.system(executable + " -sequence %s"
                           " -sformat1 fasta"
                           " -snucleotide1"
                           " -supper1"
                           " -table %i"
                           " -firstorf"
                           " -lastorf"
                           " -outseq %s"
                           " -osformat fasta"
                           " -osname %s"
                           " -reverse"
                           " -orfminsize %i"
                           " -number"
                           " -width 120"
                           " -length 0"
                           " -margin 10"
                           " -name"
                           " -description"
                           " -offset 1"
                           " -nohtml"
                           " -auto"
                           " -verbose"
                           " -mstart"  # Make optional?
              % (genome_file, codon_table, fasta_out, on_platform, min_prot_length))

    return fasta_out
