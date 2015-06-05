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

    fasta_out = "..\\" + "analysis_pipeline\\sixpack_out" + genome_file[genome_file.rfind("\\"):-len(
        ".fasta")] + "_six_frame.fasta"

    # TODO perhaps add the "features format -fformat1" parameter for custom fasta descriptors?
    os.system(executable + " -sequence %s"
                           " -sformat1 fasta"  # Input sequence format
                           " -snucleotide1"  # Sequence is nucleotide
                           " -supper1"  # Make upper case
                           " -table %i"  # Genetics code used for the translation
                           " -firstorf"  # Count the beginning of a sequence as a possible ORF, even if it's inferior to the minimal ORF size.
                           " -lastorf"  # Count the end of a sequence as a possible ORF, even if it's not finishing with a STOP, or inferior to the minimal ORF size.
                           " -outseq %s"
                           " -osformat fasta"  # Output seq format
                           " -osname %s"  # Base file name
                           " -reverse"  # Display also the translation of the DNA sequence in the 3 reverse frames
                           " -orfminsize %i"  # Minimum size of Open Reading Frames (ORFs) to display in the translations.
                           " -auto"  # Turn off prompts
                           " -mstart"  # Displays only ORFs starting with an M
              % (genome_file, codon_table, fasta_out, on_platform, min_prot_length))

    return fasta_out
