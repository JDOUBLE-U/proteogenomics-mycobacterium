__author__ = 'Jeroen'

import os
import sys

#
# def choose_min_prot_length():
#     print("What is the minimum protein length allowed?")
#     return int(input(">>> "))
#
#
# def main(genome):
#     out_file = genome[:-len("_genome.fasta")] + "proteome_six_frame.fasta"
#     visual_out = out_file[:-len(".fasta")] + ".txt"
#
#     on_platform = sys.platform
#     if on_platform == 'win32':
#         sixpack_command = "windows_sixpack.exe"
#     else:
#         sixpack_command = "./sixpack"
#
#     min_protein_len = choose_min_prot_length()
#
#     # TODO perhaps add the "features format -fformat1" parameter for custom fasta descriptors?
#     os.system(sixpack_command + " -sequence %s"
#                                 " -sformat1 fasta "
#                                 " -snucleotide1"
#                                 " -supper1 "
#                                 " -table %i"
#                                 " -firstorf "
#                                 " -lastorf "
#                                 " -outfile %s"
#                                 " -outseq %s"
#                                 " -osformat fasta "
#                                 " -osname %s"
#                                 " -reverse "
#                                 " -orfminsize %i"
#                                 " -number "
#                                 " -width 120 "
#                                 " -length 0 "
#                                 " -margin 10 "
#                                 " -name "
#                                 " -description "
#                                 " -offset 1 "
#                                 " -nohtml "
#                                 " -auto"
#                                 " -verbose"
#                                 " -mstart"  # Make optional?
#               % (genome, CODON_TABLE, visual_out, out_file, on_platform, min_protein_len))
#
#
# main("../" + "GitHub_test_files/Mycobacterium_marinum_genome.fasta")