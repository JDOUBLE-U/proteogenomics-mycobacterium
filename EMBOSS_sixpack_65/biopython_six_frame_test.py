from Bio import SeqIO

__author__ = 'Jeroen'

out_file = open("test", "w+")
handle = open("C:/Users/Jeroen/IdeaProjects/BPSYS/Mycobacterium_marinum_genome.fasta", "rU")
# handle = open("C:/Users/Jeroen/IdeaProjects/BPSYS/GCA_000018345.1_ASM1834v1_genomic2.fna", "rU")
#
# # for record in SeqIO.parse(handle, "fasta"):
# coding_dna_frame_1 = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", generic_dna)
# coding_dna_frame_2 = coding_dna_frame_1[1:]
# coding_dna_frame_3 = coding_dna_frame_1[2:]
#
# print(str(coding_dna_frame_1))
# print(str(coding_dna_frame_2))
# print(str(coding_dna_frame_3))
#
# print(coding_dna_frame_1.translate(11, cds=True))
# print(coding_dna_frame_2.translate(11, cds=True))
# print(coding_dna_frame_3.translate(11, cds=True))

# for frames in range(3):

for record in SeqIO.parse(handle, "fasta"):
    # print(record.seq)
    # print(record.seq.reverse_complement())

    for strand, base in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(4):

            length = 3 * ((len(record) - frame) // 3)  # Multiple of three
            for prot in base[frame:frame + length].translate(11).split("*"):
                if len(prot) >= 5:
                    print("%s - length %i, strand %i, frame %i" % (prot, len(prot), strand, frame))