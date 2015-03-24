#!/usr/bin/python

"""
Jeroen Merks
Jan Willem Wijnands
24-3-2015
Python WebLogo
Example: python WebLogo.py seq.fasta Logo.jpeg
"""

#Import all the modules
from weblogolib import *
import sys
import os

seq_file_1 = sys.argv[1]
logo_file = sys.argv[2]

sequence = read_seq_data(open(seq_file_1))
data = LogoData.from_seqs(sequence)
logo_options = LogoOptions(logo_start=1, logo_end=8)
logo_options.resolution = 300
logo_options.fineprint = "LUMC"
logo_options.logo_title = "MatAP activity"

my_Format = LogoFormat(data, logo_options)
my_SeqLogo = jpeg_formatter(data, my_Format)
Seq_Logo = open(logo_file, 'wb')
Seq_Logo.write(my_SeqLogo)
Seq_Logo.close()
