#!/bin/bash

# Additional information:
# =======================
#
# Just type /bin/bash proteo_mycoT_pipeline.sh and the pipeline starts.
#

# Show usage information:
if [ "$1" == "--h" ] || [ "$1" == "--help" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ]
then
	echo "" 
	echo "Way of usage:" 
	echo "" 
	echo "The user can use this script for starting the pipeline."
	echo ""
	echo "Example:" 
	echo ""
	echo "/bin/bash proteo_mycoT_pipeline.sh"
	echo ""
	
	exit  
fi

#fasta-input
seq_file_1="/home/BPSYS/motifs_test.fasta"
#Seq-jpeg
seqLogo="/home/BPSYS/seqLogo.jpeg"

#This command is for using Comet. This program analyze the mzXML-files.
/home/janwillem/Desktop/comet_binaries_2015011/comet.2015011.linux.exe MM.mzXML -Pcomet.params.new -DMycobacterium_marinum_strain_ATCC_BAA_535.fasta

#This command is PeptideProphet
runprophet interact.htm

#This command is ProteinProphet.
#usage:	ProteinProphet.pl '<interact pep prob html file1><interact pep prob html file2>....' <outfile> (ICAT) (GLYC) (XPRESS) (ASAP_PROPHET) (ACCURACY) (ASAP) (REFRESH) (DELUDE) (NOOCCAM)
perl ProteinProphet.pl



#This python-script makes sequence logo's.
python "/home/BPSYS/s1071701/WebLogo.py" ${seq_file_1} ${seqLogo}

