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


#This python-script makes sequence logo's.
python "/home/BPSYS/s1071701/WebLogo.py" ${seq_file_1} ${seqLogo}
