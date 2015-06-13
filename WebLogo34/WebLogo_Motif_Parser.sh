#Adjust 'avium' in the proteome that you're using.
#Proteome file.
proteome_fasta="/Users/janwillem/Documents/BPSYS/M_avium_proteome.fasta"
proteome_one_line_fasta="/Users/janwillem/Documents/BPSYS/M_avium_proteome_one_line.fasta"
proteome_one_line_without_first_M="/Users/janwillem/Documents/BPSYS/M_avium_proteome_one_line_without_first_M.fasta"
proteome_one_line_without_first_M_five_amino="/Users/janwillem/Documents/BPSYS/M_avium_proteome_one_line_without_first_M_5_amino.fasta"
proteome_weblogo="/Users/janwillem/Documents/BPSYS/M_avium_proteome_5_amino_after_M_Weblogo.fasta"


#Per header and sequence one line
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${proteome_fasta} > ${proteome_one_line_fasta}


#Cleave all first Methionine's
sed 's/^[A-Z]\{1\}//' ${proteome_one_line_fasta} > ${proteome_one_line_without_first_M}



#The first 5 amino acids after the methionine.
perl -e 'while($id = <>){$seq = <>; $seq = substr($seq,0,5); print $id.$seq."\n"}' ${proteome_one_line_without_first_M} > ${proteome_one_line_without_first_M_five_amino}


#Optional for use in Weblogo.
sed '/^>/d' ${proteome_one_line_without_first_M_five_amino} > ${proteome_weblogo}
