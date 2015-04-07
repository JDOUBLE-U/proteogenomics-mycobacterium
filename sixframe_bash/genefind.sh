#!/bin/bash
# Erik Walinga
# finds all genes in the genome through regex, determines the reading frames 
# of the genes and translates the gene to a protein sequence. The stopcodons
# can be defined by the user.

function translate {
  declare -A codontabel=([AAA]="K" [AAC]="N" [AAG]="K" [AAT]="N" [ACA]="T"
                         [ACC]="T" [ACG]="T" [ACT]="T" [AGA]="R" [AGC]="S"
                         [AGG]="R" [AGT]="S" [ATA]="I" [ATC]="I" [ATG]="M"
                         [ATT]="I" [CAA]="Q" [CAC]="H" [CAG]="Q" [CAT]="H"
                         [CCA]="P" [CCC]="P" [CCG]="P" [CCT]="P" [CGA]="R"
                         [CGC]="R" [CGG]="R" [CGT]="R" [CTA]="L" [CTC]="L"
                         [CTG]="L" [CTT]="L" [GAA]="E" [GAC]="D" [GAG]="E"
                         [GAT]="D" [GCA]="A" [GCC]="A" [GCG]="A" [GCT]="A"
                         [GGA]="G" [GGC]="G" [GGG]="G" [GGT]="G" [GTA]="V"
                         [GTC]="V" [GTG]="V" [GTT]="V" [TAC]="Y" 
                         [TAT]="Y" [TCA]="S" [TCC]="S" [TCG]="S"
                         [TCT]="S" [TGC]="C" [TGG]="W" [TGT]="C"
                         [TTA]="L" [TTC]="F" [TTG]="L" [TTT]="F")
  IFS=',' read -a stops <<< ${2}
  for codon in "${stops[@]}"
  do
    codontabel[${codon}]="*"
  done 
  prot=""
  for ((pos=0; pos<${#1}-2; pos+=3))
  do
    codon=${1:${pos}:3}
    if [[ -n "${codontabel[${codon}]}" ]]; then
      prot+=${codontabel[${codon}]}
    else
      frame+="X"
    fi
  done
  echo "$prot"
}

#  genoomfile, stops
echo "" > result.txt
echo ${1}
bla="ATG[ACGT]*?(("
IFS=',' read -a stops <<< ${1}
for codon in "${stops[@]:0:$((${#stops}-1))}"
  do
    bla+=${codon}"|"
  done
bla+="${stops[@]:$((${#stops}-1)):1}))"
matches=($(cat testfile.fasta | tr -d '\n' | tr -d '\r' | grep -boP ${bla} ))
count=1
for match in ${matches[@]}
do
  IFS=':' read -a dna <<< ${match}
  echo ${match}
  echo "gen ${count}, frame +$((${dna[0]}%3 + 1))" >> result.txt
  echo "$(translate ${dna[1]} ${1} )" >> result.txt
  
  count=$((${count} + 1))
done



