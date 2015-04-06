#!/bin/bash
#Erik Walinga, 130-03-2015

# Dit script maakt een six-frame translatie van een door de gebruiker
# ingevoerde streng DNA. Het resultaat wordt op het scherm getoond, in delen
# van 50 basen per regel.

#Deze functie stelt de complementaire streng op van de ingevoerde DNA streng en returnt deze.
#input:
# $1: een DNA streng
#
#output:
# De complementaire streng.
function complementair {
  comp_seq=""
  for ((pos=${#1}-1; pos>=0; pos--))
  do
    case ${1:${pos}:1} in
     "A")
	comp_seq+="T";;
     "T")
	comp_seq+="A";;
     "C")
	comp_seq+="G";;
     "G")
	comp_seq+="C";;
     *)
        comp_seq+="X";;
     esac
  done
  echo ${comp_seq}
}

#Deze functie stelt een reading frame op. Het frame dat opgesteld wordt is
#afhankelijk van de parameters die aan de functie mee worden gegeven.
#input:
# $1: de te gebruiken DNA streng.
# $2: het nummer van het frame.
# $3: het aantal spaties dat vooraf gaat aan het eerste aminozuur.
#output:
# Het reading frame.
function leesframe {
  frame=${3}
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
  IFS=',' read -a stops <<< ${4}
  for codon in "${stops[@]}"
  do
    codontabel[${codon}]="*"
  done
  for ((pos=${2}; pos<${#1}-2; pos+=3))
  do
    codon=${1:${pos}:3}
    if [[ -n "${codontabel[${codon}]}" ]]; then
      frame+=${codontabel[${codon}]}
      frame+="  "
    else
      frame+="X  "
    fi
  done
  echo "$frame"
}

#Deze functie zorgt ervoor een deel van de output op het scherm getoond wordt.
#input:
#  $1: De beginpositie van het te tonen deel in de arrays en strings
function printer {
echo "+3 ${plusarray[3]:${1}:50}"
echo "+2 ${plusarray[2]:${1}:50}"
echo "+1 ${plusarray[1]:${1}:50}"
echo "   ${inDNA:${1}:50}"
echo "   ${compDNA:${1}:50}"
echo "-1 ${minarray[1]:${1}:50}"
echo "-2 ${minarray[2]:${1}:50}"
echo "-3 ${minarray[3]:${1}:50}"
echo ""
}

#dna, stops
inDNA=${1}
compDNA=$(complementair ${inDNA})
plusarray[1]=$(leesframe ${inDNA} 0 "" ${2} )
plusarray[2]=$(leesframe ${inDNA} 1 " " ${2} )
plusarray[3]=$(leesframe ${inDNA} 2 "  " ${2} )
minarray[1]=$(leesframe ${compDNA} 0 "" ${2} )
minarray[2]=$(leesframe ${compDNA} 1 " " ${2} )
minarray[3]=$(leesframe ${compDNA} 2 "  " ${2} )
for ((num=0; num<=${#inDNA}/50; num++))
do
  printer $((${num}*50))
done
