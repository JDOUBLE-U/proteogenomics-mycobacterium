#!/bin/bash
#Erik Walinga, 30-03-2015

# this script extracts the DNA sequence between the two given positions
# from the given file.

pos1=${2}
pos2=${3}
curpos=0
target=""
go=0
while read line
do
  endpos=$((${curpos}+${#line}-1))
  if ((${pos1} > ${curpos} && ${pos1} < ${endpos})); then
    target+=${line:$((${pos1}-${curpos}-1)):$((${endpos}-${pos1}+2))}
    go=1
  fi
  if ((${pos2} > ${curpos} && ${pos2} < ${endpos})); then
    target+=${line:0:$((${pos2}-${curpos}))}
    go=0
  fi
  if ((${go} == 1 && ${curpos} > ${pos1} && ${endpos} < ${pos2})); then
    target+=${line}
  fi
  curpos=$((${curpos} + ${#line}))
done < ${1}

echo ${target}

