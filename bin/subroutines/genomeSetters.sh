#!/bin/bash

##########################################################################
# Copyright 2017, Jelena Telenius (jelena.telenius@imm.ox.ac.uk)         #
#                                                                        #
# This file is part of CCseqBasic5 .                                     #
#                                                                        #
# CCseqBasic5 is free software: you can redistribute it and/or modify    #
# it under the terms of the MIT license.
#
#
#                                                                        #
# CCseqBasic5 is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# MIT license for more details.
#                                                                        #
# You should have received a copy of the MIT license
# along with CCseqBasic5.  
##########################################################################

setBOWTIEgenomeSizes(){
    
BowtieGenome="UNDETERMINED"

#-----------Genome-sizes-for-bowtie-commands----------------------------------------------  
    
if [ "${BOWTIE}" -eq 1 ]
then
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${GENOME}" ]; then
    BowtieGenome="${BOWTIE1[$g]}"
fi

done

fi

#------------------------------------------------

if [ "${BOWTIE}" -eq 2 ]
then
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${GENOME}" ]; then
    BowtieGenome="${BOWTIE2[$g]}"
fi

done  

fi  

#------------------------------------------------

# Check that it got set ..

if [ "${BowtieGenome}" == "UNDETERMINED" ]; then 
  echo "Genome build " ${GENOME} " is not supported - aborting !"  >&2
  exit 1 
fi    

# Check that at least one index file exists..

TEMPcount=$(($( ls -1 ${BowtieGenome}* | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]; then

  echo "Bowtie genome build for ${GENOME} : no index files ${BowtieGenome} found - aborting !"  >&2
  exit 1     
fi

echo
echo "Genome ${GENOME} .  Set BOWTIE index directory and basename : ${BowtieGenome}"

}

setGenomeFasta(){
    
GenomeFasta="UNDETERMINED"

#-----------Genome-sizes-for-bowtie-commands----------------------------------------------  

# echo "_${GENOME}_"
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${GENOME}" ]; then
    GenomeFasta="${WholeGenomeFASTA[$g]}"
fi

done  

#------------------------------------------------

# Check that it got set ..

if [ "${GenomeFasta}" == "UNDETERMINED" ]; then 
  echo "Genome build " ${GENOME} " is not supported -- aborting !"  >&2
  exit 1 
fi    

# Check that the index file exists..

if [ ! -s "${GenomeFasta}" ]; then

  echo "Whole genome fasta for ${GENOME} : file not found : ${GenomeFasta} - aborting !"  >&2
  exit 1     
fi

echo
echo "Genome ${GENOME} .  Set whole genome fasta file : ${GenomeFasta}"


}


setUCSCgenomeSizes(){
    
ucscBuild="UNDETERMINED"
    
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
# echo ${supportedGenomes[$g]}

if [ "${supportedGenomes[$g]}" == "${GENOME}" ]; then
    ucscBuild="${UCSC[$g]}"
fi

done 
    
if [ "${ucscBuild}" == "UNDETERMINED" ]; then 
  echo "Genome build " ${GENOME} " is not supported --- aborting !"  >&2
  exit 1 
fi

# Check that the file exists..
if [ ! -e "${ucscBuild}" ] || [ ! -r "${ucscBuild}" ] || [ ! -s "${ucscBuild}" ]; then
  echo "Genome build ${GENOME} file ${ucscBuild} not found or empty file - aborting !"  >&2
  exit 1     
fi

echo
echo "Genome ${GENOME} . Set UCSC genome sizes file : ${ucscBuild}"
echo

}



