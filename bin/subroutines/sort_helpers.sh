#!/bin/bash

##########################################################################
# Copyright 2017, Jelena Telenius (jelena.telenius@imm.ox.ac.uk)         #
#                                                                        #
# This file is part of CCseqBasic5 .                                     #
#                                                                        #
# CCseqBasic5 is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# CCseqBasic5 is distributed in the hope that it will be useful,         #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with CCseqBasic5.  If not, see <http://www.gnu.org/licenses/>.   #
##########################################################################

sortIn1E6bunches(){
  # needs these to be set :
  # thisIsWhereIam=$( pwd )
  # sortParams="-k1,1 -k2,2n"  or sortParams="-n" etc
  # input in intoSorting.txt
  # outputs TEMPsortedMerged.txt
  
  rm -f preSorting.txt
  cp intoSorting.txt preSorting.txt
  
   # Here we make the files in 1 000 000 reads bunches
    
   head -n  1000000 preSorting.txt > forSortingfile1.txt
   tail -n +1000001 preSorting.txt > temp.txt
   mv -f temp.txt preSorting.txt
   
   TEMPPcounter=1
   while [ -s preSorting.txt ]
   do
        TEMPPcounter=$(( ${TEMPPcounter}+1 ))
        head -n  1000000 preSorting.txt > forSortingfile${TEMPPcounter}.txt
        tail -n +1000001 preSorting.txt > temp.txt
        mv -f temp.txt preSorting.txt
        
   done
   rm -f preSorting.txt
   
   echo "made this many files for sorting (each upto 1 million reads) :"  >&2
   # ls -lht | grep forSortingfile >&2
   ls forSortingfile* >&2
   
   for file in forSortingfile*
   do
       newwwName=$( echo ${file} | sed 's/forSortingfile/preSortedfile/' )
       # echo "sort -S ${BOWTIEMEMORY}M ${sortParams} -T ${thisIsWhereIam} ${file} > ${newwwName}" >&2
       sort -S ${BOWTIEMEMORY}M ${sortParams} -T ${thisIsWhereIam} ${file} > ${newwwName}
       rm -f ${file}
       
   done
   
   sort -m -S ${BOWTIEMEMORY}M ${sortParams} -T ${thisIsWhereIam} preSortedfile* > TEMPsortedMerged.txt
   rm -f preSortedfile*
   
}

sortResultTester(){
    
     # Check if all went well.. 

    countBefore=$(($( cat intoSorting.txt | grep -c "" )))
    countAfter=$(($( cat TEMPsortedMerged.txt | grep -c "" )))
    
    if [ "${countBefore}" -ne "${countAfter}" ]
    then
    
    echo "BIG TIME Sorting FAILED. " >&2
    echo "Original file had ${countBefore} data lines - sorted file had only ${countAfter} lines." >&2
    echo "EXITING!!" >&2
    exit 1
    
    else

    echo "BIG TIME Sorting SUCCEEDED ! " >&2

    fi   
    
}

sortResultInfo(){
    
    # Check if all went well.. 

    countBefore=$(($( cat intoSorting.txt | grep -c "" )))
    countAfter=$(($( cat TEMPsortedMerged.txt | grep -c "" )))
    
    if [ "${countBefore}" -ne "${countAfter}" ]
    then
    
    echo "BIG TIME Sorting FAILED." >&2
    echo "Original file had ${countBefore} data lines - sorted file had only ${countAfter} lines." >&2
    echo "CONTINUING - but may produce nonsense data!!" >&2
    
    else

    echo "BIG TIME Sorting SUCCEEDED ! " >&2

    fi
}

