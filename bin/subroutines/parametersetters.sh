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

setMparameter(){
   
mParameter=""

if [ "${CAPITAL_M}" -eq 0 ] ; then
    mParameter="-m ${LOWERCASE_M}"
else
    mParameter="-M ${CAPITAL_M}"
fi 
    
}

setParameters(){

#----------------------------------------------
# Listing current limitations, exiting if needed :

if [ "${LOWERCASE_M}" -ne 0 ] && [ "${CAPITAL_M}" -ne 0 ];
then
    printThis="Both -m and -M parameters cannot be set at the same time\nEXITING"
    printToLogFile
  exit 1
fi

#----------------------------------------------

if [ "${LOWERCASE_V}" -ne -1 ] && [ "${bowtie1MismatchBehavior}" != "" ]
then
    printThis="Bowtie1 does not allow setting -v with any other mismatch-reporting altering parameters ( --seedmms --seedlen --maqerr ) \nUse only -v, or (any) combination of --seedmms --seedlen --maqerr\nEXITING"
    printToLogFile
  exit 1
fi

if [ "${bowtie1MismatchBehavior}" != "" ]
then
    otherBowtie1Parameters="${otherBowtie1Parameters} ${bowtie1MismatchBehavior}"
fi

if [ "${LOWERCASE_V}" -ne -1 ]
then
    otherBowtie1Parameters="${otherBowtie1Parameters} -v ${LOWERCASE_V}"
fi

otherBowtie2Parameters="${otherBowtie2Parameters} ${bowtie2MismatchBehavior}"

#----------------------------------------------
#Setting the m and M parameters..

if [ "${LOWERCASE_M}" -ne 0 ] ;
then
   CAPITAL_M=0 
fi

if [ "${CAPITAL_M}" -ne 0 ];
then
   LOWERCASE_M=0
fi

if [ "${LOWERCASE_M}" -eq 0 ] && [ "${CAPITAL_M}" -eq 0 ];
then
    LOWERCASE_M=2
fi

#------------------------------------------------
# Custom adapter sequences..

if [ "${ADA31}" != "no"  ] || [ "${ADA32}" != "no" ] || [ "${ADA51}" != "no" ] || [ "${ADA52}" != "no" ] 
    then
    CUSTOMAD=1
fi

}

