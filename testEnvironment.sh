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


function finish {
echo
echo
echo "###########################################"
echo
echo "Finished with tests !"
echo
echo "Check that you got no errors in above listings !"
echo
echo "###########################################"
echo
echo
}
trap finish EXIT

# Set the default return value
exitCode=0

echo
echo "This is test for CCseqBasic configuration setup ! "
echo
echo "( For automated testing : Return value of the script is '0' if all clear or only warnings, and '1' if fatal errors encountered. )"
echo
sleep 2
echo "Running test script $0"
echo
echo "###########################################"
echo
echo "1) Testing that the UNIX basic tools (sed, awk, etc) are found"
echo "2) Testing that the needed scripts are found 'near by' the main script "
echo "3) Setting the environment - running the conf/config.sh , to set the user-defined parameters"
echo "4) Listing the set allowed genomes, and blacklist files"
echo "5) Testing that all toolkits (bowtie etc) are found in the user-defined locations"
echo "6) Testing that the user-defined public server exists"
echo
sleep 5

##########################################################################

echo "###########################################"
echo
echo "1) Testing that the UNIX basic tools (sed, awk, grep, et cetera) are found"
echo

echo "Calling sed .."
echo
sed --version | head -n 1
sed --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Calling awk .."
echo
awk --version | head -n 1
awk --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Calling grep .."
echo
grep --version | head -n 1
grep --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Calling GNU coreutils .."
echo

cat    --version | head -n 1
cat    --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
chmod  --version | head -n 1
chmod  --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
cp     --version | head -n 1
cp     --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
cut    --version | head -n 1
cut    --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
date   --version | head -n 1
date   --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
ln     --version | head -n 1
ln     --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
ls     --version | head -n 1
ls     --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
mkdir  --version | head -n 1
mkdir  --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
mv     --version | head -n 1
mv     --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
paste  --version | head -n 1
paste  --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
rm     --version | head -n 1
rm     --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
rmdir  --version | head -n 1
rmdir  --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
sort   --version | head -n 1
sort   --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
tail   --version | head -n 1
tail   --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
tr     --version | head -n 1
tr     --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
unlink --version | head -n 1
unlink --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
dirname --version | head -n 1
dirname --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 4

echo "Calling 'which'  .."
echo
which --version | head -n 1
which --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

diffLoadFailed=0
echo "Calling 'diff' (optional - needed only in this tester script) .."
echo
diff --version | head -n 1
diff --version >> /dev/null
diffLoadFailed=$?
echo

sleep 2

echo "Calling 'hostname' (optional - it is only used to print out the name of the computer) .."
echo
hostname --version 2>&1 
echo

sleep 2

echo "Calling 'module' (optional - only needed if you set your  conf/loadNeededTools.sh   to use the module environment) .."
echo
module --version 2>&1 | head -n 2
echo

sleep 3

##########################################################################

# Test that the script files exist ..

echo "###########################################"
echo
echo "2) Testing that the needed scripts are found 'near by' the main script .."
echo

sleep 2

PipeTopPath="$( dirname $0 )"
dirname $0 >> /dev/null
exitCode=$(( ${exitCode} + $? ))

# From where to call the CONFIGURATION script..

confFolder="${PipeTopPath}/conf"
mainScriptFolder="${PipeTopPath}/bin/runscripts"
helperScriptFolder="${PipeTopPath}/bin/subroutines"

echo
echo "This is where they should be ( will soon see if they actually are there ) :"
echo
echo "PipeTopPath        ${PipeTopPath}"
echo "confFolder         ${confFolder}"
echo "mainScriptFolder   ${mainScriptFolder}"
echo "helperScriptFolder ${helperScriptFolder}"
echo

sleep 4


scriptFilesMissing=0

# Check that it can find the scripts .. ( not checking all - believing that if these exist, the rest exist too )

echo
echo "Master script and its tester script :"
echo
ls ${PipeTopPath}/CCseqBasic5.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${PipeTopPath}/testEnvironment.sh 
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3
echo "Main analysis scripts :"
echo
ls ${mainScriptFolder}/QC_and_Trimming.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/analyseMappedReads.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3
echo "In-silico cutter scripts - for each Restriction Enzyme :"
echo
ls ${mainScriptFolder}/dpnIIcutGenome4.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/dpnIIcutReads4.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))

ls ${mainScriptFolder}/nlaIIIcutGenome4.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/nlaIIIcutReads4.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))

# ls ${mainScriptFolder}/hindIIIcutGenome4.pl
# scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
# ls ${mainScriptFolder}/hindIIIcutReads4.pl
# scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))

echo "UCSC downloaded tools for visualisation (downloaded 06May2014 from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ ) :"
echo
ls ${confFolder}/ucsctools/bedToBigBed
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${confFolder}/ucsctools/wigToBigWig
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3

echo
sleep 3
echo "Helper subroutines :"
echo

ls ${helperScriptFolder}/usageAndVersion.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/testers_and_loggers.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/parametersetters.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/blacklistSetters.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/genomeSetters.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/runtools.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/hubbers.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/cleaners.sh 
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))

ls ${mainScriptFolder}/intersectappend.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/fastq_scores_bowtie1.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/fastq_scores_bowtie2.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))

echo
sleep 3
echo "Scripts for read filtering :"
echo
ls ${mainScriptFolder}/filterArtifactMappers/filter.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/filterArtifactMappers/1_blat.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/filterArtifactMappers/2_psl_parser.pl
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3
echo "Scripts for drawing summary figure :"
echo
ls ${mainScriptFolder}/drawFigure/countsFromCCanalyserOutput.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/drawFigure/generatePercentages.py
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${mainScriptFolder}/drawFigure/drawFigure.py
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3
echo "Configuration setters :"
echo
ls ${confFolder}/genomeBuildSetup.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${confFolder}/loadNeededTools.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${confFolder}/serverAddressAndPublicDiskSetup.sh
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3
echo "Configuration tester helpers :"
echo
ls ${helperScriptFolder}/validateSetup/g.txt
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/validateSetup/l.txt
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
ls ${helperScriptFolder}/validateSetup/s.txt
scriptFilesMissing=$(( ${scriptFilesMissing} + $? ))
echo
sleep 3

if [ "${scriptFilesMissing}" -ne 0 ]
then
echo
echo "###########################################"
echo
echo "ERROR !   The scripts CCseqBasic5.sh is dependent on, are not found in their correct relative paths !"
echo "          Maybe your tar archive was corrupted, or you meddled with the folder structure after unpacking ?"
echo
echo "###########################################"
echo
echo "This is what you SHOULD see if you run 'tree' command in your CCseqBasic folder :"
echo
echo ' |-- CCseqBasic5.sh'
echo ' |-- testEnvironment.sh'
echo ' |'
echo ' `-- bin'
echo '     |-- runscripts'
echo '     |   |-- QC_and_Trimming.sh analyseMappedReads.pl '
echo '     |   |-- fastq_scores_bowtie1.pl fastq_scores_bowtie2.pl intersectappend.pl'
echo '     |   |-- dpnIIcutGenome4.pl dpnIIcutReads4.pl nlaIIIcutGenome4.pl nlaIIIcutReads4.pl'
# echo '     |   |   hindIIIcutGenome4.pl hindIIIcutReads4.pl'
echo '     |   |-- filterArtifactMappers'
echo '     |   |    `-- filter.sh 1_blat.sh 2_psl_parser.pl '
echo '     |   `-- drawFigure'
echo '     |        `-- countsFromCCanalyserOutput.sh drawFigure.py generatePercentages.py '
echo '     |-- subroutines'
echo '     |   |-- hubbers.sh runtools.sh '
echo '     |   |-- blacklistSetters.sh genomeSetters.sh cleaners.sh '
echo '     |   |-- testers_and_loggers.sh parameterSetters.sh usageAndVersion.sh'
echo '     |   `-- validateSetup'
echo '     |       `-- g.txt l.txt s.txt'
echo '     `-- perlHelpers'
echo '         |-- data2gff.pl fastq_scores_bowtie1/2.pl reverse_seq.pl '
echo '         `-- sam2fastq.pl trim1base3prime.pl windowingScript.pl'
echo ''
echo '`-- conf'
echo '    |-- BLACKLIST'
echo '    |   `-- hg18.bed hg19.bed mm10.bed mm9.bed'
echo '    |-- genomeBuildSetup.sh'
echo '    |-- loadNeededTools.sh'
echo '    |-- serverAddressAndPublicDiskSetup.sh'
echo '    |-- ucsctools'
echo '    |   |-- wigToBigWig '
echo '    |   `-- bedToBigBed'
echo '    `-- UCSCgenomeSizes'
echo '        |-- danRer10/7.chrom.sizes dm3.chrom.sizes galGal4.chrom.sizes'
echo '        `-- hg18/19/38.chrom.sizes mm10/9.chrom.sizes'
echo ''

sleep 4

# Return the value : 0 if only warnings, 1 if fatal problems.
exit 1

fi

exitCode=$(( ${exitCode} + ${scriptFilesMissing} ))
sleep 5

##########################################################################

# Test that user has made at least SOME changes to them (otherwise they are running with the WIMM configuration .. )


echo
echo "###########################################"
echo
echo "3) Setting the environment - running the conf/(setupscripts).sh , to set the user-defined parameters"
echo

sleep 6


# We have 3 files - if all of them are in WIMM setup, we end up with "0" as the value in the end ..
setupMade=3
genomeSetupMade=1
toolsSetupMade=1
serverSetupMade=1

TEMPcount=$(($( diff ${helperScriptFolder}/validateSetup/g.txt ${confFolder}/genomeBuildSetup.sh | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]
then
setupMade=$((${setupMade}-1))
genomeSetupMade=0
echo
echo "WARNING ! It seems you haven't set up your Bowtie Genome indices !"
echo "          Add your Bowtie indices to this file : "
echo "          ${confFolder}/genomeBuildSetup.sh "
echo
sleep 6
fi

TEMPcount=$(($( diff ${helperScriptFolder}/validateSetup/l.txt ${confFolder}/loadNeededTools.sh | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]
then
setupMade=$((${setupMade}-1))
toolsSetupMade=0
echo
echo "WARNING ! It seems you haven't set up the loading of your Needed Toolkits !"
echo "          Add your toolkit paths to this file : "
echo "          ${confFolder}/loadNeededTools.sh "
echo
echo "NOTE !!   You need to edit this file ALSO if you want to disable loading the toolkits via the above script."
echo "          To disable the loading of the tools, set : setToolLocations=0 "
echo
sleep 8
fi

TEMPcount=$(($( diff ${helperScriptFolder}/validateSetup/s.txt ${confFolder}/serverAddressAndPublicDiskSetup.sh | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]
then
setupMade=$((${setupMade}-1))
serverSetupMade=0
echo
echo "WARNING ! It seems you haven't set up your Server address and Public Disk Area !"
echo "          Add your Server address to this file : "
echo "          ${confFolder}/serverAddressAndPublicDiskSetup.sh "
echo
sleep 4
fi

# Only continue to the rest of the script, if there is some changes in the above listings ..

if [ "${setupMade}" -eq 0 ]
then
echo 
echo
echo "Could not finish testing, as you hadn't set up your environment !"
echo
echo "Set up your files according to instructions in :"
echo "http://sara.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/instructionsGeneral.html"
echo
sleep 4

exitCode=1

fi

if [ "${exitCode}" -gt 0 ]
then
exit 1
else
exit 0
fi

##########################################################################

# These have been checked earlier. Should exist now.
. ${confFolder}/loadNeededTools.sh
. ${confFolder}/genomeBuildSetup.sh
. ${confFolder}/serverAddressAndPublicDiskSetup.sh

##########################################################################

if [ "${genomeSetupMade}" -eq 1 ]; then 

setGenomeLocations 1>/dev/null

supportedGenomes=()
BOWTIE1=()
BOWTIE2=()
UCSC=()
genomesWhichHaveBlacklist=()
BLACKLIST=()

echo
sleep 4

echo "###########################################"
echo
echo "4) Listing the setup of allowed genomes, and blacklist files"
echo


echo "Supported genomes : "
echo
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do    
 echo "${supportedGenomes[$g]}"
done

echo
sleep 2

##########################################################################
echo
echo "Bowtie 1 indices : "
echo
echo -e "GENOME\tBOWTIE 1 index"
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do    

 echo -en "${supportedGenomes[$g]}\t${BOWTIE1[$g]}"

TEMPcount=$(($( ls -1 ${BOWTIE1[$g]}* | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]; then
    echo -e "\tINDICES DO NOT EXIST in the given location !!"
    exitCode=$(( ${exitCode} +1 ))
 else
    echo ""
 fi
 
done

echo
sleep 2

##########################################################################
echo
echo "Bowtie 2 indices : "
echo
echo -e "GENOME\tBOWTIE 2 index"
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do    

 echo -en "${supportedGenomes[$g]}\t${BOWTIE2[$g]}"

TEMPcount=$(($( ls -1 ${BOWTIE2[$g]}* | grep -c "" )))

if [ "${TEMPcount}" -eq 0 ]; then
    echo -e "\tINDICES DO NOT EXIST in the given location !!"
    exitCode=$(( ${exitCode} +1 ))
 else
    echo ""
 fi
 
done

echo
sleep 2

##########################################################################
echo
echo "UCSC genome size files : "
echo
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do
    
 echo -en "${supportedGenomes[$g]}\t${UCSC[$g]}"
    
 if [ ! -e "${UCSC[$g]}" ] || [ ! -r "${UCSC[$g]}" ] || [ ! -s "${UCSC[$g]}" ]; then
    echo -e "\tFILE DOES NOT EXIST in the given location !!"
    exitCode=$(( ${exitCode} +1 ))
 else
    echo ""
 fi

done

echo
sleep 5

##########################################################################
echo
echo "Genomes which have blacklist : "
echo
for g in $( seq 0 $((${#genomesWhichHaveBlacklist[@]}-1)) ); do
    
 echo -en "${supportedGenomes[$g]}\t${BLACKLIST[$g]}"

 if [ ! -e "${BLACKLIST[$g]}" ] || [ ! -r "${BLACKLIST[$g]}" ] || [ ! -s "${BLACKLIST[$g]}" ]; then
    echo -e "\tFILE DOES NOT EXIST in the given location !!"
    exitCode=$(( ${exitCode} +1 ))
 else
    echo ""
 fi
 
done

echo
sleep 2

else
 
echo "###########################################"
echo
echo "4) Skipping genome setup testing ( as genomBuildSetup.sh was not filled with genome locations )"
echo   

sleep 3
 
fi

##########################################################################

if [ "${toolsSetupMade}" -eq 1 ]; then

echo "###########################################"
echo
echo "5) Testing that all toolkits (bowtie etc) are found in the user-defined locations"
echo

setPathsForPipe 1>/dev/null

echo "Ucsctools .."
echo
wigToBigWig 2>&1 | head -n 1
which wigToBigWig >> /dev/null
exitCode=$(( ${exitCode} + $? ))
bedToBigBed      2>&1 | head -n 1
which bedToBigBed      >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 3

echo "Samtools .."
echo
samtools 2>&1 | head -n 3 | grep -v "^\s*$"
which samtools >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Bedtools .."
echo
bedtools --version
bedtools --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))

echo

sleep 2

echo "Bowtie 1 .."
echo
bowtie --version | head -n 5
bowtie --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Bowtie 2 .."
echo
bowtie2 --version | head -n 5
bowtie2 --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Flash .."
echo
flash --version | head -n 1
flash --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Trim_galore .."
echo
trim_galore --version | sed 's/^\s*//' | grep -v "^\s*$"
trim_galore --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Cutadapt .."
echo
cutadapt --version
cutadapt --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Flash .."
echo
flash --version
flash --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Blat .."
echo
blat 2>&1 | head -n 1
which blat >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2


echo "FastqQC .. "
echo "(series 0.10.x is NOT supported. Check below, that you have 0.11.x )"
echo
fastqc --version
fastqc --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 4

echo "Perl .."
echo
perl --version | head -n 5 | grep -v "^\s*$"
perl --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 2

echo "Python .."
echo "(Python 3.* is NOT supported. Check below, that you have Python 2.*)"
echo
python --version
python --version >> /dev/null
exitCode=$(( ${exitCode} + $? ))
echo

sleep 4

echo "Matplotlib - within python .."
echo "(matplotlib 1.3.* or older are NOT supported. Check below, that you have at least matplotlib 1.4.1 )"
echo
python -c "import matplotlib as mpl; print mpl.__version__" 
exitCode=$(( ${exitCode} + $? ))
echo

sleep 4

else
 
echo "###########################################"
echo
echo "5) Skipping toolkit availability testing ( as loadNeededTools.sh was not set up to find the proper tool locations )"
echo   

sleep 3
 
fi

##########################################################################

if [ "${serverSetupMade}" -eq 1 ]; then

echo "###########################################"
echo
echo "6) Testing that the user-defined public server exists"
echo

echo
echo "Public area settings : "
echo
echo "SERVERTYPE ${SERVERTYPE}"
echo "SERVERADDRESS ${SERVERADDRESS}"
echo "ADDtoPUBLICFILEPATH ${ADDtoPUBLICFILEPATH}"
echo "REMOVEfromPUBLICFILEPATH ${REMOVEfromPUBLICFILEPATH}"
echo "tobeREPLACEDinPUBLICFILEPATH ${tobeREPLACEDinPUBLICFILEPATH}"
echo "REPLACEwithThisInPUBLICFILEPATH ${REPLACEwithThisInPUBLICFILEPATH}"
echo

sleep 3

echo
echo "--------------------------------------------------------------------------"
echo "Testing the existence of public server ${SERVERTYPE}://${SERVERADDRESS} "
echo
echo "Any curl errors such as 'cannot resolve host' in the listing below, mean that the server is not available (possible typos in server name above ? ) "
echo
echo "curl --head ${SERVERTYPE}://${SERVERADDRESS} "
echo

curl --head ${SERVERTYPE}://${SERVERADDRESS}
exitCode=$(( ${exitCode} + $? ))

sleep 5

else
 
echo "###########################################"
echo
echo "6) Skipping public server existence testing ( as serverAddressAndPublicDiskSetup.sh was not filled with the server address )"
echo   

sleep 3
 
fi

# ##################################

if [ "${setupMade}" -ne 3 ]
then
echo 
echo
echo "ERROR : Could not finish testing, as you hadn't set up all of your environment !"
echo
echo "One or more of these need to be still set up :"
echo "1) genome locations"
echo "2) tool locations"
echo "3) server address"
echo 
echo "Set up your files according to instructions in :"
echo "http://sara.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/instructionsGeneral.html"
echo

exitCode=1

fi

# Return the value : 0 if only warnings, 1 if fatal problems.
if [ "${exitCode}" -gt 0 ]
then
exit 1
else
exit 0
fi

