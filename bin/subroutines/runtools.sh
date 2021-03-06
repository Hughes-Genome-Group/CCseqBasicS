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

################################################################
# Running whole genome blacklist dpnII generation

generateReBlacklist(){
    
printThis="Preparing 'too far from RE cut sites' blacklist file for ${REenzyme} cut genome) .."
printToLogFile
    
rm -f genome_${REenzyme}_blacklist.bed

# Plus/minus 300 bases to both directions
cat ${fullPathDpnGenome} | sed 's/:/\t/' | sed 's/-/\t/' | awk '{if(($3-$2)>(2*'${ampliconSize}')){print "chr"$1"\t"$2+'${ampliconSize}'"\t"$3-'${ampliconSize}'}}' > genome_${REenzyme}_blacklist.bed

testedFile="genome_${REenzyme}_blacklist.bed"
doTempFileTesting

doQuotaTesting

fullPathDpnBlacklist=$(pwd)"/genome_${REenzyme}_blacklist.bed"


}

generateReDigest(){

################################################################
# Running whole genome fasta dpnII digestion..

rm -f genome_${REenzyme}_coordinates.txt

if [ -s ${CaptureDigestPath}/${GENOME}.txt ] 
then
    
ln -s ${CaptureDigestPath}/${GENOME}.txt genome_${REenzyme}_coordinates.txt
    
else
    
    
# Running the digestion ..
# dpnIIcutGenome.pl
# nlaIIIcutGenome.pl   

printThis="Running whole genome fasta ${REenzyme} digestion.."
printToLogFile

printThis="perl ${RunScriptsPath}/${REenzyme}cutGenome4.pl ${GenomeFasta}"
printToLogFile

perl ${RunScriptsPath}/${REenzyme}cutGenome4.pl "${GenomeFasta}"

testedFile="genome_${REenzyme}_coordinates.txt"
doTempFileTesting

doQuotaTesting

fi

ls -lht

dpnGenomeName=$( echo "${GenomeFasta}" | sed 's/.*\///' | sed 's/\..*//' )
# output file :
# ${GenomeFasta}_dpnII_coordinates.txt

fullPathDpnGenome=$(pwd)"/genome_${REenzyme}_coordinates.txt"

}

runFlash(){
    
    echo
    echo "Running flash with parameters :"
    echo " m (minimum overlap) : ${flashOverlap}"
    echo " x (sum-of-mismatches/overlap-lenght) : ${flashErrorTolerance}"
    echo " p phred score min (33 or 64) : ${intQuals}"
    echo
    printThis="flash --interleaved-output -p ${intQuals} -m ${flashOverlap} -x ${flashErrorTolerance} READ1.fastq READ2.fastq > flashing.log"
    printToLogFile
    
    # flash --interleaved-output -p "${intQuals}" READ1.fastq READ2.fastq > flashing.log
    flash --interleaved-output -p "${intQuals}" -m "${flashOverlap}" -x "${flashErrorTolerance}" READ1.fastq READ2.fastq > flashing.log
    
    ls | grep out*fastq
    
    # This outputs these files :
    # flashing.log  FLASHED.fastq  out.hist  out.histogram  NONFLASHED.fastq  o

    echo "Read counts after flash :"
    
    flashedCount=0
    nonflashedCount=0
    
    if [ -s "out.extendedFrags.fastq" ] ; then
        flashedCount=$(( $( grep -c "" out.extendedFrags.fastq )/4 ))
    fi
    if [ -s "out.notCombined.fastq" ] ; then
        nonflashedCount=$(( $( grep -c "" out.notCombined.fastq )/4 ))
    fi
    
    mv -f out.extendedFrags.fastq FLASHED.fastq
    mv -f out.notCombined.fastq NONFLASHED.fastq
    
    echo "FLASHED.fastq (count of read pairs combined in flash) : ${flashedCount}"
    echo "NONFLASHED.fastq (not extendable via flash) : ${nonflashedCount}"
    
}

runCCanalyser(){
    
################################################################
# Running CAPTURE-C analyser for the aligned file..

#sampleForCCanalyser="RAW_${Sample}"
#samForCCanalyser="Combined_reads_REdig.sam"
#runDir=$( pwd )
#samDirForCCanalyser=${runDir}
#publicPathForCCanalyser="${PublicPath}/RAW"
#JamesUrlForCCanalyser="${JamesUrl}/RAW"


printThis="Running CAPTURE-C analyser for the aligned file.."
printToLogFile

testedFile="${CapturesiteFile}"
doTempFileTesting

mkdir -p "${publicPathForCCanalyser}"

printThis="perl ${RunScriptsPath}/${CCscriptname} -f ${samDirForCCanalyser}/${samForCCanalyser} -o ${CapturesiteFile} -r ${fullPathDpnGenome} --pf ${publicPathForCCanalyser} --pu ${JamesUrlForCCanalyser} -s ${sampleForCCanalyser} --genome ${GENOME} --ucscsizes ${ucscBuild} -w ${WINDOW} -i ${INCREMENT} --flashed ${FLASHED} --duplfilter ${DUPLFILTER} --maxfrags ${FRAGSPERREAD} ${otherParameters} "
printToLogFile

echo "-f Input filename "
echo "-r Restriction coordinates filename "
echo "-o Capturesitenucleotide position filename "
echo "--CCversion CS3 or CS4 or CS5 (which version of the duplicate filtering we will perform)"
echo "--pf Your public folder"
echo "--pu Your public url"
echo "-s Sample name (and the name of the folder it goes into)"
echo "-w Window size (default = 2kb)"
echo "-i Window increment (default = 200bp)"
echo "--dump Print file of unaligned reads (sam format)"
echo "--snp Force all capture points to contain a particular SNP"
echo "--limit Limit the analysis to the first n reads of the file"
echo "--genome Specify the genome (mm9 / hg18)"
echo "--ucscsizes Chromosome sizes file path"
echo "--globin Combines the two captures from the gene duplicates (HbA1 and HbA2)"
echo "--maxfrags : report at most this many fragments per read"
echo "--flashed	1 or 0 (are the reads in input sam combined via flash or not ? - run out.extended with 1 and out.not_combined with 0)"
echo "--duplfilter 1 or 0 (will the reads be duplicate filtered)\n"
echo "--parp Filter artificial chromosome chrPARP out before visualisation"
echo "--stringent enforces additional stringency - forces all reported subfragments to be unique"
echo "--stranded To replicate the strand-specific (i.e. wrong) duplicate filter of CB3a/CC3 and CB4a/CC4"
# echo "--umi Run contains UMI indices - alter the duplicate filter accordingly : ask Damien Downes how to prepare your files for pipeline, if you are interested in doing this"
echo "--wobble Wobble bin width. default 1(turned off). UMI runs recommendation 20, i.e. +/- 10bases wobble. to turn this off, set it to 1 base."


runDir=$( pwd )

# Copy used capture-site (REfragment) file for archiving purposes..
cp ${CapturesiteFile} usedCapturesiteFile.txt

# remove parameter file from possible earlier run..
rm -f parameters_for_normalisation.log

perl ${RunScriptsPath}/${CCscriptname} --CCversion "${CCversion}" -f "${samDirForCCanalyser}/${samForCCanalyser}" -o "${CapturesiteFile}" -r "${fullPathDpnGenome}" --pf "${publicPathForCCanalyser}" --pu "${JamesUrlForCCanalyser}" -s "${sampleForCCanalyser}" --genome "${GENOME}" --ucscsizes "${ucscBuild}" -w "${WINDOW}" -i "${INCREMENT}" --flashed "${FLASHED}" --duplfilter "${DUPLFILTER}" ${otherParameters}

echo "Contents of run folder :"
ls -lht

echo
echo "Contents of CCanalyser output folder ( ${sampleForCCanalyser}_${CCversion} ) "
ls -lht ${sampleForCCanalyser}_${CCversion}

echo
echo "Counts of output files - by file type :"

count=$( ls -1 ${publicPathForCCanalyser} | grep -c '.bw' )
echo
echo "${count} bigwig files (should be x2 the amount of capture-site (REfragment)s, if all had captures)"

count=$( ls -1 ${sampleForCCanalyser}_${CCversion} | grep -c '.wig' )
echo
echo "${count} wig files (should be x2 the amount of capture-site (REfragment)s, if all had captures)"

count=$( ls -1 ${sampleForCCanalyser}_${CCversion} | grep -c '.gff')
echo
echo "${count} gff files (should be x1 the amount of capture-site (REfragment)s, if all had captures)"

echo
echo "Output log files :"
ls -1 ${sampleForCCanalyser}_${CCversion} | grep '.txt'

echo
echo "Bed files :"
ls -1 ${sampleForCCanalyser}_${CCversion} | grep '.bed'

echo
echo "Sam files :"
ls -1 ${sampleForCCanalyser}_${CCversion} | grep '.sam'

echo
echo "Fastq files :"
ls -1 ${sampleForCCanalyser}_${CCversion} | grep '.fastq'   
    
}

runCCanalyserOnlyBlat(){
    
################################################################
# Running CAPTURE-C analyser for the aligned file..

printThis="Running onlyBlat CAPTURE-C analyser .."
printToLogFile

testedFile="${CapturesiteFile}"
doTempFileTesting

printThis="perl ${RunScriptsPath}/${CCscriptname} --onlyparamsforfiltering --CCversion ${CCversion} -o ${CapturesiteFile} --genome ${GENOME} -r "${fullPathDpnGenome}" --ucscsizes ${ucscBuild} ${otherParameters}"
printToLogFile

runDir=$( pwd )

# Copy used capture-site (REfragment) file for archiving purposes..
cp ${CapturesiteFile} usedCapturesiteFile.txt

# remove parameter file from possible earlier run..
rm -f parameters_for_filtering.log

perl ${RunScriptsPath}/${CCscriptname} --onlyparamsforfiltering --CCversion "${CCversion}" -o "${CapturesiteFile}" --genome "${GENOME}" -r "${fullPathDpnGenome}" --ucscsizes "${ucscBuild}" ${otherParameters}

if [ "$?" -ne 0 ];then
    printThis="Filtering parameter generation run reported error !"
    printNewChapterToLogFile
    paramGenerationRunFineOK=0
else
    printThis="Filtering parameter generation succeeded !"
    printToLogFile    
fi


echo "Contents of run folder :"
ls -lht

cat parameters_for_filtering.log
   
}


runTric(){

################################################################
# Running tri-c script

if [ "${BINNED_TRIC}" -eq "1" ]; then

printThis="Running 3-way contact matrix table generation (tri-C) in RE-fragment based AND binned modes .."
printToLogFile

printThis="perl ${RunScriptsPath}/TriC_MO_bin.pl -b ${TRIC_BIN} -sam F6_greenGraphs_combined_${Sample}_${CCversion}/COMBINED_reported_capture_reads_${CCversion}.sam -o ${CapturesiteFile} -r ${fullPathDpnGenome} --genome ${GENOME} --ucscsizes ${ucscBuild} --server ${SERVERTYPE}://${SERVERADDRESS} -pf ${diskFolder} -sf ${serverFolder} --name ${Sample} ${otherTricParameters}"
printToLogFile

perl ${RunScriptsPath}/TriC_MO_bin.pl -b ${TRIC_BIN} -sam "F6_greenGraphs_combined_${Sample}_${CCversion}/COMBINED_reported_capture_reads_${CCversion}.sam" -o ${CapturesiteFile} -r ${fullPathDpnGenome} --genome ${GENOME} --ucscsizes ${ucscBuild} --server "${SERVERTYPE}"'://'"${SERVERADDRESS}" -pf ${diskFolder} -sf ${serverFolder} --name ${Sample} ${otherTricParameters}

else

printThis="Running 3-way contact matrix table generation (tri-C) in RE-fragment based mode .."
printToLogFile

printThis="perl ${RunScriptsPath}/TriC_MO.pl -sam F6_greenGraphs_combined_${Sample}_${CCversion}/COMBINED_reported_capture_reads_${CCversion}.sam -o ${CapturesiteFile} -r ${fullPathDpnGenome} --genome ${GENOME} --ucscsizes ${ucscBuild} --server ${SERVERTYPE}://${SERVERADDRESS} -pf ${diskFolder} -sf ${serverFolder} --name ${Sample} ${otherTricParameters}"
printToLogFile

perl ${RunScriptsPath}/TriC_MO.pl -sam "F6_greenGraphs_combined_${Sample}_${CCversion}/COMBINED_reported_capture_reads_${CCversion}.sam" -o ${CapturesiteFile} -r ${fullPathDpnGenome} --genome ${GENOME} --ucscsizes ${ucscBuild} --server "${SERVERTYPE}"'://'"${SERVERADDRESS}" -pf ${diskFolder} -sf ${serverFolder} --name ${Sample} ${otherTricParameters}

    
fi

printThis="Running 3-way contact matrix visualisation (tri-C) .."
printToLogFile

TEMP_subfolder="F6_greenGraphs_combined_${Sample}_${CCversion}/${Sample}_TriC"
TEMP_tricParams1="-b ${TRIC_BIN} -o ${TEMP_subfolder} -t ${TRIC_MAX} ${TRIC_OTHER_VISUAL_PARAMS}"

rm -f USED_TriC_plottingParameters.txt

if [ ! -s "PIPE_TriC_plottingParameters.txt" ]; then

printThis="Generating plotting parameters on the fly (no PIPE_TriC_plottingParameters.txt file in run directory)"
printToLogFile
printThis="Finetuning Tri-C plotting can be done by re-naming the generated USED_TriC_plottingParameters.txt to PIPE_TriC_plottingParameters.txt \n and editing the plotting window (columns 2-3) for each capture site,\n and restarting the run with --onlyTriC to only repeat TriC analysis."
printToLogFile

else
    
printThis="Reading plotting parameters from PIPE_TriC_plottingParameters.txt file (given in the run directory)"
printToLogFile    

fi

for file in ${TEMP_subfolder}/${Sample}_*_TriC_interactions.txt
do

TEMP_tricParams2=""
TEMPcaptureSite=""


if [ ! -s "PIPE_TriC_plottingParameters.txt" ]; then

    TEMPcaptureSite=$(basename $file | sed 's/^'${Sample}'_//' | sed 's/_TriC_interactions.txt//')
    TEMPchr=$(cat ${CapturesiteFile} | grep '^'${TEMPcaptureSite}'\s' | cut -f 2)
    
    TEMPmiddle=$(cat ${CapturesiteFile} | grep '^'${TEMPcaptureSite}'\s' | cut -f 3,4 | awk '{ print int(($1+$2)/20000)*10000}')
    
    TEMPstr=$(echo ${TEMPmiddle} | awk '{if($1<100000)print 1; else print $1-100000}')
    TEMPstp=$(echo ${TEMPmiddle} | awk '{ print $1+100000}')
    TEMPchrSize=$(cat ${ucscBuild} | grep '^[Cc]hr'${TEMPchr}'\s' | sed 's/.*\s//')
    
    TEMPstpClip=$(echo -e "${TEMPstp}\t${TEMPchrSize}" | awk '{if($1>$2)print $2; else print $1}')
    
    echo -e "${TEMPcaptureSite}\t${TEMPchr}\t${TEMPstr}\t${TEMPstpClip}" >> USED_TriC_plottingParameters.txt
    
    TEMP_tricParams2="-c chr${TEMPchr} --str ${TEMPstr} --stp ${TEMPstpClip}"


else
    
    TEMPcaptureSite=$(basename $file | sed 's/^'${Sample}'_//' | sed 's/_TriC_interactions.txt//')
    TEMPchr=$(cat PIPE_TriC_plottingParameters.txt | grep '^'${TEMPcaptureSite}'\s' | cut -f 2)
    TEMPstr=$(cat PIPE_TriC_plottingParameters.txt | grep '^'${TEMPcaptureSite}'\s' | cut -f 3)
    TEMPstp=$(cat PIPE_TriC_plottingParameters.txt | grep '^'${TEMPcaptureSite}'\s' | cut -f 4)
    
    echo -e "${TEMPcaptureSite}\t${TEMPchr}\t${TEMPstr}\t${TEMPstp}" >> USED_TriC_plottingParameters.txt
    
    TEMP_tricParams2="-c chr${TEMPchr} --str ${TEMPstr} --stp ${TEMPstp}"
    
fi

printThis="Visualising capture site ${TEMPcaptureSite}"
printToLogFile
printThis="python ${RunScriptsPath}/TriC_matrix_MO.py -f ${file} ${TEMP_tricParams1} ${TEMP_tricParams2} > ${TEMP_subfolder}/${Sample}_${TEMPcaptureSite}_TriC_visualisationRun.log"
printToLogFile

python ${RunScriptsPath}/TriC_matrix_MO.py -f ${file} ${TEMP_tricParams1} ${TEMP_tricParams2} > ${TEMP_subfolder}/${Sample}_${TEMPcaptureSite}_TriC_visualisationRun.log

copyTricPdfToPublic

done

ls -lht F6_greenGraphs_combined_${Sample}_${CCversion}/${Sample}_TriC

# public html update for pdfs ..

updateTricHub

}
