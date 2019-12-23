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


doTrackExist(){
    # NEEDS THESE TO BE SET BEFORE CALL :
    #trackName=""
     
    if [ -s "${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt" ]; then
    
    echo -e "grep bigDataUrl ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt | grep -c \"${fileName}\$\" " > temp.command
    chmod u=rwx temp.command
    trackExists=$(( $(./temp.command) ))
    rm -f temp.command
    
    else
    trackExists=0
    
    fi
}

doMultiWigParent(){
    
    # NEEDS THESE TO BE SET BEFORE CALL :
    #longLabel=""
    #trackName=""
    #overlayType=""
    #windowingFunction=""
    #visibility=""
    
    echo "" >> TEMP2_tracks.txt
    echo "#--------------------------------------" >> TEMP2_tracks.txt
    echo "" >> TEMP2_tracks.txt
    
    echo "track ${trackName}" >> TEMP2_tracks.txt
    echo "container multiWig" >> TEMP2_tracks.txt
    echo "shortLabel ${trackName}" >> TEMP2_tracks.txt
    echo "longLabel ${longLabel}" >> TEMP2_tracks.txt
    echo "type bigWig" >> TEMP2_tracks.txt
    echo "visibility ${visibility}" >> TEMP2_tracks.txt
    echo "aggregate ${overlayType}" >> TEMP2_tracks.txt
    echo "showSubtrackColorOnUi on" >> TEMP2_tracks.txt
    #echo "windowingFunction maximum" >> TEMP2_tracks.txt
    #echo "windowingFunction mean" >> TEMP2_tracks.txt
    echo "windowingFunction ${windowingFunction}" >> TEMP2_tracks.txt
    echo "configurable on" >> TEMP2_tracks.txt
    echo "dragAndDrop subtracks" >> TEMP2_tracks.txt
    echo "autoScale on" >> TEMP2_tracks.txt
    echo "alwaysZero on" >> TEMP2_tracks.txt
    echo "" >> TEMP2_tracks.txt
    
}

doMultiWigChild(){
    
    # NEEDS THESE TO BE SET BEFORE CALL
    # parentTrack=""
    # trackName=""
    # fileName=".bw"
    # trackColor=""
    # trackPriority=""
    # bigWigSubfolder="${PublicPath}/FILTERED"
    
    # Does this track have data file which has non-zero size?
    if [ -s "${PublicPath}/${bigWigSubfolder}/${fileName}" ]; then
    
    echo "track ${trackName}" >> TEMP2_tracks.txt
    echo "parent ${parentTrack}" >> TEMP2_tracks.txt
    echo "bigDataUrl ${ServerAndPath}/${bigWigSubfolder}/${fileName}" >> TEMP2_tracks.txt
    # These are super long paths. using relative paths instead !
    #echo "bigDataUrl ${fileName}" >> TEMP2_tracks.txt
    echo "shortLabel ${trackName}" >> TEMP2_tracks.txt
    echo "longLabel ${trackName}" >> TEMP2_tracks.txt
    echo "type bigWig" >> TEMP2_tracks.txt
    echo "color ${trackColor}" >> TEMP2_tracks.txt
    echo "html http://${JamesUrl}/${Sample}_${CCversion}_description" >> TEMP2_tracks.txt
    echo "priority ${trackPriority}" >> TEMP2_tracks.txt
    echo "" >> TEMP2_tracks.txt
    
    else
    
    echo "Cannot find track ${PublicPath}/${bigWigSubfolder}/${fileName} - not writing it into ${parentTrack} track "  >> "/dev/stderr"
    
    fi
}

doRegularTrack(){
    
    # NEEDS THESE TO BE SET BEFORE CALL
    # trackName=""
    # longLabel=""
    # fileName=".bw"
    # trackColor=""
    # trackPriority=""
    # visibility=""
    # trackType="bb" "bw"
   
    # Is this track already written to the tracks.txt file?
    doTrackExist
    if [ "${trackExists}" -eq 0 ] ; then
   
    # Does this track have data file which has non-zero size?
    if [ -s "${publicPathForCCanalyser}/${fileName}" ] ; then

    echo "" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "#--------------------------------------" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    
    echo "track ${trackName}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    # These are super long paths. using relative paths instead !
    #echo "bigDataUrl ${ServerAndPath}/${bigWigSubfolder}/${fileName}" | sed 's/\/\//\//g' >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "bigDataUrl ${fileName}"  >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "shortLabel ${trackName}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "longLabel ${longLabel}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    
    if [ "${trackType}" = "bb" ] ; then
    echo "type bigBed" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    # If we want item name as well.."
    elif [ "${trackType}" = "bb4" ] ; then
    echo "type bigBed 4" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    # Defaults to "bw"
    else
    echo "type bigWig" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    #echo "color ${trackColor}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    fi
    
    echo "color ${trackColor}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "visibility ${visibility}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "priority ${trackPriority}" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "windowingFunction maximum" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "autoScale on" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "alwaysZero on" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    echo "" >> ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
    
    else
        echo "TRACK DESCRIPTION NOT CREATED - track ${trackName} does not have size in ${publicPathForCCanalyser}/${fileName}" >> "/dev/stderr"
    fi
    else
        echo -n ""
        # echo "TRACK DESCRIPTION NOT CREATED - track ${trackName} already exists in ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt" >> "/dev/stderr"
    fi
    
}

copyPreCCanalyserLogFilesToPublic(){
    
# Making a public folder for log files
printThis="Making a public folder for log files"
printToLogFile
mkdir -p "${PublicPath}/${Sample}_logFiles"

# Copying log files
printThis="Copying log files to public folder"
printToLogFile
cp -rf READ1_fastqc_ORIGINAL "${PublicPath}/${Sample}_logFiles"
cp -rf READ2_fastqc_ORIGINAL "${PublicPath}/${Sample}_logFiles"
cp -rf READ1_fastqc_TRIMMED "${PublicPath}/${Sample}_logFiles"
cp -rf READ2_fastqc_TRIMMED "${PublicPath}/${Sample}_logFiles"
cp -rf FLASHED_fastqc "${PublicPath}/${Sample}_logFiles"
cp -rf NONFLASHED_fastqc "${PublicPath}/${Sample}_logFiles"
cp -rf FLASHED_REdig_fastqc "${PublicPath}/${Sample}_logFiles"
cp -rf NONFLASHED_REdig_fastqc "${PublicPath}/${Sample}_logFiles"

cp -f ./read_trimming.log "${PublicPath}/${Sample}_logFiles/${Sample}_read_trimming.log"
cp -f ./flashing.log "${PublicPath}/${Sample}_logFiles/${Sample}_flashing.log"
cp -f ./out.hist "${PublicPath}/${Sample}_logFiles/${Sample}_flash.hist"
cp -f ./NONFLASHED_${REenzyme}digestion.log "${PublicPath}/${Sample}_logFiles/${Sample}_NONFLASHED_${REenzyme}digestion.log"
cp -f ./FLASHED_${REenzyme}digestion.log "${PublicPath}/${Sample}_logFiles/${Sample}_FLASHED_${REenzyme}digestion.log"

}

copyTricPdfToPublic(){
printThis="Copying pdf to public folder"
printToLogFile
rm -f ${diskFolder}/${Sample}_TriC/${Sample}_${TEMPcaptureSite}_TriC_interactions_${TRIC_BIN}_${TRIC_MAX}.pdf
cp ${TEMP_subfolder}/${Sample}_${TEMPcaptureSite}_TriC_interactions_${TRIC_BIN}_${TRIC_MAX}.pdf ${diskFolder}/${Sample}_TriC/.

}

updateTricHub(){
    
    printThis="Updating TriC html file"
    printToLogFile
    
    echo $(date) > temp_description.html
    echo "<h4>TriC visualisation files here : </h4>" >> temp_description.html
    
    for publicPdf in ${diskFolder}/${Sample}_TriC/${Sample}_*_TriC_interactions_${TRIC_BIN}_${TRIC_MAX}.pdf
    do
        TEMPbasename=$(basename ${publicPdf})
        TEMPcapture=$(echo ${TEMPbasename}| sed 's/^'${Sample}'_//' | sed 's/_TriC_interactions_'${TRIC_BIN}_${TRIC_MAX}'.pdf//')
        echo "<li><a target="_blank" href=\"${ServerAndPath}/${Sample}_TriC/${TEMPbasename}\" >${TEMPbasename}</a>  </li>" >> temp_description.html
    done

    echo "<hr />" >> temp_description.html
    
    echo "<h4>TriC run counters (before visualisation)  : </h4>" >> temp_description.html
    
    cat ${diskFolder}/${Sample}_TriC/description.html >> temp_description.html    
    mv -f temp_description.html ${diskFolder}/${Sample}_TriC/description.html
    
}

updateHub_part2c(){
    
    # Link the file to each of the existing tracks..
    seddedUrl=$( echo ${JamesUrl} | sed 's/\//\\\//g' )
    echo "sed -i 's/alwaysZero on/alwaysZero on\nhtml http\:\/\/${seddedUrl}\/${Sample}_description/' ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt " > temp.command

    chmod u=rwx temp.command
    # cat temp.command
    ./temp.command
    rm -f temp.command
    
}

updateHub_part3(){
    TEMPname=$( echo ${sampleForCCanalyser} | sed 's/_.*//' )
    echo
    if [ "${TEMPname}" == "RAW" ] || [ "${TEMPname}" == "PREfiltered" ] || [ "${TEMPname}" == "FILTERED" ] || [ "${TEMPname}" == "COMBINED" ] ; then
    echo "Generated a data hub in : ${ServerAndPath}/${TEMPname}/${sampleForCCanalyser}_${CCversion}_hub.txt"
    else
    echo "Generated a data hub in : ${ServerAndPath}/${sampleForCCanalyser}_${CCversion}_hub.txt"
    fi
    echo 'How to load this hub to UCSC : http://userweb.molbiol.ox.ac.uk/public/telenius/CaptureCompendium/CCseqBasic/DOCS/HUBtutorial_AllGroups_160813.pdf'    

}

updateHub_part3final(){
    echo
    echo "Generated a data hub for RAW data in : ${ServerAndPath}/RAW/RAW_${Sample}_${CCversion}_hub.txt"
    echo "( pre-filtered data for DEBUGGING purposes is here : ${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_${CCversion}_hub.txt )"
    echo "Generated a data hub for FILTERED data in : ${ServerAndPath}/FILTERED/FILTERED_${Sample}_${CCversion}_hub.txt"
    echo "Generated a data hub for FILTERED flashed+nonflashed combined data in : ${ServerAndPath}/COMBINED/COMBINED_${Sample}_${CCversion}_hub.txt"
    echo
    echo "Generated a COMBINED data hub (of all the above) in : ${ServerAndPath}/${Sample}_${CCversion}_hub.txt"

}

updateHub_part3finalTric(){
    echo
    echo "Generated a tri-C data hub in : ${ServerAndPath}/${Sample}_TriC/hub.txt"

}

updateHub_part3finalTricBinned(){
    echo "Generated a     BINNED tri-C data hub in : ${ServerAndPath}/${Sample}_TriC/hub_bin.txt"
    
}

updateHub_part3finalTricDescription(){
    echo "And Tri-C visualisation html page in : ${ServerAndPath}/${Sample}_TriC/description.html"

}

updateHub_part3final2(){
    echo 'How to load this hub to UCSC : http://userweb.molbiol.ox.ac.uk/public/telenius/CaptureCompendium/CCseqBasic/DOCS/HUBtutorial_AllGroups_160813.pdf'    
}

updateCCanalyserDataHub(){
    
printThis="Updating the public folder with analysis log files.."
printToLogFile

temptime=$( date +%d%m%y )

mkdir -p ${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles

#samForCCanalyser="F1_${Sample}_pre${CCversion}/Combined_reads_REdig.sam"
#samBasename=$( echo ${samForCCanalyser} | sed 's/.*\///' | sed 's/\_FLASHED.sam$//' | sed 's/\_NONFLASHED.sam$//' )

if [ -s ${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt ]
then
cp -f "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s ${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt ]
then
cp -f "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s ${sampleForCCanalyser}_${CCversion}/COMBINED_report_${CCversion}.txt ]
then
cp -f "${sampleForCCanalyser}_${CCversion}/COMBINED_report_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi


if [ -s "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report2_${CCversion}.txt" ]
then
    cp -f "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report2_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s  "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report2_${CCversion}.txt" ]
then
    cp -f "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report2_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s ${sampleForCCanalyser}_${CCversion}/COMBINED_report2_${CCversion}.txt ]
then
cp -f "${sampleForCCanalyser}_${CCversion}/COMBINED_report2_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi

if [ -s "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report3_${CCversion}.txt" ]
then
    cp -f "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report3_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s  "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report3_${CCversion}.txt" ]
then
    cp -f "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report3_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s ${sampleForCCanalyser}_${CCversion}/COMBINED_report3_${CCversion}.txt ]
then
cp -f "${sampleForCCanalyser}_${CCversion}/COMBINED_report3_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi


if [ -s "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report4_${CCversion}.txt" ]
then
    cp -f "${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report4_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
if [ -s  "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report4_${CCversion}.txt" ]
then
    cp -f "${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report4_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
 if [ -s ${sampleForCCanalyser}_${CCversion}/COMBINED_report4_${CCversion}.txt ]
then
cp -f "${sampleForCCanalyser}_${CCversion}/COMBINED_report4_${CCversion}.txt" "${publicPathForCCanalyser}/${sampleForCCanalyser}_logFiles/."
fi
    
# Make the bigbed file from the bed file of capture-site (REfragment) coordinates and used exlusions ..

tail -n +2 "${CapturesiteFile}" | sort -T $(pwd) -k1,1 -k2,2n > tempBed.bed
bedOrigName=$( echo "${CapturesiteFile}" | sed 's/\..*//' | sed 's/.*\///' )
bedname=$( echo "${CapturesiteFile}" | sed 's/\..*//' | sed 's/.*\///' | sed 's/^/'${Sample}'_/' )

# Capturesite coordinates 
tail -n +2 "${sampleForCCanalyser}_${CCversion}/${bedOrigName}.bed" | awk 'NR%2==1' | sort -T $(pwd) -k1,1 -k2,2n > tempBed.bed
bedToBigBed -type=bed9 tempBed.bed ${ucscBuild} "${sampleForCCanalyser}_${CCversion}/${bedname}_capturesite.bb"
rm -f tempBed.bed

# Exclusion fragments
tail -n +2 "${sampleForCCanalyser}_${CCversion}/${bedOrigName}.bed" | awk 'NR%2==0' | sort -T $(pwd) -k1,1 -k2,2n > tempBed.bed
bedToBigBed -type=bed9 tempBed.bed ${ucscBuild} "${sampleForCCanalyser}_${CCversion}/${bedname}_exclusion.bb"
rm -f tempBed.bed

mv -f "${sampleForCCanalyser}_${CCversion}/${bedname}_capturesite.bb" ${publicPathForCCanalyser}
mv -f "${sampleForCCanalyser}_${CCversion}/${bedname}_exclusion.bb" ${publicPathForCCanalyser}

    fileName=$( echo ${publicPathForCCanalyser}/${bedname}_capturesite.bb | sed 's/^.*\///' )
    trackName=$( echo ${fileName} | sed 's/\.bb$//' )
    longLabel="${trackName}_coordinates"
    trackColor="133,0,122"
    trackPriority="1"
    visibility="full"
    trackType="bb"
    
    doRegularTrack
    
    fileName=$( echo ${publicPathForCCanalyser}/${bedname}_exclusion.bb | sed 's/^.*\///' )
    trackName=$( echo ${fileName} | sed 's/\.bb$//' )
    longLabel="${trackName}_coordinates"
    trackColor="133,0,0"
    trackPriority="2"
    visibility="full"
    trackType="bb"
    
    doRegularTrack

    
# Add the missing tracks - if the hub was not generated properly in the perl..
    
for file in ${publicPathForCCanalyser}/*.bw
do
    fileName=$( echo ${file} | sed 's/^.*\///' )
    trackName=$( echo ${fileName} | sed 's/\.bw$//' )
    longLabel=${trackName}
    trackColor="0,0,0"
    trackPriority="200"
    visibility="hide"
    trackType="bw"
    bigWigSubfolder=${bigWigSubfolder}
    
    doRegularTrack
    
done

    updateHub_part2c

if [ -s ${runDir}/${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt ]
then

echo "##############################################"
echo "Report of FLASHED reads CCanalyser results :"
echo
cat "${runDir}/${sampleForCCanalyser}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt"
echo

fi


if [ -s ${runDir}/${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt ]
then
    
echo "##############################################"
echo "Report of NONFLASHED reads CCanalyser results :"
echo
cat "${runDir}/${sampleForCCanalyser}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt"
echo

fi

if [ -s ${runDir}/${sampleForCCanalyser}_${CCversion}/COMBINED_report_${CCversion}.txt ]
then
    
echo "##############################################"
echo "Report of COMBINED reads CCanalyser results :"
echo
cat "${runDir}/${sampleForCCanalyser}_${CCversion}/COMBINED_report_${CCversion}.txt"
echo


fi

    updateHub_part3 
    
}

alterCustomGenomeNames(){
    
printThis="Changing the genome names (custom genome ${GENOME} was used - changing the visualising genome to be ${ucscBuildName} ) .."
printToLogFile

# genomes.txt name : common part
genomesTxt="${Sample}_${CCversion}_genomes.txt"

sed -i 's/^genome\s.*/genome'${ucscBuildName}'/' ${PublicPath}/${genomesTxt}
sed -i 's/^genome\s.*/genome'${ucscBuildName}'/' ${PublicPath}/RAW/RAW_${genomesTxt}
sed -i 's/^genome\s.*/genome'${ucscBuildName}'/' ${PublicPath}/PREfiltered/PREfiltered_${genomesTxt}
sed -i 's/^genome\s.*/genome'${ucscBuildName}'/' ${PublicPath}/FILTERED/FILTERED_${genomesTxt}


}

generateCombinedDataHub(){
    
printThis="Generating summary figure for the data.."
printToLogFile


# Setting track colors : normal or redgreen blindness ?

# PINK GREEN (default)
redcolor="255,74,179"
orangecolor="255,140,0"
greencolor="62,176,145"

# RED GREEN (old default - very redgreen blindness-unfriendly colors)

if [ "${REDGREEN}" -eq 1 ];then
    
redcolor="255,0,0"
orangecolor="255,140,0"
greencolor="0,200,0"

fi


runDir=$( pwd )
mkdir F7_summaryFigure_${Sample}_${CCversion}
cd F7_summaryFigure_${Sample}_${CCversion}
echo "summaryFigure is folder where the statistics summary figure is generated" > a_folder_containing_statistics_summary_figure

echo "${CapturePlotPath}/countsFromCCanalyserOutput.sh ${runDir} ${Sample} ${CCversion} ${REenzyme}" >> "/dev/stderr"
${CapturePlotPath}/countsFromCCanalyserOutput.sh ${runDir} ${Sample} ${CCversion} ${REenzyme} > counts.py

echo "python ${CapturePlotPath}/generatePercentages.py" >> "/dev/stderr"
python ${CapturePlotPath}/generatePercentages.py > percentages.txt

echo "python ${CapturePlotPath}/drawFigure.py" >> "/dev/stderr"
python ${CapturePlotPath}/drawFigure.py

figureDimensions=$( file summary.png | sed 's/\s//g' | tr ',' '\t' | cut -f 2 | sed 's/^/width=\"/' | sed 's/x/\" height=\"/' | sed 's/$/\"/' )
cp summary.png ${publicPathForCCanalyser}/.
cp counts.py ${publicPathForCCanalyser}/counts.txt

cd ..
    
printThis="Updating the public folder with analysis log files.."
printToLogFile

temptime=$( date +%d%m%y )

# Here add :

# Generate the hub itself, as well as genomes.txt

#${publicPathForCCanalyser}/${Sample}_${CCversion}_hub.txt

# cat MES0_CC2_hub.txt
#hub MES0_CC2
#shortLabel MES0_CC2
#longLabel MES0_CC2_CaptureC
#genomesFile http://userweb.molbiol.ox.ac.uk/public/mgosden/Dilut_Cap/MES0_CC2_genomes.txt
#email james.davies@trinity.ox.ac.uk

echo "hub ${Sample}_${CCversion}" > ${PublicPath}/${Sample}_${CCversion}_hub.txt
echo "shortLabel ${Sample}_${CCversion}" >> ${PublicPath}/${Sample}_${CCversion}_hub.txt
echo "longLabel ${Sample}_${CCversion}_CaptureC" >> ${PublicPath}/${Sample}_${CCversion}_hub.txt
echo "genomesFile ${ServerAndPath}/${Sample}_${CCversion}_genomes.txt" >> ${PublicPath}/${Sample}_${CCversion}_hub.txt
echo "email jelena.telenius@gmail.com" >> ${PublicPath}/${Sample}_${CCversion}_hub.txt

#${publicPathForCCanalyser}/${Sample}_${CCversion}_genomes.txt

# cat MES0_CC2_genomes.txt 
#genome mm9
#trackDb http://sara.molbiol.ox.ac.uk/public/mgosden/Dilut_Cap//MES0_CC2_tracks.txt

#echo "genome ${GENOME}" > ${ServerAndPath}/${Sample}_${CCversion}_genomes.txt
#echo "trackDb ${ServerAndPath}/${Sample}_${CCversion}_tracks.txt" >> ${ServerAndPath}/${Sample}_${CCversion}_genomes.txt

echo "genome ${GENOME}" > TEMP_genomes.txt
echo "trackDb ${ServerAndPath}/${Sample}_${CCversion}_tracks.txt" >> TEMP_genomes.txt
mv -f TEMP_genomes.txt ${PublicPath}/${Sample}_${CCversion}_genomes.txt

# Catenate the tracks.txt files to form new tracks.txt
tracksTxt="${Sample}_${CCversion}_tracks.txt"

# [telenius@deva run15]$ cat /public/telenius/capturetests/test_040316_O/test_040316_O/CC4/RAW/RAW_test_040316_O_CC4_tracks.txt | grep track
# [telenius@deva run15]$ cat /public/telenius/capturetests/test_040316_O/test_040316_O/CC4/PREfiltered/PREfiltered_test_040316_O_CC4_tracks.txt | grep track
# [telenius@deva run15]$ cat /public/telenius/capturetests/test_040316_O/test_040316_O/CC4/FILTERED/FILTERED_test_040316_O_CC4_tracks.txt | grep track


cat ${PublicPath}/RAW/RAW_${tracksTxt} ${PublicPath}/PREfiltered/PREfiltered_${tracksTxt} ${PublicPath}/FILTERED/FILTERED_${tracksTxt} > TEMP_tracks.txt

# echo
# echo
# cat TEMP_tracks.txt
# echo
# echo

# Make proper redgreen tracks based on the RAW and FILTERED tracks..

#doMultiWigParent    
    # NEEDS THESE TO BE SET BEFORE CALL :
    #longLabel=""
    #trackName=""
    #overlayType=""
    #windowingFunction=""
    #visibility=""
    
#doMultiWigChild   
    # NEEDS THESE TO BE SET BEFORE CALL
    # parentTrack=""
    # trackName=""
    # fileName=".bw"
    # trackColor=""
    # trackPriority=""
    
  rm -f TEMP2_tracks.txt
    
 trackList=$( cat TEMP_tracks.txt | grep track | grep RAW | sed 's/^track RAW_//' )
 filenameList=$( cat TEMP_tracks.txt | grep bigDataUrl | grep RAW | sed 's/^bigDataUrl .*RAW\///' )
 
 cat TEMP_tracks.txt | grep track | grep RAW | sed 's/^track RAW_//' | sed 's/^track win_RAW_/win_/' > TEMP_trackList.txt
 cat TEMP_tracks.txt | grep bigDataUrl | grep RAW | sed 's/^bigDataUrl .*RAW\///' > TEMP_bigDataUrlList.txt
 
 list=$( paste TEMP_trackList.txt TEMP_bigDataUrlList.txt | sed 's/\s/,/' )
 echo
 echo RAW track list
 paste TEMP_trackList.txt TEMP_bigDataUrlList.txt
 echo
 
 rm -f TEMP_trackList.txt TEMP_bigDataUrlList.txt
 
  
 for track in $list
 do
    echo $track
    
    trackname=$( echo $track | sed 's/,.*//' )
    filename=$( echo $track | sed 's/.*,//')
    
    longLabel="CC_${trackname} all mapped RED, duplicate-filtered ORANGE, dupl+ploidy+blat-filtered GREEN"
    trackName="${trackname}"
    overlayType="solidOverlay"
    windowingFunction="maximum"
    visibility="hide"
    doMultiWigParent
    
    parentTrack="${trackname}"
    trackName="${trackname}_raw"
    fileName="${filename}"
    bigWigSubfolder="RAW"
    trackColor="${redcolor}"
    trackPriority="100"
    doMultiWigChild
    
 done

 
 cat TEMP_tracks.txt | grep track | grep PREfiltered  | sed 's/^track PREfiltered_//' | sed 's/^track win_PREfiltered_/win_/' > TEMP_trackList.txt
 cat TEMP_tracks.txt | grep bigDataUrl | grep PREfiltered| sed 's/^bigDataUrl .*PREfiltered\///' > TEMP_bigDataUrlList.txt
 
 list=$( paste TEMP_trackList.txt TEMP_bigDataUrlList.txt | sed 's/\s/,/' )
 
 echo
 echo PREfiltered track list
  paste TEMP_trackList.txt TEMP_bigDataUrlList.txt
 echo
 
 rm -f TEMP_trackList.txt TEMP_bigDataUrlList.txt

 
 for track in $list
 do
    echo $track
    
    trackname=$( echo $track | sed 's/,.*//')
    filename=$( echo $track | sed 's/.*,//')
    
    parentTrack="${trackname}"
    trackName="${trackname}_PREfiltered"
    fileName="${filename}"
    bigWigSubfolder="PREfiltered"
    trackColor="${orangecolor}"
    trackPriority="110"
    doMultiWigChild
    
 done
 
 cat TEMP_tracks.txt | grep track | grep FILTERED  | sed 's/^track FILTERED_//' | sed 's/^track win_FILTERED_/win_/' > TEMP_trackList.txt
 cat TEMP_tracks.txt | grep bigDataUrl | grep FILTERED | sed 's/^bigDataUrl .*FILTERED\///' > TEMP_bigDataUrlList.txt
 
 list=$( paste TEMP_trackList.txt TEMP_bigDataUrlList.txt | sed 's/\s/,/' )
 
 echo
 echo FILTERED track list
  paste TEMP_trackList.txt TEMP_bigDataUrlList.txt
 echo
 
 rm -f TEMP_trackList.txt TEMP_bigDataUrlList.txt

 
 for track in $list
 do
    echo $track
    
    trackname=$( echo $track | sed 's/,.*//')
    filename=$( echo $track | sed 's/.*,//')
    
    parentTrack="${trackname}"
    trackName="${trackname}_filtered"
    fileName="${filename}"
    bigWigSubfolder="FILTERED"
    trackColor="${greencolor}"
    trackPriority="120"
    doMultiWigChild
    
 done

 rm -f TEMP_tracks.txt
 
    # Adding the combined files and the capture-site (REfragment) tracks
    
    # Here used to be also sed 's/visibility hide/visibility full/' : to set only the COMBINED tracks visible.
    # As multi-capture samples grep more frequent, this was taken out of the commands below.
    cat ${PublicPath}/COMBINED/COMBINED_${tracksTxt} | sed 's/color 0,0,0/color '"${greencolor}"'/' \
    | sed 's/priority 200/windowingFunction maximum\npriority 10/' \
    | sed 's/bigDataUrl .*COMBINED\//bigDataUrl COMBINED\//' | grep -v "^html" \
    > TEMP3_tracks.txt
   
    cp ${PublicPath}/COMBINED/*.bb ${PublicPath}/.
   
    cat TEMP2_tracks.txt TEMP3_tracks.txt > TEMP4_tracks.txt
    rm -f TEMP2_tracks.txt TEMP3_tracks.txt 
    
    
    # Move over..
    
    mv -f TEMP4_tracks.txt ${PublicPath}/${tracksTxt}
    
    
    
    # Adding the bigbed track for BLAT-filter-marked RE-fragments :

    cat filteringLogFor_PREfiltered_${Sample}_${CCversion}/BlatPloidyFilterRun/BLAT_PLOIDY_FILTERED_OUTPUT/blatFilterMarkedREfragments.bed | sort -T $(pwd) -k1,1 -k2,2n > tempBed.bed
    bedToBigBed -type=bed4 tempBed.bed ${ucscBuild} ${sampleForCCanalyser}_${CCversion}_blatFilterMarkedREfragments.bb
    rm -f tempBed.bed
    
    mv -f ${sampleForCCanalyser}_${CCversion}_blatFilterMarkedREfragments.bb ${publicPathForCCanalyser}

    fileName=$( echo ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_blatFilterMarkedREfragments.bb | sed 's/^.*\///' )
    trackName=$( echo ${fileName} | sed 's/\.bb$//' )
    longLabel="${trackName}"
    trackColor="133,0,122"
    trackPriority="1"
    visibility="pack"
    trackType="bb4"
    bigWigSubfolder=""
    
    doRegularTrack
    
    printThis="Writing the description html-document"
    printToLogFile

# Write the beginning of the html file

    echo "<!DOCTYPE HTML PUBLIC -//W3C//DTD HTML 4.01//EN" > begin.html
    echo "http://www.w3.org/TR/html4/strict.dtd" >> begin.html
    echo ">" >> begin.html
    echo " <html lang=en>" >> begin.html
    echo " <head>" >> begin.html
    echo " <title> ${hubNameList[0]} data hub in ${genomeName} </title>" >> begin.html
    echo " </head>" >> begin.html
    echo " <body>" >> begin.html

    # Generating TimeStamp 
    TimeStamp=($( date | sed 's/[: ]/_/g' ))
    DateTime="$(date)"
    
    echo "<p>Data produced ${DateTime} with CapC pipeline (coded by James Davies, pipelined by Jelena Telenius, located in ${CapturePipePath} )</p>" > temp_description.html
    
    echo "<hr />" >> temp_description.html
    
    echo "Capturesite coordinates given to the run :" >> temp_description.html
    echo "<pre>" >> temp_description.html
    cat ${CapturesiteFile} | cut -f 1-4 | awk '{print $1"\tchr"$2"\t"$3"\t"$4}' >> temp_description.html
    echo "</pre>" >> temp_description.html
    
    echo "<hr />" >> temp_description.html
    
#    echo "<p>User manual - to understand the pipeline and the output :  <a target="_blank" href=\"http://sara.molbiol.ox.ac.uk/public/jdavies/MANUAL_for_pipe/PipeUserManual.pdf\" >CapturePipeUserManual.pdf</a></p>" >> temp_description.html
    
    echo "<hr />" >> temp_description.html
    
    echo "<p>Data located in : $(pwd)</p>" >> temp_description.html
    echo "<p>Sample name : ${Sample}, containing fastq files : ${Read1} and ${Read2}</p>" >> temp_description.html

# These have to be listed for all runs (RAW and PREfiltered and FILTERED)


# Summary figure and report ------------------------------------------------------------

    echo "<hr />" >> temp_description.html
    echo "<img" >> temp_description.html
    echo "src=${ServerAndPath}/summary.png" >> temp_description.html
    echo "${figureDimensions}" >> temp_description.html
    echo "alt=\"Summary of the analysis\"/>" >> temp_description.html
    echo "<hr />" >> temp_description.html
    echo "<li>Above figure as READ COUNTS in text file :<a target="_blank" href=\"${ServerAndPath}/counts.txt\" >readCounts.txt</a></li>"  >> temp_description.html
  
    echo "<hr />" >> temp_description.html
    
    echo "<li><b>Final counts summary files</b> :" >> temp_description.html
    echo " ( <a target="_blank" href=\"${ServerAndPath}/COMBINED/COMBINED_${Sample}_logFiles/COMBINED_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/COMBINED/COMBINED_${Sample}_logFiles/COMBINED_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html
    
    echo "<hr />" >> temp_description.html
    
# [telenius@deva CC4]$ ls RAW/RAW_test_040316_J_logFiles/
# FLASHED_REdig_report_CC4.txt  NONFLASHED_REdig_report_CC4.txt

# [telenius@deva CC4]$ ls PREfiltered/PREfiltered_test_040316_J_logFiles/
# FLASHED_REdig_report2_CC4.txt  FLASHED_REdig_report_CC4.txt  NONFLASHED_REdig_report2_CC4.txt  NONFLASHED_REdig_report_CC4.txt

    
# FASTQ reports ------------------------------------------------------------
   
    echo "<h4>FASTQC results here : </h4>" >> temp_description.html

    echo "<li>FastQC results (untrimmed) : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/READ1_fastqc_ORIGINAL/fastqc_report.html\" >READ1_fastqc_ORIGINAL/fastqc_report.html</a>   , and " >> temp_description.html
    echo " <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/READ2_fastqc_ORIGINAL/fastqc_report.html\" >READ2_fastqc_ORIGINAL/fastqc_report.html</a>  </li>" >> temp_description.html
   
    echo "<li>FastQC results (trimmed) : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/READ1_fastqc_TRIMMED/fastqc_report.html\" >READ1_fastqc_TRIMMED/fastqc_report.html</a>   , and " >> temp_description.html
    echo " <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/READ2_fastqc_TRIMMED/fastqc_report.html\" >READ2_fastqc_TRIMMED/fastqc_report.html</a>  </li>" >> temp_description.html 
  
    echo "<li>FastQC results (flash-combined) - before RE digestion: <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/FLASHED_fastqc/fastqc_report.html\" >FLASHED_fastqc/fastqc_report.html</a> </li>" >> temp_description.html 
    echo "<li>FastQC results (non-flash-combined) - before RE digestion: <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/NONFLASHED_fastqc/fastqc_report.html\" >NONFLASHED_fastqc/fastqc_report.html</a> </li>" >> temp_description.html 

    echo "<li>FastQC results (flash-combined) - after RE digestion: <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/FLASHED_REdig_fastqc/fastqc_report.html\" >FLASHED_REdig_fastqc/fastqc_report.html</a> </li>" >> temp_description.html 
    echo "<li>FastQC results (non-flash-combined) - after RE digestion: <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/NONFLASHED_REdig_fastqc/fastqc_report.html\" >NONFLASHED_REdig_fastqc/fastqc_report.html</a> </li>" >> temp_description.html 

   
    echo "<hr />" >> temp_description.html

# TRIMMING FLASHING RE-DIGESTION reports ------------------------------------------------------------
   
    echo "<h4>Trimming/flashing/RE-digestion log files here : </h4>" >> temp_description.html
    echo "<li>Harsh trim_galore trim : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_read_trimming.log\" >read_trimming.log</a>  </li>" >> temp_description.html
    echo "<li>Flashing : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_flashing.log\" >flashing.log</a>  </li>" >> temp_description.html
    echo "<li>Histogram of flashed reads : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_flash.hist\" >flash.hist</a>  </li>" >> temp_description.html
    echo "<li>RE digestion for of flash-combined reads : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_FLASHED_${REenzyme}digestion.log\" >FLASHED_${REenzyme}digestion.log</a>  </li>" >> temp_description.html
    echo "<li>RE digestion for of non-flash-combined reads : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_NONFLASHED_${REenzyme}digestion.log\" >NONFLASHED_${REenzyme}digestion.log</a>  </li>" >> temp_description.html

    echo "<hr />" >> temp_description.html

# CC-analyser reports ------------------------------------------------------------

    echo "<h4>Step-by-step analysis reports (in \"chronological order\" ) : </h4>" >> temp_description.html
    
    echo "Run log files available in : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_${QSUBOUTFILE}\" >${QSUBOUTFILE}</a> , and <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_${QSUBERRFILE}\" >${QSUBERRFILE}</a>"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    
    echo "<b>RE digestion of flash-combined reads</b> - throws out all reads with no RE cut : <a target="_blank" href=\"${ServerAndPath}/${Sample}_logFiles/${Sample}_FLASHED_${REenzyme}digestion.log\" >FLASHED_${REenzyme}digestion.log</a>  " >> temp_description.html
    echo "<br>"  >> temp_description.html

    echo "<br>"  >> temp_description.html
    echo "<b>Red graph</b> - CCanalyser results with duplicate filter turned OFF (before blat+ploidy filter)"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    echo "<li><b>Red graph, flashed reads</b> Capture script log files :"
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/FLASHED_REdig_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/FLASHED_REdig_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/FLASHED_REdig_report2_${CCversion}.txt\"  >RE cut report</a> ) "  >> temp_description.html
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/FLASHED_REdig_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html
    echo "<li><b>Red graph, non-flashed reads</b> Capture script log files :"  >> temp_description.html
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/NONFLASHED_REdig_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/NONFLASHED_REdig_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/NONFLASHED_REdig_report2_${CCversion}.txt\" >RE cut report</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/RAW/RAW_${Sample}_logFiles/NONFLASHED_REdig_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html  
    echo "<br>"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    
    
    echo "<b>Orange graph</b> - CCanalyser results with duplicate filter turned ON (before blat+ploidy filter)"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    echo "<li><b>Orange graph, flashed reads</b> Capture script log files :" >> temp_description.html   
    echo " ( <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/FLASHED_REdig_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/FLASHED_REdig_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/FLASHED_REdig_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html
    echo "<li><b>Orange graph, non-flashed reads</b> Capture script log files :" >> temp_description.html   
    echo " ( <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/NONFLASHED_REdig_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/NONFLASHED_REdig_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/NONFLASHED_REdig_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html

    echo "<li><b>Non-flashed reads coordinate adjustment in duplicate filter</b> : <a target="_blank" href=\"${ServerAndPath}/PREfiltered/PREfiltered_${Sample}_logFiles/NONFLASHED_REdig_report2_${CCversion}.txt\" >dupl_filtered_nonflashed_report2_${CCversion}.txt</a> </li>"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    
    
    echo "Ploidy+Blat-filtering log file : <a target="_blank" href=\"${ServerAndPath}/PREfiltered/filtering.log\" >filtering.log</a>"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    
    
    echo "<b>Green graph</b> - CCanalyser results for filtered data (duplicate, blat, ploidy filtered)"  >> temp_description.html
    echo "<br>"  >> temp_description.html
    echo "<li><b>Green graph in flashed reads overlay track</b> Capture script log files :"  >> temp_description.html   
    echo " ( <a target="_blank" href=\"${ServerAndPath}/FILTERED/FILTERED_${Sample}_logFiles/FLASHED_REdig_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/FILTERED/FILTERED_${Sample}_logFiles/FLASHED_REdig_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/FILTERED/FILTERED_${Sample}_logFiles/FLASHED_REdig_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html
    
    echo "<li><b>Green graph in non-flashed reads overlay track</b> Capture script log files :"  >> temp_description.html   
    echo " ( <a target="_blank" href=\"${ServerAndPath}/FILTERED/FILTERED_${Sample}_logFiles/NONFLASHED_REdig_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/FILTERED/FILTERED_${Sample}_logFiles/NONFLASHED_REdig_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/FILTERED/FILTERED_${Sample}_logFiles/NONFLASHED_REdig_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html

    echo "<br>"  >> temp_description.html
    echo "<li><b>Green graph, combined flashed and nonflashed filtered reads</b> Capture script log files :" >> temp_description.html   
    echo " ( <a target="_blank" href=\"${ServerAndPath}/COMBINED/COMBINED_${Sample}_logFiles/COMBINED_report4_${CCversion}.txt\" >Final REPORTER counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/COMBINED/COMBINED_${Sample}_logFiles/COMBINED_report3_${CCversion}.txt\" >Final counts</a> ) "  >> temp_description.html    
    echo " ( <a target="_blank" href=\"${ServerAndPath}/COMBINED/COMBINED_${Sample}_logFiles/COMBINED_report_${CCversion}.txt\"  >FULL report</a> ) "  >> temp_description.html    
    echo "</li>"  >> temp_description.html
  
    
    echo "<hr />" >> temp_description.html

    echo "</body>" > end.html
    echo "</html>"  >> end.html
    
    cat begin.html temp_description.html end.html > "${sampleForCCanalyser}_description.html"
    rm -f begin.html temp_description.html end.html
    
    # Moving the description file
    mv -f "${sampleForCCanalyser}_description.html" "${publicPathForCCanalyser}/."
    
    updateHub_part2c
    
    updateHub_part3 
    
}

