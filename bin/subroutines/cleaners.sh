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

cleanCCfolder(){
rm -f *_coordstring_${CCversion}.txt

TEMPweHaveSams=$(($(ls -1 *.sam | grep -c "")))
if [ "${TEMPweHaveSams}" -ne 0 ]; then
for file in *.sam
do
    TEMPreturnvalue=0
    bamname=$( echo $file | sed 's/.sam/.bam/' )
    if [ -s ${file} ]
    then
    samtools view -bh ${file} > ${bamname}
    TEMPreturnvalue=$?
        if [ ! -s ${bamname} ]
        then
            rm -f ${bamname}
        fi    
    fi
    
    if [ "${TEMPreturnvalue}" -eq 0 ] ; then
        rm -f $file
    else
        printThis="Couldn't sam-->bam transform file $file . Leaving it as SAM file."
        printToLogFile
        ls -lht ${file}
    fi
    ls -lht ${bamname}
done
fi
}

cleanUpRunFolder(){
    
# We want to leave somewhat explore-able structure to the output folder..

# |-- F1_test_040316_G_preCC4
# |-- F2_RAW_test_040316_G_CC4
# |-- F3_PREfiltered_test_040316_G_CC4
# |-- F4_filteringLogFor_test_040316_G_CC4
# `-- F5_FILTERED_test_040316_G_CC4

printThis="Cleaning up after ourselves - renaming folders and packing files.."
printToLogFile

# This one is already ready :
# |-- F1_test_040316_G_preCC4

# Changing names of these three :
# |-- F2_RAW_test_040316_G_CC4
# |-- F3_PREfiltered_test_040316_G_CC4
# |-- F4_filteringLogFor_test_040316_G_CC4

mv -f filteringLogFor_PREfiltered_${Sample}_${CCversion} F4_blatPloidyFilteringLog_${Sample}_${CCversion}

# Packing files.. 

cd F1_beforeCCanalyser_${Sample}_${CCversion}
echo F1_beforeCCanalyser_${Sample}_${CCversion}
echo "beforeCCanalyser folder contains trimming, flashing, REdigesting, bowtie mapping of the sample" > a_trim_flash_REdigest_bowtie_containing_folder

flashstatus="FLASHED"
echo "samtools view -hb ${flashstatus}_REdig.sam > ${flashstatus}_REdig.bam"
samtools view -hb ${flashstatus}_REdig.sam > ${flashstatus}_REdig.bam
flashstatus="NONFLASHED"
echo "samtools view -hb ${flashstatus}_REdig.sam > ${flashstatus}_REdig.bam"
samtools view -hb ${flashstatus}_REdig.sam > ${flashstatus}_REdig.bam

ls -lht FLASHED_REdig.bam
ls -lht NONFLASHED_REdig.bam
rm -f FLASHED_REdig.sam NONFLASHED_REdig.sam
cd ..


cd F2_redGraphs_${Sample}_${CCversion}
echo F2_redGraphs_${Sample}_${CCversion}
cleanCCfolder
echo "redGraphs folder is a CCanalyser run where duplicate filter is switched ON" > a_CCanalyser_run_with_duplicate_filter_switched_OFF
cd ..

cd F3_orangeGraphs_${Sample}_${CCversion}
echo F3_orangeGraphs_${Sample}_${CCversion}
echo "orangeGraphs folder is a CCanalyser run where duplicate filter is switched ON" > a_CCanalyser_run_with_duplicate_filter_switched_ON
cleanCCfolder
cd ..

thisIsWhereIam=$(pwd)
cd F4_blatPloidyFilteringLog_${Sample}_${CCversion}/BlatPloidyFilterRun/BLAT_PLOIDY_FILTERED_OUTPUT
cleanCCfolder
cd ${thisIsWhereIam}

cd F5_greenGraphs_separate_${Sample}_${CCversion}
echo F5_greenGraphs_separate_${Sample}_${CCversion}
echo "greenGraphs_separate is a CCanalyser re-run for blat+ploidy filtered data" > a_CCanalyser_run_for_blatPloidy_filtered_data
cleanCCfolder
cd ..

# Don't clean F6 if doing tri-c right after ..

if [ "${TRIC}" -eq "0" ] && [ "${TRIC_EXCL}" -eq "0" ]; then
    
cd F6_greenGraphs_combined_${Sample}_${CCversion}
echo F6_greenGraphs_combined_${Sample}_${CCversion}
echo "greenGraphs_combined is a CCanalyser run to combine flashed and nonflashed data" > a_CCanalyser_run_to_generate_final_results
cleanCCfolder
cd ..

fi

echo
echo "Output folders generated :"

ls -lht
    
}

cleanUpAfterTric(){

printThis="Cleaning up after tri-c : packing files.."
printToLogFile
   
cd F6_greenGraphs_combined_${Sample}_${CCversion}
echo F6_greenGraphs_combined_${Sample}_${CCversion}

if [ "${ONLY_TRIC}" -eq "0" ]; then
   cleanCCfolder
else
   rm -f COMBINED_reported_capture_reads_${CCversion}.sam
fi

cd ..
    
}


oneFolderSymLinks(){

# ${symLinkFolderToBe}

cd ${symLinkFolderToBe}

ls | grep ".bw$"
echo "mv -f *.bw ${BigwigsAreHere}/${symLinkFolderToBe}/."
mv -f *.bw ${BigwigsAreHere}/${symLinkFolderToBe}/.

echo "ln -s ${BigwigsAreHere}/${symLinkFolderToBe}/*.bw ."
ln -s ${BigwigsAreHere}/${symLinkFolderToBe}/*.bw .

cd ..   
    
}

makeSymbolicLinks(){
    
# Move bigwigs and generate symbolic links

printThis="Moving bigwigs and generating symbolic links.."
printToLogFile

RunDir=$( pwd )

mkdir PERMANENT_BIGWIGS_do_not_move
cd PERMANENT_BIGWIGS_do_not_move
mkdir RAW PREfiltered FILTERED
BigwigsAreHere=$( pwd )

cd ${PublicPath}

symLinkFolderToBe="RAW"
oneFolderSymLinks

symLinkFolderToBe="PREfiltered"
oneFolderSymLinks

symLinkFolderToBe="FILTERED"
oneFolderSymLinks

cd ${RunDir}
    
    
}
