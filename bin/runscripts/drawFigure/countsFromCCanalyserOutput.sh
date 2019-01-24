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

# targetDir="/full-path/to/the/run/dir"
# sample="samplename"

targetDir="$1"
Sample="$2"
CCversion="$3"
restrictionEnzyme="$4"


# All reads
# Old versions of Flash say "Total reads" when counting total read pairs. New versions of Flash say "Total pairs".
# all=$(($( cat ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/flashing.log | grep "Total reads:" | sed 's/.*:\s*//' )))
# all=$(($( cat ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/flashing.log | grep "Total pairs:" | sed 's/.*:\s*//' )))
all=$(($( cat ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/flashing.log | grep "Total [rp][ea][ai][dr]s:" | sed 's/.*:\s*//' )))

# Flashed (see notes for "All reads" above)

allflashed=$(($( cat  ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/flashing.log | grep "Combined [rp][ea][ai][dr]s:" | sed 's/.*:\s*//' )))
allnonflashed=$(($( cat ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/flashing.log | grep "Uncombined [rp][ea][ai][dr]s:" | sed 's/.*:\s*//' )))

# RE_site 	
REflashed=$(($( cat ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/FLASHED_${restrictionEnzyme}digestion.log | grep "had at least one ${restrictionEnzyme} site in them" | sed 's/\s.*//' )))
REnonflashed=$(($( cat ${targetDir}/F1_beforeCCanalyser_${Sample}_${CCversion}/NONFLASHED_${restrictionEnzyme}digestion.log | grep "had at least one ${restrictionEnzyme} site in them" | sed 's/\s.*//' )))

# continues 	
continuesToMappingFlashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^11\s" | sed 's/.*:\s*//' )))
continuesToMappingNonflashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^11\s" | sed 's/.*:\s*//' )))

# Cont_cap 	

containsCaptureFlashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^11b\s" | sed 's/.*:\s*//' )))
containsCaptureNonflashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^11b\s" | sed 's/.*:\s*//' )))

# Cap_and_rep 	

containsCapAndRepFlashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^11c\s" | sed 's/.*:\s*//' )))
containsCapAndRepNonflashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^11c\s" | sed 's/.*:\s*//' )))

# multicap 	

singleCapFlashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^11f\s" | sed 's/.*:\s*//' )))
singleCapNonflashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^11f\s" | sed 's/.*:\s*//' )))
multiCapFlashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^11g\s" | sed 's/.*:\s*//' )))
multiCapNonflashed=$(($( cat ${targetDir}/F2_redGraphs_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^11g\s" | sed 's/.*:\s*//' )))

# duplicate 

nonduplicateFlashed=$(($( cat ${targetDir}/F3_orangeGraphs_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^16\s" | sed 's/.*:\s*//' )))
nonduplicateNonflashed=$(($( cat ${targetDir}/F3_orangeGraphs_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^16\s" | sed 's/.*:\s*//' )))	

# blatploidy 	

blatploidyFlashed=$(($( cat ${targetDir}/F5_greenGraphs_separate_${Sample}_${CCversion}/FLASHED_REdig_report_${CCversion}.txt | grep "^11d\s" | sed 's/.*:\s*//' )))
blatploidyNonflashed=$(($( cat ${targetDir}/F5_greenGraphs_separate_${Sample}_${CCversion}/NONFLASHED_REdig_report_${CCversion}.txt | grep "^11d\s" | sed 's/.*:\s*//' )))

echo "all=${all}"
echo "allflashed=${allflashed}"
echo "allnonflashed=${allnonflashed}"
echo "REflashed=${REflashed}"
echo "REnonflashed=${REnonflashed}"
echo "continuesToMappingFlashed=${continuesToMappingFlashed}"
echo "continuesToMappingNonflashed=${continuesToMappingNonflashed}"
echo "containsCaptureFlashed=${containsCaptureFlashed}"
echo "containsCaptureNonflashed=${containsCaptureNonflashed}"
echo "containsCapAndRepFlashed=${containsCapAndRepFlashed}"
echo "containsCapAndRepNonflashed=${containsCapAndRepNonflashed}"
echo "singleCapFlashed=${singleCapFlashed}"
echo "singleCapNonflashed=${singleCapNonflashed}"
echo "multiCapFlashed=${multiCapFlashed}"
echo "multiCapNonflashed=${multiCapNonflashed}"
echo "nonduplicateFlashed=${nonduplicateFlashed}"
echo "nonduplicateNonflashed=${nonduplicateNonflashed}"
echo "blatploidyFlashed=${blatploidyFlashed}"
echo "blatploidyNonflashed=${blatploidyNonflashed}"



