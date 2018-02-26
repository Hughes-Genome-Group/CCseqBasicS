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

#------------------------------------------
# The codes of the pipeline 
#------------------------------------------
#
# CCseqBasic5/
#
# |
# |-- CCseqBasic5.sh
# |
# `-- bin
#     |
#     |-- runscripts
#     |   |
#     |   |-- analyseMappedReads.pl
#     |   |-- dpnIIcutGenome.pl
#     |   |-- nlaIIIcutGenome.pl
#     |   |-- dpnIIcutReads.pl
#     |   |-- nlaIIIcutReads.pl
#     |   |
#     |   |-- filterArtifactMappers
#     |   |   |
#     |   |   |-- 1_blat.sh
#     |   |   |-- 2_psl_parser.pl
#     |   |   `-- filter.sh
#     |   |
#     |   `-- drawFigure
#     |       |
#     |       |-- countsFromCCanalyserOutput.sh
#     |       |-- drawFigure.py
#     |       `-- generatePercentages.py
#     |   
#     `-- subroutines
#         |-- cleaners.sh
#         |-- hubbers.sh
#         |-- parametersetters.sh
#         |-- runtools.sh
#         |-- testers_and_loggers.sh
#         `-- usageAndVersion.sh

#------------------------------------------

function finish {
if [ $? != "0" ]; then
echo
echo "RUN CRASHED ! - check qsub.err to see why !"
echo
echo "If your run passed folder1 (F1) succesfully - i.e. you have F2 or later folders formed correctly - you can restart in same folder, same run.sh :"
echo "Just add --onlyCCanalyser to the end of run command in run.sh, and start the run normally, in the same folder you crashed now (this will overrwrite your run from bowtie output onwards)."
echo
echo "If you are going to rerun a crashed run without using --onlyCCanalyser , copy your run script to a NEW EMPTY FOLDER,"
echo "and remember to delete your malformed /public/ hub-folders (especially the tracks.txt files) to avoid wrongly generated data hubs (if you are going to use same SAMPLE NAME as in the crashed run)" 
echo

else
echo
echo "Analysis complete !"
date
fi
}
trap finish EXIT

#------------------------------------------

QSUBOUTFILE="qsub.out"
QSUBERRFILE="qsub.err"

OligoFile=""
TRIM=1
GENOME=""
WINDOW=200
INCREMENT=20
CAPITAL_M=0
LOWERCASE_M=0
BOWTIEMEMORY="256"
Sample="sample"
Read1=""
Read2=""

CUSTOMAD=-1
ADA31="no"
ADA32="no"

# trimgalore default
QMIN=20

# flash defaults
flashOverlap=10
flashErrorTolerance=0.25

saveDpnGenome=0

ucscBuild=""
otherBowtieParameters=""
bowtieMismatchBehavior=""

otherParameters=""
PublicPath="UNDETERMINED"

ploidyFilter=""
extend=20000

# Blat flags
stepSize=5 # Jon default - James used blat default, which is "tileSize", in this case thus 11 (this was the setting in CC2 and CC3 - i.e filter versions VS101 and VS102)
tileSize=11 # Jon, James default
minScore=10 # Jon new default. Jon default before2016 and CC4 default until 080916 minScore=30 - James used minScore=30 (this was the setting in CC2 and CC3 - i.e filter versions VS101 and VS102)
maxIntron=4000 # blat default 750000- James used maxIntron=4000 (this was the setting in CC2 and CC3 - i.e filter versions VS101 and VS102)
oneOff=0 # oneOff=1 would allow 1 mismatch in tile (blat default = 0 - that is also CC3 and CC2 default)

# Whether we reuse blat results from earlier run ..
# Having this as "." will search from the run dir when blat is ran - so file will not be found, and thus BLAT will be ran normally.
reuseBLATpath="."

REenzyme="dpnII"

# Skip other stages - assume input from this run has been ran earlier - to construct to THIS SAME FOLDER everything else
# but as the captureC analyser naturally crashed - this will jump right to the beginning of that part..
ONLY_CC_ANALYSER=0
# Rerun public folder generation and filling. Will not delete existing folder, but will overwrite all files (and start tracks.txt from scratch).
ONLY_HUB=0

strandSpecificDuplicates=0

#------------------------------------------

CCversion="CS5"
captureScript="analyseMappedReads"
CCseqBasicVersion="CCseqBasic5"

echo "${CCseqBasicVersion}.sh - by Jelena Telenius, 05/01/2016"
echo
timepoint=$( date )
echo "run started : ${timepoint}"
echo
echo "Script located at"
echo "$0"
echo

echo "RUNNING IN MACHINE : "
hostname --long

echo "run called with parameters :"
echo "${CCseqBasicVersion}.sh" $@
echo

#------------------------------------------

# Loading subroutines in ..

echo "Loading subroutines in .."

CaptureTopPath="$( echo $0 | sed 's/\/'${CCseqBasicVersion}'.sh$//' )"

CapturePipePath="${CaptureTopPath}/bin/subroutines"

# HUBBING subroutines
. ${CapturePipePath}/hubbers.sh

# SETTING parameter values - subroutines
. ${CapturePipePath}/parametersetters.sh

# CLEANING folders and organising structures
. ${CapturePipePath}/cleaners.sh

# TESTING file existence, log file output general messages
. ${CapturePipePath}/testers_and_loggers.sh

# RUNNING the main tools (flash, ccanalyser, etc..)
. ${CapturePipePath}/runtools.sh

# SETTING THE GENOME BUILD PARAMETERS
. ${CapturePipePath}/genomeSetters.sh

# SETTING THE BLACKLIST GENOME LIST PARAMETERS
. ${CapturePipePath}/blacklistSetters.sh

# PRINTING HELP AND VERSION MESSAGES
. ${CapturePipePath}/usageAndVersion.sh

#------------------------------------------

# From where to call the main scripts operating from the runscripts folder..

RunScriptsPath="${CaptureTopPath}/bin/runscripts"

#------------------------------------------

# From where to call the filtering scripts..
# (blacklisting regions with BLACKLIST pre-made region list, as well as on-the-fly BLAT-hit based "false positive" hits) 

CaptureFilterPath="${RunScriptsPath}/filterArtifactMappers"

#------------------------------------------

# From where to call the python plots..
# (blacklisting regions with BLACKLIST pre-made region list, as well as on-the-fly BLAT-hit based "false positive" hits) 

CapturePlotPath="${RunScriptsPath}/drawFigure"

#------------------------------------------

# From where to call the CONFIGURATION script..

confFolder="${CaptureTopPath}/conf"

#------------------------------------------

echo
echo "CaptureTopPath ${CaptureTopPath}"
echo "CapturePipePath ${CapturePipePath}"
echo "confFolder ${confFolder}"
echo "RunScriptsPath ${RunScriptsPath}"
echo "CaptureFilterPath ${CaptureFilterPath}"
echo

#------------------------------------------

# Calling in the CONFIGURA`TION script and its default setup :

# Defaulting this to "not in use" - if it is not set in the config file.
CaptureDigestPath="NOT_IN_USE"

#------------------------------------------

# Calling in the CONFIGURATION script and its default setup :

echo "Calling in the conf/config.sh script and its default setup .."

CaptureDigestPath="NOT_IN_USE"
supportedGenomes=()
BOWTIE1=()
UCSC=()
BLACKLIST=()
genomesWhichHaveBlacklist=()


# . ${confFolder}/config.sh
. ${confFolder}/genomeBuildSetup.sh
. ${confFolder}/loadNeededTools.sh
. ${confFolder}/serverAddressAndPublicDiskSetup.sh

# setConfigLocations
setPathsForPipe
setGenomeLocations

echo 
echo "Supported genomes : "
for g in $( seq 0 $((${#supportedGenomes[@]}-1)) ); do echo -n "${supportedGenomes[$g]} "; done
echo 
echo

echo 
echo "Blacklist filtering available for these genomes : "
for g in $( seq 0 $((${#genomesWhichHaveBlacklist[@]}-1)) ); do echo -n "${genomesWhichHaveBlacklist[$g]} "; done
echo 
echo

echo "Calling in the conf/serverAddressAndPublicDiskSetup.sh script and its default setup .."

SERVERTYPE="UNDEFINED"
SERVERADDRESS="UNDEFINED"
REMOVEfromPUBLICFILEPATH="NOTHING"
ADDtoPUBLICFILEPATH="NOTHING"
tobeREPLACEDinPUBLICFILEPATH="NOTHING"
REPLACEwithThisInPUBLICFILEPATH="NOTHING"

. ${confFolder}/serverAddressAndPublicDiskSetup.sh

setPublicLocations

echo
echo "SERVERTYPE ${SERVERTYPE}"
echo "SERVERADDRESS ${SERVERADDRESS}"
echo "ADDtoPUBLICFILEPATH ${ADDtoPUBLICFILEPATH}"
echo "REMOVEfromPUBLICFILEPATH ${REMOVEfromPUBLICFILEPATH}"
echo "tobeREPLACEDinPUBLICFILEPATH ${tobeREPLACEDinPUBLICFILEPATH}"
echo "REPLACEwithThisInPUBLICFILEPATH ${REPLACEwithThisInPUBLICFILEPATH}"
echo

#------------------------------------------

OPTS=`getopt -o h,m:,M:,o:,s:,w:,i:,v: --long help,dump,snp,dpn,nla,strandSpecificDuplicates,UMI,CCversion:,BLATforREUSEfolderPath:,globin:,outfile:,errfile:,limit:,pf:,genome:,R1:,R2:,saveGenomeDigest,dontSaveGenomeDigest,trim,noTrim,chunkmb:,window:,increment:,ada3read1:,ada3read2:,extend:,onlyCCanalyser,onlyHub,noPloidyFilter:,qmin:,flashBases:,flashMismatch:,stringent,trim3:,trim5:,seedmms:,seedlen:,maqerr:,stepSize:,tileSize:,minScore:,maxIntron:,oneOff:,wobblyEndBinWidth: -- "$@"`
if [ $? != 0 ]
then
    exit 1
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -h) usage ; shift;;
        -m) LOWERCASE_M=$2 ; shift 2;;
        -M) CAPITAL_M=$2 ; shift 2;;
        -o) OligoFile=$2 ; shift 2;;
        -w) WINDOW=$2 ; shift 2;;
        -i) INCREMENT=$2 ; shift 2;;
        -s) Sample=$2 ; shift 2;;
        -v) LOWERCASE_V="$2"; shift 2;;
        --help) usage ; shift;;
        --UMI) otherParameters="$otherParameters --umi" ; shift;;
        --CCversion) CCversion="$2"; shift 2;;       
        --dpn) REenzyme="dpnII" ; shift;;
        --nla) REenzyme="nlaIII" ; shift;;
        --onlyCCanalyser) ONLY_CC_ANALYSER=1 ; shift;;
        --onlyHub) ONLY_HUB=1 ; shift;;
        --R1) Read1=$2 ; shift 2;;
        --R2) Read2=$2 ; shift 2;;
        --chunkmb) BOWTIEMEMORY=$2 ; shift 2;;
        --saveGenomeDigest) saveDpnGenome=1 ; shift;;
        --dontSaveGenomeDigest) saveDpnGenome=0 ; shift;;
        --trim) TRIM=1 ; shift;;
        --noTrim) TRIM=0 ; shift;;
        --window) WINDOW=$2 ; shift 2;;
        --increment) INCREMENT=$2 ; shift 2;;
        --genome) GENOME=$2 ; shift 2;;
        --ada3read1) ADA31=$2 ; shift 2;;
        --ada3read2) ADA32=$2 ; shift 2;;
        --extend) extend=$2 ; shift 2;;
        --noPloidyFilter) ploidyFilter="--noploidyfilter " ; shift;;
        --strandSpecificDuplicates) otherParameters="$otherParameters --stranded"; strandSpecificDuplicates=1 ; shift;;
        --dump) otherParameters="$otherParameters --dump" ; shift;;
        --snp) otherParameters="$otherParameters --snp" ; shift;;
        --globin) otherParameters="$otherParameters --globin $2" ; shift 2;;
        --limit) otherParameters="$otherParameters --limit $2" ; shift 2;;
        --stringent) otherParameters="$otherParameters --stringent" ; shift 1;;
        --pf) PublicPath="$2" ; shift 2;;
        --qmin) QMIN="$2" ; shift 2;;
        --BLATforREUSEfolderPath) reuseBLATpath="$2" ; shift 2;;
        --flashBases) flashOverlap="$2" ; shift 2;;
        --flashMismatch) flashErrorTolerance="$2" ; shift 2;;
        --trim3) otherBowtieParameters="${otherBowtieParameters} --trim3 $2 " ; shift 2;;
        --trim5) otherBowtieParameters="${otherBowtieParameters} --trim5 $2 " ; shift 2;;
        --seedmms) bowtieMismatchBehavior="${bowtieMismatchBehavior} --seedmms $2 " ; shift 2;;
        --seedlen) bowtieMismatchBehavior="${bowtieMismatchBehavior} --seedlen $2 " ; shift 2;;
        --maqerr) bowtieMismatchBehavior="${bowtieMismatchBehavior} --maqerr $2 " ; shift 2;;
        --stepSize) stepSize=$2 ; shift 2;;
        --tileSize) tileSize==$2 ; shift 2;;
        --minScore) minScore=$2 ; shift 2;;
        --maxIntron) maxIntron=$2 ; shift 2;;
        --oneOff) oneOff=$2 ; shift 2;;
        --outfile) QSUBOUTFILE=$2 ; shift 2;;
        --errfile) QSUBERRFILE=$2 ; shift 2;;
        --wobblyEndBinWidth) otherParameters="$otherParameters --wobble $2" ; shift 2;;
        --) shift; break;;
    esac
done

# ----------------------------------------------

# Setting the duplicate filter style !



if [ ${CCversion} == "CS5" ] ; then
    echo
    echo "Duplicate filtering style CS5 selected ! "
    echo
elif [ ${CCversion} == "CS3" ] ; then
    echo
    echo "Duplicate filtering style CS3 selected ! "
    echo
elif [ ${CCversion} == "CS4" ] ; then
    echo
    echo "Duplicate filtering style CS4 selected ! "
    echo
else
   # Crashing here !
    printThis="Duplicate filtering style given wrong ! Give either --CCversion CS3 or --CCversion CS4 ( or default --CCversion CS5 )"
    printToLogFile
    printThis="You gave --CCversion ${CCversion}"
    printToLogFile
    printThis="EXITING ! "
    printToLogFile
    exit 1
fi


# ----------------------------------------------

echo "Parsing the data area and server locations .."

PublicPath="${PublicPath}/${Sample}/${CCversion}_${REenzyme}"

# Here, parsing the data area location, to reach the public are address..
diskFolder=${PublicPath}
serverFolder=""   
echo
parsePublicLocations
echo

tempJamesUrl="${SERVERADDRESS}/${serverFolder}"
JamesUrl=$( echo ${tempJamesUrl} | sed 's/\/\//\//g' )
ServerAndPath="${SERVERTYPE}://${JamesUrl}"

# ----------------------------------------------

# Setting artificial chromosome on, if we have it .

if [ ${GENOME} == "mm9PARP" ] ; then

# Whether we have artificial chromosome chrPARP or not, to feed to analyseMappedReads.pl (to be filtered out before visualisation)
# Will be turned on based on genome name, to become :
otherParameters="$otherParameters --parp"

fi

# ----------------------------------------------

# Modifying and adjusting parameter values, based on run flags

setBOWTIEgenomeSizes
setGenomeFasta

echo "GenomeFasta ${GenomeFasta}" >> parameters_capc.log
echo "BowtieGenome ${BowtieGenome}" >> parameters_capc.log

setUCSCgenomeSizes

echo "ucscBuild ${ucscBuild}" >> parameters_capc.log

#------------------------------------------

CaptureDigestPath="${CaptureDigestPath}/${REenzyme}"

setParameters

# ----------------------------------------------

# Loading the environment - either with module system or setting them into path.
# This subroutine comes from conf/config.sh file

printThis="LOADING RUNNING ENVIRONMENT"
printToLogFile

setPathsForPipe

#---------------------------------------------------------

# Check that the requested RE actually exists ..

if [ ! -s ${RunScriptsPath}/${REenzyme}cutReads4.pl ] || [ ! -s ${RunScriptsPath}/${REenzyme}cutGenome4.pl ] ; then

printThis="EXITING ! - Restriction enzyme ${REenzyme} is not supported (check your spelling)"
exit 1
   
fi

#---------------------------------------------------------

echo "Run with parameters :"
echo
echo "Output log file ${QSUBOUTFILE}" > parameters_capc.log
echo "Output error log file ${QSUBERRFILE}" >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "CaptureTopPath ${CaptureTopPath}" >> parameters_capc.log
echo "CapturePipePath ${CapturePipePath}" >> parameters_capc.log
echo "confFolder ${confFolder}" >> parameters_capc.log
echo "RunScriptsPath ${RunScriptsPath}" >> parameters_capc.log
echo "CaptureFilterPath ${CaptureFilterPath}" >> parameters_capc.log
echo "CaptureDigestPath ${CaptureDigestPath}" >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "Sample ${Sample}" >> parameters_capc.log
echo "Read1 ${Read1}" >> parameters_capc.log
echo "Read2 ${Read2}" >> parameters_capc.log
echo "GENOME ${GENOME}" >> parameters_capc.log
echo "GenomeIndex ${GenomeIndex}" >> parameters_capc.log
echo "OligoFile ${OligoFile}" >> parameters_capc.log
echo "REenzyme ${REenzyme}" >> parameters_capc.log
echo "ONLY_CC_ANALYSER ${ONLY_CC_ANALYSER}" >> parameters_capc.log

echo "------------------------------" >> parameters_capc.log
echo "TRIM ${TRIM}  (TRUE=1, FALSE=0)" >> parameters_capc.log
echo "QMIN ${QMIN}  (default 20)" >> parameters_capc.log

echo "CUSTOMAD ${CUSTOMAD}   (TRUE=1, FALSE= -1)"  >> parameters_capc.log

if [ "${CUSTOMAD}" -ne -1 ]; then

echo "ADA31 ${ADA31}"  >> parameters_capc.log
echo "ADA32 ${ADA32}"  >> parameters_capc.log
   
fi

echo "------------------------------" >> parameters_capc.log
echo "flashOverlap ${flashOverlap} (default 10)"  >> parameters_capc.log
echo "flashErrorTolerance ${flashErrorTolerance} (default 0.25)"  >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "saveDpnGenome ${saveDpnGenome}  (TRUE=1, FALSE=0)" >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "BOWTIEMEMORY ${BOWTIEMEMORY}"  >> parameters_capc.log
echo "CAPITAL_M ${CAPITAL_M}" >> parameters_capc.log
echo "LOWERCASE_M ${LOWERCASE_M}" >> parameters_capc.log
echo "otherBowtieParameters ${otherBowtieParameters}"  >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "reuseBLATpath ${reuseBLATpath}" >> parameters_capc.log
echo "stepSize ${stepSize}" >> parameters_capc.log
echo "tileSize ${tileSize}" >> parameters_capc.log
echo "minScore ${minScore}" >> parameters_capc.log
echo "maxIntron ${maxIntron}" >> parameters_capc.log
echo "oneOff ${oneOff}" >> parameters_capc.log
echo "extend ${extend}"  >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "ploidyFilter ${ploidyFilter}"  >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "WINDOW ${WINDOW}" >> parameters_capc.log
echo "INCREMENT ${INCREMENT}" >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "PublicPath ${PublicPath}" >> parameters_capc.log
echo "ServerUrl ${SERVERADDRESS}" >> parameters_capc.log
echo "JamesUrl ${JamesUrl}" >> parameters_capc.log
echo "ServerAndPath ${ServerAndPath}" >> parameters_capc.log
echo "otherParameters ${otherParameters}" >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "GenomeFasta ${GenomeFasta}" >> parameters_capc.log
echo "BowtieGenome ${BowtieGenome}" >> parameters_capc.log
echo "ucscBuild ${ucscBuild}" >> parameters_capc.log


cat parameters_capc.log
echo

echo "Whole genome fasta file path : ${GenomeFasta}"
echo "Bowtie genome index path : ${BowtieGenome}"
echo "Chromosome sizes for UCSC bigBed generation will be red from : ${ucscBuild}"


testedFile="${OligoFile}"
doInputFileTesting

# Making output folder.. (and crashing run if found it existing from a previous crashed run)
if [[ ${ONLY_HUB} -eq "0" ]]; then
if [[ ${ONLY_CC_ANALYSER} -eq "0" ]]; then

if [ -d F1_beforeCCanalyser_${Sample}_${CCversion} ] ; then
  # Crashing here !
  printThis="EXITING ! Previous run data found in run folder ! - delete data of previous run (or define rerun with --onlyCCanalyser )"
  printToLogFile
  exit 1
  
fi
    
mkdir F1_beforeCCanalyser_${Sample}_${CCversion}   
fi
fi

# Here crashing if public folder exists (and this is not --onlyCCanalyser run ..

if [ -d ${PublicPath} ] && [ ${ONLY_CC_ANALYSER} -eq "0" ] ; then
    # Allows to remove if it is empty..
    rmdir ${PublicPath}

if [ -d ${PublicPath} ] ; then
   # Crashing here !
  printThis="EXITING ! Existing public data found in folder ${PublicPath} "
  printToLogFile
  printThis="Delete the data before restarting the script (refusing to overwrite) "
  printToLogFile
  exit 1
fi

fi


if [[ ${ONLY_HUB} -eq "0" ]]; then
if [[ ${ONLY_CC_ANALYSER} -eq "0" ]]; then

# Copy files over..

testedFile="${Read1}"
doInputFileTesting
testedFile="${Read2}"
doInputFileTesting

if [ "${Read1}" != "READ1.fastq" ] ; then
printThis="Copying input file R1.."
printToLogFile
cp "${Read1}" F1_beforeCCanalyser_${Sample}_${CCversion}/READ1.fastq
else
printThis="Making safety copy of the original READ1.fastq : READ1.fastq_original.."
printToLogFile
cp "${Read1}" F1_beforeCCanalyser_${Sample}_${CCversion}/READ1.fastq_original
fi
doQuotaTesting

if [ "${Read2}" != "READ2.fastq" ] ; then
printThis="Copying input file R2.."
printToLogFile
cp "${Read2}" F1_beforeCCanalyser_${Sample}_${CCversion}/READ2.fastq
else
printThis="Making safety copy of the original READ2.fastq : READ2.fastq_original.."
printToLogFile
cp "${Read2}" F1_beforeCCanalyser_${Sample}_${CCversion}/READ2.fastq_original
fi
doQuotaTesting

testedFile="F1_beforeCCanalyser_${Sample}_${CCversion}/READ1.fastq"
doTempFileTesting
testedFile="F1_beforeCCanalyser_${Sample}_${CCversion}/READ2.fastq"
doTempFileTesting

# Save oligo file full path (to not to lose the file when we cd into the folder, if we used relative paths ! )
TEMPdoWeStartWithSlash=$(($( echo ${OligoFile} | awk '{print substr($1,1,1)}' | grep -c '/' )))
if [ "${TEMPdoWeStartWithSlash}" -eq 0 ]
then
 OligoFile=$(pwd)"/"${OligoFile}
fi

testedFile="${OligoFile}"
doInputFileTesting

fi
fi

# Go into output folder..
cd F1_beforeCCanalyser_${Sample}_${CCversion}

if [[ ${ONLY_HUB} -eq "0" ]]; then
if [[ ${ONLY_CC_ANALYSER} -eq "0" ]]; then

################################################################
#Check BOWTIE quality scores..

printThis="Checking the quality score scheme of the fastq files.."
printToLogFile
    
    bowtieQuals=""
    LineCount=$(($( grep -c "" READ1.fastq )/4))
    if [ "${LineCount}" -gt 100000 ] ; then
        bowtieQuals=$( perl ${RunScriptsPath}/fastq_scores_bowtie1.pl -i READ1.fastq -r 90000 )
    else
        rounds=$((${LineCount}-10))
        bowtieQuals=$( perl ${RunScriptsPath}/fastq_scores_bowtie1.pl -i READ1.fastq -r ${rounds} )
    fi
    
    echo "Flash, Trim_galore and Bowtie will be ran in quality score scheme : ${bowtieQuals}"

    # The location of "zero" for the filtering/trimming programs cutadapt, trim_galore, flash    
    intQuals=""
    if [ "${bowtieQuals}" == "--phred33-quals" ] ; then
        intQuals="33"
    else
        # Both solexa and illumina phred64 have their "zero point" in 64
        intQuals="64"
    fi

################################################################
# Fastq for original files..
printThis="Running fastQC for input files.."
printToLogFile

printThis="${RunScriptsPath}/QC_and_Trimming.sh --fastqc"
printToLogFile

${RunScriptsPath}/QC_and_Trimming.sh --fastqc

    # Changing names of fastqc folders to be "ORIGINAL"
    mv -f READ1_fastqc READ1_fastqc_ORIGINAL
    mv -f READ2_fastqc READ2_fastqc_ORIGINAL
    mv -f READ1_fastqc.zip READ1_fastqc_ORIGINAL.zip
    mv -f READ2_fastqc.zip READ2_fastqc_ORIGINAL.zip
   
    ls -lht

################################################################
# Trimgalore for the reads..

if [[ ${TRIM} -eq "1" ]]; then

printThis="Running trim_galore for the reads.."
printToLogFile

printThis="${RunScriptsPath}/QC_and_Trimming.sh -q ${intQuals} --filter 3 --qmin ${QMIN}"
printToLogFile

${RunScriptsPath}/QC_and_Trimming.sh -q "${intQuals}" --filter 3 --qmin ${QMIN}

doQuotaTesting
ls -lht

testedFile="READ1.fastq"
doTempFileTesting
testedFile="READ2.fastq"
doTempFileTesting

################################################################
# Fastq for trimmed files..
printThis="Running fastQC for trimmed files.."
printToLogFile

printThis="${RunScriptsPath}/QC_and_Trimming.sh --fastqc"
printToLogFile

${RunScriptsPath}/QC_and_Trimming.sh --fastqc

    # Changing names of fastqc folders to be "TRIMMED"
    mv -f READ1_fastqc READ1_fastqc_TRIMMED
    mv -f READ2_fastqc READ2_fastqc_TRIMMED
    mv -f READ1_fastqc.zip READ1_fastqc_TRIMMED.zip
    mv -f READ2_fastqc.zip READ2_fastqc_TRIMMED.zip
    
fi
    
################################################################
# FLASH for trimmed files..
printThis="Running FLASH for trimmed files.."
printToLogFile

runFlash

ls -lht
doQuotaTesting

rm -f READ1.fastq READ2.fastq

################################################################
# Fastq for flashed files..
printThis="Running fastQC for FLASHed and nonflashed files.."
printToLogFile

printThis="fastqc --quiet -f fastq FLASHED.fastq"
printToLogFile

fastqc --quiet -f fastq FLASHED.fastq


printThis="fastqc --quiet -f fastq NONFLASHED.fastq"
printToLogFile

fastqc --quiet -f fastq NONFLASHED.fastq


################################################################

# Running dpnII digestion for flashed file..
printThis="Running ${REenzyme} digestion for flashed file.."
printToLogFile

printThis="perl ${RunScriptsPath}/${REenzyme}cutReads4.pl FLASHED.fastq FLASHED"
printToLogFile

perl ${RunScriptsPath}/${REenzyme}cutReads4.pl FLASHED.fastq FLASHED > FLASHED_${REenzyme}digestion.log
cat FLASHED_${REenzyme}digestion.log
ls -lht
doQuotaTesting

testedFile="FLASHED_REdig.fastq"
doTempFileTesting
rm -f FLASHED.fastq

# Running dpnII digestion for non-flashed file..
printThis="Running ${REenzyme} digestion for non-flashed file.."
printToLogFile

printThis="perl ${RunScriptsPath}/${REenzyme}cutReads4.pl NONFLASHED.fastq NONFLASHED"
printToLogFile

perl ${RunScriptsPath}/${REenzyme}cutReads4.pl NONFLASHED.fastq NONFLASHED > NONFLASHED_${REenzyme}digestion.log
cat NONFLASHED_${REenzyme}digestion.log

 ls -lht
 doQuotaTesting
 
testedFile="NONFLASHED_REdig.fastq"
doTempFileTesting
rm -f NONFLASHED.fastq

################################################################
# Fastq for flashed files..
printThis="Running fastQC for RE-digested files.."
printToLogFile

printThis="fastqc --quiet -f fastq FLASHED_REdig.fastq"
printToLogFile

fastqc --quiet -f fastq FLASHED_REdig.fastq


printThis="fastqc --quiet -f fastq NONFLASHED_REdig.fastq"
printToLogFile

fastqc --quiet -f fastq NONFLASHED_REdig.fastq

################################################################
# Running Bowtie for the digested file..
printThis="Running Bowtie for the digested files.."
printToLogFile


printThis="Flashed reads Bowtie .."
printToLogFile

echo "Beginning bowtie run (outputting run command after completion) .."
setMparameter
bowtie -p 1 --chunkmb "${BOWTIEMEMORY}" ${otherBowtieParameters} ${bowtieQuals} ${mParameter} --best --strata --sam "${BowtieGenome}" FLASHED_REdig.fastq > FLASHED_REdig.sam
#bowtie -p 1 -m 2 --best --strata --sam --chunkmb 256 ${bowtieQuals} "${BowtieGenome}" Combined_reads_REdig.fastq Combined_reads_REdig.sam

testedFile="FLASHED_REdig.sam"
doTempFileTesting

doQuotaTesting

samtools view -SH FLASHED_REdig.sam | grep bowtie



printThis="Non-flashed reads Bowtie .."
printToLogFile

echo "Beginning bowtie run (outputting run command after completion) .."
setMparameter
bowtie -p 1 --chunkmb "${BOWTIEMEMORY}" ${otherBowtieParameters} ${bowtieQuals} ${mParameter} --best --strata --sam "${BowtieGenome}" NONFLASHED_REdig.fastq > NONFLASHED_REdig.sam
#bowtie -p 1 -m 2 --best --strata --sam --chunkmb 256 ${bowtieQuals} "${BowtieGenome}" Combined_reads_REdig.fastq Combined_reads_REdig.sam

testedFile="NONFLASHED_REdig.sam"
doTempFileTesting

doQuotaTesting

samtools view -SH NONFLASHED_REdig.sam | grep bowtie

# Cleaning up after ourselves ..

printThis="Cleaning up the run folder.."
printToLogFile

#ls -lht Combined_reads_REdig.bam
ls -lht FLASHED_REdig.sam
ls -lht NONFLASHED_REdig.sam
#rm -f  Combined_reads_REdig.fastq
rm -f FLASHED_REdig.fastq NONFLASHED_REdig.fastq

else
# This is the "ONLY_CC_ANALYSER" end fi - if testrun, skipped everything before this point :
# assuming existing output on the above mentioned files - all correctly formed except captureC output !
echo
echo "RE-RUN ! - running only capC analyser script, and filtering (assuming previous pipeline output in the run folder)"
echo

# Here deleting the existing - and failed - capturec analysis directory. not touching public files.

    rm -rf "../F2_redGraphs_${Sample}_${CCversion}"
    rm -rf "../F3_orangeGraphs_${Sample}_${CCversion}"
    rm -rf "../F4_blatPloidyFilteringLog_${Sample}_${CCversion}"
    rm -rf "../F5_greenGraphs_separate_${Sample}_${CCversion}"
    rm -rf "../F6_greenGraphs_combined_${Sample}_${CCversion}"
    rm -rf "../F7_summaryFigure_${Sample}_${CCversion}"
    
    rm -rf ../filteringLogFor_PREfiltered_${Sample}_${CCversion} ../RAW_${Sample}_${CCversion} ../PREfiltered_${Sample}_${CCversion} ../FILTERED_${Sample}_${CCversion} ../COMBINED_${Sample}_${CCversion}

    
# Remove the malformed public folder for a new try..
    rm -rf ${PublicPath}
    
# Restoring the input sam files..

# Run crash : we will have SAM instead of bam - if we don't check existence here, we will overwrite (due to funny glitch in samtools 1.1 )
if [ ! -s FLASHED_REdig.sam ]
then
    samtools view -h FLASHED_REdig.bam > TEMP.sam
    mv -f TEMP.sam FLASHED_REdig.sam
    if [ -s FLASHED_REdig.sam ]; then
        rm -f FLASHED_REdig.bam
    else
        echo "EXITING ! : Couldn't make FLASHED_REdig.sam from FLASHED_REdig.bam" >> "/dev/stderr"
        exit 1
    fi
fi

# Run crash : we will have SAM instead of bam - if we don't check existence here, we will overwrite (due to funny glitch in samtools 1.1 )
if [ ! -s NONFLASHED_REdig.sam ]
then
    samtools view -h NONFLASHED_REdig.bam > TEMP.sam
    mv -f TEMP.sam NONFLASHED_REdig.sam
    if [ -s NONFLASHED_REdig.sam ]; then
        rm -f NONFLASHED_REdig.bam
    else
        echo "EXITING ! : Couldn't make NONFLASHED_REdig.sam from NONFLASHED_REdig.bam" >> "/dev/stderr"
        exit 1
    fi
fi
    
fi

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

################################################################
# Store the pre-CCanalyser log files for metadata html

printThis="Store the pre-CCanalyser log files for metadata html.."
printToLogFile

copyPreCCanalyserLogFilesToPublic


dpnGenomeName=$( echo "${GenomeFasta}" | sed 's/.*\///' | sed 's/\..*//' )
# output file :
# ${GenomeFasta}_dpnII_coordinates.txt

fullPathDpnGenome=$(pwd)"/genome_dpnII_coordinates.txt"

cd ..

################################################################
# Running CAPTURE-C analyser for the aligned file..

printThis="##################################"
printToLogFile
printThis="Running CCanalyser without filtering - generating the RED graphs.."
printToLogFile
printThis="##################################"
printToLogFile

runDir=$( pwd )
dirForQuotaAsking=${runDir}
samDirForCCanalyser=${runDir}

publicPathForCCanalyser="${PublicPath}/RAW"
JamesUrlForCCanalyser="${JamesUrl}/RAW"

CCscriptname="${captureScript}.pl"


################################

printThis="Flashed reads.."
printToLogFile

sampleForCCanalyser="RAW_${Sample}"

samForCCanalyser="F1_beforeCCanalyser_${Sample}_${CCversion}/FLASHED_REdig.sam"
testedFile="${samForCCanalyser}"
doTempFileTesting

rm -f parameters_for_filtering.log

# For testing purposes..
# otherParameters="${otherParameters} --dump"

FLASHED=1
DUPLFILTER=0
runCCanalyser
doQuotaTesting

printThis="##################################"
printToLogFile

printThis="Non-flashed reads.."
printToLogFile

sampleForCCanalyser="RAW_${Sample}"

samForCCanalyser="F1_beforeCCanalyser_${Sample}_${CCversion}/NONFLASHED_REdig.sam"
testedFile="${samForCCanalyser}"
doTempFileTesting

rm -f parameters_for_filtering.log

# For testing purposes..
# otherParameters="${otherParameters} --dump"

FLASHED=0
DUPLFILTER=0
runCCanalyser
doQuotaTesting


else
# This is the "ONLY_HUB" end fi - if only hubbing, skipped everything before this point :
# assuming existing output on the above mentioned files - all correctly formed except the public folder (assumes correctly generated bigwigs, however) !
echo
echo "RE-HUB ! - running only public tracks.txt file update (assumes existing bigwig files and other hub structure)."
echo "If your bigwig files are missing (you see no .bw files in ${publicPathForCCanalyser}, or you wish to RE-LOCATE your data hub, run with --onlyCCanalyser parameter (instead of the --onlyHub parameter)"
echo "This is because parts of the hub generation are done inside captureC analyser script, and this assumes only tracks.txt generation failed."
echo

# Remove the malformed tracks.txt for a new try..
#rm -f ${publicPathForCCanalyser}/${sampleForCCanalyser}_${CCversion}_tracks.txt
rm -f ${PublicPath}/RAW/RAW_${Sample}_${CCversion}_tracks.txt
rm -f ${PublicPath}/FILTERED/FILTERED_${Sample}_${CCversion}_tracks.txt
rm -f ${PublicPath}/${Sample}_${CCversion}_tracks.txt

fi

################################################################
# Updating the public folder with analysis log files..

# to create file named ${Sample}_description.html - and link it to each of the tracks.

subfolder="RAW"
updateCCanalyserDataHub

mv -f RAW_${Sample}_${CCversion} F2_redGraphs_${Sample}_${CCversion}

#################################################################

# Running again - to make the otherwise filtered-but-not-blat-and-ploidy-filtered

printThis="##################################"
printToLogFile
printThis="Re-running CCanalyser with filtering - generating data to enter blat and ploidy filters.."
printToLogFile
printThis="##################################"
printToLogFile

runDir=$( pwd )
samDirForCCanalyser=${runDir}

publicPathForCCanalyser="${PublicPath}/PREfiltered"
JamesUrlForCCanalyser="${JamesUrl}/PREfiltered"

CCscriptname="${captureScript}.pl"

################################

printThis="Flashed reads.."
printToLogFile

sampleForCCanalyser="PREfiltered_${Sample}"

samForCCanalyser="F1_beforeCCanalyser_${Sample}_${CCversion}/FLASHED_REdig.sam"
testedFile="${samForCCanalyser}"
doTempFileTesting

rm -f parameters_for_filtering.log

FLASHED=1
DUPLFILTER=1
runCCanalyser
doQuotaTesting

# Adding the flashed filename, but not forgetting the common prefix either..
# cat parameters_for_filtering.log | grep dataprefix > prefixline
# sed -i 's/^dataprefix\s/dataprefix_FLASHED\t/' parameters_for_filtering.log
# cat parameters_for_filtering.log prefixline > FLASHED_parameters_for_filtering.log
# rm -f parameters_for_filtering.log

# Adding the flashed filename
mv -f parameters_for_filtering.log FLASHED_parameters_for_filtering.log
sed -i 's/^dataprefix\s/dataprefix_FLASHED\t/' FLASHED_parameters_for_filtering.log

printThis="##################################"
printToLogFile

printThis="Non-flashed reads.."
printToLogFile

sampleForCCanalyser="PREfiltered_${Sample}"

samForCCanalyser="F1_beforeCCanalyser_${Sample}_${CCversion}/NONFLASHED_REdig.sam"
testedFile="${samForCCanalyser}"
doTempFileTesting

rm -f parameters_for_filtering.log

FLASHED=0
DUPLFILTER=1
runCCanalyser
doQuotaTesting

# Adding the nonflashed filename
mv -f parameters_for_filtering.log NONFLASHED_parameters_for_filtering.log
sed -i 's/^dataprefix\s/dataprefix_NONFLASHED\t/' NONFLASHED_parameters_for_filtering.log

#################

# Combining parameter files..

cat FLASHED_parameters_for_filtering.log NONFLASHED_parameters_for_filtering.log | sort | uniq > parameters_for_filtering.log
rm -f FLASHED_parameters_for_filtering.log NONFLASHED_parameters_for_filtering.log


##################################
# Filtering the data..
printThis="##################################"
printToLogFile
printThis="Ploidy filtering and blat-filtering the data.."
printToLogFile
printThis="##################################"
printToLogFile

# ${CaptureFilterPath}
# /home/molhaem2/telenius/CC2/filter/VS101/filter.sh -p parameters.txt --outputToRunfolder --extend 30000
#
#        -p) parameterfile=$2 ; shift 2;;
#        --parameterfile) parameterfile=$2 ; shift 2;;
#        --noploidyfilter) ploidyfilter=0 ; shift 1;;
#        --pipelinecall) pipelinecall=1 ; shift 1;;
#        --extend) extend=$2 ; shift 2;;

echo "${CaptureFilterPath}/filter.sh -p parameters_for_filtering.log -s ${CaptureFilterPath} --pipelinecall ${ploidyFilter} --extend ${extend}"
echo "${CaptureFilterPath}/filter.sh -p parameters_for_filtering.log -s ${CaptureFilterPath} --pipelinecall ${ploidyFilter} --extend ${extend}"  >> "/dev/stderr"

#        --stepSize) stepSize=$2 ; shift 2;;
#        --tileSize) tileSize==$2 ; shift 2;;
#        --minScore) minScore=$2 ; shift 2;;
#        --maxIntron) maxIntron=$2 ; shift 2;;
#        --oneOff) oneOff=$2 ; shift 2;;

echo "--stepSize ${stepSize} --minScore ${minScore} --maxIntron=${maxIntron} --tileSize=${tileSize} --oneOff=${oneOff}"
echo "--stepSize ${stepSize} --minScore ${minScore} --maxIntron=${maxIntron} --tileSize=${tileSize} --oneOff=${oneOff}" >> "/dev/stderr"

echo "--reuseBLAT ${reuseBLATpath}"
echo "--reuseBLAT ${reuseBLATpath}" >> "/dev/stderr"

mkdir filteringLogFor_${sampleForCCanalyser}_${CCversion}
mv parameters_for_filtering.log filteringLogFor_${sampleForCCanalyser}_${CCversion}/.
cd filteringLogFor_${sampleForCCanalyser}_${CCversion}

TEMPreturnvalue=0
TEMPreturnvalue=$( ${CaptureFilterPath}/filter.sh --reuseBLAT ${reuseBLATpath} -p parameters_for_filtering.log --pipelinecall ${ploidyFilter} --extend ${extend} --stepSize ${stepSize} --minScore ${minScore} --maxIntron=${maxIntron} --tileSize=${tileSize} --oneOff=${oneOff} > filtering.log )
cat filtering.log
rm -f ${publicPathForCCanalyser}/filtering.log
cp filtering.log ${publicPathForCCanalyser}/.

if [ "${TEMPreturnvalue}" -ne 0 ]; then
    
    printThis="Filtering after BLAT was crashed !"
    printToLogFile
    printThis="( If this was BLAT-generation run, this is what you wanted ) "
    printToLogFile
    
    printThis="Your psl-files for BLAT-filtering can be found in folder :\n $( pwd )/BlatPloidyFilterRun/REUSE_blat/"
    printToLogFile

    printThis="EXITING !"
    printToLogFile
    
    exit 1
    
fi


cd ..

# By default the output of this will go to :
# ${Sample}_${CCversion}/BLAT_PLOIDY_FILTERED_OUTPUT
# because the parameter file line for data location is
# ${Sample}_${CCversion}

################################################################
# Updating the public folder with analysis PREfiltered log files..

# to create file named ${Sample}_description.html - and link it to each of the tracks.

subfolder="PREfiltered"
updateCCanalyserDataHub

mv -f PREfiltered_${Sample}_${CCversion} F3_orangeGraphs_${Sample}_${CCversion}

################################################################


printThis="##################################"
printToLogFile
printThis="Re-running CCanalyser for the filtered data.."
printToLogFile
printThis="##################################"
printToLogFile

runDir=$( pwd )
samDirForCCanalyser="${runDir}"

publicPathForCCanalyser="${PublicPath}/FILTERED"
JamesUrlForCCanalyser="${JamesUrl}/FILTERED"

CCscriptname="${captureScript}.pl"

PREVsampleForCCanalyser="${sampleForCCanalyser}"

# FLASHED

printThis="------------------------------"
printToLogFile
printThis="FLASHED file.."
printToLogFile

# keeping the "RAW" in the file name - as this part (input folder location) still needs that
ln -s filteringLogFor_${PREVsampleForCCanalyser}_${CCversion}/BlatPloidyFilterRun/BLAT_PLOIDY_FILTERED_OUTPUT/FLASHED_REdig_${CCversion}_filtered_combined.sam FLASHED_REdig.sam
samForCCanalyser="FLASHED_REdig.sam"

FILTEREDsamBasename=$( echo ${samForCCanalyser} | sed 's/.*\///' | sed 's/\.sam$//' )
testedFile="${samForCCanalyser}"
doTempFileTesting

# Now changing the identifier from "RAW" to "FILTERED" - to set the output folder

sampleForCCanalyser="FILTERED_${Sample}"

FLASHED=1
DUPLFILTER=0
runCCanalyser
doQuotaTesting

# Remove symlink
rm -f FLASHED_REdig.sam

# NONFLASHED

printThis="------------------------------"
printToLogFile
printThis="NONFLASHED file.."
printToLogFile


# keeping the "RAW" in the file name - as this part (input folder location) still needs that
ln -s filteringLogFor_${PREVsampleForCCanalyser}_${CCversion}/BlatPloidyFilterRun/BLAT_PLOIDY_FILTERED_OUTPUT/NONFLASHED_REdig_${CCversion}_filtered_combined.sam NONFLASHED_REdig.sam
samForCCanalyser="NONFLASHED_REdig.sam"

FILTEREDsamBasename=$( echo ${samForCCanalyser} | sed 's/.*\///' | sed 's/\.sam$//' )
testedFile="${samForCCanalyser}"
doTempFileTesting

# Now changing the identifier from "RAW" to "FILTERED" - to set the output folder
sampleForCCanalyser="FILTERED_${Sample}"

FLASHED=0
DUPLFILTER=0
runCCanalyser
doQuotaTesting

# Remove symlink
rm -f NONFLASHED_REdig.sam

################################################################
# Updating the public folder with analysis log files..

# to create file named ${Sample}_description.html - and link it to each of the tracks.

subfolder="FILTERED"
updateCCanalyserDataHub

mv -f FILTERED_${Sample}_${CCversion} F5_greenGraphs_separate_${Sample}_${CCversion}

################################################################

printThis="##################################"
printToLogFile
printThis="Combining FLASHED and NONFLASHED CCanalyser filtered data .."
printToLogFile
printThis="##################################"
printToLogFile

printThis="Combining sam files.."
printToLogFile

cat filteringLogFor_PREfiltered_${Sample}_${CCversion}/BlatPloidyFilterRun/BLAT_PLOIDY_FILTERED_OUTPUT/NONFLASHED_REdig_${CCversion}_filtered_combined.sam | grep -v "^@" | \
cat filteringLogFor_PREfiltered_${Sample}_${CCversion}/BlatPloidyFilterRun/BLAT_PLOIDY_FILTERED_OUTPUT/FLASHED_REdig_${CCversion}_filtered_combined.sam - > COMBINED.sam

COMBINEDsamBasename=$( echo ${samForCCanalyser} | sed 's/.*\///' | sed 's/\.sam$//' )
samForCCanalyser="COMBINED.sam"
COMBINEDsamBasename=$( echo ${samForCCanalyser} | sed 's/.*\///' | sed 's/\.sam$//' )
testedFile="${samForCCanalyser}"
doTempFileTesting

printThis="------------------------------"
printToLogFile
printThis="Running CCanalyser.."
printToLogFile


runDir=$( pwd )
samDirForCCanalyser="${runDir}"

publicPathForCCanalyser="${PublicPath}/COMBINED"
JamesUrlForCCanalyser="${JamesUrl}/COMBINED"

CCscriptname="${captureScript}.pl"

sampleForCCanalyser="COMBINED_${Sample}"

# This means : flashing is "NOT IN USE" - and marks the output tracks with name "" instead of "FLASHED" or "NONFLASHED"
FLASHED=-1
DUPLFILTER=0
runCCanalyser
doQuotaTesting

# Remove input file
rm -f COMBINED.sam

################################################################
# Updating the public folder with analysis log files..

# to create file named ${Sample}_description.html - and link it to each of the tracks.

subfolder="COMBINED"
updateCCanalyserDataHub

mv -f COMBINED_${Sample}_${CCversion} F6_greenGraphs_combined_${Sample}_${CCversion}

################################################################


if [[ ${saveDpnGenome} -eq "0" ]] ; then
  rm -f "genome_${REenzyme}_coordinates.txt"  
fi

# Generating combined data hub

sampleForCCanalyser="${Sample}"
publicPathForCCanalyser="${PublicPath}"
JamesUrlForCCanalyser="${JamesUrl}"

generateCombinedDataHub

# Cleaning up after ourselves ..

cleanUpRunFolder
# makeSymbolicLinks

# Data hub address (print to stdout) ..
updateHub_part3final

echo
echo "All done !"
echo  >> "/dev/stderr"
echo "All done !" >> "/dev/stderr"

# Copying log files

echo "Copying run log files.." >> "/dev/stderr"

cp -f ./qsub.out "${PublicPath}/${Sample}_logFiles/${Sample}_qsub.out"
cp -f ./qsub.err "${PublicPath}/${Sample}_logFiles/${Sample}_qsub.err"

echo "Log files copied !" >> "/dev/stderr"

exit 0