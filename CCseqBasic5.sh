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
    if [ "${thisWashelpRequest}" != 1 ];then
    echo
    echo "Analysis complete !"
    date
    fi
fi
}
trap finish EXIT

#------------------------------------------

# script top path setup , usage script load.

CCversion="CS5"
captureScript="analyseMappedReads"
CCseqBasicVersion="CCseqBasic5"

CaptureTopPath="$( which $0 | sed 's/\/'${CCseqBasicVersion}'.sh$//' )"
CapturePipePath="${CaptureTopPath}/bin/subroutines"

. ${CapturePipePath}/usageAndVersion.sh

#------------------------------------------

# help user cases .

thisWashelpRequest=0
if [ "$1" == '-h' ] || [ "$1" == '--help' ]
then
    thisWashelpRequest=1
    usage ;
fi


#------------------------------------------

QSUBOUTFILE="qsub.out"
QSUBERRFILE="qsub.err"

CapturesiteFile=""
TRIM=1
GENOME=""
WINDOW=200
INCREMENT=20
CAPITAL_M=0
LOWERCASE_M=0
LOWERCASE_V=-1
BOWTIEMEMORY="256"
Sample="sample"
Read1=""
Read2=""

LANES=1
GZIP=0
SINGLE_END=0

CUSTOMAD=-1
ADA31="no"
ADA32="no"

# trimgalore default
QMIN=20

# bowtie default
BOWTIE=1

# flash defaults
flashOverlap=10
flashErrorTolerance=0.25

saveDpnGenome=0

ucscBuild=""
REDGREEN=0

otherBowtie1Parameters=""
otherBowtie2Parameters=""
bowtie1MismatchBehavior=""
bowtie2MismatchBehavior=""

otherParameters=""
PublicPath="UNDETERMINED"

ploidyFilter=""
extend=20000

ampliconSize=300

# If we have many capture-site (REfragment)s, the stuff can be eased up by analysing only in cis.
onlyCis=0

# Blat flags
stepSize=5 # Jon default - James used blat default, which is "tileSize", in this case thus 11 (this was the setting in CC2 and CC3 - i.e filter versions VS101 and VS102)
tileSize=11 # Jon, James default
minScore=10 # Jon new default. Jon default before2016 and CC4 default until 080916 minScore=30 - James used minScore=30 (this was the setting in CC2 and CC3 - i.e filter versions VS101 and VS102)
minIdentity=70 # Jon, James default.
minMatch=2 # blat default
maxIntron=4000 # blat default 750000- James used maxIntron=4000 (this was the setting in CC2 and CC3 - i.e filter versions VS101 and VS102)
oneOff=0 # oneOff=1 would allow 1 mismatch in tile (blat default = 0 - that is also CC3 and CC2 default)

# Whether we reuse blat results from earlier run ..
# Having this as "." will search from the run dir when blat is ran - so file will not be found, and thus BLAT will be ran normally.
reuseBLATpath="."
ONLY_BLAT=0

REenzyme="dpnII"

# Skip other stages - assume input from this run has been ran earlier - to construct to THIS SAME FOLDER everything else
# but as the captureC analyser naturally crashed - this will jump right to the beginning of that part..
ONLY_CC_ANALYSER=0
# Rerun public folder generation and filling. Will not delete existing folder, but will overwrite all files (and start tracks.txt from scratch).
ONLY_HUB=0

strandSpecificDuplicates=0

srrFastq=0

#------------------------------------------

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

# INPUT the fastq files
. ${CapturePipePath}/inputFastqs.sh

# SETTING THE GENOME BUILD PARAMETERS
. ${CapturePipePath}/genomeSetters.sh

# SETTING THE BLACKLIST GENOME LIST PARAMETERS
. ${CapturePipePath}/blacklistSetters.sh

# SORTING HELPER SUBROUTINES
. ${CapturePipePath}/sort_helpers.sh

# PRINTING HELP AND VERSION MESSAGES
# . ${CapturePipePath}/usageAndVersion.sh
# (loaded already above - where the help user case is done)

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

# Calling in the CONFIGURATION script and its default setup :

# Defaulting this to "not in use" - if it is not set in the config file.
CaptureDigestPath="NOT_IN_USE"

#------------------------------------------

# Calling in the CONFIGURATION script and its default setup :

echo "Calling in the conf/config.sh script and its default setup .."

CaptureDigestPath="NOT_IN_USE"
supportedGenomes=()
BOWTIE1=()
BOWTIE2=()
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

OPTS=`getopt -o h,m:,M:,o:,c:,s:,w:,i:,v: --long help,dump,snp,dpn,nla,hind,gz,strandSpecificDuplicates,redGreen,onlyCis,onlyBlat,UMI,useSymbolicLinks,SRR,CCversion:,BLATforREUSEfolderPath:,globin:,outfile:,errfile:,lanes:,limit:,pf:,genome:,R1:,R2:,saveGenomeDigest,dontSaveGenomeDigest,trim,noTrim,chunkmb:,bowtie1,bowtie2,window:,increment:,ada3read1:,ada3read2:,extend:,onlyCCanalyser,onlyHub,noPloidyFilter:,qmin:,flashBases:,flashMismatch:,stringent,trim3:,trim5:,seedmms:,seedlen:,maqerr:,stepSize:,tileSize:,minScore:,minIdentity:,minMatch:,maxIntron:,oneOff:,wobblyEndBinWidth:,ampliconSize:,sonicationSize: -- "$@"`
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
        -o) CapturesiteFile=$2 ; shift 2;;
        -c) CapturesiteFile=$2 ; shift 2;;
        -w) WINDOW=$2 ; shift 2;;
        -i) INCREMENT=$2 ; shift 2;;
        -s) Sample=$2 ; shift 2;;
        -v) LOWERCASE_V=$2; shift 2;;
        --help) usage ; shift;;
        --UMI) printThis="UMI flag temporarily broken 01Nov2018\nEXITING";printToLogFile;exit 1;otherParameters="$otherParameters --umi" ; shift;;
        --useSymbolicLinks) otherParameters="$otherParameters --symlinks" ; shift;;
        --CCversion) CCversion="$2"; shift 2;;       
        --dpn) REenzyme="dpnII" ; shift;;
        --nla) REenzyme="nlaIII" ; shift;;
        --hind) REenzyme="hindIII" ; shift;;
        --onlyCCanalyser) ONLY_CC_ANALYSER=1 ; shift;;
        --onlyHub) ONLY_HUB=1 ; shift;;
        --onlyCis) onlyCis=1;otherParameters="$otherParameters --onlycis"; shift;;
        --onlyBlat) ONLY_BLAT=1 ; shift;;
        --R1) Read1=$2 ; shift 2;;
        --R2) Read2=$2 ; shift 2;;
        --lanes) LANES=$2 ; shift 2;;
        --gz) GZIP=1 ; shift;;
        --SRR) srrFastq=1 ; shift;;
        --bowtie1) BOWTIE=1 ; shift;;
        --bowtie2) BOWTIE=2 ; shift;;
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
        --ampliconSize) ampliconSize=$2 ; shift 2;;
        --sonicationSize) ampliconSize=$2 ; shift 2;;
        --strandSpecificDuplicates) otherParameters="$otherParameters --stranded"; strandSpecificDuplicates=1 ; shift;;
        --dump) otherParameters="$otherParameters --dump" ; shift;;
        --snp) otherParameters="$otherParameters --snp" ; shift;;
        --globin) otherParameters="$otherParameters --globin $2" ; shift 2;;
        --limit) otherParameters="$otherParameters --limit $2" ; shift 2;;
        --stringent) otherParameters="$otherParameters --stringent" ; shift 1;;
        --pf) PublicPath="$2" ; shift 2;;
        --redGreen) REDGREEN=1 ; shift;;
        --qmin) QMIN="$2" ; shift 2;;
        --BLATforREUSEfolderPath) reuseBLATpath="$2" ; shift 2;;
        --flashBases) flashOverlap="$2" ; shift 2;;
        --flashMismatch) flashErrorTolerance="$2" ; shift 2;;
        --trim3) otherBowtieParameters="${otherBowtieParameters} --trim3 $2 " ; shift 2;;
        --trim5) otherBowtieParameters="${otherBowtieParameters} --trim5 $2 " ; shift 2;;
        --seedmms) bowtie1MismatchBehavior="${bowtie1MismatchBehavior} --seedmms $2 " ; ${bowtie2MismatchBehavior}="${bowtie2MismatchBehavior} -N $2 "  ; shift 2;;
        --seedlen) bowtie1MismatchBehavior="${bowtie1MismatchBehavior} --seedlen $2 " ; ${bowtie2MismatchBehavior}="${bowtie2MismatchBehavior} -L $2 " ; shift 2;;
        --maqerr) bowtie1MismatchBehavior="${bowtie1MismatchBehavior} --maqerr $2 " ; shift 2;;
        --stepSize) stepSize=$2 ; shift 2;;
        --tileSize) tileSize=$2 ; shift 2;;
        --minScore) minScore=$2 ; shift 2;;
        --minMatch) minMatch=$2 ; shift 2;;
        --minIdentity) minIdentity=$2 ; shift 2;;
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
# Here parsing the parameter files - if they are not purely tab-limited, but partially space-limited, or multiple-tab limited, this fixes it.

echo
echo "PARAMETER FILES GIVEN IN RUN FOLDER (if any) :"
echo

if [ $(ls PIPE*.txt 2>/dev/null| grep -c "") -gt 0 ]; then

for file in ./PIPE*.txt
    do
        echo ${file}
        sed -i 's/\s\s*/\t/g' ${file}
    done
    
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
echo "------------------------------" >> parameters_capc.log
echo "GZIP ${GZIP} (TRUE=1, FALSE=0) - if fastq input files are gzipped "
echo "Read1 ${Read1}" >> parameters_capc.log
echo "Read2 ${Read2}" >> parameters_capc.log
echo "srrFastq ${srrFastq} (TRUE=1, FALSE=0) - if fastq files are SRR format (instead of Illumina) " >> parameters_capc.log
echo "------------------------------" >> parameters_capc.log
echo "GENOME ${GENOME}" >> parameters_capc.log
echo "GenomeIndex ${GenomeIndex}" >> parameters_capc.log
echo "CapturesiteFile ${CapturesiteFile}" >> parameters_capc.log
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
echo "ampliconSize ${ampliconSize}"  >> parameters_capc.log
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


testedFile="${CapturesiteFile}"
doInputFileTesting

#---------------------------------------------------------

# Doing the ONLY_BLAT  user case first - they doesn't need existing input files (except the capture-site (REfragment) file) - so we shouldn't enter any testing of parameters here.

if [[ "${ONLY_BLAT}" -eq "1" ]]; then
{

  paramGenerationRunFineOK=0
  
  printThis="Running ONLY BLATS (user given --onlyBlat flag, or parallel run first step)"
  printToLogFile

  # --------------------------
  
  # RE enzyme digestion (if needed .. )

  dpnGenomeName=""
  fullPathDpnGenome=""
  generateReDigest

  CCscriptname="${captureScript}.pl"
  runCCanalyserOnlyBlat


  # Return information to log file if we are parallel ..
  if [ "${paramGenerationRunFineOK}" -ne 0 ];then {

    printThis="CCanalyser to prepare BLAT runs failed."
    printToLogFile
    printThis="EXITING !"
    printToLogFile
  
  exit 1
  }
  fi
  
  # --------------------------

  ${CaptureFilterPath}/filter.sh --onlyBlat ${ONLY_BLAT} --reuseBLAT ${reuseBLATpath} -p parameters_for_filtering.log --pipelinecall --extend ${extend} --onlyCis ${onlyCis} --stepSize ${stepSize} --minScore ${minScore} --minIdentity=${minIdentity} --minMatch=${minMatch} --maxIntron=${maxIntron} --tileSize=${tileSize} --oneOff=${oneOff} --bowtieMemory ${BOWTIEMEMORY} > filtering.log
  # cat filtering.log

  if [ "$?" -ne 0 ]; then {
      printThis="Running filter.sh crashed - BLAT filtering failed !"
      printToLogFile
      printThis="EXITING !"
      printToLogFile    
      exit 1
  }  
  fi

  # --------------------------
  
  rm -rf blat_run_params
  mkdir blat_run_params
  mv blatParams.txt blat_run_params/.
  mv -f parameters_*.log blat_run_params/.
    
  echo  > How_to_use_these_BLAT_files.txt   
  echo "Use the generated BLAT filtering .psl files by adding this to your run command : " >> How_to_use_these_BLAT_files.txt
  echo  >> How_to_use_these_BLAT_files.txt 
  echo '--BLATforREUSEfolderPath '$( pwd )/BlatPloidyFilterRun/REUSE_blat/ >> How_to_use_these_BLAT_files.txt
  
  echo  > How_to_use_these_BLAT_files.txt
  echo "Your psl-files for BLAT-filtering can be found in folder :\n $( pwd )/BlatPloidyFilterRun/REUSE_blat/" >> How_to_use_these_BLAT_files.txt
  echo  >> How_to_use_these_BLAT_files.txt
  echo "Use the generated BLAT filtering .psl files by adding this to your run command : " >> How_to_use_these_BLAT_files.txt
  echo  >> How_to_use_these_BLAT_files.txt 
  echo '--BLATforREUSEfolderPath '$( pwd )/BlatPloidyFilterRun/REUSE_blat/ >> How_to_use_these_BLAT_files.txt
  echo  >> How_to_use_these_BLAT_files.txt
  echo "Here full list of generated files : " >> How_to_use_these_BLAT_files.txt
  echo >> How_to_use_these_BLAT_files.txt
  echo "ls -lht $( pwd )/BlatPloidyFilterRun/REUSE_blat/" >> How_to_use_these_BLAT_files.txt
  ls -lht $( pwd )/BlatPloidyFilterRun/REUSE_blat/ >> How_to_use_these_BLAT_files.txt
  echo >> How_to_use_these_BLAT_files.txt
  
  printThis="Your psl-files for BLAT-filtering can be found in folder :\n $( pwd )/BlatPloidyFilterRun/REUSE_blat/"
  printToLogFile

  echo "Use the generated BLAT filtering .psl files in CCseqBasic by adding this to your run command : "
  echo '--BLATforREUSEfolderPath '$( pwd )/BlatPloidyFilterRun/REUSE_blat/ 
  
  printThis="Details of this in the ouput file 'How_to_use_these_BLAT_files.txt' "
  printToLogFile

  
  echo
  echo "All done !"
  echo  >> "/dev/stderr"
  echo "All done !" >> "/dev/stderr"
  
  exit 0
  
}
fi


# ---------------------------------------

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

#--------Test if fastq paths given correctly ------------------------------------------------------

if [ -r "./PIPE_fastqPaths.txt" ] && [ "${Read1}" != "" ]; then
    printThis="PIPE_fastqPaths.txt and --R1 parameter given at the same time. Only one at a time is allowed !"
    printToLogFile
    printThis="EXITING ! "
    printToLogFile
    exit 1
fi

if [ -r "./PIPE_fastqPaths.txt" ] && [ "${Read2}" != "" ]; then
    printThis="PIPE_fastqPaths.txt and --R2 parameter given at the same time. Only one at a time is allowed !"
    printToLogFile
    printThis="EXITING ! "
    printToLogFile
    exit 1
fi
    
#--------THE-LOOP-over-all-FASTQ-file-based-data-sets------------------------------------------------------

if [ -r "./PIPE_fastqPaths.txt" ] ; then

    # Check how many columns we have.
    test=0
    test=$( cut -f 3 ./PIPE_fastqPaths.txt | grep -vc "^\s*$" )

    # If we have 2 columns paired end :
    if [ "${test}" -eq "0" ]; then

    # Fake list - only first element gets filled (as this is multi-lane support, no multi-sample support)
    fileList1=($( cat ./PIPE_fastqPaths.txt | grep -v '^\s*$' | cut -f 1 | tr '\n' ',' | sed 's/,$//' ))
    fileList2=($( cat ./PIPE_fastqPaths.txt | grep -v '^\s*$' | cut -f 2 | grep -v '^\s*$' | tr '\n' ',' | sed 's/,$//' ))
    
    # If we have 3 columns paired end :
    else
    
    cat ./PIPE_fastqPaths.txt | grep -v '^\s*$' | cut -f 1,3 | awk '{ print $2"/"$1 }'| sed 's/\/\//\//' > forRead1.txt
    cat ./PIPE_fastqPaths.txt | grep -v '^\s*$' | cut -f 2,3 | awk '{ print $2"/"$1 }'| sed 's/\/\//\//' > forRead2.txt
    
    # Fake list - only first element gets filled (as this is multi-lane support, no multi-sample support)
    fileList1=($( cat ./forRead1.txt | tr '\n' ',' | sed 's/,$//' ))
    fileList2=($( cat ./forRead2.txt | tr '\n' ',' | sed 's/,$//' ))
    
    rm -f forRead1.txt forRead2.txt
    
    fi
    
    LANES=$(($( cat ./PIPE_fastqPaths.txt | grep -v '^\s*$' | grep -c "" )))
    
else

#---------------------------------------
    
# If we didn't have PIPE_fastqPaths.txt , we have Read1 and Read2 as input parameters instead.

# Copy files over..

    # Fake list - only first element gets filled (as this is multi-lane support, no multi-sample support)
    fileList1=(${Read1})
    fileList2=(${Read2})
    
    LANES=1

#---------------------------------------

fi

# Fetch the files  ..

printRunStartArraysFastq

echo "LANES ${LANES}" >> parameters_capc.log
echo "LANES ${LANES}"

# This is a fake list - only having one element (see above)
for (( i=0; i<=$(( ${#fileList1[@]} -1 )); i++ ))
do
    
    # If we have single lane sequencing.
    if [ "$LANES" -eq 1 ] ; then 
    
    #Fetch FASTQ :
    fetchFastq
    
    else

    # If we have MULTIPLE lanes from sequencing.
    
    fetchFastqMultilane

    fi
    
done

#---------------------------------------

# Check that we have the files ..

testedFile="READ1.fastq"
doInputFileTesting
testedFile="READ2.fastq"
doInputFileTesting

mv -f READ1.fastq F1_beforeCCanalyser_${Sample}_${CCversion}/READ1.fastq
mv -f READ2.fastq F1_beforeCCanalyser_${Sample}_${CCversion}/READ2.fastq


echo "Generated fastqs :"
ls -lh F1_beforeCCanalyser_${Sample}_${CCversion} | cut -d " " -f 1,2,3,4 --complement
echo "In folder :"
pwd 

testedFile="F1_beforeCCanalyser_${Sample}_${CCversion}/READ1.fastq"
doTempFileTesting
testedFile="F1_beforeCCanalyser_${Sample}_${CCversion}/READ2.fastq"
doTempFileTesting

#---------------------------------------

# Save capture-site (REfragment) file full path (to not to lose the file when we cd into the folder, if we used relative paths ! )
TEMPdoWeStartWithSlash=$(($( echo ${CapturesiteFile} | awk '{print substr($1,1,1)}' | grep -c '/' )))
if [ "${TEMPdoWeStartWithSlash}" -eq 0 ]
then
 CapturesiteFile=$(pwd)"/"${CapturesiteFile}
fi

testedFile="${CapturesiteFile}"
doInputFileTesting

fi
fi

# Go into output folder..
cd F1_beforeCCanalyser_${Sample}_${CCversion}

if [[ ${ONLY_HUB} -eq "0" ]]; then
if [[ ${ONLY_CC_ANALYSER} -eq "0" ]]; then
    
#---------------------------------------

# SRR to Illumina - if we have --SRR flag  ..

if [ "${srrFastq}" -eq 1 ]; then

printThis="Transforming read names from SRR format to Illumina format .."
printToLogFile

printThis="${RunScriptsPath}/srr_to_illumina.pl READ1.fastq 1 "
printToLogFile

${RunScriptsPath}/srr_to_illumina.pl READ1.fastq 1

ls -lh READ1*fastq | cut -d " " -f 1,2,3,4 --complement
mv -f READ1_asIllumina.fastq READ1.fastq

printThis="${RunScriptsPath}/srr_to_illumina.pl READ1.fastq 2 "
printToLogFile

${RunScriptsPath}/srr_to_illumina.pl READ2.fastq 2

ls -lh READ2*fastq | cut -d " " -f 1,2,3,4 --complement
mv -f READ2_asIllumina.fastq READ2.fastq

testedFile="READ1.fastq"
doInputFileTesting
testedFile="READ2.fastq"
doInputFileTesting

fi

#---------------------------------------


################################################################
#Check BOWTIE quality scores..

printThis="Checking the quality score scheme of the fastq files.."
printToLogFile
    
    bowtieQuals=""
    LineCount=$(($( grep -c "" READ1.fastq )/4))
    if [ "${LineCount}" -gt 100000 ] ; then
        bowtieQuals=$( perl ${RunScriptsPath}/fastq_scores_bowtie${BOWTIE}.pl -i READ1.fastq -r 90000 )
    else
        rounds=$((${LineCount}-10))
        bowtieQuals=$( perl ${RunScriptsPath}/fastq_scores_bowtie${BOWTIE}.pl -i READ1.fastq -r ${rounds} )
    fi
    
    echo "Flash, Trim_galore and Bowtie will be ran in quality score scheme : ${bowtieQuals}"

    # The location of "zero" for the filtering/trimming programs cutadapt, trim_galore, flash    
    intQuals=""
    if [ "${bowtieQuals}" == "--phred33-quals" ] || [ "${bowtieQuals}" == "--phred33" ]; then
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
if [ "$?" -ne 0 ]; then
printThis="FastQC run failed ! Possible reasons : \n 1) did you maybe use .gz packed files without adding --gz to the run parameters ? \n 2) did you try to run with corrupted input fastq files ? \n EXITING !! "
printToLogFile
exit 1
fi

    # Changing names of fastqc folders to be "ORIGINAL"
    
    rm -rf READ1_fastqc_ORIGINAL
    rm -rf READ2_fastqc_ORIGINAL
    
    mkdir READ1_fastqc_ORIGINAL
    mkdir READ2_fastqc_ORIGINAL
    
    mv -f READ1_fastqc.html READ1_fastqc_ORIGINAL/fastqc_report.html
    mv -f READ2_fastqc.html READ2_fastqc_ORIGINAL/fastqc_report.html 
    mv -f READ1_fastqc.zip  READ1_fastqc_ORIGINAL.zip
    mv -f READ2_fastqc.zip  READ2_fastqc_ORIGINAL.zip
   
    ls -lht

################################################################
# Trimgalore for the reads..

if [[ ${TRIM} -eq "1" ]]; then

printThis="Running trim_galore for the reads.."
printToLogFile

printThis="${RunScriptsPath}/QC_and_Trimming.sh -q ${intQuals} --filter 3 --qmin ${QMIN}"
printToLogFile

${RunScriptsPath}/QC_and_Trimming.sh -q "${intQuals}" --filter 3 --qmin ${QMIN}
if [ "$?" -ne 0 ]; then
printThis="TrimGalore run failed ! Possible reasons : \n 1) did you maybe use .gz packed files without adding --gz to the run parameters ? \n 2) did you try to run with corrupted input fastq files ? \n EXITING !! "
printToLogFile
exit 1
fi

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
    
    rm -rf READ1_fastqc_TRIMMED
    rm -rf READ2_fastqc_TRIMMED
    
    mkdir READ1_fastqc_TRIMMED
    mkdir READ2_fastqc_TRIMMED
    
    mv -f READ1_fastqc.html READ1_fastqc_TRIMMED/fastqc_report.html
    mv -f READ2_fastqc.html READ2_fastqc_TRIMMED/fastqc_report.html 
    
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

rm -rf FLASHED_fastqc
rm -rf NONFLASHED_fastqc
    
mkdir FLASHED_fastqc
mkdir NONFLASHED_fastqc
    
printThis="fastqc --quiet -f fastq FLASHED.fastq"
printToLogFile

fastqc --quiet -f fastq FLASHED.fastq
mv -f FLASHED_fastqc.html FLASHED_fastqc/fastqc_report.html


printThis="fastqc --quiet -f fastq NONFLASHED.fastq"
printToLogFile

fastqc --quiet -f fastq NONFLASHED.fastq
mv -f NONFLASHED_fastqc.html NONFLASHED_fastqc/fastqc_report.html


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

rm -rf FLASHED_REdig_fastqc
rm -rf NONFLASHED_REdig_fastqc
    
mkdir FLASHED_REdig_fastqc
mkdir NONFLASHED_REdig_fastqc
    
printThis="fastqc --quiet -f fastq FLASHED_REdig.fastq"
printToLogFile

fastqc --quiet -f fastq FLASHED_REdig.fastq
if [ "$?" -ne 0 ]; then
printThis="FastqQC run failed ! Possible reasons : \n 1) did you maybe use SRR archive format fastq files without adding --SRR to the run parameters ? \n 2) Fastq files not in Illumina format ? (also here you can rescue with --SRR if your rolling read ID is the first field in the fastq @name line) \n EXITING !! "
printToLogFile
exit 1
fi

mv -f FLASHED_REdig_fastqc.html FLASHED_REdig_fastqc/fastqc_report.html


printThis="fastqc --quiet -f fastq NONFLASHED_REdig.fastq"
printToLogFile

fastqc --quiet -f fastq NONFLASHED_REdig.fastq
mv -f NONFLASHED_REdig_fastqc.html NONFLASHED_REdig_fastqc/fastqc_report.html

################################################################
# Running Bowtie for the digested file..
printThis="Running Bowtie for the digested files.."
printToLogFile


printThis="Flashed reads Bowtie .."
printToLogFile

echo "Beginning bowtie run (outputting run command after completion) .."
setMparameter

if [ "${BOWTIE}" -eq 2 ] ; then
bowtie2 -p 1 ${otherBowtie2Parameters} ${bowtieQuals} -x ${BowtieGenome} -U FLASHED_REdig.fastq > FLASHED_REdig_unfiltered.sam
echo "bowtie2 -p 1 ${otherBowtie2Parameters} ${bowtieQuals} -x ${BowtieGenome} -U FLASHED_REdig.fastq"
else
bowtie -p 1 --chunkmb "${BOWTIEMEMORY}" ${otherBowtie1Parameters} ${bowtieQuals} ${mParameter} --best --strata --sam "${BowtieGenome}" FLASHED_REdig.fastq > FLASHED_REdig_unfiltered.sam
fi

#bowtie -p 1 -m 2 --best --strata --sam --chunkmb 256 ${bowtieQuals} "${BowtieGenome}" Combined_reads_REdig.fastq Combined_reads_REdig.sam

testedFile="FLASHED_REdig_unfiltered.sam"
doTempFileTesting

doQuotaTesting

samtools view -SH FLASHED_REdig_unfiltered.sam | grep bowtie

echo
echo "Read count - in bowtie output sam file : "
echo
flashstatus="FLASHED"
echo ${flashstatus}_REdig_unfiltered.sam
cat  ${flashstatus}_REdig_unfiltered.sam | grep -cv '^@'
echo

printThis="Sam to bam transform .."
printToLogFile

flashstatus="FLASHED"
rm -f ${flashstatus}FLASHED_REdig.fastq
echo "samtools view -hb ${flashstatus}_REdig_unfiltered.sam > ${flashstatus}_REdig_unfiltered.bam"
samtools view -hb ${flashstatus}_REdig_unfiltered.sam > ${flashstatus}_REdig_unfiltered.bam
ls -lht ${flashstatus}_REdig_unfiltered.bam
rm -f ${flashstatus}_REdig_unfiltered.sam

printThis="Non-flashed reads Bowtie .."
printToLogFile

echo "Beginning bowtie run (outputting run command after completion) .."
setMparameter
if [ "${BOWTIE}" -eq 2 ] ; then
bowtie2 -p 1 ${otherBowtie2Parameters} ${bowtieQuals} -x ${BowtieGenome} -U NONFLASHED_REdig.fastq > NONFLASHED_REdig_unfiltered.sam
echo "bowtie2 -p 1 ${otherBowtie2Parameters} ${bowtieQuals} -x ${BowtieGenome} -U NONFLASHED_REdig.fastq"
else
bowtie -p 1 --chunkmb "${BOWTIEMEMORY}" ${otherBowtie1Parameters} ${bowtieQuals} ${mParameter} --best --strata --sam "${BowtieGenome}" NONFLASHED_REdig.fastq > NONFLASHED_REdig_unfiltered.sam
fi
#bowtie -p 1 -m 2 --best --strata --sam --chunkmb 256 ${bowtieQuals} "${BowtieGenome}" Combined_reads_REdig.fastq Combined_reads_REdig.sam

testedFile="NONFLASHED_REdig_unfiltered.sam"
doTempFileTesting

doQuotaTesting

samtools view -SH NONFLASHED_REdig_unfiltered.sam | grep bowtie

echo
echo "Read count - in bowtie output sam file : "
echo
flashstatus="NONFLASHED"
echo ${flashstatus}_REdig_unfiltered.sam
cat  ${flashstatus}_REdig_unfiltered.sam | grep -cv '^@'
echo

printThis="Sam to bam transform .."
printToLogFile

flashstatus="NONFLASHED"
rm -f ${flashstatus}FLASHED_REdig.fastq
echo "samtools view -hb ${flashstatus}_REdig_unfiltered.sam > ${flashstatus}_REdig_unfiltered.bam"
samtools view -hb ${flashstatus}_REdig_unfiltered.sam > ${flashstatus}_REdig_unfiltered.bam
ls -lht ${flashstatus}_REdig_unfiltered.bam
rm -f ${flashstatus}_REdig_unfiltered.sam

TEMPweAreHere=$(pwd)
cd ..

# RE enzyme digestion (if needed .. )

dpnGenomeName=""
fullPathDpnGenome=""
generateReDigest

cd ${TEMPweAreHere}

# RE enzyme genome blacklist generation (regions farther than amplicon lenght from the cut site)

fullPathDpnBlacklist=""
generateReBlacklist

# Filtering based on RE enzyme genome blacklist ..
printThis="Filtering out reads which are farther away from cut sites than amplicon size ${ampliconSize} .."
printToLogFile

flashstatus="FLASHED"
echo "bedtools intersect -v -abam ${flashstatus}_REdig_unfiltered.bam -b ${fullPathDpnBlacklist} > ${flashstatus}_REdig.bam"
bedtools intersect -v -abam ${flashstatus}_REdig_unfiltered.bam -b ${fullPathDpnBlacklist} > ${flashstatus}_REdig.bam
echo "samtools view -h ${flashstatus}_REdig.bam > ${flashstatus}_REdig.sam"
samtools view -h ${flashstatus}_REdig.bam > ${flashstatus}_REdig.sam
ls -lht ${flashstatus}_REdig.sam
rm -f ${flashstatus}_REdig.bam

flashstatus="NONFLASHED"
echo "bedtools intersect -v -abam ${flashstatus}_REdig_unfiltered.bam -b ${fullPathDpnBlacklist} > ${flashstatus}_REdig.bam"
bedtools intersect -v -abam ${flashstatus}_REdig_unfiltered.bam -b ${fullPathDpnBlacklist} > ${flashstatus}_REdig.bam
echo "samtools view -h ${flashstatus}_REdig.bam > ${flashstatus}_REdig.sam"
samtools view -h ${flashstatus}_REdig.bam > ${flashstatus}_REdig.sam
ls -lht ${flashstatus}_REdig.sam
rm -f ${flashstatus}_REdig.bam

# Cleaning up after ourselves ..

printThis="Finishing up the F1 run folder.."
printToLogFile

#ls -lht Combined_reads_REdig.bam
ls -lht FLASHED_REdig.sam
ls -lht NONFLASHED_REdig.sam

echo
echo "Read counts - in PCR-amplicon-size filtered sam files : "
echo
flashstatus="FLASHED"
echo ${flashstatus}_REdig.sam
cat  ${flashstatus}_REdig.sam | grep -cv '^@'
echo
flashstatus="NONFLASHED"
echo ${flashstatus}_REdig.sam
cat  ${flashstatus}_REdig.sam | grep -cv '^@'
echo

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

    # These are temp symlinks to filteringLogFor_PREfiltered_${Sample}_${CCversion} files
    rm -f ../FLASHED_REdig.sam ../NONFLASHED_REdig.sam
    
# Remove the malformed public folder for a new try..
    rm -rf ${PublicPath} ../PERMANENT_BIGWIGS_do_not_move
    
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

TEMPweAreHere=$(pwd)
cd ..
    
# RE enzyme digestion (if needed .. )

dpnGenomeName=""
fullPathDpnGenome=""
generateReDigest

# RE enzyme genome blacklist generation (regions farther than amplicon lenght from the cut site)

fullPathDpnBlacklist=""
generateReBlacklist

cd ${TEMPweAreHere}
    
fi

################################################################
# Store the pre-CCanalyser log files for metadata html

printThis="Store the pre-CCanalyser log files for metadata html.."
printToLogFile

copyPreCCanalyserLogFilesToPublic

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

echo "${CaptureFilterPath}/filter.sh -p parameters_for_filtering.log -s ${CaptureFilterPath} --pipelinecall ${ploidyFilter} --extend ${extend} "
echo "${CaptureFilterPath}/filter.sh -p parameters_for_filtering.log -s ${CaptureFilterPath} --pipelinecall ${ploidyFilter} --extend ${extend} "  >> "/dev/stderr"

#        --stepSize) stepSize=$2 ; shift 2;;
#        --tileSize) tileSize==$2 ; shift 2;;
#        --minScore) minScore=$2 ; shift 2;;
#        --maxIntron) maxIntron=$2 ; shift 2;;
#        --oneOff) oneOff=$2 ; shift 2;;

echo "--stepSize ${stepSize} --minScore ${minScore} --maxIntron ${maxIntron} --tileSize ${tileSize} --minIdentity=${minIdentity} --minMatch=${minMatch} --oneOff ${oneOff}"
echo "--stepSize ${stepSize} --minScore ${minScore} --maxIntron ${maxIntron} --tileSize ${tileSize} --minIdentity=${minIdentity} --minMatch=${minMatch} --oneOff ${oneOff}" >> "/dev/stderr"

echo "--reuseBLAT ${reuseBLATpath}"
echo "--reuseBLAT ${reuseBLATpath}" >> "/dev/stderr"
echo "--onlyCis ${onlyCis}"
echo "--onlyCis ${onlyCis}" >> "/dev/stderr"

mkdir filteringLogFor_${sampleForCCanalyser}_${CCversion}
mv parameters_for_filtering.log filteringLogFor_${sampleForCCanalyser}_${CCversion}/.
cd filteringLogFor_${sampleForCCanalyser}_${CCversion}

TEMPreturnvalue=0
${CaptureFilterPath}/filter.sh --reuseBLAT ${reuseBLATpath} -p parameters_for_filtering.log --pipelinecall ${ploidyFilter} --extend ${extend} --onlyCis ${onlyCis} --stepSize ${stepSize} --minScore ${minScore} --maxIntron ${maxIntron} --tileSize ${tileSize} --minIdentity=${minIdentity} --minMatch=${minMatch} --oneOff ${oneOff} > filtering.log
TEMPreturnvalue=$?
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