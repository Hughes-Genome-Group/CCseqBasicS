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


version(){
    pipesiteAddress="http://userweb.molbiol.ox.ac.uk/public/telenius/NGseqBasicManual/outHouseUsers/"
    
    versionInfo="\n${CCversion} pipeline, running ${CCseqBasicVersion}.sh \nUser manual, updates, bug fixes, etc in : ${pipesiteAddress}\n"

}


usage(){
    
    version
    echo -e ${versionInfo}
    
echo 
echo "FOR PAIRED END SEQUENCING DATA"
echo
echo "fastq --> fastqc --> trimming  --> fastqc again --> flashing --> fastqc again --> in silico RE fragments --> bowtie1 -->  in silico genome digestion --> captureC analyser --> data hub generation"
echo
echo "Is to be ran in the command line like this :"
echo "qsub -cwd -o qsub.out -e qsub.err -N MyPipeRun < ./run.sh"
echo "Where run.sh is oneliner like : '${DNasePipePath}/${nameOfScript} RunOptions' "
echo "where RunOptions are the options given to the pipe - listed below"
echo
echo "Run the script in an empty folder - it will generate all the files and folders it needs."
echo
echo "OBLIGATORY FLAGS FOR THE PIPE RUN :"
echo
echo "-c /path/to/capturefragment/file.txt : the file containing the RE-fragments within which the BIOTINYLATED CAPTURESITES reside, and their proximity exclusions (standard practise : 1000bases both directions), and possible SNP sites (see pipeline manual how to construct this file : ${manualURLpath} )"
echo "-o (synonymous to -c above)"
echo "--R1 /path/to/read1.fastq : fastq file from miseq or hiseq run (if this is .gz packed, add also flag --gz)"
echo "--R2 /path/to/read2.fastq : fastq file from miseq or hiseq run (if this is .gz packed, add also flag --gz)"
echo "--genome mm9 : genome to use to map and analyse the sample (supports most WIMM genomes - mm9,mm10,hg18,hg19 - report to Jelena if some genomes don't seem to work ! )"
echo "--pf /public/username/location/for/UCSC/visualisation/files : path to a folder where you want the data hub to be generated. Does not need to exist - pipeline will generate it for you."
echo "-s SampleName : the name of the sample - no weird characters, no whitespace, no starting with number, only characters azAz_09 allowed (letters, numbers, underscore)"
echo
echo "THE MINIMAL RUN COMMAND :"
echo "${RunScriptsPath}/${nameOfScript} -o /path/to/capturesite/file.txt --R1 /path/to/read1.fastq --R2 /path/to/read2.fastq --genome mm9 --pf /public/username/location/for/UCSC/visualisation/files"
echo "where genome (f.ex 'mm9')is the genome build the data is to be mapped to."
echo
echo "The above command expanded (to show the default settings) :"
echo "${RunScriptsPath}/${nameOfScript} -o /path/to/capturesite/file.txt --R1 /path/to/read1.fastq --R2 /path/to/read2.fastq --genome mm9 --pf /public/username/location/for/UCSC/visualisation/files"
echo "   -s sample --maxins 250 -m 2 --chunkmb 256 --trim -w 200 -i 20"
echo
echo "OPTIONAL FLAGS FOR TUNING THE PIPE RUN :"
echo
echo "HELP"
echo "-h, --help : prints this help"
echo
echo "FASTQ FILE INPUT OPTIONS"
echo "One fastq file pair can be given with --R1 and --R2 parameters. If these are gzipped, add parameter --gz"
echo "Multiple fastq file pairs can be given in PIPE_fastqPaths.txt which has file format :"
echo
echo "R1a.fastq.gz  R2a.fastq.gz  /path/to/these/files"
echo "R1b.fastq.gz  R2b.fastq.gz  /path/to/these/files"
echo "or"
echo "/path/to/this/file/R1a.fastq.gz  /path/to/this/file/R2a.fastq.gz"
echo "/path/to/this/file/R1b.fastq.gz  /path/to/this/file/R2b.fastq.gz"
echo
echo "If these files are gzipped, add parameter --gz"
echo 
echo "If these files are SRR archive format (instead of Illumina format), add parameter --SRR"
echo 
echo "OUTPUT LOG FILE NAMES"
echo "--outfile qsub.out (the STDOUT log file name in your RUN COMMAND - see above )"
echo "--errfile qsub.err (the STDERR log file name in your RUN COMMAND - see above )"
echo ""
echo "RE-RUNNING PARTS OF THE PIPE (in the case something went wrong)"
echo
echo "--onlyCCanalyser : Start the analysis from the CCanalyser script. "
echo "Using this flag will assume that the beginning of the pipe was successfull, and the run crashed at 'captureC analyser' step."
echo "The flow of the pipeline is :"
echo "fastq --> fastqc --> trimming  --> fastqc again --> flashing --> fastqc again --> in silico RE fragments --> bowtie1 -->  in silico genome digestion --> CAPTUREC ANALYSER --> DATA HUB GENERATION"
echo "This run will thus repeat only the last 2 steps - it will delete the previous output folder, and repeat the analysis from CaptureC analyser onwards, and generate the data hub."
echo "This re-run will take only ~30 minutes, when the whole pipe takes ~5h to run !"
echo "Using this flag is good idea, when your sam file (Combined_reads_REdig.sam) is OK (it looks like a BIG file full of reads),"
echo "but your capture run went somehow wonky (missing hub, missing bigwigs, typos in capture-site (REfragment) file etc)."
echo "Fix your files (capture-site (REfragment) file typos, wrong paths etc), and start the run in the same folder you ran it the first time"
echo "- it will delete the previous captureC analyser output folder, and tries to regenerate it, using the corrected parameters you gave."
echo
echo "--BLATforREUSEfolderPath /full/path/to/previous/F4_blatPloidyFilteringLog_CC4/BlatPloidyFilterRun/REUSE_blat folder"
echo "    if run crashes during or after blatting : point to the crashed folder's blat results, and avoid running blat again ! "
echo "    Remember to check (if you crashed during blat) - that all your blat output (psl) files are intact. Safest is to delete the last one of them (check which, with ls -lht)"
echo "    NOTE !!! if you combine this with --onlyCCanalyser : As the blat results live inside folder F4, --onlyCCanalyser will delete these files before re-starting the run (all folders but F1 get deleted)."
echo "        you need to COPY the REUSE_blat folder outside the original run folder, and point to that copied folder here. "
echo
# echo "--onlyHub  : runs the tracks.txt generation again. Assumes that the whole pipe ran correctly, but the tracks were written wrongly to the output tracks.txt for the data hub "
# echo "This is mainly for debugging (only Jelena should need to use to this one)"
# echo
echo "RESTRICTION ENZYME SETTINGS"
echo "--dpn  (default) : dpnII is the RE of the experiment"
echo "--nla  : nlaIII is the RE of the experiment"
echo "--hind : hindIII is the RE of the experiment"
echo
echo "BOWTIE SETTINGS"
echo "--bowtie1 / --bowtie2 (default is bowtie1 - decide if bowtie1 or bowtie2 is to be used. bowtie2 is better to long reads - read lenght more than 70b, library fragment size at PCR-amplification step more than 350b)"
echo "--chunkmb "${BOWTIEMEMORY}" - memory allocated to Bowtie, defaults to 256mb "
echo "-M 2 run with bowtie parameter M=2 (if maps more than M times, report one alignment in random) - only affects bowtie1 run"
echo "-m 2 run with bowtie parameter m=2 (if maps more than m times, do not report any alignments) - only affects bowtie1 run"
echo "-m and -M are mutually exclusive."
echo "--trim3 0 : trim the reads this many bases from the 3' end when mapping in bowtie"
echo "--trim5 0 : trim the reads this many bases from the 5' end when mapping in bowtie"
echo "-v 3 : allow up-to-this-many total mismatches per read (ignore base qualities for these mismatches). "
echo "       cannot be combined to --seedlen, --seedmms or --maqerr (below)."
echo "--seedlen 28 - alignment seed lenght (minimum 5 bases) . Seed is the high-quality bases in the 5' end of the read. Default 28 (bowtie1), 20 (bowtie2)."
echo "--seedmms 2 - max mismatches within the seed (see the 'seed' above). Allowed 0,1,2,3 mismatches in bowtie1 - default 2, allowed 0,1 in bowtie2 (per each multi-seed alignment) - default 0. "
echo "--maqerr 70 - only in Bowtie1 - max total quality values at all mismatched read positions throughout the entire alignment (not just in seed)"
echo ""
echo "ADAPTER TRIMMING SETTINGS"
echo "--trim/noTrim** (run/do-not-run TrimGalore for the data - Illumina PE standard adapter filter, trims on 3' end)"
echo "**) NOTE : use --noTrim with caution - the data will go through FLASH : this can result in combining reads on the sites of ADAPTERS instead of the reads themselves."
echo "--ada3read1 SEQUENCE --ada3read2 SEQUENCE  : custom adapters 3' trimming, R1 and R2 (give both) - these adapters will be used instead of Illumina default / atac adapters. SEQUENCE has to be in CAPITAL letters ATCG"
echo "--ada5read1 SEQUENCE --ada5read2 SEQUENCE  : custom adapters 5' trimming, R1 and R2 (give both) - these adapters will be used instead of Illumina default / atac adapters. SEQUENCE has to be in CAPITAL letters ATCG"
echo ""
echo "QUALITY TRIMMING SETTINGS"
echo "--qmin 20 (trim low quality reads up-to this phred score) - sometimes you may want to rise this to 30 or so"
echo ""
echo "FLASH SETTINGS"
echo "--flashBases 10 (when flashing, has to overlap at least this many bases to combine)"
echo "--flashMismatch 0.25 (when flashing, max this proportion of the overlapped bases are allowed to be MISMATCHES - defaults to one in four allowed to be mismatch, i.e. 0.25 )"
echo "                      sometimes you may want to lower this to 0.1 (one in ten) or 0.125 (one in eight) or so"
echo ""
echo "SAVE IN-SILICO DIGESTED WHOLE-GENOME FILE"
echo "--saveGenomeDigest"
echo "(to save time in your runs : add the output file to your digests folder, and update your conf/config.sh  !)"
echo ""
echo "PCR AMPLICON SIZE"
echo "--ampliconSize 300 (how far from RE enzyme cut site reach the 'valid fragments'. This is the max expected library fragment size, at the PCR amplification step.)"
echo "--sonicationSize 300 (synonymous to --ampliconSize)"
echo
echo "CAPTURE-C ANALYSER OPTIONS"
echo "--useSymbolicLinks (use symbolic links between run directory and public directory to store bigwig files, "
echo "   instead of storing the bigwigs straight in the public directory)"
echo "--onlyCis (to analyse only cis-mapping reporters : this flag also affects BLAT OUTPUT FILES, see below)"
echo "-s Sample name (and the name of the folder it goes into)"
echo "--snp : snp-specific run (check your capture-site (REfragment) coordinates file that you have defined the SNPs there)"
echo "--globin : combine captures globin capture sites :"
echo "  To combine ONLY alpha globin  :  --globin 1 (name your globin capture sites Hba-1 and Hba-2)"
echo "-s Sample name (and the name of the folder it goes into)"
echo
echo "WINDOWING IN CAPTUREC ANALYSER SCRIPT"
echo "Default windowing is 200b window and 20b increment"
echo "-w 200   or   --window 200  : custom window size (instead of default 200b)."
echo "-i 20    or   --increment 20 : custom window increment (instead of default 20b). "
echo
echo "BLAT FILTERING - blat parameters :"
echo "blat -oneOff=0 -minScore=10 -maxIntron=4000 -tileSize=11 -stepSize=5 -minIdentity=70 -repMatch=999999"
echo
echo "BLAT OPTIONS :"
echo "--onlyCis (to generate blat-filtering files for only cis chromosomes : this flag also affects CAPTURE-C ANALYSER, see above)"
echo "--stepSize 5 (spacing between tiles). if you want your blat run faster, set this to 11."
echo "--tileSize 11 (the size of match that triggers an alignment)"
echo "--minScore 10 (minimum match score)"
echo "--maxIntron 4000 (to make blat run quicker) (blat default value is 750000) - max intron size"
echo "--oneOff 0 (set this to 1, if you want to allow one mismatch per tile. Setting this to 1 will make blat slow.)"
echo "--BLATforREUSEfolderPath /full/path/to/previous/F4_blatPloidyFilteringLog_CC4/BlatPloidyFilterRun/REUSE_blat folder"
echo "   (enables previously ran BLAT for the same capture-site (REfragment)s, to be re-used in the run)"
echo
echo "CAPTURE-C BLAT + PLOIDY FILTERING OPTIONS"
echo "--extend 20000  Extend the Blat-filter 20000bases both directions from the psl-file regions outwards. (default 20 000)"
echo "--noPloidyFilter  Do not filter for ploidy regions (Hughes lab peak call for mm9/mm10, Duke Uni blacklisted for hg18/hg19, other genomes don't have ploidy track provided in pipeline)"
echo
echo "CAPTURE-C DUPLICATE FILTERING OPTIONS"
echo "--CCversion CS5 : Which duplicate filtering is to be done : CS3 (for short sequencing reads), CS4 (long reads), CS5 (any reads). "
echo "--strandSpecificDuplicates : To replicate the strand-specific (i.e. wrong) duplicate filter of CB3a/CC3 and CB4a/CC4"
echo "--UMI : Use UMI-style duplicate filtering."
echo "--wobblyEndBinWidth 1 : to set bin width for non-exact fragment-end coordinates."
echo "   To turn wobbly ends off, set this to 1 ( --wobblyEndBinWidth 1 ) - this is the default."
echo "   If using --UMI , --wobblyEndBinWidth 20 is recommended."
echo "   For example : --wobblyEndBinWidth 20 means : bin of 20 bases for duplicate filter :"
echo "   if all fragment coordinates are the same +/- 10 bases, ( and if --UMI is used : UMI is the same), reads are duplicates."
echo
echo "CAPTURE-C ANALYSER DEVELOPER OPTIONS"
echo "--dump : Print file of unaligned reads (sam format)"
echo "--limit n  : only analyse the first 'n' reads - for testing purposes "
echo

#echo "More info : hands-on tutorial : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DnaseCHIPpipe_TUTORIAL.pdf, comprehensive user manual : http://userweb.molbiol.ox.ac.uk/public/telenius/MANUAL_for_pipe_030214/DNasePipeUserManual_VS_100_180215.pdf , and comment lines (after the subroutine descriptions) in the script ${DNasePipePath}/DnaseAndChip_pipe_1.sh"
echo 
 

 
 exit 0
    
}

