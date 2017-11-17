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

runTempCommands(){
chmod u=rwx TEMP_commands.sh
echo
cat TEMP_commands.sh
echo >> "/dev/stderr"
echo "Running following commands :" >> "/dev/stderr"
cat TEMP_commands.sh >> "/dev/stderr"
echo "Possible error messages while running abovelisted commands :" >> "/dev/stderr"
echo >> "/dev/stderr"
./TEMP_commands.sh
rm -f TEMP_commands.sh    
}

# This is the first blat-filter automaton - test code built 11/Sep/2015 by Jelena

#for oligoName in oligoFile
#this will be some kind of awk thing - commands generated via reading in the oligo file..

#CapturePipePath="/home/molhaem2/telenius/CC2/norm/VS101"

PathForReuseBlatResults="."

CapturePipePath="UNDEFINED"

oligofile="UNDEFINED"
genomefasta="UNDEFINED"
recoordinatefile="UNDEFINED"
ucscBuild="UNDEFINED"
extend="20000"

# full path to parameter file
blatparams="UNDEFINED"

OPTS=`getopt -o o:,f:,r:,p:,u:,e: --long genomefasta:,oligofile:,refile:,pipepath:,ucscbuild:,extend:,blatparams:,reusefile: -- "$@"`
if [ $? != 0 ]
then
    exit 1
fi

eval set -- "$OPTS"

while true ; do
    case "$1" in
        -o) oligofile=$2 ; shift 2;;
        -f) genomefasta=$2 ; shift 2;;
        -r) recoordinatefile=$2 ; shift 2;;
        -p) CapturePipePath=$2 ; shift 2;;
        -u) ucscBuild=$2 ; shift 2;;
        -e) extend=$2 ; shift 2;;        
        --oligofile) oligofile=$2 ; shift 2;;
        --genomefasta) genomefasta=$2 ; shift 2;;
        --refile) recoordinatefile=$2 ; shift 2;;
        --pipepath) CapturePipePath=$2 ; shift 2;;
        --ucscbuild) ucscBuild=$2 ; shift 2;;
        --extend) extend=$2 ; shift 2;;
        --blatparams) blatparams=$2 ; shift 2;;
        --reusefile) PathForReuseBlatResults=$2 ; shift 2;;
        --) shift; break;;
    esac
done

echo
echo "blatting with 1_blat.sh"
echo
echo "Starting run with parameters :"
echo

echo "oligofile ${oligofile}"
echo "genomefasta ${genomefasta}"
echo "recoordinatefile ${recoordinatefile}"
echo "CapturePipePath ${CapturePipePath}"
echo "ucscBuild ${ucscBuild}"
echo

#module unload bedtools
#module unload blat
#module load bedtools/2.17.0
#module load blat/35
module list

echo
echo "Dividing the oligo file to one-liners.."
echo "Dividing the oligo file to one-liners.."  >> "/dev/stderr"
# Dividing the oligo file to one-liners like :
# chr     str     stp (for the exclusion fragments in the oligo coordinate file)

echo '#!/bin/bash ' > TEMP_commands.sh
echo "oligofile=${oligofile}"  >> TEMP_commands.sh
cat ${oligofile} | awk '{print "cat ${oligofile} | grep \"^"$1"\\s\" | awk * > TEMP_"$1"_coordinate.bed"}' | sed 's/*/*\{print \"chr\"$5\"\\t\"$6\"\\t\"$7\}*/' | tr "*" "'" >> TEMP_commands.sh
echo "" >> TEMP_commands.sh

# above results in lines like :
# cat ${oligofile} | grep Hba-1 | awk '{print "chr"$5"\t"$6"\t"$7}' > TEMP_Hba-1_coordinate.bed

runTempCommands

# Run fasta generation to all these exlusion fragment coordinates
echo
echo "Run fasta generation to all these exlusion fragment coordinates.."
echo "Run fasta generation to all these exlusion fragment coordinates.."  >> "/dev/stderr"

for file in TEMP*coordinate.bed
do

echo "bedtools getfasta -fi ${genomefasta} -bed ${file} -fo ${file}.fa"
bedtools getfasta -fi ${genomefasta} -bed ${file} -fo ${file}.fa

done

rm -f TEMP*coordinate.bed

for file in TEMP*coordinate.bed.fa
do

sed -i 's/a/A/g' $file
sed -i 's/t/T/g' $file
sed -i 's/c/C/g' $file
sed -i 's/g/G/g' $file

done

# Running blat for these sequences (exclusion fragment sequences)
echo 
echo "Running blat for these sequences (exclusion fragment sequences).."
echo "Running blat for these sequences (exclusion fragment sequences).."  >> "/dev/stderr"

blatParams=$( cat ${blatparams} )

for file in TEMP*coordinate.bed.fa
do

basename=$( echo $file | sed 's/_coordinate.bed.fa//' )



# If cannot find the file - i.e. if this is the first run for these oligos,
#   or if the earlier run didn't find any blat hits (this is unintended consequence of the safety feature :
#   - we don't want to SKIP blat filter for all oligos BECAUSE OF FILE ADDRESS TYPO here.. )

weRunBLAT=1
if [ -e ${PathForReuseBlatResults}/${basename}_blat.psl ]
then
if [ -s ${PathForReuseBlatResults}/${basename}_blat.psl ]
then   
weRunBLAT=0
fi
fi

if [ "${weRunBLAT}" -eq "1" ]
then

echo "blat ${blatParams} ${genomefasta} ${file} ${basename}_blat.psl"
echo "blat ${blatParams} ${genomefasta} ${file} ${basename}_blat.psl" >> "/dev/stderr"

blat ${blatParams} ${genomefasta} ${file} ${basename}_blat.psl
else

echo "Found file ${PathForReuseBlatResults}/${basename}_blat.psl - skipping blat for that oligo, using the found file instead !"
echo "Found file ${PathForReuseBlatResults}/${basename}_blat.psl - skipping blat and that oligo, using the found file instead !" >> "/dev/stderr"

cp ${PathForReuseBlatResults}/${basename}_blat.psl .

fi

done

rm -f TEMP*coordinate.bed.fa

# Works upto here.. (10/09/15)

# Remove the files which didn't get any results in blat - these exclusion regions won't thus need a blat filter!
echo 
echo "Removing the files which didn't get any results in blat - these exclusion regions won't thus need a blat filter!.."
echo "Removing the files which didn't get any results in blat - these exclusion regions won't thus need a blat filter!.."  >> "/dev/stderr"

for file in TEMP*blat.psl
do

if [ ! -s $file ]
then

rm -f $file

fi

done

echo
echo "Files for blat filtering :"
ls -lht |  grep blat.psl

# Preparing the input files for the RE-fragment generation perl script
echo 
echo "Preparing the input files for the RE-fragment generation perl script.."
echo "Preparing the input files for the RE-fragment generation perl script.."  >> "/dev/stderr"


echo '#!/bin/bash ' > TEMP_commands.sh
echo "oligofile=${oligofile}"  >> TEMP_commands.sh
cat ${oligofile} | awk '{print "cat ${oligofile} | grep \"^"$1"\\s\" > TEMP_"$1"_oligocoordinate.txt"}' >> TEMP_commands.sh
echo "" >> TEMP_commands.sh

runTempCommands

#for file in TEMP*blat.psl
#do   

#cat ${file} | grep chr | sed 's/\s\s*/\t/' > TEMPFILE
#mv -f TEMPFILE ${file}

#done



# Running the psl parser - change blat output to RE fragment coordinates to be removed from the files..
echo 
echo "Running the psl parser - change blat output to RE fragment coordinates to be removed from the files.."
echo "Running the psl parser - change blat output to RE fragment coordinates to be removed from the files.."  >> "/dev/stderr"

for file in TEMP*blat.psl
do   

basename=$( echo $file | sed 's/_blat.psl//' )
echo
echo "-----------------------"
echo "perl ${CapturePipePath}/2_psl_parser.pl -f ${file} -o ${basename}_oligocoordinate.txt -a ${oligofile} -r ${recoordinatefile}"
echo "perl ${CapturePipePath}/2_psl_parser.pl -f ${file} -o ${basename}_oligocoordinate.txt -a ${oligofile} -r ${recoordinatefile}" >> "/dev/stderr"
echo
perl ${CapturePipePath}/2_psl_parser.pl -f ${file} -o ${basename}_oligocoordinate.txt -a ${oligofile} -r ${recoordinatefile}

done

rm -f TEMP*_oligocoordinate.txt

# To reuse the blat coordinates..
rm -rf REUSE_blat
mkdir REUSE_blat
mv -f TEMP*blat.psl REUSE_blat/.

# Saving only files which have content :

for file in *_blat_filter.gfc
do

if [ ! -s $file ]
then

rm -f $file

fi

done

# Change names to not contain TEMP

for file in *_blat_filter.gfc
do

newname=$( echo $file | sed 's/TEMP_//' )
mv -f $file ${newname}

done

# Here we bedtools slop it - to extend 20 000 bases both directions..
echo
echo "Extending blat filter by ${extend} bases to both directions (default 20000 bases) .."
echo "Extending blat filter by ${extend} bases to both directions.."  >> "/dev/stderr"
echo

# Here we also transform our blat filter to gff file :
#[telenius@deva captureVS1_VS2comparison_030915]$ head testblatfilter/7thtry/TEMP_SOX2_blat_filter.gfc | sed 's/:/\tcol2\tcol3\t/' | sed 's/-/\t/' | sed 's/$/\t6\t7\t8\tfacet9=0/'
#chr8    col2    col3    12396275        12396506        6       7       8       facet9=0
#chr8    col2    col3    12396507        12396601        6       7       8       facet9=0

for file in *_blat_filter.gfc
do

newname=$( echo $file | sed 's/.gfc$//' )
cat ${file} | sed 's/:/\tcol2\tcol3\t/' | sed 's/-/\t/' | sed 's/$/\t.\t.\t.\tfacet=0/' > ${newname}.gff

bedtools slop -i ${newname}.gff -g ${ucscBuild} -b ${extend} > ${newname}_extendedBy${extend}b.gff

done

echo
echo "Files for blat filtering :"
ls -lht |  grep "_extendedBy${extend}b.gff"





