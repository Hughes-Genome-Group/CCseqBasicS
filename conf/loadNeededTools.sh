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


setPathsForPipe(){

# #############################################################################

# This is the CONFIGURATION FILE to load in the needed toolkits ( conf/loadNeededTools.sh )

# #############################################################################

# Setting the needed programs to path.

# This can be done EITHER via module system, or via EXPORTING them to the path.
# If exporting to the path - the script does not check already existing conflicting programs (which may contain executable with same names as these)

# If neither useModuleSystem or setPathsHere : the script assumes all toolkits are already in path !

# If you are using module system
useModuleSystem=1
# useModuleSystem=1 : load via module system
# useModuleSystem=0 : don't use module system

# If you are adding to path (using the script below)
setPathsHere=0
# setPathsHere=1 : set tools to path using the bottom of this script
# setPathsHere=0 : dset tools to path using the bottom of this script

# If neither useModuleSystem or setPathsHere : the script assumes all toolkits are already in path !

# #############################################################################

# PATHS_LOADED_VIA_MODULES

if [ "${useModuleSystem}" -eq 1 ]; then

module purge
# Removing all already-loaded modules to start from clean table

module load samtools/1.3
# Supports samtools 1.* versions 1.3 or newer. Does not support samtools/0.* or samtools 1.1.

module load bowtie/1.1.2
# Supports all bowtie1 versions 1.* and 0.*
module load bowtie2/2.1.0
# Supports all bowtie2 versions

module load bedtools/2.17.0
# Supports bedtools versions 2.1* . Does not support bedtools versions 2.2*

module load flash/1.2.8
# Not known if would support other flash versions. Most probably will support.

module load fastqc/0.11.4
# Will not support fastqc versions 0.10.* or older

module load trim_galore/0.3.1
# Not known if would support other trim_galore versions. Most probably will support.

# module load cutadapt/1.2.1
# If your trim_galore module does not automatically load the needed cutadapt,
# uncomment this line

# Not known if would support other cutadapt versions. Most probably will support.

module load blat/35
# Not known if would support other blat versions. Most probably will support.

module load perl/5.18.1
# Most probably will run with any perl

module load python/2.7.5
# Most probably will run with any Python

module list

# #############################################################################

# EXPORT_PATHS_IN_THIS_SCRIPT

elif [ "${setPathsHere}" -eq 1 ]; then

echo
echo "Adding tools to PATH .."
echo
    
# Note !!!!!
# - the script does not check already existing conflicting programs within $PATH (which may contain executable with same names as these)

export PATH=$PATH:/package/samtools/1.3/bin
export PATH=$PATH:/package/bowtie/1.1.2/bin
export PATH=$PATH:/package/bowtie2/2.1.0/bin
export PATH=$PATH:/package/bedtools/2.17.0/bin
export PATH=$PATH:/package/flash/1.2.8/bin
export PATH=$PATH:/package/fastqc/0.11.4/bin
export PATH=$PATH:/package/trim_galore/0.3.1/bin
export PATH=$PATH:/package/cutadapt/1.2.1/bin
export PATH=$PATH:/package/perl/5.18.1/bin
export PATH=$PATH:/package/blat/35/bin
export PATH=$PATH:/package/python/2.7.5/bin


# See notes of SUPPORTED VERSIONS above !

echo $PATH

# #############################################################################

# EXPORT_NOTHING_i.e._ASSUMING_USER_HAS_TOOLS_LOADED_VIA_OTHER_MEANS

else
    
echo
echo "Tools should already be available in PATH - not loading anything .."
echo

fi

# #########################################

# UCSCtools are taken from install directory, in any case :

export PATH=$PATH:/${confFolder}/ucsctools



}

