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


Installation instructions :

1) If you downloaded from GitHub - jump to step (2) below
   
   If you asked  Jelena (jelena __ telenius __ at __ gmail __ com) to send you the tar.gz file including the codes :
   
   Unpack with this command :
    tar --preserve-permissions -xzf CCseqBasic4.tar.gz

2) Fill in the locations (or modules) of the needed tools (bowtie, fastqc etc), and the genome builds, to the config files :
    nano CCseqBasic4/conf/loadNeededTools.sh        # Instructions as comment lines in the  loadNeededTools.sh file
    nano CCseqBasic4/conf/genomeBuildSetup.sh       # Instructions as comment lines in the genomeBuildSetup.sh file

3) Fill in your server address to the conf/serverAddressAndPublicDiskSetup.sh file
    nano CCseqBasic4/conf/serverAddressAndPublicDiskSetup.sh       # Instructions as comment lines in the file

4) Add the main script CCseqBasic4.sh to your path or BASH profile (optional), f.ex :
    export PATH:${PATH}:/where/you/unpacked/it/CCseqBasic4/CCseqBasic4.sh

5) Start using the pipe ! (no installation needed)

6) Good place to start is the pipeline's help :
    CCseqBasic4.sh --help

7) Below web site provides a test data set, hands-on tutorials, full manual, and other documentation !
   http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/
   
8) Direct link to the test data set :
   http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/testdata/exampledata.html
   