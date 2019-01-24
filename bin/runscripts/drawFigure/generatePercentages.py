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
import sys
# print  >> sys.stderr, "Preparing to run - check where we are and which version we have.."
# print  >> sys.stderr, ""
print  >> sys.stderr, "We are running in Python and machine :"
print  >> sys.stderr, sys.platform
print  >> sys.stderr, sys.version
# print  >> sys.stderr, ""
# print  >> sys.stderr, "----------------------------------------"
# print  >> sys.stderr, "Run directory :"
import os
# print  >> sys.stderr, os.getcwd()
# print  >> sys.stderr, "----------------------------------------"
# print  >> sys.stderr, "Enabling log file output.."
import syslog
# print  >> sys.stderr, "Enabling resource use and run time statistics monitoring.."
import stat

# print  >> sys.stderr, "Importing regular expressions"
import re

# print >> sys.stderr,  "----------------------------------------"
print >> sys.stderr,  "Imported (and auto-loaded) modules :"
print >> sys.stderr, (globals())

# print >> sys.stderr, "----------------------------------------"
# print >> sys.stderr, "Reading in the subroutines.."

# Making the comments above the lines..
def printList(lista) :
    for x in lista :
        print x,
    print

# print  >> sys.stderr, "----------------------------------------"
# print  >> sys.stderr, "Starting the run.."
# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Reading the input.."
# print >> sys.stderr, ""

# counts.py is a text file, in which the counted values like this :

# This reads the in-pwd-located above counts.py file in :
sys.path.append(".")
from counts import *

# print  >> sys.stderr, "all ",all
# print  >> sys.stderr, "allflashed",allflashed
# print  >> sys.stderr, "allnonflashed",allnonflashed
# print  >> sys.stderr, "REflashed",REflashed
# print  >> sys.stderr, "REnonflashed",REnonflashed
# print  >> sys.stderr, "continuesToMappingFlashed",continuesToMappingFlashed
# print  >> sys.stderr, "continuesToMappingNonflashed",continuesToMappingNonflashed
# print  >> sys.stderr, "containsCaptureFlashed",containsCaptureFlashed
# print  >> sys.stderr, "containsCaptureNonflashed",containsCaptureNonflashed
# print  >> sys.stderr, "containsCapAndRepFlashed",containsCapAndRepFlashed
# print  >> sys.stderr, "containsCapAndRepNonflashed",containsCapAndRepNonflashed
# print  >> sys.stderr, "singleCapFlashed",singleCapFlashed
# print  >> sys.stderr, "singleCapNonflashed",singleCapNonflashed
# print  >> sys.stderr, "multiCapFlashed",multiCapFlashed
# print  >> sys.stderr, "multiCapNonflashed",multiCapNonflashed
# print  >> sys.stderr, "nonduplicateFlashed",nonduplicateFlashed
# print  >> sys.stderr, "nonduplicateNonflashed",nonduplicateNonflashed
# print  >> sys.stderr, "blatploidyFlashed",blatploidyFlashed
# print  >> sys.stderr, "blatploidyNonflashed",blatploidyNonflashed

# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
# print >> sys.stderr, "Generating arrays.."
# print >> sys.stderr, ""

# Generating the lists..
r0_All_reads=100
r1_flashed=[0,1]
r2_RE_site=[0,1,2,3]
r3_continues=[0,1]
r4_Cont_cap=[0,1,2,3]
r5_Cap_and_rep=[0,1,2,3,4,5,6]
r6_multicap=[0,1,2,3]
r7_duplicate=[0,1,2,3]
r8_blatploidy=[0,1,2,3]

# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Setting values.."
# print >> sys.stderr, ""

# All flashed
r1_flashed[0]=((allflashed*1.0)/all)*100
# All nonflashed
r1_flashed[1]=((allnonflashed*1.0)/all)*100

# Flashed having RE site
r2_RE_site[0]=((REflashed*1.0)/all)*100
# Flashed not having RE site
r2_RE_site[1]=(((allflashed-REflashed)*1.0)/all)*100
# Not flashed having RE site
r2_RE_site[2]=((REnonflashed*1.0)/all)*100
# Not flashed not having RE site
r2_RE_site[3]=(((allnonflashed-REnonflashed)*1.0)/all)*100

# Continue to mapping : Flashed having RE site
r3_continues[0]=r2_RE_site[0]
# Continue to mapping : All nonflashed
r3_continues[1]=r1_flashed[1]

# Contains capture Flashed
r4_Cont_cap[0]=((containsCaptureFlashed*1.0)/all)*100
# No-contain capture Flashed
r4_Cont_cap[1]=(((REflashed-containsCaptureFlashed)*1.0)/all)*100
# Contains capture Nonflashed
r4_Cont_cap[2]=((containsCaptureNonflashed*1.0)/all)*100
# No-contain capture Nonflashed
r4_Cont_cap[3]=(((allnonflashed-containsCaptureNonflashed)*1.0)/all)*100

# Contains cap+rep Flashed
r5_Cap_and_rep[0]=((containsCapAndRepFlashed*1.0)/all)*100
# Contains cap+excl Flashed
r5_Cap_and_rep[1]=0
# Contains only cap Flashed
r5_Cap_and_rep[2]=(((containsCaptureFlashed-containsCapAndRepFlashed)*1.0)/all)*100
# Contains cap+excl Nonflashed
r5_Cap_and_rep[3]=((containsCapAndRepNonflashed*1.0)/all)*100
# Contains cap+excl Nonflashed
r5_Cap_and_rep[4]=0
# Contains only cap Nonflashed
r5_Cap_and_rep[5]=(((containsCaptureNonflashed-containsCapAndRepNonflashed)*1.0)/all)*100

# Single cap Flashed
r6_multicap[0]=((singleCapFlashed*1.0)/all)*100
# Multicap Flashed
r6_multicap[1]=((multiCapFlashed*1.0)/all)*100
# Single cap Nonflashed
r6_multicap[2]=((singleCapNonflashed*1.0)/all)*100
# Multicap Nonflashed
r6_multicap[3]=((multiCapNonflashed*1.0)/all)*100

# NonDuplicates Flashed
r7_duplicate[0]=((nonduplicateFlashed*1.0)/all)*100
# Duplicates Flashed
r7_duplicate[1]=(((singleCapFlashed-nonduplicateFlashed)*1.0)/all)*100
# NonDuplicates Nonflashed
r7_duplicate[2]=((nonduplicateNonflashed*1.0)/all)*100
# Duplicates Nonflashed
r7_duplicate[3]=(((singleCapNonflashed-nonduplicateNonflashed)*1.0)/all)*100

# Platploidy-nonfiltered Flashed
r8_blatploidy[0]=(((nonduplicateFlashed-blatploidyFlashed)*1.0)/all)*100
# Platploidy-filtered Flashed
r8_blatploidy[1]=((blatploidyFlashed*1.0)/all)*100
# Platploidy-nonfiltered Nonflashed
r8_blatploidy[2]=(((nonduplicateNonflashed-blatploidyNonflashed)*1.0)/all)*100
# Platploidy-filtered Nonflashed
r8_blatploidy[3]=((blatploidyNonflashed*1.0)/all)*100


print "r0_All_reads",r0_All_reads
print "r1_flashed",
printList(r1_flashed)
print "r2_RE_site",
printList(r2_RE_site)
print "r3_continues",
printList(r3_continues)
print "r4_Cont_cap",
printList(r4_Cont_cap)
print "r5_Cap_and_rep",
printList(r5_Cap_and_rep)
print "r6_multicap",
printList(r6_multicap)
print "r7_duplicate",
printList(r7_duplicate)
print "r8_blatploidy",
printList(r8_blatploidy)

