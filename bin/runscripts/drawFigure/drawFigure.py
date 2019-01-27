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
# print "Preparing to run - check where we are and which version we have.."
import sys
# print ""
print "We are running in machine :"
print sys.platform
print "We are running in Python version :"
print sys.version
print ""
# print "We can load from these paths :"
# print sys.path
# print "We have these auto-loaded modules"
# print sys.modules
# print "We have these built-ins :"
# print dir(__builtins__)

# print "----------------------------------------"
# print "Run directory :"
import os
# print os.getcwd()
# print "----------------------------------------"
# print "Enabling log file output.."
import syslog
# print "Enabling resource use and run time statistics monitoring.."
import stat
# print "----------------------------------------"
# print "Importing script-specific libraries.."
# print "----------------------------------------"
# 
# print "Importing regular expressions"
import re
# 
# print "Importing matplotlib"
import matplotlib as mpl
print "We are running matplotlib version :"
print mpl.__version__
# 
# print "Available back ends (instead of X windows) :"
# print (mpl.rcsetup.non_interactive_bk)
# print "Now loading the back end : "
# print "mpl.use('pdf')"
mpl.use('pdf')

# print "Importing pyplot "
import matplotlib.pyplot as plt
# print "Importing patches "
import matplotlib.patches as patch

# print "----------------------------------------"
print "Imported (and auto-loaded) modules :"
print(globals())

# print "----------------------------------------"
# print "Reading in the subroutines.."
# print >> sys.stderr, "----------------------------------------"
# print >> sys.stderr, "Reading in the subroutines.."

# Making the comments above the lines..
def writeComment(row) :
    layer1.annotate(myComment[row], (2, 88-(10*row))
    )

# Subroutines to make RED-GREEN drawing :
# row = from the top, what is the row number for this data ?
def drawTwoColorsFlashed(row) :
    print >> sys.stderr, 'Flashed'
    print >> sys.stderr, row
    # myPercentages[1]=[70,30]
    layer1.broken_barh(
                   [(0, myFlashedPercentages[row][0]), (myFlashedPercentages[row][0], myFlashedPercentages[row][1])], # X (start, width)
                   (83-(10*row), 4),             # Y (start, height) 
                   facecolors=myColors[row]
    )

def drawTwoColorsNonFlashed(row) :
    print >> sys.stderr, 'NonFlashed'
    print >> sys.stderr, row
    # myPercentages[1]=[70,30]
    layer1.broken_barh(
                   [(myPercentages[1][0], myNonFlashedPercentages[row][0]), (myPercentages[1][0]+myNonFlashedPercentages[row][0], myNonFlashedPercentages[row][1])], # X (start, width)
                   (83-(10*row), 4),             # Y (start, height) 
                   facecolors=myColors[row]
    )

# Subroutines to make RED-ORANGE-GREEN drawing :
# row = from the top, what is the row number for this data ?
def drawThreeColorsFlashed(row) :
    print >> sys.stderr, 'Flashed'
    print >> sys.stderr, row
    # myPercentages[1]=[70,30]
    layer1.broken_barh(
                   [(0, myFlashedPercentages[row][0]), (myFlashedPercentages[row][0], myFlashedPercentages[row][1]), (myFlashedPercentages[row][0]+myFlashedPercentages[row][1], myFlashedPercentages[row][2])], # X (start, width)
                   (83-(10*row), 4),             # Y (start, height) 
                   facecolors=myColors[row]
    )

def drawThreeColorsNonFlashed(row) :
    print >> sys.stderr, 'NonFlashed'
    print >> sys.stderr, row
    # myPercentages[1]=[70,30]
    layer1.broken_barh(
                   [(myPercentages[1][0], myNonFlashedPercentages[row][0]), (myPercentages[1][0]+myNonFlashedPercentages[row][0], myNonFlashedPercentages[row][1]), (myPercentages[1][0]+myNonFlashedPercentages[row][0]+myNonFlashedPercentages[row][1], myNonFlashedPercentages[row][2])], # X (start, width)
                   (83-(10*row), 4),             # Y (start, height) 
                   facecolors=myColors[row]
    )

# Subroutines to make GREEN-GREEN drawing :
# row = from the top, what is the row number for this data ?
def drawOneColorFlashed(row) :
    print >> sys.stderr, 'Flashed'
    print >> sys.stderr, row
    # myPercentages[1]=[70,30]
    layer1.broken_barh(
                   [(0, myFlashedPercentages[row])], # X (start, width)
                   (83-(10*row), 4),             # Y (start, height) 
                   facecolors=myColors[row]
    )

def drawOneColorNonFlashed(row) :
    print >> sys.stderr, 'NonFlashed'
    print >> sys.stderr, row
    # myPercentages[1]=[70,30]
    layer1.broken_barh(
                   [(myPercentages[1][0], myNonFlashedPercentages[row])], # X (start, width)
                   (83-(10*row), 4),             # Y (start, height) 
                   facecolors=myColors[row]
    )

# print "----------------------------------------"
# print "Starting the run.."
# print "----------------------------------------" 
# print ""
print "Reading the input.."
# print ""
# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Reading the input.."
# print >> sys.stderr, ""

names=[]
values = []
colors = []
valuesAsFloats = []

with open('percentages.txt') as f:
    for line in f:
        data = line.split()
        names.append(re.sub(r'_', ' ', data[0]))
        values.append(data[1:])
        temp = []
        for i, value in enumerate(data[1:]):
            temp.append(float(value))
        valuesAsFloats.append(temp)

# print "names :"
# print names
# print "values :"
# print values
# print "valuesAsFloats :"
# print valuesAsFloats
# print "valuesAsFloats[0] :"
# print valuesAsFloats[0]
# print "valuesAsFloats[1] :"
# print valuesAsFloats[1]
# print "valuesAsFloats[2] :"
# print valuesAsFloats[2]
# print "valuesAsFloats[3] :"
# print valuesAsFloats[3]
# print "valuesAsFloats[4] :"
# print valuesAsFloats[4]


# print "----------------------------------------" 
# print ""
print "Setting values.."
# print ""
# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Setting values.."
# print >> sys.stderr, ""

# Generating the lists..
myLabel=['0','1','2','3','4','5','6','7','8']
myComment=['0','1','2','3','4','5','6','7','8']
myPercentages=[0,1]
myFlashedPercentages=[0,1,2,3,4,5,6,7,8]
myNonFlashedPercentages=[0,1,2,3,4,5,6,7,8]
myColors=['0',['1','1'],['2','2'],['3','3'],['4','4'],['5','5','5'],['6','6'],['7','7'],['8','8']]

# Default colors (for color blindness support)
#
# PINK GREEN (default)
#         RGB          HEX  
# red     255,74,179   #FF4ABE
# orange  255,140,0    #FF8C00
# green    62,176,145  #3EB091
#
red='#FF4ABE'
orange='#FF8C00'
green='#3EB091'

# Setting the values.. (most of them have four values - those are set here.)
for x in range(2, 9):
    if (x != 3 and x !=5 ):
        myFlashedPercentages[x]=[valuesAsFloats[x][0],valuesAsFloats[x][1]]
        myNonFlashedPercentages[x]=[valuesAsFloats[x][2],valuesAsFloats[x][3]]

myLabel[0]='Total reads (input fastq)'
myComment[0]='Total reads (input fastq)'
myPercentages[0]=valuesAsFloats[0][0]
myColors[0]='blue'

myLabel[1]='Flashed / nonflashed'
myComment[1]='Flash-combined (light blue), non-combined (yellow)'
myPercentages[1]=[valuesAsFloats[1][0],valuesAsFloats[1][1]]
myColors[1]=['dodgerblue','gold']

myLabel[2]='Do/don\'t have RE site'
myComment[2]='With RE site(green), no RE site(red)'
myColors[2]=[green,red]

myLabel[3]='Continue to mapping'
myComment[3]='Continues to mapping :'
myFlashedPercentages[3]=valuesAsFloats[3][0]
myNonFlashedPercentages[3]=valuesAsFloats[3][1]
myColors[3]=green

myLabel[4]='Contains capture'
myComment[4]='cap(green), nocap(red)'
myColors[4]=[green,red]

myLabel[5]='Capture and/or reporter'
# myComment[5]='cap+rep(green), cap+excl(orange), only cap(red)'
myComment[5]='cap+rep(green),only cap(red) - cap+excl also red'
myFlashedPercentages[5]=[valuesAsFloats[5][0],valuesAsFloats[5][1],valuesAsFloats[5][2]]
myNonFlashedPercentages[5]=[valuesAsFloats[5][3],valuesAsFloats[5][4],valuesAsFloats[5][5]]
myColors[5]=[green,orange,red]

myLabel[6]='Multiple (different) captures'
myComment[6]='single cap (green), multicap (red)'
myColors[6]=[green,red]

myLabel[7]='Duplicate filtered'
myComment[7]='non-duplicate (green), duplicate (red)'
myColors[7]=[green,red]

myLabel[8]='Blat/ploidy filtered'
myComment[8]='no-blat-no-ploidy (green), blat and/or ploidy (red)'
myColors[8]=[green,red]

# print >> sys.stderr,"----------------------------------------" 
# print >> sys.stderr,""
# print >> sys.stderr,"Checking that the labels are not in wonky order :"
# print >> sys.stderr,""

# for x in range(0, 9):
#     print >> sys.stderr,"Label here :", myLabel[x]
#     print >> sys.stderr,"Same line in input : ", names[x]

# print >> sys.stderr,""
# for x in range(0, 2):
#     print >> sys.stderr,"Label here :", myLabel[x]
#     print >> sys.stderr,"myPercentages : ", myPercentages[x]

# print >> sys.stderr,""
# for x in range(2, 9):
#     print >> sys.stderr,"Label here :", myLabel[x]
#     print >> sys.stderr,"myFlashedPercentages : ", myFlashedPercentages[x]
#     print >> sys.stderr,"myNonFlashedPercentages : ", myNonFlashedPercentages[x]


# print "----------------------------------------" 
# print ""
print "Drawing axes and tick marks (general overlay).."
# print ""

# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Drawing axes and tick marks (general overlay).."
# print >> sys.stderr, ""

# class matplotlib.figure.Figure(figsize=None, dpi=None, facecolor=None, edgecolor=None, linewidth=0.0, frameon=None, subplotpars=None, tight_layout=None)
# matplotlib.pyplot.subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None, **fig_kw)

fig1, layer1 = plt.subplots()

# Set the overall settings here ..

# Grid on (dotted lines)
layer1.grid(True)

# Where (in whole canvas) we want to put our y-range and x-range
# 0,0 is as normally, left hand down.

# Set x-axis to be from 0 to 100
layer1.set_xlim(0, 100)
layer1.set_xticks([ 0,10,20,30,40,50,60,70,80,90,100])
layer1.set_xlabel('Percentage of input reads')

# Set y-axis to be contain all the reads..
layer1.set_ylim(0, 100)

# From bottom up (as the coordinates go that direction) :
# Copy and reverse the list..
myReverseLabels=myLabel[:]
myReverseLabels.reverse()

layer1.set_yticks([5,15,25,35,45,55,65,75,85])
layer1.set_yticklabels(myReverseLabels)

# print "----------------------------------------" 
# print ""
print "Drawing boxes and their labels.."
# print ""

# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Drawing boxes and their labels.."
# print >> sys.stderr, ""

# matplotlib.pyplot.broken_barh(xranges, yrange, hold=None, data=None, **kwargs)
# Plot horizontal bars.

myFlashedPercentages[1]=myPercentages[1]
layer1.broken_barh(
               [(0, myFlashedPercentages[1][0]), (myFlashedPercentages[1][0], myFlashedPercentages[1][1])], # X (start, width)
               (0, 75),             # Y (start, height) 
               facecolors=['lightcyan','lemonchiffon'], edgecolor = "none"
)

# Total reads (input fastq)
writeComment(0)
myFlashedPercentages[0]=myPercentages[0]
drawOneColorFlashed(0)

# Flashed / nonflashed
writeComment(1)
myFlashedPercentages[1]=myPercentages[1]
drawTwoColorsFlashed(1)

# Do/don\'t have RE site
writeComment(2)
drawTwoColorsFlashed(2)
drawTwoColorsNonFlashed(2)

# Continue to mapping
writeComment(3)
drawOneColorFlashed(3)
drawOneColorNonFlashed(3)

# Contains capture
writeComment(4)
drawTwoColorsFlashed(4)
drawTwoColorsNonFlashed(4)

# Capture and/or reporter
writeComment(5)
drawThreeColorsFlashed(5)
drawThreeColorsNonFlashed(5)

# Multiple (different) captures
writeComment(6)
drawTwoColorsFlashed(6)
drawTwoColorsNonFlashed(6)

# Duplicate filtered
writeComment(7)
drawTwoColorsFlashed(7)
drawTwoColorsNonFlashed(7)

# Blat/ploidy filtered
writeComment(8)
drawTwoColorsFlashed(8)
drawTwoColorsNonFlashed(8)

# print "----------------------------------------" 
# print ""
print "Saving figure.."
# print ""
# print >> sys.stderr, "----------------------------------------" 
# print >> sys.stderr, ""
print >> sys.stderr, "Saving figure.."
# print >> sys.stderr, ""

fig1.savefig('summary.pdf', dpi=90, bbox_inches='tight')
fig1.savefig('summary.png', dpi=90, bbox_inches='tight')







