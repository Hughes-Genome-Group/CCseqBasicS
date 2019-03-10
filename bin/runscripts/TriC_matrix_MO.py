# -*- coding: utf-8 -*-

##########################################################################
# Copyright 2019, Marieke O Oudelaar (marieke.oudelaar@univ.ox.ac.uk)
# portability fixes Jelena Telenius (jelena.telenius@imm.ox.ac.uk)       #
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


import numpy as np
import matplotlib

# To allow using cluster "non-interactive plotting" - so not to complain when X-windows is not present.
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import re as re
from collections import defaultdict
plt.rcParams['pdf.fonttype'] = 42

# Command line arguments, help function
import getopt

# Reading input params + less error prone file path parsing (instead of manual parsing) ..
import os

# The stuff to bring in clock time to report 'the script was ran at this time'
import datetime
wallclock = datetime.datetime.now()
date_and_time = wallclock.strftime("%Y-%m-%d %H:%M")

# To plot integer-based grid, we need square root, as we are turning a square 45 degrees,
# and thus have hypothenusa on the x-axis
from math import sqrt

# The stuff to rotate and clip the printed matrix
import matplotlib.transforms as mtransforms
# The stuff to scale the colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

# To be able to ask which version of python we are running ..
import sys

# The rotation and masking subroutine fetched 7March2019 from :
# https://matplotlib.org/gallery/images_contours_and_fields/affine_image.html#sphx-glr-gallery-images-contours-and-fields-affine-image-py
# And modified to the below form by Jelena Telenius
def do_plot(ax, Z, transform):
    # Select the subplot in the main figure (we do have only one -  but I tested with several, that's why)
    plt.sca(ax)
    # Plot the matrix, using 'viridis' color scale, and dimensions of sqrt(2) to allow nice transformation of 45 degrees
    im = ax.imshow(Z, interpolation='nearest',
                   origin='upper',vmin = 0.001, vmax = threshold, cmap=plt.cm.viridis,
                   extent=[-0, 2*sqrt(2), -0, 2*sqrt(2)], clip_on=True)
    
    # Set title
    title_text = "Tri-C matrix plot  \n(C) Marieke Oudelaar 2019\n\n" + my_sample + "\nplot printed " + date_and_time + "\n"
    heading = plt.title(title_text, loc='left')

    # Set x label
    label_for_plot = str(chrom) + ":" + str(start) + "-" + str(stop) + "\nBin size " + str(bin_size) + "b  :  Max value " + str(threshold) + " RPM/bin"   
    plt.xlabel(label_for_plot)

    # Turn the figure according to the 'transform' given as parameter
    trans_data = transform + ax.transData
    im.set_transform(trans_data)

    # display intended extent of the image
    x1, x2, y1, y2 = im.get_extent()
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "y  ",
            transform=trans_data)
    
    # Put the plot to match the area - we have diagonal lenght of 4 units, and
    # that makes our x-axis 4 units wide and y-axis 2 units wide (we want to hide the bottom half of the triangle)
    ax.set_xlim(-2,2)
    ax.set_ylim(2,4)
    
    # Turn off ticks, tick labels, figure frame
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_frame_on(False)
    ax.set_yticks([])
    ax.set_xticks([])
    
    # Make color bar to scale with the plot, not with the whole page
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    plt.colorbar(im, cax=cax)
    
    # Adjust marigins to look balanced
    plt.subplots_adjust(left=0.05, right=0.90, top=0.85, bottom=0.05)
    
def print_help(exitcode):
    
    print ("\nUsage :\n")
    
    print(os.path.basename(sys.argv[0]) + "\n")
          
    print ("-f --file \t OBLIGATORY \t Input file (relative or absolute path) ")
    print ("\t\t\t\t This is the   *_TriC_interactions-*.tab   output file of TriC_MO.pl script. \n")
           
    print ("-c --chr \t optional \t Chromosome (only for printing out the name in plot title)")
    print ("-l --str \t OBLIGATORY \t Start position (left)  of the visualisation plot")
    print ("-r --stp \t OBLIGATORY \t Stop  position (right) of the visualisation plot\n")
    
    print ("-b --bin \t optional \t Bin width (default 1000 bases)")
    print ("-t --threshold \t optional \t Max value of normalised interactions for a bin - values above this will be capped to value [threshold] (default 20 RPM/bin)")
    print ("\t\t\t\t The RPM bin-wise normalisation follows formula (for interactions between bins i,j) :")
    print ("\t\t\t\t (10^6/all_interactions_of_whole_matrix)*(sum_of_interactions_for_bin_ij/(REfragments_in_i+REfragments_in_j)\n")
    
    print ("-s --sample \t optional \t Sample name (defaults to basename of input file)")
    print ("-o --outdir \t optional \t Output folder (defaults to run folder)\n")
 
    print ("-h --help \t \t \t Prints this help.\n")
    
    exit(exitcode)
    
###############################################################################################################################################################

# Hardcoded parameters :
email = 'marieke.oudelaar@univ.ox.ac.uk';
version = "1.0.0";

###############################################################################################################################################################

print ( "\n" )
print ( '------------------------------------------------' )
print ( "TriC analyser matrix plotter - version " + version )
print ( '------------------------------------------------' )
print ( "Developer email " + email )
print ( "\n" )

###############################################################################################################################################################

# Help user case

if len(sys.argv) == 2 :
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print_help(0)

###############################################################################################################################################################

print ( "Run started at " +  date_and_time + "\n\n" )

print ( "Running with Python/Matplotlib/Numpy versions :" )

# The setup user is running it
temp_versionline = sys.version
temp_versionline=temp_versionline.replace('\n', ' ').replace('\r', '')
print("Python " + temp_versionline)
temp_versionline = matplotlib.__version__
temp_versionline=temp_versionline.replace('\n', ' ').replace('\r', '')
temp_versionline2 = np.version.version
temp_versionline2=temp_versionline2.replace('\n', ' ').replace('\r', '')
print("Matplotlib " + temp_versionline + " Numpy " + temp_versionline2)

# The setup it was tested in
print ("\n( the script was developed and tested in 2 setups :")
print (" Python 3.5.3 [GCC 4.4.7 20120313 (Red Hat 4.4.7-17)] Matplotlib 2.0.2 Numpy 1.15.4\n and ")
print (" Python 2.7.5 [GCC 4.4.7 20120313 (Red Hat 4.4.7-4) ] Matplotlib 2.2.2 Numpy 1.14.2 )\n\n")
    
###############################################################################################################################################################

# Read in paramters ..

# Default values ..

full_file_name_ery="/UNDEFINED/INPUT/FILE.txt" # This is full path to the input file (can be relative path i.e. FILE.txt or ../FILE.txt etc)
my_sample = ""
# this is the old "my_ery_file" : the BASENAME of the input file - without file ending if the file has one ("FILE" if parsing the above line). will be parsed from the input params
# if samplename is given in command line params, the above parsed value is not used, but rather the user-given name for the sample is used .

# my_dir = "/UNDEFINED/FOLDER/WHERE/THE/INPUT/FILE/CAN/BE/FOUND" # this is the DIRNAME of the input file.
# this is not needed any more.

my_dir_out = os.getcwd() # this is the output dir name. does not have to exist, will be generated on the fly, defaults to "here"

# specify coordinates and resoluton for plot

chrom = 'chr?'
start = -1
stop = -1
bin_size = 1000
threshold = 20

############################################

# Parameters in ..

# the following is based on tutorial in site :
# https://stackabuse.com/command-line-arguments-in-python/
# using getopt module (not using argparse module - to keep python2 support)

# read commandline arguments,
fullCmdArguments = sys.argv
# script name is the 0th argument, the actual params are 1-->
argumentsList = fullCmdArguments[1:]
# if no arguments given, crashing.
if not argumentsList:
    sys.stderr.write("ERROR : No parameters given.\nEXITING!\n")
    print_help(1)

# List valid parameters
unixOptions = "hc:vl:vr:vb:vs:vo:vf:vt:v"
gnuOptions = ["help", "chr=", "str=", "stp=", "bin=", "sample=", "outdir=", "file=", "threshold="]  

# Start parsing from command line ..
try:  
    arguments, values = getopt.getopt(argumentsList, unixOptions, gnuOptions)
except getopt.error as err:
    sys.stderr.write("ERROR : " + str(err))
    print_help(1)
    
# Now, going through the above-red stuff (if something was wrong in the parameter names , we already crashed above)

# Start - obligatory, has to be int
# Stop - obligatory, has to be int
# Bin - optional, but if given has to be int
# Sample - optional, but if given has to be string
# Outdir - optional, but if given has to be string
# Input file - obligatory (existence is tested later)

for currentArgument, currentValue in arguments:  
    if currentArgument in ("-h", "--help"):
        print_help(0)
    elif currentArgument in ('-c', '--chr'):
        if currentValue:
            chrom = currentValue
    elif currentArgument in ('-l', '--str'):
        if currentValue:
            start = int(currentValue)
    elif currentArgument in ('-r', '--stp'):
        if currentValue:
            stop = int(currentValue)
    elif currentArgument in ('-b', '--bin'):
        if currentValue:
            bin_size = int(currentValue)
    elif currentArgument in ('-s', '--sample'):
        if currentValue:
            my_sample = currentValue
    elif currentArgument in ('-o', '--outdir'):
        if currentValue:
            my_dir_out = currentValue
    elif currentArgument in ('-f', '--file'):
        if currentValue:
            full_file_name_ery = currentValue
    elif currentArgument in ('-t', '--threshold'):
        if currentValue:
            threshold = int(currentValue)

# Checking that obligatory ones were given ..

obligatory_fine = 1

if start == -1 :
    sys.stderr.write("ERROR : Start coordinate missing : has to be given with parameter -l or --str '\n")
    obligatory_fine = 0
if stop == -1 :
    sys.stderr.write("ERROR :  Stop coordinate missing : has to be given with parameter -r or --stp '\n")
    obligatory_fine = 0
if full_file_name_ery == "/UNDEFINED/INPUT/FILE.txt" :
    sys.stderr.write("ERROR : Input file missing : has to be given with parameter -f or --file '\n")
    obligatory_fine = 0    

if not obligatory_fine :
    sys.stderr.write("EXITING !! \n")
    print_help(1)
    
############################################

# For testing purposes - overwriting the above

# obligatory params
# full_file_name_ery="../triC_ery1_R2_TriC_interactions.txt" # string, no test needed, if wonky, gets crashed afterwards nicely.
# start = 32030000 # int, has to be int.
# stop = 32250000 # int, has to be int.

# optional params 
# my_sample = "R2" # no test needed, if empty the name gets parsed from full_file_name_ery
# my_dir_out = "figging" # string, no test needed - crazy whitespace containing paths should be caught in python's own subroutines
# chrom = 'chr11' # string, no test needed (is only needed for printing out purposes)
# bin_size = 1000 # int, has to be int.

###############################################################################################################################################################

# Printing out what we just set ..

print ("Starting run with parameters : \n")

print ("Genomic coordinates to plot : " + chrom + ":" + str(start) + "-" + str(stop) )
print ("Bin size : " + str(bin_size) + "b\t Treshold : max " + str(threshold) + "RPM/bin \n")

print ("Input file : '" + full_file_name_ery + "'\n")

###############################################################################################################################################################


# Setting the starting setup, with the input parameters ..

# specify file name and directories

# Ask if it exists and is a normal file (not a dir) otherwise parsing this doesn't make much sense ..
if not os.path.isfile(full_file_name_ery) :
   sys.stderr.write("Input file not found : '" + full_file_name_ery + "'\n")
   exit(1)

# specify file name and directories
    
# if my_sample is still empty, parsing full_file_name_ery to fill it
if not my_sample:
    filename, file_extension = os.path.splitext(os.path.basename(full_file_name_ery))
    my_sample=filename


# make the output dir

print ("\nMaking output directory .. \n")


if not os.path.isdir(my_dir_out):
    os.makedirs(my_dir_out)
    print ("(Output directory '"  + my_dir_out + "' created)")
else:    
    print ("(Existing output directory '" + my_dir_out + "' found)")

# Note that the above does NOT specifically catch the situation,
# where the path exists, but points to a FILE not a FOLDER.
# Whatever happens in that case, is up to Python's basic mode in these situations
# Could be no problem (as some OS's allow paths and folders with same name)


# generate the output file name start using the above info
full_file_name_out_ery = os.path.join(my_dir_out, my_sample)

###############################################################################################################################################################

# Printing out what we just set ..

print ("\nOutput directory set : '" + my_dir_out )
print ("\nOutput file names set to start with : '" + full_file_name_out_ery + "*'\n\n")

###############################################################################################################################################################

# initialise matrix
print ("Initialise matrix ..")

bin_start = start / bin_size
bin_stop = stop / bin_size
n_bins = int(bin_stop - bin_start)
matrix = np.zeros((n_bins, n_bins))

# generate matrix
print ("Generate matrix ..")

matrix_interaction_count = 0
int_dic = {}
RF_dic = defaultdict(list)

with open(full_file_name_ery) as f_ery:
    for line in f_ery:
        RF1, RF2, count = line.split()
        chr1, start_coord1, stop_coord1 = re.split(":|-", RF1)
        chr2, start_coord2, stop_coord2 = re.split(":|-", RF2)
        mid1 = (int(start_coord1) + int(stop_coord1)) / 2
        mid2 = (int(start_coord2) + int(stop_coord2)) / 2
        if mid1 >= start and mid1 < stop and mid2 >= start and mid2 < stop:
            matrix_interaction_count += int(count)      # count number on interactions contributing to matrix
            bin1 = int(mid1 / bin_size)
            bin2 = int(mid2 / bin_size)
            int_dic[mid1, mid2] = int(count)
            if mid1 not in RF_dic[bin1]:                # count number of detected restriction fragmnets per bin for normalisation
                RF_dic[bin1].append(mid1)
            if mid2 not in RF_dic[bin2]:
                RF_dic[bin2].append(mid2)

for key in int_dic:
    bin1 = int(key[0] / bin_size)
    bin2 = int(key[1] / bin_size)
    if bin1 >= bin_start and bin1 < bin_stop and bin2 >= bin_start and bin2 < bin_stop:
        corr = len(RF_dic[bin1]) * len(RF_dic[bin2])    # calculate number of restriction fragments contributing to each bin 
        # normalise for number of restriction fragments per bin and total number of interactions contributing to matrix
        matrix[int(bin1 - bin_start), int(bin2 - bin_start)] +=  1.0*int_dic[key] / corr * 1000000 / matrix_interaction_count
        

print ("Generate plot ..")

# make the subplot object(s) and fig object to call later
fig, ax = plt.subplots(1,1)
# fig, axs = plt.subplots(2,1)

# Ask for the transform and plotting - 45 degrees rotation as parameter
do_plot(ax, matrix, mtransforms.Affine2D().rotate_deg(45))

# For testing purposes - just show the plot, don't print it to a file
# plt.show()


print ("Save figure ..")

# Save file as pdf
plt.savefig(full_file_name_out_ery + "_" + str(bin_size) + "_" + str(threshold) + ".pdf", format='pdf', dpi=1000)

# Print the values of the normalised matrix (if people want to visualise it the same way)

print ("Save normalised matrix (tab-separated) ..")

# Raw data, without threshold :
np.savetxt(full_file_name_out_ery + "_" + str(bin_size) + "_RAW.tab", matrix, delimiter='\t', fmt='%2.2f')

# Raw data, clipped to the threshold used in the figure plotting :
np.clip(matrix, 0, 20)
np.savetxt(full_file_name_out_ery + "_" + str(bin_size) + "_" + str(threshold) + ".tab", matrix, delimiter='\t', fmt='%2.2f')

# Print that we finished ..
date_and_time = wallclock.strftime("%Y-%m-%d %H:%M")
print ( "\n\nRun completed at " +  date_and_time + "\n\n" )




