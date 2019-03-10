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

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

######################################################################################################################################################

# This script performs the analysis of a Tri-C experiment. It depends on the combined bam file outputted by the CCseqBasic pipeline "COMBINED_reported_capture_reads_CS5.bam" (F6 folder), ran using an oligo file without proximity exclusion (see instructions here: http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/PipeSite.html), and converted to sam format. This files contains all mapped reads containing 1 capture fragment post PCR duplicate, blat and ploidy filtering. The scripts also requires the oligo file and a file with the coordinates of the restriction fragments (of correspondong enzyme) in the genome. 

# The script selects the reporter reads in cis, performs a proximity exclusion, maps the reads to the restriction fragments and counts multi-way interactions. These are outputted in a text file with suffix "_TriC_interactions.txt". This file can be used as input for the TriC_matrix.py script to generate a basic interaction matrix.

# Output:
# 1. Report with read and interaction counts
# 2. Text file containing reported multi-way interactions per viewpoint in format: bin1 \t bin2 \t count \n

# Example of run command:
# nohup perl /t1-data/user/hugheslab/oudelaar/scripts/TriC/TriC_MO_bin.pl -sam /t1-data/user/hugheslab/oudelaar/TriC_KO_CTCF/samfiles_pipe/bowtie1/C57-1A_Tri_KO_CS5.sam -o /t1-data/user/hugheslab/oudelaar/Tri-C/pipe/tri-c_oligo_file_noprex.txt -r /t1-data1/WTSA_Dev/oudelaar/scripts/mm9_nlaIII_coordinates.txt -name C57-1A -b 1000 &

######################################################################################################################################################

# Hardcoded parameters :
my $email = 'marieke.oudelaar@univ.ox.ac.uk';
my $version = "1.0.0";

# Obligatory parameters :
my $oligo_file = "UNDEFINED";
my $sam_file = "UNDEFINED";

my $name = "UNDEFINED";

my $ucscsizes="UNDEFINED";
my $genome = "UNDEFINED";

my $dig_genome="UNDEFINED"; #Specifies filename of restriction coordinates (file made by dpngenome2.pl)

# Optional parameters :
my $store_bigwigs_here_folder="UNDEFINED";

# Parameters with hardcoded default value :
my $public_folder = "DATA_FOR_PUBLIC_FOLDER";
my $server_folder = "UNDEFINED";
my $server = 'http://sara.molbiol.ox.ac.uk';
my $use_symlinks=0; #whether to use symlinks to store bigwigs or not
my $bin_size=1000;

print STDOUT "\n" ;
print STDOUT "TriC analyser using a fixed BIN size - version $version !\n" ;
print STDOUT "Developer email $email\n" ;
print STDOUT "\n" ;

&GetOptions
(
    "sam=s" =>\ $sam_file, 	    # -sam          sam file
    "o=s" => \ $oligo_file,       # -o            oligo file (CCanalyser format)
    "c=s" => \ $oligo_file,       # -c            oligo file (CCanalyser format) - compatibility flag : CCseqBasic nowadays talks about "capture coordinate file" rather than "oligo coordinate file"
    "r=s" =>\ $dig_genome,       # -r            file with restriction coordinates in genome 
    "name=s" =>\ $name,          # -name         name of the experiment, can be whatever you like, will be used for names of output files
    "pf=s" =>\ $public_folder,   # -pf	    your public folder (e.g. /public/username)
    "sf=s" =>\ $server_folder,   # -pf	    your public folder (e.g. /public/username - how it is seen WITHIN the web address http://sara.molbiol.ox.ac.uk : in WIMM system the same as -pf)
    "server=s"=>\ $server,	    # -pu	    Your server setup (e.g. http://sara.molbiol.ox.ac.uk)
    "genome=s"=>\ $genome,	    # -genome	    Specify the genome (mm9 / hg18)
    "ucscsizes=s"=>\ $ucscsizes, # -ucscsizes    Genome sizes file 
    "symlinks"=>\ $use_symlinks,    # -symlinks     To make symlinks to bigwigs into the public area, instead of actually storing the bigwigs there.
    "b=i" => $bin_size,         # -b            bin size
);

# If symlinks were requested, setting stuff ..
if ($use_symlinks)
{
$store_bigwigs_here_folder="PERMANENT_BIGWIGS_do_not_move";
}
else
{
$store_bigwigs_here_folder="$public_folder/$name\_TriC/";
}

# If --sf is not given, assuming intra-WIMM user case (pf and sf being the same)
if ( $server_folder eq "UNDEFINED" )
{
$server_folder=$public_folder
}



# Any case updating public folders
$public_folder="$public_folder/$name\_TriC";
$server_folder="$server_folder/$name\_TriC";

# Printing out the parameters in the beginning of run

print STDOUT "Starting run with parameters :\n" ;
print STDOUT "\n" ;

print STDOUT "capturesite filename $oligo_file\n";
print STDOUT "sam file $sam_file\n";
print STDOUT "restriction enzyme coords file $dig_genome \n";
print STDOUT "sample $name\n" ;
print STDOUT "sample $bin_size\n" ;

print STDOUT "store_bigwigs_here_folder $store_bigwigs_here_folder \n" ;
print STDOUT "public_folder $public_folder \n" ;
print STDOUT "server_folder $server_folder \n" ;
print STDOUT "server $server \n";
print STDOUT "ucscsizes $ucscsizes\n" ;
print STDOUT "genome $genome\n" ;
print STDOUT "email $email\n" ;

# Crash if some important ones are not given ..

my $parameters_missing=0;
if ( $oligo_file eq "UNDEFINED" ){$parameters_missing=1}
if (   $sam_file eq "UNDEFINED" ){$parameters_missing=1}
if ( $dig_genome eq "UNDEFINED" ){$parameters_missing=1}
if ($parameters_missing) {
    die "Incomplete run command : you need to give all these parameters :\n --sam input.sam -o/-c capturesites.txt -r REcoordinates.txt\n\nAborting triC run !\n\n";
}

if ( $name eq "UNDEFINED" ){$parameters_missing=1}
if ($parameters_missing) {
    print STDERR "WARNING : triC run does not have sample name. Using default name 'UNDEFINED'\n\n";
    print STDOUT "WARNING : triC run does not have sample name. Using default name 'UNDEFINED'\n\n"
}

my $public_parameters_missing=0;
if ( $ucscsizes eq "UNDEFINED" ){$public_parameters_missing=1}
if ($public_parameters_missing) {
    print STDERR "WARNING : triC run does not have UCSC genome sizes. Bigwig tracks will not be generated.\n\n";
    print STDOUT "WARNING : triC run does not have UCSC genome sizes. Bigwig tracks will not be generated.\n\n"
}

if ( $public_folder eq "UNDEFINED" ){$public_parameters_missing=1}
if (        $genome eq "UNDEFINED" ){$public_parameters_missing=1}
if ($public_parameters_missing) {
    print STDERR "WARNING : triC run does not have complete public server area setup. Won't be able to generate a functional UCSC data hub.\n\n";
    print STDOUT "WARNING : triC run does not have complete public server area setup. Won't be able to generate a functional UCSC data hub.\n\n"
}

# Open filehandles

open (SAMFH, $sam_file) or die "can't open sam file $sam_file";
open (OLIFH, $oligo_file) or die "can't open oligo file $oligo_file";
open (GENFH, $dig_genome) or die "can't open genome file with restriction coordinates $dig_genome";

my $path = "undefined"; 
if ($sam_file =~ /(.*\/)(\V++)/) {
    $path = $1;
    }
my $dir = "$path/$name\_TriC_binned/";
unless (-d $dir) {
    mkdir $dir;
}

my $report = "$dir/$name\_TriC_analysis_report.txt";
open (REP, ">$report") or die "can't open output report file $report";

# Store oligo coordinates in oligos hash

my %oligos;
while (my $line = <OLIFH>) {
    chomp $line;
    my ($id, $chr1, $start1, $stop1, $chr2, $start2, $stop2, $SNP_loc, $SNP_base) = split (/\t/, $line);
    $oligos{$id}{"chr"}=$chr1;
    $oligos{$id}{"start"}=$start1;
    $oligos{$id}{"stop"}=$stop1;
}

# Store restriction fragment coordinates in RE_hash and sort in ascending order; to use for binary search

my %RE_hash;
while (my $line = <GENFH>) {
    chomp $line;
    my ($chr, $start, $stop) = split (/\W/, $line);
    push @{$RE_hash{$chr}}, $start;
}

my @chr_index = keys %RE_hash;
foreach my $chr (@chr_index) {
    @{$RE_hash{$chr}} = sort {$a <=> $b} @{$RE_hash{$chr}};
}

# Map reporter reads in cis to restriction fragments, remove fragments in proximity to the viewpoint, and store in data_hash; both as keys to ensure restriction fragments are only reported once per read, and in array to keep track of order and facilitate multi-way interaction analysis

my %data_hash; my %counter; 
while (my $line = <SAMFH>) {
    chomp $line;
    unless ($line =~ /^@/) {
        $counter{"01 Number of data lines in sam file"}++;
        my ($name, $flag, $chr, $start, $map_qual, $cigar, $mate_ref, $mate_start, $mate_insert, $seq, $seq_qual, $options) = split (/\t/, $line, 12);
        my $read_name = "undefined";
        if ($name =~ /(.*):PE(.*)$/) {
            $read_name = $1;
        }
        $chr =~ s/chr//; 
        my $stop = 0;
        if ($cigar =~/(\d++)(.*)/) {
            $stop = $1 + $start;                # Use length of mapped read until first indel for proximity exclusion and mapping to restriction fragment
        }
        if ($options =~ /(.*)CO:Z:(.*)(_CISREP)/) {             # Select reporters in cis
            my $oligo = $2;    
            my ($start_frag, $end_frag) = binary_search(\@{$RE_hash{$chr}}, $start, $stop, \%counter);              # Map cis reporters to restriction fragment
            unless ($start_frag =~ /error/) {
                if ($stop <= $oligos{$oligo}{"start"} - 1000 or $start >= $oligos{$oligo}{"stop"} + 1000) {             # Proximity exclusion
                    unless (exists $data_hash{$oligo}{$read_name}{"reporters"}{"$chr:$start_frag-$end_frag"}) {         # Check mapped restriction fragments are only reported once per read
                        push (@{$data_hash{$oligo}{$read_name}{"reporter_array"}}, "$chr:$start_frag-$end_frag");
                        $data_hash{$oligo}{$read_name}{"reporters"}{"$chr:$start_frag-$end_frag"} = 1;
                        $counter{"02 Total number of unique reporters in cis for $oligo post proximity exclusion"}++;
                        }
                        my $mid_frag = ($end_frag + $start_frag) / 2;
                        my $bin_start = int($mid_frag / $bin_size) * $bin_size;
                        #unless (exists $data_hash{$oligo}{$read_name}{"bins"}{"$bin_start"}) {     # Check bins are only reported once per read
                            push (@{$data_hash{$oligo}{$read_name}{"bin_array"}}, "$bin_start");
                            #$data_hash{$oligo}{$read_name}{"bins"}{"$bin_start"} = 1;
                        #}
                    }
                else {
                      $counter{"03 Number of proximity excluded fragments for $oligo"}++;
                }
            }
        }
    }
}

# Count multi-way interactions in all reads stored in data_hash and store these counts in frag_hash

my %frag_hash; 
foreach my $oligo (keys %data_hash) {
    foreach my $read_name (keys %{$data_hash{$oligo}}) {
        if (exists $data_hash{$oligo}{$read_name}{"reporter_array"}) {
            my @frags_sorted = sort @{$data_hash{$oligo}{$read_name}{"reporter_array"}};
            if ($#frags_sorted == 0) {                                              # Duplets
                $counter{"04a Number of duplets in cis with $oligo"}++;
            }
            if ($#frags_sorted == 1) {                                              # Triplets
                $counter{"04b Number of triplets in cis with $oligo"}++;
                $counter{"04 Total number of multi-way interactions in cis with $oligo"}++;
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]}) {         
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} += 1;                                
                }
            }
            if ($#frags_sorted == 2) {                                              # Quadruplets 
                $counter{"04c Number of quadruplets in cis with $oligo"}++;
                $counter{"04 Total number of multi-way interactions in cis with $oligo"} += 3;
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]}) {         
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]}) {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]} += 1;
                }
            }
            if ($#frags_sorted == 3) {                                              # Quintuplets
                $counter{"04d Number of quintuplets in cis with $oligo"}++;
                $counter{"04 Total number of multi-way interactions in cis with $oligo"} += 6;
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[3]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[3]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[3]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]}) {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[3]}) {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[3]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[3]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[3]}) {
                    $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[3]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[3]} += 1;
                }
            }
            if ($#frags_sorted == 4) {                                              # Sextuplets
                $counter{"04e Number of sextuplets in cis with $oligo"}++;
                $counter{"04 Total number of multi-way interactions in cis with $oligo"} += 10;
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[1]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[2]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[3]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[3]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[3]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[4]}) {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[4]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[0]}{$frags_sorted[4]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]}) {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[2]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[3]}) {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[3]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[3]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[4]}) {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[4]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[1]}{$frags_sorted[4]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[3]}) {
                    $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[3]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[3]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[4]}) {
                    $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[4]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[2]}{$frags_sorted[4]} += 1;
                }
                unless (exists $frag_hash{$oligo}{$frags_sorted[3]}{$frags_sorted[4]}) {
                    $frag_hash{$oligo}{$frags_sorted[3]}{$frags_sorted[4]} = 1;
                }
                else {
                    $frag_hash{$oligo}{$frags_sorted[3]}{$frags_sorted[4]} += 1;
                }
            }
        }
    }
}

my %bin_hash; 
foreach my $oligo (keys %data_hash) {
    foreach my $read_name (keys %{$data_hash{$oligo}}) {
        if (exists $data_hash{$oligo}{$read_name}{"bin_array"}) {
            my @bins_sorted = sort @{$data_hash{$oligo}{$read_name}{"bin_array"}};
            if ($#bins_sorted == 0) {                                              # Duplets
                $counter{"05a Number of binned duplets in cis with $oligo"}++;
            }
            if ($#bins_sorted == 1) {                                              # Triplets
                $counter{"05b Number of binned triplets in cis with $oligo"}++;
                $counter{"05 Total number of binned multi-way interactions in cis with $oligo"}++;
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]}) {         
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} += 1;                                
                }
            }
            if ($#bins_sorted == 2) {                                              # Quadruplets 
                $counter{"05c Number of binned quadruplets in cis with $oligo"}++;
                $counter{"05 Total number of binned multi-way interactions in cis with $oligo"} += 3;
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]}) {         
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]}) {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]} += 1;
                }
            }
            if ($#bins_sorted == 3) {                                              # Quintuplets
                $counter{"05d Number of binned quintuplets in cis with $oligo"}++;
                $counter{"05 Total number of binned multi-way interactions in cis with $oligo"} += 6;
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[3]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[3]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[3]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]}) {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[3]}) {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[3]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[3]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[3]}) {
                    $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[3]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[3]} += 1;
                }
            }
            if ($#bins_sorted == 4) {                                              # Sextuplets
                $counter{"05e Number of binned sextuplets in cis with $oligo"}++;
                $counter{"05 Total number of binned multi-way interactions in cis with $oligo"} += 10;
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[1]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[2]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[3]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[3]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[3]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[4]}) {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[4]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[0]}{$bins_sorted[4]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]}) {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[2]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[3]}) {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[3]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[3]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[4]}) {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[4]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[1]}{$bins_sorted[4]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[3]}) {
                    $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[3]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[3]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[4]}) {
                    $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[4]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[2]}{$bins_sorted[4]} += 1;
                }
                unless (exists $bin_hash{$oligo}{$bins_sorted[3]}{$bins_sorted[4]}) {
                    $bin_hash{$oligo}{$bins_sorted[3]}{$bins_sorted[4]} = 1;
                }
                else {
                    $bin_hash{$oligo}{$bins_sorted[3]}{$bins_sorted[4]} += 1;
                }
            }
        }
    }
}


# Print multi-way interactions to output file

foreach my $oligo (sort keys %frag_hash) {
    my $output_file = "$dir/$name\_$oligo\_TriC_interactions.txt";
    open (OUT, ">$output_file") or die "can't open output report file $output_file";
    foreach my $X (sort keys %{$frag_hash{$oligo}}) {
        foreach my $Y (sort keys %{$frag_hash{$oligo}{$X}}) {
            print OUT "$X\_$Y\t$frag_hash{$oligo}{$X}{$Y}\n";
        }
    }
}

foreach my $oligo (sort keys %bin_hash) {
    my $bin_output_file = "$dir/$name\_$oligo\_binned$bin_size\_TriC_interactions.txt";
    open (OUT, ">$bin_output_file") or die "can't open output report file $bin_output_file";
    foreach my $X (sort keys %{$bin_hash{$oligo}}) {
        foreach my $Y (sort keys %{$bin_hash{$oligo}{$X}}) {
            print OUT "$X\_$Y\t$bin_hash{$oligo}{$X}{$Y}\n";
        }
    }
}

# Print counts to report file

foreach my $key (sort keys %counter) {
    printf REP "%-8s %s\n", $key, $counter{$key};
}

# Generate Capture-C-like tracks containing all mapped reporters in cis for each viewpoint

print REP "\nCapture-C like tracks:\n";

unless (open(TRACKHUBC, ">$public_folder/tracks.txt")){die "Cannot open file $public_folder/tracks.txt, stopped ";}
frag_to_wig(\%frag_hash,"",'178,65,244');
close TRACKHUBC;

unless (open(TRACKHUBC, ">$public_folder/tracks_bin.txt")){die "Cannot open file $public_folder/tracks_binned.txt, stopped ";}
frag_to_wig(\%bin_hash,"BIN",'66, 215, 244');
close TRACKHUBC;

#########################################################################################################################################################################
# Creates a track hub (non-binned data)
unless (open(TRACKHUBA, ">$public_folder/hub.txt")){die "Cannot open file $public_folder/hub.txt, stopped ";}
unless (open(TRACKHUBB, ">$public_folder/genomes.txt")){die "Cannot open file $public_folder/genomes.txt, stopped ";}

print TRACKHUBA "hub $name\_TriC
shortLabel $name\_TriC
longLabel $name\_TriC
genomesFile $server/$server_folder/genomes.txt
email $email";

print TRACKHUBB "genome $genome
trackDb $server/$server_folder/tracks.txt";

close TRACKHUBA;
close TRACKHUBB;

print REP "\nThe track hub (non-binned data) can be found at:
$server/$server_folder/hub.txt
This URL just needs to be pasted into the UCSC genome browser\n\n";

close REP;

print STDOUT "\nThe track hub (non-binned data) can be found at:
$server/$server_folder/hub.txt
This URL just needs to be pasted into the UCSC genome browser\n\n";

# Creates a track hub (binned data)
unless (open(TRACKHUBA, ">$public_folder/hub_bin.txt")){die "Cannot open file $public_folder/hub_bin.txt, stopped ";}
unless (open(TRACKHUBB, ">$public_folder/genomes_bin.txt")){die "Cannot open file $public_folder/genomes_bin.txt, stopped ";}

print TRACKHUBA "hub $name\_TriC_binned
shortLabel $name\_TriC_binned
longLabel $name\_TriC_binned
genomesFile $server/$server_folder/genomes_bin.txt
email $email";

print TRACKHUBB "genome $genome
trackDb $server/$server_folder/tracks_bin.txt";

close TRACKHUBA;
close TRACKHUBB;

print REP "\nThe track hub (binned data) can be found at:
$server/$server_folder/hub_bin.txt
This URL just needs to be pasted into the UCSC genome browser\n\n";

close REP;

print STDOUT "\nThe track hub (binned data) can be found at:
$server/$server_folder/hub_bin.txt
This URL just needs to be pasted into the UCSC genome browser\n\n";


# Make html copy of the REPORT file ..

# system ("cp $report $store_bigwigs_here_folder\/description.html") == 0 or die "couldn't copy public log file\n";

my $reporthtml = "$dir/description.html";
open (REPHTML, ">$reporthtml") or die "can't open output report html file $reporthtml";
print REPHTML "<pre>\n";
open(REP, "<$report") or die "can't open output report file $report";
while (my $line = <REP>) {
    print REPHTML $line 
}
print REPHTML "</pre>\n";
close REPHTML;
close REP;

system ("mv $reporthtml $store_bigwigs_here_folder") == 0 or die "couldn't move log html file\n";

######################################################################################################################################################

sub binary_search {
    my ($chr_array, $start, $stop, $counter_hash) = @_;
    my $mid_value = ($start + $stop) / 2;
    my $first_pos = 0;
    my $last_pos = scalar @$chr_array - 1; 
    my $counter = 0;
    if (($mid_value < $$chr_array[$first_pos]) or ($mid_value > $$chr_array[$last_pos])) {
        #$$counter_hash{"Binary search error: search outside range of restriction enzyme coordinates"}++;
        return ("error", "error")
    }
    for (my $i = 0; $i < 99; $i++) {
        my $mid_search = int(($first_pos + $last_pos) / 2);
        if ($$chr_array[$mid_search] > $$chr_array[$mid_search + 1]) {
            #$$counter_hash{"Binary search error: restriction enzyme array coordinates not in ascending order"}++;
            return ("error", "error")
        }
        if (($$chr_array[$mid_search] <= $mid_value) and ($$chr_array[$mid_search + 1] > $mid_value)) {    
            if (($$chr_array[$mid_search] <= $start + 2) and ($$chr_array[$mid_search + 1] >= $stop - 2)) {    
                return ($$chr_array[$mid_search], $$chr_array[$mid_search + 1] - 1)
            }
            else {
                #$$counter_hash{"Binary search error: fragment overlaps multiple restriction sites"}++;
                return ("error", "error");
            }
        }       
        elsif ($$chr_array[$mid_search] > $mid_value) {
            $last_pos = $mid_search - 1;
        }    
        elsif ($$chr_array[$mid_search] < $mid_value) {
            $first_pos = $mid_search + 1;
        }
        else {
            #$$counter_hash{"Binary search error: end of loop reached"}++;
        }
    }
    #$$counter_hash{"Binary search error: couldn't map read to fragments"}++;
    return ("error", "error")
}

sub frag_to_wig {                                                                       
    my ($hashref, $file_name, $color) = @_;
    my %coord_hash;
    foreach my $oligo (keys %$hashref) {
        my $tracks_out = "$name\_$oligo$file_name";
        foreach my $read (keys %{$$hashref{$oligo}}) {
            if (exists $$hashref{$oligo}{$read}{"reporter_array"}) {
                my $length = scalar @{$$hashref{$oligo}{$read}{"reporter_array"}};
                for (my $i = 0; $i < $length; $i++ ) {
                    my ($chr, $start, $stop) = split(/\W/, $$hashref{$oligo}{$read}{"reporter_array"}[$i]);                   
                    my @range = ($start..$stop);
                    for(my $j = 0; $j < $#range; $j++) {
                        ${$coord_hash{$oligo}{$chr}{$range[$j]}}++;
                    }
                }
            }
        }
    open (TRACKS_OUT, ">$dir/$tracks_out.wig");
    foreach my $chr (sort {$coord_hash{$oligo}{$a} <=> $coord_hash{$oligo}{$b}} keys %{$coord_hash{$oligo}}) {
        print TRACKS_OUT "variableStep  chrom=chr$chr\n";
        foreach my $coord (sort {$a <=> $b} keys %{$coord_hash{$oligo}{$chr}}) {
            my $count = ${$coord_hash{$oligo}{$chr}{$coord}};
            print TRACKS_OUT "$coord\t$count\n";
        }
    }
    system ("wigToBigWig -clip $dir/$tracks_out.wig $ucscsizes $dir/$tracks_out.bw") == 0 or die "couldn't bigwig files\n";
    system ("mv $dir/$tracks_out.bw $store_bigwigs_here_folder") == 0 or die "couldn't move files\n";		
    system ("chmod 755 $store_bigwigs_here_folder/$tracks_out.bw") == 0 or die "couldn't chmod files\n";   
    print REP "track type=bigWig name=\"$tracks_out\" description=\"c-trap $tracks_out\" bigDataUrl=$server/$server_folder/$tracks_out.bw\n";
    
    print TRACKHUBC
"track $tracks_out
type bigWig
longLabel $tracks_out
shortLabel $tracks_out
bigDataUrl $server/$server_folder/$tracks_out.bw
visibility full
priority 200
color $color
autoScale on
alwaysZero on
windowingFunction maximum
html description

";
    
    
    }
}

