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
use strict;

my $min_length = 19;
# Specifications
my $line_counter = 0; my $f_counter = 0; my $gatc_counter=0; my %hash; my $flag; my $k;
my $gatc_0_counter=0; my $gatc_1_counter=0; my $gatc_2_counter=0; my $gatc_3_counter=0; my $fragment_counter=0; my $printed_counter=0;
my @line_labels = qw(name seq spare qscore);

my $filename = $ARGV[0];
my $filetype = $ARGV[1];

unless ($filename =~ /(.*)\.(fastq|fq)/) {die"filename does not match fq"};
my $filename_out= $1."_REdig.$2"; my $filename_dump = $1."_merge_dump.$2";
my $file_name=$1;
my $file_path="";
if ($file_name =~ /(.*\/)(\V++)/) {$file_path = $1; $file_name = $2};

if ( ($filetype ne "NONFLASHED") and ($filetype ne "FLASHED") ) {print "Second parameter has to be file type : either FLASHED or NONFLASHED\n"; exit;}

unless (open(FH, $filename)) {print "Cannot open file $filename\n"; exit;}

# opens a file in append modes for the output of the data
open FHOUT, ">$filename_out" or die $!;   
  

while ($hash{$line_labels[$f_counter]}=<FH>)  #assigns each fq line to the hash in batches of 4
{
chomp $hash{$line_labels[$f_counter]};
$f_counter++; $line_counter++;


if ($f_counter==4)
    {
    # name formats @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 1:N:0:GAGTTAGT run1
    # name formats @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 2:N:0:GAGTTAGT run2
    # /(.*):(.*):(.*):(.*):(.*):(.*):(\d++) (\d):(.*):(.*):(.*)/
    if ($hash{"name"} =~ /(.*:.*:.*:.*:.*:.*:\d++) (\d):(.*:.*:.*)/)
        {
        $hash{"PE"}=$2;
        $hash{"new name"} = $1;
        #$hash{"new name end"} = " ".$2.":".$3
        }
      # Other option : illumina old format :
     # @HWUSI-EAS100R:6:73:941:1973#0/1
     # @HWUSI-EAS100R:6:73:941:1973#0/2
     elsif($hash{"name"} =~ /(.*:.*:.*:.*:.*#\d++)\/(\d)/)
     {
        $hash{"PE"}=$2;
        $hash{"new name"} = $1;
        #$hash{"new name end"} = " ".$2.":".$3
        }    

# codes for marking
# :0 no AAGCTT cut
# :1 AAGCTT in LEFT
# :2 AAGCTT in RIGHT
# :3 AAGCTT in RIGHT and LEFT

#  : The cut sequence for HindIII is A^AGCTT
        
if ($hash{"seq"} =~ /AAGCTT/)
        {

         my @gatc_splits = split/AAGCTT/, $hash{"seq"};

             for (my $i=0; $i<$#gatc_splits+1;$i++)
             {
             $fragment_counter++;
             
             # ------------------------------------------
             # First fragment of the red-in line
             if ($i==0) #first fragment
                {
                
                # ----------------------------------------------
                # Setting the sseq and sqscore for the found fragment - to be updated in late if clauses
                if ($i==$#gatc_splits) 
                {
                #last fragment if there is just one fragment    
                $hash{"split$i"}{"sseq"}= "$gatc_splits[$i]";
                }
                else 
                {
                # other cases (normal cases)
                $hash{"split$i"}{"sseq"}= "$gatc_splits[$i]AAGCTT";
                }
                $hash{"split$i"}{"sqscore"}= substr ($hash{"qscore"},0,length $hash{"split$i"}{"sseq"});
                # ----------------------------------------------
                
                 # If starts with AAGCTT or AGCTT or not
                if ($hash{"split$i"}{"sseq"} =~ /^AAGCTT/)
                {
                # for testing :    
                # $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3a";
                # this seems to never happen - see testing results in workingDiary66 .
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3";
                if (length $hash{"split$i"}{"sseq"}>$min_length){$gatc_3_counter++;}
                
                }
                elsif ($hash{"split$i"}{"sseq"} =~ /^AGCTT/)
                {
                $hash{"split$i"}{"sseq"}= "A" . $hash{"split$i"}{"sseq"};
                $hash{"split$i"}{"sqscore"}= substr ($hash{"split$i"}{"sqscore"},0,1) . $hash{"split$i"}{"sqscore"};
                
                # for testing :
                # $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3b";
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3";                
                if (length $hash{"split$i"}{"sseq"}>$min_length){$gatc_3_counter++;}
                
                }
                else
                {
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":2";
                if (length $hash{"split$i"}{"sseq"}>$min_length){$gatc_2_counter++;}
                }
                
                
                }
                
             # ------------------------------------------
             # Middle fragments of the red-in line (not-first, not-last)
             if ($i!=0 and $i != $#gatc_splits) #middle fragment
                {
                $hash{"split$i"}{"sseq"}= "AAGCTT$gatc_splits[$i]AAGCTT";
                my $offset=0;
                for (my$j=0; $j<$i; $j++){$offset = $offset -6 + length $hash{"split$j"}{"sseq"};} # calculates the offset by looping through the lengths of the left hand side fragments and summing them
                $hash{"split$i"}{"sqscore"}= substr ($hash{"qscore"},$offset,length $hash{"split$i"}{"sseq"});
                # for testing :
                # $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3c";               
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3";
                if (length $hash{"split$i"}{"sseq"}>$min_length){$gatc_3_counter++;}

                }
                
             # ------------------------------------------
             # Last fragment of the red-in line (if there is more than one fragment : if only one fragment, it is handled in "first fragment" above)
             if ($i==$#gatc_splits and $i!=0) #last fragment if there is more than one fragment
                {
                $hash{"split$i"}{"sseq"}= "AAGCTT$gatc_splits[$i]";
                $hash{"split$i"}{"sqscore"}= substr ($hash{"qscore"},-length $hash{"split$i"}{"sseq"},length $hash{"split$i"}{"sseq"});
                
                # If ends with AAGCTT or not
                # (we cannot really know if we had cut site here, though, as if there is A in the end, that may be co-incidence..
                #  so, the reads with A in the end get marked like they would not be cut sites - as we cannot be sure)
                if ($hash{"split$i"}{"sseq"} =~ /AAGCTT$/)
                {
                # for testing :
                # $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3d";
                # this seems to never happen - see testing results in workingDiary66 .
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":3";
                if (length $hash{"split$i"}{"sseq"}>$min_length){$gatc_3_counter++;}
                }
                else
                {
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.":1";                  
                if (length $hash{"split$i"}{"sseq"}>$min_length){$gatc_1_counter++;}
                }
                }         
             # ------------------------------------------
             # end of the first-middle-last fragment read : now  'sseq'  'sqscore'  and  'sname'  have been set for all cases.
             # ------------------------------------------
             
             # If the lenght of the current fragment is long enough to enable unique mapping, proceed to printing ..
             if (length $hash{"split$i"}{"sseq"}>$min_length)
                {
                print FHOUT $hash{"split$i"}{"sname"}."\n".$hash{"split$i"}{"sseq"}."\n+\n".$hash{"split$i"}{"sqscore"}."\n";
                $printed_counter++;
                }
             }

         $gatc_counter++;
        }
    else
    {
         # We print these only if we are in NONFLASHED - otherwise we just skip over..
         if ( $filetype eq "NONFLASHED" )
         {
             print FHOUT $hash{"new name"}.":PE".$hash{"PE"}.":0".":0\n";
             #print FHOUT $hash{"new name"}.":PE".$hash{"PE"}.":0".$hash{"new name end"}."\n";
             print FHOUT $hash{"seq"}."\n";  #error check ."\t".length $hash{"split$i"}{"sseq"};
             print FHOUT "+\n";
             print FHOUT $hash{"qscore"}."\n";   #error check "\t".length $hash{"split$i"}{"sqscore"};
             $printed_counter++;
         }
             $gatc_0_counter++;
             $fragment_counter++;
             
    }
        
    #prints the data in the hash: for (my $i=0; $i<4; $i++){print $i.$hash{$line_labels[$i]}."\n"}print "\n";

    $f_counter=0

    }

#if ($line_counter>1000000){print "$line_counter lines reviewed\n$gatc_counter DpnII sites found\n";exit #}   
}


my $readcount=int( $line_counter/4 );

print "hindIII.pl command run on file: $filename\n$readcount reads reviewed, of which\n$gatc_counter had at least one hindIII site in them\n";
if ( $filetype eq "NONFLASHED" )
{
print "Here READ means either R1 or R2 part of any read - both are called read, and the combined entity does not exist, as flashing was not succesfull for these.\n";
print "In detail, \n$fragment_counter fragments was found, and of these \n$printed_counter fragments were printed - as they were longer than the set threshold $min_length\n";
print "Of the printed fragments (in FASTQ coordinates):\n$gatc_3_counter fragments had LEFT and RIGHT AAGCTT,\n$gatc_1_counter fragments had only LEFT AAGCTT,\n$gatc_2_counter fragments had only RIGHT AAGCTT site,\n$gatc_0_counter fragments had no AAGCTT.\n";
}
else {
print "In detail, \n$fragment_counter fragments was found, and of these \n$printed_counter fragments were printed - as they were longer than the set threshold $min_length , and the read contained at least one AAGCTT\n";
print "Of the printed fragments (in FASTQ coordinates):\n$gatc_3_counter fragments had LEFT and RIGHT AAGCTT,\n$gatc_1_counter fragments had only LEFT AAGCTT,\n$gatc_2_counter fragments had only RIGHT AAGCTT site.\n";
print "$gatc_0_counter reads were longer htan the set treshold $min_length, but had no AAGCTT and thus were discarded .\n";
}


# codes for marking (above)
# :0 no AAGCTT cut
# :1 AAGCTT in LEFT
# :2 AAGCTT in RIGHT
# :3 AAGCTT in RIGHT and LEFT

