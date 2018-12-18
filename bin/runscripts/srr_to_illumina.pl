#!/usr/bin/perl -w
use strict;

# Jelena Telenius 2014 - Clipping the last base from 3' end
# James Davies 2014 - Reading in the fastq in batches of four lines (from script dpnII2E.pl )

my $line_counter = 0; my $f_counter = 0; my %hash; my $flag;
my @line_labels = qw(name seq spare qscore);

my $filename = $ARGV[0];
my $r1_or_2 = $ARGV[1];
$filename =~ /(.*)\.(fastq|fq)/; my $filename_out= $1."_asIllumina.$2";

unless (open(FH, $filename)) {print "Cannot open file $filename\n"; exit;}

# opens a file in append modes for the output of the data
open FHOUT, ">$filename_out" or die $!;   
  
while ($hash{$line_labels[$f_counter]}=<FH>)  #assigns each fq line to the hash in batches of 4
{
chomp $hash{$line_labels[$f_counter]};
$f_counter++; $line_counter++;

if ($f_counter==4)
    {
     if ($hash{"name"} =~ /^(.*) /)
        {
        $hash{"name"}=$1.":1:1:1:1:1:1 ".$r1_or_2.":1:1:1";
        # @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 1:N:0:GAGTTAGT
     }
        
         print FHOUT $hash{"name"}."\n";
         print FHOUT $hash{"seq"}."\n";  
         print FHOUT "+\n";
         print FHOUT $hash{"qscore"}."\n";  


    $f_counter=0

    }

}
