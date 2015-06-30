#!/usr/bin/perl -w
####################################################################################################
### Get GC Content, CpG Content                                                                  ###
### Usage: get_gc_content.pl <fasta file>                                                        ###
### This program takes a fasta file as it's first (and only) parameter.                          ###
###                                                                                              ###
### It returns a tab delimited file (gc_out.txt): column 1 = header ID (everything between ">"   ###
### and the first space in the header), column 2 = gc content for the fasta entry,               ###
### column 3 = CpG content.                                                                      ###
###                                                                                              ###
### Fabio Ticconi                                                                                ###
### June 30, 2015                                                                                ###
###                                                                                              ###
### Jennifer Meneghin                                                                            ###
### July 23, 2009                                                                                ###
###                                                                                              ###
### This script now works properly with sequences that contain spaces.                           ###
### September 20, 2010                                                                           ###
###                                                                                              ###
### This script now also returns the total nucleotide count, along with the number of of         ###
### A's, G's, C's and T's for each fasta record.                                                 ###
### September 21, 2010                                                                           ###
####################################################################################################
#---------------------------------------------------------------------------------------------------------------------------
#Deal with passed parameters
#---------------------------------------------------------------------------------------------------------------------------

use strict;
use warnings;

if ($#ARGV == -1) {
    usage();
    exit;
}

my $fasta_file = $ARGV[0];
my $out_file = "gc_out.txt";

unless ( open(IN, "$fasta_file") ) {    
    print "Got a bad fasta file: $fasta_file\n\n";
    exit;
}
unless ( open(OUT, ">$out_file") ) {
    print "Couldn't create $out_file\n";
    exit;
}

print "Parameters:\nfasta file = $fasta_file\noutput file = $out_file\n\n";

#---------------------------------------------------------------------------------------------------------------------------
#The main event
#---------------------------------------------------------------------------------------------------------------------------

print OUT "ID\t% GC\t% CpG\tObs/Exp CpG\tTotal Count\tG Count\tC Count\tA Count\tT Count\n";

my $seq = "";
my $id;

while (<IN>)
{
    chomp;

    if (/^>/)
    {
	    #finish up previous line.
	    if (length($seq) > 0)
        {
            &process_it;
        }

	   #start new line.
	   $id = $_;
	   $id =~ s/^>(.+?)\s.+$/$1/g;
	   print OUT "$id\t";
    }
    else
    {
	   $seq = $seq . $_;
    }
}

#finish up last line.
&process_it;

close(IN);
close(OUT);

sub usage {
    print "Get GC Content, CpG Content\n";
    print "Usage: get_gc_content.pl <fasta file>\n";
    print "This program takes a fasta file as it's first (and only) parameter.\n\n";
    print "It returns a tab delimited file (gc_out.txt): column 1 = header ID (everything between \">\"\n";
    print "and the first space in the header), and column 2 = gc content for the fasta entry.\n\n";
    print "Jennifer Meneghin\n";
    print "July 23, 2009\n\n";
    print "Updated September 20, 2010:\n";
    print "This script now works properly with sequences that contain spaces.\n\n";
    print "Updated September 21, 2010:\n";
    print "This script now also returns the total nucleotide count, along with the number of of A's, G's, C's and T's for each fasta record.\n\n";
}

sub process_it {
    my @letters = split(//, $seq);

    my $length = @letters;

    my $l = "";

    my $totalcount = 0;
    my $acount = 0;
    my $tcount = 0;
    my $gcount = 0;
    my $ccount = 0;

    my $prev = "x";
    my $cpgcount = 0;

    my $cpgcontent = 0;
    my $gccontent = 0;
    my $expcpgcont = 0;

    foreach my $i (@letters) {
     if (lc($i) =~ /[a-z]/) {
         $totalcount++;
         $l = lc($i);
     }

     if ($l eq "a") {
         $acount++;
     }
     if ($l eq "t") {
         $tcount++;
     }
     if ($l eq "g") {
         $gcount++;
     }
     if ($l eq "c") {
         $ccount++;
     }
     if ($prev eq "c" && $l eq "g") {
        $cpgcount++;
    }

    $prev = $l;
}

if ($totalcount > 0) {
    $gccontent = 100 * (($gcount + $ccount) / $totalcount);
    $cpgcontent = 100 * ($cpgcount / ($totalcount - 1));
    $expcpgcont = ($cpgcount * $length) / ($ccount * $gcount);
}

print OUT "$gccontent\t$cpgcontent\t$expcpgcont\t$totalcount\t$gcount\t$ccount\t$acount\t$tcount\n";
$seq = "";
}
