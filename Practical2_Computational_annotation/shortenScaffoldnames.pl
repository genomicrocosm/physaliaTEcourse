#!/usr/bin/env/perl

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# perl shortenScaffoldnames.pl <infile> <outfile>
# =====================================================
# Shortens FASTA headers at first occurrence of " " 
# =====================================================
# Alexander Suh				13 Feb 2015
# =====================================================
# Example:
# Old FASTA header: ">DF091493.1 Bombyx mori DNA, 
# scaffold: Bm_scaf1566, strain: p50T/Dazao"
# New FASTA header: ">DF091493.1"
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use warnings;

my $file = $ARGV[0];
my $out = $ARGV[1];

open (OUT, ">$out");

open (IN, $file);
while (my $line = <IN>) {
    if ($line =~ m/^>/) {
        chomp $line;
	my @parts = split(/ /, $line);
	my $id = $parts[0];
	print OUT $id."\n";
    } else {
        printf OUT $line;
    }
}
close (IN);
close (OUT); 