#!/usr/bin/perl -w
# Convert Illumina _export.txt to a FastQ
# (c) Aaron Darling 2010
use strict;

while( my $line = <STDIN> ){
	my @vals = split( /\s+/, $line );
	my $hdr = $vals[0]."_".$vals[1].":".$vals[2].":".$vals[3].":".$vals[4].":".$vals[5].":".$vals[6]."/".$vals[7];
	print ">$hdr\n";
	print $vals[8]."\n";
	print "+$hdr\n";
	print $vals[9]."\n";
}

