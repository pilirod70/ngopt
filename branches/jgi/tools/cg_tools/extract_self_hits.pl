#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
if (scalar(@ARGV) != 1) {
	print "Usage: ".basename($0)." <m8_blast_output>\n";
	exit;
}
open (IN,"<",$ARGV[0]);
while(<IN>){
	chomp;
	my @ar = split;
	if ($ar[0] eq $ar[1]){
		print STDOUT $_."\n";

	}
}




