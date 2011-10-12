#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

open (IN,"<",$ARGV[0]);
while(<IN>){
	chomp;
	my @ar = split;
	if ($ar[0] eq $ar[1]){
		if ($ar[6] == $ar[9] && $ar[7] == $ar[8]){
			print STDOUT $_."\n";
		}
	}
}




