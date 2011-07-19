#!/usr/bin/perl -w
use strict;
use warnings;

my %seen = ();

while(<>){
	chomp;
	my ($hdr, $seq, $qual) = split /\t/ ;
	if (defined($seen{$hdr})){
		print "\@$hdr/2\n$seq\n+$hdr/2\n$qual\n";
	} else {
		print "\@$hdr/1\n$seq\n+$hdr/1\n$qual\n";
		$seen{$hdr} = 1;
	}	
}
