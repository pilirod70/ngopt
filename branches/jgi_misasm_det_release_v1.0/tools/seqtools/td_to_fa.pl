#!/usr/bin/perl -w
use strict;
use warnings;

my %seen = ();

while(<>){
	chomp;
	my ($hdr, $seq) = split /\t/ ;
	if (defined($seen{$hdr})){
		print "\>$hdr/2\n$seq\n";
	} else {
		print "\>$hdr/1\n$seq\n";
		$seen{$hdr} = 1;
	}	
}
