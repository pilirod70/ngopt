#!/usr/bin/perl -w
use strict;
use warnings;

my $cutoff = shift;
my $base = shift;

my $seq = "";

open (ABOVE,">$base.gt$cutoff.fasta");
open (BELOW,">$base.lt$cutoff.fasta");

my $cov = -1;
my $hdr = "";
while (<>) {
	if ($_ =~ />NODE_\d+_length_\d+_cov_(\d+)/){
		if ($cov > 0) {
			if ($cov >= $cutoff) {
				print ABOVE $hdr.$seq;
			} else {
				print BELOW $hdr.$seq;
			}
		}
		$cov = $1;
		$hdr = $_;
		$seq = "";
	} else {
		$seq .= $_;
	}
}

if ($cov > 0) {
	if ($cov >= $cutoff) {
		print ABOVE $hdr.$seq;
	} else {
		print BELOW $hdr.$seq;
	}
}

close ABOVE;
close BELOW;
