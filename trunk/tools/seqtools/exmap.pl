#!/usr/bin/perl -w

use strict;
use warnings;

my %keepers;

open (SAM ,"<", $ARGV[0]);

while (<SAM>) {
	chomp;
	if (substr($_,0,1) eq "@") {
		continue;
	} else {
		my @tmp = split;
		$keepers{$tmp[0]} = 1;
	}
}

open (FQ, "<", $ARGV[1]);

while (<FQ>) {
	chomp;
	my $line;
	my $hdr = substr($_,1,length($_)-3);
	if ($keepers{$hdr}) {
		print STDOUT $_."\n";
		print STDOUT <FQ>;
		print STDOUT <FQ>;
		print STDOUT <FQ>;
	} else {
		$line = <FQ>;
		$line = <FQ>;
		$line = <FQ>;
	}
}


