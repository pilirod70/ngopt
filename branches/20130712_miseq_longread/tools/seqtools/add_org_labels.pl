#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) == 0) {
	print "Usage: ".basename($0)." <org> <fasta_file>\n";
	exit 1;
}

my $org=shift;

while (<>) {
	if ($_=~m/^>/){
		chomp;
		my $ctg = (split(/\|/,substr($_,1)))[0];
		print ">$org-$ctg\n";
	} else {
		print $_;
	}
}
