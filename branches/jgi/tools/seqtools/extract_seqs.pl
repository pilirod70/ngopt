#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 2) {
	print "Usage: ".basename($0)." <list_file> <fasta_file>\n";
	exit 1;
}
my %hash = ();
open(LIST,"<",$ARGV[0]);
while (<LIST>) {
	chomp;
	$hash{$_} = 1;
}
my $keep = 0;
open(FA,"<",$ARGV[1]);
while (<FA>) {
	chomp;
	if ($_ =~ m/^>/){
		if (defined($hash{substr($_,1)})){
			$keep = 1;
			print $_."\n";
		} else {
			$keep = 0;
		}
	} else {
		if ($keep) {
			print $_."\n";
		} else {
			next;
		}	
	}
}

