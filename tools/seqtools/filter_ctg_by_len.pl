#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV != 2 && @ARGV != 3){
	print "Usage: filter_fasta_by_len.pl <min_length> <in.fasta|stdin>\n";
	exit;
}

my $cut = shift;

my $hdr = "";
my $seq = "";
my $len = 0;
while (<>){
	if ($_ =~ /^>/){
		if (length($hdr)){
			if ($len >= $cut) {
				print  $hdr.$seq;
			} 
			$len = 0;
			$seq = "";
		}
		$hdr = $_;	
	} else {
		$seq .= $_;
		$len += length($_)-1;
	}
}

