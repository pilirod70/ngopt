#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 2) {
	print STDOUT "Usage: ".basename($0)." <maf_file> <output_fasta>\n";	
	exit 1;
}

my $maf_file = $ARGV[0];
my $out_file = $ARGV[1];

my %seqs;
my $block_num=0;
open (MAF,"<",$maf_file);
my $first = 1;
while (<MAF>) {
	chomp;
	if (length($_) == 0) {
		next;
	}
	my $c = substr($_,0,1);
	if ($c eq "s") {
		my @fields = split(' ',$_);
		if (!defined($seqs{$fields[1]})) {
			$seqs{$fields[1]} = $fields[6];
		} else {
			$seqs{$fields[1]}.=$fields[6];
		}
	} else {
		next;	
	} 
}
open(FA,">",$out_file);
while ((my $key, my $value) = each(%seqs)){
	print FA ">".$key."\n";
	while (length($value) > 80) {
		print FA substr($value,0,80)."\n";
		$value = substr($value,80);
	}
	print FA $value."\n";
}
close FA;



