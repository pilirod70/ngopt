#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 2) {
	print STDOUT "Usage: ".basename($0)." <maf_file> <output_prefix>\n";	
	exit 1;
}

my $maf_file = $ARGV[0];
my $out_base = $ARGV[1];

my $block_num=0;
open (MAF,"<",$maf_file);
my $first = 1;
while (<MAF>) {
	chomp;
	if (length($_) == 0) {
		next;
	}
	my $c = substr($_,0,1);
	if ($c eq "a") {
		if ($first){
			$first = 0;
		} else {
			close FA;	
		}
		open (FA,">",$out_base."-".$block_num.".fasta");
		$block_num++;
		next;
	} elsif ($c eq "s") {
		my @fields = split(' ',$_);
		shift(@fields);
		my $seq = pop @fields;
		my $hdr = join('|',@fields);
		print FA ">".$hdr."\n";
		while (length($seq) > 80) {
			print FA substr($seq,0,80)."\n";
			$seq = substr($seq,80);
		}
		print FA $seq."\n";
	} else {
		next;	
	} 
}
close FA;



