#!/usr/bin/perl -w

use strict;
use warnings;

if (scalar(@ARGV) != 2) {
	print STDOUT "Splits a set of reads into mapped and unmapped files\n";
	print STDOUT "Usage: exmap.pl <sam_file> <reads_file>\n";
	print STDOUT "Reads can be in either fasta or fastq format\n";
	print STDOUT "Two output files will be created: reads_base.mapped.fasta/q and reads_base.unmapped.fasta/q\n";
	exit 1;
}	

my %keepers;

open (SAM ,"<", $ARGV[0]);
open (FILT, ">", $ARGV[0].".filtered");
while (<SAM>) {
	chomp;
	if ( substr($_,0,1) ne "@") {
		my @tmp = split;
		my $cigar = $tmp[5];
		if ($tmp[5] eq  "49M") {
			$keepers{$tmp[0]} = 1;
			print FILT join('\t',@tmp)."\n";
		}
	}
}

open (FQ, "<", $ARGV[1]);
my $base = substr($ARGV[1],0,rindex($ARGV[1],".f[a,q]"));
my $map_file = $base.".mapped.fasta";
my $umap_file = $base.".unmapped.fasta";
my $first = 1;
my $fastq = 0;
if (getc(FQ) eq '@') {
	$fastq = 1;
	$map_file = $base.".mapped.fastq";
	$umap_file = $base.".unmapped.fastq";
	
}

my $mapped = 0;
my $unmapped = 0;
open (MAP, ">", $map_file);
open (UMAP,">", $umap_file);
while (<FQ>) {
	chomp;
	my $line;
	my $hdr;
	if ($first){
		$hdr = $_;
		$first = 0;
		print STDOUT $hdr."\n";
	} else {
		$hdr = substr($_,1,length($_)-1);
	}
	if ($keepers{$hdr}) {
		++$mapped; 
		print STDOUT "PLUS ONE";
		print MAP $_."\n";
		print MAP <FQ>;
		if ($fastq) {
			print MAP <FQ>;
			print MAP <FQ>;
		}
	} else {
		++$unmapped;
		print UMAP $_."\n";
		print UMAP <FQ>;
		if ($fastq) {
			print UMAP <FQ>;
			print UMAP <FQ>;
		}
	}
}

close MAP;
close UMAP;

print STDOUT "Wrote ".$mapped." to ".$map_file." and ".$unmapped." to ".$umap_file."\n"; 

