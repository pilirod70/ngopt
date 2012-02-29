#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV < 1) {
	print "Usage: unshuffle.pl <output_basename> <in.fastq|stdin>\n";
	exit;
}

my $base = shift;

open (P1,">$base\_1.fastq");
open (P2,">$base\_2.fastq");

while (my $hdr = <>) {
	my $seq = <>;
	my $qhdr = <>;
	my $qseq = <>;
	print P1 $hdr.$seq.$qhdr.$qseq;
	$hdr = <>;
	$seq = <>;
	$qhdr = <>;
	$qseq = <>;
	print P2 $hdr.$seq.$qhdr.$qseq;
}

close P1;
close P2;
