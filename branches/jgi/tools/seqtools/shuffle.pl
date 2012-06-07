#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV != 2) {
	print "Usage: shuffle.pl <read1.fq> <read2.fq>\nOutput is printed to stdout.\n";
	exit;
}

my $r1 = shift;
my $r2 = shift;

open (R1,"<$r1");
open (R2,"<$r2");

while (my $hdr = <R1>){
	my $seq = <R1>;
	my $qhdr = <R1>;
	my $qseq = <R1>;
	print $hdr.$seq.$qhdr.$qseq;
	$hdr = <R2>;
	$seq = <R2>;
	$qhdr = <R2>;
	$qseq = <R2>;
	print $hdr.$seq.$qhdr.$qseq;
}
