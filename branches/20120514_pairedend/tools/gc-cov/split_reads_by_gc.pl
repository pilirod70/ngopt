#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV != 2 && @ARGV != 3){
	print "Usage: split_reads_by_gc.pl <cutoff> <output_base> <in.fastq|stdin>\n";
	exit;
}

my $cut = shift;
my $base = shift;
my $GC = 0;
my $LEN = 0;

open (ABOVE,">$base.gt$cut.fastq");
open (BELOW,">$base.lt$cut.fastq");

my $hdr = "";

my $seq = "";
while (my $hdr = <>){
	my $seq = <>;
	my $qhdr = <>;
	my $qseq = <>;
	my $gc = gc($seq);
	print STDERR $gc."\n";
	if (gc($seq) >= $cut){
		print ABOVE $hdr.$seq.$qhdr.$qseq;
	} else {
		print BELOW $hdr.$seq.$qhdr.$qseq;
	}
}

close ABOVE;
close BELOW;

sub gc {
	my $seq = shift;
	chomp $seq;
	my $count = 0;
	my $len = 0;
	for (my $i = 0; $i < length($seq); $i++) {
		my $char = substr($seq,$i,1);
		$len++ if ($char =~ /[ATCGatcg]/);
		$count++ if ($char =~ /[GCgc]/);
	}
	return $count/$len;
}

