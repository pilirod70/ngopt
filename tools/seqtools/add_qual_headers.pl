#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV)!=1){
	print STDOUT "Usage: ".basename($0)." in.fastq > out.fixed.fastq\n";
	exit 1;
}

my $fq_file = $ARGV[0];
open(FQ,"<",$fq_file);

while (<FQ>) {
	my $hdr = substr($_,1);
	my $seq = <FQ>;
	<FQ>;
	my $qual = <FQ>;
	print STDOUT "@".$hdr;
	print STDOUT $seq;
	print STDOUT "+".$hdr;
	print STDOUT $qual;
}

