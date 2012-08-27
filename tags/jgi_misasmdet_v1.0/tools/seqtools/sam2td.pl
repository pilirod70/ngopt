#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV < 2) {
	print "Usage: sam2td.pl <contig1> <contig2> <sam|stdin>.\nOutput is printed to stdout\n";
	exit;
}

my $ctg1 = shift;
my $ctg2 = shift;

while(my $line1 = <>){
	next if ($line1 =~ m/^@/);
	my $line2 = <>;
	my @ar1 = split(/\t/,$line1);
	my @ar2 = split(/\t/,$line2);
	next unless ($ar1[3] =~ /\d+/ && $ar2[3] =~ /\d+/);
	if ($ar1[2] eq $ctg1 && $ar2[2] eq $ctg2) {
		print $ar1[3]."\t".$ar2[3]."\n";
	} elsif ($ar2[2] eq $ctg1 && $ar2[2] eq $ctg2){
		print $ar2[3]."\t".$ar1[3]."\n";
	}
}

