#!/usr/bin/perl -w
use strict;
use warnings;
my $ins = shift;
my $rdlen = shift;

my $line = <>;
chomp $line;
my $genome_length=0;
while($line =~ m/^@.*LN:(\d+)/) {
	$genome_length += $1;
	$line = <>;
	chomp $line;
}

my $nreads = 1;
while($line = <>){
	$nreads++;
}

my $minlink = $nreads*$ins/$genome_length;

