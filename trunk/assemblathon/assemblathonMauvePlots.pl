#!/usr/bin/env perl
# (c) Aaron Darling 2011
# Licensed under the GPL
# Program to generate plots on a directory containing assemblathon results
# Usage: <assemblathon directory>
#

use strict;
use warnings;


sub getsomestuff
{
	my $type = shift;
	my $finder_cl = "find $ARGV[0] -name \"*__sum.txt\" \| grep $type \| sort |";
	open(FINDER, $finder_cl);
	my $alldat="";
	while(my $line=<FINDER>){
		chomp($line);
		$alldat .= " $line";
	}
	return $alldat;
}

sub getstuff
{
	my $type = shift;
	my $alldat = getsomestuff($type);
	$alldat =~ s/__sum.txt//g;
	my $plot_cl = "assemblathonMauvePlots.R $alldat";
	print "$plot_cl\n";
	`$plot_cl`;
	`rmdir -rf $type`;
	`mkdir -p $type`;
	`mv *.pdf $type`;
}

sub getsummary
{
	my $alldat = getsomestuff("c");
	`cat $alldat > summaries.txt`;
	`perl -p -i -e "s/NumContigs.*//g" summaries.txt`;
}

getsummary();
getstuff("contigs");
getstuff("scaffolds");

