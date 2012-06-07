#!/usr/bin/env perl 
# Program to submit assemblathon scoring jobs
# to a cluster managed with Sun Grid Engine
# (c) 2011 Aaron Darling
# Licensed under the GPL

use warnings;
use strict;

my $scratchDir = "/tmp";

if(@ARGV!=3){
	print STDERR "Usage: assemblathon.pl <reference genome> <assembly directory> <output directory>\n";
	exit -1;
}

my $reference = $ARGV[0];
my $assemblyDir = $ARGV[1];
my $outputDir = $ARGV[2];

my @assemblies;
open( ASSLIST, "ls -1 $assemblyDir |" );
while( my $line = <ASSLIST> ){
	chomp $line;
	if($line =~ /\.fa$/){
		push(@assemblies, $line);
	}
}

foreach my $assembly(@assemblies){
	my $scoreDir = "$outputDir"."$assembly"."_scoring";
	open (LSSER, "find $scoreDir -wholename \"*/*__sum.txt\" |");
	my $line2 = <LSSER>;
	print $line2 if defined($line2);
	next if defined($line2);
	open (LSSER, "find $scoreDir -size +1b -wholename \"*/*/alignment4/alignment4\" |");
	my $line = <LSSER>;
	my $score_qsub_cl="";
	if( defined($line) ){
		# still need to run scoring
#		chomp $line;
#		$score_qsub_cl = "qsub -q all.q -q eisen.q /home/koadman/bin/assemblathon_sgescoreonly.sh $line $scoreDir";
	}else{
		`rm -rf $scoreDir`;
		`mkdir -p $scoreDir`; 
		$score_qsub_cl = "qsub -q all.q -q eisen.q /home/koadman/bin/assemblathon_sgescoring.sh $reference $assemblyDir/$assembly $scoreDir";
	}
	print $score_qsub_cl."\n";
	`$score_qsub_cl`;
}

