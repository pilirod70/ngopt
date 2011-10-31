#!/usr/bin/env perl
# (c) Aaron Darling 2011
# Licensed under the GPL
# Program to generate plots on a directory containing assemblathon results
# Usage: <scored assembly directory>
#

use strict;
use warnings;
use File::Basename;
use File::Spec;

if(@ARGV != 1){
	print "Usage: <directory containing scored assemblies>\n";
	exit(-1);
}

get_chromosomes();
get_summary();
make_plots();
exit(0);



#
# find all the scoring summary files
#
sub get_summary_files
{
	my $finder_cl = "find $ARGV[0] -name \"*__sum.txt\" \| sort |";
	open(FINDER, $finder_cl);
	my $alldat="";
	while(my $line=<FINDER>){
		chomp($line);
		$alldat .= "$line ";
	}
	return $alldat;
}

#
# run the R script to plot assembly metrics
#
sub make_plots
{
	my $alldat = get_summary_files();
	$alldat =~ s/__sum.txt//g;
	
	# find the plotting script
	my $prefix = "";
	my $script_fullpath = File::Spec->rel2abs( $0 );
	my $script_dir = dirname($script_fullpath);
	my $plotscript = "assemblathonMauvePlots.R";
	my $testout = `$plotscript 2> /dev/null`;
	unless($testout =~ /Error/){
		$plotscript = "$script_dir/$plotscript" if( -e "$script_dir/$plotscript" );
		$plotscript = "./$plotscript" if( -e "./$plotscript" );
		die "Unable to find $plotscript\n" if($plotscript eq "assemblathonMauvePlots.R");
	}
	
	# run the plots
	my $plot_cl = "$plotscript $alldat";
	print "$plot_cl\n";
	`$plot_cl`;
}

#
# create a summaries.txt file
#
sub get_summary
{
	my $alldat = get_summary_files();
	my @files = split( /\s+/, $alldat );
	`rm -f summaries.txt`;
	my $firstfile = $files[0];
	`head -n 1 $firstfile > summaries.txt`;
	foreach my $file( @files ){
		`tail -n 1 $file >> summaries.txt`;
	};
}

sub get_chromosomes
{
	my $finder_cl = "find $ARGV[0] -name chromosomes.txt |";
	open(FINDER, $finder_cl);
	my $chr = <FINDER>;
	chomp $chr;
	`cp $chr .`;
}
