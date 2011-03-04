#!/usr/bin/perl -w 

use strict;
use warnings;

if (scalar(@ARGV) != 2) {
	print "Usage: splitCV <core_size> <ortholog_td_file>\n";
	exit 1;
}

my $num_core;
my $ortho_file;

$num_core = $ARGV[0];
$ortho_file = $ARGV[1];

open (OF,"<",$ortho_file);
open (CORE,">",$ortho_file.".core");
open (VAR,">",$ortho_file.".var");
while (<OF>) {
	chomp;
	my @orthos = split;
	if (scalar(@orthos) == $num_core) {
		print CORE $_."\n";		
	} else {
		print VAR $_."\n";
	}
}

close CORE;
close VAR;

