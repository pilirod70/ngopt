#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV < 1) {
	print "Usage: subsample_paired_sam.pl <fraction> <sam|stdin>\nOutput is printed to stdout\n";
	exit;
}

my $frac = shift;
if (!($frac > 0 && $frac < 1)){
	print "<fraction> must be greater than 0 and less than 1\n";
	exit;
}

my $r1;
my $r2;
while($r1 = <>){
	if ($r1 =~ /^\@SQ/) {
		print $r1;
	} else {
		$r2 = <>;
		if (rand() < $frac){
			print $r1.$r2;
		}
	}
}
