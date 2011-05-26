#!/usr/bin/env perl
use strict;
use warnings;

my $r1fq = $ARGV[0];
my $r2fq = $ARGV[1];
my $barcodes = $ARGV[2];
my $prefix = $ARGV[3];

#`splitBCpe $r1fq $r2fq --bcfile $barcodes --bol --prefix $prefix`;

open( BC, $barcodes );
while( my $line = <BC> ){
	my ($name, $oligo) = split( /\t/, $line );
	if($r2fq =~ /_3_sequence.txt/){
		# need to replace /3 with /2
		my $pie = "perl -p -i -e \"s/\\\/3\$/\\\/2/g\" $prefix$name"."_p2";
		print $pie."\n";
#		`$pie`;
	}
	my $assemble = "qsub -q all.q -q eisen.q /home/koadman/bin/assemble_sga_idba_sspace.sh $prefix$name"."_p1 $prefix$name"."_p2 $prefix$name"."_up $prefix$name";
	print $assemble."\n";
	`$assemble`;
}

