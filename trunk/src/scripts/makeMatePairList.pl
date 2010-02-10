#!/usr/bin/perl -w
use strict;
die "Usage: <output directory> <FastQ reads> <insert size> <insert sd> <total coverage> <read len> <mapping chr 1> <chr 1 len> ... <mapping chr n> <chr n len>\n" unless @ARGV>=6;

my $outputdir = $ARGV[0];
my $fastq = $ARGV[1];
my $totalcov = $ARGV[4];
exit if $totalcov == 0;  # nothing to see here, move along
my $readlen = $ARGV[5];
my $totallen = 0;
for( my $i = 7; $i < @ARGV; $i+=2 ){
	$totallen += $ARGV[$i];
}
my $allreads=0;
my @chr_reads;
for( my $i = 6; $i < @ARGV; $i+=2 ){
	my $prev="";
	$chr_reads[$i]=0;
	open( CURMAP, $ARGV[$i] );
	while( my $line = <CURMAP> ){
		my @f = split( /\t/, $line );
		next if( $prev eq $f[0] );
		$prev = $f[0];
		$chr_reads[$i]++;
	}
	close CURMAP;
	$allreads += $chr_reads[$i];
}

print "Found $allreads reads to estimate expected relative coverage\n";
`mkdir -p $outputdir`;
my $readcount = $totallen*$totalcov/($readlen*2);
for( my $i = 6; $i < @ARGV; $i+=2 ){
	my $chr_rc = int($readcount * ($chr_reads[$i]/$allreads));
	my $cl = "createMatePairs.pl $ARGV[$i] $ARGV[$i+1] $chr_rc $ARGV[2] $ARGV[3] > $outputdir/$i.pairs";
	print "$cl\n";
	`$cl`;
	my $cl2 = "makeFastQpairs.pl $fastq $outputdir/$i.pairs $outputdir/simfastq";
	`$cl2`;
}

