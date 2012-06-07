#!/usr/bin/perl -w
use strict;

die "Usage <FastQ reads> <mate pair list> <output Paired FastQ>\n" unless @ARGV==3;

open( FASTQ, $ARGV[0] ) || die "Unable to read $ARGV[0]\n";
open( PAIRLIST, $ARGV[1] ) || die "Unable to read mate pair list $ARGV[1]\n";
open( QPAIRS, ">>$ARGV[2]" ) || die "Unable to write $ARGV[2]\n";

my %fastq;
my $curread = "";
my $curname = "";
print "Indexing FastQ\n";
while( my $line=<FASTQ> )
{
	if($line =~ /\@SOLEXA/){
		$fastq{$curname} = $curread;
		$curread = $line;
		chop $line;
		$line =~ s/\@SOLEXA/SOLEXA/g;
		$curname = $line;
	}else{
		$curread .= $line;
	}
}

print "Writing paired reads\n";
while( my $line = <PAIRLIST> ){
	$line =~ s/>//g;
	chop $line;
	my($r1,$r2) = split( /\t/, $line );
	print QPAIRS $fastq{$r1};
	print QPAIRS $fastq{$r2};
}
