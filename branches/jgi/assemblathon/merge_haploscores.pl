#!/usr/bin/env perl
# warning: this code has not been thoroughly tested

use strict;
use warnings;
my $haplo1file = $ARGV[0];
my $haplo2file = $ARGV[1];
open( H1FILE, $haplo1file ) || die "Unable to read \"$haplo1file\"\n";
open( H2FILE, $haplo2file ) || die "Unable to read \"$haplo2file\"\n";

my %found1;
my %found2;
while( my $line = <H1FILE> ){
	my @blob = split( /\s+/, $line );
	next if $blob[0] eq "SNP_Pattern";
	$found1{$blob[6]}=$line;
}

while( my $line = <H2FILE> ){
	my @blob = split( /\s+/, $line );
	$found2{$blob[6]}=$line;
}

foreach my $key(keys(%found1)){
	next unless defined( $found2{$key} );
	# this base was miscalled against both haplotypes
	print $found1{$key};
}

