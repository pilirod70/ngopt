#!/usr/bin/perl -w
# program to remove reads from a map that aren't in a FastQ file

die "Usage: <FastQ> <Input SHRiMP map> <Output filtered SHRiMP map>\n" unless @ARGV==3;

open(FASTQ, "$ARGV[0]") || die "Unable to open $ARGV[0]\n";
open(SHRIMPIN, "$ARGV[1]") || die "Unable to open $ARGV[1]\n";
open(SHRIMPOUT, ">$ARGV[2]") || die "Unable to open $ARGV[2]\n";
my %fqreads;


#read in all FastQ entries
print "Parsing FastQ file\n";
while( my $line = <FASTQ> ){
	next unless $line =~ /^\@SOLEXA/;
	chop $line;
	$line =~ s/^\@//g;
	$fqreads{$line}=1;
}
print "Writing filtered read map\n";

while( my $line2 = <SHRIMPIN> )
{
	my @l = split( /\t/, $line2 );
	$l[0] =~ s/>//g;
	print SHRIMPOUT $line2 if( defined( $fqreads{$l[0]} ) );
}

