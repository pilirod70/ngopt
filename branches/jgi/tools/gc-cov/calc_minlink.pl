#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
if (scalar(@ARGV)<2){
	print "Usage: ".basename($0)." <insert> <read_length> <sam_file>\nif <sam_file> missing, read from stdin\n";
	exit;
}
my $ins = shift;
my $rdlen = shift;

my $line = <>;
chomp $line;
my $genome_length=0;
while($line =~ m/^@.*LN:(\d+)/) {
	$genome_length += $1;
	$line = <>;
	chomp $line;
}
print STDERR "genome length: $genome_length\n";
my $nreads = 1;
while($line = <>){
	$nreads++;
}
print STDERR "# reads: $nreads\n";
my $minlink = $nreads*$ins/(2*$genome_length);
printf("%.0f\n", $minlink);

