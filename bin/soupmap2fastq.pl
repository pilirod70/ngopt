#!/usr/bin/perl -w
# @author Aaron Darling
# Usage: soupmap2fastq.pl <read 1 mapping> <read 2 mapping> <output filtered map> <output filtered read pairs FastQ>
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 4) {
	print STDOUT " Usage: ".basename($0)." <read 1 mapping> <read 2 mapping> <output filtered map> <base output filtered read pairs FastQ>\n";
	exit 1;

}

open( READMAP, $ARGV[0] ) || die "Unable to open \"$ARGV[0]\" for reading\n";
open( READMAP2, $ARGV[1] ) || die "Unable to open \"$ARGV[1]\" for reading\n";
open( FILTEREDMAP, ">".$ARGV[2] ) || die "Unable to open \"$ARGV[2]\" for writing\n";
open( FILTEREDREADS1, ">".$ARGV[3]."_p1.fastq" ) || die "Unable to open \"$ARGV[3]\" for writing\n";
open( FILTEREDREADS2, ">".$ARGV[3]."_p2.fastq" ) || die "Unable to open \"$ARGV[3]\" for writing\n";

# parse pairs of high-quality mapped reads
my $keepstring = "49M";  # reads are 49nt long so 49M indicates complete mapping without indels
my %keepreads;
while( my $line = <READMAP> )
{
	if ($line =~ /^\@/) {
		print FILTEREDMAP $line;
		next;
	}
	my @mapline = split( /\s+/, $line );
	next unless $mapline[5] eq $keepstring;
	$keepreads{$mapline[0]} = $line;
}

while( my $line = <READMAP2> )
{
        next if ($line =~ /^\@/);
	my @mapline = split( /\s+/, $line );
	next unless $mapline[5] eq $keepstring;
	next unless defined($keepreads{$mapline[0]});
	my @mapline1 = split( /\s+/, $keepreads{$mapline[0]} );
        print FILTEREDREADS1 "\@$mapline1[0]/1\n";
	my $seq = $mapline1[9];
#	$seq = revcomp($seq) if($mapline1[1] == 16);
        print FILTEREDREADS1 "$seq\n";
        print FILTEREDREADS1 "+$mapline1[0]/1\n";
        print FILTEREDREADS1 "$mapline1[10]\n";
	$seq = $mapline[9];
#	$seq = revcomp($seq) if($mapline[1] == 16);
        print FILTEREDREADS2 "\@$mapline[0]/2\n";
        print FILTEREDREADS2 "$seq\n";
        print FILTEREDREADS2 "+$mapline[0]/2\n";
        print FILTEREDREADS2 "$mapline[10]\n";

	print FILTEREDMAP $keepreads{$mapline[0]};
        print FILTEREDMAP $line;
}
close FILTEREDREADS1;
close FILTEREDREADS2;


sub revcomp {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
} 
