#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my $winlen = 20;
my $uniq = 0;
my $help = 0;
if (!GetOptions( 
	"win=i"		=> \$winlen,
	"uniq"		=> \$uniq,
	"h"			=> \$help
    )
   ){
	printHelp()
}
printHelp() if ($help);
if (@ARGV == 1 && ! -f $ARGV[0]){
	print $ARGV[0]." not a valid file\n";
	printHelp();
}

my ($line1,$line2,$pos1,$pos2,$str1,$str2);

$line1 = <>;
my $len = substr($line1,index($line1,"LN:")+3);
my $num_win = (int$len/$winlen)+1;

my @bins = ();
for (my $i = 0; $i < $num_win; $i++) {
	$bins[$i][0] = ($i+1)*$winlen;
	$bins[$i][1] = 0;
	$bins[$i][2] = 0;
}

while($line1 = <>){
	next if ($line1 =~ /^@/);
	$line2 = <>;
	if ($uniq) {
		next unless $line1 =~ /XT:A:U/;
		next unless $line2 =~ /XT:A:U/;
	}
	$pos1 = (split(/\t/,$line1))[3];
	$pos2 = (split(/\t/,$line2))[3];
	tally($pos1,$pos2) if ($pos1 != 0 && $pos2 != 0);
}

my $mean;
my $MoM; # the mean of means
my $SoM; # the stdev of means
for my $ar (@bins) {
	if ($ar->[2] == 0) {
		$mean = 0;
	} else {
		$mean = $ar->[1]/$ar->[2];
	}
	print $ar->[0]."\t".$mean."\t".$ar->[2]."\n";
	$MoM += $mean;
	$SoM += $mean**2;
}
$MoM = $MoM/scalar(@bins);
$SoM = ($SoM/scalar(@bins) - $MoM**2)**(1/2);
printf STDERR "%.3f,%.3f\n", $MoM, $SoM; 

sub tally {
	my $r1 = shift;
	my $r2 = shift;
	my $dist = abs($r1-$r2);
	my $idx = getIndex($r1);
	$bins[$idx][1] += $dist;
	$bins[$idx][2]++;
	$idx = getIndex($r2);
	$bins[$idx][1] += $dist;
	$bins[$idx][2]++;
}

sub getIndex {
	my $sorted_idx = 0;
	my $q = shift;
	my $max = scalar(@bins)-1;
	my $min = 0;
	my $mid;
	$mid = int(($max+$min)/2);
	while ($max > $min) {
		if ($q < $bins[$mid][$sorted_idx]) {
			$max = $mid;
		} elsif ($q > $bins[$mid][$sorted_idx]) {
			$min = $mid+1;
		} else {
			 return $mid;
		}
		$mid = int(($max+$min)/2);
	}
	return $mid;
}

sub printHelp {
	print "Usage: getDistHist.pl [options] <sam|stdin>\n".
	      " where options are:\n\n".
	      "    -win       the window length to use [20]\n".
	      "    -uniq      keep only uniquely mapping reads\n".
	      "    -h         print this message and exit\n".
	      "\n".
	      "Compute the average distance to mates of reads mapping\n".
	      "in windows across a genome\n";
	exit;
}
