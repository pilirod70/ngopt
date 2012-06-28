#!/usr/bin/perl -w
use strict;
use warnings;

my $winlen = 20;

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
	next unless $line1 =~ /XT:A:U/;
	next unless $line2 =~ /XT:A:U/;
	$pos1 = (split(/\t/,$line1))[3];
	$pos2 = (split(/\t/,$line2))[3];
	tally($pos1,$pos2) if ($pos1 != 0 && $pos2 != 0);
}

for my $ar (@bins) {
	print $ar->[0]."\t".$ar->[1]."\t".$ar->[2]."\n";
}

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
