#!/usr/bin/perl -w
use strict;
use warnings;
sub muvar(\@);
print "Usage: $0 <tabd_file(s)>\n" if (scalar(@ARGV) < 2);
exit  if (scalar(@ARGV)<2);

while (@ARGV) {
	my $file = shift;
	open(IN,"<",$file);
	my (@gc, @cov, @len);
	while(<IN>){
		chomp;
		my @ctg = split /\t/;
		push (@gc,$ctg[1]);
		push (@cov,$ctg[2]);
		push (@len,$ctg[3]);
	}
	print "$file\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", (muvar(@gc),muvar(@cov),muvar(@len));
}	
sub muvar(\@) {
	my $vec = shift;
	my $sum = 0;
	for my $i (@$vec) {
		$sum += $i;
	}
	my $mean = $sum/scalar(@$vec);
	$sum = 0;
	for my $i (@$vec) {
		$sum += ($mean-$i)**2;
	}
	return ($mean, $sum/scalar(@$vec));
}
