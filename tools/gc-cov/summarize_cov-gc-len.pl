#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
sub muvar (\@);
print "Usage: ".basename($0)." <tabd_file(s)>\n" if (scalar(@ARGV) < 2);
exit  if (scalar(@ARGV)<2);

while (@ARGV) {
	my $file = shift;
	open(IN,"<",$file);
	<IN>;
	my (@gc, @cov, @len);
	while(<IN>){
		chomp;
		my @ctg = split /\t/;
		push (@gc,$ctg[1]);
		push (@cov,$ctg[2]);
		push (@len,$ctg[3]);
	}
	printf "$file\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", (muvar(@gc),muvar(@cov),muvar(@len));
}	
sub muvar (\@){
	my $vec = shift;
	my $sum = 0;
	#print STDERR "length: ".scalar(@$vec)."\n";
	print STDERR @$vec."\n" if (scalar(@$vec) == 1);
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
