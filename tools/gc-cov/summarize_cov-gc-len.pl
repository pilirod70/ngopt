#!/usr/bin/perl -w
use strict;
use warnings;
sub mean(\@);
sub var(\@$);
print "Usage: $0 <tabd_file>\n" if (scalar(@ARGV)!=1);
exit  if (scalar(@ARGV)!=1);

my $file = shift;
open(IN,"<",$file);
my @gc;
my @cov;
my @len;

while(<IN>){
	chomp;
	my @ctg = split /\t/;
	push (@gc,$ctg[1]);
	push (@cov,$ctg[2]);
	push (@len,$ctg[3]);
}
my $gc_ = mean(@gc);
my $cov_ = mean(@cov);
my $len_ = mean(@len);

print $file."\t".$gc_."\t".var(@gc,$gc_)."\t".$cov_."\t".var(@cov,$cov_)."\t".$len_."\t".var(@len,$len_)."\n";

sub mean(\@) {
	my $vec = shift;
	my $sum = 0;
	for my $i (@$vec) {
		$sum += $i;
	}
	return $sum/scalar(@$vec);
}

sub var(\@$) {
	my $vec = shift;
	my $mean = shift;
	my $sum = 0;
	for my $i (@$vec) {
		$sum += ($mean-$i)**2;
	}
	return $sum/scalar(@$vec);
}
