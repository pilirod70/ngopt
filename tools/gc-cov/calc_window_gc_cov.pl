#!/usr/bin/perl -w 
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 2) {
	print "Usage: ".basename($0)." <window_size> <pileup_file> \n";
	print "set <window_size> to -1 to calculate per contig\n";
	exit 1;
}

my $win = shift;
my $plp_file = shift;

my $curr_len = 0;
my $curr_cov = 0;
my $curr_gc = 0;
open (PLP,"<",$plp_file);

#node0_0_0_93995	12	C	1	^FT	~
#node0_0_0_93995	13	G	1	.	~
my $tmp = <PLP>;
my @line = split(/\t/,$tmp);
my $curr_ctg = $line[0];
$curr_cov += $line[1];
$curr_gc++ if (lc($line[2]) eq "c" || lc($line[2]) eq "g");
$curr_len++;
print "name\t" if ($win == -1);
print "cov\tgc\tlen\n";
while (<PLP>) {
	chomp;
	@line = split /\t/;
	if ($line[0] eq $curr_ctg) {
		$curr_gc++ if (lc($line[2]) eq "c" || lc($line[2]) eq "g");
		$curr_len++;
		$curr_cov+= $line[3];
		if ($curr_len == $win) {
			$curr_gc /= $curr_len;
			$curr_cov /= $curr_len;
			print "$curr_cov\t$curr_gc\t$curr_len\n";
			$curr_len = 0.0;
			$curr_gc = 0.0;
			$curr_cov = 0.0;
		}
	} else {
		if ($curr_len > 0.0) {
			$curr_gc /= $curr_len;
			$curr_cov /= $curr_len;
			print "$curr_ctg\t"  if ($win == -1);
			print "$curr_cov\t$curr_gc\t$curr_len\n";
		}
		$curr_ctg = $line[0];
		$curr_len = 1.0;
		$curr_cov = 1.0;
		if (lc($line[2]) eq "c" || lc($line[2]) eq "g"){
			$curr_gc = 1.0;
		} else {
			$curr_gc = 0.0;
		}
	}
}

if ($curr_len != 0.0) {
	$curr_gc /= $curr_len;
	$curr_cov /= $curr_len;
	print "$curr_ctg\t"  if ($win == -1);
	print "$curr_cov\t$curr_gc\t$curr_len\n";
}
