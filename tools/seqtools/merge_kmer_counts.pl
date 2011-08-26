#!/usr/bin/perl 
use strict;
use warnings;

die "Usage: merge_kmer_counts.pl <kmer_counts1> <kmer_counts2>\n" if (@ARGV != 2);

my %count1 = ();
my %count2 = ();

my $file = shift;
open(IN,"<",$file);

my $total1 = 0;
while(<IN>){
	chomp;
	my ($kmer, $cnt) = split ' ';
	$count1{$kmer} = $cnt;
	$total1 += $cnt;
}

$file = shift;
open(IN,"<",$file);

my $total2 = 0;
while(<IN>){
	chomp;
	my ($kmer, $cnt) = split ' ';
	$count2{$kmer} = $cnt;
	$total2 += $cnt;
}

for my $kmer (keys %count1) {
	if (defined($count2{$kmer})){
		$count2{$kmer} /= $total2;
		$count1{$kmer} /= $total1;
		print "$kmer\t".$count1{$kmer}."\t".$count2{$kmer}."\n"; 
		delete($count2{$kmer});
	} else {
		$count1{$kmer} /= $total1;
		print "$kmer\t".$count1{$kmer}."\t0\n"; 
	}
	delete($count1{$kmer});
}

for my $kmer (keys %count2){
	print "$kmer\t0\t".$count2{$kmer}."\n"; 
}
