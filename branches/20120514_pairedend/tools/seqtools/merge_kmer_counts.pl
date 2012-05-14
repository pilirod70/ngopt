#!/usr/bin/perl 
use strict;
use warnings;

die "Usage: merge_kmer_counts.pl <kmer_counts1> <kmer_counts2> <out_base>\n" if (@ARGV != 3);

my %count1 = ();
my %count2 = ();

my $file = shift;
open(IN,"<",$file);

my $total1 = 0;
my $nk1 = 0;
while(<IN>){
	chomp;
	my ($kmer, $cnt) = split ' ';
	$count1{$kmer} = $cnt;
	$total1 += $cnt;
	$nk1++;
}

$file = shift;
open(IN,"<",$file);

my $total2 = 0;
my $nk2 = 0;
while(<IN>){
	chomp;
	my ($kmer, $cnt) = split ' ';
	$count2{$kmer} = $cnt;
	$total2 += $cnt;
	$nk2++;
}

my $out_base = shift;
open(UNIQ,">","$out_base.uniq1");
open(SHARE,">","$out_base.shared");

my $uniq1 = 0;
my $uniq2 = 0;
for my $kmer (keys %count1) {
	if (defined($count2{$kmer})){
		$count2{$kmer} /= $total2;
		$count1{$kmer} /= $total1;
		print SHARE "$kmer\t".$count1{$kmer}."\t".$count2{$kmer}."\n"; 
	} else {
		$count1{$kmer} /= $total1;
		$uniq1++;
		print UNIQ "$kmer\t".$count1{$kmer}."\t0\n"; 
	}
}
close SHARE;
close UNIQ;
open(UNIQ,">","$out_base.uniq2");
for my $kmer (keys %count2){
	next if defined($count1{$kmer});
	$uniq2++;
	print UNIQ "$kmer\t0\t".$count2{$kmer}/$total2."\n"; 
}
close UNIQ;
print STDERR "set1: total=$nk1 unique=$uniq1\nset2: total=$nk2 unique=$uniq2\n";


