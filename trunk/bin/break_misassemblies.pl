#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV)!=3){
	print "Usage: ".basename($0)." <blocks_file> <contig_label_file> <fasta_file>\n";
	exit;
}


my $min_blocks = 10;

my $blocks_file = shift;
my $ctgLbl_file = shift;
my $fasta_file = shift;

my %blocks = ();

# read in our blocks
open(IN,"<",$blocks_file);
<IN>;
while(<IN>){
	chomp;
	my @ar = split(/\t/,$_);
	my $block = $ar[2];
	push(@{$blocks{$block}},\@ar);
}

my %breaks = ();
print STDERR scalar(keys %blocks)." blocks \n";

# filter out short blocks
for my $block (keys %blocks) {
	my @ar = @{$blocks{$block}};
	print STDERR "block $block has ".scalar(@ar)." points\n";
	if (scalar(@ar) < $min_blocks) {
		delete($blocks{$block});
	} else {
		# calcate the boundaries of the segment we need to remove	
		my $ctg1;
		my $ctg2;
		my $ctg1_l = 9**9**9;
		my $ctg1_r = 0;
		my $ctg2_l = 9**9**9;
		my $ctg2_r = 0;
		for my $point (@ar) {
			$ctg1 = $point->[0];
			$ctg2 = $point->[1];
			my $pt1 = $point->[5];
			if ($pt1 =~ m/c(\d+)p(\d+)/) {
				$pt1 = $2;
			} 
			my $pt2 = $point->[6];
			if ($pt2 =~ m/c(\d+)p(\d+)/) {
				$pt2 = $2;
			} 
			$ctg1_l = $pt1 if ($pt1 < $ctg1_l);
			$ctg1_r = $pt1 if ($pt1 > $ctg1_r);
			$ctg2_l = $pt2 if ($pt2 < $ctg2_l);
			$ctg2_r = $pt2 if ($pt2 > $ctg2_r);
		}
		my @tmp1 = ($ctg1_l,$ctg1_r);
		push(@{$breaks{$ctg1}},\@tmp1);

		my @tmp2 = ($ctg2_l,$ctg2_r);
		push(@{$breaks{$ctg2}},\@tmp2);
	}
}
#print STDERR scalar(keys %breaks)."\n";
#print STDERR join("\n",sort keys %breaks)."\n";

open(IN,"<",$ctgLbl_file);
my %contigs = ();

while(<IN>) {
	chomp;
	my ($id,$name) = split(/\t/,$_);
	$contigs{$name} = $id;
}

open(IN,"<",$fasta_file);
my $ctg = "";
my $seq = "";
while (<IN>){
	chomp;
	if ($_ =~ m/^>(.*)/){
#		print STDERR ">".$ctg."<".$contigs{$ctg}."\n";
		unless (length($ctg) == 0 || !defined($breaks{$contigs{$ctg}})) {
			my @ar = @{$breaks{$contigs{$ctg}}};
			@ar = sort {$a->[0] <=> $b->[0]} @ar;
			my $offset = 0;
			my $sub_seq = 0;
			for my $coord (@ar){
				print STDERR "\$coord length = ".scalar(@$coord)."\n";
				$sub_seq++;
				my $l = $coord->[0];
				my $r = $coord->[1];
				print ">$ctg\_$sub_seq\n";
				print substr($seq,$offset,$l-100)."\n";
				$sub_seq++;
				print ">$ctg\_$sub_seq\n";
				print substr($seq,$l,$r-$l+1)."\n";
				$offset = $r + 100; 
			}
		}
		$ctg = $1;
		print $_."\n" unless(defined($breaks{$contigs{$ctg}}));
	} else {
		if (defined($breaks{$contigs{$ctg}})){
			$seq .= $_; 
		} else {
			print $_."\n";
		}
	}
	
}


