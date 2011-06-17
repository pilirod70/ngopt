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

my %blkbnds = ();
my $totblks = scalar(keys %blocks); 
my $sigblks=0;
# filter out short blocks
for my $block (keys %blocks) {
	my @ar = @{$blocks{$block}};
	#print STDERR "block $block has ".scalar(@ar)." points\n";
	if (scalar(@ar) < $min_blocks) {
		delete($blocks{$block});
	} else {
		$sigblks++;
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
			# find the left-most and right-most boundaries in the two contigs
			$ctg1_l = $pt1 if ($pt1 < $ctg1_l);
			$ctg1_r = $pt1 if ($pt1 > $ctg1_r);
			$ctg2_l = $pt2 if ($pt2 < $ctg2_l);
			$ctg2_r = $pt2 if ($pt2 > $ctg2_r);
		}
		my @tmp1 = ($ctg1_l,$ctg1_r);
		push(@{$blkbnds{$ctg1}},\@tmp1);

		my @tmp2 = ($ctg2_l,$ctg2_r);
		push(@{$blkbnds{$ctg2}},\@tmp2);
	}
}



print STDERR "Using $sigblks of $totblks blocks to break contigs\n";
#print STDERR scalar(keys %blkbnds)."\n";
#print STDERR join("\n",sort keys %blkbnds)."\n";

open(IN,"<",$ctgLbl_file);
my %contigs = ();

while(<IN>) {
	chomp;
	my ($id,$name) = split(/\t/,$_);
	$contigs{$name} = $id;
}

my %seqs = ();
open(IN,"<",$fasta_file);
my $ctg = "";
my $seq = "";
while (<IN>){
	chomp;
	if ($_ =~ m/^>(.*)/){
		#print STDERR ">".$ctg."<".$contigs{$ctg}."\n";
		unless (length($ctg) == 0) {
			$seqs{$ctg} = $seq;
			$seq = "";
		}
		$ctg = $1;
	} else {
		$seq .= $_; 
	}
	
}
$seqs{$ctg} = $seq;

#my %breaks = ();
#for $ctg 

for $ctg (keys %seqs) {
	if (defined($blkbnds{$contigs{$ctg}})){
		my @ar = @{$blkbnds{$contigs{$ctg}}};
		@ar = sort {$a->[0] <=> $b->[0]} @ar;
		my $start = 0;
		my $sub_seq = 0;
		$seq = $seqs{$ctg};
		for my $coord (@ar){
			my $l = $coord->[0];
			my $r = $coord->[1];
			print STDERR "Breaking $ctg at ($l,$r)\n";
			if ($l-100-$start > 100) {
				$sub_seq++;
				print ">$ctg\_$sub_seq\n";
				print substr($seq,$start,$l-100-$start)."\n";
			} 
			if ($r-$l > 100) {
				$sub_seq++;
				print ">$ctg\_$sub_seq\n";
				print substr($seq,$l,$r-$l)."\n";
			}
			$start = $r + 100; 
		}
	} else {
		print ">$ctg\n".$seqs{$ctg}."\n";
	}
}


#unless (length($ctg) == 0 || !defined($blkbnds{$contigs{$ctg}})) {
#	my @ar = @{$blkbnds{$contigs{$ctg}}};
#	@ar = sort {$a->[0] <=> $b->[0]} @ar;
#	my $start = 0;
#	my $sub_seq = 0;
#	for my $coord (@ar){
#		$sub_seq++;
#		my $l = $coord->[0];
#		my $r = $coord->[1];
#		print ">$ctg\_$sub_seq\n";
#		print substr($seq,$start,$l-100)."\n";
#		$sub_seq++;
#		print ">$ctg\_$sub_seq\n";
#		print substr($seq,$l,$r-$l+1)."\n";
#		$start = $r + 100; 
#	}
#}
#$ctg = $1;
#print $_."\n" unless(defined($blkbnds{$contigs{$ctg}}));

