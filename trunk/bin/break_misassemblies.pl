#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV)!=3){
	print "Usage: ".basename($0)." <blocks_file> <contig_label_file> <fasta_file>\n";
	exit;
}



my $blocks_file = shift;
my $ctgLbl_file = shift;
my $fasta_file = shift;

my %blocks = ();

# read in our blocks
open(IN,"<",$blocks_file);
my $line = <IN>;
$line = <IN>;
while(!($line =~ m/^-by size/)){
	chomp $line;
	my @ar = split(/\t/,$line);
	my $block = $ar[2];
	push(@{$blocks{$block}},\@ar);
	$line = <IN>;
}
$line = <IN>;
my @count = ();
my @tmp = ();
my $total = 0;
my $min_pts = 0;
while($line = <IN>){
	chomp $line;
	@tmp = split(' ',$line);
	push(@count,[ @tmp ]);
	$total += $tmp[1];
	if ($tmp[3] == 0 && !$min_pts){
		$min_pts = $tmp[0];
	}
}


#my $p = 0;
#my $alpha = 1 - $tail_prob;
#for (my $i = 0; $i < scalar(@count); $i++){
#	$count[$i][1] /= $total;
#	$p+=$count[$i][1];
#	if ($p > $alpha){
#		$min_pts = $count[$i][0];
#		last;
#	}
#}
print STDERR "[a5_break] Breaking contings on blocks with $min_pts or more points\n";


my %blkbnds = ();
my $totblks = scalar(keys %blocks); 
my $sigblks=0;

my %cncts = ();

# filter out short blocks
for my $block ( sort {$a <=> $b} (keys %blocks) ) {
	my @ar = @{$blocks{$block}};
	#print STDERR "block $block has ".scalar(@ar)." points\n";
	if (scalar(@ar) < $min_pts) {
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
		my $first =  1;
		for my $point (@ar) {
			$ctg1 = $point->[0];
			$ctg2 = $point->[1];
			if ($first) {
				$cncts{$ctg1} = () unless defined($cncts{$ctg1});
				$cncts{$ctg1}{$ctg2} = 0 unless defined($cncts{$ctg1}{$ctg2});
				$cncts{$ctg1}{$ctg2}++;	
				$first = 0;
			}
			my @pt1 = split(/,/,$point->[5]);
			for my $pt (@pt1){
				if ($pt =~ m/c(\d+)p(\d+)/) {
					$pt = $2;
				} 
				# find the left-most and right-most boundaries in the first contig
				$ctg1_l = $pt if ($pt < $ctg1_l);
				$ctg1_r = $pt if ($pt > $ctg1_r);
			}
			my @pt2 = split(/,/,$point->[6]);
			for my $pt (@pt2){
				if ($pt =~ m/c(\d+)p(\d+)/) {
					$pt = $2;
				} 
				# find the left-most and right-most boundaries in the second contig
				$ctg2_l = $pt if ($pt < $ctg2_l);
				$ctg2_r = $pt if ($pt > $ctg2_r);
			}
		}
		print STDERR "Block $block: $ctg1|$ctg1_l-$ctg1_r  $ctg2|$ctg2_l-$ctg2_r\n";
		my @tmp1 = ($ctg1_l,$ctg1_r);
		push(@{$blkbnds{$ctg1}},@tmp1);

		my @tmp2 = ($ctg2_l,$ctg2_r);
		push(@{$blkbnds{$ctg2}},@tmp2);
	}
}

for my $c1 (sort keys %cncts){
	for my $c2 (sort keys %{$cncts{$c1}}) {
		print STDERR "$c1\t$c2\t".$cncts{$c1}{$c2}."\n";
	}

}

for my $ctg (keys %blkbnds){
	my $aref = $blkbnds{$ctg};
	my @sorted = sort {$a <=> $b} @$aref;
	delete($blkbnds{$ctg});
	push(@{$blkbnds{$ctg}},@sorted);
}


print STDERR "Using $sigblks of $totblks blocks to break contigs\n";
print STDERR scalar(keys %blkbnds)." contigs have blocks\n";

# make a map of contig-names to contig-ids
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
my $nbases = 0;
while (<IN>){
	chomp;
	if ($_ =~ m/^>(.*)/){
		unless (length($ctg) == 0) {
			$seqs{$ctg} = $seq;
			$seq = "";
		}
		$ctg = $1;
	} else {
		$nbases += length($_);
		$seq .= $_; 
	}
	
}
$seqs{$ctg} = $seq;

print STDERR "Found $nbases total bases\n";
my $currbases = 0;
for $ctg (keys %seqs) {
	if (!defined ($contigs{$ctg}) || !defined($blkbnds{$contigs{$ctg}})){
		print ">$ctg\n".$seqs{$ctg}."\n";
		$currbases += length($seqs{$ctg});
		next;
	}
	if (defined($blkbnds{$contigs{$ctg}})){
		my @ar = @{$blkbnds{$contigs{$ctg}}};
		my $sub_seq = 0;
		$seq = $seqs{$ctg};
		if ($ar[0] >= 100) {
			$sub_seq++;
			print STDERR "Breaking $ctg at $ar[0]\n";
			print ">$ctg|$sub_seq|1-$ar[0]\n";
			print substr($seq,0,$ar[0])."\n";
			$currbases += $ar[0];
		} else {
			print STDERR "Discarding region >$ctg|1-$ar[0]\n";
		}	
		for (my $i = 1; $i < scalar(@ar); $i++){
			my $l = $ar[$i-1]+1;
			my $r = $ar[$i];
			if ($r - $l + 1 >= 100){
				$sub_seq++;
				print STDERR "Breaking $ctg at $r\n";
				print ">$ctg|$sub_seq|$l-$r\n";
				print substr($seq,$l-1,$r-$l+1)."\n";
				$currbases += $r-$l+1;
			} else {
				print STDERR "Discarding region >$ctg|$l-$r\n";
			}
		}
		my $l = $ar[scalar(@ar)-1]+1;
		my $r = length($seq);
		if ($r - $l + 1 >= 100){
			$sub_seq++;
			print STDERR "Breaking $ctg at $l\n";
			print ">$ctg|$sub_seq|$l-$r\n";
			print substr($seq,$l-1,$r-$l+1)."\n";
			$currbases += $r-$l+1;
		} else {
			print STDERR "Discarding region >$ctg|$l-$r\n";
		}
		
	}
}
print STDERR "Printed $currbases total bases \n";

