#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV)!=7){
	print "Usage: ".basename($0)." <blocks_file> <contig_label_file> <fasta_file> <out_file> <min_pts> <min_len> <max_len>\n";
	exit;
}



my $blocks_file = shift;
my $ctgLbl_file = shift;
my $fasta_file = shift;
my $out_file = shift;
my $min_pts = shift;
my $min_len = shift;
my $max_len = shift;

my %blocks = ();

# read in our blocks
open(IN,"<",$blocks_file);
my $line = <IN>;
$line = <IN>;
while(!($line =~ m/^-by size/)){
	chomp $line;
	my @ar = split(/\t/,$line);
	my $block = $ar[2];
	if ($ar[5] =~ m/c(\d+)p(\d+)/){
		$ar[5] = $2;
	}
	if ($ar[6] =~ m/c(\d+)p(\d+)/){
		$ar[6] = $2;
	}
	push(@{$blocks{$block}},\@ar);
	$line = <IN>;
}
print STDERR "[a5_break] Breaking contings on blocks with $min_pts or more points\n";


my %blkbnds = ();
my $totblks = scalar(keys %blocks); 
my $sigblks=0;


for my $blk ( sort {$a <=> $b} (keys %blocks) ) {
	my $ctg1 = $blocks{$blk}[0][0];
	my $ctg2 = $blocks{$blk}[0][1];
	my $ctg1_r = 0;
	my $ctg1_l = 9**9**9;
	my $ctg2_r = 0;
	my $ctg2_l = 9**9**9;
	my $i = 0; 
	for my $pt (sort { $a->[4] <=> $b->[4]  } @{$blocks{$blk}}) {
		$i++;
		push (@$pt,$i);
		$ctg1_r = $pt->[5] if $pt->[5] > $ctg1_r;
		$ctg1_l = $pt->[5] if $pt->[5] < $ctg1_l;
	}
	$i = 0; 
	for my $pt (sort { $a->[5] <=> $b->[5]  } @{$blocks{$blk}}) {
		$i++;
		push (@$pt,$i);
		$ctg2_r = $pt->[6] if $pt->[6] > $ctg2_r;
		$ctg2_l = $pt->[6] if $pt->[6] < $ctg2_l;
	}
	my @x = ();
	my @y = ();
	for my $pt (@{$blocks{$blk}}){
		push(@x,$pt->[8]);
		push(@y,$pt->[9]);
	}
	my ($k,$z_k) = kendall(\@x,\@y);
	my $ctg1_len = ($ctg1_r-$ctg1_l);
	my $ctg2_len = ($ctg2_r-$ctg2_l);
	# skip block if stretch is too short, points aren't orderly enough (z_k), or not enough points
	next unless ( #abs($z_k) >= 2 &&
                  ($ctg1_len<=$max_len && $ctg1_len>=$min_len) && 
                  ($ctg2_len<=$max_len && $ctg2_len>=$min_len) && 
                  scalar(@x) >= $min_pts ); 
	$sigblks++;
	#printf STDERR "[a5_break] Block $blk: $ctg1|$ctg1_l-$ctg1_r  $ctg2|$ctg2_l-$ctg2_r t=%.4f Z_t=%.4f n=$min_pts\n",$k,$z_k;
	my @tmp1 = ($ctg1_l,$ctg1_r);
	push(@{$blkbnds{$ctg1}},@tmp1);

	my @tmp2 = ($ctg2_l,$ctg2_r);
	push(@{$blkbnds{$ctg2}},@tmp2);
}

for my $ctg (keys %blkbnds){
	my $aref = $blkbnds{$ctg};
	my @sorted = sort {$a <=> $b} @$aref;
	delete($blkbnds{$ctg});
	push(@{$blkbnds{$ctg}},@sorted);
}


print STDERR "[a5_break] Using $sigblks of $totblks blocks to break contigs\n";
if ($sigblks == 0) {
	print STDERR "[a5_break] No significant blocks found. Exiting and not breaking contigs.\n"; 
	print "$sigblks\n";
	exit;
}
print STDERR "[a5_break] ".scalar(keys %blkbnds)." contigs have blocks\n";

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

print STDERR "[a5_break] Found $nbases total bases\n";
my $currbases = 0;
open(OUT,">",$out_file);
for $ctg (keys %seqs) {
	if (!defined ($contigs{$ctg}) || !defined($blkbnds{$contigs{$ctg}})){
		print OUT ">$ctg\n".$seqs{$ctg}."\n";
		$currbases += length($seqs{$ctg});
		next;
	}
	if (defined($blkbnds{$contigs{$ctg}})){
		my @ar = @{$blkbnds{$contigs{$ctg}}};
		my $sub_seq = 0;
		$seq = $seqs{$ctg};
		if ($ar[0] >= 100) {
			$sub_seq++;
			print STDERR "[a5_break] Breaking $ctg at $ar[0]\n";
			print OUT ">$ctg|$sub_seq|1-$ar[0]\n";
			print OUT substr($seq,0,$ar[0])."\n";
			$currbases += $ar[0];
		} else {
			print STDERR "[a5_break] Discarding region >$ctg|1-$ar[0]\n";
		}	
		for (my $i = 1; $i < scalar(@ar); $i++){
			my $l = $ar[$i-1]+1;
			my $r = $ar[$i];
			if ($r - $l + 1 >= 100){
				$sub_seq++;
				print STDERR "[a5_break] Breaking $ctg at $r\n";
				print OUT ">$ctg|$sub_seq|$l-$r\n";
				print OUT substr($seq,$l-1,$r-$l+1)."\n";
				$currbases += $r-$l+1;
			} else {
				print STDERR "[a5_break] Discarding region >$ctg|$l-$r\n";
			}
		}
		my $l = $ar[scalar(@ar)-1]+1;
		my $r = length($seq);
		if ($r - $l + 1 >= 100){
			$sub_seq++;
			print STDERR "[a5_break] Breaking $ctg at $l\n";
			print OUT ">$ctg|$sub_seq|$l-$r\n";
			print OUT substr($seq,$l-1,$r-$l+1)."\n";
			$currbases += $r-$l+1;
		} else {
			print STDERR "[a5_break] Discarding region >$ctg|$l-$r\n";
		}
		
	}
}

close OUT;

print STDERR "[a5_break] Printed $currbases total bases \n";
print "$sigblks\n";

sub kendall {
	my $x = shift;
	my $y = shift;
	my $n = scalar(@$x);
	my $numer = 0;

	for (my $i = 0; $i < $n; $i++){
		for (my $j = 0; $j < $i; $j++){
			$numer += ($x->[$i]-$x->[$j])*($y->[$i]-$y->[$j]) > 0 ? 1 : -1;
		}
	}
	my $z_A = 3*($numer)/(($n*($n-1)*(2*$n+5)/2)**(1/2));
	return (2*($numer)/($n*($n-1)), $z_A);
}
