#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 3){
	print "Usage ".basename($0)." <min_len> <min_%id> <tabd_blast_file> \n";
	exit 1;
}
my $num_lines = 0;
my $min_aln = shift;
my $min_pid = shift;
my $file = shift;
open(IN,"<",$file);
print STDERR "Removing blast hits less than $min_aln bp\n";
#    node0_0_0_93995	node0_0_0_93995	100.00	3815	0	0	3190	7004	3190	7004	0.0	7563
#    node0_0_0_93995	node0_0_0_93995	100.00	2437	0	0	1	2437	1	2437	0.0	4831
my $num_disc = 0;
my $num_reflexive = 0;
my $num_too_short = 0;
my $num_too_diff = 0;
my %hits = ();
while(<IN>) {
	chomp;
	my @hit = split;
	$num_lines++;
	if ($hit[3] >= $min_aln) {
		if ($hit[0] eq $hit[1]) {
			if ($hit[6] != $hit[8] || $hit[7] != $hit[9]) {
			#	print $_."\n";
				add_hit(\@hit);
			} else {
				$num_reflexive++;
			}	
		} else {
			if ($hit[2] >= $min_pid) {
			#	print $_."\n";
				add_hit(\@hit);
			} else {
				$num_too_diff++;
			}
		}
	} else { 
		$num_too_short++;
	}
}

foreach my $ctg1 (keys %hits) {
	foreach my $ctg2 (keys %{$hits{$ctg1}}){
		foreach my $hit (@{$hits{$ctg1}{$ctg2}}){
			print STDOUT $ctg1."\t".$ctg2."\t".join("\t",@$hit)."\n";
		}
	}
} 

print STDERR "Removed:\n";
print STDERR "\t$num_reflexive reflexive hits\n";
print STDERR "\t$num_too_short hits with less than $min_aln aligned base pairs\n";
print STDERR "\t$num_too_diff hits with less than $min_pid percent identity\n";

sub add_hit {
	my ($hit) = shift;
	print STDERR "size = ".scalar(@$hit)." at line $num_lines\n";
	my $ctg1 = shift @$hit;
	my $ctg2 = shift @$hit;
	if ($ctg2 lt $ctg1) {
		my $tmp = $ctg2;
		$ctg2 = $ctg1;
		$ctg1 = $ctg2;
		# swap 4,5 with 6,7
		$tmp = $hit->[4];
		$hit->[4] = $hit->[6];
		$hit->[6] = $tmp;
		$tmp = $hit->[5];
		$hit->[5] = $hit->[7];
		$hit->[7] = $tmp;
	}
	if (defined($hits{$ctg1})){
		if (defined($hits{$ctg1}{$ctg2})){
#			print STDERR "$ctg1 and $ctg2 have ".scalar(@{$hits{$ctg1}{$ctg2}})." matches so far\n";	
			my $dup = 0;
			foreach my $tmp_hit (@{$hits{$ctg1}{$ctg2}}) {
#				print STDERR "1:length of hit: ".scalar(@$tmp_hit)."\n"; 
#				print STDERR "2:length of hit: ".scalar($tmp_hit)."\n"; 
#				print STDERR "hit = $tmp_hit\n"; 
				if ($hit->[4] == $tmp_hit->[4] && $hit->[6] == $tmp_hit->[6]){
					$dup = 1;
					last;
				}
			}
			if (!$dup){
				push (@{$hits{$ctg1}{$ctg2}}, $hit);
			}
		} else {
			$hits{$ctg1}{$ctg2} = ();
			push (@{$hits{$ctg1}{$ctg2}}, $hit);
		}
	} else {
		$hits{$ctg1} = ();
		$hits{$ctg1}{$ctg2} = ();
		push (@{$hits{$ctg1}{$ctg2}}, $hit);
	}
}


