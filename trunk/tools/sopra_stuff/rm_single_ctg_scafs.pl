#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename; 

if (@ARGV != 2) {
	print STDOUT "Usage: ".basename($0)." <ctginscaf.txt> <contigs.fasta>\n";
	exit 1;
}
my $debug = 0;
my $cis_file = $ARGV[0];
my $fasta_file = $ARGV[1];
my %ctg_in_scaf = ();
my %keepers = ();

open (CIS,"<",$cis_file);
while (<CIS>) {
	my ($ctg,$scaf) = split;
	print STDOUT ">$ctg< : >$scaf<\n" if $debug;
	push(@{$ctg_in_scaf{$scaf}}, $ctg); 
}
my $num_singletons = 0;
foreach my $scaf (keys %ctg_in_scaf) {	
	my $keep = (scalar($ctg_in_scaf{$scaf}) > 1); 
	foreach my $ctg (@{$ctg_in_scaf{$scaf}}){
		print STDERR "1: ".$ctg."\n" if $debug;
		$keepers{$ctg} = scalar(@{$ctg_in_scaf{$scaf}});
		if ($keepers{$ctg} == 1){
			$num_singletons++;
		}
	}	
	print STDERR "$scaf has ".scalar(@{$ctg_in_scaf{$scaf}})." contigs\n" if $debug;
}
print STDOUT "Found $num_singletons singletons.\n";

open (FA,"<",$fasta_file);
my $base = basename($fasta_file,".fasta");
my $dir = dirname($fasta_file);
open (KEEP,">",$dir."/".$base.".scaffolded.fasta");
open (SING,">",$dir."/".$base.".singletons.fasta");
my $keep = 0;
my $num_rm = 0;
my $num_kept = 0;
my $num_extra = 0;
my $num_ctgs = 0;
my @tmp; 
while (<FA>) {
	chomp;
    if (/^>(\S+)/) {
		@tmp = split /\|/,$1;
		if (defined $keepers{$tmp[0]}) {
			print STDOUT "Found contig ".$tmp[0]."\n" if $debug;
			if ($keepers{$tmp[0]}>1) {
				$keep = 1;
				print KEEP ">$1\n";
				delete($keepers{$tmp[0]});
				$num_kept++;
			} else {
				$keep = 0;
				print SING ">$1\n";
				delete($keepers{$tmp[0]});
				$num_rm++;
			}
		} else {
			print STDERR "contig not found: ".$tmp[0]."\n";
			$keep = 0;
			$num_rm++;
			$num_extra++;
		}
		$num_ctgs++;
    } else {
		if ($keep) {
			print KEEP "$_\n";
		} else {
			print SING "$_\n";
		}
    }
}
print STDOUT "Wrote $num_kept of $num_ctgs contigs to scaffolded file\n";
print STDOUT "Discarded $num_rm contigs in total, $num_extra of which were not found in <ctginscaf.txt> file\n";
$num_extra = scalar(keys %keepers);
print STDOUT "Found $num_extra extra contigs in <ctginscaf.txt> file\n";
close KEEP;
close SING;

