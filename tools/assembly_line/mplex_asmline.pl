#!/usr/bin/perl -w

use strict;
use warnings;


my $conf_file;
my $adptrs;
my $dust = 1;
my $paired = 0;

my @base;
my $base_idx=0;
my @qual_cut;
my $qual_idx=2;
my @mink;
my $mink_idx=3;
my @maxk;
my $maxk_idx=4;
my @kfreq_cut=5;
my $kcut_idx=5;
my @ins;
my $ins_idx=6;
my @min_pairs;
my $mnpr_idx=7
my @max_maphits;
my $mxhit_idx=8;
my @min_ctglen;
my $mnlen_idx=9
my @outdir;
my $out_idx=10;

my $CONF;
open CONF, "< $conf_file" or die "Can't open conf file\n";

my @args;
my $i=0;
while (<CONF>) {
	$args = split(' ',$_);
	$base[$i] = $args[$base_idx];
	$qual_cut[$i] = $args[$qual_idx];
	$mink[$i] = $args[$mink_idx];
	$maxk[$i] = $args[$maxk_idx];
	$kfreq_cut[$i] = $args[$kcut_idx];
	$ins[$i] = $args[$ins_idx];
	$min_pairs[$i] = $args[$mnpr_idx];
	$max_maphits[$i] = $args[$mxhit_idx];
	$min_ctglen[$i] = $args[$mnlen_idx];
	$outdir[$i] = $args[$out_idx];	
}	

if ($dust){
	if ($paired){ 
		`tagdust -o $base.cln.fastq -a $base.aft.fastq -fdr 0.05 $adptrs $base\_p1.fastq $base\_p2.fastq`;
		`splitBC $base.cln.fastq --bcfile $conf_file --prefix ./ --suffix .fastq --mismatches 1 --boli`;
	} else {
		`tagdust -o $base.cln.fastq -a $base.aft.fastq -fdr 0.05 $adptrs $base.fastq`;
		foreach my $org (@base) {
			`repair `
		}
	}
}else {
	if ()
	`splitBC $base.cln.fastq --bcfile $conf_file --prefix ./ --suffix .fastq --mismatches 1 --bol`;
}
