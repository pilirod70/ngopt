#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::Perl;
use Getopt::Std;
sub usage();

our ($opt_h,$opt_p,$opt_l,$opt_m,$opt_g,$opt_e,$opt_b);
my $optstr = "hp:l:m:g:e:b:";

my $min_pid = 0;
my $min_len = 1;
my $max_eval = 9**9**9;
my $min_bit = 0;
my $max_dif = 9**9**9;
my $max_gap = 9**9**9;
my $dir="/share/eisen-d6/halophile_illumina_data/scaffold/";
getopts($optstr);
usage() if ($opt_h || scalar(@ARGV)==0);
$min_pid = $opt_p if ($opt_p);
$min_len = $opt_l if ($opt_l);
$max_eval = $opt_e if ($opt_e);
$min_bit = $opt_b if ($opt_b);
$max_dif = $opt_m if ($opt_m);
$max_gap = $opt_g if ($opt_g);

while (<>) {
	chomp;
	my @line = split /\t/;
	print join("\t",@line)."\n" if ($line[4] < $max_dif && $line[5] < $max_gap && 
	                                $line[10] < $max_eval && $line[3] > $min_len && 
	                                $line[2] > $min_pid && $line[11] > $min_bit);
}


sub usage(){
	print "Usage: filter_blast.pl [options] <m8_blast_output>\n";
	print "    where options are:\n";
	print "        -h           : print this help message\n";
	print "        -p <float>   : minimum percent id\n";
	print "        -l <int>     : minimum alignment length\n";
	print "        -m <int>     : maximum mismatches\n";
	print "        -g <int>     : maximum gap-openings\n";
	print "        -e <float>   : maximum e-value\n";
	print "        -b <int>     : minimum bit-score\n";
	print "    Note: if <m8_blast_output> missing, read from stdin\n";
	exit;
}
