#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
sub printNwide($$);
sub usage();
our ($opt_h, $opt_q, $opt_n, $opt_x);
my $optstr = "hqn:x:";
getopts($optstr);

usage() if ($opt_h || scalar(@ARGV)==0);

my $src = 1;
my $min = 1;
my $max = 9**9**9;
$min = $opt_n if ($opt_n);
$max = $opt_x if ($opt_x);
$src = 0 if ($opt_q);

my %orf;
my $blast_file = shift;
open (IN,"<",$blast_file);
while (<IN>){
	chomp;
	my @line = split /\t/;
	$orf{$line[$src]} = 1;
#	print STDERR "keeping ".$line[$src]."\n";
}
print STDERR "keeping ".scalar(keys %orf)." orfs\n";
while (my $line = <>) {
	if ($line =~ m/^>/){
LINE:{	chomp $line;
		my $name = (split(' ',substr($line,1)))[0];
		if (defined ($orf{$name})){
			my $hdr = $line;
			$line = <>;
			chomp $line;
			my $seq;
			while ($line = <>){
				chomp $line;
				if ($line !~ m/^>/) {
					$seq .= $line;
				} else {
					if (length($seq) < $max && length($seq) > $min) {
						print $hdr."\n";
						printNwide(80,$seq);
					} else {
						print STDERR "$name is not within length limit\n";
					}
					redo LINE;
				}
			} 
		}
	}
	} 
}

sub printNwide($$) {
	my $width = shift;
	my $seq = shift;
	while (length($seq) > $width){
		print substr($seq,0,$width)."\n";
		$seq = substr($seq,$width);
	}
	print $seq."\n";
}

sub usage() {
	print "Usage: filter_seqs_by_blast.pl [options] <m8_blast_outpu> <fasta_file(s)>\n";
	print "    where options are:\n";
	print "        -h        : print this message\n";
	print "        -q        : keep query sequences rather than subjects\n";
	print "                    i.e. look at column 1 rather than 2\n";
	print "        -n <int>  : minimum sequence length\n";
	print "        -x <int>  : maximum sequence length\n";
	print "    Note: if <fasta_file(s)> are missing, read from stdin\n";
	exit;
}
