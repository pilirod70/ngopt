#!/usr/bin/env perl
#
# (c) 2011 Aaron Darling
# Licensed under the GPLv3
# 
# SGE commands for use on UCD cluster
#$ -S /usr/bin/perl
#$ -q eisen.q
#$ -q all.q
#$ -V

use strict;
use warnings;
use File::Basename;

# variables that may need to be changed:
# path to node-local temporary storage
my $tmpdir = "/state/partition1/koadman/tnpclean.".$ENV{"JOB_ID"};
# path to the e. coli reference db
my $ecoli_db = "/home/koadman/data/ecoli_bl21.fasta";
# minimum number of matching sites to consider the read as e. coli.
my $coli_match_threshold = 40;

die("Usage: colifilter.pl <read1file> <read1 suffix> <output dir>") unless @ARGV==3;

my $read1file = $ARGV[0];
my $read1suffix = $ARGV[1];
my $outdir = $ARGV[2];

`mkdir -p $tmpdir`;
chdir($tmpdir);

# stage the data
my $sourcedir = dirname($read1file);
get_copy_lock($sourcedir);
`cp $ecoli_db .`;
my $ecoli = basename($ecoli_db);
my $bname=basename($read1file, $read1suffix);
my $read2file = $read1file;
my $read2suffix = $read1suffix;
$read2suffix =~ s/1/2/g;
$read2file =~ s/$read1suffix/$read2suffix/g;
`cp $read1file $read2file .`;
release_copy_lock($sourcedir);

print "Processing $bname";
`bwa index $ecoli`;
`bwa aln $ecoli $bname$read1suffix > $bname.r1.sai`;
`bwa aln $ecoli $bname$read2suffix > $bname.r2.sai`;

bwasamse("$bname.r1.sai", "$bname$read1suffix", "$bname.r1.pure.fastq");
bwasamse("$bname.r2.sai", "$bname$read2suffix", "$bname.r2.pure.fastq");

`repair $bname $bname.r1.pure.fastq $bname.r2.pure.fastq`;
`gzip $bname*p?`;
my $mover = "mv $bname"."_p1.gz $outdir";
`$mover`;
$mover = "mv $bname"."_p2.gz $outdir";
`$mover`;
#`rm -rf $tmpdir`;

sub bwasamse
{
	my $sai = shift;
	my $reads = shift;
	my $out = shift;	
	open(BWA, "bwa samse $ecoli $sai $reads |");
	open(OUT, ">$out");
	while( my $line = <BWA> ){
		next if $line =~ /^@/;
		chomp $line;
		my @samline = split(/\t/, $line);
		my $rname = $samline[0];		
		my $cigar = $samline[5];
		my $matches = 0;
		my $indels = -1;
		while($cigar =~ s/(\d+)M//g){
			$matches += $1;
			$indels++;
		}
		next unless $matches >= $coli_match_threshold;
		$rname =~ s/\#0\/3$/\#0\/2/g;
		$rname =~ s/\#0$/\#0\/1/g;
		print OUT "@"."$rname\n";
		print OUT $samline[9]."\n";
		print OUT "+$rname\n";
		print OUT $samline[10]."\n";
	}
}

sub get_copy_lock
{
	my $dir = shift;
        my $notlocked = 1;
        while( $notlocked == 1 ){
                sleep(5);
                `mktemp $dir/copylock`;
                $notlocked = $?;
        }
}

sub release_copy_lock
{
	my $dir = shift;
	`rm $dir/copylock`;
}
