#!/usr/bin/env perl
use strict;

die "Usage: assemble_sga_idba_sspace.pl <read 1 fastq> <read 2 fastq> <unpaired fastq> <output base>\n" if(@ARGV<4);
idba_assemble($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);
#sga_assemble($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

sub sga_assemble {
	my $r1fq = shift;
	my $r2fq = shift;
	my $rupfq = shift;
	my $outbase = shift;
	`sga-static preprocess -p 1 -q 10 -f 20 -m 30 --phred64 $r1fq $r2fq > $outbase.pp.fastq`;
	`sga-static index -d 4000000 -t 4  $outbase.pp.fastq`;
	`sga-static correct -k 31 -i 10 -t 4  $outbase.pp.fastq`;
	`sga-static index -d 2000000 -t 4 $outbase.pp.ec.fa`;
	`sga-static qc -x 2 -t 4 $outbase.pp.ec.fa`;
	`sga-static rmdup -t 4 $outbase.pp.ec.qcpass.fa`;
	`sga-static overlap -m 30 -t 4 $outbase.pp.ec.qcpass.rmdup.fa`;
	`sga-static assemble -x 10 -b 5 -r 20 $outbase.pp.ec.qcpass.rmdup.asqg.gz`;
}

sub idba_assemble {
	my $r1fq = shift;
	my $r2fq = shift;
	my $rupfq = shift;
	my $outbase = shift;
	`sga-static preprocess -q 10 -f 20 -m 30 --phred64 $r1fq $r2fq $rupfq > $outbase.pp.fastq`;
	`sga-static index -d 4000000 -t 4  $outbase.pp.fastq`;
	`sga-static correct -k 31 -i 10 -t 4  $outbase.pp.fastq`;
# this seems to cut out way too many reads
#	`sga-static index -d 6000000 -t 4 $outbase.pp.ec.fa`;
#	`sga-static qc -x 2 -t 4 $outbase.pp.ec.fa`;
	`fastq2fasta.pl $outbase.pp.ec.fa > $outbase.clean.fa`;
	`rm $outbase.pp.*`;
	`idba -r $outbase.clean.fa -o $outbase --mink 33 --maxk 84`;
	`rm $outbase.clean.fa`;
	open( LIBRARY, ">library.txt" );
	print LIBRARY "lib1 $r1fq $r2fq 300 0.75 0";
	close LIBRARY;
	`/home/koadman/software/SSPACE-1.1_linux-x86_64/SSPACE_v1-1.pl -k 10 -a 0.2 -o 10 -l library.txt -s $outbase-contig.fa -b $outbase.sspace`;
}


