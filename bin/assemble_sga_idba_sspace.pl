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
	# get the read length
	open(R1FQ, $r1fq);
	my $maxrdlen = 0;
	my $line = <R1FQ>;
	chomp $line;
	$maxrdlen = length($line);
	close R1FQ;
	`sga-static preprocess -q 10 -f 20 -m 30 --phred64 $r1fq $r2fq $rupfq > $outbase.pp.fastq`;
	`sga-static index -d 4000000 -t 4  $outbase.pp.fastq`;
	`sga-static correct -k 31 -i 10 -t 4  $outbase.pp.fastq`;
	`fastq2fasta.pl $outbase.pp.ec.fa > $outbase.clean.fa`;
	`rm $outbase.pp.*`;
	`idba -r $outbase.clean.fa -o $outbase --mink 33 --maxk $maxrdlen`;
	`rm $outbase.clean.fa`;
	# estimate the library insert size with bwa
	# just use a subsample of 40k reads
	`bwa index -a bwtsw $outbase-contig.fa`;
	`head -n 20000 $r1fq > $r1fq.sub`;
	`tail -n 20000 $r1fq >> $r1fq.sub`;
	`head -n 20000 $r2fq > $r2fq.sub`;
	`tail -n 20000 $r2fq >> $r2fq.sub`;
	`bwa aln $outbase-contig.fa $r1fq.sub > $r1fq.sub.sai`;
	`bwa aln $outbase-contig.fa $r2fq.sub > $r2fq.sub.sai`;
	my $ins_mean = 0;
	my $ins_error = 0;
	# bwa will print the estimated insert size, let's capture it then kill the job
	open(SAMPE, "bwa sampe -P -f $outbase-pe.sam $outbase-contig.fa $r1fq.sub.sai $r2fq.sub.sai $r1fq.sub $r2fq.sub |");
	while( my $line = <SAMPE> ){
		if($line =~ /inferred external isize from \d+ pairs: (\d+) \+\/\- (\d+)/){
			$ins_mean = $1;
			my $ins_sd = $2;
			$ins_error = $ins_error*4 / $ins_mean;
			close SAMPE;
			last;
		}
	}
	`rm $r1fq.sub* $r2fq.sub*`;

	open( LIBRARY, ">library.txt" );
	print LIBRARY "lib1 $r1fq $r2fq $ins_mean $ins_error 0";
	close LIBRARY;
	`/home/koadman/software/SSPACE-1.1_linux-x86_64/SSPACE_v1-1.pl -k 10 -a 0.2 -o 10 -l library.txt -s $outbase-contig.fa -b $outbase.sspace`;
}


