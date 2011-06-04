#!/usr/bin/env perl
use strict;

die "Usage: andrews_assembly_line.pl <library file> <output base>\n" if(@ARGV!=2);
idba_assemble($ARGV[0], $ARGV[1]);
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
	my $libfile = shift;
	my $outbase = shift;
	open(LIB,"<",$libfile);
	my $lib_count = 0;
	my %hash = ();
	my %libs = ();
	my $maxrdlen = 0;
	my $files = "";
	while(<LIB>){
		if ($_ =~ m/\[LIB\]/){
			if ($lib_count > 0) {
				for my $key (sort keys %hash){
					push(@{$libs{"lib$lib_count"}},$hash{$key});
				}
			} 
			$lib_count++;
		} elsif ($_ =~ m/(p[1,2])=([\w\/\-\.]+)/) { 
			$hash{$1} = $2;
			$files .= "$2 ";
			my $len = get_rdlen($2);
			$maxrdlen = $len if ($len > $maxrdlen);
		} elsif ($_ =~ m/up=([\w\/\-\.]+)/) { 
			$hash{"up"} = $1;
			$files .= "$1 ";
			my $len = get_rdlen($1);
			$maxrdlen = $len if ($len > $maxrdlen);
		} elsif ($_ =~ m/ins=([\w\/\-\.]+)/) { 
			$hash{"ins"} = $1;
		}
	} 
	for my $key (sort keys %hash){
		push(@{$libs{"lib$lib_count"}},$hash{$key});
	}
	`sga-static preprocess -q 10 -f 20 -m 30 --phred64 $files > $outbase.pp.fastq`;
	`sga-static index -d 4000000 -t 4  $outbase.pp.fastq > index.out`;
	`sga-static correct -k 31 -i 10 -t 4  $outbase.pp.fastq > correct.out`;
	open(FQ,"<$outbase.pp.ec.fa");
	open(FA,">$outbase.clean.fa");
	while(my $hdr = <FQ>) {
		my $seq = <FQ>;
		my $qhdr = <FQ>;
		my $qseq = <FQ>;
		$hdr =~ s/^@/>/g;
		print FA $hdr.$seq;
	}
	close FA;
	my $idba_cmd = "idba -r $outbase.clean.fa -o $outbase --mink 33 --maxk $maxrdlen";
	print STDERR $idba_cmd."\n";
	`$idba_cmd > idba.out`;
	`gzip -f $outbase.clean.fa`;
	`rm $outbase.pp.* $outbase.kmer $outbase.graph`;

	# build library file
	my @library_file = ();
	my $lib_files;
	`rm $outbase.unpaired.fastq` if (-f "$outbase.unpaired.fastq");
	for my $lib (keys %libs) {
	#	print STDERR "$lib\t".join(" ",@{$libs{$lib}})."\n";
		$lib_files = \@{$libs{$lib}};
		# if we have at least two files in this library, assume the first two are paired
		if (scalar(@$lib_files) >= 2) {
			my $ins = shift @$lib_files;
			my $fq1 = shift @$lib_files;
			my $fq2 = shift @$lib_files;
			my ($ins_mean, $ins_err) = get_insert($fq1,$fq2,$outbase,$lib);
			if ($ins_mean == -1) {
				$ins_mean = $ins;
				$ins_err = 0.7;
			}
			push (@library_file, "$lib $fq1 $fq2 $ins_mean $ins_err 0\n");
		}
		# if we have one file left, treat it as unpaired. 
		if (scalar(@$lib_files) == 1){
			my $up = shift @$lib_files;
			`cat $up >> $outbase.unpaired.fastq`;
		}
	}
	open( LIBRARY, ">library.txt" );
	# sort library.txt to so that we scaffold with smaller insert libraries first
	for my $line (sort {(split(' ',$a))[3] <=>  (split(' ',$b))[3]} @library_file) {
		print LIBRARY $line;	
	}
	close LIBRARY;
	my $sspace_cmd = "SSPACE_v1-1.pl -x 1 -k 10 -a 0.2 -o 10 -l library.txt -s $outbase-contig.fa -b $outbase.sspace";
	if (-f "$outbase.unpaired.fastq") {
		print STDERR "Running SSPACE with unpaired reads\n";
		$sspace_cmd .= " -u $outbase.unpaired.fastq";
	}
	print STDERR "$sspace_cmd\n";
	`$sspace_cmd > sspace.out`;
	`rm -rf bowtieoutput/ reads/`
}


sub get_insert($$$$) {
	my $r1fq = shift;
	my $r2fq = shift;
	my $outbase = shift;
	my $lib = shift;
	# estimate the library insert size with bwa
	# just use a subsample of 40k reads
	`bwa index -a is $outbase-contig.fa`;
	`head -n 40000 $r1fq > $r1fq.sub`;
	`tail -n 40000 $r1fq >> $r1fq.sub`;
	`head -n 40000 $r2fq > $r2fq.sub`;
	`tail -n 40000 $r2fq >> $r2fq.sub`;
	`bwa aln $outbase-contig.fa $r1fq.sub > $r1fq.sub.sai`;
	`bwa aln $outbase-contig.fa $r2fq.sub > $r2fq.sub.sai`;
	# bwa will print the estimated insert size, let's capture it then kill the job
	open(SAMPE, "bwa sampe -P -f $outbase-pe.sam $outbase-contig.fa $r1fq.sub.sai $r2fq.sub.sai $r1fq.sub $r2fq.sub 2>&1 | tee $outbase.$lib.sampe.out |");
	my $ins_mean = 0;
	my $ins_error = 0;
	my $min;
	my $max;
	my $ins_sd;
	my $ins_n;
	while( my $line = <SAMPE> ){
		if($line =~ m/^\[infer_isize\] inferred external isize from (\d+) pairs: ([\d\.]+) \+\/\- ([\d\.]+)$/){
			$ins_n = $1;
			$ins_mean = $2;
			$ins_sd = $3;
			close SAMPE;
			last;
		} elsif ($line =~ m/^\[infer_isize\] low and high boundaries: (\d+) and (\d+) for estimating avg and std$/){
			$min = $1;
			$max = $2;
		}
	}
#	$min = ($ins_mean - $min)
#	$max = ($max - $ins_mean);
#	$ins_error = (($min,$max)[$min < $max])/$ins_mean;
	`rm $r1fq.sub* $r2fq.sub*`;
	if ($ins_n > 20000) {		
		$ins_error = $ins_sd*6 / $ins_mean;
		$ins_mean = sprintf("%.0f",$ins_mean);
		$ins_error = sprintf("%.3f",$ins_error);
		$ins_error =~ s/0+$//g;
	} else {
		$ins_mean = -1;
		$ins_error = 0;
	}
	return ($ins_mean, $ins_error);
}

sub get_rdlen($$){
	open(FILE, shift);
	my $line = <FILE>;
	$line = <FILE>;
	chomp $line;
	close FILE;
	return length($line);
}

