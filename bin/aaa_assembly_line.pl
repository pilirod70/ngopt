#!/usr/bin/env perl
use strict;

die "Usage: andrews_assembly_line.pl <library file> <output base>\n" if(@ARGV!=2);
my $libfile = $ARGV[0];
my $outbase = $ARGV[1];

my %libs = ();
sga_clean($libfile, $outbase, \%libs);
my $maxrdlen = fastq_to_fasta("$outbase.pp.ec.fa", "$outbase.clean.fa");
idba_assemble($outbase, $maxrdlen);
scaffold_sspace($libfile, $outbase, \%libs);

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

sub sga_clean {
	my $libfile = shift;
	my $outbase = shift;
	my $libs = shift;
	open(LIB,"<",$libfile);
	my $lib_count = 0;
	my %hash = ();
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
	system("sga-static preprocess -q 10 -f 20 -m 30 --phred64 $files > $outbase.pp.fastq");
	die "Error preprocessing reads with SGA\n" if( $? != 0 );
	my $sga_ind = "";
	my $sga_ind_kb = 4000000;
	do{
		$sga_ind = `sga-static index -d $sga_ind_kb -t 4  $outbase.pp.fastq > index.out 2>&1`;
		$sga_ind_kb = int($sga_ind_kb/2);
	}while($sga_ind =~ /bad_alloc/ || $? != 0);
	system("rm -f core*");
	die "Error indexing reads with SGA\n" if( $? != 0 );
	system("sga-static correct -k 31 -i 10 -t 4  $outbase.pp.fastq > correct.out");
	die "Error correcting reads with SGA\n" if( $? != 0 );
}


sub fastq_to_fasta {
	my $fastq = shift;
	my $fasta = shift;
	my $maxrdlen = 0;
	open(FQ,"<$outbase.pp.ec.fa");
	open(FA,">$outbase.clean.fa");
	while(my $hdr = <FQ>) {
		my $seq = <FQ>;
		my $qhdr = <FQ>;
		my $qseq = <FQ>;
		$hdr =~ s/^@/>/g;
		print FA $hdr.$seq;
		$maxrdlen = length($seq) if length($seq) > $maxrdlen;
	}
	close FA;
	return $maxrdlen - 2;	# -1 removes newline char
}

# expects a file called $outbase.clean.fa in the current working directory
sub idba_assemble {
	my $outbase = shift;
	my $maxrdlen = shift;
	my $idba_cmd = "idba -r $outbase.clean.fa -o $outbase --mink 29 --maxk $maxrdlen";
	print STDERR $idba_cmd."\n";
	`$idba_cmd > idba.out`;
	die "Error building contigs with IDBA\n" if ($? != 0);
	`gzip -f $outbase.clean.fa`;
	`rm $outbase.pp.* $outbase.kmer $outbase.graph`;
}


sub scaffold_sspace {
	my $libfile = shift;
	my $outbase = shift;
	# build library file
	my $libs = shift;
	my @library_file = ();
	my $lib_files;
	my $min_insert_size = -1;
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
			$min_insert_size = $ins_mean if($min_insert_size == -1 || $min_insert_size > $ins_mean);
			my $outtie = 0;
			$outtie = 1 if $ins_mean > 1500;	# assume that anything > 1500 used mate-pairing
			push (@library_file, "$lib $fq1 $fq2 $ins_mean $ins_err $outtie\n");
		}
		# if we have one file left, treat it as unpaired. 
		if (scalar(@$lib_files) == 1){
			my $up = shift @$lib_files;
			`cat $up >> $outbase.unpaired.fastq`;
		}
	}
	my $genome_size = get_genome_size("$outbase-contig.fa");
	print STDERR "Total contig length $genome_size\n";
	my $libraryI=1;
	open( LIBRARY, ">library_$libraryI.txt" );
	my $prev_ins = -1;
	my $prev_reads = "";
	# sort library.txt to so that we scaffold with smaller insert libraries first
	for my $line (sort {(split(' ',$a))[3] <=>  (split(' ',$b))[3]} @library_file) {
		# if we've hit a substantially different insert size, do a separate
		# round of scaffolding
		my $cur_ins = (split(' ',$line))[3];
		if($prev_ins > 0 && ($prev_ins * 2) < $cur_ins){
			close LIBRARY;
			my $exp_link = calc_explinks( $genome_size, $prev_ins, $prev_reads ); 
			print STDERR "Insert $prev_ins, expected links $exp_link\n";
			run_sspace($genome_size, $prev_ins, $exp_link, $libraryI);
			$libraryI++;
			open( LIBRARY, ">library_$libraryI.txt" );
		}
		$prev_reads = (split(' ',$line))[1];
		$prev_ins = $cur_ins;
		print LIBRARY $line;	
	}
	close LIBRARY;
	my $exp_link = calc_explinks( $genome_size, $prev_ins, $prev_reads ); 
	print STDERR "Insert $prev_ins, expected links $exp_link\n";
	run_sspace($genome_size, $prev_ins, $exp_link, $libraryI);
	`mv $outbase.lib$libraryI.sspace.final.scaffolds.fasta $outbase.sspace.final.scaffolds.fasta`;
}



sub fish_break_misasms {
	my $ctgs = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $sai = "sai"
	my $sam = "sam"
	my $broken_ctgs = "$ctgs.broken"
	`bwa index -a is $ctgs`
	`cat $fq1 $fq2 | bwa aln $ctgs - > $sai`
	`cat $fq1 $fq2 | bwa samse $ctgs $sai - > $sam`
	`get_fish_input.sh $sam . > get_fish_input.out 2> get_fish_input.err`	
	`fish -b blocks.txt > fish.out 2> fish.err`
	`break_misassemblies.pl blocks.txt contig_labels.txt $ctgs > $broken_ctgs 2> break_misasm.err` 
	return $broken_ctgs;
}

# calculate the expected number of read pairs to span a point in the assembly
sub calc_explinks {
	my $genome_size = shift;
	my $ins_len = shift;
	my $read_file = shift;
	my $read_count = 0;
	my $maxrdlen = -1;
	open( READFILE, $read_file );
	while( my $line = <READFILE> ){
		$read_count++;
		$maxrdlen = length($line) if $read_count % 4 == 2 && $maxrdlen < length($line);
	}
	$read_count /= 4;
	my $cov = $maxrdlen * $read_count / $genome_size;
	print STDERR "Lib $read_file coverage $cov\n";
	my $exp_link = $cov * $ins_len / $maxrdlen;
	return $exp_link;
}

sub run_sspace {
	my $genome_size = shift;
	my $insert_size = shift;
	my $exp_links = shift;
	my $libraryI = shift;
	my $sspace_m = int(log2($genome_size)+3.99);
	my $sspace_n = int(log2($insert_size)*1.25+.99);
	my $sspace_k = int(log($exp_links)/log(1.4)-9.5);
	# require at least 2 links to preclude chimeras
	$sspace_k = $sspace_k < 2 ? 2 : $sspace_k;	
#	$sspace_k = 5;	# gives best results on mediterranei
	my $input_fa = "$outbase-contig.fa";
	$input_fa = "$outbase.lib".($libraryI-1).".sspace.final.scaffolds.fasta" if $libraryI>1;
	my $sspace_cmd = "SSPACE_v1-1.pl -m $sspace_m -n $sspace_n -k $sspace_k -a 0.2 -o 1 -l library_$libraryI.txt -s $input_fa -b $outbase.lib$libraryI.sspace";
	if (-f "$outbase.unpaired.fastq") {
		print STDERR "Running SSPACE with unpaired reads\n";
		$sspace_cmd .= " -u $outbase.unpaired.fastq";
	}
	print STDERR "$sspace_cmd\n";
	`$sspace_cmd > sspace_lib$libraryI.out`;
	`rm -rf bowtieoutput/ reads/`
}

sub log2 {
	my $n = shift;
	return (log($n)/ log(2));
}

sub get_genome_size {
	my $fasta = shift;
	open( FASTA, $fasta );
	my $len = 0;
	while( my $line = <FASTA> ){
		next if $line =~ /^>/;
		chomp $line;
		$len += length($line);
	}
	return $len;
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

