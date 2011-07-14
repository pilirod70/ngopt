#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

die "Usage: ".basename($0)." [--begin=1-5] <library file> <output base>\n" if(@ARGV!=2 && @ARGV!=3 || (@ARGV==3 && $ARGV[0] !~ m/--begin=([1-5])/));
my $start = 1;
if (@ARGV==3) {
	$start = $1;
	shift;
}

my $libfile = $ARGV[0];
my $outbase = $ARGV[1];


my $DIR = dirname(abs_path($0));

my %LIBS = read_lib_file($libfile);
print STDERR "Found ".scalar(keys %LIBS)." libraries\n";
my $maxrdlen = -1;
my $scafs;

print "$libfile $start\n";

if ($start <= 1) {
	print "CHECKPOINT: Cleaning reads with SGA\n";
	print STDERR "CHECKPOINT: Cleaning reads with SGA\n";
	sga_clean($outbase, \%LIBS); 
} 
if ($start <= 2) {
	print "CHECKPOINT: Building contigs with IDBA\n";
	print STDERR "CHECKPOINT: Building contigs with IDBA\n";
	$maxrdlen = fastq_to_fasta("$outbase.pp.ec.fa", "$outbase.clean.fa");
	idba_assemble($outbase, $maxrdlen); 
} 
if ($start <= 3) {
	print "CHECKPOINT: Scaffolding contigs with SSPACE\n";
	print STDERR "CHECKPOINT: Scaffolding contigs with SSPACE\n";
	$scafs = scaffold_sspace($libfile, $outbase, \%LIBS, "$outbase-contig.fa");
	`mv $scafs $outbase.sspace.scaffolds.fasta`; 
} 
if ($start <= 4) {
	print "CHECKPOINT: Detecting and breaking misassemblies with FISH\n";
	print STDERR "CHECKPOINT: Detecting and breaking misassemblies with FISH\n";
	$scafs = "$outbase.sspace.scaffolds.fasta";
	$scafs = break_all_misasms($scafs,\%LIBS,"$outbase.fish"); 
	`mv $scafs $outbase.fish.broken.fasta`; 
} 
if ($start <= 5) {
	print "CHECKPOINT: Scaffolding broken contigs with SSPACE\n";
	print STDERR "CHECKPOINT: Scaffolding broken contigs with SSPACE\n";
	$scafs = "$outbase.fish.broken.fasta"; 
	$scafs = scaffold_sspace($libfile,"$outbase.rescaf",\%LIBS,$scafs);
	`mv $scafs $outbase.a5.final.fasta`;
	print "Final assembly in $outbase.a5scafs.fasta\n"; 
} 

sub sga_assemble {
	my $r1fq = shift;
	my $r2fq = shift;
	my $rupfq = shift;
	my $outbase = shift;
	`$DIR/sga preprocess -p 1 -q 10 -f 20 -m 30 --phred64 $r1fq $r2fq > $outbase.pp.fastq`;
	`$DIR/sga index -d 4000000 -t 4  $outbase.pp.fastq`;
	`$DIR/sga correct -k 31 -i 10 -t 4  $outbase.pp.fastq`;
	`$DIR/sga index -d 2000000 -t 4 $outbase.pp.ec.fa`;
	`$DIR/sga qc -x 2 -t 4 $outbase.pp.ec.fa`;
	`$DIR/sga rmdup -t 4 $outbase.pp.ec.qcpass.fa`;
	`$DIR/sga overlap -m 30 -t 4 $outbase.pp.ec.qcpass.rmdup.fa`;
	`$DIR/sga assemble -x 10 -b 5 -r 20 $outbase.pp.ec.qcpass.rmdup.asqg.gz`;
}

sub sga_clean {
	my $outbase = shift;
	my $libs = shift;
	# figure out which files we need to pass to SGA
	my $files = "";
	for my $lib (keys %$libs) {
		my @lib_files = @{$libs->{$lib}};
		if (scalar(@lib_files) >= 2) {
			shift @lib_files;
			my $fq1 = shift @lib_files;
			my $fq2 = shift @lib_files;
			$files .= "$fq1 $fq2 ";
		}
		if (scalar(@lib_files) == 1){
			my $up = shift @lib_files;
			$files .= "$up ";
		}
	}
	print STDERR "sga preprocess -q 10 -f 20 -m 30 --phred64 $files > $outbase.pp.fastq\n";
	system("$DIR/sga preprocess -q 10 -f 20 -m 30 --phred64 $files > $outbase.pp.fastq");
	die "Error preprocessing reads with SGA\n" if( $? != 0 );
	my $sga_ind = "";
	my $sga_ind_kb = 4000000;
	do{
		$sga_ind = `$DIR/sga index -d $sga_ind_kb -t 4  $outbase.pp.fastq > index.out 2> index.err`;
		$sga_ind = read_file('index.err');
		$sga_ind_kb = int($sga_ind_kb/2);
	}while(($sga_ind =~ /bad_alloc/ || $? != 0) && $sga_ind_kb > 0);
	system("rm -f core*") if (-f "core*");
	die "Error indexing reads with SGA\n" if( $? != 0 );
	system("$DIR/sga correct -k 31 -i 10 -t 4  $outbase.pp.fastq > correct.out");
	die "Error correcting reads with SGA\n" if( $? != 0 );
}

sub read_lib_file {
	my $libfile = shift;
	my %libs = ();
	open(LIB,"<",$libfile);
	my $lib_count = 0;
	my %hash = ();
	while(<LIB>){
		chomp;
		if ($_ =~ m/\[LIB\]/){
			if ($lib_count > 0) {
				for my $key (sort keys %hash){
					push(@{$libs{"lib$lib_count"}},$hash{$key});
				}
			} 
			$lib_count++;
		} elsif ($_ =~ m/shuf=([\w\/\-\.]+)/) { 
			my ($fq1, $fq2) = split_shuf($1,"$outbase.lib$lib_count");
			$hash{"p1"} = $fq1;
			$hash{"p2"} = $fq2;
		} elsif ($_ =~ m/(p[1,2])=([\w\/\-\.]+)/) { 
			$hash{$1} = $2;
		} elsif ($_ =~ m/up=([\w\/\-\.]+)/) { 
			$hash{"up"} = $1;
		} elsif ($_ =~ m/ins=([\w\/\-\.]+)/) { 
			$hash{"ins"} = $1;
		} else {
			die "Unrecognizable line in library file: >$_<\n";
		}
	} 
	for my $key (sort keys %hash){
		push(@{$libs{"lib$lib_count"}},$hash{$key});
	}
	return %libs;
}

sub read_file {
	my $file = shift;
	local $/=undef;
	open(FILE,"<",$file);
	my $str = <FILE>;
	return $str;
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
	# why are we removing 2 here?
	return $maxrdlen - 2;	# -1 removes newline char
}

sub split_shuf {
	my $shuf = shift;
	my $outbase = shift;
	my $fq1 = "$outbase\_p1.fastq";
	my $fq2 = "$outbase\_p2.fastq";
	open(FQ1,">$fq1");
	open(FQ2,">$fq2");
	open(IN,"<$shuf");
	while(<IN>){
		print FQ1 $_;
		print FQ1 <IN>;
		print FQ1 <IN>;
		print FQ1 <IN>;
		print FQ2 <IN>;
		print FQ2 <IN>;
		print FQ2 <IN>;
		print FQ2 <IN>;
	}
	close FQ1;
	close FQ2;
	return ($fq1, $fq2);
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
	my $curr_ctgs = shift;
	my @library_file = ();
	my @lib_files;
	my $min_insert_size = -1;
	# remove the unpaired file if present 
	`rm $outbase.unpaired.fastq` if (-f "$outbase.unpaired.fastq");
	for my $lib (keys %$libs) {
	#	print STDERR "$lib\t".join(" ",@{$libs{$lib}})."\n";
		# need to make a copy of this little array, so we can preserve the original for break_all_misasms
		@lib_files = @{$libs->{$lib}};
		#$lib_files = \@{$libs{$lib}};
		# if we have at least two files in this library, assume the first two are paired
		if (scalar(@lib_files) >= 2) {
			my $ins = shift @lib_files;
			my $fq1 = shift @lib_files;
			my $fq2 = shift @lib_files;
			my ($ins_mean, $ins_err) = get_insert($fq1,$fq2,$outbase,$lib,$curr_ctgs);
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
		if (scalar(@lib_files) == 1){
			my $up = shift @lib_files;
			`cat $up >> $outbase.unpaired.fastq`;
		}
	}
	my $genome_size = get_genome_size($curr_ctgs);
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
			$curr_ctgs = run_sspace($genome_size, $prev_ins, $exp_link, $libraryI, $curr_ctgs);
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
	my $fin_scafs = run_sspace($genome_size, $prev_ins, $exp_link, $libraryI,$curr_ctgs);
	`mv $fin_scafs $outbase.sspace.final.scaffolds.fasta`;
	return "$outbase.sspace.final.scaffolds.fasta";
}

sub break_all_misasms {
	my $ctgs = shift;
	my $libs = shift;
	my $outbase = shift;
	my @lib_files;
	my %tmp = %$libs;
	for my $lib (keys %$libs){
		$tmp{$lib} = @{$libs->{$lib}} if (scalar(@{$libs->{$lib}}) >= 2);
	}
	# sort libraries by insert so we break with larger inserts first.
	for my $lib (sort {$tmp{$b}[0] <=> $tmp{$a}[0]} keys %tmp){
		@lib_files = @{$libs->{$lib}};
		next unless (scalar(@lib_files) >= 2);
		my $ins = shift @lib_files;
		my $fq1 = shift @lib_files;
		my $fq2 = shift @lib_files;
		$ctgs = fish_break_misasms($ctgs,$fq1,$fq2,"$outbase.lib$lib");
	}	
	return $ctgs;
}

sub fish_break_misasms {
	my $ctgs = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $outbase = shift;
	my $sai = "$outbase.sai";
	my $sam = "$outbase.sam";
	`$DIR/bwa index -a is $ctgs > $outbase.index.out`;
	`cat $fq1 $fq2 | $DIR/bwa aln $ctgs - > $sai`;
	`cat $fq1 $fq2 | $DIR/bwa samse $ctgs $sai - > $sam`;
	`java -jar $DIR/GetFishInput.jar $sam $outbase > $outbase.fie.out`;
	`$DIR/fish -f $outbase.control.txt -b $outbase.blocks.txt > $outbase.fish.out`;
	die "Error getting blocks with FISH for $outbase\n" if ($? != 0);
	`$DIR/break_misassemblies.pl $outbase.blocks.txt contig_labels.txt $ctgs > $outbase.broken.fasta 2> $outbase.break.out`;
	die "Error getting breaking contigs after running FISH\n" if ($? != 0);
	#`rm $outbase.blocks.txt contig_labels.txt fish.* get_fish_input.* break_misasm.err $sam $sai`;
	return "$outbase.broken.fasta";
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
	#my $input_fa = "$outbase-contig.fa";
	my $input_fa = shift;
	$input_fa = "$outbase.lib".($libraryI-1).".sspace.final.scaffolds.fasta" if $libraryI>1;
	my $sspace_cmd = "$DIR/SSPACE/SSPACE -m $sspace_m -n $sspace_n -k $sspace_k -a 0.2 -o 1 -l library_$libraryI.txt -s $input_fa -b $outbase.lib$libraryI.sspace";
	if (-f "$outbase.unpaired.fastq") {
		print STDERR "Running SSPACE with unpaired reads\n";
		$sspace_cmd .= " -u $outbase.unpaired.fastq";
	}
	print STDERR "$sspace_cmd\n";
	`$sspace_cmd > sspace_lib$libraryI.out`;
	`rm -rf bowtieoutput/ reads/`;
	return "$outbase.lib$libraryI.sspace.final.scaffolds.fasta";
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

sub get_insert($$$$$) {
	my $r1fq = shift;
	my $r2fq = shift;
	my $outbase = shift;
	my $lib = shift;
	my $ctgs = shift;

	my $estimate_pair_count = 20000;
	my $require_fraction = 0.25;
	my $fq_linecount = $estimate_pair_count*2;
	# estimate the library insert size with bwa
	# just use a subsample of 40k reads
	`$DIR/bwa index -a is $ctgs`;
	`head -n $fq_linecount $r1fq > $r1fq.sub`;
	`tail -n $fq_linecount $r1fq >> $r1fq.sub`;
	`head -n $fq_linecount $r2fq > $r2fq.sub`;
	`tail -n $fq_linecount $r2fq >> $r2fq.sub`;
	`$DIR/bwa aln $ctgs $r1fq.sub > $r1fq.sub.sai`;
	`$DIR/bwa aln $ctgs $r2fq.sub > $r2fq.sub.sai`;
	# bwa will print the estimated insert size, let's capture it then kill the job
	open(SAMPE, "$DIR/bwa sampe -P -f $outbase-pe.sam $ctgs $r1fq.sub.sai $r2fq.sub.sai $r1fq.sub $r2fq.sub 2>&1 | tee $outbase.$lib.sampe.out |");
	my $ins_mean = 0;
	my $ins_error = 0;
	my $min;
	my $max;
	my $ins_sd;
	my $ins_n;
	my $found = 0;
	while( my $line = <SAMPE> ){
		if(!$found  && $line =~ m/^\[infer_isize\] inferred external isize from (\d+) pairs: ([\d\.]+) \+\/\- ([\d\.]+)$/){
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
	if ($ins_n > $require_fraction * $estimate_pair_count) {		
		$ins_error = $ins_sd*6 / $ins_mean;
		$ins_mean = sprintf("%.0f",$ins_mean);
		$ins_error = sprintf("%.3f",$ins_error);
		$ins_error =~ s/0+$//g;
	} else {
		print STDERR "Discarding estimate. Not enough data points: $ins_n\n";
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

