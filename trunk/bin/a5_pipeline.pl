#!/usr/bin/env perl
#
# a5 pipeline
# (c) 2011 Andrew Tritt and Aaron Darling
# This is free software licensed under the GPL
#
# Usage: a5_pipeline.pl <library file> <output base>
#
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
my $OUTBASE = $ARGV[1];

my $DIR = dirname(abs_path($0));
my %RAW_LIBS = read_lib_file($libfile);
print STDERR "[a5] Found the following libraries:\n";
for my $lib (sort keys %RAW_LIBS) {
	my $notfirst = 0;
	print STDERR "     $lib: ";
	for my $att (sort keys %{$RAW_LIBS{$lib}}){
		print STDERR " $att=".$RAW_LIBS{$lib}{$att};
	}
	print STDERR "\n";
}

die "[a5] No libraries found in $libfile\n" unless(keys %RAW_LIBS);
print "[a5] Found ".scalar(keys %RAW_LIBS)." libraries\n";
my $maxrdlen = -1;
my $scafs;
my $ctgs;
my $reads;
my $WD="";

print "[a5] Starting pipeline at step $start\n";

if ($start <= 1) {
	print "[a5] Cleaning reads with SGA\n";
	print STDERR "[a5] Cleaning reads with SGA\n";
	$WD = "s1";
	mkdir($WD) if ! -d $WD;
	$reads = sga_clean($OUTBASE, \%RAW_LIBS);
	$reads = tagdust($OUTBASE, $reads);
	`mv $reads $OUTBASE.ec.fastq`;
} 
if ($start <= 2) {
	my $fq_reads = "$OUTBASE.ec.fastq";
	`gunzip $fq_reads.gz` if -f "$fq_reads.gz";
	$reads = $fq_reads;
	die "[a5_s2] Can't find error corrected reads $reads" unless -f $reads;
	print "[a5_s2] Building contigs from $reads with IDBA\n";
	print STDERR "[a5_s2] Building contigs from $reads with IDBA\n";
	$WD="s2";
	mkdir($WD) if ! -d $WD;
	($reads, $maxrdlen) = fastq_to_fasta($reads,"$WD/$OUTBASE.ec.fasta");
#	`gzip -f $fq_reads`;

	print STDERR "$reads exists\n" if -f $reads;
	print STDERR "$reads does not exist\n" if ! -f $reads;
	$ctgs = idba_assemble($OUTBASE, $reads, $maxrdlen); 
	`rm $reads`;
	`mv $ctgs $OUTBASE.contigs.fasta`;
	
} 
$WD="s3";
mkdir($WD) if ! -d $WD;
print STDERR "[a5] Preprocess libraries for scaffolding with SSPACE\n";
my %PAIR_LIBS = preprocess_libs(\%RAW_LIBS,"$OUTBASE.contigs.fasta");
if ($start <= 3) {
	$ctgs = "$OUTBASE.contigs.fasta";
	die "[a5_s3] Can't find starting contigs $ctgs.\n" unless -f $ctgs;
	print "[a5_s3] Scaffolding contigs from $ctgs with SSPACE\n";
	print STDERR "[a5_s3] Scaffolding contigs from $ctgs with SSPACE\n";
	$scafs = scaffold_sspace($libfile, $OUTBASE, \%PAIR_LIBS, $ctgs);
	`mv $scafs $OUTBASE.crude.scaffolds.fasta`; 
} 
if ($start <= 4) {
	$scafs = "$OUTBASE.crude.scaffolds.fasta";
	die "[a5_s4] Can't find starting crude scaffolds $scafs\n" unless -f $scafs;
	print "[a5_s4] Detecting and breaking misassemblies in $scafs with FISH\n";
	print STDERR "[a5_s4] Detecting and breaking misassemblies in scafs with FISH\n";
	$WD="s4";
	mkdir($WD) if ! -d $WD;
	$scafs = break_all_misasms($scafs,\%PAIR_LIBS,"$OUTBASE.fish"); 
	`mv $scafs $OUTBASE.broken.scaffolds.fasta`; 
} 
if ($start <= 5) {
	$scafs = "$OUTBASE.broken.scaffolds.fasta";
	if (-z $scafs) {
		print "[a5_s5] No misassemblies found.\n";
		`cp $OUTBASE.crude.scaffolds.fasta $OUTBASE.final.scaffolds.fasta`;
	} else {
		die "[a5_s5] Can't find starting broken scaffolds $scafs\n" unless -f $scafs;
		print "[a5_s5] Scaffolding broken contigs with SSPACE\n";
		print STDERR "[a5_s5] Scaffolding broken contigs with SSPACE\n";
		$WD="s5";
		mkdir($WD) if ! -d $WD;
		$scafs = scaffold_sspace($libfile,"$OUTBASE.rescaf",\%PAIR_LIBS,$scafs);
		`mv $scafs $OUTBASE.final.scaffolds.fasta`;
	}
	print "[a5] Final assembly in $OUTBASE.final.scaffolds.fasta\n"; 
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
	my $libsref = shift;
	my %libs = %$libsref;
	# figure out which files we need to pass to SGA
	my $files = "";
	for my $lib (keys %libs) {
		if (defined($libs{$lib}{"p1"})) {
			my $fq1 = $libs{$lib}{"p1"}; 
			my $fq2 = $libs{$lib}{"p2"}; 
			$files .= "$fq1 $fq2 ";
		}
		if (defined($libs{$lib}{"up"})){
			my $up = $libs{$lib}{"up"}; 
			$files .= "$up ";
		}
	}
	my $cmd = "sga preprocess -q 10 -f 20 -m 30 --phred64 $files > $WD/$outbase.pp.fastq";
	print STDERR "[a5] $cmd\n";
	system("$DIR/$cmd");
	die "[a5] Error preprocessing reads with SGA\n" if( $? != 0 );
	my $sga_ind = "";
	my $sga_ind_kb = 4000000;
	my $err_file = "$WD/index.err";
	do{
		$cmd = "sga index -d $sga_ind_kb -t 4 $WD/$outbase.pp.fastq > $WD/index.out 2> $err_file";
		print STDERR "[a5] $cmd\n";
		$sga_ind = `$DIR/$cmd`;
		$sga_ind = read_file($err_file);
		$sga_ind_kb = int($sga_ind_kb/2);
	}while(($sga_ind =~ /bad_alloc/ || $? != 0) && $sga_ind_kb > 0);
	my $ec_file = "$WD/$outbase.pp.ec.fa";
	system("rm -f core*") if (-f "core*");
	die "[a5] Error indexing reads with SGA\n" if( $? != 0 );
	$cmd = "sga correct -k 31 -i 10 -t 4 -o $ec_file $WD/$outbase.pp.fastq > $WD/correct.out";
	print STDERR "[a5] $cmd\n";
	`rm $WD/$OUTBASE.pp.*` if (-f "$WD/$OUTBASE.pp.*");
	`rm $OUTBASE.pp.*` if (-f "$OUTBASE.pp.*");
	system("$DIR/$cmd");
	
	die "[a5] Error correcting reads with SGA\n" if( $? != 0 );
	return $ec_file;
}

#
# return a hash of hashes index as such: $hash{$library}{$component}
# where $component = 'up' || 'p1' || 'p2' || 'ins' 
#
sub read_lib_file {
	my $libfile = shift;
	my %libs = ();
	open(LIB,"<",$libfile);
	my $lib_count = 0;
	my $id = "";
	my %hash = ();
	while(<LIB>){
		chomp;
		if ($_ =~ m/\[LIB\]/){
			if ($lib_count > 0) {
				$id = "raw$lib_count" unless (length($id)); # give a generic name unless user has specified one
				for my $key (keys %hash){
					$libs{$id}{$key} = $hash{$key};
					delete($hash{$key});
				}
			} 
			$lib_count++;
			$id = "raw$lib_count";
		} elsif ($_ =~ m/shuf=([\w\/\-\.]+)/) { 
			my ($fq1, $fq2) = split_shuf($1,"$OUTBASE.lib$lib_count");
			$hash{"p1"} = $fq1;
			$hash{"p2"} = $fq2;
		} elsif ($_ =~ m/(p[1,2])=([\w\/\-\.]+)/) { 
			$hash{$1} = $2;
		} elsif ($_ =~ m/up=([\w\/\-\.]+)/) { 
			$hash{"up"} = $1;
		} elsif ($_ =~ m/ins=([\w\/\-\.]+)/) { 
			$hash{"ins"} = $1;
		} elsif ($_ =~ m/id=([\w\/\-\.]+)/) { 
			$id = $1;
			$hash{"id"} = $1;
		} else {
			die "[a5] Unrecognizable line in library file: >$_<\n";
		}
	} 
	$id = "raw$lib_count" unless (length($id)); # give a generic name unless user has specified one
	for my $key (keys %hash){
		$libs{$id}{$key} = $hash{$key};
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
	open(FQ,"<$fastq");
	open(FA,">$fasta");
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
	return ($fasta, $maxrdlen - 2); # -1 removes newline char
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

sub tagdust {
	my $outbase = shift;
	my $readfile = shift;
	my $tagdust_cmd = "tagdust -s -o $WD/$outbase.dusted.fq $DIR/../adapter.fasta $readfile";
	print STDERR "[a5] $tagdust_cmd\n";
	system("$DIR/$tagdust_cmd");
	return "$WD/$outbase.dusted.fq";
}

# expects a file called $outbase.clean.fa in the current working directory
sub idba_assemble {
	my $outbase = shift;
	my $reads = shift;
	my $maxrdlen = shift;
	$maxrdlen = 90 if $maxrdlen > 90;	# idba seems to break if the max k gets too big
	my $idba_cmd = "idba -r $reads -o $WD/$outbase --mink 29 --maxk $maxrdlen";
	print STDERR "[a5] $idba_cmd\n";
	`$DIR/$idba_cmd > $WD/idba.out`;
	die "[a5] Error building contigs with IDBA\n" if ($? != 0);
	`rm $WD/$outbase.kmer $WD/$outbase.graph`;
	return "$WD/$outbase-contig.fa";
}



sub scaffold_sspace {
	my $libfile = shift;
	my $outbase = shift;
	# build library file
	my $libsref = shift;
	my %libs = %$libsref;
	my $curr_ctgs = shift;
	my @library_file = ();
	my @lib_files;

	my $genome_size = get_genome_size($curr_ctgs);
	print STDERR "[a5] Total contig length $genome_size\n";
	my $libraryI=1;
	my @curr_lib_file = ();
	my $curr_ins = -1;
	my $prev_reads = "";
	my %run_lib;
	# sort library.txt to so that we scaffold with smaller insert libraries first
	for my $lib (sort { $libs{$a}{"ins"} <=> $libs{$b}{"ins"} } keys %libs) {
		my ($exp_link, $cov) = calc_explinks( $genome_size, $libs{$lib}{"ins"}, $libs{$lib}{"p1"} ); 
		printf STDERR "[a5] %s\: Insert %.0f, coverage %.2f, expected links %.0f\n", $libs{$lib}{"id"}, $libs{$lib}{"ins"}, $cov, $exp_link;
		if (defined($libs{$lib}{"up"})) { # run sspace with unpaired library if we have one
			$curr_ctgs = run_sspace($genome_size, $libs{$lib}{"ins"}, $exp_link, $outbase.".".$libs{$lib}{"id"},
                                                  $libs{$lib}{"libfile"}, $curr_ctgs, $libs{$lib}{"up"});
		} else {
			$curr_ctgs = run_sspace($genome_size, $libs{$lib}{"ins"}, $exp_link, $outbase.".".$libs{$lib}{"id"},
                                                  $libs{$lib}{"libfile"}, $curr_ctgs);
		}
	}

	return $curr_ctgs;
}

# Merge libraries with similar inserts to avoid gaps from complementary coverage
sub preprocess_libs {
	my $libsref = shift;
	my %libs = %$libsref;
	my $ctgs = shift;
	if (-f "$OUTBASE.unpaired.fastq"){
		print STDERR "[a5] Removing unpaired reads $OUTBASE.unpaired.fastq\n";
		`rm $OUTBASE.unpaired.fastq`;
	}
#
#       MAKE OUR INITIAL ESTIMATES OF INSERT SIZE HERE AND GET ERROR ESTIMATES.
#       ALSO CONCATENATE UNPAIRED LIBRARIES AND REMOVE FROM HASH
#
	print STDERR "[a5] Making initial estimates of insert size\n";
	for my $libid (keys %libs) {
		if (defined($libs{$libid}{"p1"})) {
			my $fq1 = $libs{$libid}{"p1"};
			my $fq2 = $libs{$libid}{"p2"};
			my ($ins_mean, $ins_err, $outtie) = get_insert($fq1,$fq2,"$OUTBASE.$libid",$ctgs);
			if ($ins_mean > 0) {
				$libs{$libid}{"ins"} = $ins_mean;
				$libs{$libid}{"err"} = $ins_err;
				$libs{$libid}{"rc"} = $outtie;
			} 
		} else {
			my $up = $libs{$libid}{"up"};
			`cat $up >> $OUTBASE.unpaired.fastq`;
			delete($libs{$libid}); # isn't paired, don't include in aggregated libraries
		}
	}
#
#       NOW MERGE SIMILIAR LIBRARIES
#
	my @curr_lib_file = ();
	my $curr_lib = "";
	my $prev_lib;
	my $prev_reads = "";
	my $libraryI = 1;
	my %processed = ();
	my $run_lib;
	print STDERR "[a5] Will merge libraries if similar enough\n";
	# sort library.txt to so that we scaffold with smaller insert libraries first
	for my $lib (sort { $libs{$a}{"ins"} <=> $libs{$b}{"ins"} } keys %libs) {
		# if we've hit a substantially different insert size, do a separate
		# round of scaffolding
		my $curr_ins = $libs{$lib}{"ins"};
		my $prev_ins = defined($prev_lib) ? $libs{$prev_lib}{"ins"} : -1;
		my $prev_min = -1;
		my $prev_min = -1;
		my $prev_min = -1;
		my $prev_min = -1;
		if (defined($prev_lib)) {
			
		}

		print STDERR "[a5] \$prev_ins = $prev_ins, \$curr_ins = $curr_ins\n";
		if($prev_ins > 0 && abs(log($prev_ins)-log($curr_ins)) > 0.1){
			print STDERR "[a5] abs(log($prev_ins)-log($curr_ins)) = ".abs(log($prev_ins)-log($curr_ins))."\n";
			# scaffold with the previous insert....
			# combine libraries if necessary, and return just one library hash
			$run_lib = aggregate_libs(\@curr_lib_file,$curr_lib,$ctgs);
			$run_lib->{"libfile"} = print_libfile("library_$libraryI.txt", $run_lib);
			for my $key (keys %$run_lib){
				$processed{$run_lib->{"id"}}{$key} = $run_lib->{$key};
			}
			# now move on to the next library...
			$libraryI++;
			@curr_lib_file = ();
			$curr_lib = "";
		}
		$prev_reads = $libs{$lib}{"p1"};
		$prev_ins = $curr_ins;
		push (@curr_lib_file,$libs{$lib});
		$curr_lib .= $lib;
		$prev_lib = $lib;
	}
	$run_lib = aggregate_libs(\@curr_lib_file,$curr_lib,$ctgs);
	$run_lib->{"libfile"} = print_libfile("library_$libraryI.txt", $run_lib);
	for my $key (keys %$run_lib){
		$processed{$run_lib->{"id"}}{$key} = $run_lib->{$key};
	}
	return %processed;
}

sub aggregate_libs {
	my $curr_lib_file = shift;
	my $curr_lib = shift;
	my $curr_ctgs = shift;
	my ($fq1, $fq2,$up);
	my %fin_lib = ();
	if (scalar(@$curr_lib_file) > 1) { # if this is an aggregate of libraries, combine them into one file
		($fq1, $fq2, $up) = merge_libraries($curr_lib_file,$curr_lib);
	} else {
		($fq1, $fq2) = ($curr_lib_file->[0]{"p1"}, $curr_lib_file->[0]{"p2"});
		$up = $curr_lib_file->[0]{"up"} if (defined($curr_lib_file->[0]{"up"}));
	}
	$fin_lib{"id"} = $curr_lib;
	$fin_lib{"p1"} = $fq1;
	$fin_lib{"p2"} = $fq2;
	if (defined($up)){
		$fin_lib{"up"} = $up;
		# add unpaired libraries to unpaired reads from this library aggregate
		`cat $OUTBASE.unpaired.fastq >> $fin_lib{"up"}` if -f "$OUTBASE.unpaired.fastq";
	} elsif (-f "$OUTBASE.unpaired.fastq") {
		$up = $fin_lib{"id"}.".up.fastq";	
		`cat $OUTBASE.unpaired.fastq > $up`;
		$fin_lib{"up"} = $up;
	}
	$fin_lib{"rc"} = 0;
	# re-estimate insert sizes
	my ($ins_mean, $ins_err, $outtie) = get_insert($fq1,$fq2,"$OUTBASE.$curr_lib",$curr_ctgs);
	$fin_lib{"ins"} = abs($ins_mean);
	$fin_lib{"err"} = abs($ins_err);	
	$fin_lib{"nlibs"} = scalar(@$curr_lib_file);
	return \%fin_lib;		
}

sub print_libfile {
	my $file = shift;
	my $libref = shift;
	#open( LIBRARY, ">$WD/library_$libraryI.txt" );
	open( LIBRARY, ">$file" );
	#print LIBRARY "$curr_lib $fq1 $fq2 $ins_mean $ins_err $outtie\n";
	print LIBRARY $libref->{"id"}." ".$libref->{"p1"}." ".$libref->{"p2"}." ".
                  $libref->{"ins"}." ".$libref->{"err"}." ".$libref->{"rc"}."\n";
	close LIBRARY;
	return $file;
}

sub merge_libraries {
	my $curr_lib_file = shift; # array of hashes
	my $curr_lib = shift;
	my $fq1 = "$curr_lib\_p1.fastq";
	my $fq2 = "$curr_lib\_p2.fastq";
	my $up = "$curr_lib\_up.fastq";
	open(my $fq1h,">","$fq1");
	open(my $fq2h,">","$fq2");
	my $uph;
	# merge files for scaffolding, reverse complementing as necessary
	print STDERR "[a5] Merging libraries";
	for my $sublib (@$curr_lib_file){
		my @ar  = split(' ',$sublib);
		print STDERR " ".$sublib->{"id"};
		open(my $fq,"<",$sublib->{"p1"});
		pipe_fastq($fq,$fq1h,$sublib->{"rc"});
		open($fq,"<",$sublib->{"p2"});
		pipe_fastq($fq,$fq2h,$sublib->{"rc"});
		if (defined($sublib->{"up"})){  # build and unpaired file if there are unpaired reads
			if (defined($uph)){
				open($uph,">","$up");
			}
			open($fq,"<",$sublib->{"up"});
			pipe_fastq($fq,$uph,$sublib->{"rc"});
		}
	}	
	close $fq1h;
	close $fq2h;
	close $uph;
	return ($fq1, $fq2, $up);
}

sub break_all_misasms {
	my $ctgs = shift;
	my $libsref = shift;
	my %libs = %$libsref;
	my $outbase = shift;
	my @lib_files;
	# sort libraries by insert so we break with larger inserts first.
	for my $lib (sort { $libs{$a}{"ins"} <=> $libs{$b}{"ins"} } keys %libs) {
		my $ins = $libs{$lib}{"ins"};
		my $fq1 = $libs{$lib}{"p1"};
		my $fq2 = $libs{$lib}{"p2"};
		$ctgs = fish_break_misasms($ctgs,$fq1,$fq2,"$outbase.lib$lib");
	}	
	return $ctgs;
}

sub fish_break_misasms {
	my $ctgs = shift;
	my $fq1 = shift;
	my $fq2 = shift;
	my $outbase = shift;
	print STDERR "[a5] Identifying misassemblies in $ctgs with $outbase\n";
	my $sai = "$WD/$outbase.sai";
	my $sam = "$WD/$outbase.sam";
#	print STDERR "[a5] Attempting to run bwa from $DIR\n";	
	`$DIR/bwa index -a is $ctgs > $WD/$outbase.index.out`;
	`cat $fq1 $fq2 | $DIR/bwa aln $ctgs - > $sai`;
	`cat $fq1 $fq2 | $DIR/bwa samse $ctgs $sai - > $sam`;
	`rm $ctgs.* $sai`;
	my $mem = "7000m";
	my $cmd = "GetFishInput.jar $sam $outbase $WD > $WD/$outbase.fie.out\n";
	print STDERR "[a5] java -Xmx$mem -jar $cmd\n"; 
	`java -Xmx$mem -jar $DIR/$cmd`;
	`gzip -f $sam`;
	$cmd = "fish -off -f $WD/$outbase.control.txt -b $WD/$outbase.blocks.txt > $WD/$outbase.fish.out";
	print STDERR "[a5] $cmd\n"; 
	`$DIR/$cmd`;
	die "[a5] Error getting blocks with FISH for $outbase\n" if ($? != 0);

	$cmd = "break_misassemblies.pl $WD/$outbase.blocks.txt $WD/$outbase.contig_labels.txt $ctgs > $WD/$outbase.broken.fasta 2> $WD/$outbase.break.out";
	print STDERR "[a5] $cmd\n"; 
	`$DIR/$cmd`;
	die "[a5] Error getting breaking contigs after running FISH\n" if ($? != 0);
	return "$WD/$outbase.broken.fasta";
}

sub pipe_fastq {
	my $from = shift;
	my $to = shift;
	my $rc = shift;
	while (my $line = <$from>) {
		print $to $line;
		$line = <$from>;
		if ($rc){
			chomp $line;
			$line =~ tr/ACGTacgt/TGCAtgca/; # complement
			$line = reverse($line);         # reverse
			print $to $line."\n";
			print $to readline($from);
			$line = <$from>;
			chomp $line;
			$line = reverse($line);         # reverse
			print $to $line."\n";
		} else {
			print $to readline($from);
			print $to readline($from);
			print $to readline($from);
		}
	}
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
	my $exp_link = $cov * $ins_len / $maxrdlen;
	return ($exp_link, $cov);
}
#run_sspace $genome_size $insert_size $exp_links $libfile $input_fa $unpaired_optional
sub run_sspace {
	my $genome_size = shift;
	my $insert_size = shift;
	my $exp_links = shift;
	my $outbase = shift;
	my $libfile = shift;
	my $input_fa = shift;

	my $sspace_m = int(log2($genome_size)+3.99);
	$sspace_m = 15 if $sspace_m < 15;
	my $sspace_n = int(log2($insert_size)*1.25+.99);
	my $sspace_k = int(log($exp_links)/log(1.4)-9.5);
	# require at least 2 links to preclude chimeras
	$sspace_k = $sspace_k < 2 ? 2 : $sspace_k;	
#	$sspace_k = 5;	# gives best results on mediterranei

	my $sspace_cmd = "SSPACE -m $sspace_m -n $sspace_n -k $sspace_k -a 0.2 -o 1 ".
                                        "-l $libfile -s $input_fa -b $outbase -d $WD";
	if (@_) {
		print STDERR "[a5] Running SSPACE with unpaired reads\n";
		my $up = shift;
		$sspace_cmd .= " -u $up";
	}
	print STDERR "[a5] $sspace_cmd > $WD/$outbase.out\n";
	`$DIR/SSPACE/$sspace_cmd > $WD/$outbase.out`;
	`rm -rf $WD/bowtieoutput/ $WD/reads/`;
	return "$WD/$outbase.final.scaffolds.fasta";
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
	my $cmd = "$DIR/bwa sampe -P $ctgs $r1fq.sub.sai $r2fq.sub.sai $r1fq.sub $r2fq.sub ".
												"> $WD/$outbase.sub.pe.sam 2> $WD/$outbase.sampe.out";
	`$cmd`;
	open(SAMPE, "<","$WD/$outbase.sampe.out");
	my $ins_mean = 0;
	my $ins_error = 0;
	my $min;
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
		} 
	}
	`rm $r1fq.sub* $r2fq.sub* $ctgs.*`;
	if ($ins_n > $require_fraction * $estimate_pair_count) {		
		$ins_error = $ins_sd*6 / $ins_mean;
		$ins_mean = sprintf("%.0f",$ins_mean);
		$ins_error = sprintf("%.3f",$ins_error);
		$ins_error =~ s/0+$//g;
	} else {
		print STDERR "[a5] Discarding estimate. Not enough data points: $ins_n\n";
		$ins_mean *= -1; 
		$ins_error *= -1;
	}
	my $orient = get_orientation("$WD/$outbase.sub.pe.sam");
	return ($ins_mean, $ins_error, $orient);
}

sub get_orientation {
	my $samfile = shift;
	open(SAM,"<$samfile");
	my %reads = ();
	my $in = 0;
	my $in_ins = 0;
	my $out = 0;
	my $out_ins = 0;
	my $tot = 0;
	while(<SAM>){
		next if ($_ =~ /^@/);
		$tot++;
		my @r1 = split /\t/;
		# skip this read pair if either the mate or the query is unmapped
		next if (get_bit($r1[1],2));
		next if (get_bit($r1[1],3));
		if (defined($reads{$r1[0]})){
			my $r2 = delete($reads{$r1[0]});
			next unless($r1[2] eq $r2->[2]);
		#	print STDERR "[a5] found read pair $r1[0] on $r1[2]\n";
			if ($r1[3] < $r2->[3]){
		#		print STDERR "[a5] r1 is upstream of r2\n";
				my $s1 = get_bit($r1[1],4);
				my $s2 = get_bit($r2->[1],4);
				# outie orientation if upstream read maps to reverse strand
				if($s1 && !$s2){
					$out++;
					$out_ins += abs($r1[1]);
				} elsif (!$s1 && $s2){
					$in++;
					$in_ins += abs($r1[1]);
				}
			} else {
		#		print STDERR "[a5] r1 is upstream of r2\n";
				my $s1 = get_bit($r1[1],4);
				my $s2 = get_bit($r2->[1],4);
				# outie orientation if upstream read maps to reverse strand
				if(!$s1 && $s2){
					$out++;
					$out_ins += abs($r1[1]);
				} elsif ($s1 && !$s2){
					$in++;
					$in_ins += abs($r1[1]);
				}
			}
		} else {
			$reads{$r1[0]} = \@r1;
		}
	}
	$tot /= 2;
	$out_ins /= $out if ($out);
	$in_ins /= $in if ($in);
	print STDERR "[a5] get_orientation: $samfile n_tot = $tot n_out = $out ins_out = $out_ins n_in = $in ins_in = $in_ins\n";
	# if majority outies or we have  
	if (($out_ins > 1500 && $out/$tot > 0.1) || $in < $out){
		return 1;
	} else {
		return 0;
	}	
}

sub get_bit {
	my $mod = 0;
	my $flag = shift;
	my $dig = shift;
	my $i = 0;
	while($flag != 0){
		$mod = $flag % 2;
		$flag = int($flag/2);
		return $mod if ($i == $dig); 
		$i++;
	}
	return 0;
}

sub get_rdlen($$){
	open(FILE, shift);
	my $line = <FILE>;
	$line = <FILE>;
	chomp $line;
	close FILE;
	return length($line);
}

