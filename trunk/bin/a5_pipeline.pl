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
my $outbase = $ARGV[1];

my $DIR = dirname(abs_path($0));
my %LIBS = read_lib_file($libfile);
die "[a5] No libraries found in $libfile\n" unless(keys %LIBS);
print "[a5] Found ".scalar(keys %LIBS)." libraries\n";
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
	$reads = sga_clean($outbase, \%LIBS);
	$reads = tagdust($outbase, $reads);
	`rm $WD/$outbase.pp.*`
	`mv $reads $outbase.ec.fastq`;
} 
if ($start <= 2) {
	$reads = "$outbase.ec.fastq";
	die "[a5_s2] Can't find error corrected reads $reads" unless -f $reads;
	print "[a5_s2] Building contigs from $reads with IDBA\n";
	print STDERR "[a5_s2] Building contigs from $reads with IDBA\n";
	$WD="s2";
	mkdir($WD) if ! -d $WD;
	$maxrdlen = fastq_to_fasta($reads,"$WD/$outbase.ec.fasta");
	$ctgs = idba_assemble($outbase, "$WD/$outbase.ec.fasta", $maxrdlen); 
	`gzip -f $reads ` 
	`rm $WD/$outbase.ec.fasta`;
	`mv $ctgs $outbase.contigs.fasta`;
	
} 
if ($start <= 3) {
	$ctgs = "$outbase.contigs.fasta";
	die "[a5_s3] Can't find starting contigs $ctgs.\n" unless -f $ctgs;
	print "[a5_s3] Scaffolding contigs from $ctgs with SSPACE\n";
	print STDERR "[a5_s3] Scaffolding contigs from $ctgs with SSPACE\n";
	$WD="s3";
	mkdir($WD) if ! -d $WD;
	$scafs = scaffold_sspace($libfile, $outbase, \%LIBS, $ctgs);
	`mv $scafs $outbase.crude.scaffolds.fasta`; 
} 
if ($start <= 4) {
	$scafs = "$outbase.crude.scaffolds.fasta";
	die "[a5_s4] Can't find starting crude scaffolds $scafs\n" unless -f $scafs;
	print "[a5_s4] Detecting and breaking misassemblies in $scafs with FISH\n";
	print STDERR "[a5_s4] Detecting and breaking misassemblies in scafs with FISH\n";
	$scafs = break_all_misasms($scafs,\%LIBS,"$outbase.fish"); 
	`mv $scafs $outbase.broken.scaffolds.fasta`; 
} 
if ($start <= 5) {
	$scafs = "$outbase.broken.scaffolds.fasta";
	die "[a5_s5] Can't find starting broken scaffolds $scafs\n" unless -f $scafs;
	print "[a5_s5] Scaffolding broken contigs with SSPACE\n";
	print STDERR "[a5_s5] Scaffolding broken contigs with SSPACE\n";
	$scafs = scaffold_sspace($libfile,"$outbase.rescaf",\%LIBS,$scafs);
	`mv $scafs $outbase.final.scaffolds.fasta`;
	print "[a5] Final assembly in $outbase.final.scaffolds.fasta\n"; 
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
		my %lib_files = %{$libs->{$lib}};
		if (defined($lib_files{"p1"})) {
			#shift @lib_files;
			my $fq1 = $lib_files{"p1"}; #shift @lib_files;
			my $fq2 = $lib_files{"p2"}; #shift @lib_files;
			$files .= "$fq1 $fq2 ";
		}
		if (defined($lib_files{"up"})){
			my $up = $lib_files{"up"}; #shift @lib_files;
			$files .= "$up ";
		}
	}
	print STDERR "[a5] sga preprocess -q 10 -f 20 -m 30 --phred64 $files > $WD/$outbase.pp.fastq\n";
	system("$DIR/sga preprocess -q 10 -f 20 -m 30 --phred64 $files > $WD/$outbase.pp.fastq");
	die "[a5] Error preprocessing reads with SGA\n" if( $? != 0 );
	my $sga_ind = "";
	my $sga_ind_kb = 4000000;
	do{
		$sga_ind = `$DIR/sga index -d $sga_ind_kb -t 4  $WD/$outbase.pp.fastq > $WD/index.out 2> $WD/index.err`;
		$sga_ind = read_file('index.err');
		$sga_ind_kb = int($sga_ind_kb/2);
	}while(($sga_ind =~ /bad_alloc/ || $? != 0) && $sga_ind_kb > 0);
	system("rm -f core*") if (-f "core*");
	die "[a5] Error indexing reads with SGA\n" if( $? != 0 );
	system("$DIR/sga correct -k 31 -i 10 -t 4  $WD/$outbase.pp.fastq > $WD/correct.out");
	die "[a5] Error correcting reads with SGA\n" if( $? != 0 );
	return "$WD/$outbase.pp.ec.fa"
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
	my %hash = ();
	my $id = "";
	while(<LIB>){
		chomp;
		if ($_ =~ m/\[LIB\]/){
			if ($lib_count > 0) {
				$id = "raw$lib_count" unless (length($id));
				$libs{$id} = %hash;
#				for my $key (sort keys %hash){
#					push($libs{"raw$lib_count"}{,$hash{$key});
#				}
			} 
			$lib_count++;
			$id = "";
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
		} elsif ($_ =~ m/id=([\w\/\-\.]+)/) { 
			$id = $1;
			$hash{"id"} = $1;
		} else {
			die "[a5] Unrecognizable line in library file: >$_<\n";
		}
	} 
	$id = "raw$lib_count" unless (length($id));
	$libs{$id} = %hash;
	$libs{"raw$lib_count"} = %hash;
#	for my $key (sort keys %hash){
#		push(@{$libs{"raw$lib_count"}},$hash{$key});
#	}
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

sub tagdust {
	my $outbase = shift;
	my $readfile = shift;
	my $tagdust_cmd = "$DIR/tagdust -o $WD/$outbase.dusted.fq $DIR/../adapter.fasta $readfile";
	print STDERR "[a5] $tagdust_cmd\n";
	system($tagdust_cmd);
	return "$WD/$outbase.dusted.fq";
}

# expects a file called $outbase.clean.fa in the current working directory
sub idba_assemble {
	my $outbase = shift;
	my $reads = shift;
	my $maxrdlen = shift;
	$maxrdlen = 90 if $maxrdlen > 90;	# idba seems to break if the max k gets too big
	my $idba_cmd = "$DIR/idba -r $reads -o $WD/$outbase --mink 29 --maxk $maxrdlen";
	print STDERR "[a5] $idba_cmd\n";
	`$idba_cmd > $WD/idba.out`;
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
#vestigial?	my $min_insert_size = -1;
	# remove the unpaired file if present 
	`rm $outbase.unpaired.fastq` if (-f "$outbase.unpaired.fastq");

#
#       MAKE OUR INITIAL ESTIMATES OF INSERT SIZE HERE AND GET ERROR ESTIMATES
#
	for my $libid (keys %libs) {
		if (defined($libs{$libid}{"p1"})) {
			my $fq1 = $libs{$libid}{"p1"};
			my $fq2 = $libs{$lib}{"p2"};
			my ($ins_mean, $ins_err, $outtie) = get_insert($fq1,$fq2,"$outbase.$lib",$curr_ctgs);
			if ($ins_mean > 0) {
				$libs{$libid}{"ins"} = $ins_mean;
				$libs{$libid}{"err"} = $ins_err;
				$libs{$libid}{"rc"} = $outtie;
			} 
		}
		# if we have one file left, treat it as unpaired. 
		if (defined($libs{$libid}{"up"})) {
			my $up = $libs{$libid}{"up"};
			`cat $up >> $outbase.unpaired.fastq`;
			delete($libs{$libid});
		}
	}
#
#       NOW MERGE LIBRARIES WITH SIMILIAR ENOUGH INSERT SIZE DISTRIBUTIONS
#
	my $genome_size = get_genome_size($curr_ctgs);
	print STDERR "[a5] Total contig length $genome_size\n";
	my $libraryI=1;
	my @curr_lib_file = ();
	my $curr_lib = "";
	my $prev_ins = -1;
	my $prev_reads = "";
	my $libfile_path;
	my $up_path;
	my %run_lib;
	# sort library.txt to so that we scaffold with smaller insert libraries first
	for my $lib (sort { $libs{$a}{"ins"} <=> $libs{$b}{"ins"} } keys %libs) {
		# if we've hit a substantially different insert size, do a separate
		# round of scaffolding
		my $cur_ins = $libs{$lib}{"ins"};
		if($prev_ins > 0 && abs(log($prev_ins)-log($cur_ins)) > 0.1){
			# scaffold with the previous insert....
			# combine libraries if necessary, and return just one library hash
			%run_lib = prep_libs_sspace(\@curr_lib_file,$libraryI,$outbase,$curr_lib,$curr_ctgs);
			$libfile_path = print_libfile("$WD/library_$libraryI.txt", %run_lib);
			my ($exp_link, $cov) = calc_explinks( $genome_size, $prev_ins, $prev_reads ); 
			print STDERR "[a5] $curr_lib -  Insert $prev_ins, coverage $cov, expected links $exp_link\n";
			if (defined($run_lib{"up"})) {
				# add unpaired libraries to unpaired reads from this library aggregate
				`cat $outbase.unpaired.fastq >> $run_lib{"up"}` if -f "$outbase.unpaired.fastq";
				$curr_ctgs = run_sspace($genome_size, $prev_ins, $exp_link, $libfile_path, $curr_ctgs, $run_lib{"up"});
				 # if we combined libraries, removed the concatenated files
				`rm $run_lib{"p1"} $run_lib{"p2"} $run_lib{"up"}` if (scalar(@curr_lib_file) > 1);
			} else {
				$curr_ctgs = run_sspace($genome_size, $prev_ins, $exp_link, $libfile_path, $curr_ctgs);
				 # if we combined libraries, removed the concatenated files
				`rm $run_lib{"p1"} $run_lib{"p2"}` if (scalar(@curr_lib_file) > 1);
			}
			# now move on to the next library...
			$libraryI++;
			@curr_lib_file = ();
			$curr_lib = "";
		}
		$prev_reads = $libs{$lib}{"p1"};
		$prev_ins = $cur_ins;
		push (@curr_lib_file,$libs{$lib});
		$curr_lib .= $lib;
	}
	# combine libraries if necessary, and return just one library hash
	%run_lib = prep_libs_sspace(\@curr_lib_file,$libraryI,$outbase,$curr_lib,$curr_ctgs);
	$libfile_path = print_libfile("$WD/library_$libraryI.txt", %run_lib);
	my ($exp_link, $cov) = calc_explinks( $genome_size, $prev_ins, $prev_reads ); 
	print STDERR "[a5] $curr_lib -  Insert $prev_ins, coverage $cov, expected links $exp_link\n";
	if (defined($run_lib{"up"})) {
		# add unpaired libraries to unpaired reads from this library aggregate
		`cat $outbase.unpaired.fastq >> $run_lib{"up"}` if -f "$outbase.unpaired.fastq";
		$curr_ctgs = run_sspace($genome_size, $prev_ins, $exp_link, $libfile_path, $curr_ctgs, $run_lib{"up"});
		 # if we combined libraries, removed the concatenated files
		`rm $run_lib{"p1"} $run_lib{"p2"} $run_lib{"up"}` if (scalar(@curr_lib_file) > 1);
	} else {
		$curr_ctgs = run_sspace($genome_size, $prev_ins, $exp_link, $libfile_path, $curr_ctgs);
		 # if we combined libraries, removed the concatenated files
		`rm $run_lib{"p1"} $run_lib{"p2"}` if (scalar(@curr_lib_file) > 1);
	}

	return $curr_ctgs;
}

sub prep_libs_sspace {
	my $curr_lib_file = shift;
	my $libraryI = shift;
	my $outbase = shift;
	my $curr_lib = shift;
	my $curr_ctgs = shift;
	my ($fq1, $fq2,$up);
	my %fin_lib = ();
	if (scalar(@$curr_lib_file) > 1) { # if this is an aggregate of libraries, combine them into one file
		($fq1, $fq2, $up) = merge_libraries($curr_lib_file,$curr_lib);
	} else {
		($fq1, $fq2) = ($curr_lib_file->[0]{"p1"}, $curr_lib_file->[0]{"p2"})
		$up = $curr_lib_file->[0]{"up"} if (defined($curr_lib_file->[0]{"up"}));
	}
	$fin_lib{"id"} = $curr_lib;
	$fin_lib{"p1"} = $fq1;
	$fin_lib{"p2"} = $fq2;
	$fin_lib{"up"} = $up if defined($up);
	$fin_lib{"rc") = 0;
	my ($ins_mean, $ins_err, $outtie) = get_insert($fq1,$fq2,"$outbase.$curr_lib",$curr_ctgs);
	$fin_lib{"ins"} = abs($ins_mean);
	$fin_lib{"err"} = abs($ins_err);	
	return %fin_lib;		
}

sub print_libfile {
	my $file = shift;
	my $lib = shift;
	#open( LIBRARY, ">$WD/library_$libraryI.txt" );
	open( LIBRARY, ">$file" );
	#print LIBRARY "$curr_lib $fq1 $fq2 $ins_mean $ins_err $outtie\n";
	print LIBRARY $lib->{"id"}." ".$lib->{"p1"}." ".$lib->{"p2"}." ".
                  $lib->{"ins"}." ".$lib->{"err"}." ".$lib->{"rc"}."\n";
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
	my $libs = shift;
	my $outbase = shift;
	my @lib_files;
	my %tmp = %$libs;
	for my $lib (keys %$libs){
		@{$tmp{$lib}} = @{$libs->{$lib}} if (scalar(@{$libs->{$lib}}) >= 2);
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
	print STDERR "[a5] Breaking misassemblies with $outbase\n";
	my $sai = "$outbase.sai";
	my $sam = "$outbase.sam";
#	print STDERR "[a5] Attempting to run bwa from $DIR\n";	
	`$DIR/bwa index -a is $ctgs > $outbase.index.out`;
	`cat $fq1 $fq2 | $DIR/bwa aln $ctgs - > $sai`;
	`cat $fq1 $fq2 | $DIR/bwa samse $ctgs $sai - > $sam`;
	`rm $ctgs.*`;
	my $cmd = "java -Xmx7000m -jar $DIR/GetFishInput.jar $sam $outbase > $outbase.fie.out\n";
	print STDERR "[a5] $cmd\n"; 
	`$cmd`;
	`gzip -f $sam`;
	$cmd = "$DIR/fish -off -f $outbase.control.txt -b $outbase.blocks.txt > $outbase.fish.out";
	print STDERR "[a5] $cmd\n"; 
	`$cmd`;
	die "[a5] Error getting blocks with FISH for $outbase\n" if ($? != 0);

	$cmd = "$DIR/break_misassemblies.pl $outbase.blocks.txt contig_labels.txt $ctgs > $outbase.broken.fasta 2> $outbase.break.out";
	print STDERR "[a5] $cmd\n"; 
	`$cmd`;
	die "[a5] Error getting breaking contigs after running FISH\n" if ($? != 0);
	#`rm $outbase.blocks.txt contig_labels.txt fish.* get_fish_input.* break_misasm.err $sam $sai`;
	return "$outbase.broken.fasta";
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
	my $libfile = shift;
	my $input_fa = shift;

	my $sspace_m = int(log2($genome_size)+3.99);
	my $sspace_n = int(log2($insert_size)*1.25+.99);
	my $sspace_k = int(log($exp_links)/log(1.4)-9.5);
	# require at least 2 links to preclude chimeras
	$sspace_k = $sspace_k < 2 ? 2 : $sspace_k;	
#	$sspace_k = 5;	# gives best results on mediterranei

	#my $input_fa = "$outbase-contig.fa";
	# why was this done? 
	#$input_fa = "$outbase.lib".($libraryI-1).".sspace.final.scaffolds.fasta" if $libraryI>1;
	my $sspace_cmd = "$DIR/SSPACE/SSPACE -m $sspace_m -n $sspace_n -k $sspace_k -a 0.2 -o 1 ".
                                        "-l $libfile -s $input_fa -b $outbase -d $WD";
                                        #"-l $libfile -s $input_fa -b $outbase.lib$libraryI.sspace";
	if (@_) {
		print STDERR "[a5] Running SSPACE with unpaired reads\n";
		
		my $up = shift;
		$sspace_cmd .= " -u $up";
	}
	print STDERR "[a5] $sspace_cmd > $WD/$outbase.out\n";
	`$sspace_cmd > $WD/$outbase.out`;
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

sub get_insert($$$$$) {
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
												"> $WD/$outbase.pe.sam 2> $WD/$outbase.sampe.out";
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
	my $orient = get_orientation("$WD/$outbase.pe.sam");
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
	while(<SAM>){
		next if ($_ =~ /^@/);
		my @r1 = split /\t/;
		# skip this read pair if either the mate or the query is unmapped
		next if (get_bit($r1[4],2));
		next if (get_bit($r1[4],3));
		if (defined($reads{$r1[0]})){
			my $r2 = delete($reads{$r1[0]});
			next unless($r1[2] eq $r2->[2]);
			if ($r1[3] < $r2->[3]){
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
	$out_ins /= $out;
	$in_ins /= $in;
	my $tot = $in + $out;
	# if majority outies 
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
		return $mod if ($dig == $i); 
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

