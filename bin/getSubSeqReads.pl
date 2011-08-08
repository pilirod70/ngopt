#!/usr/bin/perl -w
use strict;
use warnings;

die "Usage: getSubSeqReads.pl <base> <libfile> <reference>\n" if (@ARGV != 3);

my $base = shift;
my $libfile = shift;
my $reffile = shift;

my %libs = read_lib_file($libfile);

`bwa index $reffile`;

my $files = "";
for my $lib (keys %libs) {
	my @lib_files = @{$libs{$lib}};
	if (scalar(@lib_files) >= 2) {
		shift @lib_files;
		my $fq1 = shift @lib_files;
		`bwa aln $reffile $fq1 > $base.$lib.p1.sai`;
		my $fq2 = shift @lib_files;
		`bwa aln $reffile $fq2 > $base.$lib.p2.sai`;
		#`bwa sampe $reffile $base.$lib.p1.sai $base.$lib.p2.sai $fq1 $fq2 > $base.$lib.pe.sam`;
		#`extract_mapped_w_mates.pl $base.$lib.pe.sam > $base.$lib.16S.sam`;
		my %reads = ();
		open(IN,"bwa sampe $reffile $base.$lib.p1.sai $base.$lib.p2.sai $fq1 $fq2 |");
		open(SAM,">$base.$lib.16S.sam");
		open(SAM,">$base.$lib.16S.sam");
		open(FQ,">$base.$lib.16S.fastq");
		while (<IN>){
			chomp;
			if ($_ =~ m/^@/) {
				print SAM $_."\n";
				next;
			}
			my @ar = split (/\t/, $_);
			if (defined($reads{$ar[0]})){
				my $r1 = $reads{$ar[0]};
				if ($r1->[5] =~ m/.*M.*/ || $ar[5] =~ m/.*M.*/){
					print SAM join("\t",@$r1)."\n".join("\t",@ar)."\n";
					print FQ "@".$r1->[0]."\n".$r1->[9]."\n+".$r1->[0]."\n".$r1->[10]."\n";
					print FQ "@".$ar[0]."\n".$ar[9]."\n+".$ar[0]."\n".$ar[10]."\n";
				} 
				delete($reads{$ar[0]});
			} else {
				$reads{$ar[0]} = \@ar;
			}
		}		
		close SAM;
		close FQ;
	}
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
