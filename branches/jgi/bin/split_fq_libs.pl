#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV != 2) {
	print "Usage: $0 <hdr_info_file> <fastq>\n";
	exit;
}


my $hdr_info_file = shift;

open (HDR,"<",$hdr_info_file);
my %output_fh = ();
my %hdr_base = ();
while(<HDR>){
	chomp;
	my @ar = split /\t/;
	$hdr_base{$ar[0]} = $ar[1];
	open ($output_fh{$ar[0]}{0},">",$hdr_base{$ar[0]}."_p1.fastq");
	open ($output_fh{$ar[0]}{1},">",$hdr_base{$ar[0]}."_p2.fastq");
}

my %libs = ();

my $fastq_file = shift;

open (IN,"<",$fastq_file);
		
while (my $line = <IN>) {
	my $hdr;
	chomp $line;
	my @hdrAr = split(/\t/,$line);
	if ($hdrAr[0] =~ m/\@SR.* (.*) length=\d*/){ # FastQ from SRA
		$hdr = $1;
	} else {
		$hdr = substr($hdrAr[0],1);
	}
	@hdrAr = split(/:/,$hdr);
	my $lib;
	my $pair_num;
	if (@hdrAr == 5) { # Illumina 1.4+
		$lib = $hdrAr[0].":".$hdrAr[1];
		$pair_num = substr($hdr,index($hdr,"/")+1);
		$hdr = substr($hdr,0,index($hdr,"/"));
	} else { # Casava 1.8+
		@hdrAr = split(/' '/,$hdr);
		$hdr = $hdrAr[0];
		$pair_num = (split(/:/,$hdrAr[1]))[0];
		@hdrAr = split(/' '/,$hdrAr[0]);
		$lib = $hdrAr[0].":".$hdrAr[1].":".$hdrAr[2];
	}
	my $lib_hash;
	if (!defined($libs{$lib})){
		$libs{$lib} = ();
	} 
	my $seq = <IN>;
	<IN>;
	my $qual = <IN>;

	if (defined($libs{$lib}{$hdr})){ # this read is paired, and we found its mate
		my $file_pair_key = $pair_num == 1 ? 0 : 1;
		print STDERR "Found match   >$lib<  >$file_pair_key<\n";
		my $fh = $output_fh{$lib}{$file_pair_key};
		print {$fh} "\@$hdr/$pair_num\n".$seq."+$hdr/$pair_num\n".$qual;

		$file_pair_key = $libs{$lib}{$hdr}[0] ? 0 : 1;
		$fh = $output_fh{$lib}{$file_pair_key};
		print {$fh} "\@$hdr/".$libs{$lib}{$hdr}[0]."\n".
                   $libs{$lib}{$hdr}[1].
                  "+$hdr/".$libs{$lib}{$hdr}[0]."\n".
                   $libs{$lib}{$hdr}[2];

		delete($libs{$lib}{$hdr});
	} else {
		print STDERR "$hdr\n";
		$libs{$lib}{$hdr} = [ $pair_num, $seq, $qual ];
	}
}

for my $lib (keys %libs){
	if (scalar(keys %{$libs{$lib}})) { # if there reads left, print them as unpaired.
		open (FILE,">",$hdr_base{$lib}."_up.fastq");
		for my $hdr (keys %{$libs{$lib}}){
			print FILE "\@$hdr/".$libs{$lib}{$hdr}[0]."\n".
                        $libs{$lib}{$hdr}[1].
                       "+$hdr/".$libs{$lib}{$hdr}[0]."\n".
                        $libs{$lib}{$hdr}[2];
		}
		close FILE;
	} 
	# close out the other filehandles
	for my $fh (values %{$output_fh{$lib}}){
		close $fh;
	}
}

