#!/usr/bin/perl -w
use strict;
use warnings;




my $hdr_info_file = shift;

open (HDR,"<",$hdr_info_file);
my %output_fh = ();
my %hdr_base = ();
while(<HDR>){
	chomp;
	my @ar = split /\t/;
	$hdr_base{$ar[0]} = $ar[1];
	open (FILE,">",$hdr_base{$ar[0]}."_p1.fastq");
	$output_fh{"p1"} = *FILE;	
	open (FILE,">",$hdr_base{$ar[0]}."_p2.fastq");
	$output_fh{"p2"} = *FILE;	
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
		$pair_num = substr($hdr,0,index($hdr,"/"));
	} else { # Casava 1.8+
		@hdrAr = split(/' '/,$hdr);
		$pair_num = (split(/:/,$hdrAr[1]))[0];
		@hdrAr = split(/' '/,$hdrAr[0]);
		$lib = $hdrAr[0].":".$hdrAr[1].":".$hdrAr[2];
	}
	my $lib_hash;
	if (!defined($libs{$lib})){
		$libs{$lib} = ();
	} 
	my $seq = <IN>;
	<IN>:
	my $qual = <IN>;

	if (defined($libs{$lib}{$hdr})){ # this read is paired, and we found its mate
		my $fh = $output_fh{$lib}{"p$pair_num"};
		print $fh "\@$hdr/$pair_num\n".$seq."+$hdr/$pair_num\n".$qual;
		$fh = $output_fh{$lib}{"p".$libs{$lib}{$hdr}[0]};
		print $fh "\@$hdr/".$libs{$lib}{$hdr}[0]."\n".
                   $libs{$lib}{$hdr}[1].
                  "+$hdr/".$libs{$lib}{$hdr}[0]."\n".
                   $libs{$lib}{$hdr}[2];
		delete($libs{$lib}{$hdr});
	} else {
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

