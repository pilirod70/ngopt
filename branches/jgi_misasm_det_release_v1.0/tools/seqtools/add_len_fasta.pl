#!/usr/bin/perl -w
use strict;
use warnings;

my $hdr = "";
my $seq = "";
my $len = 0;
while (<>){
	if ($_ =~ /^>/){
		if (length($hdr)){
			print  $hdr."_len_$len\n".$seq;
			$len = 0;
			$seq = "";
		}
		$hdr = $_;	
		chomp $hdr;
	} else {
		$seq .= $_;
		$len += length($_)-1;
	}
}

print  $hdr."_len_$len\n".$seq;
