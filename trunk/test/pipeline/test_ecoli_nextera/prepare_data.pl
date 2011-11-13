#!/usr/bin/perl -w
use strict;
use warnings;


my $count = 0;
while (my $hdr = <>){
	my $seq = <>;
	my $qhdr = <>;
	my $qual = <>;
	next unless $count++ % 10 == 0;
#	print STDERR length($seq)."\n".length($qual)."\n";
	chomp $qual;
	$qual = p33to64($qual);
	$hdr = (split(' ', $hdr))[1];
	print "\@$hdr/1\n".$seq."+$hdr/1\n".$qual."\n";	
}

sub p33to64 {
	my $qual = shift;
	my $ret = "";
	for (my $i = 0; $i < length($qual); $i++){
		my $char .= chr(ord(substr($qual,$i,1))+5);
		$ret .= $char;
		#print STDERR substr($qual,$i,1)."(".ord(substr($qual,$i,1)).") -> ".$char."(".ord($char).")\n";
	}
	return $ret;
}
