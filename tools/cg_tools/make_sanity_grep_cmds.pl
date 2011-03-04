#!/usr/bin/perl -w
use strict;
use warnings;


if (scalar(@ARGV) != 2) {
	print "Usage: make_sanity_grep_cmds.pl <org_id_file> <file_to_grep>\n";
	exit;
}

open (IN,"<",shift);
my $grep_file = shift;

my @ar1 = ();

while(<IN>){
	chomp;
	push(@ar1,$_);
}

my $len = scalar(@ar1);

for (my $i = 0; $i < $len; $i++){
	my $org = shift @ar1;
	foreach my $org2 (@ar1) {
		print "count=\`grep $org $grep_file | grep -c $org2\`\n";
		print "echo -e \"$org\\t$org2\\t\$count\"\n";  
	}
}

