#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
sub update_count($);
if (scalar(@ARGV) != 3) {
	print "Usage: ".basename($0)." <ortho_group_prefix> <fasta_file> <groups.txt>\n";
	exit 1;
}
my $last_ortho_group = 10000;
my $id; 
my $prefix = $ARGV[0];
my $fa_file = $ARGV[1];
my $groups_file = $ARGV[2];
open (FA,"<",$fa_file);
my %cds = ();

while (<FA>) {
	chomp;
	if ($_=~m/^>/){
		$id = substr($_,1,3);
		$cds{substr($_,1)}=1;
		#print STDERR "I just hashed ".substr($_,1).".\n";
	}
}
print STDOUT "Organism: $id\n";
my $num_cds = scalar((keys %cds));
print STDOUT "Found ".$num_cds." CDS\n";
open (GROUPS,"<",$groups_file);
	
#  H4100001: HAA|HAA-001077 HAA|HAA-004244
while (<GROUPS>) {
	chomp;
	my @line = split;
	my $ortho_id = shift(@line);
#	print STDERR "Processing ortholog group with ".scalar(@line)." members\n";
	update_count($ortho_id);
	for my $feat_id (@line){
		my $org = substr($feat_id,0,index($feat_id,"|"));
#		print STDERR ">>$org<<\n";
		next if $org ne $id;
		if (defined($cds{$feat_id})){
			delete($cds{$feat_id});
#			print STDERR "I just deleted $feat_id\n";
		} else {
		#	print STDERR "Found $feat_id in $groups_file, but not in $fa_file\n";
		}
	}
}

open (GROUPS,">>",$groups_file);
my $num_singles = scalar((keys %cds));
print STDOUT "Found ".$num_singles." singletons\n";
my $count = 0;
foreach  (keys %cds) {
	$last_ortho_group++;
	print GROUPS $prefix.$last_ortho_group.": ".$_."\n";
	$count++;
}
print STDOUT "Found ".$count." singletons\n";
if ($num_cds == $num_singles) {
	print STDERR "$id was not found in $groups_file\n"; 
}
close GROUPS;

sub update_count($) {
	my $ortho_group = shift;
	if ($ortho_group =~ m/\.*(\d{5}):/){
		#my $count = substr($_,index($_,$prefix)+length($prefix),rindex($_,":"));
		if (int($1) > $last_ortho_group) {
			$last_ortho_group = int ($1);
		}
	} else {
		print STDERR "Unrecognizable ortholog id $ortho_group\n";
	}
}

