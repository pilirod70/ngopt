#!/usr/bin/perl -w

my $del_file = shift;
my $id_file = shift;

open (IN,"<",$id_file);
my %jobids = ();
while(<IN>){
	chomp;
	my ($ccid, $jobid) = split(/\t/);
	$jobids{$ccid} = $jobid;
}

my $base = "svr_delete_RAST_job atritt halophile";

print $base;
open (IN,"<",$del_file);
while(<IN>) {
	chomp;
	my $jobid = $jobids{$_};
	print " $jobid";
}
print "\n";

