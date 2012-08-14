#!/usr/bin/perl -w
if (scalar(@ARGV) != 2) {
	print "Usage: $0 <cc_id_to_del_file> <cc_id_job_id_file>\n";
	exit;
}

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

