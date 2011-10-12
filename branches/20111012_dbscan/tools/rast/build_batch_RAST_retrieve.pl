#!/usr/bin/perl -w
use strict;
use warnings;

if (scalar(@ARGV) == 0) {
	print "Usage: $0 <format> <job_id_file> <output_dir>\n";
	exit;
}

my %extensions = ( genbank => 'gbk', embl => "embl", gff3 => "gff3", rast_tarball => 'tar.gz',
			       genbank_stripped => 'gbk', embl_stripped => "embl", gff3_stripped => "gff3"); 
		  	

my $user="atritt";
my $pass="halophile";
my $format=shift;
my $job_file=shift;
my $outdir=shift;
my $ext = $extensions{$format};
open(IN,"<",$job_file);
print "#!/bin/bash\n";
while(<IN>){
	chomp;
	my ($ccid,$jobid) = split /\t/;
	print "svr_retrieve_RAST_job $user $pass $jobid $format > $outdir/$ccid.$ext\n";
}


# Usage: /home/atritt/bin/RAST/plbin/svr_retrieve_RAST_job.pl username password jobid format > output-file
