#!/usr/bin/perl -w

use File::Basename;
die "Usage: ".basename($0)." <name_file> <fasta_dir> <fasta_file_prefix> <outdir>\n" if (scalar(@ARGV) != 4);
my $TPL = "--user atritt  --passwd halophile ".
          "--genetic_code 11 --gene_caller rast --determine_family --domain Archaea";  
my $names_file = shift;
my $dir = shift;
my $pfx = shift;
my $outdir = shift;

open (NAMES,"<",$names_file);
my %name = ();
while(<NAMES>){
	chomp;
	my @line = split(/\t/,$_);
	$names{$line[0]} = $line[1];
}

print "#!/bin/bash\n";
my @files = glob "$dir/*$pfx";
for my $file (@files){
	my $base = basename($file,$pfx);
	my $name = $base; 
	if (!defined($names{$base})) {
		print STDERR "no name found in $names_file for $base\n";	
	} else {
		$name = $names{$base};
	}
	print "svr_submit_RAST_job --fasta $file --bioname \"$name\" $TPL > $outdir/$base.out 2> $outdir/$base.err\n";
}
# 
#   usage: svr_submit_RAST_job.pl  --user UserName  --passwd Passwd 
#                                 (--genbank GenbankFilename | --fasta FASTAfilename) --domain (Bacteria|Archaea)  
#                                 [--taxon_ID NCBI_taxonomy_ID(def:666666)]   
#                                 [--bioname "genus species strain" (def: "Unknown sp.")]  
#                                 [--genetic_code (11|4)] 
#                                 [--gene_caller (rast|glimmer3)] 
#                                 [--determine_family] 
#                                 [--reannotate_only] [--test] [--nonActive]
#


