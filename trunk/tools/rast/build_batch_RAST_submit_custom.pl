#!/usr/bin/perl -w

use File::Basename;
die "Usage: ".basename($0)." <name_file> <fasta_dir> <fasta_file_suffix> <outdir> <keep_file>\n".
		"<name_file> should contain two columns:  CC_ID   HUMAN_SENSIBLE_NAME\n" if (scalar(@ARGV) != 5);
my $TPL = "--user darling  --passwd C94zvhvD ".
          "--genetic_code 11 --gene_caller rast --determine_family --domain Archaea";  
my $names_file = shift;
my $dir = shift;
my $sfx = shift;
my $outdir = shift;
my $spec_file = shift;

open (NAMES,"<",$names_file);
my %name = ();
while(<NAMES>){
	chomp;
	my @line = split(/\t/,$_);
	$names{$line[0]} = $line[1];
}
#my %keep = ();
#open(NAMES,"<",$spec_file);
#while(<NAMES>){
#	chomp;
#	$keep{$_} = 1;
#}
print "#!/bin/bash\n";
my @files = glob "$dir/*$sfx";
for my $file (@files){
	my $base = basename($file,$sfx);
#	next if (!defined($keep{$base}));
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


