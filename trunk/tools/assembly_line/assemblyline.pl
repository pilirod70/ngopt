#!/usr/bin/perl -w 

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
sub usage();

unless (@ARGV) {
	usage();
	exit 0;
}
my $bin="~/bin";
#my $data_file=$ARGV[0];
#my $bc_file=$ARGV[1];

my $base = $ARGV[0];
my $qual_cut=0;
my $mink=25;
my $maxk=25;
my $kfreq_cut=1;
my $ins=0;
my $min_pairs=4;
my $max_maphits=5;
my $min_ctglen=100;
my $outdir="asmline_$base";
my $result = GetOptions (
				"minqual=i" => \$qual_cut,
				"mink=i"  => \$mink,
				"maxk=i"  => \$maxk,
				"minkfreq=i" => \$kfreq_cut,
				"insert=i" => \$ins,
				"minpair=i" => \$min_pairs,
				"maxhits=i" => \$max_maphits, 
				"minlen=i" => \$min_ctglen,
#				"help" => \$help,
#				"doc" => \$doc,
				"output=s" => \$outdir
                  ) ; 
my $abs_base_path = abs_path($base."_p1.fastq");
my $indir = dirname($abs_base_path);
#if ($indir eq "./") $indir = cwd();
print STDOUT "input dir = ", $indir,"\n"; 
print STDOUT "outdir = ".$outdir."\n";
if (!(-d $outdir)){
	mkdir $outdir;
}

my $outbase = $outdir."/".$base;
my $ctg_in = $outbase.".fastq";
my @fq_files = glob($indir."/".$base."_*p*.fastq");
my $cmd = "cat @fq_files";
`$cmd > $ctg_in`;
#open my $pipe, "$cmd |" or die "Couldn't concatenate reads\n";
#open my $fh, '>',$ctg_in or die "Couldn't open $ctg_in for concatenating into\n";
#print $fh while <$pipe>; 
#close $pipe or die "Couldn't write concatenation to $ctg_in\n";
#print STDOUT "Concatenated all files with base $base to $ctg_in for assembling contigs\n";
#exit 1;
my $stdout = system(split(' ',$cmd)) == 0 or die "system $cmd failed";

if ($qual_cut>0) {
	print STDOUT "Quality trimming data...\n";
	print STDERR "Quality trimming data...\n";
	$cmd = "$bin/trimfq $ctg_in $maxk $qual_cut";
	system($cmd);
	$ctg_in = $outbase.".htqs.fastq";
}
$cmd = "$bin/repair --shuf -s .htqs.fastq $outbase $ctg_in";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
$cmd = "$bin/fq2fa $outbase.htqs.fastq $outbase.htqs.fasta";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
$ctg_in = $outbase."_shuf.htqs.fasta";
$cmd = "$bin/idba --read $ctg_in -o $outbase.htqs --scaffold --mink $mink --maxk $maxk --minCount $kfreq_cut --minPairs $min_pairs --minLength $min_ctglen";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
my $ctgs_out="$outbase.htqs-contigs.fa";
my $scafout = $outdir."/scaffold";
if (!(-d $scafout)){
	mkdir $scafout;
}
$cmd = "$bin/prep_contigAseq_v1_x1.pl -contig $ctgs_out -mate $ctg_in -a $scafout";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
my $prep_ctgs = "$scafout/contigs_sopra.fasta";
$cmd = "$bin/repair --split -s .htqs.fasta $outbase $ctg_in";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
my $prefix="$scafout/$base.bwa";
my $fq1 = "$outbase\_p1.htqs.fastq";
my $fq2 = "$outbase\_p2.htqs.fastq";
my $sai1 = "$scafout/$base\_p1.sai";
my $sai2 = "$scafout/$base\_p2.sai";
$cmd = "$bin/bwa index -a is -p $prefix $prep_ctgs";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
$cmd = "$bin/bwa aln -f $sai1 $prefix $fq1";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
$cmd = "$bin/bwa aln -f $sai2 $prefix $fq2";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
my $max_ins = $ins+100;
my $sam = "$scafout/$base.sam";
$cmd = "$bin/bwa sampe -f $sam $prefix $sai1 $sai2 $fq1 $fq2";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";

$cmd = "$bin/parse_sam_v1_x2.pl -sam $sam -a $scafout";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
$cmd = "$bin/read_parsed_sam_v1_x3.pl -parsed $sam.sam_parsed -d $ins -a $scafout";
#system($cmd); 
system(split(' ',$cmd)) == 0 or die "system $cmd failed";
$cmd = "$bin/scaf_v1_x4.pl -o $scafout/orientdistinfo_c5 -a $scafout -w $min_pairs -L $min_ctglen";
#system($cmd);
system(split(' ',$cmd)) == 0 or die "system $cmd failed";





sub usage(){
	print "Usage : $0 <base> [options] \n";
	print "  options:\n";
	print "  --minqual\n";
	print "  --mink\n";
	print "  --maxk\n";
	print "  --minkfreq\n";
	print "  --insert\n";
	print "  --minpair\n";
	print "  --maxhits\n";
	print "  --minlen\n";
#	print "  --help\n";
#	print "  --doc\n";
	print "  --output\n";
}
