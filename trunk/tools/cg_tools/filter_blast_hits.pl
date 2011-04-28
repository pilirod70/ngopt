#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

Getopt::Long::Configure(qw{no_auto_abbrev no_ignore_case_always});
my $min_aln = 0; 
my $min_pid = 0.0;
my $min_bit = 0;
my $max_eval = 100.0;
my $max_mis = 100000000;
my $max_gap = 100000000;
my $no_reflexive = 0;
my $help = 0;
my $quiet = 0;
my $numargs = scalar(@ARGV);
GetOptions( 'alnlen|a=i' => \$min_aln,
			'percid|p=f' => \$min_pid,
			'bitscore|b=i' => \$min_bit,
			'mismatch|m=i' => \$max_mis,
			'eval|E=f' => \$max_eval,
			'gap|g=i' => \$max_gap, 
			'noreflex' => \$no_reflexive,
			'quiet|q' => \$quiet,
			'help|h' => \$help );



if ($numargs == 0 || $help){
	usage();
	exit 1;
}
summary() unless $quiet;


#my $num_lines = 0;
#my $min_aln = shift;
#my $min_pid = shift;
#my $file = shift;
#open(IN,"<",$file);
#print STDERR "Removing blast hits less than $min_aln bp\n";
#    node0_0_0_93995	node0_0_0_93995	100.00	3815	0	0	3190	7004	3190	7004	0.0	7563
#    node0_0_0_93995	node0_0_0_93995	100.00	2437	0	0	1	2437	1	2437	0.0	4831
#my $num_disc = 0;
#my $num_reflexive = 0;
#my $num_too_short = 0;
#my $num_too_diff = 0;
my %hits = ();
while(<>) {
	chomp;
	my @hit = split;
#	$num_lines++;
	if (okay(\@hit)) {
		if ($no_reflexive){
			if ($hit[0] eq $hit[1]) {
				if ($hit[6] != $hit[8] || $hit[7] != $hit[9]) {
					add_hit(\@hit);
				} 
			} 	
		} else {
			add_hit(\@hit);
		}	

	}

#	if ($hit[2] >= $min_pid && $hit[3] >= $min_aln) {
#		if ($hit[0] eq $hit[1]) {
#			if ($hit[6] != $hit[8] || $hit[7] != $hit[9]) {
#			#	print $_."\n";
#				add_hit(\@hit);
#			} else {
#				$num_reflexive++;
#			}	
#		} else {
#			if ($hit[2] >= $min_pid) {
#			#	print $_."\n";
#				add_hit(\@hit);
#			} else {
#				$num_too_diff++;
#			}
#		}
#	} else { 
#		$num_too_short++;
#	}
}

foreach my $ctg1 (keys %hits) {
	foreach my $ctg2 (keys %{$hits{$ctg1}}){
		foreach my $hit (@{$hits{$ctg1}{$ctg2}}){
			print STDOUT $ctg1."\t".$ctg2."\t".join("\t",@$hit)."\n";
		}
	}
} 

#print STDERR "Removed:\n";
#print STDERR "\t$num_reflexive reflexive hits\n";
#print STDERR "\t$num_too_short hits with less than $min_aln aligned base pairs\n";
#print STDERR "\t$num_too_diff hits with less than $min_pid percent identity\n";

sub add_hit {
	my ($hit) = shift;
#	print STDERR "size = ".scalar(@$hit)." at line $num_lines\n";
	my $ctg1 = shift @$hit;
	my $ctg2 = shift @$hit;
	if ($ctg2 lt $ctg1) {
		my $tmp = $ctg2;
		$ctg2 = $ctg1;
		$ctg1 = $ctg2;
		# swap 4,5 with 6,7
		$tmp = $hit->[4];
		$hit->[4] = $hit->[6];
		$hit->[6] = $tmp;
		$tmp = $hit->[5];
		$hit->[5] = $hit->[7];
		$hit->[7] = $tmp;
	}
	if (defined($hits{$ctg1})){
		if (defined($hits{$ctg1}{$ctg2})){
#			print STDERR "$ctg1 and $ctg2 have ".scalar(@{$hits{$ctg1}{$ctg2}})." matches so far\n";	
			my $dup = 0;
			foreach my $tmp_hit (@{$hits{$ctg1}{$ctg2}}) {
#				print STDERR "1:length of hit: ".scalar(@$tmp_hit)."\n"; 
#				print STDERR "2:length of hit: ".scalar($tmp_hit)."\n"; 
#				print STDERR "hit = $tmp_hit\n"; 
				if ($hit->[4] == $tmp_hit->[4] && $hit->[6] == $tmp_hit->[6]){
					$dup = 1;
					last;
				}
			}
			if (!$dup){
				push (@{$hits{$ctg1}{$ctg2}}, $hit);
			}
		} else {
			$hits{$ctg1}{$ctg2} = ();
			push (@{$hits{$ctg1}{$ctg2}}, $hit);
		}
	} else {
		$hits{$ctg1} = ();
		$hits{$ctg1}{$ctg2} = ();
		push (@{$hits{$ctg1}{$ctg2}}, $hit);
	}
}

sub okay {
	my ($hit) = shift;
	if ($hit->[2] >= $min_pid && $hit->[3] >= $min_aln && $hit->[4] <= $max_mis &&
		$hit->[5] <= $max_gap && $hit->[10] <= $max_eval && $hit->[11] >= $min_bit) {
		return 1;
	} else {
		return 0;
	}
}


sub usage {
	print STDOUT "filter_blast_hits.pl [options] <m8_output>\n\n".
				 "  options:\n".
				 "   -a, --alnlen=<int>      minimum alignment length [0]\n".
				 "   -p, --percid=<float>    minimum percent identity [0.0]\n".
				 "   -b, --bitscore=<int>    minimum bitscore [0]\n".
				 "   -m, --mismatch=<int>    maximum number of mismatches [100000000]\n".
				 "   -E, --eval=<float>      maximum E-value [100.0]\n".
				 "   -g, --gap=<int>         maximum number of gap openings [100000000]\n".
				 "   --noreflex              discard reflexive hits\n".
				 "   -h, --help              display this message\n\n".
				 "   -q, --quiet             do not print status and summary messages\n".
				 "if <m8_output> is absent, input will be read from stdin\n\n";

}


sub summary {
	print "minimum alignment length: $min_aln\n";
	print "minimum percent identity: $min_pid\n";
	print "minimum bit score:        $min_bit\n";
	print "maximum E-val:            $max_eval\n";
	print "maximum no. mismatches:   $max_mis\n";
	print "maximum no. gap openings: $max_gap \n";
	print "keep reflexive hits?      ".($no_reflexive ? "no" : "yes")."\n";
}
