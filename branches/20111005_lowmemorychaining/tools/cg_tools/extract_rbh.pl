#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 1){
	print "Usage ".basename($0)." <tabd_blast_file> \n";
	exit 1;
}
my $num_lines = 0;
my $file = shift;
open(IN,"<",$file);
my %hits = ();
my $count = 0;
while(<IN>) {
	chomp;
	my @hit = split;
	my ($orgQ, $Q_id) = split(/\|/, $hit[0]);
	my ($orgS, $S_id) = split(/\|/, $hit[1]);
	$hits{$orgQ} = () if (!defined($hits{$orgQ}));
	$hits{$orgQ}{$Q_id} = () if (!defined($hits{$orgQ}{$Q_id}));

#	100.00	3815	0	0	3190	7004	3190	7004	0.0	7563
#	100.00	2437	0	0	1	2437	1	2437	0.0	4831
	if (!defined($hits{$orgQ}{$Q_id}{$orgS})){
#		$hits{$orgQ}{$Q_id}{$orgS} = ();
		$hits{$orgQ}{$Q_id}{$orgS} = \@hit; 
		$count++;
	}else{
		if ($hit[11] == $hits{$orgQ}{$Q_id}{$orgS}[11]){
			if ($hit[10] == $hits{$orgQ}{$Q_id}{$orgS}[10]){
				if ($hit[2] == $hits{$orgQ}{$Q_id}{$orgS}[2]){
					if ($hit[3] > $hits{$orgQ}{$Q_id}{$orgS}[3]){   # alignment length
						$hits{$orgQ}{$Q_id}{$orgS} =  \@hit; 
					}
				} elsif ($hit[2] >  $hits{$orgQ}{$Q_id}{$orgS}[2]){ # pid 
					$hits{$orgQ}{$Q_id}{$orgS} =  \@hit; 
				}
			} elsif ($hit[10] > $hits{$orgQ}{$Q_id}{$orgS}[10]){    # e-val
				$hits{$orgQ}{$Q_id}{$orgS} =  \@hit; 
			}
		} elsif ($hit[11] < $hits{$orgQ}{$Q_id}{$orgS}[11]){        # bit-score
			$hits{$orgQ}{$Q_id}{$orgS} =  \@hit; 
		}
		$count++;
	}
}

print "Finished reading BLAST file. Processed $count hits.\n";

my @bidis = ();

my @taxa = keys %hits;
my $n_taxa = scalar(@taxa);
for (my $i = 0; $i < $n_taxa; $i++) {
	my $org1 = shift (@taxa); 
	for my $org2 (@taxa) {
		for my $id1 (keys %{$hits{$org1}}) {
			next unless(defined($hits{$org1}{$id1}{$org2}));

			my $hit = $hits{$org1}{$id1}{$org2}; 
			my $id2 = substr($hit->[1],index($hit->[1],"|")+1);
			next unless (defined($hits{$org2}{$id2}{$org1}));

			my $test_id = $hits{$org2}{$id2}{$org1}->[1];
			$test_id = substr($test_id,index($test_id,"|")+1);
			if ($test_id eq $id1) { # found bidi hit!
				push(@bidis,$hit);
				delete ($hits{$org2}{$id2}{$org1});
			}
		}
	}
}
open (FILT,">",$file.".bidi_filt");
open (PAIR,">",$file.".pairwise");
for my $aref (@bidis) {
	print FILT join("\t",@$aref)."\n";
	print PAIR join("\t",@$aref[0..1])."\n";
}

close FILT;
close PAIR;


#Use of uninitialized value $hit[1] in index at ./extract_best_bidi_hits.pl line 44, <IN> line 21851197.
#Use of uninitialized value in substr at ./extract_best_bidi_hits.pl line 44, <IN> line 21851197.
#substr outside of string at ./extract_best_bidi_hits.pl line 44, <IN> line 21851197.
#Use of uninitialized value $id2 in hash element at ./extract_best_bidi_hits.pl line 45, <IN> line 21851197.
