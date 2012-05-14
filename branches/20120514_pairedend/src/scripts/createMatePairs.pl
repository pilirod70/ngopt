#!/usr/bin/perl -w
# small program to create simulated mate paired reads from
# mapped illumina data
use strict;
use Math::Trig;

my $mapData = $ARGV[0];
my $chrLen = $ARGV[1];
my $pairCount = $ARGV[2];
my $insert_size = $ARGV[3];
my $insert_sd = $ARGV[4];

die "Usage: <SHRiMP map> <chr len> <mate pair count> <insert mean> <insert sd>" if( @ARGV != 5 );


#
# Step 1: create the read map in memory
#

my @readsPlus; # array storing reads mapped to each site
my @readsMinus; # array storing reads mapped to each site
my %readsToSitePlus;
my %readsToSiteMinus;

open( MAPDATA, "$mapData" ) || die "no taspd\n";
my $cur_read = "";
my @cur_lines;
my $lI = 0;
my $plusReadCount = 0;
my $minusReadCount = 0;
while( my $line = <MAPDATA> )
{
	print STDERR "Read $lI lines\n" if( $lI++ % 10000 == 0 );
	my @l = split( /\t/, $line );
	next if $l[0] =~ /^\#/;
	if($l[0] ne $cur_read){
		my @rs = @cur_lines;
		@cur_lines = ();
		push( @cur_lines, $line );
		$cur_read = $l[0];
		next if @rs == 0; # first read set will be empty

		# finished a read set
		# pick one read at random
		my $r = $rs[int(rand(@rs))];
		my @ll = split( /\t/, $r );
		my $start = $ll[2] eq "+" ? $ll[3] : $ll[4];
		if($ll[2] eq "+")
		{
			if(!defined($readsPlus[$start])){
				$readsPlus[$start]= [ ($ll[0]) ];
				$plusReadCount++;
			}else{
				push( @{$readsPlus[$start]}, $ll[0] );
			}
			$readsToSitePlus{$ll[0]}=$start;
		}else{
			if(!defined($readsMinus[$start])){
				$readsMinus[$start]= [ ($ll[0]) ];
				$minusReadCount++;
			}else{
				push( @{$readsMinus[$start]}, $ll[0] );
			}
			$readsToSiteMinus{$ll[0]}=$start;
		}
	}else{
		push( @cur_lines, $line );
	}
}

close MAPDATA;
print STDERR "found $plusReadCount forward and $minusReadCount reverse reads\n";

#
# Step 2:  pick pairs
#
# sample from a truncated normal centered around insert size
#
print STDERR "Mating...\n";
my @read_pairs;
my @sample_reads = keys(%readsToSitePlus);
my $sample_count = scalar(@sample_reads);
my @probs;
for(my $pairI=0; $pairI<$pairCount; $pairI++)
{
	# make sure we haven't run out of reads
	last if $sample_count == 0;
	# pick a read uniform from all reads, plus and minus
	my $rI = int(rand($sample_count));
	my $read;
	my $read_l;
	{
		$read = $sample_reads[$rI];
		$read_l = $readsToSitePlus{$read};
		my $ins_min = $read_l + $insert_size - 3*$insert_sd;
		my $ins_max = $read_l + $insert_size + 3*$insert_sd;
		# sum up count of all reads in the region
		my $sum=0;
		for( my $i=$ins_min; $i<$ins_max; $i++)
		{
			$sum += @{$readsMinus[$i%$chrLen]} if defined($readsMinus[$i%$chrLen]);
			$sum += 0 if !defined($readsMinus[$i%$chrLen]);
		}
		# weight reads by the normal probability distribution
		# this creates a mixture of the coverage distribution and the normal distribution
		my $j=0;
		my $sum2=0;
		my $ts = $insert_size - 3*$insert_sd;
		for( my $ii=$ins_min; $ii<$ins_max; $ii++)
		{
			$probs[$j] = scalar(@{$readsMinus[$ii%$chrLen]}) if defined($readsMinus[$ii%$chrLen]);
			$probs[$j] = 0 if !defined($readsMinus[$ii%$chrLen]);
			$probs[$j] *= exp(-($j+$ts-$insert_size)**2/(2*$insert_sd**2))/($insert_sd*sqrt(2*pi));
			$sum2 += $probs[$j];
			$probs[$j] = $sum2;	# set this to the cumulative sum, so we can use bin search to sample later
			$j++;
		}
		# now pick a site from the truncated normal distribution
		# choose a uniform from 0 to $sum2
		my $randy = rand($sum2);
		# find the nearest index in $probs
		my $nearest = int(BinSearch ($randy, \&cmpFunc, \@probs)+0.5);
		my $p_index = ($ins_min+$nearest)%$chrLen;
		my @bindat = @{$readsMinus[$p_index]};
		# found a site, now pick a read within that site
		my $mr = int(rand(@bindat));
		my $mate_read = $bindat[$mr];
		push( @read_pairs, [($read, $mate_read)] );
		# remove the forward and reverse reads from consideration (sampling without replacement)
		$readsMinus[$p_index][$mr] = $bindat[$#bindat];
		pop(@{$readsMinus[$p_index]});
#		delete $readsToSitePlus{$read};
#		delete $readsToSiteMinus{$mate_read};
		$sample_reads[$rI] = $sample_reads[$sample_count-1];	# achieves sampling without replacement
		$sample_count--;
		print "$read\t$mate_read\n";
	}
}

sub cmpFunc {
    $_[0] <=> $_[1];
}



sub BinSearch {
    my ($target, $cmp, $array) = @_;
    my $posmin = 0;
    my $posmax = $#$array;
    
    return undef if ! $array;
    return -0.5 if $cmp->($array->[0], $target) > 0;
    return $#$array + 0.5 if $cmp->($array->[-1], $target) < 0;
    
    while (1) {
        my $mid = int (($posmin + $posmax) / 2);
        my $result = $cmp->($array->[$mid], $target);
        
        if ($result < 0) {
            $posmin = $posmax, next if $mid == $posmin && $posmax != $posmin;
            return $mid + 0.5 if $mid == $posmin;
            $posmin = $mid;
        } elsif ($result > 0) {
            $posmax = $posmin, next if $mid == $posmax && $posmax != $posmin;
            return $mid - 0.5 if $mid == $posmax;
            $posmax = $mid;
        } else {
            return $mid;
        }
    }
}


