#!/usr/bin/perl
# Author: Andrew Tritt atritt@ucdavis.edu
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

sub matches($$);

my $translate = 0;
my $upstr = 0;
my $downstr = 0;
my $help = 0;
GetOptions( 'trans|t' => \$translate,
			'up|u' => \$upstr,
			'down|d' => \$downstr,
			'help|h' => \$help );

my $usage="extract_annotated_sequence.pl [options] <gbk_file> <term1> <term2> ... <termN> \n".
          "\n".
          "  options:\n".
          "    -t, --trans               \n".
          "    -u, --up=<int>            \n".
          "    -d, --down=<int>          \n".
          "\n".
          "if a search term begins with \'-\', annotations containing that search term will be excluded from the output\n".
		  "use \'extract_annotated_sequence.pl -t\' if CDS are to be translated\n".
		  "sequences are printed to standard output\n";
if (scalar(@ARGV) < 2 || $help){
	print $usage;
	exit;
}

my $gbk_file = shift;
my @srch_str = @ARGV;
unless(defined($gbk_file)){
	print $usage;
	exit;
}

my $in = new Bio::SeqIO( -file => $gbk_file, -format => 'genbank' );
my $out = new Bio::SeqIO( -fh => \*STDOUT , -format => 'fasta' );

while ( my $seq = $in->next_seq() ) {

	#$gbk_hash{'definition'} = $seq->desc();
	#$gbk_hash{'reptype'}    = &get_reptype( $gbk_hash{'definition'} );

	my @features = $seq->get_SeqFeatures();
	my $display = "";
	my $org = "";
	my $seq_name =$seq->display_id(); 
	foreach my $feat (@features) {
		my $ltag = "";
		my @ar;
		if ($feat->primary_tag() eq "source"){
				@ar = $feat->get_tag_values('organism');
				$org = $ar[0];
				$org =~ s/\ /_/g;
			#	print STDERR "Found source tag: >$display<\n";
				next;
		}
		$display = $org;
		if ($feat->has_tag('locus_tag')){
			@ar = $feat->get_tag_values('locus_tag');
			$ltag = $ar[0];
		}
		my $target;
		my $inc = 0;
		if ($feat->has_tag('product')){
			@ar = $feat->get_tag_values('product'); 
			for my $prod (@ar) {
				my $desc = $prod;
				$prod = lc($prod);
				my $num_match = 0;
				for my $query (@srch_str){
					$query = lc($query);
				#	print STDERR "$query\t$prod";
					if (matches($prod,$query)) {
				#		print STDERR " match: $prod \n";
						$num_match++;
					}
				}
				if ($num_match == scalar(@srch_str)) { 
					my $seq_str = ""; 
					$inc = 1;
					$display .= "-".$ltag if (length($ltag)>0);
					$display .= "-".$seq_name;
					if ($feat->strand() == -1) {
						$display .= "|".$feat->end."-".$feat->start;
						my $start = $feat->start-$downstr > 1 ? $feat->start-$downstr : 1;
						my $end = $feat->end+$upstr < $seq->length() ? $feat->end+$upstr : $seq->length();
						$seq_str = $seq->subseq($start,$end);
					} else {
						$display .= "|".$feat->start."-".$feat->end;
						my $start = $feat->start-$upstr > 1 ? $feat->start-$upstr : 1;
						my $end = $feat->end+$downstr < $seq->length() ? $feat->end+$downstr : $seq->length();
						$seq_str = $seq->subseq($start,$end);
					}
					if ($feat->primary_tag eq "CDS" && $translate){
						@ar = $feat->get_tag_values('translation');
						$seq_str = $ar[0];
					} 
					$target = Bio::Seq->new( -seq => $seq_str,
										    # -display_id => $ltag,
										     -display_id => $display, 
										     -desc => $desc );
					if ($feat->strand() == -1 && !$translate){
						$target = $target->revcom();	
					}
					$out->write_seq($target);
					last;
				}
			}
		}
		
	}
}

sub matches($$){
	my $prod = shift;
	my $query = shift;
	if ($query =~ m/^-(.*)/){
		if ($prod =~ /.*$1.*/) {
			return 0;
		} else {
			return 1;
		}
	} else {
		if ($prod =~ /.*$query.*/) {
			return 1;
		} else {
			return 0;
		}
	}
}
