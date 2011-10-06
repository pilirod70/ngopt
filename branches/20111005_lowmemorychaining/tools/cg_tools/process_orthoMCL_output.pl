#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;

if (scalar(@ARGV) != 3) {
	print "Usage: ".basename($0)." <org_id_file> <orthoMCL_output> <output_prefix>\n".
	      "<org_id_file> should have one organism id per line\n".
	      "Five output files will be created:\n".
	      " <output_prefix>.bin.txt \n".
	      " <output_prefix>.count.txt \n".
	      " <output_prefix>.perc-var.txt\n".
	      " <output_prefix>.pan_count.txt \n".
	      " <output_prefix>.core_count.txt \n";
	exit 1;
}

my $id_file = shift;
my $ortho_file = shift;
my $base = shift;
open (ID,"<",$id_file);
# orthos: a hash of hashes of arrays  of Organism x OrthoGroupID x GroupMembers: 
# the arrays are necessary for dealing with paralogs

my %orthos = ();
while (<ID>){
	chomp;
	$orthos{$_} = ();
}

open (ORTHO,"<",$ortho_file);
my @ortho_names;
my %ortho_groups = ();
while (<ORTHO>){
	chomp;
	my @line = split;
	my $ortho_id = shift @line;
	push(@ortho_names,$ortho_id);
	$ortho_groups{$ortho_id} = ();
	foreach my $mem (@line) {
		(my $org, my $ltag) = split(/\|/,$mem);
		push (@{$orthos{$org}{$ortho_id}},$ltag);
		if (defined($ortho_groups{$ortho_id}{$org})) {
			$ortho_groups{$ortho_id}{$org}++; 
		} else {
			$ortho_groups{$ortho_id}{$org} = 1;
		}
	}
}
my %dist = ();
my %pan_size = ();
my @orgs;
foreach (keys %orthos){
	$dist{$_} = ();
	$pan_size{$_} = ();
	push (@orgs, $_);
}
my $num_taxa = scalar(keys %orthos);
for (my $i = 0; $i < $num_taxa; $i++){
	my $org_id = shift @orgs;
	$pan_size{$org_id}{$org_id} = scalar(keys %{$orthos{$org_id}});
	$dist{$org_id}{$org_id} = 0.0;
#	print $org_id."\t".scalar(keys $orthos{$org_id})
	foreach (@orgs) {
		my ($core,$var) = calc_cv($orthos{$org_id},$orthos{$_});
#		if ($core == 0 || $var == 0){
			my $s1 = keys %{$orthos{$org_id}};
			my $s2 = keys %{$orthos{$_}};
			print STDERR "$org_id = $s1 $_ = $s2 : core = $core var = $var\n";
#		}
		$pan_size{$_}{$org_id} = $var+$core;
		$pan_size{$org_id}{$_} = $pan_size{$_}{$org_id};
		$dist{$org_id}{$_} = $var/$pan_size{$org_id}{$_};
		$dist{$_}{$org_id} = $dist{$org_id}{$_}
	}
}
open (PERC,">",$base.".perc-var.txt");
open (PAN,">",$base.".pan_count.txt"); 
open (CORE,">",$base.".core_count.txt"); 

my $hdr = join("\t",(keys %orthos));
print PERC "\t$hdr\n";
print PAN "\t$hdr\n";
print CORE "\t$hdr\n";

foreach my $org1 (keys %orthos) {
	print PERC $org1;
	print PAN $org1;
	print CORE $org1;
	foreach my $org2 (keys %orthos) {
		print PERC "\t".$dist{$org1}{$org2}; 
		print PAN "\t".$pan_size{$org1}{$org2}; 
		print CORE "\t".($pan_size{$org1}{$org2} * (1-$dist{$org1}{$org2})); 
	}
	print PERC "\n";
	print PAN "\n";
	print CORE "\n";
}

close PERC;
close PAN;
close CORE;

open (BIN,">",$base.".bin.txt");
open (COUNT,">",$base.".count.txt");
my @taxa = keys %orthos;
$hdr = join("\t",@taxa);
print BIN "OrthGrp\\Taxa\t$hdr\n";
print COUNT "OrthGrp\\Taxa\t$hdr\n";
foreach my $ortho (@ortho_names) {
	print BIN $ortho;
	print COUNT $ortho;
	for my $org (@taxa) {
		if (defined($ortho_groups{$ortho}{$org})){
			print BIN "\t1";
			my $count = $ortho_groups{$ortho}{$org}; 
			print COUNT "\t$count";
		} else {
			print BIN "\t0";
			print COUNT "\t0";
		}
	}
	print BIN "\n";
	print COUNT "\n";
}

close BIN;
close COUNT;

sub calc_cv {
    my ($list1,$list2) = @_;
    my $core=0.0;
	my $var = scalar(keys %$list1) + scalar(keys %$list2);
    foreach my $id1(keys %$list1){
		if(exists($list2->{$id1})){
		    $core = $core + 1.0;
			$var = $var - 2.0;
		}
    }
    return ($core,$var);
}

