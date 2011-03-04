#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
sub mu (\@);
sub var ($\@);
sub cor($\@$\@);
if (scalar(@ARGV) == 0 || $ARGV[0] eq "-h" || $ARGV[0] eq "--help") { 
	print "Usage: ".basename($0)." <tabd_file(s)>\n"; 
	exit;
}

print "org\tcov_mean\tcov_var\tgc_mean\tgc_var\tlen_mean\tlen_var\tcor_cov_gc\tcor_gc_len\tcor_len_cov\tn_ctgs\n";
while (@ARGV) {
	my $file = shift;
	open(IN,"<",$file);
	my $line = <IN>;
	chomp $line;
	my @hdr = split(/\t/, $line);
	next if (scalar(@hdr) != 4 || $hdr[1] ne "cov" || $hdr[2] ne "gc" || $hdr[3] ne "len");
#	if (scalar(@hdr) != 4 || $hdr[1] ne 'cov' || $hdr[2] ne 'gc' || $hdr[3] ne 'len'){
#		print "bad header: >".join("<\t>",@hdr)."<\n";
#		next;
#	}
	my (@gc, @cov, @len);
	while(<IN>){
		chomp;
		my @ctg = split /\t/;
		push (@cov,$ctg[1]);
		push (@gc,$ctg[2]);
		push (@len,$ctg[3]);
	}
	my $cov_mu = mu(@cov);
	my $gc_mu = mu(@gc);
	my $len_mu = mu(@len);
	my $cov_var = var($cov_mu,@cov);
	my $gc_var = var($gc_mu,@gc);
	my $len_var = var($len_mu,@len);
	my $cor_cov_gc = cor($cov_mu,@cov,$gc_mu,@gc);
	my $cor_gc_len = cor($gc_mu,@gc,$len_mu,@len);
	my $cor_len_cov = cor($len_mu,@len,$cov_mu,@cov);
	my $base = fileparse($file,("\.pre-scaf\.gc-cov-len\.txt","\.scaf\.gc-cov-len\.txt"));
	print $base;
	foreach my $stat ($cov_mu,$cov_var,$gc_mu,$gc_var,$len_mu,$len_var,$cor_cov_gc,$cor_gc_len,$cor_len_cov) {
		printf "\t%.3f",($stat);
	}
	
	printf "\t%.0d\n", (scalar(@cov));
}	
sub mu (\@){
	my $vec = shift;
	my$sum = 0;
	for my $i (@$vec) {
		$sum += $i;
	}
	return $sum/scalar(@$vec);
}

sub var ($\@) {
	my $mean = shift;
	my $vec = shift;
	my $sum = 0;
	for my $i (@$vec) {
		$sum += ($mean-$i)**2;
	}
	return $sum/scalar(@$vec);
}

sub cor($\@$\@) {
	my $mu_x = shift;
	my $vec_x = shift;
	my $mu_y = shift;
	my $vec_y = shift;
	my $sum = 0;
	for (my $i = 0; $i < scalar(@$vec_x); $i++) {
		$sum += (@$vec_x[$i] - $mu_x)*(@$vec_y[$i] - $mu_y);
	}
	return $sum/scalar(@$vec_x);
}

