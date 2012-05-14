#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
sub mu (\@\@);
sub var ($\@\@);
sub covar($\@$\@\@);
if (scalar(@ARGV) == 0 || $ARGV[0] eq "-h" || $ARGV[0] eq "--help") { 
	print "Usage: ".basename($0)." <tabd_file(s)>\n"; 
	exit;
}

print "org\tcov_mean\tcov_var\tgc_mean\tgc_var\tlen_mean\tlen_var\tcovar_cov_gc\tcovar_gc_len\tcovar_len_cov\tn_ctgs\n";
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
	my @dummy = (1) x length(@len);
	my $cov_mu = mu(@cov,@len);
	my $gc_mu = mu(@gc,@len);
	my $len_mu = mu(@len,@dummy);
	my $cov_var = var($cov_mu,@cov,@len);
	my $gc_var = var($gc_mu,@gc,@len);
	my $len_var = var($len_mu,@len,@dummy);
	my $covar_cov_gc = covar($cov_mu,@cov,$gc_mu,@gc,@len);
	my $covar_gc_len = covar($gc_mu,@gc,$len_mu,@len,@len);
	my $covar_len_cov = covar($len_mu,@len,$cov_mu,@cov,@len);
	my $base = fileparse($file,("\.pre-scaf\.gc-cov-len\.txt","\.scaf\.gc-cov-len\.txt"));
	print $base;
	foreach my $stat ($cov_mu,$cov_var,$gc_mu,$gc_var,$len_mu,$len_var,$covar_cov_gc,$covar_gc_len,$covar_len_cov) {
		printf "\t%.3f",($stat);
	}
	
	printf "\t%.0d\n", (scalar(@cov));
}	
sub mu (\@\@){
	my $vec = shift;
	my $w = shift;
	my $num = 0;
	my $den = 0;
	for (my $i = 0; $i < length(@$vec); $i++) {
		$num += @$vec[$i]*@$w[$i];
		$den += @$w[$i];
	}
	return $num/$den;
}

sub var ($\@\@) {
	my $mean = shift;
	my $vec = shift;
	my $w = shift;
	my $num = 0;
	my $den = 0;
	for (my $i = 0; $i < length(@$vec); $i++) {
		$num += (($mean - @$vec[$i])**2)*@$w[$i];
		$den += @$w[$i];
	}
	return $num/$den;
}

sub covar($\@$\@\@) {
	my $mu_x = shift;
	my $vec_x = shift;
	my $mu_y = shift;
	my $vec_y = shift;
	my $w = shift;
	my $num = 0;
	my $den = 0;
	for (my $i = 0; $i < scalar(@$vec_x); $i++) {
		$num += ((@$vec_x[$i] - $mu_x)*(@$vec_y[$i] - $mu_y))*@$w[$i];
		$den += @$w[$i];
	}
	return $num/$den;
}

