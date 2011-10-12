#!/usr/bin/perl -w
my $mapData = $ARGV[0];
my $chrLen = $ARGV[1];
my $thinprob = $ARGV[2];
my $readlen = $ARGV[3];
my $outter = $ARGV[4];

die "Usage: <SHRiMP map> <chr len> <thin probability> <read length> <output>" if( @ARGV != 5 );

my @readcount;	# array to store
for( my $i = 0; $i < $chrLen; $i++)
{
	$readcount[$i] = 0;
}

open( MAPDATA, "$mapData" ) || die "no taspd\n";
my $cur_read = "";
my @cur_lines;
my $lI = 0;
while( my $line = <MAPDATA> )
{
	print "Read $lI lines\n" if( $lI++ % 10000 == 0 );
	my @l = split( /\t/, $line );
	next if $l[0] =~ /^\#/;
	if($l[0] ne $cur_read){
		my @rs = @cur_lines;
		@cur_lines = ();
		push( @cur_lines, $line );
		$cur_read = $l[0];
		next if rand(1) < $thinprob;	# skip this read
		next if @rs == 0; # first read set will be empty

		# finished a read set
		# pick one read at random
		my $r = $rs[int(rand(@rs))];
		my @ll = split( /\t/, $r );
		my $start = $ll[2] eq "+" ? $ll[3] : $ll[4];
		my $end = $ll[2] eq "+" ? $ll[4] : $ll[3];
		my $dir = $ll[2] eq "+" ? 1 : -1;
		my $rlen = $readlen > 0 ? $readlen : abs($start-$end);
		# count each site hit by the read	
		for( my $i = 0; $i < $rlen; $i++){
			$readcount[ ($start+($dir*$i))%$chrLen ]++;
		}
	}
}

open( OUTTIE, ">$outter" ) || die "buabuab\n";
foreach my $k(@readcount){
	print OUTTIE "$k\n";
}
