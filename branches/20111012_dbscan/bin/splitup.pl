#!/usr/bin/perl -w
use strict;
my $reads_per_file = 500000;
my $input_file = $ARGV[0];
my $output_dir = $ARGV[1];
open( READIN, "$input_file" );
my $lineI=0;
my $fileI=1;
open( READOUT, ">$output_dir/reads.$fileI" );
while( my $line = <READIN> )
{
	if( $lineI>0 && $lineI % ($reads_per_file*2)==0 )
	{
		$fileI++;
		close(READOUT);
		open( READOUT, ">$output_dir/reads.$fileI" );
	}
	$lineI++;
	print READOUT $line;
}
close(READOUT);
