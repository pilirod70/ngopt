#!/usr/bin/perl -w
use File::Basename;
if (scalar(@ARGV)==0){ print "Usage: $0 file1 file2 .... fileN\n"; }
while(@ARGV) {
	my $file = shift;
	my $base = basename($file,".err");
	open (IN,"<",$file);
	<IN>;
	my $line = <IN>;
	chomp $line;
	<IN>;
	my $id = (split("\'",$line))[1];
	print $base."\t".$id."\n";

}

#
#Job '20876' was successfully started
#
