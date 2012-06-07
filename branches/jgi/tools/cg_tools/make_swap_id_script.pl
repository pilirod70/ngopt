#!/usr/bin/perl -w
use strict;
use warnings;
if (scalar(@ARGV) != 1) {
	print "\nUsage: make_swap_id_script.pl <identifier_file>\n\n".
	      "Halophile_identifiers.csv can be used to make or used as <identifier_file>\n".
	      "<identifier_file> must have the following columns in this order:\n".
	      "\n\thuman_readable_name    Genus    species    identifier\n\n";
	exit;
}

open (IN,"<",$ARGV[0]);
<IN>;
print "#!/bin/bash\n";

print "if [ \$# -ne 1 ]; then\n";
print "\techo \"Usage: \`basename \$0\` <file>\"\n";
print "\texit;\n";
print "fi\n";


while (<IN>) {
	my @line = split /\t/;
	$line[2] =~ s/\ /_/g;
	print "perl -p -i -e \'s/".$line[4]."/".$line[1]."_".$line[2]."/g\' \$1\n";
}
