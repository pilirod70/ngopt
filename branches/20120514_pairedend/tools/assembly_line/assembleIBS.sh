#!/bin/bash

	c="10"
	L="150"
	h="2.2"
	
if [ $# -ne 5 ]; then
	echo "Usage: `basename $0` <output_base> <kcut> <mnpair> <lib> <ins>"
	echo "<lib> must be shuffled and in fasta format."
	exit -1
fi
	base=$1
	kcut=$2
	w=$3

	a=$base
	if [ ! -d $a ]; then
		mkdir $a
	fi

	echo "Outputing to $base/"
	stdout="$a/std.out"
	stderr="$a/std.err"
	lib=$4
	ins=$5
	prefix="$a/$base.idba"
	echo "Begin contiging" > $stdout 2> $stderr 
	idba -r $lib -o $prefix --minCount $kcut >> $stdout 2>> $stderr 
	echo "Done contiging" >> $stdout 2>> $stderr 
	echo "Discarding contigs smaller than $L bp"
	rmctgs $L $prefix-contig.fa > $prefix-contig.$L.fa
	echo "Begin scaffolding" | tee -a $stdout $stderr
	scaf_sopra_bwa.sh $a $prefix-contig.$L.fa $lib $w $ins >> $stdout 2>> $stderr 
