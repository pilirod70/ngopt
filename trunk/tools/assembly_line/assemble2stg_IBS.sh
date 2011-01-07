#!/bin/bash

	c="10"
	L="150"
	h="2.2"
	
if [ $# -ne 3 ]; then
	echo "Usage: $0 <output_base> <kcut> <mnpair> <sm_lib> <sm_ins> <lg_lib> <lg_ins>"
	echo "Libraries must be shuffled."
	echo "Final assembly under <output_base>_2stg/scaf_lg/scaffolds_h${h}_L${L}_w<mnpair>.fasta"
	exit -1
fi
	base=$1
	kcut=$2
	w=$3

	a="${base}_2stg"
	if [ ! -d $a ]; then
		mkdir $a
	fi

	stdout="$a/std.out"
	stderr="$a/std.err"
	echo -n > $stdout
	echo -n > $stderr
	echo $base	
	sm_ins=$4
	sm=$5
	lg_ins=$6
	lg=$7
	fa_reads="$base.fasta"
	if [ ! -f $fa_reads ]; then
		cat $sm_ins $lg_ins > $fa_reads 
	fi
	prefix="$a/$base.idba"
	idba -r $fa_reads -o $prefix --minCount $kcut > $stdout 2> $stderr 
	echo "Done contiging\n" >> $stdout 2>> $stderr 
	rmctgs $L $prefix-contig.fa > $prefix-contig.$L.fa

	echo "Begin first stage of scaffolding" | tee -a $stdout $stderr
	scaf_sopra_bwa.sh $a $prefix-contig.$L.fa $sm_ins $w $sm >> $stdout 2>> $stderr 
	echo "Begin second stage of scaffolding" | tee -a $stdout $stderr
	scaf_sopra_bwa.sh $a/scaf_lg $a/scaffolds_h$h\L$L\_w$w.fasta $lg_ins $w $lg >> $stdout 2>> $stderr 
	
