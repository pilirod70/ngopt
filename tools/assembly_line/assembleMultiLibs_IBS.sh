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
	num_libs=$(( $# - 3 ))
	
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
		x=4
		while [ $x -le ($# - 1) ] ; do 
			cat $(( $x )) >> $fa_reads 
			x=$(( $x + 2 ))
		done
	fi
	prefix="$a/$base.idba"
	idba -r $fa_reads -o $prefix --minCount $kcut > $stdout 2> $stderr 
	echo "Done contiging\n" >> $stdout 2>> $stderr 
	rmctgs $L $prefix-contig.fa > $prefix-contig.$L.fa

	contigs=$prefix-contig.$L.fa
	x=4
	stage=1
	dir=$a/stage$stage
	while [ $x -le ($# - 1) ] ; do 
		lib=$(( $x )) >> $fa_reads 
		ins=$(( $x + 1 ))
		echo "Begin stage $stage of scaffolding" | tee -a $stdout $stderr
		scaf_sopra_bwa.sh $dir $contigs $lib $w $ins >> $stdout 2>> $stderr 
		x=$(( $x + 2 ))
		stage=$(( $stage + 1 ))
		contigs=$dir/scaffolds_h$h\L$L\_w$w.fasta
		dir=$a/stage$stage
	done

#	echo "Begin first stage of scaffolding" | tee -a $stdout $stderr
#	scaf_sopra_bwa.sh $a $prefix-contig.$L.fa $sm_ins $w $sm >> $stdout 2>> $stderr 
#	echo "Begin second stage of scaffolding" | tee -a $stdout $stderr
#	scaf_sopra_bwa.sh $a/scaf_lg $a/scaffolds_h$h\L$L\_w$w.fasta $lg_ins $w $lg >> $stdout 2>> $stderr 
	
