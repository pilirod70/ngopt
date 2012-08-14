#!/bin/bash

function run {
	base=$1
	kcut=$2
	
	if [ -e $base.out ]; then
		rm $base.out $base.err 
	fi


#	fq2fa $base.fastq $base.fasta >> $base.out 2>> $base.err
	idba -r $base.fasta -o $base-idba --minCount $kcut >> $base.out 2>> $base.err 

}

run jcm10879  10  &
run jcm13917  12  &
run dsm3757   13  &
run dsm21966  19  &
run jcm13562  10  &
run dsm3751   12  &
run jcm14089  11  &

