#!/bin/bash

if [ $# -ne 4 ] ; then
	echo "Usage: `basename $0` <ref> <pair1> <pair2> <max_ins> <output_dir>"
	exit
fi

ref=$1
p1=$2
p2=$3
ins=$4
out=$5
prefix=$out/$ref.bwa

bwa index -a is -p $prefix $ref > $out/index.out
bwa aln $prefix $p1 > $out/$p1.sai
bwa aln $prefix $p2 > $out/$p2.sai
bwa sampe -n 1 -a $ins $prefix $out/$p1.sai $out/$p2.sai $p1 $p2 > $out/$ref.sam

samtools view -bhS $out/$ref.sam | samtools sort -o - $out/$ref | samtools pileup -f $ref -  
