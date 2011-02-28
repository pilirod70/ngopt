#!/bin/bash

if [ $# -lt 2 ] ; then
	echo "Usage: `basename $0` <output_dir> <ref> <reads>"
	echo "<reads> can be read from stdin"
	echo "Final pileup file is printed to stdout"
	exit
fi



out=$1
prefix=$out/$ref_base.bwa
ref=$2
ref_base=`basename $ref`
reads=/dev/stdin
base=$ref_base
if [ $# -eq 3 ]; then
	reads=$3
	base=`basename $reads`
fi
bwa index -a is -p $prefix $ref > $out/index.out
bwa aln $prefix $reads > $out/$base.sai
bwa samse $prefix $out/$base.sai $reads > $out/$ref_base.sam

samtools view -bhS $out/$ref_base.sam | samtools sort -o - $out/$ref_base | samtools pileup -f $ref -  
