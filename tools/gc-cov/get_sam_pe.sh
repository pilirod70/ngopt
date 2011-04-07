#!/bin/bash

if [ $# -ne 5 ] ; then
	echo "Usage: `basename $0` <output_dir> <ref> <pair1> <pair2> <max_ins>"
	echo "Final pileup file is printed to stdout"
	exit
fi

out=$1
ref=$2
ref_base=`basename $ref`
p1=$3
b1=`basename $p1`
p2=$4
b2=`basename $p2`
ins=$5
prefix=$out/$ref_base.bwa

bwa index -a is -p $prefix $ref > $out/index.out
bwa aln $prefix $p1 > $out/$b1.sai
bwa aln $prefix $p2 > $out/$b2.sai
bwa sampe -n 1 -a $ins $prefix $out/$b1.sai $out/$b2.sai $p1 $p2 
