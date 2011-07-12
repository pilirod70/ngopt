#!/bin/bash

if [ $# -lt 3 ] ; then
	echo "Usage: `basename $0` <output_dir> <ref> <reads>"
	echo "Final sam file is printed to stdout"
	exit
fi



out=$1
ref=$2
ref_base=`basename $ref`
prefix=$out/`basename $ref`
reads=$3
base=`basename $reads`
bwa index -a is -p $prefix $ref > /dev/null
bwa aln $prefix $reads | bwa samse $prefix - $reads   

#bwa aln $prefix $reads 2> $out/aln.err | bwa samse $prefix - $reads 2> $out/samse.err > $out/out.sam 
#samtools view -bhS $out/out.sam 2> $out/view.err > $out/out.bam 
#samtools sort -o $out/out.bam $out/out.sort 2> $out/sort.err > $out/out.sort.bam 
#samtools mpileup -f $ref $out/out.sort.bam
