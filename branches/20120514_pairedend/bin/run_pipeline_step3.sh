#!/bin/bash
##
#$ -l mem_free=4G
#$ -pe threaded 3
#$ -V
#$ -cwd
#$ -S /bin/bash

# Authors: Andrew Tritt and Aaron Darling
# (c) 2011, Licensed under the GPL
if [ $# -ne 3 ]; then
	echo "Usage: run_pipeline_step3.sh <base> <final_output_directory> <contig_directory>"
	exit
fi

export JAVA_HOME="/home/koadman/software/jre1.6.0_20"
export PATH="/home/koadman/software/jre1.6.0_20/bin:$PATH"


base=$1
CTG_DIR=$3
#
# shared storage paths that may need to be changed
# 
# repository of fastq files
fq="/share/eisen-d6/halophile_illumina_data/allfastq/$base"
fq2="/share/eisen-d2/koadman/haloasm/allfastq/$base"
# location of library description files
conf="/share/eisen-d2/koadman/haloasm/lib_files/$base.libs"
ctg="$CTG_DIR/$base/$base.contigs.fasta"
# destination path for assemblies
#DEST="/share/eisen-d2/koadman/haloasm_aggressive/"
DEST=$2
# nothing below this line should need to be changed
stdoe="$PWD/$JOB_NAME.*$JOB_ID"
pipeline="a5_pipeline.pl"
LOCKFILE="$DEST/liunjlkjadftydsfbg908"

echo $JOB_ID `hostname` $JOB_NAME $@
which java
# create a temporary storage dir on the compute node
OUT=""
if [ -d "/data/scratch" ]; then
	OUT="/data/scratch/$USER"
elif [ -d "/state/partition1/" ]; then
	OUT="/state/partition1/$USER"
else
	OUT="/share/eisen-d6/$USER"
fi

if [ ! -d $DEST/$base ]; then
	mkdir -p $DEST/$base
fi
mkdir -p $OUT
mkdir $OUT/$base.$JOB_ID
cd $OUT/$base.$JOB_ID
# check to make sure no other process is trying to read a fastq
# lock 
while [ ! `mktemp -q $LOCKFILE` ]; do
	sleep 10s	
done
echo "Fetching data"
n=`ls $fq.* 2> /dev/null | wc -l`
if [ "$n" -ne "0" ]; then
	for gz in $fq.*.gz; do
		echo "$gz -> `basename $gz .gz`"
		zcat $gz > `basename $gz .gz`
	done
fi

n=`ls $fq2.* 2> /dev/null | wc -l`
if [ "$n" -ne "0" ]; then
	for gz in $fq2.*.gz; do
		echo "$gz -> `basename $gz .gz`"
		zcat $gz > `basename $gz .gz`
	done
fi
cp $conf .
cp $ctg .
# unlock it for others
rm $LOCKFILE
which $pipeline
$pipeline --begin=3 $conf $base > $DEST/$base/stdout.s3.$JOB_ID 2> $DEST/$base/stderr.s3.$JOB_ID

rm -rf *.fastq

mv $OUT/$base.$JOB_ID/* $DEST/$base/
mv $stdoe $DEST/$base/

