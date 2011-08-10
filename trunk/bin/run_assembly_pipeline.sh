#!/bin/bash
##
#$ -l mem_free=5G
#$ -pe threaded 3
#$ -V
#$ -cwd
#$ -S /bin/bash

# Authors: Andrew Tritt and Aaron Darling
# (c) 2011, Licensed under the GPL
if [ $# -ne 2 ]; then
	echo "Usage: run_assembly_pipeline.sh <base> <final_output_directory>"
	exit
fi

export JAVA_HOME="/home/koadman/software/jre1.6.0_20"
export PATH="/home/koadman/software/jre1.6.0_20/bin:$PATH"

base=$1

#
# shared storage paths that may need to be changed
# 
# repository of fastq files
fq="/share/eisen-d6/halophile_illumina_data/allfastq/$base"
fq2="/share/eisen-d2/koadman/haloasm/allfastq/$base"
# location of library description files
conf="/share/eisen-d2/koadman/haloasm/lib_files/$base.libs"
# destination path for assemblies
#DEST="/share/eisen-d2/koadman/haloasm_aggressive/"
DEST=$2
# nothing below this line should need to be changed
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
n=`ls $fq* 2> /dev/null | wc -l`
if [ "$n" -ne "0" ]; then
	for gz in $fq*.gz; do
		zcat $gz > `basename $gz .gz`
	done
fi

n=`ls $fq2* 2> /dev/null | wc -l`
if [ "$n" -ne "0" ]; then
	for gz in $fq2*.gz; do
		echo "basename $gz .gz  ==" `basename $gz .gz`
		zcat $gz > `basename $gz .gz`
	done
fi
cp $conf .
# unlock it for others
rm $LOCKFILE

$pipeline $conf $base > stdout.$JOB_ID 2> stderr.$JOB_ID

rm -rf *.fastq

stdoe="$PWD/$JOB_NAME.*$JOB_ID"
mv $OUT/$base.$JOB_ID/* $DEST/$base/
mv $stdoe $DEST/$base/

