#!/bin/bash
##
#$ -l mem_free=5G
#$ -pe threaded 3
#$ -V
#$ -cwd
#$ -S /bin/bash

# Authors: Andrew Tritt and Aaron Darling
# (c) 2011, Licensed under the GPL
base=$1

#
# shared storage paths that may need to be changed
# 
# repository of fastq files
fq="/share/eisen-d6/halophile_illumina_data/allfastq/$base"
# location of library description files
conf="/share/eisen-d6/halophile_illumina_data/assembly_line/lib_files/$base.libs"
# destination path for assemblies
DEST="/share/eisen-d2/koadman/haloasm/"

# nothing below this line should need to be changed
stdoe="$PWD/$JOB_NAME.*$JOB_ID"
pipeline="andrews_assembly_line.pl"
LOCKFILE="$DEST/liunjlkjadftydsfbg908"

echo $JOB_ID `hostname` $JOB_NAME $@
# create a temporary storage dir on the compute node
OUT=""
if [ -d "/data/scratch" ]; then
	OUT="/data/scratch/$USER"
elif [ -d "/state/partition1/" ]; then
	OUT="/state/partition1/$USER"
else
	OUT="/share/eisen-d6/$USER"
fi

mkdir -p $DEST/$base
mkdir -p $OUT
mkdir $OUT/$base.$JOB_ID
cd $OUT/$base.$JOB_ID
# check to make sure no other process is trying to read a fastq
# lock 
while [ ! `mktemp -q $LOCKFILE` ]; do
	sleep 10s	
done
echo "Fetching data"
for gz in $fq*.gz; do
	zcat $gz > `basename $gz .gz`
done
cp $conf .
# unlock it for others
rm $LOCKFILE

$pipeline $conf $base > $DEST/$base/stdout.$JOB_ID 2> $DEST/$base/stderr.$JOB_ID

rm -rf *.fastq

mv $OUT/$base.$JOB_ID/* $DEST/$base/
mv $stdoe $DEST/$base/

