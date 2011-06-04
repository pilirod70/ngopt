#!/bin/bash
##
#$ -l mem_free=5G
#$ -pe threaded 4
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -m e
#$ -M andrew.j.tritt@gmail
# I see that andrew tritt wants to get email every time someone runs this script :)
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
stdoe="$PWD/$JOB_NAME.{e,o,pe,po}$JOB_ID"
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
if [ ! -d $OUT ]; then
	mkdir -p $OUT
fi
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

$pipeline $conf $base > pipeline.out 2> pipeline.err

rm -rf *.fastq

mv $OUT/$base.$JOB_ID $DEST/$base
mv $stdoe $DEST/$base/

