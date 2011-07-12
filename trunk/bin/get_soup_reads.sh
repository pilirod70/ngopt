#!/bin/bash
##
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -m e
#$ -M andrew.j.tritt@gmail
echo $JOB_ID `hostname` $JOB_NAME $@
OUT=""
if [ -d "/data/scratch" ]; then
	OUT="/data/scratch/atritt"
elif [ -d "/state/partition1/" ]; then
	OUT="/state/partition1/atritt"
else
	OUT="/share/eisen-d4/atritt"
fi
if [ ! -d $OUT ]; then
	mkdir -p $OUT
fi
LOCKFILE="/share/eisen-d6/halophile_illumina_data/akdjfhakdsjfhalkdjfh8"

REFDIR="/share/eisen-d6/halophile_illumina_data/aaron_assemblies/scaffold"
FINDIR="/share/eisen-d6/halophile_illumina_data/aaron_assemblies/soup2_sam"
FQDIR="/share/eisen-d2/koadman/haloasm/allfastq"

mkdir $OUT/filt_soup.$JOB_ID
outdir=$OUT/filt_soup.$JOB_ID
cd $outdir
while [ ! `mktemp -q $LOCKFILE` ]; do
	sleep 10s	
done
cp $1 $2 $outdir
rm $LOCKFILE
SOUP1=$outdir/`basename $1`
SOUP2=$outdir/`basename $2`


function run {
	base=$1
	mkdir $outdir/$base
	cd $outdir/$base
	get_sam.sh . $REFDIR/$base.scaf.fasta $SOUP1 > soup1.sam
	get_sam.sh . $REFDIR/$base.*.fasta $SOUP2 > soup2.sam
	soupmap2fastq.pl soup1.sam soup2.sam $base.soup2.sam $base.soup2 	
	mv $base.soup2.sam $SAMDIR
	mv $base.soup2* $FQDIR
	rm $outdir/$base
}

for org in `cat $3`; do
	run $org
done


