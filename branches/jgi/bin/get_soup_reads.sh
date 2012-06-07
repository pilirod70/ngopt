#!/bin/bash
##
#$ -V
#$ -pe threaded 3
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
tar -xzvf $1
tar -xzvf $2
rm $LOCKFILE
SOUP1=$outdir/`basename $1 .tar.gz`
SOUP2=$outdir/`basename $2 .tar.gz`


function run {
	base=$1
	if [ ! -f $REFDIR/$base.scaf.fasta ] ; then
		echo "No scaffolds for $base"
		return 0
	fi	
	if [ -d $outdir/$base ]; then
		rm -rf $outdir/$base
	fi
	mkdir $outdir/$base
	cd $outdir/$base
	bwa index -a is -p ref $REFDIR/$base.scaf.fasta
	bwa aln -n 0 -t 4 ref $SOUP1 | bwa samse ref - $SOUP1 > soup1.sam
	bwa aln -n 0 -t 4 ref $SOUP2 | bwa samse ref - $SOUP2 > soup2.sam
	soupmap2fastq.pl soup1.sam soup2.sam $base.soup2.sam $base.soup2 	
	cp $base.soup2.sam $SAMDIR
	cp $base.soup2* $FQDIR
	rm -rf $outdir/$base
}

for org in `cat $3`; do
	run $org
done


