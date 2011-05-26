#!/bin/bash
##
#$ -pe threaded 2
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -m e
#$ -M andrew.j.tritt@gmail
echo $JOB_ID `hostname`
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

OUTPUT_DIR=$1

FQDIR="/share/eisen-d6/halophile_illumina_data/allfastq"

DEST_DIR="/share/eisen-d6/halophile_illumina_data/sga-ec/complete"

function get_insert {
	echo `grep insert\ size $1 | sed s/^.*\:\ //g | sed s/\ \(.*$//g`
}

function get_sd {
	echo `grep insert\ size $1 | sed s/^.*\:\ [0-9][0-9]*\ \(//g | sed s/\ sigma\)//g`
}

# I expect the path to contigs to end with -contig.fa (IDBA naming convention)
function get_sampe {
	ctgs=$1
	prefix=`basename $ctgs -contig.fa`
	fq1=$2
	fq2=$3
	bwa index -a is -p $prefix >> get_sampe.out 2>> get_sampe.err 
	bwa aln $prefix $fq1 > sai1 2>> get_sampe.err 
	bwa aln $prefix $fq2 > sai2 2>> get_sampe.err 
	bwa sampe sai1 sai2 $fq1 $fq2 > sam 2> stderr 2>> get_sampe.err 
	rm sai1 sai2 $prefix* get_sampe.out 
	echo get_sampe.err 
}

function sgaec {
	outbase=$1
	shift 1
	sga-static preprocess -q 10 -f 20 -m 30 --phred64 $@ > $outbase.pp.fastq
	sga-static index -d 4000000 -t 4  $outbase.pp.fastq
	sga-static correct -k 31 -i 10 -t 4 -o $outbase.sgaec.fasta $outbase.pp.fastq
	rm $outbase.pp.*
}

function maxrdlen {
	echo `cat $@ | egrep -v '(^@)|^\+|^>' | awk ' { if ( length > x ) { x = length } }END{ print x }'`	
}

# Usage: minlink pair1 pair2 contigs insert rdlen
function minlink {
	numRdBases=`cat $1 $2 | egrep -v '^@|^\+|^>' | tr -d '\n' | wc -m`
	genomeLen=`cat $3 | egrep -v '^>' | tr -d '\n' | wc -m`
	insert=$4
	rdlen=$5
	echo $( $(($numRdBases * $insert)) / $((2 * $genomeLen * $rdlen)) )
}

function run {
	base=$1
	outdir=$OUT/$base
	if [ -d $outdir ]; then
		rm -rf $outdir
	fi 
	mkdir $outdir
	cd $outdir
	echo "==== $base ===="
	zcat $FQDIR/$base.*.fastq.gz > $base.fastq
	sgaec $base.fastq $base
	MAXK=`maxrdlen $base.sgaec.fasta`
	idba --read $base.sgaec.fasta --output $base.idba --mink 35 --maxk $MAXK 
	contigs="$base.idba-contig.fa"
	gzip $base.sgaec.fasta
	i=2
	while [ $i -le $# ]; do
		lib=${!i}
		i=$(($i+1))
		zcat $FQDIR/$base.${lib}_p1.fastq.gz > $base.${lib}_p1.fastq
		zcat $FQDIR/$base.${lib}_p2.fastq.gz > $base.${lib}_p2.fastq
		get_sampe $contigs $base.${lib}_p[1,2].fastq
		ins=`get_insert get_sampe.err`
		sd=`get_sd get_sampe.err`	
		err=$(($sd*4/$ins))
		rdlen=`maxrdlen $base.${lib}_p[1,2].fastq`
		# Usage: minlink pair1 pair2 contigs insert rdlen
		MINLINK=`minlink $base.${lib}_p[1,2].fastq $contigs $ins $rdlen`
		# libname pair1 pair1 insert uncertainty rc?  
		echo "$lib $base.${lib}_p1.fastq $base.${lib}_p2.fastq $ins $err 0 $MINLINK" >> library.txt	
		if [ -f $FQDIR/$base.${lib}_up.fastq.gz ]; then
			zcat $FQDIR/$base.${lib}_up.fastq.gz >> $base.unpaired.fastq
		fi
	done
	SSPACE -a 0.2 -o 10 -l library.txt -s $contigs -u $base.unpaired.fastq -b $base.sspace	
	mv $outdir $OUTPUT_DIR
	rm -Rf $outdir
}

