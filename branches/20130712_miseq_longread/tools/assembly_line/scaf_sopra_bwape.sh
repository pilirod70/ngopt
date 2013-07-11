#!/bin/bash

if [ $# -ne 5 ]; then
	echo "Usage: `basename $0` <outdir> <contigs> <reads.fasta> <minpair> <insert>"
	echo "Reads file must be shuffled and end with \".fasta\" to conform with SOPRA"
	echo "Final product will be at <outdir>/scaffolds_h2.2_L150_w<minpair>.fasta"
	exit -1
fi

	outdir=$1
	contigs=$2
	fasta=$3
	w=$4
	ins=$5

	prefix="$outdir/contigs_sopra"	
	c="10"
	L="150"
	h="2.2"
	if [ ! -d $outdir ]; then
		mkdir $outdir
	fi
	echo "Scaffolding contigs from $contigs with paired reads from $fasta and writing output to $outdir.";
	s1_prep_contigAseq_v1.4.0.pl -contig $contigs -mate $fasta -a $outdir 
	bwa index -a is -p $prefix $outdir/contigs_sopra.fasta 
	base=`basename $fasta .fasta`
	fasta="$outdir/${base}_sopra.fasta"	
	repair --paired -p $outdir/ -s .fasta ${base}_sopra $fasta
	fas1="$outdir/${base}_sopra_p1.fasta"
	fas2="$outdir/${base}_sopra_p2.fasta"
	sai1="$outdir/${base}_sopra_p1.sai"
	sai2="$outdir/${base}_sopra_p2.sai"
	bwa aln -f $sai1 $prefix $fas1
	bwa aln -f $sai2 $prefix $fas2
	sam="$outdir/$base.sam"	
	bwa sampe -a $ins -f $sam $prefix $sai1 $sai2 $fas1 $fas2  
	s2_parse_sam_v1.4.0.pl -sam $sam -a $outdir 
	sam="${sam}_parsed"
	s3_read_parsed_sam_v1.4.0.pl -parsed $sam -d $ins -c $c -a $outdir 
	s4_scaf_v1.4.0.pl -o $outdir/orientdistinfo_c$c -w $w -L $L -a $outdir -h $h 
	
