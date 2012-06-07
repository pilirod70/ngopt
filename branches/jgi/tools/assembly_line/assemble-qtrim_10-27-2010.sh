#!/bin/bash

SUFFIX="idba-contig.fa"

function run {

	base=$1
	kcut=$2
	ins=$3
	minPairs=$4
	if [ ! -d $base ]; then
		 mkdir $base
	fi

	if [ -e $base.out ]; then
		rm $base.out $base.err 
	fi

	echo "Quality trimming" >> $base.out
	echo "Quality trimming" >> $base.err 
	
	trimfq $base.fastq 25 13 >> $base.out 2>> $base.err 
 
	echo "Contiging" >> $base.out
	echo "Contiging" >> $base.err
	
	fa_reads=$base.fasta.htqs
	
	fq2fa $base.fastq.htqs $fa_reads >> $base.out 2>> $base.err
	 
	idba -r $fa_reads -o $base/$base.idba --minCount $kcut >> $base.out 2>> $base.err
 
	echo "Preping contigs and reads" >> $base.out
	echo "Preping contigs and reads" >> $base.err
 
	repair --shuf -s .fasta.htqs $base $fa_reads >> $base.out 2>> $base.err
	s1_prep_contigAseq_v1.4.0.pl -contig $base/$base.$SUFFIX -mate ${base}_shuf.fasta.htqs -a $base  >> $base.out 2>> $base.err
 
	echo "Mapping" >> $base.out  
	echo "Mapping" >> $base.err  
	bwa index -a is -p $base/$base $base/contigs_sopra.fasta  >> $base.out 2>> $base.err
	bwa aln -q 13 -f $base/${base}_shuf.sai $base/$base $base/${base}_shuf.fasta.htqs_sopra.fasta >> $base.out 2>> $base.err
	bwa samse -f $base/$base.sam $base/$base $base/${base}_shuf.sai $base/${base}_shuf.fasta.htqs_sopra.fasta >> $base.out 2>> $base.err
 
	echo "Scaffolding" >> $base.out
	echo "Scaffolding" >> $base.err
	s2_parse_sam_v1.4.0.pl -sam $base/$base.sam -a $base  >> $base.out 2>> $base.err
	s3_read_parsed_sam_v1.4.0.pl -parsed $base/$base.sam_parsed -d $ins -a $base  >> $base.out 2>> $base.err
	s4_scaf_v1.4.0.pl -o $base/orientdistinfo_c5 -w $minPairs -L 150 -a $base   >> $base.out 2>> $base.err
 
}
 
min_pair=5
 
run atcc700873   8 437 $min_pair &
run atccBAA_652  8 446 $min_pair &
run dsm10524	 9 425 $min_pair &
run dsm1307 396 $min_pair &
run dsm13077	 9 469 $min_pair &
run dsm3757	  7 392 $min_pair &
run dsm8989	  5 448 $min_pair &
run jcm12892	 7 456 $min_pair &
run jcm13561	 7 423 $min_pair &
run jcm13563    10 449 $min_pair & 
run jcm13891     6 421 $min_pair &
run jcm14978     6 487 $min_pair &
