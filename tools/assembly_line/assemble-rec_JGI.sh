#!/bin/bash

function run {
	base=$1
	kcut=$2
	minpair=$3

	if [ ! -d $base ]; then
		mkdir $base
	fi
	echo $base	

	sm_ins="small_ins/$base.sm.rec.fasta"
	lg_ins="large_ins/$base.sm.rec.fasta"
	fa_reads="$base.rec.fasta"
	cat $sm_ins $lg_ins > $fa_reads 

	idba -r $fa_reads -o $base/$base.rec.idba --minCount $kcut > $base.out 2> $base.err 
	echo "Done contiging\n" >> $base.out 2>> $base.err 

	s1_prep_contigAseq_v1.4.0.pl -contig $base/$base.rec.idba-contig.fa -mate $sm_ins $lg_ins -a $base >> $base.out 2>> $base.err 

	prefix="$base/$base.rec.idba"
	
	bwa index -a is -p $prefix $base/contigs_sopra.fasta >> $base.out 2>> $base.err 

	lg_ins="$base/$base.lg.rec_sopra.fasta"
	sm_ins="$base/$base.sm.rec_sopra.fasta"	
	lg_sai="$base/$base.lg.rec.sai"
	sm_sai="$base/$base.sm.rec.sai"

	bwa aln -f $lg_sai $prefix $lg_ins >> $base.out 2>> $base.err 
	bwa aln -f $sm_sai $prefix $sm_ins >> $base.out 2>> $base.err 

	lg_sam="$base/$base.sm.rec.sam"
	sm_sam="$base/$base.sm.rec.sam"	

	bwa samse -f $lg_sam $prefix $lg_sai $lg_ins >> $base.out 2>> $base.err
	bwa samse -f $sm_sam $prefix $sm_sai $sm_ins >> $base.out 2>> $base.err

	s2_parse_sam_v1.4.0.pl -sam $lg_sam -a $base >> $base.out 2>> $base.err 
	s2_parse_sam_v1.4.0.pl -sam $sm_sam -a $base >> $base.out 2>> $base.err 

	lg_sam="${lg_sam}_parsed"
	sm_sam="${sm_sam}_parsed"
		
	c="10"
	s3_read_parsed_sam_v1.4.0.pl -parsed $lg_sam -d 6000 -parsed $sm_sam -d 500 -c $c -a $base >> $base.out 2>> $base.err 
	s4_scaf_v1.4.0.pl -o $base/orientdistinfo_c$c -w $minpair -L 150 -a $base -h 3.89 >> $base.out 2>> $base.err 

}

run 4084251  10 &
run 4085576  10 &
run 4090482  10 &
run 4084251 80 20 &
run 4085576 120 20 &
run 4090482 150 20 &
