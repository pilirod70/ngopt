#!/bin/bash

function run {
	base=$1
	kcut=$2
	minpair=$3

	if [ ! -d $base ]; then
		mkdir $base
	fi
	echo $base	

	ins500="500_insert/$base/$base-500bp_shuf.rec.fasta"
	ins6kb="6kb_insert/$base/$base-6K_shuf.rec.rc2.fasta"
	fa_reads="$base.rec.fasta"
#	cat $ins500 $ins6kb > $fa_reads 

#	idba -r $fa_reads -o $base/$base.rec.idba --minCount $kcut > $base.out 2> $base.err 
#	mv $base.rec.idba* $base/
	echo "Done contiging\n" >> $base.out 2>> $base.err 
	s1_prep_contigAseq_v1.4.0.pl -contig $base/$base.rec.idba-contig.fa -mate $ins500 $ins6kb -a $base >> $base.out 2>> $base.err 

	prefix="$base/$base.rec.idba"
	
	bwa index -a is -p $prefix $base/contigs_sopra.fasta >> $base.out 2>> $base.err 

	ins6kb="$base/$base-6K_shuf.rec.rc2_sopra.fasta"
	ins500="$base/$base-500bp_shuf.rec_sopra.fasta"	
	sai6kb="$base/$base-6K_shuf.rec.sai"
	sai500="$base/$base-500bp_shuf.rec.sai"

	bwa aln -f $sai6kb $prefix $ins6kb >> $base.out 2>> $base.err 
	bwa aln -f $sai500 $prefix $ins500 >> $base.out 2>> $base.err 

	sam6kb="$base/$base-6K_shuf.rec.sam"
	sam500="$base/$base-500bp_shuf.rec.sam"	

	bwa samse -f $sam6kb $prefix $sai6kb $ins6kb >> $base.out 2>> $base.err
	bwa samse -f $sam500 $prefix $sai500 $ins500 >> $base.out 2>> $base.err

	s2_parse_sam_v1.4.0.pl -sam $sam6kb -a $base >> $base.out 2>> $base.err 
	s2_parse_sam_v1.4.0.pl -sam $sam500 -a $base >> $base.out 2>> $base.err 

	sam6kb="${sam6kb}_parsed"
	sam500="${sam500}_parsed"
		
	c="10"
	s3_read_parsed_sam_v1.4.0.pl -parsed $sam6kb -d 6000 -parsed $sam500 -d 500 -c $c -a $base >> $base.out 2>> $base.err 
	s4_scaf_v1.4.0.pl -o $base/orientdistinfo_c$c -w $minpair -L 150 -a $base -h 3.89 >> $base.out 2>> $base.err 

}

#run jcm10879     6 15 &
#run jcm13917     7 15 &
#run dsm3757      9 15 &
#run dsm21966    17 15 &
#run jcm13562     7 15 &
#run dsm3751      8 15 &
run jcm14089     7 15 &
