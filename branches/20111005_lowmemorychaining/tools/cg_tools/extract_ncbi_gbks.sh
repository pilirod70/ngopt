#!/bin/bash
if [ $# -ne 3 ]; then
	echo "Usage: `basename $0` <list_file> <NCBI_Genomes_dir> <output_dir>"
	exit
fi
list=$1
NCBI_DIR=$2
OUT_DIR=$3
for bug in `cat $list`; do 
 echo $bug
 head -q -n 1 $NCBI_DIR/$bug/*.gbk | perl -p -i -e "s/\ +/\t/g" | cut -f 2 > $NCBI_DIR/$bug.2
 head -q -n 1 $NCBI_DIR/$bug/*.gbk | perl -p -i -e "s/\ +/\t/g" | cut -f 3 > $NCBI_DIR/$bug.3 
 paste $bug.3 $bug.2 | sort -g -r | cut -f2 > $NCBI_DIR/$bug.gbks
 rm $NCBI_DIR/$bug.[2,3] 
 for gbk in `cat $NCBI_DIR/$bug.gbks`; do 
  cat $NCBI_DIR/$bug/$gbk.gbk >> $NCBI_DIR/$bug.gbk
 done
 rm $NCBI_DIR/$bug.gbks
done
