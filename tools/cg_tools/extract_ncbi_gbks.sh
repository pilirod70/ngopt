#!/bin/bash
if [ $# -ne 1 ]; then
	echo "Usage: `basename $0` <list_file>"
	exit
fi

list=$1
for bug in `cat $list`; do 
 head -q -n 1 $bug/*.gbk | perl -p -i -e "s/\ +/\t/g" | cut -f 2 > $bug.2
 head -q -n 1 $bug/*.gbk | perl -p -i -e "s/\ +/\t/g" | cut -f 3 > $bug.3 
 paste $bug.3 $bug.2 | sort -g -r | cut -f2 > $bug.gbks
 rm $bug.[2,3] 
done

for gbks in *.gbks; do
 org=`basename $gbks .gbks`
 for gbk in `cat $gbks`; do 
  cat $org/$gbk.gbk >> $org.gbk
 done
 rm $gbks
done
