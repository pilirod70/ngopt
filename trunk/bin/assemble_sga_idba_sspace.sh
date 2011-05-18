#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 2

. $HOME/.bashrc
. $HOME/.bash_profile

CURDIR=`pwd`
WORKDIR=/state/partition1/koadman/assembly/$JOB_ID
mkdir -p $WORKDIR
cp $1 $WORKDIR
cp $2 $WORKDIR
cp $3 $WORKDIR
cd $WORKDIR
assemble_sga_idba_sspace.pl $1 $2 $3 $4

rm $WORKDIR/$1 
rm $WORKDIR/$2 
rm $WORKDIR/$3
cp -r $WORKDIR $CURDIR
rm -rf $WORKDIR

