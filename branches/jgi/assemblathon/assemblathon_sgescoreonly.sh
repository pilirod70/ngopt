#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -l mem_free=3G
#$ -pe threaded 2
#
# Usage: <reference> <assembly> <scoring dir>
MAUVEDIR="/home/koadman/software/mauve_2.4.0_beta/"
MAUVE="/home/koadman/software/jre1.6.0_20/bin/java -Xmx2000m -cp Mauve.jar "
WORKDIR="/state/partition1/koadman/assemblathon/$JOB_ID"
export TMPDIR="/tmp"

mkdir -p $WORKDIR
export DISPLAY="merlot.genomecenter.ucdavis.edu:1"
export LD_LIBRARY_PATH=$HOME/lib64:$HOME/lib
cd $MAUVEDIR
$MAUVE org.gel.mauve.assembly.ScoreAssembly -alignment $1 -outputDir $WORKDIR
rm $WORKDIR/*/*.sslist
mv $WORKDIR $3

