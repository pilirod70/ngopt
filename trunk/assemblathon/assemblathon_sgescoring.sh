#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -l mem_free=3G
#$ -pe threaded 2
#
# Usage: <reference> <assembly> <scoring dir>
#
# Script to score an assembly in SGE using Mauve
# Depends on existence of a VNC or real X server
#
MAUVEDIR="/home/koadman/software/mauve_2.4.0_beta/"
MAUVE="/home/koadman/software/jre1.6.0_20/bin/java -Xmx2000m -cp Mauve.jar "
WORKDIR="/state/partition1/$USER/assemblathon/$JOB_ID"
export TMP="$WORKDIR/tmp"
export TMPDIR="$WORKDIR/tmp"

mkdir -p $WORKDIR
mkdir -p $TMPDIR
# the following will send useless display requests to an unattended Vnc
export DISPLAY="merlot.genomecenter.ucdavis.edu:1"
# the following is necessary on merlot to pick up a missing X-windows lib libXtst.so
export LD_LIBRARY_PATH=/home/koadman/lib64:/home/koadman/lib
cd $MAUVEDIR
#  MAUVE org.gel.mauve.contigs.ContigOrderer -output $WORKDIR -ref $1 -draft $2
$MAUVE org.gel.mauve.assembly.ScoreAssembly -reference $1 -assembly $2 -reorder $WORKDIR -outputDir $WORKDIR
rm $WORKDIR/*/*.sslist
mv $WORKDIR $3
rm -rf /tmp/mauve*
rm -rf $TMPDIR
