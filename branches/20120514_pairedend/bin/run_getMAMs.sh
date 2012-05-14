#!/bin/bash
##
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -m e
#$ -M andrew.j.tritt@gmail.com
echo $JOB_ID `hostname`
OUT=""
if [ -d "/data/scratch" ]; then
	OUT="/data/scratch/$USER"
elif [ -d "/state/partition1/" ]; then
	OUT="/state/partition1/$USER"
else
	OUT="/share/eisen-d4/$USER"
fi
if [ ! -d $OUT ]; then
	mkdir -p $OUT
fi
MAUVE="/home/koadman/software/jre1.6.0_20/bin/java -Xmx2000m -cp /home/atritt/src/mauve_snapshot/Mauve.jar"
export DISPLAY="merlot.genomecenter.ucdavis.edu:9"
DEST_DIR="/share/eisen-d6/halophile_illumina_data/aaron_assemblies/MAMs"
REF_DIR="/home/$USER/Halophiles/ncbi_halos"
SCAF_DIR="/share/eisen-d6/halophile_illumina_data/aaron_assemblies/scaffold"
script="/home/$USER/ngopt/trunk/tools/mauve_assembly_metrics/mauveAssemblyMetrics.pl"

function run {
	base=$1
	draft=$SCAF_DIR/$base.scaf.fasta
	outdir="$OUT/$base"
	if [ ! -d $outdir ]; then
		mkdir $outdir
	fi
	cd $outdir
	$MAUVE org.gel.mauve.assembly.ScoreAssembly -reference $REF -assembly $draft -reorder -outputDir reorder > mauve.out 2> mauve.err
	$script ./ > mam.out 2> mam.err
	mv $outdir $DEST_DIR
}

# Har
REF=$REF_DIR/Haloarcula_marismortui_ATCC_43049_uid57719.gbk
run atcc43049
run jcm10024
run jcm13557
run dsm12282
run dsm6131
#Har	sp GUBF-8	atccBAA-651
#Har	sp GUBF-9	atccBAA-652

# Hbt-Hsx-Nmn-Hmc
REF=$REF_DIR/Halobacterium_NRC_1_uid57769.gbk
run jcm11222
run dsm15987

# Hcc-Hac 
REF=$REF_DIR/Halalkalicoccus_jeotgali_B3_uid50305.gbk
run dsm18796
run jcm12892
run dsm1307
run dsm5350
run dsm8989
run jcm13552

# Hfx
REF=$REF_DIR/Haloferax_volcanii_DS2_uid46845.gbk
run jcm10717
run atccBAA-1513
run atcc33959
run jcm13917
run dsm14919
run dsm18310
#Hfx	sp GUBF-1	atccBAA-644
#Hfx	sp GUBF-2	atccBAA-645
#Hfx	sp GUBF-3	atccBAA-646

# Hqr-Hsn-Hgm
REF=$REF_DIR/Halogeometricum_borinquense_DSM_11551_uid54919.gbk
run jcm14848

# Hrr
REF=$REF_DIR/Halorubrum_lacusprofundi_ATCC_49239_uid58807.gbk
run jcm13560
run jcm13916
run dsm19288
run dsm10284
run jcm10118
run jcm9100
run jcm14265
run atcc700873
run jcm14978
run dsm21995
run jcm13561
run dsm1137
run dsm14210
run jcm10247
run dsm21966

# Htg
REF=$REF_DIR/Haloterrigena_turkmenica_DSM_5511_uid43501.gbk
run dsm5511
run jcm13563
run jcm13562
run jcm13891
run dsm11522
run jcm12890
run jcm13890
run jcm14663
run dsm3751
run dsm15624
run jcm10478

# Hvx-Ncc-Nbt-Hbf
run jcm12983
run jcm10879
run dsm3393
run dsm10524
run dsm18795

# Nab
REF=$REF_DIR/Natrialba_magadii_ATCC_43099_uid46245.gbk
run dsm13077
run dsm12278
run jcm10990
run jcm10989
run atcc43099
run dsm3394
run dsm12281
run jcm10635
run jcm14089
run jcm10636

