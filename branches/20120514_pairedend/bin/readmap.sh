#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -t 1-45
SCRATCH=/state/partition1/eisenlab/
mkdir -p $SCRATCH
/home/koadman/software/SHRiMP_1_2_1/bin/rmapper-ls -B reads.$SGE_TASK_ID $1 > $SCRATCH/reads.$SGE_TASK_ID.$1.mapping
mv $SCRATCH/reads.$SGE_TASK_ID.$1.mapping .
