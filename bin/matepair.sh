#!/bin/sh
#$ -cwd
#$ -S /bin/bash

../createMatePairs.pl H_volcanii_DS2_chr1.fasta.mapping.filtered 2847757 1423878 8000 1000 > chr1.pairs
../createMatePairs.pl H_volcanii_DS2_chr2.fasta.mapping.filtered 437906 218953 8000 1000 > chr2.pairs
../createMatePairs.pl H_volcanii_DS2_chr3.fasta.mapping.filtered 6359 6359 8000 1000 > chr3.pairs
../createMatePairs.pl H_volcanii_DS2_chr4.fasta.mapping.filtered 85092 42546 8000 1000 > chr4.pairs
../createMatePairs.pl H_volcanii_DS2_chr5.fasta.mapping.filtered 635786 317893 8000 1000 > chr5.pairs

