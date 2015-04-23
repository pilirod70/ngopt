# Introduction #
This page will describe how to preprocess custom-barcode multiplexed Illumina data.
# Prerequisites #
You will need the following software in your executable path.
  * [TagDust](http://genome.gsc.riken.jp/osc/english/dataresource/) - eliminate artifactual sequence from NGS data
  * [splitBC](http://code.google.com/p/ngopt/source/browse/trunk/tools/splitbc/splitBC?r=70) - a script to split barcoded reads. Found in ngopt repository under trunk/tools/splitbc
If your reads are paired, you will also need the following program:
  * [repair](http://code.google.com/p/ngopt/source/browse/trunk/tools/pair_reads/repair.cpp?r=85) - a program for pairing and shuffling Illumina reads. Found in ngopt repository under trunk/tools/pair\_reads

You will need the following sequence files:
  * [adapter.fasta](http://code.google.com/p/ngopt/source/browse/trunk/adapter.fasta?r=85) - Illumina adapter sequences.
  * a tab-delimited file with your barcode sequences. This file should look something like this:
```
 BC1    ATGCTA
 BC2    TTGACT
 BC3    GTAGGT
```
  * a single file containing all of your reads. They do not need to be shuffled if you have paired reads
# Preprocess Multiplexed Data #
## Remove Artifactual DNA ##
Often times during library preparations, DNA fragments consisting only of adapter sequence are inserted and sequenced. To get rid of reads containing these sequences, we will use the program `tagdust`.
```
tagdust -o reads.clean.fastq -a reads.artifact.fastq adapters.fasta reads.fastq
```
where `reads.fastq` is our Illumina reads file. Reads without adapter sequence will be output to `reads.clean.fastq`
  * [TagDust](http://genome.gsc.riken.jp/osc/english/dataresource/) has one key parameter, the false discovery rate (FDR). By default this is set to 0.01. This is the false discovery rate of classifying a read as coming from **contaminant** library sequences. (i.e. lower FDR, fewer discarded reads)
## Split Barcoded Reads ##
To split our barcoded data set, we will use the `splitBC` from ngopt.
```
splitBC reads.clean.fastq --bcfile barcodes.txt --prefix ./ --suffix .fastq --bol 
```
A file will be created for each line in `barcodes.txt`. For example, if we used the example barcode file from above, we would now have **three** files, `BC1.fastq`, `BC2.fastq`, and `BC3.fastq`
## Re-pairing of Paired Reads ##
If you were processing paired-reads, your fastq files will likely contain broken pairings even if you shuffled them ahead of time. To fix this, you can use the `repair` program from ngopt.
```
repair -p ./ -s .fastq BC1 BC1.fastq
```
This will produce **three** files in the current working directory: `BC1_p1.fastq`, `BC1_p2.fastq`, `BC1_up.fastq`.
The third file, `BC1_up.fastq`, will contain reads whose pair was discarded due to adapter contamination or barcode degeneracy.

If shuffled reads are desired, you can use the `--shuf` option for `repair`. `BC1_shuf.fastq` will be created instead of `BC1_p1.fastq` and `BC1_p2.fastq`

**To correct errors with Reptile, consider using the ReptileErrorCorrection tutorial**
