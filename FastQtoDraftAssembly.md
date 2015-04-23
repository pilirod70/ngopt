# ~~DEPRECATED~~ #
This pipeline has been deprecated. A new pipeline is in the works. Check back soon for the new and improved pipeline.

# Introduction #

This page will describe how create draft genome assemblies from paired-reads in Fastq format.

# Prerequisites #
You will need the following packages/programs. They will need to be in your executable path:
  * [Reptile](http://aluru-sun.ece.iastate.edu/doku.php?id=reptile) - error-correction software
  * [IDBA](http://i.cs.hku.hk/~alse/hkubrg/projects/idba/) - Iterative De Bruijn graph de novo assembler.  See our [IDBA build notes](http://code.google.com/p/ngopt/wiki/IdbaBuildNotes)
  * [BWA](http://bio-bwa.sourceforge.net/) - short-read aligner using Burrows-Wheeler transform
  * [SOPRA](http://www.physics.rutgers.edu/~anirvans/SOPRA/) - statistical-optimization scaffolding
  * [repair](http://code.google.com/p/ngopt/source/browse/trunk/tools/pair_reads/repair.cpp?r=70) - pairs reads named using Illumina format sequence identifiers. Found in ngopt repository under trunk/tools/pair\_reads
  * [rmctgs](http://code.google.com/p/ngopt/source/browse/trunk/tools/seqtools/rmctgs.cpp?r=70) - removes small contigs from Fasta file. Found in ngopt repository under trunk/tools/seqtools.
  * [scaf\_sopra\_bwa.sh](http://code.google.com/p/ngopt/source/browse/trunk/tools/assembly_line/scaf_sopra_bwa.sh?r=70) - wrapper script for scaffolding with BWA and SOPRA. Found in ngopt repository under trunk/tools/assembly\_line.
  * [assembleIBS.sh](http://code.google.com/p/ngopt/source/browse/trunk/tools/assembly_line/assembleIBS.sh?r=70) - runs IDBA-BWA-SOPRA pipeline. Found in ngopt repository under trunk/tools/assembly\_line

# Assembling Paired Reads in FastQ format #
For this tutorial, we will assume you are using paired-end library in Illumina1.3+ Fastq format. We will refer to this Fastq file as my\_reads.fastq
## Reptile Error Correction ##
See ReptileErrorCorrection for instructions on correcting reads using Reptile and NGOpt tools.

## Format Reads ##
If you did not do error-correction with Reptile, you will need to convert your Fastq to Fasta format. We will refer to the converted file as `my_reads.fasta`.

You will need to shuffle your paired reads if they are not already shuffled. To do so, run **`repair`** like so:
```
repair --shuf -s .fastq my_reads my_reads.fasta
```
where `pe_shuf` is the output basename for our shuffled reads. This will produce two files named `my_reads_shuf.fasta` and `my_reads_up.fasta` in the current working directory.  The `-s` option or `repair`  can be used to output to a different directory. `my_reads_up.fasta` will contain reads with no pair.

## Assemble Corrected Reads ##
```
assembleIBS.sh ibs_asm 10 5 my_reads_shuf.fasta 500
```
  * `ibs_asm` is the basename for output files. It is also the name of the directory where all output files generated in the assembly line are stored.
  * the first integer argument, `10`, is the minimum frequency of a k-mer in the De Bruijn graph
  * the second integer argument, `5`, is the minimum number of pairs required to scaffold two contigs
  * `pe_shuf.rec.fasta` is our shuffled file of error corrected reads
  * the third integer argument, `500`, is the estimated insert size of your library.
> Please note that the values used here are by no means suggested values. These three values need to be tailored to your library

> Final scaffolds will be in `ibs_asm/scaffolds_h2.2_L150_w5.fasta`