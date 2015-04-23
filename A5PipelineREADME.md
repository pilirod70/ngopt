## Quickstart ##

```
a5_pipeline.pl read_1.fastq read_2.fastq mygenome
```

This will assemble the paired reads in the two files `read_1.fastq` and `read_2.fastq`.
The final scaffolded assembly will be saved in `mygenome.final.scaffolds.fasta`.


## Copyright and license ##
The A5 pipeline is (c) 2011, 2012 Andrew Tritt and Aaron Darling.
This A5 pipeline is free, open source software licensed under the GPLv3, please see the file LICENSE for details.
The A5 pipeline includes open-source components developed and copyright by 3rd parties: `tagdust`, `bwa`, `SGA`, `bowtie`, `IDBA`, `SSPACE`, and `samtools`. Please see their license agreements for further details.

## What is A5? ##
_A5_ is a pipeline for assembling DNA sequence data generated on the Illumina sequencing platform. This README will take you through the steps necessary for running _A5_.

## What A5 can't do ##
There are many situations where A5 is not the right tool for the job. In order to produce accurate results, A5 requires Illumina data with certain characteristics. A5 will likely not work well with Illumina reads shorter than around 80nt, or reads where the base qualities are low in all or most reads before 60nt. A5 assumes it is assembling homozygous haploid genomes. Use a different assembler for metagenomes and heterozygous diploid or polyploid organisms. Use a different assembler if a tool like FastQC reports your data quality is dubious. You've been warned!

## Requirements ##
A5 requires 64-bit Linux (kernel 2.6.15 or later) or Mac OS X 10.6 or later. A Java Runtime Environment is also required. Mac OS X includes Java. On Linux, check with your distribution provider for details on installing Java.

## Installation ##
Once you have downloaded and extracted the pipeline, it is sufficient to simply add the pipeline's bin/ directory to the `$PATH` environment variable:
```
export PATH=$PATH:/path/to/ngopt/bin
```
Please change /path/to/ngopt/bin appropriately, and put this command in ~/.bashrc or ~/.profile or another script that gets run at login time.
The pipeline does not need to be copied to a specific location in order to be installed. You do not need root or superuser or administrator access to install and use the pipeline with this approach.


## Usage details ##

```
Usage: a5_pipeline.pl [--begin=1-5] [--end=1-5] [--preprocessed] <lib_file> <out_base>

Or:    a5_pipeline.pl <Read 1 FastQ> <Read 2 FastQ> <out_base>

Or:    a5_pipeline.pl <Read 1,2 Interleaved FastQ> <out_base>

<out_base> is the base file name for all output files. When assembling from 
a single library, the fastq files may be given directly on the command line.
If using more than one library, a library file must be given as <lib_file>.
The library file must contain the filenames of all read files.

If --preprocessed is used, <lib_file> is expected to be the library file
created during step 2 of the pipeline, named <out_base>.preproc.libs. Note 
that this flag only applies if beginning pipeline after step 2.
```

#### Example usage 1 ####
```
a5_pipeline.pl my_libs assembly.out
```
#### Example usage 2 ####
```
a5_pipeline.pl my_reads.1.fastq my_reads.2.fastq assembly.out
```
All output will be contained within the directory <output base>. The above example will produce the following directory and files:

Running the pipeline using either of the two examples will produce the following output:

```
assembly.out.ec.fastq.gz              // Error corrected reads
assembly.out.contigs.fasta            // Contigs
assembly.out.crude.scaffolds.fasta    // Crude scaffolds:  Not checked for misassemblies
assembly.out.broken.scaffolds.fasta   // Broken scaffolds: Checked for misassemblies, but not rescaffolded
assembly.out.final.scaffolds.fasta    // Final scaffolds:  Checked for misassemblies, and rescaffolded
```

If no misassemblies are detected, `assembly.out.broken.scaffolds.fasta` will not be present.

### If you are using a library file, please read the following ###

### Making a library file ###

The library file must contain file paths to all libraries being assembled.
Separate libraries are delimited in the library file by `[LIB]`.
For a paired-end or mate-pair library, reads can be passed as two files, one for each end, or as a single shuffled (aka interleaved) file.
Additionally, unpaired reads can be passed as a single file.
For a single-end library, reads must be passed as a single file.
Reads must be in FASTQ format. In addition to file locations, you can also give an insert size for paired reads. Please note that this is NOT necessary. The given number will only be used when _A5_ is not able to reliably calculate the insert size.

Example library file:
```
[LIB]
p1=library1_p1.fastq
p2=library1_p2.fastq
up=library1_up.fastq
ins=350
[LIB]
shuf=library2.fastq
up=library2.unpaired.fq
ins=6500
[LIB]
up=unpaired_library.fastq
```


### Other notes ###

Various stages of the A5 pipeline estimate assembly parameters from subsets of the input reads. This parameter estimation can vary from run to run, and can be dependent on the order of the input reads. Therefore it is possible that two runs of A5 on the same data will not generate the exact same assembly.