# Introduction #

This page will describe how to use NGOpt tools with Reptile to correct errors in Illumina reads.

# Prerequisites #
You will need the following packages/programs. They will need to be in your executable path:
  * [Reptile](http://aluru-sun.ece.iastate.edu/doku.php?id=reptile) - error-correction software
  * NGOpt toolkit
## Step 1.  Convert Fastq to Fasta+Qual Files. ##
In the directory containing `my_reads.fastq`, run:
```
  <ngopt_trunk>/tools/quality/reptile_fq-converter ./ ./  2 
```
> > This will create **three** files, `my_reads.fa` , `my_reads.q`, and `my_reads.map.txt`. Because Reptile changes the original sequence names to arbitrary names, I have posted a hack on NGOpt that will produce a `.map.txt` file. We will use this file later to rename our reads to the original names.
> > Note: All files ending in .fastq in the current directory will be processed.
## Step 2. Create Reptile Config-File ##
Copy the Reptile config-file template from the Reptile directory and change as needed.
```
  cp <path_to_reptile>/src/config~ ./my_reads.rec.conf
  vim my_reads.rec.conf
```
Edit `my_reads.rec.conf` to specify input and output files.
  * Modify line beginning with "`InFaFile`"  to contain `my_reads.fa`
  * Modify line beginning with "`IQFile`" to contain `my_reads.q`
  * Modify line beginning with "`OErrFile`" to contain the desired name of the output file.
    * This output file will not be corrected reads. We are going to run another command on this file to get our corrected reads. Let's call the file `reptile-output`
If you are more educated in the Reptile Error Correction process, feel free to modify any of other run configurations.
## Step 3. Run Reptile Error Correction ##
Run `reptile` with the config file from previous step.
```
  reptile my_reads.rec.conf
```
This will produce the file `reptile-output`, as mentioned in the previous step.
## Step 4. Merge Corrections ##
Now that we've identified our errors, we will merge them to get back corrected reads.
```
  reptile_merger my_reads.fa reptile-output my_reads.rec.fasta
```
The third argument to `reptile_merger` is the output file. This file will contain our corrected reads.
## Step 5. Rename Corrected Reads ##
To get back the original read names, we can use a script from NGOpt--`reptile_rename`.
```
<ngopt_trunk>/tools/quality/reptile_rename my_reads.map.txt my_reads.rec.fasta > my_reads.rec.renamed.fasta
```
We now have  corrected reads with their original names in `my_reads.rec.renamed.fasta`.
## Where to next? ##
If you have a paired-end library, you might want to re-pair reads with broken pairings. You can use the program `repair` in the NGOpt toolkit to fix broken pairings.
```
repair -p ./ -s .fasta my_reads.rec my_reads.rec.renamed.fasta
```
This will produce **three** files: `my_reads.rec_p1.fasta`, `my_reads.rec_p2.fasta`, and `my_reads.rec_up.fasta`. The third file contains reads with a missing pair.