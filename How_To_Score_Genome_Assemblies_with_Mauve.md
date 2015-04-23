# Introduction #

Most assembly algorithms provide a large variety of parameters that can be adjusted to control the assembly process.  It may be necessary to tune these parameters to fit details of the sequencing strategy and genomes being assembled.  This document will illustrate how one can compare assemblies done by different algorithms or different parameter settings to select an ideal assembly strategy.  We assume that a high quality reference genome is available, and that sequence reads for that reference genome are available, and also that the reference genome is similar enough to the genomes being targeted for assembly that parameters which work well on the reference will also work well on the other genomes.  This software is designed to work with haploid and homozygous polyploid genomes only.  Use with diploids or polyploids with divergent alleles is possible but not supported at present.

The following guide applies to computing assembly metrics on the command line for a large number of genome assemblies.  Although we do not describe it here, it is also possible to activate the assembly scoring on a single assembly from within the Mauve GUI when a pairwise alignment of  reference genome and assembly has been loaded.

To see some examples of the output generated by Mauve Assembly Metrics, please look at [this supplementary PDF](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics_supplementary.pdf)

# Software requirements #

### Mauve ###
Please obtain a copy of Mauve dated June 10th 2011 or later from here:
http://gel.ahabs.wisc.edu/mauve/snapshots/

### Mauve Assembly Metrics ###
Then obtain a copy of the Mauve assembly metrics scripts from here:
http://ngopt.googlecode.com/svn/trunk/tools/mauve_assembly_metrics

Alternatively these can be checked out using subversion:

```
svn checkout http://ngopt.googlecode.com/svn/trunk/tools/mauve_assembly_metrics mauve_assembly_metrics
```

### R ###

You will also need the R statistical software, available from here:
http://r-project.org
It can be installed on Debian and Ubuntu distributions simply by running:
```
sudo apt-get install r-recommended
```
To download and install R-2.13.0 on Mac OS X try:
```
curl -o R-2.13.0.pkg http://cran.r-project.org/bin/macosx/R-2.13.0.pkg
sudo installer -pkg R-2.13.0.pkg -target /
```

### Java ###

Most computers already have Java installed.  Windows users mayneed to download and install it from http://java.sun.com

### A Graphical User Interface ###
Even though Mauve Assembly Metrics can be run from the command-line, the software requires access to a graphical user interface.  If running the software remotely on a server, one can use the `ssh -X servername` option to connect the local display to the remote server.  In the event that X forwarding is not possible we recommend using VNC [download here](http://www.realvnc.com/products/free/4.1/index.html). It will be necessary to set the DISPLAY environment variable to the screen created by `vncserver`, e.g. `export DISPLAY=:1` if VNC created X display ID 1.


# Data requirements #

  * A reference genome, either in FastA or GenBank format (preferred).  For the tutorial we will use [Haloferax\_volcanii\_DS2.gbk](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/Haloferax_volcanii_DS2.gbk)

  * Some genome assemblies.  These can also be in FastA or GenBank formats.  For the tutorial we will use [assembly1.fasta](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly1.fasta), [assembly2.fasta](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly2.fasta), and [assembly3.fasta](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly3.fasta)


# Computing Mauve Assembly Metrics #

The following instructions assume you have downloaded and unpacked the Mauve software into a directory called mauve\_snapshot\_directory and have also downloaded the test data files [Haloferax\_volcanii\_DS2.gbk](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/Haloferax_volcanii_DS2.gbk), [assembly1.fasta](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly1.fasta), [assembly2.fasta](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly2.fasta), and [assembly3.fasta](http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly3.fasta) into that directory.  If these files exist in different locations it will be necessary to adjust the paths accordingly in the following commands.

  1. **Use Mauve to score each assembly**
```
cd mauve_snapshot_directory

java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly1.fasta -reorder -outputDir assembly1_scores

java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly2.fasta -reorder -outputDir assembly2_scores

java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly3.fasta -reorder -outputDir assembly3_scores
```
> > Two caveats exist at this stage.  First, the program must be run from within the Mauve directory in order to find dependency jar libraries.  If the libraries in the "ext" folder are added to $CLASSPATH and the directory containing the progressiveMauve binary is added to $PATH then the program may run from within a different directory.  See the Mac OS X example below for a demonstration of this. Second, access to a GUI is required.  See the comments above for getting GUI support on headless and remote servers.
  1. **Generate summary plots comparing the assemblies**
```
mauveAssemblyMetrics.pl ./
```
  1. **Visualize PDFs of accuracy plots.**  A large number of PDFs will be generated in the current directory.  These plots can be used to decide on the assembly strategy that best meets your project's objectives.  For example, you may prefer the assembly with the highest fraction of coding regions reconstructed, or the assembly with the fewest structural errors, or the assembly that captures the true genomic content most precisely.
  1. **Examine the table of assembly characteristics** In addition to the visual output, a tab-delimited table containing assembly metrics for each assembly is reported in the file `summaries.txt`.  This file could be used in conjunction with an automated method to select assembly parameters, for example.

# A complete example on the Linux command line #
The following commands can be copied and pasted into a linux terminal.

```
#
# Get the software
#
wget http://gel.ahabs.wisc.edu/mauve/snapshots/2011/2011-08-31/linux-x64/mauve_linux_snapshot_2011-08-31.tar.gz
tar xvzf mauve_linux_snapshot_2011-08-31.tar.gz
cd mauve_snapshot_2011-08-31/
wget http://ngopt.googlecode.com/svn/trunk/tools/mauve_assembly_metrics/mauveAssemblyMetrics.R
wget http://ngopt.googlecode.com/svn/trunk/tools/mauve_assembly_metrics/mauveAssemblyMetrics.pl
chmod 755 mauveAssemblyMetrics.*

#
# Get the test data
#
wget http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/Haloferax_volcanii_DS2.gbk
wget http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly1.fasta
wget http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly2.fasta
wget http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly3.fasta

#
# run the scoring
#
java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly1.fasta -reorder assembly1_scores -outputDir assembly1_scores

java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly2.fasta -reorder assembly2_scores -outputDir assembly2_scores

java -cp Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly3.fasta -reorder assembly3_scores -outputDir assembly3_scores

#
# Generate plots of metrics
# ./ indicates the parent directory containing the assembly scoring directories -- the current working directory in our example
#
./mauveAssemblyMetrics.pl ./

```

At this stage a set of PDF files named `mauve_*` will exist in the current directory.  These files can be inspected to evaluate assembly strategies.



# A complete example on the Mac OS X 10.6 command line #
The following commands can be copied and pasted into a Mac OS X 10.6 terminal window.
Note that this requires Java version 6, which comes with OS X 10.6.  On older versions of OS X, Java 6 needs to be installed manually.

```
#
# Get the software
#
curl -o Mauve-snapshot_2011-08-31.dmg http://gel.ahabs.wisc.edu/mauve/snapshots/2011/2011-08-31/MacOS/Mauve-snapshot_2011-08-31.dmg
hdiutil attach Mauve-snapshot_2011-08-31.dmg
export PATH=$PATH:/Volumes/Mauve\ snapshot_2011-08-31/Mauve.app/Contents/MacOS/
export CLASSPATH="/Volumes/Mauve\ snapshot_2011-08-31/Mauve.app/Contents/Resources/Java/ext/*"
curl -o mauveAssemblyMetrics.R http://ngopt.googlecode.com/svn/trunk/tools/mauve_assembly_metrics/mauveAssemblyMetrics.R
curl -o  mauveAssemblyMetrics.pl http://ngopt.googlecode.com/svn/trunk/tools/mauve_assembly_metrics/mauveAssemblyMetrics.pl
chmod 755 mauveAssemblyMetrics.*

#
# Get the test data
#
curl -o Haloferax_volcanii_DS2.gbk http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/Haloferax_volcanii_DS2.gbk
curl -o assembly1.fasta  http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly1.fasta
curl -o assembly2.fasta http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly2.fasta
curl -o assembly3.fasta http://edhar.genomecenter.ucdavis.edu/~koadman/mauve_assembly_metrics/assembly3.fasta

#
# run the scoring
#
java -cp /Volumes/Mauve\ snapshot_2011-08-31/Mauve.app/Contents/Resources/Java/Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly1.fasta -reorder  -outputDir assembly1_scores

java -cp /Volumes/Mauve\ snapshot_2011-08-31/Mauve.app/Contents/Resources/Java/Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly2.fasta -reorder  -outputDir assembly2_scores

java -cp /Volumes/Mauve\ snapshot_2011-08-31/Mauve.app/Contents/Resources/Java/Mauve.jar org.gel.mauve.assembly.ScoreAssembly -reference Haloferax_volcanii_DS2.gbk -assembly assembly3.fasta -reorder  -outputDir assembly3_scores

#
# Generate plots of metrics
# ./ indicates the parent directory containing the assembly scoring directories -- the current working directory in our example
#
./mauveAssemblyMetrics.pl ./

```