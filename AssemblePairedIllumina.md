# Introduction #

Assembling genomes is fun


# Software Requirements #
  * Java - Most computers already have Java installed.
  * Perl - Most computers already have Perl installed.
  * [NGOpt pipeline](http://code.google.com/p/ngopt/downloads/list)

# Data Requirements #
  * One or more paired-end dataset in FastQ format in Illumina 1.3+ format. For the tutorial, we will use [GECO.tgz](http://edhar.genomecenter.ucdavis.edu/~andrew/ngopt_pipeline/GECO.tgz)
    * Single-end datasets may be included in addition to paired-end datasets

# Installing NGOpt pipeline #
After you download the archive from the link given above, extract the archive.
```
$ tar -xzf ngopt_pipeline_PLATFORM-x64.tar.gz 
```
Now add `ngopt_pipeline_PLATFORM-x64/bin/a5_pipeline.pl` to your executable path. You can also call the pipeline directly without adding it to your executable path.

# Running NGOpt pipeline #
The following instructions will assume you have downloaded and are using the test data from the tarball [GECO.tgz](http://edhar.genomecenter.ucdavis.edu/~andrew/ngopt_pipeline/GECO.tgz) and that these files are located in the current working directory.
```
cd ngopt_a5pipeline_PLATFORM-x64
bin/a5_pipeline.pl GECO.libs test
```
When complete, the current working directory will contain the final product of each step in the pipeline.

# A complete example on the Linux command line #
The following commands can be copied and pasted into a Linux Terminal
```
#
# Get software
#
wget http://ngopt.googlecode.com/files/ngopt_a5pipeline_linux-x64.tar.gz
tar -xzvf ngopt_a5pipeline_linux-x64.tar.gz

#
# Make a temp directory for organizational purposes, as many output files will be generated
#
mkdir ngopt_a5pipeline_macOS-x64/ngopt_tutorial
cd ngopt_a5pipeline_macOS-x64/ngopt_tutorial

#
# Get data
#
wget http://edhar.genomecenter.ucdavis.edu/~andrew/ngopt_pipeline/GECO.tgz
tar -xzvf GECO.tgz

#
# Run pipeline
#
../bin/a5_pipeline.pl GECO.libs test

```

# A complete example on the MacOS command line #
The following commands can be copied and pasted into the MacOS Terminal
```
#
# Get software
#
curl -o ngopt_a5pipeline_macOS-x64.tar.gz http://ngopt.googlecode.com/files/ngopt_a5pipeline_macOS-x64.tar.gz
tar -xzvf ngopt_a5pipeline_macOS-x64.tar.gz

#
# Make a temp directory for organizational purposes, as many output files will be generated
#
mkdir ngopt_a5pipeline_macOS-x64/ngopt_tutorial
cd ngopt_a5pipeline_macOS-x64/ngopt_tutorial

#
# Get data
#
curl -o GECO.tgz http://edhar.genomecenter.ucdavis.edu/~andrew/ngopt_pipeline/GECO.tgz
tar -xzvf GECO.tgz

#
# Run pipeline
#
../bin/a5_pipeline.pl GECO.libs test

```