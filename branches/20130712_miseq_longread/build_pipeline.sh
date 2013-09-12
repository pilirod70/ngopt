#!/bin/bash


function copy_bin {
	bin="bwa sga scythe idba_ud tagdust samtools trimmomatic-0.30.jar"
	for ex in $bin; do
		cp -v $1/$ex $findir/bin
		if [ ! $? ]; then return 1; fi;
	done
}

function copy_sspace {
	echo "Copying SSPACE to $findir"
	mkdir $sspace_dir
	cp -v SSPACE/SSPACE $sspace_dir/ && \
	cp -v -r SSPACE/bin $sspace_dir/ && \
	cp -v -r SSPACE/dotlib $sspace_dir/
}

function copy_bowtie {
	mkdir $sspace_dir/bowtie
	cp -vr $1/bowtie* $sspace_dir/bowtie
}

function copy_adhoc {
	cp -v bin/a5_pipeline.pl bin/GetInsertSize.jar bin/A5qc.jar $findir/bin && \
	cp -v adapter.fasta scythe_adapter.fasta $findir/
	if [ ! $? ]; then return 1; fi
	echo "Removing unnecessary .svn directories"
	for dir in `find $findir/ -name .svn`; do 
		rm -rf $dir; 
	done
}

function copy_exdat {
	echo "Copying example data and README to archive"
	mkdir -p $findir/example
	cp -v test/sequence/phiX_p[1,2].fastq test/sequence/phiX.libs test/sequence/phiX.a5.final.scaffolds.fasta $findir/example && \
#	cp -v ngopt_a5pipeline.README $findir/README
	if [ `uname` = "Darwin" ]; then
		curl http://code.google.com/p/ngopt/wiki/A5PipelineREADME > A5PipelineREADME.html	
	else
		wget -O A5PipelineREADME.html http://code.google.com/p/ngopt/wiki/A5PipelineREADME
	fi
	cp A5PipelineREADME.html $findir/README.html
	cp test.a5.sh $findir
}

function bundle_clean {
	rm -f $findir.tar.gz
	echo "Creating archive $findir.tar.gz"
	tar -czvf $findir.tar.gz $findir
	rm -rf $findir
}

function reset {
	rm -rf $findir
	mkdir -p $findir/bin
}

findir_base="ngopt_a5pipeline"

############################# Linux Build #############################


function build_linux_x64 {

	echo "Building pipeline for Linux x64"
	findir="${findir_base}_linux-x64_`date +%Y%m%d`"

	reset && \
	echo "Copying Linux binaries to $findir" && \
	copy_bin linux-x64 && \
	sspace_dir="$findir/bin/SSPACE" && \
	copy_sspace && \
	echo "Copying Linux specific bowtie binaries to SSPACE directory" && \
	copy_bowtie linux-x64 && \

	echo "Compiling Java code and bundling into executable Jar" && \
	ant -q compile jar && \
	copy_adhoc && \
	copy_exdat && \
	bundle_clean
	return 0
}

############################# Mac Build #############################

function build_osx {

	echo -e "\nBuilding pipeline for Mac OSX"
	findir="${findir_base}_macOS-x64_`date +%Y%m%d`"

	reset && \
	echo "Copying Mac binaries to $findir" && \
	copy_bin osx && \
	sspace_dir="$findir/bin/SSPACE" && \
	copy_sspace && \
	echo "Copying Mac specific bowtie binaries to SSPACE directory" && \
	copy_bowtie osx && \
	copy_adhoc && \
	copy_exdat && \
	bundle_clean
	return 0
}

build_linux_x64
if [ ! $? ]; then echo "Error building 64-bit linux pipeline" ; fi

build_osx
if [ ! $? ]; then echo "Error building Mac OS X pipeline"; fi

