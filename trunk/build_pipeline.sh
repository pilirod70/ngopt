#!/bin/bash


function copy_bin {
	bin="bwa fish sga idba tagdust"
	for ex in $bin; do
		cp $1/$ex $findir/bin
	done
}

function copy_sspace {
	echo "Copying SSPACE to $findir"
	mkdir $sspace_dir
	cp SSPACE/SSPACE $sspace_dir/ 
	cp -r SSPACE/bin $sspace_dir/
	cp -r SSPACE/dotlib $sspace_dir/
}

function copy_bowtie {
	mkdir $sspace_dir/bowtie
	cp -r $1/bowtie* $sspace_dir/bowtie
}

function copy_adhoc {
	cp bin/a5_pipeline.pl bin/break_misassemblies.pl bin/GetInsertSize.jar bin/GetFishInput.jar $findir/bin
	cp adapter.fasta $findir/
	echo "Removing unnecessary .svn directories"
	for dir in `find $findir/ -name .svn`; do 
		rm -rf $dir; 
	done
}

function copy_exdat {
	echo "Copying example data and README to archive"
	mkdir $findir/example
	cp test/sequence/phiX_p[1,2].fastq test/sequence/phiX.libs $findir/example
	cp ngopt_a5pipeline.README $findir/README
}

function bundle_clean {
	if [ -f $findir.tar.gz ]; then
		rm $findir.tar.gz
	fi
	echo "Creating archive"
	tar -czvf $findir.tar.gz $findir
	rm -rf $findir
}

function reset {
	if [ -d $findir ]; then
		rm -rf $findir
	fi
	mkdir $findir
	mkdir $findir/bin
}

findir_base="ngopt_a5pipeline"

############################# Linux Build #############################

echo "Building pipeline for Linux x86"

findir="${findir_base}_linux-x86_64"

reset
echo "Copying Linux binaries to $findir"
copy_bin linux-x64
sspace_dir="$findir/bin/SSPACE"
copy_sspace
echo "Copying Linux specific bowtie binaries to SSPACE directory"
copy_bowtie linux-x64

echo "Compiling Java code and bundling into executable Jar"
ant -q compile jar
copy_adhoc
copy_exdat
bundle_clean

############################# Mac Build #############################

echo -e "\nBuilding pipeline for Mac OSX"

findir="${findir_base}_macOS-x86_64"

reset
echo "Copying Mac binaries to $findir"
copy_bin osx
sspace_dir="$findir/bin/SSPACE"
copy_sspace
echo "Copying Mac specific bowtie binaries to SSPACE directory"
copy_bowtie osx

copy_adhoc
copy_exdat
bundle_clean
