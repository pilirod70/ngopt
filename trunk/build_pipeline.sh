#!/bin/bash

findir_base="ngopt_a5pipeline"

############################# Linux Build #############################

echo "Building pipeline for Linux x86"

findir="${findir_base}_linux-x86_64"

if [ ! -d $findir ]; then
	mkdir $findir
else 
	rm -rf $findir/*
fi

echo "Copying Linux binaries to $findir"
bin="bwa fish sga idba"
for ex in $bin; do
	cp linux-x86/$ex $findir/
done

echo "Copying SSPACE to $findir"
sspace_dir="$findir/SSPACE"
mkdir $sspace_dir
cp SSPACE/SSPACE $sspace_dir/ 
cp -r SSPACE/bin $sspace_dir/
cp -r SSPACE/dotlib $sspace_dir/
echo "Copying Linux specific bowtie binaries to SSPACE directory"
cp -r ../vendor/SSPACE_linux-x86_64/current/bowtie $sspace_dir
echo "Removing unnecessary bowtie executables"
rm $sspace_dir/bowtie/*debug*
#chmod +x $sspace_dir/bowtie/bowtie*
echo "Compiling Java code and bundling into executable Jar"
ant -q compile jar

cp bin/aaa_assembly_line.pl bin/break_misassemblies.pl GetFishInput.jar $findir/
chmod +x $findir/aaa_assembly_line.pl

echo "Removing unnecessary .svn directories"
for dir in `find $findir/ -name .svn`; do 
	rm -rf $dir; 
done

echo "Creating archive"
tar -czvf $findir.tar.gz $findir
rm -rf $findir

############################# Mac Build #############################

echo -e "\nBuilding pipeline for Mac OSX"

findir="${findir_base}_macOS-x86_64"

if [ ! -d $findir ]; then
	mkdir $findir
else 
	rm -rf $findir/*
fi

echo "Copying Mac binaries to $findir"
bin="bwa fish sga idba"
for ex in $bin; do
	cp osx/$ex $findir/
done

echo "Copying SSPACE to $findir"
sspace_dir="$findir/SSPACE"
mkdir $sspace_dir
cp SSPACE/SSPACE $sspace_dir/ 
cp -r SSPACE/bin $sspace_dir/
cp -r SSPACE/dotlib $sspace_dir/
echo "Copying Mac specific bowtie binaries to SSPACE directory"
cp -r ../vendor/SSPACE_macOS-x86_64/current/bowtie $sspace_dir
echo "Removing unnecessary executables and giving Mac bowtie binaries executable permissions"
rm $sspace_dir/bowtie/*debug*
chmod +x $sspace_dir/bowtie/bowtie*

cp bin/aaa_assembly_line.pl bin/break_misassemblies.pl GetFishInput.jar $findir/
chmod +x $findir/aaa_assembly_line.pl

echo "Removing unnecessary .svn directories"
for dir in `find $findir/ -name .svn`; do 
	rm -rf $dir; 
done

echo "Creating archive"
tar -czvf $findir.tar.gz $findir

rm -rf $findir
