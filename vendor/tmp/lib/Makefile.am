AM_LDFLAGS = -fopenmp
AM_CXXFLAGS = -Wall -O2 -fopenmp
INCLUDES = -I$(top_srcdir)/include/common -I$(top_srcdir)/include/assembly 

libassembly_a_SOURCES = \
	$(top_srcdir)/include/common/globals.h \
	$(top_srcdir)/include/common/SequenceReader.h \
	$(top_srcdir)/include/common/Kmer.h \
	$(top_srcdir)/include/common/BitOperation.h \
	$(top_srcdir)/include/common/Contig.h \
	$(top_srcdir)/include/common/Log.h \
	$(top_srcdir)/include/common/Utils.h \
	$(top_srcdir)/include/common/SequenceWriter.h \
	$(top_srcdir)/include/common/Sequence.h \
	$(top_srcdir)/include/common/Read.h \
	$(top_srcdir)/include/common/KmerVector.h \
	$(top_srcdir)/include/common/BlatRecord.h \
	$(top_srcdir)/include/common/OptionSet.h \
	$(top_srcdir)/src/common/BlatRecord.cpp \
	$(top_srcdir)/src/common/Kmer.cpp \
	$(top_srcdir)/src/common/Log.cpp \
	$(top_srcdir)/src/common/KmerVector.cpp \
	$(top_srcdir)/src/common/BitOperation.cpp \
	$(top_srcdir)/src/common/SequenceWriter.cpp \
	$(top_srcdir)/src/common/Utils.cpp \
	$(top_srcdir)/src/common/Sequence.cpp \
	$(top_srcdir)/src/common/OptionSet.cpp \
	$(top_srcdir)/src/common/Read.cpp \
	$(top_srcdir)/src/common/Contig.cpp \
	$(top_srcdir)/src/common/globals.cpp \
	$(top_srcdir)/src/common/SequenceReader.cpp \
	$(top_srcdir)/include/assembly/AssemblyUtility.h \
	$(top_srcdir)/include/assembly/Component.h \
	$(top_srcdir)/include/assembly/HashAlign.h \
	$(top_srcdir)/include/assembly/ContigGraph.h \
	$(top_srcdir)/include/assembly/HashNode.h \
	$(top_srcdir)/include/assembly/ComponentNode.h \
	$(top_srcdir)/include/assembly/KmerNode.h \
	$(top_srcdir)/include/assembly/ContigNode.h \
	$(top_srcdir)/include/assembly/ContigPath.h \
	$(top_srcdir)/include/assembly/ContigBranchGroup.h \
	$(top_srcdir)/include/assembly/ConnectionGraph.h \
	$(top_srcdir)/include/assembly/HashGraph.h \
	$(top_srcdir)/include/assembly/Path.h \
	$(top_srcdir)/include/assembly/ComponentGraph.h \
	$(top_srcdir)/include/assembly/BranchGroup.h \
	$(top_srcdir)/include/assembly/RNAGraph.h \
	$(top_srcdir)/include/assembly/ScaffoldGraph.h \
	$(top_srcdir)/include/assembly/HashTable.h \
	$(top_srcdir)/include/assembly/AbstractNode.h \
	$(top_srcdir)/include/assembly/NodeImpl.h \
	$(top_srcdir)/src/assembly/AssemblyUtility.cpp \
	$(top_srcdir)/src/assembly/HashGraph.cpp \
	$(top_srcdir)/src/assembly/ContigGraph.cpp \
	$(top_srcdir)/src/assembly/ContigBranchGroup.cpp \
	$(top_srcdir)/src/assembly/ScaffoldGraph.cpp \
	$(top_srcdir)/src/assembly/Path.cpp \
	$(top_srcdir)/src/assembly/KmerNode.cpp \
	$(top_srcdir)/src/assembly/BranchGroup.cpp \
	$(top_srcdir)/src/assembly/AbstractNode.cpp \
	$(top_srcdir)/src/assembly/NodeImpl.cpp \
	$(top_srcdir)/src/assembly/ComponentGraph.cpp \
	$(top_srcdir)/src/assembly/Component.cpp \
	$(top_srcdir)/src/assembly/HashAlign.cpp \
	$(top_srcdir)/src/assembly/HashNode.cpp \
	$(top_srcdir)/src/assembly/ConnectionGraph.cpp \
	$(top_srcdir)/src/assembly/ContigPath.cpp \
	$(top_srcdir)/src/assembly/RNAGraph.cpp \
	$(top_srcdir)/src/assembly/ComponentNode.cpp \
	$(top_srcdir)/src/assembly/HashTable.cpp \
	$(top_srcdir)/src/assembly/ContigNode.cpp 

noinst_LIBRARIES = \
	libassembly.a 

