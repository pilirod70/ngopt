/**
 * @file tidba.cpp
 * @brief Iterative De Bruijn graph de novo short read Assembler for Transcriptome
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "config.h"

#include "globals.h"

#include "AssemblyUtility.h"
#include "Log.h"
#include "OptionSet.h"
#include "Read.h"
#include "RNAGraph.h"
#include "Sequence.h"
#include "SequenceWriter.h"
#include "Utils.h"

#include <unistd.h>
#include <sys/wait.h>

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

struct TIDBAOption
{
    std::string prefix;
    std::string readfile;
    std::string long_readfile;
    std::string graphfile;
    std::string contigfile;
    std::string align_file;
    std::string kmer_file;
    std::string long_contigfile;
    std::string path_file;
    std::string isoform_file;
    int mink;
    int modk;
    int maxk;
    int prefix_length;
    int min_count;
    int trim;
    double cover;
    int min_pairs;
    bool is_scaffold;
    uint64 mask;

    TIDBAOption()
    {
        prefix = "out";
        mink = 25;
        modk = 50;
        maxk = 50;
        prefix_length = 3;
        min_count = 2;
        trim = 0;
        cover = 0;
        min_pairs = 5;
        is_scaffold = false;
    }

    void Compute()
    {
        graphfile = prefix + ".graph";
        contigfile = prefix + "-contig.fa";
        long_contigfile = prefix + "-contig-long.fa";
        align_file = prefix + ".align";
        kmer_file = prefix + ".kmer";
        path_file = prefix + ".path";
        isoform_file = prefix + "-isoforms.fa";
    }
};

TIDBAOption option;
AssemblyUtility assembly_utitlity;

void BuildKmer()
{
    kmerLength = option.mink;

    assembly_utitlity.BuildKmerFile(option.min_count, option.prefix_length, option.kmer_file);
}

void BuildSimpleGraph()
{
    kmerLength = option.mink;
    HashGraph *hash_graph = new HashGraph();
    hash_graph->ReadFrom(option.kmer_file);
    assembly_utitlity.AddInternalKmer(*hash_graph);

    hash_graph->RefreshEdges();
    hash_graph->RemoveLowCoverageDeadEnd(2*option.maxk, sqrt(hash_graph->AverageCoverage()));
    hash_graph->Trim(2);
    hash_graph->RemoveLowCoverageContigs(sqrt(hash_graph->AverageCoverage()));
    LogMessage("total kmer %lld\n", hash_graph->NumNodes());

    //Vector<Contig> result_contigs;
    deque<Contig> result_contigs;
    hash_graph->Assemble(result_contigs);
    delete hash_graph;

    ContigGraph *contig_graph = new ContigGraph;
    contig_graph->Initialize(result_contigs);
    contig_graph->WriteTo(option.graphfile);
    delete contig_graph;
}

void Iterate()
{
    kmerLength = option.mink;
    ContigGraph *contig_graph = new ContigGraph;
    contig_graph->ReadFrom(option.graphfile);
    assembly_utitlity.IterateContigGraph(*contig_graph, option.modk, false, false);
    contig_graph->SortNodes();
    contig_graph->WriteTo(option.contigfile);
    delete contig_graph;
}

void Iterate2()
{
    kmerLength = option.modk;
    assembly_utitlity.SetLongReadFile(option.path_file);

    ContigGraph *contig_graph = new ContigGraph;
    contig_graph->ReadFrom(option.contigfile);
    assembly_utitlity.IterateContigGraph(*contig_graph, option.maxk, false, false);
    contig_graph->SortNodes();
    contig_graph->WriteTo(option.long_contigfile);
    delete contig_graph;
}

void FindIsoforms()
{
    kmerLength = option.maxk;
    assembly_utitlity.SetLongReadFile(option.path_file);

    RNAGraph *scaffold_graph = new RNAGraph;
    scaffold_graph->ReadFrom(option.long_contigfile);
    assembly_utitlity.BuildConnectionFromPath(*scaffold_graph, option.path_file);

    deque<Contig> contigs;
    scaffold_graph->FindIsoforms(contigs);

    FastaWriter writer(option.isoform_file);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        //contigs[i].Decode();
        //writer.WriteFormat(contigs[i], "isoforms%d", i);
        writer.Write(contigs[i], FormatString("isoforms%d", i));
    }
}

void Scaffold()
{
    kmerLength = option.modk;

    RNAGraph *scaffold_graph = new RNAGraph;
    scaffold_graph->ReadFrom(option.contigfile);
    assembly_utitlity.BuildConnection(*scaffold_graph, option.align_file);

    scaffold_graph->SetMinPairs(option.min_pairs);
    scaffold_graph->ComputeDistance();

    deque<Contig> contigs;
    scaffold_graph->FindValidConnection(contigs);

    FastaWriter writer(option.path_file);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        //contigs[i].Decode();
        //writer.WriteFormat(contigs[i], "contig%d_%d", i, contigs[i].sum_coverage);
        writer.Write(contigs[i], FormatString("contig%d_%d", i, contigs[i].sum_coverage));
    }
}

void AlignReads()
{
    kmerLength = option.modk;
    assembly_utitlity.AlignReads(option.contigfile, option.align_file);
}

void Usage()
{
    cout << "T-IDBA: Iterative De Bruijn graph short read Assembler for transcriptome" << endl;
    cout << "Version " << VERSION << endl;
    cout << endl;
    cout << "Usage: tidba --read read-file [--output out] [options]\n" << endl;
}

OptionSet option_set;

int main(int argc, char *argv[])
{
    option.maxk = 90;

    bool is_help = false;
    option_set.AddOption("help", "h", is_help, "produce help message");
    option_set.AddOption("read", "r", option.readfile, "read file");
    option_set.AddOption("output", "o", option.prefix, "prefix of output");
    option_set.AddOption("mink", "", option.mink, "minimum k value");
    option_set.AddOption("modk", "", option.modk, "moderate k value");
    option_set.AddOption("maxk", "", option.maxk, "maximum k value");
    option_set.AddOption("minCount", "", option.min_count, "filtering threshold for each k-mer");
    option_set.AddOption("minPairs", "", option.min_pairs, "minimum number of pair-end connections to join two contigs");
    option_set.AddOption("prefixLength", "", option.prefix_length, "length of the prefix of k-mer used to split k-mer table");

    try 
    {
        option_set.ProcessOptions(argc, argv);
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
        Usage();
        cout << "Allowed Options: " << endl;
        cout << option_set.ToString() << endl;
        //cout << OptionDescriptions() << endl;
        exit(1);
    }

    if (is_help || option.readfile == "")
    {
        Usage();
        cout << "Allowed Options: " << endl;
        cout << option_set.ToString() << endl;
        //cout << OptionDescriptions() << endl;
        exit(1);
    }

    option.Compute();
    kmerLength = option.mink;

    assembly_utitlity.SetReadFile(option.readfile);
    assembly_utitlity.SetLongReadFile(option.long_readfile);

    //if (!IsExist(option.kmer_file))
        BuildKmer();

    //if (!IsExist(option.graphfile))
        BuildSimpleGraph();

    //if (!IsExist(option.contigfile))
        Iterate();

    //if (!IsExist(option.align_file))
        AlignReads();

    //if (!IsExist(option.path_file))
        Scaffold();

    //if (!IsExist(option.long_contigfile))
        Iterate2();

    //if (!IsExist(option.isoform_file))
        FindIsoforms();

    return 0;
}

