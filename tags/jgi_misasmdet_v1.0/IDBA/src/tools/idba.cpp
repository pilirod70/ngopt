/**
 * @file idba.cpp
 * @brief Iterative De Bruijn graph de novo short read Assembler
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "config.h"

#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Read.h"
#include "Utils.h"
#include "AssemblyUtility.h"
#include "OptionSet.h"
#include "SequenceWriter.h"

#include <unistd.h>
#include <sys/wait.h>

#include <cstdio>
#include <algorithm>
#include <deque>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

struct IDBAOption
{
    std::string prefix;
    std::string readfile;
    std::string long_readfile;
    std::string graphfile;
    std::string contigfile;
    std::string scaffile;
    std::string align_file;
    std::string kmer_file;
    int mink;
    int maxk;
    int prefix_length;
    int min_count;
    int trim;
    double cover;
    int min_pairs;
    bool is_scaffold;

    IDBAOption()
    {
        prefix = "out";
        mink = 25;
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
        scaffile = prefix + "-contig-mate.fa";
        align_file = prefix + ".align";
        kmer_file = prefix + ".kmer";
    }
};

IDBAOption option;
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

    LogMessage("before operation: total kmer %lld edges %lld\n", hash_graph->NumNodes(), hash_graph->NumEdges());
    hash_graph->RemoveDeadEnd(2*option.maxk);

    if (option.cover != 0)
        hash_graph->RemoveLowCoverageContigs(option.cover);
    else
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
    assembly_utitlity.IterateContigGraph(*contig_graph, option.maxk);
    contig_graph->SortNodes();
    contig_graph->WriteTo(option.contigfile);

    deque<ContigNode> &nodes = contig_graph->GetContigNodes();
    vector<int64> lengths(nodes.size());
    for (unsigned i = 0; i < nodes.size(); ++i)
        lengths[i] = nodes[i].GetSize();

    cout << FormatString("contigs: %d N50: %lld max: %lld mean: %lld total: %lld n80 %lld",
            lengths.size(), Nxx(lengths, 0.5), Maximum(lengths), Average(lengths), Sum(lengths), Nxx(lengths, 0.8)) << endl;

    delete contig_graph;
}

void Scaffold()
{
    kmerLength = option.maxk;

    ScaffoldGraph *scaffold_graph = new ScaffoldGraph;
    scaffold_graph->ReadFrom(option.contigfile);

    assembly_utitlity.BuildConnection(*scaffold_graph, option.align_file);

    scaffold_graph->SetMinPairs(option.min_pairs);
    //scaffold_graph->SetMaxPaths(10);
    scaffold_graph->ComputeDistance();
    scaffold_graph->FindUniquePaths();

    deque<Contig> contigs;
    scaffold_graph->Scaffold(contigs);

    FastaWriter writer(option.scaffile);
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        //contigs[i].Decode();
        writer.Write(contigs[i], FormatString("contig%d", i));
        //writer.WriteFormat(contigs[i], "contig%d", i);
    }

    vector<int64> lengths(contigs.size());
    for (unsigned i = 0; i < contigs.size(); ++i)
        lengths[i] = contigs[i].Size();

    cout << FormatString("contigs: %d N50: %lld max: %lld mean: %lld total: %lld n80 %lld",
            lengths.size(), Nxx(lengths, 0.5), Maximum(lengths), Average(lengths), Sum(lengths), Nxx(lengths, 0.8)) << endl;

    delete scaffold_graph;
}

void AlignReads()
{
    kmerLength = option.maxk;
    assembly_utitlity.AlignReads(option.contigfile, option.align_file);
}

void Usage()
{
    cout << "IDBA: Iterative De Bruijn graph short read Assembler" << endl;
    cout << "Version " << VERSION << endl;
    cout << endl;
    cout << "Usage: idba --read read-file [--output out] [options]\n" << endl;
}

OptionSet option_set;

int main(int argc, char *argv[])
{
    bool is_help = false;
    option_set.AddOption("help", "h", is_help, "produce help message");
    option_set.AddOption("read", "r", option.readfile, "read file");
    option_set.AddOption("long", "l", option.long_readfile, "long read file");
    option_set.AddOption("output", "o", option.prefix, "prefix of output");
    option_set.AddOption("scaffold", "", option.is_scaffold, "use pair end information to merge contigs");
    option_set.AddOption("mink", "", option.mink, "minimum k value");
    option_set.AddOption("maxk", "", option.maxk, "maximum k value");
    option_set.AddOption("minCount", "", option.min_count, "filtering threshold for each k-mer");
    option_set.AddOption("cover", "", option.cover, "the cutting coverage for contigs");
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

    //if (!IsExist(option.align_file) && option.is_scaffold)
    if (option.is_scaffold)
        AlignReads();

    //if (option.is_scaffold && !IsExist(option.scaffile))
    if (option.is_scaffold)
        Scaffold();

    return 0;
}

