/**
 * @file metaidba.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-25
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
#include "ComponentGraph.h"

#include <unistd.h>
#include <sys/wait.h>

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>
#include <queue>

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
    std::string long_path_file;
    std::string component_file;
    std::string consensus_file;
    std::string alignment_file;
    int mink;
    int maxk;
    int prefix_length;
    int min_count;
    int trim;
    double cover;
    int min_pairs;
    bool is_scaffold;
    bool is_validate;
    bool is_connect;

    int num_species;
    int times;

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
        is_connect = false;

        num_species = 1;
        times = 1;

        is_validate = false;
        //is_validate = true;
    }

    void Compute()
    {
        graphfile = prefix + ".graph";
        contigfile = prefix + "-contig.fa";
        scaffile = prefix + "-contig-mate.fa";
        align_file = prefix + ".align";
        kmer_file = prefix + ".kmer";
        long_path_file = prefix + ".long";
        component_file = prefix + ".component";
        consensus_file = prefix + ".consensus";
        alignment_file = prefix + ".alignment";
    }
};

IDBAOption option;
AssemblyUtility assembly_utitlity;
OptionSet option_set;
//ContigGraph *contig_graph;

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
    //assembly_utitlity.AddInternalKmer(*hash_graph);

    hash_graph->RefreshEdges();
    hash_graph->RemoveDeadEnd(2*option.maxk);

    if (option.cover != 0)
        hash_graph->RemoveLowCoverageContigs(option.cover);
    else
        hash_graph->RemoveLowCoverageContigs(sqrt(hash_graph->AverageCoverage()));
    LogMessage("total kmer %lld\n", hash_graph->NumNodes());

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
}

void AlignReads()
{
    kmerLength = option.maxk;
    assembly_utitlity.AlignReads(option.contigfile, option.align_file);
}

void ConnectComponents()
{
    kmerLength = option.maxk;

    ComponentGraph *contig_graph = new ComponentGraph;

    contig_graph->ReadFrom(option.contigfile);

    deque<Component> components;

    ComponentGraph *scaffold_graph = new ComponentGraph();
    scaffold_graph->ReadFrom(option.contigfile);
    assembly_utitlity.BuildConnection(*scaffold_graph, option.align_file);

    scaffold_graph->SetMinPairs(option.min_pairs);
    scaffold_graph->ComputeDistance();
    scaffold_graph->BuildComponents();

    assembly_utitlity.BuildComponentConnection(*scaffold_graph, option.align_file);
    scaffold_graph->ConnectComponents(components);

    ofstream path_stream(option.alignment_file.c_str());
    FastaWriter fconsensus(option.consensus_file);
    FastaWriter flong_path(option.long_path_file);
    FastaWriter fcomponent(option.component_file);

    int max_num_buffer = 0;
    vector<int> num_buffers;
    vector<int64> path_length;
    vector<int64> components_length;
    for (unsigned i = 0; i < components.size(); ++i)
    {
        components_length.push_back(components[i].GetSize());

        ContigNodeAdapter begin = components[i].GetBeginContig();
        ContigNodeAdapter end = components[i].GetEndContig();

        if (!begin || !end)
            continue;

        for (unsigned j = 0; j < components[i].alignment.size(); ++j)
            path_stream << components[i].alignment[j] << endl;
        string tmp;
        tmp.append(components[i].alignment.back().size(), '-');
        path_stream << tmp << endl;

        vector<string> buffer_strings = components[i].alignment;
        if (max_num_buffer < (int)buffer_strings.size())
            max_num_buffer = buffer_strings.size();

        Contig contig = components[i].longest_path;
        //contig = components[i].consensus;

        if (contig.Size() >= 100)
            path_length.push_back(contig.Size());

        if (contig.Size() >= 100)
            flong_path.Write(contig, FormatString("path%d", i));

        if (components[i].consensus.Size() >= 100)
            fconsensus.Write(components[i].consensus, FormatString("consensus%d", i));
        deque<ContigNodeAdapter> contigs = components[i].GetContigs();
        for (unsigned j = 0; j < contigs.size(); ++j)
        {
            Contig contig;
            contigs[j].GetContig(contig);
            fcomponent.Write(contig, FormatString("component%d_%d", i, j));
        }
    }

    cout << "fine" << endl;
    cout << FormatString("n50 %lld %lld %lld %lld ", 
            Nxx(components_length, 0.5), Sum(components_length),  
            Nxx(path_length, 0.5), Sum(path_length));
    cout << FormatString("max_buf %d ", max_num_buffer);
    cout << endl;
}

void Usage()
{
    cout << "Meta-IDBA: Iterative De Bruijn graph short read Assembler for metagenomic based on graph partition" << endl;
    cout << "Version " << VERSION << endl;
    cout << endl;
    cout << "Usage: metaidba --read read-file [--output out] [options]\n" << endl;
}

int main(int argc, char *argv[])
{
    bool is_help = false;
    option_set.AddOption("help", "h", is_help, "produce help message");
    option_set.AddOption("read", "r", option.readfile, "read file");
    option_set.AddOption("long", "l", option.long_readfile, "long read file");
    option_set.AddOption("output", "o", option.prefix, "prefix of output");
    //option_set.AddOption("scaffold", "", option.is_scaffold, "use pair end information to merge contigs");
    option_set.AddOption("mink", "", option.mink, "minimum k value");
    option_set.AddOption("maxk", "", option.maxk, "maximum k value");
    option_set.AddOption("minCount", "", option.min_count, "filtering threshold for each k-mer");
    option_set.AddOption("cover", "", option.cover, "the cutting coverage for contigs");
    //option_set.AddOption("num", "n", option.num_species, "the number of species");
    //option_set.AddOption("times", "t", option.times, "the number of times of branches splitting");
    option_set.AddOption("connect", "", option.is_connect, "use paired-end reads to connect components");
    option_set.AddOption("minPairs", "", option.min_pairs, "minimum number of pair-end connections to join two components");
    option_set.AddOption("prefixLength", "", option.prefix_length, "length of the prefix of k-mer used to split k-mer table");
    //option_set.AddOption("validate", "", option.is_validate, "validate the split result");

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
        exit(1);
    }

    if (is_help || option.readfile == "")
    {
        Usage();
        cout << "Allowed Options: " << endl;
        cout << option_set.ToString() << endl;
        exit(1);
    }

    option.Compute();
    kmerLength = option.mink;

    assembly_utitlity.SetReadFile(option.readfile);
    assembly_utitlity.SetLongReadFile(option.long_readfile);

    BuildKmer();
    BuildSimpleGraph();
    Iterate();
    AlignReads();
    ConnectComponents();


    return 0;
}

