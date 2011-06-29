/**
 * @file AssemblyUtility.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"
#include "AssemblyUtility.h"
//#include "Reader.h"
#include "SequenceReader.h"
#include "Utils.h"
#include "ContigGraph.h"
#include "HashAlign.h"
#include "Log.h"

#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>

using namespace std;

struct Connection
{
    ContigNodeAdapter from, to;
};

void AssemblyUtility::BuildKmerFile(int min_count, int prefix_length, const string &kmer_file)
{
    HashGraph *hash_graph = new HashGraph();
    ofstream fkmer(kmer_file.c_str(), ios_base::out | ios_base::binary);

    GetReads();
    GetLongReads();

    uint64 mask = (1ULL << prefix_length) - 1;
    for (int prefix = 0; prefix < (1 << prefix_length); ++prefix)
    {
#pragma omp parallel for
        for (int i = 0; i < (int64)reads.size(); ++i)
        {
            if (!reads[i].IsActive() || reads[i].Size() < kmerLength)
                continue;

            Sequence seq;
            seq = reads[i];
            hash_graph->InsertSequence(seq, prefix, mask);
        }

#pragma omp parallel for
        for (int i = 0; i < (int)long_reads.size(); ++i)
            hash_graph->InsertSequence(long_reads[i], prefix, mask);

        hash_graph->RefreshVertices(min_count);
        hash_graph->WriteTo(fkmer);
        hash_graph->Clear();
    }

    delete hash_graph;
}

void AssemblyUtility::AddInternalKmer(HashGraph &hash_graph)
{
    GetReads();

#pragma omp parallel for
    for (int i = 0; i < (int64)reads.size(); ++i)
    {
        if (!reads[i].IsActive() || reads[i].Size() < kmerLength)
            continue;

        Sequence seq;
        seq = reads[i];
        hash_graph.AddEdgesAndInternalKmers(seq);
    }

    LogMessage("after addback total kmer %lld\n", hash_graph.NumNodes());
}

void AssemblyUtility::IterateContigGraph(ContigGraph &contig_graph, int maxk,
        bool is_remove_deadend, bool is_merge_bubble)
{
//    Vector<Contig> contigs;
//    Vector<Kmer> branches;

    deque<Contig> contigs;
    deque<Kmer> branches;
    contig_graph.Assemble(contigs);

    GetReads();
    GetLongReads();

    cout << reads.size() << " " << long_reads.size() << endl;

    LogMessage("start iteration from mink = %d to maxk = %d\n", kmerLength, maxk);
    while (true)
    {
        contig_graph.Initialize(contigs, branches);

#pragma omp parallel for
        for (int64 i = 0; i < (int64)reads.size(); ++i)
        {
            if (reads[i].IsActive() || reads[i].Size() < kmerLength)
            {
                Sequence seq;
                seq = reads[i];
                if (!contig_graph.BuildEdgesFromSequence(seq))
                    reads[i].Inactivate();
            }
        }

#pragma omp parallel for
        for (int i = 0; i < (int)long_reads.size(); ++i)
            contig_graph.BuildEdgesFromSequence(long_reads[i]);

        contig_graph.Refresh();

        int minLength = kmerLength * 2;
        if (minLength < 75)
            minLength = 75;

        int deadend = 0;

        if (is_remove_deadend)
            deadend = contig_graph.RemoveDeadEnd(minLength);
        else
            deadend = contig_graph.RemoveDeadEnd(2);

        int stand_lone = contig_graph.RemoveStandAlone(minLength);

        if (kmerLength == maxk && is_merge_bubble)
        {
            int bubbles = contig_graph.RemoveBubble();
            LogMessage("remove %d bubbles\n", bubbles);
        }

        contig_graph.MergeContigs();

        int num = contig_graph.Assemble(contigs, branches);
        LogMessage("k = %d, remove %d dead end, %d stand alone, remain %d branches\n", 
                kmerLength, deadend, stand_lone, num);

        if (kmerLength == maxk)
            break;
        else
            ++kmerLength;
    }

#pragma omp parallel
    for (int64 i = 0; i < (int64)reads.size(); ++i)
        reads[i].Activate();
}

void AssemblyUtility::AlignReads(const string &contig_file, const string &align_file)
{
    GetReads();

//    Vector<Sequence> contigs;
//    FastAReader contig_reader(contig_file);
//    ((Reader &)contig_reader).Read(contigs);

    deque<Sequence> contigs;
    FastaReader contig_reader(contig_file);
    contig_reader.Read(contigs);

    HashAlign aligner;
    aligner.Initialize(contigs);
    cout << "aligner initialize" << endl;

    FILE *falign = OpenFile(align_file, "wb");

    Alignment aux_alignments[1000];
    const int64 BufferSize = 100000;
    vector<vector<Alignment> > alignments(BufferSize);
    for (unsigned offset = 0; offset < reads.size(); offset += BufferSize)
    {
        int buffer_size = min(BufferSize, (int64)reads.size() - offset);
#pragma omp parallel for
        for (int64 i = 0; i < buffer_size; ++i)
        {
            Sequence seq = reads[i + offset];
            aligner.AlignSequence(&seq, alignments[i]);
        }

        for (int i = 0; i < buffer_size; ++i)
        {
            int size = alignments[i].size();
            copy(alignments[i].begin(), alignments[i].end(), aux_alignments);
            fwrite(&size, sizeof(int), 1, falign);
            fwrite(aux_alignments, sizeof(Alignment), size, falign);
        }
    }

    fclose(falign);
}

//void AssemblyUtility::BuildScaffoldGraph(const string &contig_file, const string &align_file,
//        ScaffoldGraph &scaffold_graph)
void AssemblyUtility::BuildConnection(ConnectionGraph &connection_graph, const string &align_file)
{
    //scaffold_graph.ReadFrom(contig_file);
    Alignment aux_alignments[1000];
    FILE *falign = OpenFile(align_file, "rb");
    vector<Alignment> alignments1, alignments2;

    while (true)
    {
        int size;
        if (fread(&size, sizeof(int), 1, falign) == 0)
            break;

        alignments1.resize(size);
        fread(aux_alignments, sizeof(Alignment), size, falign);
        copy(aux_alignments, aux_alignments + size, alignments1.begin());

        if (fread(&size, sizeof(int), 1, falign) == 0)
            break;
        alignments2.resize(size);
        fread(aux_alignments, sizeof(Alignment), size, falign);
        copy(aux_alignments, aux_alignments + size, alignments2.begin());

        connection_graph.AddPair(alignments1, alignments2);
    }
}

void AssemblyUtility::BuildComponentConnection(ComponentGraph &connection_graph, const string &align_file)
{
    //scaffold_graph.ReadFrom(contig_file);
    Alignment aux_alignments[1000];
    FILE *falign = OpenFile(align_file, "rb");
    vector<Alignment> alignments1, alignments2;

    while (true)
    {
        int size;
        if (fread(&size, sizeof(int), 1, falign) == 0)
            break;

        alignments1.resize(size);
        fread(aux_alignments, sizeof(Alignment), size, falign);
        copy(aux_alignments, aux_alignments + size, alignments1.begin());

        if (fread(&size, sizeof(int), 1, falign) == 0)
            break;
        alignments2.resize(size);
        fread(aux_alignments, sizeof(Alignment), size, falign);
        copy(aux_alignments, aux_alignments + size, alignments2.begin());

        connection_graph.AddComponentPair(alignments1, alignments2);
    }
}

void AssemblyUtility::BuildConnectionFromPath(ConnectionGraph &connection_graph, const string &path_file)
{
    map<Kmer, ContigNodeAdapter> dict;
    map<Kmer, ContigNodeAdapter> in_dict;
    map<Kmer, ContigNodeAdapter> out_dict;

    vector<Sequence> sequences;
    vector<int> weights;

    FastaReader path_reader(path_file);

    //FastAReader pair_reader(path_file);
    //sequences.Resize(pair_reader.NumReads());
    //weights.Resize(sequences.Size());
    string comment;
    string tmp;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        //pair_reader.Read(sequences[i], comment);
        sequences.resize(sequences.size() + 1);
        weights.resize(weights.size() + 1);
        path_reader.Read(sequences[i], comment);
        Replace(comment, '_', ' ');
        stringstream ss(comment);
        ss >> tmp >> weights[i];
        //sequences[i].Encode();
    }
    

    //Vector<ContigNode> &contig_nodes = connection_graph.GetNodes();
    deque<ContigNode> &contig_nodes = connection_graph.GetContigNodes();
    for (unsigned i = 0; i < contig_nodes.size(); ++i)
    {
        ContigNodeAdapter curr(&contig_nodes[i]);

        for (int strand = 0; strand < 2; ++strand)
        {
            Kmer kmer = curr.GetBeginKmer();
            dict[kmer] = curr;
            in_dict[kmer] = curr;

            kmer = curr.GetEndKmer();
            dict[kmer] = curr;
            out_dict[kmer] = curr;

            curr.ReverseComplement();
        }
    }

    cout << "contigs node num " << contig_nodes.size() << endl;
    cout << "initialized dict " << dict.size() << " " << in_dict.size() << " " << out_dict.size() << endl;

    vector<vector<Connection> > paths(sequences.size());
#pragma omp parallel for schedule (static, 1)
    for (int i = 0; i < sequences.size(); ++i)
    {
        if (sequences[i].Size() <= kmerLength)
            continue;

        Kmer kmer;
        for (int j = 0; j < kmerLength-1; ++j)
            kmer.AddRight(sequences[i][j]);

        vector<ContigNodeAdapter> prev_nodes;
        for (int j = kmerLength-1; j < sequences[i].Size(); ++j)
        {
            kmer.AddRight(sequences[i][j]);

            if (in_dict.find(kmer) != in_dict.end())
            {
                for (unsigned x = 0; x < prev_nodes.size(); ++x)
                {
                    Connection conn;
                    conn.from = prev_nodes[x];
                    conn.to = in_dict[kmer];
                    paths[i].push_back(conn);
                }
            }

            if (out_dict.find(kmer) != out_dict.end())
                prev_nodes.push_back(out_dict[kmer]);
        }
    }

    int num_paths = 0;
    for (int i = 0; i < sequences.size(); ++i)
    {
        if (paths[i].size() >= 1)
            ++num_paths;
    }

    cout << "found path " << sequences.size() << " " << num_paths << endl;

    int num_connections = 0;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        for (unsigned j = 0; j < paths[i].size(); ++j)
                num_connections += 1;
    }

    cout << "connections " << num_connections << endl;

    num_connections = 0;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        for (unsigned j = 0; j < paths[i].size(); ++j)
        {
            num_connections += weights[i];
            connection_graph.AddConnectionWeight(paths[i][j].from, paths[i][j].to, weights[i]);
        }
    }

    cout << "added connection " << num_connections << endl;
}

//Vector<Read> &AssemblyUtility::GetReads()
deque<Read> &AssemblyUtility::GetReads()
{
    if (reads.size() == 0 && read_file != "")
    {
        FastaReader reader(read_file);
        reader.Read(reads);
//        Reader *reader = new FastAReader(read_file); 
//        reader->Read(reads);
//        delete reader;
    }

    return reads;
}

//Vector<Sequence> &AssemblyUtility::GetLongReads()
deque<Sequence> &AssemblyUtility::GetLongReads()
{
    if (long_reads.size() == 0 && long_read_file != "")
    {
        FastaReader reader(long_read_file);
        reader.Read(long_reads);
//        Reader *reader = new FastAReader(long_read_file); 
//        reader->Read(long_reads);
//        delete reader;
    }

    return long_reads;
}

