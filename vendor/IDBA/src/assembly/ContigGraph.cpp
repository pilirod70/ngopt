/**
 * @file ContigGraph.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Kmer.h"
#include "Sequence.h"
#include "Utils.h"
#include "HashNode.h"
#include "HashGraph.h"
#include "ContigGraph.h"
#include "ContigNode.h"
#include "Log.h"
#include "ContigBranchGroup.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <queue>
#include <sstream>
#include <map>

using namespace std;

void ContigGraph::ClearStatus()
{
#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
        nodes[i].ClearStatus();
}

void ContigGraph::ReadFrom(const string &filename)
{
    FastaReader reader(filename);
    ReadFrom(reader);
}

void ContigGraph::WriteTo(const string &filename)
{
    FastaWriter writer(filename);
    WriteTo(writer);
}

void ContigGraph::ReadFrom(SequenceReader &reader)
{
    Contig contig;
    string comment;
    string tmp;
    int in, out;

    nodes.clear();
    //nodes.resize(reader.NumReads());
    //for (int i = 0; i < nodes.size(); ++i)
    int index = 0;
    while (reader.Read(contig, comment))
    {
        //reader.Read(contig, comment);
        //contig.Encode();

        nodes.resize(nodes.size() + 1);
        nodes[index].SetContent(contig);

        if (nodes[index].GetSize() < kmerLength)
        {
            nodes[index].SetDeadFlag();
            continue;
        }

        Replace(comment, '_', ' ');
        stringstream ss(comment);
        ss >> tmp >> in >> out;

        nodes[index].SetInEdges(in);
        nodes[index].SetOutEdges(out);

        index++;
    }

    //cout << nodes.size() << endl;

    BuildVertices();
    Refresh();
}

void ContigGraph::WriteTo(SequenceWriter&writer)
{
    Contig contig;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        contig = nodes[i].GetContig();
        //contig.Decode();
        writer.Write(contig, FormatString("node%d_%d_%d", i, nodes[i].InEdges(), nodes[i].OutEdges()));
        //writer.WriteFormat(contig, "node%d_%d_%d", i, nodes[i].InEdges(), nodes[i].OutEdges());
    }
}

void ContigGraph::Refresh(int minCount)
{
    hash_graph->SetDeadFlag();
#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        nodes[i].ClearLock();

        KmerNode *p = hash_graph->GetNode(nodes[i].GetBeginKmer());
        if (p == NULL)
        {
            LogError("ContigGraph::Refresh() inconsistent\n");
            exit(1);
        }
        p->Data() = i;
        p->ResetDeadFlag();

        p = hash_graph->GetNode(nodes[i].GetEndKmer());
        if (p == NULL)
        {
            LogError("ContigGraph::Refresh() inconsistent\n");
            exit(1);
        }
        p->Data() = i;
        p->ResetDeadFlag();
    }

    hash_graph->Refresh();

#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
        {
            nodes[i].SetInEdges(0);
            nodes[i].SetOutEdges(0);
            continue;
        }

        ContigNodeAdapter contigAdp(&nodes[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            KmerNodeAdapter adp = hash_graph->GetNodeAdapter(contigAdp.GetEndKmer());
            contigAdp.SetOutEdges(adp.OutEdges());

            for (int x = 0; x < 4; ++x)
            {
                if ((1 << x) & adp.OutEdges())
                {
                    ContigNodeAdapter neighbor = GetNeighbor(contigAdp, x);
                    if (neighbor.IsNull())
                    {
                        adp.RemoveOutEdge(x);
                        contigAdp.RemoveOutEdge(x);
                    }
                }
            }

            contigAdp.ReverseComplement();
        }
    }
}

double ContigGraph::AverageCoverage()
{
    double sum = 0;
    int64 count = 0;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        sum += nodes[i].GetContig().sum_coverage;
        count += nodes[i].GetContig().Size() - kmerLength + 1;
    }
    return sum / count;
}

int ContigGraph::RemovePalindrome()
{
    int deadend = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Size() == kmerLength && nodes[i].GetContig().GetBeginKmer().IsPalindrome())
        {
            nodes[i].SetDeadFlag();

#pragma omp atomic
            ++deadend;
        }
    }

    Refresh();

    return deadend;
}

int ContigGraph::RemoveDeadEnd(int minLength)
{
    int deadend = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Size() < kmerLength + minLength - 1 && 
                (nodes[i].InDegree() == 0 || nodes[i].OutDegree() == 0))
        {
            nodes[i].SetDeadFlag();

#pragma omp atomic
            ++deadend;
        }
    }

    Refresh();

    return deadend;
}

int ContigGraph::RemoveStandAlone(int minLength)
{
    int deadend = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Size() < kmerLength + minLength - 1 && 
                (nodes[i].InDegree() == 0 && nodes[i].OutDegree() == 0))
        {
            nodes[i].SetDeadFlag();

#pragma omp atomic
            ++deadend;
        }
    }

    Refresh();

    return deadend;
}

int ContigGraph::RemoveLowCoverageContigs(double c)
{
    int remove = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Coverage() < c)
        {
            nodes[i].SetDeadFlag();

#pragma omp atomic
            ++remove;
        }
    }

    Refresh();

    return remove;
}

int ContigGraph::RemoveBubble()
{
    int bubbles = 0;
    for (unsigned k = 0; k < nodes.size(); ++k)
    {
        if (nodes[k].IsDead())
            continue;

        ContigNodeAdapter begin(&nodes[k]);
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigBranchGroup branch_group(this, begin);
            if (branch_group.Search())
            {
                branch_group.Merge();
                ++bubbles;
            }

            begin.ReverseComplement();
        }
    }
    
    Refresh();

    return bubbles;
}

void ContigGraph::MergeContigs()
{
    //Vector<Contig> new_contigs;
    deque<Contig> new_contigs;
    Assemble(new_contigs);
    SetContigs(new_contigs);
    Refresh();
}

void ContigGraph::Decomposite()
{
    int64 last = 0;
    for (int i = 0; i < 100; ++i)
    {
        int64 split = SplitBranches();
        cout << split << " " << 2*GetContigNodes().size() << endl;
//        int64 deadend = SplitTips(kmerLength * 2);
//        cout << deadend << " " << endl;

        if (last == split)
            break;
        last = split;
    }
}

int64 ContigGraph::SplitTips(int minLength)
{
    int deadend = 0;
#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (nodes[i].GetContig().Size() < kmerLength + minLength - 1 && 
                ((nodes[i].InDegree() == 0 && nodes[i].OutDegree() == 1)
                 || (nodes[i].OutDegree() == 0 && nodes[i].InDegree() == 1)))
        {
            Kmer begin = nodes[i].GetBeginKmer();
            Kmer end = nodes[i].GetEndKmer();

            hash_graph->GetNode(begin)->SetDeadFlag();
            hash_graph->GetNode(end)->SetDeadFlag();

#pragma omp atomic
            ++deadend;
        }
    }

    hash_graph->RefreshEdges();
    hash_graph->ResetDeadFlag();

    Refresh();

    return deadend;
}

int64 ContigGraph::SplitBranches()
{
    int count = 0;
#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
        {
#pragma omp atomic
            count += 2;
            continue; 
        }

        ContigNodeAdapter curr(&nodes[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            if (!IsConverge(curr))
            {
#pragma omp atomic
                ++count;

                set<ContigNodeAdapter> sources;
                queue<ContigNodeAdapter> qu;
                qu.push(curr);
                sources.insert(curr);

                while (!qu.empty())
                {
                    ContigNodeAdapter u = qu.front();
                    qu.pop();

                    for (int x = 0; x < 4; ++x)
                    {
                        if (u.OutEdges() & (1 << x))
                        {
                            RemoveEdge(u, x);

                            ContigNodeAdapter v = GetNeighbor(u, x);
                            v.ReverseComplement();
                            if (sources.find(v) == sources.end())
                            {
                                sources.insert(v);
                                qu.push(v);
                            }
                        }
                    }
                }
//                for (int x = 0; x < 4; ++x)
//                {
//                    if (curr.OutEdges() & (1 << x))
//                        RemoveEdge(curr, x);
//                }
            }

            curr.ReverseComplement();
        }
    }

    //cout << count << " " << nodes.size() * 2 << endl;

    Refresh();

    return count;
}

struct SearchNode
{
    ContigNodeAdapter node;
    int distance;
    int label;
};

bool ContigGraph::IsConverge(ContigNodeAdapter &curr)
{
    int TimeLimit = 10000;
    int DistanceLimit = 300;
    map<ContigNodeAdapter, int> reachable;
    queue<SearchNode> qu;

    for (int x = 0; x < 4; ++x)
    {
        if (curr.OutEdges() & (1 << x))
        {
            SearchNode search_node;
            search_node.node = GetNeighbor(curr, x);
            search_node.distance = -kmerLength + 1;
            search_node.label = x;

            if (!search_node.node.IsDead())
                qu.push(search_node);
        }
    }

    int time = 0;
    while (!qu.empty())
    {
        if (time++ == TimeLimit)
            break;

        SearchNode search_node = qu.front();
        qu.pop();

        reachable[search_node.node] |= (1 << search_node.label);
        
        if (reachable[search_node.node] == curr.OutEdges())
        {
            //cout << search_node.distance << endl;
            return true;
        }

        if (search_node.distance + search_node.node.GetSize() - kmerLength + 1 > DistanceLimit)
            continue;

        for (int x = 0; x < 4; ++x)
        {
            if (search_node.node.OutEdges() & (1 << x))
            {
                ContigNodeAdapter next = GetNeighbor(search_node.node, x);

                SearchNode new_search_node;
                new_search_node.node = next;
                new_search_node.distance = search_node.distance + search_node.node.GetSize() - kmerLength + 1;
                new_search_node.label = search_node.label;

//                if (new_search_node.node == curr)
//                    continue;

                if (reachable[new_search_node.node] & (1 << new_search_node.label))
                    continue;

                if (!new_search_node.node.IsDead())
                    qu.push(new_search_node);
            }
        }
    }

//    for (map<ContigNodeAdapter, int>::iterator iter = reachable.begin(); iter != reachable.end(); ++iter)
//    {
//        if (iter->second == curr.OutEdges())
//            return true;
//    }

    return false;
}

void ContigGraph::SplitBranches2()
{
    cout << hash_graph->NumEdges() << endl;
    deque<Kmer> edges;

    int count = 0;
//#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        ContigNodeAdapter curr(&nodes[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            if (!IsConverge(curr, edges))
            {
//#pragma omp atomic
//                ++count;
//                for (int x = 0; x < 4; ++x)
//                {
//                    if (curr.OutEdges() & (1 << x))
//                        RemoveEdge(curr, x);
//                }
            }

            curr.ReverseComplement();
        }
    }

    cout << count << " " << nodes.size() * 2 << " " << edges.size() << endl;

    hash_graph->ClearGraph();
    for (unsigned i = 0; i < edges.size(); ++i)
    {
        Kmer kmer = edges[i];
        int x = kmer.GetBase(kmerLength);
        kmer.SetBase(kmerLength, 0);
        hash_graph->GetNodeAdapter(kmer).AddOutEdge(x);

        int y = kmer.GetBase(0);
        kmer.AddRight(x);
        kmer.ReverseComplement();
        hash_graph->GetNodeAdapter(kmer).AddOutEdge(y);
    }
    Refresh();

    cout << hash_graph->NumEdges() << endl;
}

bool ContigGraph::IsConverge(ContigNodeAdapter &curr, deque<Kmer> &edges)
{
    //edges.resize(0);

    int TimeLimit = 10000;
    int DistanceLimit = 300;
    map<ContigNodeAdapter, int> reachable;
    map<ContigNodeAdapter, deque<ContigNodeAdapter> > prev;
    queue<SearchNode> qu;

    for (int x = 0; x < 4; ++x)
    {
        if (curr.OutEdges() & (1 << x))
        {
            SearchNode search_node;
            search_node.node = GetNeighbor(curr, x);
            search_node.distance = -kmerLength + 1;
            search_node.label = x;

            prev[search_node.node].push_back(curr);

            qu.push(search_node);
        }
    }

    int time = 0;
    while (!qu.empty())
    {
        if (time++ == TimeLimit)
            break;

        SearchNode search_node = qu.front();
        qu.pop();

        reachable[search_node.node] |= (1 << search_node.label);

        if (reachable[search_node.node] == curr.OutEdges())
        {
            BackTraceEdges(search_node.node, prev, edges);
            return true;
        }

        if (search_node.distance + search_node.node.GetSize() - kmerLength + 1 > DistanceLimit)
            continue;

        for (int x = 0; x < 4; ++x)
        {
            if (search_node.node.OutEdges() & (1 << x))
            {
                ContigNodeAdapter next = GetNeighbor(search_node.node, x);

                SearchNode new_search_node;
                new_search_node.node = next;
                new_search_node.distance = search_node.distance + search_node.node.GetSize() - kmerLength + 1;
                new_search_node.label = search_node.label;

                prev[new_search_node.node].push_back(search_node.node);

                qu.push(new_search_node);
            }
        }
    }

    return false;
}

void ContigGraph::BackTraceEdges(ContigNodeAdapter &curr, map<ContigNodeAdapter, deque<ContigNodeAdapter> > &prev, deque<Kmer> &edges)
{
    for (unsigned i = prev[curr].size()-1; i < prev[curr].size(); ++i)
    {
        ContigNodeAdapter p = prev[curr][i];

        Kmer x = p.GetEndKmer();
        Kmer y = curr.GetBeginKmer();

        for (int a = 0; a < 4; ++a)
        {
            Kmer tmp = x;
            tmp.AddRight(a);
            if (tmp == y)
            {
                x.SetBase(kmerLength, a);
                edges.push_back(x);
                break;
            }
        }
    }

    deque<ContigNodeAdapter> v;
    v.swap(prev[curr]);

    for (unsigned i = 0; i < v.size(); ++i)
        BackTraceEdges(v[i], prev, edges);
}

int64 ContigGraph::GetComponents(std::deque<std::set<ContigNodeAdapter> > &component_sets, 
        std::deque<string> &component_strings)
{
    component_sets.clear();
    component_strings.clear();

    int node_index = 0;

    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ContigNodeAdapter curr(&nodes[i]);

        if (curr.IsUsed() || curr.IsDead())
            continue;

        component_sets.push_back(set<ContigNodeAdapter>());
        stringstream ss;

        queue<ContigNodeAdapter> qu;
        qu.push(curr);
        curr.SetUsedFlag();

        curr.Data() = node_index++;

        ContigNodeAdapter begin;
        ContigNodeAdapter end;

        int num_nodes = 1;
        int64 sum_length = curr.GetSize();
        while (!qu.empty())
        {
            curr = qu.front();
            qu.pop();

            component_sets.back().insert(curr);

            for (int strand = 0; strand < 2; ++strand)
            {
                vector<ContigNodeAdapter> neighbors;
                GetNeighbors(curr, neighbors);
                for (unsigned j = 0; j < neighbors.size(); ++j)
                {
                    ContigNodeAdapter next = neighbors[j];

                    if (!next.IsUsed())
                    {
                        sum_length += next.GetSize();
                        next.SetUsedFlag();
                        next.Data() = node_index++;

                        if (strand == 1)
                            next.ReverseComplement();

                        qu.push(next);
                    }

                    if (strand == 1)
                        ss << next.Data() << " " << curr.Data() << endl;
                    else
                        ss << curr.Data() << " " << next.Data() << endl;
                }

                curr.ReverseComplement();
            }
        }

        component_strings.push_back(ss.str());
    }

    ClearStatus();
    for (unsigned i = 0; i < nodes.size(); ++i)
        nodes[i].Data() = i;

    return component_sets.size();
}

int ContigGraph::Assemble(deque<Contig> &result_contigs)
{
    result_contigs.clear();

    omp_lock_t lockContigs;
    omp_init_lock(&lockContigs);

#pragma omp parallel for
    for (int i = 0; i < (int)nodes.size(); ++i)
    {
        if (nodes[i].GetContig().Size() == kmerLength 
                && nodes[i].GetBeginKmer().IsPalindrome()
                && !nodes[i].IsDead())
        {
            nodes[i].SetDeadFlag();

            omp_set_lock(&lockContigs);
            result_contigs.push_back(nodes[i].GetContig());
            omp_unset_lock(&lockContigs);
        }
    }

    int total = 0;
#pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].IsDead())
            continue;

        if (!nodes[i].Lock(omp_get_thread_num()))
            continue;

        ++total;
        ContigNodeAdapter current(&nodes[i]);
        ContigPath path;
        path.Append(current);

        Contig contig;
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                current = path.GetEndNodeAdapter();
                ContigNodeAdapter next;

                if (!GetNextNodeAdapter(current, next))
                    break;

                if (next.IsDead())
                    break;

                if (next.GetLockID() == omp_get_thread_num() && IsLoop(path, next))
                    break;

                if (!next.LockPreempt(omp_get_thread_num()))
                    goto FAIL;

                path.Append(next);
            }

            path.ReverseComplement();
        }

        path.Merge(contig);

        omp_set_lock(&lockContigs);
        result_contigs.push_back(contig);
        omp_unset_lock(&lockContigs);
FAIL:
        ;
    }

    omp_destroy_lock(&lockContigs);

    return result_contigs.size();
}

int ContigGraph::Assemble(deque<Contig> &result_contigs, deque<Kmer> &branches)
{
    result_contigs.clear();
    branches.clear();

    omp_lock_t lockKmers;
    omp_init_lock(&lockKmers);

    int64 size = nodes.size();
#pragma omp parallel for
    for (int64 i = 0; i < size; ++i)
    {
        if (!nodes[i].IsDead() && nodes[i].GetSize() >= kmerLength)
        {
            ContigNodeAdapter contig_adp(&nodes[i]);
            for (int strand = 0; strand < 2; ++strand)
            {
                Kmer kmer = contig_adp.GetEndKmer();
                unsigned edges = contig_adp.OutEdges();
                
                for (int j = 0; j < 4; ++j)
                {
                    if (edges & (1 << j))
                    {
                        Kmer x = kmer;
                        if (x.IsPalindrome() && contig_adp.GetSize() > kmerLength 
                                && contig_adp.GetNucleotide(contig_adp.GetSize() - kmerLength - 1) == (3U - j))
                            continue;

                        x.SetBase(kmerLength, j);
                        omp_set_lock(&lockKmers);
                        branches.push_back(x);
                        omp_unset_lock(&lockKmers);
                    }
                }

                contig_adp.ReverseComplement();
            }
        }
    }

    ++kmerLength;

#pragma omp parallel for
    for (int64 i = 0; i < (int64)branches.size(); ++i)
    {
        Kmer revComp = branches[i];
        revComp.ReverseComplement();
        if (revComp < branches[i])
            branches[i] = revComp;
    }

    sort(branches.begin(), branches.end());
    branches.resize(unique(branches.begin(), branches.end()) - branches.begin());

    --kmerLength;

    return Assemble(result_contigs);
}

ContigNodeAdapter ContigGraph::GetNeighbor(const ContigNodeAdapter &node, int x)
{
    Kmer kmer = node.GetEndKmer();
    kmer.AddRight(x);

    KmerNode *p = hash_graph->GetNode(kmer);
    if (p == NULL)
        return ContigNodeAdapter(NULL);

    ContigNodeAdapter neighbor(&nodes[p->Data()]);
    if (neighbor.GetBeginKmer() != kmer)
        neighbor.ReverseComplement();

    if (neighbor.GetBeginKmer() != kmer)
        return ContigNodeAdapter(NULL);
    else
        return neighbor;
} 

void ContigGraph::GetNeighbors(const ContigNodeAdapter &node, 
        std::vector<ContigNodeAdapter> &neighbors)
{
    neighbors.resize(0);
    for (int x = 0; x < 4; ++x)
    {
        if ((1 << x) & node.OutEdges())
        {
            Kmer kmer = node.GetEndKmer();
            kmer.AddRight(x);

            int j = hash_graph->GetNode(kmer)->Data();
            ContigNodeAdapter neighbor(&nodes[j]);
            if (neighbor.GetBeginKmer() != kmer)
                neighbor.ReverseComplement();

            neighbors.push_back(neighbor);
        }
    }
}

void ContigGraph::BuildVertices()
{
    hash_graph->Clear();
    hash_graph->Reserve(nodes.size() * 2);

#pragma omp parallel for
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        if (nodes[i].GetContig().Size() < kmerLength)
        {
            nodes[i].SetDeadFlag();
            continue;
        }

        nodes[i].Data() = i;
        hash_graph->InsertKmer(nodes[i].GetBeginKmer());
        hash_graph->InsertKmer(nodes[i].GetEndKmer());

        hash_graph->GetNodeAdapter(nodes[i].GetEndKmer()).SetOutEdges(nodes[i].OutEdges());
        hash_graph->GetNodeAdapter(nodes[i].GetBeginKmer()).SetInEdges(nodes[i].InEdges());
    }

    hash_graph->Refresh();
}

void ContigGraph::SetContigs(deque<Contig> &contigs)
{
    nodes.clear();

    this->nodes.resize(contigs.size());
#pragma omp parallel for
    for (int i = 0; i < (int)contigs.size(); ++i)
    {
        this->nodes[i].SetContent(contigs[i]);
        this->nodes[i].Data() = i;
    }
}

void ContigGraph::AddBranches(deque<Kmer> &branches)
{
    int old_size = nodes.size();
    nodes.resize(nodes.size() + branches.size());
#pragma omp parallel for
    for (int i = 0; i < (int)branches.size(); ++i)
    {
        Contig contig(branches[i]);
        nodes[old_size + i].SetContent(contig);
        nodes[old_size + i].Data() = old_size + i;
    }
}

bool ContigGraph::GetNextNodeAdapter(ContigNodeAdapter &current, ContigNodeAdapter &next)
{
    if (current.OutDegree() != 1)
        return false;

    Kmer kmer = current.GetEndKmer();
    kmer.AddRight(BitOperation::bitToIndex[current.OutEdges()]);

    next.SetNode(&nodes[hash_graph->GetNode(kmer)->Data()]);
    if (next.GetBeginKmer() != kmer)
        next.ReverseComplement();

    return next.InDegree() == 1;
}

