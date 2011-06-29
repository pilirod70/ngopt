/**
 * @file ContigGraph.h
 * @brief ContigGraph Class which represent a compact version of de Bruijn Graph
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __CONTIG_GRAPH_H_

#define __CONTIG_GRAPH_H_

#include "globals.h"
#include "Kmer.h"
#include "Contig.h"
#include "HashGraph.h"
#include "ContigNode.h"
#include "ContigBranchGroup.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"

#include <string>
#include <set>
#include <map>
#include <deque>

class ContigGraph
{
public:
    ContigGraph() { hash_graph = new HashGraph(); }
    ContigGraph(std::deque<Contig> &contigs) 
    { hash_graph = new HashGraph(); Initialize(contigs); }
    virtual ~ContigGraph() { delete hash_graph; }

    void Clear() { hash_graph->Clear(); }
    void ClearStatus();

    virtual void Initialize(std::deque<Contig> &contigs)
    { SetContigs(contigs); BuildVertices(); }

    void Initialize(std::deque<Contig> &contigs, std::deque<Kmer> &branches)
    { SetContigs(contigs), AddBranches(branches); BuildVertices(); }

    void SortNodes()
    {
        //std::sort(nodes.Begin(), nodes.End(), Compare);
        std::sort(nodes.begin(), nodes.end(), Compare);
        BuildVertices();
        Refresh();
    }

    void ReadFrom(const std::string &filename); 
    void WriteTo(const std::string &filename);

    void ReadFrom(SequenceReader &reader);
    void WriteTo(SequenceWriter &writer);

    void BuildAllEdges() { hash_graph->AddAllEdges(); Refresh(); }
    bool BuildEdgesFromSequence(const Sequence &seq)
    { return hash_graph->AddEdgesFromSequence(seq); }
    void RemoveEdge(ContigNodeAdapter &curr, int x)
    {
        //curr.RemoveOutEdge(x);
        Kmer end = curr.GetEndKmer();
        KmerNodeAdapter end_node = hash_graph->GetNodeAdapter(end);
        end_node.RemoveOutEdge(x);

        int y = 3 - end.GetBase(0);
        Kmer next = end;
        next.AddRight(x);
        next.ReverseComplement();
        KmerNodeAdapter next_node = hash_graph->GetNodeAdapter(next);
        next_node.RemoveOutEdge(y);

        if (end.IsPalindrome())
        {
            end_node.ReverseComplement();
            end_node.RemoveOutEdge(x);
            end_node.ReverseComplement();

//            if (curr.GetSize() == kmerLength)
//            {
//                curr.ReverseComplement();
//                curr.RemoveOutEdge(x);
//                curr.ReverseComplement();
//            }
        }

        if (next.IsPalindrome())
        {
            next_node.ReverseComplement();
            next_node.RemoveOutEdge(y);
            next_node.ReverseComplement();
        }
    }

    void Refresh(int minCount = 0);
    bool Check();

    double AverageCoverage();
    int RemovePalindrome();
    int RemoveDeadEnd(int minLength);
    int RemoveStandAlone(int minLength);
    int RemoveLowCoverageContigs(double c);
    int RemoveBubble();
    void MergeContigs();

    void Decomposite();
    int64 SplitTips(int minLength);
    int64 SplitBranches();
    bool IsConverge(ContigNodeAdapter &curr);

    void SplitBranches2();
    bool IsConverge(ContigNodeAdapter &curr, std::deque<Kmer> &edges);
    void BackTraceEdges(ContigNodeAdapter &curr, std::map<ContigNodeAdapter, std::deque<ContigNodeAdapter> > &prev, std::deque<Kmer> &edges);

    int64 GetComponents(std::deque<std::set<ContigNodeAdapter> > &component_sets, std::deque<std::string> &component_strings);

    int Assemble(std::deque<Contig> &contigs);
    int Assemble(std::deque<Contig> &contigs, std::deque<Kmer> &branches);

    ContigNodeAdapter GetNeighbor(const ContigNodeAdapter &node, int x);
    void GetNeighbors(const ContigNodeAdapter &node, std::vector<ContigNodeAdapter> &neighbors);

    std::deque<ContigNode> &GetContigNodes() { return nodes; }
    int64 NumEdges() { return hash_graph->NumEdges(); }

    void ReadEdges(std::istream &is) { hash_graph->ReadFrom(is); }
    void WriteEdges(std::ostream &os) { hash_graph->WriteTo(os); }

protected:
    static bool Compare(const ContigNode &x, const ContigNode &y)
    { return x.GetSize() > y.GetSize(); }

    std::deque<ContigNode> nodes;
    //Vector<ContigNode> nodes;

private:
    ContigGraph(const ContigGraph&);
    const ContigGraph &operator =(const ContigGraph &);

    void BuildVertices();

    void SetContigs(std::deque<Contig> &contigs);
    void AddBranches(std::deque<Kmer> &branches);

    bool GetNextNodeAdapter(ContigNodeAdapter &current, ContigNodeAdapter &next);

    bool IsLoop(ContigPath &path, ContigNodeAdapter &next)
    {
        return path.GetBeginNodeAdapter().GetNode() == next.GetNode()
            || path.GetEndNodeAdapter().GetNode() == next.GetNode();
    }

    HashGraph *hash_graph;
};

#endif

