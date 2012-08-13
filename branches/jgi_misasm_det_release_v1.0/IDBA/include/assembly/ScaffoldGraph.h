/**
 * @file ScaffoldGraph.h
 * @brief ScaffoldGraph Class which represents a compact de Bruijn Graph for scaffolding
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __SCAFFOLD_GRAPH_H_

#define __SCAFFOLD_GRAPH_H_

#include "globals.h"
#include "ContigNode.h"
#include "ContigGraph.h"
#include "ContigBranchGroup.h"
#include "HashAlign.h"
#include "ConnectionGraph.h"

#include <omp.h>

#include <algorithm>
#include <vector>
#include <iostream>
#include <map>
#include <deque>

class ScaffoldGraph: public ConnectionGraph
{
public:
    explicit ScaffoldGraph(int min_pairs = 5) 
    { SetMinPairs(min_pairs); }
    explicit ScaffoldGraph(std::deque<Contig> &contigs, int min_pairs = 5) 
    { Initialize(contigs); SetMinPairs(min_pairs); }
    ~ScaffoldGraph() {}

    void Initialize(std::deque<Contig> &contigs);

    void ReadFrom(const std::string &filename)
    {
        ConnectionGraph::ReadFrom(filename);
        in_paths.resize(nodes.size());
        out_paths.resize(nodes.size());

        double sum = 0;
        int64 count = 0;
        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            sum += nodes[i].GetContig().sum_coverage;
            count += nodes[i].GetContig().Size() - kmerLength + 1;
        }

        avg_coverage = sum / count;
        delta_coverage = sum / count;
    }

    void FindUniquePaths();
    int64 Scaffold(std::deque<Contig> &contigs);
    int64 ScaffoldWithGap(std::vector<Contig> &contigs);

    std::vector<ContigPath> &GetPaths(const ContigNodeAdapter &node)
    { return (!node.IsReverse() ? out_paths[node.Data()] : in_paths[node.Data()]); }

    const std::vector<ContigPath> &GetPaths(const ContigNodeAdapter &node) const
    { return (!node.IsReverse() ? out_paths[node.Data()] : in_paths[node.Data()]); }

private:
    ScaffoldGraph(const ScaffoldGraph &);
    const ScaffoldGraph &operator =(const ScaffoldGraph &);

    void ProcessPaths(std::vector<ContigPath> &paths);
    bool FindPath(ContigNodeAdapter &node, ContigPath &path);
    
    bool IsConsistance(const ContigNodeAdapter &node, const ContigNodeAdapter &end);
    int GetDistance(const ContigNodeAdapter &node, const ContigNodeAdapter &next);
    bool IsSeed(const ContigNodeAdapter &node) 
    { return node.GetSize() >= kmerLength * 2; }

    std::vector<std::vector<ContigPath> > out_paths;
    std::vector<std::vector<ContigPath> > in_paths;

    double avg_coverage;
    double delta_coverage;
};

#endif

