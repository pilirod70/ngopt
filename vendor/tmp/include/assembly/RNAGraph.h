/**
 * @file RNAGraph.h
 * @brief RNAGraph Class which represents a compact de Bruijn Graph for Transcriptome Assembly
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __RNA_GRAPH_H_

#define __RNA_GRAPH_H_

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
#include <deque>
#include <map>
#include <set>


class RNAGraph: public ConnectionGraph
{
public:
    explicit RNAGraph(int min_pairs = 5) 
    { SetMinPairs(min_pairs); }
    explicit RNAGraph(std::deque<Contig> &contigs, int min_pairs = 5) 
    { Initialize(contigs); SetMinPairs(min_pairs); }
    ~RNAGraph() {}

    void Initialize(std::deque<Contig> &contigs)
    { ConnectionGraph::Initialize(contigs); }

    void FindValidConnection(std::deque<Contig> &contigs);

    void FindIsoforms(std::deque<Contig> &contigs);
    void FindIsoforms(ContigPath &partial_path, std::set<ContigNodeAdapter> &partial_set, std::deque<Contig> &contigs);
    void FindIsoformsPaths(ContigPath &partial_path, std::set<ContigNodeAdapter> &partial_set, std::deque<ContigPath> &paths);
    void FindIsoformsSimple(ContigPath &partial_path, std::set<ContigNodeAdapter> &partial_set, std::deque<Contig> &contigs);

    bool Check();

    int used_time;
    int num_isoforms;

protected:
    RNAGraph(const RNAGraph &);
    const RNAGraph &operator =(const RNAGraph &);

    static const int IsoformsLimit = 3;
};

#endif

