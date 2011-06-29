/**
 * @file ContigBranchGroup.h
 * @brief ContigBranchGroup Class for managing branches in ContigGraph
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __CONTIG_BRANCH_GROUP_H_

#define __CONTIG_BRANCH_GROUP_H_

#include "globals.h"
#include "Contig.h"
#include "ContigNode.h"
#include "ContigPath.h"

#include <vector>
#include <algorithm>

class ContigGraph;

class ContigBranchGroup
{
public:
    ContigBranchGroup(ContigGraph *graph, ContigNodeAdapter &begin, int max_branches = 2, int max_length = kmerLength + 2)
    {
        this->graph = graph;
        this->begin = begin;
        this->max_branches = max_branches;
        this->max_length = max_length;
    }

    bool Search();
    void Merge();

private:
    ContigGraph *graph;
    ContigNodeAdapter begin;
    std::vector<ContigPath> branches;
    int max_branches;
    int max_length;
};

#endif

