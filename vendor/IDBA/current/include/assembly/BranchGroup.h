/**
 * @file BranchGroup.h
 * @brief BranchGroup Class for managing branches in HashGraph
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __BRANCH_GROUP_H_

#define __BRANCH_GROUP_H_

#include "globals.h"
#include "Kmer.h"
#include "HashNode.h"
#include "Path.h"

#include <vector>

class HashGraph;

class BranchGroup
{
public:
    BranchGroup(HashGraph *graph, KmerNodeAdapter &begin, int max_branches = 2, int max_length = kmerLength + 2)
    {
        this->graph = graph;
        this->begin = begin;
        this->max_branches = max_branches;
        this->max_length = max_length;
    }

    bool Search();
    void Merge();


private:
    HashGraph *graph;
    KmerNodeAdapter begin;
    std::vector<Path> branches;
    int max_branches;
    int max_length;
};

#endif

