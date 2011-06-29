/**
 * @file BranchGroup.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "BranchGroup.h"
#include "HashGraph.h"
#include "HashNode.h"
#include "Kmer.h"

#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace std;


bool BranchGroup::Search()
{
    branches.reserve(max_branches);
    Path branch;
    branch.Append(begin);
    branches.push_back(branch);

    if (begin.InDegree() != 1 || begin.OutDegree() <= 1 || begin.OutDegree() > max_branches)
        return false;

//    if (!graph->LockNewContig(begin.GetNode()))
//        return false;

    for (int k = 1; k < max_length; ++k)
    {
        int num_branches = branches.size();
        for (int i = 0; i < num_branches; ++i)
        {
            KmerNodeAdapter adp = branches[i].GetEndNodeAdapter();

            if (adp.OutDegree() == 0)
                return false;

            bool first = true;
            Path current_branch = branches[i];
            for (int x = 0; x < 4; ++x)
            {
                if (adp.OutEdges() & (1 << x))
                {
                    Kmer kmer;
                    adp.GetKmer(kmer);
                    kmer.AddRight(x);
                    KmerNodeAdapter next = graph->GetNodeAdapter(kmer);

                    //if (next.GetNode()->GetStatus(HashGraphFlagDeadend))
                    if (next.GetNode()->IsDead())
                        return false;

//                    if (graph->GetThreadID(next.GetNode()) == omp_get_thread_num())
//                        return false;
//
//                    if (!graph->LockContig(next.GetNode()))
//                        return false;

                    if (first)
                    {
                        branches[i].Append(next);
                        first = false;
                    }
                    else
                    {
                        if ((int)branches.size() == max_branches)
                            return false;

                        Path new_branch = current_branch;
                        new_branch.Append(next);
                        branches.push_back(new_branch);
                    }
                }
            }
        }
    }

    for (unsigned i = 0; i < branches.size(); ++i)
    {
        if (!branches[i].IsSimplePath())
            return false;
    }

    KmerNodeAdapter end = branches[0].GetEndNodeAdapter();
    if (end.InDegree() != (int)branches.size() || end.OutDegree() != 1)
        return false;

    for (unsigned i = 0; i < branches.size(); ++i)
    {
        if (branches[i].GetEndNodeAdapter() != end || branches[i].Size() != max_length)
            return false;
    }

    if (begin == end)
        return false;

    return true;
}

void BranchGroup::Merge()
{
    int best = 0;
    for (unsigned i = 0; i < branches.size(); ++i)
    {
        branches[i].Inactivate();
        if (branches[i].Weight() > branches[best].Weight())
            best = i;
    }

    branches[best].Activate();
}

