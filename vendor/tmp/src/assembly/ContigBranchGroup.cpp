/**
 * @file ContigBranchGroup.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"
#include "ContigBranchGroup.h"
#include "ContigNode.h"
#include "Contig.h"
#include "ContigGraph.h"

#include <vector>

using namespace std;

bool ContigBranchGroup::Search()
{
    ContigPath branch;
    branch.Append(begin);
    branches.push_back(branch);

    if (begin.GetSize() <= kmerLength 
            || begin.OutDegree() <= 1 || begin.OutDegree() > max_branches)
        return false;

    //for (int k = 1; k < max_length; ++k)
    while (true)
    {
        bool extend = false;

        int num_branches = branches.size();
        for (int i = 0; i < num_branches; ++i)
        {
            if (branches[i].InternalSize() >= max_length)
                continue;

            ContigNodeAdapter adp = branches[i].GetEndNodeAdapter();

            if (adp.OutDegree() == 0)
                return false;

            vector<ContigNodeAdapter> neighbors;
            graph->GetNeighbors(adp, neighbors);

            bool first = true;
            ContigPath current_branch = branches[i];
            for (unsigned j = 0; j < neighbors.size(); ++j)
            {
                ContigNodeAdapter next = neighbors[j];
                if (next.IsDead())
                    return false;

                extend = true;

                if (first)
                {
                    branches[i].Append(next);
                    first = false;
                }
                else
                {
                    if ((int)branches.size() == max_branches)
                        return false;

                    ContigPath new_branch = current_branch;
                    new_branch.Append(next);
                    branches.push_back(new_branch);
                }
            }
        }

        if (!extend)
            break;
    }

    for (unsigned i = 0; i < branches.size(); ++i)
    {
        if (!branches[i].IsSimplePath())
            return false;
    }

    ContigNodeAdapter end = branches[0].GetEndNodeAdapter();
    if (end.InDegree() != (int)branches.size() || end.GetSize() <= kmerLength)
        return false;

    for (unsigned i = 0; i < branches.size(); ++i)
    {
        if (branches[i].GetEndNodeAdapter() != end || branches[i].InternalSize() != max_length)
            return false;
    }

    if (begin == end)
        return false;

    return true;
}

void ContigBranchGroup::Merge()
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

