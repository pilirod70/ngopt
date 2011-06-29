/**
 * @file ContigPath.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "ContigBranchGroup.h"
#include "ContigGraph.h"
#include "Contig.h"
#include "ContigNode.h"

#include <vector>

using namespace std;

int64 ContigPath::InternalSize()
{
    int64 size = 1;
    for (unsigned i = 1; i < nodes.size(); ++i)
        size += nodes[i].GetSize() - (kmerLength - 1);

    if (nodes.size() > 1)
        size -= nodes.back().GetSize() - kmerLength;

    return size;
}

int64 ContigPath::InternalDistance()
{
    int64 distance = -kmerLength + 1;
    for (unsigned i = 1; i+1 < nodes.size(); ++i)
        distance += nodes[i].GetSize() - kmerLength + 1;

    return distance;
}

bool ContigPath::IsSimplePath()
{
    for (unsigned i = 1; i+1 < nodes.size(); ++i)
    {
        if (nodes[i].InDegree() != 1 || nodes[i].OutDegree() != 1)
            return false;
    }
    return true;
}

void ContigPath::Inactivate()
{
    for (unsigned i = 0; i < nodes.size(); ++i)
        nodes[i].SetDeadFlag();
}

void ContigPath::Activate()
{
    for (unsigned i = 0; i < nodes.size(); ++i)
        nodes[i].ResetDeadFlag();
}

void ContigPath::SetUsedFlag()
{
    for (unsigned i = 0; i < nodes.size(); ++i)
        nodes[i].SetUsedFlag();
}

void ContigPath::Merge(Contig &contig)
{
    contig.Clear();
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        Contig next;
        nodes[i].GetContig(next);
//        Contig next = nodes[i].GetNode()->GetContig();
//        if (nodes[i].IsReverse())
//            next.ReverseComplement();

        if (i == 0)
            contig = next;
        else
            contig.Merge(next);
    }
}

void GappedContigPath::Merge(Contig &contig)
{
    contig.Clear();
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        Contig next;
        nodes[i].GetContig(next);

        if (i == 0)
            contig = next;
        else
            contig.Merge(next, distances[i-1]);
    }
}

