/**
 * @file Path.cpp
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
bool Path::IsSimplePath()
{
    for (unsigned i = 1; i+1 < path.size(); ++i)
    {
        if (path[i].InDegree() != 1 || path[i].OutDegree() != 1)
            return false;
    }
    return true;
}

void Path::Inactivate()
{
    for (unsigned i = 1; i+1 < path.size(); ++i)
        path[i].GetNode()->SetDeadFlag();
}

void Path::Activate()
{
    for (unsigned i = 0; i < path.size(); ++i)
        path[i].GetNode()->ResetDeadFlag();
    path[0].SetOutEdges(1 << path[1].GetNucleotide(kmerLength - 1));
    path[path.size()-1].SetInEdges(1 << (3 - path[path.size()-2].GetNucleotide(0)));
}

