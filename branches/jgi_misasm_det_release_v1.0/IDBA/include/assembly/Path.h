/**
 * @file Path.h
 * @brief Path Class which represents a path in de Bruijn Graph
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __PATH_H_

#define __PATH_H_

#include "globals.h"
#include "Kmer.h"
#include "HashNode.h"

#include <vector>

class Path
{
public:
    Path() { weight = 0; }

    void Append(KmerNodeAdapter &next)
    {
        path.push_back(next);
        weight += UnitOne / next.Count();
    }

    KmerNodeAdapter GetEndNodeAdapter() { return path.back(); }
    int Size() { return path.size(); }
    int64 Weight() { return weight; }

    bool IsSimplePath();
    void Inactivate();
    void Activate();

private:
    static const uint64 UnitOne = 100000000ULL;

    std::vector<KmerNodeAdapter> path;
    uint64 weight;
};

#endif

