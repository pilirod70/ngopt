/**
 * @file KmerVector.h
 * @brief KmerVector Class
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __KMER_VECTOR_H_

#define __KMER_VECTOR_H_

#include "globals.h"
#include "Sequence.h"
#include "Kmer.h"

#include <algorithm>
#include <vector>

using namespace std;

class KmerVector
{
public:
    explicit KmerVector(int length = 4)
    {
        this->length = length;
        this->size = (1 << (length<<1))/2 + (~length&1) * (1 << length)/2;
        values.resize(size);
        ranks.resize(size);
    }
    ~KmerVector() { }

    void Clear()
    {
        for (int i = 0; i < size; ++i)
            values[i] = 0;
    }

    KmerVector(const KmerVector &x)
    {
        length = x.length;
        size = x.size;
        values = x.values;
        ranks = x.ranks;
        weight = x.weight;
    }

    const KmerVector &operator =(const KmerVector &x)
    {
        length = x.length;
        size = x.size;
        values = x.values;
        ranks = x.ranks;
        weight = x.weight;
        return *this;
    }

    const KmerVector &operator +=(const KmerVector &x)
    {
        for (int i = 0; i < size; ++i)
            values[i] += x.values[i];
        ComputeRank();

        return *this;
    }

    const KmerVector &operator /=(double x)
    {
        for (int i = 0; i < size; ++i)
            values[i] /= x;

        return *this;
    }

    void ComputeRank();
    void Compute(Sequence *seq);
    void Compute(std::vector<Sequence> &component);

    int length;
    int size;
    int weight;
    std::vector<double> values;
    std::vector<int> ranks;
};

double DistanceSpearman(KmerVector *kv1, KmerVector *kv2);
double DistanceKendall(KmerVector *kv1, KmerVector *kv2);

double DistanceSpearmanRef(const KmerVector &kv1, const KmerVector &kv2);

#endif

