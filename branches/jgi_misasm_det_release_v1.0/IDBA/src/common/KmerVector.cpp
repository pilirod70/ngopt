/**
 * @file KmerVector.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"
#include "Kmer.h"
#include "KmerVector.h"

#include <iostream>
#include <algorithm>
#include <cstdio>

using namespace std;

static long long table[1<<20];
static int aux[MaxLine];

class Comparer
{
public:
    Comparer(KmerVector *kv) { this->kv = kv; }

    bool operator ()(int x, int y) const 
    { return kv->values[x] < kv->values[y]; }

private:
    KmerVector *kv;
};

void KmerVector::ComputeRank()
{
    for (int i = 0; i < size; ++i)
        aux[i] = i;

    sort(aux, aux + size, Comparer(this));
    for (int i = 0; i < size; ++i)
        ranks[aux[i]] = i;
}

void KmerVector::Compute(Sequence *seq)
{
    int tableSize = 1 << (length << 1);
    fill_n(table, tableSize, 0);

    uint64 kmer = 0;
    for (int i = 0; i < length-1; ++i)
    {
        kmer |= (uint64((*seq)[i]) << (i << 1));
    }

    // Count the kmers
    for (int i = length-1; i < seq->Size(); ++i)
    {
        kmer |= (uint64((*seq)[i]) << ((length-1) << 1));
        ++table[kmer];
        kmer >>= 2;
    }

    // Sum up the count of the kmer and its reverseComplement,
    // and build the kmer vector.
    int k = 0;
    for (int i = 0; i < tableSize; ++i)
    {
        uint64 j = i;
        BitOperation::ReverseComplement(j);
        j >>= (64 - (length << 1));
        if (table[i] >= 0)
            values[k++] = (table[i] + table[j])*1.0/seq->Size();
        table[i] = table[j] = -1;
    }

    ComputeRank();
}

void KmerVector::Compute(vector<Sequence> &component)
{
    int tableSize = 1 << (length << 1);
    fill_n(table, tableSize, 0);

    int total_kmer = 1;
    for (unsigned index = 0; index < (int)component.size(); ++index)
    {
        Sequence *seq = &component[index];

//        if (seq->Size() < 1000)
//            continue;

        uint64 kmer = 0;
        for (int i = 0; i < length-1; ++i)
        {
            kmer |= (uint64((*seq)[i]) << (i << 1));
        }

        // Count the kmers
        for (int i = length-1; i < seq->Size(); ++i)
        {
            kmer |= (uint64((*seq)[i]) << ((length-1) << 1));
            ++table[kmer];
            kmer >>= 2;

            ++total_kmer;
        }
    }

    // Sum up the count of the kmer and its reverseComplement,
    // and build the kmer vector.
    int k = 0;
    for (int i = 0; i < tableSize; ++i)
    {
        uint64 j = i;
        BitOperation::ReverseComplement(j);
        j >>= (64 - (length << 1));
        if (table[i] >= 0)
            values[k++] = (table[i] + table[j])*1.0/total_kmer;
        table[i] = table[j] = -1;
    }

    weight = total_kmer;

    ComputeRank();
}

double DistanceSpearman(KmerVector *kv1, KmerVector *kv2)
{
    double sum = 0;
    for (int i = 0; i < kv1->size; ++i)
        sum += abs(kv1->ranks[i] - kv2->ranks[i]);
    return sum;
}

double DistanceKendall(KmerVector *kv1, KmerVector *kv2)
{
    uint64 status[kv1->size/64 + 1];
    int v[kv1->size];
    int size = kv1->size;

    fill_n(status, kv1->size/64 + 1, 0);
    for (int i = 0; i < size; ++i)
        v[kv1->ranks[i]] = kv2->ranks[i];

    int total = size * (size-1) / 2;
    int c = 0;
    for (int i = 0; i < size; ++i)
    {
        int x = v[i];
        for (int j = 0; j < (x >> 6); ++j)
            c += __builtin_popcountll(status[j]);
        c +=__builtin_popcountll(status[x>>6] & ((1ULL << (x&63)) - 1));
        status[x>>6] |= (1ULL << (x & 63));
    }

    int d = total - c;

    return (1.0 - 1.0 * (c - d) / total) / 2;
}

double DistanceSpearmanRef(const KmerVector &kv1, const KmerVector &kv2)
{
    double sum = 0;
    for (int i = 0; i < kv1.size; ++i)
        sum += abs(kv1.ranks[i] - kv2.ranks[i]);
    return sum;
}

