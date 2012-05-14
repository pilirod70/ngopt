/**
 * @file Kmer.h
 * @brief Kmer Class
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __KMER_H_

#define __KMER_H_

#include "globals.h"
#include "BitOperation.h"

#include <algorithm>
#include <iostream>
#include <cstring>

template <int MaxKValue>
class KmerImpl
{
public:
    KmerImpl(int x = kmerLength) { std::memset(data, 0, sizeof(int64) * KmerSize); }

    bool operator <(const KmerImpl<MaxKValue> &kmer) const
    {
        for (int i = KmerSize-1; i >= 0; --i)
        {
            if (data[i] != kmer.data[i])
                return data[i] < kmer.data[i];
        }
        return false;
    }

    bool operator ==(const KmerImpl<MaxKValue> &kmer) const
    {
        for (int i = 0; i < KmerSize; ++i)
        {
            if (data[i] != kmer.data[i])
                return false;
        }
        return true;
    }

    bool operator !=(const KmerImpl<MaxKValue> &kmer) const
    {
        for (int i = 0; i < KmerSize; ++i)
        {
            if (data[i] != kmer.data[i])
                return true;
        }
        return false;
    }

    void ReverseComplement()
    {
        int size = (kmerLength + 31) >> 5;

        for (int i = 0; i < size; ++i)
            BitOperation::ReverseComplement(data[i]);

        for (int i = 0; i < (size >> 1); ++i)
            std::swap(data[i], data[size-1-i]);

        if ((kmerLength & 31) != 0)
        {
            int offset = (32 - (kmerLength & 31)) << 1;
            for (int i = 0; i+1 < size; ++i)
                data[i] = (data[i] >> offset) | data[i+1] << (64 - offset);
            data[size-1] >>= offset;
        }
    }

    void AddRight(unsigned ch)
    {
        int size = (kmerLength + 31) >> 5;
        for (int i = 0; i+1 < size; ++i)
            data[i] = (data[i] >> 2) | (data[i+1] << 62);
        data[size-1] = (data[size-1] >> 2) | (uint64(ch) << (((kmerLength - 1) & 31)<< 1));
    }

    unsigned GetBase(unsigned index) const
    {
        return (data[index>>5] >> ((index & 31) << 1)) & 3;
    }

    void SetBase(unsigned index, unsigned ch)
    {
        int offset = (index & 31) << 1;
        data[index>>5] = (data[index>>5] & ~(3ULL << offset)) | (uint64(ch) << offset);
    }

    uint64 Hash() const
    {
        uint64 key = 0;
        for (int i = 0; i < KmerSize; ++i)
            key ^= data[i];
        return (key * 1299709 + 104729) % 323780508946331ULL;
    }

    bool IsPalindrome() const
    {
        KmerImpl<MaxKValue> kmer = *this;
        kmer.ReverseComplement();
        return *this == kmer;
    }

    static const int KmerSize = (MaxKValue + 31)/32;

private:
    uint64 data[KmerSize];
};

typedef KmerImpl<MaxKValue> Kmer;

#endif

