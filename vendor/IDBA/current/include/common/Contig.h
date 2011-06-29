/**
 * @file Contig.h
 * @brief Contig Class
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __CONTIG_H_

#define __CONTIG_H_

#include "globals.h"
#include "Kmer.h"
#include "Sequence.h"

#include <algorithm>

class Contig: public Sequence
{
public:
    friend class HashGraph;
    friend class Reader;
    friend class Writer;
    friend class SequenceReader;
    friend class SequenceWriter;

    Contig()
    { sum_coverage = 0; is_tangle = false; }

    Contig(const Kmer &kmer): Sequence(kmer) 
    { sum_coverage = 1; is_tangle = false; }

    Contig(const Contig &contig): Sequence(contig)
    { sum_coverage = contig.sum_coverage; is_tangle = contig.is_tangle; }

    void Clear()
    { Sequence::Clear(); sum_coverage = 0; is_tangle = 0; }

    const Contig &operator =(const Contig &contig)
    { 
        if (this != &contig)
        {
            Assign(contig);
            sum_coverage = contig.sum_coverage; 
            is_tangle = contig.is_tangle;
        }
        return *this;
    }

    const Contig &operator =(const Sequence &seq)
    { 
        if (this != &seq)
        {
            Assign(seq);
            sum_coverage = 0;
            is_tangle = false;
        }
        return *this;
    }

    bool IsTangle() const { return is_tangle; }
    double Coverage() const { return sum_coverage * 1.0 / (Size() - kmerLength + 1); }

    Kmer GetBeginKmer() const { return GetKmer(0, kmerLength); }
    Kmer GetEndKmer() const { return GetKmer(Size() - kmerLength, kmerLength); }

    void Swap(Contig &contig)
    {
        Sequence::Swap(contig);
        std::swap(sum_coverage, contig.sum_coverage);
        std::swap(is_tangle, contig.is_tangle);
    }

    void Merge(const Contig &contig) 
    { 
        AddSequence(contig.ToCString() + kmerLength - 1, contig.Size() - kmerLength + 1);
        sum_coverage += contig.sum_coverage;
        is_tangle = false;
    }

    void Merge(const Contig &contig, int d)
    {
        if (d < 0)
        {
            AddSequence(contig.ToCString() - d, contig.Size() + d);
        }
        else
        {
            for (int i = 0; i < d; ++i)
                AddNucleotide((char)4);
            AddSequence(contig.ToCString(), contig.Size());
        }

        sum_coverage += contig.sum_coverage;
        is_tangle = false;
    }

    int sum_coverage;
private:
    bool is_tangle;
};

#endif

