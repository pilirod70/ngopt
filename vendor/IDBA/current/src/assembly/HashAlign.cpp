/**
 * @file HashAlign.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Contig.h"
#include "HashAlign.h"
#include "Utils.h"

#include <omp.h>

#include <iostream>
#include <algorithm>
#include <cstdio>

using namespace std;


void HashAlign::Initialize(deque<Sequence> &contigs)
{
    int64 total = 0;
    for (int i = 0; i < (int)contigs.size(); ++i)
        total += contigs[i].Size() - kmerLength + 1;
    Reserve(total * 1.1);

#pragma omp parallel for schedule(static, 10)
    for (int i = 0; i < (int)contigs.size(); ++i)
        InsertSequence(contigs[i]);

#pragma omp parallel for schedule(static, 10)
    for (int i = 0; i < (int)contigs.size(); ++i)
        AddSequenceInfo(&contigs[i], i);

    this->contigs.resize(contigs.size());
    for (int i = 0; i < contigs.size(); ++i)
        (Sequence&)this->contigs[i] = contigs[i];
}

void HashAlign::Initialize(vector<Contig> &contigs)
{
    int64 total = 0;
    for (int i = 0; i < (int)contigs.size(); ++i)
    {
        total += contigs[i].Size() - kmerLength + 1;
    }
    Reserve(total * 1.1);

    cout << int64(total * 1.2) << endl;
    //positions.reserve(total);

//#pragma omp parallel for
//    for (int i = 0; i < (int)contigs.size(); ++i)
//        InsertSequence((Sequence &)contigs[i]);

    //positions.resize(NumNodes());
#pragma omp parallel for schedule(static, 10)
    for (int i = 0; i < (int)contigs.size(); ++i)
    {
        InsertSequence(contigs[i]);
    }

#pragma omp parallel for schedule(static, 10)
    for (int i = 0; i < (int)contigs.size(); ++i)
    {
        //InsertSequence(&contigs[i], i);
        AddSequenceInfo(&contigs[i], i);
    }

    this->contigs.swap(contigs);
}

void HashAlign::Initialize(vector<ContigNode> &contigs)
{
    int64 total = 0;
    for (int i = 0; i < (int)contigs.size(); ++i)
        total += contigs[i].GetSize();
    positions.reserve(total);

//#pragma omp parallel for
    for (int i = 0; i < (int)contigs.size(); ++i)
        InsertSequence(&contigs[i].GetContig(), i);
}

void HashAlign::InsertSequence(const Sequence &seq)
{
    if (seq.Size() < kmerLength)
        return;

    Kmer kmer;
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight(seq[i]);
    for (int i = kmerLength-1; i < seq.Size(); ++i)
    {
        kmer.AddRight(seq[i]);
        KmerNode *p = InsertKmer(kmer);
        p->Data() = NumNodes() - 1;
    }
}

void HashAlign::InsertSequence(const Sequence *seq, int id)
{
    if (seq->Size() < kmerLength)
        return;

    Kmer kmer;
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight((*seq)[i]);
    for (int i = kmerLength-1; i < seq->Size(); ++i)
    {
        kmer.AddRight((*seq)[i]);

//        if (kmer.IsPalindrome())
//            continue;

        Kmer key = kmer;
        Kmer revComp = kmer;
        revComp.ReverseComplement();

        Position pos;
        pos.id = id;
        pos.offset = i - kmerLength + 1;
        pos.length = seq->Size();
        pos.isReverse = false;

        if (revComp < kmer)
        {
            key = revComp;
            pos.ReverseComplement();
        }

        KmerNode *p = InsertKmer(key);
        int index = positions.size();
        p->Data() = index;
        //positions[p->Data()] = pos;
        positions.push_back(pos);
        positions[index] = pos;
    }
}

void HashAlign::AddSequenceInfo(const Sequence *seq, int id)
{
    if (seq->Size() < kmerLength)
        return;

    Kmer kmer;
    for (int i = 0; i < kmerLength-1; ++i)
        kmer.AddRight((*seq)[i]);
    for (int i = kmerLength-1; i < seq->Size(); ++i)
    {
        kmer.AddRight((*seq)[i]);

//        if (kmer.IsPalindrome())
//            continue;

        Kmer key = kmer;
        Kmer revComp = kmer;
        revComp.ReverseComplement();

        Position pos;
        pos.id = id;
        pos.offset = i - kmerLength + 1;
        pos.length = seq->Size();
        pos.isReverse = false;

        if (revComp < kmer)
        {
            key = revComp;
            pos.ReverseComplement();
        }

        KmerNode *p = GetNode(key);
        p->ID() = (pos.id << 1) + pos.isReverse;
        p->Offset() = pos.offset;
    }
}

void HashAlign::BuildTable(vector<Contig> &contigs)
{
}

void HashAlign::BuildPositions(vector<Contig> &contigs)
{
}

bool HashAlign::AlignKmer(const Kmer &kmer, Position &pos)
{
    KmerNode *p = GetNode(kmer);
    if (p == NULL)
        return false;

    pos.id = p->ID() >> 1;
    pos.offset = p->Offset();
    pos.length = contigs[pos.id].Size();
    pos.isReverse = p->ID() & 1;
    //pos = positions[p->Data()];
    if (kmer != p->GetKmer())
        pos.ReverseComplement();
    
    return true;
}

void HashAlign::AlignSequence(const Sequence *seq, vector<Alignment> &alignments)
{
    alignments.resize(0);

    if (seq->Size() < kmerLength)
        return ;
    {
        vector<Alignment> condidates;

        Kmer kmer;
        for (int i = 0; i < kmerLength-1; ++i)
            kmer.AddRight((*seq)[i]);
        for (int i = kmerLength-1; i < seq->Size(); ++i)
        {
            kmer.AddRight((*seq)[i]);

            Position pos;
            bool aligned = AlignKmer(kmer, pos);

            if (!aligned)
                continue;

            {
                Kmer revComp = kmer;
                revComp.ReverseComplement();

                bool flag = false;
                for (unsigned k = 0 ; k < condidates.size(); ++k)
                {
                    Alignment &align = condidates[k];
                    if (align.contigId == pos.id && align.isReverse == pos.isReverse 
                        && align.contigOffset + align.length == pos.offset + kmerLength - 1)
                    {
                        align.length += 1;
                        flag = true;
                        break;
                    }
                }

                if (!flag)
                {
                    Alignment align;
                    align.readLength = seq->Size();
                    align.readOffset = i - kmerLength + 1;
                    align.contigId = pos.id;
                    align.contigOffset = pos.offset;
                    align.contigLength = pos.length;
                    align.length = kmerLength;
                    align.isReverse = pos.isReverse;
                    condidates.push_back(align);
                }
            }

            int l = 0;
            for (unsigned k = 0; k < condidates.size(); ++k)
            {
                Alignment &align = condidates[k];
                if (i < align.readOffset + align.length)
                    condidates[l++] = align;
                else
                    alignments.push_back(align);
            }
            condidates.resize(l);
        }

        for (unsigned k = 0; k < condidates.size(); ++k)
            alignments.push_back(condidates[k]);
    }
}

