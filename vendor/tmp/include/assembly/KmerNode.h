/**
 * @file KmerNode.h
 * @brief KmerNode Class which represents an abtract k-mer in de Bruijn Graph
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __KMER_NODE_H_

#define __KMER_NODE_H_

#include "globals.h"
#include "Kmer.h"
#include "AbstractNode.h"


class KmerNode: public AbstractNode
{
public:
    KmerNode() {}
    KmerNode(const Kmer &kmer)
    { SetContent(kmer); }

    void SetContent(const Kmer &kmer)
    {
        this->kmer = kmer;
        Clear();
    }

    const Kmer &GetKmer() const { return kmer; }

protected:

    Kmer kmer;
};

class KmerNodeAdapter: public AbstractNodeAdapter
{
public:
    KmerNodeAdapter(KmerNode *node = NULL) { SetNode(node); }
    KmerNodeAdapter(KmerNode *node, const Kmer &kmer) { SetNode(node, kmer); }

    void SetNode(KmerNode *node)
    { AbstractNodeAdapter::SetNode(node); }

    void SetNode(KmerNode *node, const Kmer &kmer)
    { AbstractNodeAdapter::SetNode(node); SetDirection(kmer); }

    void SetDirection(const Kmer &kmer) 
    { if (node != 0) is_reverse = (kmer != GetNode()->GetKmer()); }

    void GetKmer(Kmer &kmer) const
    {
        kmer = GetNode()->GetKmer();
        if (is_reverse)
            kmer.ReverseComplement();
    }

    unsigned GetNucleotide(int index) const
    {
        if (!is_reverse)
            return GetNode()->GetKmer().GetBase(index);
        else
            return 3 - GetNode()->GetKmer().GetBase(kmerLength-1 - index);
    }

//    Kmer GetKmer() const 
//    { 
//        Kmer kmer = GetNode()->GetKmer();
//        if (is_reverse)
//            kmer.ReverseComplement();
//        return kmer;
//    }

    const KmerNode *GetNode() const { return (const KmerNode *)node; }
    KmerNode *GetNode() { return (KmerNode *)node; }
};

#endif

