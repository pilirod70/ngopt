/**
 * @file Component.h
 * @brief Component Class
 * @author Yu Peng
 * @version 0.20
 * @date 2011-01-20
 */

#ifndef __COMPNONENT_H_

#define __COMPNONENT_H_

#include "Kmer.h"
#include "Contig.h"
#include "ContigNode.h"
#include "AbstractNode.h"
#include "ContigPath.h"

#include <deque>
#include <set>
#include <algorithm>

class Component
{
public:
    Component() { }
    Component(std::set<ContigNodeAdapter> &contigs)
    { SetContent(contigs); }
    Component(std::deque<ContigNodeAdapter> &contigs)
    { SetContent(contigs); }

    void SetContent(std::set<ContigNodeAdapter> &contigs)
    { this->contigs.assign(contigs.begin(), contigs.end()); Initialize(); }
    void SetContent(std::deque<ContigNodeAdapter> &contigs)
    { this->contigs.assign(contigs.begin(), contigs.end()); Initialize(); }
    void SetContent(ContigNodeAdapter &contig)
    { this->contigs.resize(1); this->contigs[0] = contig; Initialize(); }

    void Initialize()
    {
        for (unsigned i = 0; i < contigs.size(); ++i)
        {
            ContigNodeAdapter curr = contigs[i];
            if (curr.InDegree() == 0)
            {
                if (!begin)
                    begin = curr;
                else
                {
                    begin.SetNode(NULL);
                    break;
                }
            }
        }

        for (unsigned i = 0; i < contigs.size(); ++i)
        {
            ContigNodeAdapter curr = contigs[i];
            if (curr.OutDegree() == 0)
            {
                if (!end)
                    end = curr;
                else
                {
                    end.SetNode(NULL);
                    break;
                }
            }
        }

        if (contigs.size() == 1)
            begin = end = contigs[0];
    }

    void Dispose()
    { this->contigs.clear(); }

//    void Swap(ComponentNode &node)
//    {
//        if (this != &node)
//        {
//            AbstractNode::Swap(node);
//            contigs.swap(node.contigs);
//        }
//    }

    void Add(ContigNodeAdapter &contig)
    {
        contigs.push_back(contig);
        end = contig;

        Contig contig_seq;
        contig.GetContig(contig_seq);

        PrintAlignment(contig_seq.ToString());
        longest_path.Merge(contig_seq);
        consensus.Merge(contig_seq);
        CompleteAlignment();
    }

    void Add(Component &component)
    {
        for (unsigned i = 0; i < component.contigs.size(); ++i)
            contigs.push_back(component.contigs[i]);
        end = component.end;

        int length = alignment[0].size();
        for (unsigned i = 0; i < component.alignment.size(); ++i)
            PrintAlignment(component.alignment[i], length);
        longest_path.Merge(component.longest_path);
        consensus.Merge(component.consensus);
        CompleteAlignment();
    }

    ContigNodeAdapter GetBeginContig() const
    { return begin; }
    ContigNodeAdapter GetEndContig() const 
    { return end; }

    std::deque<ContigNodeAdapter> &GetContigs() { return contigs; }
    const std::deque<ContigNodeAdapter> &GetContigs() const { return contigs; }

    int GetSize() const
    {
        int64 size = 0;
        for (unsigned i = 0; i < contigs.size(); ++i)
            size += contigs[i].GetSize();
        return size;
    }

    void ReverseComplement()
    {
        longest_path.ReverseComplement();
        consensus.ReverseComplement();
        for (unsigned i = 0; i < alignment.size(); ++i)
            reverse(alignment[i].begin(), alignment[i].end());
        for (unsigned i = 0; i < contigs.size(); ++i)
            contigs[i].ReverseComplement();
        std::swap(begin, end);
        begin.ReverseComplement();
        end.ReverseComplement();
    }

    std::string graph;
    Contig longest_path;
    std::vector<std::string> alignment;
    Contig consensus;

    void PrintAlignment(const std::string &seq, int offset = -1);
    void CompleteAlignment();

//private:
    std::deque<ContigNodeAdapter> contigs;
    ContigNodeAdapter begin;
    ContigNodeAdapter end;
};

#endif

