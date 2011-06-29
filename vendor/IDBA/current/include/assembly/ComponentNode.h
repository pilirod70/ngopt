/**
 * @file ComponentNode.h
 * @brief ComponentNode Class which is node of ComponentGraph
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-26
 */

#ifndef __COMPNONENT_NODE_H_

#define __COMPNONENT_NODE_H_

#include "Kmer.h"
#include "Contig.h"
#include "ContigNode.h"
#include "AbstractNode.h"
#include "ContigPath.h"
#include "Component.h"

#include <deque>
#include <set>
#include <algorithm>

class ComponentNode: public AbstractNode
{
public:
    ComponentNode() { Clear(); }
    ComponentNode(std::set<ContigNodeAdapter> &contigs)
    { Clear(); component.SetContent(contigs); }
    ComponentNode(std::deque<ContigNodeAdapter> &contigs)
    { Clear(); component.SetContent(contigs); }

    void SetContent(std::set<ContigNodeAdapter> &contigs)
    { Clear(); component.SetContent(contigs); }
    void SetContent(std::deque<ContigNodeAdapter> &contigs)
    { Clear(); component.SetContent(contigs); }
    void SetContent(ContigNodeAdapter &contig)
    { Clear(); component.SetContent(contig); }

    void Dispose()
    { Clear(); component.Dispose(); }

//    void Swap(ComponentNode &node)
//    {
//        if (this != &node)
//        {
//            AbstractNode::Swap(node);
//            contigs.swap(node.contigs);
//        }
//    }

    void Add(ContigNodeAdapter &contig)
    { component.Add(contig); }

    void Add(ComponentNode &node)
    { this->component.Add(node.component); }

    ContigNodeAdapter GetBeginContig() const
    { return component.GetBeginContig(); }
    ContigNodeAdapter GetEndContig() const 
    { return component.GetEndContig(); }

    std::deque<ContigNodeAdapter> &GetContigs() { return component.GetContigs(); }
    const std::deque<ContigNodeAdapter> &GetContigs() const { return component.GetContigs(); }

    int GetSize() const
    { return component.GetSize(); }

    void ReverseComplement()
    { component.ReverseComplement(); }

    void PrintAlignment(const std::string &seq, int offset = -1);
    void CompleteAlignment();

    Component component;
};

class ComponentNodeAdapter: public AbstractNodeAdapter
{
public:
    ComponentNodeAdapter(ComponentNode *node = NULL, bool is_reverse = false) 
    { SetNode(node, is_reverse); }

    void SetNode(ComponentNode *node, bool is_reverse = false)
    { AbstractNodeAdapter::SetNode(node, is_reverse); }

    const ComponentNode *GetNode() const { return (const ComponentNode*)node; }
    ComponentNode *GetNode() { return (ComponentNode*)node; }

    ContigNodeAdapter GetBeginContig() const
    {
        if (!is_reverse)
        {
            return GetNode()->GetBeginContig();
        }
        else
        {
            ContigNodeAdapter tmp = GetNode()->GetEndContig();
            tmp.ReverseComplement();
            return tmp;
        }
    }

    ContigNodeAdapter GetEndContig() const
    {
        if (!is_reverse)
        {
            return GetNode()->GetEndContig();
        }
        else
        {
            ContigNodeAdapter tmp = GetNode()->GetBeginContig();
            tmp.ReverseComplement();
            return tmp;
        }
    }

    void GetContigs(std::deque<ContigNodeAdapter> &contigs) const
    {
        contigs = GetNode()->GetContigs();

        if (is_reverse)
        {
            for (unsigned i = 0; i < contigs.size(); ++i)
                contigs[i].ReverseComplement();
        }
    }

    int GetSize() const { return GetNode()->GetSize(); }
};

#endif

