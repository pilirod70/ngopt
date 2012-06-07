/**
 * @file NodeImpl.h
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-29
 */

#ifndef __NODE_IMPL_H_

#define __NODE_IMPL_H_

#include "globals.h"
#include "BitOperation.h"

#include <algorithm>
#include <iostream>


template <class T>
class NodeImpl
{
public:
    typedef T ContentType;

    NodeImpl() { Clear(); }

    ContentType &GetContent() { return content; }
    const ContentType &GetContent() const { return content; }

    unsigned Count() const { return count; }
    void SetCount(unsigned count) { this->count = count; }

    void Increase(int c = 1)
    {
        if (count + c < 65535)
#pragma omp atomic
            count += c;
        else
            count = 65535;
    }
    void Decrease(int c = 1)
    {
        if (count - c >= 0)
#pragma omp atomic
            count -= c;
        count = 0;
    }

    unsigned InEdges() const 
    { return (edges & EdgeMask); }
    unsigned OutEdges() const 
    { return ((edges >> 4) & EdgeMask); }

    void SetInEdges(unsigned in) 
    { edges = (edges & ~EdgeMask) | (in & EdgeMask); }
    void SetOutEdges(unsigned out) 
    { edges = (edges & EdgeMask) | ((out & EdgeMask) << 4); }

    void AddInEdge(unsigned x) 
    { 
#pragma omp atomic
        edges |= (1 << x);
    }
    void AddOutEdge(unsigned x) 
    { 
#pragma omp atomic
        edges |= (1 << (x+4));
    }

    void RemoveInEdge(unsigned x) 
    { 
#pragma omp atomic
        edges &= ~(1 << x);
    }
    void RemoveOutEdge(unsigned x) 
    { 
#pragma omp atomic
        edges &= ~(1 << (x+4));
    }

    int InDegree() const 
    { return BitOperation::BitCount(edges & EdgeMask); }
    int OutDegree() const 
    { return BitOperation::BitCount((edges >> 4) & EdgeMask); }

    void Clear() { edges = 0; status = 0; count = 0; }
    void ClearStatus() { status = 0; }
    void ClearCount() { count = 0; }
    void ClearLock() { ResetUsedFlag(); status = status & ~(FlagMask); }

    void SetUsedFlag() { SetStatus(FlagUsed); }
    void ResetUsedFlag() { ResetStatus(FlagUsed);}
    bool IsUsed() const { return GetStatus(FlagUsed); }

    void SetDeadFlag() { SetStatus(FlagDead); }
    void ResetDeadFlag() { ResetStatus(FlagDead); }
    bool IsDead() const { return GetStatus(FlagDead); }

    int GetLockID()
    {
        if (!(status & FlagUsed))
            return -1;
        return status & FlagMask;
    }

    bool SetLockID(int old_id, int new_id)
    {
        unsigned char old_status = status;
        if (GetLockID() != old_id)
            return false;

        if (new_id == -1)
            return __sync_bool_compare_and_swap(&status, old_status, 
                    (old_status & ~FlagMask & ~FlagUsed));
        else
            return __sync_bool_compare_and_swap(&status, old_status, 
                    (old_status & ~FlagMask) | FlagUsed | new_id);
    }

    bool Lock(int id)
    {
        if (GetLockID() != -1)
           return false;
        return SetLockID(-1, id);
    }

    bool LockPreempt(int id)
    {
        while (true)
        {
            int old_id = GetLockID();

            if (old_id >= id)
                return false;

            if (SetLockID(old_id, id))
                return true;
        }
    }

    const unsigned &Data() const { return data; }
    unsigned &Data() { return data; }

    int &ID() { return mem[0]; }
    int &Offset() { return mem[1]; }

protected:
//    void ReverseComplement()
//    { edges = ((edges & EdgeMask) << 4) | ((edges & ~EdgeMask) >> 4); }
    
//    void Swap(NodeImpl &node)
//    {
//        std::swap(count, node.count);
//        std::swap(edges, node.edges);
//        std::swap(status, node.status);
//        std::swap(data,  node.data);
//        std::swap(content, node.content);
//    }


private:
    bool GetStatus(unsigned flag) const
    { return status & flag; }

    void SetStatus(unsigned flag)
    { 
#pragma omp atomic
        status |= flag; 
    }
    
    void ResetStatus(unsigned flag)
    { 
#pragma omp atomic
        status &= ~flag; 
    }

    void ToggleStatus(unsigned flag)
    { 
#pragma omp atomic
        status ^= flag; 
    }

    static const unsigned EdgeMask = 0x000f;
    static const unsigned char FlagDead = 0x80;
    static const unsigned char FlagUsed = 0x40;
    static const unsigned char FlagMask = 0x3F;

    int mem[0];
    unsigned short count;
    unsigned char edges;
    unsigned char status;
    unsigned data;
    ContentType content;
};

template <class N>
class NodeAdapterImpl
{
public:
    typedef N NodeType;

    NodeAdapterImpl(NodeType *node = NULL, bool is_reverse = false)
    { SetNode(node, is_reverse); }

    NodeType *GetNode() { return *node; }
    const NodeType *GetNode() const { return *node; }

    void SetNode(NodeType *node, bool is_reverse = false)
    { this->node = node; this->is_reverse = is_reverse; }

    bool IsNull() const { return node == 0; }
    bool IsReverse() const { return is_reverse; }
    void ReverseComplement() { is_reverse = !is_reverse; }

    bool operator !() const { return node == 0; }

    bool operator == (const NodeAdapterImpl &adapter) const 
    { return node == adapter.node && is_reverse == adapter.is_reverse; }
    bool operator != (const NodeAdapterImpl &adapter) const
    { return node != adapter.node || is_reverse != adapter.is_reverse; }

    bool operator < (const NodeAdapterImpl &adapter) const
    { return (node != adapter.node) ? (node < adapter.node) : (is_reverse < adapter.is_reverse); }

    unsigned char InEdges() const
    { return (!is_reverse ? node->InEdges() : node->OutEdges()); }
    unsigned char OutEdges() const 
    { return (!is_reverse ? node->OutEdges() : node->InEdges()); }

    void SetInEdges(unsigned in) 
    { !is_reverse ? node->SetInEdges(in) : node->SetOutEdges(in); }
    void SetOutEdges(unsigned out) 
    { !is_reverse ? node->SetOutEdges(out) : node->SetInEdges(out); }

    void AddInEdge(unsigned x) 
    { !is_reverse ? node->AddInEdge(x): node->AddOutEdge(x); }
    void AddOutEdge(unsigned x) 
    { !is_reverse ? node->AddOutEdge(x) : node->AddInEdge(x); }

    void RemoveInEdge(unsigned x) 
    { !is_reverse ? node->RemoveInEdge(x) : node->RemoveOutEdge(x); }
    void RemoveOutEdge(unsigned x) 
    { !is_reverse ? node->RemoveOutEdge(x) : node->RemoveInEdge(x); }

    int InDegree() const 
    { return !is_reverse ? node->InDegree() : node->OutDegree(); }
    int OutDegree() const 
    { return !is_reverse ? node->OutDegree() : node->InDegree(); }

    void SetUsedFlag() { node->SetUsedFlag(); }
    void ResetUsedFlag() { node->ResetUsedFlag(); }
    bool IsUsed() const { return node->IsUsed(); }

    void SetDeadFlag() { node->SetDeadFlag(); }
    void ResetDeadFlag() { node->ResetDeadFlag(); }
    bool IsDead() const { return node->IsDead(); }

    int Count() const { return node->Count(); }
    void SetCount(unsigned count) { node->SetCount(count); }
    void Increase(int c = 1) { node->Increase(c); }
    void Decrease(int c = 1) { node->Decrease(c); }

    int GetLockID() { return node->GetLockID(); }
    bool SetLockID(int old_id, int id) { return node->SetLockID(old_id, id); }
    bool Lock(int id) { return node->Lock(id); }
    bool LockPreempt(int id) { return node->LockPreempt(id); }

    const unsigned &Data() const { return node->Data(); }
    unsigned &Data() { return node->Data(); }

protected:
    NodeType *node;
    bool is_reverse;
};

#endif

