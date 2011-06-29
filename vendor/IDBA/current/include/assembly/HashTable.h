/**
 * @file HashTable.h
 * @brief HashTable Class which is an efficient hash table of k-mers supporting parallel access
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __HASH_TABLE_H_

#define __HASH_TABLE_H_

#include "globals.h"
#include "KmerNode.h"
#include "HashNode.h"
#include "Sequence.h"

#include <omp.h>

#include <vector>
#include <cstdio>
#include <iostream>
#include <iterator>

class HashTable
{
public:
    typedef HashNode *HashNodePointer;

    class Iterator: public std::iterator<std::forward_iterator_tag, KmerNode *>
    {
    public:
        Iterator(HashTable *hashtable = NULL, int64 index = 0, HashNode *current = NULL);

        Iterator &operator ++();

        KmerNode *operator *() { return current; }
        bool operator == (const Iterator &iter)
        {
            return hashtable == iter.hashtable && index == iter.index && current == iter.current;
        }
        bool operator != (const Iterator &iter)
        {
            return !((*this) == iter);
        }

    private:
        void FindSolidPointer();

        HashTable *hashtable;
        int64 index;
        HashNode *current;
    };

    explicit HashTable(uint64 table_size = MaxHashTable, bool is_reversable = false);
    ~HashTable();

    Iterator Begin() { return Iterator(this, 0, 0); }
    Iterator End() { return Iterator(this, table_size, 0); }

    void Clear();
    void ClearStatus();
    void ClearCount();
    void SetDeadFlag();
    void ResetDeadFlag();
    void SetUsedFlag();
    void ResetUsedFlag();

    KmerNode *InsertKmer(const Kmer &kmer, int count = 1);
    KmerNode *GetNode(const Kmer &kmer);

    int64 NumNodes() const { return num_nodes; }
    int64 Capacity() const { return table_size; }
    void Reserve(int64 capacity) { if (capacity > table_size) Reallocate(capacity); }

    void ReadFrom(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary);
    void WriteTo(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary);

    void ReadFrom(std::istream &is);
    void WriteTo(std::ostream &os);

    static const uint64 HashNodeBufferSize = (1 << 12);
    static const uint64 MaxHashTable = 999997;
    static const uint64 LockSize = (1 << 20);

protected:
    HashNode *NewNode() { return NewNode(omp_get_thread_num()); }
    void FreeNode(HashNode *node) { FreeNode(node, omp_get_thread_num()); }

    HashNode *NewNode(int index);
    void FreeNode(HashNode *node, int index);

    HashNode **table;
    int64 table_size;
    int64 num_nodes;

private:
    HashTable(const HashTable &);
    const HashTable &operator =(const HashTable &);

    HashNode *AllocateNodes();
    void Reallocate(int64 new_table_size);

    bool is_reversable;
    int max_threads;
    omp_lock_t locks[LockSize];
    HashNode **heads;
    omp_lock_t *mem_locks;
    std::vector<HashNode *> backup;
    omp_lock_t lock_alloc;

};

#endif

