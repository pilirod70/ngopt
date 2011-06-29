/**
 * @file AssemblyUtility.h
 * @brief Utilities for assembly
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __ASSEMBLY_UTILITY_H_

#define __ASSEMBLY_UTILITY_H_

#include "globals.h"
#include "Read.h"
#include "Sequence.h"
#include "HashGraph.h"
#include "ContigGraph.h"
#include "ScaffoldGraph.h"
#include "ComponentGraph.h"

#include <string>
#include <deque>

class AssemblyUtility
{
public:
    void BuildKmerFile(int min_count, int prefix_length, const std::string &kmer_file);
    void AddInternalKmer(HashGraph &hash_graph);
    void IterateContigGraph(ContigGraph &contig_graph, int maxk, 
            bool is_remove_deadend = true, bool is_merge_bubble = true);
    void AlignReads(const std::string &contig_file, const std::string &align_file);

    void BuildConnection(ConnectionGraph &connection_graph, const std::string &align_file);
    void BuildComponentConnection(ComponentGraph &connection_graph, const std::string &align_file);
    void BuildConnectionFromPath(ConnectionGraph &connection_graph, const std::string &path_file);

    void SetReadFile(const std::string &read_file) 
    { this->read_file = read_file; }
    void SetLongReadFile(const std::string &long_read_file)
    { this->long_read_file = long_read_file; }

    std::deque<Read> &GetReads();
    std::deque<Sequence> &GetLongReads();

    void FreeReads() { reads.clear(); }
    void FreeLongReads() { long_reads.clear(); }

private:
    std::deque<Read> reads;
    std::deque<Sequence> long_reads;
    std::string read_file;
    std::string long_read_file;
};

#endif

