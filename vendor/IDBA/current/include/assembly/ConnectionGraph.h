/**
 * @file ConnectionGraph.h
 * @brief ConnectionGraph Class containing pair-end information
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __CONNECTION_GRAPH_H_

#define __CONNECTION_GRAPH_H_

#include "globals.h"
#include "ContigNode.h"
#include "ContigGraph.h"
#include "ContigBranchGroup.h"
//#include "Vector.h"
#include "HashAlign.h"

#include <omp.h>

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <deque>

class ContigConnection
{
public:
    ContigNodeAdapter to;
    int weight;
    int distance;
    std::vector<int> values;

    ContigConnection() { }

    ContigConnection(const ContigNodeAdapter &node, int value)
    { to = node; values.push_back(value); weight = value; }

    bool operator <(const ContigConnection &connection) const
    { return distance < connection.distance; }

    bool AddValue(const ContigNodeAdapter &node, int value);
    void ComputeDistance();
    void Translate(int offset);
    bool IsConsistant(int delta);
    void RemoveInconsistantPairs(int delta);
};

class ConnectionGraph: public ContigGraph
{
public:
    ConnectionGraph(int min_pairs = 5, int max_paths = 1) 
    { SetMinPairs(min_pairs); SetMaxPaths(max_paths); }
    ConnectionGraph(std::deque<Contig> &contigs, int min_pairs = 5, int max_paths = 1) 
    { Initialize(contigs); SetMinPairs(min_pairs); SetMaxPaths(max_paths); }
    virtual ~ConnectionGraph() {}

    void Initialize(std::deque<Contig> &contigs);

    void ReadFrom(const std::string &filename)
    {
        ContigGraph::ReadFrom(filename);
        std::fill_n(hits, MaxDistance, 0);
        in_connections.resize(nodes.size());
        out_connections.resize(nodes.size());

        for (int i = 0; i < MaxThreads; ++i)
            prev_tables[i] = new std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> >[MaxDistance];
    }

    void AddPair(std::vector<Alignment> &alignments1, std::vector<Alignment> &alignments2);
    void AddPair(Alignment a1, Alignment a2);
    void ComputeDistance();
    void ProcessConnections();

    void SetMinPairs(int min_pairs) { this->min_pairs = min_pairs; }
    void SetMaxPaths(int max_paths) { this->max_paths = max_paths; }
    //Vector<ContigNode> &GetNodes() { return nodes; }

    std::vector<ContigConnection> &GetConnections(const ContigNodeAdapter &node)
    { return (!node.IsReverse() ? out_connections[node.Data()] : in_connections[node.Data()]); }

    const std::vector<ContigConnection> &GetConnections(const ContigNodeAdapter &node) const
    { return (!node.IsReverse() ? out_connections[node.Data()] : in_connections[node.Data()]); }


    void AddConnection(ContigNodeAdapter node1, ContigNodeAdapter node2, int distance);
    void AddConnectionWeight(ContigNodeAdapter node1, ContigNodeAdapter node2, int weight);

    void AddConnection(std::vector<ContigConnection> &connections, 
            ContigNodeAdapter &node, int distance);
    void AddConnectionWeight(std::vector<ContigConnection> &connections, 
            ContigNodeAdapter &node, int weight);

    std::vector<std::vector<ContigConnection> > out_connections;
    std::vector<std::vector<ContigConnection> > in_connections;

protected:
    ConnectionGraph(const ConnectionGraph &);
    const ConnectionGraph &operator =(const ConnectionGraph &);

    void ProcessConnections(std::vector<ContigConnection> &connections);

    bool FindPath(ContigNodeAdapter &node, ContigNodeAdapter &target, ContigPath &path, int length);
    bool ConstructPath(ContigNodeAdapter &node, int length,
            std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > *prev,
            ContigPath &current_path,
            std::vector<ContigNodeAdapter> &partial_path,
            std::vector<ContigPath> &result_paths
            );
    
    static const int MaxDistance = 10000;
    static const int TimeLimit = 5000;

    std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > *prev_tables[MaxThreads];

    int64 hits[MaxDistance];
    int estimate_distance;
    int delta;
    int min_pairs;
    int max_paths;
};

#endif

