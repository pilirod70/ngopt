/**
 * @file ComponentGraph.h
 * @brief ComponentGraph class represents the relation between components in de Bruijn Graph
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-26
 */

#ifndef __COMPONENT_GRAPH_H_

#define __COMPONENT_GRAPH_H_

#include "globals.h"
#include "ContigNode.h"
#include "ContigGraph.h"
#include "ContigBranchGroup.h"
#include "HashAlign.h"
#include "ConnectionGraph.h"
#include "ComponentNode.h"

#include <omp.h>

#include <algorithm>
#include <vector>
#include <iostream>
#include <map>
#include <deque>

struct ContigAlignment
{
    Sequence seq;
    std::vector<int> pos;
};

struct QueueNode
{
    ContigNodeAdapter node;
    int distance;

    QueueNode() {}
    QueueNode(ContigNodeAdapter &n, int d): node(n), distance(d) {}

    bool operator <(const QueueNode &queue_node) const { return distance < queue_node.distance; }
};

class ComponentConnection
{
public:
    ComponentNodeAdapter to;
    int count;

    ComponentConnection() { count = 0; }
    ComponentConnection(const ComponentNodeAdapter &node) { to = node; count = 1; }

    bool operator <(const ComponentConnection &connection) const
    { return count < connection.count; }

    bool AddValue(const ComponentNodeAdapter &node)
    {
        if (to == node)
        {
            ++count;
            return true;
        }
        else
            return false;
    }
};

class ComponentGraph: public ConnectionGraph
{
public:
    explicit ComponentGraph(int min_pairs = 5) 
    { SetMinPairs(min_pairs); }
    explicit ComponentGraph(std::deque<Contig> &contigs, int min_pairs = 5) 
    { Initialize(contigs); SetMinPairs(min_pairs); }
    ~ComponentGraph() {}

    void Initialize(std::deque<Contig> &contigs);

    void ReadFrom(const std::string &filename)
    { ConnectionGraph::ReadFrom(filename); }

    void TopSortDFS(std::vector<ContigNodeAdapter> &v, ContigNodeAdapter &curr);
    void TopSort(Component &component, std::vector<QueueNode> &order, Contig &longest_path);

    void PrintComponent(Component &component);

    void RemoveCycles(ContigNodeAdapter &curr);
    void RemoveCycles();
    void BuildComponents();
    void AlignContigToLongestPath(Contig &longest, int from, int to, Contig &contig, int contig_length, ContigAlignment &contig_alignment);
    void ConnectComponents(std::deque<Component> &output_components);

    void AddComponentPair(std::vector<Alignment> &alignments1, std::vector<Alignment> &alignments2);
    void AddComponentPair(Alignment a1, Alignment a2);

    std::vector<ComponentConnection> &GetComponentConnections(const ComponentNodeAdapter &node)
    { return (!node.IsReverse() ? out_component_connections[node.Data()] : in_component_connections[node.Data()]); }

    const std::vector<ComponentConnection> &GetComponentConnections(const ComponentNodeAdapter &node) const
    { return (!node.IsReverse() ? out_component_connections[node.Data()] : in_component_connections[node.Data()]); }

    void AddComponentConnection(ComponentNodeAdapter node1, ComponentNodeAdapter node2)
    {
        AddComponentConnection(GetComponentConnections(node1), node2);
        node1.ReverseComplement();
        node2.ReverseComplement();
        AddComponentConnection(GetComponentConnections(node2), node1);
    }

    void AddComponentConnection(std::vector<ComponentConnection> &connections, ComponentNodeAdapter &node)
    {
        for (unsigned i = 0; i < connections.size(); ++i)
        {
            if (connections[i].AddValue(node))
                return;
        }

        connections.push_back(ComponentConnection(node));
    }

    std::deque<ComponentNode> &GetComponentNodes() { return components; }
    const std::deque<ComponentNode> &GetComponentNodes() const { return components; }


private:
    ComponentGraph(const ComponentGraph &);
    const ComponentGraph &operator =(const ComponentGraph &);

    std::deque<ComponentNode> components;
    std::map<ContigNodeAdapter, ComponentNodeAdapter> component_map;

    std::vector<std::vector<ComponentConnection> > out_component_connections;
    std::vector<std::vector<ComponentConnection> > in_component_connections;

    std::map<ContigNodeAdapter, int> status;
};

#endif

