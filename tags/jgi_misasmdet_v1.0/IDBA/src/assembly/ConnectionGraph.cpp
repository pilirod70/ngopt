/**
 * @file ConnectionGraph.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"
#include "ConnectionGraph.h"
#include "ContigBranchGroup.h"

#include <vector>
#include <cmath>
#include <set>

using namespace std;

bool ContigConnection::AddValue(const ContigNodeAdapter &node, int value)
{
    if (node != to)
        return false;
    values.push_back(value);
    return true;
}

void ContigConnection::ComputeDistance()
{
    std::sort(values.begin(), values.end());
    distance = values[values.size()/2];
}

void ContigConnection::Translate(int offset)
{
    for (unsigned i = 0; i < values.size(); ++i)
        values[i] += offset;
}

bool ContigConnection::IsConsistant(int delta)
{
   int ignore = values.size()/10;
   for (unsigned i = ignore; i + ignore < values.size(); ++i)
   {
       if (abs(values[i] - distance) > delta)
           return false;
   }
   return true;
}

void ContigConnection::RemoveInconsistantPairs(int delta)
{
    int index = 0;
    for (unsigned i = 0; i < values.size(); ++i)
    {
        if (abs(values[i] - distance) < delta)
            values[index++] = values[i];
    }
    values.resize(index);
}

void ConnectionGraph::Initialize(std::deque<Contig> &contigs)
{
    ContigGraph::Initialize(contigs);

    fill_n(hits, MaxDistance, 0);
    in_connections.resize(nodes.size());
    out_connections.resize(nodes.size());

    for (int i = 0; i < MaxThreads; ++i)
        prev_tables[i] = new map<ContigNodeAdapter, std::vector<ContigNodeAdapter> >[MaxDistance];
}

void ConnectionGraph::AddPair(vector<Alignment> &alignments1, vector<Alignment> &alignments2)
{
    if (alignments1.size() > 2 && alignments2.size() > 2)
        return;

    for (unsigned i = 0; i < alignments1.size(); ++i)
    {
        for (unsigned j = 0; j < alignments2.size(); ++j)
        {
            AddPair(alignments1[i], alignments2[j]);
        }
    }
}

void ConnectionGraph::AddPair(Alignment a1, Alignment a2)
{
    a2.ReverseComplement();

    if (a1.contigId == a2.contigId && a1.isReverse == a2.isReverse)
    {
        int from = a1.contigOffset - a1.readOffset;
        int to = a2.contigOffset + (a2.readLength - a2.readOffset);
        if (!(to - from < 0 || to - from >= MaxDistance))
        {
            ++hits[to - from];
        }
    }
    else
    {
        int d = - (a1.contigLength - (a1.contigOffset - a1.readOffset))
                - (a2.contigOffset + (a2.readLength - a2.readOffset));

       AddConnection(ContigNodeAdapter(&nodes[a1.contigId], a1.isReverse),
                ContigNodeAdapter(&nodes[a2.contigId], a2.isReverse), d);
    }
}

void ConnectionGraph::ComputeDistance()
{
    int64 count = 0;
    for (int i = 0; i < MaxDistance; ++i)
        count += hits[i];

    int discard = count/200;
    int from = 0;
    int sum = 0;
    while (sum + hits[from] < discard)
    {
        sum += hits[from];
        ++from;
    }

    int to = MaxDistance - 1;
    sum = 0;
    while (sum + hits[to] < discard)
    {
        sum += hits[to];
        --to;
    }

    double sum_distance = 0;
    int real_num = 0;
    for (int i = from; i <= to; ++i)
    {
        sum_distance += i * hits[i];
        real_num += hits[i];
    }

    estimate_distance = int(round(sum_distance/real_num));

    int region = real_num * 80 / 100;
    sum = hits[estimate_distance];
    int offset = 1;
    while (offset < estimate_distance && sum  < region)
    {
        sum += hits[estimate_distance - offset] + hits[estimate_distance + offset];
        ++offset;
    }

    delta = offset;

    cout << estimate_distance << " " << delta << endl;
}

void ConnectionGraph::ProcessConnections()
{
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ProcessConnections(in_connections[i]);
        ProcessConnections(out_connections[i]);
    }
}

void ConnectionGraph::AddConnection(ContigNodeAdapter node1, ContigNodeAdapter node2, int distance)
{
    AddConnection(GetConnections(node1), node2, distance);
    node1.ReverseComplement();
    node2.ReverseComplement();
    AddConnection(GetConnections(node2), node1, distance);
}

void ConnectionGraph::AddConnectionWeight(ContigNodeAdapter node1, ContigNodeAdapter node2, int weight)
{
    AddConnectionWeight(GetConnections(node1), node2, weight);
    node1.ReverseComplement();
    node2.ReverseComplement();
    AddConnectionWeight(GetConnections(node2), node1, weight);
}


void ConnectionGraph::AddConnection(std::vector<ContigConnection> &connections, 
        ContigNodeAdapter &node, int distance)
{
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        if (connections[i].AddValue(node, distance))
            return;
    }

    connections.push_back(ContigConnection(node, distance));
}

void ConnectionGraph::AddConnectionWeight(std::vector<ContigConnection> &connections, 
        ContigNodeAdapter &node, int weight)
{
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        if (connections[i].to == node)
        {
            connections[i].weight += weight;
            return;
        }
    }

    connections.push_back(ContigConnection(node, weight));
}

void ConnectionGraph::ProcessConnections(vector<ContigConnection> &connections)
{
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        connections[i].Translate(estimate_distance);
        connections[i].ComputeDistance();
    }

    int index = 0;
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        if ((int)connections[i].values.size() >= min_pairs 
                && connections[i].distance >= -2*kmerLength 
                //&& connections[i].distance <= estimate_distance[connections[i].type] * 2
                && connections[i].IsConsistant(delta))
            connections[index++] = connections[i];
    }
    connections.resize(index);

    sort(connections.begin(), connections.end());
}

bool ConnectionGraph::FindPath(ContigNodeAdapter &node, ContigNodeAdapter &target, ContigPath &path, int length)
{
    std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > *prev;
    std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > *prev_table = prev_tables[omp_get_thread_num()];
    ContigPath current_path;
    std::vector<ContigNodeAdapter> partial_path;
    std::vector<ContigPath> result_paths;

    prev = prev_table + kmerLength;

    for (int i = -kmerLength; i <= length + delta; ++i)
        prev[i].clear();

    prev[-kmerLength][node].push_back(node);

    int used_time = 0;
    for (int i = -kmerLength; i <= length + delta; ++i)
    {
        if (i >= length)
        {
            int offset = i - length;
            if (length - offset >= -kmerLength && prev[length - offset].find(target) != prev[length - offset].end())
                break;

            if (prev[i].find(target) != prev[i].end())
                break;
        }

        used_time += prev[i].size();
        if (used_time >= TimeLimit)
        {
            break;
        }

        map<ContigNodeAdapter, vector<ContigNodeAdapter> >::iterator iter = prev[i].begin();
        while (iter != prev[i].end())
        {
            ContigNodeAdapter current = iter->first;
            //Vector<ContigNodeAdapter> neighbors;
            vector<ContigNodeAdapter> neighbors;
            GetNeighbors(current, neighbors);

            int d;
            if (current == node)
                d = - kmerLength + 1;
            else
                d = i - kmerLength + 1 + current.GetSize();

            if (d <= length + delta)
            {
                for (unsigned j = 0; j < neighbors.size(); ++j)
                    prev[d][neighbors[j]].push_back(current);
            }

            ++iter;
        }
    }

    int result_length = MaxDistance;

    for (int i = 0; i <= delta; ++i)
    {
        if (length + i >= -kmerLength && prev[length + i].find(target) != prev[length + i].end())
        {
            result_length = length + i;
            break;
        }

        if (length - i >= -kmerLength && prev[length - i].find(target) != prev[length - i].end())
        {
            result_length = length - i;
            break;
        }
    }

    if (result_length == MaxDistance)
        return false;

    result_paths.resize(0);
    if (!ConstructPath(target, result_length, prev,
                current_path, partial_path, result_paths))
        return false;
    path = result_paths[0];

    return true;
}

bool ConnectionGraph::ConstructPath(ContigNodeAdapter &target, int length, 
            std::map<ContigNodeAdapter, std::vector<ContigNodeAdapter> > *prev,
            ContigPath &current_path,
            std::vector<ContigNodeAdapter> &partial_path,
            std::vector<ContigPath> &result_paths)
{
    partial_path.push_back(target);

    if (length == -target.GetSize() && length <= -kmerLength)
    {
        reverse(partial_path.begin(), partial_path.end());
        current_path.Clear();
        for (unsigned i = 0; i < partial_path.size(); ++i)
        {
            current_path.Append(partial_path[i]);
        }
        result_paths.push_back(current_path);
        reverse(partial_path.begin(), partial_path.end());

        partial_path.pop_back();

        if ((int)result_paths.size() > max_paths)
        {
            return false;
        }

        return true;
    }

    vector<ContigNodeAdapter> &v = prev[length][target];
    for (unsigned i = 0; i < v.size(); ++i)
    {
        int d = length - (-kmerLength + 1 + v[i].GetSize());
        if (!ConstructPath(v[i], d, prev,
                    current_path, partial_path, result_paths))
        {
            partial_path.pop_back();
            return false;
        }
    }

    partial_path.pop_back();

    return true;
}

