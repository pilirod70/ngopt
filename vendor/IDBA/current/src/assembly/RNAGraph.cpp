/**
 * @file RNAGraph.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "ContigBranchGroup.h"
#include "RNAGraph.h"

#include <cmath>
#include <queue>
#include <set>
#include <deque>
#include <vector>

using namespace std;

struct CandidateBranch
{
    CandidateBranch(const ContigNodeAdapter &node, int count)
    { this->node = node; this->count = count; }

    bool operator <(const CandidateBranch &c) const { return this->count < c.count; }

    ContigNodeAdapter node;
    int count;
};

void RNAGraph::FindValidConnection(std::deque<Contig> &contigs)
{
    cout << "process connection" << endl;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ProcessConnections(in_connections[i]);
        ProcessConnections(out_connections[i]);
    }
    cout << "process connection" << endl;

    int conn_count = 0;
    int unique_count = 0;

    omp_lock_t lock;
    omp_init_lock(&lock);
    vector<ContigPath> valid_paths;
    vector<int> weights;

    cout << nodes.size() << endl;
#pragma omp parallel for schedule(static, 1)
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        ContigNodeAdapter adp(&nodes[i]);

        for (int strand = 0; strand < 2; ++strand)
        {
            if (GetConnections(adp).size() > 0)
            {
                vector<ContigConnection> &connections = GetConnections(adp);
                int index = 0;
                for (unsigned j = 0; j < connections.size(); ++j)
                {
#pragma omp atomic 
                    ++conn_count;
                    ContigPath path;
                    if (FindPath(adp, connections[j].to, path, connections[j].distance))
                    {
#pragma omp atomic 
                        ++unique_count;
                        connections[index++] = connections[j];

                        omp_set_lock(&lock);
                        valid_paths.push_back(path);
                        weights.push_back(connections[j].values.size());
                        omp_unset_lock(&lock);
                    }
                }
                connections.resize(index);
            }

            adp.ReverseComplement();
        }
    }

    contigs.resize(valid_paths.size());
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        valid_paths[i].Merge(contigs[i]);
        contigs[i].sum_coverage = weights[i];
    }

    cout << unique_count << " " << conn_count << endl;
}

void RNAGraph::FindIsoforms(std::deque<Contig> &contigs)
{
    int num_components = 0;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ContigNodeAdapter curr(&nodes[i]);

        if (curr.IsUsed())
            continue;

        ++num_components;
        queue<ContigNodeAdapter> qu;
        qu.push(curr);
        curr.SetUsedFlag();

        int num_nodes = 1;
        vector<ContigNodeAdapter> nodes;
        while (!qu.empty())
        {
            curr = qu.front();
            qu.pop();

            nodes.push_back(curr);

            for (int strand = 0; strand < 2; ++strand)
            {
                vector<ContigNodeAdapter> neighbors;
                GetNeighbors(curr, neighbors);
                for (unsigned j = 0; j < neighbors.size(); ++j)
                {
                    ContigNodeAdapter next = neighbors[j];

                    if (!next.IsUsed())
                    {
                        ++num_nodes;
                        next.SetUsedFlag();

                        if (strand == 1)
                            next.ReverseComplement();

                        qu.push(next);
                    }
                }

                curr.ReverseComplement();
            }
        }

        if (nodes.size() <= 100)
        {
            std::deque<ContigPath> paths;
            for (unsigned i = 0; i < nodes.size(); ++i)
            {
                ContigNodeAdapter curr = nodes[i];

                if (curr.InDegree() == 0)
                {
                    ContigPath parital_path;
                    set<ContigNodeAdapter> partial_set;
                    parital_path.Append(curr);
                    partial_set.insert(curr);
                    used_time = 0;
                    num_isoforms = 0;

                    FindIsoformsPaths(parital_path, partial_set, paths);
                }
            }

            if (paths.size() == 0)
                continue;

            for (int i = 0; i < paths.size(); ++i)
            {
                Contig contig;
                paths[i].Merge(contig);
                contigs.push_back(contig);
            }
        }
    }
}

void RNAGraph::FindIsoforms(ContigPath &partial_path, set<ContigNodeAdapter> &partial_set, std::deque<Contig> &contigs)
{
    ++used_time;
    ContigNodeAdapter curr = partial_path.Back();

    if (used_time > TimeLimit || num_isoforms > IsoformsLimit)
        return;
    else
    {
        vector<ContigNodeAdapter> neighbors;
        GetNeighbors(curr, neighbors);
        vector<CandidateBranch> candidates;
        for (unsigned i = 0; i < neighbors.size(); ++i)
        {
            if (partial_set.find(neighbors[i]) == partial_set.end())
            {
                ContigNodeAdapter rev_cand = neighbors[i];
                rev_cand.ReverseComplement();
                vector<ContigConnection> &connections = GetConnections(rev_cand);
                int count = 0;
                for (unsigned j = 0; j < connections.size(); ++j)
                {
                    ContigNodeAdapter tmp = connections[j].to;
                    tmp.ReverseComplement();
                    if (partial_set.find(tmp) != partial_set.end())
                    {
                        count += connections[j].values.size();
                    }
                }

                {
                    candidates.push_back(CandidateBranch(neighbors[i], count));
                }
            }
        }

        sort(candidates.begin(), candidates.end());
        reverse(candidates.begin(), candidates.end());

        for (unsigned i = 0; i < candidates.size(); ++i)
        {
            partial_path.Append(candidates[i].node);
            partial_set.insert(candidates[i].node);
            FindIsoforms(partial_path, partial_set, contigs);
            partial_path.Pop();
            partial_set.erase(partial_set.find(candidates[i].node));
        }

        if (candidates.size() == 0)
        //if (GetConnections(curr).size() == 0)
        {
            ++num_isoforms;
            Contig contig;
            partial_path.Merge(contig);
            contigs.push_back(contig);
        }
    }
}

void RNAGraph::FindIsoformsPaths(ContigPath &partial_path, set<ContigNodeAdapter> &partial_set, std::deque<ContigPath> &paths)
{
    ++used_time;
    ContigNodeAdapter curr = partial_path.Back();

    if (used_time > TimeLimit || num_isoforms > IsoformsLimit)
        return;
    else
    {
        vector<ContigNodeAdapter> neighbors;
        GetNeighbors(curr, neighbors);
        vector<CandidateBranch> candidates;
        for (unsigned i = 0; i < neighbors.size(); ++i)
        {
            if (partial_set.find(neighbors[i]) == partial_set.end())
            {
                ContigNodeAdapter rev_cand = neighbors[i];
                rev_cand.ReverseComplement();
                vector<ContigConnection> &connections = GetConnections(rev_cand);
                int count = 0;
                bool flag = 0;
                for (unsigned j = 0; j < connections.size(); ++j)
                {
                    ContigNodeAdapter tmp = connections[j].to;
                    tmp.ReverseComplement();
                    if (partial_set.find(tmp) != partial_set.end())
                    {
                        count += connections[j].weight;
                    }

                    if (tmp == curr)
                        flag = 1;
                }

                {
                    candidates.push_back(CandidateBranch(neighbors[i], count));
                }
            }
        }

        sort(candidates.begin(), candidates.end());
        reverse(candidates.begin(), candidates.end());

        for (unsigned i = 0; i < candidates.size(); ++i)
        {
            partial_path.Append(candidates[i].node);
            partial_set.insert(candidates[i].node);
            FindIsoformsPaths(partial_path, partial_set, paths);
            partial_path.Pop();
            partial_set.erase(partial_set.find(candidates[i].node));
        }

        if (candidates.size() == 0)
        {
            ++num_isoforms;
            paths.push_back(partial_path);
        }
    }
}

void RNAGraph::FindIsoformsSimple(ContigPath &partial_path, set<ContigNodeAdapter> &partial_set, std::deque<Contig> &contigs)
{
    ++used_time;
    ContigNodeAdapter curr = partial_path.Back();

    if (used_time > TimeLimit || num_isoforms > IsoformsLimit)
        return;
    else
    {
        vector<ContigNodeAdapter> neighbors;
        GetNeighbors(curr, neighbors);
        for (unsigned i = 0; i < neighbors.size(); ++i)
        {
            partial_path.Append(neighbors[i]);
            //partial_set.insert(neighbors[i]);
            FindIsoformsSimple(partial_path, partial_set, contigs);
            partial_path.Pop();
            //partial_set.erase(partial_set.find(neighbors[i]));
        }

        if (neighbors.size() == 0)
        {
            ++num_isoforms;
            Contig contig;
            partial_path.Merge(contig);
            contigs.push_back(contig);
        }
    }
}

bool RNAGraph::Check()
{
    //vector<ContigNode> &contig_nodes = GetNodes();

    int is_valid = 0;

#pragma omp parallel for
    for (int i = 0; i < nodes.size(); ++i)
    {
        ContigNodeAdapter curr(&nodes[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            vector<ContigNodeAdapter> neighbors;
            GetNeighbors(curr, neighbors);
            //vector<ContigConnection> &connections = GetConnections(curr);

            for (unsigned x = 0; x < neighbors.size(); ++x)
            {
                bool flag = false;

                ContigNodeAdapter rev_comp = neighbors[x];
                rev_comp.ReverseComplement();
                vector<ContigConnection> &connections = GetConnections(rev_comp);

                for (unsigned y = 0; y < connections.size(); ++y)
                {
                    ContigNodeAdapter tmp = connections[y].to;
                    tmp.ReverseComplement();

                    if (tmp == curr && connections[y].weight >= 5)
                        flag = true;
                }

                if (!flag)
                {
#pragma omp atomic
                    is_valid++;
                    cout << i << " " << nodes[i].GetSize() << " " << is_valid << endl;
                }
            }

            curr.ReverseComplement();
        }
    }

    return is_valid == 0;
}

