/**
 * @file ScaffoldGraph.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "ContigBranchGroup.h"
#include "ScaffoldGraph.h"

#include <cmath>
#include <vector>
#include <deque>

using namespace std;

void ScaffoldGraph::Initialize(deque<Contig> &contigs)
{
    ConnectionGraph::Initialize(contigs);
    BuildAllEdges();

    in_paths.resize(nodes.size());
    out_paths.resize(nodes.size());

    double sum = 0;
    int64 count = 0;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        sum += nodes[i].GetContig().sum_coverage;
        count += nodes[i].GetContig().Size() - kmerLength + 1;
    }

    avg_coverage = sum / count;
    delta_coverage = sum / count;
    cout << "avg coverage " << avg_coverage << " delta " << delta_coverage << endl;
}

void ScaffoldGraph::FindUniquePaths()
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
    int unique_multi_degree = 0;

    cout << nodes.size() << endl;
#pragma omp parallel for schedule(static, 1)
    for (int64 i = 0; i < (int64)nodes.size(); ++i)
    {
        //cout << i << endl;
        ContigNodeAdapter adp(&nodes[i]);

        for (int strand = 0; strand < 2; ++strand)
        {
            if (GetConnections(adp).size() > 0 && IsConsistance(adp, GetConnections(adp).back().to))
                    //&& GetConnections(adp).back().to.GetSize() > estimate_distance)
            {
#pragma omp atomic 
                ++conn_count;
                ContigPath path;

                if (FindPath(adp, path))
                {
                    //omp_set_lock(&lock_path);

                    GetPaths(adp).push_back(path);
#pragma omp atomic 
                    ++unique_count;

                    if (GetConnections(adp).size() > 1)
#pragma omp atomic 
                        ++unique_multi_degree;

                    //omp_unset_lock(&lock_path);
                }

            }

            adp.ReverseComplement();
        }
    }

    cout << "process connection again" << endl;
    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        ProcessPaths(in_paths[i]);
        ProcessPaths(out_paths[i]);
    }
    cout << "process connection again" << endl;

    cout << unique_count << " " << conn_count << " " << unique_multi_degree << endl;
}

int64 ScaffoldGraph::Scaffold(std::deque<Contig> &contigs)
{
    contigs.resize(0);

    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i].IsDead() || nodes[i].IsUsed())
            continue;

        ContigNodeAdapter contig_adp(&nodes[i]);
        contig_adp.SetUsedFlag();
        ContigPath contig_path;
        contig_path.Append(contig_adp);
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                ContigNodeAdapter curr = contig_path.GetEndNodeAdapter();
                ContigPath path;

                if (GetPaths(curr).size() != 1)
                    break;

                path = GetPaths(curr)[0];

                while (true)
                {
                    ContigNodeAdapter end = path.GetEndNodeAdapter();
                    ContigNodeAdapter tmp_begin = end;
                    tmp_begin.ReverseComplement();

                    if (path.Size() > 1 &&
                            (end.IsUsed() || GetPaths(tmp_begin).size() != 1))
                        path.Pop();
                    else
                        break;
                }

                ContigPath tmp_path;
                ContigNodeAdapter tmp_end = curr;
                ContigNodeAdapter tmp_begin = path.GetEndNodeAdapter();
                tmp_begin.ReverseComplement();
                tmp_end.ReverseComplement();

                if (GetPaths(tmp_begin).size() != 1)// && GetConnections(tmp_begin).size() > 1)
                    break;

                if (path.GetEndNodeAdapter().IsUsed())
                    break;

                path.SetUsedFlag();
                contig_path.Append(path);
            }

            contig_path.ReverseComplement();
        }

        Contig contig;
        contig_path.Merge(contig);

        if (contig.Size() >= kmerLength * 3)
            contigs.push_back(contig);
    }

    return contigs.size();
}

void ScaffoldGraph::ProcessPaths(std::vector<ContigPath> &paths)
{
    int index = 0;
    for (unsigned i = 0; i < paths.size(); ++i)
    {
        bool flag = true;
        for (unsigned j = i+1; j < paths.size(); ++j)
        {
            if (paths[j].Contain(paths[i]))
            {
                flag = false;
                break;
            }
        }

        if (flag)
            paths[index++] = paths[i];
    }
    paths.resize(index);
}

bool ScaffoldGraph::FindPath(ContigNodeAdapter &node, ContigPath &path)
{
    if (GetConnections(node).size() == 0)
        return false;

    if (!IsConsistance(node, GetConnections(node).back().to))
        return false;

    //return FindPath(node, GetConnections(node).back(), path, GetConnections(node).back().distance);

    path.Append(node);

    vector<ContigConnection> &connections = GetConnections(node);
    int length = 0;
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        ContigPath middle_path;
        if (!ConnectionGraph::FindPath(path.GetEndNodeAdapter(), connections[i].to, middle_path, connections[i].distance - length))
            return false;

        length = connections[i].distance + connections[i].to.GetSize();

        path.Append(middle_path);
    }

    return true;
}

bool ScaffoldGraph::IsConsistance(const ContigNodeAdapter &node, const ContigNodeAdapter &end)
{
    return true;
//    vector<ContigConnection> &connections = GetConnections(node);
//    if (connections.back().to != end)
//        return false;
//
//    for (unsigned i = 0; i+1 < connections.size(); ++i)
//    {
//        ContigNodeAdapter middle = connections[i].to;
//        ContigNodeAdapter next = connections[i+1].to;
//
//        int d = GetDistance(node, middle) + node.GetSize() + GetDistance(middle, next);
//
//        if (abs(d - GetDistance(node, next) > delta[0]))
//        {
//            //cout << d << " " << GetDistance(node, next) << endl;
//            return false;
//        }
//    }

    return true;
}

int ScaffoldGraph::GetDistance(const ContigNodeAdapter &node, const ContigNodeAdapter &next)
{
    std::vector<ContigConnection> &connections = GetConnections(node);
    for (unsigned i = 0; i < connections.size(); ++i)
    {
        if (connections[i].to == next)
            return connections[i].distance;
    }

    return MaxDistance;
}


