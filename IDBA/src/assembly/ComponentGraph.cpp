/**
 * @file ComponentGraph.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-26
 */

#include "globals.h"

#include "ContigBranchGroup.h"
#include "ComponentGraph.h"
#include "Utils.h"

#include <cmath>
#include <vector>
#include <deque>
#include <queue>
#include <sstream>

using namespace std;


void ComponentGraph::Initialize(deque<Contig> &contigs)
{
    ConnectionGraph::Initialize(contigs);
    BuildAllEdges();
}

void ComponentGraph::TopSortDFS(vector<ContigNodeAdapter> &v, ContigNodeAdapter &curr)
{
    if (status[curr])
        return;

    vector<ContigNodeAdapter> neighbors;
    GetNeighbors(curr, neighbors);

    status[curr] = 1;
    for (int i = 0; i < neighbors.size(); ++i)
        TopSortDFS(v, neighbors[i]);
    v.push_back(curr);
}

void ComponentGraph::TopSort(Component &component, std::vector<QueueNode> &order, Contig &longest_path)
{
    order.resize(0);

    ContigNodeAdapter begin = component.GetBeginContig();
    ContigNodeAdapter end = component.GetEndContig();
    vector<ContigNodeAdapter> v;
    ContigPath path;

    TopSortDFS(v, begin);
    reverse(v.begin(), v.end());

    if (v.size() != component.GetContigs().size())
    {
        component.begin = ContigNodeAdapter(NULL);
        component.end = ContigNodeAdapter(NULL);
        return;
    }

    map<ContigNodeAdapter, int> dist;
    map<ContigNodeAdapter, ContigNodeAdapter> prev;
    dist[begin] = 0;
    prev[begin] = begin;

    for (unsigned i = 0; i < v.size(); ++i)
    {
        ContigNodeAdapter curr = v[i];
        vector<ContigNodeAdapter> neighbors;
        GetNeighbors(curr, neighbors);

        for (unsigned j = 0; j < neighbors.size(); ++j)
        {
            int tmp = dist[curr] + curr.GetSize() - kmerLength + 1;
            if (tmp > dist[neighbors[j]])
            {
                dist[neighbors[j]] = tmp;
                prev[neighbors[j]] = curr;
            }
        }

        order.push_back(QueueNode(curr, dist[curr]));
    }
    sort(order.begin(), order.end());

    vector<ContigNodeAdapter> tmp;
    ContigNodeAdapter x = end;
    while (true)
    {
        tmp.push_back(x);

        if (x == prev[x])
            break;
        else
            x = prev[x];
    }
    reverse(tmp.begin(), tmp.end());

    for (unsigned i = 0; i < tmp.size(); ++i)
        path.Append(tmp[i]);

    path.Merge(component.longest_path);
}

void ComponentGraph::PrintComponent(Component&component)
{
    vector<string> buffer_strings;
    ContigNodeAdapter begin = component.GetBeginContig();
    ContigNodeAdapter end = component.GetEndContig();

    vector<QueueNode> order;
    TopSort(component, order, component.longest_path);

    vector<ContigAlignment> contig_alignments;
    for (unsigned i = 0; i < order.size(); ++i)
    {
        ContigNodeAdapter curr = order[i].node;
        int d = order[i].distance;

        Contig contig;
        curr.GetContig(contig);
        component.PrintAlignment(contig.ToString(), (d > 0) ? d + kmerLength - 1 : 0);

        if (d + contig.Size() - kmerLength + 1 > component.longest_path.Size())
            continue;

        ContigAlignment contig_alignment;
        if (curr != end)
            AlignContigToLongestPath(component.longest_path, d, d + contig.Size() - kmerLength + 1, 
                    contig, contig.Size() - kmerLength + 1, contig_alignment);
        else
            AlignContigToLongestPath(component.longest_path, d, d + contig.Size(), 
                    contig, contig.Size(), contig_alignment);
        contig_alignments.push_back(contig_alignment);
    }

    vector<vector<string> > gap_buffers(component.longest_path.Size());
    for (unsigned i = 0; i < contig_alignments.size(); ++i)
    {
        vector<string> gap(component.longest_path.Size());
        ContigAlignment contig_alignment = contig_alignments[i];
        vector<int> &pos = contig_alignment.pos;
        Sequence &seq = contig_alignment.seq;
        for (unsigned j = 0; j < pos.size(); ++j)
        {
            if (pos[j] < 0)
                gap[-pos[j]] += seq[j];
        }

        for (unsigned j = 0; j < gap.size(); ++j)
        {
            if (gap[j] != "")
                gap_buffers[j].push_back(gap[j]);
        }
    }

    vector<string> max_gap(component.longest_path.Size());
    vector<int> offset(component.longest_path.Size() + 1);
    offset[0] = 0;
    Contig consensus;
    for (int i = 0; i < component.longest_path.Size(); ++i)
    {
        consensus.AddNucleotide(component.longest_path[i]);

        string gap;
        if (gap_buffers[i].size() > 0)
        {
            for (unsigned j = 0; j < gap_buffers[i].size(); ++j)
            {
                if (gap_buffers[i][j].size() > gap.size())
                    gap = gap_buffers[i][j];
            }
        }

        max_gap[i] = gap;
        consensus += gap;

        offset[i+1] = offset[i] + max_gap[i].size() + 1;
    }
    component.consensus = consensus;

    vector<vector<int> > votes(consensus.Size());
    for (unsigned i = 0; i < votes.size(); ++i)
        votes[i].resize(4, 0);

    vector<string> buffers;
    buffers.push_back(consensus.ToString());
    for (unsigned i = 0; i < contig_alignments.size(); ++i)
    {
        ContigAlignment contig_alignment = contig_alignments[i];
        string seq = contig_alignment.seq.ToString();
        int from = abs(contig_alignment.pos.front());
        int to = abs(contig_alignment.pos.back()) + 1;

        string result(offset[to], ' ');

        vector<int> &pos = contig_alignment.pos;
        vector<string> gaps(result.size());
        for (unsigned i = 0; i < pos.size(); ++i)
        {
            if (pos[i] >= 0)
            {
                result[offset[pos[i]]] = seq[i];
                ++votes[offset[pos[i]]][CharToCode(seq[i])];
            }
            else
            {
                gaps[-pos[i]] += contig_alignment.seq[i];
                result[offset[-pos[i]] + gaps[-pos[i]].size()] = seq[i];
                ++votes[offset[-pos[i]] + gaps[-pos[i]].size()][CharToCode(seq[i])];
            }
        }

        for (int i = offset[from]; i < offset[to]; ++i)
        {
            if (result[i] == ' ')
                result[i] = '_';
        }

        int k = -1;
        for (unsigned j = 0; j < buffers.size(); ++j)
        {
            if (buffers[j].size() <= offset[from])
            {
                k = j;
                break;
            }
        }

        if (k == -1)
        {
            buffers.push_back("");
            k = buffers.size()-1;
        }

        if (buffers[k].size() < offset[from])
            buffers[k].append(offset[from] - buffers[k].size(), ' ');

        buffers[k].append(result.substr(offset[from], offset[to] - offset[from]));
    }

    for (unsigned i = 0; i < consensus.Size(); ++i)
    {
        int x = max_element(votes[i].begin(), votes[i].end()) - votes[i].begin();
        if (votes[i][x] != 0)
            consensus[i] = x;
    }

    component.consensus = consensus;
    component.alignment = buffers;
    component.CompleteAlignment();
}

void ComponentGraph::RemoveCycles(ContigNodeAdapter &curr)
{
    if (status[curr] != 0)
        return;

    status[curr] = 1;
    for (int x = 0; x < 4; ++x)
    {
        if ((1 << x) & curr.OutEdges())
        {
            ContigNodeAdapter next = GetNeighbor(curr, x);

            if (status[next] == 0)
                RemoveCycles(next);
            else if (status[next] == 1)
            {
                RemoveEdge(curr, x);
            }
        }
    }
    status[curr] = 2;
}

void ComponentGraph::RemoveCycles()
{
    Refresh();

    status.clear();
    for (unsigned i = 0; i < components.size(); ++i)
    {
        if (!components[i].GetBeginContig() || !components[i].GetEndContig())
            continue;

        ContigNodeAdapter begin = components[i].GetBeginContig();
        RemoveCycles(begin);
    }

    Refresh();
}

void ComponentGraph::BuildComponents()
{
    RemovePalindrome();

    deque<set<ContigNodeAdapter> > component_sets;
    deque<string> component_strings;

    stringstream origin_edges;
    WriteEdges(origin_edges);
    Decomposite();
    GetComponents(component_sets, component_strings);

    components.resize(component_sets.size());
    for (unsigned i = 0; i < components.size(); ++i)
    {
        components[i].SetContent(component_sets[i]);
        components[i].component.graph = component_strings[i];
        components[i].Data() = i;

        ComponentNodeAdapter adp(&components[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            deque<ContigNodeAdapter> contigs;
            adp.GetContigs(contigs);

            for (unsigned j = 0; j < contigs.size(); ++j)
                component_map[contigs[j]] = adp;

            adp.ReverseComplement();
        }
    }
    RemoveCycles();

    status.clear();
    int not_unique = 0;
    int not_unique_length = 0;
    int num_components = components.size();
    for (int i = 0; i < num_components; ++i)
    {
        ContigNodeAdapter begin = components[i].GetBeginContig();
        ContigNodeAdapter end = components[i].GetEndContig();
        ContigPath path;

        if (!begin || !end)
        {
            ++not_unique;
            not_unique_length += components[i].GetSize();

            deque<ContigNodeAdapter> contigs = components[i].GetContigs();
            components[i].SetContent(contigs[0]);

            for (unsigned j = 0; j < contigs.size(); ++j)
            {
                ContigNodeAdapter adp = contigs[j];
                for (int strand = 0; strand < 2; ++strand)
                {
                    for (int x = 0; x < 4; ++x)
                    {
                        if ((1 << x) & adp.OutEdges())
                            RemoveEdge(adp, x);
                    }
                    adp.ReverseComplement();
                }
            }

            for (unsigned j = 1; j < contigs.size(); ++j)
            {
                ComponentNode comp;
                comp.SetContent(contigs[j]);
                components.push_back(comp);

                components.back().Data() = components.size()-1;

                ComponentNodeAdapter adp(&components.back());
                for (int strand = 0; strand < 2; ++strand)
                {
                    deque<ContigNodeAdapter> contigs;
                    adp.GetContigs(contigs);

                    for (unsigned k = 0; k < contigs.size(); ++k)
                        component_map[contigs[k]] = adp;

                    adp.ReverseComplement();
                }
            }
        }
    }

    Refresh();
    cerr << "not unique " << not_unique << " "<< components.size() << " " << not_unique_length << endl;

    for (unsigned i = 0; i < components.size(); ++i)
        PrintComponent(components[i].component);

    ReadEdges(origin_edges);
    Refresh();

    in_component_connections.resize(components.size());
    out_component_connections.resize(components.size());
}


void ComponentGraph::AlignContigToLongestPath(Contig &longest, int from, int to, Contig &contig, int contig_length, ContigAlignment &contig_alignment)
{
    Sequence a;
    Sequence b;
    a.Assign(longest, from, to - from);
    b.Assign(contig, 0, contig_length);


    contig_alignment.seq = b;
    contig_alignment.pos.resize(b.Size());

    if (a == b)
    {
        for (int i = 0; i < b.Size(); ++i)
            contig_alignment.pos[i] = i + from;
    }
    else
    {
        vector<vector<int> > table;
        table.resize(a.Size() + 1);
        for (unsigned i = 0; i < table.size(); ++i)
            table[i].resize(b.Size() + 1);

        for (int i = 0; i <= a.Size(); ++i)
            table[i][0] = i;

        for (int j = 0; j <= b.Size(); ++j)
            table[0][j] = j;

        for (int i = 1; i <= a.Size(); ++i)
        {
            for (int j = 1; j <= b.Size(); ++j)
            {
                table[i][j] = 1000000000;
                if (table[i-1][j] + 1 < table[i][j])
                    table[i][j] = table[i-1][j] + 1;
                if (table[i][j-1] + 1 < table[i][j])
                    table[i][j] = table[i][j-1] + 1;
                if (table[i-1][j-1] + (a[i-1] != b[j-1]) < table[i][j])
                    table[i][j] = table[i-1][j-1] + (a[i-1] != b[j-1]);
            }
        }

        int i = a.Size();
        int j = b.Size();

        while (true)
        {
            if (i == 0 && j == 0)
                break;

            if (i > 0 && j > 0 && table[i-1][j-1] + (a[i-1] != b[j-1]) == table[i][j])
            {
                contig_alignment.pos[j-1] = from + i-1;
                --i, --j;
            }
            else if (i > 0 && table[i-1][j] + 1 == table[i][j])
            {
                --i;
            }
            else if (j > 0 && table[i][j-1] + 1 == table[i][j])
            {
                contig_alignment.pos[j-1] = -(from + i-1);
                --j;
            }
            else
            {
                cerr << "error" << endl;
                exit(1);
            }
        }
    }
}

void ComponentGraph::ConnectComponents(deque<Component> &output_components)
{
    int valid_connections = 0;
    int unique_connections = 0;
    int add_connections = 0;
    for (unsigned i = 0; i < components.size(); ++i)
    {
        ComponentNodeAdapter curr(&components[i]);
        for (int strand = 0; strand < 2; ++strand)
        {
            vector<ComponentConnection> &connections = GetComponentConnections(curr);

            int index = 0;
            for (unsigned j = 0; j < connections.size(); ++j)
            {
                ContigNodeAdapter from = curr.GetEndContig();
                ContigNodeAdapter to = connections[j].to.GetBeginContig();
                ContigPath path;
                delta = estimate_distance;
                max_paths = 10;
                if (connections[j].count > min_pairs && !!from && !!to && FindPath(from, to, path, 0))
                {
                    bool flag = true;

                    if (flag)
                    {
                        ++valid_connections;
                        connections[index++] = connections[j];
                    }
                }
            }

            connections.resize(index);

            if (connections.size() == 1)
                ++unique_connections;

            curr.ReverseComplement();
        }
    }

    for (unsigned i = 0; i < components.size(); ++i)
    {
        ComponentNodeAdapter curr(&components[i]);
        if (curr.IsUsed() || curr.IsDead())
            continue;

        curr.SetUsedFlag();
        Component component;
        component = components[i].component;

        for (int strand = 0; strand < 2; ++strand)
        {
            curr.SetNode(&components[i], strand);

            while (true)
            {
                vector<ComponentConnection> &connections = GetComponentConnections(curr);
                ComponentNodeAdapter rev_comp = curr;
                rev_comp.ReverseComplement();

                if (connections.size() == 1)
                {
                    ComponentNodeAdapter to = connections[0].to;
                    to.ReverseComplement();
                    if (GetComponentConnections(to).size() == 1 && GetComponentConnections(to)[0].to == rev_comp
                            && !connections[0].to.IsUsed())
                    { 
                        ContigNodeAdapter from = curr.GetEndContig(); 
                        ContigNodeAdapter to = connections[0].to.GetBeginContig(); 
                        ContigPath path; 
                        FindPath(from, to, path, 0); 
                        vector<ContigNodeAdapter> &contig_nodes = path.GetNodes();
                        for (unsigned i = 1; i+1 < contig_nodes.size(); ++i)
                        {
                            component.Add(contig_nodes[i]);

                            if (component_map[contig_nodes[i]].GetNode()->component.contigs.size() == 1)
                            {
                                component_map[contig_nodes[i]].SetUsedFlag();
                            }
                        }

                        ComponentNode tmp_node = *connections[0].to.GetNode();
                        if (connections[0].to.IsReverse())
                            tmp_node.ReverseComplement();

                        component.Add(tmp_node.component);

                        curr = connections[0].to;
                        curr.SetUsedFlag();
                    }
                    else
                        break;
                }
                else
                    break;
            }

            component.ReverseComplement();
        }

        output_components.push_back(component);
    }

    cout << output_components.size() << endl;
}

void ComponentGraph::AddComponentPair(vector<Alignment> &alignments1, vector<Alignment> &alignments2)
{
    if (alignments1.size() > 2 && alignments2.size() > 2)
        return;

    for (unsigned i = 0; i < alignments1.size(); ++i)
    {
        for (unsigned j = 0; j < alignments2.size(); ++j)
        {
            AddComponentPair(alignments1[i], alignments2[j]);
        }
    }
}

void ComponentGraph::AddComponentPair(Alignment a1, Alignment a2)
{
    a2.ReverseComplement();

    ComponentNodeAdapter from = component_map[ContigNodeAdapter(&nodes[a1.contigId], a1.isReverse)];
    ComponentNodeAdapter to = component_map[ContigNodeAdapter(&nodes[a2.contigId], a2.isReverse)];

    if (from.GetNode() != to.GetNode() && from.GetNode()->component.longest_path.Size() > estimate_distance
            && to.GetNode()->component.longest_path.Size() > estimate_distance)
        AddComponentConnection(from, to);
}

