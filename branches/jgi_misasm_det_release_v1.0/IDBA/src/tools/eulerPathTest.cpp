#include "globals.h"

#include "Kmer.h"
#include "OptionSet.h"
#include "SequenceReader.h"
#include "Sequence.h"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <map>
#include <queue>
#include <stack>

using namespace std;

struct Node
{
    bool is_connected;
    Kmer kmer;
    int in[4];
    int out[4];

    int InDegree() { return in[0] + in[1] + in[2] + in[3]; }
    int OutDegree() { return out[0] + out[1] + out[2] + out[3]; }
};

map<Kmer, int> dict;
vector<Node> nodes;

void Clear()
{
    dict.clear();
    nodes.resize(0);
}

Node &GetNode(const Kmer &kmer)
{
    if (dict.find(kmer) == dict.end())
    {
        int id = dict.size();
        dict[kmer] = id;

        nodes.resize(nodes.size() + 1);
        nodes[id].kmer = kmer;
        for (int i = 0; i < 4; ++i)
            nodes[id].in[i] = nodes[id].out[i] = 0;
    }

    return nodes[dict[kmer]];
}

void AddEdge(const Kmer &kmer, int c)
{
    Node &node1 = GetNode(kmer);
    node1.out[c]++;

    int d = kmer.GetBase(0);
    Kmer kmer2 = kmer;
    kmer2.AddRight(c);
    Node &node2 = GetNode(kmer2);
    node2.in[d]++;
}

void RemoveEdge(const Kmer &kmer, int c)
{
    Node &node1 = GetNode(kmer);
    node1.out[c]--;

    int d = kmer.GetBase(0);
    Kmer kmer2 = kmer;
    kmer2.AddRight(c);
    Node &node2 = GetNode(kmer2);
    node2.in[d]--;
}

void BFS(const Kmer &kmer)
{
    queue<Kmer> qu;
    qu.push(kmer);

    GetNode(kmer).is_connected = true;

    while (!qu.empty())
    {
        Node &node = GetNode(qu.front());
        qu.pop();

        for (int i = 0; i < 4; ++i)
        {
            if (node.out[i])
            {
                Kmer kmer = node.kmer;
                kmer.AddRight(i);


                if (!GetNode(kmer).is_connected)
                {
                    GetNode(kmer).is_connected = true;
                    qu.push(kmer);
                }
            }
        }
    }
}

bool IsConnected(const Kmer &kmer)
{
    for (unsigned i = 0; i < nodes.size(); ++i)
        nodes[i].is_connected = false;

    BFS(kmer);

    for (unsigned i = 0; i < nodes.size(); ++i)
    {
        if ((nodes[i].InDegree() != 0 || nodes[i].OutDegree() != 0) && !nodes[i].is_connected)
            return false;
    }
    return true;
}

bool FindEulerPath(const Kmer &begin, const Kmer &end, vector<Kmer> &path)
{
    stack<Kmer> st;
    st.push(begin);

    while (!st.empty())
    {
        Kmer kmer = st.top();
        Node &node = GetNode(kmer);

        if (node.OutDegree() == 0)
        {
            path.push_back(kmer);
            st.pop();
        }
        else
        {
            for (int i = 0; i < 4; ++i)
            {
                if (node.out[i])
                {
                    RemoveEdge(kmer, i);
                    Kmer next = kmer;
                    next.AddRight(i);
                    st.push(next);
                    break;
                }
            }
        }
    }

    reverse(path.begin(), path.end());

    return true;

//    path.push_back(begin);
//
//    while (true)
//    {
//        Kmer kmer = path.back();
//        vector<int> candidate;
//        vector<int> cut;
//        Node &node = GetNode(kmer);
//
//        if (node.OutDegree() == 0)
//            break;
//
//        for (int i = 0; i < 4; ++i)
//        {
//            if (node.out[i])
//            {
//                RemoveEdge(kmer, i);
//                if (IsConnected(kmer))
//                    candidate.push_back(i);
//                else
//                    cut.push_back(i);
//                AddEdge(kmer, i);
//            }
//        }
//
//        if (candidate.size() > 1)
//        {
//            cout << "too many candidate" << endl;
//            return false;
//        }
//        else if (candidate.size() == 0 && cut.size() != 1)
//        {
//            cout << "impossible" << endl;
//            return false;
//        }
//        else
//        {
//            int c;
//            if (candidate.size() == 1)
//                c = candidate[0];
//            else
//                c = cut[0];
//
//            RemoveEdge(kmer, c);
//            kmer.AddRight(c);
//            path.push_back(kmer);
//        }
//    }
//
//    for (unsigned i = 0; i < nodes.size(); ++i)
//    {
//        if (nodes[i].InDegree() != 0 || nodes[i].OutDegree() != 0)
//            return false;
//    }
//
//    return true;
}


int num = 10000;
int fragment_length = 200;
int read_length = 5;

int main(int argc, char *argv[])
{
    OptionSet option_set;
    option_set.AddOption("num", "", num, "");
    option_set.AddOption("fragmentLength", "", fragment_length, "");
    option_set.AddOption("readLength", "", read_length, "");
    option_set.ProcessOptions(argc, argv);

    if (argc < 2)
    {
        fprintf(stderr, "Usage: eulerPathTest ref-file\n");
        exit(1);
    }

    kmerLength = read_length - 1;

    vector<Sequence> refs;
    Sequence seq;
    string comment;

    FastaReader reader(argv[1]);
    while (reader.Read(seq, comment))
    {
        //seq.Encode();
        refs.push_back(seq);
    }

    int sum_loop = 0;
    int num_success = 0;
    for (int round = 0; round < num; ++round)
    {
        int sum = 0;
        for (unsigned i = 0; i < refs.size(); ++i)
            sum += refs[i].Size();
        int r = rand() % sum;

        sum = 0;
        int index = 0;
        while (true)
        {
            sum += refs[index].Size();
            if (r < sum)
                break;
            else
                ++index;
        }

        while (true)
        {
            int offset = rand() % (refs[index].Size() - fragment_length + 1);
            //refs[index].GetSubSequence(seq, offset, fragment_length);
            seq.Assign(refs[index], offset, fragment_length);

            bool flag = true;
            for (int i = 0; i < fragment_length; ++i)
            {
                if (seq[i] == 4)
                {
                    flag = false;
                    break;
                }
            }

            if (flag)
                break;
        }


//    refs.resize(num);
//    for (unsigned i = 0; i < refs.size(); ++i)
//    {
//        Sequence seq = refs[i];

        Clear();
        int len = 0;
        Kmer kmer;
        for (int i = 0; i < seq.Size(); ++i)
        {
            kmer.AddRight(seq[i]);
            ++len;
            if (len >= kmerLength && i+1 < seq.Size())
                AddEdge(kmer, seq[i+1]);
        }

        Kmer begin;
        Kmer end;

        int loop = 0;
        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            if (nodes[i].InDegree() < nodes[i].OutDegree())
                begin = nodes[i].kmer;
            if (nodes[i].InDegree() > nodes[i].OutDegree())
                end = nodes[i].kmer;

            if (nodes[i].InDegree() > 1 || nodes[i].OutDegree() > 1)
                ++loop;
        }

        sum_loop += loop;

        vector<Kmer> path;

        if (FindEulerPath(begin, end, path))
        {
            Sequence result = begin;
            for (unsigned i = 1; i < path.size(); ++i)
                result += path[i].GetBase(kmerLength-1);

            if (seq == result)
            {
                cout << "success" << endl;
                ++num_success;
            }
            else
                cout << "fail" << endl;

            cout << seq << endl;
            cout << result << endl;
        }
        else
        {
            cout << "fail" << endl;
        }
    }

    cerr << sum_loop * 1.0 / refs.size() << endl;
    cerr << num_success << " " << refs.size() << endl;

    return 0;
}

