/**
 * @file compareMetaGraph.cpp
 * @brief 
 * @author Yu peng
 * @version 0.18
 * @date 2010-11-21
 */


#include "globals.h"

#include "BitOperation.h"
#include "HashGraph.h"
#include "KmerVector.h"
#include "Log.h"
#include "OptionSet.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "Utils.h"
#include "SequenceWriter.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>
#include <set>

using namespace std;

struct Node
{
    Node *parent;
    vector<Node *> childs;
    vector<string> refs;
    string taxon;
    bool isExist;
    KmerVector center;
    vector<KmerVector> points;
};

const int MaxNodes = 1000000;

char line[MaxLine];
char buf[MaxLine];
char name[MaxLine];
Node nodes[MaxNodes];
map<string, vector<int> > taxonDict;
map<string, Sequence> seq_table;
set<string> validTaxon;
vector<int> flags;
string taxonNames[] =
{
    "phylum", "class", "order", "family", "genus", "species",
};
bool is_random = false;

vector<string> refList;
vector<string> allRefs;
vector<KmerVector> allKmerVectors;
map<string, int> refRank;
string level = "genus";
int size = 2;
vector<vector<int> > tuples;
map<int, map<int, int> > graph;
vector<int> species;
int max_sub_species = 100;
int min_sub_species = 2;
string type = "meta";

HashGraph hash_graph;

void GenTuple(vector<int> &v, int k)
{
    if (v.size() == size)
    {
        tuples.push_back(v);
    }
    else
    {
        for (int i = k; i < (int)species.size(); ++i)
        {
            if (nodes[species[i]].refs.size() <= min_sub_species)
                continue;
//            if (nodes[species[i]].refs.size() <= 10)
//                continue; 

            bool flag = true;
            for (int j = 0; j < (int)v.size(); ++j)
            {
                if (graph[species[i]][v[j]] == 0)
                {
                    flag = false;
                    break;
                }
            }

            if (flag)
            {
                v.push_back(species[i]);
                GenTuple(v, i + 1);
                v.pop_back();
            }
        }
    }
}

void GetPath(Node *p, vector<Node *> &path)
{
    path.resize(0);
    while (true)
    {
        path.push_back(p);

        if (p != p->parent)
            p = p->parent;
        else
            break;
    }
    reverse(path.begin(), path.end());
}

void PrintPath(vector<Node *> &path)
{
    for (unsigned i = 0; i < path.size(); ++i)
        cout << (path[i] - nodes) << " " << path[i]->taxon << " ";
}

string LowestCommonAnsestor(Node *p, Node *q)
{
    vector<Node *> path1;
    vector<Node *> path2;

    GetPath(p, path1);
    GetPath(q, path2);

    int l = 0;
    while (l < (int)path1.size() && l < (int)path2.size())
    {
        if (path1[l] == path2[l])
            ++l;
        else
            break;
    }

    while (validTaxon.find(path1[l-1]->taxon) == validTaxon.end())
        --l;

//    if (path1[l-1]->taxon == level)
//    {
//        cout << (p - nodes) << " " << (q - nodes) << ",";
//        PrintPath(path1);
//        cout << ",";
//        PrintPath(path2);
//        cout << endl;
//    }

    return path1[l-1]->taxon;
}

void dfs(Node *node)
{
    if (node->childs.size() > 0)
    {
        int count = 0;
        node->center.Clear();
        for (unsigned i = 0; i < node->childs.size(); ++i)
        {
            dfs(node->childs[i]);
            if (node->childs[i]->refs.size() > 0)
            {
                for (unsigned j = 0; j < node->childs[i]->refs.size(); ++j)
                {
                    node->refs.push_back(node->childs[i]->refs[j]);
                    node->points.push_back(node->childs[i]->points[j]);
                }
                node->center += node->childs[i]->center;
                ++count;
            }
        }

        if (count != 0)
            node->center /= count;
    }
}

void BuildTaxonTree(char *nodes_file, char *gff_list)
{
    FILE *fnodes = OpenFile(nodes_file, "rb");
    FILE *fgff = OpenFile(gff_list, "rb");

    for (int i = 0; i < MaxNodes; ++i)
        nodes[i].isExist = false;

    while (fgets(line, MaxLine, fnodes) != NULL)
    {
        int id, pid;
        sscanf(line, "%d %s %d %s %s", &id, buf, &pid, buf, name);
        nodes[id].parent = &nodes[pid];
        nodes[id].taxon = name;
        nodes[id].isExist = true;

        if (id != pid)
            nodes[pid].childs.push_back(&nodes[id]);
    }

    fclose(fnodes);

    while (fscanf(fgff, "%s", name) != EOF)
    {
        FILE *fp = OpenFile(name, "rb");

        while (fgets(line, MaxLine, fp) != NULL)
        {
            if (line[0] == '#')
                continue;
            else
            {
                if (strstr(line, "plasmid") == NULL)
                {
                    char *p = strstr(line, "taxon");

                    if (p == NULL)
                    {
                        printf("%s\n", name);
                        fprintf(stderr, "%s\n", p);
                        fprintf(stderr, "gff format error\n");
                        exit(1);
                    }

                    int id = atoi(p + 6);

                    if (nodes[id].isExist == false)
                    {
                        fprintf(stderr, "id not exist\n");
                        exit(1);
                    }

                    strcpy(name + strlen(name)-3, "fna");

                    bool flag = false;
                    for (unsigned i = 0; i < refList.size(); ++i)
                    {
                        if (name == refList[i])
                        {
                            flag = true;
                            break;
                        }
                    }

                    if (flag)
                        break;

                    nodes[id].refs.push_back(name);

                    FILE *fref = OpenFile(name, "rb");
                    Sequence seq;
                    //ReadFasta(fref, &seq);

                    FastaReader reader(name);
                    string comment;
                    reader.Read(seq, comment);
                    //seq.Encode();
                    seq_table[name] = seq;
                    nodes[id].center.Compute(&seq);
                    nodes[id].points.push_back(nodes[id].center);
                    fclose(fref);

                    allRefs.push_back(name);
                    allKmerVectors.push_back(nodes[id].center);
                }

                break;
            }
        }

        fclose(fp);
    }

    fclose(fgff);

    //seq_table.clear();
    for (int i = 0; i < MaxNodes; ++i)
    {
        if (nodes[i].isExist && nodes[i].refs.size() > 0)
        {
            string name = nodes[i].refs[0];
            //cout << name << endl;
            name.resize(name.rfind('/'));
            name = FormatString("%s/%d.fa", name.c_str(), i);

            //cout << name << endl;
            //FastaWriter writer(name);
            Sequence seq;
            for (unsigned j = 0; j < nodes[i].refs.size(); ++j)
            {
                seq += seq_table[nodes[i].refs[j]];
                seq.AddNucleotide(4);
            }

//            if (nodes[i].refs.size() > 1)
//            {
//                for (unsigned j = 0; j < nodes[i].refs.size(); ++j)
//                    cout << nodes[i].refs[j] << endl;
//                cout << name << endl;
//                cout << endl;
//            }
//
//            writer.Write(seq, name);

            seq_table[name] = seq;
            nodes[i].refs.resize(1);
            nodes[i].refs[0] = name;
        }
    }

    for (int i = 0; i < MaxNodes; ++i)
    {
        if (nodes[i].isExist)
        {
            Node *node = nodes[i].parent;
            while (node != node->parent)
            {
                bool found = false;
                for (int j = 0; j < 6; ++j)
                {
                    if (node->taxon == taxonNames[j])
                    {
                        found = true;
                        break;
                    }
                }

                if (found)
                    break;

                node = node->parent;
            }

            nodes[i].parent = node;
        }
    }

    dfs(&nodes[1]);

    for (int i = 0; i < MaxNodes; ++i)
    {
        if (nodes[i].isExist && nodes[i].refs.size() > 0)
        {
            if (validTaxon.find(nodes[i].taxon) != validTaxon.end())
            {
                vector<Node *> path;
                GetPath(&nodes[i], path);
                bool flag = true;
                for (unsigned j = 0; j+1 < path.size(); ++j)
                {
                    if (path[j]->taxon == nodes[i].taxon)
                    {
                        flag = false;
                        break;
                    }
                }

                if (flag)
                    taxonDict[nodes[i].taxon].push_back(i);
            }
        }
    }
}

void Compare(Node *p, Node *q)
{
    string lca = LowestCommonAnsestor(p, q);

    if (lca == "root")
        return;

    if (p->refs.size() < 2 || q->refs.size() < 2)
        return;

    printf("%s, %ld, %ld, %d, %d, ", lca.c_str(), p - nodes, q - nodes, (int)p->refs.size(), (int)q->refs.size());
    fflush(NULL);

    string file1 = FormatString("/tmp/meta_control_%s1.fa", level.c_str());
    string file2 = FormatString("/tmp/meta_control_%s2.fa", level.c_str());

    FastaWriter writer1(file1);
//    writer1.Write(seq_table[p->refs[0]]);
//    writer1.Write(seq_table[p->refs[1]]);
    for (unsigned i = 0; i < min(3, int(p->refs.size())); ++i)
        writer1.Write(seq_table[p->refs[i]]);

    FastaWriter writer2(file2);
//    writer2.Write(seq_table[q->refs[0]]);
//    writer2.Write(seq_table[q->refs[1]]);
    for (unsigned i = 0; i < min(3, int(q->refs.size())); ++i)
        writer2.Write(seq_table[q->refs[i]]);

    string file_mix = FormatString("/tmp/meta_control_%s_mix.fa", level.c_str());
    string file_reads = FormatString("/tmp/meta_control_%s_reads.fa", level.c_str());

    string command = FormatString("cat %s %s > %s", file1.c_str(), file2.c_str(), file_mix.c_str());
    system(command.c_str());

    command = FormatString("bin/idba --minCount 1 -r empty -l %s --mink %d --maxk %d --cover 0.5 -o /tmp/%s_control|tail -1",  
            file_mix.c_str(), kmerLength, kmerLength, level.c_str());
    system(command.c_str());


}

void TestMeta(vector<int> &v)
{
    printf("%s, ", level.c_str());
    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", v[i]);

    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", (int)nodes[v[i]].refs.size());

    fflush(NULL);

    vector<string> ref_names;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        string filename = FormatString("/tmp/meta_%s%u.fa", level.c_str(), i);
        FastaWriter writer(filename);
        for (unsigned j = 0; j < min(max_sub_species, int(nodes[v[i]].refs.size())); ++j)
            writer.Write(seq_table[nodes[v[i]].refs[j]]);
        ref_names.push_back(filename);
    }

    string file_mix = FormatString("/tmp/meta_%s_mix.fa", level.c_str());

    stringstream ss_command;
    ss_command << "cat ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << "> " << file_mix;

    system(ss_command.str().c_str());

    ss_command.str("");

    ss_command << "bin/metaAssembler -n " << size << " -k " << kmerLength << " ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << FormatString("-o /tmp/meta_%s | tail -1", level.c_str());

    system(ss_command.str().c_str());


}

void TestControl(vector<int> &v)
{
    printf("%s, ", level.c_str());
    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", v[i]);

    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", (int)nodes[v[i]].refs.size());

    fflush(NULL);

    vector<string> ref_names;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        string filename = FormatString("/tmp/meta_control_%s%u.fa", level.c_str(), i);
        FastaWriter writer(filename);
        for (unsigned j = 0; j < min(max_sub_species, int(nodes[v[i]].refs.size())); ++j)
            writer.Write(seq_table[nodes[v[i]].refs[j]], nodes[v[i]].refs[j]);
        ref_names.push_back(filename);
    }

    string file_mix = FormatString("/tmp/meta_control_%s_mix.fa", level.c_str());

    stringstream ss_command;
    ss_command << "cat ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << "> " << file_mix;

    system(ss_command.str().c_str());

    ss_command.str("");

    string command = FormatString("bin/idba --minCount 1 -r empty -l %s --mink %d --maxk %d --cover 0.5 -o /tmp/meta_control_%s|tail -1",  
            file_mix.c_str(), kmerLength, kmerLength, level.c_str());
    system(command.c_str());

    command = FormatString("script/validateBlat /tmp/meta_control_%s-contig.fa %s 0.95 0", level.c_str(), file_mix.c_str());
    system(command.c_str());

    for (unsigned i = 0; i < ref_names.size(); ++i)
    {
        for (unsigned j = i+1; j < ref_names.size(); ++j)
        {
            string command = FormatString("bin/compareKmer -k %d %s %s", kmerLength, ref_names[i].c_str(), ref_names[j].c_str());
            system(command.c_str());
        }
    }

    //exit(1);
}

void TestMetaIDBA(vector<int> &v)
{
    printf("%s, ", level.c_str());
    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", v[i]);

    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", (int)nodes[v[i]].refs.size());

    fflush(NULL);

    vector<string> ref_names;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        string filename = FormatString("/tmp/meta_idba_%s%u.fa", level.c_str(), i);
        FastaWriter writer(filename);
        for (unsigned j = 0; j < min(max_sub_species, int(nodes[v[i]].refs.size())); ++j)
        {
            writer.Write(seq_table[nodes[v[i]].refs[j]], nodes[v[i]].refs[j]);
        }
        ref_names.push_back(filename);
    }

    string file_mix = FormatString("/tmp/meta_idba_%s_mix.fa", level.c_str());
    string file_reads = FormatString("/tmp/meta_idba_%s_reads.fa", level.c_str());

    stringstream ss_command;
    ss_command << "cat ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << "> " << file_mix;

    system(ss_command.str().c_str());

    ss_command.str("");
    ss_command << "bin/metaSimReads --readLength 75 --depth 30 --mate --num " << size << " ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << file_reads << " ";
    string command = ss_command.str();
    //string command = FormatString("bin/simReads %s %s --readLength 75 --depth 30 --mate", file_mix.c_str(), file_reads.c_str());
    if (is_random)
        command += " --type 2";
    system(command.c_str());

    ss_command.str("");

    ss_command << FormatString("time bin/metaidba -r %s --maxk %d -o /tmp/meta_idba_%s ", file_reads.c_str(), kmerLength, level.c_str());
    ss_command << " -n " << size << " ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << " | tail -1 ";

    system(ss_command.str().c_str());

    command = FormatString("bin/rawN50 /tmp/meta_idba_%s.long", level.c_str());
    system(command.c_str());
    fflush(NULL);

//    command = FormatString("script/validateBlat /tmp/meta_idba_%s.long %s 0.95 0", level.c_str(), file_mix.c_str());
//    system(command.c_str());
//
//    command = FormatString("blat -tileSize=18 -minMatch=2 %s /tmp/meta_idba_%s.component /tmp/meta_idba_%s.blat > /dev/null", 
//            file_mix.c_str(), level.c_str(), level.c_str());
//    system(command.c_str());
//
//    ss_command.str("");
//    ss_command << FormatString("bin/validateComponent -n %d /tmp/meta_idba_%s.component < /tmp/meta_idba_%s.blat ", 
//            size, level.c_str(), level.c_str());
//    for (unsigned i = 0; i < ref_names.size(); ++i)
//        ss_command << ref_names[i] << " ";
//    system(ss_command.str().c_str());
}

void TestMetaVelvet(vector<int> &v)
{
    printf("%s, ", level.c_str());
    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", v[i]);

    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", (int)nodes[v[i]].refs.size());

    fflush(NULL);

    vector<string> ref_names;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        string filename = FormatString("/tmp/meta_velvet_%s%u.fa", level.c_str(), i);
        FastaWriter writer(filename);
        for (unsigned j = 0; j < min(max_sub_species, int(nodes[v[i]].refs.size())); ++j)
        {
            writer.Write(seq_table[nodes[v[i]].refs[j]], nodes[v[i]].refs[j]);
        }
        ref_names.push_back(filename);
    }

    string file_mix = FormatString("/tmp/meta_velvet_%s_mix.fa", level.c_str());
    string file_reads = FormatString("/tmp/meta_velvet_%s_reads.fa", level.c_str());

    stringstream ss_command;
    ss_command << "cat ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << "> " << file_mix;

    system(ss_command.str().c_str());

    ss_command.str("");
    //ss_command << "bin/metaSimReads --readLength 75 --depth 30 --mate ";
    ss_command << "bin/metaSimReads --readLength 75 --depth 30 --mate --num " << size << " ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << file_reads << " ";
    string command = ss_command.str();
    //string command = FormatString("bin/simReads %s %s --readLength 75 --depth 30 --mate", file_mix.c_str(), file_reads.c_str());
    if (is_random)
        command += " --type 2";
    system(command.c_str());

    ss_command.str("");

    command = FormatString("rm -rf /tmp/meta_velvet_%s", level.c_str());
    system(command.c_str());
    fflush(NULL);

    command = FormatString("time script/velvet-pe /tmp/meta_velvet_%s %d %s | tail -1", level.c_str(), kmerLength, file_reads.c_str());
    system(command.c_str());
    fflush(NULL);

    command = FormatString("bin/splitScaffold /tmp/meta_velvet_%s/contigs.fa /tmp/meta_velvet_%s/contigs.fa.split", level.c_str(), level.c_str());
    system(command.c_str());
    fflush(NULL);

    command = FormatString("bin/rawN50 /tmp/meta_velvet_%s/contigs.fa.split", level.c_str());
    system(command.c_str());
    fflush(NULL);

//    command = FormatString("script/validateBlat /tmp/meta_velvet_%s/contigs.fa %s 0.95 0", level.c_str(), file_mix.c_str());
//    system(command.c_str());
//    fflush(NULL);
//
//    command = FormatString("script/validateBlat /tmp/meta_velvet_%s/contigs.fa.split %s 0.95 0", level.c_str(), file_mix.c_str());
//    system(command.c_str());
//    fflush(NULL);

}

void TestMetaAbyss(vector<int> &v)
{
    printf("%s, ", level.c_str());
    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", v[i]);

    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", (int)nodes[v[i]].refs.size());

    fflush(NULL);

    vector<string> ref_names;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        string filename = FormatString("/tmp/meta_abyss_%s%u.fa", level.c_str(), i);
        FastaWriter writer(filename);
        for (unsigned j = 0; j < min(max_sub_species, int(nodes[v[i]].refs.size())); ++j)
        {
            writer.Write(seq_table[nodes[v[i]].refs[j]], nodes[v[i]].refs[j]);
        }
        ref_names.push_back(filename);
    }

    string file_mix = FormatString("/tmp/meta_abyss_%s_mix.fa", level.c_str());
    string file_reads = FormatString("/tmp/meta_abyss_%s_reads.fa", level.c_str());

    stringstream ss_command;
    ss_command << "cat ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << "> " << file_mix;

    system(ss_command.str().c_str());

    ss_command.str("");
    //ss_command << "bin/metaSimReads --readLength 75 --depth 30 --mate ";
    ss_command << "bin/metaSimReads --readLength 75 --depth 30 --mate --num " << size << " ";
    for (unsigned i = 0; i < ref_names.size(); ++i)
        ss_command << ref_names[i] << " ";
    ss_command << file_reads << " ";
    string command = ss_command.str();

    //string command = FormatString("bin/simReads %s %s --readLength 75 --depth 30 --mate", file_mix.c_str(), file_reads.c_str());
    if (is_random)
        command += " --type 2";
    system(command.c_str());

    ss_command.str("");

    command = FormatString("rm -f /tmp/meta_abyss_out_%s*", level.c_str());
    system(command.c_str());

    command = FormatString("time abyss-pe k=%d n=5 in=%s name=/tmp/meta_abyss_out_%s |tail -1", kmerLength, file_reads.c_str(), level.c_str());
    system(command.c_str());

    command = FormatString("bin/rawN50 /tmp/meta_abyss_out_%s-contigs.fa", level.c_str());
    system(command.c_str());

    command = FormatString("script/validateBlat /tmp/meta_abyss_out_%s-contigs.fa %s 0.95 100", level.c_str(), file_mix.c_str());
    system(command.c_str());
}

void TestCompareKmer(vector<int> &v)
{
    printf("%s, ", level.c_str());
    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", v[i]);

    for (unsigned i = 0; i < v.size(); ++i)
        printf("%d, ", (int)nodes[v[i]].refs.size());

    if (v.size() == 1)
    {
        printf("%s, ", level.c_str());

        vector<string> x = nodes[v[0]].refs;
        vector<string> y = nodes[v[0]].refs;

        for (unsigned i = 0; i < x.size(); ++i)
        {
            for (unsigned j = i+1; j < y.size(); ++j)
            {
                string command = FormatString("bin/compareKmer %s %s", x[i].c_str(), y[j].c_str());
                //printf("%s ", command.c_str());
                fflush(NULL);
                system(command.c_str());
            }
        }
    }
    else if (v.size() == 2)
    {
        vector<string> x = nodes[v[0]].refs;
        vector<string> y = nodes[v[1]].refs;

        for (unsigned i = 0; i < x.size(); ++i)
        {
            for (unsigned j = 0; j < y.size(); ++j)
            {
                string command = FormatString("bin/compareKmer %s %s", x[i].c_str(), y[j].c_str());
                //printf("%s ", command.c_str());
                fflush(NULL);
                system(command.c_str());
            }
        }
    }
}

OptionSet option_set;
bool is_mix = false;

int main(int argc, char *argv[])
{
    //reverse(taxonNames, taxonNames + 6);

    option_set.AddOption("kmer", "", kmerLength, "k value");
    option_set.AddOption("level", "", level, "level");
    option_set.AddOption("size", "", size, "level");
    option_set.AddOption("max_sub", "", max_sub_species, "level");
    option_set.AddOption("min_sub", "", min_sub_species, "level");
    option_set.AddOption("mix", "", is_mix, "level");
    option_set.AddOption("type", "", type, "level");
    option_set.AddOption("random", "", is_random, "level");
    option_set.ProcessOptions(argc, argv);

    cout << kmerLength << " " << level << endl;

    if (argc < 3)
    {
        fprintf(stderr, "usage: cluster node.dmp gff-list [--kmer k] [--level l] [--size s] [--max_sub m] [--type t] [--mix] [--random]\n");
        exit(1);
    }

    for (int i = 0; i < 6; ++i)
        validTaxon.insert(taxonNames[i]);
    validTaxon.insert("root");
    BuildTaxonTree(argv[1], argv[2]);
    nodes[1].taxon = "root";

    species = taxonDict["species"];

    for (unsigned i = 0; i < species.size(); ++i)
    {
        for (unsigned j = i+1; j < species.size(); ++j)
        {
            if (is_mix)
            {
                string lca = LowestCommonAnsestor(&nodes[species[i]], &nodes[species[j]]);
                int x = find(taxonNames, taxonNames + 6, lca.c_str()) - taxonNames;
                if (x != 6 && x >= (find(taxonNames, taxonNames + 6, level.c_str()) - taxonNames))
                {
                    if (nodes[species[i]].refs.size() >= min_sub_species && nodes[species[j]].refs.size() >= min_sub_species)
                    {
                        graph[species[i]][species[j]] = 1;
                        graph[species[j]][species[i]] = 1;
                    }
                }
            }
            else
            {
                if (LowestCommonAnsestor(&nodes[species[i]], &nodes[species[j]]) == level)
                {
                    if (nodes[species[i]].refs.size() >= min_sub_species && nodes[species[j]].refs.size() >= min_sub_species)
                    {
                        graph[species[i]][species[j]] = 1;
                        graph[species[j]][species[i]] = 1;
                    }
                }
            }
        }
    }

    vector<int> v;
    if (level == "random")
    {
        vector<int> pool;
        for (int i = 0; i < species.size(); ++i)
        {
            if (nodes[species[i]].refs.size() >= min_sub_species)
                pool.push_back(species[i]);
        }

        vector<int> aux(pool.size());
        for (unsigned i = 0; i < aux.size(); ++i)
            aux[i] = i;

        for (int i = 0; i < 100; ++i)
        {
            v.resize(size);
            for (int j = 0; j < size; ++j)
            {
                swap(aux[j], aux[j + rand() % (aux.size() - j)]);
                v[j] = pool[aux[j]];
            }
            tuples.push_back(v);
        }
    }
    else
        GenTuple(v, 0);

    if (size >= 10)
    {
        for (unsigned i = 0; i < tuples.size(); ++i)
            swap(tuples[i], tuples[i + rand() % (tuples.size() - i)]);
    }

    for (unsigned i = 0; i < tuples.size(); ++i)
    {
//        if (i != 1)
//            continue;
//        if (i < 2)
//            continue;
//
        if (type == "meta")
            TestMeta(tuples[i]);
        else if (type == "control")
            TestControl(tuples[i]);
        else if (type == "metaidba")
            TestMetaIDBA(tuples[i]);
        else if (type == "velvet")
            TestMetaVelvet(tuples[i]);
        else if (type == "abyss")
            TestMetaAbyss(tuples[i]);
        else if (type == "tuples")
        {
            for (unsigned j = 0; j < tuples[i].size(); ++j)
                cout << tuples[i][j] << " ";
            for (unsigned j = 0; j < tuples[i].size(); ++j)
                cout << nodes[tuples[i][j]].refs.size() << " ";
            cout << endl;
        }
        else if (type == "compareKmer")
        {
            TestCompareKmer(tuples[i]);
        }
    }

    return 0;
}

