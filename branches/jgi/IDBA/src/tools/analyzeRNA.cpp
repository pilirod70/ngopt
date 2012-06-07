/**
 * @file analyzeRNA.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "BlatRecord.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "Utils.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>

using namespace std;

char buf[MaxLine];
char line[MaxLine];
map<string, int> ref_table;
map<string, int> read_table;
vector<Sequence> refs;
vector<vector<int> > matchs;
vector<Sequence> reads;
char rev_comp_table[256];
int exon_id[200000000];
//vector<int> exon_id(200000000);
vector<int> exon_from;
vector<int> exon_to;
int threshold = 10;

bool IsEnd(int length, int from, int to)
{
    if (from <= threshold || to >= length - threshold)
        return true;
    else
        return false;
}

int main(int argc, char *argv[])
{
    exon_from.resize(200000000);
    exon_to.resize(200000000);

    rev_comp_table['A'] = 'T';
    rev_comp_table['C'] = 'G';
    rev_comp_table['G'] = 'C';
    rev_comp_table['T'] = 'A';
    rev_comp_table['N'] = 'N';

    if (argc < 4)
    {
        fprintf(stderr, "analyzeRNA genome.fa read.fa read.blat\n");
        exit(1);
    }

    Sequence seq;
    string comment;
    FastaReader ref_reader(argv[1]);
    //refs.resize(ref_reader.NumReads());
    //matchs.resize(refs.size());
    int index = 0;
    while (ref_reader.Read(seq, comment))
    {
        refs.resize(refs.size() + 1);
        matchs.resize(matchs.size() + 1);
        cout << comment << " " << seq.Size() << endl;
        ref_table[comment] = index;
        refs[index] = seq;
        matchs[index].resize(seq.Size());
        ++index;
    }

    FastaReader read_reader(argv[2]);
    //reads.resize(read_reader.NumReads());
    index = 0;
    while (read_reader.Read(seq, comment))
    {
        reads.resize(reads.size() + 1);
        read_table[comment] = index;
        reads[index] = seq;
        ++index;
    }

    cout << ref_table.size() << " " << read_table.size() << endl;

    FILE *fblat = fopen(argv[3], "rb");
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);

    index = 0;
    string last_query;
    while (fgets(line, MaxLine, fblat) != NULL)
    {
        if (index % 100000 == 0)
            cout << index << endl;
        ++index;

        bool is_valid = false;
        for (int i = 0; line[i]; ++i)
        {
            if (line[i] == '\t')
                is_valid = true;
        }

        if (!is_valid)
            continue;

        //printf("%s", line);
        BlatRecord record;
        //cout << "start" << endl;
        record.Parse(line);
        //cout << "finish" << endl;
        

        if (read_table.find(record.query_name) == read_table.end())
            continue;
        if (ref_table.find(record.ref_name) == ref_table.end())
            continue;

        if (last_query == record.query_name)
            continue;
        if (record.query_to - record.query_from < record.query_length * 0.9)
            continue;
        last_query = record.query_name;

        //int query_id = read_table[record.query_name];
        int ref_id = ref_table[record.ref_name];

        for (unsigned i = 0; i < record.blocks.size(); ++i)
        {
            if (record.blocks[i].size < 25)
                continue;

            for (int j = 0; j < record.blocks[i].size; ++j)
            {
                ++matchs[ref_id][record.blocks[i].ref_from + j];
            }
        }
    }

    fill_n(exon_id, matchs[0].size(), -1);
    for (unsigned i = 0; i < matchs.size(); ++i)
    {
        int offset = 0;
        int exon_count = 0;
        int exon_num_base = 0;
        int index = 0;
        FILE *fatt = fopen(argv[5], "wb");
        while (offset < (int)matchs[i].size())
        {
            if (matchs[i][offset] > 0)
            {
                int from = offset;
                int to = from;
                int sum_matchs = 0;
                while (matchs[i][to] > 0 && to < (int)matchs[i].size())
                {
                    sum_matchs += matchs[i][to];
                    ++to;
                }

                if (to - from < 25)
                {
                    cout << to -from << endl;
                    for (int j = from; j < to; ++j)
                        matchs[i][j] = 0;
                }
                else
                {
                    if (sum_matchs >= 2 * (to - from))
                    {
                        fprintf(fatt, "%d %d\n", index, to - from);
                        ++exon_count;
                        exon_num_base += to - from;
                        for (int j = from; j < to; ++j)
                        {
                            exon_id[j] = index;
                            exon_from[j] = from;
                            exon_to[j] = to;
                        }
                        ++index;
                    }
                }

                offset = to;
            }
            else
                ++offset;
        }

        int num_base = 0;
        int num_matchs = 0;
        int num[10] = {0};

        for (unsigned j = 0; j < matchs[i].size(); ++j)
        {
            if (matchs[i][j])
            {
                ++num_base;
                int m = min(9, matchs[i][j]);
                for (int k = 1; k <= m; ++k)
                    ++num[k];
            }
            num_matchs += matchs[i][j];
        }

        cout << exon_count << " " << exon_num_base << endl;
        cout << num_base << " " << num_matchs << endl;
        for (int j = 0; j < 10; ++j)
            cout << num[j] << " ";
        cout << endl;
    }

    fseek(fblat, 0, SEEK_SET);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);

    index = 0;
    FILE *fedge = fopen(argv[4], "wb");
    int edge_count = 0;
    while (fgets(line, MaxLine, fblat) != NULL)
    {
//        if (index % 100000 == 0)
//            cout << index << endl;
        ++index;

        bool is_valid = false;
        for (int i = 0; line[i]; ++i)
        {
            if (line[i] == '\t')
                is_valid = true;
        }

        if (!is_valid)
            continue;

        BlatRecord record;
        record.Parse(line);

        if (read_table.find(record.query_name) == read_table.end())
            continue;
        if (ref_table.find(record.ref_name) == ref_table.end())
            continue;

        if (last_query == record.query_name)
            continue;
        if (record.query_to - record.query_from < record.query_length * 0.9)
            continue;
        last_query = record.query_name;

        //int query_id = read_table[record.query_name];
        //int ref_id = ref_table[record.ref_name];

        for (unsigned i = 0; i+1 < record.blocks.size(); ++i)
        {
            if (record.blocks[i].size < 25 || record.blocks[i+1].size < 25)
                continue;

            int from = record.blocks[i].ref_from;
            int to = record.blocks[i+1].ref_from;

            int from_id = exon_id[from];
            int to_id = exon_id[to];

            if (from_id != -1 && to_id != -1)
            {
                int length1 = exon_to[from] - exon_from[from];
                int from1 = from - exon_from[from];
                int to1 = from - exon_from[from] + record.blocks[i].size;

                int length2 = exon_to[to] - exon_from[to];
                int from2 = to - exon_from[to];
                int to2 = to - exon_from[to] + record.blocks[i+1].size;

                if (IsEnd(length1, from1, to1) && IsEnd(length2, from2, to2))
                    cout << length1 << " " << from1 << " " << to1 << " " << length2 << " " << from2 << " " << to2 << endl;
                else
                    continue;


                fprintf(fedge, "%d %d\n", from_id, to_id);
                ++edge_count;
            }
        }
    }
    fclose(fedge);

    cout << edge_count << endl;

    return 0;
}

