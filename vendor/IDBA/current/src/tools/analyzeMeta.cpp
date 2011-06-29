/**
 * @file analyzeMeta.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-11-07
 */

#include "globals.h"

#include "BlatRecord.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"
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
vector<string> names;
vector<Sequence> refs;
vector<vector<int> > matchs;
vector<Sequence> reads;
char rev_comp_table[256];
int exon_id[200000000];
vector<int> exon_from;
vector<int> exon_to;
int threshold = 10;
int aux[10000];
double coverage[10000];

bool Compare(int x, int y)
{
    return coverage[x] > coverage[y];
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

    if (argc < 5)
    {
        fprintf(stderr, "analyzeMeta genome.fa read.fa read.blat high-contig.fa\n");
        exit(1);
    }

    Sequence seq;
    string comment;
    FastaReader ref_reader(argv[1]);
//    refs.resize(ref_reader.NumReads());
//    matchs.resize(refs.size());
//    names.resize(refs.size());
    int index = 0;
    while (ref_reader.Read(seq, comment))
    {
        refs.resize(refs.size() + 1);
        matchs.resize(matchs.size() + 1);
        names.resize(names.size() + 1);
        sscanf(comment.c_str(), "%s", buf);

        ref_table[buf] = index;
        refs[index] = seq;
        matchs[index].resize(seq.Size());
        //names[index] = buf;
        names[index] = comment;
        ++index;
    }

//    FastAReader read_reader(argv[2]);
//    reads.resize(read_reader.NumReads());
//    index = 0;
//    while (read_reader.Read(seq, comment))
//    {
//        read_table[comment] = index;
//        reads[index] = seq;
//        ++index;
//    }

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

        BlatRecord record;
        record.Parse(line);

        if (record.query_to - record.query_from < record.query_length * 0.9)
            continue;
        last_query = record.query_name;

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

    FastaWriter writer(argv[4]);
    sprintf(buf, "%s.cover", argv[4]);
    FILE *fcover = OpenFile(buf, "wb");
    for (unsigned i = 0; i < matchs.size(); ++i)
    {
        int count = 0;
        int sum = 0;
        for (unsigned j = 0; j < matchs[i].size(); ++j)
        {
            if (matchs[i][j])
                ++count;
            sum += matchs[i][j];
        }

        //if (count > 0.9 * matchs[i].size())
        {
            coverage[i] = 1.0 * count / matchs[i].size();
            cout << names[i] << "," << count << "," << sum << "," << matchs[i].size()  << "," << 1.0 * sum / count << "," << 1.0 * count / matchs[i].size() << "," << 1.0 * sum / matchs[i].size() << endl;
        }

        int from = 0;
        count = 0;
        sum = 0;
        matchs[i].push_back(0);
        int index = 0;
        for (unsigned j = 0; j < matchs[i].size(); ++j)
        {
            if (matchs[i][j])
            {
                count++;
                sum += matchs[i][j];
            }
            else
            {
                if (count > 1000 && sum > count * 30)
                {
                    Sequence seq;
                    //refs[i].GetSubSequence(seq, from, count);
                    seq.Assign(refs[i], from, count);
                    writer.Write(seq, FormatString("%s_%d", names[i].c_str(), index++));

                    for (int k = 0; k < count; ++k)
                    {
                        fprintf(fcover, "%d,", matchs[i][k+from]);
                    }
                    fprintf(fcover, "\n");
                }

                count = 0;
                sum = 0;
                from = j + 1;
            }
        }
    }

    for (unsigned i = 0; i < matchs.size(); ++i)
        aux[i] = i;
    sort(aux, aux + matchs.size(), Compare);

    FastaWriter ref_writer("ecoli-ref.fa");
    for (int i = 0; i < 11; ++i)
    {
        Sequence seq = refs[aux[i]];
        ref_writer.Write(seq, names[aux[i]]);
    }

    return 0;
}

