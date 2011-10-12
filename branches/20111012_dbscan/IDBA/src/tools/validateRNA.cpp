/**
 * @file validateRNA.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "BlatRecord.h"
#include "Log.h"
#include "OptionSet.h"
#include "Read.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"
#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

using namespace std;

struct Interval
{
    int from, to;
};

char line[MaxLine];
char buf[MaxLine];
char tmp[MaxLine];
map<string, int> id_map;
vector<Sequence> refs;
vector<vector<int> > matchs;
vector<vector<BlatRecord> > aligned_records;
vector<int> flags;
vector<double> depths;
double similar = 0.99;
int min_cover = 1;
int dist[1000] = {0};
int min_contig = 100;
string names[100000];
double coverage[100000] = { 0 };
vector<int> valid_lengths;
double cutoff = 0;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("similar", "", similar, "similarity of alignment");
    option_set.AddOption("minCover", "", min_cover, "minimum coverage");
    option_set.AddOption("minContig", "", min_contig, "minimum length");
    option_set.AddOption("cutoff", "", cutoff, "cutoff of the expression level");
    option_set.ProcessOptions(argc, argv);

    if (argc < 4)
    {
        fprintf(stderr, "validateRNA ref-file blat-file depth-file\n");
        exit(1);
    }

    FastaReader reader(argv[1]);
    Sequence seq;
    string comment;
    int index = 0;

//    refs.resize(reader.NumReads());
//    matchs.resize(refs.size());
//    flags.resize(refs.size());

    while (reader.Read(seq, comment))
    {
        refs.resize(refs.size() + 1);
        matchs.resize(matchs.size() + 1);
        flags.resize(flags.size() + 1);
        sscanf(comment.c_str(), "%s", buf);
        names[index] = buf;
        refs[index] = seq;
        matchs[index].resize(seq.Size() + 1, 0);

//        Replace(comment, '_', ' ');
//        sscanf(comment.c_str(), "%s %s", tmp, buf);
        id_map[buf] = index++;
    }

    depths.resize(matchs.size(), 100);
    aligned_records.resize(matchs.size());

    FILE *fdepth = fopen(argv[3], "rb");
    fill(depths.begin(), depths.end(), 0);
    for (unsigned i = 0; i < depths.size(); ++i)
    {
        double n, d;
        fscanf(fdepth, "%s %lf %lf", buf, &n, &d);
        d *= 75;
        depths[id_map[buf]] = d;
    }

    cout << "initialized" << endl;

    FILE *fblat = fopen(argv[2], "rb");
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);

    int num_contigs = 0;
    int total_contigs = 0;
    string prev_name;
    while (true)
    {
        BlatRecord record;
        if (fgets(line, MaxLine, fblat) == NULL)
            break;

        record.Parse(line);
        vector<BlatRecord> records;
        records.push_back(record);
        while (fgets(line, MaxLine, fblat) != NULL)
        {
            record.Parse(line);
            if (record.query_name != records.back().query_name)
            {
                fseek(fblat, -strlen(line), SEEK_CUR);
                break;
            }
            else
                records.push_back(record);
        }

        string s = records[0].query_name;
        Replace(s, '_', ' ');
        double level;
        sscanf(s.c_str(), "%s %lf", buf, &level);
        if (level < 0)
            level = 0;
        if (level < cutoff)
            continue;

        ++total_contigs;

        int index = 0;
        for (unsigned i = 0; i < records.size(); ++i)
        {
            BlatRecord &record = records[i];

            if (record.query_length < min_contig)
                continue;
            else if (record.match_count < record.query_length * similar)
                continue;
            else
                records[index++] = record;
        }
        records.resize(index);

        if (records.size() == 0)
            continue;

        ++num_contigs;

        for (unsigned i = 0; i < records.size(); ++i)
        {
            BlatRecord &record = records[i];

            valid_lengths.push_back(record.query_length);

            int id = id_map[record.ref_name];
            if (record.match_count >= similar * record.ref_length)
            {
                //int id = id_map[record.ref_name];
                flags[id] = true;
            }

            aligned_records[id].push_back(record);

            for (unsigned i = 0; i < record.blocks.size(); ++i)
            {
                BlatBlock &block = record.blocks[i];
                int id = id_map[record.ref_name];

                for (int j = block.ref_from; j < block.ref_from + block.size && j < (int)matchs[id].size(); ++j)
                {
                    matchs[id][j] += 1;//.0 / records.size();
                }
            }
        }
    }

    int aligned = 0;
    int aligned_length = 0;
    int total = 0;
    int total_length = 0;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        int num_match = 0;
        for (unsigned j = 0; j < matchs[i].size(); ++j)
        {
            if (matchs[i][j] >= min_cover)
                num_match++;
        }

        vector<BlatRecord> &records = aligned_records[i];
        sort(records.begin(), records.end());
        int len = records.size();
        if (len > 100)
            len = 100;

        //cout << i << endl;

        int maximum_match = 0;
        for (int a = 0; a < len; ++a)
        {
            BlatRecord &record = records[a];
            for (int b = a+1; b < len; ++b)
            {
                BlatRecord &record2 = records[b];

                vector<int> tmp(matchs[i].size(), 0);

                for (unsigned x = 0; x < record.blocks.size(); ++x)
                {
                    BlatBlock &block = record.blocks[x];
                    int id = id_map[record.ref_name];

                    for (int y = block.ref_from; y < block.ref_from + block.size && y < (int)matchs[id].size(); ++y)
                    {
                        tmp[y]++;
                    }
                }

                for (unsigned x = 0; x < record2.blocks.size(); ++x)
                {
                    BlatBlock &block = record2.blocks[x];
                    int id = id_map[record2.ref_name];

                    for (int y = block.ref_from; y < block.ref_from + block.size && y < (int)matchs[id].size(); ++y)
                    {
                        tmp[y]++;
                    }
                }

                int sum = 0;
                for (unsigned x = 0; x < tmp.size(); ++x)
                {
                    if (tmp[x])
                        sum++;
                }

                if (sum > maximum_match)
                    maximum_match = sum;
            }
        }

        if (depths[i] >= 30)
        {
//            if (flags[i] == 0 && maximum_match > similar * matchs[i].size())
//                cout << names[i] << endl;

            //if (flags[i] || maximum_match > similar * matchs[i].size())
            if (flags[i])// || num_match > similar * matchs[i].size())
            {
                ++aligned;
                aligned_length += refs[i].Size();

                //cout << depths[i] << " " << endl;
            }
            total++;
            total_length += refs[i].Size();
        }
    }

    string s = argv[2];
    s.resize(s.size() - 5);

    total_contigs = 0;
    FastaReader contig_reader(s);
    while (contig_reader.Read(seq, comment))
    {
        string s = comment;
        Replace(s, '_', ' ');
        double level;
        sscanf(s.c_str(), "%s %lf", buf, &level);
        if (level < 0)
            level = 0;
        if (level < cutoff)
            continue;
        
        if (seq.Size() < min_contig)
            continue;

        ++total_contigs;
    }


    printf("sensitivity %.2f%% %d/%d %d/%d\n", aligned * 100.0 / total, aligned, total, aligned_length, total_length);
    printf("precision %.2f%% %d/%d\n", 100.0 * num_contigs / total_contigs, num_contigs, total_contigs);
    
///    int maximum = 0;
///    for (unsigned i = 0; i < aligned_records.size(); ++i)
///    {
///        cout << aligned_records[i].size() << endl;
///        if (aligned_records[i].size() > maximum)
///            maximum = aligned_records[i].size();
///    }
///    cout << maximum << endl;

//    for (unsigned i = 0; i < flags.size(); ++i)
//    {
//        cout << names[i] << " " << flags[i] << endl;
//    }

    return 0;
}


