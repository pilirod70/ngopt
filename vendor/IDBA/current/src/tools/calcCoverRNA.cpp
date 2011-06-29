/**
 * @file calcCoverRNA.cpp
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
double similar = 0.95;
int min_cover = 10;
int dist[1000] = {0};
int min_contig = 0;
string names[100000];
double coverage[100000] = { 0 };
vector<int> valid_lengths;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("similar", "", similar, "similarity of alignment");
    option_set.AddOption("minCover", "", min_cover, "minimum coverage");
    option_set.AddOption("minContig", "", min_contig, "minimum length");
    option_set.ProcessOptions(argc, argv);

    if (argc < 4)
    {
        fprintf(stderr, "analyzeAlignment ref-file blat-file cover-file\n");
        exit(1);
    }

    FILE *fp = fopen("found.list", "rb");
    vector<int> flags;
    int x;
    int num = 0;
    while (fscanf(fp, "%s %d", buf, &x) != EOF)
    {
        flags.push_back(x);
        if (x)
            ++num;
    }

    map<string, int> dict;
    FILE *fname = fopen("name.list", "rb");
    while (fscanf(fname, "%s", buf) != EOF)
    {
        dict[buf] = 1;
        //cout << buf << endl;
    }

    cerr << flags.size() << " " << num << endl;

    FastaReader reader(argv[1]);
    Sequence seq;
    string comment;
    int index = 0;

//    refs.resize(reader.NumReads());
//    matchs.resize(refs.size());

    while (reader.Read(seq, comment))
    {
        refs.resize(refs.size() + 1);
        matchs.resize(matchs.size() + 1);
        sscanf(comment.c_str(), "%s", buf);
        names[index] = buf;
        refs[index] = seq;
        matchs[index].resize(seq.Size() + 1, 0);

        id_map[buf] = index++;
    }

    cout << "initialized " << refs.size() << endl;

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

        ++total_contigs;

        if (total_contigs % 100000 == 0)
            cerr << total_contigs << endl;

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

        if (records.size() > 0)
            ++num_contigs;

        for (unsigned i = 0; i < records.size(); ++i)
        {
            BlatRecord &record = records[i];

            valid_lengths.push_back(record.query_length);

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

    cout << "parse finish" << endl;

    int64 total = 0;
    int64 total_ref_bases = 0;
    int64 num_ref = 0;
    FILE *fcover = fopen(argv[3], "wb");
    for (unsigned i = 0; i < matchs.size(); ++i)
    {
        int total_cover = 0;
        int num_match = 0;
        int last = 0;
        vector<Interval> intervals;
        for (unsigned j = 0; j < matchs[i].size(); ++j)
        {
            if (matchs[i][j] >= min_cover)
            {
                if (matchs[i][last] < min_cover)
                    last = j;
                    
                ++num_match;
                total_cover += matchs[i][j];
            }
            else
            {
                if (j - last >= 10)
                {
                    Interval inter;
                    inter.from = last;
                    inter.to = j;
                    intervals.push_back(inter);
                }

                last = j;
            }

            if (matchs[i][j] < 1000)
                dist[(int)matchs[i][j]]++;
        }

        int maximum = 0;
        int inter_index = 0;
        for (unsigned a = 0; a < intervals.size(); ++a)
        {
            int x = intervals[a].to - intervals[a].from;
            if (x > maximum)
            {
                maximum = x;
                inter_index = a;
            }
        }

        //cout << maximum << " " << matchs[i].size() << endl;

        if (flags[i])
        {
//            cout << flags[i] << " " << matchs[i].size() << " " << maximum  << " " << 1.0 * maximum / matchs[i].size()
//                << " " << num_match << " " << 1.0 * num_match / matchs[i].size() << " " << total_cover / (num_match+1) << endl;
            total += num_match;
            ++num_ref;
            total_ref_bases += matchs[i].size();

            coverage[i] = total_cover / matchs[i].size();
            coverage[i] = total_cover / (num_match + 1);
        }

        if (dict[names[i]])
        {
            cout << "mRNA" << i << " " << matchs[i].size() << endl;
            for (unsigned j = 0; j < matchs[i].size(); ++j)
                cout << matchs[i][j] << ",";
            cout << endl;
        }
    }

    sort(valid_lengths.begin(), valid_lengths.end());
    reverse(valid_lengths.begin(), valid_lengths.end());

    int64 n50 = 0;
    int64 sum = 0;
    for (unsigned i = 0; i < valid_lengths.size(); ++i)
    {
        sum += valid_lengths[i];
        if (n50 == 0 && sum > 0.5 * total_ref_bases)
            n50 = valid_lengths[i];
    }

    printf("contigs: %d/%d total basese: %lld\n", num_contigs, total_contigs, total);
    printf("N50: %lld total refs: %lld/%lld\n", n50, num_ref, total_ref_bases);
    printf("Coverage: %.4f%%\n", total*100.0 / total_ref_bases);

    for (int i = 0; i < (int)matchs.size(); ++i)
    {
        fprintf(fcover, "%s %.2f %.2f\n", names[i].c_str(), coverage[i], coverage[i]/75);
    }

    return 0;
}

