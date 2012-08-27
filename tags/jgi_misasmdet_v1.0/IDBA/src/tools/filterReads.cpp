/**
 * @file filterReads.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "BlatRecord.h"
#include "Log.h"
#include "Read.h"
#include "Sequence.h"
#include "Utils.h"
#include "OptionSet.h"
#include "SequenceReader.h"

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
double similar = 0.9;
int min_cover = 1;
int dist[1000] = {0};
int min_contig = 0;
string names[100000];
double coverage[100000] = { 0 };
vector<int> valid_lengths;
vector<string> aligned_reads;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("similar", "", similar, "similarity of alignment");
    option_set.AddOption("minCover", "", min_cover, "minimum coverage");
    option_set.AddOption("minContig", "", min_contig, "minimum length");
    option_set.ProcessOptions(argc, argv);

    if (argc < 3)
    {
        fprintf(stderr, "analyzeAlignment ref-file blat-file read-file new-read-file\n");
        exit(1);
    }

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

        if (records.size() > 0)
            ++num_contigs;

        sscanf(records[0].query_name.c_str(), "%s", buf);
        aligned_reads.push_back(buf);
    }

    
    for (unsigned i = 0; i < aligned_reads.size(); ++i)
    {
        printf("%s\n", aligned_reads[i].c_str());
    }

    //FastAReader read_reader(argv[3]);
    //FastAWriter read_writer(argv[4]);

//    Sequence seq1, seq2;
//    string comment1, comment2;;
//    index = 0;
//    while (read_reader.Read(seq1, comment1) && read_reader.Read(seq2, comment2))
//    {
//        sscanf(comment1.c_str(), "%s", buf);
//
//        if (buf == aligned_reads[index])
//        {
//            read_writer.Write(seq1, comment1);
//            read_writer.Write(seq2, comment2);
//
//            ++index;
//        }
//    }

    return 0;
}

