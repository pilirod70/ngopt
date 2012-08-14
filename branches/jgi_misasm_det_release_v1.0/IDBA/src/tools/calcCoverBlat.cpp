/**
 * @file calcCoverMummer.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "OptionSet.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"
#include "Utils.h"
#include "BlatRecord.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

char buf[MaxLine];
char line[MaxLine];
char name[MaxLine];
char ref_name[MaxLine];
vector<int64> recrod;
vector<int64> lengths;
vector<int64> valid_lengths;
set<string> valid_contigs;
map<string, int> dict;
int min_contig = 100;
double similar = 0.99;
vector<int> flags[10000];
int num_ref = 0;
OptionSet option_set;
int num_species = 1;

int main(int argc, char *argv[])
{
    option_set.AddOption("minContig", "", min_contig, "");
    option_set.AddOption("similar", "", similar, "");
    option_set.AddOption("num", "n", num_species, "");
    option_set.ProcessOptions(argc, argv);

    FastaReader ref_reader(argv[1]);
    FastaReader qry_reader(argv[2]);

    Sequence ref;
    string comment;
    int64 index = 0;
    while (ref_reader.Read(ref, comment))
    {
        sscanf(comment.c_str(), "%s", ref_name);
        dict[ref_name] = index;
        flags[index].resize(ref.Size(), 0);
        ++index;
    }
    num_ref = index;
    //cout << num_ref << endl;

    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);

    index = 0;
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        int ref_id = dict[record.ref_name];

        //if (record.match_count > similar * record.query_length && record.query_length >= min_contig)
        if (record.match_count > similar * record.query_length && record.query_length >= min_contig
                && record.match_count > similar * abs(record.ref_to - record.ref_from))
        {
            for (unsigned i = 0; i < record.blocks.size(); ++i)
            {
                BlatBlock block = record.blocks[i];
                for (unsigned j = block.ref_from; j < block.ref_from + block.size; ++j)
                    flags[ref_id][j] = 1;
            }

            valid_lengths.push_back(abs(record.query_to - record.query_from));
            valid_contigs.insert(record.query_name);
        }
    }

    int64 count = 0;
    int64 total = 0;
    for (int k = 0; k < num_ref; ++k)
    {
        for (unsigned i = 0; i < flags[k].size(); ++i)
        {
            if (flags[k][i])
                ++count;
            ++total;
        }
    }

    sort(valid_lengths.begin(), valid_lengths.end());
    reverse(valid_lengths.begin(), valid_lengths.end());
    int64 n50 = 0;
    int64 sum = 0;
    int64 n80 = 0;

    for (unsigned i = 0; i < valid_lengths.size(); ++i)
    {
        sum += valid_lengths[i];
        if (sum >= 0.5 * total && n50 == 0)
            n50 = valid_lengths[i];
        if (sum >= 0.8 * total && n80 == 0)
            n80 = valid_lengths[i];
    }

    int64 maximum = 0;
    int64 mean = 0;
    if (valid_lengths.size() > 0)
    {
        maximum = valid_lengths[0];
        mean = sum / valid_lengths.size();
    }

    FastaReader qry_reader2(argv[2]);
    int64 num_contigs = 0;
    long long sum_wrong = 0;
    int64 num_wrong = 0;
    int64 corret_contigs = 0;
    int64 sum_corret = 0;

    Sequence seq;
    while (qry_reader2.Read(seq, comment))
    {
        if (seq.Size() < min_contig)
            continue;

        sscanf(comment.c_str(), "%s", buf);
        comment = buf;
        if (valid_contigs.find(comment) == valid_contigs.end())
        {
            ++num_wrong;
            sum_wrong += seq.Size();
        }
        else
        {
            ++corret_contigs;
            sum_corret += seq.Size();
        }
    }

    printf("contigs: %lld N50: %lld coverage: %.2f%% max: %lld mean:%lld total: %lld/%lld N80: %lld\n", 
            (int64)valid_contigs.size(), n50, count * 100.0 / total,
            maximum, mean,
            count, total, n80);
    printf("error: %lld %lld correct: %lld %lld\n", num_wrong, sum_wrong, corret_contigs, sum_corret);

    return 0;
}

