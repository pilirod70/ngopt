/**
 * @file calcCoverHuman.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Sequence.h"
#include "SequenceReader.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace std;

vector<int64> lengths;
set<string> valid_contigs;
char buf[MaxLine];
char name[MaxLine];
int min_contig = 100;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "usage: caclCoverHuman contig-file < stat-file\n");
        exit(1);
    }

    int64 sum_aligned = 0;
    int64 sum_ref_length = 0;

    for (int i = 0; i < 24; ++i)
    {
        int64 aligned;
        int64 ref_length;
        scanf("%lld %lld", &aligned, &ref_length);
        sum_aligned += aligned;
        sum_ref_length += ref_length;

        int num;
        scanf("%d", &num);

        for (int j = 0;j < num; ++j)
        {
            int len;
            scanf("%d", &len);
            lengths.push_back(len);
        }


        int64 num_valid_contigs;
        scanf("%lld", &num_valid_contigs);

        for (int i = 0; i < num_valid_contigs; ++i)
        {
            scanf("%s", buf);
            valid_contigs.insert(buf);
        }
    }


    sort(lengths.begin(), lengths.end());
    reverse(lengths.begin(), lengths.end());
    int64 n50 = 0;
    int64 sum = 0;
    for (unsigned i = 0; i < lengths.size(); ++i)
    {
        sum += lengths[i];
        if (n50 == 0 && sum > sum_ref_length*0.5)
            n50 = lengths[i];
    }

    int64 num_contigs = 0;
    int64 num_wrong = 0;;
    int64 num_wrong_bases = 0;
    FastaReader reader(argv[1]);
    Sequence seq;
    string comment;
    while (reader.Read(seq, comment))
    {
        if (seq.Size() >= min_contig)
            ++num_contigs;
        else
            continue;

        sscanf(comment.c_str(), "%s", name);
        if (valid_contigs.find(name) == valid_contigs.end())
        {
            ++num_wrong;
            num_wrong_bases += seq.Size();
        }
    }


    printf("contigs: %lld/%lld rate: %.2f%% n50: %lld max: %lld mean: %lld total: %lld coverage: %.2f%% wrong contigs: %lld, %lld\n",
            (int64)valid_contigs.size(), num_contigs, 100.0 * valid_contigs.size() / num_contigs,
            n50, lengths[0], sum / lengths.size(), 
            sum_ref_length, 100.0 * sum_aligned / sum_ref_length,
            num_wrong, num_wrong_bases);

    return 0;
}

