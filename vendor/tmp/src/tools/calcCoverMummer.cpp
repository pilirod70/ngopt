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

int main(int argc, char *argv[])
{
    option_set.AddOption("minContig", "", min_contig, "");
    option_set.AddOption("similar", "", similar, "");
    option_set.ProcessOptions(argc, argv);
//    AddParameter("minContig", &min_contig, INTEGER);
//    AddParameter("similar", &similar, FLOAT);

    //ProcessParameters(argc, argv);

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
    cout << num_ref << endl;

    index = 0;
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        if (line[0] == '>')
        {
            int64 ref_len;
            int64 contig_len;
            sscanf(line + 1, "%s %s %lld %lld", ref_name, name, &ref_len, &contig_len);
            int ref_id = dict[ref_name];

            ++index;
            while (fgets(line, MaxLine, stdin) != NULL)
            {
                if (line[0] == '>')
                {
                    fseek(stdin, -strlen(line), SEEK_CUR);
                    break;
                }

                int64 ref_from, ref_to;
                int64 contig_from, contig_to;
                int64 errors;
                sscanf(line, "%lld %lld %lld %lld %lld", &ref_from, &ref_to, &contig_from, &contig_to, &errors);

                --ref_from;
                --ref_to;
                --contig_from;
                --contig_to;

                if (ref_from > ref_to)
                    swap(ref_from, ref_to);
                ++ref_to;

                if (contig_from > contig_to)
                    swap(contig_from, contig_to);
                ++contig_to;

                int64 ref_match_len = ref_to - ref_from;
                int64 contig_match_len = contig_to - contig_from;

                if (ref_match_len >= similar * contig_len
                        && contig_len >= ref_match_len * similar
                        && errors <= contig_len * (1 - similar)
                        && contig_len >= min_contig)
                {
                    for (int64 i = ref_from; i < ref_to; ++i)
                        flags[ref_id][i] = 1;

                    valid_lengths.push_back(contig_match_len);
                    valid_contigs.insert(name);
                }

                while (fgets(line, MaxLine, stdin) != NULL)
                {
                    if (line[0] == '0')
                        break;
                }
            }
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
    //qry_reader.Rewind();
    int64 num_contigs = 0;
    long long sum_wrong = 0;
    int64 num_wrong = 0;

    sprintf(buf, "%s-e", argv[2]);
    FastaWriter writer(buf);
    Sequence seq;
    while (qry_reader2.Read(seq, comment))
    {
        if (seq.Size() >= min_contig)
        {
            ++num_contigs;

            sscanf(comment.c_str(), "%s", name);
            if (valid_contigs.find(name) == valid_contigs.end())
            {
                ++num_wrong;
                sum_wrong += seq.Size();
                writer.Write(seq, comment);
            }
            else
            {
            }
        }
    }

    printf("contigs: %lld N50: %lld coverage: %.2f%% max: %lld mean:%lld total: %lld/%lld N80: %lld\n", 
            (int64)valid_contigs.size(), n50, count * 100.0 / total,
            maximum, mean,
            count, total, n80);
    printf("error: %lld %lld\n", num_wrong, sum_wrong);

    sprintf(buf, "%s.stat", argv[2]);
    FILE *fp = fopen(buf, "ab");
    fprintf(fp, "%lld %lld\n", count, total);

    fprintf(fp, "%lld\n", (int64)valid_lengths.size());
    for (unsigned i = 0; i < valid_lengths.size(); ++i)
        fprintf(fp, "%lld\n", valid_lengths[i]);

    fprintf(fp, "%lld\n", (int64)valid_contigs.size());
    set<string>::iterator iter = valid_contigs.begin();
    while (iter != valid_contigs.end())
    {
        fprintf(fp, "%s\n", iter->c_str());
        ++iter;
    }

    return 0;
}

