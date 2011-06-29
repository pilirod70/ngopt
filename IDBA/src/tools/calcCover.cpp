/**
 * @file calcCover.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Log.h"
#include "OptionSet.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"
#include "Utils.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

struct Interval
{
    int from, to;
};

int minContig = 100;

const int MaxSeqLength = 10000000;

char buf[MaxLine];
char line[MaxLine];
char readName[MaxLine];
char refName[MaxLine];
char comment[MaxLine];
vector<string> names;
vector<Sequence> sequences;
map<string, int> dict;
vector<char *> matchs;
vector<vector<int> > lengths;
vector<set<string> > founds;
set<string> allFound;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("minContig", "", minContig, "");
    option_set.ProcessOptions(argc, argv);
//    AddParameter("minContig", &minContig, INTEGER);
//
//    ProcessParameters(argc, argv);

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: calcCover contig-file ref-file\n");
        fprintf(stderr, "       [--minContig l]\n");
        exit(1);
    }

    FastaReader contigs_reader(argv[1]);
    FastaReader reference_reader(argv[2]);
    Sequence seq;
    string comment;
    while (reference_reader.Read(seq, comment))
    {
        sscanf(comment.c_str(), "%s", refName);
        dict[refName] = sequences.size();
        sequences.push_back(seq);
        char *match = new char[seq.Size()];
        fill_n(match, seq.Size(), 0);
        matchs.push_back(match);
        lengths.push_back(vector<int>());
        names.push_back(refName);
        founds.push_back(set<string>());
    }

    set<string> found;
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        if (line[0] == '-' && line[1] == '-')
            break;

        int from, to;
        sscanf(line, "%s %s %s %s %s %s %d %d",
               buf, readName, buf, buf, buf, refName, &from, &to);

        if (from > to)
            swap(from, to);

        if (to - from >= minContig)
        {
            if (dict.find(refName) == dict.end())
            {
                LogError("Input contain unknown sequence\n");
                exit(1);
            }

            int index = dict[refName];
            founds[index].insert(readName);
            allFound.insert(readName);
            lengths[index].push_back(to - from);
            for (int i = from; i < to; ++i)
                matchs[index][i] = 1;
        }
    }

    uint64 sumAlign = 0;
    uint64 sumLength = 0;
    uint64 sumCount = 0;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        vector<int> &v = lengths[i];
        int n50 = 0;
        int sum = 0;
        int n80 = 0;
        int maximum = 0;
        sort(v.begin(), v.end());
        reverse(v.begin(), v.end());
        for (unsigned j = 0; j < v.size(); ++j)
        {
            sum += v[j];
            if (v[j] > maximum)
                maximum = v[j];
            if (sum > 0.5 * sequences[i].Size() && n50 == 0)
                n50 = v[j];
            if (sum > 0.8 * sequences[i].Size() && n80 == 0)
                n80 = v[j];
        }

        int total = 0;
        int contigs = 0;
        char *match = matchs[i];
        for (int j = 0; j < sequences[i].Size(); ++j)
        {
            total += match[j];
            if (match[j] && j+1 < sequences[i].Size() && match[j+1] == 0)
                ++contigs;
        }

        {
        printf("%s\nfound: %lu n50: %d n80: %d max: %d mean: %d\n",
               names[i].c_str(), founds[i].size(), n50, n80, maximum, sum/max((int)v.size(), 1));
        printf("coverage: %.2f%% merged contigs: %d mapped length: %d %d ref length: %d\n",
               total*100.0/sequences[i].Size(), contigs, sum, total, sequences[i].Size());
        }

        sumAlign += total;
        sumLength += sequences[i].Size();
        sumCount += (founds[i].size() != 0);
    }

    printf("coverage %.2f%%\n", sumAlign * 100.0 / sumLength);

    int wrong = 0;
    int wrongLength = 0;
    while (contigs_reader.Read(seq, comment))
    {
        sscanf(comment.c_str(), "%s", buf);
        comment = buf;

        if (seq.Size() >= minContig && allFound.find(comment) == allFound.end())
        {
            ++wrong;
            wrongLength += seq.Size();
        }
    }
     printf("wrong: %d wrong length: %d\n", wrong, wrongLength);

    return 0;
}

