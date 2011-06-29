#include "globals.h"

#include "Kmer.h"
#include "OptionSet.h"
#include "SequenceReader.h"
#include "Sequence.h"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <map>
#include <queue>
#include <stack>

using namespace std;

int num = 10000;
int fragment_length = 200;
int read_length = 5;
int dist[1000] = {0};

int main(int argc, char *argv[])
{
    OptionSet option_set;
    option_set.AddOption("num", "", num, "");
    option_set.AddOption("fragmentLength", "", fragment_length, "");
    option_set.AddOption("readLength", "", read_length, "");
    option_set.ProcessOptions(argc, argv);

    if (argc < 2)
    {
        fprintf(stderr, "Usage: eulerPathTest ref-file\n");
        exit(1);
    }

    kmerLength = read_length - 1;

    vector<Sequence> refs;
    Sequence seq;
    string comment;

    FastaReader reader(argv[1]);
    while (reader.Read(seq, comment))
    {
        refs.push_back(seq);
    }

    if (num < refs.size())
        refs.resize(num);

    int64 sum_dist = 0;
    int count = 0;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        map<Kmer, int> dict;
        Sequence seq = refs[i];
        Kmer kmer;
        int len = 0;
        for (int j = 0; j < seq.Size(); ++j)
        {
            kmer.AddRight(seq[j]);
            len++;
            if (len < kmerLength)
                continue;

            if (dict.find(kmer) != dict.end())
            {
                sum_dist += j - dict[kmer];
                ++count;
                int d = j - dict[kmer];
                if (d < 1000)
                    dist[d/10]++;
            }
            dict[kmer] = j;
        }
    }

    for (int i = 0; i < 100; ++i)
        cout << i * 10 << " " << dist[i] << endl;
    cout << sum_dist << " " << count << " " << sum_dist * 1.0 / count << endl;


    return 0;
}

