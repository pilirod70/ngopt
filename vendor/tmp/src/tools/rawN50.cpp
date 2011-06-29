/**
 * @file rawN50.cpp
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
#include "Utils.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>

using namespace std;

const int MaxContigs = 100000000;

Sequence contigs[MaxContigs];
int len[MaxContigs];
char line[MaxLine];
int minContig = 100;
int refLength = 0;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("minContig", "", minContig, "");
    option_set.AddOption("refLength", "", refLength, "");
    option_set.ProcessOptions(argc, argv);
///    AddParameter("minContig", &minContig, INTEGER);
///    AddParameter("refLength", &refLength, INTEGER);
///    ProcessParameters(argc, argv);

    if (argc < 2)
    {
        cerr << "usage: rawN50 contig-file" << endl;
        throw exception();
    }

    FastaReader reader(argv[1]);

    Sequence seq;
    string comment;
    unsigned index = 0;
    while (reader.Read(seq, comment))
    {
        {
            if (seq.Size() < minContig)
                continue;

            for (int i = 0; i < seq.Size(); ++i)
                seq[i] = toupper(seq[i]);

//             if (!seq.IsChar())
//                 continue;

            //contigs[index].SetContent(seq);
            contigs[index] = seq;
            ++index;
        }
    }

    int total = index;
    int64 sum = 0;
    for (int i = 0; i < total; ++i)
    {
        len[i] = contigs[i].Size();
        sum += len[i];
    }

    sort(len, len + total);
    reverse(len, len + total);

    int n50 = 0;
    int64 s = 0;

    if (refLength > 0)
        sum = refLength;
    int n80 = 0;
    for (int i = 0; i < total; ++i)
    {
        s += len[i];
        if (s > sum*0.5 && n50 == 0)
            n50 = len[i];
        if (s > sum*0.8 && n80 == 0)
            n80 = len[i];
    }

//    for (unsigned i = 0; i < total; ++i)
//    {
//        if (len[i] >= n50)
//            cout << len[i] << endl;
//    }

    printf("contigs: %d n50: %d max: %d mean: %lld total length: %lld n80: %d\n",
           total, n50, len[0], sum/total, sum, n80);

    return 0;
}

