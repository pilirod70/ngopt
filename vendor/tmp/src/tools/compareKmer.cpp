#include "globals.h"
#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "HashGraph.h"
#include "SequenceReader.h"
#include "OptionSet.h"

#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

HashGraph hashGraph1;
HashGraph hashGraph2;
char line[MaxLine];
string comment1;
string comment2;
char buf[MaxLine];

void Usage()
{
    printf("Usage: compareKmer ref-file1 ref-file2\n");
}

int64 CommonKmer(HashGraph *graph1, HashGraph *graph2)
{
    int64 common = 0;
    for (HashGraph::Iterator iter = graph1->Begin(); iter != graph1->End(); ++iter)
    {
        KmerNode *node = *iter;
        if (graph2->GetNode(node->GetKmer()) != NULL)
            ++common;
    }
    return common;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        Usage();
        exit(1);
    }

    OptionSet option_set;
    option_set.AddOption("kmer", "k", kmerLength, "");
    option_set.ProcessOptions(argc, argv);

    FastaReader reader1(argv[1]);
    FastaReader reader2(argv[2]);

    Sequence ref1, ref2;

//    reader1.Read(ref1, comment1);
//    reader2.Read(ref2, comment2);
//
//    hashGraph1.InsertSequence(ref1);
    while (reader1.Read(ref1, comment1))
        hashGraph1.InsertSequence(ref1);

    //hashGraph2.InsertSequence(ref2);
    while (reader2.Read(ref2, comment2))
        hashGraph2.InsertSequence(ref2);

    int64 common = CommonKmer(&hashGraph1, &hashGraph2);

    printf("%.6f%%\n", 100.0 * 2 * common / (hashGraph1.NumNodes() + hashGraph2.NumNodes()));

    return 0;
}
