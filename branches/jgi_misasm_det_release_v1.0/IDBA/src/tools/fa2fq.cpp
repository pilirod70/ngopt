/**
 * @file fa2fq.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Log.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: fq2fa fq-file fa-file\n");
        exit(1);
    }

    FastaReader reader(argv[1]);
    FastqWriter writer(argv[2]);
    
    Sequence seq;
    string ss;
    int index = 0;
    while (reader.Read(seq, ss))
    {
        writer.Write(seq, ss);
        ++index;
    }
    printf("total %d reads\n", index);

    return 0;
}

