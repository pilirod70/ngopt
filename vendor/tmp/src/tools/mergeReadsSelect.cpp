/**
 * @file mergeReadsSelect.cpp
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
#include "Utils.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 4
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: mergeReads read-file1 read-file2 merge-read-file\n");
        exit(1);
    }

    FastqReader reader1(argv[1]);
    FastqReader reader2(argv[2]);
    FastaWriter writer(argv[3]);

    Sequence seq1;
    Sequence seq2;
    string comment1;
    string comment2;

    int64 index = 0;
    while (reader1.Read(seq1, comment1) && reader2.Read(seq2, comment2))
    {
//        writer.Write(seq1, comment1);
//        writer.Write(seq2, comment2);

        writer.Write(seq1, FormatString("read%d/1", index));
        writer.Write(seq2, FormatString("read%d/2", index));
//        writer.WriteFormat(seq1, "read%d/1", index);
//        writer.WriteFormat(seq2, "read%d/2", index);
        ++index;
    }

    cout << index << " pairs." << endl;

    return 0;
}

