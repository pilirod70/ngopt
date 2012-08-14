/**
 * @file filterContigs.cpp
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

using namespace std;

int min_contig = 100;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("minContig", "", min_contig, "");
    option_set.ProcessOptions(argc, argv);

//    AddParameter("minContig", &min_contig, INTEGER);
//    ProcessParameters(argc, argv);

    if (argc < 3)
    {
        fprintf(stderr, "filterContigs in-file out-file\n");
        exit(1);
    }

    FastaReader reader(argv[1]);
    FastaWriter writer(argv[2]);

    Sequence seq;
    string comment;

    int index = 0;
    while (reader.Read(seq, comment))
    {
        if (seq.Size() < min_contig)
            continue;
        else
        {
            writer.Write(seq, FormatString("contig%d %d", index++, comment.c_str()));
            //writer.WriteFormat(seq, "contig%d %s", index++, comment.c_str());
            //writer.Write(seq, comment);
        }
    }

    return 0;
}

