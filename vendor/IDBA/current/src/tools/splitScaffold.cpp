/**
 * @file splitScaffold.cpp
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

    //int index = 0;
    Sequence seq;
    string comment;
    while (reader.Read(seq, comment))
    {
        vector<Sequence> sub_seqs;

        int last = 0;
        for (int i = 0; i < seq.Size(); ++i)
        {
            if (seq[i] >= 0 && seq[i] < 4)
            {
            }
            else
            {
                if (i > last)
                {
                    Sequence sub;
                    sub.Assign(seq, last, i - last);
                    sub_seqs.push_back(sub);
                }

                last = i + 1;
            }
        }
        
        int i = seq.Size();
        if (i > last)
        {
            Sequence sub;
            sub.Assign(seq, last, i - last);
            sub_seqs.push_back(sub);
        }

        int index = 0;
        for (unsigned i = 0; i < sub_seqs.size(); ++i)
        {
            writer.Write(sub_seqs[i], FormatString("%s_%d", comment.c_str(), index++));
            //writer.WriteFormat(t, "%s_%d", comment.c_str(), index++);
        }
    }

    return 0;
}

