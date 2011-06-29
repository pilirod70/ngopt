/**
 * @file normReads.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Log.h"
#include "Sequence.h"
#include "Utils.h"
#include "OptionSet.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace std;


char line[MaxLine];
char comment[MaxLine];
char comment2[MaxLine];
FILE *fp[MaxLine];
int length = 0;
bool mate = 0;

OptionSet option_set;

int main(int argc, char *argv[])
{
    option_set.AddOption("length", "", length, "");
    option_set.AddOption("mate", "", mate, "");

    option_set.ProcessOptions(argc, argv);

    if (argc < 3
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: normReads fa-file norm-fa-file\n");
        fprintf(stderr, "       [--length l] [--mate]\n");
        exit(1);
    }

    FastaReader reader(argv[1]);
    FastaWriter writer(argv[2]);

    Sequence seq;
    Sequence seq2;
    string comment;
    string comment2;
    unsigned index = 0;

    if (mate)
    {
        while (reader.Read(seq, comment))
        {
            if (!reader.Read(seq2, comment2))
                break;

            if (!seq.IsValid() || !seq2.IsValid())
                continue;

            if (length == 0)
            {
                writer.Write(seq, FormatString("read%d/1", index));
                writer.Write(seq2, FormatString("read%d/2", index));
                index += 2;
            }
            else
            {
                if (seq.Size() >= length && seq2.Size() >= length)
                {
                    seq.Resize(length);
                    seq2.Resize(length);
                    writer.Write(seq, FormatString("read%d/1", index));
                    writer.Write(seq2, FormatString("read%d/2", index));
                    index += 2;
                }
            }
        }
    }
    else
    {
        while (reader.Read(seq, comment))
        {
            if (length == 0 && seq.IsValid())
            {
                writer.Write(seq, comment);
                ++index;
            }
            else
            {
                if (seq.Size() >= length && seq.IsValid())
                {
                    seq.Resize(length);
                    writer.Write(seq, comment);
                    ++index;
                }
            }
        }
    }

    return 0;
}

