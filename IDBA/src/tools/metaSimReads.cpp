/**
 * @file simReads.cpp
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
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;

enum { EQUAL, UNIFORM, LOG_NORMAL, LEVEL10 };

char buf[MaxLine];
bool mate = false;
int insertLength = 250;
char line[MaxLine];
int depth = 100;
double errorRate = 0.01;
int readLength = 35;
uint64 readNumber = 0;
int errorNums[1000] = {0};
bool separate = false;
vector<Sequence> refs;
vector<string> names;
vector<double> numbers;
bool is_uniform = true;
int type = 0;
const double Pi = 3.141592653589793;

double Normal()
{
    double x = rand()*1.0 / RAND_MAX;
    double y = rand()*1.0 / RAND_MAX;
    return sqrt(-2*log(x)) * cos(2*Pi*y);
}

void Simulate(Sequence &ref, uint64 total_number, FastaWriter &writer)
{
    static int ref_id = 0;
    ++ref_id;

    if (ref.Size() < readLength)
        return;

    if (mate && ref.Size() < insertLength)
        return;

    if (!mate)
    {
        Sequence seq;
        for (unsigned i = 0; i < total_number; ++i)
        {
            unsigned offset = 0;
            while (true)
            {
                offset = rand() % (ref.Size() - readLength + 1);
                seq.Assign(ref, offset, readLength);
                //ref.GetSubSequence(seq, offset, readLength);
                if (seq.IsValid())
                    break;
            }

            //seq.Encode();

            int error = 0;
            for (int j = 0; j < readLength; ++j)
            {
                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error;
                    int c = seq[j];
                    while (c == seq[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq[j] = c;
                }
            }
            //seq.Decode();

            ++errorNums[error];

            if (rand() < RAND_MAX/2)
                seq.ReverseComplement();

            writer.Write(seq, FormatString("read%d", i));
            //writer.WriteFormat(seq, "read%d", i);
        }
    }
    else
    {
        Sequence seq1;
        Sequence seq2;
        for (unsigned i = 0; i < total_number; i += 2)
        {
            unsigned offset = 0;
            while (true)
            {
                offset = rand() % (ref.Size() + 1 - insertLength);
                seq1.Assign(ref, offset, readLength);
                seq2.Assign(ref, offset + insertLength - readLength, readLength);
                //ref.GetSubSequence(seq1, offset, readLength);
                //ref.GetSubSequence(seq2, offset + insertLength - readLength, readLength);
                if (seq1.IsValid() && seq2.IsValid())
                    break;
            }

            //seq1.Encode();
            //seq2.Encode();
            int error1 = 0;
            int error2 = 0;
            for (int j = 0; j < readLength; ++j)
            {
                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error1;
                    int c = seq1[j];
                    while (c == seq1[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq1[j] = c;
                }

                if (rand()*1.0 < errorRate*RAND_MAX)
                {
                    ++error2;
                    int c = seq2[j];
                    while (c == seq2[j])
                    {
                        c = rand()/93 % 4;
                    }
                    seq2[j] = c;
                }
            }

            ++errorNums[error1];
            ++errorNums[error2];

            //seq1.Decode();
            //seq2.Decode();
            seq2.ReverseComplement();

            if (i%2 == 0)
            {
                writer.Write(seq2, FormatString("read%d_%d/1", ref_id, i));
                writer.Write(seq1, FormatString("read%d_%d/2", ref_id, i));
//                writer.WriteFormat(seq2, "read%d_%d/1", ref_id, i);
//                writer.WriteFormat(seq1, "read%d_%d/2", ref_id, i);
            }
            else
            {
                writer.Write(seq1, FormatString("read%d_%d/1", ref_id, i));
                writer.Write(seq2, FormatString("read%d_%d/2", ref_id, i));
//                writer.WriteFormat(seq1, "read%d_%d/1", ref_id, i);
//                writer.WriteFormat(seq2, "read%d_%d/2", ref_id, i);
            }
        }
    }
}

OptionSet option_set;
int num = 1;

int main(int argc, char *argv[])
{
    option_set.AddOption("depth", "", depth, "");
    option_set.AddOption("errorRate", "", errorRate, "");
    option_set.AddOption("readLength", "", readLength, "");
    option_set.AddOption("mate", "", mate, "");
    option_set.AddOption("insertLength", "", insertLength, "");
    option_set.AddOption("separate", "", separate, "");
    option_set.AddOption("type", "", type, "");
    option_set.AddOption("num", "n", num, "");

    option_set.ProcessOptions(argc, argv);

    //ProcessParameters(argc, argv);

    if (argc < 2 + num
        || strcmp(argv[1], "--help") == 0
        || strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: simRead [--num n] ref-file1 ref-file2 ... read-file\n");
        fprintf(stderr, "       [--depth d] [--mate]\n");
        fprintf(stderr, "       [--errorRate e] [--readLength l]\n");
        exit(1);
    }

    for (int i = 0; i < num; ++i)
    {
        double ratio = 1;
        if (type == EQUAL)
            ratio = 1;
        else if (type == UNIFORM)
            ratio = rand() * 1.0 / RAND_MAX;
        else if (type == LOG_NORMAL)
            ratio = exp(Normal());
        else if (type == LEVEL10)
            ratio = rand() % 1000;

        FastaReader reader(argv[1+i]);
        Sequence ref;
        string comment;
        while (reader.Read(ref, comment))
        {
            refs.push_back(ref);
            names.push_back(comment);
            numbers.push_back(ratio);
        }
    }

    FastaWriter writer(argv[1 + num]);

    int64 total_length = 0;
    for (unsigned i = 0; i < refs.size(); ++i)
        total_length += refs[i].Size();

    if (readNumber == 0)
        readNumber = (depth * total_length) / readLength;
    else
        depth = int(readLength*readNumber/total_length);

    readNumber += readNumber % 2;
    LogMessage("depth %d\treadLength %d\terrorRate %.3f\treadNumber %d\n",
        depth, readLength, errorRate, readNumber);
    LogMessage("insertLength %d\n", insertLength);
    LogMessage("mate %d\n", mate);
    LogMessage("type %d\n", type);

    numbers.resize(refs.size());

    int64 sum = 0;
    for (unsigned i = 0; i < refs.size(); ++i)
        sum += numbers[i] * refs[i].Size();

    for (unsigned i = 0; i < refs.size(); ++i)
        numbers[i] = numbers[i] * refs[i].Size() * readNumber / sum;

    for (unsigned i = 0; i < refs.size(); ++i)
        if (insertLength <= refs[i].Size())
            Simulate(refs[i], (int64)numbers[i], writer);

    sprintf(buf, "%s.list", argv[2]);
    FILE *flist = fopen(buf, "wb");
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        fprintf(flist, "%s %.6f %.6f\n", names[i].c_str(), numbers[i], numbers[i] * readLength / refs[i].Size());
    }

    return 0;
}

