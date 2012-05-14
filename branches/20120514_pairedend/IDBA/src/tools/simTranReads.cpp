/**
 * @file simTranReads.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "BlatRecord.h"
#include "Log.h"
#include "Read.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"
#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

using namespace std;

char line[MaxLine];
char buf[MaxLine];
vector<Sequence> sequences;
vector<string> names;
map<string, int> dict;
int dist[100] = {0};

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        fprintf(stderr, "simTranReads transcript.fa transcript.cover reads.fa\n");
        exit(1);
    }

    FastaReader reader(argv[1]);
    FastaWriter writer(argv[3]);

    //sequences.reserve(reader.NumReads());
    Sequence seq;
    string comment;
    while (reader.Read(seq, comment))
    {
        sequences.push_back(seq);
        names.push_back(comment);

    }

    FILE *fcover = fopen(argv[2], "rb");
    int index = 0;
    while (fgets(line, MaxLine, fcover) != NULL)
    {
        double depth = 0;
        sscanf(line, "%s %lf", buf, &depth);

        if (depth > 0)
        {
            //writer.WriteFormat(sequences[index], "%s", names[index].c_str());
            writer.Write(sequences[index], names[index]);

            comment = buf;
            Replace(comment, '_', ' ');
            sscanf(comment.c_str(), "%s", buf);
            dict[buf]++;
        }

        ++index;
    }

    cout << dict.size() << endl;
    for (map<string, int>::iterator iter = dict.begin(); iter != dict.end(); ++iter)
    {
        ++dist[iter->second];
    }
    
    cout << "tran distribution" << endl;
    for (int i = 0; i < 100; ++i)
    {
        cout << i << " " << dist[i] << endl;
    }

    return 0;
}

