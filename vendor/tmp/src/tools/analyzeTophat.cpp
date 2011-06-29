/**
 * @file analyzeTophat.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Utils.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

struct Junction
{
    string chr;
    int prev_end;
    int next_start;

    bool operator <(const Junction &junc) const
    {
        if (chr != junc.chr)
            return chr < junc.chr;
        else if (prev_end != junc.prev_end)
            return prev_end < junc.prev_end;
        else
            return next_start < junc.next_start;
    }
};

struct Record
{
    string name;
    string chr;
    vector<Junction> juncs;
    vector<int> exon_from;
    vector<int> exon_to;

    void Parse(char *line)
    {
        stringstream ss(line);
        ss >> name >> chr;

        string tmp;
        ss >> tmp >> tmp >> tmp >> tmp >> tmp;

        int size;
        ss >> size;

        string from_string;
        string to_string;
        ss >> from_string >> to_string;

        ReplaceComma(from_string);
        ReplaceComma(to_string);

        stringstream from_stream(from_string);
        stringstream to_stream(to_string);

        exon_from.resize(size);
        exon_to.resize(size);

        for (int i = 0; i < size; ++i)
        {
            from_stream >> exon_from[i];
            to_stream >> exon_to[i];
        }

        for (int i = 0; i+1 < size; ++i)
        {
            Junction junc;
            junc.chr = chr;
            junc.prev_end = exon_to[i];
            junc.next_start = exon_from[i+1];

            juncs.push_back(junc);
        }
    }
};


map<string, vector<int> > matchs;
map<Junction, int> junctions;

char buf[MaxLine];
char line[MaxLine];

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        fprintf(stderr, "analyszeTophat tophat-prefix ucsc-pos spectra.csv\n");
        exit(1);
    }

    sprintf(buf, "%s/junctions.bed", argv[1]);
    FILE *fjunc = OpenFile(buf, "rb");

    string tmp;
    while (fgets(line, MaxLine, fjunc) != NULL)
    {
        stringstream ss(line);
        Junction junc;
        ss >> junc.chr;

        int from;
        int to;
        ss >> from >> to;

        int count = 0;
        ss >> tmp >> count >> tmp >> tmp;

        int size = 0;
        ss >> tmp >> tmp >> size;
        
        string length_string;
        string offset_string;
        ss >> length_string >> offset_string;
        ReplaceComma(length_string);
        ReplaceComma(offset_string);

        stringstream length_stream(length_string);
        stringstream offset_stream(offset_string);

        vector<int> exon_from;
        vector<int> exon_to;
        for (int i = 0; i < size; ++i)
        {
            int length;
            int offset;
            length_stream >> length;
            offset_stream >> offset;

            exon_from.push_back(from + offset);
            exon_to.push_back(from + offset + length);
        }

        for (unsigned i = 0; i+1 < exon_from.size(); ++i)
        {
            junc.prev_end = exon_to[i];
            junc.next_start = exon_from[i+1];

            junctions[junc] = count;
        }
    }

//    vector<string> chr;
//    vector<int> from;
//    vector<int> to;
//    vector<int> count;
    sprintf(buf, "%s/coverage.wig", argv[1]);
    FILE *fcover = OpenFile(buf, "rb");

    int f, t, c;
    fgets(line, MaxLine, fcover);
    while (fscanf(fcover, "%s %d %d %d", buf, &f, &t, &c) != EOF)
    {
        vector<int> &match = matchs[buf];
        match.resize(t);
        for (int i = f; i < t; ++i)
            match[i] = c;
    }

    int found_isoforms = 0;
    int found_isoforms_cover = 0;
    int num_records = 0;

    FILE *fpos = OpenFile(argv[2], "rb");
    FILE *fspec = OpenFile(argv[3], "wb");

    fgets(line, MaxLine, fpos);
    while (fgets(line, MaxLine, fpos) != NULL)
    {
        ++num_records;

        Record record;
        record.Parse(line);

        int count = 0;
        for (unsigned i = 0; i < record.juncs.size(); ++i)
        {
            Junction &junc = record.juncs[i];

            if (junctions.find(junc) != junctions.end())
                ++count;
        }

        if (count == (int)record.juncs.size())
        {
            ++found_isoforms;

            int count = 0;
            int total = 0;
            int consecutive_count = 0;
            int maximum = 0;
            vector<int> &match = matchs[record.chr];
            for (unsigned i = 0; i < record.exon_from.size(); ++i)
            {
                for (int j = record.exon_from[i]; j < record.exon_to[i] && j < (int)match.size(); ++j)
                {
                    fprintf(fspec, "%d, ", match[j]);
                    if (match[j] >= 8)
                    {
                        consecutive_count++;
                        if (consecutive_count > maximum)
                            maximum = consecutive_count;
                        ++count;
                    }
                    else
                        consecutive_count = 0;

                    ++total;
                }
            }

            fprintf(fspec, "\n");

            if (count >= 0.9 * total && maximum >= 0.9 * total)
            {
                ++found_isoforms_cover;
            }

            fprintf(fspec, "%d, %d, %d, \n", count, maximum, total);
        }
    }

    cout << "found: " << found_isoforms << " " << found_isoforms_cover << " in " << num_records << endl;

    return 0;
}

