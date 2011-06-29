/**
 * @file analyzeBed.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Utils.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct Exon
{
    string chr;
    int start;
    int size;

    bool operator < (const Exon &e) const
    {
        if (chr != e.chr)
            return chr < e.chr;
        else if (start != e.start)
            return start < e.start;
        else
            return size < e.size;
    }
};

struct BedRecord
{
    string chr;
    string type;
    vector<Exon> exons;

    void Parse(char *line)
    {
        stringstream ss(line);
        ss >> chr;

        int start, end;
        ss >> start >> end;

        string tmp;
        ss >> tmp;
        type = tmp.substr(0, 2);

        ss >> tmp >> tmp >> tmp >> tmp >> tmp;

        int size;
        ss >> size;

        string size_buf;
        string start_buf;
        ss >> size_buf >> start_buf;

        ReplaceComma(size_buf);
        ReplaceComma(start_buf);

        stringstream size_ss(size_buf);
        stringstream start_ss(start_buf);

        for (int i = 0; i < size; ++i)
        {
            Exon exon;
            exon.chr = chr;

            size_ss >> exon.size;
            start_ss >> exon.start;

            exon.start += start;
            exons.push_back(exon);
        }
    }
};

char buf[MaxLine];
char line[MaxLine];
map<string, int> type_dict;
map<Exon, int> freq_dict;
map<int, int> length_dict;

int freq_dist[1000] = {0};
int length_dist[1000] = {0};
int num_dist[1000] = {0};

int main(int argc, char *argv[])
{
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        sscanf(line, "%s", buf);
        if (strcmp(buf, "track") == 0)
            break;
    }

    int total_exons = 0;
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        BedRecord record;
        record.Parse(line);

        type_dict[record.type]++;
        for (unsigned i = 0; i < record.exons.size(); ++i)
        {
            freq_dict[record.exons[i]]++;
            length_dict[record.exons[i].size]++;
        }

        total_exons += record.exons.size();
        ++num_dist[record.exons.size()];
    }

    cout << "# of different exons: " << freq_dict.size() << endl;
    cout << "# of exons: " << total_exons << endl;

    const char *types[] = { "CE", "ME", "EI", "II", "IR", };
    for (int i = 0; i < 5; ++i)
        cout << types[i] << ", " << type_dict[types[i]] << endl;

    for (map<Exon, int>::iterator iter = freq_dict.begin(); iter != freq_dict.end(); ++iter)
    {
        ++freq_dist[iter->second];
    }

    for (map<int, int>::iterator iter = length_dict.begin(); iter != length_dict.end(); ++iter)
    {
        int k = 0;
        while (iter->first >= (k+1) * 50)
            ++k;
        length_dist[k] += iter->second;
    }

    for (int i = 0; i < 100; ++i)
        cout << i << ", " << freq_dist[i] << endl;

    for (int i = 0; i < 100; ++i)
        cout << i*50 << "-" << (i+1)*50-1 << ", " << length_dist[i] << endl;

    for (int i = 0; i < 100; ++i)
        cout << i << ", " << num_dist[i] << endl;

    return 0;
}


