/**
 * @file Utils.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Log.h"
#include "Utils.h"
#include "Sequence.h"

#include <cstdio>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdarg>

#include <getopt.h>

using namespace std;

static char line[MaxLine];

string FormatString(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsprintf(line, fmt, ap);
    va_end(ap);

    return line;
}

bool IsExist(const string &filename)
{
    FILE *fp = fopen(filename.c_str(), "rb");
    if (fp == NULL)
    {
        return false;
    }
    else
    {
        fclose(fp);
        return true;
    }
}

FILE *OpenFile(const string &filename, const string &mode)
{
    FILE *fp = fopen(filename.c_str(), mode.c_str());

    if (fp == NULL)
    {
        LogError("open %s failed\n", filename.c_str());
        exit(1);
    }

    return fp;
}

int64 Sum(const vector<int64> &lengths)
{
    int64 sum = 0;
    for (unsigned i = 0; i < lengths.size(); ++i)
        sum += lengths[i];
    return sum;
}

int64 Average(const vector<int64> &lengths)
{
    if (lengths.size() == 0)
        return 0;
    else
        return Sum(lengths) / lengths.size();
}

int64 Maximum(const vector<int64> &lengths)
{
    return *max_element(lengths.begin(), lengths.end());
}

int64 Nxx(const vector<int64> &lengths, double percent)
{
    vector<int64> v(lengths);
    sort(v.begin(), v.end(), greater<int64>());

    int64 sum = Sum(v);
    int64 nxx = 0;
    int64 partial = 0;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        partial += v[i];
        if (nxx == 0 && partial > percent * sum)
        {
            nxx = v[i];
            break;
        }
    }
    return nxx;
}

void Replace(string &s, int ch1, int ch2)
{
    for (unsigned i = 0; i < s.size(); ++i)
    {
        if (s[i] == ch1)
            s[i] = ch2;
    }
}

void ReplaceComma(string &s)
{
    for (unsigned i = 0; i < s.size(); ++i)
    {
        if (s[i] == ',')
            s[i] = ' ';
    }
}

void SplitString(const string &s, vector<string> &items)
{
    items.resize(0);

    stringstream ss(s);
    string item;

    while (ss >> item)
        items.push_back(item);
}


