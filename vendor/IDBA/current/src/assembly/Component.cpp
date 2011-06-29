/**
 * @file Component.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.20
 * @date 2011-01-20
 */

#include "globals.h"

#include "Component.h"

#include <vector>

using namespace std;

void Component::PrintAlignment(const string &seq, int offset)
{
    if (offset == -1)
    {
        offset = 0;
        for (unsigned i = 0; i < alignment.size(); ++i)
        {
            if (offset < alignment[i].size())
                offset = alignment[i].size();
        }
    }

    vector<string> &buffer_strings = alignment;

    int d = offset;
    int k = -1;
    for (unsigned j = 0; j < buffer_strings.size(); ++j)
    {
        if ((int)buffer_strings[j].size() <= d)
        {
            k = j;
            break;
        }
    }
    if (k == -1)
    {
        buffer_strings.push_back("");
        k = buffer_strings.size() - 1;
    }

    if ((int)buffer_strings[k].size() < d)
        buffer_strings[k].append(d - buffer_strings[k].size(), ' ');

    if (d == 0)
        buffer_strings[k] += seq;
    else
        buffer_strings[k] += seq.substr(kmerLength - 1);
}

void Component::CompleteAlignment()
{
    vector<string> &buffer_strings = alignment;

    int length = 0;
    for (unsigned i = 0; i < buffer_strings.size(); ++i)
    {
        if ((int)buffer_strings[i].size() > length)
            length = buffer_strings[i].size();
    }

    for (unsigned i = 0; i < buffer_strings.size(); ++i)
    {
        if ((int)buffer_strings[i].size() < length)
            buffer_strings[i].append(length - buffer_strings[i].size(), ' ');
    }
}

