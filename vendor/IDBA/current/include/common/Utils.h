/**
 * @file Utils.h
 * @brief Utilities
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __UTILS_H_

#define __UTILS_H_

#include "globals.h"

#include <cstdio>
#include <vector>
#include <string>

std::string FormatString(const char *fmt, ...);

bool IsExist(const std::string &filename);
FILE *OpenFile(const std::string &filename, const std::string &mode);

inline int CodeToChar(int ch)
{
    switch (ch)
    {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        case 4: return 'N';
        default: return 'N';
    }
}

inline int CharToCode(int ch)
{
    switch (ch)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        default: return 4;
    }
} 

int64 Sum(const std::vector<int64> &lengths);
int64 Average(const std::vector<int64> &lengths);
int64 Maximum(const std::vector<int64> &lengths);
int64 Nxx(const std::vector<int64> &lengths, double percent);

void Replace(std::string &s, int ch1, int ch2);
void ReplaceComma(std::string &s);
void SplitString(const std::string &s, std::vector<std::string> &items);

#endif

