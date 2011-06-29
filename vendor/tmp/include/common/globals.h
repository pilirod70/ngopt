/**
 * @file globals.h
 * @brief Global Header
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __GLOBALS_H_

#define __GLOBALS_H_

#define DEBUG           1


typedef long long int64;
typedef unsigned long long uint64;


const int64 MaxLine = (1 << 20);
const int64 MaxThreads = (1 << 6);
const int64 MaxKValue = 96;

extern int kmerLength;


#endif

