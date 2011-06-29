/**
 * @file filterBlat.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "BlatRecord.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

using namespace std;

char buf[MaxLine];
char line[MaxLine];

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "filterBlat blat-file\n");
        exit(1);
    }

    FILE *fblat = fopen(argv[1], "rb");

    fgets(line, MaxLine, fblat);
    printf("%s", line);
    fgets(line, MaxLine, fblat);
    printf("%s", line);
    fgets(line, MaxLine, fblat);
    printf("%s", line);
    fgets(line, MaxLine, fblat);
    printf("%s", line);
    fgets(line, MaxLine, fblat);
    printf("%s", line);

    int index = 0;
    string last_query;
    while (fgets(line, MaxLine, fblat) != NULL)
    {
        BlatRecord record;
        record.Parse(line);
        
        ++index;
        if (index % 100000 == 0)
            cerr << index << endl;
        
        if (record.match_count < 0.95 * record.query_length)
            continue;

//        if (last_query == record.query_name)
//            continue;
//        last_query = record.query_name;

        printf("%s", line);
    }

    return 0;
}

