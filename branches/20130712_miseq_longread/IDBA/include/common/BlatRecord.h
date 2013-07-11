/**
 * @file BlatRecord.h
 * @brief BlatRecord Class 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __BLAT_RECORD_H_

#define __BLAT_RECORD_H_

#include "globals.h"

#include <vector>
#include <string>

struct BlatBlock
{
    int64 query_from;
    int64 ref_from;
    int64 size;
};

struct BlatRecord
{
    void Parse(const std::string &record);

    bool operator <(const BlatRecord &r) const
    { return match_count > r.match_count; }

    std::string query_name;
    std::string ref_name;
    int64 match_count;
    int64 query_from;
    int64 query_to;
    int64 query_length;
    int64 ref_from;
    int64 ref_to;
    int64 ref_length;
    bool is_reverse;

    std::vector<BlatBlock> blocks;
};

#endif

