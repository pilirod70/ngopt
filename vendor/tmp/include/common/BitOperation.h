/**
 * @file BitOperation.h
 * @brief BitOperation Class Library
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __BIT_OPERATION_H_

#define __BIT_OPERATION_H_

#include "globals.h"

class BitOperation
{
public:
    static void ReverseComplement(uint64 &value)
    {
        value = ((value & SwapMask32) << 32) | ((value & ~SwapMask32) >> 32);
        value = ((value & SwapMask16) << 16) | ((value & ~SwapMask16) >> 16);
        value = ((value & SwapMask8) << 8) | ((value & ~SwapMask8) >> 8);
        value = ((value & SwapMask4) << 4) | ((value & ~SwapMask4) >> 4);
        value = ((value & SwapMask2) << 2) | ((value & ~SwapMask2) >> 2);
        value = ~value;
    }

    static int BitCount(unsigned x)
    {
        x = (x & SwapMask1) + ((x >> 1) & SwapMask1);
        x = (x & SwapMask2) + ((x >> 2) & SwapMask2);
        x = (x & SwapMask4) + ((x >> 4) & SwapMask4);
        x = (x & SwapMask8) + ((x >> 8) & SwapMask8);
        x = (x & SwapMask16) + ((x >> 16) & SwapMask16);
        return x;
    }

    static unsigned bitToIndex[9];

private:
    static const uint64 SwapMask32 = 0x00000000FFFFFFFFULL;
    static const uint64 SwapMask16 = 0x0000FFFF0000FFFFULL;
    static const uint64 SwapMask8 = 0x00FF00FF00FF00FFULL;
    static const uint64 SwapMask4 = 0x0F0F0F0F0F0F0F0FULL;
    static const uint64 SwapMask2 = 0x3333333333333333ULL;
    static const uint64 SwapMask1 = 0x5555555555555555ULL;
};

#endif

