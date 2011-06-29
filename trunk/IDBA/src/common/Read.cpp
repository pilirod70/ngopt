/**
 * @file Read.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Sequence.h"
#include "Read.h"

#include <cstdio>
#include <algorithm>
#include <iostream>

using namespace std;

void Read::SetContent(const Sequence &seq)
{
    for (int i = 0; i < seq.Size(); ++i)
        SetNucleotide(i, seq[i]);
    length = seq.Size();
    status = 0;
}

