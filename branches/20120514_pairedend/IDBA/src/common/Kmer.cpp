/**
 * @file Kmer.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"
#include "Kmer.h"
#include "Sequence.h"

#include <iostream>

using namespace std;

ostream &operator <<(ostream &stream, const Kmer &kmer)
{
    Sequence seq;
    seq = kmer;
    //seq.Decode();
    return stream << seq;
}

