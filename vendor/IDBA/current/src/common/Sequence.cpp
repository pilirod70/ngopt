/**
 * @file Sequence.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "Log.h"
#include "Kmer.h"
#include "Sequence.h"
#include "Read.h"
#include "Utils.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

using namespace std;

istream &operator >>(istream &is, Sequence &seq)
{
    getline(is, seq.bases);
    if (!is)
        return is;

    string line;
    while (is && (isalnum(is.peek()) || is.peek() == '\n') && getline(is, line))
        seq.bases += line;
    is.clear();

    seq.Encode();

    return is;
}

ostream &operator <<(ostream &os, const Sequence &seq)
{
    Sequence tmp = seq;
    tmp.Decode();
    return os << tmp.bases;
}

bool Sequence::IsValid() const
{
    for (unsigned i = 0; i < bases.size(); ++i)
    {
        if (!IsValid(bases[i]))
            return false;
    }
    return true;
}

void Sequence::ReverseComplement()
{
    reverse(bases.begin(), bases.end());
    for (unsigned i = 0; i < bases.size(); ++i)
    {
        switch(bases[i])
        {
            case 'A':
                bases[i] = 'T';
                break;

            case 'C':
                bases[i] = 'G';
                break;

            case 'G':
                bases[i] = 'C';
                break;

            case 'T':
                bases[i] = 'A';
                break;

            case 0:
            case 1:
            case 2:
            case 3:
                bases[i] = 3 - bases[i];
                break;

            case 'N':
            case 4:
                break;

#ifdef DEBUG
            default:
                throw logic_error(FormatString("reverse complement error: unkown character %c", bases[i]));
#endif
        }
    }
}

void Sequence::TrimError()
{ 
    int length = bases.size();
    while (length > 0 && !IsValid(bases[length-1]))
        --length;
    bases.resize(length);
}

Kmer Sequence::GetKmer(int offset, int kmerLength) const
{
    Kmer kmer;
    for (int i = 0; i < kmerLength; ++i)
        kmer.AddRight(bases[offset + i]);
    return kmer;
}

void Sequence::Assign(const Read &read)
{
    bases.resize(read.Size());
    for (int i = 0; i < read.Size(); ++i)
        bases[i] = read[i];
}

void Sequence::Assign(const Kmer &kmer, int kmerLength)
{
    bases.resize(kmerLength);
    for (int i = 0; i < kmerLength; ++i)
        bases[i] = kmer.GetBase(i);
}

void Sequence::Encode()
{
    for (unsigned i = 0; i < bases.size(); ++i)
    {
        switch(bases[i])
        {
            case 'A':
                bases[i] = 0;
                break;

            case 'C':
                bases[i] = 1;
                break;

            case 'G':
                bases[i] = 2;
                break;

            case 'T':
                bases[i] = 3;
                break;

            case 'N':
                bases[i] = 4;
                break;

            default:
                bases[i] = 4;
                break;

//#ifdef DEBUG
//            default:
//                throw logic_error(FormatString("encode error: unkown character %c", bases[i]));
//#endif
        }
    }
}

void Sequence::Decode()
{
    for (unsigned i = 0; i < bases.size(); ++i)
    {
        switch(bases[i])
        {
            case 0:
                bases[i] = 'A';
                break;

            case 1:
                bases[i] = 'C';
                break;

            case 2:
                bases[i] = 'G';
                break;

            case 3:
                bases[i] = 'T';
                break;

            case 4:
                bases[i] = 'N';
                break;

#ifdef DEBUG
            default:
                throw logic_error(FormatString("decode error: unkown character %c", bases[i]));
#endif
        }
    }
}

istream &ReadFasta(istream &is, Sequence &seq, string &comment)
{
    string line;
    getline(is, line);

    if (!is)
        return is;

    comment = line.substr(1);

    return is >> seq;
}

ostream &WriteFasta(ostream &os, const Sequence &seq, const string &comment)
{
    return os << ">" << comment << "\n" << seq << "\n";
}

istream &ReadFastq(istream &is, Sequence &seq, string &comment, string &quality)
{
    string line;
    getline(is, line);
    if (!is)
        return is;

    comment = line.substr(1);
    is >> seq;
    getline(is, line);

    quality = "";
    while (is && is.peek() != '@' && getline(is, line))
        quality += line;
    is.clear();

    return is;
}

ostream &WriteFastq(ostream &os, const Sequence &seq, const string &comment, const string &quality)
{
    return os << "@" << comment << "\n" << seq << "\n" << "+" << "\n" << quality << "\n";
}

