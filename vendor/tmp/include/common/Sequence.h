/**
 * @file Sequence.h
 * @brief Sequence Class
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __SEQUENCE_H_

#define __SEQUENCE_H_

#include "globals.h"
#include "Kmer.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <istream>
#include <ostream>

class Read;

class Sequence
{
public:
    friend std::istream &operator >>(std::istream &stream, Sequence &seq);
    friend std::ostream &operator <<(std::ostream &stream, const Sequence &seq);

    Sequence() {}
    Sequence(const Sequence &seq): bases(seq.bases) {}
    Sequence(const std::string &seq): bases(seq) {}
    Sequence(const char *seq): bases(seq) {}
    Sequence(const Read &read) { Assign(read); }
    Sequence(const Kmer &kmer) { Assign(kmer, kmerLength); }

    void Clear() { bases.clear(); }

    void Swap(Sequence &seq) { std::swap(bases, seq.bases); }

    char &operator [](int index) { return bases[index]; }
    const char &operator [](int index) const { return bases[index]; }

    const Sequence &operator =(const Sequence &seq) { Assign(seq); return *this; }
    const Sequence &operator =(const std::string &seq) { Assign(seq); return *this; }
    const Sequence &operator =(const char *seq) { Assign(seq); return *this; }
    const Sequence &operator =(const Read &read) { Assign(read); return *this; }
    const Sequence &operator =(const Kmer &kmer) { Assign(kmer, kmerLength); return *this; }

    const Sequence &operator +=(const Sequence &seq) { AddSequence(seq); return *this; }
    const Sequence &operator +=(const std::string &seq) { AddSequence(seq); return *this; }
    const Sequence &operator +=(const char *seq) { AddSequence(seq); return *this; }
    const Sequence &operator +=(unsigned ch) { AddNucleotide(ch); return *this; }

    bool operator ==(const Sequence &seq) const { return bases == seq.bases; }
    bool operator !=(const Sequence &seq) const { return bases != seq.bases; }
    bool operator <(const Sequence &seq) const { return bases < seq.bases; }

    void Assign(const Sequence &seq, int offset = 0, size_t length = std::string::npos)
    { bases.assign(seq.bases, offset, length); }
    void Assign(const std::string &s, int offset = 0, size_t length = std::string::npos)
    { bases.assign(s, offset, length); }
    void Assign(const char *s, int offset = 0, size_t length = std::string::npos)
    { bases.assign(s, offset, length); }
    void Assign(const Read &read);
    void Assign(const Kmer &kmer, int kmerLength);

    std::string ToString() const { Sequence tmp(bases); tmp.Decode(); return tmp.bases; }
    const char *ToCString() const { return bases.c_str(); }

    bool IsValid() const; 

    int Size() const { return bases.size(); }
    void Resize(int l) { bases.resize(l); }

    void ReverseComplement();

    void Trim(int t) { bases.resize(std::max(0LL, (int64)bases.size() - t)); }
    void TrimError();

    //void GetSubSequence(Sequence &seq, int offset, int sub_length) const;
    Kmer GetKmer(int offset, int kmer_length) const;

    void AddNucleotide(unsigned value) { bases.append(1, value); }

protected:
    void AddSequence(const Sequence &seq) { bases.append(seq.bases); }
    void AddSequence(const std::string &seq) { bases.append(seq); }
    void AddSequence(const char *seq) { bases.append(seq); }
    void AddSequence(const char *seq, int length) { bases.append(seq, length); }

    void Encode();
    void Decode();

    bool IsValid(char ch) const
    { return ch == 'A' || ch ==  'C' || ch == 'G' || ch == 'T' || ch == 0 || ch == 1 || ch == 2 || ch == 3 || ch == 4; }

private:
    std::string bases;
};

std::istream &ReadFasta(std::istream &is, Sequence &seq, std::string &comment);
std::ostream &WriteFasta(std::ostream &os, const Sequence &seq, const std::string &comment);

std::istream &ReadFastq(std::istream &is, Sequence &seq, std::string &comment, std::string &quality);
std::ostream &WriteFastq(std::ostream &os, const Sequence &seq, const std::string &comment, const std::string &quality);

#endif

