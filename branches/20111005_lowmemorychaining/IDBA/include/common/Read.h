/**
 * @file Read.h
 * @brief Read Class whic is a compact structure for storing reads
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __READ_H_

#define __READ_H_

#include "globals.h"

class Sequence;

class Read
{
public:
    Read() { status = 0; length = 0; }
    Read(const Sequence &seq)
    { SetContent(seq); }

    const Read &operator =(const Sequence &seq)
    { SetContent(seq); ;return *this; }

    char operator[] (int index) const
    { return GetNucleotide(index); }

    int Size() const { return length; }
    void Resize(int new_length) { length = new_length; }

    bool IsActive() const { return status == 0; }
    void Activate() { status = 0; }
    void Inactivate() { status = 1; }

    static const int MaxReadLength = 100;
    static const int BufferSize = MaxReadLength/4 + 1;

private:
    unsigned GetNucleotide(unsigned index) const
    { return (compressed[index>>2] >> ((index&3) << 1)) & 3; }
    void SetNucleotide(unsigned index, unsigned char value)
    {
        compressed[index>>2] &= ~(3 << ((index&3) << 1));
        compressed[index>>2] |= value << ((index&3) << 1);
    }

    void SetContent(const Sequence &seq);

    unsigned char compressed[BufferSize];
    unsigned char status;
    short length;
};

#endif

