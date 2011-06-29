/**
 * @file SequenceReader.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-01
 */

#include "globals.h"

#include "Sequence.h"
#include "SequenceReader.h"
#include "Read.h"
#include "Contig.h"

#include <string>
#include <deque>

using namespace std;

bool SequenceReader::Read(Sequence &seq)
{
    string comment;
    return Read(seq, comment);
}

bool SequenceReader::Read(Sequence &seq, string &comment)
{
    string quality;
    return Read(seq, comment, quality);
}

bool SequenceReader::Read(Sequence &seq, string &comment, string &quality)
{
    return ReadRecord(seq, comment, quality);
}

bool SequenceReader::Read(Contig &contig)
{
    string comment;
    return Read(contig, comment);
}

bool SequenceReader::Read(Contig &contig, string &comment)
{
    string quality;
    return Read(contig, comment, quality);
}

bool SequenceReader::Read(Contig &contig, string &comment, string &quality)
{
    contig.sum_coverage = 0;
    contig.is_tangle = false;

    if (!ReadRecord(contig, comment, quality))
        return false;

    size_t index = comment.rfind('_');
    if (index != string::npos)
        contig.sum_coverage = atoi(comment.c_str() + index + 1);
    comment = comment.substr(0, index);

    return true;
}

int64 SequenceReader::Read(deque< ::Read> &reads)
{
    reads.resize(0);

    Sequence seq;
    while (Read(seq))
    {
        seq.TrimError();
        //seq.Encode();
        reads.resize(reads.size()+1);
        reads.back() = seq;
    }

    return reads.size();
}

int64 SequenceReader::Read(deque<Sequence> &sequences)
{
    sequences.resize(0);

    Sequence seq;
    while (Read(seq))
    {
        //seq.Encode();
        sequences.push_back(seq);
    }

    return sequences.size();
}

int64 SequenceReader::Read(deque<Contig> &contigs)
{
    contigs.resize(0);

    Contig contig;
    while (Read(contig))
    {
        //contig.Encode();
        contigs.push_back(contig);
    }

    return contigs.size();
}

bool FastaReader::ReadRecord(Sequence &seq, string &comment, string &quality)
{
    return ReadFasta(*is, seq, comment);
}

bool FastqReader::ReadRecord(Sequence &seq, string &comment, string &quality)
{
    return ReadFastq(*is, seq, comment, quality);
}

