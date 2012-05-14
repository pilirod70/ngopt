/**
 * @file SequenceWriter.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-01
 */


#include "globals.h"

#include "Read.h"
#include "Sequence.h"
#include "SequenceWriter.h"
#include "Contig.h"
#include "Utils.h"

#include <string>
#include <deque>

using namespace std;

bool SequenceWriter::Write(const Sequence &seq)
{
    string comment;
    return Write(seq, comment);
}

bool SequenceWriter::Write(const Sequence &seq, const string &comment)
{
    string quality;
    return Write(seq, comment, quality);
}

bool SequenceWriter::Write(const Sequence &seq, const string &comment, const string &quality)
{
    return WriteRecord(seq, comment, quality);
}

bool SequenceWriter::Write(const Contig &contig)
{
    string comment;
    return Write(contig, comment);
}

bool SequenceWriter::Write(const Contig &contig, const string &comment)
{
    string quality;
    return Write(contig, comment, quality);
}

bool SequenceWriter::Write(const Contig &contig, const string &comment, const string &quality)
{
    return WriteRecord(contig, FormatString("%s_%d", comment.c_str(), contig.sum_coverage), quality);
}

int64 SequenceWriter::Write(const deque<Sequence> &sequences)
{
    for (unsigned i = 0; i < sequences.size(); ++i)
        Write(sequences[i], FormatString("seq_%d", i));
    return sequences.size();
}

int64 SequenceWriter::Write(const deque<Contig> &contigs)
{
    for (unsigned i = 0; i < contigs.size(); ++i)
        Write(contigs[i], FormatString("contig_%d", i));
    return contigs.size();
}

bool FastaWriter::WriteRecord(const Sequence &seq, const string &comment, const string &quality)
{
    return WriteFasta(*os, seq, comment);
}

bool FastqWriter::WriteRecord(const Sequence &seq, const string &comment, const string &quality)
{
    return WriteFastq(*os, seq, comment, quality);
}

