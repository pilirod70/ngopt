/**
 * @file SequenceReader.h
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-01
 */

#ifndef __SEQUENCE_READER_H_

#define __SEQUENCE_READER_H_

#include "globals.h"

#include "Sequence.h"

#include <deque>
#include <istream>
#include <fstream>
#include <string>

class Contig;

class SequenceReader
{
public:
    SequenceReader(std::istream &is, bool is_own = false)
    { this->is = &is; this->is_own = is_own; }
    virtual ~SequenceReader() { Dispose(); }

    virtual void Dispose() { if (is_own) delete is; is = NULL; }

    bool Read(Sequence &seq);
    bool Read(Sequence &seq, std::string &comment);
    bool Read(Sequence &seq, std::string &comment, std::string &quality);

    bool Read(Contig &contig);
    bool Read(Contig &contig, std::string &comment);
    bool Read(Contig &contig, std::string &comment, std::string &quality);
    
    int64 Read(std::deque< ::Read> &reads);
    int64 Read(std::deque<Sequence> &sequences);
    int64 Read(std::deque<Contig> &contigs);

protected:
    SequenceReader(const SequenceReader &);
    const SequenceReader &operator =(const SequenceReader &);

    virtual bool ReadRecord(Sequence &seq, std::string &comment, std::string &quality) = 0;

    std::istream *is; 
    bool is_own;
};

class FastaReader: public SequenceReader
{
public:
    explicit FastaReader(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
        : SequenceReader(*(new std::ifstream(filename.c_str(), mode)), true) {}

protected:
    virtual bool ReadRecord(Sequence &seq, std::string &comment, std::string &quality);
};

class FastqReader: public SequenceReader
{
public:
    explicit FastqReader(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
        : SequenceReader(*(new std::ifstream(filename.c_str(), mode)), true) {}

protected:
    virtual bool ReadRecord(Sequence &seq, std::string &comment, std::string &quality);
};

#endif

