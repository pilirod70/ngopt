/**
 * @file SequenceWriter.h
 * @brief 
 * @author Yu Peng
 * @version 0.19
 * @date 2010-12-01
 */

#ifndef __SEQUENCE_WRITER_H_

#define __SEQUENCE_WRITER_H_

#include "globals.h"

#include "Sequence.h"

#include <deque>
#include <istream>
#include <fstream>
#include <string>

class Contig;

class SequenceWriter
{
public:
    SequenceWriter(std::ostream &os, bool is_own = false)
    { this->os = &os; this->is_own = is_own; }
    virtual ~SequenceWriter() { Dispose(); }

    virtual void Dispose() { if (is_own) delete os; os = NULL; }

    bool Write(const Sequence &seq);
    bool Write(const Sequence &seq, const std::string &comment);
    bool Write(const Sequence &seq, const std::string &comment, const std::string &quality);

    bool Write(const Contig &contig);
    bool Write(const Contig &contig, const std::string &comment);
    bool Write(const Contig &contig, const std::string &comment, const std::string &quality);
    
    int64 Write(const std::deque<Sequence> &sequences);
    int64 Write(const std::deque<Contig> &contigs);

protected:
    SequenceWriter(const SequenceWriter &);
    const SequenceWriter &operator =(const SequenceWriter &);

    virtual bool WriteRecord(const Sequence &seq, const std::string &comment, const std::string &quality) = 0;

    std::ostream *os; 
    bool is_own;
};

class FastaWriter: public SequenceWriter
{
public:
    explicit FastaWriter(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary)
        : SequenceWriter(*(new std::ofstream(filename.c_str(), mode)), true) {}

protected:
    virtual bool WriteRecord(const Sequence &seq, const std::string &comment, const std::string &quality);
};

class FastqWriter: public SequenceWriter
{
public:
    explicit FastqWriter(const std::string &filename, 
            std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary)
        : SequenceWriter(*(new std::ofstream(filename.c_str(), mode)), true) {}

protected:
    virtual bool WriteRecord(const Sequence &seq, const std::string &comment, const std::string &quality);
};

#endif
