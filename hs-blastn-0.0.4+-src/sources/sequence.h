#ifndef SEQUENCE_H_
#define SEQUENCE_H_ 1

#include "line_reader.h"
#include "utility.h"

struct SequenceTypes
{
    typedef StreamLineReader::SIZE_TYPE   SIZE_TYPE;
    typedef StreamLineReader::STRING_TYPE STRING_TYPE;
    typedef cy_utility::SimpleArray<char> SEQUENCE_TYPE;
};

class FastaReader
{
public:
    typedef SequenceTypes::SIZE_TYPE     SIZE_TYPE;
    typedef SequenceTypes::STRING_TYPE   STRING_TYPE;
    typedef SequenceTypes::SEQUENCE_TYPE SEQUENCE_TYPE;
	
	static const char*                   kClassName;
    static SIZE_TYPE                     line_number;

public:
    static bool IsNucl(unsigned char ch);
    static bool IsAmbigNucl(unsigned char ch);
    static void ParseDefLine(STRING_TYPE& line, SEQUENCE_TYPE& header);
    static void ParseDataLine(STRING_TYPE& line, SEQUENCE_TYPE& sequence);
    static void CheckDataLine(STRING_TYPE& line);
	static char GetFirstLineChar(STRING_TYPE& line);
	static bool IsEmptyLine(STRING_TYPE& line);
    static void ReadOneSeq(StreamLineReader& line_reader);
    static void IncLineNumber();
    static void SetFirstDataLine();
    static SIZE_TYPE CheckFirstDataLine();
    static SIZE_TYPE GetLineNumber();
	static void SetLineNumber(SIZE_TYPE new_line_number) { line_number = new_line_number; }
};

class Sequence
{
public:
    typedef SequenceTypes::SIZE_TYPE      SIZE_TYPE;
    typedef SequenceTypes::STRING_TYPE    STRING_TYPE;
    typedef SequenceTypes::SEQUENCE_TYPE  SEQUENCE_TYPE;

private:
    SEQUENCE_TYPE org_sequence;
    SEQUENCE_TYPE header;

public:
    void Print();
    void Clear();
    int64_t ReadOneSeq(StreamLineReader& line_reader);
    SEQUENCE_TYPE& GetHeader();
    SEQUENCE_TYPE& GetSequence();
    int64_t      GetSeqLength();
    void ToUpperCase();
};

#endif // SEQUENCE_H_
