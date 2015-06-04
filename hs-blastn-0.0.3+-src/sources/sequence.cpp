#include "sequence.h"
#include <algorithm>
#include <cctype>
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

const char* FastaReader::kClassName = "[FASTA-Reader]";
FastaReader::SIZE_TYPE FastaReader::line_number = 0;

static unsigned char BLASTNA_TABLE[256] ={
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, //15
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, // 31
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 16, 16,  // 47
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 0, 10, 1, 11, 16, 16, 2, 12, 16, 16, 7, 16, 6, 14, 16,
	16, 16, 4, 9, 3, 16, 13, 8, 16, 5, 16, 16, 16, 16, 16, 16,
	16, 0, 10, 1, 11, 16, 16, 2, 12, 16, 16, 7, 16, 6, 14, 16,
	16, 16, 4, 9, 3, 16, 13, 8, 16, 5, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
};

/////////////////////////////////////////////////////////////////
//         FastaReader
/////////////////////////////////////////////////////////////////

bool FastaReader::IsNucl(unsigned char ch)
{
    int r = BLASTNA_TABLE[ch];
    return r < 16;
}

bool FastaReader::IsAmbigNucl(unsigned char ch)
{
    int r = BLASTNA_TABLE[ch];
    return (r < 16) && (r > 3);
}

void FastaReader::CheckDataLine(FastaReader::STRING_TYPE& line)
{
    SIZE_TYPE good = 0, bad = 0, len = line.size(), pos;

    for (pos = 0; pos < len; ++pos)
    {
        unsigned char c = line[pos];
        if (IsNucl(c))
        {
            ++good;
        }
        else if (c == '-')
        {
            ++good;
        }
        else if (isspace(c) || (c >= '0' && c <= '9'))
        {

        }
        else if (c == ';')
        {
            break;
        }
        else
        {
            ++bad;
        }
    }

    if (bad >= good / 3 && (len > 3 || good == 0 || bad > good))
    {
        cerr << kClassName << " Error: Near line "
             << line_number
             << ", there's a line that doesn't look like plausible data, "
             << "but it's not marked as defline or comment."
             << endl;
        exit(1);
    }
}

void FastaReader::ParseDataLine(FastaReader::STRING_TYPE& line, FastaReader::SEQUENCE_TYPE& sequence)
{
    CheckDataLine(line);
    SIZE_TYPE len = std::min(line.size(), line.find(';'));
    SIZE_TYPE pos;

    for (pos = 0; pos < len; ++pos)
    {
        unsigned char c = line[pos];
        if (IsNucl(c) || c == '-')
        {
            sequence.push_back(c);
        }
        else if (!isspace(c))
        {
            cerr << kClassName << " error: "
                 << "There are invalid residue(s) around position "
                 << pos + 1
                 << " of line "
                 << line_number
                 << endl;
                 exit(1);
        }
    }
}

void FastaReader::ParseDefLine(FastaReader::STRING_TYPE& line, FastaReader::SEQUENCE_TYPE& header)
{
    header.clear();
    SIZE_TYPE start = 0, end = line.size() - 1;
    while (start < line.size() && isspace(line[start])) ++start;
    while (end > start && isspace(line[end])) --end;
    if (start == end)
    {
        cerr << kClassName << " error: "
             << "There is a empty header line near line "
             << line_number
             << endl;
        exit(1);
    }
    header.push_back(&line[start + 1], end - start);
    header.push_back('\0');
}

bool FastaReader::IsEmptyLine(FastaReader::STRING_TYPE& line)
{
    SIZE_TYPE len = line.size();
    SIZE_TYPE pos;
    for (pos = 0; pos < len; ++pos)
    {
        if (!isspace(line[pos])) break;
    }
    if (pos == len) return true;
    return false;
}

char FastaReader::GetFirstLineChar(FastaReader::STRING_TYPE& line)
{
    SIZE_TYPE len = line.size(), pos = 0;
    for (pos = 0; pos < len; ++pos)
    {
        if (!isspace(line[pos])) break;
    }
	return line[pos];
}

///////////////////////////////////////////////////////////////////////////
// sequence
///////////////////////////////////////////////////////////////////////////

void Sequence::Print()
{
	const SIZE_TYPE kLineLength = 60;
    char* h = &header[0];
    cout << ">" << h << endl;
    SIZE_TYPE left = org_sequence.size();
	char* seq = (char*)org_sequence.get_data();
    while (left > 0)
    {
        SIZE_TYPE len = std::min(kLineLength, left);
        left -= len;
        for (SIZE_TYPE pos = 0; pos < len; ++pos)
        {
            cout << seq[pos];
        }
		seq += len;
        cout << endl;
    }
}

void Sequence::Clear()
{
    header.clear();
    org_sequence.clear();
}

FastaReader::SEQUENCE_TYPE&
Sequence::GetHeader()
{
    return header;
}

FastaReader::SEQUENCE_TYPE&
Sequence::GetSequence()
{
    return org_sequence;
}

FastaReader::SIZE_TYPE
Sequence::GetSeqLength()
{
    return org_sequence.size();
}

void Sequence::ToUpperCase()
{
#define __IsLower(c) ((c) >= 'a' && (c) <= 'z')
#define __ToUpper(c) ((c) - 'a' + 'A')

	SIZE_TYPE len = org_sequence.size();
	SIZE_TYPE i;
	for (i = 0; i < len; ++i)
	{
		unsigned char c = org_sequence[i];
		if (__IsLower(c)) org_sequence[i] = __ToUpper(c);
	}
}

Sequence::SIZE_TYPE
Sequence::ReadOneSeq(StreamLineReader& line_reader)
{
#define __check_comment_line(c) ((c) == '#' || (c) == '!')
    Clear();
    bool need_defline = true;

    while (!line_reader.AtEof())
    {
        ++line_reader;
        FastaReader::SetLineNumber(line_reader.GetCurrentLineNumber());
        STRING_TYPE& line = line_reader.GetLine();
        if (FastaReader::IsEmptyLine(line)) continue;

        char c = FastaReader::GetFirstLineChar(line);
        if (c == '>' || c == '@')
        {
            if (need_defline)
            {
                FastaReader::ParseDefLine(line, header);
                need_defline = false;
                continue;
            }
            else
            {
                line_reader.UngetLine();
                break;
            }
        }
        else if (c == '+')
        {
            ++line_reader;
            break;
        }
        else if (__check_comment_line(c))
        {
            continue;
        }

        FastaReader::ParseDataLine(line, org_sequence);
    }

    if (org_sequence.size() > 0 && header.size() == 0)
    {
        cerr << FastaReader::kClassName << " error: Near line "
             << FastaReader::line_number
             << ". Sequence missing header."
             << endl;
        exit(1);
    }
    if (org_sequence.size() == 0 && header.size() > 0)
    {
        cerr << FastaReader::kClassName << " error: Near line "
             << FastaReader::line_number
             << ". Sequence missing data."
             << endl;
        exit(1);
    }

    return org_sequence.size();
}
