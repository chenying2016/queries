#ifndef LINE_READER_H_
#define LINE_READER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <istream>
#include <algorithm>
#include <iterator>
#include <cstdlib>

#include "buffer_line_iterator.h"

using std::ifstream;
using std::cout;
using std::cin;
using std::clog;
using std::cerr;
using std::endl;
using std::string;
using std::istream;

/*
class LineIterator
{
public:
    typedef std::input_iterator_tag iterator_tag;
    typedef string                  value_type;
    typedef std::ptrdiff_t          difference_type;
    typedef const string*           pointer;
    typedef const string&           reference;

public:
    LineIterator() : in(NULL), isValid(false) {}
    LineIterator(istream& s) : in(&s), isValid((*in) ? true : false) {}
    reference operator*() { return line; }
    pointer operator->() const { return &line; }
    LineIterator operator++()
    {
        read();
        return *this;
    }
    LineIterator operator++(int)
    {
        LineIterator tmp = *this;
        read();
        return tmp;
    }

    void reset(istream* newStream)
    {
        in = newStream;
        if (*in)
            isValid = true;
    }

    bool eof()
    {
        return !isValid;
    }

    void clear()
    {
        in = NULL;
        line.clear();
        isValid = false;
    }

    string& getString()
    {
        return line;
    }

private:
    void read()
    {
        if (*in) getline(*in, line);
        isValid = (*in) ? true : false;
    }

private:
    istream* in;
    string line;
    bool isValid;
};
*/

class StreamLineReader
{
public:
    typedef BufferLineIterator       LineIterator;
    typedef LineIterator::value_type STRING_TYPE;
    typedef STRING_TYPE::size_type   SIZE_TYPE;

    static const char* kClassName;

public:
    StreamLineReader(const char* initFileName);
    ~StreamLineReader();

    bool operator++(void);
    STRING_TYPE& GetLine();
    void UngetLine();
    bool AtEof();
    SIZE_TYPE GetCurrentLineNumber();
    void ChangeFileName(const char* newFileName);
    void OpenFile();
    void CloseFile();
    void Clear();

private:
    LineIterator lineIterator;
    STRING_TYPE* currentLine;
    SIZE_TYPE    currentLineNumber;
    const char*  fileName;
    ifstream     file;
    bool         ungetLine;
};

#endif // LINE_READER_H_
