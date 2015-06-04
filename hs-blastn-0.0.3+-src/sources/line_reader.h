#ifndef LINE_READER_H_
#define LINE_READER_H_

#include <string>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <stdarg.h>

#include "utility.h"

// for reading files, read one line at a time
class StreamLineReader
{
public:
    typedef std::string STRING_TYPE;
	typedef cy_utility::Types::SIZE_TYPE SIZE_TYPE;

private:
    static const SIZE_TYPE kBufferSize = ((SIZE_TYPE)1) << 20;
    static const char* kClassName;

	// The current line
    STRING_TYPE m_CurrentLine;
	// File buffer
    char        m_Buffer[kBufferSize];
	// Buffer size
    SIZE_TYPE   m_BufferSize;
	// Current position in buffer
    SIZE_TYPE   m_CurrentPosition;
	// How many bytes to be read in
    SIZE_TYPE   m_BytesToRead;
	// The number of current line
    SIZE_TYPE   m_CurrentLineNumber;
	// Name of file
    const char* m_FileName;
	// File
    FILE*       m_File;
	// This line is not dealt (because we come across the header of the next sequence)
    bool        m_UngetLine;

private:
    SIZE_TYPE ReadBuffer();

public:
	// Advance to the next line
    bool operator++(void);
	// Return the current line
    STRING_TYPE& GetLine();
	// This line has not been processed
    void UngetLine();
    StreamLineReader(const char* file_name);
    ~StreamLineReader();
	// End of file
    bool AtEof();
	// Return the line number
    SIZE_TYPE GetCurrentLineNumber();
	// Change the file
    void ChangeFileName(const char* new_file_name);
    void OpenFile();
    void CloseFile();
    void Clear();
};

#endif // LINE_READER_H_
