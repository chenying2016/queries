#include "line_reader.h"

const char* StreamLineReader::kClassName = "StreamLineReader";

StreamLineReader::SIZE_TYPE StreamLineReader::GetCurrentLineNumber()
{
    return this->currentLineNumber;
}

bool StreamLineReader::AtEof()
{
    if (this->ungetLine
        ||
        !this->lineIterator.eof())
        return false;
    return true;
}

bool StreamLineReader::operator++(void)
{
    if (this->ungetLine)
    {
        ungetLine = false;
        return true;
    }

    if (AtEof()) return false;

    const char kEnterKey = '\r';
    ++lineIterator;
    if (this->currentLine->size() > 0
        &&
        (*this->currentLine)[this->currentLine->size() - 1] == kEnterKey)
        this->currentLine->resize(this->currentLine->size() - 1);
    ++this->currentLineNumber;
	
    return !lineIterator.eof();
}

StreamLineReader::STRING_TYPE& StreamLineReader::GetLine()
{
    return *this->currentLine;
}

void StreamLineReader::UngetLine()
{
    this->ungetLine = true;
}

StreamLineReader::StreamLineReader(const char* initFileName)
    : lineIterator(), currentLine(&lineIterator.get_string()), currentLineNumber(0)
{
    if (initFileName != NULL)
    {
        ChangeFileName(initFileName);
        OpenFile();
    }
}

StreamLineReader::~StreamLineReader()
{
    Clear();
}

void StreamLineReader::CloseFile()
{
    if (file.is_open()) file.close();
}

void StreamLineReader::Clear()
{
    lineIterator.clear();
    CloseFile();
    currentLine->clear();
    ungetLine = false;
    currentLineNumber = 0;
}

void StreamLineReader::ChangeFileName(const char* newFileName)
{
    Clear();
    fileName = newFileName;
}

void StreamLineReader::OpenFile()
{
    Clear();
    file.open(fileName);
    if (!file)
    {
        cerr << "[" << kClassName << "] Fatal Error: Cannot open file " << fileName << " for reading." << endl;
        exit(1);
    }
    lineIterator.reset(&file);
}
