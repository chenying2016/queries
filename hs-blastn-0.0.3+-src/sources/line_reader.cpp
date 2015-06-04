#include "line_reader.h"
#include "utility.h"

const char* StreamLineReader::kClassName = "StreamLineReader";

StreamLineReader::SIZE_TYPE StreamLineReader::GetCurrentLineNumber()
{
    return m_CurrentLineNumber;
}

StreamLineReader::SIZE_TYPE
StreamLineReader::ReadBuffer()
{
    if (m_File == NULL)
    {
        OpenFile();
    }

    SIZE_TYPE s = kBufferSize;
    if (s > m_BytesToRead) s = m_BytesToRead;
    cy_utility::FileOperator::read_file(kClassName, m_FileName, m_File, m_Buffer, s);

    m_CurrentPosition = 0;
    m_BufferSize = s;
    m_BytesToRead -= s;

    return s;
}

bool StreamLineReader::AtEof()
{
    if (m_UngetLine || m_CurrentPosition < m_BufferSize || m_BytesToRead > 0) return false;
    return true;
}

bool StreamLineReader::operator++(void)
{
    if (m_UngetLine) 
    {
        m_UngetLine = false;
        return true;
    }

    m_CurrentLine.clear();

	if (AtEof()) return false;

	bool found = false;
	while (!AtEof())
	{
        SIZE_TYPE i;
        for (i = m_CurrentPosition; i < m_BufferSize; ++i)
        {
            if (m_Buffer[i] == '\n') 
			{
				found = true;
				break;
			}
        }
        m_CurrentLine.append(m_Buffer + m_CurrentPosition, i - m_CurrentPosition);
        m_CurrentPosition = i + 1;
        if (found) break;
        if (!ReadBuffer()) break;
	}

	if (m_CurrentLine.size() > 0 && m_CurrentLine[m_CurrentLine.size() - 1] == '\r') m_CurrentLine.resize(m_CurrentLine.size() - 1);
	++m_CurrentLineNumber;

	return true;
}

StreamLineReader::STRING_TYPE& StreamLineReader::GetLine()
{
    return m_CurrentLine;
}

void StreamLineReader::UngetLine()
{
    m_UngetLine = true;
}

StreamLineReader::StreamLineReader(const char* file_name)
{
    m_File = NULL;

    if (file_name == NULL)
    {
		m_File = NULL;
        Clear();
    }
    else
    {
        m_FileName = file_name;
        OpenFile();
    }
}

StreamLineReader::~StreamLineReader()
{
    Clear();
}

void StreamLineReader::CloseFile()
{
    m_File = cy_utility::FileOperator::closefile(m_File);
}

void StreamLineReader::Clear()
{
    if (m_File) CloseFile();
    m_CurrentLine.clear();
    m_BufferSize = 0;
    m_CurrentPosition = 0;
    m_CurrentLineNumber = 0;
    m_BytesToRead = 0;
    m_UngetLine = false;
}

void StreamLineReader::ChangeFileName(const char* new_file_name)
{
    Clear();
    m_FileName = new_file_name;
}

void StreamLineReader::OpenFile()
{
    Clear();
    m_File = cy_utility::FileOperator::openfile(kClassName, m_FileName, "r");

    fseek(m_File, 0ULL, SEEK_END);
    m_BytesToRead = ftell(m_File);
    fseek(m_File, 0ULL, SEEK_SET);
}
