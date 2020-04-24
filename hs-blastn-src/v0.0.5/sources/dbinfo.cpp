#include "dbinfo.h"

#include <fstream>
#include <string>

const char* DbInfo::kClassName = "DbInfoGenerator";
static const char* kDbInfo_HeaderFileSufix = ".header";
static const char* kDbInfo_DatabaseFileSufix = ".sequence";

struct HeaderLine
{
    std::string* m_Line;
    std::string::size_type pos;

    HeaderLine()
    {
        m_Line = NULL;
        pos = std::string::npos;
    }
    void ChangeLine(std::string* new_line)
    {
        m_Line = new_line;
        pos = 0;
    }
    HeaderLine& operator++(void)
    {
        ++pos;
        return *this;
    }
    bool Pass(char delim)
    {
        while (pos < m_Line->size() && (*m_Line)[pos] != delim) ++pos;
        if (pos < m_Line->size())
        {
            ++pos;
            return true;
        }
        return false;
    }
    bool IsActive()
    {
        return pos < m_Line->size();
    }
	int GetInteger(DbInfo::SIZE_TYPE& n)
    {
		using namespace cy_utility;
		
        const char* p = &(*m_Line)[pos];
        int status;
		status = NString::StringToInteger(p, n);
        return status;
    }
};

////////////////////////////////////////////////////////////////////////////////////////

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

void DbInfo::GenerateDatabaseFileName(char name[])
{
    if (m_DatabaseName == NULL)
    {
        cy_utility::Log::ErrorAndExit(kClassName, "database file name is empty.");
    }

    strcpy(name, m_DatabaseName);
    strcat(name, kDbInfo_DatabaseFileSufix);
}

void DbInfo::DumpDatabase()
{
    char name[2048];
    GenerateDatabaseFileName(name);
	
	FILE* file = cy_utility::FileOperator::openfile(kClassName, name, "w");
    const char* src = (const char*)m_Database.get_data();
	cy_utility::FileOperator::write_file(kClassName, name, file, src, m_Database.size());

    file = cy_utility::FileOperator::closefile(file);
}

#include <iostream>

void DbInfo::LoadDatabase()
{
    using std::clog;

    char name[2048];
    GenerateDatabaseFileName(name);
	
	FILE* file = cy_utility::FileOperator::openfile(kClassName, name, "r");

    fseek(file, 0ULL, SEEK_END);
    m_DatabaseLength = ftell(file);
    fseek(file, 0ULL, SEEK_SET);

    SIZE_TYPE GB = 1;
    GB = GB << 30;
    double gb = 1.0 * m_DatabaseLength / GB;
    //clog << "\tLoading " << name << ", size = " << gb << "GB\n";
    fprintf(stderr, "\tLoading %s, size = %.1gGB\n", name, gb);

    m_Database.reserve(m_DatabaseLength);
	cy_utility::FileOperator::read_file(kClassName, name, file, m_Database.get_data(), m_DatabaseLength);

    file = cy_utility::FileOperator::closefile(file);
}

void DbInfo::GenerateHeaderFileName(char name[])
{
    if (m_DatabaseName == NULL)
    {
        cy_utility::Log::ErrorAndExit(kClassName, "database file name is empty.");
    }

    strcpy(name, m_DatabaseName);
    strcat(name, kDbInfo_HeaderFileSufix);
}

void DbInfo::DumpHeaders()
{
    using std::endl;

    char name[2048];
    GenerateHeaderFileName(name);

    std::ofstream file;
    file.open(name);
    if (!file)
    {
        cy_utility::Log::ErrorAndExit(kClassName, "cannot open file %s for writing.", name);
    }

    file << "Database name=" << m_DatabaseName << endl;
    file << "Number of sequences=" << GetNumSeqs() << endl;
    file << "Number of nucls=" << GetDbLength() << endl;

    SIZE_TYPE i, n = GetNumSeqs();
    for (i = 0; i < n; ++i)
    {
        file << "subject " << i+1 << endl;
        const char* id = GetSeqHeader(i);
        file << id << endl;

        SeqInfo& seqinfo = m_SequenceInfos[i];
        file << "gi=" << seqinfo.gi;
        file << " id_offset=" << seqinfo.id_offset;
        file << " seq_offset=" << seqinfo.sequence_offset;
        file << " seq_length=" << seqinfo.sequence_length;
        file << endl;
    }
	
    // dump posted date
    file << "Posted date:\n";
    file << "\ttm_sec=" << m_PostedDate.tm_sec << endl;
    file << "\ttm_min=" << m_PostedDate.tm_min << endl;
    file << "\ttm_hour=" << m_PostedDate.tm_hour << endl;
    file << "\ttm_mday=" << m_PostedDate.tm_mday << endl;
    file << "\ttm_mon=" << m_PostedDate.tm_mon << endl;
    file << "\ttm_year=" << m_PostedDate.tm_year << endl;
    file << "\ttm_wday=" << m_PostedDate.tm_wday << endl;
    file << "\ttm_yday=" << m_PostedDate.tm_yday << endl;	

    file.close();
}

void DbInfo::Destroy()
{
	m_SequenceInfos.destroy();
	m_Headers.destroy();
	m_Database.destroy();
	m_DatabaseLength = 0;
}

void DbInfo::LoadHeaders()
{
	using namespace cy_utility;
    if (m_DatabaseName == NULL)
    {
		Log::ErrorAndExit(kClassName, "database name is empty.");
    }

    const char kDelim = '=';
    using std::string;
    char name[2048];
    char buf[64];
    GenerateHeaderFileName(name);
    StreamLineReader line_reader(name);
    SIZE_TYPE line_number = 0;
    HeaderLine header_line;
    SIZE_TYPE num_seqs;
    SeqInfo seqinfo;
    SIZE_TYPE i;
	int status;
	SIZE_TYPE x;
	
#define CHECK_FILE_END if (line_reader.AtEof()) goto unexpected_end_of_file;	
#define READ_ONE_LINE_AND_CHECK if (!(++line_reader)) goto unexpected_end_of_file;

    Clear();
	
	std::string* tline;

    // database name
    READ_ONE_LINE_AND_CHECK

    // number of sequence
	READ_ONE_LINE_AND_CHECK
    line_number = line_reader.GetCurrentLineNumber();
    header_line.ChangeLine(&line_reader.GetLine());
    if (!header_line.Pass(kDelim) || !header_line.IsActive()) goto parsing_error;
	status = header_line.GetInteger(num_seqs);
	if (status) goto parsing_error;

    // database size
	READ_ONE_LINE_AND_CHECK
    line_number = line_reader.GetCurrentLineNumber();
    header_line.ChangeLine(&line_reader.GetLine());
    if (!header_line.Pass(kDelim) || !header_line.IsActive()) goto parsing_error;
	status = header_line.GetInteger(m_DatabaseLength);
	if (status) goto parsing_error;

    // read in each each subject
    for (i = 0; i < num_seqs; ++i)
    {
        // subject i
		READ_ONE_LINE_AND_CHECK

        //header line
		READ_ONE_LINE_AND_CHECK
		
        std::string& line = line_reader.GetLine();
		m_Headers.push_back(line.c_str(), line.size());
		m_Headers.push_back('\0');
		
#define READ_SEQUENCE_ITEM(item) \
  if (!header_line.Pass(kDelim) || !header_line.IsActive()) goto parsing_error; \
  status = header_line.GetInteger((item)); \
  if (status) goto parsing_error;

        // sequence info
		READ_ONE_LINE_AND_CHECK
        header_line.ChangeLine(&line_reader.GetLine());
        line_number = line_reader.GetCurrentLineNumber();
        // gi
		READ_SEQUENCE_ITEM(seqinfo.gi);

        // id offset
        ++header_line;
		READ_SEQUENCE_ITEM(seqinfo.id_offset);

        // sequence offset
        ++header_line;
		READ_SEQUENCE_ITEM(seqinfo.sequence_offset);

        // sequence length
        ++header_line;	
		READ_SEQUENCE_ITEM(seqinfo.sequence_length);

        m_SequenceInfos.push_back(seqinfo);
    }
	
	READ_ONE_LINE_AND_CHECK
	
#define READ_DATE_ITEM(item) \
  READ_ONE_LINE_AND_CHECK \
  header_line.ChangeLine(&line_reader.GetLine()); \
  line_number = line_reader.GetCurrentLineNumber(); \
  if (!header_line.Pass(kDelim) || !header_line.IsActive()) goto parsing_error; \
  status = header_line.GetInteger(x); \
  if (status) goto parsing_error; \
  (item) = x;
	
	READ_DATE_ITEM(m_PostedDate.tm_sec);
	
	READ_DATE_ITEM(m_PostedDate.tm_min);	
	
	READ_DATE_ITEM(m_PostedDate.tm_hour);
	
	READ_DATE_ITEM(m_PostedDate.tm_mday);
	
	READ_DATE_ITEM(m_PostedDate.tm_mon);
	
	READ_DATE_ITEM(m_PostedDate.tm_year);	
	
	READ_DATE_ITEM(m_PostedDate.tm_wday);
	
	READ_DATE_ITEM(m_PostedDate.tm_yday);

    return;

unexpected_end_of_file:
    Log::ErrorAndExit(kClassName, "unexptected end of file: %s.", name);

parsing_error:
    NString::IntegerToString(buf, line_number);
    Log::ErrorAndExit(kClassName, "format error at line %s of file %s.", buf, name);
}

void DbInfo::MakePostedDate()
{
    time_t rawtime;
    time(&rawtime);
    struct tm* timeinfo;
    timeinfo = localtime(&rawtime);
    
#define TIME_MACRO(item) m_PostedDate.item = timeinfo->item;
    TIME_MACRO(tm_sec)
    TIME_MACRO(tm_min)
    TIME_MACRO(tm_hour)
    TIME_MACRO(tm_mday)
    TIME_MACRO(tm_mon)
    TIME_MACRO(tm_year)
    TIME_MACRO(tm_wday)
    TIME_MACRO(tm_yday)
}

void DbInfo::OrgDBToBLASTNADB()
{
    SIZE_TYPE len = GetDbLength();
    SIZE_TYPE i;

    for (i = 0; i < len; ++i)
    {
        unsigned char c = m_Database[i];
        c = BLASTNA_TABLE[c];
        if (IsAmbigNucl(c)) ++m_NumAmbigNucls;
        m_Database[i] = c;
    }
}

DbInfo::SIZE_TYPE
DbInfo::GetSeqOffset(GI_SIZE_TYPE gi)
{
    return m_SequenceInfos[gi].sequence_offset;
}

const char* DbInfo::GetSeqHeader(GI_SIZE_TYPE gi)
{
	SIZE_TYPE offset = m_SequenceInfos[gi].id_offset;
	const char* r = &m_Headers[offset];
	return r;
}

DbInfo::SIZE_TYPE
DbInfo::GetNumSeqs()
{
    return m_SequenceInfos.size();
}

DbInfo::SIZE_TYPE
DbInfo::GetDbLength()
{
	return m_DatabaseLength;
}

DbInfo::SIZE_TYPE
DbInfo::GetSeqLength(GI_SIZE_TYPE gi)
{
    return m_SequenceInfos[gi].sequence_length;
}

void DbInfo::Print()
{
    PrintSummary();
    SIZE_TYPE i, n;
    n = GetNumSeqs();

    using std::cout;
    using std::endl;

    for (i = 0; i < n; ++i)
    {
        const char* id = GetSeqHeader(i);
        cout << "Subject " << i
             << ": " << id;
        cout << endl;

        cout << "\tLength = " << GetSeqLength(i) << endl;
    }
}

void DbInfo::PrintSummary()
{
    using std::clog;
    using std::endl;

    SIZE_TYPE percent_amb = 100 * m_NumAmbigNucls / GetDbLength();

    clog << "[" << kClassName << "] "
         << "Summary: " << endl
         << "\tNumber of subjects: " 
         << GetNumSeqs() << endl
         << "\tNumber of nucls: "
         << GetDbLength() << endl
         << "\tNumber of ambiguous nucls: "
         << m_NumAmbigNucls << " (" 
         << percent_amb << "%)"
         << endl;
}

bool DbInfo::IsAmbigNucl(unsigned char c)
{
    return c > 3;
}

#include <iostream>

void DbInfo::AddOneSeq(Sequence& subject)
{
    SeqInfo seqinfo;
    SIZE_TYPE i, n;
	
	// header
	seqinfo.gi = 0;
	seqinfo.id_offset = m_Headers.size();
	Sequence::SEQUENCE_TYPE& header = subject.GetHeader();
	m_Headers.push_back((char*)header.get_data(), header.size());

    // sequence
    seqinfo.sequence_offset = m_Database.size();
    SEQUENCE_TYPE& sequence = subject.GetSequence();
    seqinfo.sequence_length = sequence.size();
    m_Database.push_back(&sequence[0], sequence.size());
    m_SequenceInfos.push_back(seqinfo);
}

void DbInfo::BuildDbInfo()
{
	using namespace cy_utility;
	Log::LogMsg(kClassName, "building DbInfo.");
	Timer timer;
	timer.start();
	
    m_NumAmbigNucls = 0;
    Clear();
    if (m_DatabaseName == NULL)
    {
        cy_utility::Log::ErrorAndExit(kClassName, "file name is empty.");
    }

    StreamLineReader line_reader(m_DatabaseName);
    Sequence subject;

    while (1)
    {
        if (subject.ReadOneSeq(line_reader) == -1) break;
		const char* h = (const char*)subject.GetHeader().get_data();
        AddOneSeq(subject);
    }

    m_DatabaseLength = m_Database.size();

    OrgDBToBLASTNADB();

    DumpHeaders();
	DumpDatabase();
	
	timer.end();
	double dur = timer.get_elapsed_time();
	Log::LogMsg(kClassName, "done. Time elapsed: %f secs.", dur);
	PrintSummary();
	std::clog << std::endl << std::endl;
}

void DbInfo::LoadDbInfo()
{
    Clear();
    LoadHeaders();
	LoadDatabase();
}

void DbInfo::Clear()
{
    m_SequenceInfos.clear();
    m_Headers.clear();
    m_Database.clear();
}

void DbInfo::SetDatabaseName(const char* new_name)
{
    m_DatabaseName = new_name;
    Clear();
}

DbInfo::DbInfo(const char* database_name)
{
    m_DatabaseName = database_name;
}

DbInfo::SIZE_TYPE
DbInfo::GetSeqId(DbInfo::SIZE_TYPE offset)
{
    SIZE_TYPE num_seqs = m_SequenceInfos.size();
    SIZE_TYPE low = 0, high = num_seqs - 1;
    SIZE_TYPE mid;
    SIZE_TYPE start, length;
    while (low <= high)
    {
        mid = (low + high) / 2;
        start = m_SequenceInfos[mid].sequence_offset;
        length = m_SequenceInfos[mid].sequence_length;
        if (offset < start) high = mid - 1;
        else if (offset >= start + length) low = mid + 1;
        else break;
    }
    return mid;
}

char* DbInfo::GetDb()
{
	return (char*)m_Database.get_data();
}
