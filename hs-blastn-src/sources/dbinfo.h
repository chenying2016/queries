#ifndef DBINFO_H_
#define DBINFO_H_ 1

#include "sequence.h"
#include "utility.h"
#include <vector>
#include <stdint.h>
#include <stdarg.h>
#include <cstdio>

struct DbInfoTypes
{
    typedef int64_t  Int8;
    typedef uint64_t Uint8;
    typedef int32_t  Int4;
    typedef uint32_t Uint4;

    typedef Uint8 SIZE_TYPE;
    typedef Uint8 GI_SIZE_TYPE;
};

struct SeqInfo
{
    typedef DbInfoTypes::SIZE_TYPE SIZE_TYPE;
    typedef DbInfoTypes::GI_SIZE_TYPE GI_SIZE_TYPE;

    GI_SIZE_TYPE gi;
    SIZE_TYPE id_offset;
    SIZE_TYPE sequence_offset;
    SIZE_TYPE sequence_length;
};

// This structure is used to store the date the database index is built.
struct DateStruct {
    int tm_sec; // seconds [0, 59]
    int tm_min; // minutes [0, 59]
    int tm_hour; // hour [0, 23]
    int tm_mday; // day of month [1, 31]
    int tm_mon; // month of year [0, 11]
    int tm_year; // years since 1900
    int tm_wday; // day of week [0, 6] (Sunday = 0)
    int tm_yday; // day of year [0, 365]
};

class DbInfo
{
public:
    typedef DbInfoTypes::SIZE_TYPE SIZE_TYPE;
    typedef DbInfoTypes::GI_SIZE_TYPE GI_SIZE_TYPE;
	typedef Sequence::SEQUENCE_TYPE SEQUENCE_TYPE;

private:
    const char*        m_DatabaseName;
    SIZE_TYPE          m_NumAmbigNucls;
    static const char* kClassName;

    cy_utility::SimpleArray<SeqInfo> m_SequenceInfos;
    SEQUENCE_TYPE                    m_Headers;
    SEQUENCE_TYPE                    m_Database;
    SIZE_TYPE                        m_DatabaseLength;
	DateStruct                       m_PostedDate;

private:
	// add one query to the database
    void AddOneSeq(Sequence& subject);
	// check if c is an ambiguous letter
    bool IsAmbigNucl(unsigned char c);
	// ASCII => BLASTNA
    void OrgDBToBLASTNADB();
    void DumpHeaders();
    void LoadHeaders();
    void DumpDatabase();
    void LoadDatabase();
    void GenerateDatabaseFileName(char name[]);
    void GenerateHeaderFileName(char name[]);

public:
    DbInfo(const char* database_name = NULL);
    void SetDatabaseName(const char* new_name);
    void Clear();
	void Print();
    void PrintSummary();
	// release the memory
	void Destroy();

    void BuildDbInfo();
    void LoadDbInfo();

	const char* GetDatabaseName() { return m_DatabaseName; }
	const char* GetSeqHeader(GI_SIZE_TYPE gi);
    SIZE_TYPE GetNumSeqs();
    SIZE_TYPE GetDbLength();
    SIZE_TYPE GetSeqLength(GI_SIZE_TYPE gi);
    SIZE_TYPE GetSeqOffset(GI_SIZE_TYPE gi);
    SIZE_TYPE GetSeqId(SIZE_TYPE offset);
    char* GetDb();
	DateStruct* GetPostedDate() { return &m_PostedDate; }
	void MakePostedDate();	
	SeqInfo* GetSeqInfos() { return (SeqInfo*)m_SequenceInfos.get_data(); }
};

#endif // DBINFO_H_
