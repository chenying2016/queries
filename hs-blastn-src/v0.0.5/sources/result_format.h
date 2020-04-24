#ifndef RESULT_FORMAT_H
#define	RESULT_FORMAT_H

#include <string>
#include <fstream>

#include "dbinfo.h"
#include "gapalign.h"
#include "stat.h"
#include "arguments.h"
#include "parameters.h"
#include "output_buffer.h"
#include "query_info.h"

using std::string;

struct SubjectAlignInfo
{
    Int4 sid;
    Int4 start;
    Int4 end;
    double best_evalue;
    double bit_score;
};

typedef HSP** ResultSetType;
typedef cy_utility::SimpleArray<SubjectAlignInfo> SubjectResultSetType;

class FormatUtil
{
public:
    /// The string containing that no hits were found
    static const char kNoHitsFound[];

    static void PrintDbReport(DbInfo* dbinfo, int line_length, OutputBuffer& out, bool top);

    static void PrintKAParameters(double lambda, double k, double h,
                                  int line_len, OutputBuffer& out, bool gapped);

    static void AcknowledgeQuery(QueryInfo& query_info, 
                                 Int4 qid,
                                 int line_len, 
                                 OutputBuffer& out, 
                                 bool tabular);

    static void AddSpace(OutputBuffer& out, int number);

    static void GetScoreString(double evalue, double bit_score, int raw_score,
                               string& evalue_str, string& bit_score_str,
                               string& raw_score_str);

    static void Int8ToStringWithComma(string& out_str, Int8 value);

    static void DisplayDefline(OutputBuffer& out, DbInfo* dbinfo,
                               ResultSetType& results,
                               SubjectResultSetType&);

    static int GetInt8NumDigits(Int8 value);

    static void DisplayAlignments(QueryInfo& query_info,
                                   DbInfo* dbinfo,
                                   HSP& hsp,
                                   OutputBuffer& out);
    
    static void PrintOneQueryFooter(Int4 qid, BlastScoreBlk* sbp, OutputBuffer& out);

    static void MakeDateString(string& data_str, DateStruct* date);

    static void DisplayAlignInfo(OutputBuffer& out, HSP& hsp);

    static void DisplayIdentityInfo(OutputBuffer& out, HSP& hsp);

    static void DisplaySeqAlign(QueryInfo& query_info, 
                                 Int4 qid,
                                 DbInfo* dbinfo, 
                                 ResultSetType& results, 
                                 Int4 num_hsps,
                                 SubjectResultSetType& sinfo, 
                                 OutputBuffer& out);

    /// tabular
    static void PrintQueryAndDbNames(QueryInfo& query_info,
                                     Int4 qid,
                                     DbInfo* dbinfo,
                                     OutputBuffer& out);

    static void PrintHeader(QueryInfo& query_info,
                             Int4 qid,
                             DbInfo* dbinfo, 
                             Int4 num_hsps,
                             OutputBuffer& out);

    static void PrintFieldNames(OutputBuffer& out);

    static void PrintFields(OutputBuffer& out, HSP& hsp, QueryInfo& query_info, DbInfo* dbinfo);

    static void SortOneResults(HSP** results, Int4 num_hsps, SubjectResultSetType& sinfo);
};

class OutputFormat
{
private:
    OutputBuffer m_Outfile;
    OutputOptions* output_options;
	BlastHitSavingOptions* hit_options;
    DbInfo* dbinfo;

public:
    OutputFormat(OutputOptions* options, BlastHitSavingOptions* bhso, DbInfo* _dbinfo);
    ~OutputFormat();
	void Destroy() {m_Outfile.Destroy();}
	void Init(const char* file_name) {m_Outfile.Init(file_name);}
	void CloseFile() {m_Outfile.CloseFile();}
	void FlushResults(FILE* file) { m_Outfile.write_results(file); }
	void clear() { m_Outfile.clear(); }

    void PrintProlog();
    void PrintEpilog(Int8 num_queries, BlastScoringOptions* score_options);

    void PrintOneResult(QueryInfo& query_info,
                                  Int4 qid,
                                  HSP** results,
                                  Int4 num_hsps,
                                  BlastScoreBlk* sbp);

    void PrintTabularReport(QueryInfo& query_info, 
                                      Int4 qid, 
                                      DbInfo* dbinfo, 
                                      HSP** results,
                                      Int4 num_hsps,
                                      SubjectResultSetType& sinfos);
};

#endif	/* RESULT_FORMAT_H */

