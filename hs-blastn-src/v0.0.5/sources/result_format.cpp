#include "result_format.h"

/// margins constant
static string kOneSpaceMargin = " ";
static string kTwoSpaceMargin = "  ";

/// const strings
static const string kHeader = "Sequences producing significant alignments:";
static const string kScore = "Score";
static const string kE = "E";
static const string kBits = "(Bits)";
static const string kEvalue = "E value";
static const string kValue = "Value";
static const string kPercent = "Percent";
static const string kHighest = "Highest";
static const string kQuery = "Query";
static const string kCoverage = "Query coverage";
static const string kEllipsis = "...";

/// The line length of pairwise blast output
static const int kFormatLineLength = 68;

const char FormatUtil::kNoHitsFound[] = "No hits found";

static inline void PrintSeqGI(OutputBuffer& out, const char* header)
{
	int i = 0;
	while ((header[i] != '\0') && (!isspace(header[i]))) ++i;
	out.Append(header, i);
}

void FormatUtil::PrintDbReport(DbInfo* dbinfo, int line_length,
                               OutputBuffer& out, bool top)
{
    if (top)
    {
        out << "Database: ";
        const char* db_titles = dbinfo->GetDatabaseName();
        Int4 tot_num_seqs = dbinfo->GetNumSeqs();
        Int8 tot_length = dbinfo->GetDbLength();

        out << db_titles << "\n";
        FormatUtil::AddSpace(out, 11);

        string num;
        FormatUtil::Int8ToStringWithComma(num, tot_num_seqs);
        out << num << " sequences; ";
        FormatUtil::Int8ToStringWithComma(num, tot_length);
        out << num << " total letters";
        out << "\n\n";
        return;
    }

    out << "  Database: ";
    out << dbinfo->GetDatabaseName();
    out << "\n";
    out << "    Posted date: ";
    string str;
    FormatUtil::MakeDateString(str, dbinfo->GetPostedDate());
    out << str;
    out << "\n";
    out << "  Number of letters in database: ";
    FormatUtil::Int8ToStringWithComma(str, dbinfo->GetDbLength());
    out << str;
    out << "\n";
    out << "  Number of sequences in database: ";
    FormatUtil::Int8ToStringWithComma(str, dbinfo->GetNumSeqs());
    out << str;
    out << "\n\n";
}

void FormatUtil::PrintKAParameters(double lambda, double k, double h,
                                  int line_len, OutputBuffer& out, bool gapped)
{
    out << "\n";
    char buffer[256];
    if (gapped) out << "Gapped\n";

    out << "Lambda      K        H\n";
    sprintf(buffer, "%#8.3g ", lambda);
    out << buffer;
    sprintf(buffer, "%#8.3g ", k);
    out << buffer;
    sprintf(buffer, "%#8.3g ", h);
    out << buffer;
    out << "\n";
}

static void WrapOutputLine(std::string str, int line_len, OutputBuffer& out)
{	
	list<string> arr;
	SIZE_TYPE pos = 0, len = str.size(), nl_pos = 0;
	SIZE_TYPE width = line_len;
	
	enum EScore // worst to best
	{
		eForced,
		ePunct,
		eComma,
		eSpace,
		eNewline
	};
	
	typedef pair<SIZE_TYPE, SIZE_TYPE> TWrapSubstr;
	
	TWrapSubstr best_link(0, 0);
	TWrapSubstr latest_link(0, 0);
	
	while (pos < len)
	{
		bool hyphen       = false; // "-" or empty
		SIZE_TYPE column  = 0;
		SIZE_TYPE column0 = column;
		SIZE_TYPE best_pos = string::npos;
		EScore best_score = eForced;
		
		// certain logic can be skipped if this part has no backspace,
		// which is, by far, the most common case
		bool thisPartHasBackspace = false;
		
		arr.push_back("");
		arr.back().reserve(width);
		
		SIZE_TYPE pos0 = pos;
		
		if (nl_pos <= pos)
		{
			nl_pos = str.find('\n', pos);
			if (nl_pos == string::npos)
			{
				nl_pos = len;
			}
		}
		if (column + (nl_pos - pos) <= width)
		{
			pos0 = nl_pos;
		}
		
		for (SIZE_TYPE pos2 = pos0; pos2 < len && column <= width; ++pos2, ++column)
		{
			EScore score = eForced;
			SIZE_TYPE score_pos = pos2;
			const char c = str[pos2];
			
			if (c == '\n')
			{
				best_pos = pos2;
				best_score = eNewline;
				best_link = latest_link;
				break;
			}
			else if (isspace((unsigned char)c))
			{
				score = eSpace;
			}
			else if (c == ',' && column < width && score_pos < len - 1)
			{
				score = eComma;
				++score_pos;
			}
			else if (c == '-')
			{
				switch (c)
				{
					case '(': case '[': case '{': case '<': case '`':
						score = ePunct;
						break;
					default:
						if (score_pos < len - 1 && column < width)
						{
							score = ePunct;
							++score_pos;
						}
						break;
				}
			}
			
			if (score >= best_score && score_pos > pos0)
			{
				best_pos = score_pos;
				best_score = score;
				best_link = latest_link;
			}
			
			while (pos2 < len - 1 && str[pos2 + 1] == '\b')
			{
				++pos2;
				if (column > column0)
				{
					--column;
				}
				thisPartHasBackspace = true;
			}
		}
		
		if (best_score != eNewline && column <= width)
		{
			if (best_pos != len)
			{
				best_pos = len;
				best_link = latest_link;
				thisPartHasBackspace = true;
			}
		}
		
		{{
				string::const_iterator begin = str.begin() + pos;
				string::const_iterator end = str.begin() + best_pos;
				if (thisPartHasBackspace)
				{
					string::const_iterator bs;
					while ((bs = find(begin, end, '\b')) != end)
					{
						if (bs != begin)
						{
							arr.back().append(begin, bs - 1);
						}
						else
						{
							SIZE_TYPE size = arr.back().size();
							if (size > 0)
							{
								arr.back().resize(size - 1);
							}
						}
						//skip over backspace
						begin = bs + 1;
					}
				}
				if (begin != end)
				{
					arr.back().append(begin, end);
				}
		}}
		
		if (hyphen) arr.back() += '-';
		pos = best_pos;
		
		if (best_score == eSpace)
		{
			while (str[pos] == ' ') ++pos;
			if (str[pos] == '\n') ++pos;
		}
		if (best_score == eNewline) ++pos;
		
		while (pos < len && str[pos] == '\b') ++pos;
	}
	
	list<string>::iterator iter = arr.begin();
	while (iter != arr.end())
	{
		out << iter->c_str();
		out << "\n";
		iter++;
	}
}

void FormatUtil::AcknowledgeQuery(QueryInfo& query_info, 
                                  Int4 qid, 
                                  int line_len, 
                                  OutputBuffer& out, 
                                  bool tabular)
{
    if (tabular)
    {
        out << "# Query: ";
    }
    else
    {
        out << "Query= ";
    }
	
	std::string header = query_info.GetSeqHeader(qid);
	if (!tabular)
		WrapOutputLine(header, kFormatLineLength, out);
	else
		out << header;
    
    if (!tabular)
    {
        out << "\nLength=" << query_info.GetSeqLength(qid * 2);
        out << "\n";
    }
}

void FormatUtil::AddSpace(OutputBuffer& out, int number)
{
    for (int i = 0; i < number; ++i)
        out << " ";
}

void FormatUtil::GetScoreString(double evalue, 
								double bit_score, 
								int raw_score,
                                string& evalue_str, 
								string& bit_score_str,
                                string& raw_score_str)
{
    char evalue_buf[100], bit_score_buf[100];
    static string kBitScoreFormat("%4.1lf");

    if (evalue < 1.0e-180)
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "0.0");
    }
    else if (evalue < 1.0e-99)
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);
    }
    else if (evalue < 0.0009)
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
    }
    else if (evalue < 0.1)
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "%4.3lf", evalue);
    }
    else if (evalue < 1.0)
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "%3.2lf", evalue);
    }
    else if (evalue < 10.0)
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "%2.1lf", evalue);
    }
    else
    {
        snprintf(evalue_buf, sizeof(evalue_buf), "%2.0lf", evalue);
    }

    if (bit_score > 99999){
        snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
    } else if (bit_score > 99.9){
        snprintf(bit_score_buf, sizeof(bit_score_buf), "%3.0ld", (long)bit_score);
    } else {
        snprintf(bit_score_buf, sizeof(bit_score_buf), kBitScoreFormat.c_str(), bit_score);
    }

    evalue_str = evalue_buf;
    bit_score_str = bit_score_buf;

    if (raw_score <= 0)
        raw_score = -1;

    sprintf(evalue_buf, "%d", raw_score);
    raw_score_str = evalue_buf;
}

void FormatUtil::Int8ToStringWithComma(string& out_str, Int8 value)
{
    const int kBufSize = 64;
    char buffer[kBufSize + 1];
    buffer[kBufSize] = '\0';
    int cnt = 0;
    int index = kBufSize - 1;
    do
    {
        char c = value % 10 + 48;
        buffer[index--] = c;
        ++cnt;
        if (cnt == 3)
        {
            buffer[index--] = ',';
            cnt = 0;
        }
        value /= 10;
    } while (value > 0);

    if (cnt != 0)
        out_str = buffer + (1 + index);
    else
        out_str = buffer + (2 + index);
}

int FormatUtil::GetInt8NumDigits(Int8 value)
{
    int r = 0;
    do
    {
        ++r;
        value /= 10;
    } while (value > 0);
    return r;
}

static const char BLASTNA_TO_IUPACNA[16] = {
    'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K',
    'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
};

static const char BLASTNA_TO_IUPAC_COMPLEMENT_NA[16] = {
    'T', 'G', 'C', 'A', 'Y', 'R', 'K', 'M',
    'W', 'S', 'V', 'H', 'D', 'B', 'N', '-'
};

void FormatUtil::PrintOneQueryFooter(Int4 qid, BlastScoreBlk* sbp, OutputBuffer& out)
{
    const Blast_KarlinBlk* kbp_ungap = NULL;
    if (sbp->kbp != NULL)
        kbp_ungap = sbp->kbp[qid * 2];
    const Blast_KarlinBlk* kbp_gap = NULL;
    if (sbp->kbp_gap != NULL)
        kbp_gap = sbp->kbp_gap[qid * 2];

    out << "\n";
    if (kbp_ungap)
        FormatUtil::PrintKAParameters(kbp_ungap->Lambda,
                                      kbp_ungap->K,
                                      kbp_ungap->H,
                                      kFormatLineLength,
                                      out ,false);
    if (kbp_gap)
        FormatUtil::PrintKAParameters(kbp_gap->Lambda,
                                      kbp_gap->K,
                                      kbp_gap->H,
                                      kFormatLineLength,
                                      out, true);
    out << "\n";
    out << "Effective search space used: "
        << sbp->eff_searchsp[qid * 2]
	    << "\n";
}

static const char* months[12] =
{
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
};

void FormatUtil::MakeDateString(string& date_str, DateStruct* date)
{
    char buffer[256];
    int idx = 0;
    memcpy(buffer + idx, months[date->tm_mon], strlen(months[date->tm_mon]));
    idx += strlen(months[date->tm_mon]);
    buffer[idx++] = ' ';
    sprintf(buffer + idx, "%d", date->tm_mday);
    idx += FormatUtil::GetInt8NumDigits(date->tm_mday);
    buffer[idx++] = ',';
    buffer[idx++] = ' ';
    sprintf(buffer + idx, "%d", date->tm_year + 1900);
    idx += FormatUtil::GetInt8NumDigits(date->tm_year + 1900);
    buffer[idx++] = ' ';
    buffer[idx++] = ' ';
    int hour = date->tm_hour;
    if (hour >= 12) hour -= 12;
    sprintf(buffer + idx, "%d", hour);
    idx += FormatUtil::GetInt8NumDigits(hour);
    buffer[idx++] = ':';
    sprintf(buffer + idx, "%d", date->tm_min);
    idx += FormatUtil::GetInt8NumDigits(date->tm_min);
    buffer[idx++] = ' ';
    if (date->tm_hour >= 12)
    {
        sprintf(buffer + idx, "PM");
    }
    else
    {
        sprintf(buffer + idx, "AM");
    }

    date_str = buffer;
}

void FormatUtil::DisplayAlignInfo(OutputBuffer& out, HSP& hsp)
{
    string evalue_buf, bit_score_buf, raw_score_buf;
    FormatUtil::GetScoreString(hsp.evalue,
                               hsp.bit_score,
                               hsp.score,
                               evalue_buf,
                               bit_score_buf,
                               raw_score_buf);
    out << "\n";
    out << " Score = " << bit_score_buf << " bits";
    out << " (" << raw_score_buf <<"), ";
    out << "Expect = " << evalue_buf << "\n";
}

void FormatUtil::DisplayIdentityInfo(OutputBuffer& out, HSP& hsp)
{
    int align_len = 0;
    int i;
    int gaps = 0;
    GapEditScript* esp = hsp.esp;
    for (i = 0; i < esp->size; ++i)
    {
        if (esp->op_type[i] == eGapAlignIns ||
            esp->op_type[i] == eGapAlignDel)
            gaps += esp->num[i];
        align_len += esp->num[i];
    }

    int ident_percent = (1.0 * hsp.num_ident / align_len + .005) * 100;
    out << " Identities = " << hsp.num_ident << "/" << align_len;
    out << " (" << ident_percent << "%)";
    out << ", Gaps = " << gaps << "/" << align_len;
    int gap_percent = (1.0 * gaps / align_len + .005) * 100;
    out << "(" << gap_percent << "%)\n";
    out << " Strand=Plus/";
    if (hsp.context & 1)
    {
        out << "Minus\n";
    }
    else
    {
        out << "Plus\n";
    }
    out << "\n";
}

void FormatUtil::PrintFieldNames(OutputBuffer& out)
{
    out << "# Fields: ";
    out << "query id, ";
    out << "subject id, ";
    out << "% identity, ";
    out << "alignment length, ";
    out << "mismatches, ";
    out << "gap opens, ";
    out << "q. start, ";
    out << "q. end, ";
    out << "s. start, ";
    out << "s. end, ";
    out << "evalue, ";
    out << "bit score";

    out << "\n";
}

int HSPCmpFunc_sbjid_bitscore(const void* a, const void* b)
{
    const HSP** hp1 = (const HSP**)a;
    const HSP** hp2 = (const HSP**)b;
    const HSP* h1 = *hp1;
    const HSP* h2 = *hp2;
    
    if (h1->subject_id < h2->subject_id) return -1;
    else if (h1->subject_id > h2->subject_id) return 1;
    
    if (h1->bit_score > h2->bit_score) return -1;
    else if (h1->bit_score < h2->bit_score) return 1;
    else return 0;
}

/** Precision to which e-values are compared. */
#define FUZZY_EVALUE_COMPARE_FACTOR 1e-6

/** Compares 2 real numbers up to a fixed precision.
 * @param evalue1 First evalue [in]
 * @param evalue2 Second evalue [in]
 */
int 
s_FuzzyEvalueComp(double evalue1, double evalue2)
{
	if (evalue1 < 1.0e-180 && evalue2 < 1.0e-180)
		return 0;

	if (evalue1 < (1-FUZZY_EVALUE_COMPARE_FACTOR)*evalue2)
		return -1;
	else if (evalue1 > (1+FUZZY_EVALUE_COMPARE_FACTOR)*evalue2)
		return 1;
	else 
		return 0;
}

int SaiCmpFun_evalue_sbjid(const void* v1, const void* v2)
{
    const SubjectAlignInfo* sa = (const SubjectAlignInfo*)v1;
    const SubjectAlignInfo* sb = (const SubjectAlignInfo*)v2;
    
    ASSERT(sa->sid != sb->sid);
    
    int retval;
    
    if (sa->bit_score > sb->bit_score) return -1;
    if (sa->bit_score < sb->bit_score) return 1;
    
    retval = s_FuzzyEvalueComp(sa->best_evalue, sb->best_evalue);
    if(retval != 0) return retval;
    
    if (sa->sid < sb->sid) return -1;
    else return 1;
}

void FormatUtil::SortOneResults(HSP** results,
                                Int4 num_hsps,
                                SubjectResultSetType& sinfos)
{
    qsort(results, num_hsps, sizeof(HSP*), HSPCmpFunc_sbjid_bitscore);
    sinfos.clear();
    Int4 i = 0, j;
    while (i < num_hsps)
    {
        SubjectAlignInfo sinfo;
        sinfo.best_evalue = results[i]->evalue;
        sinfo.bit_score = results[i]->bit_score;
        sinfo.start = i;
        sinfo.sid = results[i]->subject_id;
        for (j = i + 1; j < num_hsps; ++j)
        {
            if (results[j]->subject_id != sinfo.sid)
                break;
        }
        sinfo.end = j - 1;
        sinfos.push_back(sinfo);
        i = j;
    }

	SubjectAlignInfo* psinfos = (SubjectAlignInfo*)sinfos.get_data();
	qsort(psinfos, sinfos.size(), sizeof(SubjectAlignInfo), SaiCmpFun_evalue_sbjid);
}

void FormatUtil::PrintHeader(QueryInfo& query_info,
                             Int4 qid,
                             DbInfo* dbinfo, 
                             Int4 num_hsps,
                             OutputBuffer& out)
{
    FormatUtil::PrintQueryAndDbNames(query_info, qid, dbinfo, out);
    if (num_hsps > 0)
    {
        FormatUtil::PrintFieldNames(out);
        out << "# " << num_hsps << " hits found\n";
    }
	else
	{
		out << "# " << num_hsps << " hits found\n";
	}
}

void FormatUtil::PrintQueryAndDbNames(QueryInfo& query_info,
                                      Int4 qid,
                                      DbInfo* dbinfo,
                                      OutputBuffer& out)
{
    out << "# BLASTN 2.3.0+\n";
    
    FormatUtil::AcknowledgeQuery(query_info, qid, kFormatLineLength, out, true);
    
    out << "\n# Database: "
        << (const char*)dbinfo->GetDatabaseName()
        << "\n";
}

static void s_WrapOutputLine(OutputBuffer& out, const string& str)
{
    const int line_len = 60;
    bool do_wrap = false;
    int length = (int) str.size();
    if (length > line_len) {
        for (int i = 0; i < length; i ++){
            if(i > 0 && i % line_len == 0){
                do_wrap = true;
            }   
            out << str[i];
            if(do_wrap && isspace((unsigned char) str[i])){
                out << "\n";  
                do_wrap = false;
            }
        }
    } else {
        out << str;
    }
}

void FormatUtil::DisplaySeqAlign(QueryInfo& query_info, 
                                 Int4 qid,
                                 DbInfo* dbinfo, 
                                 ResultSetType& results, 
                                 Int4 num_hsps,
                                 SubjectResultSetType& sinfo, 
                                 OutputBuffer& out)
{
    if (num_hsps == 0) return;
    
	Int4 sinfo_size = sinfo.size();
	char buf[64];
    // begin to display
    for (Int4 i = 0; i < sinfo_size; ++i)
    {
		string header = dbinfo->GetSeqHeader(sinfo[i].sid);
		cy_utility::Types::SIZE_TYPE seq_len = dbinfo->GetSeqLength(sinfo[i].sid);
		cy_utility::NString::IntegerToString(buf, seq_len);
        out << "\n\n";
        out << "> "; 
		s_WrapOutputLine(out, header);
        out << "\nLength=" << buf;
        out << "\n";
        for (Int4 j = sinfo[i].start; j <= sinfo[i].end; ++j)
        {
            HSP* hsp = results[j];
            FormatUtil::DisplayAlignInfo(out, *hsp);
            FormatUtil::DisplayIdentityInfo(out, *hsp);
            FormatUtil::DisplayAlignments(query_info, dbinfo, *hsp, out);
        }
    }
}

static inline void
OutputOneLine(OutputBuffer& out,
              char qalign[],
              char salign[],
              Int4 align_len,
              Int4 q_begin,
              Int4 q_end,
              Int8 s_begin,
              Int8 s_end,
              Int4 max_num_digits)
{
    int digit1 = FormatUtil::GetInt8NumDigits(q_begin);
    int digit2 = FormatUtil::GetInt8NumDigits(s_begin);
    Int4 i;
    // first line, query
    out << "Query";
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    out << q_begin;
    FormatUtil::AddSpace(out, max_num_digits - digit1);
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    for (i = 0; i < align_len; ++i)
    {
        out << qalign[i];
    }
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    out << q_end << "\n";
    
    /// middle line
    FormatUtil::AddSpace(out, strlen("Query") + max_num_digits + 4);
    for (i = 0; i < align_len; ++i)
    {
        if (qalign[i] == salign[i])
            out << '|';
        else
            out << kOneSpaceMargin;
    }
    out << "\n";
    
    /// third line, subject
    out << "Sbjct";
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    out << s_begin;
    FormatUtil::AddSpace(out, max_num_digits - digit2);
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    for (i = 0; i < align_len; ++i)
    {
        out << salign[i];
    }
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    out << s_end;
    out << "\n\n";
}

void FormatUtil::DisplayAlignments(QueryInfo& query_info,
                                   DbInfo* dbinfo,
                                   HSP& hsp,
                                   OutputBuffer& out)
{
	const Uint1 BLASTNA_REVERSE_COMPLEMENT_TABLE[16] = 
		{3, 2, 1, 0, 5, 4, 7, 6, 8, 9, 13, 12, 11, 10, 14, 15};
	
    const int kAlignLineLen = 60;
    char salign[kAlignLineLen], qalign[kAlignLineLen];
    int i, j;
    const Uint1* q;
    const Uint1* s;
    
    Int1 inc = 1;
    if (hsp.context & 1) inc = -1;
    
    Int4 qidx = 0;
    Int4 sidx = 0;
    Int4 line_left;
    Int4 esp_index;
    Int4 line_index;
    
    Int4 q_begin;
    Int8 s_begin;
    Int4 max_num_digits = MAX(FormatUtil::GetInt8NumDigits(hsp.q_end),
                              FormatUtil::GetInt8NumDigits(hsp.s_end));
    
    if (inc < 0)
    {
        s_begin = hsp.s_end - 1;
		s = (Uint1*)dbinfo->GetDb();
        s = s + dbinfo->GetSeqOffset(hsp.subject_id) + s_begin;
        
        q_begin = query_info.GetSeqLength(hsp.context) - hsp.q_end;
        q = query_info.GetSequence(hsp.context-1) + q_begin;
        esp_index = hsp.esp->size - 1;
    }
    else
    {
        s = (Uint1*)dbinfo->GetDb();
		s = s + dbinfo->GetSeqOffset(hsp.subject_id) + hsp.s_off;
        q = query_info.GetSequence(hsp.context) + hsp.q_off;
        esp_index = 0;
        s_begin = hsp.s_off;
        q_begin = hsp.q_off;
    }
	
	int line_index_q, line_index_s;
    
    while (qidx < hsp.q_end - hsp.q_off)
    {
        line_left = kAlignLineLen;
        line_index = 0;
		line_index_q = 0;
		line_index_s = 0;
        
        while (line_index < kAlignLineLen && qidx < hsp.q_end - hsp.q_off)
        {
            EGapAlignOpType op_type = hsp.esp->op_type[esp_index];
            Int4 num_ops = hsp.esp->num[esp_index];
            line_left = kAlignLineLen - line_index;
            if (num_ops > line_left)
            {
                hsp.esp->num[esp_index] -= line_left;
                num_ops = line_left;
            }
            else
            {
                esp_index += inc;
            }
            
            switch(op_type)
            {
            case eGapAlignSub:
            {
				line_index_q += num_ops;
				line_index_s += num_ops;				
                for (i = 0; i < num_ops; ++i)
                {
                    ASSERT(line_index < kAlignLineLen);
                    qalign[line_index] = BLASTNA_TO_IUPACNA[*q];
                    if (inc > 0)
                        salign[line_index] = BLASTNA_TO_IUPACNA[*s];
                    else
                        salign[line_index] = BLASTNA_TO_IUPAC_COMPLEMENT_NA[*s];
                    q++;
                    s += inc;
                    qidx++;
                    sidx++;
                    line_index++;
                }
                break;
            }
            case eGapAlignIns:
            {
				line_index_q += num_ops;
                for (i = 0; i < num_ops; ++i)
                {
                    qalign[line_index] = BLASTNA_TO_IUPACNA[*q];
                    salign[line_index] = '-';
                    q++;
                    ++qidx;
                    line_index++;
                }
                break;
            }
            case eGapAlignDel:
            {
				line_index_s += num_ops;
                for (i = 0; i < num_ops; ++i)
                {
                    qalign[line_index] = '-';
                    if (inc > 0)
                        salign[line_index] = BLASTNA_TO_IUPACNA[*s];
                    else
                        salign[line_index] = BLASTNA_TO_IUPAC_COMPLEMENT_NA[*s];
                    s += inc;
                    ++sidx;
                    line_index++;
                }
                break;
            }
            default:
                break;
            }
        }
        
		if (inc > 0)
			OutputOneLine(out, qalign, salign, 
						line_index, 
						q_begin + 1, 
						q_begin + line_index_q,
						s_begin + 1,
						s_begin + line_index_s,
						max_num_digits);
		else
			OutputOneLine(out, qalign, salign, 
						  line_index, 
						  q_begin + 1, 
						  q_begin + line_index_q,
						  s_begin + 1,
						  s_begin + 2 - line_index_s,
						  max_num_digits);			
        
        q_begin += line_index_q;
        s_begin += line_index_s * inc;
    }
}

void FormatUtil::PrintFields(OutputBuffer& out, 
                             HSP& hsp, 
                             QueryInfo& query_info,
                             DbInfo* dbinfo)
{
    const char kFieldDelimiter = '\t';
	PrintSeqGI(out, query_info.GetSeqHeader(hsp.context>>1));
    out << kFieldDelimiter;
	PrintSeqGI(out, dbinfo->GetSeqHeader(hsp.subject_id));
    out << kFieldDelimiter;

    int align_length = 0;
    int i;
    int gaps = 0;
    int gap_opens = 0;
    char buffer[256];

    for (i = 0; i < hsp.esp->size; ++i)
    {
        align_length += hsp.esp->num[i];
        if (hsp.esp->op_type[i] == eGapAlignIns ||
            hsp.esp->op_type[i] == eGapAlignDel)
        {
            gaps += hsp.esp->num[i];
            ++gap_opens;
        }
    }

    double percent = (1.0 * hsp.num_ident / align_length) * 100;
    
    Int4 qoff, qend;
    Int8 soff, send;
    
    if (hsp.context & 1)
    {
        qoff = query_info.GetSeqLength(hsp.context) - hsp.q_end;
        qend = query_info.GetSeqLength(hsp.context) - hsp.q_off;
        soff = hsp.s_end;
        send = hsp.s_off + 1;
    }
    else
    {
        qoff = hsp.q_off;
        qend = hsp.q_end;
        soff = hsp.s_off + 1;
        send = hsp.s_end;
    }
    
    ++qoff;

    sprintf(buffer, "%.*f", 3, percent);
    out << buffer;
    out << kFieldDelimiter;

    out << align_length;
    out << kFieldDelimiter;

    out << align_length - hsp.num_ident - gaps;
    out << kFieldDelimiter;

    out << gap_opens;
    out << kFieldDelimiter;

    out << qoff;
    out << kFieldDelimiter;

    out << qend;
    out << kFieldDelimiter;

    out << soff;
    out << kFieldDelimiter;

    out << send;
    out << kFieldDelimiter;

    string evalue_str, bit_score_str, raw_score_str;
    FormatUtil::GetScoreString(hsp.evalue, hsp.bit_score, hsp.score,
                               evalue_str, bit_score_str, raw_score_str);
	
	if ((hsp.evalue >= 1.0e-180) && (hsp.evalue < 0.0009))
	{
		char evalue_buf[cy_utility::kMaxDoubleStringSize];
		int n = cy_utility::NString::DoubleToString(hsp.evalue, 2, evalue_buf, cy_utility::kMaxDoubleStringSize, cy_utility::fDoubleScientific);
		evalue_buf[n] = '\0';
		evalue_str = evalue_buf;
	}	

    out << evalue_str;
    out << kFieldDelimiter;

    out << bit_score_str;
    
    out << "\n";    
}

void FormatUtil::DisplayDefline(OutputBuffer& out, 
                                DbInfo* dbinfo, 
                                ResultSetType& results, 
                                SubjectResultSetType& sinfos)
{
    const int m_MaxScoreLen = kBits.size();
    const unsigned int m_LineLen = kFormatLineLength;
    const int m_MaxEvalueLen = kValue.size();

    FormatUtil::AddSpace(out, m_LineLen + kTwoSpaceMargin.size());
    out << kScore;
    FormatUtil::AddSpace(out, m_MaxScoreLen - kScore.size());
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    FormatUtil::AddSpace(out, 2);
    out << kE;
    out << "\n";
    out << kHeader;
    FormatUtil::AddSpace(out, m_LineLen - kHeader.size());
    FormatUtil::AddSpace(out, kOneSpaceMargin.size());
    out << kBits;
    // in case m_MaxScoreLen > kBits.size()
    FormatUtil::AddSpace(out, m_MaxScoreLen - kBits.size());
    FormatUtil::AddSpace(out, kTwoSpaceMargin.size());
    out << kValue;
    out << "\n\n";

    string evalue_str, bit_score_str, raw_score_str;

	int num_infos = sinfos.size();
    for (int i = 0; i < num_infos; ++i)
    {
        unsigned int line_length = 0;
        string line_component = "  ";
        HSP& hsp = *results[sinfos[i].start];
        FormatUtil::GetScoreString(hsp.evalue, hsp.bit_score, hsp.score,
                                   evalue_str, bit_score_str, raw_score_str);
		line_component += dbinfo->GetSeqHeader(hsp.subject_id);
        string actual_line_component;
        if (line_component.size() + line_length > m_LineLen)
        {
            actual_line_component = line_component.substr(0, m_LineLen -
                    line_length - kEllipsis.size());
            actual_line_component += kEllipsis;
        }
        else
        {
            actual_line_component = line_component.substr(0, m_LineLen - line_length);
        }
        out << actual_line_component;
        line_length += actual_line_component.size();
        // pad the short lines
        FormatUtil::AddSpace(out, m_LineLen - line_length);
        out << kTwoSpaceMargin;
        out << bit_score_str;
        FormatUtil::AddSpace(out, m_MaxScoreLen - bit_score_str.size());
        out << kTwoSpaceMargin;
        out << evalue_str;
        out << "\n"       ;
    }    
}

OutputFormat::OutputFormat(OutputOptions* options, BlastHitSavingOptions* bhso, DbInfo* _dbinfo)
        : m_Outfile(),
        output_options(options),
		hit_options(bhso),
        dbinfo(_dbinfo)
{ }

OutputFormat::~OutputFormat()
{
    m_Outfile.Destroy();
}

void OutputFormat::PrintOneResult(QueryInfo& query_info,
                                  Int4 qid,
                                  HSP** results,
                                  Int4 num_hsps,
                                  BlastScoreBlk* sbp)
{
    SubjectResultSetType sinfos;
    FormatUtil::SortOneResults(results, num_hsps, sinfos);
    
    if (output_options->outfmt == eTabular ||
        output_options->outfmt == eTabularWithComments)
    {
        PrintTabularReport(query_info, qid, dbinfo, results, num_hsps, sinfos);
        return;
    }
    
    m_Outfile << "\n\n";
    FormatUtil::AcknowledgeQuery(query_info, qid, kFormatLineLength, m_Outfile, false);
    
    if (num_hsps == 0)
    {
        m_Outfile << "\n\n"
                  << "***** " << FormatUtil::kNoHitsFound << " *****" 
                  << "\n\n";
        FormatUtil::PrintOneQueryFooter(qid, sbp, m_Outfile);
        return;
    }
    
    FormatUtil::DisplayDefline(m_Outfile, dbinfo, results, sinfos);
    FormatUtil::DisplaySeqAlign(query_info, qid, dbinfo, results, num_hsps, sinfos, m_Outfile);
    FormatUtil::PrintOneQueryFooter(qid, sbp, m_Outfile);
}

void OutputFormat::PrintTabularReport(QueryInfo& query_info, 
                                      Int4 qid, 
                                      DbInfo* dbinfo, 
                                      HSP** results,
                                      Int4 num_hsps,
                                      SubjectResultSetType& sinfos)
{
    if (output_options->outfmt == eTabularWithComments)
        FormatUtil::PrintHeader(query_info, qid, dbinfo, num_hsps, m_Outfile);
    
    Int4 num_sbjs = sinfos.size();
    Int4 i, s, e, j;
	Int4 num_targets = std::min(num_sbjs, hit_options->hitlist_size);
    for (i = 0; i < num_targets; ++i)
    {
        s = sinfos[i].start;
        e = sinfos[i].end;
        for (j = s; j <= e; ++j)
        {
            HSP* hsp = results[j];
            FormatUtil::PrintFields(m_Outfile, *hsp, query_info, dbinfo);
        }
    }
}

void OutputFormat::PrintProlog()
{
    // no header for some output types
    if (output_options->outfmt >= eXml)
        return;
    
    FormatUtil::PrintDbReport(dbinfo, kFormatLineLength, m_Outfile, true);
}

void OutputFormat::PrintEpilog(Int8 num_queries, BlastScoringOptions* score_options)
{
    if (output_options->outfmt == eTabularWithComments)
    {
        m_Outfile << "# BLAST processed " << num_queries << " queries\n";
        return;
    }
    else if (output_options->outfmt >= eTabular)
        return;

    m_Outfile << "\n\n";
    FormatUtil::PrintDbReport(dbinfo, kFormatLineLength, m_Outfile, false);
    m_Outfile << "\n\nMatrix: blastn matrix ";
    m_Outfile << score_options->reward << " " << score_options->penalty;
    m_Outfile << "\n";

    if (score_options->gapped_calculation)
    {
        m_Outfile << "Gap Penalties: Existence: " << score_options->gap_open;
        m_Outfile << ", Extension: ";
        if (score_options->gap_open == 0 && score_options->gap_extend == 0)
        {
            float gap_extend = score_options->penalty - .5 * score_options->reward;
            char buf[64];
			sprintf(buf, "%.1f\n", gap_extend);
			m_Outfile << buf;
            //m_Outfile.AppendWithFormat("%.1f\n", gap_extend);
        }
        else
        {
            m_Outfile << score_options->gap_extend << "\n";
        }
    }
}
