#include "tabular_format.h"

#include "../../corelib/hbn_package_version.h"
#include "../../ncbi_blast/cmdline_args/format_flags.hpp"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <map>
#include <set>
#include <sstream>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(align_format)

using namespace std;

static void
get_score_string(double evalue,
    double bit_score,
    double total_bit_score,
    int raw_score,
    string& evalue_str,
    string& bit_score_str,
    string& total_bit_score_str,
    string& raw_score_str)
{
    char evalue_buf[100], bit_score_buf[100], total_bit_score_buf[100];

    /* Facilitates comparing formatted output using diff */
    static string kBitScoreFormat("%4.1lf");
#ifdef CTOOLKIT_COMPATIBLE
    static bool ctoolkit_compatible = false;
    static bool value_set = false;
    if ( !value_set ) {
        if (getenv("CTOOLKIT_COMPATIBLE")) {
            kBitScoreFormat.assign("%4.0lf");
            ctoolkit_compatible = true;
        }
        value_set = true;
    }
#endif /* CTOOLKIT_COMPATIBLE */
    
    if (evalue < 1.0e-180) {
        snprintf(evalue_buf, sizeof(evalue_buf), "0.0");
    } else if (evalue < 1.0e-99) {
        snprintf(evalue_buf, sizeof(evalue_buf), "%2.0le", evalue);        
#ifdef CTOOLKIT_COMPATIBLE
        if (ctoolkit_compatible) {
            strncpy(evalue_buf, evalue_buf+1, sizeof(evalue_buf-1));
        }
#endif /* CTOOLKIT_COMPATIBLE */
    } else if (evalue < 0.0009) {
        snprintf(evalue_buf, sizeof(evalue_buf), "%3.0le", evalue);
    } else if (evalue < 0.1) {
        snprintf(evalue_buf, sizeof(evalue_buf), "%4.3lf", evalue);
    } else if (evalue < 1.0) { 
        snprintf(evalue_buf, sizeof(evalue_buf), "%3.2lf", evalue);
    } else if (evalue < 10.0) {
        snprintf(evalue_buf, sizeof(evalue_buf), "%2.1lf", evalue);
    } else { 
        snprintf(evalue_buf, sizeof(evalue_buf), "%2.0lf", evalue);
    }
    
    if (bit_score > 99999){
        snprintf(bit_score_buf, sizeof(bit_score_buf), "%5.3le", bit_score);
    } else if (bit_score > 99.9){
        snprintf(bit_score_buf, sizeof(bit_score_buf), "%3.0ld",
            (long)bit_score);
    } else {
        snprintf(bit_score_buf, sizeof(bit_score_buf), kBitScoreFormat.c_str(),
            bit_score);
    }
    if (total_bit_score > 99999){
        snprintf(total_bit_score_buf, sizeof(total_bit_score_buf), "%5.3le", 
            total_bit_score);
    } else if (total_bit_score > 99.9){
        snprintf(total_bit_score_buf, sizeof(total_bit_score_buf), "%3.0ld",
            (long)total_bit_score);
    } else {
        snprintf(total_bit_score_buf, sizeof(total_bit_score_buf), "%2.1lf",
            total_bit_score);
    }
    evalue_str = evalue_buf;
    bit_score_str = bit_score_buf;
    total_bit_score_str = total_bit_score_buf;
    if (raw_score <= 0)
      raw_score = -1;
    NStr::IntToString(raw_score_str, raw_score);
}

static void
print_query_and_db_names(const char* query_name, const char* db_name, ostringstream& out)
{
    out << "# "
        << HBN_PACKAGE_NAME
        << ' '
        << HBN_PACKAGE_VERSION
        << '\n';
    
    out << "# Query: " << query_name << '\n';
    out << "# Database: " << db_name << '\n';
}

void x_PrintFieldNames(ostringstream& m_Ostream, vector<ETabularField>& m_FieldsToShow)
{
    m_Ostream << "# Fields: ";

    ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
        if (iter != m_FieldsToShow.begin())
            m_Ostream << ", ";

        switch (*iter) {
        case eQuerySeqId:
            m_Ostream << "query id"; break;
        case eQueryGi:
            m_Ostream << "query gi"; break;
        case eQueryAccession:
            m_Ostream << "query acc."; break;
        case eQueryAccessionVersion:
            m_Ostream << "query acc.ver"; break;
        case eQueryLength:
            m_Ostream << "query length"; break;
        case eSubjectSeqId:
            m_Ostream << "subject id"; break;
        case eSubjectAllSeqIds:
            m_Ostream << "subject ids"; break;
        case eSubjectGi:
            m_Ostream << "subject gi"; break;
        case eSubjectAllGis:
            m_Ostream << "subject gis"; break;
        case eSubjectAccession:
            m_Ostream << "subject acc."; break;
        case eSubjAccessionVersion:
            m_Ostream << "subject acc.ver"; break;
        case eSubjectAllAccessions:
            m_Ostream << "subject accs."; break;
        case eSubjectLength:
            m_Ostream << "subject length"; break;
        case eQueryStart:
            m_Ostream << "q. start"; break;
        case eQueryEnd:
            m_Ostream << "q. end"; break;
        case eSubjectStart:
            m_Ostream << "s. start"; break;
        case eSubjectEnd:
            m_Ostream << "s. end"; break;
        case eQuerySeq:
            m_Ostream << "query seq"; break;
        case eSubjectSeq:
            m_Ostream << "subject seq"; break;
        case eEvalue:
            m_Ostream << "evalue"; break;
        case eBitScore:
            m_Ostream << "bit score"; break;
        case eScore:
            m_Ostream << "score"; break;
        case eAlignmentLength:
            m_Ostream << "alignment length"; break;
        case ePercentIdentical:
            m_Ostream << "% identity"; break;
        case eNumIdentical:
            m_Ostream << "identical"; break;
        case eMismatches:
            m_Ostream << "mismatches"; break;
        case ePositives:
            m_Ostream << "positives"; break;
        case eGapOpenings:
            m_Ostream << "gap opens"; break;
        case eGaps:
            m_Ostream << "gaps"; break;
        case ePercentPositives:
            m_Ostream << "% positives"; break;
        case eFrames:
            m_Ostream << "query/sbjct frames"; break; 
        case eQueryFrame:
            m_Ostream << "query frame"; break; 
        case eSubjFrame:
            m_Ostream << "sbjct frame"; break; 
        case eBTOP:
            m_Ostream << "BTOP"; break;        
        case eSubjectTaxIds:
            m_Ostream << "subject tax ids"; break;
        case eSubjectSciNames:
        	m_Ostream << "subject sci names"; break;
        case eSubjectCommonNames:
        	m_Ostream << "subject com names"; break;
        case eSubjectBlastNames:
        	m_Ostream << "subject blast names"; break;
        case eSubjectSuperKingdoms:
        	m_Ostream << "subject super kingdoms"; break;
        case eSubjectTaxId:
            m_Ostream << "subject tax id"; break;
        case eSubjectSciName:
            m_Ostream << "subject sci name"; break;
        case eSubjectCommonName:
            m_Ostream << "subject com names"; break;
        case eSubjectBlastName:
           	m_Ostream << "subject blast name"; break;
        case eSubjectSuperKingdom:
               	m_Ostream << "subject super kingdom"; break;
        case eSubjectTitle:
        	m_Ostream << "subject title"; break;
        case eSubjectAllTitles:
        	m_Ostream << "subject titles"; break;
        case eSubjectStrand:
        	m_Ostream << "subject strand"; break;
        case eQueryCovSubject:
        	m_Ostream << "% query coverage per subject"; break;
        case eQueryCovUniqSubject:
        	m_Ostream << "% query coverage per uniq subject"; break;
        case eQueryCovSeqalign:
        	m_Ostream << "% query coverage per hsp"; break;
        default:
            _ASSERT(false);
            break;
        }
    }

    m_Ostream << "\n";
}

void x_AddDefaultFieldsToShow(map<string, ETabularField>& m_FieldMap, vector<ETabularField>& m_FieldsToShow)
{
    vector<string> format_tokens;
    NStr::Split(kDfltArgTabularOutputFmt, " ", format_tokens);
    for (const auto& iter : format_tokens) {
        hbn_assert(m_FieldMap.count(iter) > 0);
        ETabularField field = m_FieldMap[iter];
        m_FieldsToShow.push_back(field);
    }
}

void
print_tabular_header(const char* query_name,
    const char* db_name,
    BlastHitList* hit_list,
    map<string, ETabularField>& m_FieldMap,
    vector<ETabularField>& m_DefaultFieldsToShow,
    ostringstream& out)
{
    print_query_and_db_names(query_name, db_name, out);
    int num_hits = 0;
    for (int i = 0; i < hit_list->hsplist_count; ++i) num_hits += hit_list->hsplist_array[i]->hspcnt;
    if (num_hits) x_PrintFieldNames(out, m_DefaultFieldsToShow);
    out << "# " << num_hits << " hits found" << '\n';
}

static void
set_tabular_repoart_pos(const BlastHSP* hsp, int* qstart, int* qend, int* sstart, int* send)
{
    int qs = hsp->hbn_query.offset;
    int qe = hsp->hbn_query.end;
    int ql = hsp->hbn_query.seq_size;
    int ss = hsp->hbn_subject.offset;
    int se = hsp->hbn_subject.end;

    // If query is on the reverse strand, because BLAST output always reverses
    // subject, not query.
    if (hsp->hbn_query.strand == FWD) {
        *qstart = qs + 1;
        *qend = qe;
        *sstart = ss + 1;
        *send = se;
    } else {
        int x = ql - qe;
        int y = ql - qs;
        *qstart = x + 1;
        *qend = y;
        *sstart = se;
        *send = ss + 1;
    }
}

static void
print_tabular_report_for_one_hsp(const BlastHSP* hsp,
    const char* query_name,
    const char* subject_name,
    ostringstream& out)
{
    const char delim = '\t';
    string evalue_str;
    string bit_score_str;
    string total_bit_score_str;
    string raw_score_str;
    get_score_string(hsp->evalue,
        hsp->bit_score,
        0,
        hsp->score,
        evalue_str,
        bit_score_str,
        total_bit_score_str,
        raw_score_str);
    if ((hsp->evalue >= 1.0e-180) && (hsp->evalue < 0.0009)){
    	evalue_str = NStr::DoubleToString(hsp->evalue, 2, NStr::fDoubleScientific);
    }

    int qs, qe, ss, se;
    set_tabular_repoart_pos(hsp, &qs, &qe, &ss, &se);
    
    out << query_name << delim 
        << subject_name << delim
        << NStr::DoubleToString(hsp->hsp_info.perc_identity, 3) << delim
        << hsp->hsp_info.align_len << delim
        << hsp->hsp_info.align_len - hsp->hsp_info.num_ident << delim
        << hsp->hsp_info.gap_opens << delim
        << qs << delim
        << qe << delim
        << ss << delim
        << se << delim
        << evalue_str << delim
        << bit_score_str << '\n';
}

static void
print_tabular_report_for_one_hitlist(const char* query_name,
    const char* db_name,
    const CSeqDB* db,
    BlastHitList* hit_list,
    const EOutputFormat outfmt,
    map<string, ETabularField>& m_FieldMap,
    vector<ETabularField>& m_DefaultFieldsToShow,
    ostringstream& out)
{
    if (outfmt == eTabularWithComments) print_tabular_header(query_name, db_name, hit_list, m_FieldMap, m_DefaultFieldsToShow, out);
    for (int i = 0; i < hit_list->hsplist_count; ++i) {
        BlastHSPList* hsp_list = hit_list->hsplist_array[i];
        if (!hsp_list->hspcnt) continue;
        const int subject_id = hsp_list->hsp_array[0]->hbn_subject.oid;
        const char* subject_name = seqdb_seq_name(db, subject_id);
        for (int j = 0; j < hsp_list->hspcnt; ++j) {
            //HBN_LOG("*** oid = %d, hspcnt = %d", hsp_list->oid, hsp_list->hspcnt);
            BlastHSP* hsp = hsp_list->hsp_array[j];
            print_tabular_report_for_one_hsp(hsp, query_name, subject_name, out);
        }
    }
}

static const char*
extract_query_name(BlastHitList* hit_list, const CSeqDB* queries)
{
    const char* query_name = NULL;
    //HBN_LOG("hsplist: %d", hit_list->hsplist_count);
    for (int i = 0; i < hit_list->hsplist_count; ++i) {
        BlastHSPList* hsp_list = hit_list->hsplist_array[i];
        if (!hsp_list) continue;
        for (int j = 0; j < hsp_list->hspcnt; ++j) {
            BlastHSP* hsp = hsp_list->hsp_array[j];
            if (!hsp) continue;
            int oid = hsp->hbn_query.oid;
            query_name = seqdb_seq_name(queries, oid);
            break;
        }
    }
    return query_name;
}

extern "C"
void
print_tabular_reports(HbnHSPResults* results, const char* db_name, const CSeqDB* db, const CSeqDB* queries, const EOutputFormat outfmt)
{
    map<string, ETabularField> field_map;
    for (size_t i = 0; i < kNumTabularOutputFormatSpecifiers; i++) {
        field_map.insert(make_pair(sc_FormatSpecifiers[i].name,
                                    sc_FormatSpecifiers[i].field));
    }
    vector<ETabularField> default_fields;
    x_AddDefaultFieldsToShow(field_map, default_fields);
    ostringstream out;
    //HBN_LOG("number of queries: %d", results->num_queries);
    for (int i = 0; i < results->num_queries; ++i) {
        BlastHitList* hit_list = results->hitlist_array + i;
        const char* query_name = extract_query_name(hit_list, queries);
        if (!query_name) {
            //HBN_LOG("%d fail to find query name", query_name);
            continue;
        }
        print_tabular_report_for_one_hitlist(query_name, db_name, db, hit_list, outfmt, field_map, default_fields, out);
    }
    kstring_t* out_buf = &results->output_buf;
    ks_clear(*out_buf);
    string out_str = out.str();
    kputsn(out_str.c_str(), out_str.size(), out_buf);
}

END_SCOPE(align_format)
END_NCBI_SCOPE