#ifndef QUERY_INFO_H
#define	QUERY_INFO_H

#include "def.h"
#include "line_reader.h"
#include "sequence.h"
#include "symdust.hpp"
#include "seq_masker.hpp"

#include <vector>

struct ContextInfo
{
    Int4 offset;
    Int4 length;
	Boolean is_valid;
};

struct QueryOffset
{
    Int4 header_offset;
    Int4 query_offset;
    Int4 query_length;
	Int4 result_offset;
	Int4 num_alignments;
};

struct QueryInfo
{
    Int4 first_context;
    Int4 last_context;
    Int4 num_queries;
    Int4 max_length;
	
	/// concatenated queries in BLASTNA format
	cy_utility::SimpleArray<Uint1>        blastna_query;
	// query names and comments
	cy_utility::SimpleArray<char>         query_names;
	// information for each context
	cy_utility::SimpleArray<ContextInfo>  contexts;
	// offset information
	cy_utility::SimpleArray<QueryOffset>  query_offsets;
	/// concatenated queries in ASCII format
	cy_utility::SimpleArray<char>         org_queries;
	// scan ranges, used for seeding only
	cy_utility::SimpleArray<Point2d_Int4> scan_ranges;
	// read queries from FASTA or FASTQ formated files
    StreamLineReader* line_reader;
    // read in a batch of queries from line_reader
    Int4 GetQueryBatch(int num_threads);
    QueryInfo(const char* fn);
	QueryInfo();
    ~QueryInfo() {if (line_reader) delete line_reader; line_reader = NULL;};
    void Clear();
	void Init(const char* fn);
    
	// print out query information
    void Print();
	// number of queries
    Int4 Size() { return num_queries; }
    Int4 Size() const { return num_queries; }
	// get query length
    Int4 GetSeqLength(Int4 context) { return contexts[context].length; }
    Int4 GetSeqLength(Int4 context) const { return contexts[context].length; }
	// get sequence
    const Uint1* GetSequence(Int4 context) { return &blastna_query[contexts[context].offset]; }
    const Uint1* GetSequence(Int4 context) const { return &blastna_query[contexts[context].offset]; }
	// get context offset
    Int4 GetContextOffset(Int4 context) { return contexts[context].offset; }
	// perform the windowmasker on queries, and get the unmasked regions (i.e., scan ranges)
    void GetScanRanges(Int4 seed_size, CSymDustMasker* dust_masker, CSeqMasker* window_masker);
	const char* GetSeqHeader(Int4 qid) { return &query_names[query_offsets[qid].header_offset]; }
	const char* GetSeqHeader(Int4 qid) const { return &query_names[query_offsets[qid].header_offset]; }
	// ASCII formated queries => BLASTNA formted queries
    void MakeBlastnaQuery();
    Int4 GetNumQueries() const { return num_queries; }
    
    static const Int4 kMaxNumQueries = 8000;
    static const Int4 kMaxQueriesLength = 40000000;
};

#endif	/* QUERY_INFO_H */

