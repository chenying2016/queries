#include "query_info.h"
#include "sequence.h"
#include "utility.h"
#include "mask_misc.h"

static unsigned char BLASTNA_TABLE[256] ={
    /*  000 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  016 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  032 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  048 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  064 */ 15, 0 /* 65, A*/, 10, 1 /*67, C*/, 11, 15, 15, 2 /*71, G*/, 12, 15, 15, 7, 15, 6, 14, 15,
    /*  080 */ 15, 15, 4 /*82,R*/, 9, 3 /*84, T*/, 15, 13, 8, 15, 5 /*89, Y*/, 15, 15, 15, 15, 15, 15,
    /*  096 */ 15, 0, 10, 1, 11, 15, 15, 2, 12, 15, 15, 7, 15, 6, 14, 15,
    /*  112 */ 15, 15, 4, 9, 3, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15,
    /*  128 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  144 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  160 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  176 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  192 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  208 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  224 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    /*  240 */ 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
};

static unsigned char BLASTNA_COMPLETE_TABLE[256] ={
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 3, 10, 2, 11, 15, 15, 1, 12, 15, 15, 7, 15, 6, 14, 15,
    15, 15, 4, 9, 0, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15,
    15, 3, 10, 2, 11, 15, 15, 1, 12, 15, 15, 7, 15, 6, 14, 15,
    15, 15, 4, 9, 0, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
};

 const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
    'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 
    'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
};

void QueryInfo::Clear()
{
    first_context = last_context = -1;
    num_queries = 0;
    max_length = 0;
    
    blastna_query.clear();
    query_names.clear();
    contexts.clear();
    query_offsets.clear();
    org_queries.clear();
    scan_ranges.clear();
}

QueryInfo::QueryInfo(const char* fn) :       
        blastna_query(kMaxQueriesLength<<1),
        contexts(kMaxNumQueries<<1),
        query_offsets(kMaxNumQueries),
        org_queries(kMaxQueriesLength),
        scan_ranges(kMaxNumQueries<<1)
{
    ASSERT(fn != NULL);
    line_reader = new StreamLineReader(fn);
    Clear();
}

QueryInfo::QueryInfo() :       
	blastna_query(kMaxQueriesLength<<1),
	contexts(kMaxNumQueries<<1),
	query_offsets(kMaxNumQueries),
	org_queries(kMaxQueriesLength),
	scan_ranges(kMaxNumQueries<<1),
	line_reader(NULL) 
{
	Clear();
}

void QueryInfo::Init(const char* fn)
{
	if (line_reader == NULL)
	{
		line_reader = new StreamLineReader(NULL);
	}
	line_reader->Clear();
	line_reader->ChangeFileName(fn);
	line_reader->OpenFile();
	
	Clear();
}

void QueryInfo::Print()
{
    int i, j;
    
    char rq[601];
    rq[600] = '\0';
    
    for (i = 0; i < num_queries; ++i)
    {
        QueryOffset& offset = query_offsets[i];
        printf(">%s\n", &query_names[offset.header_offset]);
        
        Int4 qoff = offset.query_offset;
        Int4 qlen = offset.query_length;
        memcpy(rq, &org_queries[qoff], qlen);
        rq[qlen] = '\0';
        
        printf("\n\n");
        printf("%s\n", rq);
        printf("\nQuery Lenght = %d\n\n", offset.query_length);
    }
}

Int4 QueryInfo::GetQueryBatch(int num_threads)
{
    Clear();

    Int4 tot_len = 0;
    Sequence query;
    QueryOffset offset;
	offset.result_offset = 0;
	offset.num_alignments = 0;
    int64_t max_l = 0;
	
	const Int4 max_num_queries = kMaxNumQueries * num_threads;
	const Int4 max_queries_length = kMaxQueriesLength * num_threads;
    
    while (num_queries < max_num_queries && tot_len < max_queries_length)
    {
        if (query.ReadOneSeq(*line_reader) == -1) break;
		if (query.GetSeqLength() == 0) continue;
		query.ToUpperCase();
        
        offset.header_offset = query_names.size();
        offset.query_offset = org_queries.size();
		
		offset.header_offset = query_names.size();
		offset.query_offset = org_queries.size();
		Sequence::SEQUENCE_TYPE& header = query.GetHeader();
		query_names.push_back((char*)header.get_data(), header.size());
         
		Sequence::SEQUENCE_TYPE& org_seq = query.GetSequence();
        offset.query_length = query.GetSeqLength();
        org_queries.push_back(&org_seq[0], query.GetSeqLength());
        
        ++num_queries;
        tot_len += query.GetSeqLength();
        if (query.GetSeqLength() > max_l)
            max_l = query.GetSeqLength();
        query_offsets.push_back(offset);
    }

    max_length = max_l;
    return num_queries;
}

void QueryInfo::MakeBlastnaQuery()
{
    blastna_query.clear();
    contexts.clear();
    blastna_query.push_back(QUERY_DELIMITER);
    ContextInfo cinfo;
	cinfo.is_valid = true;
    Int4 offset;
    Int4 len;
    Int4 i, j;
	
	const Uint1 BLASTNA_REVERSE_COMPLEMENT_TABLE[16] = 
		{3, 2, 1, 0, 5, 4, 7, 6, 8, 9, 13, 12, 11, 10, 14, 15};
    
    Uint1* buf = new Uint1[max_length];
	Uint1* rc_buf = new Uint1[max_length];
    
    for (i = 0; i < num_queries; ++i)
    {
        offset = query_offsets[i].query_offset;
        len = query_offsets[i].query_length;
		
		cinfo.is_valid = true;
		if (len == 0) cinfo.is_valid = false;
        
        for (j = 0; j < len; ++j)
        {
            buf[j] = BLASTNA_TABLE[(Uint1)org_queries[offset + j]];
        }
        cinfo.offset = blastna_query.size();
        cinfo.length = len;
        contexts.push_back(cinfo);
        blastna_query.push_back(buf, len);
        blastna_query.push_back(QUERY_DELIMITER);
        
        for (j = len - 1; j >= 0; --j)
        {
			rc_buf[len - j - 1] = BLASTNA_REVERSE_COMPLEMENT_TABLE[buf[j]];
        }
        cinfo.offset = blastna_query.size();
        contexts.push_back(cinfo);
        blastna_query.push_back(rc_buf, len);
        blastna_query.push_back(QUERY_DELIMITER);
    }
    
    delete[] buf;
	delete[] rc_buf;
}

/**
 * Some comments on 
 * if (window_masker != NULL) { ... } else if (dust_masker != NULL) { ... }
 * In NCBI-BLASTN, the two maskers: WindowMasker and DustMasker can perfrom on queries simulataneously.
 * However, the results of WindowMasker will override the results produced by DustMasker.
 * In NCBI-BLAST 2.2.31+, this can be seen from:
 * 1) line 188, file dust_filter.cpp: queries.SetMaskedRegions(), the results of DustMasker are stored.
 * 2) line 315, file winmask_filter.cpp: query.SetMaskedRegions(), the results WindowMasker override the results of DustMasker.
 * We do not know if this is a bug.
 * To make sure that our results will be the same as NCBI-BLASTN, we make the if-else if cluase mentioned above.
 */

void QueryInfo::GetScanRanges(Int4 seed_size, CSymDustMasker* dust_masker, CSeqMasker* window_masker)
{
    CSeqMasker::TMaskList masked_locs;
    CSeqMasker::TMaskList unmasked_locs;
    
    Int4 i;
    for (i = 0; i < num_queries; ++i)
    {
		if (contexts[i*2].is_valid == false) continue;
        Uint4 offset = contexts[i<<1].offset;
        Uint4 len = query_offsets[i].query_length;
        char* q = &org_queries[query_offsets[i].query_offset];
        
        masked_locs.clear();
        unmasked_locs.clear();
        
        if (window_masker != NULL)
        {
            CSeqVector csv(q, len);
            (*window_masker)(csv, masked_locs);
            ComplementMaskLocations(len, seed_size, masked_locs, unmasked_locs);
        }
        else if (dust_masker != NULL)
        {
            dust_masker->GetMaskedLocs(q, len, masked_locs);
            ComplementMaskLocations(len, seed_size, masked_locs, unmasked_locs);
        }
        else
        {
            unmasked_locs.push_back(std::make_pair<Uint4, Uint4>(0, len - 1));
        }
        
        Point2d_Int4 srange;

        /// forward strand
        std::vector<std::pair<Uint4, Uint4> >::const_iterator cmiter;
        for (cmiter = unmasked_locs.begin(); cmiter != unmasked_locs.end(); ++cmiter)
        {
            CSymDustMasker::TMaskedInterval const& range = *cmiter;
            srange[0] = range.first + offset;
            srange[1] = range.second + offset;
            scan_ranges.push_back(srange);
        } 
    }
    
	cy_utility::Log::Trace(__func__, "Number of scan ranges: %d", scan_ranges.size());
}
