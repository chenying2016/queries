#include "search_worker.h"

using namespace cy_utility;

#include <iostream>
using namespace std;

/****************** Auxiliary datastructures and functions ***********************/

/* A subject will be split into several pieces in the preleminary stage if it is too long. */
/** Split subject sequences if longer than this */
static const Int4 MAX_DBSEQ_LEN = 5000000; 
/** By how much should the chunks of a subject sequence overlap if it is 
    too long and has to be split */
static const Int4 DBSEQ_CHUNK_OVERLAP = 100;

static inline Int4
GetSeedDiag(Int4 soff, Int4 qoff)
{
	return soff - qoff;
}

static inline Int8
GetSubjectBlock(Int8 offset)
{
	const Int8 r = (MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP);
	Int8 block = offset / r;
	return block;
}
static inline Int8
GetLastBlock(Int8 slen)
{
	if (slen <= MAX_DBSEQ_LEN) return 0;
	const Int8 r = (MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP);
	Int8 nb = slen / r;
	return nb;
}
static inline Int8
GetSubjectBlockOffset(Int8 block)
{
	const Int8 r = (MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP);
	Int8 offset = r * block;
	return offset;
}
static inline Int8
GetSubjectBlockEnd(Int8 block, Int8 subject_length)
{
	const Int8 r = (MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP);
	Int8 start = GetSubjectBlockOffset(block);
	Int8 end = std::min(start + MAX_DBSEQ_LEN, subject_length);
	return end;
}

static int mem_cmp_soff(const void* p1, const void* p2)
{
    const MEM* m1 = (const MEM*)p1;
    const MEM* m2 = (const MEM*)p2;
	if (m1->soff < m2->soff) return -1;
	if (m1->soff > m2->soff) return 1;

    return 0;    
}

struct MemCmpSoff
{
	bool operator()(const MEM& m1, const MEM& m2)
	{
		if (m1.soff < m2.soff) return true;
		if (m1.soff > m2.soff) return false;

		return false; 	
	}
};

static inline void
s_SetHitValue(Hit& h, Int4 qoff, Int4 soff, Int4 ext_l,
              Int4 score_l, Int4 ext_r, Int4 score_r, Int4 context)
{
    h.qoff = qoff - ext_l;
    h.soff = soff - ext_l;
    h.score = score_l + score_r;
    h.len = ext_l + ext_r;
    h.context = context;
}

#ifndef BLAST_CMP
#define BLAST_CMP(a, b) ((a) > (b) ? 1 : ((a) < (b) ? -1 : 0))
#endif

struct HitSortScore
{
	bool operator()(const Hit& pa, const Hit& pb)
	{
		if (pa.score > pb.score) return true;
		if (pa.score < pb.score) return false;
		if (pa.qoff < pb.qoff) return true;
		if (pa.qoff > pb.qoff) return false;
		return false;
	}
};

/** Callback for sorting an array of initial HSP structures (not pointers to
 *  * structures!) by score. 
 *   */
static int score_compare_match(const void *v1, const void *v2)
{
    Hit *h1, *h2;
    int result = 0;

    h1 = (Hit*) v1;
    h2 = (Hit*) v2;

    /* Check if ungapped_data substructures are initialized. If not, move
 *        those array elements to the end. In reality this should never happen. */
    if (h1 == NULL && h2 == NULL)
        return 0;
    else if (h1 == NULL)
        return 1;
    else if (h2 == NULL)
        return -1;

    if (0 == (result = BLAST_CMP(h2->score,
                                 h1->score)) &&
        0 == (result = BLAST_CMP(h1->soff,
                                 h2->soff)) &&
        0 == (result = BLAST_CMP(h2->len,
                                 h1->len)) &&
        0 == (result = BLAST_CMP(h1->qoff,
                                 h2->qoff))) {
        result = BLAST_CMP(h2->len, h1->len);
    }

    return result;
}

struct SeedCmpContextGiBlockDiagQoff
{
	bool operator()(const MEM& ts1, const MEM& ts2)
	{
		if (ts1.context < ts2.context) return true;
		if (ts1.context > ts2.context) return false;
		if (ts1.sid < ts2.sid) return true;
		if (ts1.sid > ts2.sid) return false;
		if (ts1.block_id < ts2.block_id) return true;
		if (ts1.block_id > ts2.block_id) return false;
		if (ts1.diag < ts2.diag) return true;
		if (ts1.diag > ts2.diag) return false;
		if (ts1.qoff < ts2.qoff) return true;
		if (ts1.qoff > ts2.qoff) return false;
		
		return false;	
	}
};

void PrintMEM(MEM& m)
{
	cout << "context = " << m.context << ", diag = " << m.diag << ", qoff = " << m.qoff << ", soff = " << m.soff << endl;
}

void PrintMEMList(SimpleArray<MEM>& seeds)
{
	Int4 num_seeds = seeds.size(), i;
	for (i = 0; i < num_seeds; ++i)
	{
		PrintMEM(seeds[i]);
	}
}

void PrintMEMList(MEM* seeds, Int4 num_seeds)
{
	Int4 i;
	for (i = 0; i < num_seeds; ++i)
	{
		PrintMEM(seeds[i]);
	}
}

void PrintHit(Hit& m)
{
	cout << "context = " << m.context << ", len = " << m.len << ", qoff = " << m.qoff << ", soff = " << m.soff << ", score = " << m.score << endl;
}

/*********************************** ungapped and score-only gapped extension functions *************/

static inline bool LeftSameUpTo(const Uint1* query, const Uint1* subject, Int4 qoff, Int4 soff, Int4 upto)
{
	Int4 avail = MIN(qoff, soff);
	if (avail < upto) return false;
	for (Int4 i = 1; i <= upto; ++i)
	{
		if (query[qoff - i] != Blastna2Na2(subject[soff - i])) return false;
	}
	return true;
}

void BuildOneSubjectTmpSeed(SimpleArray<MEM>& seeds,
							Int4 start_index,
							Int4 num_seeds,
							Int8 gi,
							Int8 subject_start_offset,
						    Int8 subject_length,
						    Int4 seed_size,
						    const Uint1* subject,
							QueryInfo* query_info)
{
	Int8 block;
	Int8 block_start, block_end, rbe;
	Int4 offset_block;
	Int4 context;
	Int4 context_offset;
	Int4 qoff;
	MEM m; m.soff = 0; m.diag = 0;
	Int8 last_block = GetLastBlock(subject_length);
	
	Int4 i1 = 0, i2;
	while (i1 < num_seeds)
	{
		MEM& m1 = seeds[start_index + i1];
		block = GetSubjectBlock(m1.soff);
		block_start = GetSubjectBlockOffset(block);
		block_end = GetSubjectBlockEnd(block, subject_length);
		rbe = block_end;
		if (block < last_block)
			block_end -= DBSEQ_CHUNK_OVERLAP;
		
		for (i2 = i1; i2 < num_seeds; ++i2)
		{
			MEM& m2 = seeds[start_index + i2];
			if (m2.soff > block_end)
			{
				break;	
			}
			
			m2.block_id = block;
			m2.block_offset = m2.soff - block_start;
			m2.diag = m2.block_offset - m2.qoff;

			if (m2.soff + seed_size > rbe)
			{
				const Uint1* query = query_info->GetSequence(m2.context);
				const Uint1* sbj = subject + block_start;
				Int4 upto = m2.soff + seed_size - rbe;
				if (LeftSameUpTo(query, sbj, m2.qoff, m2.block_offset, upto))
					m2.soff = 0;
				else 
					m2.soff = 1;
			}
			else m2.soff = 0;

			if (block > 0 && m2.block_offset < DBSEQ_CHUNK_OVERLAP - seed_size)
			{
				m.context = m2.context;
				m.sid = m2.sid;
				m.qoff = m2.qoff;
				m.block_id = m2.block_id - 1;
				m.block_offset = MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP + m2.block_offset;
				seeds.push_back(m);
			}
		}
		
		i1 = i2;
	}
}

void BuildTmpSeed(SimpleArray<MEM>& seeds,
				  Int4 seed_size,
				  DbInfo* dbinfo,
				  QueryInfo* query_info)
{
	std::sort(&seeds[0], &seeds[0] + seeds.size(), MemCmpSoff());
	Int4 i = 0, j, num_seeds = seeds.size();
	Int8 gi, subject_start, subject_end, subject_length;
	const Uint1* subject = (const Uint1*)dbinfo->GetDb();
	
	while (i < num_seeds)
	{
		gi = dbinfo->GetSeqId(seeds[i].soff);
		subject_length = dbinfo->GetSeqLength(gi);
		subject_start = dbinfo->GetSeqOffset(gi);
		subject_end = subject_start + subject_length;
		
		for (j = i; j < num_seeds; ++j)
		{
			if (seeds[j].soff >= subject_end) break;
			seeds[j].soff -= subject_start;
			seeds[j].sid = gi;
		}
		
		BuildOneSubjectTmpSeed(seeds,
							   i,
							   j - i,
							   gi,
							   subject_start,
							   subject_length,
							   seed_size,
							   subject + dbinfo->GetSeqOffset(gi),
							   query_info);
		
		i = j;
	}
}

static inline
void ExtendLeft(const Uint1* q, Int4 qoff, const Uint1* s, Int4 soff, Int4 avail,
                Int4** matrix, Int4 X, Int4& ext, Int4& score)
{
    Int4 sum = 0;
    Int4 q_beg = qoff + 1;
    score = 0;
    ext = qoff;
    
    while (avail > 0)
    {
        Uint1 c = Blastna2Na2(s[soff]);
        sum += matrix[q[qoff]][c];
        if (sum > 0)
        {
            q_beg = qoff;
            score += sum;
            sum = 0;
        }
        else if (sum < X)
        {
            break;
        }
        --qoff;
        --soff;
        --avail;
    }
    ext = ext - q_beg + 1;
}

static inline
void ExtendRight(const Uint1* q, Int4 qoff, const Uint1* s, Int4 soff, Int4 avail,
				 Int4** matrix, Int4 X, Int4& ext, Int4& score)
{
    Int4 sum = 0;
    Int4 q_beg = qoff - 1;
    score = 0;
    ext = qoff;
    
    while (avail > 0)
    {
        Uint1 c = Blastna2Na2(s[soff]);
        sum += matrix[q[qoff]][c];
        if (sum > 0)
        {
            q_beg = qoff;
            score += sum;
            sum = 0;
        }
        else if (sum < X)
        {
            break;
        }
        ++qoff;
        ++soff;
        --avail;
    }
    ext = q_beg - ext + 1;
}


static inline
void OneContextOneSubjectBlockUngappedExtension(SimpleArray<MEM>& seeds,
												Int4 start_index,
												Int4 num_mems,
												SimpleArray<Hit>& ungapped_alignments,
												const Uint1* query,
												Int4 query_length,
												Int4 context,
												Int4** matrix,
												Int4 X,
												Int4 cutoff_score,
												Uint1* subject,
												Int4 subject_length,
												Int8 gi,
												Int8 start_offset,
												Int4 seed_size)
{
	Int4 avail;
	Int4 qoff, soff;
	Int4 i;
	Int4 ext_l, ext_r, score_l, score_r;
	Hit hit;
	Int4 last_qoff = -1;
	Int4 last_diag = seeds[start_index].diag;
	
	ungapped_alignments.clear();
	
	for (i = 0; i < num_mems; ++i)
	{
		if (seeds[start_index + i].soff) 
		{  
			continue;
		}
		if (seeds[start_index + i].diag != last_diag)
		{
			last_diag = seeds[start_index + i].diag;
			last_qoff = -1;
		}
		else
		{
			if (seeds[start_index + i].qoff < last_qoff) continue;
		}
		
		qoff = seeds[start_index + i].qoff;
		soff = seeds[start_index + i].block_offset;
		avail = std::min(soff, qoff);
		ExtendLeft(query, qoff - 1, subject, soff - 1, avail, matrix, X, ext_l, score_l);
		
		avail = std::min((query_length - qoff), (subject_length - soff));
		ExtendRight(query, qoff, subject, soff, avail, matrix, X, ext_r, score_r);
		
		if (score_l + score_r >= cutoff_score)
		{
			s_SetHitValue(hit, qoff, soff, ext_l, score_l, ext_r, score_r, context);
			ungapped_alignments.push_back(hit);
			last_qoff = qoff + ext_r;
		}
	}
}

static inline
void OneContextOneSubjectBlockGetGappedScore(SimpleArray<Hit>& ungapped_alignments,
											 SimpleArray<HSP*>& gapped_alignments,
											 const Uint1* query,
											 Int4 query_length,
											 Int4 context,
											 Int4 cutoff_score,
											 const Uint1* subject,
											 Int4 subject_length,
											 Int8 subject_start_offset,
											 Int8 gi,
											 IntervalTree* itree,
											 Int4 min_diag_separation,
											 GreedyAligner* gapped_aligner,
											 SmallObjAllocator& soa)
{
	itree->Reset(0, query_length + 1, 0, subject_length + 1);
	HSP hsp;
	HSP* new_hsp = NULL;
	hsp.subject_id = gi;
	hsp.context = context;
	Int4 num_hits = ungapped_alignments.size();
	//std::sort(&ungapped_alignments[0], &ungapped_alignments[0] + num_hits, score_compare_match());
	qsort(&ungapped_alignments[0], ungapped_alignments.size(), sizeof(Hit), score_compare_match);
	Int4 last_num_gapped_alignments = gapped_alignments.size();
	for (int i = 0; i < num_hits; ++i)
	{
		hsp.q_off = ungapped_alignments[i].qoff;
		hsp.q_end = hsp.q_off + ungapped_alignments[i].len;
		hsp.s_off = ungapped_alignments[i].soff;
		hsp.s_end = hsp.s_off + ungapped_alignments[i].len;
		hsp.score = ungapped_alignments[i].score;
		
		if (itree->IntervalTreeContainsHSP(&hsp, min_diag_separation)) continue;
		
		Int4 qoff = hsp.q_off + ungapped_alignments[i].len / 2;
		Int4 soff = hsp.s_off + ungapped_alignments[i].len / 2;
		
		gapped_aligner->GreedyGappedAlignment(query, 
											  subject, 
											  query_length, 
											  subject_length,
											  qoff ,
											  soff, 
											  false, 
											  0);   
		
		if (gapped_aligner->gap_align.score < cutoff_score) continue;
		
		new_hsp = HSPNew(soa);
		gapped_aligner->PackHSP(*new_hsp, 0, 0);
		new_hsp->context = context;
		new_hsp->subject_id = gi;
		gapped_alignments.push_back(new_hsp);
		itree->IntervalTreeAddHSP(new_hsp);   
		new_hsp = NULL;
	}
	
	Int4 num_gapped_alignments = gapped_alignments.size() - last_num_gapped_alignments;
	Blast_HSPListPurgeHSPsWithCommonEndpoints(&gapped_alignments[last_num_gapped_alignments], num_gapped_alignments, TRUE, soa);
	gapped_alignments.set_array_size(last_num_gapped_alignments + num_gapped_alignments);

	for (int i = 0; i < num_gapped_alignments; ++i)
	{
		new_hsp = gapped_alignments[last_num_gapped_alignments + i];
		new_hsp->s_off += subject_start_offset;
		new_hsp->s_end += subject_start_offset;
		new_hsp->subject_gapped_start += subject_start_offset;
	}
}


void OneContextOneSubjectPrelimSearch(SimpleArray<MEM>& tmp_seeds,
									  Int4 start_index,
									  Int4 num_mems,
									  Int4 seed_size,
									  const Uint1* query,
									  Int4 query_length,
									  Int4 context,
									  Int4** matrix,
									  Int4 ungapped_xdrop,
									  Int4 ungapped_cutoff_score,
									  Int4 gapped_cutoff_score,
									  Uint1* subject,
									  Int8 subject_length,
									  Int8 gi,
									  SimpleArray<Hit>& ungapped_alignments,
									  SimpleArray<HSP*>& gapped_alignments,
									  IntervalTree* itree,
									  Int4 min_diag_separation,
									  GreedyAligner* gapped_aligner,
									  SmallObjAllocator& soa)
{
	Int4 i = 0, j;
	while (i < num_mems)
	{
		Int8 block = tmp_seeds[start_index + i].block_id;
		j = i + 1;
		while (j < num_mems)
		{
			if (tmp_seeds[start_index + j].block_id != block) break;
			++j;
		}
		ungapped_alignments.clear();
		Int8 block_start = GetSubjectBlockOffset(block);
		Int8 subject_left = subject_length - block_start;
		Int4 block_length = MAX_DBSEQ_LEN;
		if (block_length > subject_left)
			block_length = subject_left;
		OneContextOneSubjectBlockUngappedExtension(tmp_seeds,
													start_index + i,
													j - i,
													ungapped_alignments,
													query,
													query_length,
													context,
													matrix,
													ungapped_xdrop,
													ungapped_cutoff_score,
													subject + block_start,
													block_length,
													gi,
													block_start,
													seed_size);
		
		OneContextOneSubjectBlockGetGappedScore(ungapped_alignments,
												gapped_alignments,
												query,
												query_length,
												context,
												gapped_cutoff_score,
												subject + block_start,
												block_length,
												block_start,
												gi,
												itree,
												min_diag_separation,
												gapped_aligner,
												soa);
		
		i = j;
	}
}

void PrelimSearchStage(SimpleArray<MEM>& seeds,
					   QueryInfo* query_info,
					   DbInfo* dbinfo,
					   Int4 seed_size,
					   SimpleArray<Hit>& ungapped_alignments,
					   SimpleArray<HSP*>& gapped_alignments,
					   BlastInitialWordParameters* word_params,
					   BlastScoreBlk* sbp,
					   BlastHitSavingParameters* hit_params,
					   IntervalTree* itree,
					   Int4 min_diag_separation,
					   GreedyAligner* gapped_aligner,
					   SmallObjAllocator& soa)
{
	gapped_alignments.clear();

	BuildTmpSeed(seeds, seed_size, dbinfo, query_info);
	std::sort(&seeds[0], &seeds[0] + seeds.size(), SeedCmpContextGiBlockDiagQoff());

	const Uint1* query;
	Int4 query_length;
	Int4** score_matrix = sbp->matrix->data;
	Int4 ungapped_xdrop;
	Int4 ungapped_cutoff_score;
	Int4 gapped_cutoff_score;
	
	Int4 num_seeds = seeds.size(), i = 0, j;
	Int4 last_num_gapped_alignments = 0;
	Int4 last_context = seeds[0].context;
	
	while (i < num_seeds)
	{
		MEM& seed = seeds[i];
		Int4 context = seed.context;
		Int8 diag = seed.diag;
		
		if (context != last_context)
		{
			Int4 num_context_gapped_alignments = gapped_alignments.size() - last_num_gapped_alignments;
			if (num_context_gapped_alignments > 0)
			{
				gapped_alignments[last_num_gapped_alignments]->num_ident = num_context_gapped_alignments;
			}
			last_num_gapped_alignments = gapped_alignments.size();
			last_context = context;
		}
		
		j = i + 1;
		Int8 gi = seed.sid;
		while (j < num_seeds)
		{
			if (seeds[j].context != context) break;
			if (seeds[j].sid != gi) break;
			++j;
		}
		
		query = query_info->GetSequence(context);
		query_length = query_info->GetSeqLength(context);
		ungapped_xdrop = -word_params->cutoffs[context].x_dropoff;
		ungapped_cutoff_score = word_params->cutoffs[context].cutoff_score;
		gapped_cutoff_score = hit_params->cutoffs[context].cutoff_score;
		Uint1* subject = (Uint1*)dbinfo->GetDb();
		subject += dbinfo->GetSeqOffset(gi);
		Int8 subject_length = dbinfo->GetSeqLength(gi);
		
		OneContextOneSubjectPrelimSearch(seeds,
										 i,
										 j - i,
										 seed_size,
										 query,
										 query_length,
										 context,
										 score_matrix,
										 ungapped_xdrop,
										 ungapped_cutoff_score,
										 gapped_cutoff_score,
										 subject,
										 subject_length,
										 gi,
										 ungapped_alignments,
										 gapped_alignments,
										 itree,
										 min_diag_separation,
										 gapped_aligner,
										 soa);
		
		i = j;
	}
	
	Int4 num_context_gapped_alignments = gapped_alignments.size() - last_num_gapped_alignments;
	if (num_context_gapped_alignments > 0)
	{
		gapped_alignments[last_num_gapped_alignments]->num_ident = num_context_gapped_alignments;
	}
}
