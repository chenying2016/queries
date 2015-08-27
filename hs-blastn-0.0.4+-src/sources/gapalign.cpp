#include "gapalign.h"
#include "math.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;
using std::clog;

void GapPrelimEditBlock::Print()
{
	int i, ins = 0, del = 0, len = 0;
	for (i = 0; i < num_ops; ++i)
	{
		if (edit_ops[i].op_type == eGapAlignSub)
		{
			cout << "Subs: " << edit_ops[i].num << endl;
		}
		else if (edit_ops[i].op_type == eGapAlignIns)
		{
			cout << "Ins: " << edit_ops[i].num << endl;
			ins += edit_ops[i].num;
		}
		else if (edit_ops[i].op_type == eGapAlignDel)
		{
			cout << "Del: " << edit_ops[i].num << endl;
			del += edit_ops[i].num;
		}
		len += edit_ops[i].num;
	}
	cout << endl 
		 << "total Ins: " << ins << endl
		 << "total Del: " << del << endl
		 << "total length: " << len << endl;
}

void GapEditScript::Destroy()
{
  if (size <= 0) return;
  free(op_type); op_type = NULL;
  free(num); num = NULL;
  size = 0;
}

GapEditScript*
GapEditScriptNew(SmallObjAllocator& soa, Int4 size)
{
    GapEditScript* esp;
    if (size <= 0)
        return NULL;
    
    esp = (GapEditScript*)soa.Allocate(sizeof(GapEditScript));
    
    if (esp)
    {
        esp->size = size;
        esp->num_alloc = size;
        esp->op_type = (EGapAlignOpType*)soa.Allocate(sizeof(EGapAlignOpType) * size);
        esp->num = (Int4*)soa.Allocate(sizeof(Int4) * size);
    }
    
    return esp;
}

GapEditScript*
GapEditScriptDelete(SmallObjAllocator& soa, GapEditScript* esp)
{
    if (esp)
    {
        soa.Deallocate(esp->op_type, sizeof(EGapAlignOpType) * esp->num_alloc);
        soa.Deallocate(esp->num, sizeof(Int4) * esp->num_alloc);
        soa.Deallocate(esp, sizeof(GapEditScript));
    }
    return NULL;
}

void GapEditScriptPrint(GapEditScript* esp)
{
	fprintf(stdout, "number of ops: %d\n", esp->size);
	for (int i = 0; i < esp->size; ++i)
	{
		fprintf(stdout, "op_%d: %d, num: %d.\n", i+1, esp->op_type[i], esp->num[i]);
	}
}

Int2 GapEditScriptPartialCopy(GapEditScript* new_esp, 
                              int offset, 
                              const GapEditScript* old_esp,
                              int start,
                              int stop)
{
    int size = stop - start + 1;
    int new_index = 0;
    int old_index = 0;
    
    if (old_esp == NULL || new_esp == NULL || stop - start + 1 > new_esp->size)
        return -1;
    
    old_index = start;
    for (new_index = offset; new_index < size + offset; ++new_index)
    {
        new_esp->num[new_index] = old_esp->num[old_index];
        new_esp->op_type[new_index] = old_esp->op_type[old_index];
        old_index++;
    }
    
    return 0;
}

Int2 GapPrelimEditBlock::Realloc(Int4 total_ops)
{
  if (num_ops_allocated <= total_ops) {
    Int4 new_size = total_ops * 2;
    GapPrelimEditScript* new_ops;
    new_ops = (GapPrelimEditScript*)realloc(edit_ops, new_size * sizeof(GapPrelimEditScript));
    if (new_ops == NULL) return -1;
    
    edit_ops = new_ops;
    num_ops_allocated = new_size;
  }
  
  return 0;
}

Int2 GapPrelimEditBlock::AddNew(EGapAlignOpType op_type, Int4 num)
{
  if (Realloc(num_ops + 2) != 0)
    return -1;
  last_op = op_type;
  edit_ops[num_ops].op_type = op_type;
  edit_ops[num_ops].num = num;
  num_ops++;
  
  return 0;
}

void GapPrelimEditBlock::Add(EGapAlignOpType op_type, Int4 num)
{
  if (num == 0) return;
  
  if (last_op == op_type)
    edit_ops[num_ops-1].num += num;
  else
    AddNew(op_type, num);
}

GapPrelimEditBlock::GapPrelimEditBlock()
{
  edit_ops = NULL;
  num_ops_allocated = 0;
  num_ops = 0;
  last_op = eGapAlignInvalid;
  Realloc(100);
}

void GapPrelimEditBlock::Destroy()
{
  free(edit_ops); edit_ops = NULL;
  num_ops = num_ops_allocated = 0;
  last_op = eGapAlignInvalid;
}

GapPrelimEditBlock::~GapPrelimEditBlock()
{
  Destroy();
}

void GapPrelimEditBlock::Reset()
{
  num_ops = 0;
  last_op = eGapAlignInvalid;
}

GapAlignStruct::GapAlignStruct(BlastScoringParameters* score_params,
		                       BlastExtensionParameters* ext_params)
{
  assert(score_params != NULL && ext_params != NULL);
  edit_srcipt = NULL;
  //greedy_align_mem = new SGreedyAlignMem(score_params, ext_params);
  
  Int4 max_subject_length = GREEDY_MAX_COST;
  greedy_align_mem = new SGreedyAlignMem();
  greedy_align_mem->MemoryAlloc(score_params, ext_params, max_subject_length, 0);
  
  fwd_prelim_tback = new GapPrelimEditBlock();
  rev_prelim_tback = new GapPrelimEditBlock();

  query_start = 0;
  query_stop = 0;
  subject_start = 0;
  subject_stop = 0;
  greedy_query_seed_start = 0;
  greedy_subject_seed_start = 0;
  score = 0;

  bsp = score_params;
  bep = ext_params;
}

void GapAlignStruct::GapAlignStructFree()
{
  if (edit_srcipt) delete edit_srcipt; edit_srcipt = NULL;
  if (fwd_prelim_tback) delete fwd_prelim_tback; fwd_prelim_tback = NULL;
  if (rev_prelim_tback) delete rev_prelim_tback; rev_prelim_tback = NULL;
  if (greedy_align_mem) delete greedy_align_mem; greedy_align_mem = NULL;
}

GapAlignStruct::~GapAlignStruct()
{
  GapAlignStructFree();
}

/* hsp member functions */

/** Calculate number of identities in a regular HSP.
 * @param query The query sequence [in]
 * @param subject The uncompressed subject sequence [in]
 * @param hsp All information about the HSP [in]
 * @param num_ident_ptr Number of identities [out]
 * @param align_length_ptr The alignment length, including gaps [out]
 * @param sbp Blast score blk [in]
 * @param num_pos_ptr Number of Positives [out]
 * @return 0 on success, -1 on invalid parameters or error
 */

Int2 HSP::s_Blast_HSPGetNumIdentitesAndPositives(const Uint1* query, 
               const Uint1* subject, Int4* num_ident_ptr, 
               Int4* align_length_ptr,  
               Int4* num_pos_ptr)
{
    Int4 i, num_ident, align_length;
    Uint1 *q, *s;
    Int4 q_length = q_end - q_off;
    Int4 s_length = s_end - s_off;
    Int4** matrix = NULL;
    Int4 num_pos = 0;
    
    if (!subject || !query) 
        return -1;
    
    q = (Uint1*)&query[q_off];
    s = (Uint1*)&subject[s_off];
    
    num_ident = 0;
    align_length = 0;
    
    if (esp == NULL)
    {
        /* ungapped case. Check that lengths are the same in query and subject,
         then count the number of matches.*/
        if (q_length != s_length) return -1;
        align_length = q_length;
        for (i = 0; i < align_length; ++i)
        {
            if (*q == *s)
                num_ident++;

            else if (NULL != matrix)
            {
                if (matrix[*q][*s] > 0)
                    num_pos++;
            }
            ++q;
            ++s;
        }
    } else
    {
        Int4 index;
        for (index = 0; index < esp->size; ++index)
        {
            align_length += esp->num[index];
            switch (esp->op_type[index])
            {
            case eGapAlignSub:
                for (i = 0; i < esp->num[index]; ++i)
                {
                    if (*q == *s) ++num_ident;
                    else if (NULL != matrix)
                        if (matrix[*q][*s] > 0) ++num_pos;
                    ++q;
                    ++s;
                }
                break;
            case eGapAlignDel:
                s += esp->num[index];
                break;
            case eGapAlignIns:
                q += esp->num[index];
                break;
            default:
                s += esp->num[index];
                q += esp->num[index];
                break;
            }
        }
    }
    
    if (align_length_ptr)
        *align_length_ptr = align_length;
    *num_ident_ptr = num_ident;
    if (NULL != matrix)
        *num_pos_ptr = num_pos + num_ident;
    return 0;
}

HSP* HSPNew(SmallObjAllocator& soa)
{
    HSP* hsp = (HSP*)soa.Allocate(sizeof(HSP));
    hsp->esp = NULL;
    return hsp;
}

HSP* HSPDelete(SmallObjAllocator& soa, HSP* hsp)
{
    if (hsp == NULL) return NULL;
    
    if (hsp->esp != NULL)
    {
        GapEditScriptDelete(soa, hsp->esp);
        hsp->esp = NULL;
    }
    soa.Deallocate(hsp, sizeof(HSP));
    return NULL;
}

int GetNumIndels(GapEditScript* esp)
{
    int ret = 0;
    int i;
    for (i = 0; i < esp->size; ++i)
    {
        if (esp->op_type[i] == eGapAlignIns || esp->op_type[i] == eGapAlignDel)
        {
            ret += esp->num[i];
        }
    }
    return ret;
}

/* greedy gapped alignment */

#include "greedy_align_impl.h"

Int2 GreedyAligner::GreedyGappedAlignment(const Uint1* query, const Uint1* subject, 
					  Int4 query_length, Int4 subject_length, 
					  Int4 q_off, Int4 s_off, 
					  Boolean do_traceback, int mark)
{
  const Uint1* q;
  const Uint1* s;
  Int4 score;
  Int4 X;
  Int4 q_avail, s_avail;
  Int4 q_ext_l, q_ext_r, s_ext_l, s_ext_r;
  SGreedySeed fwd_start_point;
  SGreedySeed rev_start_point;
  Uint1 rem = 0;
  if (do_traceback) rem = 4;
  Int4 q_seed_start = q_off;
  Int4 s_seed_start = s_off;
  
  GapPrelimEditBlock *fwd_prelim_tback = NULL;
  GapPrelimEditBlock *rev_prelim_tback = NULL;
  
  q_avail = query_length - q_off;
  s_avail = subject_length - s_off;
  
  q = query + q_off;
  s = subject + s_off;
  
  //X = gap_align.bep->gap_x_dropoff;
  if (mark)
      X = gap_align.bep->gap_x_dropoff_final;
  else
      X = gap_align.bep->gap_x_dropoff;
  
  if (do_traceback) {
    fwd_prelim_tback = gap_align.fwd_prelim_tback;
    rev_prelim_tback = gap_align.rev_prelim_tback;
    fwd_prelim_tback->Reset();
    rev_prelim_tback->Reset();
  }
  
  /* extend to right */
  while (TRUE)
  {
        score = BLAST_AffineGreedyAlign(q, q_avail, s, s_avail, FALSE, X,
                                  gap_align.bsp->reward, -gap_align.bsp->penalty,
                                  gap_align.bsp->gap_open, gap_align.bsp->gap_extend,
                                  &q_ext_r, &s_ext_r, gap_align.greedy_align_mem,
                                  fwd_prelim_tback, rem, &fwd_start_point);

        if (score >= 0) break;

        /* double the max distance */
        Int4 new_dist, xdrop;
        new_dist = gap_align.greedy_align_mem->max_dist * 2;
        xdrop = gap_align.greedy_align_mem->xdrop;
        gap_align.greedy_align_mem->MemoryAlloc(gap_align.bsp, NULL, new_dist, xdrop);
  }
  
  /* extend to left */
  while (TRUE)
    {
        Int4 new_dist, xdrop, score1;

        score1 = BLAST_AffineGreedyAlign(query, q_off, subject, s_off, TRUE, X,
                                   gap_align.bsp->reward, -gap_align.bsp->penalty,
                                   gap_align.bsp->gap_open, gap_align.bsp->gap_extend,
                                   &q_ext_l, &s_ext_l, gap_align.greedy_align_mem,
                                   rev_prelim_tback, rem, &rev_start_point);
        if (score1 >= 0)
        {
            score += score1;
            break;
        }
        /* double the max distance */
        new_dist = gap_align.greedy_align_mem->max_dist * 2;
        xdrop = gap_align.greedy_align_mem->xdrop;
        gap_align.greedy_align_mem->MemoryAlloc(gap_align.bsp, NULL, new_dist, xdrop);        
    }
  
  /* In basic case the greedy algorithm returns number of differences, 
   * hence we need to convert it to score
   */
  if (gap_align.bsp->gap_open == 0 && gap_align.bsp->gap_extend == 0) {
    score = (q_ext_r + s_ext_r + q_ext_l + s_ext_l) * gap_align.bsp->reward/2
	     - score * (gap_align.bsp->reward - gap_align.bsp->penalty);
  } else if (gap_align.bsp->reward %2 == 1) {
    score /= 2;
  }
  
  if (do_traceback) {
    PrelimEditBlockToGapEditScript(rev_prelim_tback, fwd_prelim_tback);
    //ASSERT(!compressed_subject);
    
    // TODO check for possible gap elimination
//    if (esp)
//    {
//    	ReduceGaps(query+q_off-q_ext_l, subject+s_off - s_ext_l);
//    }
    if (esp) s_ReduceGaps(esp, query+q_off-q_ext_l, subject+s_off-s_ext_l,
                          query+q_off+q_ext_r, subject+s_off+s_ext_r);
  }
  else {
       /* estimate the best alignment start point. This is the middle
          of the longest run of exact matches found by the aligner,
          subject to the constraint that the start point lies in the
          box formed by the optimal alignment */

       Int4 q_box_l = q_off - q_ext_l;
       Int4 s_box_l = s_off - s_ext_l;
       Int4 q_box_r = q_off + q_ext_r;
       Int4 s_box_r = s_off + s_ext_r;
       Int4 q_seed_start_l = q_off - rev_start_point.start_q;
       Int4 s_seed_start_l = s_off - rev_start_point.start_s;
       Int4 q_seed_start_r = q_off + fwd_start_point.start_q;
       Int4 s_seed_start_r = s_off + fwd_start_point.start_s;
       Int4 valid_seed_len_l = 0;
       Int4 valid_seed_len_r = 0;

       if (q_seed_start_r < q_box_r && s_seed_start_r < s_box_r) {
           valid_seed_len_r = MIN(q_box_r - q_seed_start_r,
                                  s_box_r - s_seed_start_r);
           valid_seed_len_r = MIN(valid_seed_len_r,
                                  fwd_start_point.match_length) / 2;
       } else {
           q_seed_start_r = q_off;
           s_seed_start_r = s_off;
       }

       if (q_seed_start_l > q_box_l && s_seed_start_l > s_box_l) {
           valid_seed_len_l = MIN(q_seed_start_l - q_box_l,
                                  s_seed_start_l - s_box_l);
           valid_seed_len_l = MIN(valid_seed_len_l,
                                  rev_start_point.match_length) / 2;
       } else {
           q_seed_start_l = q_off;
           s_seed_start_l = s_off;
       }

       if (valid_seed_len_r > valid_seed_len_l) {
           q_seed_start = q_seed_start_r + valid_seed_len_r;
           s_seed_start = s_seed_start_r + valid_seed_len_r;
       }
       else {
           q_seed_start = q_seed_start_l - valid_seed_len_l;
           s_seed_start = s_seed_start_l - valid_seed_len_l;
       }
  }
  GreedyGapAlignStructFill(q_off-q_ext_l, s_off-s_ext_l,
                           q_off+q_ext_r, s_off+s_ext_r,
                           q_seed_start, s_seed_start,
                           score);	
  return 0;
}

GreedyAligner::GreedyAligner(BlastScoringParameters* score_params, 
        BlastExtensionParameters* ext_params,
        SmallObjAllocator& alloc) :
        gap_align(score_params, ext_params), esp(NULL), soa(alloc) {}

GreedyAligner::~GreedyAligner()
{
    if (esp != NULL) delete esp;
}

void
GreedyAligner::PackHSP(HSP& hsp, Int4 q_off, Int8 s_off)
{ 
    hsp.q_off = gap_align.query_start + q_off;
    hsp.q_end = gap_align.query_stop + q_off;
    hsp.s_off = gap_align.subject_start + s_off;
    hsp.s_end = gap_align.subject_stop + s_off;
    hsp.score = gap_align.score;
    hsp.query_gapped_start = gap_align.greedy_query_seed_start + q_off;
    hsp.subject_gapped_start = gap_align.greedy_subject_seed_start + s_off;
    
    if (esp) hsp.esp = esp;
    esp = NULL;
}

void
GreedyAligner::PrelimEditBlockToGapEditScript(GapPrelimEditBlock* rev_prelim_tback, 
        GapPrelimEditBlock* fwd_prelim_tback)
{
  Boolean merge_ops = FALSE;
  GapPrelimEditScript* op;
  Int4 i;
  Int4 index = 0;
  Int4 size = 0;
  
  if (rev_prelim_tback == NULL || fwd_prelim_tback == NULL) return ;
  
  if (fwd_prelim_tback->num_ops > 0 && rev_prelim_tback->num_ops > 0 &&
      fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].op_type ==
      rev_prelim_tback->edit_ops[(rev_prelim_tback->num_ops)-1].op_type)
    merge_ops = TRUE;
  
  size = fwd_prelim_tback->num_ops + rev_prelim_tback->num_ops;
  if (merge_ops) size--;
  
  //if (esp != NULL) delete esp;
  //esp = new GapEditScript(size);
  
  if (esp != NULL) esp = GapEditScriptDelete(soa, esp);
  esp = GapEditScriptNew(soa, size);
  
  index = 0;
  for(i = 0; i < rev_prelim_tback->num_ops; ++i) {
    op = rev_prelim_tback->edit_ops + i;
    esp->op_type[index] = op->op_type;
    esp->num[index] = op->num;
    index++;
  }
  
  if (fwd_prelim_tback->num_ops == 0) return;
  
  if (merge_ops)
    esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;
  
  if (merge_ops)
    i = fwd_prelim_tback->num_ops - 2;
  else
    i = fwd_prelim_tback->num_ops - 1;
  
  for( ; i >= 0; --i) {
    op = fwd_prelim_tback->edit_ops + i;
    esp->op_type[index] = op->op_type;
    esp->num[index] = op->num;
    index++;
  }   
}

Int2
GreedyAligner::GreedyGapAlignStructFill(Int4 q_start, Int4 s_start, Int4 q_end, 
                        Int4 s_end, Int4 q_seed_start, Int4 s_seed_start, Int4 score)
{
  gap_align.query_start = q_start;
  gap_align.query_stop = q_end;
  gap_align.subject_start = s_start;
  gap_align.subject_stop = s_end;
  gap_align.greedy_query_seed_start = q_seed_start;
  gap_align.greedy_subject_seed_start = s_seed_start;
  gap_align.score = score;
  
  return 0;    
}

/*
 hsp sort functions
 */
int
ScoreCompareHSPs(const void* h1, const void* h2)

{
   HSP* hsp1,* hsp2;   /* the HSPs to be compared */
   int result = 0;      /* the result of the comparison */
   
   hsp1 = *((HSP**) h1);
   hsp2 = *((HSP**) h2);

   /* Null HSPs are "greater" than any non-null ones, so they go to the end
      of a sorted list. */
   if (!hsp1 && !hsp2)
       return 0;
   else if (!hsp1)
       return 1;
   else if (!hsp2)
       return -1;

   if (0 == (result = BLAST_CMP(hsp2->score, hsp1->score)) &&
       0 == (result = BLAST_CMP(hsp1->s_off, hsp2->s_off)) &&
       0 == (result = BLAST_CMP(hsp2->s_end, hsp1->s_end)) &&
       0 == (result = BLAST_CMP(hsp1->q_off, hsp2->q_off))) {
       /* if all other test can't distinguish the HSPs, then the final
          test is the result */
       result = BLAST_CMP(hsp2->q_end, hsp1->q_end);
   }
   return result;
}

extern int s_FuzzyEvalueComp(double, double);

int 
EvalueCompareHSPs(const void* v1, const void* v2)
{
    const HSP** hp1;
    const HSP** hp2;
    const HSP* h1;
    const HSP* h2;
    hp1 = (const HSP**)v1;
    hp2 = (const HSP**)v2;
    
    h1 = *hp1;
    h2 = *hp2;
    
    if (h1 == NULL && h2 == NULL)
        return 0;
    else if (h1 == NULL)
        return 1;
    else if (h2 == NULL)
        return -1;
    
    int retval = 0;
    
    retval = s_FuzzyEvalueComp(h1->evalue, h2->evalue);
    if (retval != 0) return retval;
    
    return ScoreCompareHSPs(hp1, hp2);
}

/** Callback for sorting HSPs by starting offset in query. Sorting is by
 * increasing context, then increasing query start offset, then increasing
 * subject start offset, then decreasing score, then increasing query end 
 * offset, then increasing subject end offset. Null HSPs are moved to the 
 * end of the array.
 * @param v1 pointer to first HSP [in]
 * @param v2 pointer to second HSP [in]
 * @return Result of comparison.
 */
static int
s_QueryOffsetCompareHSPs(const void* v1, const void* v2)
{
   HSP* h1,* h2;
   HSP** hp1,** hp2;

   hp1 = (HSP**) v1;
   hp2 = (HSP**) v2;
   h1 = *hp1;
   h2 = *hp2;

   if (!h1 && !h2)
      return 0;
   else if (!h1) 
      return 1;
   else if (!h2)
      return -1;

   /* If these are from different contexts, don't compare offsets */
   if (h1->context < h2->context) 
      return -1;
   if (h1->context > h2->context)
      return 1;

   if (h1->q_off < h2->q_off)
      return -1;
   if (h1->q_off > h2->q_off)
      return 1;

   if (h1->s_off < h2->s_off)
      return -1;
   if (h1->s_off > h2->s_off)
      return 1;

   /* tie breakers: sort by decreasing score, then 
      by increasing size of query range, then by
      increasing subject range. */

   if (h1->score < h2->score)
      return 1;
   if (h1->score > h2->score)
      return -1;

   if (h1->q_end < h2->q_end)
      return 1;
   if (h1->q_end > h2->q_end)
      return -1;

   if (h1->s_end < h2->s_end)
      return 1;
   if (h1->s_end > h2->s_end)
      return -1;

   return 0;
}

/** Callback for sorting HSPs by ending offset in query. Sorting is by
 * increasing context, then increasing query end offset, then increasing
 * subject end offset, then decreasing score, then decreasing query start
 * offset, then decreasing subject start offset. Null HSPs are moved to the 
 * end of the array.
 * @param v1 pointer to first HSP [in]
 * @param v2 pointer to second HSP [in]
 * @return Result of comparison.
 */
static int
s_QueryEndCompareHSPs(const void* v1, const void* v2)
{
   HSP* h1,* h2;
   HSP** hp1,** hp2;

   hp1 = (HSP**) v1;
   hp2 = (HSP**) v2;
   h1 = *hp1;
   h2 = *hp2;

   if (!h1 && !h2)
      return 0;
   else if (!h1) 
      return 1;
   else if (!h2)
      return -1;

   /* If these are from different contexts, don't compare offsets */
   if (h1->context < h2->context) 
      return -1;
   if (h1->context > h2->context)
      return 1;

   if (h1->q_end < h2->q_end)
      return -1;
   if (h1->q_end > h2->q_end)
      return 1;

   if (h1->s_end < h2->s_end)
      return -1;
   if (h1->s_end > h2->s_end)
      return 1;

   /* tie breakers: sort by decreasing score, then 
      by increasing size of query range, then by
      increasing size of subject range. The shortest range 
      means the *largest* sequence offset must come 
      first */
   if (h1->score < h2->score)
      return 1;
   if (h1->score > h2->score)
      return -1;

   if (h1->q_off < h2->q_off)
      return 1;
   if (h1->q_off > h2->q_off)
      return -1;

   if (h1->s_off < h2->s_off)
      return 1;
   if (h1->s_off > h2->s_off)
      return -1;

   return 0;
}

/* cut off the GapEditScript according to hsp offset and end */
static void 
s_CutOffGapEditScript(HSP* hsp, Int4 q_cut, Int4 s_cut, Boolean cut_begin)
{   
   int index, opid, qid, sid;
   Boolean found = FALSE;
   GapEditScript *esp = hsp->esp;
   qid = 0;
   sid = 0;
   opid = 0;
   q_cut -= hsp->q_off;
   s_cut -= hsp->s_off;
   for (index=0; index < esp->size; index++) {
       for(opid=0; opid < esp->num[index];){  
          if (esp->op_type[index] == eGapAlignSub) {
             qid++;
             sid++;
             opid++;
          } else if (esp->op_type[index] == eGapAlignDel) {
             sid+=esp->num[index];
             opid+=esp->num[index]; 
          } else if (esp->op_type[index] == eGapAlignIns) {
             qid+=esp->num[index];
             opid+=esp->num[index];
          }
          if (qid >= q_cut && sid >= s_cut) found = TRUE;
          if (found) break;
       }
       if (found) break;
   }
   
   // RMH: Unless both cut sites where found the following
   //      block would access memory outside the GapEditScript
   //      array.
   if ( found )
   {
     if (cut_begin) {
         int new_index = 0;
         if (opid < esp->num[index]) {
            ASSERT(esp->op_type[index] == eGapAlignSub);
            esp->op_type[0] = esp->op_type[index];
            esp->num[0] = esp->num[index] - opid;
            new_index++;
         } 
         ++index;
         for (; index < esp->size; index++, new_index++) {
            esp->op_type[new_index] = esp->op_type[index];
            esp->num[new_index] = esp->num[index];
         }
         esp->size = new_index;
         hsp->q_off += qid;
         hsp->s_off += sid;
     } else {
         if (opid < esp->num[index]) {
            ASSERT(esp->op_type[index] == eGapAlignSub);
            esp->num[index] = opid;
         } 
         esp->size = index+1;
         hsp->q_end = hsp->q_off + qid;
         hsp->s_end = hsp->s_off + sid;
     }
   }
}

Int4
Blast_HSPListPurgeHSPsWithCommonEndpoints(HSP** hsp_list,
                                          Int4& num_hsps,
                                          Boolean purge,
                                          SmallObjAllocator& soa)

{
   HSP** hsp_array;  /* hsp_array to purge. */
   HSP* hsp;
   Int4 i, j, k;
   Int4 hsp_count;
   
   /* If HSP list is empty, return immediately. */
   if (hsp_list == NULL || num_hsps == 0)
       return 0;

   hsp_array = hsp_list;
   hsp_count = num_hsps;

   qsort(hsp_array, hsp_count, sizeof(HSP*), s_QueryOffsetCompareHSPs);
   i = 0;
   while (i < hsp_count) {
      j = 1;
      while (i+j < hsp_count &&
             hsp_array[i] && hsp_array[i+j] &&
             hsp_array[i]->context == hsp_array[i+j]->context &&
             hsp_array[i]->q_off == hsp_array[i+j]->q_off &&
             hsp_array[i]->s_off == hsp_array[i+j]->s_off) {
         hsp_count--;
         hsp = hsp_array[i+j];
         if (!purge && (hsp->q_end > hsp_array[i]->q_end)) {
             s_CutOffGapEditScript(hsp, hsp_array[i]->q_end,
                                        hsp_array[i]->s_end, TRUE);
         } else {
             //hsp = Blast_HSPFree(hsp);
             HSPDelete(soa, hsp);
             hsp = NULL;
         }
         for (k=i+j; k<hsp_count; k++) {
             hsp_array[k] = hsp_array[k+1];
         }
         hsp_array[hsp_count] = hsp;
      }
      i += j;
   }
   
   qsort(hsp_array, hsp_count, sizeof(HSP*), s_QueryEndCompareHSPs);
   i = 0;
   while (i < hsp_count) {
      j = 1;
      while (i+j < hsp_count &&
             hsp_array[i] && hsp_array[i+j] &&
             hsp_array[i]->context == hsp_array[i+j]->context &&
             hsp_array[i]->q_end == hsp_array[i+j]->q_end &&
             hsp_array[i]->s_end == hsp_array[i+j]->s_end) {
         hsp_count--;
         hsp = hsp_array[i+j];
         if (!purge && (hsp->q_off < hsp_array[i]->q_off)) {
             s_CutOffGapEditScript(hsp, hsp_array[i]->q_off,
                                        hsp_array[i]->s_off, FALSE);
         } else {
             //hsp = Blast_HSPFree(hsp);
             HSPDelete(soa, hsp);
             hsp = NULL;
         }
         for (k=i+j; k<hsp_count; k++) {
             hsp_array[k] = hsp_array[k+1];
         }
         hsp_array[hsp_count] = hsp;
      }
      i += j;
   }

   if (purge) {
      Blast_HSPListPurgeNullHSPs(hsp_list, num_hsps);
   }

   return hsp_count;
}

Int2 Blast_HSPListPurgeNullHSPs(HSP** hsp_list, Int4& num_hsps)
{
    if (hsp_list == NULL || num_hsps == 0) return 0;
    
    Int4 index, index1;
    index1 = 0;
    for (index = 0; index < num_hsps; ++index)
    {
        if (hsp_list[index] != NULL)
        {
            hsp_list[index1] = hsp_list[index];
            ++index1;
        }
    }
    
    for (index = index1; index < num_hsps; ++index)
    {
        hsp_list[index] = NULL;
    }
    
    num_hsps = index1;
    return 0;
}
