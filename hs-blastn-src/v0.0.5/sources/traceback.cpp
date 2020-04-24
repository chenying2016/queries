/* 
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 *
 */

#include "traceback.h"
#include "stat_functions.h"

#include <algorithm>

void 
BlastGetStartForGappedAlignmentNucl (const Uint1* query, const Uint1* subject, HSP* hsp)
{
    /* We will stop when the identity count reaches to this number */
    const Int4 HSP_MAX_IDENT_RUN = 20; 
    const Uint1 *q, *s;
    Int4 index, max_offset, score, max_score, q_start, q_len;
    Int8 s_start;
    Boolean match, prev_match;
    Int4 offset = MIN(hsp->subject_gapped_start - hsp->s_off,
                      hsp->query_gapped_start - hsp->q_off);
    /* first check if the old value is ok */
    q_start = hsp->query_gapped_start;
    s_start = hsp->subject_gapped_start;
    score = -1;
    q = query + q_start;
    s = subject + s_start;
    q_len = hsp->q_end;
    while ((q-query < q_len) && (*q++ == *s++)) {
        score++;
        if (score > HSP_MAX_IDENT_RUN) 
		{
			return;
		}
    }
    q = query + q_start;
    s = subject + s_start;
    while ((q-query >= 0) && (*q-- == *s--)) {
        score++;
        if (score > HSP_MAX_IDENT_RUN) 
		{
			return;
		}
    }
    /* if the old value is not ok, try to find a better point */
    q_start = hsp->query_gapped_start - offset;
    s_start = hsp->subject_gapped_start - offset;
    q_len = MIN(hsp->s_end - s_start, hsp->q_end - q_start);
    q = query + q_start;
    s = subject + s_start;
    max_score = 0;
    max_offset = q_start;
    score = 0;
    match = FALSE;
    prev_match = FALSE; 
    for (index = q_start; index < q_start + q_len; index++) {
        match = (*q++ == *s++);
        if (match != prev_match) {
            prev_match = match;
            if (match) {
                score = 1;
            } else if (score > max_score) {
                max_score = score;
                max_offset = index - score/2;
            }
        } else if (match) {
            ++score;
            if (score > HSP_MAX_IDENT_RUN) {
                max_offset = index - HSP_MAX_IDENT_RUN/2;
                hsp->query_gapped_start = max_offset;
                hsp->subject_gapped_start = max_offset + s_start - q_start;
                return;
            } 
        }
    }
    if (match && score > max_score) {
        max_score = score;
        max_offset = index - score/2;
    }
    if (max_score > 0) {
        hsp->query_gapped_start = max_offset;
        hsp->subject_gapped_start = max_offset + s_start - q_start;
    }
}

/** Maximal subject length after which the offsets are adjusted to a 
 * subsequence.
 */
#define MAX_SUBJECT_OFFSET 2000000000
/** Approximate upper bound on a number of gaps in an HSP, needed to determine
 * the length of the subject subsequence to be retrieved for alignment with
 * traceback. 
 */
#define MAX_TOTAL_GAPS 3000

void
AdjustSubjectRange(Int8* subject_offset_ptr, 
                   Int8* subject_length_ptr,
                   Int4 query_offset, 
                   Int4 query_length, 
                   Int8* start_shift)
{
    Int8 s_offset;
    Int8 subject_length = *subject_length_ptr;
    Int4 max_extension_left, max_extension_right;
    
    /* if subject sequence is not too long, leave everything as is */
    if (subject_length < MAX_SUBJECT_OFFSET)
    {
        *start_shift = 0;
        return;
    }
    
    s_offset = *subject_offset_ptr;
    
    /* maximal extension length is the remaining length in the query,
     plus an estimate of a maximal total number of gaps.*/
    max_extension_left = query_offset + MAX_TOTAL_GAPS;
    max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;
    
    if (s_offset <= max_extension_left) *start_shift = 0;
    else
    {
        *start_shift = s_offset - max_extension_left;
        *subject_offset_ptr = max_extension_left;
    }
    
    *subject_length_ptr = 
            MIN(subject_length, s_offset + max_extension_right) - *start_shift;
}

Int2
s_Blast_HSPGetNumIdentitiesAndPositives(const Uint1* query, 
                                        const Uint1* subject,
                                        const HSP* hsp, 
                                        Int4* num_ident_ptr,
                                        Int4* align_length_ptr,
                                        const BlastScoreBlk* sbp,
                                        Int4* num_pos_ptr)
{
    Int4 i, num_ident, align_length, q_off;
    Int8 s_off;
    Uint1 *q, *s;
    
    q_off = hsp->q_off;
    s_off = hsp->s_off;
    
    if (!subject || !query || !hsp) return -1;
    
    q = (Uint1*)&query[q_off];
    s = (Uint1*)&subject[s_off];
    
    num_ident = 0;
    align_length = 0;
    
    Int4 index;
    GapEditScript* esp = hsp->esp;
    for (index = 0; index < esp->size; ++index)
    {
        align_length += esp->num[index];
        switch (esp->op_type[index])
        {
        case eGapAlignSub:
            for (i = 0; i < esp->num[i]; ++i)
                if (*s++ == *q++)++num_ident;
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
    
    if (align_length_ptr)
        *align_length_ptr = align_length;
    *num_ident_ptr = num_ident;
    
    return 0;
}

Boolean
s_UpdateReevaluatedHSP(HSP* hsp,
                       Boolean gapped,
                       Int4 cutoff_score,
                       Int4 score,
                       const Uint1* query_start,
                       const Uint1* subject_start,
                       const Uint1* best_q_start,
                       const Uint1* best_q_end,
                       const Uint1* best_s_start,
                       const Uint1* best_s_end,
                       int best_start_esp_index,
                       int best_end_esp_index,
                       int best_end_esp_num,
                       SmallObjAllocator& soa)
{
    Boolean delete_hsp = TRUE;
    
    hsp->score = score;
    
    if (hsp->score >= cutoff_score)
    {
        hsp->q_off = best_q_start - query_start;
        hsp->q_end = hsp->q_off + best_q_end - best_q_start;
        hsp->s_off = best_s_start - subject_start;
        hsp->s_end = hsp->s_off + best_s_end - best_s_start;
        
        if (gapped)
        {
            int last_num = hsp->esp->size - 1;
            if (best_end_esp_index != last_num || best_start_esp_index > 0)
            {               
                GapEditScript* esp_temp = GapEditScriptNew(soa, best_end_esp_index - best_start_esp_index + 1);
                GapEditScriptPartialCopy(esp_temp, 0, hsp->esp, best_start_esp_index, best_end_esp_index);
                GapEditScriptDelete(soa, hsp->esp);
                hsp->esp = esp_temp;
            }
            last_num = hsp->esp->size - 1;
            hsp->esp->num[last_num] = best_end_esp_num;
            ASSERT(best_end_esp_num >= 0);
        }
        delete_hsp = FALSE;
    }
    
    return delete_hsp;
}

Boolean
Blast_HSPReevaluateWithAmbiguitiesGappd(HSP* hsp,
                                        const Uint1* q,
                                        Int4 qlen,
                                        const Uint1* s,
                                        Int8 slen,
                                        const BlastHitSavingParameters* hit_params,
                                        const BlastScoringParameters* score_params,
                                        BlastScoreBlk* sbp,
                                        SmallObjAllocator& soa)
{
    Int4 sum, score, gap_open, gap_extend;
    Int4 index;
    Int4 qp, ext;
    Int8 sp;
    
    int best_start_esp_index = 0;
    int best_end_esp_index = 0;
    int current_start_esp_index = 0;
    int best_end_esp_num = 0;
    GapEditScript* esp;
    
    const Uint1* best_q_start;
    const Uint1* best_s_start;
    const Uint1* best_q_end;
    const Uint1* best_s_end;
    
    const Uint1* current_q_start;
    const Uint1* current_s_start;
    
    const Uint1 *query, *subject;
    
    Int4** matrix;
    Int2 factor = 1;
    const Uint1 kResidueMask = 0x0f;
    Int4 cutoff_score = hit_params->cutoffs[hsp->context].cutoff_score;
    
    matrix = sbp->matrix->data;
    
    if (score_params->gap_open == 0 && score_params->gap_extend == 0)
    {
        if (score_params->reward % 2 == 1)
            factor = 2;
        gap_open = 0;
        gap_extend = 
                (score_params->reward - 2 * score_params->penalty) * factor / 2;
    }
    else
    {
        gap_open = score_params->gap_open;
        gap_extend = score_params->gap_extend;
    }
    
    query = q + hsp->q_off;
    subject = s + hsp->s_off;
    score = 0;
    sum = 0;

    /* Point all pointers to the beginning of the alignment. */
    best_q_start = best_q_end = current_q_start = query;
    best_s_start = best_s_end = current_s_start = subject;
    /* There are no previous edit scripts at the beginning. */    
    
    best_end_esp_num = -1;
    esp = hsp->esp;
    if (!esp) return TRUE;
    for (index = 0; index < esp->size; ++index)
    {
        int op_index = 0;
        for (op_index = 0; op_index < esp->num[index]; )
        {
            if (esp->op_type[index] == eGapAlignSub)
            {
                sum += factor * matrix[*query & kResidueMask][*subject];
                query++;
                subject++;
                op_index++;
            }
            else if (esp->op_type[index] == eGapAlignDel)
            {
                sum -= gap_open + gap_extend * esp->num[index];
                subject += esp->num[index];
                op_index += esp->num[index];
            }
            else if (esp->op_type[index] == eGapAlignIns)
            {
                sum -= gap_open + gap_extend * esp->num[index];
                query += esp->num[index];
                op_index += esp->num[index];
            }
            
            if (sum < 0)
            {
                if (op_index < esp->num[index])
                {
                    esp->num[index] -= op_index;
                    current_start_esp_index = index;
                    op_index = 0;
                }
                else
                {
                    current_start_esp_index = index + 1;
                }
                
                sum = 0;
                
                current_q_start = query;
                current_s_start = subject;
                
                if (score < cutoff_score)
                {
                    best_q_start = query;
                    best_s_start = subject;
                    score = 0;
                    
                    best_start_esp_index = current_start_esp_index;
                    best_end_esp_index = current_start_esp_index;
                }
            }
            else if (sum > score)
            {
                score = sum;
                
                best_q_start = current_q_start;
                best_s_start = current_s_start;
                best_q_end = query;
                best_s_end = subject;
                
                best_start_esp_index = current_start_esp_index;
                best_end_esp_index = index;
                best_end_esp_num = op_index;
            }
        }
    }
    
    score /= factor;
    
    if (best_start_esp_index < esp->size && best_end_esp_index < esp->size)
    {
        ASSERT(esp->op_type[best_start_esp_index] == eGapAlignSub);
        ASSERT(esp->op_type[best_end_esp_index] == eGapAlignSub);     
        
        qp = best_q_start - q;
        sp = best_s_start - s;
        ext = 0;
        while (qp > 0 && sp > 0 && (q[--qp] == s[--sp]) && q[qp] < 4) ext++;
        best_q_start -= ext;
        best_s_start -= ext;
        esp->num[best_start_esp_index] += ext;
        if (best_end_esp_index == best_start_esp_index) best_end_esp_num += ext;
        score += ext * score_params->reward;
        
        qp = best_q_end - q;
        sp = best_s_end - s;
        ext = 0;
        while (qp < qlen && sp < slen && q[qp] < 4 && (q[qp++] == s[sp++])) ext++;
        best_q_end += ext;
        best_s_end += ext;
        esp->num[best_end_esp_index] += ext;
        best_end_esp_num += ext;
        score += ext * score_params->reward;
    }
    
    return s_UpdateReevaluatedHSP(hsp, TRUE, cutoff_score, score, q, s, 
                                  best_q_start, best_q_end, 
                                  best_s_start, best_s_end,
                                  best_start_esp_index,
                                  best_end_esp_index,
                                  best_end_esp_num, soa);
}

static Boolean s_HSPTest(HSP* hsp,
                         const BlastHitSavingOptions* hit_options,
                         Int4 align_length)
{
    return ((hsp->num_ident * 100.0 <
            align_length * hit_options->percent_identity) ||
            align_length < hit_options->min_hit_length);
}

Boolean
Blast_HSPTestIdentityAndLength(HSP* hsp,
                               const Uint1* query,
                               const Uint1* subject,
                               const BlastScoringOptions* score_options,
                               const BlastHitSavingOptions* hit_options)
{
    Int4 align_length = 0;
    Boolean delete_hsp = FALSE;
    Int2 status = 0;
    
    ASSERT(hsp && query && subject && score_options && hit_options);
    
    status = hsp->s_Blast_HSPGetNumIdentitesAndPositives(query, subject,
                                                         &hsp->num_ident,
                                                         &align_length,
                                                         NULL);
    
    ASSERT(status == 0);
    
    delete_hsp = s_HSPTest(hsp, hit_options, align_length);
    
    return delete_hsp;
}

bool HSPCmpFunc_score(const HSP& a, const HSP& b)
{
    if (a.score > b.score) return true;
    return false;
}

void Blast_HSPListAdjustOddBlastnScores(HSP** hsp_list,
                                        Int4 num_hsps,
                                        Boolean gapped_calculation,
                                        const BlastScoreBlk* sbp)
{
    int index;
    if (hsp_list == NULL ||
        num_hsps == 0 ||
        gapped_calculation == FALSE ||
        sbp->round_down == FALSE)
        return;
    
    for (index = 0; index < num_hsps; ++index)
        hsp_list[index]->score &= ~1;
    
    qsort(hsp_list, num_hsps, sizeof(HSP*), ScoreCompareHSPs);
}

Int2 Blast_HSPListGetEvalues(const QueryInfo* query_info,
                             HSP** hsp_list,
                             Int4 num_hsps,
                             Boolean gapped_calculation,
                             const BlastScoreBlk* sbp,
                             double gap_decay_rate,
                             double scaling_factor)
{
    HSP* hsp;
    Blast_KarlinBlk** kbp;
    Int4 hsp_cnt;
    Int4 index;
    Int4 kbp_context;
    double gap_decay_divisor = 1.0;
    if (hsp_list == NULL || num_hsps == 0) return 0;
    
    kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    hsp_cnt = num_hsps;
    
    if (gap_decay_rate != 0.)
        gap_decay_divisor = BLAST_GapDecayDivisor(gap_decay_rate, 1);
    
    for (index = 0; index < hsp_cnt; ++index)
    {
        hsp = hsp_list[index];
        ASSERT(hsp != NULL);
        ASSERT(scaling_factor != 0.0);
        ASSERT(sbp->round_down == FALSE || (hsp->score & 1) == 0);
        
        kbp_context = hsp->context;
        kbp[kbp_context]->Lambda /= scaling_factor;
        
        hsp->evalue = BLAST_KarlinStoE_simple(hsp->score, kbp[kbp_context],
                                             sbp->eff_searchsp[kbp_context]);
        
        hsp->evalue /= gap_decay_divisor;
        kbp[kbp_context]->Lambda *= scaling_factor;
    }
    
    return 0;
}

Int2 Blast_HSPListGetBitScores(HSP** hsp_list,
                               Int4 num_hsps,
                               Boolean gapped_calculation,
                               const BlastScoreBlk* sbp)
{
    Blast_KarlinBlk** kbp;
    Int4 index;
    HSP* hsp;
    
    if (hsp_list == NULL || num_hsps == 0) return 1;
    
    kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    
    for (index = 0; index < num_hsps; ++index)
    {
        hsp = hsp_list[index];
        ASSERT(hsp != NULL);
        ASSERT(sbp->round_down == FALSE || (hsp->score & 1) == 0);
        hsp->bit_score = 
                (hsp->score * kbp[hsp->context]->Lambda - kbp[hsp->context]->logK) / NCBIMATH_LN2;
    }
    return 0;
}

Int2 Blast_HSPListReapByEvalues(HSP** hsp_list,
                                Int4& num_hsps,
                                SmallObjAllocator& soa,
                                const BlastHitSavingOptions* hit_options)
{
    HSP* hsp;
    Int4 index = 0, i;
    double cutoff;
    
    if (hsp_list == NULL || num_hsps == 0) return 0;
    
    cutoff = hit_options->expect_value;
    
    for (i = 0; i < num_hsps; ++i)
    {
        hsp = hsp_list[i];
        ASSERT(hsp != NULL);
        
        if (hsp->evalue > cutoff)
        {
            hsp_list[i] = HSPDelete(soa, hsp_list[i]);
        }
        else
        {
            if (i > index)
                hsp_list[index] = hsp_list[i];
            ++index;
        }
    }
    
    num_hsps = index;
    
    return 0;
}

static void 
s_HSPListRescaleScores(HSP** hsp_list,
                       Int4 num_hsps,
                       double scale_factor)
{
   Int4 index;

   for (index = 0; index < num_hsps; ++index) {
      HSP* hsp = hsp_list[index];
      
      /* Remove any scaling of the calculated score */
      hsp->score = 
         (Int4) ((hsp->score+(0.5*scale_factor)) / scale_factor);
   }

   /* Sort HSPs by score again because after the loop above scores that
    * were previously different can become equal, and then the order of HSPs
    * should be determined by the tie-breaking criteria 
    * (e.g.: subject offsets, ...) */
   qsort(hsp_list, num_hsps, sizeof(HSP*), ScoreCompareHSPs);
}

Int2 s_HSPListPostTracebackUpdate(HSP** hsp_list,
                                  Int4& num_hsps,
                                  SmallObjAllocator& soa,
                                  QueryInfo* query_info,
                                  const BlastScoringParameters* score_params,
                                  const BlastHitSavingParameters* hit_params,
                                  BlastScoreBlk* sbp)
{
    BlastScoringOptions* score_options = score_params->options;
    const Boolean kGapped = score_options->gapped_calculation;
    
    double scale_factor = 1.0;
    
    Blast_HSPListAdjustOddBlastnScores(hsp_list, num_hsps, kGapped, sbp);
    
    Blast_HSPListGetEvalues(query_info, hsp_list, num_hsps, kGapped, sbp,
                            0, scale_factor);
    
    Blast_HSPListReapByEvalues(hsp_list, num_hsps, soa, hit_params->options);
    
    s_HSPListRescaleScores(hsp_list, num_hsps, score_params->scale_factor);
    
    Blast_HSPListGetBitScores(hsp_list, num_hsps, kGapped, sbp);
    
    return 0;
}

int GetNumIndels(HSP* hsp)
{
    GapEditScript* esp = hsp->esp;
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


Int2
Blast_TraceBackFromOneContextHSPList(HSP** hsp_list,
                                     Int4& num_hsps,
                                     SmallObjAllocator& soa,
                                     QueryInfo* query_info,
                                     BlastScoreBlk* sbp,
                                     DbInfo* dbinfo,
                                     IntervalTree& itree,
                                     GreedyAligner* gapped_aligner,
                                     const BlastScoringParameters* score_params,
                                     const BlastExtensionParameters* ext_params,
                                     const BlastHitSavingParameters* hit_params)
{
    if (num_hsps == 0) return 0;

    qsort(hsp_list, num_hsps, sizeof(HSP*), ScoreCompareHSPs);
    
    Int4 context = hsp_list[0]->context;
    Int4 i, index;
    const Uint1* q = query_info->GetSequence(context);
    Int4 q_length = query_info->contexts[context].length;
    
    Int8 s_start, s_length, s_shift, s_offset;
    Int4 sid;
    const Uint1* db = (const Uint1*)dbinfo->GetDb(), *ss;
    HSP* hsp;
    
    itree.Reset(0, q_length + 1, 0, dbinfo->GetDbLength() + 1);
    
    for (i = 0; i < num_hsps; ++i)
    {
        hsp = hsp_list[i];         
		sid = hsp->subject_id; 
		s_start = dbinfo->GetSeqOffset(sid);
		hsp->s_off += s_start;
		hsp->s_end += s_start;
		hsp->subject_gapped_start += s_start;
        
        if (itree.IntervalTreeContainsHSP(hsp, hit_params->options->min_diag_separation))
        {
            hsp_list[i] = HSPDelete(soa, hsp_list[i]);
            continue;
        }       
        
        BlastGetStartForGappedAlignmentNucl(q, db, hsp);           

        s_length = dbinfo->GetSeqLength(sid);
        s_offset = hsp->subject_gapped_start - s_start;
        AdjustSubjectRange(&s_offset, &s_length, hsp->query_gapped_start,
                           q_length, &s_shift);
        ss = db + s_start + s_shift;
        
        Int4 q_offset = hsp->query_gapped_start;
        gapped_aligner->GreedyGappedAlignment(q, ss, q_length, s_length,
                                              q_offset,
                                              s_offset, TRUE, 1);
        
        gapped_aligner->PackHSP(*hsp, 0, s_start + s_shift);       
        hsp->subject_id = sid; 
        
        itree.IntervalTreeAddHSP(hsp);
    }
    
    Blast_HSPListPurgeNullHSPs(hsp_list, num_hsps);  
       
    Blast_HSPListPurgeHSPsWithCommonEndpoints(hsp_list,
                                          num_hsps,
                                          FALSE,
                                          soa);

    for (index = 0; index < num_hsps; ++index)
    {
        hsp = hsp_list[index];
        if (!hsp) continue;
        sid = hsp->subject_id;
        s_start = dbinfo->GetSeqOffset(sid);
        hsp->s_off -= s_start;
        hsp->s_end -= s_start;
        hsp->subject_gapped_start -= s_start;         
    }
    
    for (index = 0; index < num_hsps; ++index)
    {
        Boolean delete_hsp = FALSE;
        hsp = hsp_list[index];
        if (!hsp) continue;
        
        ss = db + dbinfo->GetSeqOffset(hsp->subject_id);
        s_length = dbinfo->GetSeqLength(hsp->subject_id);
        
        delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGappd(hsp,
                                                             q,
                                                             q_length,
                                                             ss,
                                                             s_length,
                                                             hit_params,
                                                             score_params,
                                                             sbp,
                                                             soa);
        
        if (!delete_hsp)
            delete_hsp = Blast_HSPTestIdentityAndLength(hsp, 
                                                        q,
                                                        ss,
                                                        score_params->options,
                                                        hit_params->options);
        
        if (delete_hsp)
        {
            hsp_list[index] = HSPDelete(soa, hsp_list[index]);
        }
    }
    
    Blast_HSPListPurgeNullHSPs(hsp_list, num_hsps); 
    
    if (num_hsps > 1)
        qsort(hsp_list, num_hsps, sizeof(HSP*), ScoreCompareHSPs);
    
    itree.Reset(0, q_length + 1, 0, dbinfo->GetDbLength() + 1);
    for (index = 0; index < num_hsps; ++index)
    {
        hsp = hsp_list[index];      
        
        if (itree.IntervalTreeContainsHSP(hsp, hit_params->options->min_diag_separation))
        {
            hsp_list[index] = HSPDelete(soa, hsp_list[index]);
        }
        else
        {
            itree.IntervalTreeAddHSP(hsp);
        }
    }
    
    Blast_HSPListPurgeNullHSPs(hsp_list, num_hsps);    
    
    s_HSPListPostTracebackUpdate(hsp_list,
                                 num_hsps,
                                 soa,
                                 query_info,
                                 score_params,
                                 hit_params,
                                 sbp);  
    return 0;
}
