#include "traceback_stage.h"

#include "../../ncbi_blast/setup/blast_encoding.h"
#include "../../ncbi_blast/setup/hsp2string.h"

/** TRUE if c is between a and b; f between d and e.  Determines if the
 * coordinates are already in an HSP that has been evaluated. 
*/
#define CONTAINED_IN_HSP(a,b,c,d,e,f) \
    (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)

/** Are the two HSPs within a given number of diagonals from each other? */
#define MB_HSP_CLOSE(q1, s1, q2, s2, c) \
    (ABS(((q1)-(s1)) - ((q2)-(s2))) < c)

/** Determine whether an HSP is contained within another HSP.
 *  @param in_hsp The input HSP 
 *  @param in_q_start The start offset of the strand of the query
 *                    sequence containing in_hsp [in]
 *  @param tree_hsp An HSP from the interval tree [in]
 *  @param tree_q_start The start offset of the strand of the query
 *                      sequence containing tree_hsp [in]
 *  @param min_diag_separation Number of diagonals separating 
 *                             nonoverlapping hits (only nonzero 
 *                             for megablast) [in]
 *  @return TRUE if the second HSP envelops the first, FALSE otherwise
 */
static Boolean
s_HSPIsContained(const BlastHSP *hsp1,
                 const BlastHSP *hsp2,
                 Int4 min_diag_separation)
{
    const BlastHSP* tree_hsp = NULL;
    const BlastHSP* in_hsp = NULL;
    if (hsp1->score > hsp2->score) {
        tree_hsp = hsp1;
        in_hsp = hsp2;
    } else {
        tree_hsp = hsp2;
        in_hsp = hsp1;
    }
    /* check if alignments are from different query sequences 
       or query strands */

    if (tree_hsp->context != in_hsp->context) return FALSE;

    //if (in_q_start != tree_q_start)
    //    return FALSE;
       
    if (in_hsp->score <= tree_hsp->score &&
        SIGN(in_hsp->subject.frame) == SIGN(tree_hsp->subject.frame) &&
        CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end, 
                              in_hsp->query.offset,
                              tree_hsp->subject.offset, tree_hsp->subject.end, 
                              in_hsp->subject.offset) &&
        CONTAINED_IN_HSP(tree_hsp->query.offset, tree_hsp->query.end, 
                             in_hsp->query.end,
                             tree_hsp->subject.offset, tree_hsp->subject.end, 
                             in_hsp->subject.end)) {

        if (min_diag_separation == 0)
            return TRUE;

        if (MB_HSP_CLOSE(tree_hsp->query.offset, tree_hsp->subject.offset,
                         in_hsp->query.offset, in_hsp->subject.offset,
                         min_diag_separation) ||
            MB_HSP_CLOSE(tree_hsp->query.end, tree_hsp->subject.end,
                         in_hsp->query.end, in_hsp->subject.end,
                         min_diag_separation)) {
            return TRUE;
        }
    }

    return FALSE;
}

void
purge_contained_hsps(BlastHSPList* hsp_list, const int min_diag_seperation)
{
    for (int i = 0; i < hsp_list->hspcnt; ++i) {
        BlastHSP* hsp1 = hsp_list->hsp_array[i];
        if (!hsp1) continue;
        for (int j = i + 1; j < hsp_list->hspcnt; ++j) {
            BlastHSP* hsp2 = hsp_list->hsp_array[j];
            if (!hsp2) continue;
            if (s_HSPIsContained(hsp1, hsp2, min_diag_seperation)) hsp_list->hsp_array[j] = Blast_HSPFree(hsp2);
        }
    }
    Blast_HSPListPurgeNullHSPs(hsp_list);
}

/** Swaps insertions and deletions in an edit script for RPS BLAST search.
 * This is necessary because query and subject are switched during the
 * traceback alignment, and must be switched back.
 * @param hsp HSP structure to fix. [in] [out]
 */
static void
s_BlastHSPRPSUpdate(BlastHSP *hsp)
{
   GapEditScript *esp = hsp->gap_info;
   Int4 index;

   if (hsp->gap_info == NULL)
      return;

   for (index=0; index<esp->size; index++)
   {
      if (esp->op_type[index] == eGapAlignIns)
          esp->op_type[index] = eGapAlignDel;
      else if (esp->op_type[index] == eGapAlignDel)
          esp->op_type[index] = eGapAlignIns;
   }
}

/** Convert translation frame or strand into a context number suitable for 
 * indexing into the BlastQueryInfo::contexts array
 * @param frame Frame (allowed values: 1,2,3,-1,-2,-3, 0) [in]
 * @param program Type of BLAST program [in]
 * @return context number: 0 or 1 for nucleotide query/subjects, 
 * a value between 0 and 5 (inclusive) for translated query/subjects, and 0 for 
 * protein query/subjects.
 */
NCBI_XBLAST_EXPORT
Int4 BLAST_FrameToContext(Int2 frame, EBlastProgramType program);

/** Remove scaling from scores previously calculated on the hsp_list.
 * @param hsp_list list of HPSs with the score field calculated [in|out]
 * @param scale_factor factor by which scores are scaled, for everything other
 * than RPS-BLAST this should be 1 [in]
 * @todo rename to something which is more intention revealing, merge with
 * function of the same name in blast_kappa.c
 */
static void
s_HSPListRescaleScores(BlastHSPList* hsp_list, double scale_factor)
{
   Int4 index;

   for (index = 0; index < hsp_list->hspcnt; ++index) {
      BlastHSP* hsp = hsp_list->hsp_array[index];

      /* Remove any scaling of the calculated score */
      hsp->score =
         (Int4) ((hsp->score+(0.5*scale_factor)) / scale_factor);
   }

   /* Sort HSPs by score again because after the loop above scores that
    * were previously different can become equal, and then the order of HSPs
    * should be determined by the tie-breaking criteria
    * (e.g.: subject offsets, ...) */
   Blast_HSPListSortByScore(hsp_list);
}

/** Switches back the query and subject in all HSPs in an HSP list; also
 * reassigns contexts to indicate query context, needed to pick correct
 * Karlin block later in the code.
 * @param program Program type: RPS or RPS tblastn [in]
 * @param hsplist List of HSPs [in] [out]
 */
static void
s_BlastHSPListRPSUpdate(EBlastProgramType program, BlastHSPList *hsplist)
{
   Int4 i;
   BlastHSP **hsp;
   BlastSeg tmp;

   /* If this is not an RPS BLAST search, do not do anything. */
   if ( !Blast_ProgramIsRpsBlast(program))
      return;

   hsp = hsplist->hsp_array;
   for (i = 0; i < hsplist->hspcnt; i++) {

      /* switch query and subject (which are already in local coordinates) */
      tmp = hsp[i]->query;
      hsp[i]->query = hsp[i]->subject;
      hsp[i]->subject = tmp;

      /* Change the traceback information to reflect the query and subject
         sequences getting switched */
      s_BlastHSPRPSUpdate(hsp[i]);

      /* If query was nucleotide, set context, because it is needed in order
         to pick correct Karlin block for calculating bit scores. There are
         6 contexts corresponding to each nucleotide query sequence. */
      if (program == eBlastTypeRpsTblastn) {
          hsp[i]->context = BLAST_FrameToContext(hsp[i]->query.frame, program);
      }
   }
   Blast_HSPListSortByScore(hsplist);
}

/** Updates the e-values after the traceback alignment. Also includes relinking
 * of HSPs in case of sum statistics and calculation of bit scores.
 * @param program_number Type of BLAST program [in]
 * @param hsp_list HSPList obtained after a traceback alignment [in] [out]
 * @param query_info Query information structure [in]
 * @param score_params Scoring parameters [in]
 * @param hit_params Hit saving parameters [in]
 * @param sbp Scoring block [in]
 * @param subject_length Length of the subject sequence - needed for linking
 *                       HSPs [in]
 */
static Int2
s_HSPListPostTracebackUpdate(EBlastProgramType program_number,
   BlastHSPList* hsp_list, const BlastQueryInfo* query_info,
   const BlastScoringParameters* score_params,
   const BlastHitSavingParameters* hit_params,
   const BlastScoreBlk* sbp, Int4 subject_length)
{
   BlastScoringOptions* score_options = score_params->options;
   const Boolean kGapped = score_options->gapped_calculation;

   /* Revert query and subject to their traditional meanings.
      This involves switching the offsets around and reversing
      any traceback information */
   s_BlastHSPListRPSUpdate(program_number, hsp_list);

   /* Relink and rereap the HSP list, if needed. */
   if (hit_params->link_hsp_params) {
      //BLAST_LinkHsps(program_number, hsp_list, query_info, subject_length,
      //               sbp, hit_params->link_hsp_params, kGapped);
   } else {
      /* Only use the scaling factor from parameters structure for RPS BLAST,
       * because for other programs either there is no scaling at all, or, in
       * case of composition based statistics, Lambda is scaled as well as
       * scores, and hence scaling factor should not be used for e-value
       * computations.
       */
      double scale_factor =
         (Blast_ProgramIsRpsBlast(program_number) ?
         score_params->scale_factor : 1.0);

      /* For nucleotide search, if match score is = 2, the odd scores
         are rounded down to the nearest even number. */
#if 0
      Blast_HSPListAdjustOddBlastnScores(hsp_list, kGapped, sbp);
#endif

      Blast_HSPListGetEvalues(program_number, query_info, subject_length,
                              hsp_list, kGapped, FALSE, sbp, 0,
                              scale_factor);
   }

   Blast_HSPListReapByEvalue(hsp_list, hit_params->options);

   /* Rescale the scores by scaling factor, if necessary. This rescaling
    * should be done for all programs where scaling factor is not 1.
    */
   s_HSPListRescaleScores(hsp_list, score_params->scale_factor);

   /** Calculate and fill the bit scores. @todo: This is not the only time
    * when they are calculated, s_HSPListRescaleScores also does this in
    * blast_kappa.c.
    */
   Blast_HSPListGetBitScores(hsp_list, kGapped, sbp);

   return 0;
}

void
add_align_string(BlastHSP* hsp, const u8* query, const u8* subject, kstring_t* aligned_string)
{
    const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
    'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 
    'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
    };
    hsp->hsp_info.query_align_offset = ks_size(*aligned_string);
    int qi = hsp->hbn_query.offset;
    const u8* q = query + qi;
    for (int i = 0; i < hsp->gap_info->size; ++i) {
        EGapAlignOpType type = hsp->gap_info->op_type[i];
        int num = hsp->gap_info->num[i];
        if (type == eGapAlignSub || type == eGapAlignIns) {
            for (int p = 0; p < num; ++p, ++q, ++qi) {
                hbn_assert(qi < hsp->hbn_query.seq_size);
                u8 c = *q;
                hbn_assert(c >= 0 && c < BLASTNA_SIZE, 
                    "i = %d, nop = %d, p = %d, num = %d, c = %d, qi = %d, qsize = %d, qid = %d",
                    i, hsp->gap_info->size, p, num, c, qi, hsp->hbn_query.seq_size, hsp->hbn_query.oid);
                int dc = BLASTNA_TO_IUPACNA[c];
                kputc(dc, aligned_string);
            }
        } else if (type == eGapAlignDel) {
            for (int p = 0; p < num; ++p) {
                kputc(GAP_CHAR, aligned_string);
            }
        }
    }
    hbn_assert(qi == hsp->hbn_query.end);

    hsp->hsp_info.subject_align_offset = ks_size(*aligned_string);
    int si = hsp->hbn_subject.offset;
    const u8* s = subject + si;
    for (int i = 0; i < hsp->gap_info->size; ++i) {
        EGapAlignOpType type = hsp->gap_info->op_type[i];
        int num = hsp->gap_info->num[i];
        if (type == eGapAlignSub || type == eGapAlignDel) {
            for (int p = 0; p < num; ++p, ++s, ++si) {
                u8 c = *s;
                int dc = BLASTNA_TO_IUPACNA[c];
                kputc(dc, aligned_string);
            }
        } else if (type == eGapAlignIns) {
            for (int p = 0; p < num; ++p) {
                kputc(GAP_CHAR, aligned_string);
            }
        }
    }
    hbn_assert(si == hsp->hbn_subject.end);
}

static void
update_traceback_hsp_list_info(BlastHSPList* hsp_list, const BLAST_SequenceBlk* query_blk, const BlastQueryInfo* query_info, const u8* subject, Int4** matrix, kstring_t* aligned_string)
{
    for (int i = 0; i < hsp_list->hspcnt; ++i) {
        BlastHSP* hsp = hsp_list->hsp_array[i];
        hsp->hbn_query.offset = hsp->query.offset;
        hsp->hbn_query.end = hsp->query.end;
        hsp->hbn_subject.offset = hsp->subject.offset;
        hsp->hbn_subject.end = hsp->subject.end;
        const u8* query = query_blk->sequence_nomask + query_info->contexts[hsp->context].query_offset;
        const u8* q = hsp->hbn_query.offset + query;
        const u8* s = hsp->hbn_subject.offset + subject;
        int qi = hsp->hbn_query.offset;
        int si = hsp->hbn_subject.offset;
        int align_len = 0;
        int num_ident = 0;
        int num_positives = 0;
        int gaps = 0;
        int gap_opens = 0;
        for (int k = 0; k < hsp->gap_info->size; ++k) {
            EGapAlignOpType type = hsp->gap_info->op_type[k];
            int num = hsp->gap_info->num[k];
            align_len += num;
            switch (type)
            {
            case eGapAlignSub:
                for (int p = 0; p < num; ++p, ++q, ++s, ++qi, ++si) {
                    u8 qc = *q;
                    u8 sc = *s;
                    if (qc == sc) ++num_ident;
                    if (matrix[qc][sc] > 0) ++num_positives;
                }
                break;
            case eGapAlignDel:
                for (int p = 0; p < num; ++p, ++s, ++si) {
                    u8 qc = 15;
                    u8 sc = *s;
                    if (matrix[qc][sc] > 0) ++num_positives;
                }
                ++gap_opens;
                gaps += num;
                break;
            case eGapAlignIns:
                for (int p = 0; p < num; ++p, ++q, ++qi) {
                    u8 qc = *q;
                    u8 sc = 15;
                    if (matrix[qc][sc] > 0) ++num_positives;
                }
                ++gap_opens;
                gaps += num;
            default:
                break;
            }
        }
        hbn_assert(qi == hsp->hbn_query.end);
        hbn_assert(si == hsp->hbn_subject.end);
        hsp->hsp_info.align_len = align_len;
        hsp->hsp_info.num_ident = num_ident;
        hsp->hsp_info.num_positives = num_positives;
        hsp->hsp_info.gap_opens = gap_opens;
        hsp->hsp_info.gaps = gaps;
        hsp->num_ident = num_ident;
        hsp->num_positives = num_positives;
        if (align_len > 0) hsp->hsp_info.perc_identity = 100.0 * num_ident / align_len;
        add_align_string(hsp, query, subject, aligned_string);
    }
}

int
compute_traceback_from_hsplist(EBlastProgramType program_number,
    BlastHSPList* hsp_list,
    const BLAST_SequenceBlk* query_blk,
    const BlastQueryInfo* query_info,
    CSeqDB* subject_blk,
    const BlastScoreBlk* sbp,
    const BlastScoringParameters* score_params,
    const BlastExtensionOptions* ext_options,
    const BlastHitSavingParameters* hit_params,
    kstring_t* aligned_string)
{
    if (!hsp_list->hspcnt) return 0;
    BlastHitSavingOptions* hit_options = hit_params->options;
    BlastScoringOptions* score_options = score_params->options;
    Blast_HSPListSortByScore(hsp_list);
    BlastHSP** hsp_array = hsp_list->hsp_array;
    const int subject_id = hsp_array[0]->hbn_subject.oid;
    const int subject_length = seqdb_seq_size(subject_blk, subject_id);
    const u8* subject = subject_blk->unpacked_seq + seqdb_seq_offset(subject_blk, subject_id);
    Int4 extra_start = Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, FALSE);
    extra_start = 0;
    for (int index = extra_start; index < hsp_list->hspcnt; ++index) {
        Boolean delete_hsp = FALSE;
        BlastHSP* hsp = hsp_array[index];
        if (!hsp) continue;
        const u8* query = query_blk->sequence_nomask + query_info->contexts[hsp->context].query_offset;
        const int query_length = query_info->contexts[hsp->context].query_length;
        delete_hsp = Blast_HSPReevaluateWithAmbiguitiesGapped(hsp, query, query_length, 
                        subject, subject_length, hit_params, score_params, sbp);
        if (!delete_hsp) delete_hsp = Blast_HSPTestIdentityAndLength(program_number, hsp, query, subject, score_options, hit_options);
        if (delete_hsp) hsp_array[index] = Blast_HSPFree(hsp);
    }
    Blast_HSPListPurgeNullHSPs(hsp_list);
    if (program_number == eBlastTypeBlastn) Blast_HSPListPurgeHSPsWithCommonEndpoints(program_number, hsp_list, TRUE);
    Blast_HSPListSortByScore(hsp_list);
    purge_contained_hsps(hsp_list, hit_options->min_diag_separation);
    s_HSPListPostTracebackUpdate(program_number, hsp_list, query_info, score_params, hit_params, sbp, subject_length);
    update_traceback_hsp_list_info(hsp_list, query_blk, query_info, subject, sbp->matrix->data, aligned_string);
    return 0;
}