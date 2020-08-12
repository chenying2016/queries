#include "blast_gapalign.h"

/* Documented in blast_gapalign.h */
BlastGapAlignStruct*
BLAST_GapAlignStructFree(BlastGapAlignStruct* gap_align)
{
   if (!gap_align)
      return NULL;

   GapEditScriptDelete(gap_align->edit_script);
   GapPrelimEditBlockFree(gap_align->fwd_prelim_tback);
   GapPrelimEditBlockFree(gap_align->rev_prelim_tback);
   //if (gap_align->greedy_align_mem)
   //   s_BlastGreedyAlignsFree(gap_align->greedy_align_mem);
   GapStateFree(gap_align->state_struct);
   sfree(gap_align->dp_mem);
   //JumperGapAlignFree(gap_align->jumper);
   sfree(gap_align);
   return NULL;
}

void ksw_gen_simple_mat(int *score_array, int** score_matrix, int reward, int penalty)
{
	int i, j;
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			score_array[i * 4 + j] = i == j? reward : penalty;
	}
    
    for (int i = 0; i < 4; ++i) score_matrix[i] = score_array + i * 4;
}

/* Documented in blast_gapalign.h */
Int2
BLAST_GapAlignStructNew(const BlastScoringParameters* score_params,
   const BlastExtensionParameters* ext_params,
   Uint4 max_subject_length,
   BlastScoreBlk* sbp, BlastGapAlignStruct** gap_align_ptr)
{
   Int2 status = 0;
   BlastGapAlignStruct* gap_align;

   gap_align = (BlastGapAlignStruct*) calloc(1, sizeof(BlastGapAlignStruct));

   *gap_align_ptr = gap_align;

   gap_align->sbp = sbp;

   gap_align->gap_x_dropoff = ext_params->gap_x_dropoff;
   gap_align->max_mismatches = ext_params->options->max_mismatches;
   gap_align->mismatch_window = ext_params->options->mismatch_window;

   /* allocate memory either for dynamic programming or jumper */
   if (ext_params->options->ePrelimGapExt != eJumperWithTraceback) {
       if (ext_params->options->ePrelimGapExt == eDynProgScoreOnly) {
	   /* allocate structures for ordinary dynamic programming */
	   gap_align->dp_mem_alloc = 1000;
	   gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
						    sizeof(BlastGapDP));
	   if (!gap_align->dp_mem)
	       gap_align = BLAST_GapAlignStructFree(gap_align);
       }
       else {
#if 0
	   /* allocate structures for greedy dynamic programming */
	   max_subject_length = MIN(max_subject_length, MAX_DBSEQ_LEN);
	   max_subject_length = MIN(GREEDY_MAX_COST,
			    max_subject_length / GREEDY_MAX_COST_FRACTION + 1);
	   gap_align->greedy_align_mem =
	       s_BlastGreedyAlignMemAlloc(score_params, ext_params,
					  max_subject_length, 0);
	   if (!gap_align->greedy_align_mem)
	       gap_align = BLAST_GapAlignStructFree(gap_align);
#endif
       }
   }
   else {
#if 0
       gap_align->jumper = JumperGapAlignNew(200);
       if (ext_params->gap_x_dropoff == 0) {
           gap_align->gap_x_dropoff = 3 * MAX(-score_params->penalty,
                                              score_params->gap_open +
                                              score_params->gap_extend);
       }
#endif
   }

   if (!gap_align)
      return -1;

   gap_align->positionBased = FALSE;

   gap_align->fwd_prelim_tback = GapPrelimEditBlockNew();
   gap_align->rev_prelim_tback = GapPrelimEditBlockNew();

   ksw_gen_simple_mat(gap_align->score_array, gap_align->score_matrix, score_params->reward, score_params->penalty);

   return status;
}

/** Lower bound for scores. Divide by two to prevent underflows. */
#define MININT INT4_MIN/2

/** Minimal size of a chunk for state array allocation. */
#define	CHUNKSIZE	2097152

/** Retrieve the state structure corresponding to a given length
 * @param head Pointer to the first element of the state structures
 *        array [in]
 * @param length The length for which the state structure has to be
 *        found [in]
 * @return The found or created state structure
 */
static GapStateArrayStruct*
s_GapGetState(GapStateArrayStruct** head, Int4 length)

{
   GapStateArrayStruct*	retval,* var,* last;
   Int4	chunksize = MAX(CHUNKSIZE, length + length/3);

   length += length/3;	/* Add on about 30% so the end will get reused. */
   retval = NULL;
   if (*head == NULL) {
      retval = (GapStateArrayStruct*)
         malloc(sizeof(GapStateArrayStruct));
      retval->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
      retval->length = chunksize;
      retval->used = 0;
      retval->next = NULL;
      *head = retval;
   } else {
      var = *head;
      last = *head;
      while (var) {
         if (length < (var->length - var->used)) {
            retval = var;
            break;
         } else if (var->used == 0) {
            /* If it's empty and too small, replace. */
            sfree(var->state_array);
            var->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
            var->length = chunksize;
            retval = var;
            break;
         }
         last = var;
         var = var->next;
      }

      if (var == NULL)
      {
         retval = (GapStateArrayStruct*) malloc(sizeof(GapStateArrayStruct));
         retval->state_array = (Uint1*) malloc(chunksize*sizeof(Uint1));
         retval->length = chunksize;
         retval->used = 0;
         retval->next = NULL;
         last->next = retval;
      }
   }

#ifdef ERR_POST_EX_DEFINED
   if (retval->state_array == NULL)
      ErrPostEx(SEV_ERROR, 0, 0, "state array is NULL");
#endif

   return retval;

}

/** Remove a state from a state structure */
static Boolean
s_GapPurgeState(GapStateArrayStruct* state_struct)
{
   while (state_struct)
   {
      /*
	memset(state_struct->state_array, 0, state_struct->used);
      */
      state_struct->used = 0;
      state_struct = state_struct->next;
   }

   return TRUE;
}

/** Values for the editing script operations in traceback */
enum {
    SCRIPT_SUB           = eGapAlignSub,     /**< Substitution */
    SCRIPT_GAP_IN_A      = eGapAlignDel,     /**< Deletion */
    SCRIPT_GAP_IN_B      = eGapAlignIns,     /**< Insertion */
    SCRIPT_OP_MASK       = 0x07, /**< Mask for edit script operations */

    SCRIPT_EXTEND_GAP_A  = 0x10, /**< continue a gap in A */
    SCRIPT_EXTEND_GAP_B  = 0x40  /**< continue a gap in B */
};

Int4
ALIGN_EX(const Uint1* A, const Uint1* B, Int4 M, Int4 N, Int4* a_offset,
	Int4* b_offset, GapPrelimEditBlock *edit_block,
        BlastGapAlignStruct* gap_align,
        const BlastScoringParameters* score_params, Int4 query_offset,
        Boolean reversed, Boolean reverse_sequence,
        Boolean * fence_hit)
{
    /* See Blast_SemiGappedAlign for more general comments on
       what this code is doing; comments in this function
       only apply to the traceback computations */

    Int4 i;
    Int4 a_index;
    Int4 b_index, b_size, first_b_index, last_b_index, b_increment;
    const Uint1* b_ptr;

    BlastGapDP* score_array;

    Int4 gap_open;
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 x_dropoff;
    Int4 best_score;

    Int4** matrix = NULL;
    Int4** pssm = NULL;
    Int4* matrix_row = NULL;

    Int4 score;
    Int4 score_gap_row;
    Int4 score_gap_col;
    Int4 next_score;

    GapStateArrayStruct* state_struct;
    Uint1* edit_script_row;
    Uint1** edit_script;
    Int4 *edit_start_offset;
    Int4 edit_script_num_rows;
    Int4 orig_b_index;
    Uint1 script, next_script, script_row, script_col;
    Int4 num_extra_cells;

    matrix = gap_align->score_matrix;
    if (gap_align->positionBased) {
        pssm = gap_align->sbp->psi_matrix->pssm->data;
    }

    *a_offset = 0;
    *b_offset = 0;
    gap_open = score_params->gap_open;
    gap_extend = score_params->gap_extend;
    gap_open_extend = gap_open + gap_extend;
    x_dropoff = gap_align->gap_x_dropoff;

    if (x_dropoff < gap_open_extend)
        x_dropoff = gap_open_extend;

    if(N <= 0 || M <= 0)
        return 0;

    /* Initialize traceback information. edit_script[] is
       a 2-D array which is filled in row by row as the
       dynamic programming takes place */

    s_GapPurgeState(gap_align->state_struct);

    /* Conceptually, traceback requires a byte array of size
       MxN. To save on memory, each row of the array is allocated
       separately. edit_script[i] points to the storage reserved
       for row i, and edit_start_offset[i] gives the offset into
       the B sequence corresponding to element 0 of edit_script[i].

       Also make the number of edit script rows grow dynamically */

    edit_script_num_rows = 100;
    edit_script = (Uint1**) malloc(sizeof(Uint1*) * edit_script_num_rows);
    edit_start_offset = (Int4*) malloc(sizeof(Int4) * edit_script_num_rows);

    /* allocate storage for the first row of the traceback
       array. Because row elements correspond to gaps in A,
       the alignment can only go x_dropoff/gap_extend positions
       at most before failing the X dropoff criterion */

    if (gap_extend > 0)
        num_extra_cells = x_dropoff / gap_extend + 3;
    else
        num_extra_cells = N + 3;

    if (num_extra_cells > gap_align->dp_mem_alloc) {
        gap_align->dp_mem_alloc = MAX(num_extra_cells + 100,
                                      2 * gap_align->dp_mem_alloc);
        sfree(gap_align->dp_mem);
        gap_align->dp_mem = (BlastGapDP *)malloc(gap_align->dp_mem_alloc *
                                                  sizeof(BlastGapDP));
    }

    state_struct = s_GapGetState(&gap_align->state_struct, num_extra_cells);

    edit_script[0] = state_struct->state_array;
    edit_start_offset[0] = 0;
    edit_script_row = state_struct->state_array;

    score = -gap_open_extend;
    score_array = gap_align->dp_mem;
    score_array[0].best = 0;
    score_array[0].best_gap = -gap_open_extend;

    for (i = 1; i <= N; i++) {
        if (score < -x_dropoff)
            break;

        score_array[i].best = score;
        score_array[i].best_gap = score - gap_open_extend;
        score -= gap_extend;
        edit_script_row[i] = SCRIPT_GAP_IN_A;
    }
    state_struct->used = i + 1;

    b_size = MIN(i, N);
    best_score = 0;
    first_b_index = 0;
    if (reverse_sequence)
        b_increment = -1;
    else
        b_increment = 1;

    for (a_index = 1; a_index <= M; a_index++) {

        /* Set up the next row of the edit script; this involves
           allocating memory for the row, then pointing to it.
           It is not necessary to allocate space for offsets less
           than first_b_index (since the inner loop will never
           look at them);

           It is unknown at this point how far to the right the
           current traceback row will extend; all that's known for
           sure is that the previous row fails the X-dropoff test
           after b_size cells, and that the current row can go at
           most num_extra_cells beyond that before failing the test */

        if (gap_extend > 0)
            state_struct = s_GapGetState(&gap_align->state_struct,
                           b_size - first_b_index + num_extra_cells);
        else
            state_struct = s_GapGetState(&gap_align->state_struct,
                                        N + 3 - first_b_index);

        if (a_index == edit_script_num_rows) {
            edit_script_num_rows = edit_script_num_rows * 2;
            edit_script = (Uint1 **)realloc(edit_script,
                                            edit_script_num_rows *
                                            sizeof(Uint1 *));
            edit_start_offset = (Int4 *)realloc(edit_start_offset,
                                                edit_script_num_rows *
                                                sizeof(Int4));
        }

        edit_script[a_index] = state_struct->state_array +
                                state_struct->used + 1;
        edit_start_offset[a_index] = first_b_index;

        /* the inner loop assumes the current traceback
           row begins at offset zero of B */

        edit_script_row = edit_script[a_index] - first_b_index;
        orig_b_index = first_b_index;

        if (!(gap_align->positionBased)) {
            if(reverse_sequence)
                matrix_row = matrix[ A[ M - a_index ] ];
            else
                matrix_row = matrix[ A[ a_index ] ];
        }
        else {
            if(reversed || reverse_sequence)
                matrix_row = pssm[M - a_index];
            else
                matrix_row = pssm[a_index + query_offset];
        }

        if(reverse_sequence)
            b_ptr = &B[N - first_b_index];
        else
            b_ptr = &B[first_b_index];

        score = MININT;
        score_gap_row = MININT;
        last_b_index = first_b_index;

        for (b_index = first_b_index; b_index < b_size; b_index++) {
            int matrix_index = 0;

            b_ptr += b_increment;
            score_gap_col = score_array[b_index].best_gap;

            matrix_index = *b_ptr;

            if (matrix_index == FENCE_SENTRY) {
                if (fence_hit) {
                    *fence_hit = 1;
                }
                break;
            }

            next_score = score_array[b_index].best + matrix_row[ *b_ptr ];

            /* script, script_row and script_col contain the
               actions specified by the dynamic programming.
               when the inner loop has finished, 'script' con-
               tains all of the actions to perform, and is
               written to edit_script[a_index][b_index]. Otherwise,
               this inner loop is exactly the same as the one
               in Blast_SemiGappedAlign() */

            script = SCRIPT_SUB;
            script_col = SCRIPT_EXTEND_GAP_B;
            script_row = SCRIPT_EXTEND_GAP_A;

            if (score < score_gap_col) {
                script = SCRIPT_GAP_IN_B;
                score = score_gap_col;
            }
            if (score < score_gap_row) {
                script = SCRIPT_GAP_IN_A;
                score = score_gap_row;
            }

            if (best_score - score > x_dropoff) {

                if (first_b_index == b_index)
                    first_b_index++;
                else
                    score_array[b_index].best = MININT;
            }
            else {
                last_b_index = b_index;
                if (score > best_score) {
                    best_score = score;
                    *a_offset = a_index;
                    *b_offset = b_index;
                }

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                if (score_gap_col < (score - gap_open_extend)) {
                    score_array[b_index].best_gap = score - gap_open_extend;
                }
                else {
                    score_array[b_index].best_gap = score_gap_col;
                    script += script_col;
                }

                if (score_gap_row < (score - gap_open_extend))
                    score_gap_row = score - gap_open_extend;
                else
                    script += script_row;

                score_array[b_index].best = score;
            }

            score = next_score;
            edit_script_row[b_index] = script;
        }

        if (first_b_index == b_size || (fence_hit && *fence_hit))
            break;

        if (last_b_index + num_extra_cells + 3 >= gap_align->dp_mem_alloc) {

            gap_align->dp_mem_alloc = MAX(last_b_index + num_extra_cells + 100,
                                          2 * gap_align->dp_mem_alloc);
            score_array = (BlastGapDP *)realloc(score_array,
                                               gap_align->dp_mem_alloc *
                                               sizeof(BlastGapDP));
            gap_align->dp_mem = score_array;
        }


        if (last_b_index < b_size - 1) {
            b_size = last_b_index + 1;
        }
        else {
            while (score_gap_row >= (best_score - x_dropoff) && b_size < N) {

                score_array[b_size].best = score_gap_row;
                score_array[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                edit_script_row[b_size] = SCRIPT_GAP_IN_A;
                b_size++;
            }
        }

        /* update the memory allocator to reflect the exact number
           of traceback cells this row needed */

        state_struct->used += MAX(b_index, b_size) - orig_b_index + 1;

        if (b_size < N) {
            score_array[b_size].best = MININT;
            score_array[b_size].best_gap = MININT;
            b_size++;
        }
    }

    /* Pick the optimal path through the now complete
       edit_script[][]. This is equivalent to flattening
       the 2-D array into a 1-D list of actions. */

    a_index = *a_offset;
    b_index = *b_offset;
    script = SCRIPT_SUB;

    if (fence_hit && *fence_hit)
        goto done;

    while (a_index > 0 || b_index > 0) {
        /* Retrieve the next action to perform. Rows of
           the traceback array do not necessarily start
           at offset zero of B, so a correction is needed
           to point to the correct position */

        next_script =
            edit_script[a_index][b_index - edit_start_offset[a_index]];

        switch(script) {
        case SCRIPT_GAP_IN_A:
            script = next_script & SCRIPT_OP_MASK;
            if (next_script & SCRIPT_EXTEND_GAP_A)
                script = SCRIPT_GAP_IN_A;
            break;

        case SCRIPT_GAP_IN_B:
            script = next_script & SCRIPT_OP_MASK;
            if (next_script & SCRIPT_EXTEND_GAP_B)
                script = SCRIPT_GAP_IN_B;
            break;

        default:
            script = next_script & SCRIPT_OP_MASK;
            break;
        }

        if (script == SCRIPT_GAP_IN_A) {
            b_index--;
        }
        else if (script == SCRIPT_GAP_IN_B) {
            a_index--;
        }
        else {
            a_index--;
            b_index--;
        }
        GapPrelimEditBlockAdd(edit_block, (EGapAlignOpType)script, 1);
    }

 done:
    sfree(edit_start_offset);
    sfree(edit_script);
    return best_score;
}

/** Convert the initial list of traceback actions from a non-OOF
 *  gapped alignment into a blast edit script. Note that this routine
 *  assumes the input edit blocks have not been reversed or rearranged
 *  by calling code
 *  @param rev_prelim_tback Traceback from extension to the left [in]
 *  @param fwd_prelim_tback Traceback from extension to the right [in]
 *  @return Pointer to the resulting edit script, or NULL if there
 *          are no traceback actions specified
 */
GapEditScript*
Blast_PrelimEditBlockToGapEditScript (GapPrelimEditBlock* rev_prelim_tback,
                                      GapPrelimEditBlock* fwd_prelim_tback)
{
   Boolean merge_ops = FALSE;
   GapEditScript* esp;
   GapPrelimEditScript *op;
   Int4 i;
   Int4 index=0;
   Int4 size = 0;

   if (rev_prelim_tback == NULL || fwd_prelim_tback == NULL)
      return NULL;

   /* The fwd_prelim_tback script will get reversed here as the traceback started from the highest scoring point
     and worked backwards. The rev_prelim_tback script does NOT get reversed.  Since it was reversed when the
     traceback was produced it's already "forward" */

   if (fwd_prelim_tback->num_ops > 0 && rev_prelim_tback->num_ops > 0 &&
       fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].op_type ==
         rev_prelim_tback->edit_ops[(rev_prelim_tback->num_ops)-1].op_type)
     merge_ops = TRUE;

   size = fwd_prelim_tback->num_ops+rev_prelim_tback->num_ops;
   if (merge_ops)
     size--;

   esp = GapEditScriptNew(size);

   index = 0;
   for (i=0; i < rev_prelim_tback->num_ops; i++) {
      op = rev_prelim_tback->edit_ops + i;
      esp->op_type[index] = op->op_type;
      esp->num[index] = op->num;
      index++;
   }

   if (fwd_prelim_tback->num_ops == 0)
      return esp;

   if (merge_ops)
       esp->num[index-1] += fwd_prelim_tback->edit_ops[(fwd_prelim_tback->num_ops)-1].num;

   /* If we merge, then we skip the first one. */
   if (merge_ops)
      i = fwd_prelim_tback->num_ops - 2;
   else
      i = fwd_prelim_tback->num_ops - 1;

   for (; i >= 0; i--) {
      op = fwd_prelim_tback->edit_ops + i;
      esp->op_type[index] = op->op_type;
      esp->num[index] = op->num;
      index++;
   }

   return esp;
}

Int2 BLAST_GappedAlignmentWithTraceback(EBlastProgramType program,
        const Uint1* query, const Uint1* subject, BlastGapAlignStruct* gap_align,
        const BlastScoringParameters* score_params,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length,
        Boolean * fence_hit)
{
    Boolean found_end;
    Int4 score_right, score_left, private_q_length, private_s_length;
    Int4 q_length, s_length;
    Boolean is_ooframe = score_params->options->is_ooframe;
    Int2 status = 0;
    Boolean switch_seq = FALSE;
    GapPrelimEditBlock *fwd_prelim_tback;
    GapPrelimEditBlock *rev_prelim_tback;

    if (gap_align == NULL)
        return -1;

    fwd_prelim_tback = gap_align->fwd_prelim_tback;
    rev_prelim_tback = gap_align->rev_prelim_tback;
    GapPrelimEditBlockReset(fwd_prelim_tback);
    GapPrelimEditBlockReset(rev_prelim_tback);

    found_end = FALSE;

    q_length = query_length;
    s_length = subject_length;
    if (is_ooframe) {
       /* The mixed frame sequence is shifted to the 3rd position, so its
          maximal available length for extension is less by 2 than its
          total length. */
       switch_seq = (program == eBlastTypeBlastx);
       if (switch_seq) {
          q_length -= CODON_LENGTH - 1;
       } else {
          s_length -= CODON_LENGTH - 1;
       }
    }

    score_left = 0;

    if(is_ooframe) {
#if 0
       /* NB: Left extension does not include starting point corresponding
          to the offset pair; the right extension does. */
       score_left =
          s_OutOfFrameSemiGappedAlignWrap(query+q_start, subject+s_start,
             q_start, s_start, &private_q_length, &private_s_length,
             FALSE, rev_prelim_tback, gap_align, score_params, q_start,
             TRUE, switch_seq);
       gap_align->query_start = q_start - private_q_length;
       gap_align->subject_start = s_start - private_s_length;
#endif
    } else {
       /* NB: The left extension includes the starting point
          [q_start,s_start]; the right extension does not. */
       // <AG> score_left =
       //    Blast_SemiGappedAlign(query, subject, q_start+1, s_start+1,
       //       &private_q_length, &private_s_length, FALSE, rev_prelim_tback,
       //       gap_align, score_params, q_start, FALSE, TRUE,
       //       fence_hit);

        score_left = ALIGN_EX(query, subject, q_start+1, s_start+1,  &private_q_length, &private_s_length, rev_prelim_tback,
                              gap_align,
                              score_params, q_start, FALSE /*reversed*/, TRUE /*reverse_sequence*/,
             fence_hit);

       gap_align->query_start = q_start - private_q_length + 1;
       gap_align->subject_start = s_start - private_s_length + 1;
    }

    score_right = 0;

    if ((! (fence_hit && *fence_hit)) &&
        (q_start < q_length) &&
        (s_start < s_length)) {

       found_end = TRUE;
       if(is_ooframe) {
#if 0
          score_right =
             s_OutOfFrameSemiGappedAlignWrap(query+q_start-1, subject+s_start-1,
                q_length-q_start, s_length-s_start,
                &private_q_length, &private_s_length, FALSE, fwd_prelim_tback,
                gap_align, score_params, q_start, FALSE, switch_seq);
#endif
        } else {
            // <AG> score_right =
            //    Blast_SemiGappedAlign(query+q_start, subject+s_start,
            //       q_length-q_start-1, s_length-s_start-1, &private_q_length,
            //       &private_s_length, FALSE, fwd_prelim_tback, gap_align,
            //       score_params, q_start, FALSE, FALSE,
            //       fence_hit);
            score_right =
               ALIGN_EX(query+q_start, subject+s_start,
                  q_length-q_start-1, s_length-s_start-1, &private_q_length,
                  &private_s_length, fwd_prelim_tback, gap_align,
                  score_params, q_start, FALSE, FALSE,
                  fence_hit);
        }

        gap_align->query_stop = q_start + private_q_length + 1;
        gap_align->subject_stop = s_start + private_s_length + 1;
    }

    if (found_end == FALSE) {
        gap_align->query_stop = q_start - 1;
        gap_align->subject_stop = s_start - 1;
    }

    if(is_ooframe) {
#if 0
        Int4 nucl_align_length;
        if (program == eBlastTypeBlastx) {
            nucl_align_length = gap_align->query_stop -
                                gap_align->query_start;
        }
        else {
            nucl_align_length = gap_align->subject_stop -
                                gap_align->subject_start;
        }
        status = s_BlastOOFTracebackToGapEditScript(rev_prelim_tback,
                       fwd_prelim_tback, nucl_align_length,
                       &gap_align->edit_script);
#endif
    } else {
        Int4 i;
        GapEditScript *esp = Blast_PrelimEditBlockToGapEditScript(
                                                        rev_prelim_tback,
                                                        fwd_prelim_tback);

        /* rarely (typically when the scoring system changes between
           the score-only and traceback stages, as happens with
           composition-based statistics) it is possible to compute an
           optimal alignment with a leading or trailing gap. Prune
           these unneeded gaps here and update the score and
           alignment boundaries */

        gap_align->edit_script = esp;
        if (esp) {
            while (esp->size && esp->op_type[0] != eGapAlignSub) {
                score_left += score_params->gap_open +
                             esp->num[0] * score_params->gap_extend;

                if (esp->op_type[0] == eGapAlignDel)
                    gap_align->subject_start += esp->num[0];
                else
                    gap_align->query_start += esp->num[0];

                for (i = 1; i < esp->size; i++) {
                    esp->op_type[i-1] = esp->op_type[i];
                    esp->num[i-1] = esp->num[i];
                }
                esp->size--;
            }
            i = esp->size;
            while (esp->size && esp->op_type[i-1] != eGapAlignSub) {
                score_right += score_params->gap_open +
                             esp->num[i-1] * score_params->gap_extend;

                if (esp->op_type[i-1] == eGapAlignDel)
                    gap_align->subject_stop -= esp->num[i-1];
                else
                    gap_align->query_stop -= esp->num[i-1];

                esp->size--;
                i--;
                ASSERT(esp->size == i);
            }
        }
    }

    gap_align->score = score_right + score_left;
    return status;
}