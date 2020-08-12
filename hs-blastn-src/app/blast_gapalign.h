#ifndef ALGO_BLAST_CORE__BLAST_GAPALIGN__H
#define ALGO_BLAST_CORE__BLAST_GAPALIGN__H

#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/setup/gapinfo.h"
#include "../../ncbi_blast/setup/blast_parameters.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Auxiliary structure for dynamic programming gapped extension */
typedef struct {
  Int4 best;            /**< score of best path that ends in a match
                             at this position */
  Int4 best_gap;        /**< score of best path that ends in a gap
                             at this position */
} BlastGapDP;

/** Structure supporting the gapped alignment */
typedef struct BlastGapAlignStruct {
   Boolean positionBased; /**< Is this PSI-BLAST? */
   GapStateArrayStruct* state_struct; /**< Structure to keep extension 
                                                state information */
   GapEditScript* edit_script; /**< The traceback (gap) information */
   GapPrelimEditBlock *fwd_prelim_tback; /**< traceback from right extensions */
   GapPrelimEditBlock *rev_prelim_tback; /**< traceback from left extensions */
   ///SGreedyAlignMem* greedy_align_mem;/**< Preallocated memory for the greedy 
   ///                                      gapped extension */
   BlastGapDP* dp_mem; /**< scratch structures for dynamic programming */
   Int4 dp_mem_alloc;  /**< current number of structures allocated */
   BlastScoreBlk* sbp; /**< Pointer to the scoring information block */
   Int4 gap_x_dropoff; /**< X-dropoff parameter to use */
   Int4 max_mismatches;  /**< Max number of mismatches for jumper */
   Int4 mismatch_window; /**< Window sie for mismatches for jumper */
   Int4 query_start; /**< query start offset of current alignment */
   Int4 query_stop; /**< query end offseet of current alignment */
   Int4 subject_start;  /**< subject start offset current alignment */
   Int4 subject_stop; /**< subject end offset of current alignment */
   Int4 greedy_query_seed_start;  /**< for greedy alignments, the query 
                                       offset of the gapped start point */
   Int4 greedy_subject_seed_start;  /**< for greedy alignments, the subject
                                         offset of the gapped start point */
   Int4 score;   /**< Return value: alignment score */
   
   int score_array[16];
   int* score_matrix[4];

   ///JumperGapAlign* jumper;   /**< data for jumper alignment */
} BlastGapAlignStruct;

/** Initializes the BlastGapAlignStruct structure 
 * @param score_params Parameters related to scoring alignments [in]
 * @param ext_params parameters related to gapped extension [in]
 * @param max_subject_length Maximum length of any subject sequence (needed 
 *        for greedy extension allocation only) [in]
 * @param sbp The scoring information block [in]
 * @param gap_align_ptr The BlastGapAlignStruct structure [out]
*/
NCBI_XBLAST_EXPORT
Int2
BLAST_GapAlignStructNew(const BlastScoringParameters* score_params, 
   const BlastExtensionParameters* ext_params, 
   Uint4 max_subject_length, BlastScoreBlk* sbp, 
   BlastGapAlignStruct** gap_align_ptr);

/** Deallocates memory in the BlastGapAlignStruct structure */
NCBI_XBLAST_EXPORT
BlastGapAlignStruct* 
BLAST_GapAlignStructFree(BlastGapAlignStruct* gap_align);

/** Perform a gapped alignment with traceback
 * @param program Type of BLAST program [in]
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param gap_align The gapped alignment structure [in] [out]
 * @param score_params Scoring parameters [in]
 * @param q_start Offset in query where to start alignment [in]
 * @param s_start Offset in subject where to start alignment [in]
 * @param query_length Maximal allowed extension in query [in]
 * @param subject_length Maximal allowed extension in subject [in]
 * @param fence_hit True is returned here if overrun is detected. [in]
 */
NCBI_XBLAST_EXPORT
Int2 BLAST_GappedAlignmentWithTraceback(EBlastProgramType program, 
        const Uint1* query, const Uint1* subject, 
        BlastGapAlignStruct* gap_align, 
        const BlastScoringParameters* score_params,
        Int4 q_start, Int4 s_start, Int4 query_length, Int4 subject_length,
        Boolean * fence_hit);

#ifdef __cplusplus
}
#endif

#endif // ALGO_BLAST_CORE__BLAST_GAPALIGN__H