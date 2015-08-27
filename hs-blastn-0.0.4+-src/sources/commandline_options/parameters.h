#ifndef BLAST_PARAMETERS_H
#define	BLAST_PARAMETERS_H

#include "options.h"
#include "stat.h"
#include "sequence.h"

/** Computed values used as parameters for gapped alignments.
 *  Note that currently the value of the X-dropoff parameter
 *  is fixed for all search contexts
 */
struct BlastExtensionParameters {
   BlastExtensionOptions* options; /**< The original (unparsed) options. */
   Int4 gap_x_dropoff; /**< X-dropoff value for gapped extension (raw) */
   Int4 gap_x_dropoff_final;/**< X-dropoff value for the final gapped
                               extension (raw) */

   BlastExtensionParameters(BlastExtensionOptions* ex_options,
                            BlastScoreBlk* sbp,
                            QueryInfo* query_info);
   ~BlastExtensionParameters();
};

/** All the gapped cutoff values that can change
 *  from context to context
 */
struct BlastGappedCutoffs {
   Int4 cutoff_score; /**< Raw cutoff score corresponding to the e-value
                         provided by the user if no sum stats, the lowest score
                         to attempt linking on if sum stats are used. */
   Int4 cutoff_score_max; /**< Raw cutoff score corresponding to the e-value
                         provided by user, cutoff_score must be <= this. */
};

/** Parameter block that contains a pointer to BlastHitSavingOptions
 * and the values derived from it.
 */
struct BlastHitSavingParameters {
   BlastHitSavingOptions* options; /**< The original (unparsed) options. */
   Int4 cutoff_score_min; /**< smallest cutoff score across all contexts */
   BlastGappedCutoffs *cutoffs; /**< per-context gapped cutoff information */
   /*BlastLinkHSPParameters* link_hsp_params;*/ /**< Parameters for linking HSPs
                                               with sum statistics; linking
                                               is not done if NULL. */
   Boolean restricted_align; /**< TRUE if approximate score-only gapped
                                  alignment is used */
   Boolean do_sum_stats;  /**< TRUE if sum stats will be used.  Can override the
                                  do_sum_stats Boolean in the options if criteria for
                                  doing sum stats are not met.  */
   Int4 mask_level; /**< Only keep the highest scoring HSP when more than
                          one HSP overlaps the same region of the query by
                          more than or equal to mask_level %. -RMH- */
   Int4* low_score;             /**< lowest ungapped score that can trigger a
				gapped alignment if the histlist is already full.
                                One value for each query. */

    BlastHitSavingParameters(BlastHitSavingOptions* hs_options,
                             const BlastScoreBlk* sbp,
                             const QueryInfo* query_info,
                             Int8 avg_subj_length);
    Int2 BlastHitSavingParametersUpdate(const BlastScoreBlk* sbp,
                                        const QueryInfo* query_info,
                                        Int8 avg_subj_length);
    ~BlastHitSavingParameters();
};

/** All the ungapped cutoff values that can change
 *  from context to context
 */
struct BlastUngappedCutoffs {
    Int4 x_dropoff_init; /**< Raw X-dropoff value specified by the
                              bit score in BlastInitialWordOptions */
    Int4 x_dropoff;   /**< Raw X-dropoff value used in the ungapped extension */
    Int4 cutoff_score; /**< Cutoff score for saving ungapped hits. */
    Int4 reduced_nucl_cutoff_score; /**< for blastn, a reduced cutoff score
                                      for use with approximate ungapped
                                      alignments */
};

/** Parameter block that contains a pointer to BlastInitialWordOptions
 * and the values derived from it.
 */
struct BlastInitialWordParameters {
   const BlastInitialWordOptions* options; /**< The original (unparsed) options. */

   Int4 x_dropoff_max; /**< largest X-drop cutoff across all contexts */
   Int4 cutoff_score_min; /**< smallest cutoff score across all contexts */
   BlastUngappedCutoffs *cutoffs;   /**< cutoff values (one per context) */

   Int4 nucl_score_table[256]; /**< the combined score of all match/mismatch
                                    combinations for aligning four bases */
   Boolean ungapped_extension; /**< Should an ungapped extension be
                                  performed? */

   BlastInitialWordParameters(const BlastInitialWordOptions* word_options,
                              const BlastHitSavingParameters* hit_params,
                              const BlastScoreBlk* sbp,
                              QueryInfo* query_info,
                              Int8 subject_length);
   Int2 BlastInitialWordParamtersUpdate(const BlastHitSavingParameters* hit_params,
                                        const BlastScoreBlk* sbp,
                                        QueryInfo* query_info,
                                        Int8 subj_length);
   ~BlastInitialWordParameters();
};

/** Parameters for setting up effective lengths and search spaces.
 * The real database size values to be used for statistical calculations, if
 * there are no overriding values in options.
 */
struct BlastEffectiveLengthsParameters {
   BlastEffectiveLengthsOptions* options; /**< User provided values for these
                                             parameters */
   Int8 real_db_length; /**< Total database length to use in search space
                           calculations. */
   Int4 real_num_seqs;  /**< Number of subject sequences to use for search
                           space calculations */

   BlastEffectiveLengthsParameters(BlastEffectiveLengthsOptions* eff_options,
		   Int8 db_length, Int4 num_seqs) :
		   options(eff_options),
		   real_db_length(db_length),
		   real_num_seqs(num_seqs) {}
   ~BlastEffectiveLengthsParameters() {};
};

/** Scoring parameters block
 *  Contains scoring-related information that is actually used
 *  for the blast search
 */
struct BlastScoringParameters {
   BlastScoringOptions *options; /**< User-provided values for these params */
   Int2 reward;      /**< Reward for a match */
   Int2 penalty;     /**< Penalty for a mismatch */
   Int4 gap_open;    /**< Extra penalty for starting a gap (scaled version) */
   Int4 gap_extend;  /**< Penalty for each gap residue  (scaled version) */
   Int4 shift_pen;   /**< Penalty for shifting a frame in out-of-frame
                        gapping (scaled version) */
   double scale_factor; /**< multiplier for all cutoff scores */

   BlastScoringParameters(BlastScoringOptions* score_options,
		   BlastScoreBlk* sbp)
   {
       ASSERT(score_options != NULL);
       options = score_options;
       scale_factor = sbp->scale_factor;
       reward = sbp->reward;
       penalty = sbp->penalty;
       gap_open = score_options->gap_open * (Int4)scale_factor;
       gap_extend = score_options->gap_extend * (Int4)scale_factor;
   }
};

#endif	/* BLAST_PARAMETERS_H */

