#ifndef __SEARCH_SETUP_H
#define __SEARCH_SETUP_H

#include "../../ncbi_blast/setup/blast_parameters.h"
#include "../../ncbi_blast/setup/blast_stat.h"
#include "../../corelib/seqdb.h"
#include "hbn_options.h"
#include "hbn_options_handle.h"

#ifdef __cplusplus
extern "C" {
#endif

BlastScoreBlk*
CSetupFactory__CreateScoreBlock(const HbnOptionsHandle* opts_memento,
                                BLAST_SequenceBlk* queries,
                                BlastQueryInfo* query_info);

Int2 
BLAST_GapAlignSetUp(EBlastProgramType program_number,
    const text_t* db,
    const BlastScoringOptions* scoring_options,
    const BlastEffectiveLengthsOptions* eff_len_options,
    const BlastExtensionOptions* ext_options,
    const BlastHitSavingOptions* hit_options,
    const BlastInitialWordOptions* word_options,
    BlastQueryInfo* query_info, 
    BlastScoreBlk* sbp, 
    BlastScoringParameters** score_params,
    BlastExtensionParameters** ext_params,
    BlastHitSavingParameters** hit_params,
    BlastEffectiveLengthsParameters** eff_len_params,
    BlastInitialWordParameters** word_params);                                

#ifdef __cplusplus
}
#endif

#endif // __SEARCH_SETUP_H