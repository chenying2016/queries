#ifndef __TRACEBACK_STAGE_H
#define __TRACEBACK_STAGE_H

#include "search_setup.h"
#include "../../ncbi_blast/setup/blast_hits.h"

#ifdef __cplusplus
extern "C" {
#endif

void
add_align_string(BlastHSP* hsp, const u8* query, const u8* subject, kstring_t* aligned_string);

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
    kstring_t* aligned_string);

#ifdef __cplusplus
}
#endif

#endif