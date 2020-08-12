#ifndef __HBN_EXTEND_SUBSEQ_HIT_H
#define __HBN_EXTEND_SUBSEQ_HIT_H

#include "hbn_options.h"
#include "subseq_hit.h"
#include "../../algo/hbn_traceback.h"
#include "../../algo/init_hit_finder.h"
#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/setup/blast_hits.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    vec_subseq_hit fwd_sbjct_subseq_list;
    vec_subseq_hit rev_sbjct_subseq_list;
    vec_subseq_hit sbjct_subseq_list;
    InitHitFindData* hit_finder;
    HbnTracebackData* traceback_data;
    vec_chain_seed chain_seed_list;
} HbnSubseqHitExtnData;

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataNew(const HbnProgramOptions* opts);

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataFree(HbnSubseqHitExtnData* data);

void
hbn_extend_query_subseq_hit_list(HbnSubseqHit* subseq_hit_array,
    int subseq_hit_count,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_id,
    const int fwd_query_context,
    const int query_length,
    vec_u8* subject_v,
    const text_t* db,
    const HbnProgramOptions* opts,
    HbnSubseqHitExtnData* data,
    BlastHitList* hit_list,
    HbnHSPResults* results);

#ifdef __cplusplus
}
#endif

#endif // __HBN_EXTEND_SUBSEQ_HIT_H