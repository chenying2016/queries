#ifndef __PRIMER_MAP_CHAIN_DP_H
#define __PRIMER_MAP_CHAIN_DP_H

#include "../../algo/chain_dp.h"

#ifdef __cplusplus
extern "C" {
#endif

ChainWorkData*
ChainWorkDataNew_PrimerMap(int min_seed_cnt, int min_can_score);

int chaining_find_candidates_primer_map(ChainWorkData* data,
        ChainSeed* chain_seed_array,
        int chain_seed_count,
        const int is_maximal_exact_match,
        const int subject_strand,
        vec_init_hit* init_hit_list,
        vec_chain_seed* chain_seed_list);

#ifdef __cplusplus
}
#endif

#endif // __PRIMER_MAP_CHAIN_DP_H