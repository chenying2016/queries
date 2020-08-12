#ifndef __CHAIN_DP_H
#define __CHAIN_DP_H

#include "../corelib/hbn_aux.h"
#include "../corelib/gapped_candidate.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int length;
    idx qoff;
    idx soff;
    int sdir;
    u64 hash;
} ChainSeed;

typedef kvec_t(ChainSeed) vec_chain_seed;

void
validate_mem(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const ChainSeed* cdpsa,
    const int cdpsc);

typedef struct {
    int max_dist_ref;
    int max_dist_qry;
    int max_band_width;
    int max_skip;
    int min_cnt;
    int min_score;
    int max_join_long;
    int max_join_short;
    int min_join_flank_sc;
    double min_join_flank_ratio;
    int dump_info;

    vec_int     f;
    vec_int     p;
    vec_int     t;
    vec_int     v;
    vec_int_pair u;
    vec_chain_seed seeds;
    vec_chain_seed fwd_seeds;
    vec_chain_seed rev_seeds;
} ChainWorkData;

ChainWorkData*
ChainWorkDataNew(int min_seed_cnt, int min_can_score);

ChainWorkData*
ChainWorkDataFree(ChainWorkData* data);

int
find_best_kmer_match(ChainWorkData* data,
    int* best_kmer_match_index,
    int* best_kmer_match_score);

int chaining_find_candidates(ChainWorkData* data,
        ChainSeed* chain_seed_array,
        int chain_seed_count,
        const int is_maximal_exact_match,
        const int subject_strand,
        vec_init_hit* init_hit_list,
        vec_chain_seed* chain_seed_list);

void
join_chains(ChainWorkData* data, 
    HbnInitHit* hit_array,
    int *hit_count,
    ChainSeed* chain_seed_array);

void
join_adjacent_chains(ChainWorkData* data, 
    vec_init_hit* init_hit_list, 
    vec_chain_seed* chain_seed_list);

#ifdef __cplusplus
}
#endif

#endif // __CHAIN_DP_H