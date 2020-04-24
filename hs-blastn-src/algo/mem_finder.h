#ifndef __MEM_FINDER_H
#define __MEM_FINDER_H

#include "../corelib/gapped_candidate.h"
#include "../corelib/seqdb.h"
#include "hbn_lookup_table.h"
#include "chain_dp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    u64 hash;
    int offset;
    int occ;
    int strand;
} KmerInfo;

typedef kvec_t(KmerInfo) vec_kmif;

typedef struct {
    u8* pv;
    vec_kmif kmif_list;
    u8* cv;
    int* seed_count;
    int* kmer_stats;
    int kmer_size;
} SmallLookupTable;

typedef struct {
    vec_kmif qry_kmif_list;
    vec_kmif ref_kmif_list;
    vec_chain_seed fwd_chain_seed_list;
    vec_init_hit fwd_init_hit_list;
    vec_chain_seed rev_chain_seed_list;
    vec_init_hit rev_init_hit_list;
    vec_init_hit init_hit_list;
    const u8* fwd_ref;
    const u8* rev_ref;
    SmallLookupTable* slt;
    ChainWorkData* chain_data;
    int ref_size;
    const char* ref_name;
    int kmer_size;
    int window_size;
    int mem_size;
} MaximalExactMatchWorkData;

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataNew(int kmer_size, int window_size, int mem_size, int chain_score);

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataFree(MaximalExactMatchWorkData* data);

int
MaximalExactMatchWorkData_FindCandidates(
    MaximalExactMatchWorkData* data, 
    const u8* read, 
    const int read_size,
    const int read_strand);

void 
MaximalExactMatchWorkData_Init(MaximalExactMatchWorkData* data, 
    const u8* fwd_ref,
    const u8* rev_ref, 
    const int ref_size);

void
validate_mem(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const ChainSeed* cdpsa,
    const int cdpsc);

    

#ifdef __cplusplus
}
#endif

#endif // __MEM_FINDER_H