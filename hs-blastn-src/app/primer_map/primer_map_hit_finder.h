#ifndef __PRIMER_MAP_HIT_FINDER_H
#define __PRIMER_MAP_HIT_FINDER_H

#include "primer_map_chain_dp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    u32 hash;
    int offset;
    int context;
    u64 key;
} PrimerMapWord;

typedef kvec_t(PrimerMapWord) vec_pm_word;

typedef struct {
    int query_context;
    int query_offset;
    int primer_context;
    int primer_offset;
    int length;
    u64 key;
} PrimerMapSeed;

typedef kvec_t(PrimerMapSeed) vec_pm_seed;

typedef struct {
    int seq_index;
    int context;
    const char* name;
    int size;
    const u8* sequence;
} PrimerMapContextInfo;

typedef kvec_t(PrimerMapContextInfo) vec_pm_ctxi;

typedef struct {
    int word_size;
    int word_stride;
    ChainWorkData* chain;
    vec_pm_seed seed_list;
    int num_query_contexts;
    vec_pm_ctxi query_context_list;
    vec_pm_word query_word_list;
    int num_primer_contexts;
    vec_pm_ctxi primer_context_list;
    vec_pm_word primer_word_list;
    vec_init_hit hit_list;
    vec_chain_seed hit_seed_list;
} PrimerMapHitFindData;

PrimerMapHitFindData*
PrimerMapHitFindDataNew(int word_size, int word_stride, int chain_score);

PrimerMapHitFindData*
PrimerMapHitFindDataFree(PrimerMapHitFindData* data);

const u8*
PrimerMapHitFindData_ExtractPrimer(PrimerMapHitFindData* data, int context);

int
PrimerMapHitFindData_PrimerSize(PrimerMapHitFindData* data, int context);

const char*
PrimerMapHitFindData_PrimerName(PrimerMapHitFindData* data, int context);

void
PrimerMapHitFindData_CleanPrimerData(PrimerMapHitFindData* data);

void
PrimerMapHitFindData_AddOnePrimer(PrimerMapHitFindData* data,
    const int primer_index,
    const u8* primer,
    const int primer_size,
    const char* primer_name);

const u8*
PrimerMapHitFindData_ExtractQuery(PrimerMapHitFindData* data, int context);

int
PrimerMapHitFindData_QuerySize(PrimerMapHitFindData* data, int context);

const char*
PrimerMapHitFindData_QueryName(PrimerMapHitFindData* data, int context);

void
PrimerMapHitFindData_CleanQueryData(PrimerMapHitFindData* data);

void
PrimerMapHitFindData_AddOneQuery(PrimerMapHitFindData* data,
    const int query_index,
    const u8* query,
    const int query_size,
    const char* query_name);

void
PrimerMapHitFindData_BuildPrimerWordList(PrimerMapHitFindData* data);

void
PrimerMapHitFindData_BuildQueryWordList(PrimerMapHitFindData* data);

void
PrimerMapHitFindData_FindHits(PrimerMapHitFindData* data);

#ifdef __cplusplus
}
#endif

#endif // __PRIMER_MAP_HIT_FINDER_H