#ifndef __INIT_HIT_FINDER_H
#define __INIT_HIT_FINDER_H

#include "../corelib/hbn_aux.h"
#include "chain_dp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    u32 hash;
    int offset;
    int strand;
} InitHitFindWord;

typedef kvec_t(InitHitFindWord) vec_ihf_word;

typedef struct {
    int query_offset;
    int subject_offset;
} InitHitFindSeed;

typedef kvec_t(InitHitFindSeed) vec_ihf_seed;

typedef struct {
    int word_size;
    int word_stride;

    int query_oid;
    const char* query_name;
    const u8* fwd_query;
    const u8* rev_query;
    int query_size;

    int subject_oid;
    const char* subject_name;
    const u8* fwd_subject;
    const u8* rev_subject;
    int subject_size;

    vec_ihf_word query_word_list;
    vec_ihf_word subject_word_list;
    vec_ihf_seed fwd_subject_seed_list;
    vec_ihf_seed rev_subject_seed_list;

    ChainWorkData* chain;
    vec_init_hit hit_list;
    vec_chain_seed hit_seed_list;
} InitHitFindData;

InitHitFindData*
InitHitFindDataNew(int word_size, int word_stride, int chain_score);

InitHitFindData*
InitHitFindDataFree(InitHitFindData* data);

void
InitHitFindData_SetupCnsAlignInfoFromHit(InitHitFindData* data,
    HbnInitHit* hit,
    vec_chain_seed* seed_list);

void
InitHitFindData_SetupMapAlignInfoFromHit(InitHitFindData* data,
    HbnInitHit* hit,
    vec_chain_seed* seed_list);

void
InitHitFindData_FindHits(InitHitFindData* data, int subject_dir);

void
InitHitFindData_Init(InitHitFindData* data,
    const int subject_oid,
    const char* subject_name,
    const u8* fwd_subject,
    const u8* rev_subject,
    const int subject_size);

void
InitHitFindData_AddQuery(InitHitFindData* data,
    const int query_oid,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size);

void
print_init_hit_range(HbnInitHit* hit);

#ifdef __cplusplus
}
#endif

#endif // __INIT_HIT_FINDER_H