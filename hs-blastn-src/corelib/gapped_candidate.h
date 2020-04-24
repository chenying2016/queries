#ifndef __GAPPED_CANDIDATE_H
#define __GAPPED_CANDIDATE_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int qdir;
    int qoff;
    int qbeg;
    int qend;
    int qsize;
    int sid;
    int sdir;
    size_t soff;
    size_t sbeg;
    size_t send;
    size_t ssize;
    int score;
    void* chain_seed_array;
    int chain_seed_offset;
    int chain_seed_count;
} HbnInitHit;

typedef kvec_t(HbnInitHit) vec_init_hit;

#define dump_init_hit(output_func, stream, hit) \
    output_func(stream, "%d\t%d\t%d\t%d\t%d\t%d\t%zu\t%zu\t%d\n", \
        (hit).qid, \
        (hit).qdir, \
        (hit).qoff, \
        (hit).qsize, \
        (hit).sid, \
        (hit).sdir, \
        (hit).soff, \
        (hit).ssize, \
        (hit).score)

void ks_introsort_init_hit_score_gt(size_t n, HbnInitHit* a);
void ks_introsort_init_hit_qid_lt(size_t n, HbnInitHit* a);
void ks_introsort_init_hit_sid_lt(size_t n, HbnInitHit* a);

typedef struct {
    int qid;
    int qdir;
    int sid;
    int sdir;
    int score;
} HbnConsensusInitHit;

typedef kvec_t(HbnConsensusInitHit) vec_cns_hit;

#define dump_cns_hit(output_func, stream, cns_hit) \
    output_func(stream, "%d\t%d\t%d\t%d\t%d\n", \
        (cns_hit).qid, \
        (cns_hit).qdir, \
        (cns_hit).sid, \
        (cns_hit).sdir, \
        (cns_hit).score)

void ks_introsort_cns_hit_sid_lt(size_t n, HbnConsensusInitHit* a);
void ks_introsort_cns_hit_score_gt(size_t n, HbnConsensusInitHit* a);

#ifdef __cplusplus
}
#endif

#endif // __GAPPED_CANDIDATE_H