#ifndef __SUBSEQ_HIT_H
#define __SUBSEQ_HIT_H

#include "../../corelib/kvec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int qdir;
    int qoff;
    int sid;
    int soff;
    size_t sfrom;
    size_t sto;
    int score;
} HbnSubseqHit;

#define dump_subseq_hit(output_func, stream, hit) \
    output_func(stream, "qid = %d, qdir = %d, qoff = %d, sid = %d, soff = %d, sfrom = %zu, sto = %zu, score = %d\n", \
        (hit).qid, \
        (hit).qdir, \
        (hit).qoff, \
        (hit).sid, \
        (hit).soff, \
        (hit).sfrom, \
        (hit).sto, \
        (hit).score)

typedef kvec_t(HbnSubseqHit) vec_subseq_hit;

void sort_subseq_hit_score_gt(size_t n, HbnSubseqHit* a);

void sort_subseq_hit_sid_lt(size_t n, HbnSubseqHit* a);

void sort_subseq_hit_sfrom_lt(size_t n, HbnSubseqHit* a);

#ifdef __cplusplus
}
#endif

#endif // __SUBSEQ_HIT_H