#ifndef __BACKUP_RESULTS_H
#define __BACKUP_RESULTS_H

#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/setup/blast_hits.h"
#include "hbn_options.h"
#include "hbn_results.h"
#include "tabular_format.h"
#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

void
dump_m4_hits(const text_t* query_vol,
    const text_t* subject_vol,
    HbnHSPResults* results,
    const HbnProgramOptions* opts);

void
dump_one_result_set(const CSeqDB* queries,
    const CSeqDB* db,
    HbnHSPResults* results, 
    const HbnProgramOptions* opts,
    FILE* out, 
    FILE* backup_out, 
    pthread_mutex_t* out_lock);

void
recover_qi_vs_sj_results(const CSeqDB* queries, const CSeqDB* db, const HbnProgramOptions* opts, FILE* in, FILE* out);

#ifdef __cplusplus
}
#endif

#endif // __BACKUP_RESULTS_H